#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <complex.h>
#include <math.h>

#include <omp.h>

#include "fastsphere.h"
#include "config.h"
#include "init.h"
#include "util.h"
#include "scatmat.h"
#include "farfield.h"
#include "spreflect.h"
#include "ptsrc.h"

#ifndef _OPENMP
int omp_get_max_threads () { return 1; }
#endif /* _OPENMP */

void usage () {
	fprintf (stderr, "USAGE: fastsphere [-h] [-n bounces] [-m dist] [-i input] [-o output] [-t samples]\n");
	fprintf (stderr, "\t-h\tPrint this message and exit\n");
	fprintf (stderr, "\t-i\tSpecify configuration file path (default: stdin)\n");
	fprintf (stderr, "\t-o\tSpecify output radiation file path (default: stdout)\n");
	fprintf (stderr, "\t-m\tSpecify measurement distance (default: infinite)\n");
	fprintf (stderr, "\t-n\tSpecify maximum number of bounces (default: 10)\n");
	fprintf (stderr, "\t-t\tSpecify number of theta samples in scattering pattern (default: optimized)\n");

	exit (EXIT_FAILURE);
}

int writebvec (FILE *out, complex double *vec, int n, int m) {
	int size[2], len;

	size[0] = n; size[1] = m;
	len = n * m;

	/* Write the size of the vector. */
	fwrite (size, sizeof(int), 2, out);
	/* Write the vector itself. */
	fwrite (vec, sizeof(complex double), len, out);

	return len;
}

int main (int argc, char **argv) {
	int nspheres, nsptype, n, nterm, i, j, nit, nbounce = 0;
	int ntbg, ntheta = 0, mxthd = 1, ndig;
	sptype *sparms, bgspt;
	spscat *slist;
	bgtype bg;
	exctparm exct;
	shdata shtr, shroot, shinterp;
	itconf itc;
	trdesc *trans;
	double mslen = -1, err;

	complex double *rhs, *sol, *radpat, *finc,
		*oinc = NULL, *lastwv = NULL, *rpinterp = NULL;

	FILE *fptr = NULL;
	char *inname = NULL, *outname = NULL, ch;

	while ((ch = getopt (argc, argv, "hi:o:m:n:t:")) != -1) {
		switch (ch) {
		case 'i':
			inname = optarg;
			break;
		case 'o':
			outname = optarg;
			break;
		case 'm':
			mslen = atof (optarg);
			break;
		case 'n':
			nbounce = atoi (optarg);
			break;
		case 't':
			ntheta = atoi (optarg);
			break;
		case 'h': default:
			usage ();
		}
	}

	if (!inname) fptr = stdin;
	else fptr = critopen (inname, "r");

	mxthd = omp_get_max_threads ();

	/* Only initialize FFTW threads if OpenMP is used. */
	if (mxthd > 1) fftw_init_threads ();

	if (nbounce > 0)
		readcfg (fptr, &nspheres, &nsptype, &sparms, 
				&bgspt, &slist, &bg, &exct, &itc, &ndig);
	else {
		readcfg (fptr, &nspheres, &nsptype, &sparms, 
				NULL, &slist, &bg, &exct, &itc, &ndig);
		/* There is no enclosing sphere; the k-vector and density should
		 * not contrast from background material. */
		bgspt.k = bg.k;
		bgspt.rho = 1.0;
	}
	fprintf (stderr, "Parsed configuration for %d spheres at %g MHz\n", nspheres, exct.f / 1e6);
	sphinit (sparms, nsptype, bgspt.k, bgspt.rho, &shtr, -1, ndig);
	fprintf (stderr, "Initialized spherical harmonic data for degree %d (%d digits accurate)\n", shtr.deg, ndig);
	trans = sphbldfmm (slist, nspheres, bgspt.k, &shtr);
	fprintf (stderr, "Built FMM translators for all spheres\n");

	/* Multithreaded FFT for the root-level transform, since that is
	 * a serial point in the code. All other FFTs are serialized, because
	 * parallelization occurs at the sphere level. */
	if (mxthd > 1) {
		fftw_plan_with_nthreads (mxthd);
		fprintf (stderr, "Using %d threads for root-level FFT\n", mxthd);
	}

	if (nbounce > 0) {
		/* Use the lowest-possible sampling rate for the enclosing sphere. */
		esbdinit (&bgspt, bg.k, 1.0, &shroot, 0, ndig);
		fprintf (stderr, "Built data for enclosing sphere\n");
	} else {
		/* Use the greater of the estimated root bandwidth or that requested. */
		i = rootorder (slist, nspheres, bg.k, ndig);
		ntheta = MAX (i + (i % 2) + 1, ntheta);
		fshtinit (&shroot, i, ntheta, 2 * ntheta, 1);
		fprintf (stderr, "Built spherical harmonic data for far-field\n");
	}

	nterm = shtr.ntheta * shtr.nphi;
	n = nspheres * nterm;
	rhs = calloc (2 * n, sizeof(complex double));
	sol = rhs + n;
	
	ntbg = shroot.ntheta * shroot.nphi;

	radpat = calloc (2 * ntbg, sizeof(complex double));
	finc = radpat + ntbg;

	if (nbounce > 0) {
		oinc = calloc (2 * ntbg, sizeof(complex double));
		lastwv = oinc + ntbg;
	}

	/* Build the RHS contribution from point sources. */
#pragma omp parallel for private(i) default(shared)
	for (i = 0; i < exct.nps; ++i)
		ptsrcexp (finc, bg.k, &shroot, exct.psmag[i],
				exct.psloc + 3 * i, exct.psax + 3 * i, exct.alpha[i]);

	if (exct.npw || nbounce > 0) {
		/* Convert to a harmonic expansion. */
		ffsht (finc, &shroot, 0);
		shscale (finc, &shroot, 1); 
		
		/* Include contribution of incident plane waves. */
#pragma omp parallel for private(i) default(shared)
		for (i = 0; i < exct.npw; ++i)
			shincident (shroot.deg, finc, shroot.nphi,
					exct.pwmag[i], exct.theta[i], exct.phi[i]); 
		
		/* Include the appropriate scaling factor. */
		shscale (finc, &shroot, -1); 
		
		if (nbounce > 0) {
			spreflect (finc, finc, bgspt.transmit, shroot.deg, shroot.nphi, 0, 1);
			memcpy (oinc, finc, ntbg * sizeof(complex double));
		}
		
		ifsht (finc, &shroot, 0);
	}

	fartonear (rhs, finc, slist, nspheres, bgspt.k, &shtr, &shroot);
	sprflpw (rhs, slist, nspheres, &shtr);

	fprintf (stderr, "Built RHS vector\n");

	/* The RHS is the initial guess. */
	memcpy (sol, rhs, n * sizeof(complex double));

	for (j = 0, nit = 1; j < itc.restart && nit > 0; ++j)
		nit = gmres (sol, rhs, j, slist, nspheres, trans, &shtr, &itc);

	neartofar (radpat, sol, slist, nspheres, bgspt.k, &shroot, &shtr);

	if (nbounce > 0) {
		/* Reset the incident field. */
		memcpy (finc, oinc, ntbg * sizeof(complex double));
		/* Add in the reflection of the inner spheres. */
		ffsht (radpat, &shroot, 0);
		spreflect (finc, radpat, bgspt.reflect, shroot.deg, shroot.nphi, 1, -1);
		/* Store a copy of this field for comparison in the next cycle. */
		memcpy (lastwv, finc, ntbg * sizeof(complex double));
	}

	for (i = 0, --nbounce; i < nbounce; ++i) {
		/* Distribute standing wave to internal spheres. */
		ifsht (finc, &shroot, 0);
		fartonear (rhs, finc, slist, nspheres, bgspt.k, &shtr, &shroot);

		/* Reflect incoming wave from internal spheres. */
		sprflpw (rhs, slist, nspheres, &shtr);

		/* Solve for the fields scattered by inner spheres. */
		for (j = 0, nit = 1; j < itc.restart && nit > 0; ++j)
			nit = gmres (sol, rhs, j, slist, nspheres, trans, &shtr, &itc);

		/* Compute the far-field pattern of the internal spheres. */
		neartofar (radpat, sol, slist, nspheres, bgspt.k, &shroot, &shtr);

		/* Reset the standing-wave pattern to the incident field. */
		memcpy (finc, oinc, ntbg * sizeof(complex double));

		/* Reflect the scattered field into the standing wave. */
		ffsht (radpat, &shroot, 0);
		spreflect (finc, radpat, bgspt.reflect, shroot.deg, shroot.nphi, 1, -1);

		/* Find the relative change in the standing-wave pattern. */
		err = rmserror (finc, lastwv, ntbg);
		memcpy (lastwv, finc, ntbg * sizeof(complex double));

		fprintf (stderr, "Iteration %d complete (error: %g)\n", i + 1, err);
		fflush (stdout);

		if (err < itc.eps) break;
	}

	if (nbounce > 0) {
		/* Transmit the scattered field through the outer boundary. */
		spreflect (radpat, radpat, bgspt.transmit + shroot.deg, shroot.deg, shroot.nphi, 0, 1);
		
		/* Now add in the reflected standing-wave coefficients. */
		spreflect (radpat, finc, bgspt.reflect + shroot.deg, shroot.deg, shroot.nphi, 1, 1);

		/* Back to a far-field signature if the harmonic
		 * representation won't be needed again. */
		if (mslen <= 0 && ntheta <= shroot.ntheta)
			ifsht (radpat, &shroot, 0);
	}

	/* Measure at a finite distance, if desired. */
	if (mslen > 0) {
		/* Scale the measurement distance to wavelengths. */
		mslen *= exct.f / bg.cabs;
		fprintf (stderr, "Computing scattered field at radius %g wavelengths\n", mslen);

		/* Forward transform if necessary. */
		if (nbounce <= 0) ffsht (radpat, &shroot, 0);

		/* Scale the SH coefficients appropriately. */
		shscale (radpat, &shroot, 1);

		/* Now include the radial factor in the SH coefficients. */
		shradial (shroot.deg, radpat, shroot.nphi, bg.k, mslen);

		/* To angular representation unless interpolation is needed. */
		if (nbounce <= 0 || ntheta <= shroot.ntheta)
			ifsht (radpat, &shroot, 0);
	}

	/* Interpolate the far-field pattern, if the root transform
	 * contains fewer points than those requested. In such cases,
	 * the field representation is already harmonic. */
	if (nbounce > 0 && ntheta > shroot.ntheta) {
		fprintf (stderr, "Interpolating for bandwidth %d\n", ntheta);
		/* Initialize the spherical harmonic transform for interpolation. */
		fshtinit (&shinterp, shroot.deg, ntheta, 2 * ntheta, 0);

		/* Allocate the expanded far-field pattern. */
		rpinterp = calloc (shinterp.ntheta * shinterp.nphi, sizeof(complex double));

		/* Copy the harmonic expansion for interpolation. */
		copysh (shroot.deg, rpinterp, shinterp.nphi, radpat, shroot.nphi);

		/* Perform the inverse harmonic transform to do the interpolation. */
		ifsht (rpinterp, &shinterp, 0);
	}

	if (!outname) fptr = stdout;
	else fptr = critopen (outname, "w");
	if (rpinterp) {
		fprintf (stderr, "Writing solution for %d theta and %d phi samples\n", shinterp.ntheta, shinterp.nphi);
		writebvec (fptr, rpinterp, shinterp.nphi, shinterp.ntheta);
	} else {
		fprintf (stderr, "Writing solution for %d theta and %d phi samples\n", shroot.ntheta, shroot.nphi);
		writebvec (fptr, radpat, shroot.nphi, shroot.ntheta);
	}
	fclose (fptr);

	clrspheres (sparms, nsptype);
	sphclrfmm (trans, nspheres * nspheres);
	fshtfree (&shroot);
	fshtfree (&shtr);

	free (trans);
	if (exct.pwmag) free (exct.pwmag);
	if (exct.theta) free (exct.theta);
	if (exct.psmag) free (exct.psmag);
	if (exct.psloc) free (exct.psloc);

	free (rhs);
	free (radpat);

	if (rpinterp) {
		fshtfree (&shinterp);
		free (rpinterp);
	}

	if (nbounce > 0) free (oinc);

	if (mxthd > 1) fftw_cleanup_threads ();

	fftw_cleanup ();

	return EXIT_SUCCESS;
}

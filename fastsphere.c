#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <complex.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif /* _OPENMP */

#include "fastsphere.h"
#include "config.h"
#include "init.h"
#include "util.h"
#include "scatmat.h"
#include "farfield.h"
#include "spreflect.h"

void usage () {
	fprintf (stderr, "USAGE: fastsphere [-h] [-n bounces] [-m dist] [-i input] [-o output] [-r rhsfile]\n");
	fprintf (stderr, "\t-h\tPrint this message and exit\n");
	fprintf (stderr, "\t-i\tSpecify configuration file path (default: stdin)\n");
	fprintf (stderr, "\t-o\tSpecify output radiation file path (default: stdout)\n");
	fprintf (stderr, "\t-r\tSpecify output RHS file path (default: none)\n");
	fprintf (stderr, "\t-m\tSpecify measurement distance (default: infinite)\n");
	fprintf (stderr, "\t-n\tSpecify maximum number of bounces (default: 10)\n");

	exit (EXIT_FAILURE);
}

double rmserror (complex double *v, complex double *r, int n) {
	double err = 0, errd = 0, e;
	int i;

	for (i = 0; i < n; ++i) {
		e = cabs (v[i] - r[i]);
		err += e * e;
		e = cabs (r[i]);
		errd += e * e;
	}

	return sqrt (err / errd);
}

FILE *critopen (char *fname, char *mode) {
	FILE *fptr;

	fptr = fopen (fname, mode);

	if (!fptr) {
		fprintf (stderr, "ERROR: Could not open input file %s.\n", fname);
		exit (EXIT_FAILURE);

		return NULL;
	}

	return fptr;
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
	int nspheres, nsptype, n, nterm, i, nbounce = 0, ntbg;
	sptype *sparms, bgspt;
	spscat *slist;
	bgtype bg;
	exctparm exct;
	shdata shtr, shroot;
	itconf itc;
	trdesc *trans;
	double mslen = -1, err;

	complex double *rhs, *sol, *radpat, *finc,
		*oinc = NULL, *lastwv = NULL;

	FILE *fptr = NULL;
	char *inname = NULL, *outname = NULL, *rhsname = NULL, ch;

	while ((ch = getopt (argc, argv, "hi:o:r:m:n:")) != -1) {
		switch (ch) {
		case 'i':
			inname = optarg;
			break;
		case 'o':
			outname = optarg;
			break;
		case 'r':
			rhsname = optarg;
			break;
		case 'm':
			mslen = atof (optarg);
			break;
		case 'n':
			nbounce = atoi (optarg);
			break;
		case 'h': default:
			usage ();
		}
	}

	if (!inname) fptr = stdin;
	else fptr = critopen (inname, "r");

#ifdef _OPENMP
	/* Only initialize FFTW threads if OpenMP is used. */
	fftw_init_threads ();
#endif

	if (nbounce > 0)
		readcfg (fptr, &nspheres, &nsptype, &sparms, &bgspt, &slist, &bg, &exct, &itc);
	else {
		readcfg (fptr, &nspheres, &nsptype, &sparms, NULL, &slist, &bg, &exct, &itc);
		/* There is no enclosing sphere; the k-vector and density should
		 * not contrast from background material. */
		bgspt.k = bg.k;
		bgspt.rho = 1.0;
	}
	fprintf (stderr, "Parsed configuration for %d spheres at %g MHz\n", nspheres, exct.f / 1e6);
	sphinit (sparms, nsptype, bgspt.k, bgspt.rho, &shtr, -1);
	fprintf (stderr, "Initialized spherical harmonic data for degree %d\n", shtr.deg);
	trans = sphbldfmm (slist, nspheres, bgspt.k, &shtr);
	fprintf (stderr, "Built FMM translators for all spheres\n");

#ifdef _OPENMP
	/* Multithreaded FFT for the root-level transform, since that is
	 * a serial point in the code. All other FFTs are serialized, because
	 * parallelization occurs at the sphere level. */
#pragma omp parallel default(shared)
	if (omp_get_thread_num () == 0) n = omp_get_num_threads ();
	fftw_plan_with_nthreads (n);
	fprintf (stderr, "Using %d threads for root-level FFT\n", n);
#endif /* _OPENMP */

	if (nbounce > 0) {
		esbdinit (&bgspt, bg.k, 1.0, &shroot);
		fprintf (stderr, "Built data for enclosing sphere\n");
	} else {
		i = rootorder (slist, nspheres, bg.k);
		n = 2 * i - 1; /* The number of angular samples (in each direction). */
		fshtinit (&shroot, i, n, 2 * n);
		fprintf (stderr, "Built spherical harmonic data for far-field\n");
	}

	nterm = shtr.ntheta * shtr.nphi + 2;
	n = nspheres * nterm;
	rhs = calloc (2 * n, sizeof(complex double));
	sol = rhs + n;
	
	ntbg = shroot.ntheta * shroot.nphi + 2;

	radpat = calloc (2 * ntbg, sizeof(complex double));
	finc = radpat + ntbg;

	if (nbounce > 0) {
		oinc = calloc (2 * ntbg, sizeof(complex double));
		lastwv = oinc + ntbg;
	}

	/* Build the RHS vector consisting of multiple plane waves. */
	for (i = 0; i < exct.npw; ++i)
		shincident (shroot.deg, finc, shroot.nphi,
				(exct.mag)[i], (exct.theta)[i], (exct.phi)[i]);

	/* Include the appropriate scaling factor. */
	shscale (finc, &shroot, -1);

	if (nbounce > 0) {
		spreflect (finc, finc, bgspt.transmit, shroot.deg, shroot.nphi, 0, 1);
		memcpy (oinc, finc, ntbg * sizeof(complex double));
	}

	ifsht (finc, &shroot);
	fartonear (rhs, finc, slist, nspheres, bgspt.k, &shtr, &shroot);
	sprflpw (rhs, slist, nspheres, &shtr);

	fprintf (stderr, "Built RHS vector\n");

	if (rhsname) {
		fptr = critopen (rhsname, "w");
		writebvec (fptr, rhs, n, 1);
		fclose (fptr);
	}

	/* The RHS is the initial guess. */
	memcpy (sol, rhs, n * sizeof(complex double));

	itsolve (sol, rhs, slist, nspheres, trans, &shtr, &itc);

	neartofar (radpat, sol, slist, nspheres, bgspt.k, &shroot, &shtr);

	if (nbounce > 0) {
		/* Reset the incident field. */
		memcpy (finc, oinc, ntbg * sizeof(complex double));
		/* Add in the reflection of the inner spheres. */
		ffsht (radpat, &shroot);
		spreflect (finc, radpat, bgspt.reflect, shroot.deg, shroot.nphi, 1, -1);
		/* Store a copy of this field for comparison in the next cycle. */
		memcpy (lastwv, finc, ntbg * sizeof(complex double));
	}

	for (i = 0, --nbounce; i < nbounce; ++i) {
		/* Distribute standing wave to internal spheres. */
		ifsht (finc, &shroot);
		fartonear (rhs, finc, slist, nspheres, bgspt.k, &shtr, &shroot);

		/* Reflect incoming wave from internal spheres. */
		sprflpw (rhs, slist, nspheres, &shtr);

		/* Solve for the fields scattered by inner spheres. */
		itsolve (sol, rhs, slist, nspheres, trans, &shtr, &itc);

		/* Compute the far-field pattern of the internal spheres. */
		neartofar (radpat, sol, slist, nspheres, bgspt.k, &shroot, &shtr);

		/* Reset the standing-wave pattern to the incident field. */
		memcpy (finc, oinc, ntbg * sizeof(complex double));

		/* Reflect the scattered field into the standing wave. */
		ffsht (radpat, &shroot);
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

		/* Back to a far-field signature. */
		ifsht (radpat, &shroot);
	}

	/* Measure at a finite distance, if desired. */
	if (mslen > 0) {
		/* Scale the measurement distance to wavelengths. */
		mslen *= exct.f / bg.cabs;
		fprintf (stderr, "Computing scattered field at radius %g wavelengths\n", mslen);

		ffsht (radpat, &shroot);

		/* Scale the SH coefficients appropriately. */
		shscale (radpat, &shroot, 1);

		/* Now include the radial factor in the SH coefficients. */
		shradial (shroot.deg, radpat, shroot.nphi, bg.k, mslen);

		/* Back to angular representation. */
		ifsht (radpat, &shroot);
	}

	if (!outname) fptr = stdout;
	else fptr = critopen (outname, "w");
	fprintf (stderr, "Writing solution for %d theta and %d phi samples\n", shroot.ntheta, shroot.nphi);
	writebvec (fptr, radpat, shroot.nphi, shroot.ntheta);
	fclose (fptr);

	clrspheres (sparms, nsptype);
	sphclrfmm (trans, nspheres * nspheres);
	fshtfree (&shroot);
	fshtfree (&shtr);

	free (trans);
	free (exct.mag);
	free (exct.theta);

	free (rhs);
	free (radpat);
	if (nbounce > 0) free (oinc);

#ifdef _OPENMP
	fftw_cleanup_threads ();
#endif /* _OPENMP */

	fftw_cleanup ();

	return EXIT_SUCCESS;
}

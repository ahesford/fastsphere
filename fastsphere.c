#include <stdio.h>
#include <stdlib.h>
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
	fprintf (stderr, "USAGE: fastsphere [-h] [-i input] [-o output] [-r rhsfile]\n");
	fprintf (stderr, "\t-h\tPrint this message and exit\n");
	fprintf (stderr, "\t-i\tSpecify configuration file path (default: stdin)\n");
	fprintf (stderr, "\t-o\tSpecify output radiation file path (default: stdout)\n");
	fprintf (stderr, "\t-r\tSpecify output RHS file path (default: none)\n");

	exit (EXIT_FAILURE);
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
	int nspheres, nsptype, n, i, nterm, trunc;
	sptype *sparms, bgspt;
	spscat *slist;
	bgtype bg;
	exctparm exct;
	shdata shtr, shroot;
	itconf itc;
	trdesc *trans;
	complex double *rhs, *sol, *radpat;

	FILE *fptr = NULL;
	char *inname = NULL, *outname = NULL, *rhsname = NULL, ch;

	while ((ch = getopt (argc, argv, "hi:o:r:")) != -1) {
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

	readcfg (fptr, &nspheres, &nsptype, &sparms, &bgspt, &slist, &bg, &exct, &itc);
	fprintf (stderr, "Parsed configuration for %d spheres at %g MHz\n", nspheres, exct.f / 1e6);
	if (bgspt.r < 0) sphinit (sparms, nsptype, bg.k, 1.0, &shtr);
	else sphinit (sparms, nsptype, bgspt.k, bgspt.rho, &shtr);
	fprintf (stderr, "Initialized spherical harmonic data for degree %d\n", shtr.deg);
	if (bgspt.r < 0) trans = sphbldfmm (slist, nspheres, bg.k, &shtr);
	else trans = sphbldfmm (slist, nspheres, bgspt.k, &shtr);
	fprintf (stderr, "Built FMM translators for all spheres\n");

#ifdef _OPENMP
	/* Multithreaded FFT for the root-level transform, since that is
	 * a serial point in the code. All other FFTs are serialized, because
	 * parallelization occurs at the sphere level. */
	fftw_plan_with_nthreads (omp_get_num_threads ());
	fprintf (stderr, "Using %d threads for root-level FFT\n", omp_get_num_threads());
#endif /* _OPENMP */

	/* Set up far-field (or enclosing sphere) transform data. */
	if (bgspt.r < 0) {
		n = rootorder (slist, nspheres, bg.k);
		n = 2 * n - 1; /* The number of angular samples (in each direction). */
		fshtinit (&shroot, shtr.deg, n, 2 * n);
		fprintf (stderr, "Built spherical harmonic data for far-field\n");
	} else {
		esbdinit (&bgspt, bg.k, 1.0, &shroot);
		fprintf (stderr, "Built data for enclosing sphere\n");
	}

	nterm = shtr.ntheta * shtr.nphi;
	n = nspheres * nterm;
	if (bgspt.r >= 0) n += shroot.ntheta * shroot.nphi;
	rhs = calloc (2 * n, sizeof(complex double));
	sol = rhs + n;

	if (bgspt.r < 0) {
		/* Build a translator from the source location to each sphere.
		 * This is the incoming incident field for each sphere. */
		trunc = 2 * shtr.deg - 1; 

#pragma omp parallel private(i) default(shared)
{
		/* These variables are private. */
		double rsrc;
		trdesc trinc;
		trinc.trunc = trunc;
		trinc.type = TRPLANE;

		for (i = 0; i < nspheres; ++i) {
			trinc.sdir[0] = slist[i].cen[0] - exct.cen[0];
			trinc.sdir[1] = slist[i].cen[1] - exct.cen[1];
			trinc.sdir[2] = slist[i].cen[2] - exct.cen[2];
			
			rsrc = sqrt (DVDOT(trinc.sdir,trinc.sdir));
			trinc.sdir[0] /= rsrc;
			trinc.sdir[1] /= rsrc;
			trinc.sdir[2] /= rsrc;
			
			trinc.kr = bg.k * rsrc;
			trinc.trdata = rhs + i * nterm; 
			
			translator (&trinc, shtr.ntheta, shtr.nphi, shtr.theta);
		}
}

		/* Convert the incoming incident field to the reflected incident
		 * field, which is the RHS for this problem. */
		sprflpw (rhs, slist, nspheres, &shtr);
	} else {
		double rsrc;
		trdesc trinc;

		trinc.trunc = 2 * shroot.deg - 1;
		trinc.type = TRPLANE;

		trinc.sdir[0] = -exct.cen[0];
		trinc.sdir[0] = -exct.cen[0];
		trinc.sdir[0] = -exct.cen[0];

		rsrc = sqrt (DVDOT(trinc.sdir,trinc.sdir));
		trinc.sdir[0] /= rsrc;
		trinc.sdir[1] /= rsrc;
		trinc.sdir[2] /= rsrc;

		trinc.kr = bg.k * rsrc;
		trinc.trdata = rhs + nspheres * nterm;

		/* The plane-wave incident field. */
		translator (&trinc, shroot.ntheta, shroot.nphi, shroot.theta);

		/* Take the incident field to SH coefficients. */
		ffsht (trinc.trdata, &shroot);
	}

	fprintf (stderr, "Built RHS vector\n");

	if (rhsname) {
		fptr = critopen (rhsname, "w");
		writebvec (fptr, rhs, nterm, nspheres);
		fclose (fptr);
	}

	itsolve (sol, rhs, slist, nspheres, &bgspt, trans, &shtr, &shroot, &itc);
	fprintf (stderr, "Iteration complete, computing radiation pattern\n");

	fflush (stdout);

	/* Compute the far-field radiation pattern of the object. */
	if (bgspt.r < 0) {
		n = shroot.ntheta * shroot.nphi;
		radpat = malloc (n * sizeof(complex double));
		neartofar (radpat, sol, slist, nspheres, bg.k, &shroot, &shtr);
	} else {
		complex double *sptr;
		n = nspheres * nterm;
		radpat = rhs + n;
		sptr = sol + n;

		/* Compute the far-field pattern of the internal spheres. */
		neartofar (radpat, sol, slist, nspheres, bgspt.k, &shroot, &shtr);
		ffsht (radpat, &shroot);

		/* Transmit this field through the outer boundary. */
		spreflect (radpat, radpat, bgspt.transmit + shroot.deg, shroot.deg, shroot.nphi, 0);

		/* Now add in the reflected standing-wave coefficients. */
		spreflect (radpat, sptr, bgspt.reflect + shroot.deg, shroot.deg, shroot.nphi, 1);

		ifsht (radpat, &shroot);
	}

	if (!outname) fptr = stdout;
	else fptr = critopen (outname, "w");
	fprintf (stderr, "Writing solution for %d theta and %d phi samples\n", shroot.ntheta, shroot.nphi);
	writebvec (fptr, radpat, shroot.nphi, shroot.ntheta);
	fclose (fptr);

	clrspheres (sparms, nsptype);
	sphclrfmm (trans, nspheres * nspheres);
	free (trans);
	free (rhs);
	if (bgspt.r < 0) free (radpat);
	fshtfree (&shroot);
	fshtfree (&shtr);

#ifdef _OPENMP
	fftw_cleanup_threads ();
#endif /* _OPENMP */

	fftw_cleanup ();

	return EXIT_SUCCESS;
}

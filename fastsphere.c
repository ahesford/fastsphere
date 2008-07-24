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
	int nspheres, nsptype, n, nterm;
	sptype *sparms, bgspt;
	spscat *slist;
	bgtype bg;
	exctparm exct;
	shdata shtr, shroot;
	itconf itc;
	trdesc *trans, trinc;
	double rsrc;

	complex double *rhs, *sol, *radpat, *sptr;

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
	sphinit (sparms, nsptype, bgspt.k, bgspt.rho, &shtr);
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

	/* Set up far-field (or enclosing sphere) transform data. */
	esbdinit (&bgspt, bg.k, 1.0, &shroot);
	fprintf (stderr, "Built data for enclosing sphere\n");

	nterm = shtr.ntheta * shtr.nphi;
	n = nspheres * nterm + shroot.ntheta * shroot.nphi;
	rhs = calloc (2 * n, sizeof(complex double));
	sol = rhs + n;

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

	fprintf (stderr, "Built RHS vector\n");

	if (rhsname) {
		fptr = critopen (rhsname, "w");
		writebvec (fptr, rhs, n, 1);
		fclose (fptr);
	}

	itsolve (sol, rhs, slist, nspheres, &bgspt, trans, &shtr, &shroot, &itc);
	fprintf (stderr, "Iteration complete, computing radiation pattern\n");

	fflush (stdout);

	/* Compute the far-field radiation pattern of the object. */
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

	if (!outname) fptr = stdout;
	else fptr = critopen (outname, "w");
	fprintf (stderr, "Writing solution for %d theta and %d phi samples\n", shroot.ntheta, shroot.nphi);
	writebvec (fptr, radpat, shroot.nphi, shroot.ntheta);
	fclose (fptr);

	clrspheres (sparms, nsptype);
	sphclrfmm (trans, nspheres * nspheres);
	free (trans);
	free (rhs);
	fshtfree (&shroot);
	fshtfree (&shtr);

#ifdef _OPENMP
	fftw_cleanup_threads ();
#endif /* _OPENMP */

	fftw_cleanup ();

	return EXIT_SUCCESS;
}

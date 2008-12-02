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
	fprintf (stderr, "USAGE: fastsphere [-h] [-m dist] [-i input] [-o output] [-r rhsfile]\n");
	fprintf (stderr, "\t-h\tPrint this message and exit\n");
	fprintf (stderr, "\t-i\tSpecify configuration file path (default: stdin)\n");
	fprintf (stderr, "\t-o\tSpecify output radiation file path (default: stdout)\n");
	fprintf (stderr, "\t-r\tSpecify output RHS file path (default: none)\n");
	fprintf (stderr, "\t-m\tSpecify measurement distance (default: infinite)\n");

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
	int nspheres, nsptype, n, nterm, i;
	sptype *sparms, bgspt;
	spscat *slist;
	bgtype bg;
	exctparm exct;
	shdata shtr, shroot;
	itconf itc;
	trdesc *trans;
	double mslen = -1;

	complex double *rhs, *sol, *radpat, *sptr;

	FILE *fptr = NULL;
	char *inname = NULL, *outname = NULL, *rhsname = NULL, ch;

	while ((ch = getopt (argc, argv, "hi:o:r:m:")) != -1) {
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

	esbdinit (&bgspt, bg.k, 1.0, &shroot);
	fprintf (stderr, "Built data for enclosing sphere\n");

	nterm = shtr.ntheta * shtr.nphi + 2;
	n = nspheres * nterm + shroot.ntheta * shroot.nphi + 2;
	rhs = calloc (2 * n, sizeof(complex double));
	sol = rhs + n;

	n = nspheres * nterm;
	radpat = rhs + n;
	sptr = sol + n;

	/* Build the RHS vector consisting of multiple plane waves. */
	for (i = 0; i < exct.npw; ++i)
		shincident (shroot.deg, radpat, shroot.nphi,
				(exct.mag)[i], (exct.theta)[i], (exct.phi)[i]);

	/* Include the appropriate scaling factor. */
	shscale (radpat, &shroot, -1);

	/* Compute the transmitted component of the incident field. */
	spreflect (radpat, radpat, bgspt.transmit, shroot.deg, shroot.nphi, 0, 1);
	ifsht (radpat, &shroot);

	fprintf (stderr, "Built RHS vector\n");

	if (rhsname) {
		fptr = critopen (rhsname, "w");
		writebvec (fptr, radpat, shroot.nphi, shroot.ntheta);
		fclose (fptr);
	}

	itsolve (sol, rhs, slist, nspheres, &bgspt, trans, &shtr, &shroot, &itc);
	fprintf (stderr, "Iteration complete, computing radiation pattern\n");

	fflush (stdout);

	/* Compute the far-field pattern of the internal spheres. */
	neartofar (radpat, sol, slist, nspheres, bgspt.k, &shroot, &shtr);

	/* Convert the boundary fields to SH coefficients for reflections. */
	ffsht (radpat, &shroot);
	ffsht (sptr, &shroot);

	/* Transmit this field through the outer boundary. */
	spreflect (radpat, radpat, bgspt.transmit + shroot.deg, shroot.deg, shroot.nphi, 0, 1);

	/* Now add in the reflected incident wave coefficients. */
	spreflect (radpat, sptr, bgspt.reflect + shroot.deg, shroot.deg, shroot.nphi, 1, 1);

	/* Measure at a finite distance, if desired. */
	if (mslen > 0) {
		/* Scale the measurement distance to wavelengths. */
		mslen *= exct.f / bg.cabs;
		fprintf (stderr, "Computing scattered field at radius %g wavelengths\n", mslen);
		/* Scale the SH coefficients appropriately. */
		shscale (radpat, &shroot, 1);
		/* Now include the radial factor in the SH coefficients. */
		shradial (shroot.deg, radpat, shroot.nphi, bg.k, mslen);
	}

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

	free (exct.mag);
	free (exct.theta);

#ifdef _OPENMP
	fftw_cleanup_threads ();
#endif /* _OPENMP */

	fftw_cleanup ();

	return EXIT_SUCCESS;
}

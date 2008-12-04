#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <complex.h>
#include <math.h>

#include "fastsphere.h"
#include "config.h"
#include "init.h"
#include "util.h"
#include "scatmat.h"
#include "farfield.h"

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
	sptype *sparms;
	spscat *slist;
	bgtype bg;
	exctparm exct;
	shdata shtr, shroot;
	itconf itc;
	trdesc *trans;
	complex double *rhs, *sol, *radpat;

	double mslen = -1;

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

	readcfg (fptr, &nspheres, &nsptype, &sparms, &slist, &bg, &exct, &itc);
	fprintf (stderr, "Parsed configuration for %d spheres at %g MHz\n", nspheres, exct.f / 1e6);
	sphinit (sparms, nsptype, &bg, &shtr);
	fprintf (stderr, "Initialized spherical harmonic data for degree %d\n", shtr.deg);
	trans = sphbldfmm (slist, nspheres, &bg, &shtr);
	fprintf (stderr, "Built FMM translators for all spheres\n");

	i = rootorder (slist, nspheres, bg.k);
	n = 2 * i - 1; /* The number of angular samples (in each direction). */
	fshtinit (&shroot, i, n, 2 * n);
	fprintf (stderr, "Built spherical harmonic data for far-field\n");

	nterm = shtr.ntheta * shtr.nphi + 2;
	n = nspheres * nterm;
	rhs = calloc (2 * n, sizeof(complex double));
	sol = rhs + n;
	
	n = shroot.ntheta * shroot.nphi + 2;
	radpat = calloc (n, sizeof(complex double));

	/* Build the RHS vector consisting of multiple plane waves. */
	for (i = 0; i < exct.npw; ++i)
		shincident (shroot.deg, radpat, shroot.nphi,
				(exct.mag)[i], (exct.theta)[i], (exct.phi)[i]);

	/* Include the appropriate scaling factor. */
	shscale (radpat, &shroot, -1);
	ifsht (radpat, &shroot);

	fartonear (rhs, radpat, slist, nspheres, bg.k, &shtr, &shroot);

	/* Convert the incoming incident field to the reflected incident
	 * field, which is the RHS for this problem. */
	buildrhs (rhs, slist, nspheres, &shtr);
	fprintf (stderr, "Built RHS vector\n");

	if (rhsname) {
		fptr = critopen (rhsname, "w");
		writebvec (fptr, rhs, nterm, nspheres);
		fclose (fptr);
	}

	itsolve (sol, rhs, slist, nspheres, trans, &shtr, &itc);
	fprintf (stderr, "Iteration complete, computing radiation pattern\n");

	fflush (stdout);

	/* Compute the far-field radiation pattern of the object. */
	neartofar (radpat, sol, slist, nspheres, bg.k, &shroot, &shtr);

	if (mslen > 0) {
		mslen *= exct.f / bg.cabs;
		fprintf (stderr, "Computing scattered field at radius %g wavelengths\n", mslen);
		ffsht (radpat, &shroot);
		shscale (radpat, &shroot, 1);
		shradial (shroot.deg, radpat, shroot.nphi, bg.k, mslen);
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
	free (radpat);
	fshtfree (&shroot);
	fshtfree (&shtr);

	return EXIT_SUCCESS;
}

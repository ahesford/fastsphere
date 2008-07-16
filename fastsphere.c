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

	readcfg (fptr, &nspheres, &nsptype, &sparms, &bgspt, &slist, &bg, &exct, &itc);
	fprintf (stderr, "Parsed configuration for %d spheres at %g MHz\n", nspheres, exct.f / 1e6);
	sphinit (sparms, nsptype, bg.k, 1.0, &shtr);
	fprintf (stderr, "Initialized spherical harmonic data for degree %d\n", shtr.deg);
	trans = sphbldfmm (slist, nspheres, bg.k, &shtr);
	fprintf (stderr, "Built FMM translators for all spheres\n");

	n = rootorder (slist, nspheres, bg.k);
	n = 2 * n - 1; /* The number of angular samples (in each direction). */
	fshtinit (&shroot, shtr.deg, n, 2 * n);
	fprintf (stderr, "Built spherical harmonic data for far-field\n");

	nterm = shtr.ntheta * shtr.nphi;
	n = nspheres * nterm;
	rhs = calloc (2 * n, sizeof(complex double));
	sol = rhs + n;
	
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
	n = shroot.ntheta * shroot.nphi;
	radpat = malloc (n * sizeof(complex double));
	neartofar (radpat, sol, slist, nspheres, bg.k, &shroot, &shtr);

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

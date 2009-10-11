#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <complex.h>
#include <math.h>

#include "fastsphere.h"
#include "config.h"
#include "init.h"
#include "util.h"
#include "scatmat.h"
#include "farfield.h"
#include "spreflect.h"

void usage () {
	fprintf (stderr, "USAGE: fastsphere [-h] [-e] [-m mx [my mz [Mx My Mz]]]\n\t\t[-n Nx [Ny Nz]] [-i input] [-o output]\n\n");
	fprintf (stderr, "\t-h\tPrint this message and exit\n");
	fprintf (stderr, "\t-e\tSpecify the existence of an enclosing sphere in the input file\n");
	fprintf (stderr, "\t-m\tSpecify the lower and upper box corner in wavelengths\n\t\t(default: tight)\n");
	fprintf (stderr, "\t-n\tSpecify the number of pixels in each dimension (default: 100)\n");
	fprintf (stderr, "\t-i\tSpecify configuration file path (default: stdin)\n");
	fprintf (stderr, "\t-o\tSpecify output contrast file path (default: stdout)\n");

	exit (EXIT_FAILURE);
}

int writebvec (FILE *out, complex float *vec, int ndim, int *nelts) {
	int len, i;

	len = nelts[0];
	for (i = 1; i < ndim; ++i) len *= nelts[i];

	/* Write the size of the vector. */
	fwrite (nelts, sizeof(int), ndim, out);
	/* Write the vector itself. */
	fwrite (vec, sizeof(complex float), len, out);

	return len;
}

int insphere (double *pt, double *cen, double r) {
	double dist, dx[3];

	dx[0] = pt[0] - cen[0];
	dx[1] = pt[1] - cen[1];
	dx[2] = pt[2] - cen[2];

	/* Find the distance between the point and the center. */
	dist = sqrt (dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);

	/* Return whether or not the point is in the sphere. */
	return (dist <= r);
}

int bldcontrast (complex float *ct, int *nelt, double *bmin, double *bmax,
		sptype *bgs, spscat *slist, int nsphere) {
	double cell[3], zero[3] = {0, 0, 0};
	int ntot = nelt[0] * nelt[1] * nelt[2];

	cell[0] = (bmax[0] - bmin[0]) / nelt[0];
	cell[1] = (bmax[1] - bmin[1]) / nelt[1];
	cell[2] = (bmax[2] - bmin[2]) / nelt[2];

#pragma omp parallel default(shared)
	{
	double cen[3];
	int i, j, idx[3];
	complex float ctval;

	for (i = 0; i < ntot; ++i) {
		/* Find the cell index. */
		idx[0] = i % nelt[0];
		idx[1] = (i / nelt[0]) % nelt[2];
		idx[2] = i / (nelt[0] * nelt[1]);

		/* Find the cell center. */
		cen[0] = bmin[0] + ((double)idx[0] + 0.5) * cell[0];
		cen[1] = bmin[1] + ((double)idx[1] + 0.5) * cell[1];
		cen[2] = bmin[2] + ((double)idx[2] + 0.5) * cell[2];

		ctval = 0;

		/* Check if the point is in an enclosing sphere, if it exists. */
		if (bgs && insphere (cen, zero, bgs->r))
			ctval = (complex float)(bgs->k);

		/* If the point is in an inner sphere, set the wave number. */
		for (j = 0; j < nsphere; ++j) 
			if (insphere (cen, slist[j].cen, slist[j].spdesc->r))
				ctval = (complex float)(slist[j].spdesc->k);

		/* Convert the wave number to the contrast. */
		ctval /= (2 * M_PI);
		ctval = ctval * ctval - 1;

		/* Set the contrast value in the grid. */
		ct[i] = ctval;
	}
	}

	return ntot;
}

int main (int argc, char **argv) {
	int nspheres, nsptype, n, i;
	int autobox = 1, nelt[3] = {100, 100, 100};
	double boxlim[6];

	complex float *contrast;

	FILE *fptr = NULL;
	char *inname = NULL, *outname = NULL, ch;

	sptype *sparms, *bgptr = NULL, bgspt;
	spscat *slist;
	bgtype bg;
	exctparm exct;
	itconf itc;

	while ((ch = getopt (argc, argv, "hi:o:em:n:")) != -1) {
		switch (ch) {
		case 'i':
			inname = optarg;
			break;
		case 'o':
			outname = optarg;
			break;
		case 'e':
			bgptr = &bgspt;
			break;
		case 'm':
			/* Specify the box limits. */
			autobox = sscanf (optarg, "%lf %lf %lf %lf %lf %lf",
					boxlim, boxlim + 1, boxlim + 2,
					boxlim + 3, boxlim + 4, boxlim + 5);

			if (autobox == 3) {
				/* Set symmetric bounds from one corner. */
				boxlim[0] = -(boxlim[3] = ABS(boxlim[0]));
				boxlim[1] = -(boxlim[4] = ABS(boxlim[1]));
				boxlim[2] = -(boxlim[5] = ABS(boxlim[2]));
			} else if (autobox == 1) {
				/* Set symmetric bounds from one dimension. */
				boxlim[0] = boxlim[1] = boxlim[2] = -ABS(boxlim[0]);
				boxlim[3] = boxlim[4] = boxlim[5] =  ABS(boxlim[0]);
			} else if (autobox != 6) usage ();

			/* Don't automatically specify limits. */
			autobox = 0;
			break;
		case 'n':
			i = sscanf (optarg, "%d %d %d", nelt, nelt + 1, nelt + 2);

			if (i == 1) nelt[1] = nelt[2] = nelt[0];
			else if (i != 3) usage ();
			break;
		case 'h': default:
			usage ();
		}
	}

	if (!inname) fptr = stdin;
	else fptr = critopen (inname, "r");

	readcfg (fptr, &nspheres, &nsptype, &sparms, bgptr, &slist, &bg, &exct, &itc);
	fprintf (stderr, "Parsed configuration for %d spheres at %g MHz\n", nspheres, exct.f / 1e6);

	/* Automatically set box dimensions if necessary. */
	if (autobox && bgptr) {
		boxlim[0] = boxlim[1] = boxlim[2] = -bgspt.r;
		boxlim[3] = boxlim[4] = boxlim[5] = -bgspt.r;
	} else if (autobox) {
		/* Set the initial bounds to enclose the first sphere. */
		boxlim[0] = slist->cen[0] - slist->spdesc->r;
		boxlim[1] = slist->cen[1] - slist->spdesc->r;
		boxlim[2] = slist->cen[2] - slist->spdesc->r;
		boxlim[3] = slist->cen[0] + slist->spdesc->r;
		boxlim[4] = slist->cen[1] + slist->spdesc->r;
		boxlim[5] = slist->cen[2] + slist->spdesc->r;
		for (i = 1; i < nspheres; ++i) {
			boxlim[0] = MIN(boxlim[0], slist[i].cen[0] - slist[i].spdesc->r);
			boxlim[1] = MIN(boxlim[1], slist[i].cen[1] - slist[i].spdesc->r);
			boxlim[2] = MIN(boxlim[2], slist[i].cen[2] - slist[i].spdesc->r);
			boxlim[3] = MAX(boxlim[3], slist[i].cen[0] + slist[i].spdesc->r);
			boxlim[4] = MAX(boxlim[4], slist[i].cen[1] + slist[i].spdesc->r);
			boxlim[5] = MAX(boxlim[5], slist[i].cen[2] + slist[i].spdesc->r);
		}
	}

	/* Allocate the contrast map. */
	contrast = malloc (nelt[0] * nelt[1] * nelt[2] * sizeof(complex float));

	bldcontrast (contrast, nelt, boxlim, boxlim + 3, bgptr, slist, nspheres);

	if (!outname) fptr = stdout;
	else fptr = critopen (outname, "w");
	fprintf (stderr, "Writing contrast file.\n");
	writebvec (fptr, contrast, 3, nelt);
	fclose (fptr);

	clrspheres (sparms, nsptype);
	free (exct.mag);
	free (exct.theta);
	free (contrast);

	return EXIT_SUCCESS;
}

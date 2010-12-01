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
	fprintf (stderr, "USAGE: spherepix [-h] [-e] [-m mx [my mz [Mx My Mz]]]\n\t\t[-n Nx [Ny Nz]] [-i input] [-o output]\n\n");
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

/* Determin if a point is inside a sphere. */
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

/* Compute the Laplacian of the inverse of the square root of density. The
 * density is already rooted. Watch for edges of the domain. */
float lapden (float *r, float *lr, float *nr, double *c, int pos, int *nelt) {
	float dlap, nval, pval;
	int x, y;

	x = pos % nelt[0];
	y = pos / nelt[0];

	/* Contribution of the x offsets with bounds checking. */
	if (x >= nelt[0] - 1) nval = 1.0;
	else nval = 1.0 / r[pos + 1];

	if (x <= 0) pval = 1.0;
	else pval = 1.0 / r[pos - 1];

	dlap = (pval + nval - 2.0 / r[pos]) / (c[0] * c[0]);

	/* Contribution of the y offsets with bounds checking. */
	if (y >= nelt[1] - 1) nval = 1.0;
	else nval = 1.0 / r[pos + nelt[0]];

	if (y <= 0) pval = 1.0;
	else pval = 1.0 / r[pos - nelt[0]];

	dlap += (pval + nval - 2.0 / r[pos]) / (c[1] * c[1]);

	/* Contribution of the z offsets with bounds checking. */
	if (!nr) nval = 1.0;
	else nval = 1.0 / nr[pos];

	if (!lr) pval = 1.0;
	else pval = 1.0 / lr[pos];

	dlap += (pval + nval - 2.0 / r[pos]) / (c[2] * c[2]);

	return dlap;
}

int augct (complex float *k, float *r, float *lr, float *nr, int *nelt, double *cell) {
	int i, npx = nelt[0] * nelt[1];
	float dval;

#pragma omp parallel for default(shared) private(i,dval)
	for (i = 0; i < npx; ++i) {
		/* Compute the Laplacian of the inverse square root of the density. */
		dval = lapden (r, lr, nr, cell, i, nelt);
		/* Scale by density square root and normalize by wave number. */
		dval *= r[i] / (4.0 * M_PI * M_PI);

		/* Subtract the density term from the contrast. */
		k[i] -= dval;
	}

	return npx;
}

/* Build the contrast and density maps for a slab. */
int bldct (complex float *ct, float *density, int *nelt, double *blim,
		double *cell, sptype *bgs, spscat *slist, int nsphere, int zidx) {
	double zero[3] = {0, 0, 0};
	int ntot = nelt[0] * nelt[1];

#pragma omp parallel default(shared)
{
	double cen[3];
	int i, j, idx[2];
	complex float ctval;
	float dval;
	
	cen[2] = blim[2] + ((double)zidx + 0.5) * cell[2];

	/* Build the density-free contrast and the density map. */
#pragma omp for
	for (i = 0; i < ntot; ++i) {
		/* Find the cell index. */
		idx[0] = i % nelt[0];
		idx[1] = i / nelt[0];

		/* Find the cell center. */
		cen[0] = blim[0] + ((double)idx[0] + 0.5) * cell[0];
		cen[1] = blim[1] + ((double)idx[1] + 0.5) * cell[1];

		/* Set the background contrast. */
		ctval = 2 * M_PI;
		dval = 1.0;

		/* Check if the point is in an enclosing sphere, if it exists. */
		if (bgs && insphere (cen, zero, bgs->r)) {
			ctval = (complex float)(bgs->k);
			dval = (float)(bgs->rho);
		}

		/* If the point is in an inner sphere, set the wave number. */
		for (j = 0; j < nsphere; ++j) 
			if (insphere (cen, slist[j].cen, slist[j].spdesc->r)) {
				ctval = (complex float)(slist[j].spdesc->k);
				dval = (float)(slist[j].spdesc->rho);
			}

		/* Convert the wave number to the contrast. */
		ctval /= (2 * M_PI);
		ctval = ctval * ctval - 1;

		/* Set the contrast value and density in the grid. */
		ct[i] = ctval;
		density[i] = sqrt(dval);
	}
}
	return ntot;
}

int main (int argc, char **argv) {
	int nspheres, nsptype, n, i, npx;
	int autobox = 1, nelt[3] = {100, 100, 100};
	double boxlim[6], cell[3];

	complex float *k, *nk, *kslab;
	float *density, *lr, *r, *nr;

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

			switch (autobox) {
			case 1:
				/* Set symmetric bounds from one dimension. */
				boxlim[0] = boxlim[1] = boxlim[2] = -ABS(boxlim[0]);
				boxlim[3] = boxlim[4] = boxlim[5] =  ABS(boxlim[0]);
				break;
			case 3:
				/* Set symmetric bounds from one corner. */
				boxlim[0] = -(boxlim[3] = ABS(boxlim[0]));
				boxlim[1] = -(boxlim[4] = ABS(boxlim[1]));
				boxlim[2] = -(boxlim[5] = ABS(boxlim[2]));
				break;
			case 6:
				/* Nothing to te done for fully specified box. */
				break;
			default:
				usage ();
			}

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

	/* Compute the cell dimensions. */
	cell[0] = (boxlim[3] - boxlim[0]) / nelt[0];
	cell[1] = (boxlim[4] - boxlim[1]) / nelt[1];
	cell[2] = (boxlim[5] - boxlim[2]) / nelt[2];

	npx = nelt[0] * nelt[1];

	/* Allocate the contrast and density map for a slab. */
	kslab = malloc (2 * npx * sizeof(complex float));
	density = malloc (3 * npx * sizeof(float));

	/* Point to the slab data stores. */
	k = kslab;
	nk = k + npx;
	lr = NULL;
	r = density;
	nr = r + npx;

	if (!outname) fptr = stdout;
	else fptr = critopen (outname, "w");
	fprintf (stderr, "Writing contrast file.\n");

	/* Write the header. */
	fwrite (nelt, sizeof(int), 3, fptr);

	/* Construct the first slab of data. */
	bldct (k, r, nelt, boxlim, cell, bgptr, slist, nspheres, 0);

	for (i = 1; i < nelt[2]; ++i) {
		/* Construct the next slab of data. */
		bldct (nk, nr, nelt, boxlim, cell, bgptr, slist, nspheres, i);

		/* Build and write the previous slab. */
		augct (k, r, lr, nr, nelt, cell);
		fwrite (k, sizeof(complex float), npx, fptr);

		/* Update the media pointers. */
		k = kslab + (i % 2) * npx;
		nk = kslab + ((i + 1) % 2) * npx;

		lr = density + ((i - 1) % 3) * npx;
		r = density + (i % 3) * npx;
		nr = density + ((i + 1) % 3) * npx;
	}

	/* Build and write the last slab. */
	augct (k, r, lr, NULL, nelt, cell);
	fwrite (k, sizeof(complex float), npx, fptr);

	fclose (fptr);

	clrspheres (sparms, nsptype);
	free (exct.mag);
	free (exct.theta);
	free (kslab);
	free (density);

	return EXIT_SUCCESS;
}

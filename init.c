#include <stdlib.h>
#include <complex.h>
#include <math.h>

#include "fastsphere.h"
#include "translator.h"
#include "util.h"
#include "init.h"

int sphinit (sptype *spt, int nspt, bgtype *bg, shdata *shtr) {
	int i, deg, ntheta;
	double maxrad = 0.0;
	sptype *sptr;

	/* Find the maximum sphere radius. */
	for (i = 0, sptr = spt; i < nspt; ++i, ++sptr)
		maxrad = MAX (maxrad, sptr->r);

	/* Use six digits of accuracy. */
	deg = exband (bg->k * maxrad, 6);
	ntheta = 2 * deg - 1;

	/* Initialize the SH transform data. */
	fshtinit (shtr, deg, ntheta);

	/* Initialize and populate the SH reflection coefficients. */
#pragma omp parallel for private(i,sptr) default(shared)
	for (i = 0; i < nspt; ++i) {
		sptr = spt + i;
		sptr->reflect = malloc (deg * sizeof(complex double));
		spbldrc (sptr->reflect, bg->k, sptr->k, sptr->r, deg);
	}

	return nspt;
}

void clrspheres (sptype *spt, int nspt) {
	int i;

#pragma omp parallel for private(i) default(shared)
	for (i = 0; i < nspt; ++i) free (spt[i].reflect);
}

int sphbldfmm (complex double ***trans, spscat *sph, int nsph,
		bgtype *bg, shdata *shtr) {
	int i, j, k, nsq, nterm, trunc;
	complex double *tptr;
	long ntrans;

	nsq = nsph * nsph;
	nterm = shtr->ntheta * shtr->nphi;
	trunc = 2 * shtr->deg - 1;

	*trans = malloc (nsq * sizeof(complex double *));

	/* Allocate space for all translators. */
	ntrans = (long) (nsq - nsph) * (long) nterm;

	/* The **trans pointer is never used, but keep it pointing to the
	 * backing store for convenience of freeing later on. */
	**trans = malloc (ntrans * sizeof(complex double));
	tptr = **trans;

	/* Set up the translator pointers into the backing store.
	 * This is a separate loop because the construction loop
	 * should be parallelized. */
	for (k = 1; k < nsq; ++k) {
		j = k / nsph;	/* Source sphere. */
		i = k % nsph;	/* Destination sphere. */

		if (i == j) {
			(*trans)[k] = NULL;
			continue;
		}

		(*trans)[k] = tptr;
		tptr += nterm;
	}

#pragma omp parallel private(i,j,k) default(shared)
{
	/* These variables are private. */
	complex double kr;
	double sdir[3], dist;

	/* Start with the first valid translator. */
#pragma omp for
	for (k = 1; k < nsq; ++k) {
		j = k / nsph;	/* Source sphere. */
		i = k % nsph;	/* Destination sphere. */

		/* This translation is never used. */
		if (i == j) continue;

		/* Translation direction. */
		sdir[0] = sph[i].cen[0] - sph[j].cen[0];
		sdir[1] = sph[i].cen[1] - sph[j].cen[1];
		sdir[2] = sph[i].cen[2] - sph[j].cen[2];

		/* Translation distance. */
		dist = sqrt(sdir[0] * sdir[0] + sdir[1] * sdir[1] + sdir[2] * sdir[2]);
		kr = bg->k * dist;

		/* Normalize translation direction. */
		sdir[0] /= dist; sdir[1] /= dist; sdir[2] /= dist;

		/* Build the translator. */
		translator ((*trans)[k], trunc, shtr->ntheta, shtr->nphi,
				shtr->theta, kr, sdir);
	}
}

	return ntrans;
}

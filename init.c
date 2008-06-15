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
	for (i = 0, sptr = spt; i < nspt; ++i, ++sptr) {
		sptr->reflect = malloc (deg * sizeof(complex double));
		spbldrc (sptr->reflect, bg->k, sptr->k, sptr->r, deg);
	}

	return nspt;
}

void clrspheres (sptype *spt, int nspt) {
	int i;
	sptype *sptr;

	for (i = 0, sptr = spt; i < nspt; ++i, ++sptr) free (sptr->reflect);
}

int sphbldfmm (complex double ***trans, spscat *sph, int nsph,
		bgtype *bg, shdata *shtr) {
	int i, j, k, nsq, nterm, trunc;
	complex double kr, *tptr;
	double sdir[3], dist;
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

	/* Start with the first valid translator. */
	for (k = 1; k < nsq; ++k) {
		j = k / nsph;		/* Source sphere. */
		i = mod (k, nsph);	/* Destination sphere. */

		if (i == j) {
			/* This translation is never used. */
			(*trans)[k] = NULL;
			continue;
		}

		/* Translation direction. */
		sdir[0] = sph[i].cen[0] - sph[j].cen[0];
		sdir[1] = sph[i].cen[1] - sph[j].cen[1];
		sdir[2] = sph[i].cen[2] - sph[j].cen[2];

		/* Translation distance. */
		dist = sdir[0] * sdir[0] + sdir[1] * sdir[1] + sdir[2] * sdir[2];
		dist = sqrt(dist);

		kr = bg->k * dist;

		/* Normalize translation direction. */
		sdir[0] /= dist; sdir[1] /= dist; sdir[2] /= dist;

		/* Build the translator. */
		translator (tptr, trunc, shtr->ntheta, shtr->nphi,
				shtr->theta, kr, sdir);

		/* Set the pointer to the translator. */
		(*trans)[k] = tptr;
		tptr += nterm;
	}

	return ntrans;
}

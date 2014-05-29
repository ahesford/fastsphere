#include <stdlib.h>
#include <complex.h>
#include <math.h>

#include "fastsphere.h"
#include "translator.h"
#include "spreflect.h"
#include "util.h"
#include "init.h"

int sphinit (sptype *spt, int nspt, complex double bgk,
		double bgrho, shdata *shtr, int nt, int ndig) {
	int i, deg = 0;
	sptype *sptr;

	/* Find the maximum spherical harmonic degree required. */
	for (i = 0, sptr = spt; i < nspt; ++i, ++sptr) {
		sptr->deg = exband (bgk * sptr->r, ndig);
		deg = MAX(deg, sptr->deg);
	}

	/* Use a default value for angular samples if none specified.
	 * The default value is forced to be odd. */
	if (nt < 1) nt = deg + (deg % 2) + 1;

	/* Initialize the SH transform data. */
	fshtinit (shtr, deg, nt, 2 * nt, 1);

	/* Initialize and populate the SH reflection coefficients. */
#pragma omp parallel for private(i,sptr) default(shared)
	for (i = 0; i < nspt; ++i) {
		sptr = spt + i;
		/* The reflection coefficients go out to the maximum degree,
		 * but are all zero beyond the sphere's maximum degree. */
		sptr->reflect = calloc (sptr->deg, sizeof(complex double));
		/* The background relative density is unity, and the sphere
		 * density is specified relative to the background. */
		spbldrc (sptr, bgk, bgrho);
	}

	return nspt;
}

int esbdinit (sptype *sbd, complex double bgk, 
		double bgrho, shdata *shtr, int nang, int ndig) {
	/* The spherical harmonic degree for the background sphere. */
	sbd->deg = exband (bgk * sbd->r, ndig);
	/* The number of angular samples in theta is forced to be odd,
	 * and at least equal to the degree required. */
	nang = MAX(sbd->deg + (sbd->deg % 2) + 1, nang);

	/* Initialize the SH transform data. */
	fshtinit (shtr, sbd->deg, nang, 2 * nang, 1);

	/* Allocate space for reflection and transmission coefficients. */
	sbd->reflect = calloc (4 * sbd->deg, sizeof(complex double));
	sbd->transmit = sbd->reflect + 2 * sbd->deg;

	/* Compute the reflection and transmission coefficients. */
	esbldrc (sbd, bgk, bgrho);

	return 1;
}

void clrspheres (sptype *spt, int nspt) {
	int i;

#pragma omp parallel for private(i) default(shared)
	for (i = 0; i < nspt; ++i) free (spt[i].reflect);
}

trdesc* sphbldfmm (spscat *sph, int nsph, complex double bgk, shdata *shtr) {
	int i, j, k, nsq, nterm;
	trdesc *trans;

	nsq = nsph * nsph;
	nterm = shtr->ntheta * shtr->nphi;

	trans = calloc (nsq, sizeof(trdesc));

	/* Set up the translators in parallel. */
#pragma omp parallel private(i,j,k) default(shared)
{
	double dist, *sdir;

#pragma omp for
	for (k = 1; k < nsq; ++k) {
		j = k / nsph;	/* Source sphere. */
		i = k % nsph;	/* Destination sphere. */

		if (i == j) {
			trans[k].type = TRNONE;
			continue;
		}

		/* Convenient pointer to the translation axis. */
		sdir = trans[k].sdir;

		/* Translation direction. */
		sdir[0] = sph[i].cen[0] - sph[j].cen[0];
		sdir[1] = sph[i].cen[1] - sph[j].cen[1];
		sdir[2] = sph[i].cen[2] - sph[j].cen[2];

		/* Translation distance. */
		dist = sqrt(DVDOT(sdir,sdir));
		trans[k].kr = bgk * dist;

		/* Normalize translation direction. */
		sdir[0] /= dist; sdir[1] /= dist; sdir[2] /= dist;

		trans[k].trunc = MAX(sph[j].spdesc->deg, sph[i].spdesc->deg);

		/* Set the translator type. */
		trans[k].type = TRPLANE;

		/* Allocate and build the translator array. */
		trans[k].trdata = malloc (nterm * sizeof(complex double));
		translator (trans + k, shtr->ntheta, shtr->nphi, shtr->theta);
	}
}

	return trans;
}

void sphclrfmm (trdesc *trans, int ntrans) {
	int i;

	for (i = 0; i < ntrans; ++i)
		if (trans[i].type == TRPLANE && trans[i].trdata)
			free (trans[i].trdata);
}

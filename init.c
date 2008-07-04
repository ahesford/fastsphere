#include <stdlib.h>
#include <complex.h>
#include <math.h>

#include "fastsphere.h"
#include "translator.h"
#include "spreflect.h"
#include "shrotate.h"
#include "util.h"
#include "init.h"

int sphinit (sptype *spt, int nspt, bgtype *bg, shdata *shtr) {
	int i, deg = 0, nang;
	sptype *sptr;

	/* Find the maximum spherical harmonic degree required. */
	for (i = 0, sptr = spt; i < nspt; ++i, ++sptr) {
		sptr->deg = exband (bg->k * sptr->r, 6);
		deg = MAX(deg, sptr->deg);
	}

	nang = 2 * deg - 1; /* The number of angular samples (per dimension). */

	/* Initialize the SH transform data. */
	fshtinit (shtr, deg, nang, 2 * nang);

	/* Initialize and populate the SH reflection coefficients. */
#pragma omp parallel for private(i,sptr) default(shared)
	for (i = 0; i < nspt; ++i) {
		sptr = spt + i;
		/* The reflection coefficients go out to the maximum degree,
		 * but are all zero beyond the sphere's maximum degree. */
		sptr->reflect = calloc (deg, sizeof(complex double));
		/* The background relative density is unity, and the sphere
		 * density is specified relative to the background. */
		spbldrc (sptr->reflect, bg->k, sptr->k,
				1.0, sptr->rho, sptr->r, sptr->deg);
	}

	return nspt;
}

void clrspheres (sptype *spt, int nspt) {
	int i;

#pragma omp parallel for private(i) default(shared)
	for (i = 0; i < nspt; ++i) free (spt[i].reflect);
}

trdesc* sphbldfmm (spscat *sph, int nsph, bgtype *bg, shdata *shtr) {
	int i, j, k, nsq, nterm;
	trdesc *trans;

	nsq = nsph * nsph;
	nterm = shtr->ntheta * shtr->nphi;

	trans = calloc (nsq, sizeof(trdesc));

	/* Set up the translators in parallel. */
#pragma omp parallel private(i,j,k) default(shared)
{
	double dist, *sdir, trunc, rad;

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
		trans[k].kr = bg->k * dist;

		/* Normalize translation direction. */
		sdir[0] /= dist; sdir[1] /= dist; sdir[2] /= dist;

		trunc = MAX((sph[i].spdesc)->deg, (sph[j].spdesc)->deg);
		trans[k].trunc = 2 * trunc - 1;

		/* Find the rotation angles of the translation direction. */
		getangles (&(trans[k].theta), &(trans[k].chi), &(trans[k].phi), sdir);

		/* Find the maximum radius for this translation. */
		rad = MAX((sph[i].spdesc)->r, (sph[j].spdesc)->r);

		/* If this is a dense translation, do nothing else. */
		if (dist < 2 * rad) {
			trans[k].type = TRDENSE;
			continue;
		}

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

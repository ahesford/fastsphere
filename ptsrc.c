#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <float.h>

#include "fastsphere.h"
#include "fsht.h"
#include "translator.h"
#include "util.h"

/* Compute directivity envelope dot(s, ax) * exp(alpha * (1 - dot(s, ax)**2)),
 * where s and ax are both unit vectors. If ax is NULL or has zero magnitude,
 * 1.0 is always returned to indicate an omnidirectional directivity envelope. */
double directivity (double *ax, double *s, double alpha) {
	double ca;

	/* No directivity axis means the pattern is omnidirectional. */
	if (!ax || DVDOT(ax, ax) < DBL_EPSILON) return 1.0;

	ca = DVDOT(ax, s);
	return ca * exp(alpha * (1. - ca * ca));
}

/* Compute the incoming far-field expansion of a point source embedded in a
 * material with background wave number bgk and located at coordinates loc
 * relative to the center of the insonified sphere described by the shdata
 * structure sphere. The expansion is added to the representations in spec. If
 * ax is not null and does not have zero magnitude, it represents the
 * directivity axis of the source that has a beam width parameter alpha. The
 * axis must be normalized by the caller. If ax is NULL or has zero magnitude,
 * an undirected point source is used regardless of the value of alpha. */
int ptsrcexp (complex double *spec, complex double bgk, shdata *sphere, 
		complex double mag, double *loc, double *ax, double alpha) {
	trdesc trans;
	int nterm, i, j;
	double s[3], dist, phi, dphi, st;
	complex double *tptr, *sptr, scale = mag * bgk / (4 * M_PI);

	/* Allocate space to hold the far-field expansion of a point source. */
	nterm = sphere->ntheta * sphere->nphi;
	
	/* The translation direction is the negative of the source direction. */
	trans.sdir[0] = -loc[0];
	trans.sdir[1] = -loc[1];
	trans.sdir[2] = -loc[2];

	/* Translation distance. */
	dist = sqrt(DVDOT(trans.sdir, trans.sdir));
	trans.kr = bgk * dist;

	/* Normalize translation direction. */
	trans.sdir[0] /= dist; trans.sdir[1] /= dist; trans.sdir[2] /= dist;

	trans.trunc = sphere->deg;
	trans.type = TRPLANE;

	/* Build the translator representing an undirected point-source. */
	trans.trdata = malloc (nterm * sizeof(complex double));
	translator (&trans, sphere->ntheta, sphere->nphi, sphere->theta);

	/* There is no directivity axis if it has zero magnitude. */
	if (ax && DVDOT(ax, ax) < DBL_EPSILON) ax = NULL;

	dphi = 2 * M_PI / MAX(sphere->nphi, 1);

#pragma omp critical(incaug)
	for (i = 0, tptr = trans.trdata, sptr = spec; i < sphere->ntheta; ++i) {
		s[2] = sphere->theta[i];
		st = sin(acos(s[2]));
		for (j = 0; j < sphere->nphi; ++j, ++tptr, ++sptr) {
			phi = j * dphi;
			s[0] = st * cos(phi);
			s[1] = st * sin(phi);
			*sptr += (*tptr) * scale * directivity (ax, s, alpha);
		}
	}

	free (trans.trdata);
	return nterm;
}

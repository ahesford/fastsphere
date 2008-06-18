#include <stdlib.h>

#include <math.h>
#include <complex.h>

#include "spbessel.h"
#include "translator.h"
#include "util.h"

/* Build the translator coefficient of order ord for precomputed Hankel terms
 * in the direction specified by angles theta and phi. The relative
 * translation direction is provided by sdir. lgwork is a workspace. */
complex double transang (int ord, complex double *hfn, double *lgwork,
		double *sdir, double theta, double phi) {
	double s[3], st, sds;
	complex double tsol = 0;
	int i;

	st = sin(theta);
	s[0] = st * cos(phi);
	s[1] = st * sin(phi);
	s[2] = cos(theta);

	sds = DVDOT(s,sdir);

	/* Compute the Legendre polynomials. */
	legpoly (ord, sds, lgwork);

	/* Sum the translator for the specified angle. */
	for (i = 0; i < ord; ++i) tsol += lgwork[i] * hfn[i];

	return tsol;
}

/* Compute the translator of order ord for a given argument kr, with a
 * reference direction sdir. The values are stored in trans. In the theta
 * plane, ntheta samples are chosen according to Legendre-Gauss quadrature,
 * plus one value at each pole. There are nphi equally-spaced phi samples. */
int translator (complex double *trans, int ord, int ntheta, int nphi,
		double *theta, complex double kr, double *sdir) {
	double phi, dphi, *lgwork;
	complex double *hfn, *tptr, cscale[4] = { 1.0, I, -1.0, -I };
	int i, j, iscale;

	if (ord < 0) return -1;

	/* Allocate storage space for the Legendre and Hankel functions. */
	hfn = malloc (ord * sizeof(complex double));
	lgwork = malloc (ord * sizeof(double));

	spbesh (hfn, kr, ord);

	/* The value (1i)^i wraps every fourth value of index i, so use the
	 * cscale array to exploit that. */
	/* iscale is 2 * i + 1, so fold that into the loop. */
	for (i = 0, iscale = 1; i < ord; ++i, iscale += 2)
		hfn[i] *= cscale[i % 4] * iscale;

	/* Set the phi increment and the theta samples. */
	dphi = 2 * M_PI / MAX(nphi, 1);

	/* Build the non-pole translator values. */
	for (i = 0, tptr = trans; i < ntheta; ++i) {
		for (j = 0; j < nphi; ++j, ++tptr) {
			phi = j * dphi;
			*tptr = transang (ord, hfn, lgwork, sdir, theta[i], phi);
		}
	}

	free (lgwork);
	free (hfn);
}

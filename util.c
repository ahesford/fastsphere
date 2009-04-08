#include <stdlib.h>
#include <string.h>

#include <math.h>
#include <complex.h>
#include <float.h>

#include <gsl/gsl_sf_legendre.h>

#include "util.h"
#include "spbessel.h"


/* Compute the number of harmonics required using the excess bandwidth
 * formula. */
int exband (complex double kr, double ndig) {
	double akr;
	int l;
	
	akr = creal (kr);
	ndig *= ndig;

	l = ceil (akr + 1.8 * cbrt(ndig * akr));

	return l;
}

/* Computes the Legendre polynomials up to order n for the argument x. The
 * values are stored in the array v. */
int legpoly (int n, double x, double *v) {
	int i;

	/* Don't bother computing anything for order less than zero. */
	if (n < 0) return -1;

	/* Make sure the argument is within [-1,1] as expected. */
	if (fabs(x) > 1.0 + DBL_EPSILON) return -2;

	/* The first function value. */
	v[0] = 1.0;

	/* If the order happens to be zero, there is no going further. */
	if (n < 1) return 0;

	/* The second function value. */
	v[1] = x;

	/* The recursion formula for Legendre polynomials. */
	for (i = 1; i < n - 1; ++i) 
		v[i + 1] = ((2 * i + 1) * x * v[i] - i * v[i - 1]) / (i + 1);

	return 0;
}

/* Computes the wave number from the relative sound speed and
 * unitless (dB) attenuation coefficient. */
complex double buildkvec (double cr, double alpha) {
	double kr, ki;

	ki = log(10) * alpha / 20;
	kr = 2 * M_PI / cr;

	return kr + I * ki;
}

/* Copy the spherical harmonic representation from one location to another. */
int copysh (int deg, complex double *out, int ldo, complex double *in, int ldi) {
	int i, j, offo, offi;

	/* Copy the spherical harmonic coefficients into the right place. */
	for (i = 0; i < deg; ++i) {
		offo = i * ldo;
		offi = i * ldi;
		/* Copy the zero-order coefficients. */
		out[offo] = in[offi];

		/* Copy the other coefficients. */
		for (j = 1; j <= i; ++j) {
			out[offo + j] = in[offi + j];
			out[offo + ldo - j] = in[offi + ldi - j];
		}
	}

	return deg;
}

/* Multiply in a radial component to the spherical harmonic. */
int shradial (int deg, complex double *a, int lda, complex double k, double r) {
	complex double kr, *hlkr;
	int i, j, off;

	hlkr = malloc (deg * sizeof(complex double));

	/* Compute the Hankel functions for all degrees. */
	kr = k * r;
	spbesh (hlkr, kr, deg);

	for (i = 0; i < deg; ++i) {
		off = i * lda;

		/* The zero-order coefficients. */
		a[off] *= hlkr[i];

		/* The other coefficients. */
		for (j = 1; j <= i; ++j) {
			a[off + j] *= hlkr[i];
			a[off + lda - j] *= hlkr[i];
		}
	}

	free (hlkr);

	return deg;
}

/* Compute the harmonic coefficients for an incident plane wave. */
int shincident (int deg, complex double *a, int lda,
		complex double mag, double theta, double phi) {
	double cth, *lgvals;
	int l, m, dm1, npm, off;
	complex double scale[4] = { 1, I, -1, -I }, cx, fx;

	dm1 = deg - 1;
	cth = cos(theta);

	lgvals = malloc (deg * sizeof(double));

	/* Compute the Legendre polynomials of zero order for all degrees. */
	gsl_sf_legendre_sphPlm_array (dm1, 0, cth, lgvals);

	/* Compute the zero-order coefficients for all degrees. */
	for (l = 0; l < deg; ++l) a[l * lda] += scale[l % 4] * mag * lgvals[l];

	for (m = 1; m < deg; ++m) {
		npm = lda - m;

		/* Compute the Legendre polynomials of fixed order for all degrees. */
		gsl_sf_legendre_sphPlm_array (dm1, m, cth, lgvals);

		/* Compute the phi variation. */
		cx = cexp (I * m * phi);

#pragma omp critical(incaug)
		for (l = m; l < deg; ++l) {
			off = l * lda;
			fx = scale[l % 4] * mag * lgvals[l - m];
			/* The positive-order coefficient. */
			a[off + m] += fx * conj (cx);
			/* The negative-order coefficient. */
			a[off + npm] += fx * cx;
		}
	}

	free (lgvals);

	return deg;
}

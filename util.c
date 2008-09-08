#include <string.h>

#include <math.h>
#include <complex.h>
#include <float.h>

#include "util.h"

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

/* Binary search to find an interval containing a data point. Return is the
 * index of the first value greater than or equal to the test value. */
int interval (double val, double *arr, int n) {
	int low = 0, high = n, mid;

	while (low < high) {
		mid = low + ((high - low) / 2);

		if (arr[mid] < val) low = mid + 1;
		else high = mid;
	}

	return low;
}

/* Return the Lagrange interpolation coefficients. */
int lgpoly (double *vals, double *nd, double sloc, int n) {
	int i, j;
	double den;

	for (i = 0; i < n; ++i) {
		den = vals[i] = 1;

		for (j = 0; j < i; ++j) {
			vals[i] *= sloc - nd[j];
			den *= nd[i] - nd[j];
		}

		for (j = i + 1; j < n; ++j) {
			vals[i] *= sloc - nd[j];
			den *= nd[i] - nd[j];
		}

		vals[i] /= den;
	}

	return n;
}

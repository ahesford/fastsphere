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

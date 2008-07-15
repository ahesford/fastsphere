#include <stdio.h>
#include <stdlib.h>

#include <complex.h>
#include <math.h>

#include "spbessel.h"
#include "util.h"

complex double dspbesj0 (complex double);
complex double dspbesh0 (complex double);
void spbesjinit (complex double *, complex double, int);
void spbeshinit (complex double *, complex double, int);

int spbesh (complex double *vals, complex double z, int n) {
	int i, nmo;

	if (n < 0) return -1;

	/* Build the first orders of the Hankel function. */
	spbeshinit (vals, z, n);

	/* Recursion fills out the rest of the orders. */
	nmo = n - 1;
	for (i = 1; i < nmo; ++i)
		vals[i+1] = (2 * i + 1) * vals[i] / z - vals[i - 1];

	return 0;
}

int spbesj (complex double *vals, complex double z, int n) {
	int i, nmo;

	if (n < 0) return -1;

	/* Build the first orders of the Hankel function. */
	spbesjinit (vals, z, n);

	/* Recursion fills out the rest of the orders. */
	nmo = n - 1;
	for (i = 1; i < nmo; ++i)
		vals[i+1] = (2 * i + 1) * vals[i] / z - vals[i - 1];

	return 0;
}

/* Compute the zero-order and unity-order spherical Bessel function and store
 * in jf[0] and jf[1], respectively. */
void spbesjinit (complex double *jf, complex double x, int nmax) {
	jf[0] = csin(x) / x;

	/* Don't compute the unity-order if nmax says not to. */
	if (nmax < 1) return;

	jf[1] = (jf[0] - ccos(x)) / x;
}

/* Compute the zero-order and unity-order spherical Hankel function and store
 * in hf[0] and hf[1], respectively. */
void spbeshinit (complex double *hf, complex double x, int nmax) {
	complex double cxfact;

	cxfact = -cexp (I * x) / x;

	hf[0] = I * cxfact;

	/* Don't compute the unity-order if nmax says not to. */
	if (nmax < 1) return;

	hf[1] = cxfact * (1 + I / x);
}

/* Returns the derivative of the zero-order spherical Hankel function. */
complex double dspbesh0 (complex double x) {
	complex double fact;

	fact = (x + I) / (x * x);

	return fact * cexp (I * x);
}

/* Returns the derivative of the zero-order spherical Bessel function. */
complex double dspbesj0 (complex double x) {
	complex double sx, cx;

	sx = csin (x) / (x * x);
	cx = ccos (x) / x;

	return cx - sx;
}

/* Computes the derivatives of the spherical Hankel functions, given that the
 * Hankel functions have been precomputed. */
int dspbesh (complex double *vals, complex double *hfn, complex double z, int n) {
	int i;

	/* The first value is not computed using the recursion formula. */
	vals[0] = dspbesh0 (z);

	/* Use the recursion formula. Observe that hfn must contain values of
	 * the Hankel function including order n, which means it must have
	 * n + 1 elements, for the formula to find the order n - 1
	 * derivative. */
	for (i = 1; i < n; ++i) 
		vals[i] = 0.5 * (hfn[i - 1] - hfn[i + 1]) - hfn[i] / (2 * z);

	return 0;
}

/* Computes the derivatives of the spherical Bessel functions, given that the
 * Bessel functions have been precomputed. */
int dspbesj (complex double *vals, complex double *hfn, complex double z, int n) {
	int i;

	/* The first value is not computed using the recursion formula. */
	vals[0] = dspbesj0 (z);

	/* Use the recursion formula. Observe that hfn must contain values of
	 * the Bessel function including order n, which means it must have
	 * n + 1 elements, for the formula to find the order n - 1
	 * derivative. */
	for (i = 1; i < n; ++i) 
		vals[i] = 0.5 * (hfn[i - 1] - hfn[i + 1]) - hfn[i] / (2 * z);

	return 0;
}

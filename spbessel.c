#include <stdio.h>
#include <stdlib.h>

#include <complex.h>
#include <math.h>

#include "spbessel.h"
#include "util.h"

typedef void (besinit)(complex double *, complex double, int);
typedef complex double (cfunc)(complex double);

besinit spbesjinit;
besinit spbeshinit;
cfunc dspbesj0;
cfunc dspbesh0;
cfunc ddspbesj0;
cfunc ddspbesh0;

int spbessel (complex double *, complex double, int, besinit);
int dspbessel (complex double *, complex double *, complex double, int, cfunc);
int ddspbessel (complex double *, complex double *, complex double, int, cfunc);

/* This routine uses a common recurrence formula to compute, in the array vals,
 * any spherical Bessel function with argument z, up to (but not including)
 * order n. The function initializer is responsible for populating the 0-order
 * and, when applicable, 1-order values into the vals array. */
int spbessel (complex double *vals, complex double z, int n, besinit initializer) {
	int i, nmo;

	if (n < 0) return -1;

	/* Build the first orders of the appropriate Bessel function. */
	(*initializer)(vals, z, n);

	/* Recursion fills out the rest of the orders. */
	nmo = n - 1;
	for (i = 1; i < nmo; ++i)
		vals[i+1] = (2 * i + 1) * vals[i] / z - vals[i - 1];

	return 0;
}

/* Use the recurrence relation to build regular Bessel functions. */
int spbesh (complex double *vals, complex double z, int n) {
	return spbessel(vals, z, n, spbeshinit);
}

/* Use the recurrence relation to build Hankel functions. */
int spbesj (complex double *vals, complex double z, int n) {
	return spbessel(vals, z, n, spbesjinit);
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

	fact = cexp (I * x) / x;
	return fact * (1 + I / x);
}

/* Returns the derivative of the zero-order spherical Bessel function. */
complex double dspbesj0 (complex double x) {
	complex double sx, cx;

	sx = csin (x) / x;
	cx = ccos (x);

	return (cx - sx) / x;
}

/* Returns the second derivative of the zero-order spherical Hankel function. */
complex double ddspbesh0 (complex double x) {
	complex double fact, xsq, xcu;

	xsq = x * x;
	xcu = xsq * x;

	fact = cexp (I * x);

	return fact * (I * xsq - 2 * (x + I)) / xcu;
}

/* Returns the second derivative of the zero-order spherical Bessel function. */
complex double ddspbesj0 (complex double x) {
	complex double sx, cx, xsq, xcu;

	xsq = x * x;
	xcu = xsq * x;

	sx = csin (x);
	cx = ccos (x);

	return ((2 - xsq) * sx - 2 * x * cx) / xcu;
}

/* This routine uses a common recurrence formula to compute, in the array vals,
 * the first derivative of any spherical Bessel function with argument z, up to
 * order n, given the array hfn of precomputed Bessel functions (of the same
 * kind) up to order n. The function initializer returns the complex value of
 * the derivative of the 0-order function. */
int dspbessel (complex double *vals, complex double *hfn,
		complex double z, int n, cfunc initializer) {
	int i;

	/* The first value is computed using the direct initializer. */
	vals[0] = (*initializer)(z);

	/* Use the recursion formula. */
	for (i = 1; i < n; ++i)
		vals[i] = hfn[i - 1] - (i + 1) * hfn[i] / z;

	return 0;
}

/* Use the recurrence relation to build regular Bessel function derivatives. */
int dspbesh (complex double *vals, complex double *hfn, complex double z, int n) {
	return dspbessel(vals, hfn, z, n, dspbesh0);
}

/* Use the recurrence relation to build Hankel function derivatives. */
int dspbesj (complex double *vals, complex double *hfn, complex double z, int n) {
	return dspbessel(vals, hfn, z, n, dspbesj0);
}

/* This routine uses a common recurrence formula to compute, in the array vals,
 * the second derivative of any spherical Bessel function with argument z, up
 * to order n, given the array hfn of precomputed Bessel functions (of the same
 * kind) up to order n. The function initializer returns the complex value of
 * the second derivative of the 0-order function. */
int ddspbessel (complex double *vals, complex double *hfn,
		complex double z, int n, cfunc initializer) {
	int i;
	complex double zsq = z * z;

	/* The first value is computed using the direct initializer. */
	vals[0] = (*initializer)(z);

	/* Use the recursion formula. */
	for (i = 1; i < n; ++i)
		vals[i] = ((i + 2) * (i + 1) / zsq - 1) * hfn[i] - 2 * hfn[i - 1] / z;

	return 0;
}

/* Use the recurrence relation to build regular Bessel function derivatives. */
int ddspbesh (complex double *vals, complex double *hfn, complex double z, int n) {
	return ddspbessel(vals, hfn, z, n, ddspbesh0);
}

/* Use the recurrence relation to build Hankel function derivatives. */
int ddspbesj (complex double *vals, complex double *hfn, complex double z, int n) {
	return ddspbessel(vals, hfn, z, n, ddspbesj0);
}

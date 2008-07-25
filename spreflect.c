#include <stdlib.h>

#include <complex.h>

#include "spbessel.h"
#include "spreflect.h"

int spbldrc (complex double *vals, complex double k0, complex double k1,
		double rho0, double rho1, double r, int ord) {
	complex double *jl0, *hl0, *jl1, *hl1, *djl0, *dhl0, *djl1, *dhl1,
		k0r, k1r;
	complex double gamma;
	int i;

	/* Set up the arrays for all of the Bessel function values needed. */
	jl0 = malloc (8 * ord * sizeof(complex double));
	hl0 = jl0 + ord;
	jl1 = hl0 + ord;
	hl1 = jl1 + ord;
	djl0 = hl1 + ord;
	dhl0 = djl0 + ord;
	djl1 = dhl0 + ord;
	dhl1 = djl1 + ord;

	k0r = k0 * r;
	k1r = k1 * r;

	gamma = rho0 * k1 / (rho1 * k0);

	/* Compute the four Bessel function values. */
	spbesj (jl0, k0r, ord);
	spbesh (hl0, k0r, ord);
	spbesj (jl1, k1r, ord);
	spbesh (hl1, k1r, ord);

	/* Now compute all of the (whole-argument) derivatives. */
	dspbesj (djl0, jl0, k0r, ord);
	dspbesh (dhl0, hl0, k0r, ord);
	dspbesj (djl1, jl1, k1r, ord);
	dspbesh (dhl1, hl1, k1r, ord);

	/* Now build the reflection coefficients for each order. */
	for (i = 0; i < ord; ++i) {
		k0r = jl1[i] * djl0[i] - gamma * jl0[i] * djl1[i];
		k1r = gamma * hl0[i] * djl1[i] - jl1[i] * dhl0[i];

		vals[i] = k0r / k1r;
	}

	free (jl0);

	return ord;
}

/* This should be able to handle in-place reflection, where inc and scat are
 * pointers to the same data. */
int spreflect (complex double *scat, complex double *inc,
		complex double *spr, int ord, int nphi, int osgn, int asgn) {
	int i, j, off, npj;

	for (i = 0; i < ord; ++i) {
		off = i * nphi;

		/* Handle the zero-order case. */
		scat[off] = osgn * scat[off] + asgn * inc[off] * spr[i];
		
		/* Handle the nonzero orders. */
		for (j = 1; j <= i; ++j) {
			npj = nphi - j;
			scat[off + j] = osgn * scat[off + j] + asgn * inc[off + j] * spr[i];
			scat[off + npj] = osgn * scat[off + npj] + asgn * inc[off + npj] * spr[i];
		}
	}

	return ord;
}

int esbldrc (complex double *reflect, complex double *transmit, complex double k0,
		complex double k1, double rho0, double rho1, double r, int ord) {
	complex double *jl0, *hl0, *jl1, *hl1, *djl0, *dhl0, *djl1, *dhl1,
		*outrfl, *outtr, gamma, k0r, k1r;
	int i;

	outrfl = reflect + ord;
	outtr = transmit + ord;

	/* Set up the arrays for all of the Bessel function values needed. */
	jl0 = malloc (8 * ord * sizeof(complex double));
	hl0 = jl0 + ord;
	jl1 = hl0 + ord;
	hl1 = jl1 + ord;
	djl0 = hl1 + ord;
	dhl0 = djl0 + ord;
	djl1 = dhl0 + ord;
	dhl1 = djl1 + ord;

	k0r = k0 * r;
	k1r = k1 * r;

	gamma = rho0 * k1 / (rho1 * k0);

	/* Compute the four Bessel function values. */
	spbesj (jl0, k0r, ord);
	spbesh (hl0, k0r, ord);
	spbesj (jl1, k1r, ord);
	spbesh (hl1, k1r, ord);

	/* Now compute all of the (whole-argument) derivatives. */
	dspbesj (djl0, jl0, k0r, ord);
	dspbesh (dhl0, hl0, k0r, ord);
	dspbesj (djl1, jl1, k1r, ord);
	dspbesh (dhl1, hl1, k1r, ord);

	/* Now build the reflection coefficients for each order. */
	for (i = 0; i < ord; ++i) {
		k1r = gamma * djl1[i] * hl0[i] - jl1[i] * dhl0[i];

		k0r = djl0[i] * hl0[i] - jl0[i] * dhl0[i];
		transmit[i] = k0r / k1r;

		k0r = gamma * hl0[i] * dhl1[i] - hl1[i] * dhl0[i];
		reflect[i] = k0r / k1r;

		k1r = dhl0[i] * jl0[i] - djl0[i] * hl0[i];

		k0r = gamma * djl1[i] * jl0[i] - djl0[i] * jl1[i];
		outrfl[i] = k0r / k1r;

		k0r = gamma * dhl1[i] * jl0[i] - djl0[i] * hl0[i];
		outtr[i] = k0r / k1r;
	}

	free (jl0);

	return ord;
}

#include <stdlib.h>

#include <complex.h>

#include "spbessel.h"
#include "spreflect.h"

int spbldrc (complex double *vals, complex double k0, complex double k1,
		double rho0, double rho1, double r, int ord) {
	complex double *jl0, *hl0, *jl1, *hl1, *djl0, *dhl0, *djl1, *dhl1,
		k0r, k1r;
	double gamma;
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

	gamma = rho0 * creal (k1) / (rho1 * creal (k0));

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
		k0r = k0 * jl1[i] * djl0[i] - gamma * k1 * jl0[i] * djl1[i];
		k1r = gamma * k1 * hl0[i] * djl1[i] - k0 * jl1[i] * dhl0[i];

		vals[i] = k0r / k1r;
	}

	free (jl0);

	return ord;
}

/* This should be able to handle in-place reflection, where inc and scat are
 * pointers to the same data. */
int spreflect (complex double *scat, complex double *inc,
		complex double *spr, int ord, int nphi) {
	int i, j, off, npj;

	for (i = 0; i < ord; ++i) {
		off = i * nphi;

		/* Handle the zero-order case. */
		scat[off] = inc[off] * spr[i];

		/* Handle the nonzero orders. */
		for (j = 1; j <= i; ++j) {
			npj = nphi - j;
			scat[off + j] = inc[off + j] * spr[i];
			scat[off + npj] = inc[off + npj] * spr[i];
		}
	}

	return ord;
}

int spinvrfl (complex double *scat, complex double *inc,
		complex double *spr, int ord, int nphi) {
	int i, j, off, npj;

	for (i = 0; i < ord; ++i) {
		off = i * nphi;

		/* Handle the zero-order case. */
		scat[off] = inc[off] / spr[i];

		/* Handle the nonzero orders. */
		for (j = 1; j <= i; ++j) {
			npj = nphi - j;
			scat[off + j] = inc[off + j] / spr[i];
			scat[off + npj] = inc[off + npj] / spr[i];
		}
	}

	return ord;
}

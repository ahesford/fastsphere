#include <stdlib.h>

#include <float.h>
#include <complex.h>

#include "spbessel.h"
#include "spreflect.h"
#include "util.h"

int spbldrc (sptype *sphere, complex double k0, double rho0) {
	complex double *jl0, *jl1, *jl2, *hl0, *djl0, *djl1, *djl2, *dhl0;
	complex double *ddjl1, *ddjl2, *reflect;
	complex double f, g, rs, gamma, k0r, k1r, ik2r, k1k2sq, ik2rsq, num, den;
	int i, ord, shear, nt;

	ord = sphere->deg;
	reflect = sphere->reflect;

	jl0 = malloc (10 * ord * sizeof(complex double));
	jl1 = jl0 + ord;
	jl2 = jl1 + ord;
	hl0 = jl2 + ord;
	djl0 = hl0 + ord;
	djl1 = djl0 + ord;
	djl2 = djl1 + ord;
	dhl0 = djl2 + ord;
	ddjl1 = dhl0 + ord;
	ddjl2 = ddjl1 + ord;

	/* Compute the Bessel arguments for compressive terms. */
	k0r = k0 * sphere->r;
	k1r = sphere->k * sphere->r;

	/* Compute the Bessel functions for compressive terms. */
	spbesj (jl0, k0r, ord);
	spbesj (jl1, k1r, ord);
	spbesh (hl0, k0r, ord);
	/* And whole-argument derivatives of compressive Bessel functions. */
	dspbesj (djl0, jl0, k0r, ord);
	dspbesj (djl1, jl1, k1r, ord);
	dspbesh (dhl0, hl0, k0r, ord);

	/* Compute the inverse Bessel argument for shear terms with no loss. */
	ik2r = invwavenum (sphere->csh, 0.) / sphere->r;
	if (cabs(ik2r) > DBL_EPSILON) {
		/* If there is a shear mode, compute the shear Bessel terms. */
		complex double k2r = 1 / ik2r;
		spbesj (jl2, k2r, ord);
		dspbesj (djl2, jl2, k2r, ord);
		ddspbesj (ddjl2, jl2, k2r, ord);
		ddspbesj (ddjl1, jl1, k1r, ord);
	} else {
		/* Without a shear mode, the shear terms are don't cares. */
		ik2r = 0.;
		for (i = 0; i < ord; ++i) {
			jl2[i] = djl2[i] = ddjl1[i] = 0;
			/* This appears in a denominator; set to
			 * unity to avoid indeterminate forms. */
			ddjl2[i] = 1.;
		}
	}

	k1k2sq = (k1r * ik2r) * (k1r * ik2r);
	ik2rsq = ik2r * ik2r;

	gamma = (sphere->rho / rho0) * k0r;

	for (i = 0; i < ord; ++i) {
		nt = i * (i + 1);
		f = (1 - 2 * k1k2sq) * jl1[i] - k1k2sq * ddjl1[i];
		g = ik2rsq * jl2[i] - ik2r * djl2[i];
		rs = ik2rsq * (k1r * djl1[i] - jl1[i]);
		rs /= (ik2rsq * (nt - 2) * jl2[i] + ddjl2[i]);
		num = (f - 4 * nt * rs * g) * gamma * djl0[i];
		num += (2 * nt * rs * jl2[i] - k1r * djl1[i]) * jl0[i];
		den = 2 * nt * rs * (2 * gamma * g * dhl0[i] - jl2[i] * hl0[i]);
		den -= (gamma * f * dhl0[i] - k1r * djl1[i] * hl0[i]);
		reflect[i] = num / den;
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

int esbldrc (sptype *sphere, complex double k0, double rho0) {
	complex double *jl0, *hl0, *jl1, *hl1, *djl0, *dhl0, *djl1, *dhl1;
	complex double *reflect, *transmit, *outrfl, *outtr;
	complex double gamma, k0r, k1r;
	int i, ord;

	ord = sphere->deg;

	reflect = sphere->reflect;
	transmit = sphere->transmit;
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

	k0r = k0 * sphere->r;
	k1r = sphere->k * sphere->r;

	gamma = (rho0 / sphere->rho) * (sphere->k / k0);

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

		k1r = dhl0[i] * jl0[i] - hl0[i] * djl0[i];

		k0r = gamma * djl1[i] * jl0[i] - jl1[i] * djl0[i];
		outrfl[i] = k0r / k1r;

		k0r = gamma * dhl1[i] * jl0[i] - hl1[i] * djl0[i];
		outtr[i] = k0r / k1r;
	}

	free (jl0);

	return ord;
}

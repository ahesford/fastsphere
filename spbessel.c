#include <stdlib.h>

#include <complex.h>
#include <math.h>

#include "spbessel.h"

complex double dspbesj0 (complex double);
complex double dspbesh0 (complex double);

int spbesh (complex double *vals, complex double z, int n) {
	double zr, zi, *cr, *ci, fnu = 0.5;
	int one = 1, nz, ierr;
	complex double fact;

	cr = malloc (2 * n * sizeof(double));
	ci = cr + n;

	zr = creal(z);
	zi = cimag(z);

	/* Call the complex Hankel function routine. */
	zbesh_ (&zr, &zi, &fnu, &one, &one, &n, cr, ci, &nz, &ierr);

	if (ierr) {
		free (cr);
		return ierr;
	}

	/* The spherical scaling factor. */
	fact = csqrt (M_PI_2 / z);

	/* Scale and collapse the real and imaginary parts. */
	for (nz = 0; nz < n; ++nz) 
		vals[nz] = fact * (cr[nz] + I * ci[nz]);

	free (cr);

	return 0;
}

int spbesj (complex double *vals, complex double z, int n) {
	double zr, zi, *cr, *ci, fnu = 0.5;
	int one = 1, nz, ierr;
	complex double fact;

	cr = malloc (2 * n * sizeof(double));
	ci = cr + n;

	zr = creal(z);
	zi = cimag(z);

	/* Call the complex Bessel function routine. */
	zbesj_ (&zr, &zi, &fnu, &one, &n, cr, ci, &nz, &ierr);

	if (ierr) {
		free (cr);
		return ierr;
	}

	/* The spherical scaling factor. */
	fact = csqrt (M_PI_2 / z);

	/* Scale and collapse the real and imaginary parts. */
	for (nz = 0; nz < n; ++nz) 
		vals[nz] = fact * (cr[nz] + I * ci[nz]);

	free (cr);

	return 0;
}

/* Returns the derivative of the zero-order spherical Hankel function. */
complex double dspbesh0 (complex double x) {
	complex double fact;

	fact = (x + I) / (x * x);

	return fact * cexp (1i * x);
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

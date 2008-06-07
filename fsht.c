#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <complex.h>
#include <math.h>

#include <fftw3.h>

#include <gsl/gsl_sf_legendre.h>

#include "fsht.h"

void gaqd_ (int *, double *, double *, double *, int *, int *);

int fshtinit (shdata *dat, int deg) {
	int ierr;

	/* The number of theta samples isn't always the same as the number of
	 * spherical harmonic degrees used. */
	dat->ntheta = deg;

	/* This is actually one more than the maximum degree. */
	dat->deg = deg;

	/* This many phi values are required for fast evaluation with FFT. */
	dat->nphi = 2 * deg - 1;

	/* Allocate the theta points and weights. */
	dat->theta = malloc (2 * dat->ntheta * sizeof(double));
	dat->weights = dat->theta + dat->ntheta;

	/* Allocate the Legendre polynomial values. */
	dat->lgvals = malloc (dat->deg * sizeof(double));

	/* Allocate the FFT data buffer. */
	dat->fftbuf = fftw_malloc (dat->ntheta * dat->nphi * sizeof(complex double));

	dat->fplan = fftw_plan_many_dft (1, &(dat->nphi), dat->ntheta, dat->fftbuf,
			&(dat->nphi), 1, dat->nphi, dat->fftbuf, &(dat->nphi),
			1, dat->nphi, FFTW_FORWARD, FFTW_MEASURE);

	dat->bplan = fftw_plan_many_dft (1, &(dat->nphi), dat->ntheta, dat->fftbuf,
			&(dat->nphi), 1, dat->nphi, dat->fftbuf, &(dat->nphi),
			1, dat->nphi, FFTW_BACKWARD, FFTW_MEASURE);

	/* Find the Legendre-Gauss quadrature points. */
	gaqd_ (&(dat->ntheta), dat->theta, dat->weights, NULL, NULL, &ierr);
	if (ierr) return ierr;

	return 0;
}

int fshtfree (shdata *dat) {
	if (dat->theta) free (dat->theta);
	if (dat->lgvals) free (dat->lgvals);
	if (dat->fftbuf) fftw_free (dat->fftbuf);

	fftw_destroy_plan (dat->fplan);
	fftw_destroy_plan (dat->bplan);

	dat->ntheta = dat->nphi = dat->deg = 0;

	return 0;
}

/* Scale the components of the spherical harmonic coefficients to relate the
 * far-field signature to the spherical harmonic expansion. The default (for
 * non-negative sgn) properly scales SH coefficients AFTER a forward
 * transform (angular to spherical). If sgn is negative, it scales the SH
 * coefficients BEFORE an inverse transform (spherical to angular). */
int shscale (complex double *samp, shdata *dat, int sgn) {
	complex double cscale[4] = { I, -1.0, -I, 1.0 };
	int i, j, off, idx;

	/* For negative sgn, flip the signs of the imaginary multipliers. */
	if (sgn < 0) {
		cscale[0] = -I;
		cscale[2] = I;
	}

	for (i = 0; i < dat->deg; ++i) {
		idx = i % 4;
		off = i * dat->nphi;

		/* Scale the zero order for all degrees. */
		samp[off] *= cscale[idx];

		/* Scale the nonzero orders for all degrees. */
		for (j = 1; j <= i; ++j) {
			samp[off + j] *= cscale[idx];
			samp[off + dat->nphi - j] *= cscale[idx];
		}
	}

	return 0;
}

/* Forward spherical harmonic transform: take samples of the function in theta
 * and phi into SH coefficients. */
int ffsht (complex double *samp, shdata *dat) {
	int i, j, k, aoff, npk, dm1;
	double cth, pc, scale;
	complex double *beta;

	/* Copy the samples to the FFT buffer and perform the FFT. */
	memcpy (dat->fftbuf, samp, dat->ntheta * dat->nphi * sizeof(complex double));
	fftw_execute (dat->fplan);

	/* Zero out the input samples, to prepare for storage of coefficients. */
	memset (samp, 0, dat->ntheta * dat->nphi * sizeof(complex double));

	beta = dat->fftbuf;

	dm1 = dat->deg - 1;

	for (i = 0; i < dat->ntheta; ++i) {
		/* Some scale factors that will be required. */
		scale = 2 * M_PI * dat->weights[i] / dat->nphi;
		cth = cos(dat->theta[i]);

		/* Build the Legendre polynomials that we need. */
		gsl_sf_legendre_sphPlm_array (dm1, 0, cth, dat->lgvals);

		/* Handle m = 0 for all degrees. */
		for (j = 0; j < dat->deg; ++j)
			samp[j * dat->nphi] += scale * beta[0] * dat->lgvals[j];

		/* Handle nonzero orders for all relevant degrees. */
		for (k = 1; k < dat->deg; ++k) {
			npk = dat->nphi - k;
			gsl_sf_legendre_sphPlm_array (dm1, k, cth, dat->lgvals);
			for (j = k; j < dat->deg; ++j) {
				aoff = j * dat->nphi;
				pc = scale * dat->lgvals[j - k];
				/* The positive-order coefficient. */
				samp[aoff + k] += pc * beta[k];
				/* The negative-order coefficient. */
				samp[aoff + npk] += pc * beta[npk];
			}
		}

		beta += dat->nphi;
	}

	return dat->deg;
}

/* Inverse spherical harmonic transform: take SH coefficients to sample of the
 * function in theta and phi. */
int ifsht (complex double *samp, shdata *dat) {
	int i, j, k, aoff, npk, dm1;
	double cth;
	complex double *beta;

	/* Zero out the FFT buffer to prepare for evaluation of function. */
	memset (dat->fftbuf, 0, dat->ntheta * dat->nphi * sizeof(complex double));

	beta = dat->fftbuf;

	dm1 = dat->deg - 1;

	for (i = 0; i < dat->ntheta; ++i) {
		/* Some factors that will be used several times. */
		cth = cos(dat->theta[i]);

		/* Build the Legendre polynomials that we need. */
		gsl_sf_legendre_sphPlm_array (dm1, 0, cth, dat->lgvals);

		for (j = 0; j < dat->deg; ++j)
			beta[0] += samp[j * dat->nphi] * dat->lgvals[j];

		for (k = 1; k < dat->deg; ++k) {
			npk = dat->nphi - k;
			gsl_sf_legendre_sphPlm_array (dm1, k, cth, dat->lgvals);
			for (j = k; j < dat->deg; ++j) {
				aoff = j * dat->nphi;
				/* The positive-order coefficient. */
				beta[k] += samp[aoff + k] * dat->lgvals[j - k];
				/* The negative-order coefficient. */
				beta[npk] += samp[aoff + npk] * dat->lgvals[j - k];
			}
		}

		/* Move the next theta value in the FFT array. */
		beta += dat->nphi;
	}

	/* Perform the inverse FFT and copy the values to the storage area. */
	fftw_execute (dat->bplan);
	memcpy (samp, dat->fftbuf, dat->ntheta * dat->nphi * sizeof(complex double));

	return dat->deg;
}

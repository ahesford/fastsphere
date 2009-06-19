#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <complex.h>
#include <math.h>

#include <fftw3.h>

#include <gsl/gsl_sf_legendre.h>

#include "fsht.h"
#include "util.h"

void gaqd_ (int *, double *, double *, double *, int *, int *);

int fshtinit (shdata *dat, int deg, int ntheta, int nphi) {
	int ierr;
	complex double *fftbuf;

	/* This is actually one more than the maximum degree. */
	dat->deg = deg;

	/* The number of theta samples must be at least double the degree. */
	dat->ntheta = MAX(ntheta, 2 * deg + 1);
	/* The number of phi samples should be at least twice the theta samples. */
	dat->nphi = MAX(2 * dat->ntheta, nphi);

	/* Temporarily allocate an FFT data buffer for planning. */
	fftbuf = fftw_malloc (dat->ntheta * dat->nphi * sizeof(complex double));

	/* Plan the forward and inverse transforms. */
	dat->fplan = fftw_plan_many_dft (1, &(dat->nphi), dat->ntheta, fftbuf,
			&(dat->nphi), 1, dat->nphi, fftbuf, &(dat->nphi), 1,
			dat->nphi, FFTW_FORWARD, FFTW_MEASURE);

	dat->bplan = fftw_plan_many_dft (1, &(dat->nphi), dat->ntheta, fftbuf,
			&(dat->nphi), 1, dat->nphi, fftbuf, &(dat->nphi), 1,
			dat->nphi, FFTW_BACKWARD, FFTW_MEASURE);

	/* The FFT buffer is no longer necessary, and will be reallocated
	 * on-the-fly when it is needed later. */
	fftw_free (fftbuf);

	return 0;
}

int fshtfree (shdata *dat) {
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
	double sth, cth, pc, scale, *lgvals, theta, dtheta;
	complex double *beta, *fftbuf;

	/* Copy the samples to the FFT buffer and perform the FFT. */
	fftbuf = fftw_malloc (dat->ntheta * dat->nphi * sizeof(complex double));
	memcpy (fftbuf, samp, dat->ntheta * dat->nphi * sizeof(complex double));

	lgvals = malloc (dat->deg * sizeof(double));

	/* Perform an in-place FFT of the buffer. */
	fftw_execute_dft (dat->fplan, fftbuf, fftbuf);

	/* Zero out the input samples, to prepare for storage of coefficients. */
	memset (samp, 0, dat->ntheta * dat->nphi * sizeof(complex double));

	beta = fftbuf;

	dm1 = dat->deg - 1;

	/* The integration scale factor. */
	scale = 2 * M_PI * M_PI / (dat->nphi * dat->ntheta);
	/* The step in theta. */
	dtheta = M_PI / dat->ntheta;

	for (i = 0, theta = dtheta / 2; i < dat->ntheta; ++i, theta += dtheta) {
		/* Some factors that will be required. */
		cth = cos (theta);
		sth = sin (theta);

		/* Build the Legendre polynomials that we need. */
		gsl_sf_legendre_sphPlm_array (dm1, 0, cth, lgvals);

		/* Handle m = 0 for all degrees. */
		for (j = 0; j < dat->deg; ++j)
			samp[j * dat->nphi] += scale * sth * beta[0] * lgvals[j];

		/* Handle nonzero orders for all relevant degrees. */
		for (k = 1; k < dat->deg; ++k) {
			npk = dat->nphi - k;
			gsl_sf_legendre_sphPlm_array (dm1, k, cth, lgvals);
			for (j = k; j < dat->deg; ++j) {
				aoff = j * dat->nphi;
				pc = scale * sth * lgvals[j - k];
				/* The positive-order coefficient. */
				samp[aoff + k] += pc * beta[k];
				/* The negative-order coefficient. */
				samp[aoff + npk] += pc * beta[npk];
			}
		}

		beta += dat->nphi;
	}

	/* The FFT buffer is no longer required. */
	fftw_free (fftbuf);

	return dat->deg;
}

/* Inverse spherical harmonic transform: take SH coefficients to sample of the
 * function in theta and phi. */
int ifsht (complex double *samp, shdata *dat) {
	int i, j, k, aoff, npk, dm1, n;
	double cth, *lgvals, theta, dtheta;
	complex double *beta, *fftbuf, polval[2];

	/* The non-polar samples. */
	n = dat->ntheta * dat->nphi;

	/* Zero out the FFT buffer to prepare for evaluation of function. */
	fftbuf = fftw_malloc (n * sizeof(complex double));
	memset (fftbuf, 0, n * sizeof(complex double));

	lgvals = malloc (dat->deg * sizeof(double));

	beta = fftbuf;

	dm1 = dat->deg - 1;

	dtheta = M_PI / dat->ntheta;

	for (i = 0, theta = dtheta / 2; i < dat->ntheta; ++i, theta += dtheta) {
		/* Some factors that will be used several times. */
		cth = cos(theta);

		/* Build the Legendre polynomials that we need. */
		gsl_sf_legendre_sphPlm_array (dm1, 0, cth, lgvals);

		for (j = 0; j < dat->deg; ++j)
			beta[0] += samp[j * dat->nphi] * lgvals[j];

		for (k = 1; k < dat->deg; ++k) {
			npk = dat->nphi - k;
			gsl_sf_legendre_sphPlm_array (dm1, k, cth, lgvals);
			for (j = k; j < dat->deg; ++j) {
				aoff = j * dat->nphi;
				/* The positive-order coefficient. */
				beta[k] += samp[aoff + k] * lgvals[j - k];
				/* The negative-order coefficient. */
				beta[npk] += samp[aoff + npk] * lgvals[j - k];
			}
		}

		/* Move the next theta value in the FFT array. */
		beta += dat->nphi;
	}

	polval[0] = polval[1] = 0;

	/* Handle the poles. Note that the zero-order Legendre functions are
	 * polynomials, so the evenness and oddness can be exploited. */
	gsl_sf_legendre_sphPlm_array (dm1, 0, 1, lgvals);

	for (j = 0; j < dat->deg; ++j) {
		/* Sign change for the south pole. */
		aoff = 1 - 2 * (j % 2);
		polval[0] += samp[j * dat->nphi] * lgvals[j];
		polval[1] += aoff * samp[j * dat->nphi] * lgvals[j];
	}

	/* Perform the inverse FFT and copy the values to the storage area. */
	fftw_execute_dft (dat->bplan, fftbuf, fftbuf);
	memcpy (samp, fftbuf, n * sizeof(complex double));

	/* Copy the polar values. */
	samp[n] = polval[0];
	samp[n + 1] = polval[1];

	/* Eliminate the FFT buffer. */
	fftw_free (fftbuf);
	free (lgvals);

	return dat->deg;
}

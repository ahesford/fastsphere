#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <complex.h>
#include <math.h>

#include <fftw3.h>

#include <gsl/gsl_sf_legendre.h>

#include "fsht.h"
#include "util.h"

int fshtinit (shdata *dat, int deg, int ntheta, int nphi, int measfft) {
	int ierr;
	complex double *fftbuf;

	/* This is actually one more than the maximum degree. */
	dat->deg = deg;

	/* The number of theta samples must at least equal the SH degree. */
	if (ntheta < deg) dat->ntheta = deg;
	else dat->ntheta = ntheta;

	/* At least this many phi values are required for fast FFT evaluation. */
	if (2 * deg - 1 > nphi) dat->nphi = 2 * deg - 1;
	else dat->nphi = nphi;

	/* Allocate the theta points and weights. */
	dat->theta = calloc (2 * dat->ntheta, sizeof(double));
	dat->weights = dat->theta + dat->ntheta;

	/* Find the Legendre-Gauss quadrature points. */
	gauleg (dat->ntheta, dat->theta, dat->weights);

	/* Temporarily allocate an FFT data buffer for planning. */
	fftbuf = fftw_malloc (dat->ntheta * dat->nphi * sizeof(complex double));

	/* Plan the forward and inverse transforms. */
	dat->fplan = fftw_plan_many_dft (1, &(dat->nphi), dat->ntheta, fftbuf,
			&(dat->nphi), 1, dat->nphi, fftbuf, &(dat->nphi), 1,
			dat->nphi, FFTW_FORWARD, measfft ? FFTW_MEASURE : FFTW_ESTIMATE);

	dat->bplan = fftw_plan_many_dft (1, &(dat->nphi), dat->ntheta, fftbuf,
			&(dat->nphi), 1, dat->nphi, fftbuf, &(dat->nphi), 1,
			dat->nphi, FFTW_BACKWARD, measfft ? FFTW_MEASURE : FFTW_ESTIMATE);

	/* The FFT buffer is no longer necessary, and will be reallocated
	 * on-the-fly when it is needed later. */
	fftw_free (fftbuf);

	return 0;
}

int fshtfree (shdata *dat) {
	if (dat->theta) free (dat->theta);

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
int ffsht (complex double *samp, shdata *dat, int maxdeg) {
	int i, j, k, aoff, npk, dm1, deg = maxdeg;
	long lgn, lgi;
	double pc, scale, *lgvals;
	complex double *beta, *fftbuf;

	/* Copy the samples to the FFT buffer and perform the FFT. */
	fftbuf = fftw_malloc (dat->ntheta * dat->nphi * sizeof(complex double));
	memcpy (fftbuf, samp, dat->ntheta * dat->nphi * sizeof(complex double));

	/* Perform an in-place FFT of the buffer. */
	fftw_execute_dft (dat->fplan, fftbuf, fftbuf);

	/* Zero out the input samples, to prepare for storage of coefficients. */
	memset (samp, 0, dat->ntheta * dat->nphi * sizeof(complex double));

	beta = fftbuf;

	/* If no maximum degree is specified, use the default. */
	if (maxdeg < 1) deg = dat->deg;

	dm1 = deg - 1;

	/* Create storage for all Legendre polynomials. */
	lgn = gsl_sf_legendre_array_n (dm1);
	lgvals = malloc (lgn * sizeof(double));


	for (i = 0; i < dat->ntheta; ++i) {
		/* The scale factor that will be required. */
		scale = 2 * M_PI * dat->weights[i] / dat->nphi;

		/* Build the Legendre polynomials that we need.
		 * Don't use Condon-Shortley phase factor. */
		gsl_sf_legendre_array_e (GSL_SF_LEGENDRE_SPHARM, dm1, dat->theta[i], 1., lgvals);

		/* Handle m = 0 for all degrees. */
		for (j = 0; j < deg; ++j) {
			lgi = gsl_sf_legendre_array_index(j, 0);
			samp[j * dat->nphi] += scale * beta[0] * lgvals[lgi];
		}

		/* Handle nonzero orders for all relevant degrees. */
		for (k = 1; k < deg; ++k) {
			npk = dat->nphi - k;
			for (j = k; j < deg; ++j) {
				aoff = j * dat->nphi;
				lgi = gsl_sf_legendre_array_index(j, k);
				pc = scale * lgvals[lgi];
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

	return deg;
}

/* Inverse spherical harmonic transform: take SH coefficients to sample of the
 * function in theta and phi. */
int ifsht (complex double *samp, shdata *dat, int maxdeg) {
	int i, j, k, aoff, npk, dm1, n, deg = maxdeg;
	long lgn, lgi;
	double *lgvals;
	complex double *beta, *fftbuf;

	/* The non-polar samples. */
	n = dat->ntheta * dat->nphi;

	/* Zero out the FFT buffer to prepare for evaluation of function. */
	fftbuf = fftw_malloc (n * sizeof(complex double));
	memset (fftbuf, 0, n * sizeof(complex double));

	beta = fftbuf;

	/* Use the default degree if no maximum is specified. */
	if (maxdeg < 1) deg = dat->deg;

	dm1 = deg - 1;

	/* Create storage for all Legendre polynomials. */
	lgn = gsl_sf_legendre_array_n(dm1);
	lgvals = malloc (lgn * sizeof(double));

	for (i = 0; i < dat->ntheta; ++i) {
		/* Build the Legendre polynomials that we need. */
		/* Don't use the Condon-Shortley phase factor. */
		gsl_sf_legendre_array_e (GSL_SF_LEGENDRE_SPHARM, dm1, dat->theta[i], 1., lgvals);

		for (j = 0; j < deg; ++j) {
			lgi = gsl_sf_legendre_array_index(j, 0);
			beta[0] += samp[j * dat->nphi] * lgvals[lgi];
		}

		for (k = 1; k < deg; ++k) {
			npk = dat->nphi - k;
			for (j = k; j < deg; ++j) {
				aoff = j * dat->nphi;
				lgi = gsl_sf_legendre_array_index(j, k);
				/* The positive-order coefficient. */
				beta[k] += samp[aoff + k] * lgvals[lgi];
				/* The negative-order coefficient. */
				beta[npk] += samp[aoff + npk] * lgvals[lgi];
			}
		}

		/* Move the next theta value in the FFT array. */
		beta += dat->nphi;
	}

	/* Perform the inverse FFT and copy the values to the storage area. */
	fftw_execute_dft (dat->bplan, fftbuf, fftbuf);
	memcpy (samp, fftbuf, n * sizeof(complex double));

	/* Eliminate the FFT buffer. */
	fftw_free (fftbuf);
	free (lgvals);

	return deg;
}

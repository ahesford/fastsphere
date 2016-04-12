#include <stdlib.h>
#include <string.h>

#include <math.h>
#include <complex.h>
#include <float.h>

#ifdef _MACOSX
#include <Accelerate/Accelerate.h>
#else
#include <cblas.h>
#endif

#include <gsl/gsl_sf_legendre.h>

#include "util.h"
#include "spbessel.h"

/* Use the modified Gram-Schmidt process to compute (in place) the portion of
 * the n-dimensional vector v orthogonal to each of the nv vectors s. The
 * projection * of the vector onto each of the basis vectors is stored in the
 * length-nv array c. */
int cmgs (complex double *v, complex double *c, complex double *s, int n, int nv) {
	int i, j, k;
	complex double *sv, cv;

	for (i = 0, sv = s; i < nv; ++i, sv += n) {
		c[i] = 0;
		k = 0;

		do {
			cv = pardot (sv, v, n);
			c[i] += cv;
#pragma omp parallel for default(shared) private(j)
			for (j = 0; j < n; ++j) v[j] -= cv * sv[j];
		} while (cabs(cv / c[i]) > IMGS_TOL && ++k < IMGS_ITS);
		
	}

	return n;
}

complex double pardot (complex double *x, complex double *y, int n) {
	complex double dp;

	/* Compute the local portion. */
	cblas_zdotc_sub (n, x, 1, y, 1, &dp);

	return dp;
}

/* The MSE between two vectors. */
double rmserror (complex double *v, complex double *r, int n) {
	double err = 0, errd = 0, e;
	int i;

	for (i = 0; i < n; ++i) {
		e = cabs (v[i] - r[i]);
		err += e * e;
		e = cabs (r[i]);
		errd += e * e;
	}

	return sqrt (err / errd);
}

/* Open a file or die. */
FILE *critopen (char *fname, char *mode) {
	FILE *fptr;

	fptr = fopen (fname, mode);

	if (!fptr) {
		fprintf (stderr, "ERROR: Could not open input file %s.\n", fname);
		exit (EXIT_FAILURE);

		return NULL;
	}

	return fptr;
}

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
complex double wavenum (double cr, double alpha) {
	double kr, ki;

	ki = log(10) * alpha / 20;
	kr = 2 * M_PI / cr;

	return kr + I * ki;
}

/* Computes the inverse wave number from the relative sound speed and
 * unitless (dB) attenuation coefficient. */
complex double invwavenum (double cr, double alpha) {
	double kr, ki, mag, csq;

	csq = cr * cr;

	kr = 2 * M_PI * cr;
	ki = log(10) * alpha / 20;
	mag = 4 * M_PI * M_PI + ki * ki * csq;

	return (kr / mag) - I * (ki * csq / mag);
}

/* Copy the spherical harmonic representation from one location to another. */
int copysh (int deg, complex double *out, int ldo, complex double *in, int ldi) {
	int i, j, offo, offi;

	/* Copy the spherical harmonic coefficients into the right place. */
	for (i = 0; i < deg; ++i) {
		offo = i * ldo;
		offi = i * ldi;
		/* Copy the zero-order coefficients. */
		out[offo] = in[offi];

		/* Copy the other coefficients. */
		for (j = 1; j <= i; ++j) {
			out[offo + j] = in[offi + j];
			out[offo + ldo - j] = in[offi + ldi - j];
		}
	}

	return deg;
}

/* Multiply in a radial component to the spherical harmonic. */
int shradial (int deg, complex double *a, int lda, complex double k, double r) {
	complex double kr, *hlkr;
	int i, j, off;

	hlkr = malloc (deg * sizeof(complex double));

	/* Compute the Hankel functions for all degrees. */
	kr = k * r;
	spbesh (hlkr, kr, deg);

	for (i = 0; i < deg; ++i) {
		off = i * lda;

		/* The zero-order coefficients. */
		a[off] *= hlkr[i];

		/* The other coefficients. */
		for (j = 1; j <= i; ++j) {
			a[off + j] *= hlkr[i];
			a[off + lda - j] *= hlkr[i];
		}
	}

	free (hlkr);

	return deg;
}

/* Compute the harmonic coefficients for an incident plane wave. */
int shincident (int deg, complex double *a, int lda,
		complex double mag, double theta, double phi) {
	double cth, *lgvals;
	int l, m, dm1, npm, off;
	long lgi, lgn;
	complex double scale[4] = { 1, -I, -1, I }, cx, fx;

	/* The coefficient has a 4 pi factor that must be included. */
	mag *= 4 * M_PI;

	dm1 = deg - 1;
	cth = cos(theta);

	lgn = gsl_sf_legendre_array_n (dm1);
	lgvals = malloc (lgn * sizeof(double));

	/* Compute the Legendre polynomials of zero order for all degrees. */
	gsl_sf_legendre_array_e (GSL_SF_LEGENDRE_SPHARM, dm1, cth, 1., lgvals);

	/* Compute the zero-order coefficients for all degrees. */
	for (l = 0; l < deg; ++l) {
		lgi = gsl_sf_legendre_array_index(l, 0);
		a[l * lda] += scale[l % 4] * mag * lgvals[lgi];
	}

	for (m = 1; m < deg; ++m) {
		npm = lda - m;

		/* Compute the phi variation. */
		cx = cexp (I * m * phi);

#pragma omp critical(incaug)
		for (l = m; l < deg; ++l) {
			off = l * lda;
			lgi = gsl_sf_legendre_array_index(l, m);
			fx = scale[l % 4] * mag * lgvals[lgi];
			/* The positive-order coefficient. */
			a[off + m] += fx * conj (cx);
			/* The negative-order coefficient. */
			a[off + npm] += fx * cx;
		}
	}

	free (lgvals);

	return deg;
}

static int lgval (double *p, double *dp, double t, int m) {
	double p0 = 1.0, p1 = t;
	int k;

	/* Set values explicitly for low orders. */
	if (m < 1) {
		*p = 1.0;
		*dp = 0.0;
		return m;
	} else if (m < 2) {
		*p = t;
		*dp = 1.0;
		return m;
	}

	/* Otherwise compute the values explicitly. */
	for (k = 1; k < m; ++k) {
		*p = ((2.0 * k + 1.0) * t * p1 - k * p0) / (1.0 + k);
		p0 = p1; p1 = *p;
	}

	*dp = m * (p0 - t * p1) / (1.0 - t * t);

	return m;
}

int gauleg (int m, double *nodes, double *weights) {
	int i, j, nroots = (m + 1) / 2;
	double t, p, dp, dt;
	const int maxit = 100;
	const double tol = 1.0e-14;

#pragma omp parallel for default(shared) private(i,j,t,p,dp,dt)
	for (i = 0; i < nroots; ++i) {
		t = cos (M_PI * (i + 0.75) / (m + 0.5));

		for (j = 0; j < maxit; ++j) {
			/* Compute the value of the Legendre polynomial. */
			lgval (&p, &dp, t, m);

			/* Perform a Newton-Raphson update. */
			dt = -p / dp;
			t += dt;

			/* Break if convergence has been achieved. */
			if (fabs(dt) < tol) break;
		}

		/* Update the nodes and weights. */
		nodes[i] = t;
		nodes[m - i - 1] = -t;

		weights[i] = 2.0 / (1.0 - t * t) / (dp * dp);
		weights[m - i - 1] = weights[i];
	}

	return 0;
}

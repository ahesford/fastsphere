#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <math.h>
#include <float.h>

#include <gsl/gsl_sf_legendre.h>

#include "translator.h"
#include "shrotate.h"
#include "util.h"

/* Find the polar and rotation angles of the new z-axis for translation. */
int getangles (double *theta, double *chi, double axis[3]) {
	double st;

	*theta = acos (axis[2]); /* The polar angle of the new z-axis. */

	st = sqrt (axis[0] * axis[0] + axis[1] * axis[1]);

	/* If theta is zero or pi, the rotation angles vanish. */
	if (st < DBL_MIN) {
		*chi = 0;
		return 0;
	}

	/* The rotation angle of the new z-axis. */
	*chi = acos (axis[0] / st);

	/* Rewrap the angle if it's in the lower half plane. */
	if (axis[1] < 0) *chi = 2 * M_PI - *chi;

	return 1;
}

/* Build the initial values of the Hvn function for rotation. */
int buildhvn (double theta, double *hvn, int nmax, int mmax) {
	int i, j, dm1, lda;
	double fact, *buf, arg;

	dm1 = nmax - 1;
	lda = 2 * mmax - 1;

	arg = cos (theta);

	buf = malloc (nmax * sizeof(double));

	/* Build the zero-order array. */
	gsl_sf_legendre_sphPlm_array (dm1, 0, arg, buf);

	for (j = 0; j < nmax; ++j) {
		fact = sqrt(4 * M_PI / (2 * j + 1));
		hvn[ELT(0,j,lda)] = fact * buf[j]; 
	}

	for (i = 1; i < mmax; ++i) {
		/* Build the other orders. */
		gsl_sf_legendre_sphPlm_array (dm1, i, arg, buf);

		for (j = i; j < nmax; ++j) {
			fact = sqrt (4 * M_PI / (2 * j + 1));
			/* The positive order for this degree. */
			hvn[ELT(i,j,lda)] = fact * buf[j - i]; 
			/* The negative order for this degree. */
			hvn[ELT(-i,j+1,lda)] = fact * buf[j - i]; 
		}
	}

	free (buf);

	return nmax;
}

int nexthvn (double theta, double *hvn, int m, int nmax, int lda) {
	int i, j, ilim;
	double st, omct, opct, bnm, am, bm, bp;

	st = sin(theta);
	omct = 1.0 - cos(theta);
	opct = 1.0 + cos(theta);
	
	ilim = nmax - m;

	/* Perform the recursion. */
	for (i = m + 2; i < ilim; ++i) {
		bnm = (double)((i - m) * (i - m - 1));
		bp = -0.5 * sqrt(i * (i + 1) / bnm);
		am = sqrt(i * i / bnm);


		/* The zero-order coefficient. */
		hvn[ELT(0,i-1,lda)] = am * st * hvn[ELT(0,i,lda)]
			+ bp * (opct * hvn[ELT(-1,i+1,lda)] - omct * hvn[ELT(1,i,lda)]);

		/* The non-zero orders. */
		for (j = 1; j < i; ++j) {
			bp = 0.5 * sqrt((i - j) * (i - j + 1) / bnm);
			bm = -0.5 * sqrt((i + j) * (i + j + 1) / bnm);
			am = sqrt((i + j) * (i - j) / bnm);

			/* Positive j. */
			hvn[ELT(j,i-1,lda)] = am * st * hvn[ELT(j,i,lda)]
				+ bp * opct * hvn[ELT(j-1,i,lda)]
				- bm * omct * hvn[ELT(j+1,i,lda)];

			/* Negative j. The special IDX macro must be used
			 * because the column -j+1 can be negative when j = 1.
			 * Otherwise, the faster ELT macro can be used. */
			hvn[ELT(-j,i,lda)] = am * st * hvn[ELT(-j,i+1,lda)]
				+ bm * opct * hvn[ELT(-j-1,i+1,lda)]
				- bp * omct * hvn[IDX(-j+1,i,lda)];
		}
	}

	return nmax * lda;
}

/* Rotate the SH coefficients according to rotation angles throt and chrot. */
int shrotate (complex double *vin, int deg, int lda, trdesc *trans) {
	double *hvn, theta, chi;
	complex double pfz, mfz, *avp, *avm, *buf, iscale[4] = { 1, -I, -1, I };
	int nmax, i, j, m, ldb;

	nmax = 2 * deg - 1;
	ldb = 2 * nmax - 1;

	hvn = malloc (nmax * ldb * sizeof(double));
	buf = malloc (deg * nmax * sizeof(complex double));

	/* Find the rotation angles of the new axis. If the rotation angle
	 * vanishes, adjust the scaling parameters accordingly. */
	if (!getangles (&theta, &chi, trans->sdir))
		iscale[0] = iscale[1] = iscale[2] = iscale[3] = 1.0;

	/* Build the initial values of the H translation function. */
	buildhvn (theta, hvn, nmax, nmax);

	/* Handle the m = 0 case specifically. */
	for (i = 0; i < deg; ++i) {
		avp = vin + ELT(i,0,lda);
		buf[ELT(0,i,nmax)] = (*avp) * hvn[ELT(0,i,ldb)];
		
		/* Handle the positive and negative orders for these degrees. */
		for (j = 1; j <= i; ++j) {
			buf[ELT(j,i,nmax)] = (*avp) * hvn[ELT(j,i,ldb)];
			buf[ELT(-j,i+1,nmax)] = (*avp) * hvn[ELT(-j,i+1,ldb)];
		}
	}

	/* Handle the higher-order cases. */
	for (m = 1; m < deg; ++m) {
		/* Calculate the Hvn samples for the next m. */
		nexthvn (theta, hvn, m - 1, nmax, ldb);

		/* Calculate the phase terms. */
		pfz = cexp (I * m * chi);
		mfz = cexp (-I * m * chi);
		
		for (i = m; i < deg; ++i) {
			avp = vin + ELT(m,i,lda);
			avm = vin + ELT(-m,i+1,lda);

			buf[ELT(0,i,nmax)] += hvn[ELT(0,i,ldb)] * (pfz * (*avp) + mfz * (*avm));
			
			/* Handle the positive and negative orders for these degrees. */
			for (j = 1; j <= i; ++j) {
				buf[ELT(j,i,nmax)] += pfz * (*avp) * hvn[ELT(j,i,ldb)]
					+ mfz * (*avm) * hvn[ELT(-j,i+1,ldb)];
				buf[ELT(-j,i+1,nmax)] += pfz * (*avp) * hvn[ELT(-j,i+1,ldb)]
					+ mfz * (*avm) * hvn[ELT(j,i,ldb)];
			}
		}
	}

	/* Copy the result back to the input buffer. */
	for (i = 0; i < deg; ++i) {
		vin[ELT(0,i,lda)] = buf[ELT(0,i,nmax)];

		/* Have to scale by exp(i * j * phi), where j is the order. */
		for (j = 1; j <= i; ++j) {
			vin[ELT(j,i,lda)] = iscale[j % 4] * buf[ELT(j,i,nmax)];
			vin[ELT(-j,i+1,lda)] = conj (iscale[j % 4]) * buf[ELT(-j,i+1,nmax)];
		}
	}

	free (hvn);
	free (buf);

	return deg;
}

int main (int argc, char **argv) {
	int nmax, n, m, lda, i;
	complex double *coeff;
	trdesc trans;
	double r;

	nmax = atoi (argv[1]);
	n = atoi (argv[2]);
	m = atoi (argv[3]);

	trans.sdir[0] = strtod (argv[4], NULL);
	trans.sdir[1] = strtod (argv[5], NULL);
	trans.sdir[2] = strtod (argv[6], NULL);

	r = sqrt (DVDOT (trans.sdir, trans.sdir));
	trans.sdir[0] /= r; trans.sdir[1] /= r; trans.sdir[2] /= r;

	lda = 2 * nmax - 1;
	
	coeff = calloc (nmax * lda, sizeof(complex double));
	coeff[IDX(m,n,lda)] = 1.0;

	shrotate (coeff, nmax, lda, &trans);
	
	n = lda * nmax;
	for (i = 0; i < n; ++i) {
		if (!(i % lda)) printf ("Degree %d\n", i / lda);
		printf ("%20.15g %20.15g\n", creal(coeff[i]), cimag(coeff[i]));
	}

	free (coeff);

	return EXIT_SUCCESS;
}

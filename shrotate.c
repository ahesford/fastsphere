#include <stdlib.h>
#include <string.h>

#include <math.h>
#include <float.h>

#include <gsl/gsl_sf_legendre.h>

#include "translator.h"
#include "util.h"

#define BPNM(n,m) sqrt(((n)-(m)-1)*((n)-(m))/((2*(n)-1)*(2*(n)+1)))
#define ANM(n,m) sqrt(((n)+1+(m))*((n)+1-(m))/((2*(n)+1)*(2*(n)+3)))
#define BNM(n,m) SGN(m)*BPNM(n,m)

#define IDX(n,m,lda) ((m) < 0 ? ((n)+1) * (lda) - (m) : (n) * (lda) + (m))

/* Find the polar and rotation angles of the new z-axis for translation. */
int getangles (double *theta, double *chi, double axis[3]) {
	double st;

	*theta = acos (axis[2]); /* The polar angle of the new z-axis. */

	st = sqrt (axis[0] * axis[0] + axis[1] * axis[1]);

	/* If the z-axis hasn't changed, nothing changes. */
	if (*theta < DBL_MIN) {
		*theta = *chi = 0;
		return 0;
	}

	/* The y-coordinate of the new z-axis vector is given by
	 * y = sin(theta) sin(chi). Hence, chi is easily computed. Note
	 * that the x-coordiante is x = sin(theta) cos(chi), but this
	 * doesn't work well, since cosine is an even function. */
	*chi = asin (axis[1] / st);

	return 1;
}

/* Build the initial values of the Hvn function for rotation. */
int buildhvn (double theta, double *hvn, int nmax, int mmax) {
	int i, j, dm1, lda, sgn;
	double fact, *buf, arg;

	dm1 = nmax - 1;
	lda = 2 * mmax - 1;

	arg = cos (theta);

	buf = malloc (nmax * sizeof(double));

	/* Build the zero-order array. */
	gsl_sf_legendre_sphPlm_array (dm1, 0, arg, buf);
	fact = sqrt(4 * M_PI);

	for (j = 0; j < nmax; ++j)
		hvn[IDX(j,0,lda)] = fact * buf[j]; 

	for (i = 1; i < mmax; ++i) {
		/* Build the other orders. */
		gsl_sf_legendre_sphPlm_array (dm1, i, arg, buf);
		sgn = 1 - 2 * (i % 2);

		for (j = i; j < nmax; ++j) {
			fact = sqrt (4 * M_PI / (2 * j + 1));
			/* The positive order for this degree. */
			hvn[IDX(j,i,lda)] = sgn * fact * buf[j - i]; 
			/* The negative order for this degree. */
			hvn[IDX(j,-i,lda)] = sgn * fact * buf[j - i]; 
		}
	}

	free (buf);

	return nmax;
}

int nexthvn (double theta, double *hvn, int m, int nmax, int mmax) {
	int lda, i, j;
	double st, omct, opct, bnm;

	lda = 2 * mmax - 1;

	st = sin(theta);
	omct = 1.0 - cos(theta);
	opct = 1.0 + cos(theta);

	/* Perform the recursion. */
	for (i = 1; i < nmax; ++i) {
		bnm = BNM(i,m);

		/* The zero-order coefficient. */
		hvn[IDX(i-1,0,lda)] = (0.5 * BNM(i,-1) * omct * hvn[IDX(i,1,lda)]
				- 0.5 * BNM(i,-1) * opct * hvn[IDX(i,-1,lda)]
				- ANM(i-1,0) * st * hvn[IDX(i,0,lda)]) / bnm;

		/* The non-zero orders. */
		for (j = 1; j < i; ++i) {
			/* Positive j. */
			hvn[IDX(i-1,j,lda)] = (0.5 * BNM(i,-j-1) * omct * hvn[IDX(i,j+1,lda)]
				- 0.5 * BNM(i,j-1) * opct * hvn[IDX(i,j-1,lda)]
				- ANM(i-1,j) * st * hvn[IDX(i,j,lda)]) / bnm;
			/* Negative j. */
			hvn[IDX(i-1,-j,lda)] = (0.5 * BNM(i,j-1) * omct * hvn[IDX(i,-j+1,lda)]
				- 0.5 * BNM(i,-j-1) * opct * hvn[IDX(i,-j-1,lda)]
				- ANM(i-1,-j) * st * hvn[IDX(i,-j,lda)]) / bnm;
		}
	}

	return nmax * mmax;
}

/* Rotate the SH coefficients according to rotation angles throt and chrot. */
int shrotate (complex double *vin, int deg, int lda, trdesc *trans) {
	double *hvn, theta, chi;
	complex double pfz, mfz, *avp, *avm, *buf,
		iscale[4] = { 1, -I, -1, I }, imscale[4] = { 1, I, -1, -I };
	int nmax, i, j, m, idx, nidx;

	nmax = 2 * deg - 1;

	hvn = malloc (nmax * nmax * sizeof(double));
	buf = malloc (deg * nmax * sizeof(complex double));

	/* Find the rotation angles of the new axis. */
	getangles (&theta, &chi, trans->sdir);

	/* Build the initial values of the H translation function. */
	buildhvn (theta, hvn, nmax, deg);

	/* Handle the m = 0 case specifically. */
	for (i = 0; i < deg; ++i) {
		avp = vin + IDX(i,0,lda);
		buf[IDX(i,0,nmax)] = (*avp) * hvn[IDX(i,0,nmax)];
		
		/* Handle the positive and negative orders for these degrees. */
		for (j = 1; j <= i; ++j) {
			idx = IDX(i,j,nmax);
			buf[idx] = (*avp) * hvn[idx];
			idx = IDX(i,-j,nmax);
			buf[idx] = (*avp) * hvn[idx];
		}
	}

	/* Handle the higher-order cases. */
	for (m = 1; m < deg; ++m) {
		/* Calculate the Hvn samples for the next m. */
		nexthvn (theta, hvn, m - 1, nmax, deg);
		
		/* Calculate the phase terms. */
		pfz = cexp (I * m * chi);
		mfz = cexp (-I * m * chi);
		
		for (i = m; i < deg; ++i) {
			avp = vin + IDX(i,m,lda);
			avm = vin + IDX(i,-m,lda);

			buf[IDX(i,0,nmax)] += hvn[IDX(i,0,nmax)] * (pfz * (*avp) + mfz * (*avm));
			
			/* Handle the positive and negative orders for these degrees. */
			for (j = 1; j <= i; ++j) {
				idx = IDX(i,j,nmax);
				nidx = IDX(i,-j,nmax);

				buf[idx] += pfz * (*avp) * hvn[idx]
					+ mfz * (*avm) * hvn[nidx];

				buf[nidx] += pfz * (*avp) * hvn[nidx]
					+ mfz * (*avm) * hvn[idx];
			}
		}
	}

	/* Copy the result back to the input buffer. */
	for (i = 0; i < deg; ++i) {
		vin[IDX(i,0,lda)] = buf[IDX(i,0,nmax)];

		/* Have to scale by (-i)^nu, where nu is the order. */
		for (j = 1; j <= i; ++j) {
			vin[IDX(i,j,lda)] = iscale[j % 4] * buf[IDX(i,j,nmax)];
			vin[IDX(i,-j,lda)] = imscale[j % 4] * buf[IDX(i,-j,nmax)];
		}
	}

	free (hvn);
	free (buf);

	return deg;
}

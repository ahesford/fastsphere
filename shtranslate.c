#include <stdlib.h>

#include <complex.h>
#include <math.h>

#include "util.h"

#define ANM(n,m) sqrt((double)((n)+1+(m))*((n)+1-(m))/(double)((2*(n)+1)*(2*(n)+3)))
#define BNM(n,m) SGN(m)*sqrt((double)((n)-(m)-1)*((n)-(m))/(double)((2*(n)-1)*(2*(n)+1)))

int shtransrfl (complex double *trans, int n, int m) {
	int l, llim, sgnm, sgnl, nmax;

	llim = n - m;
	sgnm = 1 - 2 * (m % 2);

	/* Mirror the fixed-n line to the fixed-l line. */
	for (l = m + 1; l < llim; ++l) {
		sgnl = 1 - 2 * (l % 2);
		trans[ELT(m,l,n)] = sgnl * sgnm * trans[ELT(l,m,n)];
	}

	return llim - m - 1;
}

int shtransfill (complex double *trans, int deg, int m) {
	int nlim, i, l, n, nmax, sgnn, sgnl;
	double anm, alm, alm1, anm1;

	/* The maximum degree we will be computing. */
	nmax = 2 * deg - 1;

	/* Handle the n = m case, with reflection. */
	nlim = nmax - m - 1;
	n = m;
	anm = ANM(n,m);
	/* This is the negative of (-1)^n since (-1)^(n+1) is required. */
	sgnn = 2 * (n % 2) - 1;
	for (l = n + 1; l < nlim; ++l) {
		sgnl = 1 - 2 * (l % 2);
		alm = ANM(l,m);
		alm1 = ANM(l-1,m);

		trans[ELT(l,n+1,nmax)] = (alm1 * trans[ELT(l-1,n,nmax)]
			- alm * trans[ELT(l+1,n,nmax)]) / anm;

		trans[ELT(n+1,l,nmax)] = sgnn * sgnl * trans[ELT(l,n+1,nmax)];
	}

	/* Now handle the subsequent n lines, with reflection. */
	for (n = m + 1; n < deg; ++n) {
		anm = ANM(n,m);
		anm1 = ANM(n-1,m);
		sgnn = 2 * (n % 2) - 1;

		nlim = nmax - n - 1;
		for (l = n + 1; l < nlim; ++l) {
			sgnl = 1 - 2 * (l % 2);
			alm = ANM(l,m);
			alm1 = ANM(l-1,m);

			trans[ELT(l,n+1,nmax)] = (alm1 * trans[ELT(l-1,n,nmax)]
					- alm * trans[ELT(l+1,n,nmax)]
					- anm1 * trans[ELT(l,n-1,nmax)]) / anm;

			trans[ELT(n+1,l,nmax)] = sgnn * sgnl * trans[ELT(l,n+1,nmax)];
		}
	}

	return deg;
}

int shtransinit (complex double *trans, int deg, complex double kr) {
	int nmax, i, sgn;

	nmax = 2 * deg - 1;

	/* Build the initial coefficients for translation for fixed n. */
	spbesh (trans, kr, nmax);

	/* Scale these values, then clone them for fixed l. */
	for (i = 1; i < nmax; ++i) {
		sgn = 1 - 2 * (i % 2);
		/* Scale the fixed-l value. */
		trans[ELT(0,i,nmax)] = sqrt(2.0 * i + 1.0) * trans[ELT(i,0,nmax)];
		/* The fixed-l value also needs the sign factored in. */
		trans[ELT(i,0,nmax)] = sgn * trans[ELT(0,i,nmax)];
	}

	shtransfill (trans, deg, 0);

	return deg;
}

int shtransnext (complex double *trans, int deg, int m) {
	int nmax, llim, l, sgnn, sgnl, mp1;
	double bnm, blm, blm1;

	nmax = 2 * deg - 1;
	llim = nmax - 1 - m;

	/* Build the n = m + 1 case for all l, and reflect. */
	mp1 = m + 1;
	bnm = BNM(mp1, -mp1);
	sgnn = 1 - 2 * (mp1 % 2);
	for (l = mp1; l < llim; ++l) {
		sgnl = 1 - 2 * (l % 2);
		blm1 = BNM(l+1,m);
		blm = BNM(l,-mp1);

		trans[ELT(l,mp1,nmax)] = (blm1 * trans[ELT(l+1,m,nmax)]
				- blm * trans[ELT(l-1,m,nmax)]) / bnm;
		trans[ELT(mp1,l,nmax)] = sgnn * sgnl * trans[ELT(l,mp1,nmax)];
	}

	/* Fill out the n > m + 1 cases. */
	shtransfill (trans, deg, mp1);

	return mp1;
}

int shtranslate (complex double *coeff, int deg, int lda, complex double kr) {
	int l, nmax, m, n;
	complex double *trans, *buf;

	nmax = 2 * deg - 1;

	trans = malloc (nmax * nmax * sizeof(complex double));
	buf = calloc (deg * nmax, sizeof(complex double));

	/* Build the translator for m = 0. */
	shtransinit (trans, deg, kr);

	/* Find the coefficients for m = 0. */
	for (l = 0; l < deg; ++l)
		for (n = 0; n < deg; ++n) 
			buf[ELT(0,l,nmax)] += trans[ELT(l,n,nmax)] * coeff[ELT(0,n,lda)];

	shtransnext (trans, deg, 0);

	/* Find the coefficients for values abs(m) > 0. */
	for (m = 1; m < deg; ++m) {
		for (l = m; l < deg; ++l) {
			for (n = m; n < deg; ++n) {
				/* The positive orders. */
				buf[ELT(m,l,nmax)] += trans[ELT(l,n,nmax)] * coeff[ELT(m,n,lda)];
				/* The negative orders. */
				buf[ELT(-m,l+1,nmax)] += trans[ELT(l,n,nmax)] * coeff[ELT(-m,n+1,lda)];
			}
		}

		/* Build the translator for the next order. */
		shtransnext (trans, deg, m);
	}

	/* Copy the new coefficients over the old values. */
	for (l = 0; l < deg; ++l) {
		/* The zero-order. */
		coeff[ELT(0,l,lda)] = buf[ELT(0,l,nmax)];
		for (m = 1; m <= l; ++m) {
			/* The positive orders. */
			coeff[ELT(m,l,lda)] = buf[ELT(m,l,nmax)];
			/* The negative orders. */
			coeff[ELT(-m,l+1,lda)] = buf[ELT(-m,l+1,nmax)];
		}
	}

	free (trans);
	free (buf);

	return deg;
}

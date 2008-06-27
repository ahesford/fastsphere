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
	int nmax, ilim, i, sgnn, sgnl;
	double bnm, blm, blm1;
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <math.h>
#include <float.h>

#include "shrotate.h"
#include "util.h"

/* Find the polar and rotation angles of the new z-axis for translation. */
int getangles (double *theta, double *chi, double *phi, double axis[3]) {
	double st;

	*theta = acos (axis[2]); /* The polar angle of the new z-axis. */

	st = sqrt (axis[0] * axis[0] + axis[1] * axis[1]);

	/* The default rotation angle of the old z-axis in the new
	 * coordinate system is pi / 2. */
	*phi = M_PI_2;

	if (st < DBL_EPSILON) {
		*chi = 0;
		return 0;
	}

	/* The rotation angle of the new z-axis. */
	*chi = acos (axis[0] / st);

	/* Rewrap the angle if it's in the lower half plane. */
	if (axis[1] < 0) *chi = 2 * M_PI - *chi;

	return 1;
}

void coordmat (complex double *dmat, double theta, double chi, double phi) {
	double st, ct, sc, cc, sp, cp;
	double r[9];

	ct = cos(theta);
	st = sin(theta);

	cc = cos(chi);
	sc = sin(chi);

	cp = cos(phi);
	sp = sin(phi);

	/* The angle phi is always pi / 2, so is ignored. */
	r[0] = -(sp * sc + ct * cp * cc);
	r[1] = cp * sc - ct * sp * cc;
	r[2] = st * cc;
	r[3] = sp * cc - ct * cp * sc;
	r[4] = -(cp * cc + ct * sp * sc);
	r[5] = st * sc;
	r[6] = st * cp;
	r[7] = st * sp;
	r[8] = ct;

	/* First column. */
	dmat[0] = 0.5 * (r[4] + r[0] + I * (r[1] - r[3]));
	dmat[1] = M_SQRT1_2 * (r[2] - I * r[5]);
	dmat[2] = 0.5 * (r[4] - r[0] + I * (r[1] + r[3]));

	/* Second column. */
	dmat[3] = M_SQRT1_2 * (r[6] + I * r[7]);
	dmat[4] = r[8];
	dmat[5] = M_SQRT1_2 * (-r[6] + I * r[7]);

	/* Third column. */
	dmat[6] = 0.5 * (r[4] - r[0] - I * (r[1] + r[3]));
	dmat[7] = -M_SQRT1_2 * (r[2] + I * r[5]);
	dmat[8] = 0.5 * (r[4] + r[0] + I * (r[3] - r[1]));
}

int nextrot (complex double *rot, int l, int lda, complex double *dmat) {
	int mp, m, cnum, dnum, cden, ldb;
	double cpm, cmm, dmm;
	complex double *buf;

	ldb = 2 * l + 1;
	buf = malloc ((l + 1) * ldb * sizeof(complex double));

	for (mp = 0; mp < l; ++mp) {
		/* First handle the m = 0 case. */
		cden = (l + mp) * (l - mp);
		cnum = l * (l - 1);
		cpm = sqrt((double)cnum / (2.0 * (double)cden));
		dmm = (double)l / sqrt((double)cden);

		buf[ELT(0,mp,ldb)] = dmm * dmat[4] * rot[ELT(0,mp,lda)]
			- dmat[5] * cpm * rot[ELT(-1,mp+1,lda)] 
			- dmat[3] * cpm * rot[ELT(1,mp,lda)];

		/* Now handle the nonzero m. */
		for (m = 1; m <= l; ++m) {
			cnum = (l + m) * (l + m - 1);
			cpm = sqrt((double)cnum / (2.0 * (double)cden));
			cnum = (l - m) * (l - m - 1);
			cmm = sqrt((double)cnum / (2.0 * (double)cden));
			dnum = (l + m) * (l - m);
			dmm = sqrt((double)dnum / (double)cden);
			
			buf[ELT(m,mp,ldb)] = dmm * dmat[4] * rot[ELT(m,mp,lda)]
				- dmat[5] * cpm * rot[ELT(m-1,mp,lda)]
				- dmat[3] * cmm * rot[ELT(m+1,mp,lda)];
			/* For the negative, -m+1 can be nonnegative, so
			 * use the IDX macro to avoid errors. */
			buf[ELT(-m,mp+1,ldb)] = dmm * dmat[4] * rot[ELT(-m,mp+1,lda)]
				- dmat[5] * cmm * rot[ELT(-m-1,mp+1,lda)]
				- dmat[3] * cpm * rot[IDX(-m+1,mp,lda)];
		}
	}

	/* Handle the final case. */
	mp = l;
	cden = (l + mp) * (l + mp - 1);
	cnum = l * (l - 1);
	cpm = sqrt((double)cnum / (double)cden);
	dmm = (double)l * sqrt(2.0 / (double)cden);

	buf[ELT(0,mp,ldb)] = dmat[8] * cpm * rot[ELT(-1,mp,lda)]
		+ dmat[6] * cpm * rot[ELT(1,mp-1,lda)]
		- dmat[7] * dmm * rot[ELT(0,mp-1,lda)];

	for (m = 1; m <= l; ++m) {
		cnum = (l + m) * (l + m - 1);
		cpm = sqrt((double)cnum / (double)cden);
		cnum = (l - m) * (l - m - 1);
		cmm = sqrt((double)cnum / (double)cden);
		dnum = 2 * (l + m) * (l - m);
		dmm = sqrt((double)dnum / (double)cden);

		buf[ELT(m,mp,ldb)] = dmat[8] * cpm * rot[ELT(m-1,mp-1,lda)]
			+ dmat[6] * cmm * rot[ELT(m+1,mp-1,lda)]
			- dmat[7] * dmm * rot[ELT(m,mp-1,lda)];

		/* Use the IDX macro to avoid issues when -m+1 is
		 * not negative, just as above. */
		buf[ELT(-m,mp+1,ldb)] = dmat[8] * cmm * rot[ELT(-m-1,mp,lda)]
			+ dmat[6] * cpm * rot[IDX(-m+1,mp-1,lda)]
			- dmat[7] * dmm * rot[ELT(-m,mp,lda)];
	}

	/* Copy the buffer into the original location. */
	for (mp = 0; mp <= l; ++mp) {
		rot[ELT(0,mp,lda)] = buf[ELT(0,mp,ldb)];

		for (m = 1; m <= l; ++m) {
			rot[ELT(m,mp,lda)] = buf[ELT(m,mp,ldb)];
			rot[ELT(-m,mp+1,lda)] = buf[ELT(-m,mp+1,ldb)];
		}
	}

	free (buf);
	return l * lda;
}

/* Rotate the SH coefficients into the new axis. */
int shrotate (complex double *vin, int deg, int lda,
		double theta, double chi, double phi) {
	complex double *rot, *buf, *vptr, dmat[9];
	int l, m, mp, ldb, sm;

	ldb = 2 * deg - 1;

	rot = calloc (deg * ldb, sizeof(complex double));
	buf = malloc (ldb * sizeof(complex double));

	/* The initial rotator doesn't rotate anything. */
	rot[0] = 1.0;
	coordmat (dmat, theta, chi, phi);

	for (l = 1; l < deg; ++l) {
		/* Compute the next rotation from the previous one. */
		nextrot (rot, l, ldb, dmat);

		/* Find the column offset. */
		vptr = vin + ELT(0,l,lda);

		/* Handle the zero-order coefficient. */
		*buf = rot[ELT(0,0,ldb)] * vptr[0];

		for (m = 1; m <= l; ++m) {
			sm = 1 - 2 * (m % 2);
			*buf += rot[ELT(m,0,ldb)] * vptr[m];
			*buf += sm * rot[ELT(-m,1,ldb)] * vptr[lda-m];
		}

		/* Handle the nonzero orders. */
		for (mp = 1; mp <= l; ++mp) {
			buf[mp] = rot[ELT(0,mp,ldb)] * vptr[0];
			buf[ldb-mp] = conj(rot[ELT(0,mp,ldb)]) * vptr[0];
			
			for (m = 1; m <= l; ++m) {
				sm = 1 - 2 * (m % 2);
				buf[mp] += rot[ELT(m,mp,ldb)] * vptr[m];
				buf[mp] += sm * rot[ELT(-m,mp+1,ldb)] * vptr[lda-m];
				buf[ldb-mp] += sm * conj(rot[ELT(-m,mp+1,ldb)]) * vptr[m];
				buf[ldb-mp] += conj(rot[ELT(m,mp,ldb)]) * vptr[lda-m];
			}
		}

		/* Copy the finished column into the input. */
		*vptr = *buf;
		for (mp = 1; mp <= l; ++mp) {
			vptr[mp] = buf[mp];
			vptr[lda-mp] = buf[ldb-mp];
		}
	}

	free (rot);
	free (buf);

	return deg;
}

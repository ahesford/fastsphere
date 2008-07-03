#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <math.h>
#include <float.h>

#include "shrotate.h"
#include "util.h"

#define SC (-I * M_SQRT1_2)

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

int nextrot (complex double *rot, int l, int lda) {
	int mp, m, cnum, dnum, cden, ldb;
	double cpm, cmm, dmm;
	complex double *buf;

	ldb = 2 * l + 1;
	buf = malloc ((l + 1) * ldb * sizeof(complex double));

	for (mp = 0; mp < l; ++mp) {
		/* First handle the m = 0 case. */
		cden = 2 * (l + mp) * (l - mp);
		cnum = l * (l - 1);
		cpm = sqrt((double)cnum / (double)cden);

		buf[ELT(0,mp,ldb)] =
			SC * cpm * (rot[ELT(-1,mp+1,lda)] + rot[ELT(1,mp,lda)]);

		/* Now handle the nonzero m. */
		for (m = 1; m <= l; ++m) {
			cnum = (l + m) * (l + m - 1);
			cpm = sqrt((double)cnum / (double)cden);
			cnum = (l - m) * (l - m - 1);
			cmm = sqrt((double)cnum / (double)cden);
			
			buf[ELT(m,mp,ldb)] = SC * (cpm * rot[ELT(m-1,mp,lda)]
					+ cmm * rot[ELT(m+1,mp,lda)]);
			/* For the negative, -m+1 can be nonnegative, so
			 * use the IDX macro to avoid errors. */
			buf[ELT(-m,mp+1,ldb)] = SC * (cmm * rot[ELT(-m-1,mp+1,lda)]
					+ cpm * rot[IDX(-m+1,mp,lda)]);
		}
	}

	/* Handle the final case. */
	mp = l;
	cden = (l + mp) * (l + mp - 1);
	cnum = l * (l - 1);
	cpm = sqrt((double)cnum / (double)cden);
	dmm = (double)l * sqrt(2.0 / (double)cden);

	buf[ELT(0,mp,ldb)] = cpm * (rot[ELT(1,mp-1,lda)] - rot[ELT(-1,mp,lda)]) / 2.0
		- SC * dmm * rot[ELT(0,mp-1,lda)];

	for (m = 1; m <= l; ++m) {
		cnum = (l + m) * (l + m - 1);
		cpm = sqrt((double)cnum / (double)cden);
		cnum = (l - m) * (l - m - 1);
		cmm = sqrt((double)cnum / (double)cden);
		dnum = 2 * (l + m) * (l - m);
		dmm = sqrt((double)dnum / (double)cden);

		buf[ELT(m,mp,ldb)] = (cmm * rot[ELT(m+1,mp-1,lda)]
				- cpm * rot[ELT(m-1,mp-1,lda)]) / 2.0
			- SC * dmm * rot[ELT(m,mp-1,lda)];

		/* Use the IDX macro to avoid issues when -m+1 is
		 * not negative, just as above. */
		buf[ELT(-m,mp+1,ldb)] = (cpm * rot[IDX(-m+1,mp-1,lda)]
				- cmm * rot[ELT(-m-1,mp,lda)]) / 2.0
			- SC * dmm * rot[ELT(-m,mp,lda)];
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

/* Rotate about the z-axis by an angle theta. */
int zrotate (complex double *vin, int deg, int lda, double theta) {
	int l, m;
	complex double fact;

	for (m = 1; m < deg; ++m) {
		fact = cexp(I * m * theta);
		for (l = m; l < deg; ++l) {
			vin[ELT(m,l,lda)] *= fact;
			vin[ELT(-m,l+1,lda)] *= conj(fact);
		}
	}

	return lda;
}

/* Rotate the SH coefficients so the old y-axis coincides with the
 * new z-axis, the old z-axis coincides with the new y-axis and the
 * old x-axis is the negative of the new x-axis. */
int xrotate (complex double *vin, int deg, int lda) {
	complex double *rot, *buf, *vptr;
	int l, m, mp, ldb, sm;

	ldb = 2 * deg - 1;

	rot = calloc (deg * ldb, sizeof(complex double));
	buf = malloc (ldb * sizeof(complex double));

	/* The initial rotator doesn't rotate anything. */
	rot[0] = 1.0;

	for (l = 1; l < deg; ++l) {
		/* Compute the next rotation from the previous one. */
		nextrot (rot, l, ldb);

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

int shrotate (complex double *vin, int deg, int lda,
		double theta, double chi, double phi) {
	int nrot = 0;

	if (fabs(chi) > DBL_EPSILON) {
		zrotate (vin, deg, lda, chi);
		++nrot;
	}

	if (fabs(theta) > DBL_EPSILON) {
		xrotate (vin, deg, lda);
		zrotate (vin, deg, lda, theta);
		xrotate (vin, deg, lda);
		++nrot;
	}

	if (fabs(phi) > DBL_EPSILON) {
		zrotate (vin, deg, lda, phi);
		++nrot;
	}

	return nrot;
}

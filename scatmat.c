#include <stdlib.h>
#include <complex.h>
#include <string.h>

#include "fastsphere.h"
#include "spreflect.h"
#include "scatmat.h"
#include "fsht.h"

void drivezgmres_ (int *, int *, int *, int *, complex double *,
		int *, int *, double *, int *, double *);
void initzgmres_ (int *, double *);
void zgemv_ (char *, int *, int *, complex double *, complex double *, int *,
		complex double *, int *, complex double *, complex double *, int *);

int buildrhs (complex double *rhs, spscat *spl, int nsph, shdata *shtr) {
	int i, nterm;
	complex double *vptr;

	nterm = shtr->ntheta * shtr->nphi;

#pragma omp parallel for private(i,vptr) default(shared)
	for (i = 0; i < nsph; ++i) {
		vptr = rhs + i * nterm;

		/* Multiply by the reflection coefficient in SH space. */
		ffsht (vptr, shtr);
		spreflect (vptr, vptr, (spl + i)->spdesc->reflect, shtr->deg, shtr->nphi);
		ifsht (vptr, shtr);
	}

	return nsph;
}

int scatmat (complex double *vout, complex double *vin, spscat *spl,
		int nsph, complex **trans, shdata *shtr) {
	int i, j, k, off, nterm;
	complex double *buf, *voptr, *viptr;

	nterm = shtr->ntheta * shtr->nphi;

	buf = malloc (nterm * sizeof(complex double));

	for (i = 0; i < nsph; ++i) {
		voptr = vout + i * nterm;

		/* Copy the outbound field at the sphere. */
		memcpy (voptr, vin + i * nterm, nterm * sizeof(complex double));
		memset (buf, 0, nterm * sizeof(complex double));

		/* Add up the translated fields. */
		for (j = 0; j < nsph; ++j) {
			if (j == i) continue;
			off = j * nsph + i;
			viptr = vin + j * nterm;

			for (k = 0; k < nterm; ++k) 
				buf[k] += trans[off][k] * viptr[k];
		}

		/* Apply the reflection coefficient in SH space. */
		ffsht (buf, shtr);
		spreflect (buf, buf, (spl + i)->spdesc->reflect, shtr->deg, shtr->nphi);
		ifsht (buf, shtr);

		/* Subtract this from the outbound field. */
		for (k = 0; k < nterm; ++k)
			voptr[k] -= buf[k];
	}

	free (buf);

	return nsph;
}

int itsolve (complex double *sol, complex double *rhs, spscat *spl, int nsph,
		complex **trans, shdata *shtr, itconf *itc) {
	int icntl[7], irc[5], lwork, info[3], i, n, nterm, one = 1;
	double rinfo[2], cntl[5];
	complex double *zwork, *tx, *ty, *tz, zone = 1.0, zzero = 0.0;

	nterm = shtr->ntheta * shtr->nphi;
	n = nterm * nsph;

	lwork = itc->restart * itc->restart + itc->restart * (n + 5) + 5 * n + 1;
	zwork = calloc (lwork, sizeof(complex double));

	initzgmres_ (icntl, cntl);

	icntl[2] = 6; /* Print information to stdout. */
	icntl[3] = 0; /* No preconditioner. */
	icntl[4] = 0; /* Use MGS for orthogonalization. */
	icntl[5] = 1; /* The incident field is the initial guess. */
	icntl[6] = itc->iter; /* Maximum iteration count. */

	cntl[0] = itc->eps; /* Convergence tolerance. */

	/* The initial guess: the RHS. */
	memcpy (zwork, rhs, n * sizeof(complex double));
	/* The unpreconditioned RHS. */
	memcpy (zwork + n, rhs, n * sizeof(complex double));

	do {
		/* The GMRES driver for double complex values. */
		drivezgmres_ (&n, &n, &(itc->restart), &lwork, zwork,
				irc, icntl, cntl, info, rinfo);
		if (!(info[0]) && !(irc[0])) break;

		switch (irc[0]) {
		case GMV:
			tx = zwork + irc[1] - 1;
			ty = zwork + irc[3] - 1;
			scatmat (ty, tx, spl, nsph, trans, shtr);
			break;
		case GDP:
			tx = zwork + irc[1] - 1;
			ty = zwork + irc[2] - 1;
			tz = zwork + irc[3] - 1;

			/* Compute the scalar products in one pass, using
			 * LAPACK for a matrix-vector product. */
			zgemv_ ("C", &n, irc + 4, &zone, tx, &n,
					ty, &one, &zzero, tz, &one);
			break;
		default: break;
		}
	} while (irc[0]);

	if (info[0]) printf ("ZGMRES: return value: %d\n", info[0]);

	memcpy (sol, zwork, n * sizeof(complex double));

	printf ("ZGMRES: %d iterations, %0.6g PBE, %0.6g BE.\n", info[1], rinfo[0], rinfo[1]);

	free (zwork);
}

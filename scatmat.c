#include <stdlib.h>
#include <complex.h>
#include <string.h>

#include "fastsphere.h"
#include "spreflect.h"
#include "translator.h"
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
		int nsph, trdesc *trans, shdata *shtr) {
	int off, k, nterm, n, nsq;
	complex double *buf, *oshc, *voptr, *viptr;

	nterm = shtr->ntheta * shtr->nphi;
	n = nterm * nsph;
	nsq = nsph * nsph;

	buf = calloc (2 * n, sizeof(complex double));
	oshc = buf + n;

	/* Copy the outgoing plane wave coefficients into a temporary vector. */
	memcpy (oshc, vin, n * sizeof(complex double));

	/* Zero out the output buffer. */
	memset (vout, 0, n * sizeof(complex double));

	/* Transform all outgoing plane waves into spherical harmonics. */
#pragma omp parallel for private(off) default(shared)
	for (off = 0; off < nsph; ++off)
		ffsht (oshc + off * nterm, shtr);

	/* Perform the translations. */
#pragma omp parallel private(off,k,voptr,viptr) default(shared)
{
	complex double *dtr;
	int i, j;

	/* Buffer for doing in-place spherical harmonic translations. */
	dtr = malloc (nterm * sizeof(complex double));

#pragma omp for
	for (off = 0; off < nsq; ++off) {
		j = off / nsph;	/* Source sphere. */
		i = off % nsph;	/* Destination sphere. */

		/* Don't bother with self-translations. */
		if (i == j) continue;

		if (trans[off].type == TRPLANE) {
			/* Do the diagonal, plane-wave translation. */
			voptr = vout + i * nterm;
			viptr = vin + j * nterm;
			/* Copy to output, but only one thread at a time. */
#pragma omp critical(outplane)
			for (k = 0; k < nterm; ++k) 
				voptr[k] += trans[off].trdata[k] * viptr[k];
		} else if (trans[off].type == TRDENSE) {
			/* Do the dense, spherical-harmonic translation. */
			memcpy (dtr, oshc + j * nterm, nterm * sizeof(complex double));

			shrotate (dtr, shtr->deg, shtr->nphi, trans[off].theta,
					trans[off].chi, trans[off].phi);
			shtranslate (dtr, shtr->deg, shtr->nphi, trans[off].kr);
			shrotate (dtr, shtr->deg, shtr->nphi, trans[off].theta,
					trans[off].phi, trans[off].chi);
			voptr = buf + i * nterm;
			/* Copy to output, but only one thread at a time. */
#pragma omp critical(outsphere)
			for (k = 0; k < nterm; ++k)
				voptr[k] += dtr[k];
		}
	}

	free (dtr);
}

#pragma omp parallel for private(off, k,voptr,viptr) default(shared)
	for (off = 0; off < nsph; ++off) {
		voptr = vout + off * nterm;
		viptr = buf + off * nterm;

		/* Convert plane-wave translations into harmonics. */
		ffsht (voptr, shtr);
		/* Add in the direct harmonic translations. */
		for (k = 0; k < nterm; ++k) voptr[k] += viptr[k]; 
		
		/* Apply the reflection coefficient in SH space. */
		spreflect (voptr, voptr, (spl + off)->spdesc->reflect, shtr->deg, shtr->nphi);

		/* Take the reflections back to plane-waves. */
		ifsht (voptr, shtr);
	}
	
	/* Subtract the incoming field from the outgoing field. */
#pragma omp parallel for private(off) default(shared)
	for (off = 0; off < n; ++off) vout[off] = vin[off] - vout[off];

	free (buf);

	return nsph;
}

int itsolve (complex double *sol, complex double *rhs, spscat *spl, int nsph,
		trdesc *trans, shdata *shtr, itconf *itc) {
	int icntl[7], irc[5], lwork, info[3], n, nterm, one = 1;
	double rinfo[2], cntl[5];
	complex double *zwork, *tx, *ty, *tz, zone = 1.0, zzero = 0.0;

	nterm = shtr->ntheta * shtr->nphi;
	n = nterm * nsph;

	lwork = itc->restart * itc->restart + itc->restart * (n + 5) + 5 * n + 2;
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

	return info[1];
}

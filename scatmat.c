#include <stdlib.h>
#include <complex.h>
#include <string.h>

#include "fastsphere.h"
#include "spreflect.h"
#include "translator.h"
#include "scatmat.h"
#include "fsht.h"
#include "farfield.h"
#include "sprsinterp.h"

void drivezgmres_ (int *, int *, int *, int *, complex double *,
		int *, int *, double *, int *, double *);
void initzgmres_ (int *, double *);
void zgemv_ (char *, int *, int *, complex double *, complex double *, int *,
		complex double *, int *, complex double *, complex double *, int *);

/* Reflect incoming plane waves from the surfaces of all spheres. */
int sprflpw (complex double *rhs, spscat *spl, int nsph, shdata *shtr) {
	int i, nterm;
	complex double *vptr;

	nterm = shtr->ntheta * shtr->nphi + 2;

#pragma omp parallel for private(i,vptr) default(shared)
	for (i = 0; i < nsph; ++i) {
		vptr = rhs + i * nterm;

		/* Multiply by the reflection coefficient in SH space. */
		ffsht (vptr, shtr);
		spreflect (vptr, vptr, (spl + i)->spdesc->reflect, shtr->deg, shtr->nphi, 0, 1);
		ifsht (vptr, shtr);
	}

	return nsph;
}

/* Compute translations between all spheres. Augments the output vector, does
 * not overwrite it. */
int sptrans (complex double *vout, complex double *vin,
		int nsph, trdesc *trans, shdata *shtr) {
	int nterm, nsq;

	nterm = shtr->ntheta * shtr->nphi + 2;
	nsq = nsph * nsph;

	/* Perform the translations. */
#pragma omp parallel default(shared)
{
	complex double *voptr, *viptr;
	int i, j, off, k;

#pragma omp for
	for (off = 0; off < nsq; ++off) {
		j = off / nsph;	/* Source sphere. */
		i = off % nsph;	/* Destination sphere. */

		/* Don't bother with self-translations. Also ignore dense
		 * translations for the moment. */
		if (i == j || trans[off].type != TRPLANE) continue;

		/* Do the diagonal, plane-wave translation. */
		voptr = vout + i * nterm;
		viptr = vin + j * nterm;
		/* Copy to output, but only one thread at a time. */
#pragma omp critical(outplane)
		for (k = 0; k < nterm; ++k) 
			voptr[k] += trans[off].trdata[k] * viptr[k];
	}
}

	return nsph;
}

/* Compute the MVP between the scattering matrix and a specified vector. */
int scatmat (complex double *vout, complex double *vin, spscat *spl, int nsph,
		sptype *bgspt, trdesc *trans, shdata *shtr, shdata *bgtr, sprow *imat) {
	int nterm, n, i, ntbg;
	complex double *voptr, *viptr;

	nterm = shtr->ntheta * shtr->nphi + 2;
	ntbg = bgtr->ntheta * bgtr->nphi + 2;
	n = nterm * nsph;

	voptr = vout + n;
	viptr = vin + n;

	/* Initialize the output bufer. */
	fartonear (vout, viptr, spl, nsph, bgspt->k, shtr, bgtr, imat);

	/* Compute the spherical translations. */
	sptrans (vout, vin, nsph, trans, shtr);

	/* Compute the reflections of plane waves at sphere surfaces. */
	sprflpw (vout, spl, nsph, shtr);

	/* Subtract the incoming field from the outgoing field. */
#pragma omp parallel for private(i) default(shared)
	for (i = 0; i < n; ++i) vout[i] = vin[i] - vout[i];

	/* Compute the far-field pattern of the internal spheres. */
	neartofar (voptr, vin, spl, nsph, bgspt->k, bgtr, shtr, imat);

	/* Compute the reflection of the far-field pattern. */
	ffsht (voptr, bgtr);
	spreflect (voptr, voptr, bgspt->reflect, bgtr->deg, bgtr->nphi, 0, 1);
	ifsht (voptr, bgtr);

	/* Add the contribution from the standing-wave pattern. */
#pragma omp parallel for private(i) default(shared)
	for (i = 0; i < ntbg; ++i) voptr[i] += viptr[i];

	return nsph;
}

int itsolve (complex double *sol, complex double *rhs, spscat *spl,
		int nsph, sptype *bgspt, trdesc *trans, shdata *shtr,
		shdata *bgtr, itconf *itc, sprow *imat) {
	int icntl[7], irc[5], lwork, info[3], n, nterm, one = 1;
	double rinfo[2], cntl[5];
	complex double *zwork, *tx, *ty, *tz, zone = 1.0, zzero = 0.0;

	nterm = shtr->ntheta * shtr->nphi + 2;
	n = nterm * nsph + bgtr->ntheta * bgtr->nphi + 2;

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
			scatmat (ty, tx, spl, nsph, bgspt, trans, shtr, bgtr, imat);
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

#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <string.h>

#include "fastsphere.h"
#include "spreflect.h"
#include "translator.h"
#include "scatmat.h"
#include "fsht.h"
#include "farfield.h"

void drivezgmres_ (int *, int *, int *, int *, complex double *,
		int *, int *, double *, int *, double *);
void initzgmres_ (int *, double *);
void zgemv_ (char *, int *, int *, complex double *, complex double *, int *,
		complex double *, int *, complex double *, complex double *, int *);

static complex double pardot (complex double *x, complex double *y, int n) {
	complex double dp;

	/* Compute the local portion. */
	cblas_zdotc_sub (n, x, 1, y, 1, &dp);

	return dp;
}

/* Reflect incoming plane waves from the surfaces of all spheres. */
int sprflpw (complex double *rhs, spscat *spl, int nsph, shdata *shtr) {
	int i, nterm;
	complex double *vptr;
	spscat *sp;

	nterm = shtr->ntheta * shtr->nphi;

#pragma omp parallel for private(i,vptr,sp) default(shared)
	for (i = 0; i < nsph; ++i) {
		sp = spl + i;
		vptr = rhs + i * nterm;

		/* Multiply by the reflection coefficient in SH space. */
		ffsht (vptr, shtr, sp->spdesc->deg);
		spreflect (vptr, vptr, (spl + i)->spdesc->reflect,
				sp->spdesc->deg, shtr->nphi, 0, 1);
		ifsht (vptr, shtr, sp->spdesc->deg);
	}

	return nsph;
}

/* Compute translations between all spheres. Augments the output vector, does
 * not overwrite it. */
int sptrans (complex double *vout, complex double *vin,
		int nsph, trdesc *trans, shdata *shtr) {
	int nterm, nsq;

	nterm = shtr->ntheta * shtr->nphi;
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
int scatmat (complex double *vout, complex double *vin, spscat *spl,
		int nsph, trdesc *trans, shdata *shtr) {
	int nterm, n, i;

	nterm = shtr->ntheta * shtr->nphi;
	n = nterm * nsph;

	/* Initialize the output bufer. */
	memset (vout, 0, n * sizeof(complex double));

	/* Compute the spherical translations. */
	sptrans (vout, vin, nsph, trans, shtr);

	/* Compute the reflections of plane waves at sphere surfaces. */
	sprflpw (vout, spl, nsph, shtr);

	/* Subtract the incoming field from the outgoing field. */
#pragma omp parallel for private(i) default(shared)
	for (i = 0; i < n; ++i) vout[i] = vin[i] - vout[i];

	return nsph;
}

int bicgstab (complex double *sol, complex double *rhs, int guess, spscat *spl,
		int nsph, trdesc *trans, shdata *shtr, itconf *itc) {
	int i, j, n, nterm;
	complex double *r, *rhat, *v, *p, *t;
	complex double rho, alpha, omega, beta;
	float err, rhn;

	nterm = shtr->ntheta * shtr->nphi;
	n = nterm * nsph;

	rho = alpha = omega = 1.;

	/* Allocate and zero the work arrays. */
	r = calloc (5 * n, sizeof(complex double));
	rhat = r + n;
	v = rhat + n;
	p = v + n;
	t = p + n;

	/* Compute the norm of the right-hand side for residual scaling. */
	rhn = sqrt(creal(pardot (rhs, rhs, n)));

	/* Compute the inital matrix-vector product for the input guess. */
	if (guess) scatmat (r, sol, spl, nsph, trans, shtr);

	/* Subtract from the RHS to form the residual. */
#pragma omp parallel for default(shared) private(j)
	for (j = 0; j < n; ++j) r[j] = rhs[j] - r[j];

	if (!guess) memset (sol, 0, n * sizeof(complex double));
		
	/* Copy the initial residual as the test vector. */
	memcpy (rhat, r, n * sizeof(complex double));

	/* Find the norm of the initial residual. */
	err = sqrt(creal(pardot (r, r, n))) / rhn;
	printf ("True residual: %g\n", err);

	/* Run iterations until convergence or the maximum is reached. */
	for (i = 0; i < itc->iter && err > itc->eps; ++i) {
		/* Pre-compute portion of beta from previous iteration. */
		beta = alpha / (rho * omega);
		/* Compute rho for this iteration. */
		rho = pardot (rhat, r, n);
		/* Include the missing factor in beta. */
		beta *= rho;

		/* Update the search vector. */
#pragma omp parallel for default(shared) private(j)
		for (j = 0; j < n; ++j)
			p[j] = r[j] + beta * (p[j] - omega * v[j]);

		/* Compute the first search step, v = A * p. */
		scatmat (v, p, spl, nsph, trans, shtr);

		/* Compute the next alpha. */
		alpha = rho / pardot (rhat, v, n);

#pragma omp parallel for default(shared) private(j)
		for (j = 0; j < n; ++j) {
			/* Update the solution vector. */
			sol[j] += alpha * p[j];
			/* Update the residual vector. */
			r[j] -= alpha * v[j];
		}

		/* Compute the scaled residual norm and stop if convergence
		 * has been achieved. */
		err = sqrt(creal(pardot (r, r, n))) / rhn;
		printf ("BiCG-STAB(%0.1f): %g\n", 0.5 + i, err);
		if (err < itc->eps) break;

		/* Compute the next search step, t = A * r. */
		scatmat (t, r, spl, nsph, trans, shtr);

		/* Compute the update direction. */
		omega = pardot (t, r, n) / pardot (t, t, n);

		/* Update both the residual and the solution guess. */
#pragma omp parallel for default(shared) private(j)
		for (j = 0; j < n; ++j) {
			/* Update the solution vector. */
			sol[j] += omega * r[j];
			/* Update the residual vector. */
			r[j] -= omega * t[j];
		}
	
		/* Compute the scaled residual norm. */
		err = sqrt(creal(pardot (r, r, n))) / rhn;
		printf ("BiCG-STAB(%d): %g\n", i + 1, err);
	}

	free (r);
	return i;
}

int itsolve (complex double *sol, complex double *rhs, spscat *spl, int nsph,
		trdesc *trans, shdata *shtr, itconf *itc) {
	int icntl[8], irc[5], lwork, info[3], n, nterm, one = 1;
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
	icntl[5] = 1; /* Specify an initial guess. */
	icntl[6] = itc->iter; /* Maximum iteration count. */

	cntl[0] = itc->eps; /* Convergence tolerance. */

	/* The initial guess: the previous solution. */
	memcpy (zwork, sol, n * sizeof(complex double));
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

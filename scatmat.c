#include <stdlib.h>
#include <complex.h>
#include <float.h>
#include <math.h>
#include <string.h>

#ifdef _MACOSX
#include <Accelerate/Accelerate.h>
#else
#include <cblas.h>
#endif

#include "fastsphere.h"
#include "spreflect.h"
#include "translator.h"
#include "scatmat.h"
#include "fsht.h"
#include "farfield.h"
#include "util.h"

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
	double err, rhn;

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
	rhn = cblas_dznrm2 (n, rhs, 1);

	/* Compute the inital matrix-vector product for the input guess. */
	if (guess) scatmat (r, sol, spl, nsph, trans, shtr);

	/* Subtract from the RHS to form the residual. */
#pragma omp parallel for default(shared) private(j)
	for (j = 0; j < n; ++j) r[j] = rhs[j] - r[j];

	if (!guess) memset (sol, 0, n * sizeof(complex double));
		
	/* Copy the initial residual as the test vector. */
	memcpy (rhat, r, n * sizeof(complex double));

	/* Find the norm of the initial residual. */
	err = cblas_dznrm2(n, r, 1) / rhn;
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
		err = cblas_dznrm2 (n, r, 1) / rhn;
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
		err = cblas_dznrm2 (n, r, 1) / rhn;
		printf ("BiCG-STAB(%d): %g\n", i + 1, err);
	}

	free (r);
	return i;
}

int gmres (complex double *sol, complex double *rhs, int guess, spscat *spl,
		int nsph, trdesc *trans, shdata *shtr, itconf *itc) {
	int nterm = shtr->ntheta * shtr->nphi, n = nterm * nsph;
	long lwork;
	int i, j, one = 1, mit = itc->iter;
	complex double *h, *v, *beta, *y;
	complex double *vp, *hp, *s, cr, cone = 1.;
	double rhn, err, *c;

	/* Allocate space for all required complex vectors. */
	lwork = (mit + 1) * (mit + n + 1) + mit;
	v = calloc (lwork, sizeof(complex double));	/* The Krylov subspace. */
	beta = v + n * (mit + 1);		/* The least-squares RHS. */
	h = beta + mit + 1;			/* The upper Hessenberg matrix. */
	s = h + (mit + 1) * mit;		/* Givens rotation sines. */

	/* Allocate space for the Givens rotation cosines. */
	c = malloc (mit * sizeof(double));

	/* Compute the norm of the RHS for residual scaling. */
	rhn = cblas_dznrm2 (n, rhs, 1);

	/* Compute the initial matrix-vector product for the input guess. */
	if (guess) scatmat (v, sol, spl, nsph, trans, shtr);

	/* Subtract from the RHS to form the residual. */
#pragma omp parallel for default(shared) private(j)
	for (j = 0; j < n; ++j) v[j] = rhs[j] - v[j];

	/* Zero the initial guess if one wasn't provided. */
	if (!guess) memset (sol, 0, n * sizeof(complex double));

	/* Find the norm of the initial residual. */
	err = cblas_dznrm2(n, v, 1);

	/* Construct the initial Arnoldi vector by normalizing the residual. */
#pragma omp parallel for default(shared) private(j)
	for (j = 0; j < n; ++j) v[j] /= err;

	/* Construct the vector beta for the minimization problem. */
	beta[0] = err;

	/* Report the RRE. */
	err /= rhn;
	printf ("True residual: %g\n", err);

	for (i = 0; i < mit && err > itc->eps; ++i) {
		/* Point to the working space for this iteration. */
		vp = v + i * n;
		hp = h + i * (mit + 1);

		/* Compute the next expansion of the Krylov space. */
		scatmat (vp + n, vp, spl, nsph, trans, shtr);
		/* Perform modified Gram-Schmidt to orthogonalize the basis. */
		/* This also builds the Hessenberg matrix column. */
		cmgs (vp + n, hp, v, n, i + 1);
		/* Compute the norm of the next basis vector. */
		hp[i + 1] = cblas_dznrm2(n, vp + n, 1);

		/* Avoid breakdown. */
		if (cabs(hp[i + 1]) <  DBL_EPSILON) {
			++i;
			break;
		}

		/* Normalize the basis vector. */
#pragma omp parallel for default(shared) private(j)
		for (j = 0; j < n; ++j) vp[n + j] /= creal(hp[i + 1]);

		/* Apply previous Givens rotations to the Hessenberg column. */
		for (j = 0; j < i; ++j) 
			zrot_ (&one, (void *)(hp + j), &one,
					(void *)(hp + j + 1), &one,
					(void *)(c + j), (void *)(s + j));

		/* Compute the Givens rotation for the current iteration. */
		zlartg_ ((void *)(hp + i), (void *)(hp + i + 1), 
				(void *)(c + i), (void *)(s + i), (void *)(&cr));
		/* Apply the current Givens rotation to the Hessenberg column. */
		hp[i] = cr;
		hp[i + 1] = 0;
		/* Perform the rotation on the vector beta. */
		zrot_ (&one, (void *)(beta + i), &one, 
				(void *)(beta + i + 1), &one, 
				(void *)(c + i), (void *)(s + i));

		/* Estimate the RRE for this iteration. */
		err = cabs(beta[i + 1]) / rhn;
		printf ("GMRES(%d): %g\n", i, err);
	}

	/* If there were any GMRES iterations, update the solution. */
	if (i > 0) {
		/* Compute the minimizer of the least-squares problem. */
		cblas_ztrsv (CblasColMajor, CblasUpper, CblasNoTrans,
				CblasNonUnit, i, h, mit + 1, beta, 1);
		
		/* Compute the update to the solution. */
		cblas_zgemv (CblasColMajor, CblasNoTrans, n, i,
				&cone, v, n, beta, 1, &cone, sol, 1);
	}

	free (v);
	free (c);

	return i;
}

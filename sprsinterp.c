#include <stdlib.h>

#include <math.h>
#include <complex.h>

#include "sprsinterp.h"
#include "util.h"
#include "fsht.h"

/* Compute the matrix-vector product for a sparse matrix. */
int smvmpy (complex double *out, int nr, sprow *rows, complex double *in) {
	int i, j;

	/* Zeroes the output vector beforehand. */
	for (i = 0; i < nr; ++i) {
		out[i] = 0;

		for (j = 0; j < rows[i].nnz; ++j) 
			out[i] += rows[i].val[j] * in[rows[i].idx[j]];
	}

	return nr;
}

/* Compute the matrix-vector product for a transposed sparse matrix. */
int smvtmpy (complex double *out, int nr, sprow *rows, complex double *in) {
	int i, j;

	/* Does not zero the output vector. */
	for (i = 0; i < nr; ++i)
		for (j = 0; j < rows[i].nnz; ++j)
			out[rows[i].idx[j]] += rows[i].val[j] * in[i];

	return nr;
}

/* Build a row of the sparse interpolation matrix. */
int bldsrow (sprow *rw, int ord, int nt, int np,
		double *tf, double *pf, int *ti, int *pi) {
	int i, j, nelt, om1 = ord - 1, ts = 0, te = ord;

	/* The number of non-polar elements. */
	nelt = nt * np;

	rw->nnz = 0;

	/* Deal the north pole if necessary. */
	if (ti[0] < 0) {
		rw->idx[rw->nnz] = nelt;
		rw->val[rw->nnz] = tf[0];
		++(rw->nnz);
		ts = 1; /* Start the grid samples after the polar value. */
	}

	/* Deal with the sourth pole if necessary. */
	if (ti[om1] >= nt) {
		rw->idx[rw->nnz] = nelt + 1;
		rw->val[rw->nnz] = tf[om1];
		++(rw->nnz);
		--te; /* End the grid samples before the polar value. */
	}

	/* Loop through the regular grid values. */
	for (i = ts; i < te; ++i) {
		for (j = 0; j < ord; ++j) {
			/* The location in the input array. */
			nelt = np * ti[i] + pi[j];

			rw->idx[rw->nnz] = nelt;
			rw->val[rw->nnz] = tf[i] * pf[j];
			++(rw->nnz);
		}
	}

	return rw->nnz;
}

/* Build the sparse interpolation matrix. */
int intpmat (sprow *rows, shdata *shout, shdata *shin, int ord) {
	int nrow, ncol, *tidx, *pidx, i, j, k, tint, pint, om1, elused = 0;
	double *tsmp, *psmp, *tc, *pc, dphin, dphout;
	sprow *crow;

	/* Allocate the sample buffers. */
	tidx = malloc (2 * ord * sizeof(int));
	pidx = tidx + ord;

	tsmp = malloc (4 * ord * sizeof(double));
	psmp = tsmp + ord;
	tc = psmp + ord;
	pc = tc + ord;

	/* Allocate the sparse matrix storage block. */
	nrow = shout->ntheta * shout->nphi + 2;
	ncol = shin->ntheta * shin->nphi + 2;

	rows->val = malloc (nrow * ord * ord * sizeof(double));
	rows->idx = malloc (nrow * ord * ord * sizeof(int));

	dphin = 2 * M_PI / MAX(shin->nphi, 1);
	dphout = 2 * M_PI / MAX(shout->nphi, 1);

	om1 = ord - 1;

	/* Build the rows of the matrix. */
	for (i = 0, crow = rows; i < shout->ntheta; ++i) {
		/* Find the near neighbors in theta. */
		tint = interval (shout->theta[i], shin->theta, shin->ntheta);
		tidx[0] = MIN (MAX (tint - ord / 2, -1), shin->ntheta - ord + 1);

		/* The first point, if it is off the end, should be zero. */
		tsmp[0] = (tidx[0] < 0) ? 0 : shin->theta[tidx[0]];

		for (k = 1; k < om1; ++k) {
			tidx[k] = tidx[0] + k;
			tsmp[k] = shin->theta[tidx[k]];
		}

		/* The last point, if it is off the end, should be pi. */
		tidx[om1] = tidx[0] + om1;
		tsmp[om1] = (tidx[om1] >= shin->ntheta) ? M_PI : shin->theta[tidx[om1]];

		/* Build the Legendre coefficients for theta. */
		lgpoly (tc, tsmp, shout->theta[i], ord);

		for (j = 0; j < shout->nphi; ++j, ++crow) {
			/* Find the near neighbors in phi. */
			pint = (j * shin->nphi) / shout->nphi + 1;

			/* Set up the first sample. */
			pidx[0] = pint - ord / 2;
			psmp[0] = dphin * pidx[0];

			/* Don't allow the index to be out-of-bounds. */
			if (pidx[0] < 0) pidx[0] += shin->nphi;
			pidx[0] %= shin->nphi;

			/* Fill out the remaining samples. */
			for (k = 1; k < ord; ++k) {
				pidx[k] = (pidx[0] + k) % shin->nphi;
				psmp[k] = psmp[k - 1] + dphin;
			}

			/* Build the Legendre coefficients for phi. */
			lgpoly (pc, psmp, j * dphout, ord);

			/* Set up the data storage. */
			crow->idx = rows[0].idx + elused;
			crow->val = rows[0].val + elused;

			/* Build the row, and increment the number of used elements. */
			elused += bldsrow (crow, ord, shin->ntheta,
					shin->nphi, tc, pc, tidx, pidx);
		}
	}

	/* Interpolation (or copying) of the north pole. */
	crow->nnz = 1;
	crow->idx = rows[0].idx + elused;
	crow->val = rows[0].val + elused;
	crow->idx[0] = ncol - 2;
	crow->val[0] = 1;
	++elused;

	/* Interpolation (or copying) of the south pole. */
	++crow;
	crow->nnz = 1;
	crow->idx = rows[0].idx + elused;
	crow->val = rows[0].val + elused;
	crow->idx[0] = ncol - 1;
	crow->val[0] = 1;
	++elused;

	/* Tighten the matrix storage requirements. */
	rows->val = realloc (rows->val, elused * sizeof(double));
	rows->idx = realloc (rows->idx, elused * sizeof(int));
	
	free (tsmp);
	free (tidx);

	return elused;
}

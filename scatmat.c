#include <stdlib.h>
#include <complex.h>
#include <string.h>

#include "fastsphere.h"
#include "spreflect.h"
#include "scatmat.h"
#include "fsht.h"

int scatmat (complex double *vout, complex double *vin, spscat *spl,
		int nsph, complex **trans, shdata *shtr) {
	int i, j, k, off, nterm;
	complex double *buf, *voptr, *viptr;

	nterm = shtr->ntheta * shtr->nphi;

	buf = malloc (nterm * sizeof(complex double));

	for (i = 0; i < nsph; ++i) {
		voptr = vout + i * nterm;

		/* Compute the inverse reflection of the local wave. */
		memcpy (voptr, vin + i * nterm, nterm * sizeof(complex double));
		ffsht (voptr, shtr);
		spinvrfl (voptr, voptr, (spl + i)->spdesc->reflect, shtr->deg, shtr->nphi);
		ifsht (voptr, shtr);

		/* Subtract off the translated fields from other spheres. */
		for (j = 0; j < nsph; ++j) {
			if (j == i) continue;
			off = j * nsph + i;
			viptr = vin + i * nterm;

			for (k = 0; k < nterm; ++k) 
				voptr[k] -= trans[off][k] * viptr[k];
		}
	}

	free (buf);

	return nsph;
}

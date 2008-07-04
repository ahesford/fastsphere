#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <string.h>

#include "fastsphere.h"
#include "farfield.h"
#include "util.h"

/* Compute the required order of the root sphere. */
int rootorder (spscat *slist, int nsph, bgtype *bg) {
	int i, l;
	double rad = 0.0, clen;

	for (i = 0; i < nsph; ++i) {
		/* The radius of the center of the sphere. */
		clen = sqrt(DVDOT(slist[i].cen, slist[i].cen));
		/* Make this a radius enclosing all of the current sphere. */
		clen += (slist[i].spdesc)->r;
		/* Find the maxium such radius. */
		rad = MAX(rad, clen);
	}

	/* Use the excess bandwidth formula to find the number of terms. */
	l = exband (bg->k * rad, 6);

	return l;
}

void expnsh (complex double *out, complex double *in, shdata *shout, shdata *shin) {
	int i, j, offo, offi;

	/* Blank the output buffer. */
	memset (out, 0, shout->ntheta * shout->nphi * sizeof(complex double));

	/* Copy the spherical harmonic coefficients into the right place. */
	for (i = 0; i < shin->deg; ++i) {
		offo = i * shout->nphi;
		offi = i * shin->nphi;
		/* Copy the zero-order coefficients. */
		out[offo] = in[offi];

		/* Copy the other coefficients. */
		for (j = 1; j <= i; ++j) {
			out[offo + j] = in[offi + j];
			out[offo + shout->nphi - j] = in[offi + shin->nphi - j];
		}
	}
}

int farfield (complex double *vout, complex double *vin, spscat *slist,
		int nsph, bgtype *bg, shdata *shout, shdata *shin) {
	int ntin, ntout;
	double dphi;

	ntin = shin->ntheta * shin->nphi;
	ntout = shout->ntheta * shout->nphi;
	dphi = 2 * M_PI / MAX(shout->nphi, 1);

	memset (vout, 0, ntout * sizeof(complex double));

#pragma omp parallel default(shared)
{
	int i, j, k, l;
	double s[3], sdc, sth, phi;
	complex double *buf, *vp, sfact;

	buf = malloc (ntout * sizeof(complex double));

#pragma omp for
	for (i = 0; i < nsph; ++i) {
		vp = vin + i * ntin;

		ffsht (vp, shin);	/* SH coefficients for the sphere. */
		expnsh (buf, vp, shout, shin); /* Expand for interpolation. */
		ifsht (vp, shin);	/* Back to angular samples. */
		ifsht (buf, shout);	/* Interpolated angular samples. */

		/* Add the phase-shifted sphere pattern to the total pattern. */
		for (j = 0, l = 0; j < shout->ntheta; ++j) {
			s[2] = cos((shout->theta)[j]);
			sth = sin((shout->theta)[j]);
			for (k = 0; k < shout->nphi; ++k, ++l) {
				phi = k * dphi;
				s[0] = sth * cos(phi);
				s[1] = sth * sin(phi);

				/* Compute the phase-shift factor. */
				sdc = DVDOT(s, slist[i].cen);
				sfact = cexp (-I * bg->k * sdc);
				/* Augment the pattern, with synchronization. */
#pragma omp critical(outrad)
				vout[l] += sfact * buf[l];
			}
		}
	}

	free (buf);
}

	return ntout;
}

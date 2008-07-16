#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <string.h>

#include "fastsphere.h"
#include "farfield.h"
#include "util.h"

/* Compute the required order of the root sphere. */
int rootorder (spscat *slist, int nsph, complex double bgk) {
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
	l = exband (bgk * rad, 6);

	return l;
}

/* Shift and combine outgoing plane waves from small spheres into one large
 * plane-wave pattern. */
int neartofar (complex double *vout, complex double *vin, spscat *slist,
		int nsph, complex double bgk, shdata *shout, shdata *shin) {
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
		/* Expand for interpolation. */
		copysh (buf, shout->deg, shout->nphi, vp, shin->deg, shin->nphi);
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
				sfact = cexp (-I * bgk * sdc);
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

/* Anterpolate and distribute an incoming field to smaller spheres. */
int fartonear (complex double *vout, complex double *vin, spscat *slist,
		int nsph, complex double bgk, shdata *shout, shdata *shin) {
	int ntin, ntout;
	double dphi;
	complex double *buf;

	ntin = shin->ntheta * shin->nphi;
	ntout = shout->ntheta * shout->nphi;
	dphi = 2 * M_PI / MAX(shout->nphi, 1);

	buf = malloc (ntout * sizeof(complex double));
	ffsht (vin, shin);	/* Incoming plane waves to SH coefficients. */
	/* Copy the low-degree coefficients for anterpolation. */
	copysh (buf, shout->deg, shout->nphi, vin, shout->deg, shout->nphi);
	ifsht (vin, shin);	/* Plane waves for the incoming field. */
	ifsht (buf, shout); 	/* Anterpolated angular samples. */

#pragma omp parallel default(shared)
{
	int i, j, k, l;
	double s[3], sdc, sth, phi;
	complex double *vp, sfact;

#pragma omp for
	for (i = 0; i < nsph; ++i) {
		vp = vout + i * ntin;

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
				sfact = cexp (-I * bgk * sdc);
				/* Perform the shift. */
				vp[l] = sfact * buf[l];
			}
		}
	}
}

	free (buf);

	return ntout;
}

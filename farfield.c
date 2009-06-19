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

	ntin = shin->ntheta * shin->nphi + 2;
	ntout = shout->ntheta * shout->nphi + 2;
	dphi = 2 * M_PI / MAX(shout->nphi, 1);

	memset (vout, 0, ntout * sizeof(complex double));

#pragma omp parallel default(shared)
{
	int i, j, k, l;
	double s[3], sdc, sth, phi, dtheta, theta;
	complex double *vp, sfact, *buf;

	buf = malloc (ntout * sizeof(complex double));
	dtheta = M_PI / shout->ntheta;

#pragma omp for
	for (i = 0; i < nsph; ++i) {
		vp = vin + i * ntin;

		/* Interpolate the spherical scattered field. */
		ffsht (vp, shin);
		memset (buf, 0, ntout * sizeof(complex double));
		copysh (shin->deg, buf, shout->nphi, vp, shin->nphi);
		ifsht (vp, shin);
		ifsht (buf, shout);

		/* Add the phase-shifted sphere pattern to the total pattern. */
		for (j = 0, l = 0, theta = dtheta / 2; j < shout->ntheta; ++j, theta += dtheta) {
			s[2] = cos(theta);
			sth = sin(theta);
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

		/* Augment the poles with synchronization. The offset should
		 * be properly computed from the for loop. */
		l = shout->ntheta * shout->nphi;
#pragma omp critical(outrad)
{
		vout[l] += cexp(-I * bgk * slist[i].cen[2]) * buf[l];
		vout[l + 1] += cexp (I * bgk * slist[i].cen[2]) * buf[l + 1];
}
	}

	free (buf);
}

	return ntout;
}

/* Anterpolate and distribute an incoming field to smaller spheres. Input is
 * a plane-wave expansion, output is plane-wave expansion. */
int fartonear (complex double *vout, complex double *vin, spscat *slist,
		int nsph, complex double bgk, shdata *shout, shdata *shin) {
	int ntin, ntout;
	double dphi;

	ntin = shin->ntheta * shin->nphi + 2;
	ntout = shout->ntheta * shout->nphi + 2;
	dphi = 2 * M_PI / MAX(shin->nphi, 1);

#pragma omp parallel default(shared)
{
	int i, j, k, l;
	double s[3], sdc, sth, phi, theta, dtheta;
	complex double *vp, *buf;

	buf = malloc (ntin * sizeof(complex double));
	dtheta = M_PI / shin->ntheta;

#pragma omp for
	for (i = 0; i < nsph; ++i) {
		vp = vout + i * ntout;

		/* Shift the phase of the sphere pattern. */
		for (j = 0, l = 0, theta = dtheta / 2; j < shin->ntheta; ++j, theta += dtheta) {
			s[2] = cos(theta);
			sth = sin(theta);
			for (k = 0; k < shin->nphi; ++k, ++l) {
				phi = k * dphi;
				s[0] = sth * cos(phi);
				s[1] = sth * sin(phi);

				/* Compute the phase-shift factor. */
				sdc = DVDOT(s, slist[i].cen);
				buf[l] = vin[l] * cexp (I * bgk * sdc);
			}
		}
		
		/* The poles don't contribute to the forward SHT, and they
		 * are re-evaluated by the inverse SHT. Don't compute them. */
		ffsht (buf, shin);
		memset (vp, 0, ntout * sizeof(complex double));
		copysh (shout->deg, vp, shout->nphi, buf, shin->nphi);
		ifsht (vp, shout);
	}

	free (buf);
}

	return ntout;
}

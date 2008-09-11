#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <string.h>

#include "fastsphere.h"
#include "farfield.h"
#include "util.h"

#include "sprsinterp.h"

int writebvec (FILE *, complex double *, int, int);

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
int neartofar (complex double *vout, complex double *vin, spscat *slist, int nsph,
		complex double bgk, shdata *shout, shdata *shin, sprow *imat) {
	int ntin, ntout;
	double dphi;

	ntin = shin->ntheta * shin->nphi + 2;
	ntout = shout->ntheta * shout->nphi + 2;
	dphi = 2 * M_PI / MAX(shout->nphi, 1);

	memset (vout, 0, ntout * sizeof(complex double));

#pragma omp parallel default(shared)
{
	int i, j, k, l;
	double s[3], sdc, sth, phi;
	complex double *vp, sfact, *buf;

	buf = malloc (ntout * sizeof(complex double));

#pragma omp for
	for (i = 0; i < nsph; ++i) {
		vp = vin + i * ntin;

		/* Interpolate the spherical scattered field. The output
		 * buffer is first zeroed by the multiplication. */
		smvmpy (buf, ntout, imat, vp);

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
int fartonear (complex double *vout, complex double *vin, spscat *slist, int nsph,
		complex double bgk, shdata *shout, shdata *shin, sprow *imat) {
	int ntin, ntout;
	double dphi;

	ntin = shin->ntheta * shin->nphi + 2;
	ntout = shout->ntheta * shout->nphi + 2;
	dphi = 2 * M_PI / MAX(shin->nphi, 1);

#pragma omp parallel default(shared)
{
	int i, j, k, l;
	double s[3], sdc, sth, phi;
	complex double *vp, *buf;

	buf = malloc (ntin * sizeof(complex double));

#pragma omp for
	for (i = 0; i < nsph; ++i) {
		vp = vout + i * ntout;

		/* Shift the phase of the sphere pattern. */
		for (j = 0, l = 0; j < shin->ntheta; ++j) {
			s[2] = cos((shin->theta)[j]);
			sth = sin((shin->theta)[j]);
			for (k = 0; k < shin->nphi; ++k, ++l) {
				phi = k * dphi;
				s[0] = sth * cos(phi);
				s[1] = sth * sin(phi);

				/* Compute the phase-shift factor. */
				sdc = DVDOT(s, slist[i].cen);
				buf[l] = vin[l] * cexp (I * bgk * sdc);
			}
		}
		
		/* Anterpolate the result. */
		ffsht (buf, shin);
		memset (vp, 0, ntout * sizeof(complex double));
		copysh (shout->deg, vp, shout->nphi, buf, shin->nphi);
		ifsht (vp, shout);
	}

	free (buf);
}

	return ntout;
}

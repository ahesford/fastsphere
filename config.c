#include <stdio.h>
#include <stdlib.h>

#include <complex.h>
#include <math.h>
#include <float.h>

#include <string.h>

#include "fastsphere.h"
#include "scatmat.h"
#include "util.h"
#include "config.h"

#define BUFLEN 1024

int nextline (FILE *input, char *buf, int maxlen) {
	int slen;
	char *cptr;

	/* Just keep reading lines until we find something we want. */
	do {
		/* Pull out the line, give up on errors. */
		if (!fgets (buf, maxlen, input)) return 0;

		/* Skip any leading tabs or spaces. */
		/* The NUL termination occurs before running off the end. */
		cptr = buf;
		while (*cptr == ' ' || *cptr == '\t' || *cptr == '\n') ++cptr;

		/* This line is either a comment or is empty. */
		if (*cptr == '#' || *cptr == '\0') continue;

		/* Remove leading whitespace. */
		slen = strlen (cptr) + 1;
		memmove (buf, cptr, slen * sizeof(char));

		/* Strings with only whitespace have already been skipped.
		 * There is no risk of running over the beginning. */
		cptr = buf + slen - 2;
		while (*cptr == ' ' || *cptr == '\t' || *cptr == '\n') {
			*cptr = '\0';
			--cptr;
		}

		return 1;
	} while (1);

	return 0;
}

int readcfg (FILE *cfgin, int *nspheres, int *nsptype, sptype **spt, sptype *encl,
		spscat **spl, bgtype *bg, exctparm *exct, itconf *itc, int *ndig) {
	int i, tp;
	double rb, ib, *loc, *ax, tracc = -1.0;
	char buf[BUFLEN];
	sptype *stptr;
	spscat *ssptr;

	if (!nextline (cfgin, buf, BUFLEN)) return 0;

	/* Background sound speed, attenuation, density. */
	if (sscanf (buf, "%lf %lf %lf", &(bg->cabs), &(bg->alpha),
				&(bg->rho0)) != 3) return 0;

	/* Convert attenuation from dB per cm * MHz to dB per wavelength. */
	bg->alpha *= 1e-4 * bg->cabs;

	/* Build the background wave number. The relative sound speed is unity. */
	bg->k = wavenum (1.0, bg->alpha);

	if (!nextline (cfgin, buf, BUFLEN)) return 0;

	/* The excitation frequency and source counts. */
	exct->npw = exct->nps = 0;
	if (sscanf (buf, "%lf %d %d", &(exct->f), &(exct->npw), &(exct->nps)) < 2)
		return 0;

	/* Convert excitation frequency into Hz. */
	exct->f *= 1e6;

	if (exct->npw > 0) {
		/* Allocate space for plane wave magnitudes and directions. */
		exct->pwmag = malloc (exct->npw * sizeof(complex double));
		exct->theta = malloc (2 * exct->npw * sizeof(double));
		exct->phi = exct->theta + exct->npw;
	} else {
		exct->pwmag = NULL;
		exct->theta = exct->phi = NULL;
	}

	if (exct->nps > 0) {
		/* Allocate space for point source configurations. */
		exct->psmag = malloc (exct->nps * sizeof(complex double));
		exct->psloc = malloc (7 * exct->nps * sizeof(double));
		exct->psax = exct->psloc + 3 * exct->nps;
		exct->alpha = exct->psax + 3 * exct->nps;
	} else {
		exct->psmag = NULL;
		exct->psloc = exct->psax = exct->alpha = NULL;
	}

	/* Read the plane-wave excitation parameters. */
	for (i = 0; i < exct->npw; ++i) {
		if (!nextline (cfgin, buf, BUFLEN)) return 0;

		/* Excitation parameters. */
		if (sscanf (buf, "%lf %lf %lf %lf", &rb, &ib,
					exct->theta + i, exct->phi + i) != 4)
			return 0;

		/* Assemble the magnitude. */
		(exct->pwmag)[i] = rb + I * ib;

		/* Convert the angles to radians. */
		(exct->theta)[i] = DEG2RAD((exct->theta)[i]);
		(exct->phi)[i] = DEG2RAD((exct->phi)[i]);
	}

	/* Read the point-source excitation parameters. */
	for (i = 0; i < exct->nps; ++i) {
		if (!nextline (cfgin, buf, BUFLEN)) return 0;

		loc = exct->psloc + 3 * i;
		ax = exct->psax + 3 * i;

		/* Excitation parameters. */
		if (sscanf (buf, "%lf %lf %lf %lf %lf %lf %lf %lf %lf",
					&rb, &ib, loc, loc + 1, loc + 2,
					ax, ax + 1, ax + 2, exct->alpha + i) != 9)
			return 0;

		/* Assemble the magnitude. */
		(exct->psmag)[i] = rb + I * ib;

		/* Convert coordinates to wavelengths. */
		loc[0] *= exct->f / bg->cabs;
		loc[1] *= exct->f / bg->cabs;
		loc[2] *= exct->f / bg->cabs;

		/* Ensure directivity axis has unit length. */
		rb = MAG(ax[0], ax[1], ax[2]);
		if (rb > DBL_EPSILON) {
			ax[0] /= rb;
			ax[1] /= rb;
			ax[2] /= rb;
		} else ax[0] = ax[1] = ax[2] = 0.0;
	}
	
	if (!nextline (cfgin, buf, BUFLEN)) return 0;

	/* The number of iterations, restart and solver tolerance. */
	if (sscanf (buf, "%d %d %lf", &(itc->iter), &(itc->restart), &(itc->eps)) != 3)
		return 0;

	if (!nextline (cfgin, buf, BUFLEN)) return 0;

	/* Read about an enclosing sphere, if there is one. */
	if (encl != NULL) {
		/* The enclosing sphere, if there is any. */
		if (sscanf (buf, "%lf %lf %lf %lf", &(encl->r), &(encl->c),
					&(encl->alpha), &(encl->rho)) != 4 && encl->r > 0)
			return 0;
		
		encl->r *= exct->f / bg->cabs; /* Radius to wavelengths. */
		encl->c /= bg->cabs; /* Sound speed to relative sound speed. */
		encl->alpha *= 1e-4 * bg->cabs; /* Attenuation to dB per wavelength. */
		encl->rho /= bg->rho0; /* Density to relative density. */
		
		/* Build the wave number for the enclosing sphere. */
		encl->k = wavenum (encl->c, encl->alpha);
		
		if (!nextline (cfgin, buf, BUFLEN)) return 0;
	}

	/* The number of spheres, and the number of unique sphere types. 
	 * The translation and harmonic accuracy argument is optional. */
	if (sscanf (buf, "%d %d %lf", nspheres, nsptype, &tracc) < 2) return 0;

	/* Compute the desired accuracy or use the default. */
	if (tracc > 0) *ndig = MAX(MIN_ACC_DIG, (int)(-floor(log10(tracc))));
	else *ndig = DEF_ACC_DIG;

	/* Allocate the sphere list and type list. */
	*spt = malloc (*nsptype * sizeof (sptype));
	*spl = malloc (*nspheres * sizeof (spscat));

	for (i = 0, stptr = *spt; i < *nsptype; ++i, ++stptr) {
		if (!nextline (cfgin, buf, BUFLEN)) return 0;

		/* By default, the shear wave speed is 0. */
		stptr->csh = 0;

		/* The radius, sound speed, attenuation, density and
		 * optional shear wave speed  of each sphere type. */
		if (sscanf (buf, "%lf %lf %lf %lf %lf", &(stptr->r),
					&(stptr->c), &(stptr->alpha),
					&(stptr->rho), &(stptr->csh)) < 4)
			return 0;

		/* Convert radius to wavelengths. */
		stptr->r *= exct->f / bg->cabs;

		/* Convert wave speeds to relative speeds. */
		stptr->c /= bg->cabs;
		stptr->csh /= bg->cabs;

		/* Convert attenuation to dB per wavelength. */
		stptr->alpha *= 1e-4 * bg->cabs;

		/* Build the wave number. */
		stptr->k = wavenum (stptr->c, stptr->alpha);

		/* Convert density to relative density. */
		stptr->rho /= bg->rho0;
	}

	for (i = 0, ssptr = *spl; i < *nspheres; ++i, ++ssptr) {
		if (!nextline (cfgin, buf, BUFLEN)) return 0;

		/* The type and center coordinates of each sphere. */
		if (sscanf (buf, "%d %lf %lf %lf", &tp, ssptr->cen,
					ssptr->cen + 1, ssptr->cen + 2) != 4)
			return 0;

		ssptr->spdesc = (*spt) + tp;

		/* Convert coordinates to wavelengths. */
		(ssptr->cen)[0] *= exct->f / bg->cabs;
		(ssptr->cen)[1] *= exct->f / bg->cabs;
		(ssptr->cen)[2] *= exct->f / bg->cabs;
	}

	return *nspheres;
}

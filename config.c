#include <stdio.h>
#include <stdlib.h>

#include <complex.h>
#include <math.h>

#include <string.h>

#include "fastsphere.h"
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

int readcfg (FILE *cfgin, int *nspheres, int *nsptype, sptype **spt,
		spscat **spl, bgtype *bg, exctparm *exct) {
	int i, tp;
	char buf[BUFLEN];
	sptype *stptr;
	spscat *ssptr;

	if (!nextline (cfgin, buf, BUFLEN)) return 0;

	/* Background sound speed and attenuation. */
	if (sscanf (buf, "%lf %lf", &(bg->cabs), &(bg->alpha)) != 2) return 0;

	/* Convert attenuation from dB per cm * MHz to dB per wavelength. */
	bg->alpha *= 1e-4 * bg->cabs;

	if (!nextline (cfgin, buf, BUFLEN)) return 0;
	
	/* The excitation frequency and direction (theta and phi, in degrees). */
	if (sscanf (buf, "%lf %lf %lf %lf", &(exct->f), exct->cen,
				exct->cen + 1, exct->cen + 2) != 4)
		return 0;

	/* Convert excitation frequency into Hz. */
	exct->f *= 1e6;

	if (!nextline (cfgin, buf, BUFLEN)) return 0;

	/* The number of spheres, and the number of unique sphere types. */
	if (sscanf (buf, "%d %d", nspheres, nsptype) != 2) return 0;

	/* Allocate the sphere list and type list. */
	*spt = malloc (*nsptype * sizeof (sptype));
	*spl = malloc (*nspheres * sizeof (spscat));

	for (i = 0, stptr = *spt; i < *nsptype; ++i, ++stptr) {
		if (!nextline (cfgin, buf, BUFLEN)) return 0;

		/* The radius, sound speed and attenuation of each sphere type. */
		if (sscanf (buf, "%lf %lf %lf", &(stptr->r), &(stptr->c),
					&(stptr->alpha)) != 3)
			return 0;

		/* Convert radius to wavelengths. */
		stptr->r *= exct->f / bg->cabs;

		/* Convert sound speed to relative sound speed. */
		stptr->c /= bg->cabs;

		/* Convert attenuation to dB per wavelength. */
		stptr->alpha *= 1e-4 * bg->cabs;
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

	/* Build the background wave number. The relative sound speed is unity. */
	bg->k = buildkvec (1.0, bg->alpha);

	/* Build the wave numbers for each type of sphere. */
	for (i = 0, stptr = *spt; i < *nsptype; ++i, ++stptr)
		stptr->k = buildkvec (stptr->c, stptr->alpha);

	return *nspheres;
}

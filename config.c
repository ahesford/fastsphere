#include <stdio.h>
#include <stdlib.h>

#include <complex.h>
#include <math.h>

#include <string.h>

#include "fastsphere.h"
#include "util.h"

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

int readcfg (FILE *cfgin, sptype **spt, spscat **spl, bgtype *bg, exctparm *exct) {
	int nsptype = 0, nspheres = 0, i, tp;
	char buf[BUFLEN];

	if (!nextline (cfgin, buf, BUFLEN)) return 0;

	/* Background sound speed and attenuation. */
	if (sscanf (buf, "%lf %lf", &(bg->c), &(bg->alpha)) != 2) return 0;

	if (!nextline (cfgin, buf, BUFLEN)) return 0;

	/* The number of spheres, and the number of unique sphere types. */
	if (sscanf (buf, "%d %d", &nspheres, &nsptype) != 2) return 0;

	/* Allocate the sphere list and type list. */
	*spt = malloc (nsptype * sizeof (sptype));
	*spl = malloc (nspheres * sizeof (spscat));

	for (i = 0; i < nsptype; ++i) {
		if (!nextline (cfgin, buf, BUFLEN)) return 0;

		/* The radius, sound speed and attenuation of each sphere. */
		if (sscanf (buf, "%lf %lf %lf", &((*spt)[i].r),
				&((*spt)[i].c), &((*spt)[i].alpha)) != 3)
			return 0;
	}

	for (i = 0; i < nspheres; ++i) {
		if (!nextline (cfgin, buf, BUFLEN)) return 0;

		/* The type and center coordinates of each sphere. */
		if (sscanf (buf, "%d %lf %lf %lf", &tp, (*spl)[i].cen,
				(*spl)[i].cen + 1, (*spl)[i].cen + 2) != 4)
			return 0;
	}

	if (!nextline (cfgin, buf, BUFLEN)) return 0;
	
	/* The excitation frequency and direction (theta and phi, in degrees). */
	if (sscanf (buf, "%lf %lf %lf", &(exct->f), &(exct->theta), &(exct->phi)) != 3)
		return 0;

	/* Convert the excitation angles to radians. */
	exct->theta *= M_PI / 180.0;
	exct->phi *= M_PI / 180;

	/* Build the background wave number. */
	bg->k = buildkvec (bg->c, bg->alpha, exct->f);

	/* Build the wave numbers for each type of sphere. */
	for (i = 0; i < nsptype; ++i)
		(*spt)[i].k = buildkvec ((*spt)[i].c, (*spt)[i].alpha, exct->f);

	return nspheres;
}

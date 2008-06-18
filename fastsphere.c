#include <stdio.h>
#include <stdlib.h>

#include <complex.h>
#include <math.h>

#include "fastsphere.h"
#include "config.h"
#include "init.h"
#include "util.h"
#include "scatmat.h"

int main (int argc, char **argv) {
	int nspheres, nsptype, n, i, nterm, trunc;
	sptype *sparms;
	spscat *slist;
	bgtype bg;
	exctparm exct;
	shdata shtr;
	itconf itc;
	complex double **trans, *rhs, *sol;

	readcfg (stdin, &nspheres, &nsptype, &sparms, &slist, &bg, &exct, &itc);
	fprintf (stderr, "Parsed configuration for %d spheres at %g MHz\n", nspheres, exct.f / 1e6);
	sphinit (sparms, nsptype, &bg, &shtr);
	fprintf (stderr, "Initialized spherical harmonic data for degree %d\n", shtr.deg);
	sphbldfmm (&trans, slist, nspheres, &bg, &shtr);
	fprintf (stderr, "Built FMM translators for all spheres\n");

	nterm = shtr.ntheta * shtr.nphi;
	n = nspheres * nterm;
	rhs = calloc (2 * n, sizeof(complex double));
	sol = rhs + n;
	
	/* Build a translator from the source location to each sphere.
	 * This is the incoming incident field for each sphere. */
	trunc = 2 * shtr.deg - 1;

#pragma omp parallel private(i) default(shared)
{
	/* These variables are private. */
	double sdir[3], rsrc;
	complex double kr;

	for (i = 0; i < nspheres; ++i) {
		sdir[0] = slist[i].cen[0] - exct.cen[0];
		sdir[1] = slist[i].cen[1] - exct.cen[1];
		sdir[2] = slist[i].cen[2] - exct.cen[2];

		rsrc = sqrt (DVDOT(sdir,sdir));
		sdir[0] /= rsrc; sdir[1] /= rsrc; sdir[2] /= rsrc;
		kr = bg.k * rsrc;

		translator (rhs + i * nterm, trunc, shtr.ntheta,
				shtr.nphi, shtr.theta, kr, sdir);
	}
}

	/* Convert the incoming incident field to the reflected incident
	 * field, which is the RHS for this problem. */
	buildrhs (rhs, slist, nspheres, &shtr);
	fprintf (stderr, "Built RHS vector\n");

	itsolve (sol, rhs, slist, nspheres, trans, &shtr, &itc);
	fprintf (stderr, "Iteration complete, writing solution\n");

	fflush (stdout);

	for (i = 0; i < n; ++i)
		printf ("%20.15f %20.15f\n", creal (sol[i]), cimag (sol[i]));

	fflush (stdout);

	clrspheres (sparms, nsptype);
	free (*trans);
	free (trans);
	free (rhs);

	return EXIT_SUCCESS;
}

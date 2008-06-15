#include <stdio.h>
#include <stdlib.h>

#include <complex.h>
#include <math.h>

#include "fastsphere.h"
#include "config.h"
#include "init.h"
#include "scatmat.h"

int main (int argc, char **argv) {
	int nspheres, nsptype, n, i, j, nterms;
	sptype *sparms;
	spscat *slist;
	bgtype bg;
	exctparm exct;
	shdata shtr;
	itconf itc;
	complex double **trans, *rhs, *sol, val;
	double sdir[3], phase;

	readcfg (stdin, &nspheres, &nsptype, &sparms, &slist, &bg, &exct);
	sphinit (sparms, nsptype, &bg, &shtr);
	sphbldfmm (&trans, slist, nspheres, &bg, &shtr);

	itc.iter = 100; itc.restart = 10; itc.eps = 1e-6;

	nterms = shtr.ntheta * shtr.nphi;
	n = nspheres * nterms;
	rhs = calloc (2 * n, sizeof(complex double));
	sol = rhs + n;
	
	sdir[0] = cos(exct.phi) * sin(exct.theta);
	sdir[1] = sin(exct.phi) * sin(exct.theta);
	sdir[2] = cos(exct.theta);

	/* Set the excitation for each sphere. */
	for (j = 0; j < nspheres; ++j) {
		phase = sdir[0] * slist[j].cen[0] + sdir[1] * slist[j].cen[1]
			+ sdir[2] * slist[j].cen[2];
		val = cexp (-I * phase) / (double)shtr.nphi;
		for (i = 0; i < shtr.nphi; ++i)
			rhs[i] = val;
	}

	buildrhs (rhs, slist, nspheres, &shtr);

	itsolve (sol, rhs, slist, nspheres, trans, &shtr, &itc);

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

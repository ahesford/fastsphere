#include <stdio.h>
#include <stdlib.h>

#include "fastsphere.h"
#include "config.h"
#include "init.h"
#include "scatmat.h"

int main (int argc, char **argv) {
	int nspheres, nsptype, n, i;
	sptype *sparms;
	spscat *slist;
	bgtype bg;
	exctparm exct;
	shdata shtr;
	itconf itc;
	complex double **trans, *rhs, *sol;

	readcfg (stdin, &nspheres, &nsptype, &sparms, &slist, &bg, &exct);
	sphinit (sparms, nsptype, &bg, &shtr);
	sphbldfmm (&trans, slist, nspheres, &bg, &shtr);

	itc.iter = 100; itc.restart = 10; itc.eps = 1e-3;

	n = nspheres * shtr.ntheta * shtr.nphi;
	rhs = calloc (2 * n, sizeof(complex double));
	sol = rhs + n;

	for (i = 0; i < shtr.nphi; ++i)
		rhs[i] = 1.0 / (double)shtr.nphi;

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

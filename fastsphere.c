#include <stdio.h>
#include <stdlib.h>

#include "fastsphere.h"
#include "config.h"
#include "init.h"

int main (int argc, char **argv) {
	int nspheres, nsptype;
	sptype *sparms;
	spscat *slist;
	bgtype bg;
	exctparm exct;
	shdata shtr;
	complex double **trans;

	readcfg (stdin, &nspheres, &nsptype, &sparms, &slist, &bg, &exct);
	sphinit (sparms, nsptype, &bg, &shtr);
	sphbldfmm (&trans, slist, nspheres, &bg, &shtr);

	clrspheres (sparms, nsptype);
	free (*trans);
	free (trans);

	return EXIT_SUCCESS;
}

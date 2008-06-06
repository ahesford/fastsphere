#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <complex.h>
#include <math.h>

#include <fftw3.h>

#include "fsht.h"

void prtdata (FILE *output, complex double *vals, int nth, int nph) {
	int i, j;
	complex double *vptr;

	for (i = 0, vptr = vals; i < nth; ++i) 
		for (j = 0; j < nph; ++j, ++vptr)
			fprintf (output, "%d %d %20.15g %20.15g\n", i, j,
					creal(*vptr), cimag(*vptr));
}

void readat (FILE *input, complex double *vals, int n) {
	double dr, di;
	int i;

	for (i = 0; i < n; ++i) {
		fscanf (input, "%lf %lf", &dr, &di);
		vals[i] = dr + I * di;
	}
}

double vnorm (complex double *vals, int n) {
	double nrm, v;
	int i;

	for (i = 0, nrm = 0; i < n; ++i) {
		v = cabs (vals[i]);
		nrm += v * v;
	}

	return sqrt (nrm);
}

int main (int argc, char **argv) {
	int l = 10, i, j;
	shdata dat;
	complex double *vals, *ref;

	if (argc > 1) l = atoi (argv[1]);

	fshtinit (&dat, l);

	j = dat.ntheta * dat.nphi;

	vals = malloc (2 * j * sizeof(complex double));
	ref = vals + j;

	if (!getenv ("SH_UNITY")) readat (stdin, vals, j);
	else {
		for (i = 0; i < j; ++i) vals[i] = 1.0;
	}

	/* Store a reference value. */
	for (i = 0; i < j; ++i) ref[i] = vals[i];

	ffsht (vals, &dat);

	if (!getenv ("SH_QUIET")) {
		printf ("FORWARD\n");
		prtdata (stdout, vals, dat.ntheta, dat.nphi);
	}

	ifsht (vals, &dat);

	if (!getenv ("SH_QUIET")) {
		printf ("INVERSE\n");
		prtdata (stdout, vals, dat.ntheta, dat.nphi);
	}

	j = dat.ntheta * dat.nphi;
	
	for (i = 0; i < j; ++i) vals[i] -= ref[i];

	printf ("NORMALIZED ERROR: %20.15g\n", vnorm (vals, j) / vnorm (ref, j));

	fshtfree (&dat);
	free (vals);

	return EXIT_SUCCESS;
}

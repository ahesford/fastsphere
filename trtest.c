#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <complex.h>
#include <math.h>

#include <fftw3.h>

#include "fsht.h"

void prtdata (FILE *output, complex double *vals, int nth, int nph, double *theta) {
	int i, j;
	double phi;
	complex double *vptr;

	/* Print indices if there is no theta array. */
	if (!theta) {
		for (i = 0, vptr = vals; i < nth; ++i)
			for (j = 0; j < nph; ++j, ++vptr)
				fprintf (output, "%d %d %20.15g %20.15g\n", i, j,
					creal(*vptr), cimag(*vptr));
		return;
	}

	for (i = 0, vptr = vals; i < nth; ++i) 
		for (j = 0; j < nph; ++j, ++vptr) {
			phi = 2 * M_PI * j / nph;
			fprintf (output, "%g %g %20.15g %20.15g\n", theta[i], phi,
					creal(*vptr), cimag(*vptr));
		}
}

void readat (FILE *input, complex double *vals, int n) {
	double dr, di;
	int i;

	for (i = 0; i < n; ++i) {
		fscanf (input, "%lf %lf", &dr, &di);
		vals[i] = dr + I * di;
	}
}

int main (int argc, char **argv) {
	int l = 10, i, j;
	shdata dat;
	complex double *vals, *trans, kr;
	double sdir[3] = { 0.0, 0.0, 1.0 }, trlen = 10;

	/* Read some parameters. */
	if (argc > 1) l = atoi (argv[1]);
	if (argc > 2) trlen = strtod (argv[2], NULL);

	fshtinit (&dat, l);

	j = dat.ntheta * dat.nphi;

	vals = malloc (2 * j * sizeof(complex double));
	trans = vals + j;

	kr = 2 * M_PI * trlen;

	/* Build the FMM translator in the z-direction. */
	translator (trans, dat.deg, dat.ntheta, dat.nphi, dat.theta, kr, sdir);

	/* Read the SH coefficients. */
	readat (stdin, vals, j);

	/* Translate from SH coefficients to angular samples. */
	ifsht (vals, &dat);

	if (!getenv ("SH_QUIET")) {
		printf ("ANGULAR SAMPLES\n");
		prtdata (stdout, vals, dat.ntheta, dat.nphi, dat.theta);
	}

	/* Perform the translation. */
	for (i = 0; i < j; ++i) vals[i] *= trans[i];

	if (!getenv ("SH_QUIET")) {
		printf ("ANGULAR SAMPLES (INCOMING)\n");
		prtdata (stdout, vals, dat.ntheta, dat.nphi, dat.theta);
	}

	ffsht (vals, &dat);

	if (!getenv ("SH_QUIET")) {
		printf ("SPHERICAL HARMONICS (INCOMING)\n");
		prtdata (stdout, vals, dat.ntheta, dat.nphi, NULL);
	}

	fshtfree (&dat);
	free (vals);

	return EXIT_SUCCESS;
}

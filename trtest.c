#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <complex.h>
#include <math.h>

#include <fftw3.h>

#include "fsht.h"
#include "translator.h"

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
	int l, i, j, m, n;
	shdata dat;
	complex double *vals, *trans, kr, ka;
	double sdir[3] = { 0.0, 0.0, 1.0 }, trlen = 10;
	FILE *output;
	char fname[1024], *ffmt;

	if (argc < 4) return EXIT_FAILURE;

	/* Read some parameters. */
	ka = 2 * M_PI * strtod (argv[1], NULL);
	trlen = strtod (argv[2], NULL);
	ffmt = argv[3];

	l = exband (ka, 6);

	fprintf (stderr, "Translation of %d harmonics.\n", l);

	fshtinit (&dat, l, 2 * l - 1, 2 * l - 1);

	j = dat.ntheta * dat.nphi;

	vals = malloc (2 * j * sizeof(complex double));
	trans = vals + j;

	kr = 2 * M_PI * trlen;

	/* Build the FMM translator in the z-direction. */
	translator (trans, 2 * l - 1, dat.ntheta, dat.nphi, dat.theta, kr, sdir);

	for (n = 0; n < l; ++n) {
		fprintf (stderr, "Translating (%d,0)\n", n);
		memset (vals, 0, j * sizeof(complex double));
		vals[n * dat.nphi] = 1.0;
		
		/* Scale the coefficients in preparation for an inverse transform. */
		shscale (vals, &dat, -1);
		
		/* Translate from SH coefficients to angular samples. */
		ifsht (vals, &dat);
		
		/* Perform the translation. */
		for (i = 0; i < j; ++i) vals[i] *= trans[i]; 
		
		/* Translate from angular samples back to SH coefficients. */
		ffsht (vals, &dat);
		
		/* Scale the coefficients to get the appropriate values. */
		shscale (vals, &dat, 1);

		sprintf (fname, ffmt, n, 0);
		output = fopen (fname, "w");
		prtdata (output, vals, dat.deg, dat.nphi, NULL);
		fclose (output);

		for (m = 1; m <= n; ++m) {
			fprintf (stderr, "Translating (%d,%d)\n", n, m);
			memset (vals, 0, j * sizeof(complex double));
			vals[n * dat.nphi + m] = 1.0;
			
			/* Scale the coefficients in preparation for an inverse transform. */
			shscale (vals, &dat, -1);
			
			/* Translate from SH coefficients to angular samples. */
			ifsht (vals, &dat);
			
			/* Perform the translation. */
			for (i = 0; i < j; ++i) vals[i] *= trans[i];
			
			/* Translate from angular samples back to SH coefficients. */
			ffsht (vals, &dat);
			
			/* Scale the coefficients to get the appropriate values. */
			shscale (vals, &dat, 1);
			sprintf (fname, ffmt, n, m);
			output = fopen (fname, "w");
			prtdata (output, vals, dat.deg, dat.nphi, NULL);
			fclose (output);

			fprintf (stderr, "Translating (%d,%d)\n", n, -m);
			memset (vals, 0, j * sizeof(complex double));
			vals[(n + 1) * dat.nphi - m] = 1.0;
			/* Scale the coefficients in preparation for an inverse transform. */
			shscale (vals, &dat, -1);
			
			/* Translate from SH coefficients to angular samples. */
			ifsht (vals, &dat);
			
			/* Perform the translation. */
			for (i = 0; i < j; ++i) vals[i] *= trans[i];
			
			/* Translate from angular samples back to SH coefficients. */
			ffsht (vals, &dat);
			
			/* Scale the coefficients to get the appropriate values. */
			shscale (vals, &dat, 1);
			sprintf (fname, ffmt, n, -m);
			output = fopen (fname, "w");
			prtdata (output, vals, dat.deg, dat.nphi, NULL);
			fclose (output);
		}
	}


	fshtfree (&dat);
	free (vals);

	return EXIT_SUCCESS;
}

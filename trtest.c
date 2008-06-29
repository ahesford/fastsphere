#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <complex.h>
#include <math.h>

#include "util.h"
#include "shtranslate.h"

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
	int l, lda, i, j, m, n;
	complex double *vals, kr, ka;
	double trlen = 10;
	FILE *output;
	char fname[1024], *ffmt;

	if (argc < 4) return EXIT_FAILURE;

	/* Read some parameters. */
	ka = 2 * M_PI * strtod (argv[1], NULL);
	trlen = strtod (argv[2], NULL);
	ffmt = argv[3];

	l = exband (ka, 6);
	lda = 2 * l - 1;

	fprintf (stderr, "Translation of %d harmonics.\n", l);

	j = l * lda;

	vals = malloc (j * sizeof(complex double));

	kr = 2 * M_PI * trlen;

	for (n = 0; n < l; ++n) {
		fprintf (stderr, "Translating (%d,0)\n", n);
		memset (vals, 0, j * sizeof(complex double));
		vals[ELT(0,n,lda)] = 1.0;

		shtranslate (vals, l, lda, kr);

		sprintf (fname, ffmt, n, 0);
		output = fopen (fname, "w");
		prtdata (output, vals, l, lda, NULL);
		fclose (output);

		for (m = 1; m <= n; ++m) {
			fprintf (stderr, "Translating (%d,%d)\n", n, m);
			memset (vals, 0, j * sizeof(complex double));
			vals[ELT(m,n,lda)] = 1.0;

			shtranslate (vals, l, lda, kr);
			
			sprintf (fname, ffmt, n, m);
			output = fopen (fname, "w");
			prtdata (output, vals, l, lda, NULL);
			fclose (output);

			fprintf (stderr, "Translating (%d,%d)\n", n, -m);
			memset (vals, 0, j * sizeof(complex double));
			vals[ELT(-m,n+1,lda)] = 1.0;

			shtranslate (vals, l, lda, kr);

			sprintf (fname, ffmt, n, -m);
			output = fopen (fname, "w");
			prtdata (output, vals, l, lda, NULL);
			fclose (output);
		}
	}

	free (vals);

	return EXIT_SUCCESS;
}

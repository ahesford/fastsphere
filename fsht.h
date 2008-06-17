#ifndef __FSHT_H_
#define __FSHT_H_

#include <complex.h>
#include <fftw3.h>

typedef struct {
	int ntheta, nphi, deg;
	double *theta, *weights;
	fftw_plan fplan, bplan;
} shdata;

int fshtinit (shdata *, int, int);
int fshtfree (shdata *);

int ffsht (complex double *, shdata *);
int ifsht (complex double *, shdata *);

int shscale (complex double *, shdata *, int);

#endif /* __FSHT_H_ */

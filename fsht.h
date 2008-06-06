#ifndef __FSHT_H_
#define __FSHT_H_

#include <complex.h>
#include <fftw3.h>

typedef struct {
	int ntheta, nphi, deg;
	double *theta, *weights, *lgvals;
	complex double *fftbuf;
	fftw_plan fplan, bplan;
} shdata;

int fshtinit (shdata *, int);
int fshtfree (shdata *);

int ffsht (complex double *, shdata *);
int ifsht (complex double *, shdata *);

#endif /* __FSHT_H_ */

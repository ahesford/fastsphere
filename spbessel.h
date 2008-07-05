#ifndef __SPBESSEL_H_
#define __SPBESSEL_H_

#include <complex.h>

/* The Amos library routines upon which the spherical routines are based. */
void zbesh_ (double *, double *, double *, int *, int *, int *, double *,
		double *, int *, int *);
void zbesj_ (double *, double *, double *, int *, int *, double *,
		double *, int *, int *);

/* Compute the spherical Bessel functions for integer orders. */
int spbesh (complex double *, complex double, int);
int spbeshl (complex long double *, complex long double, int);
int spbesj (complex double *, complex double, int);
int spbesjl (complex long double *, complex long double, int);

/* Compute the derivatives of spherical Bessel functions for integer orders. */
int dspbesh (complex double *, complex double *, complex double, int);
int dspbesj (complex double *, complex double *, complex double, int);
int dspbeshl (complex long double *, complex long double *, complex long double, int);
int dspbesjl (complex long double *, complex long double *, complex long double, int);

#endif /* __SPBESSEL_H_ */

#ifndef __TRANSLATOR_H_
#define __TRANSLATOR_H_

#include <complex.h>

int exband (complex double, double);
int legpoly (int, double, double *);
complex double transang (int, complex double *, double *, double *, double, double);
int translator (complex double *, int, int, int, double *, complex double, double *);

#endif /* __TRANSLATOR_H_ */

#ifndef __SPREFLECT_H_
#define __SPREFLECT_H_

#include <complex.h>

int spbldrc (complex double *, complex double, complex double,
		double, double, double, int);
int esbldrc (complex double *, complex double *, complex double,
		complex double, double, double, double, int);

int spreflect (complex double *, complex double *, complex double *, int, int, int);

#endif /* __SPREFLECT_H_ */

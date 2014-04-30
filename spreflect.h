#ifndef __SPREFLECT_H_
#define __SPREFLECT_H_

#include <complex.h>
#include "fastsphere.h"

int spbldrc (sptype *, complex double, double);
int esbldrc (sptype *, complex double, double);

int spreflect (complex double *, complex double *, complex double *, int, int, int, int);

#endif /* __SPREFLECT_H_ */

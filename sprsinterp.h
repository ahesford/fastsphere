#ifndef __SPRSINTERP_H_
#define __SPRSINTERP_H_

#include <complex.h>

#include "fsht.h"

typedef struct {
	int nnz, *idx;
	double *val;
} sprow;

int smvmpy (complex double *, int, sprow *, complex double *);
int smvtmpy (complex double *, int, sprow *, complex double *);

int bldsrow (sprow *, int, int, int, double *, double *, int *, int *);

int intpmat (sprow *, shdata *, shdata *, int);

#endif /* __SPRSINTERP_H_ */

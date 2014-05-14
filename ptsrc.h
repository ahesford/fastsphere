#ifndef __PTSRC_H_
#define __PTSRC_H_

#include <complex.h>

#include "fsht.h"

double directivity (double *, double *, double);
int ptsrcexp (complex double *, complex double, shdata *, 
		complex double, double *, double *, double);

#endif /* __PTSRC_H_ */

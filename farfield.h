#ifndef __FARFIELD_H_
#define __FARFIELD_H_

#include <complex.h>

#include "sprsinterp.h"
#include "fastsphere.h"
#include "fsht.h"

int rootorder (spscat *, int, complex double);
int neartofar (complex double *, complex double *, spscat *, int,
		complex double, shdata *, shdata *, sprow *);
int fartonear (complex double *, complex double *, spscat *, int,
		complex double, shdata *, shdata *, sprow *);

#endif /* __FARFIELD_H_ */

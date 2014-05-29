#ifndef __FARFIELD_H_
#define __FARFIELD_H_

#include <complex.h>

#include "fastsphere.h"
#include "fsht.h"

int rootorder (spscat *, int, complex double, int);
int neartofar (complex double *, complex double *, spscat *, int,
		complex double, shdata *, shdata *);
int fartonear (complex double *, complex double *, spscat *, int,
		complex double, shdata *, shdata *);

#endif /* __FARFIELD_H_ */

#ifndef __FARFIELD_H_
#define __FARFIELD_H_

#include <complex.h>

#include "fastsphere.h"
#include "fsht.h"

int rootorder (spscat *, int, bgtype *);
void expnsh (complex double *, complex double *, shdata *, shdata *);
int farfield (complex double *, complex double *, spscat *, int,
		bgtype *, shdata *, shdata *);

#endif /* __FARFIELD_H_ */

#ifndef __INIT_H_
#define __INIT_H_

#include <stdio.h>

#include "fsht.h"
#include "fastsphere.h"
#include "translator.h"

int sphinit (sptype *, int, complex double, double, shdata *, int, int);
int esbdinit (sptype *, complex double, double, shdata *, int, int);
void clrspheres (sptype *, int);

trdesc* sphbldfmm (spscat *, int, complex double, shdata *);
void sphclrfmm (trdesc *, int);

#endif /* __INIT_H_ */

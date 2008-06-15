#ifndef __CONFIG_H_
#define __CONFIG_H_

#include <stdio.h>

#include "fsht.h"
#include "fastsphere.h"

int sphinit (sptype *, int, bgtype *, shdata *);
void clrspheres (sptype *, int);

int sphbldfmm (complex double ***, spscat *, int, bgtype *, shdata *);

#endif /* __CONFIG_H_ */

#ifndef __INIT_H_
#define __INIT_H_

#include <stdio.h>

#include "fsht.h"
#include "fastsphere.h"
#include "translator.h"

int sphinit (sptype *, int, bgtype *, shdata *);
void clrspheres (sptype *, int);

trdesc* sphbldfmm (spscat *, int, bgtype *, shdata *);
void sphclrfmm (trdesc *, int);

#endif /* __INIT_H_ */

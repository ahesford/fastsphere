#ifndef __CONFIG_H_
#define __CONFIG_H_

#include <stdio.h>

#include "fastsphere.h"
#include "scatmat.h"

/* Minimum and default (when unspecified) number of digits of 
 * accuracy for diagonal translation and harmonic expansions. */
#define MIN_ACC_DIG 3
#define DEF_ACC_DIG 6

int nextline (FILE *, char *, int);

int readcfg (FILE *, int *, int *, sptype **, sptype *,
		spscat **, bgtype *, exctparm *, itconf *, int *);

int setspht (sptype *, int, int, int);
int clrspht (sptype *, int);

#endif /* __CONFIG_H_ */

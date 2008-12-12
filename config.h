#ifndef __CONFIG_H_
#define __CONFIG_H_

#include <stdio.h>

#include "fastsphere.h"
#include "scatmat.h"

int nextline (FILE *, char *, int);

int readcfg (FILE *, int *, int *, sptype **, sptype *,
		spscat **, bgtype *, exctparm *, itconf *);

int setspht (sptype *, int, int, int);
int clrspht (sptype *, int);

#endif /* __CONFIG_H_ */

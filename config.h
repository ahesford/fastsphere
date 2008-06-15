#ifndef __CONFIG_H_
#define __CONFIG_H_

#include <stdio.h>

#include "fastsphere.h"

int nextline (FILE *, char *, int);

int readcfg (FILE *, int *, int *, sptype **, spscat **, bgtype *, exctparm *);

int setspht (sptype *, int, int, int);
int clrspht (sptype *, int);

#endif /* __CONFIG_H_ */

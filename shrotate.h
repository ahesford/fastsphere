#ifndef __SHROTATE_H_
#define __SHROTATE_H_

#include <complex.h>

int getangles (double *, double *, double *, double *);
int zrotate (complex double *, int, int, double);
int xrotate (complex double *, int, int);
int shrotate (complex double *, int, int, double, double, double);

#endif /* __SHROTATE_H_ */

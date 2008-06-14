#ifndef __UTIL_H_
#define __UTIL_H_

#include <complex.h>

/* Some helpful macros. */
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

int exband (complex double, double);
int legpoly (int, double, double *);

complex double buildkvec (double, double, double);
#endif /* __UTIL_H_ */

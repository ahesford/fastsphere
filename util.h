#ifndef __UTIL_H_
#define __UTIL_H_

#include <complex.h>

/* Some helpful macros. */
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define SGN(a) ((a) < 0 ? -1 : 1)

#define DVDOT(a,b) ((a)[0] * (b)[0] + (a)[1] * (b)[1] + (a)[2] * (b)[2])

int exband (complex double, double);
int legpoly (int, double, double *);

complex double buildkvec (double, double);
#endif /* __UTIL_H_ */

#ifndef __UTIL_H_
#define __UTIL_H_

#include <stdio.h>
#include <complex.h>

/* Some helpful macros. */
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define SGN(a) ((a) < 0 ? -1 : 1)
#define ABS(a) ((a) < 0 ? -(a) : (a))

#define DVDOT(a,b) ((a)[0] * (b)[0] + (a)[1] * (b)[1] + (a)[2] * (b)[2])
#define ELT(n,m,lda) ((m) * (lda) + (n))
#define IDX(n,m,lda) ((n) < 0 ? ELT(n,m+1,lda) : ELT(n,m,lda))

#ifdef _FBSD
complex double cexp (complex double);
complex double csin (complex double);
complex double ccos (complex double);
#endif /* _FBSD */

double rmserror (complex double *, complex double *, int);
FILE *critopen (char *, char *);

int exband (complex double, double);
int legpoly (int, double, double *);

complex double buildkvec (double, double);

int copysh (int, complex double *, int, complex double *, int);
int shradial (int, complex double *, int, complex double, double);
int shincident (int, complex double *, int, complex double, double, double);
#endif /* __UTIL_H_ */

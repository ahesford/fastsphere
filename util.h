#ifndef __UTIL_H_
#define __UTIL_H_

#include <stdio.h>
#include <complex.h>

/* Some helpful macros. */
#ifndef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#endif
#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#endif
#define SGN(a) ((a) < 0 ? -1 : 1)
#define ABS(a) ((a) < 0 ? -(a) : (a))

#define DVDOT(a,b) ((a)[0] * (b)[0] + (a)[1] * (b)[1] + (a)[2] * (b)[2])
#define ELT(n,m,lda) ((m) * (lda) + (n))
#define IDX(n,m,lda) ((n) < 0 ? ELT(n,m+1,lda) : ELT(n,m,lda))

#define DEG2RAD(a) ((a) * M_PI / 180.0)
#define MAG(a,b,c) sqrt(((a) * (a)) + ((b) * (b)) + ((c) * (c)))

#define IMGS_TOL 1e-3
#define IMGS_ITS 2

#ifdef _FBSD
complex double cexp (complex double);
complex double csin (complex double);
complex double ccos (complex double);
#endif /* _FBSD */

int cmgs (complex double *, complex double *, complex double *, int, int);
complex double pardot (complex double *, complex double *, int);

double rmserror (complex double *, complex double *, int);
FILE *critopen (char *, char *);

int exband (complex double, double);
int legpoly (int, double, double *);

complex double wavenum (double, double);
complex double invwavenum (double, double);

int copysh (int, complex double *, int, complex double *, int);
int shradial (int, complex double *, int, complex double, double);
int shincident (int, complex double *, int, complex double, double, double);

int gauleg (int, double *, double *);
#endif /* __UTIL_H_ */

#ifndef __TRANSLATOR_H_
#define __TRANSLATOR_H_

#include <complex.h>

typedef enum { TRNONE = 0, TRPLANE, TRDENSE } trtype;

/* Describes the type of translation being performed. */
typedef struct {
	trtype type;
	complex double kr, *trdata;
	double sdir[3];
	int trunc;
} trdesc;

complex double transang (int, complex double *, double *, double *, double, double);
int translator (trdesc *, int, int, double *);

#endif /* __TRANSLATOR_H_ */

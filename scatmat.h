#ifndef __SCATMAT_H_
#define __SCATMAT_H_

#include <complex.h>

#include "fastsphere.h"
#include "translator.h"
#include "fsht.h"

typedef struct {
	int iter, restart;
	double eps;
} itconf;

enum GMRES_INFO { GEX = 0, GMV, GLP, GRP, GDP };

int buildrhs (complex double *, spscat *, int, shdata *);
int scatmat (complex double *, complex double *, spscat *, int, trdesc *, shdata *);
int itsolve (complex double *, complex double *,
		spscat *, int, trdesc *, shdata *, itconf *);

#endif /* __SCATMAT_H_ */

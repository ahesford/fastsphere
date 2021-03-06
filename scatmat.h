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

int sprflpw (complex double *, spscat *, int, shdata *);
int sptrans (complex double *, complex double *, int, trdesc *, shdata *);
int scatmat (complex double *, complex double *, spscat *, int,
		trdesc *, shdata *);
int bicgstab (complex double *, complex double *, int, spscat *, int,
		trdesc *, shdata *, itconf *);
int gmres (complex double *, complex double *, int, spscat *, int,
		trdesc *, shdata *, itconf *);

#endif /* __SCATMAT_H_ */

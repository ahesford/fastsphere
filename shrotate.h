#ifndef __SHROTATE_H_
#define __SHROTATE_H_

#include <complex.h>

#include "translator.h"

int getangles (double *, double *, double *);
int shrotate (complex double *, int, int, trdesc *);

#endif /* __SHROTATE_H_ */

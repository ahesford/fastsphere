#ifndef __SHTRANSLATE_H_
#define __SHTRANSLATE_H_

#include <complex.h>

#define ANM(n,m) sqrt((double)((n)+1+(m))*((n)+1-(m))/(double)((2*(n)+1)*(2*(n)+3)))
#define BNM(n,m) SGN(m)*sqrt((double)((n)-(m)-1)*((n)-(m))/(double)((2*(n)-1)*(2*(n)+1)))

int shtransrfl (complex double *, int, int);
int shtranslate (complex double *, int, int, complex double);

#endif /* __SHTRANSLATE_H_ */

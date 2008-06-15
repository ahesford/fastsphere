#ifndef __FASTSPHERE_H_
#define __FASTSPHERE_H_

#include <complex.h>

#include "fsht.h"

typedef struct {
	/* The complex wave number is k = kr + i * ki. The parts are defined:
	 * 	kr = 2 * pi / cr, where cr is the relative sound speed
	 * 	ki = -alpha * log (10) / 20, where alpha is attenuation in dB
	 * 	per wavelength.
	 * Thus exp(i * k * x) = exp(i * kr * x) / 10^(alpha * x / 20). */
	complex double k;

	/* Reflection coefficients for spherical harmonics. */
	complex double *reflect;

	double r; 	/* Radius in wavelengths. */
	double c;	/* Sound speed relative to background. */
	double alpha;	/* Attenuation (dB per wavelength). */
} sptype;

typedef struct {
	sptype *spdesc;	/* Pointer to sphere type structure. */
	double cen[3];	/* Location of the center of the sphere (wavelengths). */
} spscat;

typedef struct {
	/* The complex wave number. See above for description. */
	complex double k;

	double cabs;	/* Absolute sound speed (m/s). */
	double alpha;	/* Attenuation (dB per wavelength). */
} bgtype;

typedef struct {
	double f;		/* Excitation frequency. */
	double theta, phi;	/* Excitation direction. */
} exctparm;

#endif /* __FASTSPHERE_H_ */

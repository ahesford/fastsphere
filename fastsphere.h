#ifndef __FASTSPHERE_H_
#define __FASTSPHERE_H_

#include <complex.h>

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
	double rho;	/* Density relative to background. */
	int deg;	/* Maximum spherical harmonic degree. */
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
	double rho0;	/* Grams per cubic centimeter. */
} bgtype;

typedef struct {
	double f;	/* Excitation frequency. */
	int npw;	/* Number of plane waves. */
	complex double *mag;	/* Complex magnitude. */
	double *theta, *phi;	/* Angular positions. */
} exctparm;

#endif /* __FASTSPHERE_H_ */

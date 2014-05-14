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
	complex double *reflect, *transmit;

	double r; 	/* Radius in wavelengths. */
	double c;	/* Sound speed relative to background. */
	double csh;	/* Shear wave speed relative to background. */
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
	int npw, nps;	/* Number of plane waves and point sources. */
	complex double *pwmag;	/* Complex magnitudes, plane waves. */
	double *theta, *phi;	/* Angular positions, plane waves. */
	complex double *psmag;  /* Complex magnitudes, point sources. */
	double *psloc, *psax;   /* Coordinates and directivity axes, point sources. */
	double *alpha;		/* Directivity widths, point sources. */
} exctparm;

#endif /* __FASTSPHERE_H_ */

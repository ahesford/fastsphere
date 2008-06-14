#ifndef __FASTSPHERE_H_
#define __FASTSPHERE_H_

typedef struct {
	/* The complex wave number is k = kr + i * ki. The parts are defined:
	 * 	kr = 2 * pi * f / c, where c is the real sound speed
	 * 	ki = -alpha * 1e-4 * log (10) / 20, where alpha is attenuation
	 * 	in dB per cm * MHz.
	 * Thus exp(i * k * x) = exp(i * kr * x) / 10^(1e-4 * alpha * x / 20).
	 * Note the 1e-4 factor converts dB per cm * MHz to dB per m * Hz. */
	complex double k;

	/* Reflection coefficients for spherical harmonics. */
	complex double *reflect;

	double r; 	/* Radius. */
	double c;	/* Sound speed. */
	double alpha;	/* Attenuation (dB per cm * MHz). */
	double gamma;	/* Outside:inside contrast ratio, cb / c. */

	shdata *shinfo;	/* Information for spherical harmonic transforms. */
} sptype;

typedef struct {
	sptype *spdesc;	/* Pointer to sphere type structure. */
	double cen[3];	/* Location of the center of the sphere. */
} spscat;

#endif /* __FASTSPHERE_H_ */

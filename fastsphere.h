#ifndef __FASTSPHERE_H_
#define __FASTSPHERE_H_

/* 
 * The complex wave number is k = kr + i * ki. The parts are defined:
 * kr = 2 * pi * f / c, where c is the real sound speed
 * ki = -alpha * 1e-4 * log (10) / 20, where alpha is attenuation in dB per cm * MHz.
 *
 * Thus exp(i * k * x) = exp(i * kr * x) / 10^(1e-4 * alpha * x / 20).
 * Note the 1e-4 factor converts dB per cm * MHz to dB per m * Hz.
 */

typedef struct {
	/* Wave number and reflection coefficients for angular samples. */
	complex double k, *reflect;
	/* Radius, sound speed, attenuation (dB per cm * MHz). */
	double r, c, alpha;
	/* Information for spherical harmonic transforms. */
	shdata *shtinfo;
} sptype;

typedef struct {
	/* Pointer to sphere type structure. */
	sptype *spdesc;
	/* Location of the center of the sphere. */
	double center[3];
} spscat;

#endif /* __FASTSPHERE_H_ */

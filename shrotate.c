#include <stdlib.h>

#include <math.h>
#include <float.h>

/* Find the polar and rotation angles of the new z-axis for translation. */
int getangles (double *theta, double *chi, double axis[3]) {
	double st;

	*theta = acos (axis[2]); /* The polar angle of the new z-axis. */

	st = sqrt (axis[0] * axis[0] + axis[1] * axis[1]);

	/* If the z-axis hasn't changed, nothing changes. */
	if (*theta < DBL_MIN) {
		*theta = *phi = *chi = 0;
		return 0;
	}

	/* The y-coordinate of the new z-axis vector is given by
	 * y = sin(theta) sin(chi). Hence, chi is easily computed. Note
	 * that the x-coordiante is x = sin(theta) cos(chi), but this
	 * doesn't work well, since cosine is an even function. */
	*chi = asin (axis[1] / st);

	return 1;
}

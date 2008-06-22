#include <stdlib.h>

#include <math.h>
#include <float.h>

int getangles (double *theta, double *phi, double *chi, double axis[3]) {
	double st;

	*theta = acos (axis[2]); /* The polar rotation angle. */

	st = sqrt (axis[0] * axis[0] + axis[1] * axis[1]);

	/* If the z-axis hasn't changed, nothing changes. */
	if (*theta < DBL_MIN) {
		*theta = *phi = *chi = 0;
		return 0;
	}

	/* The old z-axis is always at phi = pi / 2. This is thre result of an
	 * arbitrary choice of the new x- and y-axis orientations. */
	*phi = M_PI_2;

	/* The rotation angle of the new z-axis must derive from the y-axis
	 * coordinate, since the x-axis coordinate uses the even cosine. */
	*chi = asin (axis[1] / st);

	return 1;
}

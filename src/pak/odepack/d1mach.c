#include <cdecs.h>

/*
*  double-precision machine constants
*  d1mach( 1) = b**(emin-1), the smallest positive magnitude. 
*  d1mach( 2) = b**emax*(1 - b**(-t)), the largest magnitude.
*  d1mach( 3) = b**(-t), the smallest relative spacing.
*  d1mach( 4) = b**(1-t), the largest relative spacing.
*  d1mach( 5) = log10(b) 
*/

double FORTRAN_NAME(d1mach)(int *pi)
{
    int           i;
    static int    first = 1;
    static double dmach[5];

    i = *pi;
    if (first)
    {
	first = 0;
	dmach[0] = DBL_MIN;
	dmach[1] = DBL_MAX;
	dmach[2] = 0.5*DBL_EPSILON;
	dmach[3] = DBL_EPSILON;
	dmach[4] = log10(2.0);
    }

    if (i < 1 || i > 5)
    {
        (void) printf("d1mach - i = %d out of bounds\n",i);
        return 0.0;
    }
    return dmach[i - 1];
}		/*end FORTRAN_NAME(d1mach) */

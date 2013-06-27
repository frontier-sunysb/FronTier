#include <stdio.h>
#include <math.h>

main()
{
	double x,y;
	int i;
	for (i = 0; i < 200; ++i)
	{
	    x = i*0.1 - 10.0;
	    y = erf(x);
	    printf("%f %f\n",x,y);
	}
}

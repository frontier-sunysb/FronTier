
#include "MonteCarlo.h"

extern  double gauss_newton(
	POINTER params,
        unsigned short int xsubi[3])
{
        double    t, z;
        int      num_iter;
        static  double  alpha, sqrt2, erf1, logsqrtpi;
        static  int first = YES;
	GAUSS_PARAMS *gauss_params = (GAUSS_PARAMS*)params;
	double mu = gauss_params->mu;
	double sigma = gauss_params->sigma;

        if (first == YES)
        {
            first = NO;
            alpha = 2.0/sqrt(PI);
            sqrt2 = sqrt(2.0);
            logsqrtpi = 0.5*log(PI);
            erf1 = erf(1.0);
        }
        t = 2.0*erand48(xsubi) - 1.0;
        if (t > erf1)
            z =  sqrt(-logsqrtpi - log(1.0-t));
        else if (t < -erf1)
            z = -sqrt(-logsqrtpi - log(1.0+t));
        else
            z = t/alpha;

        for (num_iter=0; (num_iter<10) && (fabs(t-erf(z))>EPSILON); ++num_iter)
        {
            z += (t - erf(z))/(alpha*exp(-z*z));
        }
        return sqrt2*sigma*z + mu;
}               /*end gauss_newton*/

extern  double gauss_box_muller(
	POINTER params,
        unsigned short int xsubi[3])
{
	double x1,x2,y;
	GAUSS_PARAMS *gauss_params = (GAUSS_PARAMS*)params;
	double mu = gauss_params->mu;
	double sigma = gauss_params->sigma;

        x1 = erand48(xsubi);
        x2 = erand48(xsubi);

	y = sqrt(-2.0*log(x1))*cos(2.0*PI*x2);
	y = sigma*y + mu;
	return y;
}	/* end gauss_box_muller */

extern  double gauss_center_limit(
	POINTER params,
        unsigned short int xsubi[3])
{
	GAUSS_PARAMS *gauss_params = (GAUSS_PARAMS*)params;
	double mu = gauss_params->mu;
	double sigma = gauss_params->sigma;
	double x;
	int i;
	for (i = 0, x = 0.0; i < 12; ++i)
	    x += erand48(xsubi);
	x -= 6.0;
	x = sigma*x + mu;
	return x;
}	/* end gauss_center_limit */

extern  double dist_exponential(
	POINTER params,
        unsigned short int xsubi[3])
{
	double x,y;
	int i;
	EXP_PARAMS *exp_params = (EXP_PARAMS*)params;
	double lambda = exp_params->lambda;
	x = erand48(xsubi);
	y = -1.0/lambda*log(x);
	return y;
}	/* end gauss_center_limit */

extern  double dist_power(
	POINTER params,
        unsigned short int xsubi[3])
{
	POWER_PARAMS *power_params = (POWER_PARAMS*)params;
	double x,y;
	int i,power = power_params->power;
	y = erand48(xsubi);
	for (i = 1; i < power; ++i)
	{
	    x = erand48(xsubi);
	    if (y < x) y = x;
	}
	return y;
}	/* end gauss_center_limit */

extern  double dist_middle(
	POINTER params,
        unsigned short int xsubi[3])
{
	double x[3],y,tmp;
	int i,j;
	for (i = 0; i < 3; ++i)
	{
	    y = erand48(xsubi);
	    for (j = 0; j < i; ++j)
	    {
		if (y < x[j])
		{
		    tmp = x[j];
		    x[j] = y;
		    y = tmp;
		}	
	    }
	    x[i] = y;
	}
	return x[1];
}	/* end gauss_center_limit */

extern  double dist_cauchy(
	POINTER params,
        unsigned short int xsubi[3])
{
	double x,y;
	GAUSS_PARAMS *gauss_params = (GAUSS_PARAMS*)params;
	double mu = gauss_params->mu;
	double sigma = gauss_params->sigma;
	x = 2.0*erand48(xsubi) - 1.0;
	y = tan(0.5*PI*x);
	y = sigma*y + mu;
	return y;
}	/* end gauss_center_limit */


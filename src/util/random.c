/***************************************************************
FronTier is a set of libraries that implements differnt types of 
Front Traking algorithms. Front Tracking is a numerical method for 
the solution of partial differential equations whose solutions have 
discontinuities.  

Copyright (C) 1999 by The University at Stony Brook. 

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
****************************************************************/

#include <cdecs.h>

EXPORT  double gauss_newton(
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

        for (num_iter=0; (num_iter<10) && (fabs(t-erf(z))>10e-14); ++num_iter)
        {
            z += (t - erf(z))/(alpha*exp(-z*z));
        }
        return sqrt2*sigma*z + mu;
}               /*end gauss_newton*/

EXPORT  double gauss_box_muller(
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

EXPORT  double gauss_center_limit(
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

EXPORT  double dist_exponential(
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

EXPORT  double dist_power(
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

EXPORT  double dist_middle(
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

EXPORT  double dist_cauchy(
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

EXPORT  double dist_uniform(
	POINTER params,
        unsigned short int xsubi[3])
{
	double x,y;
	UNIFORM_PARAMS *uniform_params = (UNIFORM_PARAMS*)params;
	double a = uniform_params->a;
	double b = uniform_params->b;
	double k = b - a;

	x = erand48(xsubi);
	y = k*x + a;
	return y;
}	/* end dist_uniform */

EXPORT  double dist_stable(
	POINTER params,
        unsigned short int xsubi[3])
{
	double u,e,x,y;
	static double piover2,b,c,s,temptan,r_alpha;
	static UNIFORM_PARAMS *uniform_params;
	static EXP_PARAMS *exp_params;
	STABLE_PARAMS *stable_params = (STABLE_PARAMS*)params;
	double alpha = stable_params->alpha;
	double beta = stable_params->beta;
	double sigma = stable_params->sigma;
	double mu = stable_params->mu;

	if (uniform_params == NULL)
	{
	    FT_ScalarMemoryAlloc((POINTER*)&uniform_params,
			sizeof(UNIFORM_PARAMS));
	    FT_ScalarMemoryAlloc((POINTER*)&exp_params,
			sizeof(EXP_PARAMS));
	    uniform_params->a = -PI/2.0;
	    uniform_params->b =  PI/2.0;
	    exp_params->lambda = 1.0;
	    piover2 = PI/2.0;
	    temptan = beta*tan(piover2*alpha);
	    s = pow(1+temptan*temptan,1.0/2.0/alpha);
	    b = atan(temptan)/alpha;
	    c = beta*sigma*log(sigma)/piover2+mu;
	    r_alpha = 1.0/alpha;
	}

	/* u is a variable with uniform distribution from -PI/2 to PI/2 */
	u = dist_uniform((POINTER)uniform_params,xsubi);
	/* e is a variable with exponential distribution: lambda = 1 */
	e = dist_exponential((POINTER)exp_params,xsubi);
	if (alpha == 1.0)
	{
	    x = ((piover2+beta*u)*tan(u)-beta*log((piover2*e*cos(u))/
				(piover2+beta*u)))/piover2;
	    y = sigma*x + c;
	}
	else
	{
	    x = s*sin(alpha*(u+b))/pow(cos(u),r_alpha)*
				pow(cos(u-alpha*(u+b))/e,r_alpha-1.0);
	    y = sigma*x+mu;
	}
	return y;
}	/* end dist_stable */

EXPORT  double dist_gig(
        POINTER params,
        unsigned short int xsubi[3])
{
	GIG_PARAMS *gig_params = (GIG_PARAMS*)params;
	static UNIFORM_PARAMS *U,*V;
	double lambda = gig_params->lambda;
	double psi = gig_params->psi;
	double chi = gig_params->chi;
	double beta = sqrt(psi*chi);
	double m,a,b,c,p,q,phi,xminus,xplus,vplus,uminus,uplus,u,v,x;
	double x0,x1,k1,k2,k3,A1,A2,A3,A,h;

	if (lambda < 0.0)
	{
	    static GIG_PARAMS *new_gig_params;
	    if (new_gig_params == NULL)
		FT_ScalarMemoryAlloc((POINTER*)&new_gig_params,
			sizeof(GIG_PARAMS));
	    new_gig_params->lambda = -lambda;
	    new_gig_params->psi = chi;
	    new_gig_params->chi = psi;
	    double gig = dist_gig((POINTER)new_gig_params,xsubi);
	    return 1.0/gig;
	}

	if (U == NULL)
        {
            FT_ScalarMemoryAlloc((POINTER*)&U,sizeof(UNIFORM_PARAMS));
            FT_ScalarMemoryAlloc((POINTER*)&V,sizeof(UNIFORM_PARAMS));
        }

	if (lambda > 1.0 || beta > 1.0)
	{
	    m = (sqrt((lambda-1.0)*(lambda-1.0)+beta*beta)+lambda-1.0)/beta;
	    a = -2.0*(lambda+1.0)/beta-m;
	    b = 2.0*(lambda-1.0)*m/beta-1.0;
	    c = m;
	    p = b-a*a/3.0;
	    q = 2.0*a*a*a/27.0-a*b/3.0+c;
	    phi = acos(0.5*q*sqrt(-27.0/p)/p);
	    xminus = sqrt(-4.0*p/3.0)*cos(phi/3.0+4.0*PI/3.0)-a/3.0;
	    xplus = sqrt(-4.0*p/3.0)*cos(phi/3.0)-a/3.0;
	    vplus = pow(m,(lambda-1.0)/2.0)*exp(-beta/4.0*(m+1.0/m));
	    uminus = (xminus-m)*pow(xminus,(lambda-1.0)/2.0)
		*exp(-beta/4.0*(xminus+1.0/xminus));
	    uplus = (xplus-m)*pow(xplus,(lambda-1.0)/2.0)
                *exp(-beta/4.0*(xplus+1.0/xplus));
	    U->a = uminus;	U->b = uplus;
	    V->a = 0.0;		V->b = vplus;
	    do
	    {
	    	u = dist_uniform((POINTER)U,xsubi);
	    	v = dist_uniform((POINTER)V,xsubi);
		x = u/v+m;
	    } while (x<0 || v*v>pow(x,lambda-1.0)*exp(-beta/2.0*(x+1.0/x)));
	}
	else if (beta >= min(1.0/2.0,2.0*sqrt(1.0-lambda)/3.0))
	{
	    m = beta/(sqrt((1.0-lambda)*(1.0-lambda)+beta*beta)+1.0-lambda);
	    xplus = (sqrt((lambda+1.0)*(lambda+1.0)+beta*beta)+lambda+1.0)/beta;
	    vplus = pow(m,(lambda-1.0)/2.0)*exp(-beta/4.0*(m+1.0/m));
	    uplus = xplus*pow(xplus,(lambda-1.0)/2.0)
		*exp(-beta/4.0*(xplus+1.0/xplus));
	    U->a = 0.0;      U->b = uplus;
            V->a = 0.0;         V->b = vplus;
	    do
	    {
		u = dist_uniform((POINTER)U,xsubi);
                v = dist_uniform((POINTER)V,xsubi);
                x = u/v;
	    } while (v*v>pow(x,lambda-1)*exp(-beta/2.0*(x+1.0/x)));
	}
	else
	{
	    m = beta/(sqrt((lambda-1.0)*(lambda-1.0)+beta*beta)-lambda+1.0);
	    x0 = beta/(1.0-lambda);
	    x1 = max(x0,2.0/beta);
	    k1 = pow(m,lambda-1.0)*exp(-beta/2.0*(m+1.0/m));
	    A1 = k1*x0;
	    if (x0 < 2.0/beta)
	    {
		k2 = exp(-beta);
		A2 = (lambda == 0.0) ? k2*log(2.0/beta*beta):
			k2*(pow(2.0/beta,lambda)-pow(x0,lambda))/lambda;
	    }
	    else	
	    {
		k2 = 0.0;
		A2 = 0.0;
	    }
	    k3 = pow(x1,lambda-1.0);
	    A3 = 2.0*k3*exp(-x1*beta/2.0)/beta;
	    A = A1 + A2 + A3;
	    U->a = 0.0;		U->b = 1.0;
            V->a = 0.0;		V->b = A;
	    do
	    {
		u = dist_uniform((POINTER)U,xsubi);
                v = dist_uniform((POINTER)V,xsubi);
		if (v <= A1)
		{
		    x = x0*v/A1;
		    h = k1;
		}
		else if (v <= A1+A2)
		{
		    v = v-A1;
		    x = (lambda == 0) ? beta*exp(v*exp(beta)):
			pow(pow(x0,lambda)+v*lambda/k2,1.0/lambda);
		    h = k2*pow(x,lambda-1.0);
		}
		else
		{
		    v = v-A1-A2;
		    x = -2.0/beta*log(exp(-x1*beta/2.0)-0.5*v*beta/k3);
		    h = k3*exp(-x*beta*0.5);
		}
	    } while (u*h>pow(x,lambda-1)*exp(-beta/2.0*(x+1.0/x)));
	}
	return sqrt(chi/psi)*x;
}

EXPORT  double dist_gh(
        POINTER params,
        unsigned short int xsubi[3])
{
	GH_PARAMS *gh_params = (GH_PARAMS*)params;
	double lambda = gh_params->lambda;
	double alpha = gh_params->alpha;
	double beta = gh_params->beta;
	double delta = gh_params->delta;
	double mu = gh_params->mu;

	static GIG_PARAMS *gig_params;
	static GAUSS_PARAMS *U;
	double gig, u;

	if (U==NULL)
	{
	    FT_ScalarMemoryAlloc((POINTER*)&gig_params,sizeof(GIG_PARAMS));
            FT_ScalarMemoryAlloc((POINTER*)&U,sizeof(GAUSS_PARAMS));
	}
	gig_params->lambda = lambda;
	gig_params->psi = alpha*alpha-beta*beta;
	gig_params->chi = delta*delta;
	U->sigma = 1.0;
	U->mu = 0.0;
	gig = dist_gig((POINTER)gig_params,xsubi);
	u = gauss_box_muller((POINTER)U,xsubi);

	return mu + beta*gig + sqrt(gig)*u;
}

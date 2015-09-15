/***************************************************************
FronTier is a set of libraries that implements different types of 
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

#include "cFluid.h"

static int BisectionFindRoot(double(*f)(double,void*),void*,double*,
			double,double,double,double);
static int SecantFindRoot(double(*f)(double,void*),void*,double*,
			double,double,double,double);
static int NewtonFindRoot(double(*f)(double,void*),double(*f_prime)(double,void*),void*,double*,
			double,double,double);
static double UlMinusUr(double,void*);
static double UlMinusUr_Prime(double,void*);
static double u_left_centered_wave(RIEM_STATE,double);
static double u_right_centered_wave(RIEM_STATE,double);
static void OutputSolution(RIEMANN_SOLN*,RIEM_STATE*,double,int);

extern boolean RiemannSolnAtXi(		/* xi: x/t */
	RIEMANN_SOLN *riem_soln,
	RIEM_STATE *soln_at_xi,
	double xi)
{
	double sl,sr,pl,rhol,pr,rhor,ul,ur,cl,cr,c_star_l,c_star_r,A;
	double p_star,u_star,d_star_l,d_star_r;
	double gamma_l = riem_soln->left_state.gamma;
	double gamma_r = riem_soln->right_state.gamma;
	
        rhol = riem_soln->left_state.d;
	pl = riem_soln->left_state.p;
	ul = riem_soln->left_state.u;
	cl = sqrt(pl*gamma_l/rhol);

	rhor = riem_soln->right_state.d;
	pr = riem_soln->right_state.p;
	ur = riem_soln->right_state.u;
	cr = sqrt(pr*gamma_r/rhor);
	
	u_star = riem_soln->left_center_state.u;
	p_star = riem_soln->left_center_state.p;
	d_star_l = riem_soln->left_center_state.d;
	d_star_r = riem_soln->right_center_state.d;
	c_star_l = sqrt(p_star*gamma_l/d_star_l);
	c_star_r = sqrt(p_star*gamma_r/d_star_r);

        if(riem_soln->left_wave.wave_type == LF_SHOCK)
	{   
	    sl = (rhol*ul - d_star_l*u_star)/(rhol - d_star_l);
	    if(xi < sl)
	        *soln_at_xi = riem_soln->left_state;
	    else if(xi > sl && u_star > xi)
	        *soln_at_xi = riem_soln->left_center_state;
	}	
	if(riem_soln->right_wave.wave_type == RF_SHOCK)
	{    
	    sr = (rhor*ur - d_star_r*u_star)/(rhor - d_star_r);
	    if(xi > sr)
	        *soln_at_xi = riem_soln->right_state;
	    else if(sr > xi && xi > u_star)
	        *soln_at_xi = riem_soln->right_center_state;
	}	
	if(riem_soln->left_wave.wave_type == GAMMA_PLUS)
	{
	    A = pl/pow(rhol,gamma_l); 
	    if(xi < ul - cl && xi > -1.0)
	        *soln_at_xi = riem_soln->left_state;
	    	
	    else if(xi > (u_star - c_star_l) && xi < u_star)
	        *soln_at_xi = riem_soln->left_center_state;
	    
	    else if(xi > (ul - cl) && xi < (u_star - c_star_l))
	    {
		RIEM_STATE inrfw;
		double c;
		inrfw.u = (2.0/(gamma_l - 1.0)*(cl + xi) + ul)/
				(2.0/(gamma_l - 1.0) + 1.0);
		c = inrfw.u - xi;
		inrfw.d = pow(pow(c,2)/(A*gamma_l),1.0/(gamma_l - 1.0));
		inrfw.p = A*pow(inrfw.d, gamma_l);
		*soln_at_xi = inrfw;
	    }
	}
	if(riem_soln->right_wave.wave_type == GAMMA_MINUS)
	{ 
	    A = pr/pow(rhor,gamma_r);
	    if(xi > ur + cr && xi < 1.0)
	        *soln_at_xi = riem_soln->right_state;
	    else if(xi < (u_star + c_star_r) && xi > u_star)
	        *soln_at_xi = riem_soln->right_center_state;
	    else if(xi > (u_star + c_star_r) && xi < (ur + cr))
	    {
	        RIEM_STATE inrfw;
		double c;
		inrfw.u = (2.0/(gamma_r - 1.0)*(xi - cr) + ur)/
				(2.0/(gamma_r - 1.0) + 1.0);
		c = xi - inrfw.u;
		inrfw.d = pow(pow(c,2)/(A*gamma_r),1.0/(gamma_r - 1.0));
		inrfw.p = A*pow(inrfw.d, gamma_r);
		*soln_at_xi = inrfw;
	    }
	}
	
	return YES;
}	/* end RiemannSolnAtXi */

extern boolean RiemannSolution(
	RIEMANN_INPUT input,
	RIEMANN_SOLN *riem_soln)
{
	void *params;
	double (*func)(double,void*);
	double (*func_prime)(double,void*);
	double epsilon,delta;
	double p_star,p0,p1;
	double gamma_l,gamma_r;

	epsilon = 0.0000001;
	delta = 0.0000001;
	gamma_l = input.left_state.gamma;
	gamma_r = input.right_state.gamma;
	if(input.left_state.p != input.right_state.p)
	{    
	    p0 = input.left_state.p;
	    p1 = input.right_state.p;
	}
	else 
	{
	   printf("\nAdjusting p0 and p1...");
	   p0 = 0.001;
	   p1 = 100;
	}
	params = (void*)&input;
	func = UlMinusUr;
	func_prime = UlMinusUr_Prime;
	//NewtonFindRoot(func,func_prime,params,&p_star,p0,delta,epsilon);
	//BisectionFindRoot(func,params,&p_star,p0,p1,delta,epsilon);
	SecantFindRoot(func,params,&p_star,p0,p1,delta,epsilon);
	//Set left and right states the same as inputted states.
	riem_soln->left_state = input.left_state;
	riem_soln->right_state = input.right_state;

	//Set the pressures (p_star) of the left and right center states. 
	riem_soln->left_center_state.p = p_star;
	riem_soln->right_center_state.p = p_star;
        
	//Use u_left_centered_wave(p_star) or u_right_centered_wave(p_star) 
	//to get u_star.
	riem_soln->left_center_state.u = 
			u_left_centered_wave(input.left_state,p_star);
	riem_soln->right_center_state.u = 
			u_right_centered_wave(input.right_state,p_star);

	//Use either the isentropic law or the Hugoniot equations to 
	//get rho star left and rho star right.
	double pl,pr,rhol,rhor;
	pl = input.left_state.p;
	rhol = input.left_state.d;
	pr = input.right_state.p;
	rhor = input.right_state.d;

	//Left center wave
	if(p_star < pl) //Rarefaction
	{   
	    riem_soln->left_center_state.d = rhol*pow(p_star/pl, 1.0/gamma_l); 
	    riem_soln->left_wave.wave_type = GAMMA_PLUS;
	}
	else //Shock wave
	{
	    riem_soln->left_center_state.d = (p_star/(gamma_l - 1.0) + 
				0.5*(p_star + pl))*rhol/(pl/(gamma_l - 1.0) + 
				0.5*(p_star + pl));
            riem_soln->left_wave.wave_type = LF_SHOCK;
	}	    
	
	//Right center wave
	if(p_star < pr) //Rarefaction
	{
	    riem_soln->right_center_state.d = rhor*pow(p_star/pr, 1.0/gamma_r); 
	    riem_soln->right_wave.wave_type = GAMMA_MINUS;
	}
	else //Shock wave
	{
            riem_soln->right_center_state.d = (p_star/(gamma_r - 1.0) + 
				0.5*(p_star + pr))*rhor/(pr/(gamma_r - 1.0) + 
				0.5*(p_star + pr));
	    riem_soln->right_wave.wave_type = RF_SHOCK;
	}

	if(fabs(riem_soln->right_center_state.d) > 0.001)
	{
	    riem_soln->contact.wave_type = CONTACT;
	    riem_soln->contact.speed_contact = riem_soln->right_center_state.u;
	}
	else
	    riem_soln->contact.wave_type = VACUUM;
	return YES;
}	/* end RiemannSolution */

/*	bisection method to solve equation f(x) = 0        */
static int BisectionFindRoot(        	/* return YES and NO */ 
	double (*f)(double,void*),	/* function f(x) */
        void *f_params,         	/* function parameters, unknown type */
        double *px,         	/* pointer to solution */
	double x0,		/* input of initial guess 0 */
        double x1,         	/* input of initial guess 1 */
	double delta,		/* tolerance for x */
        double epsilon) 		/* tolerance for f(x) */
{
        double c; /*The midpoint value*/
	int i;
	
	if(f(x0,f_params)*f(x1,f_params) > 0)
	{
	    printf("Invalid intial choices for x0 and x1. Aborting program.");
	    return NO;
	}
        
	while((fabs(x0 - x1) > delta) && 
	       fabs((f(x0,f_params) - f(x1,f_params)) > epsilon))
	{
	    //Compute the midpoint of x0 and x1.
            c = (x0 + x1)/2.0;

            //Update either x0 or x1.
	    if(f(x0,f_params)*f(c,f_params) > 0)
		x0 = c;
	    else
		x1 = c;
	}
        *px = (x0 + x1)/2.0;
	return YES;
}	/* end bisection_find_root */

/*	secant method to solve equation f(x) = 0	*/
static int SecantFindRoot(        	
	 double (*f)(double,void*),	/* function f(x) */
         void *f_params,         	/* function parameters, unknown type */
         double *px,         	/* pointer to solution */
       	 double x0,		/* input of initial guess 1 */
         double x1,         	/* input of initial guess 0 */
         double epsilon,		/* tolerance for x */
         double delta) 		/* tolerance for f(x) */
{
        double d,xn,xn_1,fxn,fxn_1;
	int count = 0;
	xn = x1;
	xn_1 = x0;
	fxn = f(xn,f_params);
	fxn_1 = f(xn_1,f_params);
        while((fabs(xn - xn_1) > delta) && fabs(fxn - fxn_1) > epsilon)
	{
	    count++;
	    d = (xn - xn_1)/(fxn - fxn_1)*fxn;
	    xn_1 = xn;
	    xn = fabs(xn - d);
	    fxn_1 = fxn;
	    fxn = f(xn,f_params);
	    if (count > 100)
	    {
		printf("Riemann solution using secant method not converge!\n");
		clean_up(ERROR);
	    }
	}
	*px = xn;
	return YES;
}	/* end SecantFindRoot */

/*	Newton method to solve equatin f(x) = 0 	*/
static int NewtonFindRoot(
	double (*f)(double,void*),      /* function f(x) */
	double (*f_prime)(double,void*),/* function f(x) prime */
        void *f_params,                 /* function parameters, unknown type */
        double *px,             /* pointer to solution */
        double x0,              /* input of initial guess */
        double epsilon,         /* tolerance for x */
        double delta)           /* tolerance for f(x) */
{
	double d,xn,xn_1,fxn,fxn_prime;
	xn = x0;
	xn_1 = x0 + 2.0*delta + 1;
	fxn = f(xn,f_params);
        fxn_prime = f_prime(xn,f_params);
	while(fabs(xn - xn_1) > delta && fabs(fxn) > epsilon)
	{
	    xn_1 = xn;
	    xn = fabs(xn - fxn/fxn_prime);
	    fxn = f(xn,f_params);
	    fxn_prime = f_prime(xn,f_params);
	}	
	*px = xn;
	return YES;
}	/* end NewtonFindRoot */

static double UlMinusUr(
	double p,
	void *params)
{
	double ul,ur;

	RIEMANN_INPUT *riem_input = (RIEMANN_INPUT*)params;
	ul = u_left_centered_wave(riem_input->left_state,p);
	ur = u_right_centered_wave(riem_input->right_state,p);
	return ul - ur;
}	/* end UlMinusUr */

static double UlMinusUr_Prime(
	double p,
	void *params)
{
	double rhol,pl,ul,c_l_star,rhol_star,M_l,UL;
	double rhor,pr,ur,c_r_star,rhor_star,M_r,UR;
	double mu_squared;
	double gamma_l,gamma_r;

	RIEMANN_INPUT *riem_input = (RIEMANN_INPUT*)params;
	gamma_l = riem_input->left_state.gamma;
	gamma_r = riem_input->right_state.gamma;
	rhol = riem_input->left_state.d;
	pl = riem_input->left_state.p;
	ul = riem_input->left_state.u;

	rhor = riem_input->right_state.d;
	pr = riem_input->right_state.p;
	ur = riem_input->right_state.u;
	
	if(p < pl)
	{
            c_l_star = sqrt(gamma_l*p/(rhol*pow(p/pl,1.0/gamma_l)));
	    rhol_star = gamma_l*p/(c_l_star*c_l_star);
	
	    UL = -1.0/(c_l_star*rhol_star);
	}
	else
	{
            mu_squared = (gamma_l - 1.0)/(gamma_l + 1.0);
            M_l = sqrt(rhol*(mu_squared/(1.0 - mu_squared)*pl + 
				1.0/(1.0 - mu_squared)*p));
	    
	    UL = (-1.0*M_l + rhol/(2.0*(1.0 - mu_squared)*M_l)*
				(p + pl))/(M_l*M_l);
	}

	if(p < pr)
	{
  	    c_r_star = sqrt(gamma_r*p/(rhol*pow(p/pr,1.0/gamma_r)));
	    rhor_star = gamma_r*p/(c_r_star*c_r_star);

	    UR = 1.0/(c_r_star*rhor_star);    
	}
	else
	{
	    mu_squared = (gamma_r - 1.0)/(gamma_r + 1.0);
	    M_r = sqrt(rhor*(mu_squared/(1.0 - mu_squared)*pr + 
				1.0/(1.0 - mu_squared)*p));

	    UR = (M_r - rhor/(2.0*(1.0 - mu_squared)*M_r)*(p + pr))/(M_r*M_r);
	}
	return UL - UR;
}

static double u_left_centered_wave(
	RIEM_STATE state,
	double p)
{
	double u_across_wave, mu_squared, c, c_star, M;
	double gamma = state.gamma;

	if (p < state.p) //For rarefaction 
	{
	    c = sqrt(gamma*state.p/state.d);
	    c_star = sqrt(gamma*p/(state.d*pow(p/state.p,1.0/gamma)));
	    u_across_wave = state.u + 2.0/(gamma - 1.0)*(c - c_star);
	}
	else //For shock 
	{
	    mu_squared = (gamma - 1.0)/(gamma + 1.0);
	    M = sqrt(state.d*(mu_squared/(1.0 - mu_squared)*state.p + 
				1.0/(1.0 - mu_squared)*p));
	    u_across_wave = state.u + (state.p - p)/M;
	}
	return u_across_wave;
}	/* end u_across_centered_wave */

static double u_right_centered_wave(
	RIEM_STATE state,
	double p)
{
	double u_across_wave, mu_squared, c, c_star, M;
	double gamma = state.gamma;
	if (p < state.p) //For rarefaction 
	{
	    c = sqrt(gamma*state.p/state.d);
	    c_star = sqrt(gamma*p/(state.d*pow(p/state.p,1.0/gamma)));
	    u_across_wave = state.u - 2.0/(gamma - 1.0)*(c - c_star); 
	}
	else //For shock 
	{
	    mu_squared = (gamma - 1.0)/(gamma + 1.0);
	    M = -1.0*sqrt(state.d*(mu_squared/(1.0 - mu_squared)*state.p + 
				1.0/(1.0 - mu_squared)*p));
	    u_across_wave = state.u + (state.p - p)/M;
	}
	return u_across_wave;
}	/* end u_across_centered_wave */

static void OutputSolution(RIEMANN_SOLN *sol,RIEM_STATE *sol_at_xi,double xi,int flag)
{
	printf("\nThe solution is...");
	printf("\nleft state:");
	printf("\ndensity: %f", sol->left_state.d);
	printf("\npressure: %f", sol->left_state.p);
	printf("\nvelocity: %f", sol->left_state.u);

	printf("\n\nleft center state:");
        printf("\ndensity: %f", sol->left_center_state.d);
        printf("\npressure: %f", sol->left_center_state.p);
        printf("\nvelocity: %f", sol->left_center_state.u);

        printf("\n\nright center state:");
        printf("\ndensity: %f", sol->right_center_state.d);
        printf("\npressure: %f", sol->right_center_state.p);
	printf("\nvelocity: %f", sol->right_center_state.u);

        printf("\n\nright state:");
	printf("\ndensity: %f", sol->right_state.d);
        printf("\npressure: %f", sol->right_state.p);
	printf("\nvelocity: %f", sol->right_state.u);

	printf("\n\nThe left center wave is a ");
	if(sol->left_wave.wave_type == 2)
	    printf("gamma plus rarefaction wave.");
	else
	    printf("left facing shock wave.");

	printf("\nThe right center wave is a ");
	if(sol->right_wave.wave_type == 3)
	    printf("gamma minus rarefaction wave.");
	else
	    printf("right facing shock wave.");
 
        printf("\nThe center wave is a ");
	if(sol->contact.wave_type == 6)
	    printf("contact discontinuity.");
	else
	    printf("vacuum.");
	
	if(flag == 1)
	{
            printf("\nThe solution at xi = %f is:", xi);
	    printf("\ndensity: %f", sol_at_xi->d);
	    printf("\npressure: %f", sol_at_xi->p);
	    printf("\nvelocity: %f", sol_at_xi->u);
	}
}	/* end OutputSolution */

#include "inverse.h"
#define eps_plus  80
#define eps_minus  2

static int tangential_direction(int i,double *T,double *N,int dim);
static int exact_Gradient(POINTER params,int D,double *G,double *P);

extern double exact_eps(
	int D)
{
    	double v;
    	v = (D==1) ? eps_plus : eps_minus;
    	return v;
}

extern double exact_solution(
	POINTER params,
	int D, 
	double *P)
{
    	double x,y,z;
	CIM_PARAMS *cim_params = (CIM_PARAMS*)params;
	int dim = cim_params->dim;
	int Run_case = cim_params->Run_case;
	double soln;
	x = P[0];
	y = P[1];
	z = P[2];

	switch(Run_case)
	{
	case 1:
	    if (dim == 2)
	    {
		soln = x;
	    }
	    else if (dim == 3)
	    {
		soln = (D == 1) ? x*x+y*y+z*z+x*y+y*z+x*z :
			x*x+2*y*y+3*z*z+4*x*y+5*y*z+6*x*z;
	    }
	    break;
	case 2:
	    if (dim == 2)
	    {
		soln = (D == 1) ? x : x-1;
	    }
	    else if (dim == 3)
	    {
		soln = (D == 1) ? x*y+x*x*x*x+y*y*y*y+x*z*z+cos(2*x+y*y+z*z*z) :
			x*x*x+x*y*y+y*y*y+z*z*z*z+sin(3*(x*x+y*y));
	    }
	    break;
	case 3:
	    if (dim == 2)
	    {
		soln = (D == 1) ? x : y;
	    }
	    break;
	case 4:
	    if (dim == 2)
	    {
		soln = (D == 1) ? x*x+y*y+x*y+x*x*x : 4*x*x+3*y*y+5*x*y;
	    }
	    break;
	case 5:
	    if (dim == 2)
	    {
		soln = (D == 1) ? x*x+2*x*x*y+3*x*y*y+y*y : 
				  x*x*x+7*x*x*y+6*x*y*y+8*y*y*y;
	    }
	    break;
	case 6:
	    if (dim == 2)
	    {
		soln = (D == 1) ? 1+log(2*sqrt(x*x+y*y)) : 1;
	    }
	    break;
	case 7:
	    if (dim == 2)
	    {
	    	if (D == 1) 
	    	    soln = (sqrt(x*x+y*y)-0.5)*(sqrt(x*x+y*y)-0.5);
	    	else 
	    	    soln = 0;
	    }
	    break;
	case 8:
	    {
		double DH = sqrt(0.106086285*0.15);
		double r0 = 2.0;
        	double e = exact_eps(D);
                double r = sqrt(x*x+y*y+z*z);
        	double C1 =(exact_eps(-1)-exact_eps(1)*(1+DH*r0))*exp(-DH*r0)/
				(exact_eps(-1)*exact_eps(1)*r0);
        	double C2 = (1+DH*r0)*exp(-DH*r0);
		if (D == 1)
                    soln = 1/C2/(e*r)*exp(-DH*r);
                else 
                    soln = 1/(e*r)+ C1/C2;
	    }
	    break;
	}
	return soln;
}	/* end exact_solution */

static int exact_Gradient(
	POINTER params,
	int D, 
	double *G, 
	double *P)
{
    	int i;
    	double x,y,z;
	CIM_PARAMS *cim_params = (CIM_PARAMS*)params;
	int Run_case = cim_params->Run_case;
	int dim = cim_params->dim;
    
    	for (i = 0; i < dim; ++i) G[i] = 0;

    	x = P[0]; y = P[1]; z = P[2];
	switch (Run_case)
	{
	case 1:
	    if (dim == 2)
	    {
	    	G[0] = 1.0;
	    	G[1] = 0.0;
	    }
	    else if (dim == 3)
	    {
		if (D == 1)
		{
		    G[0] = 2*P[0]+1*P[1]+1*P[2];
            	    G[1] = 1*P[0]+2*P[1]+1*P[2];
            	    G[2] = 1*P[0]+1*P[1]+2*P[2];
		}
		else
		{
		    G[0] = 2*P[0]+4*P[1]+6*P[2];
            	    G[1] = 4*P[0]+4*P[1]+5*P[2];
            	    G[2] = 6*P[0]+5*P[1]+6*P[2];
		}
	    }
	    break;
	case 2:
	    if (dim == 2)
	    {
	    	G[0] = 1.0;
	    	G[1] = 0.0;
	    }
	    else if (dim == 3)
	    {
		if (D == 1)
		{
		    G[0] = y+4*x*x*x+z*z-sin(2*x+y*y+z*z*z)*2;
            	    G[1] = x+4*y*y*y-sin(2*x+y*y+z*z*z)*2*y;
            	    G[2] = 2*x*z-sin(2*x+y*y+z*z*z)*3*z*z;
		}
		else
		{
		    G[0] = 3*x*x+y*y+cos(3*(x*x+y*y))*6*x;
            	    G[1] = 2*x*y+3*y*y+cos(3*(x*x+y*y))*6*y;
            	    G[2] = 4*z*z*z;
		}
	    }
	    break;
	case 3:
	    if (D == 1)
	    {
	    	G[0] = 1.0;
	    	G[1] = 0.0;
	    }
	    else
	    {
	    	G[0] = 0.0;
	    	G[1] = 1.0;
	    }
	    break;
	case 4:
	    if (D == 1)
	    {
	    	G[0] = 2*P[0]+P[1]+3*P[0]*P[0];
	    	G[1] = 2*P[1]+P[0];
	    }
	    else
	    {
	    	G[0] = 8*P[0]+5*P[1];
	    	G[1] = 6*P[1]+5*P[0];
	    }
	    break;
	case 5:
	    if (D == 1)
	    {
	    	G[0] = 2*x+4*x*y+3*y*y;
	    	G[1] = 2*x*x+6*x*y+2*y;
	    }
	    else
	    {
	    	G[0] = 3*x*x+14*x*y+6*y*y;
	    	G[1] = 7*x*x+12*x*y+24*y*y;
	    }
	    break;
	case 6:
	    if (D == 1)
	    {
	    	G[0] = x/(x*x+y*y);
	    	G[1] = y/(x*x+y*y);
	    }
	    else
	    {
	    	G[0] = 0.0;
	    	G[1] = 0.0;
	    }
	    break;
	case 7:
	    if(D == 1) 
	    {
	    	G[0] = 2*(sqrt(x*x+y*y)-0.5)*x/sqrt(x*x+y*y);
	    	G[1] = 2*(sqrt(x*x+y*y)-0.5)*y/sqrt(x*x+y*y);
	    } 
	    else 
	    {
	    	G[0] = 0;
	    	G[1] = 0;
	    }
	    break;
	case 8:
	    {
		double DH = sqrt(0.106086285*0.15);
		double rsq = x*x+y*y+z*z;
        	double r0 = 2.0;
        	double e = exact_eps(D);        
                double r = sqrt(rsq);
        	double C2 = (1+DH*r0)*exp(-DH*r0);

                if(D == -1)
		{
                    G[0] -= x/r/rsq;
                    G[1] -= y/r/rsq;
                    G[2] -= z/r/rsq;                        
                } 
		else 
		{
                    G[0] -= x/C2/rsq*(1/r+DH)*exp(-DH*r); 
                    G[1] -= y/C2/rsq*(1/r+DH)*exp(-DH*r); 
                    G[2] -= z/C2/rsq*(1/r+DH)*exp(-DH*r);
                }
	    }
	    break;
	}
    	return YES;
}	/* exact_Gradient */

extern double exact_source(
	POINTER params,
	int D, 
	double *P)
{
	CIM_PARAMS *cim_params = (CIM_PARAMS*)params;
    	double x,y,z,v;
	int dim = cim_params->dim;
	int Run_case = cim_params->Run_case;

	x = P[0];
	y = P[1];
	z = P[2];
	switch (Run_case)
	{
	case 1:
	    if (dim == 2)
	    {
		v = 0.0;
	    }
	    else if (dim == 3)
	    {
		if (D == 1)
		    v = exact_eps(D)*(-6.0);
		else
		    v = exact_eps(D)*(-12.0);
	    }
	    break;
	case 2:
	    if (dim == 2)
	    {
		v = 0.0;
	    }
	    else if (dim == 3)
	    {
		if (D == 1)
		{
		    v = exact_eps(D)*(-12*x*x-12*y*y-2*x+
				cos(2*x+y*y+z*z*z)*(4+4*y*y+9*z*z*z*z)
                 		+sin(2*x+y*y+z*z*z)*(2+6*z));
		}
		else
		{
		    v = exact_eps(D)*(-8*x-6*y-12*z*z+
				36*(x*x+y*y)*sin(3*(x*x+y*y))
                 		-12*cos(3*(x*x+y*y)));
		}
	    }
	    break;
	case 3:
	    v = 0;
	    break;
	case 4:
	    if (D == 1)
		v = exact_eps(D)*(-4-6*P[0]);
	    else
		v = exact_eps(D)*(-14);
	    break;
	case 5:
	    if (D == 1)
		v = exact_eps(D)*(-4-6*P[0]-4*P[1]);
	    else
		v = exact_eps(D)*(-18*P[0]-62*P[1]);
	    break;
	case 6:
	    v = 0;
	    break;
	case 7:
	    if (D == 1) 
	    {
	    	v = exact_eps(D)*(1-4*sqrt(x*x+y*y))/sqrt(x*x+y*y);
	    } 
	    else 
	    {
	    	v = 0;
	    }
	    break;
	case 8:
	    {
        	double *h = cim_params->h;
                v = 0;
        	if (fabs(x)<0.5*h[0] && fabs(y)<0.5*h[1] && fabs(z)<0.5*h[2])
            	    v = 1;
	    }
	    break;
	}
    	return v;
}

extern double exact_jump_eps_gradu_dot_n(
	POINTER params,
	int D, 
	double *N, 
	double *P)
{
    	int i;
    	double v1, v2, G1[MAXD], G2[MAXD];
	CIM_PARAMS *jparams = (CIM_PARAMS*)params;
	int dim = jparams->dim;

    	exact_Gradient(params,D, G1, P);
    	exact_Gradient(params,-D, G2, P);
    	v1 = v2 = 0;
    	for (i = 0; i < dim; ++i) 
	{
            v1 += G1[i]*N[i];
            v2 += G2[i]*N[i];
    	}
    	return exact_eps(-D)*v2-exact_eps(D)*v1;
}

extern double exact_jump_u(
	POINTER params,
	int D, 
	double *P)
{
    	return exact_solution(params,-D,P)-exact_solution(params,D,P);
}

extern double exact_jump_gradu_dot_t(
	POINTER params,
	int D, 
	int i, 
	double *N, 
	double *P)
{
    	int j;
	double v1, v2, t[MAXD], G1[MAXD], G2[MAXD];
	CIM_PARAMS *jparams = (CIM_PARAMS*)params;
	int dim = jparams->dim;
	
	tangential_direction(i, t, N, dim);
    	exact_Gradient(params,D, G1, P);
    	exact_Gradient(params,-D, G2, P);
    	v1 = v2 = 0;
    	for (j = 0; j < dim; ++j) 
	{
            v1 += G1[j]*t[j];
            v2 += G2[j]*t[j];
	}
    	return v2-v1;
}

static int tangential_direction(
	int i, 
	double *T, 
	double *N,
	int dim)
{
    	int j;
    	double L;
    
    	for (j = 0; j < dim; ++j) 
	    T[j] = 0.0;
	L = sqrt(1.0-N[i]*N[i]);
	if ( L > CIM_TOL ) 
	{
	    for(j=0;j<dim;++j) 
	    {
		T[j] = -N[i]*N[j];
		T[j] += (i==j) ? 1 : 0;
		T[j] /= L;
	    }
	    j = 1;
	} 
	else 
	{
    	    j = 0;
    	}
    	return j;
}

extern double intfc_func_case2(
	POINTER params,
	double *P)
{
	double theta,r2,r,v;
	double rot_angle = 0.0;
	double Rbar = 0.5;
	double Rosc = 0.5;
	double Nosc = 6.0;

	theta = atan2(P[1],P[0])+rot_angle;
        r2 = P[0]*P[0]+P[1]*P[1];
        r  = sqrt(r2);
        v = r - Rbar*(1.0+Rosc*sin(Nosc*theta-PI/4.0));
	return v;
}	/* end intfc_func_case2 */

extern double intfc_func_case3(
	POINTER params,
	double *P)
{
	double r,theta,phi,v;

	r = sqrt(P[0]*P[0]+P[1]*P[1]+P[2]*P[2]);
        theta = atan2(P[1],P[0]);
        phi = acos(P[2]/r);
        v = r - 0.5 - 0.2*sin(2*theta)*sin(phi);
	return v;
}	/* end intfc_func_case3 */

extern double intfc_func_case4(
	POINTER params,
	double *P)
{
	double r,v;
	r = sqrt(P[0]*P[0]+P[1]*P[1]);
        v = (r-0.6)*(r-0.6) + P[2]*P[2] - 0.09;
	return v;
}	/* end intfc_func_case4 */

extern double intfc_func_case5(
	POINTER params,
	double *P)
{
	double x,y,z,v;
	x = 7*P[0]+6;
        y = 7*P[1];
        z = 7*P[2];
        if (P[0] >  1 || P[0] < -1 || P[1] >  1 || 
	    P[1] < -1 || P[2] >  1 || P[2] < -1) 
	    v = 1;
        else 
	    v = pow(x, 0.4e1) + pow(y, 0.4e1) + 0.15e1*pow(z, 0.4e1) + 
	    	0.2e1*x*x*y*y + 0.2e1*x*x*z*z + 0.2e1*y*y*z*z 
		- 0.94e2*x*x + 0.78e2*y*y - 0.94e2*z*z + 0.1521e4;
	return v;
}	/* end intfc_func_case5 */

extern double intfc_func_case6(
	POINTER params,
	double *P)
{
	double r,p,rc,v,x,y,z,xi[14],yi[14],zi[14];
	int i;
	rc = 0.6;
        xi[0] = yi[0] = xi[1] = yi[1] = 0;
        zi[0] = rc;
        zi[1] = -rc;
        for (i=0; i<5; ++i) 
	{
            xi[i+2] = rc*2*cos((2*i*PI)/5.0)/sqrt(5.0);
            yi[i+2] = rc*2*sin((2*i*PI)/5.0)/sqrt(5.0);
            zi[i+2] = rc/sqrt(5.0);
            xi[i+7] = rc*2*cos((2*i-1)*PI/5.0)/sqrt(5.0);
            yi[i+7] = rc*2*sin((2*i-1)*PI/5.0)/sqrt(5.0);
            zi[i+7] = -rc/sqrt(5.0);
        }
        r = sqrt(P[0]*P[0]+P[1]*P[1]+P[2]*P[2]);
        v = r-rc;
        for(i=0; i<12; ++i) 
	{
            x = P[0]-xi[i]; y = P[1]-yi[i]; z = P[2]-zi[i];
            v -= 2*exp(-(x*x+y*y+z*z)/0.04);
        }
	return v;
}	/* end intfc_func_case6 */

extern double intfc_func_case7(
	POINTER params,
	double *P)
{
	double v,r2,XC[20],YC[20],ZC[20],RC[20];
	int k,N_protein = 6;
	v = 1e10;
        for (k=0; k<N_protein; ++k) 
	{
            r2 = (P[0]-XC[k])*(P[0]-XC[k])+
                 (P[1]-YC[k])*(P[1]-YC[k])+
                 (P[2]-ZC[k])*(P[2]-ZC[k])-RC[k]*RC[k];
            if (r2 < v) 
		v = r2;
        }
}	/* end intfc_func_case7 */

extern double intfc_func_case8(
	POINTER params,
	double *P)
{
	double v;

	v = P[0]*P[0]+P[1]*P[1]+P[2]*P[2]-2.0*2.0;// A 0.5001*0.5001;
	return v;
}	/* end intfc_func_case8 */

extern double intfc_func_case9(
	POINTER params,
	double *P)
{
	double v;
	v = P[0]*P[0]*4.0+P[1]*P[1]*9.0+P[2]*P[2]*16.0-1.001;
	return v;
}	/* end intfc_func_case9 */

extern double intfc_func_case10(
	POINTER params,
	double *P)
{
	double v,r,x,y,z;
	double Px[18] = {-2.027,-1.669,-0.453, 0.751, 0.393,-0.823,-2.284,
		-2.888,-2.527,-1.426,-0.196,-0.701, 1.007, 1.612, 0.149, 
		1.251,-1.081,-0.575};
	double Py[18] = { 0.954, 0.234,-0.687, 0.148, 0.868, 1.788, 0.208, 
		1.617,-0.368, 0.980,-1.189,-1.441, 0.894,-0.515, 0.121, 1.470, 
		2.291, 2.543};
	double Pz[18] = {-0.651, 0.665, 0.441,-0.040,-1.356,-1.132,-1.418,
		-0.483,0.997, 1.435, 1.385,-0.320, 0.727,-0.208,-2.126,-1.688,
		-2.076,-0.371};
	double Pr[18] = {1.7,1.7,1.7,1.7,1.7,1.7,1.2,1.2,1.2,1.2,1.2,1.2,
		1.2,1.2,1.2,1.2,1.2,1.2};

	int i;
	v = 8.8271e-2;
        for (i=0; i<18; i++)
	{
            x = P[0]-Px[i]; y = P[1]-Py[i]; z = P[2]-Pz[i]; r = Pr[i];
            v -= exp(-2.0*(x*x+y*y+z*z)/(r*r))/(r*r*r);
        }
	return v;
}	/* end intfc_func_case10 */

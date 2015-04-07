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


/*
*				fvola.c:
*
*	Application of FronTier to Black-Scholes Equation.
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*	
*/

#include "finance.h"

struct _C_PARAMS {
	double sigma0;
};
typedef struct _C_PARAMS C_PARAMS;

struct _L_PARAMS {
	double sigma0;
	double k;		// slope
};
typedef struct _L_PARAMS L_PARAMS;

struct _Q_PARAMS {
	double sigma0;
	double k;
};
typedef struct _Q_PARAMS Q_PARAMS;

struct _S_PARAMS {
	double sigma0;
	double A;		// amplitude
	double nu;		// frequency
	double phi0;		// initial phase
};
typedef struct _S_PARAMS S_PARAMS;

struct _H_PARAMS {
	double sigma0;
};
typedef struct _H_PARAMS H_PARAMS;

static double constant_volatility(Front*,POINTER);
static double linear_volatility(Front*,POINTER);
static double quadratic_volatility(Front*,POINTER);
static double sine_volatility(Front*,POINTER);
static double heston_volatility(Front*,POINTER);

extern void getVolatility(Front *front)
{
	PARAMS *eqn_params = (PARAMS*)front->extra1;
	POINTER vparams = eqn_params->vparams;
	switch (eqn_params->v_type)
	{
	case CONSTANT_V:
	    eqn_params->sigma[0] = constant_volatility(front,vparams);
	    return;
	case LINEAR_V:
	    eqn_params->sigma[0] = linear_volatility(front,vparams);
	    return;
	case QUADRATIC_V:
	    eqn_params->sigma[0] = quadratic_volatility(front,vparams);
	    return;
	case SINE_V:
	    eqn_params->sigma[0] = sine_volatility(front,vparams);
	    return;
	case HESTON_V:
	    eqn_params->sigma[0] = heston_volatility(front,vparams);
	    return;
	}
}	/* end getVolatility */

extern void initVolatility(
	char *inname,
	Front *front)
{
	static C_PARAMS c_params;
	static L_PARAMS l_params;
	static Q_PARAMS q_params;
	static S_PARAMS s_params;
	static H_PARAMS h_params;
	char string[100];

	PARAMS *eqn_params = (PARAMS*)front->extra1;
	FILE *infile = fopen(inname,"r");
	CursorAfterString(infile,"Enter volatility function type:");
        fscanf(infile,"%s",string);
        (void) printf("%s\n",string);
	switch (string[0])
	{
	case 'c':
	case 'C':
	    eqn_params->v_type = CONSTANT_V;
	    eqn_params->vparams = (POINTER)&c_params;
	    CursorAfterString(infile,"Enter sigma0:");
            fscanf(infile,"%lf",&c_params.sigma0);
            (void) printf("%f\n",c_params.sigma0);
	    break;
	case 'l':
	case 'L':
	    eqn_params->v_type = LINEAR_V;
	    eqn_params->vparams = (POINTER)&l_params;
	    CursorAfterString(infile,"Enter sigma0:");
            fscanf(infile,"%lf",&l_params.sigma0);
            (void) printf("%f\n",l_params.sigma0);
	    CursorAfterString(infile,"Enter the slope:");
            fscanf(infile,"%lf",&l_params.k);
            (void) printf("%f\n",l_params.k);
	    break;
	case 'q':
	case 'Q':
	    eqn_params->v_type = QUADRATIC_V;
	    eqn_params->vparams = (POINTER)&q_params;
	    CursorAfterString(infile,"Enter sigma0:");
            fscanf(infile,"%lf",&q_params.sigma0);
            (void) printf("%f\n",q_params.sigma0);
	    CursorAfterString(infile,"Enter parabolic coefficient:");
            fscanf(infile,"%lf",&q_params.k);
            (void) printf("%f\n",q_params.k);
	    break;
	case 's':
	case 'S':
	    eqn_params->v_type = SINE_V;
	    eqn_params->vparams = (POINTER)&s_params;
	    CursorAfterString(infile,"Enter sigma0:");
            fscanf(infile,"%lf",&s_params.sigma0);
            (void) printf("%f\n",s_params.sigma0);
	    CursorAfterString(infile,"Enter amplitude:");
            fscanf(infile,"%lf",&s_params.A);
            (void) printf("%f\n",s_params.A);
	    CursorAfterString(infile,"Enter frequency:");
            fscanf(infile,"%lf",&s_params.nu);
            (void) printf("%f\n",s_params.nu);
	    CursorAfterString(infile,"Enter phase:");
            fscanf(infile,"%lf",&s_params.phi0);
            (void) printf("%f\n",s_params.phi0);
	    break;
	case 'h':
	case 'H':
	    eqn_params->v_type = HESTON_V;
	    eqn_params->vparams = (POINTER)&h_params;
	    break;
	default:
	    (void) printf("Unknown type of volatility function!\n");
	    (void) clean_up(ERROR);
	}
	fclose(infile);
}	/* end initVolatility */

static double constant_volatility(
	Front *front,
	POINTER params)
{
	C_PARAMS *cparams = (C_PARAMS*)params;
	return cparams->sigma0;
}	/* end constant_volatility */

static double linear_volatility(
	Front *front,
	POINTER params)
{
	double time = front->time;
	L_PARAMS *lparams = (L_PARAMS*)params;
	return lparams->sigma0 + lparams->k*time;
}	/* end linear_volatility */

static double quadratic_volatility(
	Front *front,
	POINTER params)
{
	Q_PARAMS *qparams = (Q_PARAMS*)params;
	double time = front->time;
	return qparams->sigma0 + qparams->k*sqr(time);
}	/* end quadratic_volatility */

static double sine_volatility(
	Front *front,
	POINTER params)
{
	S_PARAMS *sparams = (S_PARAMS*)params;
	double time = front->time;
	double sigma0 = sparams->sigma0;
	double A = sparams->A;
	double nu = sparams->nu;
	double phi0 = sparams->phi0;
	return sigma0 + A*sin(2.0*nu*PI*time + phi0*PI/180.0);
}	/* end sine_volatility */

static double heston_volatility(
	Front *front,
	POINTER params)
{
}	/* end heston_volatility */

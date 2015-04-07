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
 * weno.c
 *
 *  Created on: Jun 12, 2011
 *      Author: yli
 */

#include "weno.h"

double (*flux_func)(double,double);
double (*dF)(double,double);

static double linear_flux(double,double);
static double burger_flux(double,double);
static double linear_dF(double,double);
static double burger_dF(double,double);
static double wave_speed;

void setEquationType(
	PARAMS params)
{
	switch (params.eqn_type)
	{
	case LINEAR_EQN:
	    flux_func = linear_flux;
	    dF = linear_dF;
	    wave_speed = params.a;
	    break;
	case BURGER_EQN:
	    flux_func = burger_flux;
	    dF = burger_dF;
	    break;
	default:
	    (void) printf("Unknown equation type!\n");
	    clean_up(ERROR);
	}
}	/* end setEquationType */

static double linear_flux(
	double wave_speed,
	double u)
{
        return wave_speed*u;
}

static double linear_dF(
	double wave_speed,
	double u)
{
        return wave_speed;
}

static double burger_flux(
	double wave_speed,
	double u)
{
        return 0.5*u*u;
}

static double burger_dF(
	double wave_speed,
	double u)
{
        return u;
}

void Weno5(
	int mesh_size, 
	double *u_old, 
	double *u_new,
	double dx, 
	double dt) 
{
	/* use fourth-order Runge-Kutta scheme */
	Runge_Kutta_4th(mesh_size,u_old,u_new,dx,dt);
//	TVD_RK_3th(mesh_size,u_old,u_new,dx,dt);
}

void TVD_RK_3th(
	int mesh_size, 
	double *u_old, 
	double *u_new,
	double dx, 
	double dt) 
{
	double *u1, *u2;
	double *flux;
	int i;
	int nrad = 3; /* 5th-order weno */
	int extend_size = mesh_size + 2 * nrad;
	double *u_extend;
	double lambda = -dt/dx;

	FT_VectorMemoryAlloc((POINTER*)&u1,mesh_size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&u2,mesh_size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&u_extend,extend_size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&flux,mesh_size,sizeof(double));

	/* Set the value on extended mesh */
	for (i = 0; i < mesh_size; i++) {
	    u_extend[i + nrad] = u_old[i];
	}
	for (i = 1; i <= nrad; i++) {
	    u_extend[nrad - i] = u_old[mesh_size - 1 - i];
	    u_extend[nrad + mesh_size - 1 + i] = u_old[i];
	}

	Weno5_Get_Flux(u_extend,flux,lambda,mesh_size);
	for (i = 0; i < mesh_size; i++)
		u1[i] = u_old[i] + flux[i];

	for (i = 0; i < mesh_size; i++) {
	    u_extend[i + nrad] = u1[i];
	}
	for (i = 1; i <= nrad; i++) {
	    u_extend[nrad - i] = u1[mesh_size - 1 - i];
	    u_extend[nrad + mesh_size - 1 + i] = u1[i];
	}

	Weno5_Get_Flux(u_extend,flux,lambda,mesh_size);
	for (i = 0; i < mesh_size; i++)
		u2[i] = 0.75*u_old[i] + 0.25*u1[i] + 0.25*flux[i];

	for (i = 0; i < mesh_size; i++) {
	    u_extend[i + nrad] = u2[i];
	}
	for (i = 1; i <= nrad; i++) {
	    u_extend[nrad - i] = u2[mesh_size - 1 - i];
	    u_extend[nrad + mesh_size - 1 + i] = u2[i];
	}

	Weno5_Get_Flux(u_extend,flux,lambda,mesh_size);
	for (i = 0; i < mesh_size; i++)
		u_new[i] = 1.0/3.0*u_old[i] + 2.0/3.0*u2[i] + 2.0/3.0*flux[i];
	FT_FreeThese(4,u1,u2,u_extend,flux);
}

void Runge_Kutta_4th(
	int mesh_size, 
	double *u_old, 
	double *u_new,
	double dx, 
	double dt) 
{
	double *u1, *u2, *u3;
	double *flux;
	int i;
	int nrad = 3; /* 5th-order weno */
	int extend_size = mesh_size + 2 * nrad;
	double *u_extend;
	double lambda = -dt/dx;

	FT_VectorMemoryAlloc((POINTER*)&u1,mesh_size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&u2,mesh_size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&u3,mesh_size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&u_extend,extend_size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&flux,mesh_size,sizeof(double));

	/* Set the value on extended mesh */
	for (i = 0; i < mesh_size; i++) {
	    u_extend[i + nrad] = u_old[i];
	}
	for (i = 1; i <= nrad; i++) {
	    u_extend[nrad - i] = u_old[mesh_size - 1 - i];
	    u_extend[nrad + mesh_size - 1 + i] = u_old[i];
	}

	Weno5_Get_Flux(u_extend,flux,lambda,mesh_size);
	for (i = 0; i < mesh_size; i++)
		u1[i] = u_old[i] + 0.5*flux[i];

	/* Set the value on extended mesh */
	for (i = 0; i < mesh_size; i++) {
	    u_extend[i + nrad] = u1[i];
	}
	for (i = 1; i <= nrad; i++) {
	    u_extend[nrad - i] = u1[mesh_size - 1 - i];
	    u_extend[nrad + mesh_size - 1 + i] = u1[i];
	}
	Weno5_Get_Flux(u_extend,flux,lambda,mesh_size);
	for (i = 0; i < mesh_size; i++)
		u2[i] = u_old[i] + 0.5*flux[i];

	/* Set the value on extended mesh */
	for (i = 0; i < mesh_size; i++) {
	    u_extend[i + nrad] = u2[i];
	}
	for (i = 1; i <= nrad; i++) {
	    u_extend[nrad - i] = u2[mesh_size - 1 - i];
	    u_extend[nrad + mesh_size - 1 + i] = u2[i];
	}
	Weno5_Get_Flux(u_extend,flux,lambda,mesh_size);
	for (i = 0; i < mesh_size; i++)
		u3[i] = u_old[i] + flux[i];

	/* Set the value on extended mesh */
	for (i = 0; i < mesh_size; i++) {
	    u_extend[i + nrad] = u3[i];
	}
	for (i = 1; i <= nrad; i++) {
	    u_extend[nrad - i] = u3[mesh_size - 1 - i];
	    u_extend[nrad + mesh_size - 1 + i] = u3[i];
	}
	Weno5_Get_Flux(u_extend,flux,lambda,mesh_size);
	for (i = 0; i < mesh_size; i++)
	    u_new[i] = 1.0/3*(- u_old[i] + u1[i] + 2.0*u2[i] + u3[i])
		     + 1.0/6*flux[i];
	FT_FreeThese(5,u1,u2,u3,u_extend,flux);
}

void Weno5_Get_Flux(
	double *u, 
	double *flux, 
	double lambda,
	int n)
{
	int nrad = 3; /* 5th-order weno */
	int i,j,k;
	int extend_size = n + 2*nrad;
	const double eps = 1.e-8;
	const int p = 2;
	static double *f; /* f(u(x)) */
	static double *fp, *fm;
	static int max_n = 0;
	double max_df, norm;

	/* coefficients for 2rd-order ENO interpolation stencil */
	double a[3][3] = {{1.0/3.0, -7.0/6.0, 11.0/6.0}, 
			  {-1.0/6.0, 5.0/6.0, 1.0/3.0}, 
			  {1.0/3.0, 5.0/6.0, -1.0/6.0}};

	/* Optimal weights C_k */
	double c[3] = {0.1,0.6,0.3};

	double is[3]; /* a smoothness measurement of the flux function */
	double alpha[3];
	double omega[3]; /* weights for stencils */
	double q[3]; /* ENO approximation of flux */
	double sum;

	if (max_n < n)
        {
            max_n = n;
            if (f != NULL)
            {
                FT_FreeThese(3,f,fp,fm);
            }
            FT_VectorMemoryAlloc((POINTER*)&f,extend_size,sizeof(double));
            FT_VectorMemoryAlloc((POINTER*)&fp,n+1,sizeof(double));
            FT_VectorMemoryAlloc((POINTER*)&fm,n+1,sizeof(double));
        }

	/* Find the maximum of fabs(df(u))*/
	max_df = 0;
	for (i = 0; i < extend_size; i++) 
	{
	    norm = dF(wave_speed,u[i]) > 0 ? dF(wave_speed,u[i]) : -dF(wave_speed,u[i]);
	    max_df = max_df > norm ? max_df : norm;
	}

	/* f[i] = 0.5 * (f_{i - nrad} + max_df * u[i])
	 * ensure that df[i] > 0*/
	for (i = 0; i < extend_size; i++)
	    f[i] = 0.5 * (flux_func(wave_speed,u[i]) + max_df * u[i]);

	for (j = 0; j < n + 1; j++)
	/* To approximate flux at x_{j+1/2}, we use stencil
	 S_k = {x_{j+k-2}, x_{j+k-1}, x_{j+k}} = {x[j+k-2+nrad],
	 x[j+k-1+nrad], x[j+k+nrad]} */
	{
		/* compute the weights omega[] */
	    is[0] = 13.0/12.0*(f[j] - 2.0*f[j+1] + f[j+2])*(f[j] - 2.0
		    *f[j+1] + f[j+2]) + 0.25*(f[j] - 4.0*f[j+1] + 3.0
		    *f[j+2])*(f[j] - 4.0*f[j+1] + 3.0*f[j+2]);
	    is[1] = 13.0/12.0*(f[j+1] - 2.0*f[j+2] + f[j+3])
		    *(f[j+1] - 2.0*f[j+2] + f[j+3]) + 0.25*(f[j+1]
		    -f[j+3])*(f[j+1] - f[j+3]);
	    is[2] = 13.0/12.0*(f[j+2] - 2.0*f[j+3] + f[j+4])
		    *(f[j+2] - 2.0*f[j+3] + f[j+4]) + 0.25*(3.0*f[j+2] 
		    - 4.0*f[j+3] + f[j+4])*(3.0*f[j+2] - 4.0*f[j+3] + f[j+4]);

	    for (i = 0; i < nrad; i++)
		alpha[i] = c[i]/pow(eps + is[i],p);

	    sum = alpha[0] + alpha[1] + alpha[2];
	    for (i = 0; i < nrad; i++)
		omega[i] = alpha[i]/sum;

		/* compute ENO approximation of the flux */
	    for (i = 0; i < nrad; i++) 
	    {
		q[i] = 0.0;
		for (k = 0; k < nrad; k++)
		    q[i] += a[i][k] * f[j + i + k];
	    }

	    /* compute the linear combination of the r candidate
	     stencils to get a higher order
	     approximation of the flux */
	    fp[j] = 0.0;
	    for (i = 0; i < nrad; i++)
		fp[j] += omega[i]*q[i];
	}

	/* f[i] = 0.5 * (f_{i - nrad} - max_df * u[i])
	 * ensure that df[i] < 0*/
	for (i = 0; i < extend_size; i++)
	    f[i] = 0.5*(flux_func(wave_speed,u[i]) - max_df*u[i]);
	for (j = 0; j < n + 1; j++)
	/* To approximate flux at x_{j+1/2}, we use stencil S_k =
	 {x_{j+1-k+2}, x_{j+1-k+1}, x_{j+1-k}} = {x[j+1-k+2+nrad],
	 x[j+1-k+1+nrad], x[j+1-k+nrad]} */

	{
	    /* compute the weights omega[] */
	    is[0] = 13.0/12.0*(f[j+5] - 2.0*f[j+4] + f[j+3])
		    *(f[j+5] - 2.0*f[j+4] + f[j+3]) + 0.25*(f[j+5]
		    - 4.0*f[j+4] + 3.0*f[j+3])*(f[j+5] - 4.0*f[j+4]
		    + 3.0*f[j+3]);
	    is[1] = 13.0/12.0*(f[j+4] - 2.0*f[j+3] + f[j+2])
		    *(f[j+4] - 2.0*f[j+3] + f[j+2]) + 0.25*(f[j+4]
		    - f[j+2])*(f[j+4] - f[j+2]);
	    is[2] = 13.0/12.0*(f[j+3] - 2.0*f[j+2] + f[j+1])
		    *(f[j+3] - 2.0*f[j+2] + f[j+1]) + 0.25*(3.0*f[j+3] 
		    - 4.0*f[j+2] + f[j+1])*(3.0*f[j+3] - 4.0*f[j+2] + f[j+1]);

	    for (i = 0; i < nrad; i++)
		alpha[i] = c[i]/pow(eps + is[i], p);

	    sum = alpha[0] + alpha[1] + alpha[2];
	    for (i = 0; i < nrad; i++)
		omega[i] = alpha[i]/sum;

		/* compute ENO approximation of the flux */
	    for (i = 0; i < nrad; i++) 
	    {
		q[i] = 0.0;
		for (k = 0; k < nrad; k++)
		    q[i] += a[i][k]*f[j+5-i-k];
	    }

		/* compute the linear combination of the r candidate stencils
		 to get a higher order approximation of the flux */
	    fm[j] = 0.0;
	    for (i = 0; i < nrad; i++)
		fm[j] += omega[i]*q[i];
	}

	for (j = 0; j < n; j++)
	    flux[j] = lambda*(fp[j+1] + fm[j+1] - fp[j] - fm[j]);

}

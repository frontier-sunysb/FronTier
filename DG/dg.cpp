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
 * gd.c
 *
 *  Created on: July 1, 2015
 *      Author: xlchen
 */

#include "dg.h"

double (*flux_func)(double,double);
double (*dF)(double,double);

static double linear_flux(double,double);
static double burger_flux(double,double);
static double linear_dF(double,double);
static double burger_dF(double,double);
static double wave_speed;
static double min_mod(double,double,double);

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

void Dg5(
	int mesh_size, 
	double **coef_old, 
	double **coef_new,
	double dx, 
	double dt) 
{
	/* use third-order TVD Runge-Kutta scheme */
	TVD_RK_3th(mesh_size,coef_old,coef_new,dx,dt);
}

void TVD_RK_3th(
	int mesh_size, 
	double **coef_old, 
	double **coef_new,
	double dx,
 	double dt) 
{
	double **coef1, **coef2;
	int i, j;
	int nrad = 1;
	int extend_size = mesh_size + 2 * nrad;
	double **coef_extend;
	double **rhs;

	FT_MatrixMemoryAlloc((POINTER*)&coef_extend,
				extend_size,5,sizeof(double));
	FT_MatrixMemoryAlloc((POINTER*)&coef1,mesh_size,5,sizeof(double));
	FT_MatrixMemoryAlloc((POINTER*)&coef2,mesh_size,5,sizeof(double));
	FT_MatrixMemoryAlloc((POINTER*)&rhs,mesh_size,5,sizeof(double));

	/* Set the value on extended mesh */
	for (i = 0; i < mesh_size; ++i) 
	    for (j = 0; j < 5; ++j)
		coef_extend[i + nrad][j] = coef_old[i][j];
	for (i = 1; i <= nrad; ++i) 
	{
	    for (j = 0; j < 5; ++j)
	    {
		coef_extend[nrad - i][j] = coef_old[mesh_size - i][j];
		coef_extend[nrad + mesh_size - 1 + i][j] = coef_old[i - 1][j];
	    }
	}

	Dg5_rhs_eval(mesh_size,coef_extend,rhs,dx);
	for (i = 0; i < mesh_size; ++i)
            for (j = 0; j < 5; ++j)
                coef_old[i][j] = coef_extend[i + nrad][j];
	for (i = 0; i < mesh_size; ++i)
	    for (j = 0; j < 5; ++j)
		coef1[i][j] = coef_old[i][j] + dt * rhs[i][j];

        for (i = 0; i < mesh_size; ++i)
            for (j = 0; j < 5; ++j)
                coef_extend[i + nrad][j] = coef1[i][j];
	for (i = 1; i <= nrad; ++i)
        {
            for (j = 0; j < 5; ++j)
            {   
                coef_extend[nrad - i][j] = coef1[mesh_size - i][j];
                coef_extend[nrad + mesh_size - 1 + i][j] = coef1[i - 1][j];
            }
        }

	Dg5_rhs_eval(mesh_size,coef_extend,rhs,dx);
	for (i = 0; i < mesh_size; ++i)
            for (j = 0; j < 5; ++j)
                coef1[i][j] = coef_extend[i + nrad][j];
        for (i = 0; i < mesh_size; ++i)
            for (j = 0; j < 5; ++j)
                coef2[i][j] = 0.75*coef_old[i][j] + 0.25 *coef1[i][j] + 
				 0.25 * dt * rhs[i][j];

	for (i = 0; i < mesh_size; ++i)
            for (j = 0; j < 5; ++j)
                coef_extend[i + nrad][j] = coef2[i][j];
        for (i = 1; i <= nrad; ++i)
        {
            for (j = 0; j < 5; ++j)
            {
                coef_extend[nrad - i][j] = coef2[mesh_size - i][j];
                coef_extend[nrad + mesh_size - 1 + i][j] = coef2[i - 1][j];
            }
        }

        Dg5_rhs_eval(mesh_size,coef_extend,rhs,dx);
	for (i = 0; i < mesh_size; ++i)
            for (j = 0; j < 5; ++j)
                coef2[i][j] = coef_extend[i + nrad][j];
        for (i = 0; i < mesh_size; ++i)
            for (j = 0; j < 5; ++j)
                coef_new[i][j] = 1.0/3*coef_old[i][j] + 2.0/3 *coef2[i][j] +
                                 2.0/3 * dt * rhs[i][j];

	FT_FreeThese(4,coef1,coef2,coef_extend,rhs);
}	/* TVD_RK_3rd */

void Dg5_rhs_eval(
	int mesh_size,
	double **coef, 
	double **rhs,
	double dx)
{
	int nrad = 1;
	int i,j,l;
	double ugau[5];
	double xgau[5] = {0.0,-sqrt(5.0 - 2.0*sqrt(10.0/7.0))/3.0,
                        sqrt(5.0 - 2.0*sqrt(10.0/7.0))/3.0,
                        -sqrt(5.0 + 2.0*sqrt(10.0/7.0))/3.0,
                        sqrt(5.0 + 2.0*sqrt(10.0/7.0))/3.0};
        double wgau[5] = {128.0/225.0, (322.0+13.0*sqrt(70.0))/900.0,
                        (322.0+13.0*sqrt(70.0))/900.0,
                        (322.0-13.0*sqrt(70.0))/900.0,
                        (322.0-13.0*sqrt(70.0))/900.0};
        double AI[5][5] = {{225.0/64/dx, 0.0, -525.0/8/dx, 0.0, 945.0/4/dx},
                           {0.0, 75.0/dx, 0.0, -420.0/dx, 0.0},
                           {-525.0/8/dx, 0.0, 2205.0/dx, 0.0, -9450.0/dx},
                           {0.0, -420.0/dx, 0.0, 2800.0/dx, 0.0},
                           {945.0/4/dx, 0.0, -9450.0/dx, 0.0, 44100.0/dx}};
	double **tmp;
	double *uleft, *uright,*fleft,*fright,*flux;
	double *uave;

	FT_MatrixMemoryAlloc((POINTER*)&tmp,mesh_size,5,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&uleft,mesh_size + 2,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&uright,mesh_size + 2,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&fleft,mesh_size + 2,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&fright,mesh_size + 2,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&flux,mesh_size + 1,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&uave,mesh_size + 2*nrad + 2,
							sizeof(double));

	/* calculate the volume integral terms */
	for (i = 0; i < mesh_size; ++i)
	{
	    ugau[0] = coef[i+nrad][0];
	    for (l = 1; l < 5; ++l) /* 5 gaussian points */
	    {
		ugau[l] = 0.0;
		for (j = 0; j < 5; ++j) /* 5th order */
		    ugau[l] += coef[i+nrad][j] * pow(0.5 * xgau[l], j);
	    }
	    tmp[i][0] = 0.0;
	    for (j = 1; j < 5; ++j)
	    {
		tmp[i][j] = 0.0;
		for (l = 1; l < 5; ++l)
		    tmp[i][j] += wgau[l] * flux_func(wave_speed,ugau[l]) * 
				pow(0.5 * xgau[l], j - 1) * j / dx;
		if (j == 1)
		    tmp[i][j] += wgau[0] * flux_func(wave_speed,ugau[0]) / dx;
		tmp[i][j] *= dx * 0.5;
	    }
	    /* calculate the cell averages */
	    uave[i + nrad + 1] = 0.0;
	    for (l = 0; l < 5; ++l) /* 5 gaussian points */
		uave[i + nrad + 1] += wgau[l] * ugau[l];
	    uave[i + nrad + 1] *= 0.5;
	}
	/* periodic boundary for the cell averages */
	for (i = 1; i <= nrad; ++i)
        {
	    uave[nrad - i + 1] = uave[mesh_size + nrad - i + 1];
	    uave[mesh_size + nrad + i] = uave[nrad + i];
        }
	uave[0] = uave[mesh_size];
	uave[mesh_size+2*nrad+1] = uave[2*nrad+1];
	/* calculate the cell boundary terms */
	for (i = 0; i < mesh_size + 2; ++i)
	{
	    uleft[i] = 0.0;
	    uright[i] = 0.0;
	    for (j = 0; j < 5; ++j)
	    {
		uleft[i] += coef[i+nrad-1][j] * pow(0.5, j);
		uright[i] += coef[i+nrad-1][j] * pow(-0.5, j);
	    }
	}
	/* TVB limiter */
	for (i = 0; i < mesh_size + 2; ++i)
	{
	    double umod_l, umod_r;
	    double M = 40.0;
	    umod_l = uleft[i] - uave[i + nrad];
	    umod_r = uave[i + nrad] - uright[i];
	    boolean changed = NO;
	    if (fabs(umod_l) > M * dx * dx)
	    {
		umod_l = min_mod(umod_l,uave[i + nrad + 1] - uave[i + nrad],
				uave[i + nrad] - uave[i + nrad - 1]);
		changed = YES;
	    }
	    if (fabs(umod_r) > M * dx * dx)
	    {
                umod_r = min_mod(umod_r,uave[i + nrad + 1] - uave[i + nrad],
                                uave[i + nrad] - uave[i + nrad - 1]);
		changed = YES;
	    }
	    uleft[i] = uave[i + nrad] + umod_l;
	    uright[i] = uave[i + nrad] - umod_r;
	    if (changed)
	    {
		coef[i][1] = uleft[i] - uright[i];
		coef[i][2] = 2.0 * (uleft[i] - 2.0 * coef[i][0] + uright[i]);
		coef[i][3] = 0.0;
		coef[i][4] = 0.0;
	    }
	}
	/* end of TVB limiter */
	for (i = 0; i < mesh_size + 2; ++i)
	{
	    fleft[i] = flux_func(wave_speed,uleft[i]);
	    fright[i] = flux_func(wave_speed,uright[i]);
	}
	double alpha = fabs(dF(wave_speed,coef[nrad][0]));
	for (int i = nrad + 1; i < mesh_size + nrad; ++i)
	    if (fabs(dF(wave_speed,coef[i][0])) > alpha)
		alpha = fabs(dF(wave_speed,coef[i][0]));
	for (i = 0; i < mesh_size + 1; ++i)
	    flux[i] = 0.5 *(fright[i + 1] + fleft[i] - 
			alpha * (uright[i + 1] - uleft[i]));
	for (i = 0; i < mesh_size; ++i)
	{
	    for (j = 0; j < 5; ++j)
	    {
		tmp[i][j] -= flux[i + 1] * pow(0.5, j);
		tmp[i][j] += flux[i] * pow(-0.5, j);
	    }
	}
	/* multiply the inverse matrix */
	for (i = 0; i < mesh_size; ++i)
        {
	    for (j = 0; j < 5; ++j)
	    {
		rhs[i][j] = 0.0;
		for (l = 0; l < 5; ++l)
		    rhs[i][j] += AI[j][l] * tmp[i][l];
	    }
	}
	FT_FreeThese(7,tmp,uleft,uright,fleft,fright,flux,uave);
}	/* Dg5_rhs_eval */

double min_mod(
	double a,
	double b,
	double c)
{
	double tmp;
	if (a > 0.0 && b > 0.0 && c > 0.0)
	{
	    tmp = a < b ? a : b;
	    return tmp < c ? tmp : c;
	}
	else if (a < 0.0 && b < 0.0 && c < 0.0)
	{
	    tmp = a > b ? a : b;
	    return tmp > c? tmp : c;
	}
	else
	    return 0.0;
}	/* min_mod function */

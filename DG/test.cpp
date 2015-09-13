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
 * test.c
 *
 *  Created on: July 1, 2015
 *      Author: xlchen
 */

#include "dg.h"

static void readDgParams(PARAMS*,char*);
static void setSolutionType(INIT_TYPE);

char *in_name,*out_name;

static void (*exact_soln)(double,double,double,double &);
static void wave_func(double,double,double,double &);
static void hamp_func(double,double,double,double &);
static void cosine_func(double,double,double,double &);
static void square_func(double,double,double,double &);
static void error_func(double*,double*,double*,int);
static void L2_Projection(double,double,double,double*);

int main(int argc, char **argv) 
{
	F_BASIC_DATA f_basic; /* Only some basic features of f_basic is used */
        Front front;          /* Only some basic feature of front is used */
        double *x,*u_old,*u_sol; /* mesh points and solutions */
	double **coef_old, **coef_new;
        double *x_movie,*y_movie;
        double dx,dt;
        char movie_caption[100];
        char time_label[100];
        char xg_name[200];
        char gd_name[200];
        FILE *xg_file;
        double xmin,xmax,umin,umax,height;
	double L,U;
        double CFL;     /* CFL safety factor */
        int buffer_size = 1;
	int i,imin,imax,mesh_size,expanded_mesh_size;
	int j;
	int boundary_type[2];
	PARAMS params;
	double wave_speed;
	double err[3];

	FT_Init(argc,argv,&f_basic);
        in_name         = f_basic.in_name;
        out_name        = f_basic.out_name;

	FT_ReadSpaceDomain(in_name,&f_basic);
	L = f_basic.L[0];
        U = f_basic.U[0];
        mesh_size = f_basic.gmax[0];
	boundary_type[0] = f_basic.boundary[0][0];
        boundary_type[1] = f_basic.boundary[0][1];
        imin = 0;
        imax = mesh_size;

	FT_ReadTimeControl(in_name,&front);
        CFL = Time_step_factor(&front);
        FT_ResetTime(&front);
	readDgParams(&params,in_name);
	setEquationType(params);
	setSolutionType(params.init_type);
	wave_speed = params.a;

	/* Allocating memory for arrays */
        expanded_mesh_size = mesh_size + 2*buffer_size;
        FT_VectorMemoryAlloc((POINTER*)&u_old,mesh_size,
                                        sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&u_sol,mesh_size,
                                        sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&x,mesh_size,
                                        sizeof(double));
        FT_MatrixMemoryAlloc((POINTER*)&coef_old,mesh_size,5,
                                        sizeof(double));
        FT_MatrixMemoryAlloc((POINTER*)&coef_new,mesh_size,5,
                                        sizeof(double));
	dx = (U - L)/mesh_size;
        FT_TimeControlFilter(&front);

	/* Set the initial condition */
        for (i = imin; i < imax; i++)
            x[i] = L + (i + 0.5)*dx;
	for (i = imin; i < imax; ++i)
	{
	    exact_soln(wave_speed,0.0,x[i],u_old[i]);
	    L2_Projection(wave_speed,x[i],dx,coef_old[i]);
	}
	
	/* Set frame margin for GD movie output */
        xmin = L;       xmax = U;
        umin = HUGE;    umax = -HUGE;
        for (i = 0; i < mesh_size; i++)
        {
            if (umin > u_old[i]) umin = u_old[i];
            if (umax < u_old[i]) umax = u_old[i];
        }
        height = umax - umin;
        umin -= 0.15*height;    umax += 0.15*height;

        /* Open and initialize GD movie output */
        sprintf(movie_caption,"u vs. x");
        sprintf(gd_name,"%s/soln.gif",out_name);
        gd_initplot(gd_name,movie_caption,xmin,xmax,umin,umax,2);

	/* Time loop */
	front.dt = dt = CFL*pow(dx,2.0)/wave_speed;
        for (;;)
        {
            /* Advancing numerical solution */
	    printf("dx = %f  dt = %f\n",dx,dt);
	    Dg5(mesh_size,coef_old,coef_new,dx,front.dt);

            /* Swapping solution storage */
            for (i = imin; i < imax; ++i)
		for (j = 0; j < 5; ++j)
                    coef_old[i][j] = coef_new[i][j];
	    for (i = imin; i < imax; ++i)
		u_old[i] = coef_old[i][0];

            /* Time and step control */
            FT_AddTimeStepToCounter(&front);
            (void) printf("\ntime = %20.14f   step = %5d   ",
                                front.time,front.step);

	    /* Movie frame */
            if (FT_IsDrawTime(&front))
            {
                /* Numerical solution */
                x_movie = x;
                y_movie = u_old;
                gd_plotdata(mesh_size,x_movie,y_movie);

                /* Exact solution */
		for (i = imin; i < imax; ++i)
                    exact_soln(wave_speed,front.time,x[i],u_sol[i]);
                x_movie = x;
                y_movie = u_sol;
                //gd_plotdata(mesh_size,x_movie,y_movie);

                /* Time label */
                sprintf(time_label,"Time = %6.3f",front.time);
                gd_plotframe(time_label);
            }
            /* Output date control */
            if (FT_IsSaveTime(&front))
            {
                /* Numerical solution */
		sprintf(xg_name,"%s/num_sol-%d.xg",out_name,front.ip);
                xg_file = fopen(xg_name,"w");
                for (i = imin; i < imax; i++)
                    fprintf(xg_file,"%f  %f\n",x[i],u_old[i]);
                fclose(xg_file);

                /* Exact solution */
		for (i = imin; i < imax; ++i)
                    exact_soln(wave_speed,front.time,x[i],u_sol[i]);
                sprintf(xg_name,"%s/exc-%d.xg",out_name,front.ip);
                xg_file = fopen(xg_name,"w");
                for (i = imin; i < imax; i++)
                    fprintf(xg_file,"%f  %f\n",x[i],u_sol[i]);
                fclose(xg_file);
            }

	    /* Termination control */
            if (FT_TimeLimitReached(&front))
            {
                front.dt = dt;
                FT_TimeControlFilter(&front); /* reduce time step for output */
                (void) printf("next dt = %20.14f\n",front.dt);
		for (i = imin; i < imax; ++i)
                    exact_soln(wave_speed,front.time,x[i],u_sol[i]);
		error_func(u_sol,u_old,err,mesh_size);
                break;
            }

            front.dt = dt;
            FT_TimeControlFilter(&front); /* reduce time step for output */
            (void) printf("next dt = %20.14f\n",front.dt);
        }
        gd_closeplot();
	free(x);
	free(u_old);
	free(u_sol);
	free(coef_old);
	free(coef_new);
	clean_up(0);
}

static void error_func(
	double *u_num,
	double *u_exc,
	double err[3],
	int mesh_size)
{
	int i;
	double tmp;

	err[0]=0.0;
	err[1]=0.0;
	err[2]=0.0;
	for (i = 0; i < mesh_size; i++)
	{
	    tmp = (u_exc[i]-u_num[i] > 0) ? u_exc[i]-u_num[i] : u_num[i]-u_exc[i];
	    err[0] += tmp;
	    err[1] += pow(tmp,2);
	    err[2] = (err[2] > tmp) ? err[2] : tmp;
	}
	err[0] /= mesh_size;
	err[1] = sqrt(err[1]/mesh_size);
	printf("\n%d\t%e\t%e\t%e\n",mesh_size,err[0],err[1],err[2]);
}

static void readDgParams(
	PARAMS *params,
	char *inname)
{
	char string[100];
        FILE *infile = fopen(inname,"r");

	CursorAfterString(infile,"Enter type of equation:");
        fscanf(infile,"%s",string);
        (void) printf(" %s\n",string);
        switch (string[0])
        {
        case 'l':
        case 'L':
	    params->eqn_type = LINEAR_EQN;
	    break;
        case 'b':
        case 'B':
	    params->eqn_type = BURGER_EQN;
	    break;
	default:
	    (void) printf("Unknown equation type!\n");
	    clean_up(ERROR);
	}
	CursorAfterString(infile,"Enter type of initial condition:");
        fscanf(infile,"%s",string);
        (void) printf(" %s\n",string);
        switch (string[0])
        {
        case 's':
        case 'S':
	    params->init_type = SQUARE;
	    break;
        case 'h':
        case 'H':
	    params->init_type = HAMP;
	    break;
        case 'c':
        case 'C':
	    params->init_type = COSINE;
	    break;
        case 'w':
        case 'W':
	    params->init_type = WAVE;
	    break;
	default:
	    (void) printf("Unknown initial condition type!\n");
	    clean_up(ERROR);
	}
	CursorAfterString(infile,"Enter coefficient a:");
        fscanf(infile,"%lf",&params->a);
        (void) printf(" %f\n",params->a);
	fclose(infile);
}	/* end readWenoParams */

static void wave_func(
        double a,
        double time,
        double x,
        double & u)
{
        double arg;

	arg = x - a*time;
	u = sin(2.0*PI*(arg - 1.0));
}

static void hamp_func(
        double a,
        double time,
        double x,
        double & u)
{
        double arg;

	arg = x - a*time;
	if (arg >= 1.0 && arg <= 2.0)
	    u = sin(PI*(arg - 1.0));
	else
	    u=0.0;
}

static void cosine_func(
        double a,
        double time,
        double x,
        double & u)
{
        double arg;

	arg = x - a*time;
	u=  0.5 + sin(2.0 * PI * arg);
}

static void square_func(
        double a,
        double time,
        double x,
        double & u)
{
        double arg;

	arg = x - a*time;

	if (arg >= 1.0 && arg <= 2.0)
	    u = 1.0;
	else
	    u = 0.0;
}

static void setSolutionType(
        INIT_TYPE init_type)
{
        switch (init_type)
        {
        case WAVE:
	    exact_soln = wave_func;
            break;
        case HAMP:
	    exact_soln = hamp_func;
            break;
	case COSINE:
	    exact_soln = cosine_func;
	    break;
        case SQUARE:
	    exact_soln = square_func;
            break;
        default:
            (void) printf("Unknown equation type!\n");
            clean_up(ERROR);
        }
}      /* end setSolutionType */

static void L2_Projection(
	double wave_speed,
        double x,
        double dx,
        double *local_coef)
{
        double ugau[5],LHS[5];
	int j,l;
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

	for (j = 0; j < 5; ++j)
	    exact_soln(wave_speed,0.0,x+xgau[j]*dx*0.5,ugau[j]);
	for (j = 0; j < 5; ++j)
	{
	    LHS[j] = wgau[1] * ugau[1] * pow(xgau[1] * 0.5, j) +
			wgau[2] * ugau[2] * pow(xgau[2] * 0.5, j) +
			wgau[3] * ugau[3] * pow(xgau[3] * 0.5, j) +
			wgau[4] * ugau[4] * pow(xgau[4] * 0.5, j);
	    if (j == 0)
		LHS[j] += wgau[0] * ugau[0] * 1.0;
	    LHS[j] *= dx * 0.5;
	}
	for (int j = 0; j < 5; ++j)
	{
	    local_coef[j] = 0.0;
	    for (int l = 0; l < 5; ++l)
		local_coef[j] += AI[j][l] * LHS[l];
	}
}	/* L2_Projection to obtain the initial coefficients*/

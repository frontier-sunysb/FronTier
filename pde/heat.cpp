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

/*
*			      wave.cpp:
*
*	Sample wave equation code for AMS-528, Numerical PDE.
*	Coded by Xiaolin Li
*
*/

#include <FronTier.h>

enum    _NUM_SCHEME {
        ERROR_SCHEME    =       -1,
        CENTRAL_EXPLICIT,
        CENTRAL_IMPLICIT,
        CRANK_NICOLSON,
        MORE_SCHEMES
};
typedef enum _NUM_SCHEME NUM_SCHEME;

struct _HEAT_PARAMS {
        double b;               	/* Diffusion coefficient */
        NUM_SCHEME num_scheme;  	/* Numerical scheme choice */
	boolean explicit_scheme;	/* Yes if the scheme is explicit */
};
typedef struct _HEAT_PARAMS HEAT_PARAMS;

char *in_name,*out_name;
/* for exact solution */
static double gauss_func(double,double,double);
static void exact_gauss_soln(double,double,double*,double*,int);

/* functions for different schemes */
static void explicit_step(double*,double*,int,int,double,int);
static void implicit_step(double*,double*,int,int,double,int);
static double central_explicit(double,double*);
static void tridiagonal_solver(int,double*,double*,double*,double*,double*);
static void init_heat_params(char*,HEAT_PARAMS*);

int main(int argc, char **argv)
{ 
	F_BASIC_DATA f_basic; /* Only some basic features of f_basic is used */
	Front front;	      /* Only some basic feature of front is used */
	double *x,*u_old,*u_new,*u_sol; /* mesh points and solutions */
	double *x_movie,*y_movie;
	double dx,dt;
	char movie_caption[100];
	char time_label[100];
	char xg_name[200];
	char gd_name[200];
	FILE *xg_file;
	double xmin,xmax,umin,umax,height;
	double (*num_solver)(double,double*);
	double R,*u_stencil;
	double L,U;
	double CFL;	/* CFL safety factor */
	int buffer_size = 3;   /* if needed, can be changed to other numbers */
	int i,imin,imax,mesh_size,expanded_mesh_size;
	int boundary_type[2];
	double b;		/* diffusion coefficient */
	HEAT_PARAMS heat_params;

	/* I/O names */
	FT_Init(argc,argv,&f_basic);
	in_name		= f_basic.in_name;
	out_name	= f_basic.out_name;
	b = 1.0;

	/* Get domain information */
	FT_ReadSpaceDomain(in_name,&f_basic);
	L = f_basic.L[0];
	U = f_basic.U[0];
	mesh_size = f_basic.gmax[0];
	boundary_type[0] = f_basic.boundary[0][0];
	boundary_type[1] = f_basic.boundary[0][1];
	imin = buffer_size;
	imax = imin + mesh_size;

	/* Time control of the run and I/O */
	FT_ReadTimeControl(in_name,&front);
	init_heat_params(in_name,&heat_params);

	CFL = Time_step_factor(&front);
	FT_ResetTime(&front);

	/* Allocating memory for arrays */
	expanded_mesh_size = mesh_size + 2*buffer_size;
	FT_VectorMemoryAlloc((POINTER*)&x,mesh_size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&u_old,expanded_mesh_size,
					sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&u_new,expanded_mesh_size,
					sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&u_sol,expanded_mesh_size,
					sizeof(double));

	/* Initialization of states */

	dx = (U - L)/(mesh_size - 1);

	front.dt = dt = CFL*0.5*dx*dx/b;
	FT_TimeControlFilter(&front);

	/* Set the initial condition */
	for (i = 0; i < expanded_mesh_size; i++)
	    x[i] = L + (i - imin)*dx;
	exact_gauss_soln(b,0.0,x,u_old,expanded_mesh_size);

	/* Assign numerical scheme */
	num_solver = central_explicit;

	/* Set frame margin for GD movie output */
	xmin = L;	xmax = U;
	umin = HUGE;	umax = -HUGE;
	for (i = 0; i < mesh_size; i++)
	{
	    if (umin > u_old[i]) umin = u_old[i];
	    if (umax < u_old[i]) umax = u_old[i];
	}
	height = umax - umin;
	umin -= 0.15*height;	umax += 0.15*height;

	/* Open and initialize GD movie output */
	sprintf(movie_caption,"u vs. x");
	sprintf(gd_name,"%s/soln.gif",out_name);
	gd_initplot(gd_name,movie_caption,xmin,xmax,umin,umax,2);

	/* Time loop */
	for (;;)
	{
	    R = b*front.dt/dx/dx;
	    /* Advancing numerical solution */
	    if (heat_params.explicit_scheme == YES)
	    	explicit_step(u_old,u_new,imin,imax,R,heat_params.num_scheme);
	    else
	    	implicit_step(u_old,u_new,imin,imax,R,heat_params.num_scheme);

	    /* Swapping solution storage */
	    for (i = imin; i < imax; i++)
	    {
		u_old[i] = u_new[i];
	    }

	    /* Time and step control */
	    FT_AddTimeStepToCounter(&front);
	    (void) printf("\ntime = %20.14f   step = %5d   ",
                        	front.time,front.step);

	    /* Update boundary condition */
	    for (i = 0; i < buffer_size; ++i)
	    {
	    	u_old[i] = gauss_func(b,front.time,x[i]);
	    	u_old[imax+i] = gauss_func(b,front.time,x[imax+i]);
	    }

	    /* Movie frame */
	    if (FT_IsMovieFrameTime(&front))
	    {
		/* Numerical solution */
		x_movie = x+buffer_size;
		y_movie = u_old+buffer_size;
		gd_plotdata(mesh_size,x_movie,y_movie);

		/* Exact solution */
		exact_gauss_soln(b,front.time,x,u_sol,expanded_mesh_size);
		x_movie = x+buffer_size;
		y_movie = u_sol+buffer_size;
		gd_plotdata(mesh_size,x_movie,y_movie);

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
		fprintf(xg_file,"\"u vs. x\"\n");
		for (i = imin; i < imax; ++i)
		{
		    fprintf(xg_file,"%f  %f\n",x[i],u_old[i]);
		}
		fclose(xg_file);

		/* Exact solution */
		//wave_func(a,front.time,x,u_sol,expanded_mesh_size);
		sprintf(xg_name,"%s/exc-%d.xg",out_name,front.ip);
		xg_file = fopen(xg_name,"w");
		fprintf(xg_file,"\"u vs. x\"\n");
		for (i = imin; i < imax; ++i)
		{
		    fprintf(xg_file,"%f  %f\n",x[i],u_sol[i]);
		}
		fclose(xg_file);
	    }

	    /* Termination control */
	    if (FT_TimeLimitReached(&front))
	    {
	        front.dt = dt;
	        FT_TimeControlFilter(&front); /* reduce time step for output */
	        (void) printf("next dt = %20.14f\n",front.dt);
                break;
	    }

	    front.dt = dt;
	    FT_TimeControlFilter(&front); /* reduce time step for output */
	    (void) printf("next dt = %20.14f\n",front.dt);
	}
	gd_closeplot();
}

static void exact_gauss_soln(
        double mu,
        double time,
        double *x,
        double *u,
        int mesh_size)
{
	int i;
        double arg;

        for (i = 0; i < mesh_size; ++i)
        {
	    u[i] = gauss_func(mu,time,x[i]);
	}
}

static double gauss_func(
	double mu,
	double time,
	double x)
{
	double arg;
	double u;

	u = exp(-x*x/(4.0*mu*(time+1.0)))/sqrt(4.0*PI*mu*(time+1.0));
	return u;
}

static double central_explicit(
	double R,
	double *u)
{
	double u_new;
	u_new = u[0] + R*(u[1] - 2.0*u[0] + u[-1]);
	return u_new;
}

static void tridiagonal_solver(
	int n,
	double *a,
	double *b,
	double *c,
	double *d,
	double *u)
{
	static double *A,*B;
	int i;
	double denom;

	if (A == NULL)
	{
	    FT_VectorMemoryAlloc((POINTER*)&A,n-1,sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&B,n-1,sizeof(double));
	}

	A[n-2] = -a[n-1]/b[n-1];
	B[n-2] =  d[n-1]/b[n-1];

	for (i = n-2; i > 0; --i)
	{
	    denom = A[i]*c[i] + b[i];
	    A[i-1] = -a[i]/denom;
	    B[i-1] = (d[i] - B[i]*c[i])/denom;
	}

	u[0] = (d[0] - B[0]*c[0])/(A[0]*c[0] + b[0]);
	for (i = 0; i < n-1; ++i)
	{
	    u[i+1] = A[i]*u[i] + B[i];
	}
}
	
static void init_heat_params(
	char *in_name,
	HEAT_PARAMS *heat_params)
{
	FILE *infile = fopen(in_name,"r");
	char string[100];

	CursorAfterString(infile,"Enter diffusion coefficient:");
	fscanf(infile,"%lf",&heat_params->b);
	(void) printf(" %f\n",heat_params->b);

	heat_params->explicit_scheme = YES;
	CursorAfterString(infile,"Enter numerical scheme:");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	if (string[0] == 'c' || string[0] == 'C')
	{
	    if (string[8] == 'e' || string[8] == 'E')
	    	heat_params->num_scheme = CENTRAL_EXPLICIT;
	    else if (string[8] == 'i' || string[8] == 'I')
	    {
	    	heat_params->num_scheme = CENTRAL_IMPLICIT;
		heat_params->explicit_scheme = NO;
	    }
	    else if (string[1] == 'r' || string[1] == 'R')
	    {
	    	heat_params->num_scheme = CRANK_NICOLSON;
		heat_params->explicit_scheme = NO;
	    }
	}
	else
	{
	    (void) printf("Unknown numerical scheme!\n");
	    clean_up(ERROR);
	}
}	/* end init_heat_params */

static void explicit_step(
	double *u_old,
	double *u_new,
	int imin,
	int imax,
	double R,
	int num_scheme)
{
	int i;
	double *u_stencil;
	double (*num_solver)(double,double*);

	/* Assign numerical scheme */
	switch (num_scheme)
	{
	case CENTRAL_EXPLICIT:
	    num_solver = central_explicit;
	    break;		/* more schemes to be added below */
	default:
	    (void) printf("Numerical scheme unimplemented!\n");
	    clean_up(ERROR);
	}

	for (i = imin; i < imax; i++)
	{
	    u_stencil = u_old + i;
	    u_new[i] = num_solver(R,u_stencil);
	}
}

static void implicit_step(
	double *u_old,
	double *u_new,
	int imin,
	int imax,
	double R,
	int num_scheme)
{
	int i,n = imax - imin;
	static double *a,*b,*c,*d;
	double *pu_old,*pu_new;

	if (a == NULL) /* Allocate memory for coefficients */
	{
	    FT_VectorMemoryAlloc((POINTER*)&a,n,sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&b,n,sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&c,n,sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&d,n,sizeof(double));
	}

	pu_old = u_old+imin;
	pu_new = u_new+imin;
	for (i = 0; i < n; ++i)
	{
	    switch (num_scheme)
	    {
	    case CENTRAL_IMPLICIT:
	    	a[i] = -R;
		c[i] = -R;
		b[i] = 1.0 + 2.0*R;
		d[i] = pu_old[i];
		if (i == 0)
		{
		    a[i] = 0.0;
		    d[i] += R*pu_old[i-1];
		}
		else if (i == n-1)
		{
		    c[i] = 0.0;
		    d[i] += R*pu_old[i+1];
		}
	    	break;
	    case CRANK_NICOLSON:
	    	a[i] = -0.5*R;
		c[i] = -0.5*R;
		b[i] =  1.0 + R;
		d[i] = pu_old[i] + 0.5*R*
			(pu_old[i+1] - 2.0*pu_old[i] + pu_old[i-1]);
		if (i == 0)
		{
		    a[i] = 0.0;
		    d[i] += 0.5*R*pu_old[i-1];
		}
		else if (i == n-1)
		{
		    c[i] = 0.0;
		    d[i] += 0.5*R*pu_old[i+1];
		}
	    	break;
	    }
	}
	tridiagonal_solver(n,a,b,c,d,pu_new);
}

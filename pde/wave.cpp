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
*			      wave.cpp:
*
*	Sample wave equation code for AMS-528, Numerical PDE.
*	Coded by Xiaolin Li
*
*/

#include <FronTier.h>

enum	_NUM_SCHEME {
	ERROR_SCHEME	=	-1,
	UPWIND_SCHEME	=	1,
	LAX_FRIEDRICH,
	LAX_WENDROFF,
	MAC_CORMACK,
	CENTRAL_EXPLICIT,
	CENTRAL_IMPLICIT,
	CRANK_NICOLSON,
	MORE_SCHEMES
};
typedef enum _NUM_SCHEME NUM_SCHEME;

struct _WAVE_PARAMS {
	double a;		/* Wave speed */
	NUM_SCHEME num_scheme;	/* Numerical scheme choice */
	boolean explicit_scheme;
};
typedef struct _WAVE_PARAMS WAVE_PARAMS;

char *in_name,*out_name;
static void wave_func(double,double,double*,double*,int);
static void init_wave_params(char*,WAVE_PARAMS*);
static void explicit_step(double*,double*,int,int,double,int);
static void implicit_step(double*,double*,int,int,double,int);
static void tridiagonal_solver(int,double*,double*,double*,double*,double*);

/* functions for different schemes */
static double upwind(double,double*);
static double lax_wendroff(double,double*);
static double lax_friedrich(double,double*);
static double mac_cormack(double,double*);
static double central_explicit(double,double*);

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
	double R;
	double L,U;
	double CFL;	/* CFL safety factor */
	int buffer_size = 3;   /* if needed, can be changed to other numbers */
	int i,imin,imax,mesh_size,expanded_mesh_size;
	int boundary_type[2];
	WAVE_PARAMS wave_params;
	double a;		/* wave speed */
	NUM_SCHEME num_scheme;

	/* I/O names */
	FT_Init(argc,argv,&f_basic);
	in_name		= f_basic.in_name;
	out_name	= f_basic.out_name;
	a = 1.0;

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
	CFL = Time_step_factor(&front);
	FT_ResetTime(&front);

	/* Init probblem specific parameters */
	init_wave_params(in_name,&wave_params);
	a = wave_params.a;
	num_scheme = wave_params.num_scheme;

	/* Allocating memory for arrays */
	expanded_mesh_size = mesh_size + 2*buffer_size;
	FT_VectorMemoryAlloc((POINTER*)&x,expanded_mesh_size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&u_old,expanded_mesh_size,
					sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&u_new,expanded_mesh_size,
					sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&u_sol,expanded_mesh_size,
					sizeof(double));

	/* Initialization of states */

	dx = (U - L)/(mesh_size - 1);

	front.dt = dt = CFL*dx/a;
	FT_TimeControlFilter(&front);

	/* Set the initial condition */
	for (i = 0; i < expanded_mesh_size; i++)
	    x[i] = L + (i - imin)*dx;
	wave_func(a,0.0,x,u_old,expanded_mesh_size);

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
	    R = a*front.dt/dx;

	    /* Advancing numerical solution */
	    if (wave_params.explicit_scheme == YES)
	    	explicit_step(u_old,u_new,imin,imax,R,wave_params.num_scheme);
	    else
	    	implicit_step(u_old,u_new,imin,imax,R,wave_params.num_scheme);

	    /* Swapping solution storage */
	    for (i = imin; i < imax; i++)
	    {
		u_old[i] = u_new[i];
	    }

	    /* Time and step control */
	    FT_AddTimeStepToCounter(&front);
	    FT_PrintTimeStamp(&front);

	    /* Movie frame */
	    if (FT_IsDrawTime(&front))
	    {
		/* Numerical solution */
		x_movie = x+buffer_size;
		y_movie = u_old+buffer_size;
		gd_plotdata(mesh_size,x_movie,y_movie);

		/* Exact solution */
		wave_func(a,front.time,x,u_sol,expanded_mesh_size);
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
		wave_func(a,front.time,x,u_sol,expanded_mesh_size);
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

static void wave_func(
	double a,
	double time,
	double *x,
	double *u,
	int mesh_size)
{
	int i;
	double arg;

	for (i = 0; i < mesh_size; ++i)
	{
	    arg = x[i] - a*time;
	    if (arg >= 1.0 && arg <= 2.0)
		u[i] = sin(PI*(arg - 1.0))*sin(PI*(arg - 1.0));
	    else
		u[i] = 0.0;
	}
}

static double upwind(
	double R,
	double *u)
{
	double u_new;

	u_new = (R > 0.0) ? 
		u[0] - R*(u[0] - u[-1]) :
		u[0] - R*(u[1] - u[0]);
	return u_new;
}

static double lax_friedrich(
	double R,
	double *u)
{
	double u_new;
	u_new = 0.5*(u[1] + u[-1]) - 0.5*R*(u[1] - u[-1]);
	return u_new;
	/* Scheme to be implemented */
}

static double lax_wendroff(
	double R,
	double *u)
{
	double u_new;
	u_new = u[0] - 0.5*R*(u[1] - u[-1]) + 0.5*R*R*(u[1] - 2.0*u[0] + u[-1]);
	return u_new;
	/* Scheme to be implemented */
}

static double central_explicit(
	double R,
	double *u)
{
	double u_new;
	u_new = u[0] - 0.5*R*(u[1] - u[-1]);
	return u_new;
	/* Scheme to be implemented */
}

static double mac_cormack(
	double R,
	double *u)
{
	double u_tmp[2],u_new;
	u_tmp[0] = u[-1] - R*(u[0] - u[-1]);
	u_tmp[1] = u[0] - R*(u[1] - u[0]);
	u_new = 0.5*(u[0] + u_tmp[1]) - 0.5*R*(u_tmp[1] - u_tmp[0]);
	return u_new;
}

static void init_wave_params(
	char *in_name,
	WAVE_PARAMS *wave_params)
{
	FILE *infile = fopen(in_name,"r");
	char string[100];

	CursorAfterString(infile,"Enter wave speed:");
	fscanf(infile,"%lf",&wave_params->a);
	(void) printf(" %f\n",wave_params->a);

	wave_params->explicit_scheme = YES;
	CursorAfterString(infile,"Enter numerical scheme:");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	if (string[0] == 'u' || string[0] == 'U')
	    wave_params->num_scheme = UPWIND_SCHEME;
	else if (string[0] == 'm' || string[0] == 'M')
	    wave_params->num_scheme = MAC_CORMACK;
	else if (string[0] == 'l' || string[0] == 'L')
	{
	    if (string[4] == 'w' || string[4] == 'W')
	    	wave_params->num_scheme = LAX_WENDROFF;
	    else if (string[4] == 'f' || string[4] == 'F')
	    	wave_params->num_scheme = LAX_FRIEDRICH;
	}
	else if (string[0] == 'c' || string[0] == 'C')
	{
	    if (string[8] == 'e' || string[8] == 'E')
	    	wave_params->num_scheme = CENTRAL_EXPLICIT;
	    else if (string[8] == 'i' || string[8] == 'I')
	    {
	    	wave_params->num_scheme = CENTRAL_IMPLICIT;
		wave_params->explicit_scheme = NO;
	    }
	    else if (string[1] == 'r' || string[1] == 'R')
	    {
	    	wave_params->num_scheme = CRANK_NICOLSON;
		wave_params->explicit_scheme = NO;
	    }
	}
	else
	{
	    (void) printf("Unknown numerical scheme!\n");
	    clean_up(ERROR);
	}
}	/* end init_wave_params */

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
	case UPWIND_SCHEME:
	    num_solver = upwind;
	    break;
	case LAX_WENDROFF:
	    num_solver = lax_wendroff;
	    break;
	case LAX_FRIEDRICH:
	    num_solver = lax_friedrich;
	    break;
	case MAC_CORMACK:
	    num_solver = mac_cormack;
	    break;
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
	    	a[i] = -0.5*R;
		c[i] =  0.5*R;
		b[i] =  1.0;
		d[i] = pu_old[i];
		if (i == 0)
		{
		    a[i] = 0.0;
		    d[i] += 0.5*R*pu_old[i-1];
		}
		else if (i == n-1)
		{
		    c[i] = 0.0;
		    d[i] -= 0.5*R*pu_old[i+1];
		}
	    	break;
	    case CRANK_NICOLSON:
	    	a[i] = -0.25*R;
		c[i] =  0.25*R;
		b[i] =  1.0;
		d[i] = pu_old[i] - 0.25*R*(pu_old[i+1] - pu_old[i-1]);
		if (i == 0)
		{
		    a[i] = 0.0;
		    d[i] += 0.25*R*pu_old[i-1];
		}
		else if (i == n-1)
		{
		    c[i] = 0.0;
		    d[i] -= 0.25*R*pu_old[i+1];
		}
	    	break;
	    }
	}
	tridiagonal_solver(n,a,b,c,d,pu_new);
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
	

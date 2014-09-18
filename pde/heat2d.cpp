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
        EXPLICIT,
	IMPLICIT_SPLIT,
	IMPLICIT_ADI,
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
static void explicit_step(double**,double**,int*,int*,double,double,int);
static void implicit_step(double**,double**,int*,int*,double,double,int);
static double central_explicit(double,double*);
static void tridiagonal_solver(int,double*,double*,double*,double*,double*);
static void init_heat_params(char*,HEAT_PARAMS*);

int main(int argc, char **argv)
{ 
	F_BASIC_DATA f_basic; /* Only some basic features of f_basic is used */
	Front front;	      /* Only some basic feature of front is used */
	double *x,*y,**u_old,**u_new,**u_sol; /* mesh points and solutions */
	double dx,dy,dh,dt;
	double (*num_solver)(double,double*);
	double Rx,Ry,*u_stencil;
	double *L,*U;
	int *mesh_size,expanded_mesh_size[2];
	double CFL;	/* CFL safety factor */
	int buffer_size = 3;   /* if needed, can be changed to other numbers */
	int i,j,imin[2],imax[2];
	int boundary_type[2][2];
	double b;		/* diffusion coefficient */
	HEAT_PARAMS heat_params;
        NUM_SCHEME num_scheme;  	/* Numerical scheme choice */

	/* I/O names */
	FT_Init(argc,argv,&f_basic);
	in_name		= f_basic.in_name;
	out_name	= f_basic.out_name;
	b = 1.0;

	/* Get domain information */
	FT_ReadSpaceDomain(in_name,&f_basic);
	L = f_basic.L;
	U = f_basic.U;
	mesh_size = f_basic.gmax;

	/* Time control of the run and I/O */
	FT_ReadTimeControl(in_name,&front);
	init_heat_params(in_name,&heat_params);
	num_scheme = heat_params.num_scheme;

	CFL = Time_step_factor(&front);
	FT_ResetTime(&front);

	/* Allocating memory for arrays */
	for (i = 0; i < 2; ++i)
	{
	    imin[i] = buffer_size;
	    imax[i] = imin[i] + mesh_size[i];
	    expanded_mesh_size[i] = mesh_size[i] + 2*buffer_size;
	    boundary_type[i][0] = f_basic.boundary[i][0];
	    boundary_type[i][1] = f_basic.boundary[i][1];
	}

	FT_VectorMemoryAlloc((POINTER*)&x,expanded_mesh_size[0],sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&y,expanded_mesh_size[1],sizeof(double));
	FT_MatrixMemoryAlloc((POINTER*)&u_old,expanded_mesh_size[0],
				expanded_mesh_size[1],sizeof(double));
	FT_MatrixMemoryAlloc((POINTER*)&u_new,expanded_mesh_size[0],
				expanded_mesh_size[1],sizeof(double));
	FT_MatrixMemoryAlloc((POINTER*)&u_sol,expanded_mesh_size[0],
				expanded_mesh_size[1],sizeof(double));

	/* Initialization of states */

	dx = (U[0] - L[0])/(mesh_size[0] - 1);
	dy = (U[1] - L[1])/(mesh_size[1] - 1);
	dh = std::min(dx,dy);

	front.dt = dt = CFL*0.5*sqr(dh)/b;
	FT_TimeControlFilter(&front);

	/* Set the initial condition */
	for (i = 0; i < expanded_mesh_size[0]; i++)
	    x[i] = L[0] + (i - imin[0])*dx;

	/* Assign numerical scheme */

	/* Time loop */
	for (;;)
	{
	    Rx = b*front.dt/sqr(dx);
	    Ry = b*front.dt/sqr(dy);
	    /* Advancing numerical solution */
	    if (heat_params.explicit_scheme == YES)
		explicit_step(u_old,u_new,imin,imax,Rx,Ry,num_scheme);
	    else
		implicit_step(u_old,u_new,imin,imax,Rx,Ry,num_scheme);

	    /* Swapping solution storage */
	    for (i = imin[0]; i < imax[0]; i++)
	    for (j = imin[1]; j < imax[1]; j++)
	    {
		;
	    }

	    /* Time and step control */
	    FT_AddTimeStepToCounter(&front);
	    FT_PrintTimeStamp(&front);

	    /* Update boundary condition */
	    for (i = 0; i < buffer_size; ++i)
	    {
		;
	    }

	    /* Movie frame */
	    if (FT_IsMovieFrameTime(&front))
	    {
	    }
	    /* Output date control */
	    if (FT_IsSaveTime(&front))
	    {
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
	if (string[0] == 'e' || string[0] == 'E')
	{
	    heat_params->num_scheme = EXPLICIT;
	}
	else if (string[0] == 'i' || string[0] == 'I')
	{
	    if (string[9] == 's' || string[9] == 'S')
	    	heat_params->num_scheme = IMPLICIT_SPLIT;
	    else if (string[9] == 'a' || string[9] == 'A')
	    	heat_params->num_scheme = IMPLICIT_ADI;
	}
	else
	{
	    (void) printf("Unknown numerical scheme!\n");
	    clean_up(ERROR);
	}
}	/* end init_heat_params */

static void explicit_step(
	double **u_old,
	double **u_new,
	int *imin,
	int *imax,
	double Rx,
	double Ry,
	int num_scheme)
{
	int i,j;


	for (i = imin[0]; i < imax[0]; i++)
	for (j = imin[1]; j < imax[1]; j++)
	{
		;
	}
}

static void implicit_step(
	double **u_old,
	double **u_new,
	int *imin,
	int *imax,
	double Rx,
	double Ry,
	int num_scheme)
{
	int i,n;
	static double *a,*b,*c,*d,*pu_new;

	if (a == NULL) /* Allocate memory for coefficients */
	{
	    n = 0;
	    for (i = 0; i < 2; ++i)
		if (n < imax[i] - imin[i]) n = imax[i] - imin[i];
	    FT_VectorMemoryAlloc((POINTER*)&a,n,sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&b,n,sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&c,n,sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&d,n,sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&pu_new,n,sizeof(double));
	}

	for (i = 0; i < n; ++i)
	{
	    switch (num_scheme)
	    {
	    case IMPLICIT_SPLIT:
		tridiagonal_solver(n,a,b,c,d,pu_new);
	    	break;
	    case IMPLICIT_ADI:
		tridiagonal_solver(n,a,b,c,d,pu_new);
	    	break;
	    }
	}
}

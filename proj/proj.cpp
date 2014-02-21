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
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*/

#include "proj.h"
#define MAXDIM 101



char *in_name,*out_name;
boolean binary = NO;

typedef struct{
        double *vals;
        int m; // Number of elements
} Vector;

static void rand_vector(Vector, unsigned short int*);
static void xprint_vector(const char*,const char*,double*,double*,int);
static void solve_for_p3(Vector,double*,double);
static void solve_for_p5(Vector,double*,double);

static void rand_vector(
	Vector vec, 
	unsigned short int *seed)
{
        int i;

        for(i = 0; i < vec.m; i++)
        {
            vec.vals[i] = 2.0*erand48(seed) - 1.0;
        }
}	/* end rand_vector */

static void xprint_vector(
	const char *color,
	const char *fname,
	double *x,
	double *y,
	int N)
{
        int i;
	FILE *ofile = fopen(fname,"w");

	fprintf(ofile,"Next\n");
	fprintf(ofile,"!%s vs x\n",fname);
	fprintf(ofile,"color=%s\n",color);
	fprintf(ofile,"thickness=1.5\n");
        for(i = 0; i < N; i++)
        {
            fprintf(ofile,"%f %f\n",x[i],y[i]);
        }
	fprintf(ofile,"\n");
}	/* end xprint_vector */

int main(int argc, char **argv)
{
	static Front front;
	static F_BASIC_DATA f_basic;
	unsigned short int seed[3];
        Vector u;
	double *x;

	/* Initialize basic computational data */

	FT_Init(argc,argv,&f_basic);//Read parameters from command line
	f_basic.dim = 1;
	f_basic.size_of_intfc_state = 0;
	
	//Initialize Petsc before the FT_StartUp
	PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);

        in_name                 = f_basic.in_name;
        out_name                = f_basic.out_name;

	//FT_ReadSpaceDomain(in_name,&f_basic);
	//FT_InitDebug(in_name);

        double h = 1.0/(MAXDIM-1);
        int i,j;

	u.m = MAXDIM;
	FT_VectorMemoryAlloc((POINTER*)&u.vals,MAXDIM,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&x,MAXDIM,sizeof(double));
	
	// u will be a random vector
	seed[0] = 12;
	seed[1] = 36;
	seed[2] = 64;
        rand_vector(u,seed);
	u.vals[0] = u.vals[1];
	for (i = 0; i < MAXDIM; ++i)
	    x[i] = i*h;

        xprint_vector("blue","u0.xg",x,u.vals,u.m);

	solve_for_p3(u,x,h);
	solve_for_p5(u,x,h);

	PetscFinalize();
	clean_up(0);
}

static void solve_for_p3(
        Vector u,
	double *x,
	double h)
{
	int N = u.m;
        int ilower = 0;
        int iupper = N-2;
        PETSc solver;
        PetscInt num_iter = 0;
        double rel_residual = 0.0;
	int i;
	static double *soln,*b;
	static int *i_to_I,I;
	static double *p,*du,*w;

	if (soln == NULL)
	{
	    FT_VectorMemoryAlloc((POINTER*)&du,N,sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&w,N,sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&p,N,sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&soln,iupper,sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&b,N,sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&i_to_I,N,sizeof(int));
	}

	for (i = 2; i < N-1; ++i)
	{
	    du[i] = (u.vals[i+1] - u.vals[i-1])/(2.0*h);
	}
	du[1] = (u.vals[2] - u.vals[1])/(2.0*h);
	du[0] = du[1];
	du[N-1] = du[N-2];

        xprint_vector("red","du-3.xg",x,du,N);

	// I have added here.
	solver.Create(ilower, iupper-1, 3, 3);
        solver.Reset_A();
        solver.Reset_b();
        solver.Reset_x();

	i_to_I[0] = i_to_I[N-1] = -1;
	for (i = 1; i < N-1; ++i)
	{
	    b[i] = sqr(h)*du[i];
	    i_to_I[i] = i - 1;
	}
	b[1] += -1.0;

        for(i = 1; i < N-1; i++)
        {
	    I = i_to_I[i];
	    if (i != 1)
            	solver.Set_A(I,I-1,1.0);
	    if (i != N-2)
	    {
                solver.Set_A(I,I,-2.0);
            	solver.Set_A(I,I+1,1.0);
	    }
	    else
            	solver.Set_A(I,I,-1.0);
            solver.Set_b(I,b[i]);
	    
	}
        solver.SetMaxIter(40000);
        solver.SetTol(1e-10);

        solver.Solve();
	solver.GetNumIterations(&num_iter);
	solver.GetFinalRelativeResidualNorm(&rel_residual);
	solver.Get_x(soln);

        for(i = 1; i < N-1; i++)
	{
	    I = i_to_I[i];
	    p[i] = soln[I];
	}
	p[0] = 1.0;
	p[N-1] = p[N-2];
        for(i = 1; i < N-1; i++)
	{
	    printf("LHS = %f  RHS = %f  error = %g\n",p[i+1]-2.0*p[i]+p[i-1],
			sqr(h)*du[i],p[i+1]-2.0*p[i]+p[i-1]-sqr(h)*du[i]);
	}

        xprint_vector("green","p-3.xg",x,p,N);
 
	w[0] = u.vals[0];
        for(i = 1; i < N-2; i++)
	{
	    w[i] = u.vals[i] - (p[i+1] - p[i-1])/(2.0*h);
	}
	w[N-2] = u.vals[N-2] - 
			(p[N-2] - p[N-3])/(2.0*h);
	w[N-1] = u.vals[N-1];
	w[0] = w[1];

        xprint_vector("orange","w-3.xg",x,w,N);
	du[1] = (w[2] - w[1])/(2.0*h);
	for (i = 2; i < N-1; ++i)
	{
	    du[i] = (w[i+1] - w[i-1])/(2.0*h);
	}
	du[0] = du[1];
	du[N-1] = du[N-2];
        xprint_vector("orange","du1-3.xg",x,du,N);

}	/* end solve_for_p3 */

static void solve_for_p5(
        Vector u,
	double *x,
	double h)
{
	int N = u.m;
        int ilower = 0;
        int iupper = N-2;
        PETSc solver;
        PetscInt num_iter = 0;
        double rel_residual = 0.0;
	int i;
	static double *soln,*b;
	static int *i_to_I,I;
	static double *p,*du,*w;

	if (soln == NULL)
	{
	    FT_VectorMemoryAlloc((POINTER*)&du,N,sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&w,N,sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&p,N,sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&soln,iupper,sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&b,N,sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&i_to_I,N,sizeof(int));
	}

	for (i = 2; i < N-1; ++i)
	{
	    du[i] = (u.vals[i+1] - u.vals[i-1])/(2.0*h);
	}
	du[1] = (u.vals[2] - u.vals[1])/(2.0*h);
	du[0] = du[1];
	du[N-1] = du[N-2];

        xprint_vector("red","du-3.xg",x,du,N);

	// I have added here.
	solver.Create(ilower, iupper-1, 5, 5);
        solver.Reset_A();
        solver.Reset_b();
        solver.Reset_x();

	i_to_I[0] = i_to_I[N-1] = -1;
	for (i = 1; i < N-1; ++i)
	    i_to_I[i] = i - 1;

	for (i = 1; i < N-1; ++i)
	{
	    b[i] = 4.0*sqr(h)*du[i];
	}

	b[1] += -1.0;
	I = i_to_I[1];
        solver.Set_A(I,I,-1.0);
        solver.Set_A(I,I+1,-1.0);
        solver.Set_A(I,I+2,1.0);

	b[2] += -1.0;
	I = i_to_I[2];
        solver.Set_A(I,I,-2.0);
        solver.Set_A(I,I+2,1.0);

        for(i = 3; i < N-3; i++)
        {
	    I = i_to_I[i];
            solver.Set_A(I,I,-2.0);
            solver.Set_A(I,I-2,1.0);
            solver.Set_A(I,I+2,1.0);
	}
	I = i_to_I[N-3];
        solver.Set_A(I,I+1,1.0);
        solver.Set_A(I,I,-2.0);
        solver.Set_A(I,I-2,1.0);
	I = i_to_I[N-2];
        solver.Set_A(I,I,-1.0);
        solver.Set_A(I,I-2,1.0);

        solver.SetMaxIter(40000);
        solver.SetTol(1e-10);
        for(i = 1; i < N-1; i++)
	{
	    I = i_to_I[i];
            solver.Set_b(I,b[i]);
	}

        solver.Solve();
	solver.GetNumIterations(&num_iter);
	solver.GetFinalRelativeResidualNorm(&rel_residual);
	solver.Get_x(soln);

        for(i = 1; i < N-1; i++)
	{
	    I = i_to_I[i];
	    p[i] = soln[I];
	}
	p[0] = 1.0;
	p[N-1] = p[N-2];
        for(i = 3; i < N-3; i++)
	{
	    printf("LHS = %f  RHS = %f  error = %g\n",p[i+2]-2.0*p[i]+p[i-2],
			4*sqr(h)*du[i],p[i+2]-2.0*p[i]+p[i-2]-4*sqr(h)*du[i]);
	}

        xprint_vector("green","p-5.xg",x,p,N);
 
	w[0] = u.vals[0];
        for(i = 1; i < N-2; i++)
	{
	    w[i] = u.vals[i] - (p[i+1] - p[i-1])/(2.0*h);
	}
	w[N-2] = u.vals[N-2] - 
			(p[N-2] - p[N-3])/(2.0*h);
	w[N-1] = u.vals[N-1];
	w[0] = w[1];

        xprint_vector("orange","w-5.xg",x,w,N);
	du[1] = (w[2] - w[1])/(2.0*h);
	for (i = 2; i < N-1; ++i)
	{
	    du[i] = (w[i+1] - w[i-1])/(2.0*h);
	}
	du[0] = du[1];
	du[N-1] = du[N-2];
        xprint_vector("orange","du1-5.xg",x,du,N);

}	/* end solve_for_p5 */

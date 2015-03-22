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

static void rand_vector(double*,int,unsigned short int*);
static void xprint_vector(const char*,const char*,double*,double*,int);
static void solve_for_dual(double*,int,double*,double);
static void divergence(double*,double*,double,int);

static void rand_vector(
	double *u,
	int N,
	unsigned short int *seed)
{
        int i;

        for(i = 0; i < N; i++)
        {
            u[i] = 2.0*erand48(seed) - 1.0;
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
        double *u;
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

	FT_VectorMemoryAlloc((POINTER*)&u,MAXDIM,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&x,MAXDIM,sizeof(double));
	
	// u will be a random vector
	seed[0] = 12;
	seed[1] = 36;
	seed[2] = 64;
        rand_vector(u,MAXDIM,seed);

	for (i = 0; i < MAXDIM; ++i)
	    x[i] = i*h;

        xprint_vector("blue","u.xg",x,u,MAXDIM);

	solve_for_dual(u,MAXDIM,x,h);

	PetscFinalize();
	clean_up(0);
}

static void solve_for_dual(
        double *u,
	int N,
	double *x,
	double h)
{
        PETSc solver;
        PetscInt num_iter = 0;
        double rel_residual = 0.0;
	int i;
	static double *soln,*b;
	static double *p,*du,*w;

	if (soln == NULL)
	{
	    FT_VectorMemoryAlloc((POINTER*)&du,N+1,sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&w,N,sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&p,N+1,sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&soln,N-1,sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&b,N-1,sizeof(double));
	}

	divergence(u,du,h,N);
        xprint_vector("red","du-dual.xg",x,du,N-1);

	solver.Create(0, N-2, 3, 3);
        solver.Reset_A();
        solver.Reset_b();
        solver.Reset_x();

	for (i = 0; i < N-1; ++i)
	    b[i] = sqr(h)*du[i+1];
	b[0] += -1.0;
	for (i = 0; i < N-1; ++i)
	{
	    solver.Set_b(i,b[i]);
	    solver.Set_A(i,i,-2.0);
	    if (i!=0)
		solver.Set_A(i,i-1,1.0);
	    if (i!=N-2)
		solver.Set_A(i,i+1,1.0);
	}
	solver.Set_A(N-2,N-2,-1.0);
	    
        solver.SetMaxIter(40000);
        solver.SetTol(1e-10);
        solver.Solve();
	solver.GetNumIterations(&num_iter);
	solver.GetFinalRelativeResidualNorm(&rel_residual);
	solver.Get_x(soln);

        for(i = 1; i < N; i++)
	    p[i] = soln[i-1];
	p[0] = 1.0;
	p[N] = p[N-1];
        for(i = 1; i < N; i++)
	{
	    printf("LHS = %f  RHS = %f  error = %g\n",p[i+1]-2.0*p[i]+p[i-1],
			sqr(h)*du[i-1],p[i+1]-2.0*p[i]+p[i-1]-sqr(h)*du[i-1]);
	}
        xprint_vector("green","pr-dual.xg",x,p,N);
 
        for(i = 1; i < N-1; i++)
	    w[i] = u[i] - (p[i+1] - p[i])/h;
	w[0] = w[1];
	w[N-1] = w[N-2];
        xprint_vector("orange","w.xg",x,w,N);

	divergence(w,du,h,N);
        xprint_vector("orange","dw-dual.xg",x,du,N-1);

}	/* end solve_for_dual */

static void divergence(
	double *u,
	double *du,
	double h,
	int N)
{
	int i;
	for (i = 1; i < N; ++i)
	    du[i] = (u[i] - u[i-1])/h;	
	du[N-1] = du[N-2];
}	/* end divergence */

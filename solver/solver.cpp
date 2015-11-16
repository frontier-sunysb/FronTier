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

/*		PETSc.c
 *  Only for one node.
 *      This class PETSc is created to be used a handy interface
 *  to the function calls to PETSc. For each algebric equation 
 *  Ax=b, one PETSc instance is needed, since this instance has 
 *  the storage of these three variables. 
*/ 
#include "solver.h"
#include "petscmat.h"

PETSc::PETSc()
{
	x = NULL;			/* approx solution, RHS*/
	b = NULL;
  	A = NULL;            		/* linear system matrix */
  	
  	ksp = NULL;        		/* Krylov subspace method context */
	nullsp = NULL;
	pc = NULL;
	
	KSPCreate(PETSC_COMM_WORLD,&ksp);
}

PETSc::PETSc(int ilower, int iupper, int d_nz, int o_nz)
{	
	x = NULL;      			/* approx solution, RHS*/
	b = NULL;
  	A = NULL;            		/* linear system matrix */
  	
  	ksp = NULL;          		/* Krylov subspace method context */
	nullsp = NULL;
	pc = NULL;
	Create(ilower, iupper, d_nz, o_nz);	
	KSPCreate(PETSC_COMM_WORLD,&ksp);
}

void PETSc::Create(int ilower, int iupper, int d_nz, int o_nz)
{	
	Create(PETSC_COMM_WORLD, ilower, iupper, d_nz, o_nz);	
}

void PETSc::Create(
	MPI_Comm Comm, 
	int ilower, 
	int iupper, 
	int d_nz, 
	int o_nz)
{	
	int n	= iupper - ilower +1;
	
	comm 	= Comm;
	iLower	= ilower;	
	iUpper 	= iupper;	
	
	MatCreateAIJ(Comm,n,n,PETSC_DECIDE,PETSC_DECIDE,
				d_nz,PETSC_NULL,o_nz,PETSC_NULL,&A);	
	ierr = PetscObjectSetName((PetscObject) A, "A");
	ierr = MatSetFromOptions(A);		
	
	// b
	ierr = VecCreate(PETSC_COMM_WORLD, &b);	
	ierr = PetscObjectSetName((PetscObject) b, "b");
	ierr = VecSetSizes(b, n, PETSC_DECIDE);	
	ierr = VecSetFromOptions(b);
	
	ierr = VecCreate(PETSC_COMM_WORLD,&x);
	ierr = PetscObjectSetName((PetscObject) x, "X");
	ierr = VecSetSizes(x, n, PETSC_DECIDE);	
	ierr = VecSetFromOptions(x);
}

PETSc::~PETSc()
{
	if(x!=NULL)
	{
		VecDestroy(&x);
		x = NULL;
	}
	if(b!=NULL)
	{
		VecDestroy(&b);
		b = NULL;
	}
	if(A!=NULL)
	{
		MatDestroy(&A);
		A = NULL;
	}
	if(ksp!=NULL)
	{
		KSPDestroy(&ksp);
		ksp = NULL;
	}
	if(nullsp!=NULL)
	{
		MatNullSpaceDestroy(&nullsp);
		nullsp = NULL;
	}
}

void PETSc::Reset_A()	// Reset all entries to zero ;
{
	MatZeroEntries(A);
}
void PETSc::Reset_b()  //  Reset all entries to zero ;
{
        VecZeroEntries(b);
}
void PETSc::Reset_x()
{
        VecZeroEntries(x);
}

// A
void PETSc::Set_A(PetscInt i, PetscInt j, double val)	// A[i][j]=val;
{
	ierr = MatSetValues(A,1,&i,1,&j,&val,INSERT_VALUES);
}

void PETSc::Add_A(PetscInt i, PetscInt j, double val)	// A[i][j]+=val;
{	
	ierr = MatSetValues(A,1,&i,1,&j,&val,ADD_VALUES);
}

void PETSc::Get_row_of_A(PetscInt i, PetscInt *ncol, PetscInt **cols, double **row)
{	
	ierr = MatGetRow(A,i,ncol,(const PetscInt**)cols,
			(const PetscScalar**)row);
	ierr = MatRestoreRow(A,i,ncol,(const PetscInt**)cols,
			(const PetscScalar**)row);
}

// x
void PETSc::Set_x(PetscInt i, double val)	// x[i]=val;
{
	ierr = VecSetValues(x,1,&i,&val,INSERT_VALUES);	
}

void PETSc::Add_x(PetscInt i, double val)	// x[i]+=val;
{
	ierr = VecSetValues(x,1,&i,&val,ADD_VALUES);
}

void PETSc::Set_b(PetscInt i, double val)	// x[i]=val;
{
	ierr = VecSetValues(b,1,&i,&val,INSERT_VALUES);
}

void PETSc::Add_b(
	PetscInt i, 
	double val)	// x[i]+=val;
{
	ierr = VecSetValues(b,1,&i,&val,ADD_VALUES);
}

void PETSc::Get_x(double *p)
{
	PetscScalar      *values;
	VecGetArray(x,&values);
	for(int i = 0; i < iUpper-iLower+1; i++)
		p[i] = values[i];	
        VecRestoreArray(x,&values); 
}

void PETSc::Get_b(double *p)
{
	PetscScalar      *values;
	VecGetArray(b,&values);
	for(int i = 0; i < iUpper-iLower+1; i++)
		p[i] = values[i];	
        VecRestoreArray(b,&values); 
}

void PETSc::Get_x(double *p, 
	int n, 
	int *global_index)
{
}

void PETSc::SetMaxIter(int val)
{
	PetscInt maxits;
	double rtol, atol, dtol;
	
	KSPGetTolerances(ksp, &rtol, &atol, &dtol, &maxits);
	ierr = KSPSetTolerances(ksp, rtol, atol, dtol, val);
}	/* end SetMaxIter */

void PETSc::SetTol(double val)
{
	PetscInt maxits;
	double rtol, atol, dtol;
	
	KSPGetTolerances(ksp, &rtol, &atol, &dtol, &maxits);
	ierr = KSPSetTolerances(ksp, val, atol, dtol, maxits);
}

void PETSc::SetKDim(int val)
{
	
}

void PETSc::GetNumIterations(PetscInt *num_iterations)
{
	KSPGetIterationNumber(ksp,num_iterations);        
}	/* end GetNumIterations */

void PETSc::GetFinalRelativeResidualNorm(double *rel_resid_norm)
{
	KSPGetResidualNorm(ksp,rel_resid_norm);
}	/* end GetFinalRelativeResidualNorm */

void PETSc::Solve_GMRES(void)
{
        
        start_clock("Assemble matrix and vector");
	ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
  	ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
  	
  	ierr = VecAssemblyBegin(x);
  	ierr = VecAssemblyEnd(x);
  	
  	ierr = VecAssemblyBegin(b);
  	ierr = VecAssemblyEnd(b);
	stop_clock("Assembly matrix and vector");


        KSPSetOperators(ksp,A,A);
	KSPSetType(ksp,KSPGMRES);

        KSPSetFromOptions(ksp);
        KSPSetUp(ksp);

	start_clock("KSPSolve");
        KSPSolve(ksp,b,x);
	stop_clock("KSPSolve");

}	/* end Solve_GMRES */

void PETSc::Solve(void)
{
#if defined __HYPRE__
	Solve_HYPRE();
#else // defined __HYPRE__*/
	Solve_BCGSL();
#endif // defined __HYPRE__
}	/* end Solve */

void PETSc::Solve_BCGSL(void)
{
        
        start_clock("Assemble matrix and vector");
	ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
  	ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
  	
  	ierr = VecAssemblyBegin(x);
  	ierr = VecAssemblyEnd(x);
  	
  	ierr = VecAssemblyBegin(b);
  	ierr = VecAssemblyEnd(b);
	stop_clock("Assembly matrix and vector");

        KSPSetOperators(ksp,A,A);
        KSPSetType(ksp,KSPBCGSL);
	KSPBCGSLSetEll(ksp,2);

        KSPSetFromOptions(ksp);
        KSPSetUp(ksp);

	start_clock("KSPSolve");
        KSPSolve(ksp,b,x);
	stop_clock("KSPSolve");
}

void PETSc::Solve_withPureNeumann_GMRES(void)
{
	if (debugging("trace"))
	    printf("Entering Solve_withPureNeumann_GMRES()\n");
	PC pc;
    	ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
  	ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
  	
  	ierr = VecAssemblyBegin(x);
  	ierr = VecAssemblyEnd(x);
  	
  	ierr = VecAssemblyBegin(b);
  	ierr = VecAssemblyEnd(b);
  	
	
	MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,PETSC_NULL,&nullsp);
        KSPSetNullSpace(ksp,nullsp);
	MatNullSpaceRemove(nullsp,b);

	
        KSPSetOperators(ksp,A,A);

	KSPSetType(ksp,KSPGMRES);

        KSPSetFromOptions(ksp);
	start_clock("KSPSetUp in pure neumann solver");
        KSPSetUp(ksp);
	stop_clock("KSPSetUp in pure neumann solver");
	start_clock("Petsc Solve in pure neumann solver");
        KSPSolve(ksp,b,x);
	stop_clock("Petsc Solve in pure neumann solver");
	printf("Leaving Solve_withPureNeumann_GMRES()\n");
}	/* end Solve_withPureNeumann_GMRES */

void PETSc::Solve_withPureNeumann(void)
{
	Solve_withPureNeumann_ML();
	//Solve_withPureNeumann_GMRES();
	//Solve_withPureNeumann_BCGSL();
}	/* end Solve_withPureNeumann */

void PETSc::Solve_withPureNeumann_HYPRE(void)
{
        PC pc;
	if (debugging("trace"))
        printf("Entering Solve_withPureNeumann_HYPRE()\n");
        ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
        ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);

        ierr = VecAssemblyBegin(x);
        ierr = VecAssemblyEnd(x);

        ierr = VecAssemblyBegin(b);
        ierr = VecAssemblyEnd(b);


        MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,PETSC_NULL,&nullsp);
        KSPSetNullSpace(ksp,nullsp);
        MatNullSpaceRemove(nullsp,b);

        KSPSetType(ksp,KSPBCGS);
        KSPSetOperators(ksp,A,A);
        KSPGetPC(ksp,&pc);
	start_clock("HYPRE preconditioner");
        PCSetType(pc,PCHYPRE);
        PCHYPRESetType(pc,"boomeramg");
        KSPSetFromOptions(ksp);
        KSPSetUp(ksp);
	stop_clock("HYPRE preconditioner");
        start_clock("Petsc Solve in pure neumann solver");
        KSPSolve(ksp,b,x);
        stop_clock("Petsc Solve in pure neumann solver");
	if (debugging("trace"))
	printf("Leaving Solve_withPureNeumann_HYPRE()\n");

}

void PETSc::Solve_withPureNeumann_BCGSL(void)
{
	printf("Entering Solve_withPureNeumann_BCGSL()\n");
	ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
  	ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
  	
  	ierr = VecAssemblyBegin(x);
  	ierr = VecAssemblyEnd(x);
  	
  	ierr = VecAssemblyBegin(b);
  	ierr = VecAssemblyEnd(b);
  	
	
	MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,PETSC_NULL,&nullsp);
        KSPSetNullSpace(ksp,nullsp);
	MatNullSpaceRemove(nullsp,b);
	
        start_clock("PCG in pure neumann solver");
        KSPGetPC(ksp,&pc);
        PCSetType(pc,PCGAMG);
        stop_clock("PCG in pure neumann solver");
        KSPSetOperators(ksp,A,A);
        
	KSPSetType(ksp,KSPBCGSL);
	KSPBCGSLSetEll(ksp,2);

        KSPSetFromOptions(ksp);
        KSPSetUp(ksp);

	start_clock("Petsc Solve in pure neumann solver");
        KSPSolve(ksp,b,x);
	stop_clock("Petsc Solve in pure neumann solver");
	printf("Leaving Solve_withPureNeumann_BCGSL()\n");
}	/* end Solve_withPureNeumann_BCGSL */

void PETSc::Print_A(const char *filename)
{
	PetscViewer viewer;
        ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
        ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
        PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewer);
        MatView(A, viewer);
        PetscViewerDestroy(&viewer);
}	/* end Print_A */

void PETSc::Print_b(const char *filename)
{
        ierr = VecAssemblyBegin(b);
        ierr = VecAssemblyEnd(b);
	PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD,
			PETSC_VIEWER_ASCII_MATLAB);
        VecView(b, PETSC_VIEWER_STDOUT_WORLD);
}	/* end Print_b */

extern void viewTopVariable(
	Front *front,
	double *var,
	boolean set_bounds,
	double var_min,
	double var_max,
	char *dirname,
	char *var_name)
{
	HDF_MOVIE_VAR hdf_movie_var;
	HDF_MOVIE_VAR *hdf_movie_var_save = front->hdf_movie_var;
	front->hdf_movie_var = &hdf_movie_var;
	hdf_movie_var.num_var = 1;
	FT_MatrixMemoryAlloc((POINTER*)&hdf_movie_var.var_name,1,100,
				sizeof(char));
	FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var.top_var,1,
				sizeof(double*));
	FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var.preset_bound,1,
				sizeof(boolean));
	FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var.var_min,1,
				sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var.var_max,1,
				sizeof(double));
	sprintf(hdf_movie_var.var_name[0],"%s",var_name);
	hdf_movie_var.preset_bound[0] = set_bounds;
	hdf_movie_var.var_min[0] = var_min;
	hdf_movie_var.var_max[0] = var_max;
	hdf_movie_var.top_var[0] = var;
	gview_var2d_on_top_grid(front,dirname);

	FT_FreeThese(5,hdf_movie_var.var_name,hdf_movie_var.top_var,
				hdf_movie_var.preset_bound,
				hdf_movie_var.var_min,hdf_movie_var.var_max);
	front->hdf_movie_var = hdf_movie_var_save;
}	/* end viewTopVariable */

#if defined __HYPRE__
void PETSc::Solve_HYPRE(void)
{
        PC pc;
        start_clock("Assemble matrix and vector");
        ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
        ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);

        ierr = VecAssemblyBegin(x);
        ierr = VecAssemblyEnd(x);

        ierr = VecAssemblyBegin(b);
        ierr = VecAssemblyEnd(b);
        stop_clock("Assembly matrix and vector");

	KSPSetType(ksp,KSPBCGS);
        KSPSetOperators(ksp,A,A);
        KSPGetPC(ksp,&pc);
	PCSetType(pc,PCHYPRE);
        PCHYPRESetType(pc,"boomeramg");
        KSPSetFromOptions(ksp);
        KSPSetUp(ksp);

        start_clock("KSPSolve");
        KSPSolve(ksp,b,x);
        stop_clock("KSPSolve");

}
#endif // defined __HYPRE__

void PETSc::Solve_withPureNeumann_ML(void)
{
	if (debugging("trace"))
	    printf("Entering Solve_withPureNeumann_ML()\n");
	PC pc;
	start_clock("Assemble Matrix in pure neumann solver");
    	ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
  	ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
  	
  	ierr = VecAssemblyBegin(x);
  	ierr = VecAssemblyEnd(x);
  	
  	ierr = VecAssemblyBegin(b);
  	ierr = VecAssemblyEnd(b);
	stop_clock("Assemble Matrix in pure neumann solver");
  	
	
	MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,PETSC_NULL,&nullsp);
        KSPSetNullSpace(ksp,nullsp);
	MatNullSpaceRemove(nullsp,b);
	
        KSPSetOperators(ksp,A,A);

        KSPGetPC(ksp,&pc);
        PCSetType(pc,PCML);

	KSPSetType(ksp,KSPGMRES);

        KSPSetFromOptions(ksp);
	start_clock("KSP setup in pure neumann solver");
        KSPSetUp(ksp);
	stop_clock("KSP setup in pure neumann solver");
	start_clock("Petsc Solve in pure neumann solver");
        KSPSolve(ksp,b,x);
	stop_clock("Petsc Solve in pure neumann solver");
	printf("Leaving Solve_withPureNeumann_ML()\n");
}	/* end Solve_withPureNeumann_GMRES */


void PETSc::Solve_LU(void)
{
	PC pc;
        start_clock("Assemble matrix and vector");
        ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
        ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);

        ierr = VecAssemblyBegin(x);
        ierr = VecAssemblyEnd(x);

        ierr = VecAssemblyBegin(b);
        ierr = VecAssemblyEnd(b);
        stop_clock("Assembly matrix and vector");


        KSPSetType(ksp,KSPPREONLY);
	KSPGetPC(ksp,&pc);
	PCSetType(pc,PCLU);
        KSPSetOperators(ksp,A,A);
        KSPSetFromOptions(ksp);
        KSPSetUp(ksp);

        start_clock("KSPSolve");
        KSPSolve(ksp,b,x);
        stop_clock("KSPSolve");
} /*direct solver, usually give exact solution for comparison*/

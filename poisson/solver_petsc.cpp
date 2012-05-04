/*		solver_petsc.c
 *  This class PETSc is created to be used a handy interface
 *  to the function calls to PETSc. For each algebric equation 
 *  Ax=b, one PETSc instance is needed, since this instance has 
 *  the storage of these three variables. 
*/ 
#include <stdio.h>
#include <stdlib.h>

#include <solver_petsc.h>

PETSc::PETSc()
{
	x = NULL;			/* approx solution, RHS*/
	b = NULL;
  	A = NULL;            		/* linear system matrix */
  	
  	ksp = NULL;        		/* Krylov subspace method context */
	nullsp = NULL;
	
	KSPCreate(PETSC_COMM_WORLD,&ksp);
}

PETSc::PETSc(int ilower, int iupper, int d_nz, int o_nz)
{	
	x = NULL;      			/* approx solution, RHS*/
	b = NULL;
  	A = NULL;            		/* linear system matrix */
  	
  	ksp = NULL;          		/* Krylov subspace method context */
	nullsp = NULL;
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
	
	MatCreateMPIAIJ(PETSC_COMM_WORLD,n,n,PETSC_DECIDE,PETSC_DECIDE,
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
		VecDestroy(x);
		x = NULL;
	}
	if(b!=NULL)
	{
		VecDestroy(b);
		b = NULL;
	}
	if(A!=NULL)
	{
		MatDestroy(A);
		A = NULL;
	}
	if(ksp!=NULL)
	{
		KSPDestroy(ksp);
		ksp = NULL;
	}
	if(nullsp!=NULL)
	{
		MatNullSpaceDestroy(nullsp);
		nullsp = NULL;
	}
}

void PETSc::Reset_A()	// Reset all entries to zero ;
{
	MatZeroEntries(A);
}

// A
void PETSc::Set_A(int i, int j, double val)	// A[i][j]=val;
{
	ierr = MatSetValues(A,1,&i,1,&j,&val,INSERT_VALUES);
}

void PETSc::Add_A(int i, int j, double val)	// A[i][j]+=val;
{	
	ierr = MatSetValues(A,1,&i,1,&j,&val,ADD_VALUES);
}

void PETSc::Get_row_of_A(int i, int *ncol, int **cols, double **row)
{	
	ierr = MatGetRow(A,i,ncol,(const PetscInt**)cols,
			(const PetscScalar**)row);
	ierr = MatRestoreRow(A,i,ncol,(const PetscInt**)cols,
			(const PetscScalar**)row);
}

// x
void PETSc::Set_x(int i, double val)	// x[i]=val;
{
	ierr = VecSetValues(x,1,&i,&val,INSERT_VALUES);	
}

void PETSc::Add_x(int i, double val)	// x[i]+=val;
{
	ierr = VecSetValues(x,1,&i,&val,ADD_VALUES);
}

void PETSc::Set_b(int i, double val)	// x[i]=val;
{
	ierr = VecSetValues(b,1,&i,&val,INSERT_VALUES);
}

void PETSc::Add_b(
	int i, 
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
	int maxits;
	double rtol, atol, dtol;
	
	KSPGetTolerances(ksp, &rtol, &atol, &dtol, &maxits);
	ierr = KSPSetTolerances(ksp, rtol, atol, dtol, val);
}

void PETSc::SetTol(double val)
{
	int maxits;
	double rtol, atol, dtol;
	
	KSPGetTolerances(ksp, &rtol, &atol, &dtol, &maxits);
	ierr = KSPSetTolerances(ksp, val, atol, dtol, maxits);
}

void PETSc::SetKDim(int val)
{
	
}

void PETSc::GetNumIterations(int *num_iterations)
{
	KSPGetIterationNumber(ksp,num_iterations);        
}

void PETSc::GetFinalRelativeResidualNorm(double *rel_resid_norm)
{
	KSPGetResidualNorm(ksp,rel_resid_norm);
}

void PETSc::GetExtremeSingularValues(double *max, double *min)
{
     KSPComputeExtremeSingularValues(ksp,max,min);
}

void PETSc::Solve(void)
{
	ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
  	ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
  	
  	ierr = VecAssemblyBegin(x);
  	ierr = VecAssemblyEnd(x);
  	
  	ierr = VecAssemblyBegin(b);
  	ierr = VecAssemblyEnd(b);
	
        KSPSetOperators(ksp,A,A,DIFFERENT_NONZERO_PATTERN);
        KSPSetType(ksp,KSPBCGS);

	KSPSetComputeSingularValues(ksp, PETSC_TRUE);
        KSPSetFromOptions(ksp);
        KSPSetUp(ksp);

        KSPSolve(ksp,b,x);
}

void PETSc::Solve_withPureNeumann(void)
{
	ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
  	ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
  	
  	ierr = VecAssemblyBegin(x);
  	ierr = VecAssemblyEnd(x);
  	
  	ierr = VecAssemblyBegin(b);
  	ierr = VecAssemblyEnd(b);
  	
	
	MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,PETSC_NULL,&nullsp);
        KSPSetNullSpace(ksp,nullsp);
	
        KSPSetOperators(ksp,A,A,DIFFERENT_NONZERO_PATTERN);
        
	//KSPSetType(ksp,KSPBCGS);
	KSPSetType(ksp,KSPGMRES);
        KSPSetFromOptions(ksp);
        KSPSetUp(ksp);
        KSPSolve(ksp,b,x);
}

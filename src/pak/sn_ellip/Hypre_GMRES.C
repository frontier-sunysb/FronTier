/*      This class Hypre is created to be used a handy interface
 *  to the function calls to Hypre. For each algebric equation 
 *  Ax=b, one Hypre_GMRES instance is needed, since this instance has 
 *  the storage of these three variables. 
*/ 
#include <stdio.h>
#include <stdlib.h>

#include "utilities.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_mv.h"
#include "HYPRE_IJ_mv.h"
#include "HYPRE_parcsr_ls.h"

#include <Hypre_GMRES.h>
#include <gutil.h>

Hypre_GMRES::Hypre_GMRES(int ilower, int iupper)
{	
	comm 	= MPI_COMM_WORLD;
	iLower	= ilower;	iUpper = iupper;
	// A
	HYPRE_IJMatrixCreate(comm, iLower, iUpper, iLower, iUpper, &ij_A);
	HYPRE_IJMatrixSetObjectType(ij_A, HYPRE_PARCSR);
	HYPRE_IJMatrixInitialize(ij_A);
	// b
	HYPRE_IJVectorCreate(comm, iLower, iUpper, &ij_b);
	HYPRE_IJVectorSetObjectType(ij_b, HYPRE_PARCSR);
	HYPRE_IJVectorInitialize(ij_b);
	// x
	HYPRE_IJVectorCreate(comm, iLower, iUpper, &ij_x);
	HYPRE_IJVectorSetObjectType(ij_x, HYPRE_PARCSR);
	HYPRE_IJVectorInitialize(ij_x);
	
	// solver
	HYPRE_ParCSRGMRESCreate(comm, &solver);
	HYPRE_ParCSRGMRESSetMaxIter(solver, 500);     	 //optional 
	HYPRE_ParCSRGMRESSetTol(solver, 1.0e-07);    	 //optional 
	
}

void Hypre_GMRES::Create(int ilower, int iupper)
{	
	comm 	= MPI_COMM_WORLD;
	iLower	= ilower;	iUpper = iupper;
	// A
	HYPRE_IJMatrixCreate(comm, iLower, iUpper, iLower, iUpper, &ij_A);
	HYPRE_IJMatrixSetObjectType(ij_A, HYPRE_PARCSR);
	HYPRE_IJMatrixInitialize(ij_A);
	// b
	HYPRE_IJVectorCreate(comm, iLower, iUpper, &ij_b);
	HYPRE_IJVectorSetObjectType(ij_b, HYPRE_PARCSR);
	HYPRE_IJVectorInitialize(ij_b);
	// x
	HYPRE_IJVectorCreate(comm, iLower, iUpper, &ij_x);
	HYPRE_IJVectorSetObjectType(ij_x, HYPRE_PARCSR);
	HYPRE_IJVectorInitialize(ij_x);
	
	// solver
	HYPRE_ParCSRGMRESCreate(comm, &solver);
	HYPRE_ParCSRGMRESSetMaxIter(solver, 500);     	 //optional 
	HYPRE_ParCSRGMRESSetTol(solver, 1.0e-07);    	 //optional 
	
}

void Hypre_GMRES::Create(MPI_Comm Comm, int ilower, int iupper)
{	
	comm 	= Comm;
	iLower	= ilower;	iUpper = iupper;
	// A
	HYPRE_IJMatrixCreate(comm, iLower, iUpper, iLower, iUpper, &ij_A);
	HYPRE_IJMatrixSetObjectType(ij_A, HYPRE_PARCSR);
	HYPRE_IJMatrixInitialize(ij_A);
	// b
	HYPRE_IJVectorCreate(comm, iLower, iUpper, &ij_b);
	HYPRE_IJVectorSetObjectType(ij_b, HYPRE_PARCSR);
	HYPRE_IJVectorInitialize(ij_b);
	// x
	HYPRE_IJVectorCreate(comm, iLower, iUpper, &ij_x);
	HYPRE_IJVectorSetObjectType(ij_x, HYPRE_PARCSR);
	HYPRE_IJVectorInitialize(ij_x);
	
	// solver
	HYPRE_ParCSRGMRESCreate(comm, &solver);
	HYPRE_ParCSRGMRESSetMaxIter(solver, 500);     	 //optional 
	HYPRE_ParCSRGMRESSetTol(solver, 1.0e-07);    	 //optional 
	//HYPRE_BoomerAMGSetStrongThreshold(solver, 0.8);  //strong_threshold: 0-1	
}

void Hypre_GMRES::Initialize(void)
{
	HYPRE_IJMatrixInitialize(ij_A);
	HYPRE_IJVectorInitialize(ij_b);
	HYPRE_IJVectorInitialize(ij_x);
}

Hypre_GMRES::~Hypre_GMRES()
{
	HYPRE_ParCSRGMRESDestroy(solver);		
	HYPRE_IJMatrixDestroy(ij_A);
	HYPRE_IJVectorDestroy(ij_b);
	HYPRE_IJVectorDestroy(ij_x);
}
// A
void Hypre_GMRES::Set_A(int i, int j, double val)	// A[i][j]=val;
{
	int  nrows, ncols[1], rows[1], cols[1];
	double values[1];	
	
	nrows=1;   ncols[0]=1;   rows[0]=i;   cols[0]=j;   values[0]=val;
	HYPRE_IJMatrixSetValues(ij_A, nrows, ncols, rows, cols, values);	
	//Debug_Print("\n%d %d %f", i, j, val);
}
void Hypre_GMRES::Add_A(int i, int j, double val)	// A[i][j]+=val;
{
	int  nrows, ncols[1], rows[1], cols[1];
	double values[1];	
	
	nrows=1;   ncols[0]=1;   rows[0]=i;   cols[0]=j;   values[0]=val;
	HYPRE_IJMatrixAddToValues(ij_A, nrows, ncols, rows, cols, values);					
}
// x
void Hypre_GMRES::Set_x(int i, double val)	// x[i]=val;
{
	int nvalues, indices[1];
	double values[1];
	 
	nvalues=1;   indices[0]=i;   values[0]=val;
	HYPRE_IJVectorSetValues(ij_x, nvalues, indices, values);	
}
void Hypre_GMRES::Add_x(int i, double val)	// x[i]+=val;
{
	int nvalues, indices[1];
	double values[1];
	 
	nvalues=1;   indices[0]=i;   values[0]=val;
	HYPRE_IJVectorAddToValues(ij_x, nvalues, indices, values);	
}
// b
void Hypre_GMRES::Set_b(int i, double val)	// x[i]=val;
{
	int nvalues, indices[1];
	double values[1];
	 
	nvalues=1;   indices[0]=i;   values[0]=val;
	HYPRE_IJVectorSetValues(ij_b, nvalues, indices, values);	
}
void Hypre_GMRES::Add_b(int i, double val)	// x[i]+=val;
{
	int nvalues, indices[1];
	double values[1];
	 
	nvalues=1;   indices[0]=i;   values[0]=val;
	HYPRE_IJVectorAddToValues(ij_b, nvalues, indices, values);	
}
void Hypre_GMRES::Set_b(double *b)
{
        int i, nvalues, indices[1];
        int index=0;
        for(i=iLower; i<=iUpper; i++, index++)  // index is for array b who's index starts from 0;
        {
                nvalues=1;      indices[0]=i;
	        HYPRE_IJVectorSetValues(ij_b, nvalues, indices, &b[index]);
        }
}

									
void Hypre_GMRES::Get_x(double *p)
{
	int i, nvalues, indices[1];
	int index=0;	
	for(i=iLower; i<=iUpper; i++, index++)	// index is for array p who's index starts from 0;
	{		
		nvalues=1;	indices[0]=i;    		
		HYPRE_IJVectorGetValues(ij_x, nvalues, indices, &p[index]);
	}	
	
}

void Hypre_GMRES::Set_x(double *p)
{
        int i, nvalues, indices[1];
        int index=0;
        for(i=iLower; i<=iUpper; i++, index++)  // index is for array p who's index starts from 0;
        {
                nvalues=1;      indices[0]=i;
                HYPRE_IJVectorSetValues(ij_x, nvalues, indices, &p[index]);
        }

}
									

void Hypre_GMRES::SetPrecond(void)
{
//	HYPRE_BoomerAMGCreate(&pcg_precond);
//	HYPRE_GMRESSetPrecond(solver,
 //                            (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
//                            (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup,
//                                   pcg_precond);
}

void Hypre_GMRES::SetMaxIter(int val)
{
	HYPRE_ParCSRGMRESSetMaxIter(solver, val);     	  //optional 
}
void Hypre_GMRES::SetTol(double val)
{
	HYPRE_ParCSRGMRESSetTol(solver, val);    	  //optional 
}
void Hypre_GMRES::SetKDim(int val)
{
	HYPRE_ParCSRGMRESSetKDim(solver, val);    	  //optional 
}
void Hypre_GMRES::GetNumIterations(int *num_iterations)
{
	HYPRE_ParCSRGMRESGetNumIterations(solver, num_iterations); //Return the number of iterations taken 
}
void Hypre_GMRES::GetFinalRelativeResidualNorm(double *rel_resid_norm)
{
	HYPRE_ParCSRGMRESGetFinalRelativeResidualNorm(solver, rel_resid_norm);
}


void Hypre_GMRES::Solve(void)
{
	HYPRE_IJMatrixAssemble(ij_A);
	HYPRE_IJMatrixGetObject(ij_A, (void **) &par_A);
	HYPRE_IJVectorAssemble(ij_b);
	HYPRE_IJVectorGetObject(ij_b, (void **) &par_b);
	HYPRE_IJVectorAssemble(ij_x);
	HYPRE_IJVectorGetObject(ij_x, (void **) &par_x);
	
	// HYPRE_ParCSRGMRESCreate(&solver); 	// solver is created in Hypre_GMRES::Hypre_GMRES()
	HYPRE_ParCSRGMRESSetup(solver, par_A,  par_b,  par_x);
	HYPRE_ParCSRGMRESSolve(solver, par_A,  par_b,  par_x);
}

// this function should be used after Hypre_GMRES::Solve()
void Hypre_GMRES::BoomerAMGSolve(void)
{
	// solver
	HYPRE_BoomerAMGCreate(&pcg_precond);
	HYPRE_BoomerAMGSetMaxIter(pcg_precond, 500);     	 //optional 
	HYPRE_BoomerAMGSetTol(pcg_precond, 1.0e-07);    	 //optional 
	
	HYPRE_BoomerAMGSetup(pcg_precond, par_A,  par_b,  par_x);
	HYPRE_BoomerAMGSolve(pcg_precond, par_A,  par_b,  par_x);

	
	HYPRE_BoomerAMGDestroy(pcg_precond);			
}

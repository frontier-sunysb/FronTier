/*      This Hypre class is created to be used as a handy interface
 *  to the function calls of the Algebric Multigrid Method from the 
 *  hypre library.  There are other methods in the hypre library, such 
 *  as GMRES, etc. But their interfaces are not implemented in the 
 *  Hypre class yet.
 *  
 *  To know how to use the Hypre class, please read the two example 
 * programs: 
 *	Hypre_main0.c for running on one processor;
 *	Hypre_main.c  for running on two processors.
 *
 *  To know how to compile the Hypre class, please read the makefile.
*/ 
#ifndef  CLASS_Hypre_GMRES_H
#define  CLASS_Hypre_GMRES_H

#include <stdio.h>
#include <stdlib.h>

#include "utilities.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_mv.h"
#include "HYPRE_IJ_mv.h"
#include "HYPRE_parcsr_ls.h"
/*#include "krylov.h" */

#include <Solver.h>

class Hypre_GMRES: public SOLVER
{
	MPI_Comm  comm;				/* set to be MPI_COMM_WORLD. */
	int iLower;
	int iUpper;				/* global row range */
	
	HYPRE_IJMatrix      	ij_A;
	HYPRE_ParCSRMatrix  	par_A;
	HYPRE_IJVector		ij_b;
	HYPRE_ParVector		par_b;
	HYPRE_IJVector		ij_x;
	HYPRE_ParVector		par_x;
	HYPRE_Solver		solver;
	HYPRE_Solver		pcg_precond;
		
public:
	Hypre_GMRES()
	{
	};
	Hypre_GMRES(int ilower, int iupper);		/* global row range of A, x, b on this processor */
	~Hypre_GMRES();
	void Create(int ilower, int iupper);	/* same as Hypre(int, int) */
	void Create(MPI_Comm Comm, int ilower, int iupper);	/* same as Hypre(int, int) */

	void Initialize(void);			/* Initialize A, b, x for value correction. */
	void Set_A(int i, int j, double val);	/* A[i][j]=val; */
	void Add_A(int i, int j, double val);	/* A[i][j]=A[i][j]+val; */
	void Set_x(int i, double val);		/* x[i]=val; */
	void Add_x(int i, double val);		/* x[i]=x[i]+val; */
	void Set_b(int i, double val);		/* b[i]=val; */
	void Add_b(int i, double val);		/* b[i]=b[i]+val; */
	void Set_b(double *b);
	void Get_x(double *p);		/* get the x from ij_x to p.	 */
	void Set_x(double *p);

	void SetPrecond(void);		/* set AMG preconditioner */
	void SetMaxIter(int val);	/* (optional) Set maximum number of iterations  */
	void SetTol(double val);	/* (Optional) Set the convergence tolerance  */
	void SetKDim(int k_dim);	/* (Optional) Set the maximum size of the Krylov space  */
	void GetNumIterations(int *num_iterations);	/* (Optional) Return the number of iterations taken  */
	void GetFinalRelativeResidualNorm(double *rel_resid_norm);
	void Solve(void);
	void BoomerAMGSolve(void);		/* call BoomerAMG() , this function should be used after Hypre_GMRES::Solve() */
	
	 void Read_A(char *filename)
	{
		HYPRE_IJMatrixRead(filename, comm, HYPRE_PARCSR, &ij_A);
	};
	void Print_A(char *filename)	/* (for debug) output x to a file with name 'filename' */
	{	
		HYPRE_IJMatrixAssemble(ij_A);				/* is this redundant? */
		HYPRE_IJMatrixGetObject(ij_A, (void **) &par_A);		
		HYPRE_IJMatrixPrint(ij_A, filename);
	};
	void Read_b(char *filename)
        {       
                HYPRE_IJVectorRead(filename, comm, HYPRE_PARCSR, &ij_b);
        };                              
        void Print_b(char *filename)    /* (for debug) output x to a file with name 'filename'  */
        {       
		HYPRE_IJVectorAssemble(ij_b);
		HYPRE_IJVectorGetObject(ij_b, (void **) &par_b);
		HYPRE_IJVectorPrint(ij_b, filename);    
        };                
	void Read_x(char *filename)
	{
		HYPRE_IJVectorRead(filename, comm, HYPRE_PARCSR, &ij_x);
	};
	void Print_x(char *filename)	/* (for debug) output x to a file with name 'filename' */
	{	
		HYPRE_IJVectorAssemble(ij_x);
		HYPRE_IJVectorGetObject(ij_x, (void **) &par_x);
		HYPRE_IJVectorPrint(ij_x, filename);
	};
		
	void SetIOutDat(int ioutdat)    /* GMRES has no such function; */
	{
	};
};

#endif  /* #ifndef CLASS_HYPRE_H */


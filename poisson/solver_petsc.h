/*      
 * 
 * File:   solver_petsc.h
 * Author: shuqiang (robert) wang 
 *
 * To use this file, the following function has to be used:
 *	PetscInitialize(&argc,&args,(char *)0,help)
 * to initilize the Petsc package.
 *
*/ 

#ifndef  _PETSC_H
#define  _PETSC_H

#if defined(c_plusplus) || defined(__cplusplus)
#undef c_plusplus
#undef __cplusplus
extern "C"
{
#include "petscksp.h"
}
#define c_plusplus 1
#define __cplusplus 1
#endif

#include "solver.h"


class PETSc: public SOLVER
{
public:	
	MPI_Comm  comm;			// set to be MPI_COMM_WORLD.
	int iLower;
	int iUpper;			// global row range
	
	Vec x;      			/* approx solution, RHS*/
	Vec b;
  	Mat A;          		/* linear system matrix */
  	
  	KSP   ksp;          		/* Krylov subspace method context */
	MatNullSpace	nullsp;
		
	int ierr;
	int its;			// numer of iterations;

public:
	PETSc();
	PETSc(int ilower, int iupper, int d_nz, int o_nz);		
		// global row range of A, x, b on this processor
	~PETSc();
	void Create(int ilower, int iupper, int d_nz, int o_nz);	
		// same as Hypre(int, int)
	void Create(MPI_Comm Comm, int ilower, int iupper, int d_nz, int o_nz);	
		// same as Hypre(int, int)
	
	void Reset_A();				// Set A[i][j]=0.0;
	void Set_A(int i, int j, double val);	// A[i][j]=val;
	void Add_A(int i, int j, double val);	// A[i][j]=A[i][j]+val;
	void Get_row_of_A(int i, int *ncol, int **cols, double **row);
	void Set_x(int i, double val);		// x[i]=val;
	void Add_x(int i, double val);		// x[i]=x[i]+val;
	void Set_b(int i, double val);		// b[i]=val;
	void Add_b(int i, double val);		// b[i]=b[i]+val;
	void Get_x(double *p);		// get the x from ij_x to p.	
	void Get_b(double *p);		// get the b from ij_x to p.	
	void Get_x(double *p, int n, int *global_index);
	
	void SetMaxIter(int val); 	// Set maximum number of iterations 
	void SetTol(double val);	// Set the convergence tolerance 
	void SetKDim(int k_dim);	
			// Set the maximum size of the Krylov space 
	void GetNumIterations(int *num_iterations);	
			// Return the number of iterations taken 
	void GetFinalRelativeResidualNorm(double *rel_resid_norm);
	void GetExtremeSingularValues(double *max, double *min);
	void Solve(void);
	void Solve_withPureNeumann(void);
};
#endif


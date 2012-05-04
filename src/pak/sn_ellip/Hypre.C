/*      This class Hypre is created to be used a handy interface
 *  to the function calls to Hypre. For each algebric equation 
 *  Ax=b, one Hypre instance is needed, since this instance has 
 *  the storage of these three variables. 
*/ 
#include <stdio.h>
#include <stdlib.h>

#include "utilities.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_mv.h"
#include "HYPRE_IJ_mv.h"
#include "HYPRE_parcsr_ls.h"

#include <Hypre.h>

	


Hypre::Hypre(int ilower, int iupper)
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
	HYPRE_BoomerAMGCreate(&solver);
	HYPRE_BoomerAMGSetMaxIter(solver, 500);     	 //optional 
	HYPRE_BoomerAMGSetTol(solver, 1.0e-07);    	 //optional 
	//HYPRE_BoomerAMGSetStrongThreshold(solver, 0.8);  //strong_threshold: 0-1	
	//int grid_relax_type[] ={ 3, 3, 3, 3};
	//HYPRE_BoomerAMGSetGridRelaxType(solver, grid_relax_type);
	
}

void Hypre::Create(int ilower, int iupper)
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

	grid_relax_type = hypre_CTAlloc(int,4);
	grid_relax_type[0] = 6;
	grid_relax_type[1] = 6;
	grid_relax_type[2] = 6;
	grid_relax_type[3] = 6;

	num_grid_sweeps = hypre_CTAlloc(int,4);
	
	grid_relax_points = hypre_CTAlloc(int *,4);

    // fine grid 
    num_grid_sweeps[0] = 3;
    grid_relax_points[0] = hypre_CTAlloc(int,3);
    grid_relax_points[0][0] = -2;
    grid_relax_points[0][1] = -1;
    grid_relax_points[0][2] = 1;
   
   // /* down cycle */
    num_grid_sweeps[1] = 4;
    grid_relax_points[1] = hypre_CTAlloc(int,4); 
    grid_relax_points[1][0] = -1;
    grid_relax_points[1][1] = 1;
    grid_relax_points[1][2] = -2;
    grid_relax_points[1][3] = -2;
   
    /* up cycle */
    num_grid_sweeps[2] = 4;
    grid_relax_points[2] = hypre_CTAlloc(int,4); 
    grid_relax_points[2][0] = -2;
    grid_relax_points[2][1] = -2;
    grid_relax_points[2][2] = 1;
    grid_relax_points[2][3] = -1;
      
	/* coarsest grid */
	num_grid_sweeps[3] = 1;
	grid_relax_points[3] = hypre_CTAlloc(int,1);
	grid_relax_points[3][0] = 0;

	

	//int max_levels = 25;
	//double *relax_weight      = hypre_CTAlloc(double, max_levels);
	//int	   *smooth_option	  = hypre_CTAlloc(int, max_levels);
	//int      smooth_num_sweep = 1;

	//for (int i=0; i < max_levels; i++)
	//{
	//	relax_weight[i] = 1.;
	//	smooth_option[i] = -1;
	//}


	HYPRE_BoomerAMGCreate(&solver);
	HYPRE_BoomerAMGSetCoarsenType(solver, 1);
	HYPRE_BoomerAMGSetMeasureType(solver, 1);
	HYPRE_BoomerAMGSetTol(solver, 1.0e-07);    	 //optional 
	HYPRE_BoomerAMGSetStrongThreshold(solver, 0.75);  //strong_threshold: 0-1	
	//HYPRE_BoomerAMGSetTruncFactor(solver, 0.0);
	//HYPRE_BoomerAMGSetCycleType(solver, 1);
	
	//HYPRE_BoomerAMGSetNumGridSweeps(solver, num_grid_sweeps);
	//HYPRE_BoomerAMGSetGridRelaxType(solver, grid_relax_type);
	//HYPRE_BoomerAMGSetGridRelaxPoints(solver, grid_relax_points);
	//HYPRE_BoomerAMGSetRelaxWeight(solver, relax_weight);
	//HYPRE_BoomerAMGSetSmoothOption(solver, smooth_option);
	//HYPRE_BoomerAMGSetSmoothNumSweep(solver, smooth_num_sweep);
	//HYPRE_BoomerAMGSetMaxLevels(solver, max_levels);

	//HYPRE_BoomerAMGSetMaxRowSum(solver, 1.0);
	//HYPRE_BoomerAMGSetVariant(solver, 0);
	//HYPRE_BoomerAMGSetOverlap(solver, 1);	
	HYPRE_BoomerAMGSetMaxIter(solver, 500);     	 //optional 
	
	

}

void Hypre::Create(MPI_Comm Comm, int ilower, int iupper)
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
	HYPRE_BoomerAMGCreate(&solver);
	HYPRE_BoomerAMGSetMaxIter(solver, 500);     	 //optional 
	HYPRE_BoomerAMGSetTol(solver, 1.0e-07);    	 //optional 
	//HYPRE_BoomerAMGSetStrongThreshold(solver, 0.8);  //strong_threshold: 0-1	
}


Hypre::~Hypre()
{
	HYPRE_BoomerAMGDestroy(solver);		
	HYPRE_IJMatrixDestroy(ij_A);
	HYPRE_IJVectorDestroy(ij_b);
	HYPRE_IJVectorDestroy(ij_x);
	//delete grid_relax_type;			// for BoomerAMG only
	//delete num_grid_sweeps;	
	//delete grid_relax_points[0];	
	//delete grid_relax_points[1];	
	//delete grid_relax_points[2];	
	//delete grid_relax_points[3];	
	//delete grid_relax_points;
	
}
// A
void Hypre::Set_A(int i, int j, double val)	// A[i][j]=val;
{
	int  nrows, ncols[1], rows[1], cols[1];
	double values[1];	
	
	nrows=1;   ncols[0]=1;   rows[0]=i;   cols[0]=j;   values[0]=val;
	HYPRE_IJMatrixSetValues(ij_A, nrows, ncols, rows, cols, values);					
}
void Hypre::Add_A(int i, int j, double val)	// A[i][j]+=val;
{
	int  nrows, ncols[1], rows[1], cols[1];
	double values[1];	
	
	nrows=1;   ncols[0]=1;   rows[0]=i;   cols[0]=j;   values[0]=val;
	HYPRE_IJMatrixAddToValues(ij_A, nrows, ncols, rows, cols, values);					
}
// x
void Hypre::Set_x(int i, double val)	// x[i]=val;
{
	int nvalues, indices[1];
	double values[1];
	 
	nvalues=1;   indices[0]=i;   values[0]=val;
	HYPRE_IJVectorSetValues(ij_x, nvalues, indices, values);	
}
void Hypre::Add_x(int i, double val)	// x[i]+=val;
{
	int nvalues, indices[1];
	double values[1];
	 
	nvalues=1;   indices[0]=i;   values[0]=val;
	HYPRE_IJVectorAddToValues(ij_x, nvalues, indices, values);	
}
// b
void Hypre::Set_b(int i, double val)	// x[i]=val;
{
	int nvalues, indices[1];
	double values[1];
	 
	nvalues=1;   indices[0]=i;   values[0]=val;
	HYPRE_IJVectorSetValues(ij_b, nvalues, indices, values);	
}
void Hypre::Add_b(int i, double val)	// x[i]+=val;
{
	int nvalues, indices[1];
	double values[1];
	 
	nvalues=1;   indices[0]=i;   values[0]=val;
	HYPRE_IJVectorAddToValues(ij_b, nvalues, indices, values);	
}

void Hypre::Set_b(double *b)
{
        int i, nvalues, indices[1];
        int index=0;
        for(i=iLower; i<=iUpper; i++, index++)  // index is for array b who's index starts from 0;
        {
                nvalues=1;      indices[0]=i;
                HYPRE_IJVectorSetValues(ij_b, nvalues, indices, &b[index]);
        }
}
									

void Hypre::Get_x(double *p)
{
	int i, nvalues, indices[1];
	int index=0;	
	for(i=iLower; i<=iUpper; i++, index++)	// index is for array p who's index starts from 0;
	{		
		nvalues=1;	indices[0]=i;    		
		HYPRE_IJVectorGetValues(ij_x, nvalues, indices, &p[index]);
	}	
}
void Hypre::Set_x(double *p)
{
	int i, nvalues, indices[1];
	int index=0;	
	for(i=iLower; i<=iUpper; i++, index++)	// index is for array p who's index starts from 0;
	{		
		nvalues=1;	indices[0]=i;    		
		HYPRE_IJVectorSetValues(ij_x, nvalues, indices, &p[index]);
	}	
}


void Hypre::SetMaxIter(int val)
{
	HYPRE_BoomerAMGSetMaxIter(solver, val);     	  //optional 
}
void Hypre::SetTol(double val)
{
	HYPRE_BoomerAMGSetTol(solver, val);    	  //optional 
}
void Hypre::SetStrongThreshold(double val)	// for AMG
{
	HYPRE_BoomerAMGSetStrongThreshold(solver, val);  //strong_threshold: 0-1
}
void Hypre::GetNumIterations(int *num_iterations)
{
	HYPRE_BoomerAMGGetNumIterations(solver, num_iterations); //Return the number of iterations taken 
}
void Hypre::GetFinalRelativeResidualNorm(double *rel_resid_norm)
{
	HYPRE_BoomerAMGGetFinalRelativeResidualNorm(solver, rel_resid_norm);
}

void Hypre::Solve(void)
{
	HYPRE_IJMatrixAssemble(ij_A);
	HYPRE_IJMatrixGetObject(ij_A, (void **) &par_A);
	HYPRE_IJVectorAssemble(ij_b);
	HYPRE_IJVectorGetObject(ij_b, (void **) &par_b);
	HYPRE_IJVectorAssemble(ij_x);
	HYPRE_IJVectorGetObject(ij_x, (void **) &par_x);
	
	// HYPRE_BoomerAMGCreate(&solver); 	// solver is created in Hypre::Hypre()
	HYPRE_BoomerAMGSetup(solver, par_A,  par_b,  par_x);
	HYPRE_BoomerAMGSolve(solver, par_A,  par_b,  par_x);
	
}


void Hypre::Solve2(int total_iterations)
{
	HYPRE_IJMatrixAssemble(ij_A);
	HYPRE_IJMatrixGetObject(ij_A, (void **) &par_A);
	HYPRE_IJVectorAssemble(ij_b);
	HYPRE_IJVectorGetObject(ij_b, (void **) &par_b);
	HYPRE_IJVectorAssemble(ij_x);
	HYPRE_IJVectorGetObject(ij_x, (void **) &par_x);

	// create solver, note that solver is created in Hypre()
	HYPRE_ParCSRGMRESCreate(comm, &gmres_solver);
	// set solver parameters
	HYPRE_ParCSRGMRESSetMaxIter(gmres_solver, 500);     	 //optional 
	HYPRE_ParCSRGMRESSetTol(gmres_solver, 1.0e-07);    	 //optional	
	HYPRE_BoomerAMGSetMaxIter(solver, 5);

	// set up the solver and 
	int i;
	for(i=0; i<total_iterations; i++)
	{
		// GMRES
		HYPRE_ParCSRGMRESSetup(gmres_solver, par_A,  par_b,  par_x);
		HYPRE_ParCSRGMRESSolve(gmres_solver, par_A,  par_b,  par_x);
		// BoomerAMG
		HYPRE_BoomerAMGSetup(solver, par_A,  par_b,  par_x);
		HYPRE_BoomerAMGSolve(solver, par_A,  par_b,  par_x);
		
	}	
	
	// leaving
	HYPRE_ParCSRGMRESDestroy(gmres_solver);	
	
}

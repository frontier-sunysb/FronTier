/*       
 * File:   solver.h
 * Author: shuqiang (robert) wang 
 *
 * Questions:
 * 1) Should the virtual functions be make pure virtual functions?
 * 2) The function names do not conform to any standard.
 */

#ifndef  LSOLVER_H
#define  LSOLVER_H

#include <stdio.h>
#include <stdlib.h>

class SOLVER
{
public:
	SOLVER(){};
	SOLVER(int ilower, int iupper, int d_nz, int o_nz){};
	virtual ~SOLVER(){};
	virtual void Create(int ilower, int iupper, int d_nz, int o_nz){};
	
	virtual void Set_A(int i, int j, double val){};		// A[i][j]=val;
	virtual void Add_A(int i, int j, double val){};		// A[i][j]=A[i][j]+val;
	virtual void Set_x(int i, double val){};		// x[i]=val;
	virtual void Set_x(double *p){};			// x[i]=p[i];
	virtual void Add_x(int i, double val){};		// x[i]=x[i]+val;
	virtual void Get_x(double *p){};			// get the x from ij_x to p.		
	virtual void Get_x(double *p, int n, int *global_index){};
	virtual void Set_b(int i, double val){};		// b[i]=val;
	virtual void Set_b(double *b){};	
	virtual void Add_b(int i, double val){};		// b[i]=b[i]+val;

	virtual void SetMaxIter(int val){};	
	virtual void GetFinalRelativeResidualNorm(double *rel_resid_norm){};
	virtual void GetNumIterations(int *num_iterations){};
	virtual void GetExtremeSingularValues(double *max, double *min){};

	virtual void Solve(void){};	
	virtual void Solve_withPureNeumann(void){};	
	virtual void Read_A(char *filename){};
	virtual void Print_A(char *filename){};
	virtual void Read_b(char *filename){};
	virtual void Print_b(char *filename){};
	virtual void Read_x(char *filename){};
	virtual void Print_x(char *filename){};
	virtual void test(void){};
	
};

#endif


/*******************************************************************
 * 		solver_lapack.h
 * a wrap of dgesvx.f from lapack for solving Ax=b.
 *******************************************************************/
#ifndef  CLASS_LAPACK_H
#define  CLASS_LAPACK_H

#include <stdio.h>
#include <stdlib.h>
#include "solver.h"

class LAPACK: public SOLVER {
public:
	LAPACK();
	LAPACK(int ilower, int iupper, int d_nz, int o_nz);
	~LAPACK();
	void Create(int ilower, int iupper, int d_nz, int o_nz);
	
	void Set_A(int i, int j, double val);		// A[i][j]=val;
	void Add_A(int i, int j, double val);		// A[i][j]=A[i][j]+val;
	//void Set_x(int i, double val);			// x[i]=val;
	//void Set_x(double *p);					// x[i]=p[i];
	//void Add_x(int i, double val);			// x[i]=x[i]+val;
	void Get_x(double *p);						// get the x from ij_x to p.		
	//void Get_x(double *p, int n, int *global_index);
	void Set_b(int i, double val);				// b[i]=val;
	void Set_b(double *b);	
	void Add_b(int i, double val);				// b[i]=b[i]+val;

	//void SetMaxIter(int val){};	
	void GetFinalRelativeResidualNorm(double *rel_resid_norm);		// this error norm is not a relative residual norm
	//void GetNumIterations(int *num_iterations){};

	void Solve(void);	
	//virtual void Solve_withPureNeumann(void){};	
	//virtual void Read_A(char *filename){};
	void Print_A(char *filename);
	//virtual void Read_b(char *filename){};
	void Print_b(char *filename);
	//virtual void Read_x(char *filename){};
	void Print_x(char *filename);
	//virtual void test(void){};

private:
	// scalar arguments
	char m_equed;		// output: 
	char m_fact;		// input: 'N' or "E"
	char m_trans;		// input: 'N'
	int  m_info;		// output: 0, <0, >0	
	int  m_lda;			// input: n
	int  m_ldaf;		// input: n
	int  m_ldb;			// input: n
	int  m_ldx;			// input: n
	int  m_n;			// input: n
	int  m_nrhs;		// input: 1
	double	m_rcond;	// output: 

	// array arguments
	int  *m_ipiv;		// output: m_ipiv[n]
	int	 *m_iwork;		// workspace: m_iwork[n]
	double	*m_A;		// input/output: m_A[n][n]
	double	*m_Af;		// output: m_Af[n][n]
	double	*m_b;		// input/output: m_b[n]
	double 	*m_berr;	// output: m_berr[1]
	double 	*m_c;		// output: m_c[n]
	double  *m_ferr;	// output: m_ferr[1]
	double	*m_r;		// output: m_r[n]
	double	*m_work;	// workspace/output: m_work[4*n]
	double	*m_x;		// output: m_x[n]

};

#endif  // #ifndef CLASS_SOLVER_H


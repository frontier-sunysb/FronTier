/*****************************************************************
 * 		solver_lapack.c
 * a wrap of dgesvx.f from lapack for solving Ax=b.
 *
 * b should always be a vector.
 ** ***************************************************************/
#include "solver_lapack.h"

// the following FROTRAN_NAME is copied from FronTier
#if defined(cray) || defined(_AIX)
#	define 	FORTRAN_NAME(a)		a
#else
#	define	FORTRAN_NAME(a)		a ## _
#endif

#if defined(__cplusplus)
extern "C" {
#endif
//FORTRAN_NAME(dgesvx)(&m_fact, &m_trans, &m_n, &m_nrhs, m_A, &m_lda, m_Af, &m_ldaf, m_ipiv,
//					 &m_equed, m_r, m_c, m_b, &m_ldb, m_x, &m_ldx, &m_rcond, m_ferr, 
//					 m_berr, m_work, m_iwork, &m_info)		
void FORTRAN_NAME(dgesvx)(char*,char*,int*,int*,double*,int*,double*,int*,int*,
						  char*,double*,double*,double*,int*,double*,int*,double*,double*,
						  double*,double*,int*,int*);
#if defined(__cplusplus)
}
#endif


LAPACK::LAPACK()
{
	m_A = NULL;
}
LAPACK::~LAPACK()
{
	if(m_A!=NULL)
	{
		delete [] m_ipiv;
		delete [] m_iwork;
		delete [] m_A;
		delete [] m_Af;
		delete [] m_b;
		delete [] m_berr;
		delete [] m_c;
		delete [] m_ferr;
		delete [] m_r;
		delete [] m_work;
		delete [] m_x;
	}
}
LAPACK::LAPACK(int ilower, int iupper, int d_nz, int o_nz)
{
	Create(ilower,iupper,d_nz,o_nz);	
}

void LAPACK::Create(int ilower, int iupper, int d_nz, int o_nz)
{
	m_n		= iupper - ilower + 1;
	m_fact 	= 'N';
	m_trans = 'N';
	m_lda	= m_n;
	m_ldaf	= m_n;
	m_ldb	= m_n;
	m_ldx	= m_n;
	m_nrhs	= 1;

	m_ipiv 	= new int[m_n];			// out
	m_iwork	= new int[m_n];			// work
	m_A		= new double[m_n*m_n];	// in
	m_Af	= new double[m_n*m_n];	// out
	m_b		= new double[m_n];		// in
	m_berr	= new double[1];		// out
	m_c		= new double[m_n];		// out
	m_ferr  = new double[1];		// out
	m_r		= new double[m_n];		// out
	m_work	= new double[4*m_n];	// work
	m_x		= new double[m_n];		// out

	int i;
	for(i=0; i<m_n; i++)
	{
		m_b[i] = 0;
		m_x[i] = 0;
	}
	for(i=0; i<m_n*m_n; i++)
		m_A[i] = 0;
}

void LAPACK::Set_A(int i, int j, double val)
{
	m_A[i+j*m_n] = val;
}
void LAPACK::Add_A(int i, int j, double val)
{
	m_A[i+j*m_n] += val;
}
void LAPACK::Get_x(double *x)
{
	int i;
	for(i=0; i<m_n; i++)
		x[i] = m_x[i];
}
void LAPACK::Set_b(int i, double val)
{
	m_b[i] = val;
}
void LAPACK::Set_b(double *b)
{
	int i;
	for(i=0; i<m_n; i++)
		m_b[i] = b[i];
}
void LAPACK::Add_b(int i, double val)
{
	m_b[i] += val;
}
// this error norm is not a relative residual nor
void LAPACK::GetFinalRelativeResidualNorm(double *rel_resid_norm)
{
	if(m_info==0)
		printf("LAPACK::GetFinalRelativeResidualNorm: successful. \n");
	else if(m_info<0)
		printf("LAPACK::GetFinalRelativeResidualNorm: illegal value exists. \n");
	else
		printf("LAPACK::GetFinalRelativeResidualNorm: singular matrix. \n");
}

void LAPACK::Solve(void)
{
	FORTRAN_NAME(dgesvx)(&m_fact, &m_trans, &m_n, &m_nrhs, m_A, &m_lda, m_Af, &m_ldaf, m_ipiv,
						 &m_equed, m_r, m_c, m_b, &m_ldb, m_x, &m_ldx, &m_rcond, m_ferr, 
						 m_berr, m_work, m_iwork, &m_info);
}
void LAPACK::Print_A(char *filename)
{
	int i, j;
	for(j=0; j<m_n; j++)
	for(i=0; i<m_n; i++)
		printf("A[%d][%d] = %f \n", i, j, m_A[i+j*m_n]);
}
void LAPACK::Print_b(char *filename)
{
	int i;
	for(i=0; i<m_n; i++)
		printf("b[%d] = %f \n", i, m_b[i]);
}
void LAPACK::Print_x(char *filename)
{
	int i;
	for(i=0; i<m_n; i++)
		printf("x[%d] = %f \n", i, m_x[i]);
}

/**************************************************************************
 *		vector/matrix operation functions 
 **************************************************************************/
#include "melting.h"
#include "melting_basic.h"

void VectorZero(int n, double *vector)
{
	int i;
	for (i=0; i<n; i++)
		vector[i] = 0;	
}

void VectorZero(int n, int *vector)
{
	int i;
	for (i=0; i<n; i++)
		vector[i] = 0;	
}

void VectorCopy(int n, double *vector1, double *vector2)
{
	int i;
	for (i=0; i<n; i++)
		vector1[i] = vector2[i];
}

void MatrixPrint(int my_rank, char *name, int m, int n, double **matrix)
{
	int i, j;
	printf("\n PID %d: %s(%4d,%4d)", my_rank, name, m, n);
	for (i=0; i<m; i++)
	    for (j=0; j<n; j++)
		printf("\n PID %d: %s(%4d, %4d)=\t%12.4e", 
			my_rank, name, i, j, matrix[i][j]);
	printf("\n");
}

double MatrixMax(int m, int n, double **matrix)
{
	int i, j;
	double max = matrix[0][0];
	for (i=0; i<m; i++)
		for (j=0; j<n; j++)
			if (matrix[i][j]>max)
				max = matrix[i][j];
	return max;
			
}

double MatrixMin(int m, int n, double **matrix)
{
	int i, j;
	double min = matrix[0][0];
	for (i=0; i<m; i++)
		for (j=0; j<n; j++)
			if (matrix[i][j]<min)
				min = matrix[i][j];
	return min;
			
}

void MatrixZero(int m, int n, double **matrix)
{
	int i,j;
	for (i=0; i<m; i++)
		for (j=0; j<n; j++)
			matrix[i][j] = 0;
	
}

void MatrixIdentity(int n, double **matrix)
{
	int i,j;
	for (i=0; i<n; i++)
	{
		for (j=0; j<n; j++)
			matrix[i][j] = 0;
		matrix[i][i] = 1;
	}
}

void MatrixCopy(int m, int n, double **U2, double **U1)
{
	int i,j;
	for (i=0; i<m; i++)
		for (j=0; j<n; j++)
			U2[i][j] = U1[i][j];	
}

void ArrayZero(int m, int n, int l, double ***matrix)
{
	int i, j, k;
	for (i=0; i<m; i++)
		for (j=0; j<n; j++)
			for (k=0; k<l; k++)
				matrix[i][j][k] = 0;
}

void Vector2Matrix(int m, int n, double **matrix, double *vector)
{
	int i, j;
	for (i=0; i<m; i++)
		for (j=0; j<n; j++)
			matrix[i][j] = vector[i*n+j];
}
void ArrayPrint(int my_rank, char*name, int m, int n, int l, double ***array)
{
	int i, j, k;
	printf("\n PID %d: %s(%4d,%4d,%4d)", my_rank, name, m, n, l);
	for (i=0; i<m; i++)
	    for (j=0; j<n; j++)
		for (k=0; k<l; k++)
		    printf("\n PID %d: %s(%4d, %4d, %4d)=\t%12.4e", 
		    		my_rank, name, i, j, k, array[i][j][k]);
}


// C = AB
void MatrixMultiply(int n, double **C, double **A, double **B)	
{
	int i, j, k;
	for (i=0; i<n; i++)
		for (j=0; j<n; j++)
		{
			C[i][j] = 0;
			for (k=0; k<n; k++)
				C[i][j] += A[i][k]*B[k][j];
		}
}
// C = A'A
void MatrixMultiply(int n, double **C, double **A)	
{
	int i, j, k;
	for (i=0; i<n; i++)
		for (j=0; j<n; j++)
		{
			C[i][j] = 0;
			for (k=0; k<n; k++)
				C[i][j] += A[k][i]*A[k][j];
		}
}

// C = AB
void MatrixMultiply(double C[3][3], double A[3][3], double B[3][3])
{
	int i, j, k;
	for (i=0; i<3; i++)
	    for (j=0; j<3; j++)
	    {
		C[i][j] = 0;
		for (k=0; k<3; k++)
		    C[i][j] += A[i][k]*B[k][j];
	    }
}

// c = Ab
void MatrixVectorMultiply(int n, double *c, double **A, double *b)
{
	int i, j;
	for (i = 0; i < n; i++)
	{
	    c[i] = 0;
	    for (j = 0; j < n; j++)
		c[i] += A[i][j]*b[j];
	}
}

// B = inverse(A); 
void SymmetricMatrixInverse(double B[3][3], double A[3][3])
{
	double det = A[0][0]*A[1][1]*A[2][2] + A[0][1]*A[1][2]*A[2][0] + 
		A[0][2]*A[1][0]*A[2][1] - ( A[0][0]*A[1][2]*A[2][1] + 
		A[0][1]*A[1][0]*A[2][2] + A[0][2]*A[1][1]*A[2][0] );
	B[0][0] = 1/det*(A[1][1]*A[2][2] - A[1][2]*A[2][1]);	
	B[0][1] = -1/det*(A[1][0]*A[2][2] - A[1][2]*A[2][0]);	
	B[0][2] = 1/det*(A[1][0]*A[2][1] - A[1][1]*A[2][0]);
	B[1][0] = B[0][1];					
	B[1][1] = 1/det*(A[0][0]*A[2][2] - A[0][2]*A[2][0]);	
	B[1][2] = -1/det*(A[0][0]*A[2][1] - A[0][1]*A[2][0]);
	B[2][0] = B[0][2];					
	B[2][1] = B[1][2];						
	B[2][2] = 1/det*(A[0][0]*A[1][1] - A[0][1]*A[1][0]);
}

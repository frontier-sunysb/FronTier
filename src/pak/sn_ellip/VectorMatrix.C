/**************************************************************************
 *		vector/matrix operation functions 
 **************************************************************************/
#include <stdio.h>
#include <VectorMatrix.h>

void NewVector(int n, int **vector)
{
	*vector = new int[n];
}


void ReNewVector(int n, int &sizevector, int **vector)
{
	if(n>sizevector)
	{	
		if(sizevector>0)
			DeleteVector(sizevector, vector);
		sizevector = n;	
		NewVector(n, vector);		
	}
}

void ReNewVector(int n, int &sizevector, double **vector)
{
	if(n>sizevector)
	{	
		if(sizevector>0)
			DeleteVector(sizevector, vector);
		sizevector = n;	
		NewVector(n, vector);		
	}
}

void NewVector(int n, double **vector)
{
	*vector = new double[n];
}

void NewMatrix(int n, int ***matrix)
{
	int i;
	(*matrix) = new int*[n];
	for (i=0; i<n; i++) 
		(*matrix)[i] = new int[n];
}

void NewMatrix(int n, double ***matrix)
{
	int i;
	(*matrix) = new double*[n];
	for (i=0; i<n; i++) 
		(*matrix)[i] = new double[n];
}

void NewMatrix(int m, int n, int ***matrix)
{
	int i;
	(*matrix) = new int*[m];
	for (i=0; i<m; i++) 
		(*matrix)[i] = new int[n];
}

void NewMatrix(int m, int n, char ***matrix)
{
	int i;
	(*matrix) = new char*[m];
	for (i=0; i<m; i++) 
	{
		(*matrix)[i] = new char[n];
		(*matrix)[i][0] = '\0';
	}
}

void NewMatrix(int m, int n, double ***matrix)
{
	int i;
	(*matrix) = new double*[m];
	for (i=0; i<m; i++) 
		(*matrix)[i] = new double[n];
}

// array[m][n][l]
void NewArray(int m, int n, int l, double ****array)
{
	int i, j;
	(*array) = new double**[m];
	for (i=0; i<m; i++)
	{
		(*array)[i] = new double *[n];
		for(j=0; j<n; j++)
			(*array)[i][j] = new double[l];			
	}
			
}
// array[m][n][l]
void NewArray(int m, int n, int l, int ****array)
{
	int i, j;
	(*array) = new int**[m];
	for (i=0; i<m; i++)
	{
		(*array)[i] = new int *[n];
		for(j=0; j<n; j++)
			(*array)[i][j] = new int[l];			
	}
			
}


void DeleteVector(int n, int **vector)
{
	if(n>0)
		delete [] *vector;
}

void DeleteVector(int n, double **vector)
{
	if(n>0 && *vector!=NULL)
	{
		delete [] *vector;
		*vector = NULL;
	}
}

void DeleteMatrix(int n, int ***matrix)
{
	int i;
	for(i=0; i<n; i++)
		delete [] (*matrix)[i];
	delete [] *matrix;
}

void DeleteMatrix(int n, double ***matrix)
{
	int i;
	for(i=0; i<n; i++)
		delete [] (*matrix)[i];
	delete [] *matrix;
}

void DeleteMatrix(int m, int n, int ***matrix)
{
	int i;
	for(i=0; i<m; i++)
		delete [] (*matrix)[i];
	delete [] *matrix;
}
void DeleteMatrix(int m, int n, char ***matrix)
{
	int i;
	for(i=0; i<m; i++)
		delete [] (*matrix)[i];
	delete [] *matrix;
}
void DeleteMatrix(int m, int n, double ***matrix)
{
	int i;
	for(i=0; i<m; i++)
		delete [] (*matrix)[i];
	delete [] *matrix;
}

void DeleteArray(int m, int n, int l, double ****array)
{
	int i, j;
	for(i=0; i<m; i++)
	{
		for(j=0; j<n; j++)
			delete [] (*array)[i][j];
		delete [] (*array)[i];
	}
	delete [] *array;
}
void DeleteArray(int m, int n, int l, int ****array)
{
	int i, j;
	for(i=0; i<m; i++)
	{
		for(j=0; j<n; j++)
			delete [] (*array)[i][j];
		delete [] (*array)[i];
	}
	delete [] *array;
}

void VectorIntPrint(int my_rank, char *name, int n, int *vector)
{
	int i;
	printf("\n PID %d: %s", my_rank, name);
	for(i=0; i<n; i++)
		printf("\n PID %d: %d", my_rank, vector[i]); 
}

void VectorPrint(int my_rank, char *name, int n, double *vector)
{
	int i;
	printf("\n PID %d: %s(%4d)", my_rank, name, n);
	for(i=0; i<n; i++)
		printf("\n PID %d: %s(%4d)=\t%f", my_rank, name, i, vector[i]); 
	printf("\n");
}

void VectorZero(int n, double *vector)
{
	int i;
	for(i=0; i<n; i++)
		vector[i] = 0;
	
}
void VectorCopy(int n, double *vector1, double *vector2)
{
	int i;
	for(i=0; i<n; i++)
		vector1[i] = vector2[i];
}

void MatrixPrint(int my_rank, char *name, int m, int n, double **matrix)
{
	int i, j;
	printf("\n PID %d: %s(%4d,%4d)", my_rank, name, m, n);
	for(i=0; i<m; i++)
		for(j=0; j<n; j++)
			printf("\n PID %d: %s(%4d, %4d)=\t%12.4e", my_rank, name, i, j, matrix[i][j]);
	printf("\n");
}

double MatrixMax(int m, int n, double **matrix)
{
	int i, j;
	double max = matrix[0][0];
	for(i=0; i<m; i++)
		for(j=0; j<n; j++)
			if(matrix[i][j]>max)
				max = matrix[i][j];
	return max;
			
}

double MatrixMin(int m, int n, double **matrix)
{
	int i, j;
	double min = matrix[0][0];
	for(i=0; i<m; i++)
		for(j=0; j<n; j++)
			if(matrix[i][j]<min)
				min = matrix[i][j];
	return min;
			
}


void MatrixZero(int m, int n, double **matrix)
{
	int i,j;
	for(i=0; i<m; i++)
		for(j=0; j<n; j++)
			matrix[i][j] = 0;
	
}

void MatrixIdentity(int n, double **matrix)
{
	int i,j;
	for(i=0; i<n; i++)
	{
		for(j=0; j<n; j++)
			matrix[i][j] = 0;
		matrix[i][i] = 1;
	}
}

void MatrixCopy(int m, int n, double **U2, double **U1)
{
	int i,j;
	for(i=0; i<m; i++)
		for(j=0; j<n; j++)
			U2[i][j] = U1[i][j];	
}


void Vector2Matrix(int m, int n, double **matrix, double *vector)
{
	int i, j;
	for(i=0; i<m; i++)
		for(j=0; j<n; j++)
			matrix[i][j] = vector[i*n+j];
}


void ArrayPrint(int my_rank, char*name, int m, int n, int l, double ***array)
{
	int i, j, k;
	printf("\n PID %d: %s(%4d,%4d,%4d)", my_rank, name, m, n, l);
	for(i=0; i<m; i++)
		for(j=0; j<n; j++)
			for(k=0; k<l; k++)
				printf("\n PID %d: %s(%4d, %4d, %4d)=\t%12.4e", my_rank, name, i, j, k, array[i][j][k]);
}


// C = AB
void MatrixMultiply(int n, double **C, double **A, double **B)	
{
	int i, j, k;
	for(i=0; i<n; i++)
		for(j=0; j<n; j++)
		{
			C[i][j] = 0;
			for(k=0; k<n; k++)
				C[i][j] += A[i][k]*B[k][j];
		}
}
// C = AB
void MatrixMultiply(double C[3][3], double A[3][3], double B[3][3])
{
	int i, j, k;
	for(i=0; i<3; i++)
		for(j=0; j<3; j++)
		{
			C[i][j] = 0;
			for(k=0; k<3; k++)
				C[i][j] += A[i][k]*B[k][j];
		}
}
// c = Ab
void MatrixVectorMultiply(int n, double *c, double **A, double *b)
{
	int i, j;
	for(i=0; i<n; i++)
		{
			c[i] = 0;
			for(j=0; j<n; j++)
				c[i] += A[i][j]*b[j];
		}
}

// B = inverse(A); 
void SymmetricMatrixInverse(double B[3][3], double A[3][3])
{
	double det = 	 A[0][0]*A[1][1]*A[2][2] + A[0][1]*A[1][2]*A[2][0] + A[0][2]*A[1][0]*A[2][1] 
		     - ( A[0][0]*A[1][2]*A[2][1] + A[0][1]*A[1][0]*A[2][2] + A[0][2]*A[1][1]*A[2][0] );
	B[0][0] = 1/det * (A[1][1]*A[2][2] - A[1][2]*A[2][1]);	B[0][1] = - 1/det * (A[1][0]*A[2][2] - A[1][2]*A[2][0]);	B[0][2] =   1/det * (A[1][0]*A[2][1] - A[1][1]*A[2][0]);
	B[1][0] = B[0][1];					B[1][1] =   1/det * (A[0][0]*A[2][2] - A[0][2]*A[2][0]);	B[1][2] = - 1/det * (A[0][0]*A[2][1] - A[0][1]*A[2][0]);
	B[2][0] = B[0][2];					B[2][1] = B[1][2];						B[2][2] =   1/det * (A[0][0]*A[1][1] - A[0][1]*A[1][0]);
}



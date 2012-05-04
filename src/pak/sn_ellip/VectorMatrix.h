/**************************************************************************
 *		uni_array/bi_array operation functions 
 **************************************************************************/
#ifndef VECTORMATRIX_H
#define VECTORMATRIX_H

void NewVector(int n, int **uni_array);
void ReNewVector(int n, int &sizeuni_array, int **uni_array);		/* */
void ReNewVector(int n, int &sizeuni_array, double **uni_array);		/* */
void NewVector(int n, double **uni_array);
void NewMatrix(int n, int ***bi_array);
void NewMatrix(int n, double ***bi_array);
void NewMatrix(int m, int n, int ***bi_array);
void NewMatrix(int m, int n, char ***bi_array);
void NewMatrix(int m, int n, double ***bi_array);
void NewArray(int m, int n, int l, double ****array);
void NewArray(int m, int n, int l, int ****array);
void DeleteVector(int n, int **uni_array);
void DeleteVector(int n, double **uni_array);
void DeleteMatrix(int n, int ***bi_array);
void DeleteMatrix(int n, double ***bi_array);
void DeleteMatrix(int m, int n, int ***bi_array);
void DeleteMatrix(int m, int n, char ***bi_array);
void DeleteMatrix(int m, int n, double ***bi_array);
void DeleteArray(int m, int n, int l, double ****array);
void DeleteArray(int m, int n, int l, int ****array);

void VectorIntPrint(int my_rank, char *name, int n, int *uni_array);
void VectorPrint(int my_rank, char *name, int n, double *uni_array);
void VectorZero(int n, double *uni_array);
void VectorCopy(int n, double *uni_array1, double *uni_array2);

void MatrixPrint(int my_rank, char *name, int m, int n, double **bi_array);
double MatrixMax(int m, int n, double **bi_array);
double MatrixMin(int m, int n, double **bi_array);
void MatrixZero(int m, int n, double **bi_array);
void MatrixIdentity(int n, double **bi_array);
void MatrixCopy(int m, int n, double **U2, double **U1);

void Vector2Matrix(int m, int n, double **bi_array, double *uni_array);

void ArrayPrint(int my_rank, char*name, int m, int n, int l, double ***array);

void MatrixMultiply(int n, double **C, double **A, double **B);	/* C = AB */
void MatrixMultiply(double C[3][3], double A[3][3], double B[3][3]);	/* C = AB */
void MatrixVectorMultiply(int n, double *C, double **A, double *b);	/* c = Ab */


void SymmetricMatrixInverse(double B[3][3], double A[3][3]);	/* B = inverse(A);  */



#endif	/* #ifndef VECTORMATRIX_H */

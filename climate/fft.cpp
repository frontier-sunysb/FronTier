/***************************************************************
FronTier is a set of libraries that implements differnt types of 
Front Traking algorithms. Front Tracking is a numerical method for 
the solution of partial differential equations whose solutions have 
discontinuities.  

Copyright (C) 1999 by The University at Stony Brook. 

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
****************************************************************/


#include <fftw3.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
/*
*	Perform a n-dimensional FFT
*	For real to complex transform, use dir = 1: 
	    input:  N0xN1x...xN{n-1} real array
	    output: N0xN1x...xN{n-1}/2+1 complex array
	For complex to real transform, use dir = -1:
	    input:  N0xN1x...xN{n-1}/2+1 complex array
	    output:  N0xN1x...xN{n-1} real array
*
*/
extern bool fftnd(
	fftw_complex *in,
	int dim, /*rank of array, can be 1, 2, 3 or any positive integer*/
	int *N,  /*number in each dim*/
	int dir)
{
	fftw_complex *out;
	fftw_plan p;
	int Nc, Nr, i;
	unsigned int flags = FFTW_ESTIMATE; /*or FFTW_MEASURE*/
	static fftw_complex *cplx_array;
	static double *real_array;
	static int N_max = 0;

	/*count total number for input and output*/
	Nc = 1;
	for (i = 0; i < dim-1; i++)
	    Nc *= N[i];
	Nc *= N[dim-1]/2 + 1; /*last dim is cut to (N/2)+1 for complex array*/
	Nr = 1;
	for (i = 0; i < dim; i++)
	    Nr *= N[i];	      /*no cut for real array*/
	
	if (Nr > N_max)
	{
	    N_max = Nr;
            fftw_free(cplx_array);
            fftw_free(real_array);
	    real_array = new double[Nr];
	    cplx_array = new fftw_complex[Nc];
	}

	switch (dir)
	{
	    case 1:
		p = fftw_plan_dft_r2c(dim,N,real_array,cplx_array,flags);
		for (i = 0; i < Nr; i++ )
		    real_array[i] = in[i][0];
		break;
	    case -1:
		p = fftw_plan_dft_c2r(dim,N,cplx_array,real_array,flags);
		for (i = 0; i < Nc; i++ )
		{
		    cplx_array[i][0] = in[i][0]/Nr; /*normalization*/
		    cplx_array[i][1] = in[i][1]/Nr; /*normalization*/
		}
		break;
	    default:
		printf("Dir can only be -1 and 1 in FFT: unknown %d\n",dir);
		break;  
	}	

	fftw_execute(p); /*excute FFT*/

	switch(dir)
	{
	    case 1:
		for (i = 0; i < Nc; i++)
		{
		    in[i][0] = cplx_array[i][0];
		    in[i][1] = cplx_array[i][1];	
		}
		break;
	    case -1:
		for (i = 0; i < Nr; i++)
		{		
    		    in[i][0] = real_array[i];
		    in[i][1] = 0.0;
		}
		break;
	    default:
                printf("Dir can only be -1 and 1 in FFT: unknown %d\n",dir);
                break;
	}
	/*destroy plan*/
	fftw_destroy_plan(p);
}

/*following functions are for verification purpose*/
void printMatrix(const char*,fftw_complex*,int);

int one_dim_test()
{
	const static int N = 100;
	fftw_complex mycomplex[N];
	int dim[1];
	dim[0] = N;
	
	FILE *file;
	int i;

	for (i = 0; i < N; i++)
	{
	    mycomplex[i][0] = 0.7*sin(2*M_PI*50*i/(N-1))+sin(2*M_PI*120*i/(N-1));
	    mycomplex[i][1] = 0.0;
	}

	fftnd(mycomplex,1,dim,1);
	fftnd(mycomplex,1,dim,-1);
	printMatrix("fft1d_test",mycomplex,N);
	return 1;
}

int two_dim_test()
{
	const static int M = 64, N = 64;
	int i,j,index;
	int dim[2];
	dim[0] = M; dim[1] = N;
	FILE *file;
	double wn, L = 1.0;

	fftw_complex myarray[M*N];

	for (j = 0; j < N; j++)
	for (i = 0; i < M; i++)
	{
	    wn = (2*M_PI/L)*sqrt(i*i+j*j);
	    index = j * M + i;
	    myarray[index][0] = sin(2*M_PI*i/M)*cos(2*M_PI*j/N); 
	    myarray[index][1] = 0.0;
	}
	fftnd(myarray,2,dim,1);
	fftnd(myarray,2,dim,-1);
	printMatrix("fft2d_test",myarray,M*N);
        return 1;
}

int three_dim_test()
{
	const static int Nx = 32, Ny = 64, Nz = 128;
	int i,j,k,index;
	FILE *file;
	double wn, L = 1.0;
	int dim[3];
	dim[0] = Nx; dim[1] = Ny; dim[2] = Nz;

	fftw_complex myarray[Nx*Ny*Nz];

	for (k = 0; k < Nz; k++)
	for (j = 0; j < Ny; j++)
	for (i = 0; i < Nx; i++)
	{
	    wn = (2*M_PI/L)*sqrt(i*i+j*j+k*k);
	    index = Nx*(Ny * k + j) + i;
	    myarray[index][0] = cos(2*M_PI*i/Nx)*cos(2*M_PI*j/Ny)*sin(2*M_PI*k/Nz); 
	    myarray[index][1] = 0.0;
	}
	fftnd(myarray,3,dim,1);
	fftnd(myarray,3,dim,-1);
	printMatrix("fft3d_test",myarray,Nx*Ny*Nz);
        return 1;
}

int two_dim_filter()
{
	const static int M = 64, N = 64;
	int i,j,index;
	FILE *file;
	double wn, phi, L = 1.0;

	fftw_complex myarray[M*N];
	int dim[2];
	dim [0] = M; dim[1] = N;

	srand(time(NULL));

	for (i = 0; i < M; i++)
	for (j = 0; j < N/2+1; j++)
	{
	    index = i * (N/2+1) + j;
	    if (i * i + j *j > 4)
	    {
		myarray[index][0] = myarray[index][1] = 0.0;
		continue;
	    }
	    wn = (2*M_PI/L)*sqrt(i*i+j*j);
	    phi  = (double)rand() / (RAND_MAX + 1.0);
	    myarray[index][0] = wn*wn*exp(-wn*wn/pow(2*M_PI*4.7568/L,2))
				 * cos(2*M_PI*phi); 
	    myarray[index][1] = wn*wn*exp(-wn*wn/pow(2*M_PI*4.7568/L,2))
				 * sin(2*M_PI*phi);
	}
	fftnd(myarray,2,dim,-1);
	printMatrix("fft2d_filter_test",myarray,M*N);
        return 1;
}

void printMatrix(const char* filename,fftw_complex* complex_array,int size)
{
	FILE* file;
	int i;
	file = fopen(filename,"w");
	for (i = 0; i < size; i++)
	{
	   fprintf(file,"%f %f\n",complex_array[i][0],complex_array[i][1]);
	}
	fclose(file);
}


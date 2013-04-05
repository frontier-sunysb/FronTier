#include <cuda.h>
#include <stdio.h>
#include "airfoil_gpu.cuh"

__global__ void kernel_spring(SPRING_VERTEX*,double*,double*,double**,double**,
			double**,double*,double*,int*,int);

extern void gpu_spring_solver(
	SPRING_VERTEX *sv,
	double **x_pos,
	double **v_pos,
	int size)
{
	static SPRING_VERTEX *dev_sv;
	static double *dev_x_store,*dev_v_store;
	static double *dev_k,*dev_len0;
	static int *dev_ix;
	static double **dev_x_pos,**dev_v_pos;
	static double **dev_x_nb_store;
	int i,j,n;
	double *x_test;
	int total_num_nb;
	double *data;
	int *ix;
	printf("Entering gpu_spring_solver()\n");

	x_test = (double*) malloc(3*size*sizeof(double));
	cudaMalloc((void**)&dev_sv,size*sizeof(SPRING_VERTEX));
	cudaMalloc((void**)&dev_x_store,3*size*sizeof(double));
	cudaMalloc((void**)&dev_x_pos,size*sizeof(double*));
	cudaMalloc((void**)&dev_v_pos,size*sizeof(double*));

	cudaMemcpy(dev_sv,sv,size*sizeof(SPRING_VERTEX),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_x_store,*x_pos,3*size*sizeof(double),
				cudaMemcpyHostToDevice);
	cudaMemcpy(dev_v_store,*v_pos,3*size*sizeof(double),
				cudaMemcpyHostToDevice);
	cudaMemcpy(x_test,dev_x_store,3*size*sizeof(double),
				cudaMemcpyDeviceToHost);
	total_num_nb = 0;
	for (i = 0; i < size; i++)
	    total_num_nb += sv[i].num_nb;
	data = (double*) malloc(total_num_nb*sizeof(double));
	ix = (int*) malloc(total_num_nb*sizeof(int));
	n = 0;
	for (i = 0; i < size; i++)
	for (j = 0; j < sv[i].num_nb; j++)
	{
	    data[n++] = sv[i].k[j];
	}
	cudaMalloc((void**)&dev_k,total_num_nb*sizeof(double));
	cudaMalloc((void**)&dev_x_nb_store,total_num_nb*sizeof(double*));
	cudaMemcpy(dev_k,data,total_num_nb*sizeof(double),
				cudaMemcpyHostToDevice);
	n = 0;
	for (i = 0; i < size; i++)
	for (j = 0; j < sv[i].num_nb; j++)
	{
	    data[n++] = sv[i].len0[j];
	}
	cudaMalloc((void**)&dev_len0,total_num_nb*sizeof(double));
	cudaMemcpy(dev_len0,data,total_num_nb*sizeof(double),
				cudaMemcpyHostToDevice);
	cudaMalloc((void**)&dev_ix,total_num_nb*sizeof(int));
	n = 0;
	for (i = 0; i < size; i++)
	for (j = 0; j < sv[i].num_nb; j++)
	{
	    ix[n++] = sv[i].ix_nb[j];
	}
	cudaMalloc((void**)&dev_ix,total_num_nb*sizeof(int));
	cudaMemcpy(dev_ix,ix,total_num_nb*sizeof(int),
				cudaMemcpyHostToDevice);
	kernel_spring<<<1,1>>>(dev_sv,dev_x_store,dev_v_store,dev_x_pos,
			dev_v_pos,dev_x_nb_store,dev_k,dev_len0,dev_ix,size);

}	/* end gpu_spring_solver */

__global__ void kernel_spring(
	SPRING_VERTEX *sv,
	double *dev_x_store,
	double *dev_v_store,
	double **dev_x_pos,
	double **dev_v_pos,
	double **dev_x_nb_store,
	double *dev_k,
	double *dev_len0,
	int *dev_ix,
	int size)
{
	int i,j,n;
	n = 0;
	for (i = 0; i < size; ++i)
	{
	    sv[i].k = dev_k + n; 
	    sv[i].len0 = dev_len0 + n; 
	    sv[i].ix_nb = dev_ix + n; 
	    sv[i].x_nb = dev_x_nb_store + n; 
	    n += sv[i].num_nb;
	}
	for (i = 0; i < size; ++i)
	{
	    sv[i].x = dev_x_store + i*3; 
	    sv[i].v = dev_v_store + i*3; 
	    dev_x_pos[i] = dev_x_store + i*3; 
	    dev_v_pos[i] = dev_v_store + i*3; 
	    for (j = 0; j < sv[i].num_nb; ++j)
	    {
	    	sv[i].x_nb[j] = dev_x_store + sv[i].ix_nb[j]*3; 
	    }
	}
}	/* end kernel_spring */
	

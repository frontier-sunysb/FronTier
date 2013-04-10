#include <cuda.h>
#include <stdio.h>
#include "airfoil_gpu.cuh"

__global__ void preset_kernel_spring(SPRING_VERTEX*,
			double**,double*,double*,int*,int);
__global__ void set_kernel_spring(SPRING_VERTEX*,
			double*,double*,double**,double**,
			double*,double*,double**,double**,
			double*,double*,double**,double**,
			double*,double**,int);
__global__ void pos_to_old(double**,double**,
			double**,double**,int);
__global__ void comp_spring_accel(SPRING_VERTEX*,double**,int);
__device__ void dev_comp_spring_accel(SPRING_VERTEX,double*);
__global__ void RK_1(double**,double**,
			double**,double**,double**,double**,
			double**,double,int);
__global__ void RK_2(double**,double**,
			double**,double**,double**,double**,
			double**,double,int);
__global__ void RK_3(double**,double**,
			double**,double**,double**,double**,
			double**,double,int);
__global__ void RK_4(double**,double**,
			double**,double**,double**,double**,
			double**,double,int);

extern void gpu_spring_solver(
	SPRING_VERTEX *sv,
	double **x_pos,
	double **v_pos,
	int size,
	int n_tan,
	double dt)
{
	static SPRING_VERTEX *dev_sv;
	static double *dev_x_store,*dev_v_store;
	static double *dev_x_old_store,*dev_v_old_store;
	static double *dev_x_new_store,*dev_v_new_store;
	static double *dev_accel_store;
	static double *dev_k,*dev_len0;
	static int *dev_ix;
	static double **dev_x_pos,**dev_v_pos;
	static double **dev_x_old,**dev_v_old;
	static double **dev_x_new,**dev_v_new;
	static double **dev_accel;
	static double **dev_x_nb_store;
	static int first = 1;
	int i,j,n;
	int total_num_nb;
	double *data;
	int *ix;
	int TPB = 256;
	int NB = (size+TPB-1)/TPB;

	static double *x_data;
	printf("Entering gpu_spring_solver()\n");

	if (first)
	{
	    first = 0;
	    cudaMalloc((void**)&dev_sv,size*sizeof(SPRING_VERTEX));
	    cudaMalloc((void**)&dev_x_store,3*size*sizeof(double));
	    cudaMalloc((void**)&dev_v_store,3*size*sizeof(double));
	    cudaMalloc((void**)&dev_x_pos,size*sizeof(double*));
	    cudaMalloc((void**)&dev_v_pos,size*sizeof(double*));
	    /* The folllowing are for internally used data */
	    cudaMalloc((void**)&dev_x_old_store,3*size*sizeof(double));
	    cudaMalloc((void**)&dev_v_old_store,3*size*sizeof(double));
	    cudaMalloc((void**)&dev_x_old,size*sizeof(double*));
	    cudaMalloc((void**)&dev_v_old,size*sizeof(double*));
	    cudaMalloc((void**)&dev_x_new_store,3*size*sizeof(double));
	    cudaMalloc((void**)&dev_v_new_store,3*size*sizeof(double));
	    cudaMalloc((void**)&dev_x_new,size*sizeof(double*));
	    cudaMalloc((void**)&dev_v_new,size*sizeof(double*));
	    cudaMalloc((void**)&dev_accel_store,3*size*sizeof(double));
	    cudaMalloc((void**)&dev_accel,size*sizeof(double*));

	    total_num_nb = 0;
	    for (i = 0; i < size; i++)
	    	total_num_nb += sv[i].num_nb;

	    data = (double*)malloc(total_num_nb*sizeof(double));
	    ix = (int*)malloc(total_num_nb*sizeof(int));
	    x_data = (double*) malloc(3*size*sizeof(double));

	    /* This call will copy all direct elements of sv */
	    /* including num_nb,m,lambda. */
	    cudaMemcpy(dev_sv,sv,size*sizeof(SPRING_VERTEX),
				cudaMemcpyHostToDevice);

	    n = 0;
	    for (i = 0; i < size; i++)
	    for (j = 0; j < sv[i].num_nb; j++)
	    {
	    	data[n++] = sv[i].k[j];
	    }
	    cudaMalloc((void**)&dev_k,total_num_nb*sizeof(double));
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
	    n = 0;
	    for (i = 0; i < size; i++)
	    for (j = 0; j < sv[i].num_nb; j++)
	    {
	    	ix[n++] = sv[i].ix_nb[j];
	    }
	    cudaMalloc((void**)&dev_ix,total_num_nb*sizeof(int));
	    cudaMemcpy(dev_ix,ix,total_num_nb*sizeof(int),
				cudaMemcpyHostToDevice);
	    cudaMalloc((void**)&dev_x_nb_store,total_num_nb*sizeof(double*));
	    /* The following call will set x,v,x_nb,k,len0 for sv */
	    
	    preset_kernel_spring<<<1,1>>>(dev_sv,
			dev_x_nb_store,dev_k,dev_len0,dev_ix,size);

	    set_kernel_spring<<<NB,TPB>>>(dev_sv,
			dev_x_store,dev_v_store,
                        dev_x_pos,dev_v_pos,
                        dev_x_old_store,dev_v_old_store,
			dev_x_old,dev_v_old,
			dev_x_new_store,dev_v_new_store,
			dev_x_new,dev_v_new, 
			dev_accel_store,dev_accel,
			size);

	}

	/* The following two copy calls input data for each time step */

        for (i = 0; i < size; i++)
        for (j = 0; j < 3; j++)
        {
            x_data[i*3+j] = x_pos[i][j];
        }
	cudaMemcpy(dev_x_store,x_data,3*size*sizeof(double),
				cudaMemcpyHostToDevice);
        for (i = 0; i < size; i++)
        for (j = 0; j < 3; j++)
        {
            x_data[i*3+j] = v_pos[i][j];
        }
	cudaMemcpy(dev_v_store,x_data,3*size*sizeof(double),
				cudaMemcpyHostToDevice);
	pos_to_old<<<NB,TPB>>>(dev_x_old,dev_x_pos,
				dev_v_old,dev_v_pos,size);
	comp_spring_accel<<<NB,TPB>>>(dev_sv,dev_accel,size);

	for (n = 0; n < n_tan; ++n)
        {
	    RK_1<<<NB,TPB>>>(
			dev_x_new,dev_v_new,dev_x_pos,dev_v_pos,
			dev_x_old,dev_v_old,dev_accel,dt,size);
	    comp_spring_accel<<<NB,TPB>>>(dev_sv,dev_accel,size);

	    RK_2<<<NB,TPB>>>(
			dev_x_new,dev_v_new,dev_x_pos,dev_v_pos,
			dev_x_old,dev_v_old,dev_accel,dt,size);
	    comp_spring_accel<<<NB,TPB>>>(dev_sv,dev_accel,size);

	    RK_3<<<NB,TPB>>>(
			dev_x_new,dev_v_new,dev_x_pos,dev_v_pos,
			dev_x_old,dev_v_old,dev_accel,dt,size);
	    comp_spring_accel<<<NB,TPB>>>(dev_sv,dev_accel,size);

	    RK_4<<<NB,TPB>>>(
			dev_x_new,dev_v_new,dev_x_pos,dev_v_pos,
			dev_x_old,dev_v_old,dev_accel,dt,size);
	    if (n != n_tan-1)
            {
		comp_spring_accel<<<NB,TPB>>>(dev_sv,dev_accel,size);
		
		pos_to_old<<<NB,TPB>>>(dev_x_old,
			dev_x_pos,dev_v_old,dev_v_pos,size);
            } 
	}

	cudaMemcpy(x_data,dev_x_store,3*size*sizeof(double),
				cudaMemcpyDeviceToHost);
        for (i = 0; i < size; i++)
        for (j = 0; j < 3; j++)
        {
            x_pos[i][j] = x_data[i*3+j];
        }
	cudaMemcpy(x_data,dev_v_store,3*size*sizeof(double),
				cudaMemcpyDeviceToHost);
        for (i = 0; i < size; i++)
        for (j = 0; j < 3; j++)
        {
            v_pos[i][j] = x_data[i*3+j];
        }
}	/* end gpu_spring_solver */

__global__ void preset_kernel_spring(
	SPRING_VERTEX *sv,
	double **dev_x_nb_store,
	double *dev_k,
	double *dev_len0,
	int *dev_ix,
	int size)
{
	int i,n = 0;

	for (i = 0; i < size; ++i)
	{
	    sv[i].k = dev_k + n; 
	    sv[i].len0 = dev_len0 + n; 
	    sv[i].ix_nb = dev_ix + n; 
	    sv[i].x_nb = dev_x_nb_store + n; 
	    n += sv[i].num_nb;
	}
}	/* end preset_kernel_spring */
	
__global__ void set_kernel_spring(
	SPRING_VERTEX *sv,
	double *dev_x_store,
	double *dev_v_store,
	double **dev_x_pos,
	double **dev_v_pos,
	double *dev_x_old_store,
	double *dev_v_old_store,
	double **dev_x_old,
	double **dev_v_old,
	double *dev_x_new_store,
	double *dev_v_new_store,
	double **dev_x_new,
	double **dev_v_new,
	double *dev_accel_store,
	double **dev_accel,
	int size)
{
	int i,j;
	i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i < size)
	{
	    sv[i].x = dev_x_store + i*3; 
	    sv[i].v = dev_v_store + i*3; 
	    dev_x_pos[i] = dev_x_store + i*3; 
	    dev_v_pos[i] = dev_v_store + i*3; 
	    dev_x_old[i] = dev_x_old_store + i*3; 
	    dev_v_old[i] = dev_v_old_store + i*3; 
	    dev_x_new[i] = dev_x_new_store + i*3; 
	    dev_v_new[i] = dev_v_new_store + i*3; 
	    dev_accel[i] = dev_accel_store + i*3; 
	    for (j = 0; j < sv[i].num_nb; ++j)
	    {
	    	sv[i].x_nb[j] = dev_x_store + sv[i].ix_nb[j]*3; 
	    }
	}
}	/* end set_kernel_spring */
	
__global__ void pos_to_old(
	double **dev_x_old,
	double **dev_x_pos,
	double **dev_v_old,
	double **dev_v_pos,
	int size)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i < size)
        {
            for ( int j = 0; j < 3; ++j)
            {
        	dev_x_old[i][j] = dev_x_pos[i][j];
        	dev_v_old[i][j] = dev_v_pos[i][j];
            }
	}
}	/* end pos_to_old */

__global__ void comp_spring_accel(
	SPRING_VERTEX *sv,
	double **dev_accel,
	int size)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < size)
        {
	    dev_comp_spring_accel(sv[i],dev_accel[i]);
	}
}

__device__ void dev_comp_spring_accel(
	SPRING_VERTEX sv,
	double *accel)
{
	int i,k;
	double len,vec[3];
	int dim = 3;
	for (k = 0; k < dim; ++k) 
            accel[k] = 0.0;
        for (i = 0; i < sv.num_nb; ++i)
        {
            len = 0.0;
            for (k = 0; k < dim; ++k)
            {
                vec[k] = sv.x_nb[i][k] - sv.x[k];
                len += vec[k]*vec[k];
            }
            len = sqrt(len);
            for (k = 0; k < dim; ++k)
            {
                vec[k] /= len;
                accel[k] += sv.k[i]*((len - sv.len0[i])*vec[k])/sv.m;
            }
        }
        for (k = 0; k < dim; ++k)
        {
            accel[k] += -sv.lambda*sv.v[k]/sv.m + sv.ext_accel[k];
        }
}	/* end dev_compute_spring_accel */

__global__ void	RK_1(
	double **dev_x_new,
	double **dev_v_new,
	double **dev_x_pos,
	double **dev_v_pos,
	double **dev_x_old,
	double **dev_v_old,
	double **dev_accel,
	double dt,
	int size)
{
	int i,j;
	i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < size)
        {
	    for (j = 0; j < 3; ++j)
            {
                dev_x_new[i][j] = dev_x_old[i][j] + dt*dev_v_old[i][j]/6.0;
                dev_v_new[i][j] = dev_v_old[i][j] + dt*dev_accel[i][j]/6.0;
                dev_x_pos[i][j] = dev_x_old[i][j] + 0.5*dev_v_old[i][j]*dt;
                dev_v_pos[i][j] = dev_v_old[i][j] + 0.5*dev_accel[i][j]*dt;
            }
	}
}

__global__ void	RK_2(
	double **dev_x_new,
	double **dev_v_new,
	double **dev_x_pos,
	double **dev_v_pos,
	double **dev_x_old,
	double **dev_v_old,
	double **dev_accel,
	double dt,
	int size)
{
	int i,j;
	i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < size)
        {
	    for (j = 0; j < 3; ++j)
            {
		dev_x_new[i][j] += dt*dev_v_pos[i][j]/3.0;
		dev_v_new[i][j] += dt*dev_accel[i][j]/3.0;
		dev_x_pos[i][j] = dev_x_old[i][j] + 0.5*dev_v_pos[i][j]*dt;
		dev_v_pos[i][j] = dev_v_old[i][j] + 0.5*dev_accel[i][j]*dt;
            }
	}
}

__global__ void	RK_3(
	double **dev_x_new,
	double **dev_v_new,
	double **dev_x_pos,
	double **dev_v_pos,
	double **dev_x_old,
	double **dev_v_old,
	double **dev_accel,
	double dt,
	int size)
{
	int i,j;
	i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < size)
        {
	    for (j = 0; j < 3; ++j)
            {
		dev_x_new[i][j] += dt*dev_v_pos[i][j]/3.0;
		dev_v_new[i][j] += dt*dev_accel[i][j]/3.0;
		dev_x_pos[i][j] = dev_x_old[i][j] + dev_v_pos[i][j]*dt;
		dev_v_pos[i][j] = dev_v_old[i][j] + dev_accel[i][j]*dt;
            }
	}
}

__global__ void	RK_4(
	double **dev_x_new,
	double **dev_v_new,
	double **dev_x_pos,
	double **dev_v_pos,
	double **dev_x_old,
	double **dev_v_old,
	double **dev_accel,
	double dt,
	int size)
{
	int i,j;
	i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < size)
        {
	    for (j = 0; j < 3; ++j)
            {
		dev_x_new[i][j] += dt*dev_v_pos[i][j]/6.0;
		dev_v_new[i][j] += dt*dev_accel[i][j]/6.0;
		dev_x_pos[i][j] = dev_x_new[i][j];
		dev_v_pos[i][j] = dev_v_new[i][j];
            }
	}
}

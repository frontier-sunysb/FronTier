#include <cuda.h>
#include <stdio.h>
#include "airfoil_gpu.cuh"

__global__ void kernel_1(
	double* dev_x_old,
	double* dev_x_mid,
	double* dev_x_new,
	double* dev_v_old,
	double* dev_v_mid,
	double* dev_v_new,
	double* dev_f_mid,
	double dt,
	int dev_size)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if( i < dev_size)
	{
	    dev_x_new[i] = dev_x_new[i] + dt*dev_v_mid[i]/3.0;	
	    dev_v_new[i] = dev_v_new[i] + dt*dev_f_mid[i]/3.0;
	    dev_x_mid[i] = dev_x_old[i] + dt*dev_v_mid[i];
	    dev_v_mid[i] = dev_v_old[i] + dt*dev_f_mid[i];
	}   	
}

extern void call_gputest(
	double **x_old,
	double **x_mid,
	double **x_new,
	double **v_old,
	double **v_mid,
	double **v_new,
	double **f_mid,
	double dt,
	int size)
{
	int dev_size = size*3;
	int bit_size = dev_size*sizeof(double);
	double *dev_x_old;
	double *dev_x_mid;
	double *dev_x_new;
	double *dev_v_old;
	double *dev_v_mid;
	double *dev_v_new;
	double *dev_f_mid;

	cudaMalloc((void **)&dev_x_old,bit_size);	
	cudaMalloc((void **)&dev_x_mid,bit_size);	
	cudaMalloc((void **)&dev_x_new,bit_size);	
	cudaMalloc((void **)&dev_v_old,bit_size);	
	cudaMalloc((void **)&dev_v_mid,bit_size);	
	cudaMalloc((void **)&dev_v_new,bit_size);	
	cudaMalloc((void **)&dev_f_mid,bit_size);	

	for( int i = 0; i < size; i++)
	{
	    cudaMemcpy(dev_x_old + i*3,x_old[i],3*sizeof(double),cudaMemcpyHostToDevice);
	    cudaMemcpy(dev_x_mid + i*3,x_mid[i],3*sizeof(double),cudaMemcpyHostToDevice);
	    cudaMemcpy(dev_x_new + i*3,x_new[i],3*sizeof(double),cudaMemcpyHostToDevice);
	    cudaMemcpy(dev_v_old + i*3,v_old[i],3*sizeof(double),cudaMemcpyHostToDevice);
	    cudaMemcpy(dev_v_mid + i*3,v_mid[i],3*sizeof(double),cudaMemcpyHostToDevice);
	    cudaMemcpy(dev_v_new + i*3,v_new[i],3*sizeof(double),cudaMemcpyHostToDevice);
	    cudaMemcpy(dev_f_mid + i*3,f_mid[i],3*sizeof(double),cudaMemcpyHostToDevice);
	}

	kernel_1<<<(dev_size + 255)/256, 256>>>(dev_x_old,dev_x_mid,dev_x_new,dev_v_old,
			dev_v_mid,dev_v_new,dev_f_mid,dt,dev_size);

       	for( int i = 0; i < size; i++)
        {
            cudaMemcpy(x_old[i],dev_x_old + i*3,3*sizeof(double),cudaMemcpyDeviceToHost);
            cudaMemcpy(x_mid[i],dev_x_mid + i*3,3*sizeof(double),cudaMemcpyDeviceToHost);
            cudaMemcpy(x_new[i],dev_x_new + i*3,3*sizeof(double),cudaMemcpyDeviceToHost);
            cudaMemcpy(v_old[i],dev_v_old + i*3,3*sizeof(double),cudaMemcpyDeviceToHost);
            cudaMemcpy(v_mid[i],dev_v_mid + i*3,3*sizeof(double),cudaMemcpyDeviceToHost);
            cudaMemcpy(v_new[i],dev_v_new + i*3,3*sizeof(double),cudaMemcpyDeviceToHost);
            cudaMemcpy(f_mid[i],dev_f_mid + i*3,3*sizeof(double),cudaMemcpyDeviceToHost);
	}
	printf("OK\tdev_size:\t%d\n",dev_size);
}

__global__ void cudakernel(float *buf,const long n,const int m)
{
    	int i = threadIdx.x + blockIdx.x*blockDim.x;
    	buf[i] = 1.0f*i/n;
    	for(int j = 0; j < m; j++)
            buf[i] = buf[i]*buf[i] - 0.25f;
}

extern void call_gpu_dummy()
{
	const long n = 1024*1024;
	const int m = 10000;
	float data[n]; 
	float *d_data; 
    
	cudaMalloc(&d_data,n*sizeof(float));
	cudakernel<<<n/256,256>>>(d_data,n,m);
	cudaMemcpy(data,d_data,n*sizeof(float),cudaMemcpyDeviceToHost);
	cudaFree(d_data);
}

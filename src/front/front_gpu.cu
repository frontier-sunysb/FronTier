#include <cuda.h>
#include <stdio.h>
#include "front_gpu.cuh"

__global__ void cudakernel(float *buf, const long n, const int m)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    buf[i] = 1.0f * i / n;
    for(int j = 0; j < m; j++)
        buf[i] = buf[i] * buf[i] - 0.25f;
}

extern void call_gpu_dummy()
{
        const long n = 1024 * 1024;
        const int m = 10000;
        float data[n];
        float *d_data;

        cudaMalloc(&d_data, n * sizeof(float));
        cudakernel<<<n/256, 256>>>(d_data, n, m);
        cudaMemcpy(data, d_data, n * sizeof(float), cudaMemcpyDeviceToHost);
        cudaFree(d_data);

}


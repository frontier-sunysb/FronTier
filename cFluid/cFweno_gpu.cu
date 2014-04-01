#include <cuda.h>
#include <stdio.h>
#include <cFluid.h>
#include "cFweno_gpu.cuh"

#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))

static void HandleError( cudaError_t err, const char *file, int line )
{
    if (err != cudaSuccess) {
        printf( "%s in %s at line %d\n",cudaGetErrorString(err),file,line);
        exit( EXIT_FAILURE );
    }
}

__global__ void link_kernel(double**,double**,double*,double*,double*,double*,
		double*,double*,double*,double*,double*,double*,double*);
__global__ void flux_kernel(double,double,double,double,double,int,int,
		double**,double**);
__device__ void u2f(double*,int,int);
__device__ void matmvec(double *b, double RL[5][5], double *x);
__device__ double weno5_scal(double *f);
__device__ void f2is(double *f, double *s);

extern void weno5_get_flux_gpu(
	double gamma, 
	double lambda, 
	int extend_size, 
	int ghost_size, 
	double** u_old, 
	double** flux)
{
	double **dev_u, **dev_flux;
	double *dev_u_dens, *dev_u_engy, *dev_u_pres;
	double *dev_u_momn_1, *dev_u_momn_2, *dev_u_momn_3;
	double *dev_flux_dens, *dev_flux_engy;
	double *dev_flux_momn_1, *dev_flux_momn_2, *dev_flux_momn_3;

	double a,v, maxeig[5];
	int block_size = 512;
	int shared_mem_size_Bytes=11*(block_size+2*ghost_size-1)*sizeof(double);

//startClock("Alloc_mem");
	HANDLE_ERROR(cudaMalloc((void**)&dev_u,6*sizeof(double*)));	
	HANDLE_ERROR(cudaMalloc((void**)&dev_flux,5*sizeof(double*)));	

	HANDLE_ERROR(cudaMalloc((void**)&dev_u_dens,
				extend_size*sizeof(double)));
	HANDLE_ERROR(cudaMalloc((void**)&dev_u_engy,
				extend_size*sizeof(double)));
	HANDLE_ERROR(cudaMalloc((void**)&dev_u_pres,
				extend_size*sizeof(double)));
	HANDLE_ERROR(cudaMalloc((void**)&dev_u_momn_1,
				extend_size*sizeof(double)));
	HANDLE_ERROR(cudaMalloc((void**)&dev_u_momn_2,
				extend_size*sizeof(double)));
	HANDLE_ERROR(cudaMalloc((void**)&dev_u_momn_3,
				extend_size*sizeof(double)));

	HANDLE_ERROR(cudaMalloc((void**)&dev_flux_dens,
				extend_size*sizeof(double)));
	HANDLE_ERROR(cudaMalloc((void**)&dev_flux_engy,
				extend_size*sizeof(double)));
	HANDLE_ERROR(cudaMalloc((void**)&dev_flux_momn_1,
				extend_size*sizeof(double)));
	HANDLE_ERROR(cudaMalloc((void**)&dev_flux_momn_2,
				extend_size*sizeof(double)));
	HANDLE_ERROR(cudaMalloc((void**)&dev_flux_momn_3,
				extend_size*sizeof(double)));
//stopClock("Alloc_mem");

//startClock("link_pointers");
	link_kernel<<<1,1>>>(dev_u, dev_flux, 
		dev_u_dens, dev_u_momn_1, dev_u_momn_2, dev_u_momn_3, 
		dev_u_engy, dev_u_pres, 
		dev_flux_dens, dev_flux_momn_1, dev_flux_momn_2, 
		dev_flux_momn_3, dev_flux_engy);

	HANDLE_ERROR(cudaDeviceSynchronize());
//stopClock("link_pointers");

//startClock("memcpy_H2D");
	HANDLE_ERROR(cudaMemcpy(dev_u_dens,u_old[0],extend_size*
				sizeof(double),cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(dev_u_momn_1,u_old[1],extend_size*
				sizeof(double),cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(dev_u_momn_2,u_old[2],extend_size*
				sizeof(double),cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(dev_u_momn_3,u_old[3],extend_size*
				sizeof(double),cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(dev_u_engy,u_old[4],extend_size*
				sizeof(double),cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(dev_u_pres,u_old[5],extend_size*
				sizeof(double),cudaMemcpyHostToDevice));
//stopClock("memcpy_H2D");
	
//startClock("gpu_compute");
        for(int i = 0; i < extend_size; ++i)
        {
            a = sqrt(gamma * u_old[5][i]/u_old[0][i]);
            v = u_old[1][i]/u_old[0][i];
            maxeig[0] = std::max(maxeig[0], fabs(v - a));
            maxeig[1] = std::max(maxeig[1], fabs(v));
            maxeig[4] = std::max(maxeig[4], fabs(v + a));
        }

	flux_kernel<<<(extend_size-2*ghost_size+1+block_size-1)/block_size,
		block_size,shared_mem_size_Bytes>>>(gamma, lambda, maxeig[0], 
		maxeig[1], maxeig[4],extend_size, ghost_size, dev_u, dev_flux);

	HANDLE_ERROR(cudaDeviceSynchronize());
//stopClock("gpu_compute");

//startClock("memcpy_D2H");
	HANDLE_ERROR(cudaMemcpy(flux[0],dev_flux_dens,extend_size*
				sizeof(double),cudaMemcpyDeviceToHost));
	HANDLE_ERROR(cudaMemcpy(flux[1],dev_flux_momn_1,extend_size*
				sizeof(double),cudaMemcpyDeviceToHost));
	HANDLE_ERROR(cudaMemcpy(flux[2],dev_flux_momn_2,extend_size*
				sizeof(double),cudaMemcpyDeviceToHost));
	HANDLE_ERROR(cudaMemcpy(flux[3],dev_flux_momn_3,extend_size*
				sizeof(double),cudaMemcpyDeviceToHost));
	HANDLE_ERROR(cudaMemcpy(flux[4],dev_flux_engy,extend_size*
				sizeof(double),cudaMemcpyDeviceToHost));
//stopClock("memcpy_D2H");

//startClock("Free_mem");
	HANDLE_ERROR(cudaFree(dev_u));	
	HANDLE_ERROR(cudaFree(dev_flux));	
	
	HANDLE_ERROR(cudaFree(dev_u_dens));	
	HANDLE_ERROR(cudaFree(dev_u_momn_1));	
	HANDLE_ERROR(cudaFree(dev_u_momn_2));	
	HANDLE_ERROR(cudaFree(dev_u_momn_3));	
	HANDLE_ERROR(cudaFree(dev_u_engy));	
	HANDLE_ERROR(cudaFree(dev_u_pres));	

	HANDLE_ERROR(cudaFree(dev_flux_dens));	
	HANDLE_ERROR(cudaFree(dev_flux_momn_1));	
	HANDLE_ERROR(cudaFree(dev_flux_momn_2));	
	HANDLE_ERROR(cudaFree(dev_flux_momn_3));	
	HANDLE_ERROR(cudaFree(dev_flux_engy));	
//stopClock("Free_mem");
}

__global__ void link_kernel( 
	double **dev_u, 
	double **dev_flux, 
	double *dev_u_dens, 
	double *dev_u_momn_1, 
	double *dev_u_momn_2, 
	double *dev_u_momn_3, 
	double *dev_u_engy, 
	double *dev_u_pres, 
	double *dev_flux_dens, 
	double *dev_flux_momn_1, 
	double *dev_flux_momn_2, 
	double *dev_flux_momn_3, 
	double *dev_flux_engy)
{
	dev_u[0] = dev_u_dens;	
	dev_u[1] = dev_u_momn_1;	
	dev_u[2] = dev_u_momn_2;	
	dev_u[3] = dev_u_momn_3;	
	dev_u[4] = dev_u_engy;	
	dev_u[5] = dev_u_pres;	

	dev_flux[0] = dev_flux_dens;	
	dev_flux[1] = dev_flux_momn_1;	
	dev_flux[2] = dev_flux_momn_2;	
	dev_flux[3] = dev_flux_momn_3;	
	dev_flux[4] = dev_flux_engy;	
}

__global__ void	flux_kernel(
	double gamma, 
	double lambda, 
	double dev_maxe_0,
	double dev_maxe_1,
	double dev_maxe_4,
	int extend_size, 
	int ghost_size, 
	double **dev_u, 
	double **dev_flux)
{
	extern __shared__ double uf[];

	int gindex = threadIdx.x + blockIdx.x * blockDim.x + ghost_size;
	int lindex = threadIdx.x + ghost_size;
	int block_extend_size = blockDim.x + 2 * ghost_size - 1;

	double h, gm = gamma-1.0;
	double u_mid[11];
	double RL[5][5];
	double sten_u[6][5], sten_f[6][5];
	double u[5];
	double maxeig[5] = {dev_maxe_0,dev_maxe_1,dev_maxe_1,dev_maxe_1,
				dev_maxe_4};
	double gfluxp[5][5], gfluxm[5][5];
	double f_tmp[5];
	double ff[5];

	if(gindex < extend_size-ghost_size+1 )
	{
	    for(int i = 0; i< 6; i++)
	    {
		uf[lindex + i*block_extend_size] = dev_u[i][gindex];
	    }

	    u2f(uf,lindex,block_extend_size);

	    if (threadIdx.x < ghost_size)
	    {
		for(int i = 0; i< 6; i++)
		{
		    uf[threadIdx.x + i*block_extend_size] = 
				dev_u[i][gindex - ghost_size];
		}

		u2f(uf,threadIdx.x,block_extend_size);
	    }

	    if (lindex > blockDim.x || gindex + 2*ghost_size - 1 > extend_size)
	    {
		for(int i = 0; i< 6; i++)
		{
		    uf[lindex+ghost_size-1+i*block_extend_size]=
				dev_u[i][gindex+ghost_size-1];
		}

		u2f(uf,lindex+ghost_size-1,block_extend_size);
	    }

	    __syncthreads();

	    for(int i = 0; i < 5; i++)
	    {
		u_mid[i] = 0.5*(uf[lindex + i*block_extend_size] + 
				uf[lindex-1 + i*block_extend_size]);
	    }

	    u_mid[5] = u_mid[1]/u_mid[0];
	    u_mid[6] = u_mid[2]/u_mid[0];
	    u_mid[7] = u_mid[3]/u_mid[0];
	    u_mid[9] = sqr(u_mid[5]) + sqr(u_mid[6]) + sqr(u_mid[7]);
	    u_mid[8] = (u_mid[4] - 0.5*u_mid[0]*u_mid[9])*(gamma - 1.0);
	    u_mid[10] = sqrt(gamma*u_mid[8]/u_mid[0]);

	    h = 0.5*u_mid[9] + gamma*u_mid[8]/(gamma - 1.0)/u_mid[0];

            RL[0][0] = h + u_mid[10]*(u_mid[5] - u_mid[10]) / gm;
            RL[0][1] = -1.0*(u_mid[5] + u_mid[10] / gm);
            RL[0][2] = -1.0*u_mid[6];
            RL[0][3] = -1.0*u_mid[7];
            RL[0][4] = 1.0;

            RL[1][0] = -2.0*h + 4.0*sqr(u_mid[10])/gm;
            RL[1][1] = 2.0*u_mid[5];
            RL[1][2] = 2.0*u_mid[6];
            RL[1][3] = 2.0*u_mid[7];
            RL[1][4] = -2.0;

            RL[2][0] = -2.0*u_mid[6]*sqr(u_mid[10])/gm;
            RL[2][1] = 0;
            RL[2][2] = 2.0*sqr(u_mid[10])/gm;
            RL[2][3] = 0;
            RL[2][4] = 0;

            RL[3][0] = -2.0*u_mid[7]*sqr(u_mid[10])/gm;
            RL[3][1] = 0;
            RL[3][2] = 0;
            RL[3][3] = 2.0*sqr(u_mid[10])/gm;
            RL[3][4] = 0;

            RL[4][0] = h - u_mid[10]*(u_mid[5] + u_mid[10])/gm;
            RL[4][1] = -1.0*u_mid[5] + u_mid[10]/gm;
            RL[4][2] = -1.0*u_mid[6];
            RL[4][3] = -1.0*u_mid[7];
            RL[4][4] = 1.0;

            for(int i = 0; i < 5; ++i)
            for(int j = 0; j < 5; ++j)
	    {
                RL[i][j] *= gm/(2.0*sqr(u_mid[10]));
	    }

            for(int i = 0; i < 6; ++i)
            {
                for(int j = 0; j < 5; ++j)
                    u[j] = uf[lindex - ghost_size + i + j*block_extend_size];
                matmvec(sten_u[i],RL,u);
                for(int j = 0; j < 5; ++j)
                    u[j] = uf[lindex-ghost_size + i + (6+j)*block_extend_size];
                matmvec(sten_f[i],RL,u);
            }

            for(int j = 0; j < 5; ++j)
            for(int k = 0; k < 5; ++k)
            {
                gfluxp[j][k] = 0.5*(sten_f[k][j] + maxeig[j]*sten_u[k][j]);
                gfluxm[j][k] = 0.5*(sten_f[5 - k][j] - maxeig[j]*
                                sten_u[5 - k][j]);
            }

            for(int j = 0; j < 5; ++j)
            {
                f_tmp[j] = weno5_scal(gfluxp[j]);
                f_tmp[j] += weno5_scal(gfluxm[j]);
            }
	
	    RL[0][0] = 1.0;
            RL[0][1] = 1.0;
            RL[0][2] = 0.0;
            RL[0][3] = 0.0;
            RL[0][4] = 1.0;

            RL[1][0] = u_mid[5] - u_mid[10];
            RL[1][1] = u_mid[5];
            RL[1][2] = 0.0;
            RL[1][3] = 0.0;
            RL[1][4] = u_mid[5] + u_mid[10];

            RL[2][0] = u_mid[6];
            RL[2][1] = u_mid[6];
            RL[2][2] = 1.0;
            RL[2][3] = 0.0;
            RL[2][4] = u_mid[6];

            RL[3][0] = u_mid[7];
            RL[3][1] = u_mid[7];
            RL[3][2] = 0.0;
            RL[3][3] = 1.0;
            RL[3][4] = u_mid[7];

            RL[4][0] = h - u_mid[5]*u_mid[10];
            RL[4][1] = 0.5*u_mid[9];
            RL[4][2] = u_mid[6];
            RL[4][3] = u_mid[7];
            RL[4][4] = h + u_mid[5]*u_mid[10];;

	    matmvec(ff, RL, f_tmp);
	
	    for(int j = 0; j < 5; ++j)
            {
	        dev_flux[j][gindex] = ff[j];
		if (isnan(dev_flux[j][gindex]))
                {
		    printf("extend_size:\t%d\n",extend_size);
		    printf("In weno5_flux(): f[%d][%d] = %f\n",j,gindex,
                                        dev_flux[j][gindex]);
		    for( int k = 0; k < 11; ++k)
                        printf("umid[%d] = %f\n",k,u_mid[k]);
                    printf("\n");
		}
	    }	
	}
}

__device__ void	u2f(
	double *uf, 
	int lindex, 
	int extend_size)
{
	double v = uf[lindex + 1*extend_size] / uf[lindex];
	uf[lindex + 6*extend_size] = uf[lindex + 1*extend_size];
	uf[lindex + 7*extend_size] = v * uf[lindex + 1*extend_size] + 
				uf[lindex + 5*extend_size];
	uf[lindex + 8*extend_size] = v * uf[lindex + 2*extend_size];
	uf[lindex + 9*extend_size] = v * uf[lindex + 3*extend_size];
	uf[lindex + 10*extend_size] = v * (uf[lindex + 4*extend_size] + 
				uf[lindex + 5*extend_size]);
}

__device__ void matmvec(
	double *b, 
	double L[5][5], 
	double *x)
{
        for(int i = 0; i < 5; ++i)
        {
            b[i] = 0.0;
            for(int j = 0; j < 5; ++j)
            {
                b[i] += L[i][j] * x[j];
            }
        }
}

__device__ double weno5_scal(double *f)
{
        int i, j;
        const double eps = 1.e-16;
        const int p = 2;

        double c[3] = {0.1, 0.6, 0.3}; //*** Optimal weights C_k 
        double is[3]; //*** a smoothness measurement of the flux function
        double alpha[3];
        double omega[3]; // weights for stencils
        double sum;
        double a[3][5] = {{1.0/3, -7.0/6, 11.0/6, 0, 0},
                          {0, -1.0/6, 5.0/6, 1.0/3, 0},
                          {0, 0, 1.0/3, 5.0/6, -1.0/6}};
                //*** coefficients for 2nd-order ENO interpolation stencil
        double w[5]; //weight for every point
        double flux = 0.0;

        f2is(f,is);

        sum = 0.0;

        for(i = 0; i < 3; ++i)
        {
            alpha[i] = c[i]/pow(eps + is[i],p);
            sum += alpha[i];
        }
        for(i = 0; i < 3; ++i)
        {
            omega[i] = alpha[i] / sum;
        }

        for(i = 0; i < 5; ++i)
        {
            w[i] = 0.0;
            for(j = 0; j < 3; ++j)
            {
                w[i] += omega[j] * a[j][i];
            }
        }
        for(i = 0; i < 5; ++i)
        {
            flux += w[i] * f[i];
        }
        return flux;
}

__device__ void f2is(
        double *f,
        double *s)
{
        s[0] = 13.0/12*sqr((f[0] - 2.0*f[1] + f[2])) +
                0.25*sqr((f[0] - 4.0*f[1] + 3.0*f[2]));
        s[1] = 13.0/12*sqr((f[1] - 2.0*f[2] + f[3])) +
                0.25*sqr((f[1] - f[3]));
        s[2] = 13.0/12*sqr((f[2] - 2.0*f[3] + f[4])) +
                0.25*sqr((3.0*f[2] - 4.0*f[3] + f[4]));
}

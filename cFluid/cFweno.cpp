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

#include <cFluid.h>
#include <ctime>
#include <time.h>

#if defined(__GPU__)
#include "cFweno_gpu.cuh"
#endif

static double weno5_scal(double *f);
static void matmvec(double *b, double L[5][5], double *x);
static void f2is(double *f, double *s);
static void u2f(double *u, double *f);
static void weno5_get_flux(POINTER,int,int,double**,double**);
static void arti_compression(POINTER,double*,double*,double,double,double*,int,double &c);

extern void WENO_flux(
        POINTER params,
        SWEEP *vst,
        FSWEEP *vflux,
        int n)
{
	static double **u_old;
	static double **flux;
	int i,extend_size,ghost_size;
	SCHEME_PARAMS *scheme_params = (SCHEME_PARAMS*)params;
	double lambda = scheme_params->lambda;

	if (u_old == NULL)
	{
	    FT_VectorMemoryAlloc((POINTER*)&u_old,6,sizeof(double*));
	    FT_VectorMemoryAlloc((POINTER*)&flux,5,sizeof(double*));
	}
	u_old[0] = vst->dens;
	u_old[1] = vst->momn[0];
	u_old[2] = vst->momn[1];
	u_old[3] = vst->momn[2];
	u_old[4] = vst->engy;
	u_old[5] = vst->pres;
	flux[0] = vflux->dens_flux;
	flux[1] = vflux->momn_flux[0];
	flux[2] = vflux->momn_flux[1];
	flux[3] = vflux->momn_flux[2];
	flux[4] = vflux->engy_flux;
	ghost_size = 3;
	extend_size = n + 2*ghost_size;

#if defined(__GPU__)
	//startClock("Total_time_flux_gpu");
	weno5_get_flux_gpu(scheme_params->gamma,scheme_params->lambda,
				extend_size,ghost_size,u_old,flux);
	//stopClock("Total_time_flux_gpu");
#else
	weno5_get_flux(params,extend_size,ghost_size,u_old,flux);
#endif
	for (i = ghost_size; i < n+ghost_size; ++i)
	{
	    vflux->dens_flux[i] = -lambda*(flux[0][i+1] - flux[0][i]);
	    vflux->momn_flux[0][i] = -lambda*(flux[1][i+1] - flux[1][i]);
	    vflux->momn_flux[1][i] = -lambda*(flux[2][i+1] - flux[2][i]);
	    vflux->momn_flux[2][i] = -lambda*(flux[3][i+1] - flux[3][i]);
	    vflux->engy_flux[i] = -lambda*(flux[4][i+1] - flux[4][i]);
	    if (isnan(vflux->dens_flux[i]))
	    {
		int j;
		for (j = 0; j < extend_size; ++j)
		    printf("u = %f %f %f %f %f\n",
			u_old[0][j],u_old[1][j],u_old[2][j],
			u_old[3][j],u_old[4][j]);
		clean_up(ERROR);
	    }
	}
}	/* end weno5_flux */

static void weno5_get_flux(
	POINTER params,
	int extend_size, 
	int ghost_size, 
	double **u_old,
        double **flux)
{
    	int i, j, k;
    	double u_mid[11];	// rho,mx,my,mz,e,u,v,w,p,u^2+v^2+w^2,a
    	double h; 	//total enthalpy h=0.5*(u^2+v^2+w^2) + (a^2)/(gma - 1)
    	double R[5][5], L[5][5];
    	double sten_u[6][5], sten_f[6][5]; //f_-2, f_-1, f_0, f_1, f_2, f_3
    	double gfluxp[5][5], gfluxm[5][5]; 
    	double maxeig[5] = {0 ,0 ,0, 0, 0};
	double v,u[9],f_tmp[5],ff[5];
	static double **f;
	static int max_n = 0;
	SCHEME_PARAMS *scheme_params = (SCHEME_PARAMS*)params;
	double gamma = scheme_params->gamma;
	double vecp[5][4],vecm[5][4];
    	double a,gm = gamma - 1.0;
	double c,f_ac_tmp,ac_alpha,ac_alpha_nume,ac_alpha_deno;
	double f_prevp[5],f_nowp[5],f_prevm[5],f_nowm[5];
	double gflux_tmp[5];

	if (max_n < extend_size)
        {
            if (max_n != 0)
                FT_FreeThese(1,f);
            max_n = extend_size;
            FT_MatrixMemoryAlloc((POINTER*)&f,extend_size,5,sizeof(double));
        }

    	for(i = 0; i < extend_size; ++i)
    	{
	    a = sqrt(gamma * u_old[5][i]/u_old[0][i]);
	    v = u_old[1][i]/u_old[0][i];
	    maxeig[0] = std::max(maxeig[0], fabs(v - a));
	    maxeig[1] = std::max(maxeig[1], fabs(v));
	    maxeig[4] = std::max(maxeig[4], fabs(v + a));
	    for (j = 0; j < 6; ++j)
                u[j] = u_old[j][i];
	    u2f(u,f[i]);
    	}
    	maxeig[2] = maxeig[1];
    	maxeig[3] = maxeig[1];

    	for(i = ghost_size; i < extend_size - ghost_size + 1; ++i)
    	{
	    /*** Get u_1/2 ***/

	    for(j = 0; j < 5; ++j)
	    {
	    	u_mid[j] = 0.5*(u_old[j][i-1] + u_old[j][i]);
	    }
	    u_mid[5] = u_mid[1]/u_mid[0];
	    u_mid[6] = u_mid[2]/u_mid[0];
	    u_mid[7] = u_mid[3]/u_mid[0];
	    u_mid[9] = sqr(u_mid[5]) + sqr(u_mid[6]) + sqr(u_mid[7]);
	    u_mid[8] = (u_mid[4] - 0.5*u_mid[0]*u_mid[9])*(gamma - 1.0);
	    u_mid[10] = sqrt(gamma*u_mid[8]/u_mid[0]);

	    /*** R(u_1/2) & R^-1(u_1/2) ***/

	    h = 0.5*u_mid[9] + gamma*u_mid[8]/(gamma - 1.0)/u_mid[0];

	    R[0][0] = 1.0;
	    R[0][1] = 1.0;
	    R[0][2] = 0.0;
	    R[0][3] = 0.0;
	    R[0][4] = 1.0;

	    R[1][0] = u_mid[5] - u_mid[10];
	    R[1][1] = u_mid[5];
	    R[1][2] = 0.0;
	    R[1][3] = 0.0;
	    R[1][4] = u_mid[5] + u_mid[10];

	    R[2][0] = u_mid[6];
	    R[2][1] = u_mid[6];
	    R[2][2] = 1.0;
	    R[2][3] = 0.0;
	    R[2][4] = u_mid[6];

	    R[3][0] = u_mid[7];
	    R[3][1] = u_mid[7];
	    R[3][2] = 0.0;
	    R[3][3] = 1.0;
	    R[3][4] = u_mid[7];

	    R[4][0] = h - u_mid[5]*u_mid[10];
	    R[4][1] = 0.5*u_mid[9];
	    R[4][2] = u_mid[6];
	    R[4][3] = u_mid[7];
	    R[4][4] = h + u_mid[5]*u_mid[10];;

	    L[0][0] = h + u_mid[10]*(u_mid[5] - u_mid[10]) / gm;
	    L[0][1] = -1.0*(u_mid[5] + u_mid[10] / gm);
	    L[0][2] = -1.0*u_mid[6];
	    L[0][3] = -1.0*u_mid[7];
	    L[0][4] = 1.0;

	    L[1][0] = -2.0*h + 4.0*sqr(u_mid[10])/gm;
	    L[1][1] = 2.0*u_mid[5];
	    L[1][2] = 2.0*u_mid[6];
	    L[1][3] = 2.0*u_mid[7];
	    L[1][4] = -2.0;

	    L[2][0] = -2.0*u_mid[6]*sqr(u_mid[10])/gm;
	    L[2][1] = 0;
	    L[2][2] = 2.0*sqr(u_mid[10])/gm;
	    L[2][3] = 0;
	    L[2][4] = 0;

	    L[3][0] = -2.0*u_mid[7]*sqr(u_mid[10])/gm;
	    L[3][1] = 0;
	    L[3][2] = 0;
	    L[3][3] = 2.0*sqr(u_mid[10])/gm;
	    L[3][4] = 0;

	    L[4][0] = h - u_mid[10]*(u_mid[5] + u_mid[10])/gm;
	    L[4][1] = -1.0*u_mid[5] + u_mid[10]/gm;
	    L[4][2] = -1.0*u_mid[6];
	    L[4][3] = -1.0*u_mid[7];
	    L[4][4] = 1.0;

	    for(j = 0; j < 5; ++j)
	    for(k = 0; k < 5; ++k)
		L[j][k] *= gm/(2.0*sqr(u_mid[10])); 

	    /*** Get R^-1 * u and R^-1 * F ***/	    

	    for(j = 0; j < 6; ++j)
	    {
		for (k = 0; k < 5; ++k)
                    u[k] = u_old[k][i - ghost_size + j];
	    	matmvec(sten_u[j],L,u);
	    	matmvec(sten_f[j],L,f[i - ghost_size + j]);
	    }

	    for(j = 0; j < 5; ++j)
	    for(k = 0; k < 5; ++k)
	    {
		gfluxp[j][k] = 0.5*(sten_f[k][j] + maxeig[j]*sten_u[k][j]);
		gfluxm[j][k] = 0.5*(sten_f[5 - k][j] - maxeig[j]*
				sten_u[5 - k][j]);
	    }

	    for(j = 0; j < 5; ++j)
	    {
		f_prevp[j] = f_nowp[j];
	    	f_nowp[j] = weno5_scal(gfluxp[j]);
		f_tmp[j] = f_nowp[j];
		/* artificial compression */
	    	for(k = 0; k < 5; ++k)
			gflux_tmp[k] = 0.5*(sten_f[5-k][j] + maxeig[j]*sten_u[5-k][j]);
		arti_compression(params,gfluxp[j],gflux_tmp,f_prevp[j],f_nowp[j],vecp[j],1,c);
		f_tmp[j] += c;
		/* end of artificial compression */
		f_prevm[j] = f_nowm[j];
	    	f_nowm[j] = weno5_scal(gfluxm[j]);
		f_tmp[j] += f_nowm[j];
		/* artificial compression */
	    	for(k = 0; k < 5; ++k)
			gflux_tmp[k] = 0.5*(sten_f[k][j] - maxeig[j]*sten_u[k][j]);
		arti_compression(params,gfluxm[j],gflux_tmp,f_prevm[j],f_nowm[j],vecm[j],-1,c);
		f_tmp[j] += c;
		/* end of artificial compression */
	    }

	    matmvec(ff, R, f_tmp);
	    for(j = 0; j < 5; ++j)
	    {
                flux[j][i] = ff[j];
		if (isnan(f[j][i]))
		{
		    (void) printf("In weno5_flux(): flux[%d][%d] = %f\n",j,i,
					f[j][i]);
		    for (k = 0; k < extend_size; ++k)
                        printf("u[%d] = %f %f %f %f %f\n",k,u_old[0][k],
				u_old[1][k],u_old[2][k],u_old[3][k],
				u_old[4][k]);
                    clean_up(ERROR);
		}
	    }
    	}
}

static double weno5_scal(double *f)
{
    	int i, j;
    	const double eps = 1.e-40;
    	const int p = 2;

    	double d[3] = {0.1, 0.6, 0.3}; //*** Optimal weights C_k 
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
	    alpha[i] = d[i]/pow(eps + is[i],p);
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

static void f2is(
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

static void matmvec(
	double *b, 
	double L[5][5], 
	double *x)
{
    	int i, j;

    	for(i = 0; i < 5; ++i)
    	{
	    b[i] = 0.0;
	    for(j = 0; j < 5; ++j)
	    {
	    	b[i] += L[i][j] * x[j]; 
	    }
    	}
}

static void u2f(
	double *u,
        double *f)
{
	double v = u[1]/u[0];

    	f[0] = u[1];
    	f[1] = v*u[1] + u[5];
    	f[2] = v*u[2];
    	f[3] = v*u[3];
    	f[4] = v*(u[4] + u[5]);
}

static void arti_compression(
	POINTER params,
	double *sten,
	double *sten_tmp,
	double f_prev,
	double f_now,
	double *vec,
	int sign,
	double &c)
{
	double tmp,f_ac_tmp;
	double ac_alpha,ac_alpha_nume,ac_alpha_deno;
	SCHEME_PARAMS *scheme_params = (SCHEME_PARAMS*)params; 

	c = 0.0;
	if (scheme_params->artificial_compression == NO) return;
	ac_alpha_nume = fabs(sten[3] - 2.0*sten[2] + sten[1]);
	ac_alpha_deno = fabs(sten[3] - sten[2]) + fabs(sten[2] - sten[1]);
	f_ac_tmp = weno5_scal(sten_tmp);
	if (sign == 1)
	{
	    vec[1] = vec[0];
	    vec[0] = f_ac_tmp - f_now;
	    vec[2] = sten[3] - f_now;
	    vec[3] = vec[1] + f_prev - sten[1];
	}
	else if (sign == -1)
	{
	    vec[0] = vec[1];
	    vec[1] = f_ac_tmp - f_now;
	    vec[2] = sten[1] - f_prev;
    	    vec[3] = vec[1] + f_now - sten[3];
	}
	if (ac_alpha_deno != 0.0)
	{
	    ac_alpha = 33.0 * pow(ac_alpha_nume/ac_alpha_deno,2);
	    if (ac_alpha > 1.0)
	    {
		if (vec[0] > 0.0 && vec[1] > 0.0)
		    tmp = (vec[0] < vec[1]) ? vec[0] : vec[1];
		else if (vec[0] < 0.0 && vec[1] < 0.0)
		    tmp = (vec[0] > vec[1]) ? vec[0] : vec[1];
		else
		    tmp = 0.0;
		tmp *= ac_alpha / 2.0;
		if (tmp > 0.0 && vec[2] > 0.0 && vec[3] > 0.0)
		{
		    c = (tmp < vec[2]) ? tmp : vec[2];
		    c = (c < vec[3]) ? c : vec[3];
		}
		else if (tmp < 0.0 && vec[2] < 0.0 && vec[3] < 0.0)
		{
		    c = (tmp > vec[2]) ? tmp : vec[2];
		    c = (c > vec[3]) ? c : vec[3];
		}
		else
		    c = 0.0;
	    }
	}
}

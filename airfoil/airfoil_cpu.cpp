#include "airfoil_cpu.h"

extern void call_cputest(double **x_old,double **x_mid,double **x_new,double **v_old,
		double **v_mid,double **v_new,double **f_mid,double dt,int size)
{
	for (int i = 0; i < size; ++i)
	for (int j = 0; j < 3; ++j)
	{
	    x_new[i][j] += dt*v_mid[i][j]/3.0;
	    v_new[i][j] += dt*f_mid[i][j]/3.0;
	    x_mid[i][j] = x_old[i][j] + v_mid[i][j]*dt;
	    v_mid[i][j] = v_old[i][j] + f_mid[i][j]*dt;
	}
}

extern void call_cpu_dummy()
{
	const long n = 1024*1024;
	const int m = 10000;
	float data[n];
	float *d_data;

	for(int i = 0; i < n; i++)
	{
	    data[i] = 1.0f*i/n;
	    for(int j = 0; j < m; j++)
	    {
		data[i] = data[i]*data[i] - 0.25f;
	    }
	}
}

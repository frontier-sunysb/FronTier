#ifndef _AIRFOIL_CPU_
#define _AIRFOIL_CPU_

void call_cputest(double **x_old, double **x_mid, double **x_new, double **v_old, double **v_mid, double **v_new, double **f_mid, double dt, int size);

void call_cpu_dummy();

#endif

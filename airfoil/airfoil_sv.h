#ifndef _AIRFOIL_SV_
#define _AIRFOIL_SV_

struct _SPRING_VERTEX {
        double *x;
        double *v;
	double *f;
        double *ext_impul;
	int ix;
        int num_nb;
        double m;
        double lambda;
        double **x_nb;
	int *ix_nb;
        double *k;
        double *len0;
        double ext_accel[3];
};
typedef struct _SPRING_VERTEX SPRING_VERTEX;

struct _GLOBAL_POINT {
        double x[3];
        double v[3];
	double f[3];
        double impuls[3];
        long gindex;
};
typedef struct _GLOBAL_POINT GLOBAL_POINT;

#endif

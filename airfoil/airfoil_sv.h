#ifndef _AIRFOIL_SV_
#define _AIRFOIL_SV_

struct _SPRING_VERTEX {
        double *x;
        double *v;
	double *f;
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

#endif

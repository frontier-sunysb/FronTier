/**********************************************************************
 * 		ifluid_basic.h
 **********************************************************************/

#ifndef _FT_IFLUID_BASIC_H_
#define _FT_IFLUID_BASIC_H_

enum _IF_PROB_TYPE {
        ERROR_TYPE = -1,
        TWO_FLUID_BUBBLE = 1,
        TWO_FLUID_RT,
        TWO_FLUID_KH,
        FLUID_SOLID_CIRCLE,
        BUBBLE_SURFACE,
        FLUID_RIGID_BODY,
        ROTOR_ONE_FLUID,
        ROTOR_TWO_FLUID,
	CHANNEL_FLOW,
        FLUID_CRYSTAL
};
typedef enum _IF_PROB_TYPE IF_PROB_TYPE;

struct _STATE {
	double dens;			/* Density */
        double pres;                    /* Pressure */
        double phi;                    	/* Potential */
        double vel[MAXD];               /* Velocities */
        double vort;                    /* Vorticity in 2D */
        double vort3d[MAXD];            /* Vorticity in 3D */
	double Impct[MAXD];             /* Accum impact from external force */
};
typedef struct _STATE STATE;

extern void restart_set_dirichlet_bdry_function(Front*);
extern void iF_flowThroughBoundaryState(double*,HYPER_SURF*,Front*,POINTER,
                        POINTER);
extern void iF_timeDependBoundaryState(double*,HYPER_SURF*,Front*,POINTER,
                        POINTER);
extern void ifluid_point_propagate(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
extern void ifluid_compute_force_and_torque(Front*,CURVE*,double,double*,
                        double*);
extern void setInitialIntfc(Front*,LEVEL_FUNC_PACK*,char*,IF_PROB_TYPE);
extern void init_fluid_state_func(Incompress_Solver_Smooth_Basis*,IF_PROB_TYPE);
extern void read_iFparams(char*,IF_PARAMS*);
extern void read_iF_prob_type(char*,IF_PROB_TYPE*);
extern void recordBdryEnergyFlux(Front*,char*);
#endif

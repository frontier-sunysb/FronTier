
#include <FronTier.h>

enum _PROB_TYPE {
        ERROR_TYPE = -1,
	FLUID_CRYSTAL
};
typedef enum _PROB_TYPE PROB_TYPE;

enum {
	CRYSTAL_BOUNDARY	= FIRST_PHYSICS_WAVE_TYPE
};

struct _STATE {
	double dens;			/* Density */
        double pres;                    /* Pressure */
	double phi;			/* Potential */
        double vel[MAXD];               /* Velocities */
        double vort;                    /* Vorticity */
        double vort3d[MAXD];            /* Vorticity in 3D */
	double Impct[MAXD];             /* Accum impact from external force */
        double solute;                  /* Solute concentration */
};
typedef struct _STATE STATE;

extern void read_dirichlet_bdry_data(char*,Front*);
extern void ifluid_point_propagate(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
extern void crystal_point_propagate(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
extern void read_ss_dirichlet_bdry_data(char*,Front*,F_BASIC_DATA);
extern void init_fluid_state_func(Front*,Incompress_Solver_Smooth_Basis*);
extern void ss_flowThroughBoundaryState(double*,HYPER_SURF*,Front*,POINTER,
                        POINTER);
extern void read_fluid_params(char*,IF_PARAMS*);
extern void initFrontStates(Front*);

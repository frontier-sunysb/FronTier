#include <FronTier.h>

enum _RG_PROB_TYPE {
        ERROR_TYPE = -1,
        FLUID_SOLID_CIRCLE,
        BUBBLE_SURFACE,
        FLUID_RIGID_BODY,
        ROTOR_ONE_FLUID,
        ROTOR_TWO_FLUID,
	WINDMILL_2D,
	WINDMILL_3D,
	BEE_3D,
	HELICOPTER_3D
};
typedef enum _RG_PROB_TYPE RG_PROB_TYPE;

/*	rgbody.c functions */

struct _STATE {
	double dens;			/* Density */
	double phi;                     /* Potential */
        double pres;                    /* Pressure */
        double vel[MAXD];               /* Velocities */
        double vort;                    /* Vorticity in 2D */
        double vort3d[MAXD];            /* Vorticity in 3D */
	double Impct[MAXD];             /* Accum impact from external force */
};
typedef struct _STATE STATE;

struct _RG_PARAMS {
	int dim;
	double  total_mass;             /* Total mass */
        double  moment_of_inertial;     /* Moment of inertial about the axis */
        double  center_of_mass[MAXD];   /* Center of mass */
        double  rotation_dir[MAXD];     /* Direction of rotation */
        double  rotation_cen[MAXD];     /* Center of rotation */
        double  cen_of_mass_velo[MAXD]; /* Center of mass velocity */
        double  angular_velo;           /* Angular velocity of rotation */
	MOTION_TYPE motion_type;
};
typedef struct _RG_PARAMS RG_PARAMS;

void init_moving_bodies(Front*,LEVEL_FUNC_PACK*,char*,RG_PROB_TYPE);
void ifluid_compute_force_and_torque(Front*,HYPER_SURF*,double,double*,
			double*);
void record_moving_body_data(char*,Front*);
void read_iFparams(char*,IF_PARAMS*);
void read_rg_prob_type(char*,RG_PROB_TYPE*);
void read_movie_options(char*,IF_PARAMS*);
void read_dirichlet_bdry_data(char*,Front*,F_BASIC_DATA);
void restart_set_dirichlet_bdry_function(Front*);
void ifluid_point_propagate(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
extern void init_fluid_state_func(Incompress_Solver_Smooth_Basis*,RG_PROB_TYPE);

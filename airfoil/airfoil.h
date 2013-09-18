#ifndef _AIRFOIL_
#define _AIRFOIL_

#include <FronTier.h>
#include "airfoil_sv.h"

typedef struct {
        double N[MAXD];         /* normal of the plane */
        double P[MAXD];         /* a point on the plane */
        double cen[MAXD];       /* center of the vent */
        double radius;          /* radius of the vent */
} CONSTR_PARAMS;

enum _PERTURBATION_TYPE {
	NO_PERT		=	1,
	PARALLEL_RAND_PERT,
	ORTHOGONAL_RAND_PERT,
	LINEAR_PERT,
	RADIAL_PERT,
	SINE_PERT
};
typedef enum _PERTURBATION_TYPE PERTURBATION_TYPE;

enum _STRING_NODE_TYPE {
	FIXED_END		=	1,
	FREE_END,
	LOADED_END
};
typedef enum _STRING_NODE_TYPE STRING_NODE_TYPE;

enum _AF_NODE_TYPE {
	UNKNOWN_AF_NODE = -1,
	LOAD_NODE	= 1,
	GORE_NODE,
	STRING_NODE,
	PRESET_NODE
};
typedef enum _AF_NODE_TYPE AF_NODE_TYPE;

enum _SPRING_MODEL {
	UNKNOWN_MODEL 	= -1,
	MODEL1	= 1,
	MODEL2,
	MODEL3
};
typedef enum _SPRING_MODEL SPRING_MODEL;

struct _PERT_PARAMS {
	PERTURBATION_TYPE pert_type;
	int dir;
	double x0,xl,xu;
	double pert_amp;
	double cen[MAXD];
	double pert_radius;
};
typedef struct _PERT_PARAMS PERT_PARAMS;

typedef struct {
        int dim;
        POINTER level_func_params;
        NS_SCHEME num_scheme;
        SPRING_MODEL spring_model;
	boolean no_fluid;
	boolean is_parachute_system;
	boolean attach_gores;
	boolean attach_fixer;
	boolean cut_vent;
        boolean use_total_mass;
	boolean use_gpu;
	PERT_PARAMS pert_params;
	STRING_NODE_TYPE start_type;
	STRING_NODE_TYPE end_type;
	double gore_len_fac;
        double rho1;
        double rho2;
        double mu1;
        double mu2;
        double U1[MAXD];
        double U2[MAXD];
        double gravity[MAXD];		/* gravitational force */
	double payload;
        double surf_tension;
        double smoothing_radius;
	double ks;			/* spring constant of surface */
	double kl;			/* spring constant of string curves */
	double kg;                      /* spring constant of gore curves */
	double lambda_s;		/* damping factor of surface */
	double lambda_l;		/* damping factor of string curves */
	double lambda_g;                /* damping factor of gore curves */
	double m_s;			/* point mass of surface */
	double m_l;			/* point mass of string curves */
	double m_g;                     /* point mass of gore curves */
	double total_string_mass;	/* Total mass of string chord */
	double total_canopy_mass;	/* Total mass of string chord */
        double total_gore_mass;         /* Total mass of gore */
	double gamma;			/* canopy porosity */
	double area_dens;		/* canopy area density */
	double min_len;
	int    n_tan;			/* number of sub-steps for tan prop */
	int    num_opt_round;		/* number of canopy optimizations 
							rounds*/ 
	int    num_smooth_layers;	/* number of layer to smooth high
					   frequency velocity */
        IF_FIELD *field;
} AF_PARAMS;

/*	airfoil.cpp functions */

typedef struct {
        double i1,i2;
        double cen1[2],cen2[2];
} DOUBLE_VORTEX_PARAMS;

typedef struct {
        double tcen[MAXD]; /* toroidal center */
	double R0;	/* distance between poloidal and toroidal centers */
	double v0;	/* amplitude of velocity */
	double stop_time;
} TOROIDAL_PARAMS;

typedef struct {
        double cen[MAXD]; /* parabolic center */
	double v0;	  /* velocity at center */
	double a;	  /* concavity (downward) */
	double stop_time;
} PARABOLIC_PARAMS;

typedef struct {
        double cen[MAXD]; /* parabolic center */
	double R;	  /* radius of for v0 */
	double v0;	  /* velocity at center */
	double stop_time;
} SINGULAR_PARAMS;

struct _STRING_PARAMS {
	int num_strings;
	double start_angle;
	double coords_load[MAXD];
	double cen[MAXD];
	double shift[MAXD];
	double theta;
	double phi;
	double L[MAXD],U[MAXD];
	double P[MAXD];
};
typedef struct _STRING_PARAMS STRING_PARAMS;

struct _PARALLEL_GORE_PARAMS {
        int gores_n;
        double gores_start_x;
        double gores_dis;

        double coords_load[MAXD];
};

typedef struct _PARALLEL_GORE_PARAMS PARALLEL_GORE_PARAMS;

enum _LOAD_TYPE {
	NO_LOAD 	= 0,
	FREE_LOAD,
	RIGID_LOAD
};
typedef enum _LOAD_TYPE LOAD_TYPE;

typedef struct {
        boolean lower_bdry[MAXD];
        boolean upper_bdry[MAXD];
        double L[MAXD];         /* Lower bounds of box */
        double U[MAXD];         /* Upper bounds of box */
	LOAD_TYPE lower_side[MAXD];
	LOAD_TYPE upper_side[MAXD];
	double lower_mass[MAXD];
	double upper_mass[MAXD];
	double lower_force[MAXD][MAXD];
	double upper_force[MAXD][MAXD];
} BDRY_PARAMS;

typedef struct {
	LOAD_TYPE load_type;
	int dir;
	double load_mass;
	double point_mass;
	double force[MAXD];
	double ave_accel;
} C_PARAMS;

struct _AF_NODE_EXTRA {
	AF_NODE_TYPE af_node_type;
};
typedef struct _AF_NODE_EXTRA AF_NODE_EXTRA;

struct _REGISTERED_PTS {
	int num_pts;
	int *global_ids;
};
typedef struct _REGISTERED_PTS REGISTERED_PTS;

struct _PARACHUTE_SET{
	Front *front;
        SURFACE *canopy;
        CURVE **mono_hsbdry;
        CURVE **gore_hsbdry;
        NODE **string_node;
	NODE **gore_nodes;
        NODE *load_node;
        CURVE **string_curves;
        int num_mono_hsbdry;
        int num_gore_hsbdry;
        int num_strings;
	int num_gore_nodes;
	double ks;
	double kl;
	double kg;
	double lambda_s;
	double lambda_l;
	double lambda_g;
	double m_s;
	double m_l;
	double m_g;
	double V_surf[MAXD];
	double V_load[MAXD];
	int n_cps;		/* Number of points on canopy */
	int n_sps;		/* Number of points on string */
	double dt;
};

typedef struct _PARACHUTE_SET PARACHUTE_SET;

void read_iFparams(char*,IF_PARAMS*);
void read_movie_options(char*,IF_PARAMS*);
void read_dirichlet_bdry_data(char*,Front*,F_BASIC_DATA);
void restart_set_dirichlet_bdry_function(Front*);
void liquid_point_propagate(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);

// afinit.cpp
extern void setInitialIntfcAF(Front*,LEVEL_FUNC_PACK*,char*);

/* afinit3d.cpp */
extern void initEllipticSurf(FILE*,Front*,LEVEL_FUNC_PACK*);
extern void initParabolicSurf(FILE*,Front*,LEVEL_FUNC_PACK*);
extern void initPlaneSurf(FILE*,Front*,LEVEL_FUNC_PACK*);

// afprop.cpp
extern void airfoil_point_propagate(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
extern void elastic_point_propagate(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
extern void ifluid_point_propagate(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
extern void fourth_order_elastic_set_propagate(Front*,double);
extern void airfoil_curve_propagate(Front*,POINTER,CURVE*,CURVE*,double);
extern  int numOfGoreHsbdry(INTERFACE*);
extern  int numOfMonoHsbdry(INTERFACE*);
extern  int numOfGoreNodes(INTERFACE*);
extern boolean is_load_node(NODE*);
extern boolean is_gore_node(NODE*);
extern boolean is_bdry_node(NODE*);
extern boolean is_string_node(NODE*);
extern double springCharTimeStep(Front*);	// spring characteristic time

// aftest.cpp
extern void second_order_elastic_curve_propagate(Front*,Front*,INTERFACE*,
                                CURVE*,CURVE*,double);
extern void second_order_elastic_surf_propagate(Front*,double);
extern void set_equilibrium_mesh(Front*);
extern void unsort_surf_point(SURFACE*);
extern void print_airfoil_stat(Front*,char*);
extern void fixed_length_tan_curve_propagate(Front*,Front*,INTERFACE*,
                                CURVE*,CURVE*,double);
extern void fourth_order_elastic_curve_propagate(Front*,Front*,INTERFACE*,
                                CURVE*,CURVE*,double);
extern void fourth_order_elastic_surf_propagate(Front*,double);
extern void legacy_fourth_order_elastic_curve_propagate(Front*,Front*,
				INTERFACE*,CURVE*,CURVE*,double);
extern void legacy_fourth_order_elastic_surf_propagate(Front*,double);

// afcnpy.cpp
extern void coating_mono_hyper_surf(Front*);
extern void compute_total_canopy_force(Front*,double*,double*);
extern int airfoil_velo(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*,
                                double*);
extern int af_find_state_at_crossing(Front*,int*,GRID_DIRECTION,
                        int,POINTER*,HYPER_SURF**,double*);
extern void fourth_order_parachute_propagate(Front*,PARACHUTE_SET*);
extern void assign_node_field(NODE*,double**,double**,int*);
extern void assign_curve_field(CURVE*,double**,double**,int*);
extern void assign_surf_field(SURFACE*,double**,double**,int*);
extern void compute_surf_accel1(PARACHUTE_SET*,SURFACE*,double**,double**,
				double**,int*);
extern void compute_surf_accel2(PARACHUTE_SET*,SURFACE*,double**,double**,
				double**,int*);
extern void compute_curve_accel1(PARACHUTE_SET*,CURVE*,double**,double**,
				double**,int*);
extern void compute_node_accel1(PARACHUTE_SET*,NODE*,double**,double**,double**,
				int*);
extern void compute_curve_accel2(PARACHUTE_SET*,CURVE*,double**,double**,
				double**,int*);
extern void compute_node_accel2(PARACHUTE_SET*,NODE*,double**,double**,double**,
				int*);
extern void compute_curve_accel3(PARACHUTE_SET*,CURVE*,double**,double**,
				double**,int*);
extern void compute_node_accel3(PARACHUTE_SET*,NODE*,double**,double**,double**,
				int*);
extern void propagate_surface(PARACHUTE_SET*,SURFACE*,double**,int*);
extern void propagate_curve(PARACHUTE_SET*,CURVE*,double**,int*);
extern void propagate_node(PARACHUTE_SET*,NODE*,double**,int*);
extern boolean is_registered_point(SURFACE*,POINT*);

// afsetd.cpp
extern void count_node_neighbors(NODE*,SPRING_VERTEX*,int*);
extern void count_curve_neighbors(CURVE*,SPRING_VERTEX*,int*);
extern void count_surf_neighbors(SURFACE*,SPRING_VERTEX*,int*);
extern void count_canopy_spring_neighbors(PARACHUTE_SET*,SPRING_VERTEX*);
extern void count_string_spring_neighbors(PARACHUTE_SET*,SPRING_VERTEX*);
extern void set_node_spring_vertex(PARACHUTE_SET*,NODE*,double**,double**,
				SPRING_VERTEX*,int*);
extern void set_curve_spring_vertex(PARACHUTE_SET*,CURVE*,double**,double**,
				SPRING_VERTEX*,int*);
extern void set_surf_spring_vertex(PARACHUTE_SET*,SURFACE*,double**,double**,
				SPRING_VERTEX*,int*);
extern void set_canopy_spring_vertex(PARACHUTE_SET*,double**,double**,
				SPRING_VERTEX*);
extern void set_string_spring_vertex(PARACHUTE_SET*,double**,double**,
				SPRING_VERTEX*);
extern void set_spring_vertex_memory(SPRING_VERTEX*,int);
extern void compute_spring_accel1(SPRING_VERTEX,double*,int);
extern void generic_spring_solver(SPRING_VERTEX*,double**,double**,int,int,int,
				double);

// afvelo.cpp
extern void setMotionParams(Front*);
extern void resetFrontVelocity(Front*);

// afmodule.cpp
extern void initParachuteModules(Front*);

// afdata.cpp
extern void printAfExtraDada(Front*,char*);
extern void readAfExtraDada(Front*,char*);
extern void printHyperSurfQuality(Front*);
extern void optimizeElasticMesh(Front*);
extern void modifyInitialization(Front*);

// afprops.cpp

// sprModel/modules.cpp
extern void initSpringModel(Front*);
#endif


#include <FronTier.h>

enum {
	CRYSTAL_BOUNDARY	= FIRST_PHYSICS_WAVE_TYPE
};

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

#include "FronTier.h"

#define CIM_TOL 1e-8

struct _CIM_PARAMS {
	double h[MAXD];
	int Run_case;
	int intfc_num;
	int dim;
};
typedef struct _CIM_PARAMS CIM_PARAMS;

struct _STATE {
	double u;
	double uex;
};
typedef struct _STATE STATE;

// Old interface.h

extern double exact_eps(int D);
extern double exact_solution(POINTER jparams,int D, double *P);
extern double exact_source(POINTER params, int D, double *P);
extern double exact_jump_u(POINTER jparams,int D, double *P);
extern double exact_jump_eps_gradu_dot_n(POINTER jparams,int D,
		double *N,double *P);
extern double exact_jump_gradu_dot_t(POINTER params,int D,int i,
		double *N,double *P);
extern int cim_find_state_at_crossing(Front*,int*,GRID_DIRECTION,int,
        			POINTER*,HYPER_SURF**,double*);
extern double getStateU(POINTER);
extern void initCimIntfcParams(char*,Front*,LEVEL_FUNC_PACK*);
extern double intfc_func_case2(POINTER,double*);
extern double intfc_func_case3(POINTER,double*);
extern double intfc_func_case4(POINTER,double*);
extern double intfc_func_case5(POINTER,double*);
extern double intfc_func_case6(POINTER,double*);
extern double intfc_func_case7(POINTER,double*);
extern double intfc_func_case8(POINTER,double*);
extern double intfc_func_case9(POINTER,double*);
extern double intfc_func_case10(POINTER,double*);

class C_CARTESIAN{
        Front *front;
public:
        C_CARTESIAN(Front &front);
        int dim;

        // On topological grid
        RECT_GRID *top_grid;
        double *top_L,*top_U,*top_h,hmin;
        int *top_gmax;
        COMPONENT *top_comp;

        double *array;          // for scatter states;
        double *source;         // for scatter states;
        double *diff_coeff;     // for scatter states;
        double *numer_soln;     // numerical solution;
        double *exact_soln;     // exact solution;

        int *lbuf,*ubuf;
        int *i_to_I,*I_to_i;            // Index mapping for 1D
        int **ij_to_I,**I_to_ij;        // Index mapping for 2D
        int ***ijk_to_I,**I_to_ijk;     // Index mapping for 3D
	int w_type;
	COMPONENT pos_comp,neg_comp;

        // Sweeping limites
        int imin,jmin,kmin;
        int imax,jmax,kmax;

	// constructor
        ~C_CARTESIAN();

	int size;

	// for parallel partition
        int             NLblocks,ilower,iupper;
        int             *n_dist;
	// mesh: full cells mesh
	void solve();
	void setDomain();
	void setIndexMap(void);
	void scatMeshArray();
        void setGlobalIndex();
	void initMovieVariables(void);

	void cim_solve(void);
	void cim_solve2d(void);
	void cim_solve3d(void);
	void old_solve(void);
};


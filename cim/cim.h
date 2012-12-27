#include "FronTier.h"
#include "solver.h"

#define CIM_TOL 1e-8

enum _COMPARE_METHOD {
        NO_COMPARE = 1,
        CAUCHY_COMPARE,
        EXACT_COMPARE
};
typedef enum _COMPARE_METHOD COMPARE_METHOD;

enum _JUMP_TYPE {
        EXACT_JUMP = 1,
        CONSTANT_JUMP,
	FUNCTIONAL_JUMP
};
typedef enum _JUMP_TYPE JUMP_TYPE;

struct _FIELD {
        double *u;
        double *uex;
};
typedef struct _FIELD FIELD;

struct _CIM_PARAMS {
	FIELD *field;
	int dim;
	double h[MAXD];

	COMPARE_METHOD compare_method;
	JUMP_TYPE jump_type;

	int intfc_num;

	// For constant jump condition
	double jump_u;
	double jump_eps_grad_u_dot_n;
	double jump_grad_u_dot_t;

	// For exact jump condition
	int Run_case;

	char base_dir_name[200];
	int num_step;
	int *steps;
	F_BASIC_DATA *f_basic;
};
typedef struct _CIM_PARAMS CIM_PARAMS;

struct _STATE {
	double u;
	double uex;
};
typedef struct _STATE STATE;

// Old interface.h

extern double sourceFunc(POINTER params, int D, double *P);
extern double jumpU(POINTER jparams,int D, double *P);
extern double jumpEpsGradDotNorm(POINTER jparams,int D,
		double *N,double *P);
extern double jumpGradDotTan(POINTER params,int D,int i,
		double *N,double *P);

extern double getStateU(POINTER);
extern void assignStateU(double,POINTER);
extern void initCimIntfcParams(char*,Front*,LEVEL_FUNC_PACK*);
extern void print_front_states(FILE *outfile,Front *front);
extern void read_front_states(FILE *infile,Front *front);
extern int cim_find_state_at_crossing(Front*,int*,GRID_DIRECTION,int,
        			POINTER*,HYPER_SURF**,double*);

// Exact cases
extern double intfc_func_case2(POINTER,double*);
extern double intfc_func_case3(POINTER,double*);
extern double intfc_func_case4(POINTER,double*);
extern double intfc_func_case5(POINTER,double*);
extern double intfc_func_case6(POINTER,double*);
extern double intfc_func_case7(POINTER,double*);
extern double intfc_func_case8(POINTER,double*);
extern double intfc_func_case9(POINTER,double*);
extern double intfc_func_case10(POINTER,double*);
extern double exact_eps(int D);
extern double exact_solution(POINTER jparams,int D, double *P);
extern double exact_jump_u(POINTER jparams,int D, double *P);
extern double exact_jump_gradu_dot_t(POINTER params,int D,int i,
		double *N,double *P);
extern double exact_jump_eps_gradu_dot_n(POINTER params,int D,
		double *N,double *P);
extern double exact_source(POINTER params, int D, double *P);


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
	CIM_PARAMS *params;

        double *array;          // for scatter states;
        double *source;         // for scatter states;
        double *diff_coeff;     // for scatter states;
        double *exact_soln;     // for scatter states;
	FIELD *field;

        int *lbuf,*ubuf;
        int *i_to_I,*I_to_i;            // Index mapping for 1D
        int **ij_to_I,**I_to_ij;        // Index mapping for 2D
        int ***ijk_to_I,**I_to_ijk;     // Index mapping for 3D
	int w_type;
	COMPONENT pos_comp,neg_comp;
	Front *base_front;

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
	void printFrontInteriorStates(char *out_name);
	void readFrontInteriorStates(char *restart_state_name);
	void readBaseFront(CIM_PARAMS *cim_params,int i);
	void readBaseStates(char *restart_name);

	void cim_solve(void);
	void old_solve(void);
	void compareWithExacySoln();
	void compareWithBaseSoln();
};


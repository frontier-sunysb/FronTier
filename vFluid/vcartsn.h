#include <vector>
#include "FronTier.h"

struct _STATE {
	double vel[MAXD];
	double eta;
};
typedef struct _STATE STATE;

#define         SOLID_COMP              0
#define         LIQUID_COMP             3

enum _INIT_STATE{
	ZERO_STATE = 1,
	CONST_STATE,
	RAND_STATE,
	STEP_STATE,
	SHOCK_STATE,
	DISK_STATE
};
typedef enum _INIT_STATE INIT_STATE;

struct _VF_FIELD {
	double **vel;
	double *eta;
};
typedef struct _VF_FIELD VF_FIELD;

struct _PARAMS {
        int dim;
	VF_FIELD *field;
};
typedef struct _PARAMS PARAMS;

class RECTANGLE {
public:
        int index;                      // rectangle index
        int comp;
        double area;
        double coords[MAXD];
        int icoords[MAXD];

        RECTANGLE();

        void setCoords(double*,int);
};

class CARTESIAN{
	Front *front;
public:
	CARTESIAN(Front &front);
	
	// member data: RECT_GRID
	int dim;

	// On topological grid
	RECT_GRID *top_grid;
	double *array;		// for scatter states;
	double *source;		// for source of parabolic solver;
	double *top_L,*top_U,*top_h,hmin;
	int *top_gmax;
	COMPONENT *top_comp;
	PARAMS *eqn_params;
	VF_FIELD *field;
	int comp_size;

	int *lbuf,*ubuf,*gmax;
	int *i_to_I,*I_to_i;		// Index mapping for 1D
	int **ij_to_I,**I_to_ij;	// Index mapping for 2D
	int ***ijk_to_I,**I_to_ijk;	// Index mapping for 3D

	// Sweeping limites
	int imin,jmin,kmin;
	int imax,jmax,kmax;

	// member data: mesh storage
	std::vector<RECTANGLE> cell_center;

	double m_t;                     // time
	double m_dt;			// time increment
	double min_dt;			// Minimum dt to use non-explicit
	double max_dt;			// Maximum dt to use non-explicit

	// constructor
	~CARTESIAN();

	// for parallel partition
	int             NLblocks,ilower,iupper;
        int             *n_dist;

	// mesh: full cells mesh
	void initMesh(void);		// setup the cartesian grid
	void setDomain(void);
	void setComponent(void);	// init components	
	void setAdvectionDt(void);	// compute time step 	
	void readFrontInteriorState(char*);
	void printFrontInteriorState(char*);

	void computeAdvection();

	// interface functions
	void makeGridIntfc();
	void deleteGridIntfc();

	// Extra plot functions
	void initDrawVariables();

	// parallelization related functions
	void scatMeshArray();
	void setGlobalIndex(COMPONENT);

	// physics calculation
	void initFlowState(void);

	void setIndexMap(COMPONENT);
		// for compProjWithSmoothProperty(), 
		// should be changed to use setIndexMap() only.

	// -------------------------------------------------------
	// 		incompressible solver functions

	// main step function
	void solve(double dt);		

	void getVelocity(double *p, double *U);

	// ----------------------------------------------------------
	// 		utility functions
	// ----------------------------------------------------------

	void getRectangleIndex(int indexRectangle, int &i, int &j);
	void getRectangleIndex(int indexRectangle, int &i, int &j, int &k);
	void getRectangleCenter(int index, double *coords);
	void getRectangleCenter(int index0, int index1, double *coords);
	int  getRectangleComponent(int index);	// the center component
	int  getComponent(int *icoords);	
	int  getComponent(double *coords);	
};

extern double getStateXvel(POINTER);
extern double getStateYvel(POINTER);
extern double getStateEta(POINTER);
extern void read_vF_dirichlet_bdry_data(Front*);

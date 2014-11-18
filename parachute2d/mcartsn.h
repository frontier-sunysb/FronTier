#include "iFluid.h"
/*mcartesian.h*/
/*Following is for scalar advection: temperature*/
enum _INIT_STATE{
	ZERO_STATE = 1,
	CONST_STATE,
	RAND_STATE,
	STEP_STATE,
	SHOCK_STATE,
	DISK_STATE
};
typedef enum _INIT_STATE INIT_STATE;

struct _IF_FIELD;
typedef struct _IF_FIELD IF_FIELD;

struct _PARAMS {
        int dim;
	INIT_STATE init_state;
	IF_FIELD *field;
	double x0; /*place of discontinuity*/
	double T0;/*initial temperature */
	double Tl;/*initial left temperature */
	double Tr;/*initial right temperature */
	double D;    /*molecular diffusivities of the temperature*/
	int dir; /*direction of step*/
	double center[MAXD]; /* center of the disk state case */
        double radius; /* raduis of the disk state case */
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
	IF_FIELD *field;
	int comp_size;

	int *lbuf,*ubuf,*gmax;
	int *i_to_I,*I_to_i;		// Index mapping for 1D
	int **ij_to_I,**I_to_ij;	// Index mapping for 2D
	int ***ijk_to_I,**I_to_ijk;	// Index mapping for 3D

	// Sweeping limites
	int imin,jmin,kmin;
	int imax,jmax,kmax;

	enum BC_TYPE { 									// used by ADVECTION
		BC_PERIODIC = 1,
		BC_Extrapolation0 = 2,
		BC_Extrapolation1 = 3,
		BC_InflowOutflow = 4};	
	BC_TYPE m_bc[4];								// down, right, up, left 		

	// member data: mesh storage
	std::vector<RECTANGLE> cell_center;

	double m_t;                     // time
	double m_dt;			// time increment
	double min_dt;			// Minimum dt to use non-explicit

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
	void setBoundary(void);		// set up boundary conditions 	
	void readFrontInteriorState(char*);
	void printFrontInteriorState(char*);

	void computeAdvection();

	void computeAdvectionExplicitCim(COMPONENT);
	void computeAdvectionExplicit(COMPONENT);
	void computeAdvectionImplicit(COMPONENT);
	void computeAdvectionCN(COMPONENT);
	void computeAdvectionCim();

	void pointExplicitCimSolver(int*,COMPONENT);
	void computeSource();

	// interface functions
	void makeGridIntfc();
	void deleteGridIntfc();

	// Extra plot functions
	void oneDimPlot(char*);
	void xgraphOneDimPlot(char*);
	void initMovieVariables();
	void augmentMovieVariables(const char*);
	void vtk_plot_temperature2d(char*);
        void vtk_plot3d(const char*);

	// Extra movie functions
	void temperatureMovie(char*);

	void checkStates();

	// parallelization related functions
	//
	void scatMeshArray();
	void setGlobalIndex(COMPONENT);

	// physics calculation
	void setInitialCondition(void);

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
	int getRectangleComponent(int index);	// the center component
	
	double getDistance(double *coords0, double *coords1);
	
			// incompletely implemented
	void getNearestInterfacePoint(double *q,double *p); 
		
	int  getComponent(int *icoords);	
	int  getComponent(double *coords);	
	void save(char *filename);
};


extern double getStateTemp(POINTER);
extern void read_params(char*,PARAMS*);

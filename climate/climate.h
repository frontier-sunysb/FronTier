/**********************************************************************
 * 		liquid.h
 * the code is a direct modification of the code by Leveque[1].
 *
 * References:
 * [1] R.J. Leveque, High-resolution conservative algorithms for
 *     advection in incompressible flow.
 *
 **********************************************************************/

#include <FronTier.h>
#include <vector>
#include <petscksp.h>
#include <assert.h>

#define         LIQUID_COMP2           3
#define         SOLID_COMP             1
#define		alternate_comp(comp) 					\
		(comp) == LIQUID_COMP2 ? SOLID_COMP : LIQUID_COMP2

enum _CL_PROB_TYPE {
	BUOYANCY_TEST,
	CLIMATE,
	CHANNEL_TEST,
	ENTRAINMENT,
	RANDOM_FIELD,
	PARTICLE_TRACKING
};
typedef enum _CL_PROB_TYPE CL_PROB_TYPE;

enum _NUM_SCHEME {
	UNSPLIT_EXPLICIT = 1,
	UNSPLIT_EXPLICIT_CIM,
	UNSPLIT_IMPLICIT,
	UNSPLIT_IMPLICIT_CIM,
	CRANK_NICOLSON
};
typedef enum _NUM_SCHEME NUM_SCHEME;

enum _INIT_STATE{
	ZERO_STATE = 1,
	RAND_STATE,
	TAYLOR_STATE,
	PRESET_STATE
};
typedef enum _INIT_STATE INIT_STATE;

struct _PHASE_FIELD {
        double *temperature;
	double *vapor;
	double *supersat;
	double *pres;
	double *drops;
        double **vel;
};
typedef struct _PHASE_FIELD PHASE_FIELD;

struct _MOVIE_OPTION {
        boolean plot_pres;
        boolean plot_vort;
        boolean plot_velo;
        boolean plot_temperature;
	boolean plot_vapor;
        boolean plot_cross_section[MAXD];  /* 3D 0: yz; 1: zx; 2: xy */
};
typedef struct _MOVIE_OPTION MOVIE_OPTION;

struct _PARTICLE {
	double radius;
	double center[MAXD];
	double vel[MAXD];
	double R0;
	double rho;
	int    Gindex;
	boolean flag;
};
typedef struct _PARTICLE PARTICLE;

struct _PARAMS {
        int dim;
	NUM_SCHEME num_scheme;
	INIT_STATE init_state;
	PHASE_FIELD *field;
	MOVIE_OPTION *movie_option;
	int pde_order;
	int num_phases;
	int num_drops;
        double *Ti;    /* melting temperature at the interface */
	double *T0;	/* Ambient temperature of the phase */
	double D;    /*molecular diffusivities of the temperature*/
	double min_temperature;
	double max_temperature;
	double rho_l; /*density of water droplet*/
	double K; /*coefficient for condensation*/
	double L[MAXD];
	double U[MAXD]; /*subdomain for particles*/
	boolean no_droplets;
	boolean droplets_fixed;
	CL_PROB_TYPE prob_type;
	PARTICLE *particle_array;
};
typedef struct _PARAMS PARAMS;

typedef class CARTESIAN CARTESIAN_EB;

/**************************************************************************
 *		vector/matrix operation functions 
 **************************************************************************/

void VectorIntPrint(int my_rank, char *name, int n, int *vector);
void VectorPrint(int my_rank, char *name, int n, double *vector);
void VectorZero(int n, double *vector);
void VectorZero(int n, int *vector);
void VectorCopy(int n, double *vector1, double *vector2);
void MatrixPrint(int my_rank, char *name, int m, int n, double **matrix);
void MatrixZero(int m, int n, double **matrix);
void MatrixIdentity(int n, double **matrix);
void MatrixCopy(int m, int n, double **U2, double **U1);
void ArrayZero(int m, int n, int l, double ***matrix);
void Vector2Matrix(int m, int n, double **matrix, double *vector);
void ArrayPrint(int my_rank, char*name, int m, int n, int l, double ***array);
void MatrixMultiply(int n, double **C, double **A, double **B);	// C = AB
void MatrixMultiply(int n, double **C, double **A);		// C = A'A 
void MatrixMultiply(double C[3][3], double A[3][3], double B[3][3]);  // C = AB
void MatrixVectorMultiply(int n, double *C, double **A, double *b);   // c = Ab
void SymmetricMatrixInverse(double B[3][3], double A[3][3]); // B = inverse(A); 
double MatrixMax(int m, int n, double **matrix);
double MatrixMin(int m, int n, double **matrix);

/******************************************************************************
 * 		lcartsn.h
 * Class for solving advection diffusion equation
 * dtT+u*divT = a1*laplace(T) + f
 *
 * the main function is 
 * 	CARTESIAN::solve().
 *
 * References:
 ******************************************************************************/

class SOLVER;
class CARTESIAN;

class RECTANGLE {
public:
	int index;			// rectangle index
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
	PHASE_FIELD *field;
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
	std::vector<RECTANGLE> 	cell_center;

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

#define MAX_STEP 10000

class VCARTESIAN{
	Front *front;
public:
	VCARTESIAN(Front &front);

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
	PHASE_FIELD *field;
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
	std::vector<RECTANGLE> 	cell_center;

	double m_t;                     // time
	double m_dt;			// time increment
	double min_dt;			// Minimum dt to use non-explicit

	// constructor
	~VCARTESIAN();

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
	void computeSupersat();
	void computeCondensation();
	double computeReactTimeScale(double,double,int);

	void pointExplicitCimSolver(int*,COMPONENT);
	void computeSource();

	// interface functions
	void makeGridIntfc();
	void deleteGridIntfc();

	// Extra plot functions
	void oneDimPlot(char*);
	void xgraphOneDimPlot(char*);
	void initMovieVariables();
	void checkOutput();
	void checkField();
	void recordField(char *, const char *);
	void recordRadius(char *);
	void recordMixingLine();
        void vtk_plot3d(const char*,double*);

	// Extra movie functions
	void temperatureMovie(char*);

	void checkStates();

	// parallelization related functions
	//
	void scatMeshArray();
	void setGlobalIndex(COMPONENT);

	// physics calculation
	void setInitialCondition(void);
	void setInitialVapor(void);
	void setParallelVapor(void);

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
	void recordWaterBalance();
	void recordSampleRadius();
	void recordTKE();
};

extern void readPhaseParams(Front*);
extern void initPhaseIntfc(char*,int,LEVEL_FUNC_PACK*,PARAMS*);
extern void melting_point_propagate(Front*,POINTER,POINT*,POINT*,
                    HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
extern void read_crt_dirichlet_bdry_data(char*,Front*,F_BASIC_DATA);
extern void read_melt_dirichlet_bdry_data(char*,Front*,F_BASIC_DATA);
extern void melt_flowThroughBoundaryState(double*,HYPER_SURF*,Front*,
			POINTER,POINTER);
extern void read_fluid_params(Front*);
extern void init_fluid_state_func(Front*,Incompress_Solver_Smooth_Basis*);		
extern void assignStateTemperature(double,POINTER);
extern void assignStateVapor(double,POINTER);
extern double getStateTemperature(POINTER);
extern double getStateVapor(POINTER);
extern double getStateSuper(POINTER);
extern double jumpT(POINTER,int,double*);
extern double jumpEpsGradDotNorm(POINTER,int,double*,double*);
extern double jumpGradDotTan(POINTER,int,int,double*,double*);
extern void initWaterDrops(Front*);
extern void compute_ice_particle_force(Front*,HYPER_SURF*,double, double*, double*);
extern void CondensationPreAdvance(Front*);
extern void ParticlePropagate(Front*);
extern void read_CL_prob_type(Front*);
extern void readWaterDropsParams(Front*,char*);
extern void printDropletsStates(Front*,char*);
/*plot functions*/
extern void gv_plot_scatter(Front*);
extern void vtk_plot_scatter(Front*);
extern void vtk_plot_sample_traj(Front*);
/*Statistics functions*/
extern void  Deviation(PARTICLE*,int,double&,double&);
extern double* ComputePDF(double*,int,double,int&,double&,double&);

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
#include <fftw3.h>
#include <vector>
#include <petscksp.h>
#include <assert.h>

#define         LIQUID_COMP2           3
#define         SOLID_COMP             1
#define		alternate_comp(comp) 					\
		(comp) == LIQUID_COMP2 ? SOLID_COMP : LIQUID_COMP2

enum _CL_PROB_TYPE {
	PARTICLE_TRACKING = 1,
	ENTRAINMENT = 1,
	RANDOM_FIELD
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
	CONST_STATE,
	RAND_STATE,
	TAYLOR_STATE,
	PRESET_STATE,
	FOURIER_STATE,
	LR_STATE, /*left and right state*/
	TB_STATE  /*top and bottom state*/
};
typedef enum _INIT_STATE INIT_STATE;

struct _PHASE_FIELD {
        double *temperature;
	double *vapor;
	double *cloud;
	double *supersat;
	double *pres;
	double *mrad;  /*mean radius in a cell*/
	double *drops; /*number of droplets in a cell*/
        double **vel;
	double **ext_accel;
};
typedef struct _PHASE_FIELD PHASE_FIELD;

struct _MOVIE_OPTION {
        boolean plot_pres;
        boolean plot_vort;
        boolean plot_velo;
        boolean plot_temperature;
	boolean plot_vapor;
	boolean plot_particles;
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
	INIT_STATE init_vapor_state;
	INIT_STATE init_drop_state;
	PHASE_FIELD *field;
	MOVIE_OPTION *movie_option;
	int pde_order;
	int num_phases;
	int num_drops;
	double T0;/*initial temperature */
	double qv0;/*initial vapor mixing ratio*/
	double qs; /*saturated vapor mixing ratio*/
	double D;    /*molecular diffusivities of the temperature*/
	double rho_l; /*density of water droplet*/
	double Lh;/*laten heat*/
	double Rv;/*individual gas constant for water vapor*/
	double Rd;/*individual gas constant for dry air*/
	double Kc;/*coefficient of thermal conductivity of air*/
	double Cp;/*specific heat with pressure held constant*/
	double K; /*coefficient for condensation*/
	double L[MAXD];
	double U[MAXD]; /*subdomain for particles*/
	boolean no_droplets;
	boolean if_condensation;
	boolean if_sedimentation;
	boolean if_macroscopic;
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
	void augmentMovieVariables(const char*);
	void checkOutput();
	void checkField();
	void recordField(char *, const char *);
	void recordPDF(char *, const char *);
	void recordRadius(char *);
	void recordMixingLine();
	void recordClusteringIndex();
	void recordSampleParticles();
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
	void (*getInitialState)(COMPONENT*,double*,PHASE_FIELD*,int,int,PARAMS*);
	void setParallelVapor(void);
	void initPresetParticles();

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

class MACRO{
	Front *front;
	PARAMS *params;
	double T; /*Macro temperature*/
	double Q; /*Macro vapor mixing ratio*/
	double Cd; /*Macro Condensation rate*/
	double P; /*Macro pressure*/
	double S; /*Macro supersaturation*/
	double Wm; /*Mean vertical velocity*/
	double rho0; /*density for ambient air*/
	double Gamma; /*dry-adiabatic lapse rate*/
	int dim;
public: 
	~MACRO();
	MACRO(Front &front);
	void computeCondRate();
	void computeTemp();
	void computeVapor();
	void computeSupersat();
	void computePressure();
	void solve();
	void output();
	void ParticleShift();
	void initMovieVariables();
	double getSupersat();
	double getVapor();
	double getPressure();
	double getTemp();
	void   setW(double);
};

extern void readPhaseParams(Front*);
extern void initPhaseIntfc(char*,int,LEVEL_FUNC_PACK*,PARAMS*);
extern void read_crt_dirichlet_bdry_data(char*,Front*,F_BASIC_DATA);
extern void read_melt_dirichlet_bdry_data(char*,Front*,F_BASIC_DATA);
extern void melt_flowThroughBoundaryState(double*,HYPER_SURF*,Front*,
			POINTER,POINTER);
extern void read_fluid_params(Front*);
extern void init_fluid_state_func(Front*,Incompress_Solver_Smooth_Basis*);	
extern void init_vapor_state_func(Front*,VCARTESIAN*);	
extern void assignStateTemperature(double,POINTER);
extern void assignStateVapor(double,POINTER);
extern double getStateTemperature(POINTER);
extern double getStateVapor(POINTER);
extern double getStateSuper(POINTER);
extern void initWaterDrops(Front*);
extern void compute_ice_particle_force(Front*,HYPER_SURF*,double, double*, double*);
extern void CondensationPreAdvance(Front*);
extern void ParticlePropagate(Front*,MACRO*);
extern void setParticleGlobalIndex(PARTICLE*,int);
extern void read_CL_prob_type(Front*);
extern void readWaterDropsParams(Front*,char*);
extern void printDropletsStates(Front*,char*);
/*plot functions*/
extern void gv_plot_scatter(Front*);
extern void vtk_plot_scatter(Front*);
extern void vtk_plot_sample_traj(Front*);
/*Statistics functions*/
extern void  Deviation(double*,int,double&,double&);
extern double* ComputePDF(double*,int,double&,int,double&,double&);
extern bool  fftnd(fftw_complex*, int, int*,int);

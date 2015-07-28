/**********************************************************************
 * 		crystal.h					      *
 **********************************************************************/

#ifndef _FT_CRYSTAL_H_
#define _FT_CRYSTAL_H_

#include <FronTier.h>
#if defined(IFLUID_STATE)
#include "ifluid_state.h"
#else /* defined(IFLUID_STATE) */
#include "crystal_state.h"
#endif /* defined(IFLUID_STATE) */
#include <vector>
#include <petscksp.h>
#include <assert.h>

#define         CRYSTAL_COMP            0
#define         SOLUTE_COMP             3

enum _DF_SCHEME {
	UNSPLIT_EXPLICIT = 1,
	UNSPLIT_IMPLICIT,
	CRANK_NICOLSON,
	CONSERVATIVE_EXPLICIT,
	CONSERVATIVE_IMPLICIT,
	CONSERVATIVE_CRANK_NICOLSON
};
typedef enum _DF_SCHEME DF_SCHEME;

struct _CRT_FIELD {
        double *solute;
	double **vel;
	double *source;
	CELL_PART *cell_part;
};
typedef struct _CRT_FIELD CRT_FIELD;

enum _POINT_PROP_SCHEME {
        EXPLICIT_EULER = 1,
        IMPLICIT_EULER,
        MIDDLE_POINT,
	CONSTANT_STATE
};
typedef enum _POINT_PROP_SCHEME POINT_PROP_SCHEME;

enum _REACTION_TYPE {
        DEPOSITION_ONLY = 1,
        DISSOLUTION_ONLY,
	REVERSIBLE_REACTION
};
typedef enum _REACTION_TYPE REACTION_TYPE;

struct _CRT_PARAMS {
        int dim;
	DF_SCHEME num_scheme;
	POINT_PROP_SCHEME point_prop_scheme;
	boolean add_curvature;
	REACTION_TYPE reaction_type;
	CRT_FIELD *field;	// field of solute concentration
	int pde_order;		// Order of Runge-Kutta
	double C_0;		// Ambient concentration
        double C_eq;    /* solute concentration in equilibrium with solid */
        double rho_s;   /* density of the precipitated solid phase */
        double D;       /* diffusion coefficient of the solute concentration */
        double k;       /* local reaction rate coefficient */
	double gap;	/* initial gap between ambient solute and crystal */
	double max_solute,min_solute;
	double (*func)(POINTER,double*);
        POINTER func_params;
	double (*crystal_dens_func)(POINTER,double*);
        POINTER crystal_dens_params;
};
typedef struct _CRT_PARAMS CRT_PARAMS;

/*	Crystal density functions */
struct _HALF_MOON_PARAMS {
	int dir;
	double cutline;
	double lower_dens;
	double upper_dens;
};
typedef struct _HALF_MOON_PARAMS HALF_MOON_PARAMS;

struct _FUEL_SAMPLE_PARAMS {
        Front *front;
	double rho_1;
	double rho_2;
	double rho_3;
};
typedef struct _FUEL_SAMPLE_PARAMS FUEL_SAMPLE_PARAMS;

#define		is_conservative(num_scheme) 			\
		((num_scheme) == CONSERVATIVE_EXPLICIT ||	\
		 (num_scheme) == CONSERVATIVE_IMPLICIT ||	\
		 (num_scheme) == CONSERVATIVE_CRANK_NICOLSON)

extern double   half_moon_density(POINTER,double*);
extern double   perturbed_density(POINTER,double*);
extern double   fuel_sample_density(POINTER,double*);

/*	Crystal curve functions */
extern double   crystal_seed_curve(POINTER,double*);

typedef class C_CARTESIAN C_CARTESIAN_EB;

/******************************************************************************
 * 		lcartsn.h
 * A simple incompressible flow solver using the ghost fluid method and the
 * projection method.
 *
 * the main function is 
 * 	C_CARTESIAN::solve().
 *
 * References:
 ******************************************************************************/

class SOLVER;
class C_CARTESIAN;
class PARABOLIC_SOLVER;

//------------------------------------------------------
//		MESH
//------------------------------------------------------

class C_RECTANGLE {
public:
	int index;			// rectangle index
	int comp;			 
	double area;
	double coords[MAXD];	
	int icoords[MAXD];

	C_RECTANGLE();

	void setCoords(double*,int);
};


class C_CARTESIAN{
	Front *front;
public:
	C_CARTESIAN(Front &front);

	// member data: RECT_GRID
	int dim;

	// On topological grid
	RECT_GRID *top_grid;
	double *array;		// for scatter states;
	double *top_L,*top_U,*top_h,hmin;
	int *top_gmax;
	COMPONENT *top_comp;
	COMPONENT *top_comp_old;	// Save for new time step
	CRT_PARAMS *cRparams;
	CRT_FIELD *field;

	int *lbuf,*ubuf,*gmax;
	int *i_to_I,*I_to_i;		// Index mapping for 1D
	int **ij_to_I,**I_to_ij;	// Index mapping for 2D
	int ***ijk_to_I,**I_to_ijk;	// Index mapping for 3D
	double **time_data;
	double accum_influx;
	double total_solute0;
	double total_matter0;
	double total_crystl0;

	// Sweeping limites
	int imin,jmin,kmin;
	int imax,jmax,kmax;
	int mesh_size;

	// member data: mesh storage
	std::vector<C_RECTANGLE> 	cell_center;

	double m_t;                     // time
	double m_dt;			// time increment
        double max_dt;

	// constructor
	~C_CARTESIAN();

	// for parallel partition
	int             NLblocks,ilower,iupper;
        int             *n_dist;

	// mesh: full cells mesh
	void initMesh(void);		// setup the cartesian grid
	void setComponent(void);	// init components	
	void setAdvectionDt(void);	// compute time step 	
	void setBoundary(void);		// set up boundary conditions 	
	void printFrontInteriorStates();
	void copySolute(void);

	void computeAdvection(void);

	void computeAdvectionIM1d(void);
	void computeAdvectionIM2d(void);
	void computeAdvectionIM3d(void);
	void conservativeIntfcCrossing(int);
	void computeSource(void);
	void computeCellSource(int,double*);

	// interface functions
	void setDomain();

	// Extra plot functions
	void crystalDraw();
	void oneDimPlot(char*);
	void xgraphOneDimPlot(char*);
	void vtk_plot_concentration2d(char*);
	void initSampleSolute(char *in_name);
        void sampleSolute();
        void sampleSolute2d();
        void sampleSolute3d();
        void timeStepAnalysis(boolean);
        double computeTotalMass(double*,double*,double*,double*,boolean);

	// Extra movie functions
        void initMovieVariables(void);
	void setInitialCondition(void);

	void sampleStates(SAMPLE);

	// parallelization related functions
	//
	void scatMeshArray();
	void setGlobalIndex();

	// physics calculation
	void readFrontInteriorStates(char*);

	void setIndexMap(void);
		// for compProjWithSmoothProperty(), 
		// should be changed to use setIndexMap() only.

	// -------------------------------------------------------
	// 		incompressible solver functions

	void getInitialState(C_RECTANGLE&); 

	// main step function
	void solve(double dt);		

	// ----------------------------------------------------------
	// 		utility functions
	// ----------------------------------------------------------

	void getRectangleIndex(int indexRectangle, int &i, int &j);
	void getRectangleIndex(int indexRectangle, int &i, int &j, int &k);
	void getRectangleIndex(int indexRectangle, int*);
	void getRectangleCenter(int index, double *coords);
	void getRectangleCenter(int index0, int index1, double *coords);
	
	void save(char *filename);
};

//crysub.cpp
extern void 	solute_print_front_states(FILE*,Front*);
extern void 	solute_read_front_states(FILE*,Front*);
extern double   getStateSolute(POINTER);
extern void	crystal_point_propagate(Front*,POINTER,POINT*,POINT*,	
			HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
extern void	read_crystal_params(char*,CRT_PARAMS*);
extern void	intfc_1d_data(Front*,double*,double*);
extern void 	read_seed_params(int,FILE*,SEED_PARAMS*);
extern void 	read_crt_dirichlet_bdry_data(char*,Front*,F_BASIC_DATA);
extern void 	fractal_dimension(Front*,SEED_PARAMS*,double*,double*);
extern void	setInitialIntfc(Front*,LEVEL_FUNC_PACK*,char*);
extern void	initFrontStates(Front*);
extern void 	read_restart_params(int,char*,Front*);
extern boolean  bdryReached(Front*);
#endif

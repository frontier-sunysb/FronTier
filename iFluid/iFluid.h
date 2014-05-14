/**********************************************************************
 * 		iFluid.h
 **********************************************************************/

#ifndef _FT_IFLUID_H_
#define _FT_IFLUID_H_

#include <vector>
#include <petscksp.h>
#include <assert.h>

#include "FronTier.h"
#include "solver.h"

#define         SOLID_COMP		1
#define         LIQUID_COMP1		2
#define         LIQUID_COMP2		3
#define		LIQUID_COMP		3

#define		ifluid_comp(comp)   (((comp) == LIQUID_COMP1 || 	\
		comp == LIQUID_COMP2) ? YES : NO)

enum EBM_COORD
{
    COORD_X = 0,  COORD_Y = 1,  COORD_Z = 2
};

enum _DOMAIN_STATUS {
        NOT_SOLVED      =       0,
        TO_SOLVE,
        SOLVED
};
typedef enum _DOMAIN_STATUS DOMAIN_STATUS;

struct _STATE {
        double dens;                    /* Density */
        double pres;                    /* Pressure */
        double phi;                     /* Potential */
        double vel[MAXD];               /* Velocities */
        double vort;                    /* Vorticity in 2D */
        double impulse[MAXD];            /* Accum impact from external force */
	double solute;			/* For subsurface problem */
	double temperature;             /* For melting with flow problem */
	double vapor;                   /* For climate problem */
        double supersat;		/* For climate problem */
};
typedef struct _STATE STATE;

struct _IF_FIELD {
	double **vel;			/* Velocities */
	double *temperature;            /* Temperature */
	double *phi;
	double *q;
	double *pres;			/* Pressure */
	double *vort;			/* Vorticity in 2D */
	double *mu;
	double *rho;
	double **grad_q;
	double **f_surf;		// Surface force (such as tension)

	double *div_U;
	double *d_phi;			/* Dual grid phi */
	double *nu_t;			/* Turbulent viscosity */
};
typedef struct _IF_FIELD IF_FIELD;

enum _PROJC_METHOD {
	ERROR_PROJC_SCHEME		= -1,
        SIMPLE			=  1,
        BELL_COLELLA,
        KIM_MOIN,
        PEROT_BOTELLA
};
typedef enum _PROJC_METHOD PROJC_METHOD;

enum _ADVEC_METHOD {
	ERROR_ADVEC_SCHEME		= -1,
        UPWIND			=  1,
        WENO
};
typedef enum _ADVEC_METHOD ADVEC_METHOD;

enum _ELLIP_METHOD {
	ERROR_ELLIP_SCHEME		= -1,
	SIMPLE_ELLIP		= 1,
	DUAL_ELLIP,
	DOUBLE_ELLIP,
	CIM_ELLIP
};
typedef enum _ELLIP_METHOD ELLIP_METHOD;

enum _DOMAIN_METHOD {
	BY_COMPONENT		= 1,
	BY_CROSSING,
	BY_WALL
};
typedef enum _DOMAIN_METHOD DOMAIN_METHOD;

struct _NS_SCHEME {
	PROJC_METHOD projc_method;
	ADVEC_METHOD advec_method;
	ELLIP_METHOD ellip_method;
};
typedef struct _NS_SCHEME NS_SCHEME;

typedef struct {
        int dim;
        POINTER level_func_params;
	NS_SCHEME num_scheme;
        double rho1;
        double rho2;
	double mu1;
	double mu2;
	double U1[MAXD];
	double U2[MAXD];
	double gravity[MAXD];
	double U_ambient[MAXD];
	double surf_tension;
	double smoothing_radius;
	double ub_speed;
	COMPONENT m_comp1;
	COMPONENT m_comp2;
	IF_FIELD *field;
	int adv_order;
	boolean total_div_cancellation;
	boolean buoyancy_flow;
	boolean if_buoyancy;
	double  ref_temp;
	boolean if_ref_pres;
	boolean use_eddy_visc;	/* Yes if to use eddy viscosity */
	double  ref_pres;
	double  Amplitute; /*Amplitute of velocity*/
} IF_PARAMS;

struct _FLOW_THROUGH_PARAMS {
        POINT *oldp;
        COMPONENT comp;
};
typedef struct _FLOW_THROUGH_PARAMS FLOW_THROUGH_PARAMS;

enum _TIME_FUNC_TYPE {
	CONSTANT		=  1,
	PULSE_FUNC,
	SINE_FUNC	
};
typedef enum _TIME_FUNC_TYPE TIME_FUNC_TYPE;

struct _TIME_DEPENDENT_PARAMS {
	TIME_FUNC_TYPE td_type;
	double v_base[MAXD],p_base;
	double v_peak[MAXD],p_peak;
	double v_tail[MAXD],p_tail;
	double v_amp[MAXD],p_amp;
	double omega,phase;
	double T[10];
};
typedef struct _TIME_DEPENDENT_PARAMS TIME_DEPENDENT_PARAMS;

/******************************************************************************
 * 		lcartsn.h
 * A simple incompressible flow solver using the ghost fluid method and the
 * projection method.
 *
 * the main function is 
 * 	L_CARTESIAN::solve().
 *
 * References:
 ******************************************************************************/

class SOLVER;
class Incompress_Solver_Basis;

class L_RECTANGLE {
public:
	int comp;			 
	double m_coords[MAXD];	
	int icoords[MAXD];

	L_RECTANGLE();

	void setCoords(double*,int);
};

class Incompress_Solver_Basis{
public:
       	Incompress_Solver_Basis() {}; // constructor
	virtual ~Incompress_Solver_Basis() {};

};

class Incompress_Solver_Smooth_Basis:public Incompress_Solver_Basis{
public:
        //constructor
	Incompress_Solver_Smooth_Basis(Front &front);
	virtual ~Incompress_Solver_Smooth_Basis() {};

	double m_dt;
	double accum_dt;
	double max_speed;
	double min_pressure;
        double max_pressure;
        double min_value; //for debugging
	double max_value; //for debugging
	double max_dt;
	double min_dt;
	double *top_h;
	double vmin[MAXD],vmax[MAXD];
	int dim;
	int icrds_max[MAXD];

	void initMesh(void);
	void setAdvectionDt(void); 
			//using max speed and hmin to determine max_dt, min_dt
	void readFrontInteriorStates(char *state_name);
	void printFrontInteriorStates(char *state_name);
	void initMovieVariables(void);
	void getVelocity(double *p, double *U);
	void initSampleVelocity(char *in_name);

	//Initialization of States
	void (*getInitialState) (COMPONENT,double*,IF_FIELD*,int,int,
				IF_PARAMS*);
	int (*findStateAtCrossing)(Front*,int*,GRID_DIRECTION,int,
				POINTER*,HYPER_SURF**,double*);
	int (*findStateAtCGCrossing)(Front*,int*,GRID_DIRECTION,int,
				POINTER*,HYPER_SURF**,double*);
	void applicationSetComponent();
	void applicationSetStates();

	//For debugging test
	void solveTest(const char *msg);

	//User interface
	virtual void setInitialCondition(void) = 0;
	virtual void setParallelVelocity(void) = 0;
	virtual void solve(double dt) = 0; // main step function
        virtual void vtk_plot_scalar(char*, const char*) = 0;

protected:
	Front *front;
	IF_PARAMS *iFparams;
	IF_FIELD  *field;
	
	// On dual topological grid
	RECT_GRID *top_grid;
	double *array;
	double *source;
	double *diff_coeff;
	COMPONENT *top_comp;
	int *top_gmax;
	int *lbuf, *ubuf;
	double *top_L, *top_U;
	int **ij_to_I, **I_to_ij;
	int ***ijk_to_I, **I_to_ijk;
	int *domain_status;
	int smin[MAXD],smax[MAXD];
	// Sweeping limites
	int imin, jmin, kmin;
	int imax, jmax, kmax;
	// for parallel partition
	int NLblocks, ilower, iupper;
	int *n_dist;
	// for dual/comp overlapping
	int offset[MAXD];

	// On comp topological grid
	RECT_GRID *ctop_grid;
	double *carray;
	double *csource;
	COMPONENT *ctop_comp;
	int *ctop_gmax;
	int *clbuf, *cubuf;
	double *ctop_L, *ctop_U;
	int **cij_to_I, **cI_to_ij;
	int ***cijk_to_I, **cI_to_ijk;
	int *cdomain_status;
	int csmin[MAXD],csmax[MAXD];
	// Sweeping limites
	int cimin, cjmin, ckmin;
	int cimax, cjmax, ckmax;
	// for parallel partition
	int cNLblocks, cilower, ciupper;
	int *cn_dist;

	// Index shift between dual and comp grids 
	int ishift[MAXD];

	//member data: mesh storage
	std::vector<L_RECTANGLE>   cell_center;

	//member data:
	int    m_comp[2];
	double m_mu[2];
	double m_rho[2];// two components at most
	double m_sigma; //surface tension
	double m_smoothing_radius;// used by smoothing function

	double hmin; //smallest spacing
	double mu_min; //smallest viscocity
	double rho_min;// smallest density
	double m_t;
	double m_t_old, m_t_int, m_t_new;

protected:
	void setComponent(void); //init components;
	void setDomain();
	void setDualDomain();

	// parallelization related functions
	void scatMeshArray(void);
	void setGlobalIndex(void);
	void setDualGlobalIndex(void);
	void setIndexMap(void);
	void setDualIndexMap(void);
	void paintAllGridPoint(int status);
	void paintSolvedGridPoint();
	void setReferencePressure();
	boolean paintToSolveGridPoint();
	boolean nextConnectedPoint(int*,GRID_DIRECTION,int*,int,int*,int*);

/****************  Functions related to solve() *********/

	virtual void copyMeshStates(void) = 0;
	virtual void computeGradientQ(void) = 0;
	virtual void surfaceTension(double*, HYPER_SURF_ELEMENT*,
		HYPER_SURF*, double*, double) = 0;

	
/***********************  Utility functions  *******************/

	void   getRectangleIndex(int indexRectangle, int &i, int &j);
	void   getRectangleIndex(int indexRectangle, int &i, int &j, int &k);
	int    getRectangleComponent(int index);	// the center component
	void   getRectangleCenter(int index, double *coords);
	double getDistance(double *coords0, double *coords1);
	int    getComponent(int *icoords);	
	int    getComponent(double *coords);	
	void   save(char *filename);
	double computeFieldPointDiv(int*, double**);
	double computeDualFieldPointDiv(int*, double**);
	double computeDualMu(int*, double*);
	double computeFieldPointMuTurb(int*, double);
	void   computeFieldPointGrad(int*, double*, double*);
	void   computeDualFieldPointGrad(int*, double*, double*);
	void   checkVelocityDiv(const char*);
	void   computeDualFieldPointrho(int*);
/************* TMP Functions which are not implemented or used ***********/

	void computeSubgridModel(void);    // subgrid model by Hyunkyung Lim
	void getNearestInterfacePoint(COMPONENT,double*,double*,double*,
					double*); 
};

///////////////Interface for Embedded Boundary Method////////////////////

class Incompress_Solver_EBM:public Incompress_Solver_Basis{
public:
        Incompress_Solver_EBM(Front &front) {};//constructor
	~Incompress_Solver_EBM() {};
};
/////////////////////////////////////////////////////////////////////////////////
//
class Incompress_Solver_Smooth_2D_Basis:
public 	Incompress_Solver_Smooth_Basis{
public:
        Incompress_Solver_Smooth_2D_Basis(Front &front):
	Incompress_Solver_Smooth_Basis(front) {};
	virtual ~Incompress_Solver_Smooth_2D_Basis() {};

protected:
	double getSmoothingFunction(double r);
	double getSmoothingFunctionD(double*, double*);
	double smoothedDeltaFunction(double*, double*);
	double smoothedStepFunction(double*, double*, int);
	void sampleVelocity();
	void setSmoothedProperties(void);
};


class Incompress_Solver_Smooth_3D_Basis:
public 	Incompress_Solver_Smooth_Basis{
public:
        Incompress_Solver_Smooth_3D_Basis(Front &front):
	Incompress_Solver_Smooth_Basis(front) {};
	virtual ~Incompress_Solver_Smooth_3D_Basis() {};

protected:
	double getSmoothingFunction(double r);
	double getSmoothingFunctionD(double*, double*);
	double smoothedDeltaFunction(double*, double*);
	double smoothedStepFunction(double*, double*, int);
	void sampleVelocity();
	void setSmoothedProperties(void);
};

class Incompress_Solver_Smooth_2D_Cartesian:
public 	Incompress_Solver_Smooth_2D_Basis{
public:
        Incompress_Solver_Smooth_2D_Cartesian(Front &front):
	Incompress_Solver_Smooth_2D_Basis(front) {};
	~Incompress_Solver_Smooth_2D_Cartesian() {};

	void setInitialCondition(void);
	void setParallelVelocity(void);
	void solve(double dt);
        void vtk_plot_scalar(char*, const char*);
protected:
	void copyMeshStates(void);
	void computeAdvection(void);
	void computeDiffusion(void);
	void computeDiffusionCN(void);
	void computeDiffusionExplicit(void);
	void computeDiffusionImplicit(void);
	void computeDiffusionParab(void);
	void computeProjection(void);
	void computeProjectionCim(void);
	void computeProjectionSimple(void);
	void computeProjectionDouble(void);
	void computeProjectionDual(void);
	void computePressure(void);
	void computePressurePmI(void);
	void computePressurePmII(void);
	void computePressurePmIII(void);
	void computeGradientQ(void);
	void computeNewVelocity(void);
	void computeNewVelocityDual(void);
	void extractFlowThroughVelocity(void);
	void computeSourceTerm(double *coords, double *source);
	void surfaceTension(double*, HYPER_SURF_ELEMENT*, HYPER_SURF*, 
				double*, double);
	void updateComponent(void);

	/***************   Low level computation functions  *************/
	double getVorticity(int i, int j);
};

class Incompress_Solver_Smooth_3D_Cartesian:
public 	Incompress_Solver_Smooth_3D_Basis{
public:
        Incompress_Solver_Smooth_3D_Cartesian(Front &front):
	Incompress_Solver_Smooth_3D_Basis(front) {};
	~Incompress_Solver_Smooth_3D_Cartesian() {};

	void setInitialCondition(void);
	void setParallelVelocity(void);
	void solve(double dt);
	void solveTest(const char *msg);
        void vtk_plot_scalar(char*, const char*);
protected:
	void copyMeshStates(void);
	void computeAdvection(void);
	void computeDiffusion(void);
	void computeDiffusionCN(void);
	void computeDiffusionExplicit(void);
	void computeDiffusionImplicit(void);
	void computeDiffusionParab(void);
	void computeProjection(void);
	void computeProjectionCim(void);
	void computeProjectionSimple(void);
	void computeProjectionDouble(void);
	void computeProjectionDual(void);
	void computePressure(void);
	void computePressurePmI(void);
	void computePressurePmII(void);
	void computePressurePmIII(void);
	void computeGradientQ(void);
	void computeNewVelocity(void);
	void computeNewVelocityDual(void);
	void updateComponent(void);
	boolean InsideSolid(int*);
	void extractFlowThroughVelocity(void);
	void computeSourceTerm(double *coords, double *source);
	void surfaceTension(double*, HYPER_SURF_ELEMENT*, HYPER_SURF*, 
				double*, double);
};

extern double getStatePres(POINTER);
extern double getStatePhi(POINTER);
extern double getStateVort(POINTER);
extern double getStateXvel(POINTER);
extern double getStateYvel(POINTER);
extern double getStateZvel(POINTER);
extern double getStateXimp(POINTER);
extern double getStateYimp(POINTER);
extern double getStateZimp(POINTER);
extern double getStateComp(POINTER);
extern double getPressure(Front*,double*,double*);
extern double getPhiFromPres(Front*,double);
extern double burger_flux(double,double,double);
extern double linear_flux(double,double,double,double);
extern void fluid_print_front_states(FILE*,Front*);
extern void fluid_read_front_states(FILE*,Front*);
extern void read_iF_dirichlet_bdry_data(char*,Front*,F_BASIC_DATA);
extern boolean isDirichletPresetBdry(Front*,int*,GRID_DIRECTION,COMPONENT);
extern int ifluid_find_state_at_crossing(Front*,int*,GRID_DIRECTION,
			int,POINTER*,HYPER_SURF**,double*);
extern int ifluid_find_state_at_cg_crossing(Front*,int*,GRID_DIRECTION,
			int,POINTER*,HYPER_SURF**,double*);
extern int ifluid_find_state_at_dual_crossing(Front*,int*,GRID_DIRECTION,
			int,POINTER*,HYPER_SURF**,double*);
extern double p_jump(POINTER,int,double*);
extern double grad_p_jump_n(POINTER,int,double*,double*);
extern double grad_p_jump_t(POINTER,int,int,double*,double*);
extern boolean neumann_type_bdry(int);

#endif

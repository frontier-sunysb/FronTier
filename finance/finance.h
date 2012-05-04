/**********************************************************************
 * 				finance.h			      *
 **********************************************************************/

#include <FronTier.h>
#include <vector>
#include <petscksp.h>
#include <assert.h>

#define         EXCERCISE_COMP            	2
#define         BLACK_SCHOLES_COMP              3

#define         bs_comp(comp)   ((comp) == BLACK_SCHOLES_COMP) ? YES : NO

struct _FIELD {
        double *option_price;                   /* Pressure */
};
typedef struct _FIELD FIELD;

struct _F_MOVIE_OPTION {
        boolean plot_price;
	boolean plot_cross_section[MAXD];
};
typedef struct _F_MOVIE_OPTION F_MOVIE_OPTION;

enum _F_TYPE {
	ERROR_TYPE = -1,
	EURO_CALL_OPTION = 1,
	EURO_PUT_OPTION = 2,
	AMRI_CALL_OPTION = 3,
	AMRI_PUT_OPTION = 4
};
typedef enum _F_TYPE F_TYPE;

enum _V_TYPE {
	CONSTANT_V = 0,
	LINEAR_V = 1,
	QUADRATIC_V = 2,
	SINE_V = 3,
	HESTON_V = 4
};
typedef enum _V_TYPE V_TYPE;

enum _NUM_SCHEME {
	UNSPLIT_IMPLICIT = 1,
	CRANK_NICOLSON
};
typedef enum _NUM_SCHEME NUM_SCHEME;

struct _PARAMS {
	int dim;
	F_TYPE f_type;
	V_TYPE v_type;
	POINTER vparams;
	NUM_SCHEME num_scheme;
	F_MOVIE_OPTION	*movie_option;
	double *rate_of_change;	// Option price on stock value
	double *option_price;	// Option price on stock value
	double *temp_price;     // intermediate price while looking for SF
	double E;		// Excercise price
	double sigma[MAXD];	// Volotility
	double r;		// Interest Rate
	double D;               // Dividend
	double oldroot;
	double a;
	double b;
	double c;
	double sample_S;
        int idx;                // index for current front
	bool findSF;
	FIELD *field;
};
typedef struct _PARAMS PARAMS;

enum _TD_TYPE {
	UNKNOWN_TD_TYPE = 0,
	EXP_DECAY
};
typedef enum _TD_TYPE TD_TYPE;

struct _TIME_DEPENDENT_PARAMS {
	TD_TYPE type;
	double E;
	double D;
	double r;
	double decay_rate;
};
typedef struct _TIME_DEPENDENT_PARAMS TIME_DEPENDENT_PARAMS;

typedef double STATE;

extern double   price_of_state(POINTER);
extern double	extend_from_put_exc(double,double*,double*,double*,double);
extern double	extend_from_call_exc(double,double*,double*,double*,double);

typedef class CARTESIAN CARTESIAN_EB;

/******************************************************************************
 * 		lcartsn.h
 * A simple incompressible flow solver using the ghost fluid method and the
 * projection method.
 *
 * the main function is 
 * 	CARTESIAN::solve().
 *
 * References:
 ******************************************************************************/

class SOLVER;
class CARTESIAN;

//typedef class CARTESIAN CARTESIAN_EB;
//enum VISITED_TYPE {UNVISITED, VISITED, PARTIAL_VISITED, FULL_VISITED};

//-------------------------------------------------
//		STATES
// NOTE:
//      INC_STATE/INC_STATE_RECT_EDGEshould be put into 
// lcartsn.h. However, there are some trouble
// to compile in that way.
//-------------------------------------------------
// states inside a cell

class INC_STATE{
public:
	double P;		// Option price

	double sigma[MAXD];	// Volotility 
	double r;		// Interest rate

	void setZero(void);
};
// states on edge

//------------------------------------------------------
//		MESH
//------------------------------------------------------

class RECTANGLE {
public:
	int index;			// rectangle index
	int comp;			 
	INC_STATE state;
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
	double *top_L,*top_U,*top_h;
	COMPONENT *top_comp;
	PARAMS *eqn_params;
	FIELD *field;

	int *top_gmax;
	int *lbuf,*ubuf;
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

	// constructor
	~CARTESIAN();

	// for parallel partition
	int             NLblocks,ilower,iupper;
        int             *n_dist;

	// mesh: full cells mesh
	void initMesh(void);		// setup the cartesian grid
	void setComponent(void);	// init components	
	void setAdvectionDt(void);	// compute time step 	
	void readOptionPrice(char*);
	void printOptionPrice(char*);
	void printTimeData(char*);
	void copyOptionPrice(void);
	void copyBackOptionPrice(void);
	void initMovieVariables(void);

	void computeAdvection(void);

	void computeAdvectionExplicit(void);
	void computeAdvectionExplicit1d(void);
	void computeAdvectionExplicit2d(void);
	void computeAdvectionExplicit3d(void);

	void computeAdvectionImplicit(void);
	void computeAdvectionImplicit1d(void);
	void computeAdvectionImplicit2d(void);
	void computeAdvectionImplicit3d(void);

	void checkStates();

	// parallelization related functions
	//
	void scatMeshArray();
	void setGlobalIndex();

	// physics calculation
	void setInitialCondition(void);

	void setDomain();
	void setIndexMap(void);
		// for compProjWithSmoothProperty(), 
		// should be changed to use setIndexMap() only.

	// -------------------------------------------------------
	// 		incompressible solver functions

	void computeSourceTerm(double *coords, double t, INC_STATE &state); 

	void getInitialState(double *coords, INC_STATE &state); 
	void getInitialState(double *coords, STATE *state); 

	// main step function
	void solve(double dt);		

	void getVelocity(double *p, double *U);

	// ----------------------------------------------------------
	// 		utility functions
	// ----------------------------------------------------------

	void getRectangleCenter(int index, double *coords);
	
	int  getComponent(int *icoords);	
	void save(char *filename);
};


extern void   read_movie_options(char*,PARAMS*);
extern void   excercise_point(INTERFACE*,double*,double*);
extern double getStatePrice(POINTER);
extern double linear_extension(double,double,double,double,double,double);
extern void   initVolatility(char*,Front*);
extern void   getVolatility(Front*);

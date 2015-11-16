/**********************************************************************
 * 		solver.h
 **********************************************************************/

#ifndef _FT_IFLUID_SOLVER_H_
#define _FT_IFLUID_SOLVER_H_

#include <FronTier.h>
#include <vector>
#include <petscksp.h>
#include <petscpc.h>
#include <assert.h>

enum {
        NO_PDE_BOUNDARY                 = 0,
        CONST_V_PDE_BOUNDARY            = 1,
        CONST_P_PDE_BOUNDARY,
        NEUMANN_PDE_BOUNDARY,
        DIRICHLET_PDE_BOUNDARY,
        MOVING_BOUNDARY,
        MIXED_PDE_BOUNDARY
};

class SOLVER
{
public:
	SOLVER(){};
	SOLVER(int ilower, int iupper, int d_nz, int o_nz){};
	virtual ~SOLVER(){};
	virtual void Create(int ilower, int iupper, int d_nz, int o_nz){};
	
	virtual void Set_A(PetscInt i, PetscInt j, double val){};
						// A[i][j]=val;
	virtual void Add_A(PetscInt i, PetscInt j, double val){};
						// A[i][j]=A[i][j]+val;
	virtual void Set_x(PetscInt i, double val){};	// x[i]=val;
	virtual void Set_x(double *p){};		// x[i]=p[i];
	virtual void Add_x(PetscInt i, double val){};	// x[i]=x[i]+val;
	virtual void Get_x(double *p){};	// get the x from ij_x to p.		
	virtual void Get_x(double *p, int n, int *global_index){};
	virtual void Set_b(PetscInt i, double val){};	// b[i]=val;
	virtual void Set_b(double *b){};	
	virtual void Add_b(PetscInt i, double val){};	// b[i]=b[i]+val;

	virtual void SetMaxIter(int val){};	
	virtual void GetFinalRelativeResidualNorm(double *rel_resid_norm){};
	virtual void GetNumIterations(int *num_iterations){};

	virtual void Solve(void){};	
	virtual void Solve_withPureNeumann(void){};	
	virtual void Read_A(char *filename){};
	virtual void Print_A(char *filename){};
	virtual void Read_b(char *filename){};
	virtual void Print_b(char *filename){};
	virtual void Read_x(char *filename){};
	virtual void Print_x(char *filename){};
	virtual void test(void){};
};

class PETSc: public SOLVER
{
public:	
	MPI_Comm  comm;			// set to be MPI_COMM_WORLD.
	int iLower;
	int iUpper;			// global row range
	
	Vec x;      			/* approx solution, RHS*/
	Vec b;
  	Mat A;          		/* linear system matrix */
  	
  	KSP   ksp;          		/* Krylov subspace method context */
	PC    pc;
	MatNullSpace	nullsp;
		
	PetscErrorCode ierr;
	int its;			// numer of iterations;

public:
	PETSc();
	PETSc(int ilower, int iupper, int d_nz, int o_nz);		
		// global row range of A, x, b on this processor
	~PETSc();
	void Create(int ilower, int iupper, int d_nz, int o_nz);	
		// same as Hypre(int, int)
	void Create(MPI_Comm Comm, int ilower, int iupper, int d_nz, int o_nz);
		// same as Hypre(int, int)
	
	void Reset_A();				// Set A[i][j]=0.0;
	void Reset_b();
	void Reset_x();
	void Set_A(PetscInt i, PetscInt j, double val);	// A[i][j]=val;
	void Add_A(PetscInt i, PetscInt j, double val);	// A[i][j]=A[i][j]+val;
	void Get_row_of_A(PetscInt i, PetscInt *ncol, PetscInt **cols, double **row);
	void Set_x(PetscInt i, double val);		// x[i]=val;
	void Add_x(PetscInt i, double val);		// x[i]=x[i]+val;
	void Set_b(PetscInt i, double val);		// b[i]=val;
	void Add_b(PetscInt i, double val);		// b[i]=b[i]+val;
	void Get_x(double *p);		// get the x from ij_x to p.	
	void Get_b(double *p);		// get the b from ij_x to p.
	void Get_x(double *p, int n, int *global_index);
	
	void SetMaxIter(int val); 	// Set maximum number of iterations 
	void SetTol(double val);	// Set the convergence tolerance 
	void SetKDim(int k_dim);	
			// Set the maximum size of the Krylov space 
	void GetNumIterations(PetscInt *num_iterations);	
			// Return the number of iterations taken 
	void GetFinalRelativeResidualNorm(double *rel_resid_norm);
	void Solve(void);
	void Solve_GMRES(void);
	void Solve_BCGSL(void);
	void Solve_LU(void);
#if defined(__HYPRE__)
	void Solve_HYPRE(void);
#endif
	void Solve_withPureNeumann(void);
	void Solve_withPureNeumann_GMRES(void);
	void Solve_withPureNeumann_BCGSL(void);
	void Solve_withPureNeumann_HYPRE(void);
	void Solve_withPureNeumann_ML(void);
	virtual void Print_A(const char *filename);
        virtual void Print_b(const char *filename);
};

class PARABOLIC_SOLVER{
        Front *front;
public:
        PARABOLIC_SOLVER(Front &front);

	// Time step for one call
	double dt;
	double sub_dt;

        // On topological grid
        boolean first;
	int *i_to_I;
	int **ij_to_I;
	int ***ijk_to_I;
	int ilower;
	int iupper;
	int order;		/* order of Runge-Kutta */
	CELL_PART *cell_part;	/* Cell partition structure */

	COMPONENT soln_comp;
	COMPONENT obst_comp;
	double *var;		/* field variable of old step */
	double *soln;		/* field variable of new step */
	double *source;		/* source field */
	double **a;		/* advection field */
	double D;
	double *nu;             /* Variable diffusion coefficient */
	double var_obst;	/* default solution in obst_comp */
	double bdry_influx;
	double (*getStateVarFunc)(POINTER);
	void set_solver_domain(void);

	void solveEX(void);
	void solveCN(void);
	void solve(double *var_in,double *var_out);
	void solve1d(double *var_in,double *var_out);
	void solve2d(double *var_in,double *var_out);
	void solve3d(double *var_in,double *var_out);
	void solveIM(void);

	void solveCEX(void);
	void solveCCN(void);
	void solveCIM(void);
	void Csolve(double *var_in,double *var_out);
	void Csolve1d(double *var_in,double *var_out);
	void Csolve2d(double *var_in,double *var_out);
	void Csolve3d(double *var_in,double *var_out);

	double compBdryFlux(double*,double);
	void addMeshState(double *ans,double C1,double *var1,double C2,
				double *var2);
	int (*findStateAtCrossing)(Front*,int*,GRID_DIRECTION,int,
                                POINTER*,HYPER_SURF**,double*);
	int (*findCrossingInfo)(Front*,int*,GRID_DIRECTION,int,
                                double*,	// crx_old
				double*,	// crx_new
				double*);	// flux
private:
        // Dimension
        int dim;
        COMPONENT *top_comp;
        int *top_gmax;
	int imin,jmin,kmin;
	int imax,jmax,kmax;
	int array_size;
        double *array;          // for scatter states;
	double *top_h;
	double *top_L;
	double var_min;
        double var_max;
        double checkSolver(int *icoords,boolean print_details,double *var_in);
};

class DUAL_ELLIPTIC_SOLVER{
        Front *front;
public:
        DUAL_ELLIPTIC_SOLVER(Front &front);

        // On topological grid
	int *i_to_I;
	int **ij_to_I;
	int ***ijk_to_I;
	int ilower;
	int iupper;

	double obst_comp;
	double *soln;		/* field variable of new step */
	double *source;		/* source field */
	double *D;		/* div(D*grad)phi = source */
	void set_solver_domain(void);
	void solve(double *soln);
	double (*getStateVar)(POINTER);
	int (*findStateAtCrossing)(Front*,int*,GRID_DIRECTION,int,
                                POINTER*,HYPER_SURF**,double*);
	double checkSolver(int *icoords, boolean print_details);
private:
        // Dimension
        int dim;
        COMPONENT *top_comp,*ctop_comp;
	int *top_gmax;
	double *top_h;
	double *ctop_L;
        int *ctop_gmax;
	int cimin,cjmin,ckmin;
	int cimax,cjmax,ckmax;
	int offset[MAXD];
        double *array;          // for scatter states;
	int array_size;
	double max_soln;
	double min_soln;
	void solve1d(double *soln);
	void solve2d(double *soln);
	void solve3d(double *soln);
	void get_dual_D(int*,double*);
	double dual_average_D_2d(int dir, int nb, int**,COMPONENT**);
	double dual_average_D_3d(int dir, int nb, int***,COMPONENT***);
};

class ELLIPTIC_SOLVER{
        Front *front;
public:
        ELLIPTIC_SOLVER(Front &front);

        // On topological grid
	int *i_to_I;
	int **ij_to_I;
	int ***ijk_to_I;
	int ilower;
	int iupper;

	double porosity;
	double *soln;		/* field variable of new step */
	double *source;		/* source field */
	double *D;		/* div(D*grad)phi = source */
	void set_solver_domain(void);
	void solve(double *soln);
	void dsolve(double *soln);
	double (*getStateVar)(POINTER);
	int (*findStateAtCrossing)(Front*,int*,GRID_DIRECTION,int,
                                POINTER*,HYPER_SURF**,double*);
	double checkSolver(int *icoords,boolean print_details);
	double dcheckSolver(int *icoords,boolean print_details);
	int skip_neumann_solver;
private:
        // Dimension
        int dim;
        COMPONENT *top_comp;
	double *top_h;
	double *top_L;
        int *top_gmax;
	int imin,jmin,kmin;
	int imax,jmax,kmax;
        double *array;          // for scatter states;
	int array_size;
	double max_soln;
	double min_soln;
	void solve1d(double *soln);
	void solve2d(double *soln);
	void solve3d(double *soln);
	void dsolve2d(double *soln);
	void dsolve3d(double *soln);
};

struct _SWEEP {
	double **vel;
	double *rho;
};
typedef struct _SWEEP SWEEP;

struct _FSWEEP {
	double **vel_flux;
};
typedef struct _FSWEEP FSWEEP;

class HYPERB_SOLVER{
        Front *front;
public:
        HYPERB_SOLVER(Front &front);

	// Time step for one call
	double dt;

	int order;		/* order of Runge-Kutta */
	int size;

	COMPONENT soln_comp1;
	COMPONENT soln_comp2;
	COMPONENT obst_comp;
	double porosity;
	double **var;		/* field variable of old step */
	double **soln;		/* field variable of new step */
	double **source;	/* source field */
	double *rho;
	double rho1;
	double rho2;
	double var_obst;	/* default solution in obst_comp */
	double max_speed;
	int (*findStateAtCrossing)(Front*,int*,GRID_DIRECTION,int,
                                POINTER*,HYPER_SURF**,double*);
	void (*numericalFlux)(SWEEP*,FSWEEP*,double,int,int,int);
	void solveRungeKutta();
	double (*getStateVel[3])(POINTER);
private:
        // Dimension
        int dim;
        COMPONENT *top_comp;
	double *top_L,*top_h;
        int *top_gmax;
	int imin,jmin,kmin;
	int imax,jmax,kmax;
	int nrad;
        double *array;          // for scatter states;
	SWEEP *st_field,st_tmp;
        FSWEEP *st_flux;
        double **a,*b;
	void setSolverDomain(void);
	void allocMeshVst(SWEEP*);
	void allocMeshFlux(FSWEEP*);
	void allocDirectionalVstFlux(SWEEP*,FSWEEP*);
	void resetFlux(FSWEEP*);
	void copyToMeshVst(SWEEP*);
	void copyMeshVst(SWEEP,SWEEP*);
	void computeMeshFlux(SWEEP,FSWEEP*);
	void addMeshFluxToVst(SWEEP*,FSWEEP,double);
	void copyFromMeshVst(SWEEP);
	void addFluxInDirection(int,SWEEP*,FSWEEP*);
	void addFluxInDirection1d(int,SWEEP*,FSWEEP*);
	void addFluxInDirection2d(int,SWEEP*,FSWEEP*);
	void addFluxInDirection3d(int,SWEEP*,FSWEEP*);
	void appendGhostBuffer(SWEEP*,SWEEP*,int,int*,int,int);
	void setNeumannStates(SWEEP*,SWEEP*,HYPER_SURF*,POINTER,int*,int,
				int,int,int,int);
	void setDirichletStates(SWEEP*,SWEEP*,HYPER_SURF*,POINTER,int*,int,
				int,int,int,int);
	void setElasticStates(SWEEP*,SWEEP*,HYPER_SURF*,POINTER,int*,int,
				int,int,int,int);
	void addSourceTerm(SWEEP*,FSWEEP*);
};


extern	void upwind_flux(SWEEP*,FSWEEP*,double,int,int,int);
extern	void weno5_flux(SWEEP*,FSWEEP*,double,int,int,int);

/* For CIM solver */

#define CIM_TOL 1e-8

typedef struct {
	int N[MAXD], NA;
	double *u;
} Unknown;

typedef struct {
	int	N;		/* row size of matrix */
	int 	K;		/* nonzero elements of A */
	int	*i, *j;		/* (i,j) <-> [k], 1<= k <= K	*/
	double  *a;		/* a[k] -> a(i,j) */
} MTX;

typedef struct {
	int	N;		/* grid size */
	MTX	P;		/* prolongation matrix */
	MTX	A;		/* sparse matrix */
	double	*x;		/* vector structure */
	double	*b;		/* vector structure */
} GRID;

#define tr(s) ((s) == 1 ? 0 : 1)

typedef struct {
	int s[2*MAXD];		// indicator of intersections
	double a[2*MAXD];	// relative locations of intersections
	double n[2*MAXD][MAXD];	// normal vectors at the intersections
	int c[MAXD][MAXD];	// use only in cim2, cross derivative indicator
} CIM_STRUCT;

typedef struct {
	int    indx[MAXD][5*MAXD-1];
	double coef[MAXD][5*MAXD-1];
	double J[MAXD][6*MAXD];
} CIM_COEF;

typedef struct {
    	int indx[MAXD][6*MAXD+4];
    	double coef[MAXD][6*MAXD+4];
    	double v[MAXD];
} CIM_GRAD;

typedef struct {
    	int ID, ID_S[2], ID_R, R_dir, R[6*MAXD];
    	double C[6*MAXD];
} GHOST_VALUE;

typedef struct {
	int s[MAXD];
	double a[MAXD];
	double n[MAXD][MAXD];
	int c[MAXD][MAXD];
} CIM2_STRUCT;

typedef struct {
	int    indx[MAXD][6*MAXD-1];
	double coef[MAXD][6*MAXD-1];
	double J[MAXD][3*MAXD];
} CIM2_COEF;

class CIM_ELLIPTIC_SOLVER{
        Front *front;
public:
        CIM_ELLIPTIC_SOLVER(Front &front);

        int *i_to_I,*I_to_i;            // Index mapping for 1D
        int **ij_to_I,**I_to_ij;        // Index mapping for 2D
        int ***ijk_to_I,**I_to_ijk;     // Index mapping for 3D
	int ilower;
	int iupper;
	int size;
	boolean solve_front_state;

	POINTER jparams;	/* Params for jump functions */
	double *soln;		/* field variable of new step */
	double *source;		/* source field */
	double diff_coeff[2];	/* div(diff_coeff*grad)phi = source */
	double *kappa;
	int w_type;
	COMPONENT pos_comp,neg_comp;
	void set_solver_domain(void);
	void solve(double *soln);
	void (*assignStateVar)(double,POINTER);
	double (*getStateVar)(POINTER);
	int (*findStateAtCrossing)(Front*,int*,GRID_DIRECTION,int,
                                POINTER*,HYPER_SURF**,double*);
	double (*solutionJump)(POINTER jparams,int D,double *coords);
	double (*gradJumpDotN)(POINTER jparams,int D,double *N,double *P);
	double (*gradJumpDotT)(POINTER jparams,int D,int i,double *N,double *P);
	double (*exactSolution)(POINTER jparams,int D,double *P);
private:
        // On topological grid
        int dim;
        RECT_GRID *top_grid;
        COMPONENT *top_comp;
	double *top_h;
	double *top_L;
        int *top_gmax;
	int imin,jmin,kmin;
	int imax,jmax,kmax;

	// CIM specific variables
	int  *NB[2*MAXD],Ord[5];
	char *D, *S[2*MAXD], *Order;
	MTX A;
	double *b, *V;
	Unknown u;
	int NEWCIM;

	// CIM specific functions
	int HCIM_Matrix_Generation();    
	int HCIM_Matrix_Generation_k(int k,int *n);    
	int Generate_Domain_and_NB();
	int Check_Type();
	int Check_Type_k(int k);
	void Set_Intersection();
	int Search_Ghost_Value(int k);
	void cimSolveFrontState();
	void cimIntfcPointState(double*,int,double*,double*);

	void solve1d(double *soln);
	void solve2d(double *soln);
	void solve3d(double *soln);
};

class CIM_PARAB_SOLVER{
        Front *front;
public:
        CIM_PARAB_SOLVER(Front &front);

        int *i_to_I,*I_to_i;            // Index mapping for 1D
        int **ij_to_I,**I_to_ij;        // Index mapping for 2D
        int ***ijk_to_I,**I_to_ijk;     // Index mapping for 3D
	int ilower;
	int iupper;
	int size;
	boolean solve_front_state;

	POINTER jparams;	/* Params for jump functions */
	double m_dt;		/* Time step for parabolic equation */
	double *soln;		/* field variable of new step */
	double *source;		/* source field */
	double diff_coeff[2];	/* div(diff_coeff*grad)phi = source */
	double *kappa;
	int w_type;
	COMPONENT pos_comp,neg_comp;
	void set_solver_domain(void);
	void solve(double *soln);
	void (*assignStateVar)(double,POINTER);
	double (*getStateVar)(POINTER);
	int (*findStateAtCrossing)(Front*,int*,GRID_DIRECTION,int,
                                POINTER*,HYPER_SURF**,double*);
	double (*solutionJump)(POINTER jparams,int D,double *coords);
	double (*gradJumpDotN)(POINTER jparams,int D,double *N,double *P);
	double (*gradJumpDotT)(POINTER jparams,int D,int i,double *N,double *P);
	double (*exactSolution)(POINTER jparams,int D,double *P);
private:
        // On topological grid
        int dim;
        RECT_GRID *top_grid;
        COMPONENT *top_comp;
	double *top_h;
	double *top_L;
        int *top_gmax;
	int imin,jmin,kmin;
	int imax,jmax,kmax;

	// CIM specific variables
	int  *NB[2*MAXD],Ord[5];
	char *D, *S[2*MAXD], *Order;
	MTX A;
	double *b, *V;
	Unknown u;
	int NEWCIM;

	// CIM specific functions
	int HCIM_Matrix_Generation();    
	int HCIM_Matrix_Generation_k(int k,int *n);    
	int Generate_Domain_and_NB();
	int Check_Type();
	int Check_Type_k(int k);
	void Set_Intersection();
	int Search_Ghost_Value(int k);
	void cimSolveFrontState();
	void cimIntfcPointState(double*,int,double*,double*);

	void solve1d(double *soln);
	void solve2d(double *soln);
	void solve3d(double *soln);
};

extern void viewTopVariable(Front*,double*,boolean,double,double,char*,char*);
extern double   compBdryFlux(Front*,double*,double,int,double,
		double(*)(POINTER));
#endif

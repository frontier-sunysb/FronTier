#include <vector>
#include <petscksp.h>
#include <assert.h>

#include "FronTier.h"


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


// Old interface.h

extern double exact_eps(int D);
extern double exact_solution(int D, double *P);
extern double exact_source(int D, double *P);
extern double exact_jump_u(int D, double *P);
extern double exact_jump_eps_gradu_dot_n(int D,
		double *N,double *P,int dim);
extern double exact_jump_gradu_dot_t(int D,int i,double *N,
		double *P,int dim);
extern int exact_Gradient(int D,double *G,double *P,int dim);
extern int tangential_direction(int i,double *T,double *N,int dim);

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

        int *lbuf,*ubuf;
        int *i_to_I,*I_to_i;            // Index mapping for 1D
        int **ij_to_I,**I_to_ij;        // Index mapping for 2D
        int ***ijk_to_I,**I_to_ijk;     // Index mapping for 3D

        // Sweeping limites
        int imin,jmin,kmin;
        int imax,jmax,kmax;

	/*TMP*/
	int  *NB[2*MAXD],Ord[5];
	char *D, *S[2*MAXD], *Order;
	MTX A;
	double *uex, *b, *V;
	Unknown u;
	int NEWCIM;
	int size;

	// constructor
        ~C_CARTESIAN();

	// for parallel partition
        int             NLblocks,ilower,iupper;
        int             *n_dist;
	// mesh: full cells mesh
	void solve();
        void setComponent(void);        // init components
	void setDomain();
	void setIndexMap(void);
	void scatMeshArray();
        void setGlobalIndex();

	int HCIM_Matrix_Generation();    
	int Generate_Domain_and_NB();
	int Check_Type();
	void Set_Intersection();
	int Search_Ghost_Value(int k);
	
	void cim_solve(void);
};

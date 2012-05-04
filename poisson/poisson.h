/*      	poison.h
 *
 *	PetscInitialize(&argc,&args,(char *)0,help);
 *
 *  To know how to compile the Hypre class, please read the makefile.
*/ 

#include <FronTier.h>
#include "petscksp.h"

#define		MAX_NUM_VERTEX_IN_CELL		20
#define 	Cross_prod(x0,y0,x1,y1,xm,ym)                            \
        	((x0 - xm)*(y1 - ym) - (x1 - xm)*(y0 - ym))

struct _P_CELL2D {
        double    stencil[3][3];         
	boolean     set[3][3];
        double    r_side;
};
typedef struct _P_CELL2D P_CELL2D;

struct _LAPLACE {
        double 	D;
        double 	*solute;
	void 	*sparams;
	void  	(*flux_func)(POINTER,double*,double*);
	double  (*sfunc)(POINTER,double*);
	double  (*solution)(POINTER,double*);
};
typedef struct _LAPLACE LAPLACE;

typedef double STATE;

extern  void embed_bdry_poison(POINTER,Front*);
extern	void s_hyp_solution(double*,COMPONENT,HYPER_SURF*,SIDE,Front*,
			POINTER,POINTER,POINTER);

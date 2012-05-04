/***********************************************************************************************************
 *			SN_ELLIP.h
 * this class is used to solve an elliptic problem:
 *	
 * this is based on a cell center finite difference method; 
 * 		
 * see also FEM2D.h
 ***********************************************************************************************************/

#ifndef SN_ELLIP_H
#define SN_ELLIP_H
#include <stdio.h>
#include <math.h>
/*#undef MPI_Comm  */
/*#undef MPI_Group */
#include <mpi.h>

#include <gutil.h>
#include <pgutil.h>
#include <Hypre.h>
#include <Hypre_GMRES.h>
#include <VectorMatrix.h>
#ifndef PI
#define PI 3.14159265358979323846264338327950288419716939937510
#endif	




/* see also QUADTREE */
enum INDEX_LOCALGLOBAL {INDEX_LOCAL = -1, INDEX_GLOBAL = -2};				/* all enum type should be changed to capital. */
enum RECEIVER_SENDER {RECEIVER=0, SENDER};

#define Beta_Minus	1
#define Beta_Plus	1
#define Radius		0.5
#define Elevation	2


class SN_ELLIP {

	double m_X1[2], m_X2[2];
	double m_x1[2], m_x2[2];		/* computational domain */
	double m_dx, m_dy;
	int m_nx, m_ny;				/* grid number in the x/y direction */
	int m_nfluxbuffer;			/* number of buffers for flux in each direction, default to be 2 */
	
	
	/* member for parallel partition */
	RECEIVER_SENDER m_whoami;
	MPI_Comm  m_comm;	
	int m_comm_col, m_comm_row, m_comm_size, m_comm_rank;
	/* the ranks for this node and its neighbors. */
	int m_NRank, m_WRank, m_SRank, m_ERank, m_MRank;
	
	int m_pi, m_pj;		/* the node lies on the pith col, pjth row; */
	
	
	int m_starting_index, m_ending_index;
	int *m_BoundaryVertexIndex_North;		/* north boundary index */
	int *m_BoundaryVertexIndex_West;		/* west boundary index */
	int *m_BoundaryVertexIndex_South;		/* south boundary index */
	int *m_BoundaryVertexIndex_East;		/* east boundary index */
	int *m_GlobalIndex;				/* used to show whether the vertex is a local vertex or global vertex */
	
	double *m_BoundaryPotential_North;		/* m_comm_col*m_nx */
	double *m_BoundaryPotential_East;		/* m_comm_row*m_ny */
	
	double m_G;
	double **m_Density;			
	
	
	double *m_unknown;
	double *m_unknown_North, *m_unknown_West, *m_unknown_South, *m_unknown_East;
public:
	Hypre		m_solver;	
	/*Hypre_GMRES	m_solver;	 */
	
	
	RECEIVER_SENDER whoami(void);		/* return which nodes are first to get information from their neighbors */
	int  getIndex(int i, int j);
	void setGlobalIndex(void);		/* GB_GlobalVertexIndexing */
	void setGlobalIndex_Boundary(void);	/* GB_GlobalVertexIndexing_UpdateBoundary */
	void getGlobalSolution(void);
	void getGlobalSolution_Boundary(void);	/* */
	
	double computeIntegralTheta(double r, double z, double r1, double z1);	
	void computeBoundaryPotential(void);
	double *m_buffer;
	void setMatrix(void);
	
	
	void getFlux_Boundary(double **fx, double **fy);
	
	
public:	
	SN_ELLIP();
	~SN_ELLIP();
	
	
	void setParameters(double X1[2], double X2[2], double x1[2], double x2[2], int nx, int ny, double G, double **density);
	void setParaParameters(MPI_Comm comm, int row, int col, int nrank, int wrank, int srank, int erank);
	
	
	double func_f(int i, int j);
	
	void solve(void);
	void getPotential(double **phi);
	
	void setFluxBuffer(int n);
	void getFlux(double **fx, double **fy);		/* fx[m_nx+2m_nfluxbuffer][m_ny+2m_nfluxbuffer], fy[m_nx+2m_nfluxbuffer][m_ny+2m_nfluxbuffer] */
	
	void print(char *filename);	
	void printFlux(double **fx, double **fy);	
};

#endif



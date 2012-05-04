/* 
 * File:   ebm2d.h
 * Author: shuqiang (robert) wang
 *
 * Created on September 28, 2008, 9:28 PM
 * 
 * This class is used for the Embedded Boundary Method for solving either 
 * elliptic boundary value or elliptic interface problem.
 *
 * The EBM2D is in fact an finite colume method with the unknowns defined at the
 * cell center. For the elliptic boundary value problem, only one unknown in
 * each cell is needed. However, for the elliptic interface problem, 1 or 2 or
 * 3 unknowns are needed at the cell center. 
 *
 * There are many different way of storing/referening the unknowns at the cell 
 * center. Therefore, it is necessary to hide the implementation so that other
 * modules are not dependent on the implementation. Even if the implementation
 * is changed, there is no need to change other modules.
 *
 * The equation to be solved is the following
 *      
 *      Div * beta Grad phi = rou
 *
 * References:
 * 1) Hans Johansen and Phillip Colella, A Cartesian Grid Embedded Boundary 
 *    Method for Poisson's Equation on Irregular Domains, JCP 147, 60-85 (1998).
 *
 * See also lpetsc.c.
 *
 * To Do:
 * 1) get information from FronTier using the following functions:
 *         component(double*,INTERFACE*);
 *         intersections();
 *              i_intersections2d();
 *         insert_grid_crossing2d();
 *         crossing_in_direction();
 *         nearest_interface_point();
 * 
 */

#ifndef _EBM2D_H
#define	_EBM2D_H

#include <solver.h>
#include <utility>

#include <FronTier.h>
#undef vector 
#undef assign 
#include <vector>
#include <ebm.h>


/*****************************************************************************
 *                 EBM2D_LAPLACE
 * The EBM2D_LAPACE is used to represent the elliptic interface problem:
 *
 *      Div * beta Grad phi = rou
 * see also EBM3D_LAPLACE.
 *****************************************************************************/
class EBM2D_LAPLACE {
public:
    EBM2D_LAPLACE();
    virtual double getBeta(double coords[2], int comp);             // beta
    virtual double getExactSolution(double coords[2], int comp);    // phi
    virtual double getFlux(double coords[2], int comp, double normal[2]);
    virtual double getRightHandSide(double coords[2], int comp);    // rou    
};

/******************************************************************************
 *                  EBM2D_POINT
 * It is simply a wrapper for 2D coordidate.
 * see also EBM3D_POINT.
 *****************************************************************************/
class EBM2D_POINT {
public:
     double m_coords[2];
     EBM2D_POINT();
     EBM2D_POINT(double coords[2]);
     void getCoords(double coords[2]);
     void setCoords(double coords[2]);
};

/*****************************************************************************
 *                  EBM2D_EDGE 
 * The component always point from a smaller component to a bigger component!
 * This structure is used to calculate the length/normal/center of the part 
 * of the intfc inside a grid block.
 * 
 * To Do: merge this with EBM2D_PARTIALEDGE.
 * 
 * see also EBM3D_TRIANGLE.
 *****************************************************************************/
class EBM2D_EDGE {
     int m_comp[2];                // negative and then positive component.
     int m_index[2];
public:
     double m_vertex[2][2];        // m_vertex[0] is the coords of the 1st point.
     void getComp(int comp[2]);
     void getIndex(int index[2]);
     
     int setEdge(int index0, int comp0, int index1, int comp1, double v0[2], double v1[2]);
     int setEdge(double v0[2], double v1[2]);       // used for debugging only.

     double getLength(void);
     void getNormal(double normal[2]);
     void getLengthNormalCenter(double &length, double normal[2], double center[2]);
};

/******************************************************************************
 *          EBM2D_PARTIALEDGE 
 * see also EBM3D_PARTIALFACE.
 ******************************************************************************/
class EBM2D_PARTIALEDGE {            
public:
    int m_comp;
    int m_index;    
    double m_length;   
    double m_center[2];
    
    // first function to call in a sequence of steps to get the partial face center.        
    void begin(int c, int index);
    //void begin(int c, int index, double area, double center[2]);
    void addEdge(double v0[2], double v1[2]);
    void addEdge(double v0[2], double v1[2], double v2[2]);
    void end(void);      // the last function 

    // for debugging
    std::vector<EBM2D_EDGE> m_edges;
    void debug_begin(void);
    void debug_addEdge(double v0[2], double v1[2]);
private:
    double getEdgeLength(double v0[2], double v1[2]);
    void getEdgeCenter(double v0[2], double v1[2], double center[2]);
};

/******************************************************************************
 *      EBM2D_CELL
 * This class is the struct for a single cell.
 *
 * For a given cell, the order of the corner is the following:
 *      |3 2|
 *      |0 1| 
 * For a given cell, the order of the edge is the following:
 *       2
 *      3 1
 *       0
 * For getCellConectionMatrix, the connection matrix means that  
 * if matrix[i][j]!=0 then there is an edge between vertex i and vertex j.
 * The ordering the vertices are the following:
 *      362
 *      7 5
 *      041
 * where 1,3,5,6 means the crxing point of the edge.
 ******************************************************************************/
class EBM2D_CELL {    
    // indices to the cell center unknowns, stored at the four corners of 
    // the cell 
    int m_nCellUnknowns;
    int m_index[4];
    CELL_TYPE m_cellType;
public:
    EBM2D_CELL();   
    void init(void);
    
    CELL_TYPE getCellType(void);
    void setCellType(int *comp=m_comp);
    int getCellUnknownIndex(int index[4]);
    void setCellUnknownIndex(int index[4]);
    int getCellUnknownIndex(int corner);      // corner: 0-3.        
    int getCellUnknownNumber(void);
    int getCellIntfcIndex(void);              // the first intfc index
    
    // the following functions are used for the generalized marching cubes method
    // to get the number of unknowns, the volume fractions, the area of the surfaces
    // and their normal.
    void setComp(int comp[4]);
    void setCompCoordsCrx(int comp[4], double coords[4][2], double crx[5][2]);
    
    int setCellUnknownIndex(int start_index);      // setComp must be called before this.
    void setCellUnknownIndex(int corner, int index);
    void setCellUnknownIndex(int start_index, int comp[5], int numCellUnknowns);      

    void getUnknownToCompMap(std::map<int, int> &unknownToComp);
    void getAreaAndIntfc(std::map<int, double> &areas,
			 std::map<int, std::vector<EBM2D_EDGE> > &intfc);
    void getTriAreaAndIntfc(int verticies[3], int edges[3],
			    std::map<int, double> &volumes,
			    std::map<int, std::vector<EBM2D_EDGE> > &intfc);

    void getAveragedIntfcLengthNormalCenter(std::map<int, std::vector<EBM2D_EDGE> > &intfc,
					    std::map<int, double> &intfcLength,
					    std::map<int, EBM2D_POINT> &intfcNormal,
					    std::map<int, EBM2D_POINT> &intfcCenter);
   
    void getPartialEdges(std::map<int, std::vector<EBM2D_PARTIALEDGE> > &partialedges);
    void getPartialEdges(int corner[2], int e, std::vector<EBM2D_PARTIALEDGE> &partialedges);


    // functions for parallel 
    void resetCellUnknownIndex(int ilower);    

    // the following functions are used for debugging or output
    void test(void);
    void test_checkAreas(std::map<int, double> &areas);
    void test_getPartialEdges(std::map<int, std::vector<EBM2D_PARTIALEDGE> > &partialedges);
    void saveIntfc_tecplot(const char *filename, std::map<int, std::vector<EBM2D_EDGE> > &intfc);

public:
    void assign(int corner[2], int v0, int v1);
    void assign(int corner[3], int v0, int v1, int v2);
    void assign(int corner[4], int v0, int v1, int v2, int v3);
  
    static void getCenter(double v0[2], double v1[2], double center[2]);  
    static double getDistance(double v0[2], double v1[2]);
    static void interpolate(double v0[2], double v1[2], double t, double p[2]);
    static double getEdgeLength(double v0[2], double v1[2]);
    static double getTriArea(double v0[2], double v1[2], double v2[2]);
    
    static int *m_comp;
    static double (*m_coords)[2];
    static double (*m_crx)[2];
    static int m_cellEdges[4][3];
    static int m_cellEdgePairs[5][2];
						 
};


/*
 *      EBM2D_CARTESIAN
 * This class is the struct for the whole method. It uses the EBM2D_CELL as the 
 * cell of the whole domain.
 * 
 * This implementation uses an underlying cartesian grid to partition the 
 * computational domain.   
 * 
 * 
 */
class EBM2D_CARTESIAN {    
    Front *m_pFront;    
    SOLVER *m_pSolver;
    EBM2D_LAPLACE *m_pLaplace;
    
    int m_gmax[2];            // cartesian grid.
    int m_myid;               
    int m_lbuf[2], m_ubuf[2]; // used for parallel
    int m_ilower, m_iupper;   // matrix range

    EBM2D_CELL **m_cells;
    int m_nComps;             // size of m_comps[].
    int *m_comps;             // components for the vertices of the grid cells.    
    int m_nLocalUnknowns;     // number of local unknowns.
    double **m_unknowns[2];    // unknowns for the bigger/smaller component.
    int **m_compList[2];

    int m_i;           // used for member function default parameters.
    int m_j;           // for the current cell index.

    double m_dx;
    double m_dy;
    double m_mindd;           // min, max od m_dx, m_dy.
    double m_maxdd;                 
    double m_cellArea;

    double m_mindd4;          // for getCrx().
    int m_MAX_ITER;           // for getCrx().
    
    // grids
    void getCellCenter(double center[2], int i, int j);
    void getCellEdgeCenter(int face, double center[2], int i, int j);
    double getCellEdgeCenter(int dir, int i);
    double getCellEdgeLength(int face, int i, int j);

    int getCellCornerIndex(int corner, int i, int j);  // corner: 0,1,2,3.
    void getCellCornerCoords(int corner, double coords[2], int i, int j);
    void getCellCornerCoords(double coords[4][2], int i, int j);

    void getCurrentCell(int &i, int &j);
    void setCurrentCell(int i, int j);
    void locateCell(double coords[2], int &i, int &j);
    
    // Components
    int getCellCornerComp(int index);  
    int getCellCornerComp(int corner, int i, int j);
    int getCellEdgeComp(int face, int comp[2], int i, int j);
    
    // unknown indices
    int getCellUnknownIndex(int neighbor, int comp, int i, int j);          
    int getCellUnknownIndex(int comp, int i, int j);    // use the current cell index (m_i, m_j)
    int getCellUnknownIndexByTrial(int comp, int i, int j);
    int getCellNeighborUnknownIndex(int face, int comp, int i, int j);
    int getNeighborByPlane(EBM_PLANE plane, int x, int comp, int i, int j);
    int getNeighborAcrossPlane(EBM_PLANE plane, int x, int comp, int i, int j);
    
    void getCellCompList(int comp[2], int i, int j);
    void setCellCompList(int comp[2], int i, int j);
    double getCellUnknown(int comp, int i, int j);
    double getCellUnknownAcrossPlane(EBM_PLANE plane, int x, int comp, int i, int j);
    void getCellUnknown(double unknown[2], int i, int j);
    void setCellUnknown(double unknown[2], int i, int j);
    
    void setCellUnknown(double *x);
   
    // crossing & shapes
    EBM2D_CELL *getCell(int i, int j);  
    CELL_TYPE getCellType(int i, int j);       
    int getCellUnknownNumber(int comp[5], int i, int j);     // 0,1,2,3    

    boolean getCellEdgeCrx(int face, double crx[2], int i, int j);  // edge: 0,1,2,3.        
    void getCellComp(int &c0, int &c1, int i, int j);
    void getCellComp(int comp[4], int i, int j);
    
    void getCellCompAndCrx(int comp[5], double coords[8][2], int i, int j);
    void getCellCompCoordsCrx(int comp[4], double coords[4][2], double crx[5][2], int i, int j);
 
    // to be implemented using FronTier functions
    int getComp(double *coords);
    void getCrx(double p0[2], double p1[2], double crx[2]);

    void setEquivalentComponents(void);        // modify components.

    // matrix setup    
    void setMatrixForInternalCell(int i, int j);   
    void setMatrixForPartialCell_2Unknowns(int i, int j);    
    void setMatrixForBoundaryCell(int i, int j);

    void setMatrixForPartialCell_3Unknowns(int i, int j);             // not implemented.

    void setMatrixForPartialCell_2Unknowns_2ndOrder(int i, int j);    
    void setMatrixForPartialCell_2Unknowns_Case0111(int comp[5], int index[5], double coords[8][2], int i, int j);
    void setMatrixForPartialCell_2Unknowns_Case0100(int comp[5], int index[5], double coords[8][2], int i, int j);
    void setMatrixForPartialCell_2Unknowns_Case0010(int comp[5], int index[5], double coords[8][2], int i, int j);
    void setMatrixForPartialCell_2Unknowns_Case0001(int comp[5], int index[5], double coords[8][2], int i, int j);
    void setMatrixForPartialCell_2Unknowns_Case0110(int comp[5], int index[5], double coords[8][2], int i, int j);
    void setMatrixForPartialCell_2Unknowns_Case0011(int comp[5], int index[5], double coords[8][2], int i, int j);
    void setMatrixForPartialCell_2Unknowns_byRotation(int rotation, int comp[5], int index[5], double coords[8][2], int i, int j);
    
    void getRotatedCellsIndex(int rotation, int IIndex[3][3], int JIndex[3][3], int i, int j);
    void getRotatedCoordsIndex(int rotation, int CIndex[8]);    
    void getRotatedDirection(int rotation, int DIndex[2]);   

    // matrix setup for intfc/boundary.
    void setMatrixForInterface(int intfcIndex[3], double intfcCenter[2], double intfcNormal[2], 
			       int intfcComp[2], double intfcLength, double areaFraction[2], int i, int j);
    void setMatrixForBoundary_Dirichlet(int intfcIndex[2], double intfcCenter[2], double intfcNormal[2], 
			      int intfcComp[2], double intfcLength, double areaFraction[2], int i, int j);
    void setMatrixForBoundary_Neumann(int intfcIndex[2], double intfcCenter[2], double intfcNormal[2], 
			      int intfcComp[2], double intfcLength, double areaFraction[2], int i, int j);
    
    // flux stencil for the partial edges
    void getFluxStencilForPartialEdge(int face, double center[2], double length, int comp,
				    std::vector< std::pair<int,double> > & stencil,
				    int i, int j);
    void getFluxStencilForPartialEdge2(int face, double center[2], double length, int comp,
				     std::vector< std::pair<int,double> > & stencil,
				     int i, int j);

    double getAreaFraction(int index, std::map<int, double> &areas);


    // flux stencil for the boundary/intfc
    void getSolutionInterpStencil(double edgeCenter[2], double normal[2], int comp,
				  std::vector<int> &indices, 
				  std::vector<double> &coefficients, 
				  std::vector<double> &distances, 
				  int i, int j);    
    
    double getDerivativeAtCellCenter(EBM_PLANE plane, int comp, int i, int j);    // get derivative perpendicular to plane
    //double getDerivative(double coords[2], double normal[2], int comp, int i, int j);

    void get3PointsSolutionInterpCoeffs(double coords[3], double x, double coeffs[3]);
    void get3PointsDerivativeInterpCoeffs(double coords[3], double x, double coeffs[3]);
    
    // functions for parallel
    void pp_setup(void);
                          
    void pp_resetCellUnknownIndex(void);                                   
    void pp_resetBufferCellUnknownIndex(void);              
    void pp_packCellUnknownIndex(int dir, int size, int **data);
    void pp_unpackCellUnknownIndex(int dir, int size, int **data);
   
    void pp_resetBufferCellCompList(void);
    void pp_packCellCompList(int dir, int size, int **data);
    void pp_unpackCellCompList(int dir, int size, int **data);

    void pp_resetBufferCellUnknown(void);
    void pp_packCellUnknown(int dir, int size, double **data);
    void pp_unpackCellUnknown(int dir, int size, double **data);
            
    // Setup the data structure.        
    void setup(void);       
    
    // The following functions are used for the EBM2D methods.     
    void setMatrix(void);

public:
    EBM2D_CARTESIAN(Front *front, SOLVER *solver, EBM2D_LAPLACE *laplace);
    ~EBM2D_CARTESIAN();
    void solve(void); 

    // output functions
    void saveInterface_Tecplot(const char *filename);
    void saveReconstructedIntfc_Tecplot(const char *filename, int i, int j);
    void saveInterfaceEdges_Tecplot(const char *filename);
    void savePartialCells_Tecplot(const char *filename, 
				  std::map<int, std::vector<EBM2D_PARTIALEDGE> > &mapPartialEdges,
				  std::map<int, std::vector<EBM2D_EDGE> > &mapIntfc,
				  int i, int j);
    
    void saveStates_Tecplot(void);
    void saveStates_VTK(void);
    void saveDerivatives_Tecplot(void);
    void saveComponent_Tecplot(void);
    
    // debugging functions
    void debug(int i, int j);    
    void debug(std::vector< std::pair<int,double> > &stencil);
    void debug_printNeighborComp(int i, int j);
    void debug_printCellComp(int comp[4], int i, int j);
    void debug_printCellUnknownIndex(int i, int j);
    void debug_saveUnknownIndex(void);
    void debug_saveCompList(const char*filename);
    

    double m_minCellArea;
    double m_minPartialEdgeLength;
    double m_minIntfcLength;
	 
    void debug_getMinStatistics(int i, int j);

};



#endif	/* _EBM2D_H */


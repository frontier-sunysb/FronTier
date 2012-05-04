/*
 * File:   ebm3d.h
 * Author: shuqiang (robert) wang
 *
 * Created on Oct 15, 2008.
 *
 * This class is used for the Embedded Boundary Method for solving either
 * elliptic boundary value or elliptic interface problem.
 *
 * The EBM is in fact an finite column method with the unknowns defined at the
 * cell center. For the elliptic boundary value problem, only one unknown in
 * each cell is needed. However, for the elliptic interface problem, 1 or 2 or
 * 3 unknowns are needed at the cell center.
 *
 * There are many different way of storing/referencing the unknowns at the cell
 * center. Therefore, it is necessary to hide the implementation so that other
 * modules are not dependent on the implementation. Even if the implementation
 * is changed, there is no need to change other modules.
 *
 * The equation to be solved is the following
 *
 *      Div * beta Grad phi = rou
 *
 * This is a 3D version implementation of the EBM metthod.
 * References:
 * 1) Hans Johansen and Phillip Colella, A Cartesian Grid Embedded Boundary
 *    Method for Poisson's Equation on Irregular Domains, JCP 147, 60-85 (1998).
 *
 * See also lpetsc.c.
 * See also EBM.h/EBM.c
 */

#ifndef _EBM3D_H
#define	_EBM3D_H

#include <solver.h>
#include <utility>

#include <FronTier.h>
#undef vector
#undef assign
#include <vector>
#include <map>
#include <ebm.h>


/*
 * The EBM3D_LAPACE is used to represent the elliptic interface problem:
 *
 *      Div * beta Grad phi = rou
 *
 */
class EBM3D_LAPLACE {
public:
    EBM3D_LAPLACE();
    virtual double getBeta(double coords[3], int comp);             // beta
    virtual double getExactSolution(double coords[3], int comp);    // phi
    virtual double getFlux(double coords[3], int comp, double normal[3]);
    virtual double getRightHandSide(double coords[3], int comp);    // rou
};

/*****************************************************************************
 *          EBM3D_INDEX
 *****************************************************************************/
class EBM3D_INDEX {
public: 
     int m_index[3];
     void getIndex(int index[3]);
     void setIndex(int index[3]);
};

/******************************************************************************
 *          EBM3D_POINT
 ******************************************************************************/
class EBM3D_POINT {
public:
    double m_coords[3];
    EBM3D_POINT();
    EBM3D_POINT(double coords[3]);    
    void getCoords(double coords[3]);
    void setCoords(double coords[3]);        
};

/******************************************************************************
 *          EBM3D_TRIANGLE
 * always point from a smaller component to a bigger component!
 ******************************************************************************/
class EBM3D_TRIANGLE {
    int m_comp[2];      // negative and then positive component
    int m_index[2];          
public:   
    double m_vertex[3][3];
    void getComp(int comp[2]);
    void getIndex(int index[2]);
    
    int setTriangle(int index0, int comp0, int index1, int comp1, double v0[3], double v1[3], double v2[3]);    
    int setTriangle(double v0[3], double v1[3], double v2[3]);           // used for debugging only.

    double getArea(void);
    void getNormal(double normal[3]);  
    void getAreaNormalCenter(double &area, double normal[3], double center[3]);     // ??? center ???
};

/******************************************************************************
 *          EBM3D_INTFC
 * used to calculate the area, normal, center of the intfc patch inside a grid
 * block.
 ******************************************************************************/
//class EBM3D_INTFC {
//    int m_comp[2];
//    int m_index[2];
//    
//    double m_area;
//    double m_normal[3];
//    double m_center[3];
//    
//public:   
//    void begin(EBM3D_TRIANGLE &triangle);
//    void addTri(EBM3D_TRIANGLE &triangle);    
//    void end(void);
//    
//    void getComp(int comp[2]);
//    void getIndex(int index[2]);
//    double getArea(void);
//    void getNormal(double normal[3]);
//    void getAreaNormalCenter(double &area, double normal[3], double center[3]);    
//};

/******************************************************************************
 *          EBM3D_PARTIALFACE 
 ******************************************************************************/
class EBM3D_PARTIALFACE {            
public:
    int m_comp;
    int m_index;    
    double m_area;   
    double m_center[3];
    
    // first function to call in a sequence of steps to get the partial face center.        
    void begin(int c, int index);
    //void begin(int c, int index, double area, double center[3]);
    void addTri(double v0[3], double v1[3], double v2[3]);
    void addQuad(double v0[3], double v1[3], double v2[3], double v3[3]);
    void end(void);      // the last function 

    // for debugging
    std::vector<EBM3D_TRIANGLE> m_triangles;
    void debug_begin(void);
    void debug_addTri(double v0[3], double v1[3], double v2[3]);
private:
    double getTriArea(double v0[3], double v1[3], double v2[3]);
    void getTriCenter(double v0[3], double v1[3], double v2[3], double center[3]);
};


/******************************************************************************
 *          EBM3D_CELL
 ******************************************************************************/
class EBM3D_CELL {
    int m_nCellUnknowns;
    int m_index[8];     
    CELL_TYPE  m_cellType;       
public:
    EBM3D_CELL();  
    void init(void);
    
    CELL_TYPE getCellType(void); 
    void setCellType(int *comp=m_comp);       // comp[8]
    int getCellUnknownIndex(int index[8]);
    void setCellUnknownIndex(int index[8]);
    int getCellUnknownIndex(int corner);      // corner: 0-7.        
    int getCellUnknownNumber(void);   
    int getCellIntfcIndex(void);                // the first intfc index    
    
    // the following functions are used for the generalized marching cubes 
    // method to get the number of unkowns, the volume fractions, the area of 
    // the surfaces and their normal.
    void setComp(int comp[8]);
    void setCompCoordsCrx(int comp[8], double coords[8][3], double crx[19][3]);         
    void getCellFaceCompCoords(int face, int comp[4], double coords[4][3]);
    int getCellFaceNearestCorner(int face, int comp, double p[3]);

    int setCellUnknownIndex(int start_index);       // setComp must be called before this function.
    void setCellUnknownIndex(int corner, int index);
    void getUnknownToCompMap(std::map<int, int> &unknownToComp);    
    
    void getVolumeAndIntfc(std::map<int, double> &volumes,
			   std::map<int, std::vector<EBM3D_TRIANGLE> > &intfc);    
    void getTetraVolumeAndIntfc(int vertices[4], int edges[6], 
                        std::map<int, double> &volumes,
                        std::map<int, std::vector<EBM3D_TRIANGLE> > &intfc);
   
    void getAveragedIntfcAreaNormalCenter(std::map<int, std::vector<EBM3D_TRIANGLE> > &intfc,
                        std::map<int, double> &intfcArea,
                        std::map<int, EBM3D_POINT> &intfcNormal,
                        std::map<int, EBM3D_POINT> &intfcCenter);
    
    void getPartialFaces(std::map<int, std::vector<EBM3D_PARTIALFACE> > &partialfaces);
    void getPartialFaces(int corner[4], int edge[5], std::vector<EBM3D_PARTIALFACE> &partialfaces);

    // functions for parallel 
    void resetCellUnknownIndex(int ilower);    
        
    // the following functions are used for debugging or output.
    void test(void);
    static void test_checkVolumes(std::map<int, double> &volumes);
    void test_getPartialFaces(std::map<int, std::vector<EBM3D_PARTIALFACE> > &partialfaces);
    void saveIntfc_tecplot(const char *filename, std::map<int, std::vector<EBM3D_TRIANGLE> > &intfc);
    
public:
    
    void setFaceVertices(int corner[4], int v0, int v1, int v2, int v3);
    void setFaceEdges(int edge[5], int e0, int e1, int e2, int e3, int e4);
    void setTetraVertices(int vertices[4], int v0, int v1, int v2, int v3);
    void setTetraEdges(int edges[6], int e0, int e1, int e2, int e3, int e4, int e5);
        
    void getCenter(double a[3], double b[3], double center[3]);
    void getCenter(double a[3], double b[3], double c[3], double center[3]);
    void getCenter(double a[3], double b[3], double c[3], double d[3], double center[3]);

    void getComp(int comp[4], int v0, int v1, int v2, int v3);
    void getCoords(double coords[4][3], int v0, int v1, int v2, int v3);

    static double getDistance(double v0[3], double v1[3]);
    static double getDistance2(double v0[3], double v1[3]);
    static void interpolate(double v0[3], double v1[3], double t, double p[3]);
    static double getTriArea(double v0[3], double v1[3], double v2[3]);
    double getTriArea(int corner, int crx0, int crx1);
    static double getTetraVolume(double v0[3], double v1[3], double v2[3], double v3[3]);
    double getTriangularPrismVolume(double v0[3], double v1[3], double v2[3], 
                                    double v3[3], double v4[3], double v5[3]);
    
    static int *m_comp;
    static double (*m_coords)[3];
    static double (*m_crx)[3];
    
public:    
    static int m_cubeEdges[8][7];      
    static int m_cubeEdgePairs[19][2];  
};


/*******************************************************************************
 *          EBM3D_STATE
 * In each cell, there are two different states, which are stored using only one
 * EBM3D_STATE.
 *******************************************************************************/
class EBM3D_STATE {
public:
     int m_compList[2];
     int m_unknownIndex[2];            // global unknown
     double m_unknowns[2];
};

/******************************************************************************
 *          EBM3D_CARTESIAN
 * This class is the struct for the whole method. It uses the EBM_CELL as the
 * cell of the whole domain.
 *
 * This implementation uses an underlying cartesian grid to partition the
 * computational domain.
 * Note that the implementations are different than that of 2D EBM.
 *
 *****************************************************************************/
class EBM3D_CARTESIAN {
    Front *m_pFront;
    SOLVER *m_pSolver;
    EBM3D_LAPLACE *m_pLaplace;

    int m_gmax[3];            // cartesian grid.
    int m_myid;               
    int m_lbuf[3], m_ubuf[3]; // used for parallel
    int m_ilower, m_iupper;    // matrix range
    
    EBM3D_CELL ***m_cells;
    int m_nComps;             // size of m_comps[]
    int *m_comps;             // components for the vertices of the grid cells.
    int m_nLocalUnknowns;
    //double *m_unknowns;       // storages for the unknowns.    
    double ***m_unknowns[2]; 
    int ***m_compList[2];

    static int m_i;           // used for member function default parameters.
    static int m_j;           // for the current cell index.
    static int m_k;           
    double m_cellVolume;   // the minimum cell volume.
    double m_dx;
    double m_dy;
    double m_dz;
    double m_mindd;
    double m_maxdd;           // min, max of m_dx, m_dy, m_dz.

    // grids
    void getCellCenter(double center[3], int i, int j, int k);        
    void getCellFaceCenter(int face, double center[3], int i, int j, int k);      
    double getCellEdgeCenter(int dir, int i);
    double getCellEdgeLength(int dir, int i, int j, int k);        

    int getCellCornerIndex(int corner, int i, int j, int k);  
    void getCellCornerCoords(int corner, double coords[3], int i, int j, int k);
    void getCellCornerCoords(double coords[8][3], int i, int j, int k);

    void getCurrentCell(int &i, int &j, int &k);
    void setCurrentCell(int i, int j, int k);      // make sure that this is called before getCellUnknownIndex
    void locateCell(double coords[3], int &i, int &j, int &k);

    // Components
    int getCellCornerComp(int index);
    int getCellCornerComp(int corner, int i, int j, int k);

    // unknown indices    
    int getCellUnknownIndex(int comp, int i, int j, int k);    
    int getCellUnknownIndexByTrial(int comp, int i, int j, int k);      // try to get the unknown index
    int getCellNeighborUnknownIndex(int face, int comp, int i, int j, int k);
    int getNeighborByPlane(EBM_PLANE plane, int x, int y, int comp, int i, int j, int k);
    int getNeighborByPlane(EBM_PLANE plane, int x, int y, int comp, int index[3], int i, int j, int k);

    double getCellUnknown(int comp, int i, int j, int k);
    
    void getCellCompList(int comp[2], int i, int j, int k);
    void setCellCompList(int comp[2], int i, int j, int k);
    void getCellUnknown(double unknown[2], int i, int j, int k);
    void setCellUnknown(double unknown[2], int i, int j, int k);
    
    void setCellUnknown(double *x);

    // crossing & shapes
    EBM3D_CELL *getCell(int i, int j, int k);  
    CELL_TYPE getCellType(int i, int j, int k);     // ???
    
    //boolean getCellEdgeCrx(int edge, double crx[2], int i=m_i, int j=m_j, int k=m_k);  // ???
    //void getCellUnknownBoundary(int comp[9],
    //            int numCellUnknowns, int corner, std::vector< std::pair<double,double> > &points, int i=m_i, int j=m_j, int k=m_k); // ???
    void getCellComp(int &c0, int &c1, int i, int j, int k);
    void getCellComp(int comp[8], int i, int j, int k);
    
    void getCellCompCoordsCrx(int comp[8], double coords[8][3], double crx[19][3], int i, int j, int k);  // ???
    void getCellFaceCompCoords(int face, int comp[4], double coords[4][3]);
    
    // to be implemented using FronTier functions
    int getComp(double *coords);   
    void getCrx(double p0[3], double p1[3], double crx[3]);
    void setEquivalentComponents(void);        // modify components.
    
    // matrix setup
    void setMatrixForInternalCell(int i, int j, int k);
    void setMatrixForPartialCell_2Unknowns(int i, int j, int k);
    void setMatrixForBoundaryCell(int i, int j, int k);
    
    // matrix setup for intfc/boundary.
    void setMatrixForInterface(int intfcIndex[3], double intfcCenter[3], double intfcNormal[3], 
			       int intfcComp[2], double intfcArea, double volumeFraction[2], int i, int j, int k);
    void setMatrixForInterface2(int intfcIndex[3], double intfcCenter[3], double intfcNormal[3], 
			       int intfcComp[2], double intfcArea, double volumeFraction[2], int i, int j, int k);
    void setMatrixForBoundary(int intfcIndex[2], double intfcCenter[3], double intfcNormal[3], 
			      int intfcComp[2], double intfcArea, double volumeFraction[2], int i, int j, int k);
    
    // flux stencil for the partial faces
    void getFluxStencilForPartialFace(int face, double center[3], int comp,
				      std::vector< std::pair<int,double> > &stencil,
				      int i, int j, int k);
    void getFluxStencilForPartialFace2(int face, double center[3], int comp, 
				       std::vector< std::pair<int,double> > & stencil,
				       int i, int j, int k);    
    
    double getVolumeFraction(int index, std::map<int, double> &volumes);

    // flux stencil for the boundary/interface
    void getSolutionInterpStencil(double intfcCenter[3], double intfcNormal[3], int comp,
				  std::vector<int> &indices, std::vector<double> &coefficients, 
				  std::vector<double> &distances,
				  int i, int j, int k);            
    void getSolutionInterpStencil2(double intfcCenter[3], double intfcNormal[3], int comp,
				  std::vector<int> &indices, std::vector<double> &coefficients, 
				  std::vector<double> &distances,
				  int i, int j, int k);            
    void getStencilPointsAndDistance(double crx[3], double normal[3], int comp, double coords[3], double &distance, int i, int j, int k);      
    void get4PointsStencil(double crx[3], int comp, 
			   std::vector<int> &indices, std::vector<double> &coefficients);
    void getCrossingAlongNormal(double coords[3], double normal[3], int dir, int sign, double crx[3], int i, int j, int k);
        
    // interpolation
    void get2PointsSolutionInterpCoeffs(double coords[2], double x, double coeffs[2]);              // linear
    void get2PointsDerivativeInterpCoeffs(double coords[2], double x, double coeffs[2]);            // linear
    void get3PointsSolutionInterpCoeffs(double coords[3], double x, double coeffs[3]);              // quadratic
    void get3PointsDerivativeInterpCoeffs(double coords[3], double x, double coeffs[3]);            // quadratic
    void get3PointsSolutionInterpCoeffs_Plane(double coords[3][2], double x[2], double coeffs[3]);  // bilinear

    void get6PointsSolutionInterpCoeffs_Plane(EBM_PLANE plane, int comp, double crx[3], int index[6], double coeffs[6], int i, int j, int k);
    void get6PointsSolutionInterpCoeffs_Plane(double coords[6][2], double x[2], double coeffs[6]);  // biquadratic

    void get4PointsSolutionInterpCoeffs(double coords[4][3], double x[3], double coeffs[4]);   
    void get10PointsSolutionInterpCoeffs(double coords[10][3], double x[3], double coeffs[10]);     // quadratic
    
    void getBilinearInterpolationCoeff(double x1,double y1,double x2,double y2,double x,double y,double coeff[4]);

    void getCoordsByPlane(int plane, double x[3], double y[2]);

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

    // The following functions are used for the EBM methods.
    void setMatrix(void);    
public:
    EBM3D_CARTESIAN(Front *front, SOLVER *solver, EBM3D_LAPLACE *laplace);
    ~EBM3D_CARTESIAN();
    
    void solve(void);    

    // output functions
    void saveInterface_Tecplot(const char *filename);
    void saveReconstructedInterface_Tecplot(const char *filename=NULL, int i=-1, int j=-1, int k=-1); // use m_hfile in default.
    void saveInterfaceEdges_Tecplot(const char *filename);
    void savePartialCells_Tecplot(const char *filename, 
				  std::map<int, std::vector<EBM3D_PARTIALFACE> > &mapPartialFaces,
				  std::map<int, std::vector<EBM3D_TRIANGLE> > &mapIntfc,
				  int i, int j, int k);
    void saveStates_Tecplot(int i=-1, int j=-1, int k=-1);
    void saveStates_VTK(int i=-1, int j=-1, int k=-1);

    void saveComponent_Tecplot(const char *filename);
    
    // debugging functions
    void debug(std::vector< std::pair<int,double> > &stencil);
    void debug_printNeighborComp(int i, int j, int k);
    void debug_printNeighborComp(EBM_PLANE plane, int i, int j, int k);
    void debug_printCellComp(int comp[8], int i, int j, int k);
    void debug_get4PointsSolutionInterpCoeffs(void); 
    void debug_get6PointsSolutionInterpCoeffs_Plane(void);
    void debug_get10PointsSolutionInterpCoeffs(void);
    
    
private:
    void debug_printNormal(FILE *hfile, double center[3], double normal[3]);
    void setPlotRange(int range[3][2], int i=-1, int j=-1, int k=-1);
    void debug_saveUnknownIndex(void);
};


#endif	/* _EBM3D_H */


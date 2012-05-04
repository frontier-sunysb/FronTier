/*
 *          ebm2d.c
 * See also ebm2d.h
 */
#include <stdio.h>

#include "solver.h"
#include <stdlib.h>
#include <math.h>

#include <FronTier.h>
#include <FronTier/intfc/geom.h>
//#include "int.h"
//#include "fdecs.h"

#undef assign 
#undef vector 
#include <ebm2d.h>
#include <frontier_ppgrid.h>
#include <geometry.h>


/*****************************************************************************
 *      EBM2D_LAPLACE
 * The EBM2D_LAPACE is used to represent the elliptic interface problem. 
 *
 *      Div * beta Grad phi = rou
 *
 *****************************************************************************/
EBM2D_LAPLACE::EBM2D_LAPLACE()
{       
}
// beta
double EBM2D_LAPLACE::getBeta(double coords[2], int comp)
{
    if(comp==1)
        return 1;
    else if(comp==2)
        return 10;
    else
        printf("EBM2D_LAPLACE::getBeta: unknown case comp=%d\n", comp);
}
// phi
double EBM2D_LAPLACE::getExactSolution(double coords[2], int comp)
{
    double r = sqrt(coords[0]*coords[0]+coords[1]*coords[1]);    
    if(comp==1)
        return r*r*r;
    else if(comp==2)
        return r*r*r/10 + (1-1.0/10)*(1.0/8);        
    else
	 printf("EBM2D_LAPLACE::getExactSolution: unknown case comp=%d\n", comp);
}

double EBM2D_LAPLACE::getFlux(double coords[2], int comp, double normal[2])
{
    double r = sqrt(coords[0]*coords[0]+coords[1]*coords[1]);    
    return 3*coords[0]*r*normal[0]+3*coords[1]*r*normal[1];
}

// rou    
double EBM2D_LAPLACE::getRightHandSide(double coords[2], int comp)
{
    //return 0;
    double r = sqrt(coords[0]*coords[0]+coords[1]*coords[1]); 
    return 9*r;
}

/******************************************************************************
 *                  EBM2D_POINT
 *****************************************************************************/
EBM2D_POINT::EBM2D_POINT()
{
}

EBM2D_POINT::EBM2D_POINT(double coords[2])
{
     setCoords(coords);
}
void EBM2D_POINT::getCoords(double coords[2])
{
     coords[0] = m_coords[0];
     coords[1] = m_coords[1];
}
void EBM2D_POINT::setCoords(double coords[2])
{
     m_coords[0] = coords[0];
     m_coords[1] = coords[1];
}

/******************************************************************************
 *      EBM2D_EDGE
 ******************************************************************************/

/*
 * is there problem with the following two functions when used with setMatrixForBoundary()?????
 */
void EBM2D_EDGE::getComp(int comp[2])
{
     comp[0] = m_comp[0];
     comp[1] = m_comp[1];
}
void EBM2D_EDGE::getIndex(int index[2])
{
     index[0] = m_index[0];
     index[1] = m_index[1];
}
/*
 * always point from a smaller component to a bigger component!
 * return the smaller component.
 */
int EBM2D_EDGE::setEdge(int index0, int comp0, 
			  int index1, int comp1, 
			  double v0[2], double v1[2])
{
    if(comp0<comp1)
    {
        m_comp[0] = comp0;
        m_comp[1] = comp1;
        m_index[0] = index0;
        m_index[1] = index1;

        for(int i=0; i<2; i++)    
        {
            m_vertex[0][i] = v0[i];
            m_vertex[1][i] = v1[i];
	}
    }
    else    
    {
        m_comp[0] = comp1;
        m_comp[1] = comp0;
        m_index[0] = index1;
        m_index[1] = index0;

        for(int i=0; i<2; i++)    
        {
            m_vertex[0][i] = v1[i];
	    m_vertex[1][i] = v0[i];
            
        } 
    }
    return m_comp[0];   // m_comp[0] is always the smaller component.
}
/*
 * used for debugging with EBM2D_PARTIALEDGE.
 */
int EBM2D_EDGE::setEdge(double v0[2], double v1[2])
{
     for(int i=0; i<2; i++)    
     {
	  m_vertex[0][i] = v0[i];
	  m_vertex[1][i] = v1[i];         
     }
}

double EBM2D_EDGE::getLength(void)
{
    double a[2];
    difference(m_vertex[1],m_vertex[0],a,2);
    return sqrt(a[0]*a[0]+a[1]*a[1]);    
}
void EBM2D_EDGE::getNormal(double normal[2])
{
    double a[2];
    difference(m_vertex[1],m_vertex[0],a,2);
    double r = sqrt(a[0]*a[0]+a[1]*a[1]);
    
    normal[0] =  a[1]/r;
    normal[1] = -a[0]/r;  
}

void EBM2D_EDGE::getLengthNormalCenter(double &length, double normal[2], double center[2])     // ??? center ???
{
    double a[2];
    difference(m_vertex[1],m_vertex[0],a,2);
    
    length = sqrt(a[0]*a[0]+a[1]*a[1]);
    
    normal[0] =  a[1]/length;
    normal[1] = -a[0]/length;
    
    // center: what kind of center is needed?
    center[0] = (m_vertex[0][0] + m_vertex[1][0])/2;
    center[1] = (m_vertex[0][1] + m_vertex[1][1])/2;
}

/******************************************************************************
 *          EBM2D_PARTIALEDGE
 ******************************************************************************/
/*
 * first function to call in a sequence of steps to get the partial face center.
 */
void EBM2D_PARTIALEDGE::begin(int c, int index)
{
    m_comp = c;
    m_index = index;
    m_length = 0;
    m_center[0] = 0;
    m_center[1] = 0;

    debug_begin();
}
/*void EBM2D_PARTIALEDGE::begin(int c, int index, double area, double center[2])
{
    m_comp = c;
    m_index = index;
    m_area = area;
    m_center[0] = center[0] * area;
    m_center[1] = center[1] * area;
}*/
void EBM2D_PARTIALEDGE::addEdge(double v0[2], double v1[2])
{
    double tmp_length;
    double center[2];
    tmp_length = getEdgeLength(v0,v1);
    getEdgeCenter(v0,v1,center);
    
    m_length += tmp_length;
    for(int i=0; i<2; i++)
        m_center[i] += tmp_length*center[i];
    
    debug_addEdge(v0,v1);
}
void EBM2D_PARTIALEDGE::addEdge(double v0[2], double v1[2], double v2[2])
{
    addEdge(v0,v1);
    addEdge(v1,v2);
}
void EBM2D_PARTIALEDGE::end(void)
{
    for(int i=0; i<2; i++)
        m_center[i] /= m_length;
}

void EBM2D_PARTIALEDGE::debug_begin(void)
{
     m_edges.clear();
}
void EBM2D_PARTIALEDGE::debug_addEdge(double v0[2], double v1[2])
{
     EBM2D_EDGE edge;
     edge.setEdge(v0,v1);
     m_edges.push_back(edge);
}

double EBM2D_PARTIALEDGE::getEdgeLength(double v0[2], double v1[2])
{
    double a[2];
    difference(v1,v0,a,2);

    return sqrt(a[0]*a[0]+a[1]*a[1]);
}
void EBM2D_PARTIALEDGE::getEdgeCenter(double v0[2], double v1[2], double center[2])
{
    int i;
    for(i=0; i<2; i++)
        center[i] = 1.0/2*(v0[i]+v1[i]);
}


/*****************************************************************************
 *          EBM2D_CELL
 *****************************************************************************/

EBM2D_CELL::EBM2D_CELL()
{
     init();
}

void EBM2D_CELL::init(void)
{
     m_nCellUnknowns = -1;
     m_index[0]=m_index[1]=m_index[2]=m_index[3]=-1;    
}

CELL_TYPE EBM2D_CELL::getCellType(void)
{
     return m_cellType;
}

/* there are at most two different components in comp[4] for now.
 * calculate the min & max of comp[4].
 *        EXTERNAL  if max<0;
 *        BOUNDARY if min<0;
 *        PARTIAL    if min!=max;
 *        INTERNAL otherwise.
 */ 
void EBM2D_CELL::setCellType(int *comp)
{
     int min, max;
     min = max = comp[0];
     for(int i=1; i<4; i++)
	  if(comp[i]<min)
	       min = comp[i];
	  else if(comp[i]>max)
	       max = comp[i];
     if(max<0)
	  m_cellType = EXTERNAL;
     else if(min<0)
	  m_cellType = BOUNDARY;
     else if(min!=max)
	  m_cellType = PARTIAL;
     else
	  m_cellType = INTERNAL;
}


/*
 * Return the unknown indices in index[4] and also find the maximum of them.
 * It returns the maximun of all indices.
 */
int EBM2D_CELL::getCellUnknownIndex(int index[4])
{
    int max = index[0] = m_index[0];    
    for(int i=1;i<4;i++)
    {
        index[i] = m_index[i];
        if(max<m_index[i])
            max = m_index[i];
    }
    return max;
}
void EBM2D_CELL::setCellUnknownIndex(int index[4])
{
     for(int i=0; i<4; i++)
	  m_index[i] = index[i];
}

int EBM2D_CELL::getCellUnknownIndex(int corner)
{
    return m_index[corner];
}

/*
 * return the unknown numbers in the current cell.
 * These are for the unknowns stored at the cell center only.
 */
int EBM2D_CELL::getCellUnknownNumber(void)
{
    return m_nCellUnknowns;
}

/*
 * the first intfc index
 * 3D code to be changed!!!
 */
int EBM2D_CELL::getCellIntfcIndex(void)
{
     int max = m_index[0];
     for(int i=1; i<4; i++)
	  if(m_index[i]>max)
	       max = m_index[i];
    return max + 1;
}

/*
 * the following functions are used for the generalized marching cubes 
 * method to get the number of unkowns, the volume fractions, the area of 
 * the surfaces and their normal.
 */
void EBM2D_CELL::setComp(int comp[4])
{
    m_comp = comp;
}

/*
 * Comp[4]:      the components of the vertices;
 * coords[4][2]: the coordinates of the vertices;
 * crx[5][2]: the crossing on the 5 edges/diagonals.
 */
void EBM2D_CELL::setCompCoordsCrx(int comp[4], double coords[4][2], double crx[5][2])
{
    m_comp = comp;
    m_coords = coords;
    m_crx = crx;   
}


/*
 * use marching tetrahedra method. There are 2 triangle in each cell.
 * there are at most 2 components in each cell.
 * need to change the 3D code.
 * the unknown with a bigger component has smaller unknown index.
 */
int EBM2D_CELL::setCellUnknownIndex(int start_index)
{
    int i, p, index, end_index;
    index = end_index = start_index;

    /* // find the biggest non-negative component to start with
    int max_c = m_comp[0];
    int max_i = 0;
    for(i=1; i<4; i++)
	 if(m_comp[i]>=0&&max_c<m_comp[i])           // at least 1 positive component.
	 {
	      max_c = m_comp[i];
	      max_i = i;
	 }
	   
    // loop around the corner. 
    for(p=0, i=max_i; p<4; p++, i= (++i)%4)   */
    for(i=0; i<4; i++)
    {       
	 if(m_index[i]<0 && m_comp[i]>=0)
	 {
	      m_index[i] = index;
	      setCellUnknownIndex(i, index);
	      end_index = index;
	      index++;
	 }        
    }

     m_nCellUnknowns = end_index - start_index + 1;      // num of unknonws in this cell.         
     return m_nCellUnknowns;
 }

 void EBM2D_CELL::setCellUnknownIndex(int corner, int index)
 {    
     int pa, pb, j;
     for(j=0; j<3; j++)
     {
	 pa = corner;
	 pb = m_cellEdges[pa][j];
	 if(pb<0)
	     break;
	 if(m_index[pb]>=0)
	     continue;
	 if(m_comp[pa]==m_comp[pb])
	 {
	     m_index[pb] = index;
	     setCellUnknownIndex(pb, index);
	 }
     }    
 }

 /*
  * The comp[5] stores the comps at the four corners of the cell and the comp of
  * the center of the cell.  
  */
 void EBM2D_CELL::setCellUnknownIndex(int start_index, int comp[5], int numCellUnknowns)
 {    
      printf("EBM2D_CELL::setCellUnknownIndex: deprecated. \n");
     switch(numCellUnknowns)
     {
	 case 1:
	     m_index[0]=m_index[1]=m_index[2]=m_index[3] = start_index;            
	     break;
	 case 2:
	     //boolean first = TRUE;            
	     m_index[0]= start_index;
	     for(int k=1; k<4; k++)            
		 if(comp[k]!=comp[0])
		     m_index[k] = start_index+1;
		 else
		     m_index[k] = start_index;            
	     break;                        
	 case 3:     // two cases             
	     // comp[] has the following pattern
	     // 2 1
	     //  1
	     // 1 2
	     if(comp[0]==comp[5])
	     {
		 m_index[0] = m_index[2] = start_index;
		 m_index[1] = start_index + 1;
		 m_index[3] = start_index + 2;
	     }
	     // comp[] has the following pattern
	     // 2 1
	     //  2
	     // 1 2
	     else
	     {
		 m_index[0] = start_index;
		 m_index[1] = m_index[3] = start_index + 1;
		 m_index[2] = start_index + 2;
	     }            
	     break;          
     }

     m_nCellUnknowns = numCellUnknowns;
 }

 void EBM2D_CELL::getUnknownToCompMap(std::map<int, int> &unknownToComp)
 {
     for(int i=0; i<4; i++)    
	 unknownToComp[m_index[i]] = m_comp[i];    
 }


 /*
  *  get the areas and the intfcs inside the current grid block.
  *  to be tested,
  */
 void EBM2D_CELL::getAreaAndIntfc(std::map<int, double> &areas,                        
			 std::map<int, std::vector<EBM2D_EDGE> > &intfc)
 {
     int vertices[3], edges[3];

     // tri 0
     assign(vertices, 0,1,2);
     //assign(edges, 0,1,4);
     assign(edges, 2,1,4);
     getTriAreaAndIntfc(vertices,edges, areas, intfc);
     // tri 1 
     assign(vertices, 0,2,3);
     //assign(edges, 4,2,3);
     assign(edges, 4,3,0);
     getTriAreaAndIntfc(vertices,edges, areas, intfc);
 }

 void EBM2D_CELL::getTriAreaAndIntfc(int vertices[3], int edges[3],
			 std::map<int, double> &areas,
			 std::map<int, std::vector<EBM2D_EDGE> > &intfc)
 {
     int c0 = m_index[vertices[0]];
     int c1 = m_index[vertices[1]];
     int c2 = m_index[vertices[2]];
     int min_comp;

     EBM2D_EDGE edge;

     double triArea = getTriArea(m_coords[vertices[0]], m_coords[vertices[1]], m_coords[vertices[2]]);

     // 4 cases    
     if(c0==c1&&c0==c2)   // 000
     {        
	 areas[c0] += triArea;

     }
     else if(c0!=c1&&c0!=c2) // 011
     {
	 double v = getTriArea(m_coords[vertices[0]], m_crx[edges[0]], m_crx[edges[2]]);        
	 areas[c0] += v;
	 areas[c1] += triArea - v;
	 min_comp = edge.setEdge(c0,m_comp[vertices[0]],c1,m_comp[vertices[1]],m_crx[edges[0]], m_crx[edges[2]]);
	 intfc[min_comp].push_back(edge);
     }
     else if(c1!=c0&&c1!=c2) // 101
     {
	 double v = getTriArea(m_coords[vertices[1]], m_crx[edges[1]], m_crx[edges[0]]);       
	 areas[c1] += v;
	 areas[c0] += triArea - v;
	 min_comp = edge.setEdge(c1,m_comp[vertices[1]],c0,m_comp[vertices[0]],m_crx[edges[1]], m_crx[edges[0]]);
	 intfc[min_comp].push_back(edge);
     }
     else if(c2!=c0&&c2!=c1) // 110
     {
	 double v = getTriArea(m_coords[vertices[2]], m_crx[edges[2]], m_crx[edges[1]]);       
	 areas[c2] += v;
	 areas[c0] += triArea - v;
	 min_comp = edge.setEdge(c2,m_comp[vertices[2]],c0,m_comp[vertices[0]],m_crx[edges[2]], m_crx[edges[1]]);
	 intfc[min_comp].push_back(edge);
     }
     else
     {
	 printf("EBM2D_CELL::getTriAreaIntfc: something is wrong!\n");
     }
 }


/*
 * this function needs to be modified and tested.
 */
 void EBM2D_CELL::getAveragedIntfcLengthNormalCenter(std::map<int, std::vector<EBM2D_EDGE> > &intfc,
			 std::map<int, double> &intfcLength,
			 std::map<int, EBM2D_POINT> &intfcNormal,
			 std::map<int, EBM2D_POINT> &intfcCenter)
 {
     int i;
     double length, normal[2], center[2], tmp;
     double sumLength, sumNormal[2], sumCenter[2];
     EBM2D_POINT point;
     std::map<int, std::vector<EBM2D_EDGE> >::iterator iter;
     std::vector<EBM2D_EDGE>::iterator iterEdge;

     iter = intfc.begin();
     while(iter!=intfc.end())
     {        
	 sumLength = 0;
	 sumNormal[0] = sumNormal[1] = 0;
	 sumCenter[0] = sumCenter[1] = 0;
	 iterEdge = (*iter).second.begin();
	 while(iterEdge!=(*iter).second.end())
	 {
	     (*iterEdge).getLengthNormalCenter(length,normal,center);

	     sumLength += length;
	     for(int j=0; j<2; j++)            
	     {
		 sumNormal[j] += normal[j] * length;
		 sumCenter[j] += center[j] * length;
	     }
	     iterEdge++;
	 }   

	 i = (*iter).first;

	 tmp = sqrt(sumNormal[0]*sumNormal[0] + sumNormal[1]*sumNormal[1]);
	 for(int j=0; j<2; j++)
	 {
	     sumNormal[j] /= tmp;        // ??? should this be sumArea ???
	     sumCenter[j] /= sumLength;
	 }

	 intfcLength[i] = sumLength;         // wrongly used "area" before!


	 point.setCoords(sumNormal);
	 intfcNormal[i] = EBM2D_POINT(point);
	 point.setCoords(sumCenter);
	 intfcCenter[i] = EBM2D_POINT(point);

	 iter++;
     }
 }

 /*
  * determine the partial faces of a cell boundary.
  * do we need to consider the orentation of the face?
  */
 void EBM2D_CELL::getPartialEdges(std::map<int, std::vector<EBM2D_PARTIALEDGE> > &partialedges)
 {
      int corner[2];
      std::vector<EBM2D_PARTIALEDGE> edges;
      partialedges.clear();
      // south
      assign(corner, 0,1);
      edges.clear();
      //getPartialEdges(corner, 0, edges);
      getPartialEdges(corner, 2, edges);
      partialedges[FACE_SOUTH] = edges;

      // east
      assign(corner, 1,2);
      edges.clear();
      getPartialEdges(corner, 1, edges);
      partialedges[FACE_EAST] = edges;

      // north
      assign(corner, 3,2);
      edges.clear();
      //getPartialEdges(corner, 2, edges);
      getPartialEdges(corner, 3, edges);
      partialedges[FACE_NORTH] = edges;

      // west
      assign(corner, 0,3);
      edges.clear();
      //getPartialEdges(corner, 3, edges);
      getPartialEdges(corner, 0, edges);
      partialedges[FACE_WEST] = edges;
 }

 void EBM2D_CELL::getPartialEdges(int corner[2], int e, std::vector<EBM2D_PARTIALEDGE> &partialedges)
 {   
     int comp[] = {m_comp[corner[0]], m_comp[corner[1]]};
     double dx = getDistance(m_coords[corner[0]], m_coords[corner[1]]);

     EBM2D_PARTIALEDGE edge;
     double length;

     if(comp[0]==comp[1])    // 00
     {
	 edge.begin(comp[0], m_index[corner[0]]);
	 edge.addEdge(m_coords[corner[0]], m_coords[corner[1]]);
	 edge.end();        
	 partialedges.push_back(edge);
     }    
     else                 // 01
     {
	 edge.begin(comp[0], m_index[corner[0]]);
	 edge.addEdge(m_coords[corner[0]], m_crx[e]);
	 edge.end();
	 partialedges.push_back(edge);

	 edge.begin(comp[1], m_index[corner[1]]);
	 edge.addEdge(m_crx[e], m_coords[corner[1]]);
	 edge.end();
	 partialedges.push_back(edge);
     }    
 }

void EBM2D_CELL::resetCellUnknownIndex(int ilower)
{
     int i;
     for(i=0; i<4; i++)
	  if(m_index[i]>=0)
	       m_index[i] += ilower;
}

 void EBM2D_CELL::test(void)
 {
 }

 void EBM2D_CELL::test_checkAreas(std::map<int, double> &areas)
 {
     std::map<int, double>::iterator iter;
     iter = areas.begin();
     double total = 0;
     while(iter!=areas.end())
     {
	 total += (*iter).second;
	 iter++;
     }
     std::cout << "total areas is " << total << std::endl;
 }

/*
 * calculate total length for each edges.
 */
 void EBM2D_CELL::test_getPartialEdges(std::map<int, std::vector<EBM2D_PARTIALEDGE> > &partialedges)
 {
     double length = 0;
     std::map<int, std::vector<EBM2D_PARTIALEDGE> >::iterator mapIter;

     mapIter = partialedges.begin();
     printf("edge lengths:");
     while(mapIter!=partialedges.end())
     {
	 length = 0;
	 for(int i=0; i<(*mapIter).second.size(); i++)
	     length += (*mapIter).second[i].m_length;
	 printf("%f  ", length);
	 mapIter++;
     }
     printf("\n");
 }

 void EBM2D_CELL::saveIntfc_tecplot(const char *filename, std::map<int, std::vector<EBM2D_EDGE> > &intfc)
 {
     FILE *hfile;
     if ((hfile = fopen(filename,"w")) == NULL)
     {
	 (void) printf("EBM2D_CELL::saveIntfc_tecplot: "
		       "can't open %s\n",filename);
	 return;
     }

     fprintf(hfile, "VARIABLES = x y z\n");

     std::map<int, std::vector<EBM2D_EDGE> >::iterator iterMap;

     EBM2D_EDGE edge;

     int first;
     iterMap = intfc.begin();        
     while(iterMap!=intfc.end())
     {
	 first = (*iterMap).first;
	 fprintf(hfile, "# intfc patch %d \n", first);
	 for(int i=0; i<(*iterMap).second.size(); i++)
	 {
	     fprintf(hfile, "ZONE N=%d, E=%d, DATAPACKING=POINT, ZONETYPE=FELINESEG\n",
		     2, 1);
	     edge = (*iterMap).second[i];
	     for(int j=0; j<2; j++)
		 fprintf(hfile,"%-9g %-9g %-9g\n", edge.m_vertex[j][0], edge.m_vertex[j][1]);            
	     fprintf(hfile, "1 2\n");
	 }        

	 iterMap++;
     }   
     (void) fclose(hfile);
 }

 void EBM2D_CELL::assign(int corner[2], int v0, int v1)
 {
      corner[0] = v0;
      corner[1] = v1;
 }
 void EBM2D_CELL::assign(int corner[3], int v0, int v1, int v2)
 {
      corner[0] = v0;
      corner[1] = v1;
      corner[2] = v2;
 }
 void EBM2D_CELL::assign(int corner[4], int v0, int v1, int v2, int v3)
 {
      corner[0] = v0;
      corner[1] = v1;
      corner[2] = v2;
      corner[3] = v3;
 }

 void EBM2D_CELL::getCenter(double v0[2], double v1[2], double center[2])
 {
     center[0] = 1.0/2*(v0[0]+v1[0]);
     center[1] = 1.0/2*(v0[1]+v1[1]);
 }
 double EBM2D_CELL::getDistance(double v0[2], double v1[2])
 {
      double a[2];
      difference(v0,v1,a,2);
      return sqrt(a[0]*a[0]+a[1]*a[1]);
 } 
 void EBM2D_CELL::interpolate(double v0[2], double v1[2], double t, double p[2])
 {
      p[0] = (1-t)*v0[0] + t*v1[0];
      p[1] = (1-t)*v0[1] + t*v1[1];
 }
 double EBM2D_CELL::getEdgeLength(double v0[2], double v1[2])
 {
      return getDistance(v0,v1);
 }
 double EBM2D_CELL::getTriArea(double v0[2], double v1[2], double v2[2])
 {
      double a[2], b[2], area;
      difference(v1,v0,a,2);
      difference(v2,v0,b,2);

      area = 1.0/2 * fabs( (a[0]*b[1]-a[1]*b[0]) ); 
 }

 int *EBM2D_CELL::m_comp;
 double (*EBM2D_CELL::m_coords)[2];
 double (*EBM2D_CELL::m_crx)[2];

 /*
  * m_cellEdges[i]: the vertices from which there exists edges connecting the ith
  *                 vertices directly.
  */
 int EBM2D_CELL::m_cellEdges[4][3] = {{1,2,3},
				      {0,2,-1},
				      {0,1,3},
				      {0,2,-1}
				     };

 /*
  * m_cellEdgePairs[i]: the two vertices for the ith edges/diagonals.
  */
// int EBM2D_CELL::m_cellEdgePairs[5][2] = {{0,1}, {1,2}, {2,3}, {3,0}, {0,2}};
 int EBM2D_CELL::m_cellEdgePairs[5][2] = {{0,3}, {1,2}, {0,1}, {2,3}, {0,2}};

 /*****************************************************************************
  *          EBM2D_CARTESIAN
  *****************************************************************************/
// int EBM2D_CARTESIAN::m_i = -1;
// int EBM2D_CARTESIAN::m_j = -1;

 void EBM2D_CARTESIAN::getCellCenter(double center[2],int i, int j)
 {
      i -= m_lbuf[0];
      j -= m_lbuf[1];
     RECT_GRID *rect_grid = m_pFront->rect_grid;
     center[0] = cell_center(i,0,rect_grid);
     center[1] = cell_center(j,1,rect_grid);
 }

 /*
  *
  * see also EBM3D_CARTESIAN::getCellFaceCenter(). note their difference.
  */
 void EBM2D_CARTESIAN::getCellEdgeCenter(int face, double center[2], int i, int j)
 {
 #ifdef EBM_DEBUG
      if(face<0||face>3)
	   printf("EBM2D_CARTESIAN::getCellEdgeCenter: face = %d.\n", face);
 #endif
      
      i -= m_lbuf[0];
      j -= m_lbuf[1];
     RECT_GRID *rect_grid = m_pFront->rect_grid;

     switch(face)
     {
     case FACE_SOUTH:    //0:
	  center[0] = cell_center(i,0,rect_grid);
	  center[1] = cell_edge(j,1,rect_grid);
	  break;
     case FACE_NORTH:    //2:
	  center[0] = cell_center(i,0,rect_grid);
	  center[1] = cell_edge(j+1,1,rect_grid);
	  break;
     case FACE_EAST:     //1
	  center[0] = cell_edge(i+1,0,rect_grid);
	  center[1] = cell_center(j,1,rect_grid);
	  break;
     case FACE_WEST:    //3:
	  center[0] = cell_edge(i,0,rect_grid);
	  center[1] = cell_center(j,1,rect_grid);
	  break;
     }
 }

 double EBM2D_CARTESIAN::getCellEdgeCenter(int dir, int i)
 {
 #ifdef EBM_DEBUG
      if(dir<0||dir>1)
	   printf("EBM2D_CARTESIAN::getCellEdgeCenter: dir = %d.\n", dir);
 #endif
      i -= m_lbuf[dir];

     RECT_GRID *rect_grid = m_pFront->rect_grid;
     return cell_center(i,dir,rect_grid);

     /*if(dir==0)    
	 return cell_center(i,0,rect_grid);
     else
     return cell_center(i,1,rect_grid);*/
 }

 /*
  */

 double EBM2D_CARTESIAN::getCellEdgeLength(int face, int i, int j)
 {
 #ifdef EBM_DEBUG
      if(face<0||face>3)
	   printf("EBM2D_CARTESIAN::getCellEdgeLength: face = %d.\n", face);
 #endif

      i -= m_lbuf[0];
      j -= m_lbuf[1];
      
     RECT_GRID *rect_grid = m_pFront->rect_grid;
     switch(face)
     {
	 case FACE_SOUTH:
	 case FACE_NORTH:
	     return cell_width(i, 0, rect_grid);
	     break;
	 case FACE_EAST:
	 case FACE_WEST:
	     return cell_width(j, 1, rect_grid);
     }
 }

 /*
  * the order of the cell corner index is the following
  *  ....
  *  6789...
  *  012345
  */
 int EBM2D_CARTESIAN::getCellCornerIndex(int corner, int i, int j)  // corner: 0,1,2,3.
 {
 #ifdef EBM_DEBUG
      if(corner<0||corner>3)
	   printf("EBM2D_CARTESIAN::getCellCornerIndex: corner = %d \n", corner);
 #endif

      int length = m_lbuf[0]+m_gmax[0]+m_ubuf[0];
     //int base = i + j*(m_gmax[0]+1); 
     int base = i + j * (length+1);
     int additional = 0;
     switch(corner)
     {
     case 0:
	  additional = 0;
	  break;
     case 1:
	  additional = 1;
	  break;
     case 2:
	  additional = length + 2;     //m_gmax[0]+2;
	  break;
     case 3:
	  additional = length + 1;     //m_gmax[0]+1;
	  break;            
     }
     return base + additional;
 }

 void EBM2D_CARTESIAN::getCellCornerCoords(int corner, double coords[2], int i, int j)
 {
      i -= m_lbuf[0];
      j -= m_lbuf[1];

     RECT_GRID *rect_grid = m_pFront->rect_grid;
     double x,y,dx,dy;
     dx = cell_width(i,0,rect_grid);    
     dy = cell_width(j,1,rect_grid);    
     x = cell_edge(i,0,rect_grid);
     y = cell_edge(j,1,rect_grid);

      switch(corner)
      {
      case 0:
	   coords[0] = x;
	   coords[1] = y;
	   break;
      case 1:
	   coords[0] = x + dx;
	   coords[1] = y;
	   break;
      case 2:
	   coords[0] = x + dx;
	   coords[1] = y + dy;
	   break;
      case 3:
	   coords[0] = x;
	   coords[1] = y + dy;
	   break;
      }
 }
 void EBM2D_CARTESIAN::getCellCornerCoords(double coords[4][2], int i, int j)
 {
      i -= m_lbuf[0];
      j -= m_lbuf[1];

     RECT_GRID *rect_grid = m_pFront->rect_grid;
     double x,y,dx,dy;
     dx = cell_width(i,0,rect_grid);    
     dy = cell_width(j,1,rect_grid);    

     // 0
     coords[0][0] = x = cell_edge(i,0,rect_grid);
     coords[0][1] = y = cell_edge(j,1,rect_grid);
     // 1 
     coords[1][0] = x + dx;
     coords[1][1] = y;
     // 2
     coords[2][0] = x + dx;
     coords[2][1] = y + dy;
     // 3
     coords[3][0] = x;
     coords[3][1] = y + dy;
 }


 void EBM2D_CARTESIAN::getCurrentCell(int &i, int &j)
 {
      i = m_i;
      j = m_j;
 }
/*
 * make sure that this is called before getCellUnknownIndex
 */
 void EBM2D_CARTESIAN::setCurrentCell(int i, int j)
 {
      m_i = i;
      m_j = j;
 }
 void EBM2D_CARTESIAN::locateCell(double coords[2], int &i, int &j)
 {
     RECT_GRID *rect_grid = m_pFront->rect_grid;
     i = cell_index(coords[0],0,rect_grid);
     j = cell_index(coords[1],1,rect_grid);

     i += m_lbuf[0];
     j += m_lbuf[1];
 }


 int EBM2D_CARTESIAN::getCellCornerComp(int index)
 {
      return m_comps[index];
 }

 int EBM2D_CARTESIAN::getCellCornerComp(int corner, int i, int j)
 {
 #ifdef EBM_DEBUG
      if(corner<0||corner>3)
	   printf("EBM2D_CARTESIAN::getCellCornerComp: corner = %d.\n", corner);
 #endif

      int index = getCellCornerIndex(corner, i,j);
      if(index<0||index>=m_nComps)
	   return -1;
      else
	   return getCellCornerComp(index);
 }

/*
 * get the two comps for the face of the cell[i][j].
 *       | 
 *   cell| 
 *       |
 * the face lies to the right of the cell. With counter-clockwise order, we have
 * comp[0], comp[1].
 */ 
int EBM2D_CARTESIAN::getCellEdgeComp(int face, int comp[2], int i, int j)
{
     int corner[2];
     switch(face)
     {
     case FACE_SOUTH:
	  corner[0] = 0;
	  corner[1] = 1;
	  break;
     case FACE_EAST:
	  corner[0] = 1;
	  corner[1] = 2;
	  break;
     case FACE_NORTH:
	  corner[0] = 2;
	  corner[1] = 3;
	  break;
     case FACE_WEST:
	  corner[0] = 3;
	  corner[1] = 0;
	  break;
     }
     comp[0] = getCellCornerComp(corner[0],i,j);
     comp[1] = getCellCornerComp(corner[1],i,j);
}

 /*
  * neighbor: 
  *      0: south
  *      1: east
  *      2: north
  *      3: west
  * These are for the unknowns stored at the cell center only
  */
 int EBM2D_CARTESIAN::getCellUnknownIndex(int neighbor, int comp, int i, int j)
 {
 #ifdef EBM_DEBUG
      if(neighbor<0||neighbor>3)
	   printf("EBM2D_CARTESIAN::getCellUnknownIndex: neighbor = %d.\n", neighbor);
 #endif

      printf("EBM2D_CARTESIAN::getCellUnknownIndex: deprecated. \n"); 

     int corner0 = neighbor;
     int corner1 = (neighbor+1) % 4;
     int cornerIndex0 = getCellCornerIndex(corner0,i,j);    
     int cornerComp0 = getCellCornerComp(cornerIndex0);
     if(cornerComp0==comp)
	 return getCell(i,j)->getCellUnknownIndex(corner0);
     else if(getCellCornerComp(getCellCornerIndex(corner1,i,j))==comp)
	 return getCell(i,j)->getCellUnknownIndex(corner1);    
     else
     {
	 for(int k=0; k<4; k++)     
	      if(getCellCornerComp(getCellCornerIndex(k,i,j))==comp)
	     {
		 printf("EBM2D_CARTESIAN::getCellUnknownIndex: something is wrong 1!");            
		 return getCell(i,j)->getCellUnknownIndex(k) ;
	     }        
     }
     printf("EBM2D_CARTESIAN::getCellUnknownIndex: something is wrong 2!");            
 }

/*
 * use the current cell index (m_i, m_j)
 * simple logic and potential errors???
 */
 int EBM2D_CARTESIAN::getCellUnknownIndex(int comp, int i, int j)
 {   
 #ifdef EBM_DEBUG
      if(comp<0||comp>3)
	   printf("ID %d: EBM2D_CARTESIAN::getCellUnknownIndex: negative component, comp = %d.\n", 
		  m_myid, comp);
 #endif
      int index = getCellUnknownIndexByTrial(comp,i,j);
      return index;

/*      if(index>=0)
	   return index;

      printf("ID %d: EBM2D_CARTESIAN::getCellUnknownIndex: can't found comp=%d "
	     " with (%d,%d)\n", m_myid, comp, i,j);
*/    

 /*    int neighbor;
     if(i+1==m_i)
	 neighbor = 1;
     else if(i-1==m_i)
	 neighbor = 3;
     else if(j+1==m_j)
	 neighbor = 2;
     else if(j-1==m_j)
	 neighbor = 0;
     else if(i<m_i)
	 neighbor = 1;
     else if(i>m_i)
	 neighbor = 3;
     else if(j<m_j)
	 neighbor = 2;
     else if(j>m_j)
	 neighbor = 0;
	 return getCellUnknownIndex(neighbor,comp,i,j); */

 }

 /*
  * try to get the unknown index.
  * return 
  *         >=0 if succeed;
  *         <0  if failed.
  * This function is the same as EBM3D_CARTESIAN::getCellUnknownIndex except for the return if it failed.
  */
 int EBM2D_CARTESIAN::getCellUnknownIndexByTrial(int comp, int i, int j)
 {
      //return getCell(i,j,k)->getCellUnknownIndex(0);
      if(i<0 || j<0 || 
	 i>=m_lbuf[0]+m_gmax[0]+m_ubuf[0] || 
	 j>=m_lbuf[1]+m_gmax[1]+m_ubuf[1])
	   return -1;

      int corner, index, debug[4];
      for(corner=0; corner<4; corner++)
      {
	   index = getCellCornerIndex(corner, i,j);
	   debug[corner] = getCellCornerComp(index);
	   if(debug[corner]==comp)
		return getCell(i,j)->getCellUnknownIndex(corner);
      }
      return -1;
 }


 int EBM2D_CARTESIAN::getCellNeighborUnknownIndex(int face, int comp, int i, int j)
 {
     int I, J;
     I = i;
     J = j;

     switch(face-face%2)
     {
	 case PLANE_X:
	     if(face==FACE_WEST)     // west      
		 I = i-1;
	     else                    // east
		 I = i+1;            
	     break;
	 case PLANE_Y:
	     if(face==FACE_SOUTH)    // south
		 J = j-1;
	     else                    // north
		 J = j+1;
	     break;
     }
     return getCellUnknownIndex(comp, I,J);
 }

 /*
  * find the neighbor along one of the two coordinate planes: x, y.
  * parameter:
  *      plane: PLANE_X, PLANE_Y.
  *      x: -1,0,1, the difference along the two coordinates on the plane.
  *      comp: the current component.
  */ 
 int EBM2D_CARTESIAN::getNeighborByPlane(EBM_PLANE plane, int x, int comp, int i, int j)
 {   
     switch(plane)
     {
	 case PLANE_X:   // east, west
	     j += x;
	     break;
	 case PLANE_Y:   // south, north
	     i += x;
	     break;
     }
     return getCellUnknownIndex(comp, i,j);
 }

int EBM2D_CARTESIAN::getNeighborAcrossPlane(EBM_PLANE plane, int x, int comp, int i, int j)
{
     switch(plane)
     {
	 case PLANE_X:   // east, west
	     i += x;
	     break;
	 case PLANE_Y:   // south, north
	     j += x;
	     break;
     }
     return getCellUnknownIndex(comp, i,j);
}


void EBM2D_CARTESIAN::getCellCompList(int comp[2], int i, int j)
{
     comp[0] = m_compList[0][i][j];
     comp[1] = m_compList[1][i][j];
}
void EBM2D_CARTESIAN::setCellCompList(int comp[2], int i, int j)
{
     m_compList[0][i][j] = comp[0];
     m_compList[1][i][j] = comp[1];
}

/*
 * used to get the knowns at cell (i,j)
 */
double EBM2D_CARTESIAN::getCellUnknown(int comp, int i, int j)
{
     if(comp==m_compList[0][i][j])
	  return m_unknowns[0][i][j];
     else if(comp==m_compList[1][i][j])
	  return m_unknowns[1][i][j];
     
     printf("ID %d: EBM2D_CARTESIAN::getCellUnknown: comp=%d, i=%d,j=%d.\n",
	    m_myid, comp, i,j);
     printf("ID %d:\tm_compList[*][%d][%d]={%d,%d}\n",
	    m_myid,i,j,m_compList[0][i][j],m_compList[1][i][j]);
	    
     return HUGE_VAL;
}

double EBM2D_CARTESIAN::getCellUnknownAcrossPlane(EBM_PLANE plane, 
						  int x, 
						  int comp, 
						  int i, int j)
{
     switch(plane)
     {
	 case PLANE_X:   // east, west
	     i += x;
	     break;
	 case PLANE_Y:   // south, north
	     j += x;
	     break;
     }
     return getCellUnknown(comp, i,j);
}

void EBM2D_CARTESIAN::getCellUnknown(double unknown[2], int i, int j)
{
     unknown[0] = m_unknowns[0][i][j];
     unknown[1] = m_unknowns[1][i][j];
}
void EBM2D_CARTESIAN::setCellUnknown(double unknown[2], int i, int j)
{
     m_unknowns[0][i][j] = unknown[0];
     m_unknowns[1][i][j] = unknown[1];
}

/*
 * should this function be moved into EBM2D_CELL?
 */
void EBM2D_CARTESIAN::setCellUnknown(double *x)
{
     int index, c0, c1, numCellUnknowns;
     CELL_TYPE type;

     index = 0;
     for(int j=m_lbuf[1]; j<m_lbuf[1]+m_gmax[1]; j++)
     for(int i=m_lbuf[0]; i<m_lbuf[0]+m_gmax[0]; i++)
     {   
	  type = getCell(i,j)->getCellType();
	  numCellUnknowns = getCell(i,j)->getCellUnknownNumber();
	  getCellComp(c0, c1, i, j);
	  m_compList[0][i][j] = c0;
	  m_compList[1][i][j] = c1;
	  switch(type)
	  {
	  case INTERNAL:   // numCellUnknowns == 1
	       m_unknowns[0][i][j] = x[index];
	       break;
	  case PARTIAL:
	       m_unknowns[0][i][j] = x[index];
	       m_unknowns[1][i][j] = x[index+1];
	       break;
	  case BOUNDARY:
	       if(c0>=0)
		    m_unknowns[0][i][j] = x[index];
	       else
		    m_unknowns[1][i][j] = x[index];
	       break;
	  }

	  if(numCellUnknowns>0)
	  {
	       numCellUnknowns = numCellUnknowns * 2 - 1;
	       index += numCellUnknowns;
	  }
     }         
}


 EBM2D_CELL *EBM2D_CARTESIAN::getCell(int i, int j)
 {
     return &m_cells[i][j];
 }
 CELL_TYPE EBM2D_CARTESIAN::getCellType(int i, int j)
 {
      getCell(i,j)->getCellType();

      /* Table *T = table_of_interface(m_pFront->emb_grid_intfc);
     BLK_EDGE	**blk_edge = T->blk_edge;
     return blk_edge[i][j].ctype;*/
 }

 /*
  * Get the number of unknowns in the cell(i,j) stored at the cell center only.
  *
  */
 int EBM2D_CARTESIAN::getCellUnknownNumber(int comp[5], int i, int j)
 {    
      static boolean bPrinted = FALSE;
      if(!bPrinted)
      {
	   bPrinted = TRUE;
	   printf("EBM2D_CARTESIAN::getCellUnknownNumber: deprecated.\n");
      }
     int k;
     for(int k=0; k<4; k++)
	  comp[k] = getCellCornerComp(getCellCornerIndex(k,i,j));
     int compChange = 0;
     for(k=1; k<4; k++)
	 if(comp[k]!=comp[k-1])
	     compChange++;
     if(comp[3]!=comp[0])
	 compChange++;
     // compChange: 0,2,4
     compChange = compChange/2+1;       // return 1,2,3

     if(compChange==3)
     {
	 double center[2];
	 getCellCenter(center,i,j);
	 comp[4] = getComp(center);
     }

     return compChange;      
 }


 /*
  * _BLK_EDGE::edge[][] has the following ordering
  *            (1,1)
  *      (0,0)       (0,1)
  *            (1,0)
  * In this function, (i,j) denotes the current cell index; (edgeI,edgeJ) denotes
  * the index of the edges in cell (i,j).
  */ 
 boolean EBM2D_CARTESIAN::getCellEdgeCrx(int face, double crx[2], int i, int j)  
 {
      i -= m_lbuf[0];
      j -= m_lbuf[1];


     int edgeI, edgeJ;
     switch(face)
     {
     case FACE_SOUTH:
	  edgeI = 1;  edgeJ = 0;  break;            
     case FACE_EAST:
	  edgeI = 0;  edgeJ = 1;  break;
     case FACE_NORTH:
	  edgeI = 1;  edgeJ = 1;  break;
     case FACE_WEST:
	  edgeI = 0;  edgeJ = 0;  break;            
     }

     Table *T = table_of_interface(m_pFront->emb_grid_intfc);
     BLK_EDGE	**blk_edge = T->blk_edge;

     // check whether we have a valid crxing.
     if(blk_edge[i][j].edge[edgeI][edgeJ].etype==FULL_EDGE)
	 return FALSE;

     // return the crxing since it is valid.
     double *cen, length;
     cen = blk_edge[i][j].edge[edgeI][edgeJ].cen;    
     length = blk_edge[i][j].edge[edgeI][edgeJ].length;

     RECT_GRID *rect_grid = m_pFront->rect_grid;
     switch(face)
     {
     case FACE_SOUTH:
     case FACE_NORTH:
	  crx[0] = cen[0] + (cen[0]<cell_center(i,0,rect_grid)? 1 : -1)*length/2;
	  crx[1] = cen[1];
	  break;
     case FACE_EAST:
     case FACE_WEST:
	  crx[0] = cen[0];
	  crx[1] = cen[1] + (cen[1]<cell_center(j,1,rect_grid)? 1 : -1)*length/2;
	  break;
     }    
     return TRUE;    // yes, we have a crxing.
 }

/* 
 * get the two distinct components stored at the corners of the cell. The comps
 * can be negative.
 */ 
void EBM2D_CARTESIAN::getCellComp(int &c0, int &c1, int i, int j)
{
     int index = getCellCornerIndex(0, i,j);
     c0 = m_comps[index];
     for(int corner=1; corner<4; corner++)
     {
	  index = getCellCornerIndex(corner,i,j);
	  if(index<0||index>m_nComps)
	       c1 = -1;
	  else
	       c1 = m_comps[index];
	  if(c1!=c0)
	       return;
     }
}
 void EBM2D_CARTESIAN::getCellComp(int comp[4], int i, int j)
 {
     int index[4];     
     for(int corner=0; corner<4; corner++)
     {
	 index[corner] = getCellCornerIndex(corner,i,j);
	 if(index[corner]<0||index[corner]>=m_nComps)
	      comp[corner] = -1;
	 else
	      comp[corner] = m_comps[index[corner]];
     }
 }

/*
 * 372
 * 4 5
 * 061
 */
 void EBM2D_CARTESIAN::getCellCompAndCrx(int comp[5], double coords[8][2], int i, int j)
 { 
      static boolean bPrinted = FALSE;
      if(!bPrinted)
      {
	   bPrinted = TRUE;
	   printf("EBM2D_CARTESIAN::getCellCompAndCrx: deprecated.\n");
      }

      getCellUnknownNumber(comp,i,j);                                
      getCellCornerCoords(coords, i, j);
     
      /*for(int k=4; k<8; k++)
	  if(!getCellEdgeCrx(k-4,coords[k],i,j))
	 {
	     coords[k][0] = HUGE_VAL;
	     coords[k][1] = HUGE_VAL;
	     }*/
      getCrx(coords[0], coords[1], coords[6]);
      getCrx(coords[1], coords[2], coords[5]);
      getCrx(coords[2], coords[3], coords[7]);
      getCrx(coords[3], coords[0], coords[4]);
 }
 void EBM2D_CARTESIAN::getCellCompCoordsCrx(int comp[4], 
					    double coords[4][2], 
					    double crx[5][2], 
					    int i, int j)
 {
      int index[4];  

     // comp
     getCellComp(comp,i,j);

     // coords
     getCellCornerCoords(coords,i,j);

     // two method to get the edge crxing.
     // the first method use information from FronTier
     // the second method use bipartion method.

     // get the crx using FronTier information.
     
     boolean bCrx[4] = {FALSE,FALSE,FALSE,FALSE};
     if(comp[0]!=comp[3])
     {
	  getCellEdgeCrx(FACE_WEST,crx[0],i,j);
	  bCrx[0] = TRUE;
     }
     if(comp[1]!=comp[2])
     {
	  getCellEdgeCrx(FACE_EAST,crx[1],i,j);
	  bCrx[1] = TRUE;
     }
     if(comp[0]!=comp[1])
     {
	  getCellEdgeCrx(FACE_SOUTH,crx[2],i,j);
	  bCrx[2] = TRUE;
     }
     if(comp[3]!=comp[2])
     {
	  getCellEdgeCrx(FACE_NORTH,crx[3],i,j);
	  bCrx[3] = TRUE;
     }

     if(bCrx[0]&&bCrx[2])
     	  GEOMETRY::intersectLineLine2D(crx[0],crx[2],coords[0],coords[2],crx[4]);
     else if(bCrx[0]&&bCrx[1])
	  GEOMETRY::intersectLineLine2D(crx[0],crx[1],coords[0],coords[2],crx[4]);
     else if(bCrx[3]&&bCrx[2])
	  GEOMETRY::intersectLineLine2D(crx[3],crx[2],coords[0],coords[2],crx[4]);
     else if(bCrx[3]&&bCrx[1])
	  GEOMETRY::intersectLineLine2D(crx[3],crx[1],coords[0],coords[2],crx[4]);

     return;
     

     // get crx using bipartition
     int pa, pb;
     for(int m=0; m<5; m++)
     {        
	 pa = EBM2D_CELL::m_cellEdgePairs[m][0];
	 pb = EBM2D_CELL::m_cellEdgePairs[m][1];

	 if(comp[pa]!=comp[pb])        
	     getCrx(coords[pa], coords[pb], crx[m]);
	 else
	 {
	      crx[m][0] = crx[m][1] = HUGE_VAL;
	 }
	 } 
 }


/*
 * to be implemented using FronTier functions
 */
 int EBM2D_CARTESIAN::getComp(double *coords)
 {
      int comp = component(coords, m_pFront->interf);
      return comp;

      /*double min = -0.7, max = 0.7;
      if( coords[0]<min || coords[0]>max || 
	  coords[1]<min || coords[1]>max)
	  return -1;
      else
      return 1;*/
      
      //     return 1;

      const double R = 0.5;
     const double Center[2] = {0, 0};
     double distance = sqrt(sqr(coords[0] - Center[0]) + sqr(coords[1] - Center[1])) 
	  - R;

     if(distance<m_mindd*m_mindd)
	 return 1;
     else
	 return -2;
 }
 void EBM2D_CARTESIAN::getCrx(double p0[2], double p1[2], double crx[2])
 {
     double a[2], b[2];

     a[0] = p0[0];   a[1] = p0[1];
     b[0] = p1[0];   b[1] = p1[1];

     int ca, cb, cc;
     ca = getComp(a);
     cb = getComp(b);
     if(ca==cb)
     {
	 crx[0] = HUGE_VAL;
	 crx[1] = HUGE_VAL;
	 return;
     }

     for(int i=0; i<m_MAX_ITER; i++)
     {
	 crx[0] = 0.5*(a[0]+b[0]);
	 crx[1] = 0.5*(a[1]+b[1]);
	 cc = getComp(crx);
	 if(cc==ca)
	 {
	     a[0] = crx[0];
	     a[1] = crx[1];
	 }
	 else
	 {
	     b[0] = crx[0];
	     b[1] = crx[1];
	 }
     }

     double length = EBM2D_CELL::getDistance(p0,p1);
     double length0 = EBM2D_CELL::getDistance(p0,crx);
     double t = length0/length;
     //double T = 1.0/100;
     double T = m_mindd4; 
     if(t<T)
     {
	  EBM2D_CELL::interpolate(p0,p1,T,crx);
	  printf("EBM2D_CARTESIAN::getCrx: t<T.\n");
     }
     if((1-t)<T)
     {
	  EBM2D_CELL::interpolate(p0,p1,1-T,crx);
	  printf("EBM2D_CARTESIAN::getCrx: 1-t<T.\n");
     }
 }

/*
 * modify components.
 */
void EBM2D_CARTESIAN::setEquivalentComponents(void)
{
    for(int j=0; j<=m_lbuf[1]+m_gmax[1]+m_ubuf[1]; j++)
        for(int i=0; i<=m_lbuf[0]+m_gmax[0]+m_ubuf[0]; i++)   
        {
	     if(m_comps[getCellCornerIndex(0,i,j)]==0)
		  m_comps[getCellCornerIndex(0,i,j)] = 2;
        }    
}

/*
 * the following also deals with the Dirichlet boundary condition.
 */
 void EBM2D_CARTESIAN::setMatrixForInternalCell(int i, int j)
 {      
      //printf("EBM2D_CARTESIAN::setMatrixForInternalCell: (%d,%d)\n", i,j);
//      double fraction = m_cellArea;

     int comp = getCellCornerComp(getCellCornerIndex(0,i,j));

     int index0, index1, index2, index3, index;
     double cellCenter[2], edgeCenter[4][2], beta;
     double dx, dy, rhs, value;

     dx = getCellEdgeLength(0,i,j);      // uniform edge length only
     dy = getCellEdgeLength(1,i,j);      // uniform edge length only    

     // indices for the cells
     //index  = getCell(i,j)->getCellUnknownIndex(0);
     //index0 = (j!=0) ? getCellUnknownIndex(comp, i,j-1) : -1;
     //index1 = (i!=(m_gmax[0]-1)) ? getCellUnknownIndex(comp, i+1,j) : -1;
     //index2 = (j!=(m_gmax[1]-1)) ? getCellUnknownIndex(comp, i,j+1) : -1;
     //index3 = (i!=0) ? getCellUnknownIndex(comp, i-1,j) : -1;      
     index  = getCell(i,j)->getCellUnknownIndex(0);
     index2 = getCellUnknownIndex(comp, i,j-1);
     index1 = getCellUnknownIndex(comp, i+1,j);
     index3 = getCellUnknownIndex(comp, i,j+1);
     index0 = getCellUnknownIndex(comp, i-1,j);

     // cell center
     getCellCenter(cellCenter,i,j);
     // edge centers
     for(int k=0; k<4; k++)
	  getCellEdgeCenter(k,edgeCenter[k],i,j);             


     // now for each edges, setup the system of equations
     // edge 0, FACE_WEST
     beta = m_pLaplace->getBeta(edgeCenter[0], comp);
     if(index0>=0)
     {
	 m_pSolver->Add_A(index, index0, beta*dy/dx);
	 m_pSolver->Add_A(index, index, -beta*dy/dx);
     }
     else
     {
	 value = m_pLaplace->getExactSolution(edgeCenter[0], comp);
	 m_pSolver->Add_A(index, index, -2*beta*dy/dx);
	 m_pSolver->Add_b(index, -2*value*beta*dy/dx); 
     }

     // edge 1, FACE_EAST
     beta = m_pLaplace->getBeta(edgeCenter[1], comp);
     if(index1>=0)
     {            
	 m_pSolver->Add_A(index, index1, beta*dy/dx);
	 m_pSolver->Add_A(index, index, -beta*dy/dx);
     }
     else
     {
	 value = m_pLaplace->getExactSolution(edgeCenter[1], comp);
	 m_pSolver->Add_A(index, index, -2*beta*dy/dx);
	 m_pSolver->Add_b(index, -2*value*beta*dy/dx); 
     }
     // edge 2, FACE_SOUTH
     beta = m_pLaplace->getBeta(edgeCenter[2], comp);
     if(index2>=0)
     {
	 m_pSolver->Add_A(index, index2, beta*dx/dy);
	 m_pSolver->Add_A(index, index, -beta*dx/dy);
     }
     else
     {
	 value = m_pLaplace->getExactSolution(edgeCenter[2], comp);
	 m_pSolver->Add_A(index, index, -2*beta*dx/dy);
	 m_pSolver->Add_b(index, -2*value*beta*dx/dy); 
     }
     // edge 3, FACE_NORTH
     beta = m_pLaplace->getBeta(edgeCenter[3], comp);
     if(index3>=0)
     {
	 m_pSolver->Add_A(index, index3, beta*dx/dy);
	 m_pSolver->Add_A(index, index, -beta*dx/dy);
     }
     else
     {
	 value = m_pLaplace->getExactSolution(edgeCenter[3], comp);
	 m_pSolver->Add_A(index, index, -2*beta*dx/dy);
	 m_pSolver->Add_b(index, -2*value*beta*dx/dy); 
     }
     // rhs
     rhs = m_pLaplace->getRightHandSide(cellCenter, comp);
     m_pSolver->Add_b(index, rhs*dx*dy);           

 }


 /*
  * setup the matrix when there are two unknowns stored at the cell center.
  */
 void EBM2D_CARTESIAN::setMatrixForPartialCell_2Unknowns(int i, int j)
 {
      //printf("EBM2D_CARTESIAN::setMatrixForPartialCell_2Unknowns: (%d,%d)\n",i,j);
     int comp[4], size, e;
     double coords[4][2], crx[5][2], beta;   

     getCellCompCoordsCrx(comp,coords,crx,i,j);    
     getCell(i,j)->setCompCoordsCrx(comp,coords,crx);

     // intfc
     std::map<int, double> areas;
     std::map<int, std::vector<EBM2D_EDGE> > mapIntfc;  
     std::map<int, double> mapIntfcLength;
     std::map<int, EBM2D_POINT> mapIntfcNormal;
     std::map<int, EBM2D_POINT> mapIntfcCenter;

     getCell(i,j)->getAreaAndIntfc(areas, mapIntfc);
     getCell(i,j)->getAveragedIntfcLengthNormalCenter(mapIntfc, mapIntfcLength, 
						      mapIntfcNormal, mapIntfcCenter);

     if(mapIntfcLength.size()>1)
     {
	 printf("EBM2D_CARTESIAN::setMatrixForPartialCell_2Unknowns: "
		"intfcLength.size()>1 is undelt case!");
	 exit(0);
     }

     std::map<int, std::vector<EBM2D_EDGE> >::iterator iterMapIntfc;
     EBM2D_EDGE edge;
     EBM2D_POINT point;

     int index = getCell(i,j)->getCellIntfcIndex();    // the first intfc index.
     int intfcIndex[3], intfcComp[2], p;
     double intfcCenter[2], intfcLength, intfcNormal[2], areaFraction[2];

     iterMapIntfc = mapIntfc.begin();
     while(iterMapIntfc!=mapIntfc.end())
     {        
	 edge = (*iterMapIntfc).second[0];
	 edge.getIndex(intfcIndex);        
	 intfcIndex[2] = index++;
	 edge.getComp(intfcComp);

	 p = (*iterMapIntfc).first;
	 intfcLength = mapIntfcLength[p];
	 mapIntfcNormal[p].getCoords(intfcNormal);
	 mapIntfcCenter[p].getCoords(intfcCenter);

	 areaFraction[0] = getAreaFraction(intfcIndex[0], areas);
	 areaFraction[1] = getAreaFraction(intfcIndex[1], areas);
	
	 setMatrixForInterface(intfcIndex, intfcCenter, intfcNormal, 
			       intfcComp, intfcLength, areaFraction, i,j);

	 iterMapIntfc++;
     }    

     // the two partial areas    
     std::map<int, double>::iterator iterAreas;
     std::map<int, int> unknownToComp;    
     getCell(i,j)->getUnknownToCompMap(unknownToComp);

     int c;
     double vol, rhs, cellCenter[2];
     getCellCenter(cellCenter,i,j);

     iterAreas = areas.begin();
     while(iterAreas!=areas.end())
     {
	 index = (*iterAreas).first;
	 c = unknownToComp[index];
	 vol = (*iterAreas).second;
	 rhs = m_pLaplace->getRightHandSide(cellCenter, c);
	 
	 areaFraction[0] = getAreaFraction(index, areas);
	 m_pSolver->Add_b(index, rhs*vol*areaFraction[0]);   

	 iterAreas++;
     }


     // partial edges
     std::map<int, std::vector<EBM2D_PARTIALEDGE> > mapPartialEdges;
     std::map<int, std::vector<EBM2D_PARTIALEDGE> >::iterator iterMapPartialEdges;
     std::vector<EBM2D_PARTIALEDGE>::iterator iterPartialEdge;
     std::vector< std::pair<int,double> > stencil;
     std::vector< std::pair<int,double> >::iterator iterStencil;
     EBM2D_PARTIALEDGE partialedge;

     getCell(i,j)->getPartialEdges(mapPartialEdges);

     iterMapPartialEdges = mapPartialEdges.begin();
     while(iterMapPartialEdges!=mapPartialEdges.end())
     {
	 e = (*iterMapPartialEdges).first;
	 iterPartialEdge = (*iterMapPartialEdges).second.begin();
	 while(iterPartialEdge!=(*iterMapPartialEdges).second.end())
	 {
	     partialedge = (*iterPartialEdge);            
	     stencil.clear(); 
	     if((*iterMapPartialEdges).second.size()==1)
		  getFluxStencilForPartialEdge(e, partialedge.m_center, partialedge.m_length, 
					       partialedge.m_comp, stencil, i,j);            
	     else
		  getFluxStencilForPartialEdge2(e, partialedge.m_center, partialedge.m_length, 
						partialedge.m_comp, stencil, i,j);            

	     beta = m_pLaplace->getBeta(partialedge.m_center, partialedge.m_comp);

	     areaFraction[0] = getAreaFraction(partialedge.m_index, areas);
	     iterStencil = stencil.begin();            
	     while(iterStencil!=stencil.end())
	     {
		 m_pSolver->Add_A(partialedge.m_index, (*iterStencil).first, 
				  beta*partialedge.m_length * (*iterStencil).second * areaFraction[0]);        
		 iterStencil++;
	     }             
	     iterPartialEdge++;
	 }
	 iterMapPartialEdges++;        
     }


 }

 /*
  * now for simplicity, only Neumann boundary condition (given) is delt with.
  */
 void EBM2D_CARTESIAN::setMatrixForBoundaryCell(int i, int j)
 {
      //printf("EBM2D_CARTESIAN::setMatrixForBoundaryCell: (%d,%d)\n", i,j);

     int comp[4], size, e;
     double coords[4][2], crx[5][2], beta;   

     getCellCompCoordsCrx(comp,coords,crx,i,j);    
     getCell(i,j)->setCompCoordsCrx(comp,coords,crx);
     
     // boundary
     std::map<int, double> areas;
     std::map<int, std::vector<EBM2D_EDGE> > mapIntfc;  
     std::map<int, double> mapIntfcLength;
     std::map<int, EBM2D_POINT> mapIntfcNormal;
     std::map<int, EBM2D_POINT> mapIntfcCenter;

     getCell(i,j)->getAreaAndIntfc(areas, mapIntfc);
     getCell(i,j)->getAveragedIntfcLengthNormalCenter(mapIntfc, mapIntfcLength, mapIntfcNormal, mapIntfcCenter);

     if(mapIntfcLength.size()>1)
     {
	 printf("EBM2D_CARTESIAN::setMatrixForBoundaryCell: intfcArea.size()>1 is undelt case!");
	 exit(0);
     }

     std::map<int, std::vector<EBM2D_EDGE> >::iterator iterMapIntfc;
     EBM2D_EDGE edge;
     EBM2D_POINT point;

     int index;
     int intfcIndex[2], intfcComp[2], p;
     double intfcCenter[3], intfcLength, intfcNormal[3], areaFraction[2];

     iterMapIntfc = mapIntfc.begin();
     while(iterMapIntfc!=mapIntfc.end())
     {        
	 edge = (*iterMapIntfc).second[0];
	 edge.getIndex(intfcIndex);          
	 edge.getComp(intfcComp);

	 p = (*iterMapIntfc).first;
	 intfcLength = mapIntfcLength[p];
	 mapIntfcNormal[p].getCoords(intfcNormal);
	 mapIntfcCenter[p].getCoords(intfcCenter);

	 areaFraction[0] = getAreaFraction(intfcIndex[0], areas);
	 areaFraction[1] = getAreaFraction(intfcIndex[1], areas);

	 setMatrixForBoundary_Dirichlet(intfcIndex, intfcCenter, intfcNormal, intfcComp, intfcLength, areaFraction, i,j);
	 //setMatrixForBoundary_Neumann(intfcIndex, intfcCenter, intfcNormal, intfcComp, intfcLength, areaFraction, i,j);

	 iterMapIntfc++;
     }    

     // the partial areas
     std::map<int, double>::iterator iterAreas;
     std::map<int, int> unknownToComp;    
     getCell(i,j)->getUnknownToCompMap(unknownToComp);

     int c;
     double vol, rhs, cellCenter[2];
     getCellCenter(cellCenter,i,j);

     iterAreas = areas.begin();
     while(iterAreas!=areas.end())
     {
	 index = (*iterAreas).first;
	 c = unknownToComp[index];
	 if(c>=0)               // inside the domain
	 {     
	      vol = (*iterAreas).second;
	      rhs = m_pLaplace->getRightHandSide(cellCenter, c);

	      areaFraction[0] = getAreaFraction(index, areas);
	      m_pSolver->Add_b(index, rhs*vol*areaFraction[0]);
	 }
	 iterAreas++;
     } 

     // partial edges
     std::map<int, std::vector<EBM2D_PARTIALEDGE> > mapPartialEdges;
     std::map<int, std::vector<EBM2D_PARTIALEDGE> >::iterator iterMapPartialEdges;
     std::vector<EBM2D_PARTIALEDGE>::iterator iterPartialEdge;
     std::vector< std::pair<int,double> > stencil;
     std::vector< std::pair<int,double> >::iterator iterStencil;
     EBM2D_PARTIALEDGE partialedge;

     getCell(i,j)->getPartialEdges(mapPartialEdges);

     iterMapPartialEdges = mapPartialEdges.begin();
     while(iterMapPartialEdges!=mapPartialEdges.end())
     {
	 e = (*iterMapPartialEdges).first;
	 iterPartialEdge = (*iterMapPartialEdges).second.begin();
	 while(iterPartialEdge!=(*iterMapPartialEdges).second.end())
	 {
	     partialedge = (*iterPartialEdge);            
	     if(partialedge.m_comp>=0)
	     {
		  stencil.clear(); 
		  
		  if((*iterMapPartialEdges).second.size()==1)
		       getFluxStencilForPartialEdge(e, partialedge.m_center, partialedge.m_length, 
						    partialedge.m_comp, stencil, i,j);            
		  else
		       getFluxStencilForPartialEdge2(e, partialedge.m_center, partialedge.m_length, 
						     partialedge.m_comp, stencil, i,j);            
		  beta = m_pLaplace->getBeta(partialedge.m_center, partialedge.m_comp);
		  
		  areaFraction[0] = getAreaFraction(partialedge.m_index, areas);
		  iterStencil = stencil.begin();            
		  while(iterStencil!=stencil.end())
		  {
		       m_pSolver->Add_A(partialedge.m_index, (*iterStencil).first, 
					beta*partialedge.m_length * (*iterStencil).second * areaFraction[0]);                
		       iterStencil++;
		  }             
	     }
	     iterPartialEdge++;
	 }
	 iterMapPartialEdges++;        
     }

     
 }


/*
 * setup the matrix when there are three unknowns stored at the cell center.
 */
void EBM2D_CARTESIAN::setMatrixForPartialCell_3Unknowns(int i, int j)
{
     printf("EBM2D_CARTESIAN::setMatrixForPartialCell_3Unknowns: undefined!\n");
}

void EBM2D_CARTESIAN::setMatrixForPartialCell_2Unknowns_2ndOrder(int i, int j)
{
     int comp[5], index[5];
     double coords[8][2];
     getCellCompAndCrx(comp,coords,i,j);
     // index[4] stores the unknown on the interface.
     index[4] = getCell(i,j)->getCellUnknownIndex(index) + 1;    
     //for(int k=0; k<5; k++)
     //  if(index[k]==25)
     //       printf(" EBM2D_CARTESIAN::setMatrixForPartialCell_2Unknowns_2ndOrder: debug here!\n");

     if(comp[1]==comp[2]&&comp[1]==comp[3])  // 0111
     {
	  //setMatrixForPartialCell_2Unknowns_Case0111(i,j,comp,index,coords);
	  //printf("0111\n");
	  setMatrixForPartialCell_2Unknowns_byRotation(1,comp,index,coords,i,j);
     }    
     else if(comp[0]==comp[2]&&comp[0]==comp[3])   // 0100
     {
	  //setMatrixForPartialCell_2Unknowns_Case0100(i,j,comp,index,coords);
	  //printf("0100\n");
	  setMatrixForPartialCell_2Unknowns_byRotation(2,comp,index,coords,i,j);
     }
     else if(comp[0]==comp[1]&&comp[0]==comp[3])   // 0010
     {
	  //setMatrixForPartialCell_2Unknowns_Case0010(i,j,comp,index,coords);
	  //printf("0010\n");
	  setMatrixForPartialCell_2Unknowns_byRotation(3,comp,index,coords,i,j);
     }
     else if(comp[0]==comp[1]&&comp[0]==comp[2])   // 0001
     {
	  //setMatrixForPartialCell_2Unknowns_Case0001(i,j,comp,index,coords);
	  setMatrixForPartialCell_2Unknowns_byRotation(0,comp,index,coords,i,j);
     }
     else if(comp[0]==comp[3]&&comp[1]==comp[2])   // 0110
     {
	  setMatrixForPartialCell_2Unknowns_Case0110(comp,index,coords,i,j);
     }
     else if(comp[0]==comp[1]&&comp[2]==comp[3])       // 0011
     {
	  setMatrixForPartialCell_2Unknowns_Case0011(comp,index,coords,i,j);
     }
     else
	 printf("EBM2D_CARTESIAN::setMatrixForPartialCell_2Unknowns. ?\n");
}

/* 
 * Going along the edges of the cells counter-clockwise.
 * the intfc edge inside the cell has index intfcIndex.
 */
 void EBM2D_CARTESIAN::setMatrixForPartialCell_2Unknowns_Case0111(int comp[5], int index[5], double coords[8][2], int i, int j)
 {      
     double cellCenter[2], edgeCenter[2];

     m_pSolver->Add_A(index[1],index[1], 1);
     m_pSolver->Add_A(index[0],index[0], 1);
     m_pSolver->Add_A(index[4],index[4], 1);

     getCellCenter(cellCenter,i,j);    
     m_pSolver->Add_b(index[1], m_pLaplace->getExactSolution(cellCenter, comp[1]));    
     m_pSolver->Add_b(index[0], m_pLaplace->getExactSolution(cellCenter, comp[0]));
     edgeCenter[0] = 0.5*(coords[4][0]+coords[7][0]);
     edgeCenter[1] = 0.5*(coords[4][1]+coords[7][1]);
     m_pSolver->Add_b(index[4], m_pLaplace->getExactSolution(edgeCenter, comp[0]));
 }
 void EBM2D_CARTESIAN::setMatrixForPartialCell_2Unknowns_Case0100(int comp[5], int index[5], double coords[8][2], int i, int j)
 {
     double cellCenter[2], edgeCenter[2];

     m_pSolver->Add_A(index[2],index[2], 1);
     m_pSolver->Add_A(index[1],index[1], 1);
     m_pSolver->Add_A(index[4],index[4], 1);

     getCellCenter(cellCenter,i,j);    
     m_pSolver->Add_b(index[2], m_pLaplace->getExactSolution(cellCenter, comp[2]));
     m_pSolver->Add_b(index[1], m_pLaplace->getExactSolution(cellCenter, comp[1]));
     edgeCenter[0] = 0.5*(coords[5][0]+coords[4][0]);
     edgeCenter[1] = 0.5*(coords[5][1]+coords[4][1]);
     m_pSolver->Add_b(index[4], m_pLaplace->getExactSolution(edgeCenter, comp[0]));
 }
 void EBM2D_CARTESIAN::setMatrixForPartialCell_2Unknowns_Case0010(int comp[5], int index[5], double coords[8][2], int i, int j)
 {
     double cellCenter[2], edgeCenter[2];

     m_pSolver->Add_A(index[3],index[3], 1);
     m_pSolver->Add_A(index[2],index[2], 1);
     m_pSolver->Add_A(index[4],index[4], 1);

     getCellCenter(cellCenter,i,j);
     m_pSolver->Add_b(index[3], m_pLaplace->getExactSolution(cellCenter, comp[3]));
     m_pSolver->Add_b(index[2], m_pLaplace->getExactSolution(cellCenter, comp[2]));
     edgeCenter[0] = 0.5*(coords[6][0]+coords[5][0]);
     edgeCenter[1] = 0.5*(coords[6][1]+coords[5][1]);
     m_pSolver->Add_b(index[4], m_pLaplace->getExactSolution(edgeCenter, comp[0]));
 }

 void EBM2D_CARTESIAN::setMatrixForPartialCell_2Unknowns_Case0001(int comp[5], int index[5], double coords[8][2], int i, int j)
 {
     //printf("EBM2D_CARTESIAN::setMatrixForPartialCell_2Unknowns_Case0001(): "
     //        "coords[0] ={%f,%f}, cell=(%d,%d)\n", coords[0][0],coords[0][1],i,j);    
     /*double cellCenter[2], edgeCenter[2];

     m_pSolver->Add_A(index[0],index[0], 1);
     m_pSolver->Add_A(index[3],index[3], 1);
     m_pSolver->Add_A(index[4],index[4], 1);

     getCellCenter(i,j,cellCenter);
     m_pSolver->Add_b(index[0], m_pLaplace->getExactSolution(cellCenter, comp[0]));
     m_pSolver->Add_b(index[3], m_pLaplace->getExactSolution(cellCenter, comp[3]));
     edgeCenter[0] = 0.5*(coords[7][0]+coords[6][0]);
     edgeCenter[1] = 0.5*(coords[7][1]+coords[6][1]);
     m_pSolver->Add_b(index[4], m_pLaplace->getExactSolution(edgeCenter, comp[0]));*/

     double edgeLength, beta, edgeCenter[2];
     double dx = getCellEdgeLength(0,i,j);
     double dy = getCellEdgeLength(1,i,j);
     double ratio0, ratio1;


     // south    01   
     getCellEdgeCenter(FACE_SOUTH,edgeCenter,i,j);
     beta = m_pLaplace->getBeta(edgeCenter, comp[0]);
     m_pSolver->Add_A(index[0], getCellUnknownIndex(comp[0],i,j-1), beta*dx/dy);
     m_pSolver->Add_A(index[0], index[0], -beta*dx/dy);

     // east     11
     getCellEdgeCenter(FACE_EAST,edgeCenter,i,j);
     beta = m_pLaplace->getBeta(edgeCenter, comp[0]);
     m_pSolver->Add_A(index[0], getCellUnknownIndex(comp[0],i+1,j), beta*dy/dx);
     m_pSolver->Add_A(index[0], index[0], -beta*dy/dx);    

     // north    11
     edgeLength    = coords[2][0]-coords[6][0];
     edgeCenter[0] = coords[2][0] - edgeLength/2;
     edgeCenter[1] = coords[2][1];    
     beta = m_pLaplace->getBeta(edgeCenter, comp[0]);
     ratio0 = (dx/2+edgeLength/2)/dx;
     ratio1 = 1-ratio0;
     m_pSolver->Add_A(index[0], getCellUnknownIndex(comp[0],i,j+1),   ratio0*beta*edgeLength/dy);
     m_pSolver->Add_A(index[0], index[0],                            -ratio0*beta*edgeLength/dy);
     m_pSolver->Add_A(index[0], getCellUnknownIndex(comp[0],i+1,j+1), ratio1*beta*edgeLength/dy);
     m_pSolver->Add_A(index[0], getCellUnknownIndex(comp[0],i+1,j),  -ratio1*beta*edgeLength/dy);

     edgeLength    = dx - edgeLength;
     edgeCenter[0] = coords[6][0] - edgeLength/2;
     beta = m_pLaplace->getBeta(edgeCenter, comp[3]);
     ratio0 = (dx/2+edgeLength/2)/dx;    
     ratio1 = 1-ratio0;
     m_pSolver->Add_A(index[3], getCellUnknownIndex(comp[3],i,j+1),   ratio0*beta*edgeLength/dy);
     m_pSolver->Add_A(index[3], index[3],                            -ratio0*beta*edgeLength/dy);
     m_pSolver->Add_A(index[3], getCellUnknownIndex(comp[3],i-1,j+1), ratio1*beta*edgeLength/dy);
     m_pSolver->Add_A(index[3], getCellUnknownIndex(comp[3],i-1,j),  -ratio1*beta*edgeLength/dy);    

     // west     10
     edgeLength    = coords[7][1]-coords[0][1];
     edgeCenter[0] = coords[3][0];        
     edgeCenter[1] = coords[7][1] - edgeLength/2;    
     beta = m_pLaplace->getBeta(edgeCenter, comp[0]);
     ratio0 = (dy/2+edgeLength/2)/dy;
     ratio1 = 1-ratio0;
     m_pSolver->Add_A(index[0], getCellUnknownIndex(comp[0],i-1,j),   ratio0*beta*edgeLength/dx);
     m_pSolver->Add_A(index[0], index[0],                            -ratio0*beta*edgeLength/dx);
     m_pSolver->Add_A(index[0], getCellUnknownIndex(comp[0],i-1,j-1), ratio1*beta*edgeLength/dx);
     m_pSolver->Add_A(index[0], getCellUnknownIndex(comp[0],i,j-1),  -ratio1*beta*edgeLength/dx);        // wrongly used (i-1,j).

     edgeLength    = dy - edgeLength;
     edgeCenter[1] = coords[3][1] - edgeLength/2;
     beta = m_pLaplace->getBeta(edgeCenter, comp[3]);
     ratio0 = (dy/2+edgeLength/2)/dy;
     ratio1 = 1-ratio0;
     m_pSolver->Add_A(index[3], getCellUnknownIndex(comp[3],i-1,j),   ratio0*beta*edgeLength/dx);
     m_pSolver->Add_A(index[3], index[3],                            -ratio0*beta*edgeLength/dx);
     m_pSolver->Add_A(index[3], getCellUnknownIndex(comp[3],i-1,j+1), ratio1*beta*edgeLength/dx);
     m_pSolver->Add_A(index[3], getCellUnknownIndex(comp[3],i,j+1),  -ratio1*beta*edgeLength/dx);       


     // intfc    
     edgeLength = sqrt( (coords[6][0]-coords[7][0])*(coords[6][0]-coords[7][0])
		       +(coords[6][1]-coords[7][1])*(coords[6][1]-coords[7][1]));
     edgeCenter[0] = 0.5*(coords[6][0]+coords[7][0]);    //wrongly used "0.5/..."
     edgeCenter[1] = 0.5*(coords[6][1]+coords[7][1]);
     double normal[2];
     normal[0] = - (coords[6][1]-coords[7][1])/edgeLength;
     normal[1] =   (coords[6][0]-coords[7][0])/edgeLength;
     int intfcIndex[3], intfcComp[2];
     intfcIndex[0] = index[0];
     intfcIndex[1] = index[3];
     intfcIndex[2] = index[4];
     intfcComp[0] = comp[0];
     intfcComp[1] = comp[3];
     double areaFraction[] = {1,1};
     setMatrixForInterface(intfcIndex,edgeCenter,normal, intfcComp, edgeLength,areaFraction,i,j);

     // two areas
     // the smaller area first
     double area, rhs, cellCenter[2];
     getCellCenter(cellCenter,i,j);

     area = 0.5*(coords[6][0]-coords[3][0])*(coords[3][1]-coords[7][1]);
     rhs = m_pLaplace->getRightHandSide(cellCenter, comp[3]);
     m_pSolver->Add_b(index[3], rhs*area);

     area = dx*dy-area;
     rhs = m_pLaplace->getRightHandSide(cellCenter, comp[0]);
     m_pSolver->Add_b(index[0], rhs*area); 

 }
 void EBM2D_CARTESIAN::setMatrixForPartialCell_2Unknowns_Case0110(int comp[5], int index[5], double coords[8][2], int i, int j)
 {
     double edgeLength, beta, edgeCenter[2];
     double dx = getCellEdgeLength(0,i,j);
     double dy = getCellEdgeLength(1,i,j);
     double ratio0, ratio1;


     // east      
     getCellEdgeCenter(FACE_EAST,edgeCenter,i,j);
     beta = m_pLaplace->getBeta(edgeCenter, comp[1]);
     m_pSolver->Add_A(index[1], getCellUnknownIndex(comp[1],i+1,j), beta*dy/dx);
     m_pSolver->Add_A(index[1], index[1], -beta*dy/dx);

     // west     
     getCellEdgeCenter(FACE_WEST,edgeCenter,i,j);
     beta = m_pLaplace->getBeta(edgeCenter, comp[0]);
     m_pSolver->Add_A(index[0], getCellUnknownIndex(comp[0],i-1,j), beta*dy/dx);
     m_pSolver->Add_A(index[0], index[0], -beta*dy/dx);

     // north    
     //edgeLength    = coords[2][0]-coords[6][0];
     edgeLength    = coords[2][0]-coords[7][0];
     edgeCenter[0] = coords[2][0] - edgeLength/2;
     edgeCenter[1] = coords[2][1];    
     beta = m_pLaplace->getBeta(edgeCenter, comp[1]);
     ratio0 = (dx/2+edgeLength/2)/dx;
     ratio1 = 1-ratio0;
     m_pSolver->Add_A(index[1], getCellUnknownIndex(comp[1],i,j+1),   ratio0*beta*edgeLength/dy);
     m_pSolver->Add_A(index[1], index[1],                            -ratio0*beta*edgeLength/dy);
     m_pSolver->Add_A(index[1], getCellUnknownIndex(comp[1],i+1,j+1), ratio1*beta*edgeLength/dy);
     m_pSolver->Add_A(index[1], getCellUnknownIndex(comp[1],i+1,j),  -ratio1*beta*edgeLength/dy);

     edgeLength    = dx - edgeLength;
     //edgeCenter[0] = coords[6][0] - edgeLength/2;
     edgeCenter[0] = coords[7][0] - edgeLength/2;
     beta = m_pLaplace->getBeta(edgeCenter, comp[0]);
     ratio0 = (dx/2+edgeLength/2)/dx;    
     ratio1 = 1-ratio0;
     m_pSolver->Add_A(index[0], getCellUnknownIndex(comp[0],i,j+1),   ratio0*beta*edgeLength/dy);
     m_pSolver->Add_A(index[0], index[0],                            -ratio0*beta*edgeLength/dy);
     m_pSolver->Add_A(index[0], getCellUnknownIndex(comp[0],i-1,j+1), ratio1*beta*edgeLength/dy);
     m_pSolver->Add_A(index[0], getCellUnknownIndex(comp[0],i-1,j),  -ratio1*beta*edgeLength/dy);    

     // south     
     //edgeLength    = coords[1][0]-coords[4][0];
     edgeLength    = coords[1][0]-coords[6][0];
     edgeCenter[0] = coords[1][0] - edgeLength/2;
     edgeCenter[1] = coords[1][1];    
     beta = m_pLaplace->getBeta(edgeCenter, comp[1]);
     ratio0 = (dx/2+edgeLength/2)/dx;
     ratio1 = 1-ratio0;
     m_pSolver->Add_A(index[1], getCellUnknownIndex(comp[1],i,j-1),   ratio0*beta*edgeLength/dy);
     m_pSolver->Add_A(index[1], index[1],                            -ratio0*beta*edgeLength/dy);
     m_pSolver->Add_A(index[1], getCellUnknownIndex(comp[1],i+1,j-1), ratio1*beta*edgeLength/dy);
     m_pSolver->Add_A(index[1], getCellUnknownIndex(comp[1],i+1,j),  -ratio1*beta*edgeLength/dy);

     edgeLength    = dx - edgeLength;
     //edgeCenter[0] = coords[4][0] - edgeLength/2;
     edgeCenter[0] = coords[6][0] - edgeLength/2;
     beta = m_pLaplace->getBeta(edgeCenter, comp[0]);
     ratio0 = (dx/2+edgeLength/2)/dx;    
     ratio1 = 1-ratio0;
     m_pSolver->Add_A(index[0], getCellUnknownIndex(comp[0],i,j-1),   ratio0*beta*edgeLength/dy);
     m_pSolver->Add_A(index[0], index[0],                            -ratio0*beta*edgeLength/dy);
     m_pSolver->Add_A(index[0], getCellUnknownIndex(comp[0],i-1,j-1), ratio1*beta*edgeLength/dy);
     m_pSolver->Add_A(index[0], getCellUnknownIndex(comp[0],i-1,j),  -ratio1*beta*edgeLength/dy);      


     // intfc    
     /*edgeLength = sqrt( (coords[6][0]-coords[4][0])*(coords[6][0]-coords[4][0])
		       +(coords[6][1]-coords[4][1])*(coords[6][1]-coords[4][1]));
     edgeCenter[0] = 0.5*(coords[6][0]+coords[4][0]);    //wrongly used "0.5/..."
     edgeCenter[1] = 0.5*(coords[6][1]+coords[4][1]);*/
     edgeLength = sqrt( (coords[6][0]-coords[7][0])*(coords[6][0]-coords[7][0])
		       +(coords[6][1]-coords[7][1])*(coords[6][1]-coords[7][1]));
     edgeCenter[0] = 0.5*(coords[6][0]+coords[7][0]);    //wrongly used "0.5/..."
     edgeCenter[1] = 0.5*(coords[6][1]+coords[7][1]);
     double normal[2];
     /*normal[0] = - (coords[6][1]-coords[4][1])/edgeLength;
     normal[1] =   (coords[6][0]-coords[4][0])/edgeLength;*/
     normal[0] = - (coords[7][1]-coords[6][1])/edgeLength;
     normal[1] =   (coords[7][0]-coords[6][0])/edgeLength;
     int intfcIndex[3], intfcComp[2];
     intfcIndex[0] = index[1];
     intfcIndex[1] = index[0];
     intfcIndex[2] = index[4];
     intfcComp[0] = comp[1];
     intfcComp[1] = comp[0];
     
     double areaFraction[] = {1,1};
     setMatrixForInterface(intfcIndex,edgeCenter,normal, intfcComp, edgeLength,areaFraction,i,j);

     // two areas
     // the left area first
     double area, rhs, cellCenter[2];
     getCellCenter(cellCenter,i,j);

     //area = 0.5*((coords[6][0]-coords[3][0])+(coords[4][0]-coords[0][0]))*dy;
     area = 0.5*((coords[7][0]-coords[3][0])+(coords[6][0]-coords[0][0]))*dy;
     rhs = m_pLaplace->getRightHandSide(cellCenter, comp[0]);
     m_pSolver->Add_b(index[0], rhs*area);

     area = dx*dy-area;
     rhs = m_pLaplace->getRightHandSide(cellCenter, comp[1]);
     m_pSolver->Add_b(index[1], rhs*area); 



 }
 void EBM2D_CARTESIAN::setMatrixForPartialCell_2Unknowns_Case0011(int comp[5], int index[5], double coords[8][2], int i, int j)
 {
     /*double cellCenter[2], edgeCenter[2];

     m_pSolver->Add_A(index[2],index[2], 1);
     m_pSolver->Add_A(index[1],index[1], 1);
     m_pSolver->Add_A(index[4],index[4], 1);

     getCellCenter(i,j,cellCenter);
     m_pSolver->Add_b(index[2], m_pLaplace->getExactSolution(cellCenter, comp[2]));
     m_pSolver->Add_b(index[1], m_pLaplace->getExactSolution(cellCenter, comp[1]));
     //edgeCenter[0] = 0.5*(coords[5][0]+coords[7][0]);
     //edgeCenter[1] = 0.5*(coords[5][1]+coords[7][1]);
      edgeCenter[0] = 0.5*(coords[5][0]+coords[4][0]);
      edgeCenter[1] = 0.5*(coords[5][1]+coords[4][1]);
     m_pSolver->Add_b(index[4], m_pLaplace->getExactSolution(edgeCenter, comp[0]));*/

     double edgeLength, beta, edgeCenter[2];
     double dx = getCellEdgeLength(0,i,j);
     double dy = getCellEdgeLength(1,i,j);
     double ratio0, ratio1;


     // south    
     getCellEdgeCenter(FACE_SOUTH,edgeCenter,i,j);
     beta = m_pLaplace->getBeta(edgeCenter, comp[0]);
     m_pSolver->Add_A(index[0], getCellUnknownIndex(comp[0],i,j-1), beta*dx/dy);
     m_pSolver->Add_A(index[0], index[0], -beta*dx/dy);

     // north     
     getCellEdgeCenter(FACE_NORTH,edgeCenter,i,j);
     beta = m_pLaplace->getBeta(edgeCenter, comp[3]);
     m_pSolver->Add_A(index[3], getCellUnknownIndex(comp[3],i,j+1), beta*dx/dy);
     m_pSolver->Add_A(index[3], index[3], -beta*dx/dy);

     // east    
     edgeLength    = coords[5][1]-coords[0][1];
     edgeCenter[0] = coords[2][0];        
     edgeCenter[1] = coords[5][1] - edgeLength/2;    
     beta = m_pLaplace->getBeta(edgeCenter, comp[0]);
     ratio0 = (dy/2+edgeLength/2)/dy;
     ratio1 = 1-ratio0;
     m_pSolver->Add_A(index[0], getCellUnknownIndex(comp[0],i+1,j),   ratio0*beta*edgeLength/dx);
     m_pSolver->Add_A(index[0], index[0],                            -ratio0*beta*edgeLength/dx);
     m_pSolver->Add_A(index[0], getCellUnknownIndex(comp[0],i+1,j-1), ratio1*beta*edgeLength/dx);
     m_pSolver->Add_A(index[0], getCellUnknownIndex(comp[0],i,j-1),  -ratio1*beta*edgeLength/dx);        // wrongly used (i-1,j).

     edgeLength    = dy - edgeLength;
     edgeCenter[1] = coords[2][1] - edgeLength/2;
     beta = m_pLaplace->getBeta(edgeCenter, comp[3]);
     ratio0 = (dy/2+edgeLength/2)/dy;
     ratio1 = 1-ratio0;
     m_pSolver->Add_A(index[3], getCellUnknownIndex(comp[3],i+1,j),   ratio0*beta*edgeLength/dx);
     m_pSolver->Add_A(index[3], index[3],                            -ratio0*beta*edgeLength/dx);
     m_pSolver->Add_A(index[3], getCellUnknownIndex(comp[3],i+1,j+1), ratio1*beta*edgeLength/dx);
     m_pSolver->Add_A(index[3], getCellUnknownIndex(comp[3],i,j+1),  -ratio1*beta*edgeLength/dx);      

     // west     
     //edgeLength    = coords[7][1]-coords[0][1];
     edgeLength    = coords[4][1]-coords[0][1];
     edgeCenter[0] = coords[3][0];        
     //edgeCenter[1] = coords[7][1] - edgeLength/2;    
     edgeCenter[1] = coords[4][1] - edgeLength/2;    
     beta = m_pLaplace->getBeta(edgeCenter, comp[0]);
     ratio0 = (dy/2+edgeLength/2)/dy;
     ratio1 = 1-ratio0;
     m_pSolver->Add_A(index[0], getCellUnknownIndex(comp[0],i-1,j),   ratio0*beta*edgeLength/dx);
     m_pSolver->Add_A(index[0], index[0],                            -ratio0*beta*edgeLength/dx);
     m_pSolver->Add_A(index[0], getCellUnknownIndex(comp[0],i-1,j-1), ratio1*beta*edgeLength/dx);
     m_pSolver->Add_A(index[0], getCellUnknownIndex(comp[0],i,j-1),  -ratio1*beta*edgeLength/dx);        // wrongly used (i-1,j).

     edgeLength    = dy - edgeLength;
     edgeCenter[1] = coords[3][1] - edgeLength/2;
     beta = m_pLaplace->getBeta(edgeCenter, comp[3]);
     ratio0 = (dy/2+edgeLength/2)/dy;
     ratio1 = 1-ratio0;
     m_pSolver->Add_A(index[3], getCellUnknownIndex(comp[3],i-1,j),   ratio0*beta*edgeLength/dx);
     m_pSolver->Add_A(index[3], index[3],                            -ratio0*beta*edgeLength/dx);
     m_pSolver->Add_A(index[3], getCellUnknownIndex(comp[3],i-1,j+1), ratio1*beta*edgeLength/dx);
     m_pSolver->Add_A(index[3], getCellUnknownIndex(comp[3],i,j+1),  -ratio1*beta*edgeLength/dx);       


     // intfc    
     /*edgeLength = sqrt( (coords[5][0]-coords[7][0])*(coords[5][0]-coords[7][0])
		       +(coords[5][1]-coords[7][1])*(coords[5][1]-coords[7][1]));
     edgeCenter[0] = 0.5*(coords[5][0]+coords[7][0]);    //wrongly used "0.5/..."
     edgeCenter[1] = 0.5*(coords[5][1]+coords[7][1]);*/
     edgeLength = sqrt( (coords[5][0]-coords[4][0])*(coords[5][0]-coords[4][0])
		       +(coords[5][1]-coords[4][1])*(coords[5][1]-coords[4][1]));
     edgeCenter[0] = 0.5*(coords[5][0]+coords[4][0]);    //wrongly used "0.5/..."
     edgeCenter[1] = 0.5*(coords[5][1]+coords[4][1]);
     double normal[2];
     /*normal[0] = - (coords[5][1]-coords[7][1])/edgeLength;
       normal[1] =   (coords[5][0]-coords[7][0])/edgeLength;*/
     normal[0] = - (coords[5][1]-coords[4][1])/edgeLength;
     normal[1] =   (coords[5][0]-coords[4][0])/edgeLength;
     int intfcIndex[3], intfcComp[2];
     intfcIndex[0] = index[0];
     intfcIndex[1] = index[3];
     intfcIndex[2] = index[4];
     intfcComp[0] = comp[0];
     intfcComp[1] = comp[3];
     
     double areaFraction[] = {1,1};
     setMatrixForInterface(intfcIndex,edgeCenter,normal, intfcComp, edgeLength,areaFraction,i,j);

     // two areas
     // the smaller area first
     double area, rhs, cellCenter[2];
     getCellCenter(cellCenter,i,j);

     //area = 0.5*((coords[3][1]-coords[7][1])+(coords[2][1]-coords[5][1]))*dx;
     area = 0.5*((coords[3][1]-coords[4][1])+(coords[2][1]-coords[5][1]))*dx;
     rhs = m_pLaplace->getRightHandSide(cellCenter, comp[3]);
     m_pSolver->Add_b(index[3], rhs*area);

     area = dx*dy-area;
     rhs = m_pLaplace->getRightHandSide(cellCenter, comp[0]);
     m_pSolver->Add_b(index[0], rhs*area);  
 }


 void EBM2D_CARTESIAN::setMatrixForPartialCell_2Unknowns_byRotation(int rotation, int comp[5], int index[5], double coords[8][2], int i, int j)
 {   
     // get the rotation indices.
     int IIndex[3][3], JIndex[3][3], CIndex[8], DIndex[2], sign = 1;
     getRotatedCellsIndex(rotation,IIndex,JIndex,i,j);
     getRotatedCoordsIndex(rotation,CIndex);    
     getRotatedDirection(rotation,DIndex);

     double edgeLength, beta, edgeCenter[2];
     double dx = getCellEdgeLength(DIndex[0],i,j);
     double dy = getCellEdgeLength(DIndex[1],i,j);   
     double ratio0, ratio1;    

     // south    01   
     getCellEdgeCenter(CIndex[0],edgeCenter,i,j);
     beta = m_pLaplace->getBeta(edgeCenter, comp[CIndex[0]]);
     m_pSolver->Add_A(index[CIndex[0]], getCellUnknownIndex(comp[CIndex[0]], /*i,j-1*/IIndex[1][0], JIndex[1][0]), beta*dx/dy);
     m_pSolver->Add_A(index[CIndex[0]], index[CIndex[0]], -beta*dx/dy);

     // east     11
     getCellEdgeCenter(CIndex[1],edgeCenter,i,j);
     beta = m_pLaplace->getBeta(edgeCenter, comp[CIndex[0]]);
     m_pSolver->Add_A(index[CIndex[0]], getCellUnknownIndex(comp[CIndex[0]], /*i+1,j*/IIndex[2][1],JIndex[2][1]), beta*dy/dx);
     m_pSolver->Add_A(index[CIndex[0]], index[CIndex[0]], -beta*dy/dx);    

     // north    11
     //edgeLength    = coords[CIndex[2]][DIndex[0]]-coords[CIndex[6]][DIndex[0]];
     edgeLength    = coords[CIndex[2]][DIndex[0]]-coords[CIndex[7]][DIndex[0]];
     sign = edgeLength>0 ? 1: -1;
     edgeLength *= sign;
     edgeCenter[DIndex[0]] = coords[CIndex[2]][DIndex[0]] - sign*edgeLength/2;
     edgeCenter[DIndex[1]] = coords[CIndex[2]][DIndex[1]];    
     beta = m_pLaplace->getBeta(edgeCenter, comp[CIndex[0]]);
     ratio0 = (dx/2+edgeLength/2)/dx;
     ratio1 = 1-ratio0;
     m_pSolver->Add_A(index[CIndex[0]], getCellUnknownIndex(comp[CIndex[0]], /*i,j+1*/IIndex[1][2],JIndex[1][2]),   ratio0*beta*edgeLength/dy);
     m_pSolver->Add_A(index[CIndex[0]], index[CIndex[0]],                            -ratio0*beta*edgeLength/dy);
     m_pSolver->Add_A(index[CIndex[0]], getCellUnknownIndex(comp[CIndex[0]], /*i+1,j+1*/IIndex[2][2],JIndex[2][2]), ratio1*beta*edgeLength/dy);
     m_pSolver->Add_A(index[CIndex[0]], getCellUnknownIndex(comp[CIndex[0]], /*i+1,j*/IIndex[2][1],JIndex[2][1]),  -ratio1*beta*edgeLength/dy);

     edgeLength    = dx - edgeLength;
     //edgeCenter[DIndex[0]] = coords[CIndex[6]][DIndex[0]] - sign*edgeLength/2;
     edgeCenter[DIndex[0]] = coords[CIndex[7]][DIndex[0]] - sign*edgeLength/2;
     beta = m_pLaplace->getBeta(edgeCenter, comp[CIndex[3]]);
     ratio0 = (dx/2+edgeLength/2)/dx;    
     ratio1 = 1-ratio0;
     m_pSolver->Add_A(index[CIndex[3]], getCellUnknownIndex(comp[CIndex[3]], /*i,j+1*/IIndex[1][2],JIndex[1][2]),   ratio0*beta*edgeLength/dy);
     m_pSolver->Add_A(index[CIndex[3]], index[CIndex[3]],                            -ratio0*beta*edgeLength/dy);
     m_pSolver->Add_A(index[CIndex[3]], getCellUnknownIndex(comp[CIndex[3]], /*i-1,j+1*/IIndex[0][2],JIndex[0][2]), ratio1*beta*edgeLength/dy);
     m_pSolver->Add_A(index[CIndex[3]], getCellUnknownIndex(comp[CIndex[3]], /*i-1,j*/IIndex[0][1],JIndex[0][1]),  -ratio1*beta*edgeLength/dy);    

     // west     10
     //edgeLength    = coords[CIndex[7]][DIndex[1]]-coords[CIndex[0]][DIndex[1]];
     edgeLength    = coords[CIndex[4]][DIndex[1]]-coords[CIndex[0]][DIndex[1]];
     sign = edgeLength>0 ? 1 : -1;
     edgeLength *= sign;
     edgeCenter[DIndex[0]] = coords[CIndex[3]][DIndex[0]];        
     //edgeCenter[DIndex[1]] = coords[CIndex[7]][DIndex[1]] - sign*edgeLength/2;    
     edgeCenter[DIndex[1]] = coords[CIndex[4]][DIndex[1]] - sign*edgeLength/2;    
     beta = m_pLaplace->getBeta(edgeCenter, comp[CIndex[0]]);
     ratio0 = (dy/2+edgeLength/2)/dy;
     ratio1 = 1-ratio0;
     m_pSolver->Add_A(index[CIndex[0]], getCellUnknownIndex(comp[CIndex[0]], /*i-1,j*/IIndex[0][1],JIndex[0][1]),   ratio0*beta*edgeLength/dx);
     m_pSolver->Add_A(index[CIndex[0]], index[CIndex[0]],                            -ratio0*beta*edgeLength/dx);
     m_pSolver->Add_A(index[CIndex[0]], getCellUnknownIndex(comp[CIndex[0]], /*i-1,j-1*/IIndex[0][0],JIndex[0][0]), ratio1*beta*edgeLength/dx);
     m_pSolver->Add_A(index[CIndex[0]], getCellUnknownIndex(comp[CIndex[0]], /*i,j-1*/IIndex[1][0], JIndex[1][0]),  -ratio1*beta*edgeLength/dx);        // wrongly used (i-1,j).

     edgeLength    = dy - edgeLength;
     edgeCenter[DIndex[1]] = coords[CIndex[3]][DIndex[1]] - sign*edgeLength/2;
     beta = m_pLaplace->getBeta(edgeCenter, comp[CIndex[3]]);
     ratio0 = (dy/2+edgeLength/2)/dy;
     ratio1 = 1-ratio0;
     m_pSolver->Add_A(index[CIndex[3]], getCellUnknownIndex(comp[CIndex[3]], /*i-1,j*/IIndex[0][1],JIndex[0][1]),   ratio0*beta*edgeLength/dx);
     m_pSolver->Add_A(index[CIndex[3]], index[CIndex[3]],                            -ratio0*beta*edgeLength/dx);
     m_pSolver->Add_A(index[CIndex[3]], getCellUnknownIndex(comp[CIndex[3]], /*i-1,j+1*/IIndex[0][2],JIndex[0][2]), ratio1*beta*edgeLength/dx);
     m_pSolver->Add_A(index[CIndex[3]], getCellUnknownIndex(comp[CIndex[3]], /*i,j+1*/IIndex[1][2],JIndex[1][2]),  -ratio1*beta*edgeLength/dx);       


     // intfc    
     //edgeLength = sqrt( (coords[CIndex[6]][0]-coords[CIndex[7]][0])*(coords[CIndex[6]][0]-coords[CIndex[7]][0])
     //	       +(coords[CIndex[6]][1]-coords[CIndex[7]][1])*(coords[CIndex[6]][1]-coords[CIndex[7]][1]));
     edgeLength = sqrt( (coords[CIndex[4]][0]-coords[CIndex[7]][0])*(coords[CIndex[4]][0]-coords[CIndex[7]][0])
		   +(coords[CIndex[4]][1]-coords[CIndex[7]][1])*(coords[CIndex[4]][1]-coords[CIndex[7]][1]));
     //edgeCenter[0] = 0.5*(coords[CIndex[6]][0]+coords[CIndex[7]][0]);    //wrongly used "0.5/..."
     //edgeCenter[1] = 0.5*(coords[CIndex[6]][1]+coords[CIndex[7]][1]);    // no need to use DIndex for the above 4 lines.
     edgeCenter[0] = 0.5*(coords[CIndex[4]][0]+coords[CIndex[7]][0]);    //wrongly used "0.5/..."
     edgeCenter[1] = 0.5*(coords[CIndex[4]][1]+coords[CIndex[7]][1]);    // no need to use DIndex for the above 4 lines.
     double normal[2],tmp;    
     //normal[DIndex[0]] = (coords[CIndex[6]][DIndex[0]]-coords[CIndex[7]][DIndex[0]])/edgeLength;
     //normal[DIndex[1]] = (coords[CIndex[6]][DIndex[1]]-coords[CIndex[7]][DIndex[1]])/edgeLength;
     normal[DIndex[0]] = (coords[CIndex[7]][DIndex[0]]-coords[CIndex[4]][DIndex[0]])/edgeLength;
     normal[DIndex[1]] = (coords[CIndex[7]][DIndex[1]]-coords[CIndex[4]][DIndex[1]])/edgeLength;
     tmp = normal[0];
     normal[0] = -normal[1];
     normal[1] = tmp;

     int intfcIndex[3], intfcComp[2];
     intfcIndex[0] = index[CIndex[0]];
     intfcIndex[1] = index[CIndex[3]];
     intfcIndex[2] = index[4];   //wrongly used index[CIndex[4]];
     intfcComp[0] = comp[CIndex[0]];
     intfcComp[1] = comp[CIndex[3]];
     
     double areaFraction[] = {1,1};
     setMatrixForInterface(intfcIndex,edgeCenter,normal, intfcComp, edgeLength,areaFraction,i,j);

     // two areas
     // the smaller area first
     double area, rhs, cellCenter[2];
     getCellCenter(cellCenter,i,j);

     //area = 0.5*(coords[CIndex[6]][DIndex[0]]-coords[CIndex[3]][DIndex[0]])
     //       *(coords[CIndex[3]][DIndex[1]]-coords[CIndex[7]][DIndex[1]]);
     area = 0.5*(coords[CIndex[7]][DIndex[0]]-coords[CIndex[3]][DIndex[0]])
	       *(coords[CIndex[3]][DIndex[1]]-coords[CIndex[4]][DIndex[1]]);
     rhs = m_pLaplace->getRightHandSide(cellCenter, comp[CIndex[3]]);
     m_pSolver->Add_b(index[CIndex[3]], rhs*area);

     area = dx*dy-area;
     rhs = m_pLaplace->getRightHandSide(cellCenter, comp[CIndex[0]]);
     m_pSolver->Add_b(index[CIndex[0]], rhs*area); 
 }

/*
  * This function generates two index matrix for i and j as the following:
  * IIndex[][]
  *      i-1 i i+1
  *      i-1 i i+1
  *      i-1 i i+1
  * JIndex[][]
  *      j+1 j+1 j+1
  *      j   j   j
  *      j-1 j-1 j-1
  * Then this function rotate the two matrix clockwise by rotation*pi/2 around
  * the middle of the two matrix: IIndex[1][1] & JIndex[1][1].
  * For example, if rotation==1, we have
  * IIndex[][]
  *      i-1 i-1 i-1
  *      i   i   i
  *      i+1 i+1 i+1.
  * Similarly for JIndex[][].   
  */ 
 void EBM2D_CARTESIAN::getRotatedCellsIndex(int rotation, int IIndex[3][3], int JIndex[3][3], int i, int j)
 {
     switch(rotation%4)
     {
	 case 0:
	     IIndex[0][2] = i-1;     IIndex[1][2] = i;       IIndex[2][2] = i+1;
	     IIndex[0][1] = i-1;     IIndex[1][1] = i;       IIndex[2][1] = i+1;
	     IIndex[0][0] = i-1;     IIndex[1][0] = i;       IIndex[2][0] = i+1;

	     JIndex[0][2] = j+1;     JIndex[1][2] = j+1;     JIndex[2][2] = j+1;
	     JIndex[0][1] = j;       JIndex[1][1] = j;       JIndex[2][1] = j;
	     JIndex[0][0] = j-1;     JIndex[1][0] = j-1;     JIndex[2][0] = j-1;            
	     break;
	 case 1:
	     IIndex[0][2] = i-1;     IIndex[1][2] = i-1;     IIndex[2][2] = i-1;
	     IIndex[0][1] = i;       IIndex[1][1] = i;       IIndex[2][1] = i;
	     IIndex[0][0] = i+1;     IIndex[1][0] = i+1;     IIndex[2][0] = i+1;

	     JIndex[0][2] = j-1;     JIndex[1][2] = j;       JIndex[2][2] = j+1;
	     JIndex[0][1] = j-1;     JIndex[1][1] = j;       JIndex[2][1] = j+1;
	     JIndex[0][0] = j-1;     JIndex[1][0] = j;       JIndex[2][0] = j+1;  
	     break;
	 case 2:
	     IIndex[0][2] = i+1;     IIndex[1][2] = i;       IIndex[2][2] = i-1;
	     IIndex[0][1] = i+1;     IIndex[1][1] = i;       IIndex[2][1] = i-1;
	     IIndex[0][0] = i+1;     IIndex[1][0] = i;       IIndex[2][0] = i-1;

	     JIndex[0][2] = j-1;     JIndex[1][2] = j-1;     JIndex[2][2] = j-1;
	     JIndex[0][1] = j;       JIndex[1][1] = j;       JIndex[2][1] = j;
	     JIndex[0][0] = j+1;     JIndex[1][0] = j+1;     JIndex[2][0] = j+1;
	     break;
	 case 3:
	     IIndex[0][2] = i+1;     IIndex[1][2] = i+1;     IIndex[2][2] = i+1;
	     IIndex[0][1] = i;       IIndex[1][1] = i;       IIndex[2][1] = i;
	     IIndex[0][0] = i-1;     IIndex[1][0] = i-1;     IIndex[2][0] = i-1;

	     JIndex[0][2] = j+1;     JIndex[1][2] = j;       JIndex[2][2] = j-1;
	     JIndex[0][1] = j+1;     JIndex[1][1] = j;       JIndex[2][1] = j-1;
	     JIndex[0][0] = j+1;     JIndex[1][0] = j;       JIndex[2][0] = j-1; 
	     break;
     }
 }
 /*
  * The function generates index[] according to the following ordering
  *      362
  *      7 5
  *      041
  * index[] = {0,1,2,3,4,5,6,7}
  * Then it rotates the index[] clockwise by rotation*pi/4. For example, if 
  * rotation==1, then we have
  *      073
  *      4 6
  *      152
  * index[] = {1,2,3,0,5,6,7,4}
  *
  * the above ordering of index has been changed to the following:
  *      372
  *      4 5
  *      061
  */
 void EBM2D_CARTESIAN::getRotatedCoordsIndex(int rotation, int CIndex[8])
 {
     switch(rotation%4)
     {
	 case 0:
	     CIndex[0] = 0;   CIndex[1] = 1;   CIndex[2] = 2;   CIndex[3] = 3;
	     //CIndex[4] = 4;   CIndex[5] = 5;   CIndex[6] = 6;   CIndex[7] = 7;            
	     CIndex[4] = 4;   CIndex[5] = 5;   CIndex[6] = 6;   CIndex[7] = 7;            
	     break;
	 case 1:
	     CIndex[0] = 1;   CIndex[1] = 2;   CIndex[2] = 3;   CIndex[3] = 0;
	     //CIndex[4] = 5;   CIndex[5] = 6;   CIndex[6] = 7;   CIndex[7] = 4;
	     CIndex[4] = 6;   CIndex[5] = 7;   CIndex[6] = 5;   CIndex[7] = 4;
	     break;
	 case 2:
	     CIndex[0] = 2;   CIndex[1] = 3;   CIndex[2] = 0;   CIndex[3] = 1;
	     //CIndex[4] = 6;   CIndex[5] = 7;   CIndex[6] = 4;   CIndex[7] = 5;
	     CIndex[4] = 5;   CIndex[5] = 4;   CIndex[6] = 7;   CIndex[7] = 6;
	     break;
	 case 3:
	     CIndex[0] = 3;   CIndex[1] = 0;   CIndex[2] = 1;   CIndex[3] = 2;
	     //CIndex[4] = 7;   CIndex[5] = 4;   CIndex[6] = 5;   CIndex[7] = 6;
	     CIndex[4] = 7;   CIndex[5] = 6;   CIndex[6] = 4;   CIndex[7] = 5;
	     break;
     }
 }

 /*
  * Get the coordinates direction in DIndex[2].
  */
 void EBM2D_CARTESIAN::getRotatedDirection(int rotation, int DIndex[2])
 {
     switch(rotation)
     {
	 case 0:
	 case 2:
	     DIndex[0] = 0;
	     DIndex[1] = 1;
	     break;
	 case 1:        
	 case 3:
	     DIndex[0] = 1;
	     DIndex[1] = 0;
	     break;
     }
 }



 /*
  * index0 for the unknown in the negative side
  * index1 for the unknown in the positive side
  * index2 for the unknown on the interface
  * the normal points from the negative side to the positive side
  * comp[0] for negative component, comp[1] for positive component.
  */
 void EBM2D_CARTESIAN::setMatrixForInterface(int intfcIndex[3],double intfcCenter[2], 
					     double intfcNormal[2], int intfcComp[2], 
					     double intfcLength, double areaFraction[2], 
					     int i, int j)
 {
     std::vector<int> indices[2];
     std::vector<double> coefficients[2];
     std::vector<double> distances[2];
     double beta[2];

     // for component 1
     beta[1] = m_pLaplace->getBeta(intfcCenter, intfcComp[1]);
     getSolutionInterpStencil(intfcCenter, intfcNormal, intfcComp[1], indices[1], 
			      coefficients[1], distances[1],i,j);    

     // for component 0;    
     intfcNormal[0] = -intfcNormal[0];     intfcNormal[1] = -intfcNormal[1];
     beta[0] = m_pLaplace->getBeta(intfcCenter, intfcComp[0]);
     getSolutionInterpStencil(intfcCenter, intfcNormal, intfcComp[0], indices[0],
			      coefficients[0], distances[0],i,j);

     double coords[3], coeffs[2][3];
     int k, kk;

     // component 0 flux
     coords[0] = 0;
     coords[1] = distances[0][0];
     coords[2] = distances[0][1];
     get3PointsDerivativeInterpCoeffs(coords, 0, coeffs[0]);    

     m_pSolver->Add_A(intfcIndex[0], intfcIndex[2], 
		      - beta[0]*intfcLength*coeffs[0][0] * areaFraction[0]);
     for(k=0; k<2; k++)    
	 for(kk=0; kk<3; kk++)
	     m_pSolver->Add_A(intfcIndex[0], indices[0][k*3+kk], 
			      - beta[0]*intfcLength*coeffs[0][k+1]  
			      * coefficients[0][k*3+kk] * areaFraction[0]);    

     // component 1 flux
     coords[0] = 0;
     coords[1] = distances[1][0];
     coords[2] = distances[1][1];
     get3PointsDerivativeInterpCoeffs(coords, 0, coeffs[1]);    

     m_pSolver->Add_A(intfcIndex[1], intfcIndex[2], 
		      - beta[1]*intfcLength*coeffs[1][0] * areaFraction[1]);
     for(k=0; k<2; k++)    
	 for(kk=0; kk<3; kk++)
	     m_pSolver->Add_A(intfcIndex[1], indices[1][k*3+kk], 
			      - beta[1]*intfcLength*coeffs[1][k+1]
			      * coefficients[1][k*3+kk] * areaFraction[1]);
     
     // interface equation
     double fraction = std::min(areaFraction[0],areaFraction[1]);

     m_pSolver->Add_A(intfcIndex[2], intfcIndex[2], 
		      + beta[0]*intfcLength*coeffs[0][0] * fraction);
     for(k=0; k<2; k++)    
	 for(kk=0; kk<3; kk++)
	     m_pSolver->Add_A(intfcIndex[2], indices[0][k*3+kk], 
			      + beta[0]*intfcLength*coeffs[0][k+1]
			      * coefficients[0][k*3+kk] * fraction);  

     m_pSolver->Add_A(intfcIndex[2], intfcIndex[2], 
		      + beta[1]*intfcLength*coeffs[1][0] * fraction);
     for(k=0; k<2; k++)    
	 for(kk=0; kk<3; kk++)
	     m_pSolver->Add_A(intfcIndex[2], indices[1][k*3+kk], 
			      + beta[1]*intfcLength*coeffs[1][k+1]
			      * coefficients[1][k*3+kk] * fraction);   
 }

/*
 * Deal with Dirichlet boundary conditions.
 * Modified from EBM2D_CARTESIAN::setMatrixForInterface().
 */
void EBM2D_CARTESIAN::setMatrixForBoundary_Dirichlet(int intfcIndex[2], double intfcCenter[2], double intfcNormal[2], 
			      int intfcComp[2], double intfcLength, double areaFraction[2], int i, int j)
{
     if(intfcIndex[0]>=0||intfcIndex[1]<0)
	  printf("EBM2D_CARTESIAN::setMatrixForBoundary_Dirichlet: logic error!\n");     

     std::vector<int> indices;
     std::vector<double> coefficients;
     std::vector<double> distances;
     double beta;

     // for component 1
     beta = m_pLaplace->getBeta(intfcCenter, intfcComp[1]);
     getSolutionInterpStencil(intfcCenter, intfcNormal, intfcComp[1], indices, coefficients, distances,i,j);    

     double coords[3], coeffs[3];
     int k, kk;

     // component 1 flux
     coords[0] = 0;
     coords[1] = distances[0];
     coords[2] = distances[1];
     get3PointsDerivativeInterpCoeffs(coords, 0, coeffs);    

     double potential = m_pLaplace->getExactSolution(intfcCenter, intfcComp[1]);

     //m_pSolver->Add_A(intfcIndex[1], intfcIndex[2], - beta*intfcLength*coeffs[0] * areaFraction[1]);
     m_pSolver->Add_b(intfcIndex[1],  beta*intfcLength*coeffs[0] * areaFraction[1] * potential);
     for(k=0; k<2; k++)    
	 for(kk=0; kk<3; kk++)
	     m_pSolver->Add_A(intfcIndex[1], indices[k*3+kk], - beta*intfcLength*coeffs[k+1]*coefficients[k*3+kk] * areaFraction[1]);

}


 /*
  * Deal with Neumann boundary condition.
  */
 void EBM2D_CARTESIAN::setMatrixForBoundary_Neumann(int intfcIndex[2], double intfcCenter[2], double normal[2], 
					    int intfcComp[2], double intfcArea, double areaFraction[2], int i, int j)
 {
      double flux, value, beta;   
     if(intfcIndex[0]>=0)
     {
	  printf("EBM2D_CARTESIAN::setMatrixForBoundary: logic error!\n");
	  flux = m_pLaplace->getFlux(intfcCenter, intfcComp[0], normal);    
	  beta = m_pLaplace->getBeta(intfcCenter, intfcComp[0]);
	  m_pSolver->Add_b(intfcIndex[0], -flux*beta*intfcArea * areaFraction[0]);    // negative because it is moved to the rhs.
     }
     if(intfcIndex[1]>=0)
     {
	  flux = - m_pLaplace->getFlux(intfcCenter, intfcComp[1], normal); // negative because it has negative normal
	  beta = m_pLaplace->getBeta(intfcCenter, intfcComp[1]);
	  m_pSolver->Add_b(intfcIndex[1], -flux*beta*intfcArea * areaFraction[1]);    // negative because it is moved to the rhs.
	  //printf("EBM2D_CARTESIAN::setMatrixForBoundary: %d, %f\n", intfcIndex[1], -flux*beta*intfcArea);
	  //printf("\t intfcindex[1=%d, -flux*beta*intfcArea = %e\n", intfcIndex[1], -flux*beta*intfcArea);
     }
     //value = m_pLaplace->getExactSolution(intfcCenter, intfcComp[0]);  // assuming the solution is continuous at the interface.
     //m_pSolver->Add_A(intfcIndex[2], intfcIndex[2], 1);
     //m_pSolver->Add_b(intfcIndex[2], value); 
 }

/*
 * get the flux stencil for the partial edge using 1st order accurate finite difference.
 */
void EBM2D_CARTESIAN::getFluxStencilForPartialEdge(int face, double center[2], 
						   double length, int comp,
						   std::vector< std::pair<int,double> > & stencil,
						   int i, int j)
{
     double dh;
     int I, J;
     I = i;
     J = j;

     switch(face-face%2)
     {
	 case PLANE_X:
	     if(face==FACE_WEST)     // west      
		 I = i-1;
	     else                    // east
		 I = i+1;            
	     dh = getCellEdgeLength(0,i,j);
	     break;
	 case PLANE_Y:
	     if(face==FACE_SOUTH)    // south
		 J = j-1;
	     else                    // north
		 J = j+1;
	     dh = getCellEdgeLength(1,i,j);
	     break;
     }

     stencil.push_back(std::pair<int,double>(getCellUnknownIndex(comp,i,j), -1/dh));
     stencil.push_back(std::pair<int,double>(getCellUnknownIndex(comp,I,J),  1/dh));
}

/*
 * get the flux stencil for the partial edge using 2rd order accurate method.
 *
 * It is assumed that the edge is a real partial edge. The ordering of the four
 * neighbor cell is the following:
 *           |
 *      i3,j3| i2,j2           / \
 *     ------|-----             |
 *       i,j | i1,j1            |
 *           |                  dy
 *      dx --->
 */
void EBM2D_CARTESIAN::getFluxStencilForPartialEdge2(int face, double center[2], 
						    double length, int comp,
						    std::vector< std::pair<int,double> > & stencil,
						    int i, int j)
{
     double dx,dy;
     int i1,j1,i2,j2,i3,j3,c[2];
     
     getCellEdgeComp(face, c, i,j);
     
     switch(face)
     {
     case FACE_SOUTH:
	  i1 = i;	    j1 = j-1;
	  if(c[0]==comp)
	  {
	       i2 = i-1;    j2 = j-1;
	       i3 = i-1;    j3 = j;
	  }
	  else
	  {
	       i2 = i+1;    j2 = j-1;
	       i3 = i+1;    j3 = j;
	  }
	  dx = getCellEdgeLength(1,i,j);
	  dy = getCellEdgeLength(0,i,j);
	  break;
     case FACE_EAST:
	  i1 = i+1;         j1 = j;
	  if(c[0]==comp)
	  {
	       i2 = i+1;    j2 = j-1;
	       i3 = i;      j3 = j-1;
	  }
	  else
	  {
	       i2 = i+1;    j2 = j+1;
	       i3 = i;      j3 = j+1;
	  }
	  dx = getCellEdgeLength(0,i,j);
	  dy = getCellEdgeLength(1,i,j);
	  break;
     case FACE_NORTH:
	  i1 = i;           j1 = j+1;
	  if(c[0]==comp)
	  {
	       i2 = i+1;    j2 = j+1;
	       i3 = i+1;    j3 = j;
	  }
	  else
	  {
	       i2 = i-1;    j2 = j+1;
	       i3 = i-1;    j3 = j;
	  }	  
	  dx = getCellEdgeLength(1,i,j);
	  dy = getCellEdgeLength(0,i,j);
	  break;
     case FACE_WEST:
	  i1 = i-1;         j1 = j;
	  if(c[0]==comp)
	  {
	       i2 = i-1;    j2 = j+1;
	       i3 = i;      j3 = j+1;
	  }
	  else
	  {
	       i2 = i-1;    j2 = j-1;
	       i3 = i;      j3 = j-1;
	  }	  
	  dx = getCellEdgeLength(0,i,j);
	  dy = getCellEdgeLength(1,i,j);
	  break;
     }
     
     double t = (dy/2 - length/2)/dy;
     int index[4];
     index[0] = getCellUnknownIndex(comp,i,j);      
     index[1] = getCellUnknownIndex(comp,i1,j1);
     index[2] = getCellUnknownIndex(comp,i2,j2);
     index[3] = getCellUnknownIndex(comp,i3,j3);
     
     if(index[2]<0||index[3]<0)
     {
	  stencil.push_back(std::pair<int,double>(index[0],-1/dx));
	  stencil.push_back(std::pair<int,double>(index[1], 1/dx));
	  printf("ID %d: EBM2D_CARTESIAN::getFluxStencilForPartialEdge2: 1st order stencil is used!\n");
     }
     else
     {
	  //printf("ID %d: EBM2D_CARTESIAN::getFluxStencilForPartialEdge2: 2st order stencil is used!\n");
	  stencil.push_back(std::pair<int,double>(index[0],-1/dx * (1-t)));
	  stencil.push_back(std::pair<int,double>(index[1], 1/dx * (1-t)));
	  stencil.push_back(std::pair<int,double>(index[3],-1/dx * t));
	  stencil.push_back(std::pair<int,double>(index[2], 1/dx * t));
     }
}

double EBM2D_CARTESIAN::getAreaFraction(int index, std::map<int, double> &areas)
{
     return 1;
     return m_cellArea;

     std::map<int, double>::iterator iter;
     iter = areas.begin();
     while(iter!=areas.end())
     {
	  if(index==(*iter).first)
	       return (*iter).second/m_cellArea;
	  iter++;
     }
     return 1;
}

 /*
  * In the positive normal[] direction, find crossing of the coordinates lines; 
  * then find the 3 points stencil to interpolate the unknowns at the crossing.
  */
void EBM2D_CARTESIAN::getSolutionInterpStencil(double edgeCenter[2], 
					       double normal[2], 
					       int comp,
					       std::vector<int> &indices, 
					       std::vector<double> &coefficients, 
					       std::vector<double> &distances, 
					       int i, int j)
{    
     // first get the direction of the 3 points stencil used for the interpolation
     int plane, dir, sign;
     plane = PLANE_X;
     if(fabs(normal[0])<fabs(normal[1]))
     {
	  plane = PLANE_Y;          //dir = 0;     // PLANE_Y
	  //dir = 0;
     }

     if(normal[plane/2]>0)          //if(normal[1-dir]>0)
	 sign = 1;
     else
	 sign = -1;

     dir = plane==PLANE_X ? 1 : 0;
     double d = getCellEdgeLength(dir,i,j);     // assuming uniform grid       

     int targetI,targetJ;
     targetI = targetJ = -1;
     double t, cellCenter[2], crossing[2];
     double coords[3], coeffs[3];

     indices.clear();
     coefficients.clear();
     distances.clear();    

     for(int k=1; k<=2; k++)
     {
	  if(plane==PLANE_X)
	  {
	       targetI = i + sign * k;
	       crossing[0] = getCellEdgeCenter(0, targetI);
	       crossing[1] = (crossing[0]-edgeCenter[0])*normal[1]/normal[0]+edgeCenter[1];
	  }
	  else      //(plane==PLANE_Y)         //if(dir==0)
	  {
	       targetJ = j + sign * k;
	       crossing[1] = getCellEdgeCenter(1, targetJ);
	       crossing[0] = (crossing[1]-edgeCenter[1])*normal[0]/normal[1]+edgeCenter[0];
	  }
	 
	  locateCell(crossing, targetI, targetJ);
	  getCellCenter(cellCenter, targetI, targetJ);        
	  
	  coords[0] = cellCenter[dir] - d;
	  coords[1] = cellCenter[dir];
	  coords[2] = cellCenter[dir] + d;

	  get3PointsSolutionInterpCoeffs(coords, crossing[dir], coeffs);
	  if(plane==PLANE_X)
	  {
	       indices.push_back(getCellUnknownIndex(comp, targetI,targetJ-1));
	       indices.push_back(getCellUnknownIndex(comp, targetI,targetJ));
	       indices.push_back(getCellUnknownIndex(comp, targetI,targetJ+1));
	  }
	  else //if(plane==PLANE_Y)           //if(dir==0)
	  {
	       indices.push_back(getCellUnknownIndex(comp, targetI-1,targetJ));
	       indices.push_back(getCellUnknownIndex(comp, targetI,targetJ));
	       indices.push_back(getCellUnknownIndex(comp, targetI+1,targetJ));            
	  }
	 
	 coefficients.push_back(coeffs[0]);
	 coefficients.push_back(coeffs[1]);
	 coefficients.push_back(coeffs[2]);
	 distances.push_back(
	      sqrt( (crossing[0]-edgeCenter[0])*(crossing[0]-edgeCenter[0])+
		    +(crossing[1]-edgeCenter[1])*(crossing[1]-edgeCenter[1])
		   ));
     }
 }


/*
 * get derivative perpendicular to plane (positive side)
 */
double EBM2D_CARTESIAN::getDerivativeAtCellCenter(EBM_PLANE plane, int comp, int i, int j)
{
     // 5 neighbors across the plane [-2,-1,0,1,2]
     boolean good[5] = {FALSE,FALSE,FALSE,FALSE,FALSE};

     double dx = getCellEdgeLength(plane/2,i,j);
     double values[3], coords[3], coeffs[3];
     
     good[2] = getNeighborAcrossPlane(plane,0,comp,i,j)>=0 ? TRUE : FALSE;
     if(!good[2])
     {
	  printf("EBM2D_CARTESIAN::getDerivativeAtCellCenter: no comp %d in the current cell \n", comp);
	  exit(0);
     }
     
     // first case: -1,1
     good[1] = getNeighborAcrossPlane(plane,-1,comp,i,j)>=0 ? TRUE : FALSE;
     good[3] = getNeighborAcrossPlane(plane, 1,comp,i,j)>=0 ? TRUE : FALSE;
     if(good[1]&&good[3])
     {
	  values[0] = getCellUnknownAcrossPlane(plane,-1,comp,i,j);
	  values[1] = getCellUnknownAcrossPlane(plane, 1,comp,i,j);
	  return (values[1]-values[0])/(2*dx);
     }

     // other case
     
     if(!good[1])             // right side
     {
	  values[0] = getCellUnknownAcrossPlane(plane, 0,comp,i,j);
	  values[1] = getCellUnknownAcrossPlane(plane, 1,comp,i,j);

	  good[4] = getNeighborAcrossPlane(plane, 2,comp,i,j)>=0 ? TRUE : FALSE;
	  if(good[4])
	  {
	       coords[0] = 0;
	       coords[1] = dx;
	       coords[2] = 2*dx;
	       get3PointsDerivativeInterpCoeffs(coords,0,coeffs);
	       values[2] = getCellUnknownAcrossPlane(plane, 2,comp,i,j);
	       return values[0]*coeffs[0]+values[1]*coeffs[1]+values[2]*coeffs[2];
	  }
	  else
	  {
	       printf("EBM2D_CARTESIAN::getDerivativeAtCellCenter: "
		      "1st order derivative used for cell (%d,%d).\n",i,j);
	       return (values[1]-values[0])/dx;
	  }
     }
     else                     // left side
     {
	  values[1] = getCellUnknownAcrossPlane(plane, -1,comp,i,j);
	  values[2] = getCellUnknownAcrossPlane(plane,  0,comp,i,j);

	  good[0] = getNeighborAcrossPlane(plane, -2,comp,i,j)>=0 ? TRUE : FALSE;
	  if(good[0])
	  {
	       coords[0] =-2*dx;
	       coords[1] =-dx;
	       coords[2] = 0;
	       get3PointsDerivativeInterpCoeffs(coords,0,coeffs);
	       values[0] = getCellUnknownAcrossPlane(plane, -2,comp,i,j);
	       return values[0]*coeffs[0]+values[1]*coeffs[1]+values[2]*coeffs[2];
	  }
	  else
	  {
	       printf("EBM2D_CARTESIAN::getDerivativeAtCellCenter: "
		      "1st order derivative used for cell (%d,%d).\n",i,j);
	       return (values[2]-values[1])/dx;
	  }
     }
}
/*double EBM2D_CARTESIAN::getDerivative(double coords[2], double normal[2], int comp, int i, int j)
{   
}*/

 /*
  * Give 3 points on the x coordinates coords[0], coords[1], coords[2], get the
  * interpolation coefficients for the point with coordinate x.
  */ 
 void EBM2D_CARTESIAN::get3PointsSolutionInterpCoeffs(double coords[3], double x, double coeffs[3])
 {
     double y0,y1,y;
     y0 = coords[1] - coords[0];
     y1 = coords[2] - coords[0];
     y  = x - coords[0];
     coeffs[0] = (y-y0)*(y-y1)/(y0*y1);    
     coeffs[1] = y*(y-y1)/(y0*(y0-y1));
     coeffs[2] = y*(y0-y)/(y1*(y0-y1));
 }

 /*
  * Get the interpolation coefficients for the first order derivative.
  */ 
 void EBM2D_CARTESIAN::get3PointsDerivativeInterpCoeffs(double coords[3], double x, double coeffs[3])
 {
     double y0,y1,y, denominator;
     y0 = coords[1] - coords[0];
     y1 = coords[2] - coords[0];
     y  = x - coords[0];
     coeffs[0] = -(-2*y+y0+y1)/(y0*y1);    
     coeffs[1] = (2*y-y1)/(y0*(y0-y1));
     coeffs[2] = (y0-2*y)/(y1*(y0-y1));
 }    

 
/*
 * set m_lbuf[], m_ubuf[].
 * It is implicitly assumed that 
 *          m_lbuf[0]==m_ubuf[0] && m_lbuf[1]==m_ubuf[1]. 
 */
void EBM2D_CARTESIAN::pp_setup(void)
{
     m_lbuf[0] = m_ubuf[0] = 2;
     m_lbuf[1] = m_ubuf[1] = 2;
}
/*
 * reset the unknown index at each cell to global index, excluding the buffered cell.
 *  1) calculate the range of unknown indices: [m_ilower, m_iupper].
 *  2) reset each cell's unknown index, excluding buffered cells.
 *
 * this function should be called after m_nLocalUnknowns is set.
 */                          
void EBM2D_CARTESIAN::pp_resetCellUnknownIndex(void)                              
{
     int i,j, *n_dist, num_nodes;

     // 1) calculate the range of unknown indices: [m_ilower, m_iupper];
     num_nodes = pp_numnodes();
     FT_VectorMemoryAlloc((POINTER*)&n_dist, num_nodes, sizeof(int));
     n_dist[m_myid] = m_nLocalUnknowns;
     pp_global_imax(n_dist, num_nodes);
     
     //for(i=0; i<num_nodes; i++)
     // printf("ID %d: n_dist[%d]=%d \n", m_myid, i, n_dist[i]);

     m_ilower = 0;
     for(i=1; i<=m_myid; i++)
	  m_ilower += n_dist[i-1];
     m_iupper = m_ilower + m_nLocalUnknowns - 1;
     //printf("ID %d: m_ilower=%d, m_iupper=%d\n", m_myid, m_ilower, m_iupper); 

     // 2) reset each cell's unknown index.
     for(j=m_lbuf[1]; j<m_lbuf[1]+m_gmax[1]; j++)      
     for(i=m_lbuf[0]; i<m_lbuf[0]+m_gmax[0]; i++)
	  getCell(i,j)->resetCellUnknownIndex(m_ilower);
     
     FT_FreeThese(1, n_dist);
}
/*
 * Reset the unknown index at each buffered cell to global index using 
 * information from the node's neighbors.
 */
void EBM2D_CARTESIAN::pp_resetBufferCellUnknownIndex(void)      
{
     int myid = pp_mynode();
     int **sendbuff, **recvbuff, size;

     // in each cell, there are 4 index.
     int length[] = {m_lbuf[0]+m_gmax[0]+m_ubuf[0],
		     m_lbuf[1]+m_gmax[1]+m_ubuf[1]};
     int maxBuffSize = std::max( length[1]*m_lbuf[0],
				 length[0]*m_lbuf[1] );
     size = maxBuffSize * 4;
          
     FT_MatrixMemoryAlloc((POINTER*)&sendbuff, 2, size, sizeof(int));
     FT_MatrixMemoryAlloc((POINTER*)&recvbuff, 2, size, sizeof(int));
     
     //printf("ID %d: size of sendbuff/recvbuff is (%d,%d). \n", pp_mynode(),2,size);

     for(int dir=0; dir<2; dir++)        // for X direction, then Y direction.
     {
	  switch(dir)
	  {
	  case 0:
	       size = length[1]*m_lbuf[0] * 4;
	       break;
	  case 1:
	       size = length[0]*m_lbuf[1] * 4;
	       break;
	  }
	  
	  //printf("ID %d: pp_packCellUnknownIndex(), dir=%d, size=%d\n", myid,dir,size);
	  pp_packCellUnknownIndex(dir, size, sendbuff);

	  //printf("ID %d: sendbuff[0]=0x%x, sendbuff[1]=0x%x\n", pp_mynode(),sendbuff[0],sendbuff[1]);
	  //printf("ID %d: recvbuff[0]=0x%x, recvbuff[1]=0x%x\n", pp_mynode(),recvbuff[0],recvbuff[1]);

	  for(int i=0; i<size; i++)
	  {
	       recvbuff[0][i] = -1;
	       recvbuff[1][i] = -1;
	  }
	  //printf("ID %d: FronTierPPGrid::sendrecv()\n", myid);
	  FronTierPPGrid::sendrecv(dir, m_pFront, (void**)sendbuff, (void**)recvbuff, size, MPI_INT);

	  //printf("ID %d: pp_unpackCellUnknownIndex(), dir = %d\n", myid, dir);
	  pp_unpackCellUnknownIndex(dir, size, recvbuff);
     }
     FT_FreeThese(2, sendbuff, recvbuff);
}
/* 
 * size is the size of **data: 
 *      data[2][size];
 */
void EBM2D_CARTESIAN::pp_packCellUnknownIndex(int dir, int size, int **data)
{
     int i, j, p;
     if(dir==0)        // WEST, EAST    
     {
	  p = 0;
	  for(j=0; j<m_lbuf[1]+m_gmax[1]+m_ubuf[1]; j++)
	       for(i=m_lbuf[0]; i<m_lbuf[0]+m_lbuf[0]; i++)
	       {
		    getCell(i,j)->getCellUnknownIndex(&(data[0][p]));
		    p += 4;
	       }
	  //printf("ID %d: size=%d, p=%d\n", pp_mynode(),size,p);
	  p = 0;
	  for(j=0; j<m_lbuf[1]+m_gmax[1]+m_ubuf[1]; j++)
	       for(i=m_lbuf[0]+m_gmax[0]-m_ubuf[0]; i<m_lbuf[0]+m_gmax[0]; i++)
	       {
		    getCell(i,j)->getCellUnknownIndex(&(data[1][p]));
		    p += 4;
	       }
	  //printf("ID %d: size=%d, p=%d\n", pp_mynode(),size,p);
     }
     else              // SOUTH, NORTH
     {
	  p = 0;
	  for(j=m_lbuf[1]; j<m_lbuf[1]+m_lbuf[1]; j++)
	       for(i=0; i<m_lbuf[0]+m_gmax[0]+m_ubuf[0]; i++)
	       {
		    getCell(i,j)->getCellUnknownIndex(&(data[0][p]));
		    p += 4;
	       }
	  //printf("ID %d: size=%d, p=%d\n", pp_mynode(),size,p);
	  p = 0;
	  for(j=m_lbuf[1]+m_gmax[1]-m_ubuf[1]; j<m_lbuf[1]+m_gmax[1]; j++)
	       for(i=0; i<m_lbuf[0]+m_gmax[0]+m_ubuf[0]; i++)
	       {
		    getCell(i,j)->getCellUnknownIndex(&(data[1][p]));
		    p += 4;
	       }
	  //printf("ID %d: size=%d, p=%d\n", pp_mynode(),size,p);
     }
}

void EBM2D_CARTESIAN::pp_unpackCellUnknownIndex(int dir, int size, int **data)
{
     
     int i, j, p, *index;
     if(dir==0)        // WEST, EAST    
     {
	  p = 0;
	  for(j=0; j<m_lbuf[1]+m_gmax[1]+m_ubuf[1]; j++)
	       for(i=0; i<m_lbuf[0]; i++)
	       {
		    getCell(i,j)->setCellUnknownIndex(&data[0][p]);
		    p += 4;
	       }
	  p = 0;
	  for(j=0; j<m_lbuf[1]+m_gmax[1]+m_ubuf[1]; j++)
	       for(i=m_lbuf[0]+m_gmax[0]; i<m_lbuf[0]+m_gmax[0]+m_ubuf[0]; i++)
	       {
		    getCell(i,j)->setCellUnknownIndex(&data[1][p]);
		    p += 4;
	       }
     }
     else              // SOUTH, NORTH
     {
	  p = 0;
	  for(j=0; j<m_lbuf[1]; j++)
	       for(i=0; i<m_lbuf[0]+m_gmax[0]+m_ubuf[0]; i++)
	       {
		    getCell(i,j)->setCellUnknownIndex(&data[0][p]);
		    p += 4;
	       }
	  p = 0;
	  for(j=m_lbuf[1]+m_gmax[1]; j<m_lbuf[1]+m_gmax[1]+m_ubuf[1]; j++)
	       for(i=0; i<m_lbuf[0]+m_gmax[0]+m_ubuf[0]; i++)
	       {
		    getCell(i,j)->setCellUnknownIndex(&data[1][p]);
		    p += 4;
	       }
     }
}					  

void EBM2D_CARTESIAN::pp_resetBufferCellCompList(void)
{
     int **sendbuff, **recvbuff, size;

     // in each cell, there are 2 unknowns.
     int maxBuffSize = std::max( (m_lbuf[0]+m_gmax[0]+m_ubuf[0])*m_lbuf[1],
				 (m_lbuf[1]+m_gmax[1]+m_ubuf[1])*m_lbuf[0] );
     size = maxBuffSize * 2;
     FT_MatrixMemoryAlloc((POINTER*)&sendbuff, 2, size, sizeof(int));
     FT_MatrixMemoryAlloc((POINTER*)&recvbuff, 2, size, sizeof(int));
     
     for(int dir=0; dir<2; dir++)
     {
	  if(dir==0)
	       size = (m_lbuf[1]+m_gmax[1]+m_ubuf[1])*m_lbuf[0] * 2;
	  else
	       size = (m_lbuf[0]+m_gmax[0]+m_ubuf[0])*m_lbuf[1] * 2;
	  
	  pp_packCellCompList(dir, size, sendbuff);
	  //printf("ID %d: dir=%d\n", m_myid, dir);
	  //printf("ID %d: size=%d\n", m_myid, size);
	  //for(int i=0; i<size; i++)
	  //     printf("ID %d: sendbuf[*][%d]={%d,%d}\n", m_myid,i,sendbuff[0][i],sendbuff[1][i]);
	  for(int i=0; i<size; i++)
	  {
	       recvbuff[0][i] = -1;
	       recvbuff[1][i] = -1;
	  }
	  FronTierPPGrid::sendrecv(dir, m_pFront, (void**)sendbuff, (void**)recvbuff, size, MPI_INT);
	  pp_unpackCellCompList(dir, size, recvbuff);

	  //for(int i=0; i<size; i++)
	  //     printf("ID %d: recvbuf[*][%d]={%d,%d}\n", m_myid,i,recvbuff[0][i],recvbuff[1][i]);
     }
     FT_FreeThese(2, sendbuff, recvbuff);
}
void EBM2D_CARTESIAN::pp_packCellCompList(int dir, int size, int **data)
{
     int i, j, p;

     if(dir==0)        // WEST, EAST    
     {
	  p = 0;
	  for(j=0; j<m_lbuf[1]+m_gmax[1]+m_ubuf[1]; j++)
	       for(i=m_lbuf[0]; i<m_lbuf[0]+m_lbuf[0]; i++)
	       {
		    getCellCompList(&data[0][p],i,j);
		    p += 2;
	       }
	  p = 0;
	  for(j=0; j<m_lbuf[1]+m_gmax[1]+m_ubuf[1]; j++)
	       for(i=m_lbuf[0]+m_gmax[0]-m_ubuf[0]; i<m_lbuf[0]+m_gmax[0]; i++)
	       {
		    getCellCompList(&data[1][p],i,j);
		    p += 2;
	       }
     }
     else              // SOUTH, NORTH
     {
	  p = 0;
	  for(j=m_lbuf[1]; j<m_lbuf[1]+m_lbuf[1]; j++)
	       for(i=0; i<m_lbuf[0]+m_gmax[0]+m_ubuf[0]; i++)
	       {
		    getCellCompList(&data[0][p],i,j);
		    p += 2;
	       }
	  p = 0;
	  for(j=m_lbuf[1]+m_gmax[1]-m_ubuf[1]; j<m_lbuf[1]+m_gmax[1]; j++)
	       for(i=0; i<m_lbuf[0]+m_gmax[0]+m_ubuf[0]; i++)
	       {
		    getCellCompList(&data[1][p],i,j);
		    p += 2;
	       }
     }
}
void EBM2D_CARTESIAN::pp_unpackCellCompList(int dir, int size, int **data)
{
     int i, j, p;
     if(dir==0)        // WEST, EAST    
     {
	  p = 0;
	  for(j=0; j<m_lbuf[1]+m_gmax[1]+m_ubuf[1]; j++)
	       for(i=0; i<m_lbuf[0]; i++)
	       {
		    setCellCompList(&data[0][p],i,j);
		    p += 2;
	       }
	  p = 0;
	  for(j=0; j<m_lbuf[1]+m_gmax[1]+m_ubuf[1]; j++)
	       for(i=m_lbuf[0]+m_gmax[0]; i<m_lbuf[0]+m_gmax[0]+m_ubuf[0]; i++)
	       {
		    setCellCompList(&data[1][p],i,j);
		    p += 2;
	       }
     }
     else              // SOUTH, NORTH
     {
	  p = 0;
	  for(j=0; j<m_lbuf[1]; j++)
	       for(i=0; i<m_lbuf[0]+m_gmax[0]+m_ubuf[0]; i++)
	       {
		    setCellCompList(&data[0][p],i,j);
		    p += 2;
	       }
	  p = 0;
	  for(j=m_lbuf[1]+m_gmax[1]; j<m_lbuf[1]+m_gmax[1]+m_ubuf[1]; j++)
	       for(i=0; i<m_lbuf[0]+m_gmax[0]+m_ubuf[0]; i++)
	       {
		    setCellCompList(&data[1][p],i,j);
		    p += 2;
	       }
     }
}


/*
 * reset the unknown for the buffered cells.
 */
void EBM2D_CARTESIAN::pp_resetBufferCellUnknown(void)
{
     double **sendbuff, **recvbuff, size;

     // in each cell, there are 2 unknowns.
     int maxBuffSize = std::max( (m_lbuf[0]+m_gmax[0]+m_ubuf[0])*m_lbuf[1],
				 (m_lbuf[1]+m_gmax[1]+m_ubuf[1])*m_lbuf[0] );
     size = maxBuffSize * 2;
     FT_MatrixMemoryAlloc((POINTER*)&sendbuff, 2, size, sizeof(double));
     FT_MatrixMemoryAlloc((POINTER*)&recvbuff, 2, size, sizeof(double));
     
     for(int dir=0; dir<2; dir++)
     {
	  if(dir==0)
	       size = (m_lbuf[1]+m_gmax[1]+m_ubuf[1])*m_lbuf[0] * 2;
	  else
	       size = (m_lbuf[0]+m_gmax[0]+m_ubuf[0])*m_lbuf[1] * 2;

	  pp_packCellUnknown(dir, size, sendbuff);
	  FronTierPPGrid::sendrecv(dir, m_pFront, (void**)sendbuff, (void**)recvbuff, size, MPI_DOUBLE);
	  pp_unpackCellUnknown(dir, size, recvbuff);
     }
     FT_FreeThese(2, sendbuff, recvbuff);
}
void EBM2D_CARTESIAN::pp_packCellUnknown(int dir, int size, double **data)
{
     int i, j, p;

     if(dir==0)        // WEST, EAST    
     {
	  p = 0;
	  for(j=0; j<m_lbuf[1]+m_gmax[1]+m_ubuf[1]; j++)
	       for(i=m_lbuf[0]; i<m_lbuf[0]+m_lbuf[0]; i++)
	       {
		    getCellUnknown(&data[0][p],i,j);
		    p += 2;
	       }
	  p = 0;
	  for(j=0; j<m_lbuf[1]+m_gmax[1]+m_ubuf[1]; j++)
	       for(i=m_lbuf[0]+m_gmax[0]-m_ubuf[0]; i<m_lbuf[0]+m_gmax[0]; i++)
	       {
		    getCellUnknown(&data[1][p],i,j);
		    p += 2;
	       }
     }
     else              // SOUTH, NORTH
     {
	  p = 0;
	  for(j=m_lbuf[1]; j<m_lbuf[1]+m_lbuf[1]; j++)
	       for(i=0; i<m_lbuf[0]+m_gmax[0]+m_ubuf[0]; i++)
	       {
		    getCellUnknown(&data[0][p],i,j);
		    p += 2;
	       }
	  p = 0;
	  for(j=m_lbuf[1]+m_gmax[1]-m_ubuf[1]; j<m_lbuf[1]+m_gmax[1]; j++)
	       for(i=0; i<m_lbuf[0]+m_gmax[0]+m_ubuf[0]; i++)
	       {
		    getCellUnknown(&data[1][p],i,j);
		    p += 2;
	       }
     }
}
void EBM2D_CARTESIAN::pp_unpackCellUnknown(int dir, int size, double **data)
{
     int i, j, p;
     if(dir==0)        // WEST, EAST    
     {
	  p = 0;
	  for(j=0; j<m_lbuf[1]+m_gmax[1]+m_ubuf[1]; j++)
	       for(i=0; i<m_lbuf[0]; i++)
	       {
		    setCellUnknown(&data[0][p],i,j);
		    p += 2;
	       }
	  p = 0;
	  for(j=0; j<m_lbuf[1]+m_gmax[1]+m_ubuf[1]; j++)
	       for(i=m_lbuf[0]+m_gmax[0]; i<m_lbuf[0]+m_gmax[0]+m_ubuf[0]; i++)
	       {
		    setCellUnknown(&data[1][p],i,j);
		    p += 2;
	       }
     }
     else              // SOUTH, NORTH
     {
	  p = 0;
	  for(j=0; j<m_lbuf[1]; j++)
	       for(i=0; i<m_lbuf[0]+m_gmax[0]+m_ubuf[0]; i++)
	       {
		    setCellUnknown(&data[0][p],i,j);
		    p += 2;
	       }
	  p = 0;
	  for(j=m_lbuf[1]+m_gmax[1]; j<m_lbuf[1]+m_gmax[1]+m_ubuf[1]; j++)
	       for(i=0; i<m_lbuf[0]+m_gmax[0]+m_ubuf[0]; i++)
	       {
		    setCellUnknown(&data[1][p],i,j);
		    p += 2;
	       }
     }
}


 /*
  * Setup the data structure.        
  * 
  */
 void EBM2D_CARTESIAN::setup(void)
 {
      m_myid = pp_mynode();
      printf("ID %d: EBM2D_CARTESIAN::setup: \n", m_myid);

      int i, j;

     // the following init is moved here 
     // because EBM2D_CARTESIAN::getComp() need m_mindd.
     m_i = m_j = -1;
     m_dx = getCellEdgeLength(0,0,0);
     m_dy = getCellEdgeLength(1,0,0);
     m_mindd = std::min(m_dx,m_dy);
     m_maxdd = std::max(m_dx,m_dy);
     m_cellArea = m_dx*m_dy;

     m_mindd4 = m_mindd*m_mindd*m_mindd*m_mindd;
     m_MAX_ITER = std::max(30, (int) (-4*log(m_mindd)+1) );
     
     printf("ID %d: \tm_MAX_ITER = %d\n", m_myid, m_MAX_ITER);

     // m_gmax[]
     RECT_GRID *rect_grid = m_pFront->rect_grid;
     m_gmax[0] = rect_grid->gmax[0];
     m_gmax[1] = rect_grid->gmax[1];
     pp_setup();          // m_lbuff, m_ubuff.

     // m_cells
     FT_MatrixMemoryAlloc((POINTER*)&m_cells, 
	      m_lbuf[0]+m_gmax[0]+m_ubuf[0], 
	      m_lbuf[1]+m_gmax[1]+m_ubuf[1], 
	      sizeof(EBM2D_CELL));
     for(j=0; j<m_lbuf[1]+m_gmax[1]+m_ubuf[1]; j++)
	  for(i=0; i<m_lbuf[0]+m_gmax[0]+m_ubuf[0]; i++)
	       m_cells[i][j].init();    // call the constructor

     // m_comps
     double p[2];
     m_nComps = (m_lbuf[0]+m_gmax[0]+m_ubuf[0]+1)*(m_lbuf[1]+m_gmax[1]+m_ubuf[1]+1);
     FT_VectorMemoryAlloc((POINTER*)&m_comps, m_nComps, sizeof(int));
     for(j=0; j<=m_lbuf[1]+m_gmax[1]+m_ubuf[1]; j++)
	 for(i=0; i<=m_lbuf[0]+m_gmax[0]+m_ubuf[0]; i++)   
	 {
	      getCellCornerCoords(0,p,i,j);
	      m_comps[getCellCornerIndex(0,i,j)] = getComp(p);
	 }

     // the following is used to modify the components from 0 to 2
     // so that we only have components 1, 2 left.
     setEquivalentComponents();  

     // m_values
     m_nLocalUnknowns = 0;
     int numCellUnknowns;
     int comp[5];        // comp[4] for the comp of the center;
     
     CELL_TYPE type;

     for(j=m_lbuf[1]; j<m_lbuf[1]+m_gmax[1]; j++)
	 for(i=m_lbuf[0]; i<m_lbuf[0]+m_gmax[0]; i++)
	 {            
	      getCellComp(comp,i,j);
	      getCell(i,j)->setComp(comp);
	      getCell(i,j)->setCellType(comp);
	      type = getCell(i,j)->getCellType();
	      
	      if(type==EXTERNAL)
		   continue;
	      else                           // INTERNAL, PARTIAL, BOUNDARY
	      {
		   numCellUnknowns = getCell(i,j)->setCellUnknownIndex(m_nLocalUnknowns);
		   numCellUnknowns = numCellUnknowns * 2 - 1;
		   //if(m_nLocalUnknowns<=23&&23<m_nLocalUnknowns+numCellUnknowns)
		   //printf("debug here!\n");
		   m_nLocalUnknowns += numCellUnknowns;
	      }
        }    
    
     for(int d=0; d<2; d++)
     {
	  FT_MatrixMemoryAlloc((POINTER*)&m_unknowns[d], 
		   (m_lbuf[0]+m_gmax[0]+m_ubuf[0]), (m_lbuf[1]+m_gmax[1]+m_ubuf[1]),
		   sizeof(double));
	  FT_MatrixMemoryAlloc((POINTER*)&m_compList[d], 
		   (m_lbuf[0]+m_gmax[0]+m_ubuf[0]), (m_lbuf[1]+m_gmax[1]+m_ubuf[1]),
		   sizeof(int));

	  for(j=0; j<m_lbuf[1]+m_gmax[1]+m_ubuf[1]; j++)
	       for(i=0; i<m_lbuf[0]+m_gmax[0]+m_ubuf[0]; i++)
	       {
		    m_unknowns[d][i][j] = 0;
		    m_compList[d][i][j] = -1;
	       }
     }
    
}

/*
 * The following functions are used for the EBM2D methods.     
 */
void EBM2D_CARTESIAN::setMatrix(void)
{    
     // this is very important to speed up the matrix setup!
     //m_pSolver->Create(0, m_nLocalUnknowns-1, 12, 0);
     m_pSolver->Create(m_ilower, m_iupper, 12, 0);
        
    for(int j=m_lbuf[1]; j<m_lbuf[1]+m_gmax[1]; j++)
        for(int i=m_lbuf[0]; i<m_lbuf[0]+m_gmax[0]; i++)   
        {   
	    setCurrentCell(i,j);           // this is good for debugging.

	    switch(getCell(i,j)->getCellType())
	    {
	    case INTERNAL:
		 setMatrixForInternalCell(i,j);
		 break;
	    case PARTIAL:
		 debug_getMinStatistics(i,j);
		 setMatrixForPartialCell_2Unknowns(i,j);
		 //setMatrixForPartialCell_2Unknowns_2ndOrder(i,j);
		 break;
	    case BOUNDARY:
		 setMatrixForBoundaryCell(i,j);
		 break;
	    }	    
        }    
}


 EBM2D_CARTESIAN::EBM2D_CARTESIAN(Front *front, SOLVER *solver, EBM2D_LAPLACE *laplace)
 {
      m_pFront = front;
      m_pSolver = solver;
      m_pLaplace = laplace;
 }

 EBM2D_CARTESIAN::~EBM2D_CARTESIAN()
 {
      FT_FreeThese(2, m_cells, m_comps);
      FT_FreeThese(4, m_unknowns[0], m_unknowns[1], m_compList[0], m_compList[1]);
 }

void EBM2D_CARTESIAN::solve(void)
{
     setup();
     printf("ID %d: EBM2D_CARTESIAN::solve: \n", m_myid);
     printf("ID %d: \tm_gmax = {%d,%d}\n", m_myid, m_gmax[0],m_gmax[1]);
     printf("ID %d: \tm_dx = %e, m_dy = %e\n", m_myid, m_dx, m_dy);
     printf("ID %d: \tm_nLocalUnknowns = %d\n", m_myid, m_nLocalUnknowns);
   
     printf("ID %d: \tpp_resetCellUnknownIndex().\n", m_myid);
     pp_resetCellUnknownIndex();
     printf("ID %d: \tpp_resetBufferCellUnknownIndex().\n", m_myid);
     pp_resetBufferCellUnknownIndex();
     
     debug_saveUnknownIndex();
     //return;

     printf("ID %d: \tsetMatrix()\n", m_myid);
     setMatrix();
     printf("ID %d: \tm_minCellArea/m_cellArea = %f\n", 
	    m_myid, m_minCellArea/m_cellArea);
     printf("ID %d: \tm_minPartialEdgeLength/m_dx = %f\n", 
	    m_myid, m_minPartialEdgeLength/m_dx);
     printf("ID %d: \tm_minIntfcLength/m_dx = %f\n", 
	    m_myid, m_minIntfcLength/m_dx);

     printf("ID %d: \tsolve()\n", m_myid);
     m_pSolver->Solve();    
     //m_pSolver->Solve_withPureNeumann();

     int iter;
     double residual, maxSV, minSV;
     maxSV = minSV = 0;
     m_pSolver->GetNumIterations(&iter);
     m_pSolver->GetFinalRelativeResidualNorm(&residual);            
     m_pSolver->GetExtremeSingularValues(&maxSV,&minSV);
     printf("ID %d: \tnumber of iterations = %d\n", m_myid, iter);
     printf("ID %d: \trel residual norm = %e\n", m_myid, residual);
     printf("ID %d: \textreme singular values: maxSV=%f, minSV=%f\n", 
	    m_myid, maxSV, minSV);
     
     double *x = new double[m_nLocalUnknowns];
     m_pSolver->Get_x(x);
     //for(int i=0; i<m_nLocalUnknowns; i++)
     //  printf("\t x[%d]=%e\n", i, x[i]);
     
     //debug_saveCompList("comp0");
     setCellUnknown(x);
     delete [] x;

     //debug_saveCompList("comp1");
     pp_resetBufferCellCompList();
     
     //debug_saveCompList("comp2");

     pp_resetBufferCellUnknown();

     
}

///////////////////////////////////////////////////////////////////////////
// The following functions are used to debug the program.    
///////////////////////////////////////////////////////////////////////////
/*
 * Save the interface into a tecplot file.
 */
void EBM2D_CARTESIAN::saveInterface_Tecplot(const char *filename)
{
    RECT_GRID *rect_grid = m_pFront->rect_grid;
    INTERFACE *intfc    = m_pFront->interf;

    int i, j;
    int xmax = rect_grid->gmax[0];
    int ymax = rect_grid->gmax[1];
    double x, y;

    FILE *hfile = fopen(filename, "w");
    if(hfile==NULL)
    {
        printf("EBM2D_CARTESIAN::saveInterface_Tecplot: can't open %s \n", 
                filename);
        exit(0);
    }

    // secondly print out the interface

    if(exists_interface(intfc))
    {
        CURVE		**curs;
        CURVE		*curve;
        BOND		*bond;
        
        fprintf(hfile, "VARIABLES = X, Y, VALUES \n");
        
        for(curs=intfc->curves; curs && *curs; curs++)	
        {
            curve = *curs;
            
            fprintf(hfile, "ZONE I=%d J=%d F=POINT \n", curve->num_points, 1);
            bond=curve->first;
            fprintf(hfile, "%.4f %.4f 1 \n", 
                    bond->start->_coords[0], bond->start->_coords[1]);
            for(bond=curve->first; bond!=NULL; bond=bond->next)
                fprintf(hfile, "%.4f %.4f 1 \n", 
                        bond->end->_coords[0], bond->end->_coords[1]);	
        }					
    }
    
    
    fprintf(hfile, "VARIABLES = X, Y, VALUES \n");
    fprintf(hfile, "ZONE I=%d, J=%d, F=POINT\n", xmax+1, ymax+1);
    
    
    for(int j=0; j<=ymax; j++)
        for(int i=0; i<=xmax; i++)
        {
            fprintf(hfile, "%f %f %d \n",
                    cell_edge(i,0,rect_grid), cell_edge(j,1,rect_grid),
                    1);
        }
    
    fclose(hfile);
}

void EBM2D_CARTESIAN::saveReconstructedIntfc_Tecplot(const char *filename, int i, int j)
{
    int comp[4], index, size;
    double coords[4][2], crx[5][2];
    
    std::map<int, double> areas;
    std::map<int, std::vector<EBM2D_EDGE> > intfc;    
    std::map<int, std::vector<EBM2D_EDGE> >::iterator iterMap;
    
    EBM2D_EDGE edge;
    
    FILE *hfile;
    if((hfile = fopen(filename,"w")) == NULL )
    {
	 (void) printf("EBM2D_CARTESIAN::saveReconstructedIntfc_Tecplot: "
		       "can't open %s\n",filename);
	 return;
    }
    
    fprintf(hfile, "VARIABLES = X Y\n");
    
    for(int pj=m_lbuf[1]; pj<m_lbuf[1]+m_gmax[1]; pj++)
	 for(int pi=m_lbuf[0]; pi<m_lbuf[0]+m_gmax[0]; pi++)
	 {
	      getCellCompCoordsCrx(comp,coords,crx,pi,pj);
	      //getCell(pi,pj)->init();
	      getCell(pi,pj)->setCompCoordsCrx(comp,coords,crx);
              
	      areas.clear();
	      intfc.clear();
	      getCell(pi,pj)->getAreaAndIntfc(areas, intfc);
              
	      iterMap = intfc.begin();                    
	      while(iterMap!=intfc.end())
	      {
		   size = (*iterMap).second.size();
		   fprintf(hfile, "ZONE N=%d, E=%d, DATAPACKING=POINT, ZONETYPE=FELINESEG\n",
			   2*size, size);
		   for(int m=0; m<(*iterMap).second.size(); m++)
		   {            
                        edge = (*iterMap).second[m];
                        for(int n=0; n<2; n++)
			     fprintf(hfile,"%-9g %-9g \n", 
				     edge.m_vertex[n][0], edge.m_vertex[n][1]);                     
		   }        
		   for(int m=0; m<(*iterMap).second.size(); m++)
                        fprintf(hfile, "%d %d\n", 2*m+1, 2*m+2);
		   
		   iterMap++;
	      }   
	 }
    
    (void) fclose(hfile);    
}


/*
 * Save the edges of the emb_grid_intfc.            
 */
void EBM2D_CARTESIAN::saveInterfaceEdges_Tecplot(const char *filename)
{   
    Table *T = table_of_interface(m_pFront->emb_grid_intfc);
    BLK_EDGE	**blk_edge = T->blk_edge;
    
    RECT_GRID *rect_grid = m_pFront->rect_grid;
    
    FILE *hfile=fopen(filename, "w");
    if(hfile==NULL)
    {
        printf("EBM2D_CARTESIAN::saveInterfaceEdges_Tecplot: can't open %s\n",
                filename);
        exit(0);
    }
    
    fprintf(hfile, "VARIABLES = X, Y, VALUES \n");
    fprintf(hfile, "ZONE I=%d, J=%d, F=POINT\n", m_gmax[0], m_gmax[1]);
    
    
    for(int j=m_lbuf[1]; j<m_lbuf[1]+m_gmax[1]; j++)
        for(int i=m_lbuf[0]; i<m_lbuf[0]+m_gmax[0]; i++)
        {
            fprintf(hfile, "%f %f %d \n",
                    cell_center(i,0,rect_grid), cell_center(j,1,rect_grid),
                    blk_edge[i][j].ctype);
        }
    
    fclose(hfile);    
}



/*
 * save boundaries of the partial cells including the partialfaces.
 * mapPartialFaces contains that boundaries of 6 faces of the cells;
 * mapIntfc contains the intfc inside the cells and the key is the component at its left side.
 *
 * EBM3D_CELL::getPartialFaces() generates mapPartialFaces;
 * EBM3D_CELL::getVolumeAndIntfc() generates mapIntfc;
 *
 */
void EBM2D_CARTESIAN::savePartialCells_Tecplot(const char *filename,
					       std::map<int, std::vector<EBM2D_PARTIALEDGE> > &mapPartialEdges,
					       std::map<int, std::vector<EBM2D_EDGE> > &mapIntfc,
					       int i, int j)
{
     int nedge;
     
     FILE *hfile = fopen(filename, "w");
     if(hfile==NULL)
     {
	  printf("EBM2D_CARTESIAN::savePartialCells_Tecplot: can't open %s for writing!\n", filename);
	  exit(0);
     }

     fprintf(hfile, "VARIABLES = X Y \n");
     fprintf(hfile, "# (%d,%d) \n", i,j);
     std::map<int, std::vector<EBM2D_PARTIALEDGE> >::iterator iterPartialEdges;
     iterPartialEdges = mapPartialEdges.begin();
     while(iterPartialEdges!=mapPartialEdges.end())
     {
	  fprintf(hfile, "# edge = %d \n", (*iterPartialEdges).first);
	  
	  std::vector<EBM2D_PARTIALEDGE> &partialedges = (*iterPartialEdges).second;
	  for(int m = 0; m<partialedges.size(); m++)
	  {
	       std::vector<EBM2D_EDGE> &edges = partialedges[m].m_edges;
	       nedge = edges.size();
	       if(nedge>0)
	       {
		    fprintf(hfile, "# comp = %d\n", partialedges[m].m_comp);
		    fprintf(hfile, "ZONE N=%d, E=%d, DATAPACKING=POINT, ZONETYPE=FELINESEG \n", nedge*2, nedge);
		    for(int n=0; n<nedge; n++)
		    {
			 for(int nn=0; nn<2; nn++)
			      fprintf(hfile, "%e %e\n", 
				      edges[n].m_vertex[nn][0], 
				      edges[n].m_vertex[nn][1]);
				    
		    }
		    for(int n=0; n<nedge; n++)
			 fprintf(hfile, "%d %d\n", 2*n+1, 2*n+2);
	       }
	       
	  }
	  iterPartialEdges++;
     }

     fprintf(hfile, "# intfc\n");
     std::map<int, std::vector<EBM2D_EDGE> >::iterator iterIntfc;
     iterIntfc = mapIntfc.begin();
     while(iterIntfc!=mapIntfc.end())
     {
	  std::vector<EBM2D_EDGE> &edges = (*iterIntfc).second;
	  nedge = edges.size();
	  fprintf(hfile, "ZONE N=%d, E=%d, DATAPACKING=POINT, ZONETYPE=FELINESEG \n", nedge*2, nedge);
	  for(int n=0; n<nedge; n++)
	  {
	       for(int nn=0; nn<2; nn++)
		    fprintf(hfile, "%e %e\n", 
			    edges[n].m_vertex[nn][0], 
			    edges[n].m_vertex[nn][1]);			  
	  }
	  for(int n=0; n<nedge; n++)
	       fprintf(hfile, "%d %d\n", 2*n+1, 2*n+2);
	       iterIntfc++;
     }

     fclose(hfile);
    
/*    //printf("cell = 0x%x\n", getCell(20,15));
    RECT_GRID *rect_grid = m_pFront->rect_grid;    
    
    FILE *hfile=fopen(filename, "w");
    if(hfile==NULL)
    {
        printf("EBM2D_CARTESIAN::savePartialCells_Tecplot: can't open %s \n", 
                filename);
        exit(0);
    }
      
    int numCellUnknowns, comp[5];
    std::vector< std::pair<double,double> > points;
    
    for(int j=0; j<m_gmax[1]; j++)
        for(int i=0; i<m_gmax[0]; i++)
        {   
	     numCellUnknowns = getCellUnknownNumber(comp,i,j);
            
            if(numCellUnknowns==1)
                continue;           
            
            //if(!(j==15&&i==20))
            //    continue;
            
            for(int corner=0; corner<4; corner++)
            {
                getCellUnknownBoundary(comp, numCellUnknowns, 
				       corner, points,i,j);
                
                fprintf(hfile, "VARIABLES = X, Y, VALUES \n");
                fprintf(hfile, "ZONE I=%d, J=%d, F=POINT\n", points.size()+1, 1);
                for(int k=0; k<points.size(); k++)
                    fprintf(hfile, "%f %f %f\n", 
                            points[k].first, points[k].second,
                            0.0);
                fprintf(hfile, "%f %f %f\n", 
                            points[0].first, points[0].second,
                            0.0);
            }
        }            
	fclose(hfile); */
}
void EBM2D_CARTESIAN::saveStates_Tecplot(void)
{
     char filename[100];
     sprintf(filename, "states_%d.plt", m_myid);

    RECT_GRID *rect_grid = m_pFront->rect_grid;
    INTERFACE *intfc    = m_pFront->interf;

    FILE *hfile=fopen(filename, "w");
    if(hfile==NULL)
    {
        printf("EBM2D_CARTESIAN::saveStates_Tecplot: can't open %s \n", filename);
        exit(0);
    }
    
    fprintf(hfile, "VARIABLES = X, Y, Approximate, Exact, MaxError \n");
    fprintf(hfile, "ZONE I=%d, J=%d, F=POINT\n", m_gmax[0], m_gmax[1]);
    //fprintf(hfile, "ZONE I=%d, J=%d, F=POINT\n", m_lbuf[0]+m_gmax[0]+m_ubuf[0], m_lbuf[1]+m_gmax[1]+m_ubuf[1]);
    
    int comp[5], max_i, max_j;
    double max_error = 0, approximate, exact, error, tmp;
    double cellCenter[2];
    max_i = max_j = 0;
    //for(int j=0; j<m_lbuf[1]+m_gmax[1]+m_ubuf[1]; j++)
    //    for(int i=0; i<m_lbuf[0]+m_gmax[0]+m_ubuf[0]; i++)
    for(int j=m_lbuf[1]; j<m_lbuf[1]+m_gmax[1]; j++)
        for(int i=m_lbuf[0]; i<m_lbuf[0]+m_gmax[0]; i++)
        {    
	     getCellComp(comp,i,j);
	     getCellCenter(cellCenter,i,j);
	     approximate = exact = error = 0;
	     for(int k=0; k<4; k++)
	     {
		  if(comp[k]<0)
		       continue;
		  //tmp = m_unknowns[getCell(i,j)->getCellUnknownIndex(k)];                            
		  tmp = getCellUnknown(comp[k],i,j);
		  tmp = fabs( tmp - m_pLaplace->getExactSolution(cellCenter, comp[k]) );
		  if(tmp>error)
		  {
		       //approximate = m_unknowns[getCell(i,j)->getCellUnknownIndex(k)];
		       approximate = getCellUnknown(comp[k],i,j);
		       exact = m_pLaplace->getExactSolution(cellCenter, comp[k]);
		       error = tmp;
		  }
		  if(error>max_error)
		  {
		       max_i = i;
		       max_j = j;
		       max_error = error;
		  }
		  //exact = m_pLaplace->getExactSolution(cellCenter, comp[k]);
	     }     

	     fprintf(hfile, "%f %f %e %e %e\n", 
		     cellCenter[0], cellCenter[1],
		     approximate, exact, error);                        
        }
            
    fclose(hfile);    
    
    printf("ID %d: EBM2D_CARTESIAN::saveStates_Tecplot: \n", m_myid);
    printf("ID %d: \tm_gmax[]={%d,%d}\n", m_myid, m_gmax[0],m_gmax[1]);
    printf("ID %d: \tthe max error is %f in cell (%d,%d)\n",m_myid, max_error,max_i,max_j);
}

void EBM2D_CARTESIAN::saveStates_VTK(void)
{
     char filename[100];
     sprintf(filename, "states_%d.vtk", m_myid);

    RECT_GRID *rect_grid = m_pFront->rect_grid;
    INTERFACE *intfc    = m_pFront->interf;

    FILE *hfile=fopen(filename, "w");
    if(hfile==NULL)
    {
        printf("EBM2D_CARTESIAN::saveStates_VTK: can't open %s \n", filename);
        exit(0);
    }

    double cellCenter[2];    
    double origin[2];
    getCellCornerCoords(0,origin,m_lbuf[0],m_lbuf[1]);
    
    fprintf(hfile,"# vtk DataFile Version 3.0\n");
    fprintf(hfile,"Max_Error\n");
    fprintf(hfile,"ASCII\n");
    fprintf(hfile,"DATASET STRUCTURED_POINTS\n");
    fprintf(hfile,"DIMENSIONS %d %d 1\n", m_gmax[0]+1,m_gmax[1]+1);
    fprintf(hfile,"SPACING %f %f 0\n", m_dx, m_dy);
    fprintf(hfile,"ORIGIN %f %f 0\n", origin[0], origin[1]);
    fprintf(hfile,"CELL_DATA %d\n", m_gmax[0]*m_gmax[1]);
    fprintf(hfile,"SCALARS Max_Error double 1\n");
    fprintf(hfile,"LOOKUP_TABLE DEFAULT\n");
    
    int comp[5], max_i, max_j;
    double max_error = 0, approximate, exact, error, tmp;

    max_i = max_j = 0;
    //for(int j=0; j<m_lbuf[1]+m_gmax[1]+m_ubuf[1]; j++)
    //    for(int i=0; i<m_lbuf[0]+m_gmax[0]+m_ubuf[0]; i++)
    for(int j=m_lbuf[1]; j<m_lbuf[1]+m_gmax[1]; j++)
        for(int i=m_lbuf[0]; i<m_lbuf[0]+m_gmax[0]; i++)
        {    
	     getCellComp(comp,i,j);
	     getCellCenter(cellCenter,i,j);
	     approximate = exact = error = 0;
	     for(int k=0; k<4; k++)
	     {
		  if(comp[k]<0)
		       continue;
		  //tmp = m_unknowns[getCell(i,j)->getCellUnknownIndex(k)];                            
		  tmp = getCellUnknown(comp[k],i,j);
		  tmp = fabs( tmp - m_pLaplace->getExactSolution(cellCenter, comp[k]) );
		  if(tmp>error)
		  {
		       //approximate = m_unknowns[getCell(i,j)->getCellUnknownIndex(k)];
		       approximate = getCellUnknown(comp[k],i,j);
		       exact = m_pLaplace->getExactSolution(cellCenter, comp[k]);
		       error = tmp;
		  }
		  if(error>max_error)
		  {
		       max_i = i;
		       max_j = j;
		       max_error = error;
		  }
		  //exact = m_pLaplace->getExactSolution(cellCenter, comp[k]);
	     }     

	     //fprintf(hfile, "%f %f %e %e %e\n", 
	     //     cellCenter[0], cellCenter[1],
	     //     approximate, exact, error);
             fprintf(hfile, "%f\n", error);
        }
            
    fclose(hfile);    
    
    printf("ID %d: EBM2D_CARTESIAN::saveStates_VTK: \n", m_myid);
    printf("ID %d: \tm_gmax[]={%d,%d}\n", m_myid, m_gmax[0],m_gmax[1]);
    printf("ID %d: \tthe max error is %f in cell (%d,%d)\n",m_myid, max_error,max_i,max_j);
}

void EBM2D_CARTESIAN::saveDerivatives_Tecplot(void)
{
     char filename[100];
     sprintf(filename, "derivatives_%d.plt", m_myid);

    RECT_GRID *rect_grid = m_pFront->rect_grid;
    INTERFACE *intfc    = m_pFront->interf;

    FILE *hfile=fopen(filename, "w");
    if(hfile==NULL)
    {
        printf("EBM2D_CARTESIAN::saveDerivatives_Tecplot: can't open %s \n", filename);
        exit(0);
    }
    
    fprintf(hfile, "VARIABLES = X, Y, App_dx, Exact_dx, MaxError_dx, App_dy, Exact_dy, MaxError_dy \n");
    fprintf(hfile, "ZONE I=%d, J=%d, F=POINT\n", m_gmax[0], m_gmax[1]);
    //fprintf(hfile, "ZONE I=%d, J=%d, F=POINT\n", m_lbuf[0]+m_gmax[0]+m_ubuf[0], m_lbuf[1]+m_gmax[1]+m_ubuf[1]);
    
    int comp[5], max_i[2], max_j[2];
    double max_error[2], approximate, exact, error, tmp;
    double cellCenter[2], beta;
    double normal[2][2] = {{1,0},{0,1}};

    max_i[0] = max_i[1] = max_j[0] = max_j[1] = 0;
    max_error[0] = max_error[1] = 0;
    //for(int j=0; j<m_lbuf[1]+m_gmax[1]+m_ubuf[1]; j++)
    //    for(int i=0; i<m_lbuf[0]+m_gmax[0]+m_ubuf[0]; i++)
    for(int j=m_lbuf[1]; j<m_lbuf[1]+m_gmax[1]; j++)
        for(int i=m_lbuf[0]; i<m_lbuf[0]+m_gmax[0]; i++)
        {    
	     getCellComp(comp,i,j);
	     getCellCenter(cellCenter,i,j);

	     approximate = exact = error = 0;
	     // x derivative
	     for(int k=0; k<4; k++)
	     {
		  if(comp[k]<0)
		       continue;
		  beta = m_pLaplace->getBeta(cellCenter,comp[k]);
		  tmp = beta*getDerivativeAtCellCenter(PLANE_X, comp[k],i,j);
		  tmp = fabs( tmp - m_pLaplace->getFlux(cellCenter,comp[k],normal[0]) );
		  if(tmp>error)
		  {
		       approximate = beta*getDerivativeAtCellCenter(PLANE_X, comp[k],i,j);
		       exact = m_pLaplace->getFlux(cellCenter,comp[k],normal[0]);
		       error = tmp;
		  }
		  if(error>max_error[0])
		  {
		       max_i[0] = i;
		       max_j[0] = j;
		       max_error[0] = error;
		  }
	     }     
	     fprintf(hfile, "%f %f %e %e %e", 
		     cellCenter[0], cellCenter[1],
		     approximate, exact, error);                        

	     // y derivative
	     approximate = exact = error = 0;
	     for(int k=0; k<4; k++)
	     {
		  if(comp[k]<0)
		       continue;
		  beta = m_pLaplace->getBeta(cellCenter,comp[k]);
		  tmp = beta*getDerivativeAtCellCenter(PLANE_Y, comp[k],i,j);
		  tmp = fabs( tmp - m_pLaplace->getFlux(cellCenter,comp[k],normal[1]) );
		  if(tmp>error)
		  {
		       approximate = beta*getDerivativeAtCellCenter(PLANE_Y, comp[k],i,j);
		       exact = m_pLaplace->getFlux(cellCenter,comp[k],normal[1]);
		       error = tmp;
		  }
		  if(error>max_error[1])
		  {
		       max_i[1] = i;
		       max_j[1] = j;
		       max_error[1] = error;
		  }
	     }     
	     fprintf(hfile, " %e %e %e\n", 
		     approximate, exact, error);                        
        }
            
    fclose(hfile);    
    
    printf("ID %d: EBM2D_CARTESIAN::saveDerivatives_Tecplot: \n", m_myid);
    printf("ID %d: \tm_gmax[]={%d,%d}\n", m_myid, m_gmax[0],m_gmax[1]);
    printf("ID %d: \tthe max error for dx is %f in cell (%d,%d)\n",m_myid, max_error[0],max_i[0],max_j[0]);
    printf("ID %d: \tthe max error for dy is %f in cell (%d,%d)\n",m_myid, max_error[1],max_i[1],max_j[1]);
}

void EBM2D_CARTESIAN::saveComponent_Tecplot(void)
{
     char filename[100];
     sprintf(filename, "comps_%d.plt", m_myid);

    RECT_GRID *rect_grid = m_pFront->rect_grid;
    INTERFACE *intfc    = m_pFront->interf;

    FILE *hfile=fopen(filename, "w");
    if(hfile==NULL)
    {
        printf("EBM3D_CARTESIAN::saveComponent_Tecplot: can't open %s \n", filename);
        exit(0);
    }
    
    fprintf(hfile, "VARIABLES = X, Y, VALUES \n");
    fprintf(hfile, "ZONE I=%d, J=%d, F=POINT\n", m_lbuf[0]+m_gmax[0]+m_ubuf[0]+1, m_lbuf[1]+m_gmax[1]+m_ubuf[1]+1);
    
    double comp, coords[2];    
    for(int j=0; j<=m_lbuf[1]+m_gmax[1]+m_ubuf[1]; j++)
    for(int i=0; i<=m_lbuf[0]+m_gmax[0]+m_ubuf[0]; i++)
    {            
        comp = getCellCornerComp(getCellCornerIndex(0,i,j));
	getCellCornerCoords(0, coords, i,j);
        fprintf(hfile, "%f %f %e\n", 
                coords[0], coords[1], comp);                        
    }
            
    fclose(hfile);        
}


void EBM2D_CARTESIAN::debug(int i, int j)
{
    printf("EBM2D_CARTESIAN::debug(%d,%d)\n",i,j);
    double cellCenter[2];
    getCellCenter(cellCenter,i,j);
    printf("\t cellCenter = {%f,%f}", cellCenter[0], cellCenter[1]);
    int index[4];
    getCell(i,j)->getCellUnknownIndex(index);
    printf("\t unknown index = {%d,%d,%d,%d}\n", index[0],index[1],index[2],index[3]);
}

void EBM2D_CARTESIAN::debug(std::vector< std::pair<int,double> > &stencil)
{
    printf("stencil:\n");
    for(int i=0; i<stencil.size(); i++)    
        printf("\t(%d,%f)\n",stencil[i].first, stencil[i].second);    
}

void EBM2D_CARTESIAN::debug_printNeighborComp(int i, int j)
{
    int comp[4];    
    for(int pi=i-1; pi<=i+1; pi++)
        for(int pj=j-1; pj<=j+1; pj++)
	{
	     getCellComp(comp,pi,pj); 
	     debug_printCellComp(comp,pi,pj);
	}    
}

void EBM2D_CARTESIAN::debug_printCellComp(int comp[4], int i, int j)
{
     printf("(%d,%d) has comp = {%d,%d,%d,%d}\n",
	    i,j, comp[0], comp[1], comp[2], comp[3]);
}

void EBM2D_CARTESIAN::debug_printCellUnknownIndex(int i, int j)
{
     int index[4];
     getCell(i,j)->getCellUnknownIndex(index);
     printf("m_cells[%d][%d].m_index = {%d,%d,%d,%d}\n", i,j,
	    index[0], index[1], index[2], index[3]);
}

void EBM2D_CARTESIAN::debug_saveUnknownIndex(void)
{
     char filename[100];
     sprintf(filename, "unknownindex_%d.plt", pp_mynode());
    
    RECT_GRID *rect_grid = m_pFront->rect_grid;
    INTERFACE *intfc    = m_pFront->interf;

    FILE *hfile=fopen(filename, "w");
    if(hfile==NULL)
    {
        printf("EBM2D_CARTESIAN::debug_saveUnknownIndex: can't open %s \n", filename);
        exit(0);
    }
    
    fprintf(hfile, "VARIABLES = X, Y, I, J, Index0, Index1, Index2, Index3 \n");
    fprintf(hfile, "ZONE I=%d, J=%d, F=POINT\n", 
	    m_lbuf[0]+m_gmax[0]+m_ubuf[0], m_lbuf[1]+m_gmax[1]+m_ubuf[1]);

    int index[4];
    double center[2];
    for(int j=0; j<m_lbuf[1]+m_gmax[1]+m_ubuf[1]; j++)
        for(int i=0; i<m_lbuf[0]+m_gmax[0]+m_ubuf[0]; i++)
        {    
	     getCellCenter(center,i,j);
	     getCell(i,j)->getCellUnknownIndex(index);
	     fprintf(hfile, "%e %e %d %d %d %d %d %d\n", 
		     center[0],center[1],
		     i,j,index[0],index[1],index[2],index[3]);
        }
            
    fclose(hfile);    
}

void EBM2D_CARTESIAN::debug_saveCompList(const char*filename)
{
     char full_filename[100];
     sprintf(full_filename, "%s_%d.plt", filename, pp_mynode());
    
    RECT_GRID *rect_grid = m_pFront->rect_grid;
    INTERFACE *intfc    = m_pFront->interf;

    FILE *hfile=fopen(full_filename, "w");
    if(hfile==NULL)
    {
        printf("EBM2D_CARTESIAN::debug_saveCompList: can't open %s \n", full_filename);
        exit(0);
    }
    
    fprintf(hfile, "VARIABLES = X, Y, I, J, c0, c1\n");
    fprintf(hfile, "ZONE I=%d, J=%d, F=POINT\n", 
	    m_lbuf[0]+m_gmax[0]+m_ubuf[0], m_lbuf[1]+m_gmax[1]+m_ubuf[1]);

    int c[2];
    double center[2];
    for(int j=0; j<m_lbuf[1]+m_gmax[1]+m_ubuf[1]; j++)
        for(int i=0; i<m_lbuf[0]+m_gmax[0]+m_ubuf[0]; i++)
        {    
	     getCellCenter(center,i,j);
	     getCellCompList(c,i,j);
	     fprintf(hfile, "%e %e %d %d %d %d\n", 
		     center[0],center[1],
		     i,j,c[0],c[1]);
        }
            
    fclose(hfile);    
}


void EBM2D_CARTESIAN::debug_getMinStatistics(int i, int j)
{
     static boolean bInited = FALSE;
     if(!bInited)
     {
	  bInited = TRUE;
	  m_minCellArea = HUGE_VAL;
	  m_minPartialEdgeLength = HUGE_VAL;
	  m_minIntfcLength = HUGE_VAL;
     }

     int comp[4], size, e;
     double coords[4][2], crx[5][2], beta;   

     getCellCompCoordsCrx(comp,coords,crx,i,j);    
     getCell(i,j)->setCompCoordsCrx(comp,coords,crx);

     // intfc
     std::map<int, double> areas;
     std::map<int, std::vector<EBM2D_EDGE> > mapIntfc;  
     std::map<int, double> mapIntfcLength;
     std::map<int, EBM2D_POINT> mapIntfcNormal;
     std::map<int, EBM2D_POINT> mapIntfcCenter;

     getCell(i,j)->getAreaAndIntfc(areas, mapIntfc);
     getCell(i,j)->getAveragedIntfcLengthNormalCenter(mapIntfc, mapIntfcLength, 
						      mapIntfcNormal, mapIntfcCenter);

     std::map<int, std::vector<EBM2D_EDGE> >::iterator iterMapIntfc;
     EBM2D_EDGE edge;
     EBM2D_POINT point;

     int index = getCell(i,j)->getCellIntfcIndex();    // the first intfc index.
     int intfcIndex[3], intfcComp[2], p;
     double intfcCenter[2], intfcLength, intfcNormal[2], areaFraction[2];

     iterMapIntfc = mapIntfc.begin();
     while(iterMapIntfc!=mapIntfc.end())
     {        
	  //edge = (*iterMapIntfc).second[0];
	  //edge.getIndex(intfcIndex);        
	  //intfcIndex[2] = index++;
	 //edge.getComp(intfcComp);

	 p = (*iterMapIntfc).first;
	 intfcLength = mapIntfcLength[p];
	 if(intfcLength<m_minIntfcLength)
	      m_minIntfcLength = intfcLength;

	 //mapIntfcNormal[p].getCoords(intfcNormal);
	 //mapIntfcCenter[p].getCoords(intfcCenter);

	 iterMapIntfc++;
     }    

     // the two partial areas    
     std::map<int, double>::iterator iterAreas;
     std::map<int, int> unknownToComp;    
     getCell(i,j)->getUnknownToCompMap(unknownToComp);

     int c;
     double vol, rhs, cellCenter[2];
     getCellCenter(cellCenter,i,j);

     iterAreas = areas.begin();
     while(iterAreas!=areas.end())
     {
	  //index = (*iterAreas).first;
	  //c = unknownToComp[index];
	 vol = (*iterAreas).second;

	 if(vol<m_minCellArea)
	      m_minCellArea = vol;

	 //rhs = m_pLaplace->getRightHandSide(cellCenter, c);
	 
	 //areaFraction[0] = getAreaFraction(index, areas);
	 //m_pSolver->Add_b(index, rhs*vol*areaFraction[0]);   

	 iterAreas++;
     }


     // partial edges
     std::map<int, std::vector<EBM2D_PARTIALEDGE> > mapPartialEdges;
     std::map<int, std::vector<EBM2D_PARTIALEDGE> >::iterator iterMapPartialEdges;
     std::vector<EBM2D_PARTIALEDGE>::iterator iterPartialEdge;
     std::vector< std::pair<int,double> > stencil;
     std::vector< std::pair<int,double> >::iterator iterStencil;
     EBM2D_PARTIALEDGE partialedge;

     getCell(i,j)->getPartialEdges(mapPartialEdges);

     iterMapPartialEdges = mapPartialEdges.begin();
     while(iterMapPartialEdges!=mapPartialEdges.end())
     {
	 e = (*iterMapPartialEdges).first;
	 iterPartialEdge = (*iterMapPartialEdges).second.begin();
	 while(iterPartialEdge!=(*iterMapPartialEdges).second.end())
	 {
	     partialedge = (*iterPartialEdge);            
	     if(partialedge.m_length<m_minPartialEdgeLength)
		  m_minPartialEdgeLength = partialedge.m_length;
	     //stencil.clear(); 
	     //if((*iterMapPartialEdges).second.size()==1)
	     //  getFluxStencilForPartialEdge(e, partialedge.m_center, partialedge.m_length, 
	     //			       partialedge.m_comp, stencil, i,j);            
	     //else
	     //  getFluxStencilForPartialEdge2(e, partialedge.m_center, partialedge.m_length, 
	     //				partialedge.m_comp, stencil, i,j);            

	     //beta = m_pLaplace->getBeta(partialedge.m_center, partialedge.m_comp);

	     //areaFraction[0] = getAreaFraction(partialedge.m_index, areas);
	     //iterStencil = stencil.begin();            
	     //while(iterStencil!=stencil.end())
	     //{
	     //	 m_pSolver->Add_A(partialedge.m_index, (*iterStencil).first, 
	     //		  beta*partialedge.m_length * (*iterStencil).second * areaFraction[0]);        
	     // iterStencil++;
	     //}             
	     iterPartialEdge++;
	 }
	 iterMapPartialEdges++;        
     }

}


/******************************************************************************
 * The following is used purely as a scratch paper region for some common codes. 
 
 * for enum of non-buffered cell index:
    
     for(j=m_lbuf[1]; j<m_lbuf[1]+m_gmax[1]; j++)
          for(i=m_lbuf[0]; i<m_lbuf[0]+m_gmax[0]; i++)

 * for enum of buffered cell index:

 *   WEST
     for(j=0; j<m_lbuf[1]+m_gmax[1]+m_ubuf[1]; j++)
          for(i=0; i<m_lbuf[0]; i++)
 *   EAST
     for(j=0; j<m_lbuf[1]+m_gmax[1]+m_ubuf[1]; j++)
          for(i=m_lbuf[0]+m_gmax[0]; i<m_lbuf[0]+m_gmax[0]+m_ubuf[0]; i++)
 *   SOUTH
     for(j=0; j<m_lbuf[1]; j++)
          for(i=0; i<m_lbuf[0]+m_gmax[0]+m_ubuf[0]; i++)
 *   NORTH
     for(j=m_lbuf[1]+m_gmax[1]; j<m_lbuf[1]+m_gmax[1]+m_ubuf[1]; j++)
          for(i=0; i<m_lbuf[0]+m_gmax[0]+m_ubuf[0]; i++)

 * for enum of all cell index:
     
     for(j=0; j<m_lbuf[1]+m_gmax[1]+m_ubuf[1]; j++)
          for(i=0; i<m_lbuf[0]+m_gmax[0]+m_ubuf[0]; i++)

 *****************************************************************************/
 

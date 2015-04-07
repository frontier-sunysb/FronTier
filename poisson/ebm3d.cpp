/***************************************************************
FronTier is a set of libraries that implements different types of 
Front Traking algorithms. Front Tracking is a numerical method for 
the solution of partial differential equations whose solutions have 
discontinuities.  

Copyright (C) 1999 by The University at Stony Brook. 

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
****************************************************************/

/*
 *          ebm3d.c
 * See also ebm3d.h
 */

#include <stdio.h>
#include <stdlib.h>
#include "solver.h"
#include "solver_lapack.h"

#include "ebm3d.h"

#include <FronTier.h>
//#include <FronTier/intfc/geom.h>
//#include <FronTier/intfc/int.h>
//#include "fdecs.h"

#undef assign 
#undef vector 
//#include <ebm3d.h>
#include <math.h>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>

#include <frontier_ppgrid.h>

/*****************************************************************************
 *      EBM3D_LAPLACE
 * The EBM_LAPACE is used to represent the elliptic interface problem. 
 *
 *      Div * beta Grad phi = rou
 *
 *****************************************************************************/
EBM3D_LAPLACE::EBM3D_LAPLACE()
{       
}
// beta
double EBM3D_LAPLACE::getBeta(double coords[3], int comp)
{
    // return 1;         // for debug.
    if(comp==1)
        return 1;
    else if(comp==2)
	 return 10;
    else
        printf("EBM_LAPLACE::getBeta: unknown case comp=%d", comp);
}
// phi
double EBM3D_LAPLACE::getExactSolution(double coords[3], int comp)
{    
    double r = sqrt(coords[0]*coords[0]+coords[1]*coords[1]+coords[2]*coords[2]);    
    //return coords[0] + coords[1] + coords[2];         // for debug.

    if(comp==1)
	 return r*r*r;
    else if(comp==2)
	 return r*r*r/10 + (1-1.0/10)*(1.0/8);        
    else
	 printf("EBM3D_LAPLACE::getExactSolution: unknown case comp=%d\n", comp);
}

double EBM3D_LAPLACE::getFlux(double coords[3], int comp, double normal[3])
{
     //return normal[0] + normal[1] + normal[2];              // for debug.
     
     double r = sqrt(coords[0]*coords[0]+coords[1]*coords[1]+coords[2]*coords[2]);    
     return 3*coords[0]*r*normal[0]+3*coords[1]*r*normal[1]+3*coords[2]*r*normal[2];
}

// rou    
double EBM3D_LAPLACE::getRightHandSide(double coords[3], int comp)
{
     //return 0;               // for debug.
     double r = sqrt(coords[0]*coords[0]+coords[1]*coords[1]+coords[2]*coords[2]); 
     return 12*r;
}

/*****************************************************************************
 *          EBM3D_INDEX
 * a class for holding 3 integers.
 *****************************************************************************/
void EBM3D_INDEX::getIndex(int index[3])
{
     index[0] = m_index[0];
     index[1] = m_index[1];
     index[2] = m_index[2];
}
void EBM3D_INDEX::setIndex(int index[3])
{
     m_index[0] = index[0];
     m_index[1] = index[1];
     m_index[2] = index[2];
}

/******************************************************************************
 *          EBM3D_POINT
 * simple wrap for 3D point. 
 *****************************************************************************/
EBM3D_POINT::EBM3D_POINT()
{
}
EBM3D_POINT::EBM3D_POINT(double coords[3])
{
    setCoords(coords);
}
void EBM3D_POINT::getCoords(double coords[3])
{
    coords[0] = m_coords[0];
    coords[1] = m_coords[1];
    coords[2] = m_coords[2];
}
void EBM3D_POINT::setCoords(double coords[3])
{
    m_coords[0] = coords[0];
    m_coords[1] = coords[1];
    m_coords[2] = coords[2];
}

/******************************************************************************
 *      EBM3D_TRIANGLE
 ******************************************************************************/

/*
 * is there problem with the following two functions when used with setMatrixForBoundary()?????
 */
void EBM3D_TRIANGLE::getComp(int comp[2])
{
     comp[0] = m_comp[0];
     comp[1] = m_comp[1];
}
void EBM3D_TRIANGLE::getIndex(int index[2])
{
     index[0] = m_index[0];
     index[1] = m_index[1];
}
// always point from a smaller component to a bigger component!
// return the smaller component.
int EBM3D_TRIANGLE::setTriangle(int index0, int comp0, int index1, int comp1, double v0[3], double v1[3], double v2[3])
{
    if(comp0<comp1)
    {
        m_comp[0] = comp0;
        m_comp[1] = comp1;
        m_index[0] = index0;
        m_index[1] = index1;

        for(int i=0; i<3; i++)    
        {
            m_vertex[0][i] = v0[i];
            m_vertex[1][i] = v1[i];
            m_vertex[2][i] = v2[i];
        }
    }
    else    
    {
        m_comp[0] = comp1;
        m_comp[1] = comp0;
        m_index[0] = index1;
        m_index[1] = index0;

        for(int i=0; i<3; i++)    
        {
            m_vertex[0][i] = v0[i];
            m_vertex[1][i] = v2[i];         
            m_vertex[2][i] = v1[i];
        } 
    }
    return m_comp[0];   // m_comp[0] is always the smaller component.
}
/*
 * used for debugging with EBM3D_PARTIALFACE.
 */
int EBM3D_TRIANGLE::setTriangle(double v0[3], double v1[3], double v2[3])
{
     for(int i=0; i<3; i++)    
     {
	  m_vertex[0][i] = v0[i];
	  m_vertex[1][i] = v1[i];           // changed from v2
	  m_vertex[2][i] = v2[i];           // changed from v1
     }
}

double EBM3D_TRIANGLE::getArea(void)
{
    double a[3], b[3], normal[3];
    difference(m_vertex[1],m_vertex[0],a,3);
    difference(m_vertex[2],m_vertex[0],b,3);   
    
    Cross3d(a,b,normal);
    double r = Mag3d(normal);
    return r/2;    
}
void EBM3D_TRIANGLE::getNormal(double normal[3])
{
    double a[3], b[3];
    difference(m_vertex[1],m_vertex[0],a,3);
    difference(m_vertex[2],m_vertex[0],b,3);
    
    Cross3d(a,b,normal);    
    double r = Mag3d(normal);
    
    normal[0] /= r;
    normal[1] /= r;
    normal[2] /= r;    
}

void EBM3D_TRIANGLE::getAreaNormalCenter(double &area, double normal[3], double center[3])     // ??? center ???
{
    double a[3], b[3];
    difference(m_vertex[1],m_vertex[0],a,3);
    difference(m_vertex[2],m_vertex[0],b,3);
    
    Cross3d(a,b,normal);    
    double r = Mag3d(normal);
    
    area = r/2;
    
    normal[0] /= r;
    normal[1] /= r;
    normal[2] /= r;     
    
    // center: what kind of center is needed?
    center[0] = (m_vertex[0][0] + m_vertex[1][0])/2;
    center[1] = (m_vertex[0][1] + m_vertex[1][1])/2;
    center[2] = (m_vertex[0][2] + m_vertex[1][2])/2;
    
    center[0] = 2.0/3*center[0] + 1.0/3*m_vertex[2][0];
    center[1] = 2.0/3*center[1] + 1.0/3*m_vertex[2][1];
    center[2] = 2.0/3*center[2] + 1.0/3*m_vertex[2][2];        
}
/******************************************************************************
 *          EBM3D_INTFC
 ******************************************************************************/
/*
void EBM3D_INTFC::begin(EBM3D_TRIANGLE & triangle)
{
    triangle.getComp(m_comp);
    triangle.getIndex(m_index);
    triangle.getAreaNormalCenter(m_area, m_normal, m_center);
    for(int i=0; i<3; i++)
    {
        m_normal[i] *= m_area;
        m_center[i] *= m_area;
    }
}
void EBM3D_INTFC::addTri(EBM3D_TRIANGLE &triangle)
{
    double area, normal[3], center[3];
    triangle.getAreaNormalCenter(area, normal, center);
    m_area += area;
    for(int i=0; i<3; i++)
    {
        m_normal[i] += area*normal[i];
        m_center[i] += area*center[i];
    }
}

void EBM3D_INTFC::end(void)
{
    double tmp = Mag3d(m_normal);
    for(int i=0; i<3; i++)
    {
        m_normal[i] /= tmp;
        m_center[i] /= m_area;
    }    
}
void EBM3D_INTFC::getComp(int comp[2])
{
    comp[0] = m_comp[0];
    comp[1] = m_comp[1];
}
void EBM3D_INTFC::getIndex(int index[2])
{
    index[0] = m_index[0];
    index[1] = m_index[1];
}
double EBM3D_INTFC::getArea(void)
{
    return m_area;
} 
void EBM3D_INTFC::getNormal(double normal[3])
{
    normal[0] = m_normal[0];
    normal[1] = m_normal[1];
    normal[2] = m_normal[2];
}
void EBM3D_INTFC::getAreaNormalCenter(double &area, double normal[3], double center[3])
{
    area = m_area;
    normal[0] = m_normal[0];
    normal[1] = m_normal[1];
    normal[2] = m_normal[2];
    center[0] = m_center[0];
    center[1] = m_center[1];
    center[2] = m_center[2];    
}
*/

/******************************************************************************
 *          EBM3D_PARTIALFACE
 ******************************************************************************/
// first function to call in a sequence of steps to get the partial face center.
void EBM3D_PARTIALFACE::begin(int c, int index)
{
    m_comp = c;
    m_index = index;
    m_area = 0;
    m_center[0] = 0;
    m_center[1] = 0;
    m_center[2] = 0;

    debug_begin();
}
/*void EBM3D_PARTIALFACE::begin(int c, int index, double area, double center[3])
{
    m_comp = c;
    m_index = index;
    m_area = area;
    m_center[0] = center[0] * area;
    m_center[1] = center[1] * area;
    m_center[2] = center[2] * area;    
    }*/
void EBM3D_PARTIALFACE::addTri(double v0[3], double v1[3], double v2[3])
{
    double tmp_area;
    double center[3];
    tmp_area = getTriArea(v0,v1,v2);
    getTriCenter(v0,v1,v2,center);
    
    m_area += tmp_area;
    for(int i=0; i<3; i++)
        m_center[i] += tmp_area*center[i];
    
    debug_addTri(v0,v1,v2);
}
void EBM3D_PARTIALFACE::addQuad(double v0[3], double v1[3], double v2[3], double v3[3])
{
    addTri(v0,v1,v2);
    addTri(v0,v2,v3);
}
void EBM3D_PARTIALFACE::end(void)
{
    for(int i=0; i<3; i++)
        m_center[i] /= m_area;
}

void EBM3D_PARTIALFACE::debug_begin(void)
{
     m_triangles.clear();
}
void EBM3D_PARTIALFACE::debug_addTri(double v0[3], double v1[3], double v2[3])
{
     EBM3D_TRIANGLE tri;
     tri.setTriangle(v0,v1,v2);
     m_triangles.push_back(tri);
}

double EBM3D_PARTIALFACE::getTriArea(double v0[3], double v1[3], double v2[3])
{
    double a[3], b[3], normal[3];
    difference(v1,v0,a,3);
    difference(v2,v0,b,3);
    
    Cross3d(a,b,normal);
    double r = Mag3d(normal);
    return r/2;    
}
void EBM3D_PARTIALFACE::getTriCenter(double v0[3], double v1[3], double v2[3], double center[3])
{
    int i;
    for(i=0; i<3; i++)
        center[i] = 1.0/3*(v0[i]+v1[i]+v2[i]);
}


/*******************************************************************************
 *          EBM3D_CELL 
 ******************************************************************************/
/*
 * Return the unknown indices in index[8] and also find the maximum of them.
 * It returns the maximun of all indices.
 */
EBM3D_CELL::EBM3D_CELL()
{
    init();
}

void EBM3D_CELL::init(void)
{
    m_nCellUnknowns = -1;
    for(int i=0; i<8; i++)
        m_index[i] = -1;
}

CELL_TYPE EBM3D_CELL::getCellType(void)
{
     return m_cellType;
}

/* there are at most two different components in comp[8] for now.
 * calculate the min & max of comp[8].
 *        EXTERNAL  if max<0;
 *        BOUNDARY if min<0;
 *        PARTIAL    if min!=max;
 *        INTERNAL otherwise.
 */ 
void EBM3D_CELL::setCellType(int *comp)
{
     int min, max;
     min = max = comp[0];
     for(int i=1; i<8; i++)
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
 * Return the unknown indices in index[8] and return the maximum of them.
 */
int EBM3D_CELL::getCellUnknownIndex(int index[8])
{
    int maximun = index[0] = m_index[0];    
    for(int i=1;i<8;i++)
    {
        index[i] = m_index[i];
        if(maximun<m_index[i])
            maximun = m_index[i];
    }
    return maximun;
}

void EBM3D_CELL::setCellUnknownIndex(int index[8])
{
     for(int i=0; i<8; i++)
	  m_index[i] = index[i];
}

int EBM3D_CELL::getCellUnknownIndex(int corner)
{
    return m_index[corner];
}

// return the unknown numbers in the current cell.
// These are for the unknowns stored at the cell center only.
int EBM3D_CELL::getCellUnknownNumber(void)
{
    return m_nCellUnknowns;
}

// the first intfc index
int EBM3D_CELL::getCellIntfcIndex(void)
{
     //return m_index[0] + m_nCellUnknowns;
     int max = m_index[0];
     for(int i=1; i<8; i++)
	  if(m_index[i]>max)
	       max = m_index[i];
     return max + 1;
}

// the following functions are used for the generalized marching cubes 
// method to get the number of unkowns, the volume fractions, the area of 
// the surfaces and their normal.
void EBM3D_CELL::setComp(int comp[8])
{
    m_comp = comp;
}

// Comp[8]:      the components of the vertices;
// coords[8][3]: the coordinates of the vertices;
// crx[19][3]: the crossing on the 19 edges/diagonals.
void EBM3D_CELL::setCompCoordsCrx(int comp[8], double coords[8][3], double crx[19][3])
{
    m_comp = comp;
    m_coords = coords;
    m_crx = crx;   
}

/*
 * this function should only be called after EBM3D_CELL::setCompCoordsCrx() is called.
 */
void EBM3D_CELL::getCellFaceCompCoords(int face, int comp[4], double coords[4][3])
{
     switch(face)
     {
     case FACE_WEST:
	  getComp(comp,0,3,7,4);
	  getCoords(coords,0,3,7,4);
	  break;
     case FACE_EAST:
	  getComp(comp,1,2,6,5);
	  getCoords(coords,1,2,6,5);
	  break;
     case FACE_SOUTH:
	  getComp(comp,0,1,5,4);
	  getCoords(coords,0,1,5,4);
	  break;
     case FACE_NORTH:
	  getComp(comp,3,2,6,7);
	  getCoords(coords,3,2,6,7);
	  break;
     case FACE_UP:
	  getComp(comp,4,5,6,7);
	  getCoords(coords,4,5,6,7);
	  break;
     case FACE_DOWN:
	  getComp(comp,0,1,2,3);
	  getCoords(coords,0,1,2,3);
	  break;
     }
}

/*
 * find the nearest corner among the 4 corners of the cell face, with the same components.
 * return 0,1,2,3      if successful,
 *        -1           if not.
 */ 
int EBM3D_CELL::getCellFaceNearestCorner(int face, int comp, double p[3])
{
     int Comps[4];
     double Coords[4][3];

     getCellFaceCompCoords(face,Comps,Coords);

     // find first 
     int i, index;
     for(index=0; index<4; index++)
	  if(comp==Comps[index])
	       break;
     // is there a valid corner?
     if(index>=4)
	  return -1;

     // the first valid corner is index
     double dist[4], max = HUGE_VAL;
     for(i=index; i<4; i++)
	  if(comp==Comps[i])
	  {
	       dist[i] = getDistance2(Coords[i],p);
	       if(dist[i]<max)
	       {
		    max = dist[i];
		    index = i;
	       }
	  }
     return index;
}

// use marching tetrahedra method. There are 6 tetrahedra in each cube.
int EBM3D_CELL::setCellUnknownIndex(int start_index)
{
    int i, index, end_index;
    index = end_index = start_index;
    for(i=0; i<8; i++)
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

void EBM3D_CELL::setCellUnknownIndex(int corner, int index)
{    
    int pa, pb, j;
    for(j=0; j<7; j++)
    {
        pa = corner;
        pb = m_cubeEdges[pa][j];
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

void EBM3D_CELL::getUnknownToCompMap(std::map<int, int> &unknownToComp)
{
    for(int i=0; i<8; i++)    
        unknownToComp[m_index[i]] = m_comp[i];    
} 

/*
 *  get the volumes and the intfcs inside the current grid block.
 *  to be tested,
 */
void EBM3D_CELL::getVolumeAndIntfc(std::map<int, double> &volumes,                        
                        std::map<int, std::vector<EBM3D_TRIANGLE> > &intfc)
{
    int vertices[4], edges[6];
    
    // tetra 0
    setTetraVertices(vertices, 0,6,4,5);
    setTetraEdges(edges, 18,17,8,12,5,4);    
    getTetraVolumeAndIntfc(vertices,edges,volumes,intfc);
    // tetra 1
    setTetraVertices(vertices, 0,1,6,5);
    setTetraEdges(edges, 0,13,18,12,9,5);
    getTetraVolumeAndIntfc(vertices,edges,volumes,intfc);
    // tetra 2
    setTetraVertices(vertices, 0,1,2,6);
    setTetraEdges(edges, 0,1,16,18,13,10);
    getTetraVolumeAndIntfc(vertices,edges,volumes,intfc);
    // tetra 3
    setTetraVertices(vertices, 0,2,3,6);
    setTetraEdges(edges, 16,2,3,18,10,14);
    getTetraVolumeAndIntfc(vertices,edges,volumes,intfc);
    // tetra 4
    setTetraVertices(vertices, 0,6,3,7);
    setTetraEdges(edges, 18,14,3,15,6,11);
    getTetraVolumeAndIntfc(vertices,edges,volumes,intfc);
    // tetra 5
    setTetraVertices(vertices, 0,6,7,4);
    setTetraEdges(edges, 18,6,15,8,17,7);        
    getTetraVolumeAndIntfc(vertices,edges,volumes,intfc);
    
    //test_checkVolumes(volumes);
}


void EBM3D_CELL::getTetraVolumeAndIntfc(int vertices[4], int edges[6],
                        std::map<int, double> &volumes,
                        std::map<int, std::vector<EBM3D_TRIANGLE> > &intfc)
{
    int c0 = m_index[vertices[0]];
    int c1 = m_index[vertices[1]];
    int c2 = m_index[vertices[2]];
    int c3 = m_index[vertices[3]];    
    int min_comp;
    
    EBM3D_TRIANGLE triangle;
    
    double tetraVolume = getTetraVolume(m_coords[vertices[0]], m_coords[vertices[1]], m_coords[vertices[2]], m_coords[vertices[3]]);
    
    // 8 cases    
    if(c0==c1&&c0==c2&&c0==c3)   // 0000
    {        
        volumes[c0] += tetraVolume;
        
    }
    else if(c0!=c1&&c0!=c2&&c0!=c3) // 0111
    {
        double v = getTetraVolume(m_coords[vertices[0]], m_crx[edges[0]], m_crx[edges[2]], m_crx[edges[3]]);        
        volumes[c0] += v;
        volumes[c1] += tetraVolume - v;
        min_comp = triangle.setTriangle(c0,m_comp[vertices[0]],c1,m_comp[vertices[1]],m_crx[edges[0]], m_crx[edges[2]], m_crx[edges[3]]);
        intfc[min_comp].push_back(triangle);
    }
    else if(c1!=c0&&c1!=c2&&c1!=c3) // 1011
    {
        double v = getTetraVolume(m_coords[vertices[1]], m_crx[edges[1]], m_crx[edges[0]], m_crx[edges[4]]);       
        volumes[c1] += v;
        volumes[c0] += tetraVolume - v;
        min_comp = triangle.setTriangle(c1,m_comp[vertices[1]],c0,m_comp[vertices[0]],m_crx[edges[1]], m_crx[edges[0]], m_crx[edges[4]]);
        intfc[min_comp].push_back(triangle);
    }
    else if(c2!=c0&&c2!=c1&&c2!=c3) // 1101
    {
        double v = getTetraVolume(m_coords[vertices[2]], m_crx[edges[2]], m_crx[edges[1]], m_crx[edges[5]]);       
        volumes[c2] += v;
        volumes[c0] += tetraVolume - v;
        min_comp = triangle.setTriangle(c2,m_comp[vertices[2]],c0,m_comp[vertices[0]],m_crx[edges[2]], m_crx[edges[1]], m_crx[edges[5]]);
        intfc[min_comp].push_back(triangle);
    }
    else if(c3!=c0&&c3!=c1&&c3!=c2) // 1110
    {
        double v = getTetraVolume(m_coords[vertices[3]], m_crx[edges[4]], m_crx[edges[3]], m_crx[edges[5]]);
        volumes[c3] += v;
        volumes[c0] += tetraVolume - v;
        min_comp = triangle.setTriangle(c3,m_comp[vertices[3]],c0,m_comp[vertices[0]],m_crx[edges[4]], m_crx[edges[3]], m_crx[edges[5]]);
        intfc[min_comp].push_back(triangle);
    }
    else if(c0==c1&&c0!=c2&&c2==c3) // 0011
    {
        double v = getTriangularPrismVolume(m_coords[vertices[0]], m_crx[edges[2]], m_crx[edges[3]],
                                            m_coords[vertices[1]], m_crx[edges[1]], m_crx[edges[4]]);
        volumes[c0] += v;
        volumes[c3] += tetraVolume - v;
        min_comp = triangle.setTriangle(c0,m_comp[vertices[0]],c3,m_comp[vertices[3]],m_crx[edges[2]], m_crx[edges[3]], m_crx[edges[4]]);
        intfc[min_comp].push_back(triangle);
        min_comp = triangle.setTriangle(c0,m_comp[vertices[0]],c3,m_comp[vertices[3]],m_crx[edges[1]], m_crx[edges[2]], m_crx[edges[4]]);
        intfc[min_comp].push_back(triangle);
    }
    else if(c0==c2&&c0!=c1&&c1==c3) // 0101
    {
        double v = getTriangularPrismVolume(m_coords[vertices[0]], m_crx[edges[3]], m_crx[edges[0]],
                                            m_coords[vertices[2]], m_crx[edges[5]], m_crx[edges[1]]);
        volumes[c0] += v;
        volumes[c1] += tetraVolume - v;
        min_comp = triangle.setTriangle(c0,m_comp[vertices[0]],c1,m_comp[vertices[1]],m_crx[edges[0]], m_crx[edges[1]], m_crx[edges[3]]);
        intfc[min_comp].push_back(triangle);
        min_comp = triangle.setTriangle(c0,m_comp[vertices[0]],c1,m_comp[vertices[1]],m_crx[edges[1]], m_crx[edges[5]], m_crx[edges[3]]);
        intfc[min_comp].push_back(triangle);
    }
    else if(c0==c3&&c0!=c1&&c1==c2) // 0110
    {
        double v = getTriangularPrismVolume(m_coords[vertices[0]], m_crx[edges[0]], m_crx[edges[2]],
                                            m_coords[vertices[3]], m_crx[edges[4]], m_crx[edges[5]]);
        volumes[c0] += v;
        volumes[c1] += tetraVolume - v;
        min_comp = triangle.setTriangle(c0,m_comp[vertices[0]],c1,m_comp[vertices[1]],m_crx[edges[0]], m_crx[edges[2]], m_crx[edges[4]]);
        intfc[min_comp].push_back(triangle);
        min_comp = triangle.setTriangle(c0,m_comp[vertices[0]],c1,m_comp[vertices[1]],m_crx[edges[2]], m_crx[edges[5]], m_crx[edges[4]]);
        intfc[min_comp].push_back(triangle);
    }
    else
    {
        printf("EBM3D_CELL::getTetraVolumeIntfc: something is wrong!\n");
    }
}


// this function needs to be modified and tested.
// to be tested,
void EBM3D_CELL::getAveragedIntfcAreaNormalCenter(std::map<int, std::vector<EBM3D_TRIANGLE> > &intfc,
                        std::map<int, double> &intfcArea,
                        std::map<int, EBM3D_POINT> &intfcNormal,
                        std::map<int, EBM3D_POINT> &intfcCenter)
{
    int i;
    double area, normal[3], center[3], tmp;
    double sumArea, sumNormal[3], sumCenter[3];
    EBM3D_POINT point;
    std::map<int, std::vector<EBM3D_TRIANGLE> >::iterator iter;
    std::vector<EBM3D_TRIANGLE>::iterator iterTri;
    
    iter = intfc.begin();
    while(iter!=intfc.end())
    {        
        sumArea = 0;
        sumNormal[0] = sumNormal[1] = sumNormal[2] = 0;
        sumCenter[0] = sumCenter[1] = sumCenter[2] = 0;
        iterTri = (*iter).second.begin();
        while(iterTri!=(*iter).second.end())
        {
            (*iterTri).getAreaNormalCenter(area,normal,center);
            
            sumArea += area;
            for(int j=0; j<3; j++)            
            {
                sumNormal[j] += normal[j] * area;
                sumCenter[j] += center[j] * area;
            }
            
            iterTri++;
        }   
     
        i = (*iter).first;
        
        tmp = Mag3d(sumNormal);
        for(int j=0; j<3; j++)
        {
            sumNormal[j] /= tmp;        // ??? should this be sumArea ???
            sumCenter[j] /= sumArea;
        }
        
        intfcArea[i] = sumArea;         // wrongly used "area" before!
        
        
        point.setCoords(sumNormal);
        intfcNormal[i] = EBM3D_POINT(point);
        point.setCoords(sumCenter);
        intfcCenter[i] = EBM3D_POINT(point);
        
        iter++;
    }
}

/*
 * determine the partial faces of a cell boundary.
 * do we need to consider the orentation of the face?
 */
void EBM3D_CELL::getPartialFaces(std::map<int, std::vector<EBM3D_PARTIALFACE> > &partialfaces)
{
    int corner[4], edge[5];
    std::vector<EBM3D_PARTIALFACE> faces;
    partialfaces.clear();
    // south
    setFaceVertices(corner, 0,1,5,4);
    setFaceEdges(edge, 0,9,4,8,12);
    faces.clear();
    getPartialFaces(corner,edge,faces);
    partialfaces[FACE_SOUTH] = faces;
    // east
    setFaceVertices(corner, 1,2,6,5);
    setFaceEdges(edge, 1,10,5,9,13);
    faces.clear();
    getPartialFaces(corner,edge,faces);
    partialfaces[FACE_EAST] = faces;
    // north
    setFaceVertices(corner, 3,2,6,7);
    setFaceEdges(edge, 2,10,6,11,14);
    faces.clear();
    getPartialFaces(corner,edge,faces);
    partialfaces[FACE_NORTH] = faces;
    // west
    setFaceVertices(corner, 0,3,7,4);
    setFaceEdges(edge, 3,11,7,8,15);
    faces.clear();
    getPartialFaces(corner,edge,faces);
    partialfaces[FACE_WEST] = faces;
    // down
    setFaceVertices(corner, 0,1,2,3);
    setFaceEdges(edge, 0,1,2,3,16);    
    faces.clear();
    getPartialFaces(corner,edge,faces);
    partialfaces[FACE_DOWN] = faces;
    // up
    setFaceVertices(corner, 4,5,6,7);
    setFaceEdges(edge, 4,5,6,7,17);
    faces.clear();
    getPartialFaces(corner,edge,faces);
    partialfaces[FACE_UP] = faces;
    
    //test_getPartialFaces(partialfaces);
}

/*
 * the ordering of the corner[] is the following
 *      32
 *      01
 * the ordering of the edge[] is the following
 *      2
 *     341
 *      0
 * This is used with the generalized marching cubes method. Therefore, there is
 * always the diagonal connecting corner[0] and corner[2].
 * seems to be fine.
 */
void EBM3D_CELL::getPartialFaces(int corner[4], int edge[5], std::vector<EBM3D_PARTIALFACE> &partialfaces)
{   
    int comp[] = {m_comp[corner[0]], m_comp[corner[1]], m_comp[corner[2]], m_comp[corner[3]]};
    double dx = getDistance(m_coords[corner[0]], m_coords[corner[1]]);
    double dy = getDistance(m_coords[corner[1]], m_coords[corner[2]]);
    
    EBM3D_PARTIALFACE face;
    double area, cubeArea, center[3], center0[3], center1[3];
    
    cubeArea = dx*dy;
    
    if(comp[0]==comp[1]&&comp[0]==comp[2]&&comp[0]==comp[3])    // 0000
    {
        getCenter(m_coords[corner[0]],m_coords[corner[2]],center);
        //face.begin(comp[0], m_index[corner[0]], cubeArea, center);             
	face.begin(comp[0], m_index[corner[0]]);
	face.addQuad(m_coords[corner[0]], m_coords[corner[1]], m_coords[corner[2]], m_coords[corner[3]]);
        face.end();        
        partialfaces.push_back(face);
    }    
    else if(comp[1]==comp[2]&&comp[1]==comp[3])                 // 0111
    {
        //area  = getTriArea(corner[0], edge[0], edge[4]);
        //area += getTriArea(corner[0], edge[4], edge[3]);
        //getCenter(m_coords[corner[0]], m_crx[edge[0]], m_crx[edge[4]], m_crx[edge[3]], center);
        //face.setPartialFace(comp[0], m_index[corner[0]], area, center);                 
        face.begin(comp[0], m_index[corner[0]]);
        face.addQuad(m_coords[corner[0]], m_crx[edge[0]], m_crx[edge[4]], m_crx[edge[3]]);
        face.end();
        partialfaces.push_back(face);
        //getCenter(m_coords[corner[2]], m_coords[corner[3]], m_crx[edge[3]], m_crx[edge[4]], center);
        //getCenter(m_coords[corner[1]], m_coords[corner[2]], m_crx[edge[4]], m_crx[edge[0]], center0);
        //getCenter(center, center0, center);
        //face.setPartialFace(comp[1], m_index[corner[1]], cubeArea - area, center);
        face.begin(comp[1], m_index[corner[1]]);
        face.addQuad(m_coords[corner[2]], m_coords[corner[3]], m_crx[edge[3]], m_crx[edge[4]]);
        face.addQuad(m_coords[corner[1]], m_coords[corner[2]], m_crx[edge[4]], m_crx[edge[0]]);
        face.end();
        partialfaces.push_back(face);
    }    
    else if(comp[0]==comp[2]&&comp[0]==comp[3])                 // 0100
    {
        //area  = getTriArea(corner[1], edge[1], edge[0]);        
        //getCenter(m_coords[corner[1]], m_crx[edge[1]], m_crx[edge[0]], center);
        //face.setPartialFace(comp[1], m_index[corner[1]], area, center);        
        face.begin(comp[1], m_index[corner[1]]);
        face.addTri(m_coords[corner[1]], m_crx[edge[1]], m_crx[edge[0]]);
        face.end();
        partialfaces.push_back(face);
        //getCenter(m_coords[corner[0]], m_coords[corner[2]], m_coords[corner[3]], center);
        //getCenter(m_coords[corner[0]], m_crx[edge[0]], m_crx[edge[1]], m_coords[corner[2]], center0);
        //getCenter(center, center0, center);
        //face.setPartialFace(comp[0], m_index[corner[0]], cubeArea - area, center);
        face.begin(comp[0], m_index[corner[0]]);
        face.addTri(m_coords[corner[0]], m_coords[corner[2]], m_coords[corner[3]]);
        face.addQuad(m_coords[corner[0]], m_crx[edge[0]], m_crx[edge[1]], m_coords[corner[2]]);
        face.end();
        partialfaces.push_back(face);
    }
    else if(comp[0]==comp[1]&&comp[0]==comp[3])                 // 0010
    {
        //area  = getTriArea(corner[2], edge[2], edge[4]);
        //area += getTriArea(corner[2], edge[4], edge[1]);
        //getCenter(m_coords[corner[2]], m_crx[edge[2]], m_crx[edge[4]], m_crx[edge[1]], center);
        //face.setPartialFace(comp[2], m_index[corner[2]], area, center);
        face.begin(comp[2], m_index[corner[2]]);
        face.addQuad(m_coords[corner[2]], m_crx[edge[2]], m_crx[edge[4]], m_crx[edge[1]]);
        face.end();
        partialfaces.push_back(face);
        //getCenter(m_coords[corner[3]], m_coords[corner[0]], m_crx[edge[4]], m_crx[edge[2]], center);
        //getCenter(m_coords[corner[0]], m_coords[corner[1]], m_crx[edge[1]], m_crx[edge[4]], center0);
        //getCenter(center, center0, center);
        //face.setPartialFace(comp[0], m_index[corner[0]], cubeArea - area, center);
        face.begin(comp[0], m_index[corner[0]]);
        face.addQuad(m_coords[corner[3]], m_coords[corner[0]], m_crx[edge[4]], m_crx[edge[2]]);
        face.addQuad(m_coords[corner[0]], m_coords[corner[1]], m_crx[edge[1]], m_crx[edge[4]]);
        face.end();
        partialfaces.push_back(face);
    }
    else if(comp[0]==comp[1]&&comp[0]==comp[2])                 // 0001
    {
        //area  = getTriArea(corner[3], edge[3], edge[2]);       
        //getCenter(m_coords[corner[3]], m_crx[edge[3]], m_crx[edge[2]], center);
        //face.setPartialFace(comp[3], m_index[corner[3]], area, center);        
        face.begin(comp[3], m_index[corner[3]]);
        face.addTri(m_coords[corner[3]], m_crx[edge[3]], m_crx[edge[2]]);
        face.end();
        partialfaces.push_back(face);        
        //getCenter(m_coords[corner[0]], m_coords[corner[1]], m_coords[corner[2]], center);
        //getCenter(m_coords[corner[2]], m_crx[edge[2]], m_crx[edge[3]], m_coords[corner[0]], center0);
        //getCenter(center, center0, center);
        //face.setPartialFace(comp[0], m_index[corner[0]], cubeArea - area, center);
        face.begin(comp[0], m_index[corner[0]]);
        face.addTri(m_coords[corner[0]], m_coords[corner[1]], m_coords[corner[2]]);
        face.addQuad(m_coords[corner[2]], m_crx[edge[2]], m_crx[edge[3]], m_coords[corner[0]]);
        face.end();
        partialfaces.push_back(face);
    }
    else if(comp[0]==comp[3]&&comp[1]==comp[2])                 // 0110
    {
        //area  = getTriArea(corner[0], edge[0], edge[4]);
        //area += 0.5*cubeArea - getTriArea(corner[2], edge[2], edge[4]);
        //getCenter(m_coords[corner[0]], m_crx[edge[0]], m_crx[edge[4]], center);
        //getCenter(m_coords[corner[0]], m_crx[edge[4]], m_crx[edge[2]], m_coords[corner[3]], center0);
        //getCenter(center, center0, center);
        //face.setPartialFace(comp[0], m_index[corner[0]], area, center);
        face.begin(comp[0], m_index[corner[0]]);
        face.addTri(m_coords[corner[0]], m_crx[edge[0]], m_crx[edge[4]]);
        face.addQuad(m_coords[corner[0]], m_crx[edge[4]], m_crx[edge[2]], m_coords[corner[3]]);
        face.end();
        partialfaces.push_back(face);
        //getCenter(m_coords[corner[2]], m_crx[edge[2]], m_crx[edge[4]], center);
        //getCenter(m_coords[corner[2]], m_crx[edge[4]], m_crx[edge[0]], m_coords[corner[1]], center0);
        //getCenter(center, center0, center);
        //face.setPartialFace(comp[1], m_index[corner[1]], cubeArea - area, center);
        face.begin(comp[1], m_index[corner[1]]);
        face.addTri(m_coords[corner[2]], m_crx[edge[2]], m_crx[edge[4]]);
        face.addQuad(m_coords[corner[2]], m_crx[edge[4]], m_crx[edge[0]], m_coords[corner[1]]);
        face.end();
        partialfaces.push_back(face);
    }
    else if(comp[0]==comp[1]&&comp[2]==comp[3])                 // 0011
    {
        //area  = getTriArea(corner[0], edge[4], edge[3]);
        //area += 0.5*cubeArea - getTriArea(corner[2], edge[4], edge[1]);
        //getCenter(m_coords[corner[0]], m_crx[edge[4]], m_crx[edge[3]], center);
        //getCenter(m_coords[corner[1]], m_crx[edge[1]], m_crx[edge[4]], m_coords[corner[0]], center0);
        //getCenter(center, center0, center);
        //face.setPartialFace(comp[0], m_index[corner[0]], area, center);
        face.begin(comp[0], m_index[corner[0]]);
        face.addTri(m_coords[corner[0]], m_crx[edge[4]], m_crx[edge[3]]);        
        face.addQuad(m_coords[corner[1]], m_crx[edge[1]], m_crx[edge[4]], m_coords[corner[0]]);
        face.end();
        partialfaces.push_back(face);        
        //getCenter(m_coords[corner[2]], m_crx[edge[4]], m_crx[edge[1]], center);
        //getCenter(m_coords[corner[3]], m_crx[edge[3]], m_crx[edge[4]], m_coords[corner[2]], center0);
        //getCenter(center, center0, center);
        //face.setPartialFace(comp[3], m_index[corner[3]], cubeArea - area, center);
        face.begin(comp[3], m_index[corner[3]]);
        face.addTri(m_coords[corner[2]], m_crx[edge[4]], m_crx[edge[1]]);
        face.addQuad(m_coords[corner[3]], m_crx[edge[3]], m_crx[edge[4]], m_coords[corner[2]]);
        face.end();
        partialfaces.push_back(face);
    }
    // two more cases not delt with
}


void EBM3D_CELL::resetCellUnknownIndex(int ilower)
{
     int i;
     for(i=0; i<8; i++)
	  if(m_index[i]>=0)
	       m_index[i] += ilower;
}

void EBM3D_CELL::test(void)
{
    int comp[][8] = {
                     {2,2,2,2,2,2,1,1},
                     //{0,0,0,0,0,0,0,0}, // case 0
                     
/*
                     {1,0,0,0,0,0,0,0}, // case 1
                     {0,1,0,0,0,0,0,0},
                     {0,0,1,0,0,0,0,0},
                     {0,0,0,1,0,0,0,0},
                     {0,0,0,0,1,0,0,0},
                     {0,0,0,0,0,1,0,0},
                     {0,0,0,0,0,0,1,0},
                     {0,0,0,0,0,0,0,1},
*/
/*
                     {1,1,0,0,0,0,0,0}, // case 2
                     {0,1,1,0,0,0,0,0},
                     {0,0,1,1,0,0,0,0},
                     {0,0,0,1,1,0,0,0},
                     {0,0,0,0,1,1,0,0},
                     {0,0,0,0,0,1,1,0},
                     {0,0,0,0,0,0,1,1},
                     {1,0,0,0,0,0,0,1},
*/
/*
                     {1,0,0,0,0,1,0,0}, // case 3
                     {0,1,0,0,0,0,1,0},
                     {0,0,1,0,0,0,0,1},
                     {1,0,0,1,0,0,0,0},
                     {0,1,0,0,1,0,0,0},
                     {0,0,1,0,0,1,0,0},
                     {0,0,0,1,0,0,1,0},
                     {0,0,0,0,1,0,0,1},
*/
/*
                     {1,0,0,0,0,0,1,0}, // case 4
                     {0,1,1,1,0,0,0,0}, // case 5
                     {1,1,0,0,0,0,1,0}, // case 6
                     {0,1,0,0,1,0,1,0}, // case 7
                     {1,1,1,1,0,0,0,0}, // case 8
                     {1,0,1,1,0,0,0,1}, // case 9
                     {1,0,1,0,1,0,1,0}, // case 10
                     {1,0,1,1,0,0,1,0}, // case 11
                     {0,1,1,1,1,0,0,0}, // case 12
                     {1,0,1,0,0,1,0,1}, // case 13
                     {0,1,1,1,0,0,0,1}, // case 14                     
*/
                     };
    
    double coords[][3] = {{0,0,0}, {1,0,0}, {1,1,0}, {0,1,0},
                          {0,0,1}, {1,0,1}, {1,1,1}, {0,1,1}};
    double crx[19][3] =  {{0.5,0,0}, {1,0.5,0}, {0.5,1,0}, {0,0.5,0},   
                          {0.5,0,1}, {1,0.5,1}, {0.5,1,1}, {0,0.5,1},   
                          {0,0,0.5}, {1,0,0.5}, {1,1,0.5}, {0,1,0.5},  
                          {0.5,0,0.5}, {1,0.5,0.5}, {0.5,1,0.5}, {0,0.5,0.5}, {0.5,0.5,0}, {0.5,0.5,1},
                          {0.5,0.5,0.5}};   
    
    std::map<int, double> volumes;
    std::map<int, std::vector<EBM3D_TRIANGLE> > intfc;
    
       
    for(int i=0; i<sizeof(comp)/sizeof(comp[0]); i++)
    {
        std::ostringstream ossfilename;
        const char *filename;
        
        std::cout << "checking case " << i << " ..." << std::endl;
        
        volumes.clear();
        intfc.clear();
        init();
        
        setCompCoordsCrx(comp[i],coords,crx);        
        setCellUnknownIndex(0);
        getVolumeAndIntfc(volumes, intfc);
        
        std::cout << "intfc.size() = " << intfc.size() << std::endl;
        std::cout << "volumes.size() = " << volumes.size() << std::endl;
        ossfilename << "ebm_cell_" << i << ".plt";
        filename = ossfilename.str().c_str();
        saveIntfc_tecplot(filename, intfc);
        test_checkVolumes(volumes);
        
        std::cout << std::endl;
    }    
}

void EBM3D_CELL::test_checkVolumes(std::map<int, double> &volumes)
{
    std::map<int, double>::iterator iter;
    iter = volumes.begin();
    double total = 0;
    while(iter!=volumes.end())
    {
        total += (*iter).second;
        iter++;
    }
    std::cout << "total volumes is " << total << std::endl;
}

// calculate total area for each face.
void EBM3D_CELL::test_getPartialFaces(std::map<int, std::vector<EBM3D_PARTIALFACE> > &partialfaces)
{
    double area = 0;
    std::map<int, std::vector<EBM3D_PARTIALFACE> >::iterator mapIter;
    
    mapIter = partialfaces.begin();
    printf("face areas:");
    while(mapIter!=partialfaces.end())
    {
        area = 0;
        for(int i=0; i<(*mapIter).second.size(); i++)
            area += (*mapIter).second[i].m_area;
        printf("%f  ", area);
        mapIter++;
    }
    printf("\n");
}

void EBM3D_CELL::saveIntfc_tecplot(const char *filename, std::map<int, std::vector<EBM3D_TRIANGLE> > &intfc)
{
    FILE *hfile;
    if ((hfile = fopen(filename,"w")) == NULL)
    {
        (void) printf("EBM3D_CELL::saveIntfc_tecplot: "
                      "can't open %s\n",filename);
        return;
    }

    fprintf(hfile, "VARIABLES = \"x\", \"y\", \"z\"\n");
    
    std::map<int, std::vector<EBM3D_TRIANGLE> >::iterator iterMap;
    
    EBM3D_TRIANGLE triangle;
    
    int first;
    iterMap = intfc.begin();        
    while(iterMap!=intfc.end())
    {
/*
        fprintf(hfile, "ZONE N=%d, E=%d, DATAPACKING=POINT, ZONETYPE=FETRIANGLE\n",
                (*iterMap).second.size()*3, (*iterMap).second.size());
        
        for(int i=0; i<(*iterMap).second.size(); i++)
        {
            triangle = (*iterMap).second[i];
            for(int j=0; j<3; j++)
                fprintf(hfile,"%-9g %-9g %-9g\n", triangle.m_vertex[j][0], triangle.m_vertex[j][1], triangle.m_vertex[j][2]);            
        }
        
        for(int i=0; i<(*iterMap).second.size(); i++)
            fprintf(hfile,"%d %d %d\n", 1+i*3+0, 1+i*3+1, 1+i*3+2);        
*/
        first = (*iterMap).first;
        fprintf(hfile, "# intfc patch %d \n", first);
        for(int i=0; i<(*iterMap).second.size(); i++)
        {
            fprintf(hfile, "ZONE N=%d, E=%d, DATAPACKING=POINT, ZONETYPE=FETRIANGLE\n",
                                3,1);
            triangle = (*iterMap).second[i];
            for(int j=0; j<3; j++)
                fprintf(hfile,"%-9g %-9g %-9g\n", triangle.m_vertex[j][0], triangle.m_vertex[j][1], triangle.m_vertex[j][2]);            
            fprintf(hfile, "1 2 3\n");
        }        
        
        iterMap++;
    }   
    (void) fclose(hfile);
}

void EBM3D_CELL::setFaceVertices(int corner[4], int v0, int v1, int v2, int v3)
{
    corner[0] = v0;
    corner[1] = v1;
    corner[2] = v2;
    corner[3] = v3;
}
void EBM3D_CELL::setFaceEdges(int edge[5], int e0, int e1, int e2, int e3, int e4)
{
    edge[0] = e0;
    edge[1] = e1;
    edge[2] = e2;
    edge[3] = e3;
    edge[4] = e4;
}
    

void EBM3D_CELL::setTetraVertices(int vertices[4], int v0, int v1, int v2, int v3)
{
    vertices[0] = v0;
    vertices[1] = v1;
    vertices[2] = v2;
    vertices[3] = v3;    
}
void EBM3D_CELL::setTetraEdges(int edges[6], int e0, int e1, int e2, int e3, int e4, int e5)
{
    edges[0] = e0;
    edges[1] = e1;
    edges[2] = e2;
    edges[3] = e3;
    edges[4] = e4;
    edges[5] = e5;   
}


void EBM3D_CELL::getCenter(double a[3], double b[3], double center[3])
{
    center[0] = 1.0/2*(a[0]+b[0]);
    center[1] = 1.0/2*(a[1]+b[1]);
    center[2] = 1.0/2*(a[2]+b[2]);
}
void EBM3D_CELL::getCenter(double a[3], double b[3], double c[3], double center[3])
{
    center[0] = 1.0/3*(a[0]+b[0]+c[0]);
    center[1] = 1.0/3*(a[1]+b[1]+c[1]);
    center[2] = 1.0/3*(a[2]+b[2]+c[2]);
}
void EBM3D_CELL::getCenter(double a[3], double b[3], double c[3], double d[3], double center[3])
{
    center[0] = 1.0/4*(a[0]+b[0]+c[0]+d[0]);
    center[1] = 1.0/4*(a[1]+b[1]+c[1]+d[1]);
    center[2] = 1.0/4*(a[2]+b[2]+c[2]+d[2]);
}

void EBM3D_CELL::getComp(int comp[4], int v0, int v1, int v2, int v3)
{
     comp[0] = m_comp[v0];
     comp[1] = m_comp[v1];
     comp[2] = m_comp[v2];
     comp[3] = m_comp[v3];
}
void EBM3D_CELL::getCoords(double coords[4][3], int v0, int v1, int v2, int v3)
{
     for(int i=0; i<3; i++)
     {
	  coords[0][i] = m_coords[v0][i];
	  coords[1][i] = m_coords[v1][i];
	  coords[2][i] = m_coords[v2][i];
	  coords[3][i] = m_coords[v3][i];
     }
}
double EBM3D_CELL::getDistance(double v0[3], double v1[3])
{
    double a[3];
    difference(v0,v1,a,3);
    return Mag3d(a);
}
double EBM3D_CELL::getDistance2(double v0[3], double v1[3])
{
     double a[3];
     difference(v0,v1,a,3);
     return a[0]*a[0]+a[1]*a[1]+a[2]*a[2];
}
void EBM3D_CELL::interpolate(double v0[3], double v1[3], double t, double p[3])
{
     p[0] = (1-t)*v0[0] + t*v1[0];
     p[1] = (1-t)*v0[1] + t*v1[1];
     p[2] = (1-t)*v0[2] + t*v1[2];
}

double EBM3D_CELL::getTriArea(double v0[3], double v1[3], double v2[3])
{
    double a[3], b[3], normal[3];
    difference(v1,v0,a,3);
    difference(v2,v0,b,3);
    
    Cross3d(a,b,normal);
    double r = Mag3d(normal);
    return r/2;    
}

// corner must be one of the 8 corners of the cube.
// crx0, crx1 must be the crossing of the cube.
double EBM3D_CELL::getTriArea(int corner, int crx0, int crx1)
{
    double area = getTriArea(m_coords[corner], m_crx[crx0], m_crx[crx1]);
    return area;
}

double EBM3D_CELL::getTetraVolume(double v0[3], double v1[3], double v2[3], double v3[3])
{
    double a[3], b[3], c[3], volume;
    difference(v1,v0,a,3);
    difference(v2,v0,b,3);
    difference(v3,v0,c,3);
    volume = 1.0/6 * Det3d(a,b,c);    
    //if(volume<0)
    //    printf("EBM3D_CELL::getTetraVolume: something is wrong! volume = %f\n", volume);
    return volume;
}

double EBM3D_CELL::getTriangularPrismVolume(double v0[3], double v1[3], double v2[3], 
                                    double v3[3], double v4[3], double v5[3])
{
    double v = getTetraVolume(v0,v1,v2,v4);
    v += getTetraVolume(v0,v4,v2,v5);
    v += getTetraVolume(v0,v4,v5,v3);
    return v;
}

int *EBM3D_CELL::m_comp;
double (*EBM3D_CELL::m_coords)[3];
double (*EBM3D_CELL::m_crx)[3];

/*
 * m_cubeEdges[i]: the vertices with which there exists edges connecting the ith 
 *                 vertices directly.
 */
int EBM3D_CELL::m_cubeEdges[8][7] = {{1,2,3,4,5,6,7},
                                    {0,2,5,6,-1,-1,-1},
                                    {0,1,3,6,-1,-1,-1},
                                    {0,2,6,7,-1,-1,-1},
                                    {0,5,6,7,-1,-1,-1},
                                    {0,1,4,6,-1,-1,-1},
                                    {0,1,2,3,4,5,7},
                                    {0,3,4,6,-1,-1,-1}   
                                    };
/*
 * m_cubeEdgePairs[i]: the two vertices for the ith edges/diagonals.
 */
int EBM3D_CELL::m_cubeEdgePairs[19][2] = {{0,1}, {1,2}, {2,3}, {3,0},
					  {4,5}, {5,6}, {6,7}, {7,4},
					  {0,4}, {1,5}, {2,6}, {3,7},
					  {0,5}, {1,6}, {3,6}, {0,7}, {0,2}, {4,6}, {0,6}
                                          };


/*****************************************************************************
 *          EBM3D_CARTESIAN
 *****************************************************************************/
/*
 * Here a different method for naming the faces is used. 
 *      0,1 for the face perpendicular to the x coordinate,
 *      2,3 for the face perpendicular to the y coordinate,
 *      4,5 for the face perpendicular to the z coordinate. 
 */

int EBM3D_CARTESIAN::m_i = -1;
int EBM3D_CARTESIAN::m_j = -1;
int EBM3D_CARTESIAN::m_k = -1;


void EBM3D_CARTESIAN::getCellCenter(double center[3], int i, int j, int k)
{
     i -= m_lbuf[0];
     j -= m_lbuf[1];
     k -= m_lbuf[2];

    RECT_GRID *rect_grid = m_pFront->rect_grid;
    center[0] = cell_center(i,0,rect_grid);
    center[1] = cell_center(j,1,rect_grid);
    center[2] = cell_center(k,2,rect_grid);
}

/*
 * face: 0, west
 *       1, east
 *       2, south
 *       3, north
 *       4, down
 *       5, up.
 */
void EBM3D_CARTESIAN::getCellFaceCenter(int face, double center[3], int i, int j, int k)
{

     i -= m_lbuf[0];
     j -= m_lbuf[1];
     k -= m_lbuf[2];

    RECT_GRID *rect_grid = m_pFront->rect_grid;

    switch((face-face%2))
    {
        case PLANE_X:   // 0,1              
            if(face==0) 
                center[0] = cell_edge(i,0,rect_grid);
            else
                center[0] = cell_edge(i+1,0,rect_grid);            
            center[1] = cell_center(j,1,rect_grid);
            center[2] = cell_center(k,2,rect_grid);
            break;
        case PLANE_Y:   // 2,3            
            if(face==2) 
                center[1] = cell_edge(j,1,rect_grid);
            else
                center[1] = cell_edge(j+1,1,rect_grid);            
            center[0] = cell_center(i,0,rect_grid);
            center[2] = cell_center(k,2,rect_grid);
            break;
        case PLANE_Z:   // 4,5
            if(face==4) 
                center[2] = cell_edge(k,2,rect_grid);
            else
                center[2] = cell_edge(k+1,2,rect_grid);            
            center[0] = cell_center(i,0,rect_grid);
            center[1] = cell_center(j,1,rect_grid);
            break;        
    }
}


/*
 * dir: 0,1,2 for the coordinates.
 */
double EBM3D_CARTESIAN::getCellEdgeCenter(int dir, int i)
{
     i -= m_lbuf[dir];

    RECT_GRID *rect_grid = m_pFront->rect_grid;
    return cell_center(i,dir,rect_grid);
}
double EBM3D_CARTESIAN::getCellEdgeLength(int dir, int i, int j, int k)
{
     i -= m_lbuf[0];
     j -= m_lbuf[1];
     k -= m_lbuf[2];

    RECT_GRID *rect_grid = m_pFront->rect_grid;
    switch(dir%3)
    {
        case 0:
            return cell_width(i, 0, rect_grid);
            break;
        case 1:
            return cell_width(j, 1, rect_grid);
            break;
        case 2:            
            return cell_width(k, 2, rect_grid);
            break;
    }
}
// corner: 0,1,2,3,4,5,6,7
int EBM3D_CARTESIAN::getCellCornerIndex(int corner, int i, int j, int k)  
{
     int length[] = {m_lbuf[0]+m_gmax[0]+m_ubuf[0],
		     m_lbuf[1]+m_gmax[1]+m_ubuf[1]};
     //int base = i + j*(m_gmax[0]+1) + k*(m_gmax[0]+1)*(m_gmax[1]+1);    
     int base = i + j*(length[0]+1) + k*(length[0]+1)*(length[1]+1);
     int additional = 0;
    
     if(corner>=4)
     {
	  additional = (length[0]+1)*(length[1]+1);
	  corner -= 4;
     }
     switch(corner)
     {
     case 0:
	  additional += 0;
	  break;
     case 1:
	  additional += 1;
	  break;
     case 2:
	  additional += length[0]+2;
	  break;
     case 3:
	  additional += length[0]+1;
	  break;            
     }
     return base + additional;
}

void EBM3D_CARTESIAN::getCellCornerCoords(int corner, double coords[3], 
					  int i, int j, int k)
 {
      i -= m_lbuf[0];
      j -= m_lbuf[1];
      k -= m_lbuf[2];

     RECT_GRID *rect_grid = m_pFront->rect_grid;
     double x,y,z,dx,dy,dz;
     dx = cell_width(i,0,rect_grid);    
     dy = cell_width(j,1,rect_grid);    
     dz = cell_width(k,2,rect_grid);
     x = cell_edge(i,0,rect_grid);
     y = cell_edge(j,1,rect_grid);
     z = cell_edge(k,2,rect_grid);

     if(corner<4)
	  coords[2] = z;
     else
     {
	  coords[2] = z + dz;
	  corner -= 4;
     }

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
void EBM3D_CARTESIAN::getCellCornerCoords(double coords[8][3], int i, int j, int k)
 {
      i -= m_lbuf[0];
      j -= m_lbuf[1];
      k -= m_lbuf[2];

     RECT_GRID *rect_grid = m_pFront->rect_grid;
     double x,y,z,dx,dy,dz;
     dx = cell_width(i,0,rect_grid);    
     dy = cell_width(j,1,rect_grid);    
     dz = cell_width(k,2,rect_grid);

     // 0
     coords[0][0] = x = cell_edge(i,0,rect_grid);
     coords[0][1] = y = cell_edge(j,1,rect_grid);
     coords[0][2] = z = cell_edge(k,2,rect_grid);
     // 1 
     coords[1][0] = x + dx;
     coords[1][1] = y;
     coords[1][2] = z;
     // 2
     coords[2][0] = x + dx;
     coords[2][1] = y + dy;
     coords[2][2] = z;
     // 3
     coords[3][0] = x;
     coords[3][1] = y + dy;
     coords[3][2] = z;
     // 4
     coords[4][0] = x;
     coords[4][1] = y;
     coords[4][2] = z + dz;
     // 5 
     coords[5][0] = x + dx;
     coords[5][1] = y;
     coords[5][2] = z + dz;
     // 6
     coords[6][0] = x + dx;
     coords[6][1] = y + dy;
     coords[6][2] = z + dz;
     // 7
     coords[7][0] = x;
     coords[7][1] = y + dy;
     coords[7][2] = z + dz;
 }


void EBM3D_CARTESIAN::getCurrentCell(int &i, int &j, int &k)
{
    i = m_i;
    j = m_j;
    k = m_k;
}
// make sure that this is called before getCellUnknownIndex
void EBM3D_CARTESIAN::setCurrentCell(int i, int j, int k)
{
    m_i = i;
    m_j = j;
    m_k = k;
}

void EBM3D_CARTESIAN::locateCell(double coords[3], int &i, int &j, int &k)
{
    RECT_GRID *rect_grid = m_pFront->rect_grid;
    i = cell_index(coords[0],0,rect_grid);
    j = cell_index(coords[1],1,rect_grid);
    k = cell_index(coords[2],2,rect_grid);

    i += m_lbuf[0];
    j += m_lbuf[1];
    k += m_lbuf[2];
}


int EBM3D_CARTESIAN::getCellCornerComp(int index)
{
    return m_comps[index];
}

int EBM3D_CARTESIAN::getCellCornerComp(int corner, int i, int j, int k)
{
     int index = getCellCornerIndex(corner, i,j,k);
     if(index<0||index>=m_nComps)
	  return -1;
     else
	  return getCellCornerComp(index);
}

int EBM3D_CARTESIAN::getCellUnknownIndex(int comp, int i, int j, int k)
{   
    int index = getCellUnknownIndexByTrial(comp,i,j,k);
    return index;
/*
    if(index>=0)
	 return index;
    printf("EBM3D_CARTESIAN::getCellUnknownIndex: can't found comp=%d "
            " with (%d,%d,%d)\n", comp, i,j,k);
    printf(" (m_i,m_j,m_k) is (%d,%d,%d)\n", m_i,m_j,m_k);
    debug_printNeighborComp();    
    exit(0);*/
}

/*
 * try to get the unknown index.
 * return 
 *         >=0 if succeed;
 *         <0  if failed.
 * This function is the same as EBM3D_CARTESIAN::getCellUnknownIndex except for the return if it failed.
 */
int EBM3D_CARTESIAN::getCellUnknownIndexByTrial(int comp, int i, int j, int k)
{
     //return getCell(i,j,k)->getCellUnknownIndex(0);
     if(i<0||j<0||k<0 || 
	i>=m_lbuf[0]+m_gmax[0]+m_ubuf[0] ||
	j>=m_lbuf[1]+m_gmax[1]+m_ubuf[1] ||
	k>=m_lbuf[2]+m_gmax[2]+m_ubuf[2])
	  return -1;

     int corner, index, debug[8];
     for(corner=0; corner<8; corner++)
     {
	  index = getCellCornerIndex(corner, i,j,k);
	  debug[corner] = getCellCornerComp(index);
	  if(debug[corner]==comp)
	       return getCell(i,j,k)->getCellUnknownIndex(corner);
     }
     return -1;
}

int EBM3D_CARTESIAN::getCellNeighborUnknownIndex(int face, int comp, 
						 int i, int j, int k)
{
    int I, J, K;
    I = i;
    J = j;
    K = k;
    
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
        case PLANE_Z:
            if(face==FACE_DOWN)     // down
                K = k-1;
            else                    // up                
                K = k+1;            
            break;           
    }
    return getCellUnknownIndex(comp, I,J,K);
}

/*
 * find the neighbor along one of the three coordinate planes: x, y, z.
 * parameter:
 *      plane: PLANE_X, PLANE_Y, PLANE_Z.
 *      x,y: -1,0,1, the difference along the two coordinates on the plane.
 *      comp: the current component.
 */ 
int EBM3D_CARTESIAN::getNeighborByPlane(EBM_PLANE plane, int x, int y, int comp, int i, int j, int k)
{   
    switch(plane)
    {
        case PLANE_X:   // east, west
            j += x;
            k += y;
            break;
        case PLANE_Y:   // south, north
            i += x;
            k += y;
            break;
        case PLANE_Z:   // down, up
            i += x;
            j += y;
            break;
    }
    return getCellUnknownIndex(comp, i,j,k);
}
int EBM3D_CARTESIAN::getNeighborByPlane(EBM_PLANE plane, int x, int y, int comp, int index[3], int i, int j, int k)
{   
    
    switch(plane)
    {
        case PLANE_X:   // east, west
	    index[0] = i;
	    index[1] = j + x;
	    index[2] = k + y;
            break;
        case PLANE_Y:   // south, north
	    index[0] = i + x;
	    index[1] = j;
	    index[2] = k + y;
            break;
        case PLANE_Z:   // down, up
	    index[0] = i + x;
	    index[1] = j + y;
	    index[2] = k;
            break;
    }
    return getCellUnknownIndex(comp, index[0],index[1],index[2]);
}

/*
 * used to get the knowns at cell (i,j)
 */
double EBM3D_CARTESIAN::getCellUnknown(int comp, int i, int j, int k)
{
     if(comp==m_compList[0][i][j][k])
	  return m_unknowns[0][i][j][k];
     else if(comp==m_compList[1][i][j][k])
	  return m_unknowns[1][i][j][k];
     
     printf("ID %d: EBM3D_CARTESIAN::getCellUnknown: comp=%d, i=%d,j=%d,k=%d.\n",
	    m_myid, comp, i,j,k);
     printf("ID %d:\tm_compList[*][%d][%d][%d]={%d,%d}\n",
	    m_myid,i,j,k,m_compList[0][i][j][k],m_compList[1][i][j][k]);
	    
     return HUGE_VAL;
}

void EBM3D_CARTESIAN::getCellCompList(int comp[2], int i, int j, int k)
{
     comp[0] = m_compList[0][i][j][k];
     comp[1] = m_compList[1][i][j][k];
}
void EBM3D_CARTESIAN::setCellCompList(int comp[2], int i, int j, int k)
{
     m_compList[0][i][j][k] = comp[0];
     m_compList[1][i][j][k] = comp[1];
}

void EBM3D_CARTESIAN::getCellUnknown(double unknown[2], int i, int j, int k)
{
     unknown[0] = m_unknowns[0][i][j][k];
     unknown[1] = m_unknowns[1][i][j][k];
}
void EBM3D_CARTESIAN::setCellUnknown(double unknown[2], int i, int j, int k)
{
     m_unknowns[0][i][j][k] = unknown[0];
     m_unknowns[1][i][j][k] = unknown[1];
}

/*
 * should this function be moved into EBM2D_CELL?
 */
void EBM3D_CARTESIAN::setCellUnknown(double *x)
{
     int index, c0, c1, numCellUnknowns;
     CELL_TYPE type;

     index = 0;
     
     for(int k=m_lbuf[2]; k<m_lbuf[2]+m_gmax[2]; k++)
     for(int j=m_lbuf[1]; j<m_lbuf[1]+m_gmax[1]; j++)
     for(int i=m_lbuf[0]; i<m_lbuf[0]+m_gmax[0]; i++)
     {   
	  type = getCell(i,j,k)->getCellType();
	  numCellUnknowns = getCell(i,j,k)->getCellUnknownNumber();
	  getCellComp(c0, c1, i, j, k);
	  m_compList[0][i][j][k] = c0;
	  m_compList[1][i][j][k] = c1;
	  switch(type)
	  {
	  case INTERNAL:   // numCellUnknowns == 1
	       m_unknowns[0][i][j][k] = x[index];
	       break;
	  case PARTIAL:
	       m_unknowns[0][i][j][k] = x[index];
	       m_unknowns[1][i][j][k] = x[index+1];
	       break;
	  case BOUNDARY:
	       if(c0>=0)
		    m_unknowns[0][i][j][k] = x[index];
	       else
		    m_unknowns[1][i][j][k] = x[index];
	       break;
	  }
	  if(numCellUnknowns>0)
	  {
	       numCellUnknowns = numCellUnknowns * 2 - 1;
	       index += numCellUnknowns;
	  }
     }         
}


EBM3D_CELL *EBM3D_CARTESIAN::getCell(int i, int j, int k)
{
    return &m_cells[i][j][k];
}

CELL_TYPE EBM3D_CARTESIAN::getCellType(int i, int j, int k)
{
     getCell(i,j,k)->getCellType();
}

/* 
 * get the two distinct components stored at the corners of the cell. The comps
 * can be negative.
 */ 
void EBM3D_CARTESIAN::getCellComp(int &c0, int &c1, int i, int j, int k)
{
     int index = getCellCornerIndex(0, i,j,k);
     c0 = m_comps[index];
     for(int corner=1; corner<8; corner++)
     {
	  index = getCellCornerIndex(corner,i,j,k);
	  if(index<0||index>m_nComps)
	       c1 = -1;
	  else
	       c1 = m_comps[index];
	  if(c1!=c0)
	       return;
     }
}

void EBM3D_CARTESIAN::getCellComp(int comp[8], int i, int j, int k)
{
    int index[8];     
    for(int corner=0; corner<8; corner++)
    {
        index[corner] = getCellCornerIndex(corner,i,j,k);
	if(index[corner]<0||index[corner]>=m_nComps)
	     comp[corner] = -1;
	else
	     comp[corner] = m_comps[index[corner]];
    }
}

void EBM3D_CARTESIAN::getCellCompCoordsCrx(int comp[8], double coords[8][3], double crx[19][3], int i, int j, int k)
{
    int index[8];  
    
    // comp
    getCellComp(comp,i,j,k);
    
    // coords
    getCellCornerCoords(coords,i,j,k);

    // crx
    int pa, pb;
    for(int m=0; m<19; m++)
    {        
        pa = EBM3D_CELL::m_cubeEdgePairs[m][0];
        pb = EBM3D_CELL::m_cubeEdgePairs[m][1];
        
        if(comp[pa]!=comp[pb])        
            getCrx(coords[pa], coords[pb], crx[m]);
	else
	{
	     crx[m][0] = crx[m][1] = crx[m][2] = HUGE_VAL;
	}
    }    
}


int EBM3D_CARTESIAN::getComp(double *coords)
{
     /*double min = -0.7, max = 0.7;
     if( coords[0]<min || coords[0]>max || 
	 coords[1]<min || coords[1]>max ||
	 coords[2]<min || coords[2]>max)
	  return -1;
     else
     return 1;*/
     
     //return 1;

    const double R = 0.5;
    const double Center[3] = {0, 0, 0};
    double distance = sqrt(sqr(coords[0] - Center[0]) + sqr(coords[1] - Center[1]) +
			sqr(coords[2] - Center[2])) - R;
    
    if(distance<m_mindd*m_mindd)
        return 1;
    else
        return 2;
}

void EBM3D_CARTESIAN::getCrx(double p0[3], double p1[3], double crx[3])
{
     //for(int i=0; i<3; i++)
     //	  crx[i] = 1.0/2 * (p0[i] + p1[i]);

    const int MAX_ITER = 40;
    double a[3], b[3];

    a[0] = p0[0];   a[1] = p0[1];   a[2] = p0[2];   
    b[0] = p1[0];   b[1] = p1[1];   b[2] = p1[2];   

    int ca, cb, cc;
    ca = getComp(a);
    cb = getComp(b);
    if(ca==cb)
    {
        crx[0] = HUGE_VAL;
        crx[1] = HUGE_VAL;
        crx[2] = HUGE_VAL;
        return;
    }
    
    for(int i=0; i<MAX_ITER; i++)
    {
        crx[0] = 0.5*(a[0]+b[0]);
        crx[1] = 0.5*(a[1]+b[1]);
        crx[2] = 0.5*(a[2]+b[2]);
        cc = getComp(crx);
        if(cc==ca)
        {
            a[0] = crx[0];
            a[1] = crx[1];
            a[2] = crx[2];
        }
        else
        {
            b[0] = crx[0];
            b[1] = crx[1];
            b[2] = crx[2];
        }
    }
    
    
    double length = EBM3D_CELL::getDistance(p0,p1);
    double length0 = EBM3D_CELL::getDistance(p0,crx);
    
    double t = length0/length;
    //double T = 1.0/100;
    double T = m_mindd*m_mindd*m_mindd*m_mindd;
    if(t<T)
    {
	 EBM3D_CELL::interpolate(p0,p1,T,crx);
	 //printf("EBM3D_CARTESIAN::getCrx: t<T.\n");
    }
    if((1-t)<T)
    {
	 EBM3D_CELL::interpolate(p0,p1,1-T,crx);
	 //printf("EBM3D_CARTESIAN::getCrx: 1-t<T.\n");
    }
}

// modify components.
void EBM3D_CARTESIAN::setEquivalentComponents(void)
{
    for(int k=0; k<=m_gmax[2]; k++)
    for(int j=0; j<=m_gmax[1]; j++)
    for(int i=0; i<=m_gmax[0]; i++)   
    {
        if(m_comps[getCellCornerIndex(0,i,j,k)]==0)
            m_comps[getCellCornerIndex(0,i,j,k)] = 2;
    }    
}

/* 
 * the following also deals with the Dirichlet boundary condition.
 * face ordering:
 *      0,1 for faces perpendicular to x coordinates
 *      2,3 for faces perpendicular to y coordinates
 *      4,5 for faces perpendicular to z coordinates
 */
void EBM3D_CARTESIAN::setMatrixForInternalCell(int i, int j, int k)
{
    int comp = getCellCornerComp(getCellCornerIndex(0,i,j,k));
    
    int index0, index1, index2, index3, index4, index5, index;
    double cellCenter[3], faceCenter[6][3], beta;
    double dx, dy, dz, rhs, value;
    
    dx = getCellEdgeLength(0,i,j,k);      // uniform edge length only
    dy = getCellEdgeLength(1,i,j,k);      // uniform edge length only    
    dz = getCellEdgeLength(2,i,j,k);      // uniform edge length only 
        
    // indices for the cells
    //index  = getCell(i,j,k)->getCellUnknownIndex(0);    
    //index0 = (i!=0) ? getCellUnknownIndex(comp,i-1,j,k) : -1;      
    //index1 = (i!=(m_gmax[0]-1)) ? getCellUnknownIndex(comp,i+1,j,k) : -1;
    //index2 = (j!=0) ? getCellUnknownIndex(comp,i,j-1,k) : -1;
    //index3 = (j!=(m_gmax[1]-1)) ? getCellUnknownIndex(comp,i,j+1,k) : -1;
    //index4 = (k!=0) ? getCellUnknownIndex(comp,i,j,k-1) : -1;
    //index5 = (k!=(m_gmax[2]-1)) ? getCellUnknownIndex(comp,i,j,k+1) : -1;
    index  = getCell(i,j,k)->getCellUnknownIndex(0);    
    index0 = getCellUnknownIndex(comp,i-1,j,k);      
    index1 = getCellUnknownIndex(comp,i+1,j,k);
    index2 = getCellUnknownIndex(comp,i,j-1,k);
    index3 = getCellUnknownIndex(comp,i,j+1,k);
    index4 = getCellUnknownIndex(comp,i,j,k-1);
    index5 = getCellUnknownIndex(comp,i,j,k+1);
    
    // cell center
    getCellCenter(cellCenter,i,j,k);
    // face centers
    for(int face=0; face<6; face++)
        getCellFaceCenter(face,faceCenter[face],i,j,k);  
    
    // now for each faces, setup the system of equations    

    // face 0 
    beta = m_pLaplace->getBeta(faceCenter[0], comp);
    if(index0>=0)
    {            
        m_pSolver->Add_A(index, index0, beta*dy*dz/dx);
        m_pSolver->Add_A(index, index, -beta*dy*dz/dx);
    }
    else
    {
        value = m_pLaplace->getExactSolution(faceCenter[0], comp);
        m_pSolver->Add_A(index, index, -2*beta*dy*dz/dx);
        m_pSolver->Add_b(index, -2*value*beta*dy*dz/dx); 
    }    
    // face 1
    beta = m_pLaplace->getBeta(faceCenter[1], comp);
    if(index1>=0)
    {
        m_pSolver->Add_A(index, index1, beta*dy*dz/dx);
        m_pSolver->Add_A(index, index, -beta*dy*dz/dx);
    }
    else
    {
        value = m_pLaplace->getExactSolution(faceCenter[1], comp);
        m_pSolver->Add_A(index, index, -2*beta*dy*dz/dx);
        m_pSolver->Add_b(index, -2*value*beta*dy*dz/dx); 
    }
    
    // face 2
    beta = m_pLaplace->getBeta(faceCenter[2], comp);
    if(index2>=0)
    {
        m_pSolver->Add_A(index, index2, beta*dx*dz/dy);
        m_pSolver->Add_A(index, index, -beta*dx*dz/dy);
    }
    else
    {
        value = m_pLaplace->getExactSolution(faceCenter[2], comp);
        m_pSolver->Add_A(index, index, -2*beta*dx*dz/dy);
        m_pSolver->Add_b(index, -2*value*beta*dx*dz/dy); 
    }
    // face 3
    beta = m_pLaplace->getBeta(faceCenter[3], comp);
    if(index3>=0)
    {
        m_pSolver->Add_A(index, index3, beta*dx*dz/dy);
        m_pSolver->Add_A(index, index, -beta*dx*dz/dy);
    }
    else
    {
        value = m_pLaplace->getExactSolution(faceCenter[3], comp);
        m_pSolver->Add_A(index, index, -2*beta*dx*dz/dy);
        m_pSolver->Add_b(index, -2*value*beta*dx*dz/dy); 
    }
    
    // face 4
    beta = m_pLaplace->getBeta(faceCenter[4], comp);
    if(index4>=0)
    {
        m_pSolver->Add_A(index, index4, beta*dx*dy/dz);
        m_pSolver->Add_A(index, index, -beta*dx*dy/dz);
    }
    else
    {
        value = m_pLaplace->getExactSolution(faceCenter[4], comp);
        m_pSolver->Add_A(index, index, -2*beta*dx*dy/dz);
        m_pSolver->Add_b(index, -2*value*beta*dx*dy/dz); 
    }
    // face 5
    beta = m_pLaplace->getBeta(faceCenter[5], comp);
    if(index5>=0)
    {
        m_pSolver->Add_A(index, index5, beta*dx*dy/dz);
        m_pSolver->Add_A(index, index, -beta*dx*dy/dz);
    }
    else
    {
        value = m_pLaplace->getExactSolution(faceCenter[5], comp);
        m_pSolver->Add_A(index, index, -2*beta*dx*dy/dz);
        m_pSolver->Add_b(index, -2*value*beta*dx*dy/dz); 
    }
    
    // rhs
    rhs = m_pLaplace->getRightHandSide(cellCenter, comp);
    m_pSolver->Add_b(index, rhs*dx*dy*dz);           
       
}

void EBM3D_CARTESIAN::setMatrixForPartialCell_2Unknowns(int i, int j, int k)
{
     //printf("EBM3D_CARTESIAN::setMatrixForPartialCell_2Unknowns: (%d,%d,%d)\n",i,j,k);
    int comp[8], size, face;
    double coords[8][3], crx[19][3], beta;   
        
    getCellCompCoordsCrx(comp,coords,crx,i,j,k);    
    getCell(i,j,k)->setCompCoordsCrx(comp,coords,crx);

    // intfc
    std::map<int, double> volumes;
    std::map<int, std::vector<EBM3D_TRIANGLE> > mapIntfc;  
    std::map<int, double> mapIntfcArea;
    std::map<int, EBM3D_POINT> mapIntfcNormal;
    std::map<int, EBM3D_POINT> mapIntfcCenter;
    
    getCell(i,j,k)->getVolumeAndIntfc(volumes, mapIntfc);
    getCell(i,j,k)->getAveragedIntfcAreaNormalCenter(mapIntfc, mapIntfcArea, mapIntfcNormal, mapIntfcCenter);
    
    if(mapIntfcArea.size()>1)
    {
        printf("EBM3D_CARTESIAN::setMatrixForPartialCell_2Unknowns: intfcArea.size()>1 is undelt case!");
        exit(0);
    }
    
    std::map<int, std::vector<EBM3D_TRIANGLE> >::iterator iterMapIntfc;
    EBM3D_TRIANGLE triangle;
    EBM3D_POINT point;
    
    int index = getCell(i,j,k)->getCellIntfcIndex();    // the first intfc index.
    int intfcIndex[3], intfcComp[2], p;
    double intfcCenter[3], intfcArea, intfcNormal[3], volumeFraction[2];
    
    iterMapIntfc = mapIntfc.begin();
    while(iterMapIntfc!=mapIntfc.end())
    {        
        triangle = (*iterMapIntfc).second[0];
        triangle.getIndex(intfcIndex);        
        intfcIndex[2] = index++;
        triangle.getComp(intfcComp);
        
        p = (*iterMapIntfc).first;
        intfcArea = mapIntfcArea[p];
        mapIntfcNormal[p].getCoords(intfcNormal);
        mapIntfcCenter[p].getCoords(intfcCenter);
        
	volumeFraction[0] = getVolumeFraction(intfcIndex[0], volumes);
	volumeFraction[1] = getVolumeFraction(intfcIndex[1], volumes);
	
        setMatrixForInterface(intfcIndex, intfcCenter, intfcNormal, intfcComp, intfcArea, volumeFraction, i,j,k);
	//setMatrixForInterface2(intfcIndex, intfcCenter, intfcNormal, intfcComp, intfcArea, volumeFraction, i,j,k);
        
        iterMapIntfc++;
    }    
    
    // the two partial volumes    
    std::map<int, double>::iterator iterVolumes;
    std::map<int, int> unknownToComp;    
    getCell(i,j,k)->getUnknownToCompMap(unknownToComp);
    
    int c;
    double vol, rhs, cellCenter[3];
    getCellCenter(cellCenter,i,j,k);
    
    iterVolumes = volumes.begin();
    while(iterVolumes!=volumes.end())
    {
        index = (*iterVolumes).first;
        c = unknownToComp[index];
        vol = (*iterVolumes).second;
        rhs = m_pLaplace->getRightHandSide(cellCenter, c);
	
	volumeFraction[0] = getVolumeFraction(index, volumes);
        m_pSolver->Add_b(index, rhs*vol*volumeFraction[0]);   
        
        iterVolumes++;
    }

        
    // faces
    std::map<int, std::vector<EBM3D_PARTIALFACE> > mapPartialFaces;
    std::map<int, std::vector<EBM3D_PARTIALFACE> >::iterator iterFace;
    std::vector<EBM3D_PARTIALFACE>::iterator iterPartialFace;
    std::vector< std::pair<int,double> > stencil;
    std::vector< std::pair<int,double> >::iterator iterStencil;
    EBM3D_PARTIALFACE partialface;
    
    getCell(i,j,k)->getPartialFaces(mapPartialFaces);
    
    iterFace = mapPartialFaces.begin();
    while(iterFace!=mapPartialFaces.end())
    {
        face = (*iterFace).first;
        iterPartialFace = (*iterFace).second.begin();
        while(iterPartialFace!=(*iterFace).second.end())
        {
            partialface = (*iterPartialFace);            
            stencil.clear(); 
	    //getFluxStencilForPartialFace(face, partialface.m_center, partialface.m_comp, stencil, i,j,k);            
	    getFluxStencilForPartialFace2(face, partialface.m_center, partialface.m_comp, stencil, i,j,k);            
	    
	    beta = m_pLaplace->getBeta(partialface.m_center, partialface.m_comp);
     
	    volumeFraction[0] = getVolumeFraction(partialface.m_index, volumes);
	    iterStencil = stencil.begin();            
            while(iterStencil!=stencil.end())
            {
                m_pSolver->Add_A(partialface.m_index, (*iterStencil).first, 
				 beta*partialface.m_area * (*iterStencil).second * volumeFraction[0]);                
                iterStencil++;
            }             
            iterPartialFace++;
        }
        iterFace++;        
    }
        
    
}

/*
 * now for simplicity, only Neumann boundary condition (given) is delt with.
 */
void EBM3D_CARTESIAN::setMatrixForBoundaryCell(int i, int j, int k)
{
     printf("EBM3D_CARTESIAN::setMatrixForBoundaryCell: (%d,%d,%d)\n", i,j,k);
     //printf("\t unknown index = %d\n", getCellUnknownIndex(1,i,j,k));
    
    int comp[8], size, face;
    double coords[8][3], crx[19][3], beta;   
        
    getCellCompCoordsCrx(comp,coords,crx,i,j,k);    
    getCell(i,j,k)->setCompCoordsCrx(comp,coords,crx);

    // boundary
    std::map<int, double> volumes;
    std::map<int, std::vector<EBM3D_TRIANGLE> > mapIntfc;  
    std::map<int, double> mapIntfcArea;
    std::map<int, EBM3D_POINT> mapIntfcNormal;
    std::map<int, EBM3D_POINT> mapIntfcCenter;
    
    getCell(i,j,k)->getVolumeAndIntfc(volumes, mapIntfc);
    //EBM3D_CELL::test_checkVolumes(volumes);
    getCell(i,j,k)->getAveragedIntfcAreaNormalCenter(mapIntfc, mapIntfcArea, mapIntfcNormal, mapIntfcCenter);

    if(mapIntfcArea.size()>1)
    {
        printf("EBM3D_CARTESIAN::setMatrixForBoundaryCell: intfcArea.size()>1 is undelt case!");
        exit(0);
    }
    
    std::map<int, std::vector<EBM3D_TRIANGLE> >::iterator iterMapIntfc;
    EBM3D_TRIANGLE triangle;
    EBM3D_POINT point;
    
    int index;
    int intfcIndex[2], intfcComp[2], p;
    double intfcCenter[3], intfcArea, intfcNormal[3], volumeFraction[2];
    
    iterMapIntfc = mapIntfc.begin();
    while(iterMapIntfc!=mapIntfc.end())
    {        
        triangle = (*iterMapIntfc).second[0];
        triangle.getIndex(intfcIndex);          
        triangle.getComp(intfcComp);
        
        p = (*iterMapIntfc).first;
        intfcArea = mapIntfcArea[p];
        mapIntfcNormal[p].getCoords(intfcNormal);
        mapIntfcCenter[p].getCoords(intfcCenter);
        
	volumeFraction[0] = getVolumeFraction(intfcIndex[0], volumes);
	volumeFraction[1] = getVolumeFraction(intfcIndex[1], volumes);
	
        setMatrixForBoundary(intfcIndex, intfcCenter, intfcNormal, intfcComp, intfcArea, volumeFraction, i,j,k);
        //printf("\t intfcCenter = {%e,%e,%e}, intfcArea = %e \n", intfcCenter[0], intfcCenter[1], intfcCenter[2], intfcArea);
        iterMapIntfc++;
    }    
    
    // the partial volumes    
    std::map<int, double>::iterator iterVolumes;
    std::map<int, int> unknownToComp;    
    getCell(i,j,k)->getUnknownToCompMap(unknownToComp);
    
    int c;
    double vol, rhs, cellCenter[3];
    getCellCenter(cellCenter,i,j,k);
    
    iterVolumes = volumes.begin();
    while(iterVolumes!=volumes.end())
    {
        index = (*iterVolumes).first;
        c = unknownToComp[index];
	if(c>=0)               // inside the domain
	{     
	     vol = (*iterVolumes).second;
	     rhs = m_pLaplace->getRightHandSide(cellCenter, c);

	     volumeFraction[0] = getVolumeFraction(index, volumes);
	     m_pSolver->Add_b(index, rhs*vol*volumeFraction[0]);
	     //printf("%d: %e\n", index, rhs*vol);
	     //printf("\tvol = %e\n", vol);
        }
        iterVolumes++;
    } 


    // faces
    std::map<int, std::vector<EBM3D_PARTIALFACE> > mapPartialFaces;
    std::map<int, std::vector<EBM3D_PARTIALFACE> >::iterator iterFace;
    std::vector<EBM3D_PARTIALFACE>::iterator iterPartialFace;
    std::vector< std::pair<int,double> > stencil;
    std::vector< std::pair<int,double> >::iterator iterStencil;
    EBM3D_PARTIALFACE partialface;
    
    getCell(i,j,k)->getPartialFaces(mapPartialFaces);
    
    iterFace = mapPartialFaces.begin();
    while(iterFace!=mapPartialFaces.end())
    {
        face = (*iterFace).first;
        iterPartialFace = (*iterFace).second.begin();
        while(iterPartialFace!=(*iterFace).second.end())
        {
            partialface = (*iterPartialFace); 
	    if(partialface.m_comp>=0)             // check whether it is outside of the domain.
	    {
		 //printf("\tface = %d, m_center = {%e,%e,%e} \n", face, 
		 //partialface.m_center[0], partialface.m_center[1], partialface.m_center[2]);
		 stencil.clear(); 
		 getFluxStencilForPartialFace(face, partialface.m_center, 
					    partialface.m_comp, stencil, i,j,k);		 
		 beta = m_pLaplace->getBeta(partialface.m_center, partialface.m_comp);
		 
		 volumeFraction[0] = getVolumeFraction(partialface.m_index, volumes);
		 iterStencil = stencil.begin();            
		 while(iterStencil!=stencil.end())
		 {
		      m_pSolver->Add_A(partialface.m_index, (*iterStencil).first, 
				       beta*partialface.m_area * (*iterStencil).second * volumeFraction[0]);                
		      iterStencil++;
		 }             
	    }
            iterPartialFace++;
        }
        iterFace++;        
    }
        
    
    // test_EBM3D -s 40
    //if(i==20&&j==6&&k==6)
    //{
    // printf("\tsavePartialCells_Tecplot(...)\n");
    // savePartialCells_Tecplot("partialcell_20_6_6.plt", mapPartialFaces, mapIntfc, i,j,k);
    //}

    /*char filename[200];
    sprintf(filename, "partialcell_%d_%d_%d.plt", i,j,k);
    savePartialCells_Tecplot(filename, mapPartialFaces, mapIntfc, i,j,k);*/
}


void EBM3D_CARTESIAN::setMatrixForInterface(int intfcIndex[3], double intfcCenter[3], 
					    double normal[3], int intfcComp[2], 
					    double intfcArea, double volumeFraction[2], 
					    int i, int j, int k)
{
     //debug_printNormal(intfcCenter, normal);

     /*double flux, value;            
    flux = m_pLaplace->getFlux(intfcCenter, intfcComp[0], normal);  
    beta = m_pLaplace->getBeta(intfcCenter, intfcComp[0]);
    m_pSolver->Add_b(intfcIndex[0], -flux*beta*intfcArea);    // negative because it is moved to the rhs.
    
    flux = -m_pLaplace->getFlux(intfcCenter, intfcComp[1], normal); // negative because it has negative normal
    beta = m_pLaplace->getBeta(intfcCenter, intfcComp[1]);
    m_pSolver->Add_b(intfcIndex[1], -flux*beta*intfcArea);    // negative because it is moved to the rhs.
    
    value = m_pLaplace->getExactSolution(intfcCenter, intfcComp[0]);  // assuming the solution is continuous at the interface.
    m_pSolver->Add_A(intfcIndex[2], intfcIndex[2], 1);
    m_pSolver->Add_b(intfcIndex[2], value); 

    return; */

    std::vector<int> indices[2];
    std::vector<double> coefficients[2];
    std::vector<double> distances[2];
    double beta[2];
    
    // for component 1
    beta[1] = m_pLaplace->getBeta(intfcCenter, intfcComp[1]);
    getSolutionInterpStencil(intfcCenter, normal, intfcComp[1], indices[1], coefficients[1], distances[1], i,j,k);    
        
    // for component 0;    
    normal[0] = -normal[0];     normal[1] = -normal[1];     normal[2] = -normal[2];
    beta[0] = m_pLaplace->getBeta(intfcCenter, intfcComp[0]);
    getSolutionInterpStencil(intfcCenter, normal, intfcComp[0], indices[0], coefficients[0], distances[0], i,j,k);
    
    double coords[2], coeffs[2][2];
    int n;
    
    // component 0 flux
    coords[0] = 0;    
    coords[1] = distances[0][0];
    get2PointsDerivativeInterpCoeffs(coords, 0, coeffs[0]);    

    m_pSolver->Add_A(intfcIndex[0], intfcIndex[2], - beta[0]*intfcArea*coeffs[0][0] * volumeFraction[0]);
    for(n=0; n<coefficients[0].size(); n++)
        m_pSolver->Add_A(intfcIndex[0], indices[0][n], -beta[0]*intfcArea*coeffs[0][1]*coefficients[0][n] * volumeFraction[0]);
    
    // component 1 flux
    coords[0] = 0;    coords[1] = distances[1][0];
    get2PointsDerivativeInterpCoeffs(coords, 0, coeffs[1]);   
    m_pSolver->Add_A(intfcIndex[1], intfcIndex[2], - beta[1]*intfcArea*coeffs[1][0] * volumeFraction[1]);
    for(n=0; n<coefficients[1].size(); n++)
        m_pSolver->Add_A(intfcIndex[1], indices[1][n], -beta[1]*intfcArea*coeffs[1][1]*coefficients[1][n] * volumeFraction[1]);
    
    // interface equation
    m_pSolver->Add_A(intfcIndex[2], intfcIndex[2], + beta[0]*intfcArea*coeffs[0][0]);
    for(n=0; n<coefficients[0].size(); n++)
        m_pSolver->Add_A(intfcIndex[2], indices[0][n], + beta[0]*intfcArea*coeffs[0][1]*coefficients[0][n]);
    
    m_pSolver->Add_A(intfcIndex[2], intfcIndex[2], + beta[1]*intfcArea*coeffs[1][0]);
    for(n=0; n<coefficients[1].size(); n++)
        m_pSolver->Add_A(intfcIndex[2], indices[1][n], + beta[1]*intfcArea*coeffs[1][1]*coefficients[1][n]);    
    
}

void EBM3D_CARTESIAN::setMatrixForInterface2(int intfcIndex[3],double intfcCenter[3], 
					     double intfcNormal[3], int intfcComp[2], 
					     double intfcArea, double volumeFraction[2], 
					     int i, int j, int k)
{

     std::vector<int> indices[2];
     std::vector<double> coefficients[2];
     std::vector<double> distances[2];
     double beta[2];

     // for component 1
     beta[1] = m_pLaplace->getBeta(intfcCenter, intfcComp[1]);
     getSolutionInterpStencil2(intfcCenter, intfcNormal, intfcComp[1], indices[1], coefficients[1], distances[1],i,j,k);    

     // for component 0;    
     intfcNormal[0] = -intfcNormal[0];     
     intfcNormal[1] = -intfcNormal[1]; 
     intfcNormal[2] = -intfcNormal[2]; 
     beta[0] = m_pLaplace->getBeta(intfcCenter, intfcComp[0]);
     getSolutionInterpStencil2(intfcCenter, intfcNormal, intfcComp[0], indices[0], coefficients[0], distances[0],i,j,k);

     double coords[3], coeffs[2][3];
     int m, n;

     // component 0 flux
     coords[0] = 0;
     coords[1] = distances[0][0];
     coords[2] = distances[0][1];
     get3PointsDerivativeInterpCoeffs(coords, 0, coeffs[0]);    

     m_pSolver->Add_A(intfcIndex[0], intfcIndex[2], - beta[0]*intfcArea*coeffs[0][0] * volumeFraction[0]);
     for(m=0; m<2; m++)    
	 for(n=0; n<6; n++)
	     m_pSolver->Add_A(intfcIndex[0], indices[0][m*6+n], - beta[0]*intfcArea*coeffs[0][m+1]*coefficients[0][m*6+n] * volumeFraction[0]);    

     // component 1 flux
     coords[0] = 0;
     coords[1] = distances[1][0];
     coords[2] = distances[1][1];
     get3PointsDerivativeInterpCoeffs(coords, 0, coeffs[1]);    

     m_pSolver->Add_A(intfcIndex[1], intfcIndex[2], - beta[1]*intfcArea*coeffs[1][0] * volumeFraction[1]);
     for(m=0; m<2; m++)    
	  for(n=0; n<6; n++)
	       m_pSolver->Add_A(intfcIndex[1], indices[1][m*6+n], - beta[1]*intfcArea*coeffs[1][m+1]*coefficients[1][m*6+n] * volumeFraction[1]);

     // interface equation
     m_pSolver->Add_A(intfcIndex[2], intfcIndex[2], + beta[0]*intfcArea*coeffs[0][0]);
     for(m=0; m<2; m++)    
	  for(n=0; n<6; n++)
	       m_pSolver->Add_A(intfcIndex[2], indices[0][m*6+n], + beta[0]*intfcArea*coeffs[0][m+1]*coefficients[0][m*6+n]);  

     m_pSolver->Add_A(intfcIndex[2], intfcIndex[2], + beta[1]*intfcArea*coeffs[1][0]);
     for(m=0; m<2; m++)    
	  for(n=0; n<6; n++)
	       m_pSolver->Add_A(intfcIndex[2], indices[1][m*6+n], + beta[1]*intfcArea*coeffs[1][m+1]*coefficients[1][m*6+n]);   
 }


/*
 * Currently, only Neumann boundary condition (give) is delt with.
 */
void EBM3D_CARTESIAN::setMatrixForBoundary(int intfcIndex[2], double intfcCenter[3], 
					   double normal[3], int intfcComp[2], 
					   double intfcArea, double volumeFraction[2], 
					   int i, int j, int k)
{
     double flux, value, beta;   
    if(intfcIndex[0]>=0)
    {
	 printf("EBM3D_CARTESIAN::setMatrixForBoundary: ????\n");
	 flux = m_pLaplace->getFlux(intfcCenter, intfcComp[0], normal);    
	 beta = m_pLaplace->getBeta(intfcCenter, intfcComp[0]);
	 m_pSolver->Add_b(intfcIndex[0], -flux*beta*intfcArea * volumeFraction[0]);    // negative because it is moved to the rhs.
    }
    if(intfcIndex[1]>=0)
    {
	 flux = - m_pLaplace->getFlux(intfcCenter, intfcComp[1], normal); // negative because it has negative normal
	 beta = m_pLaplace->getBeta(intfcCenter, intfcComp[1]);
	 m_pSolver->Add_b(intfcIndex[1], -flux*beta*intfcArea * volumeFraction[1]);    // negative because it is moved to the rhs.
    }
    //value = m_pLaplace->getExactSolution(intfcCenter, intfcComp[0]);  // assuming the solution is continuous at the interface.
    //m_pSolver->Add_A(intfcIndex[2], intfcIndex[2], 1);
    //m_pSolver->Add_b(intfcIndex[2], value); 
}

/*
 * get the flux through the partial face of a cartesian cell.
 * This function is only 1st order accurate.
 * see also EBM3D_CARTESIAN::getFluxStencilOffsetByPlane_partialFace().
 */
void EBM3D_CARTESIAN::getFluxStencilForPartialFace(int face, double center[3], int comp,
				    std::vector< std::pair<int,double> > & stencil,
				    int i, int j, int k)
{
    double dh;
    int I, J, K;
    I = i;
    J = j;
    K = k;
    
    switch(face-face%2)
    {
        case PLANE_X:
            if(face==FACE_WEST)     // west      
                I = i-1;
            else                    // east
                I = i+1;            
	    dh = getCellEdgeLength(0,i,j,k);
            break;
        case PLANE_Y:
            if(face==FACE_SOUTH)    // south
                J = j-1;
            else                    // north
                J = j+1;
	    dh = getCellEdgeLength(1,i,j,k);
            break;
        case PLANE_Z:
            if(face==FACE_DOWN)     // down
                K = k-1;
            else                    // up                
                K = k+1;            
	    dh = getCellEdgeLength(2,i,j,k);
            break;           
    }
    
    stencil.push_back(std::pair<int,double>(getCellUnknownIndex(comp,i,j,k), -1/dh));
    stencil.push_back(std::pair<int,double>(getCellUnknownIndex(comp,I,J,K),  1/dh));
}


/*
 * find the interpolation stencil offset for the flux on the plane.
 *
 * input parameter:
 *      plane: PLANE_X, PLANE_Y, PLANE_Z.
 *      center[3]:  partial face center.
 *      comp: comp of the partial face.
 *      i,j,k: cell index.
 * output parameter:
 *      stencil:
 *         each element is a pair of (index, coeff) to be used for the flux interpolation on the plane.
 */
void EBM3D_CARTESIAN::getFluxStencilForPartialFace2(
     int face, double center[3], int comp, 
     std::vector< std::pair<int,double> > &stencil ,
     int i, int j, int k)
{
     double faceCenter[3];
     double x, y, dx, dy, x1, y1, x2, y2, dh;
    
     getCellFaceCenter(face, faceCenter, i, j, k);    

     EBM_PLANE plane;
     int I, J, K;
     I = i;
     J = j;
     K = k;
    
    switch(face-face%2)
    {
        case PLANE_X:
            plane = PLANE_X;
            if(face==FACE_WEST)     // west      
                I = i-1;
            else                    // east
                I = i+1;            
            
            x1 = faceCenter[1];
            y1 = faceCenter[2];
            x = center[1];
            y = center[2];
            dx = getCellEdgeLength(1,i,j,k);
            dy = getCellEdgeLength(2,i,j,k);
            dh = getCellEdgeLength(0,i,j,k);
            break;
        case PLANE_Y:
            plane = PLANE_Y;
            if(face==FACE_SOUTH)    // south
                J = j-1;
            else                    // north
                J = j+1;
            
            x1 = faceCenter[0];
            y1 = faceCenter[2];
            x = center[0];
            y = center[2];
            dx = getCellEdgeLength(0,i,j,k);
            dy = getCellEdgeLength(2,i,j,k);
            dh = getCellEdgeLength(1,i,j,k);
            break;
        case PLANE_Z:
            plane = PLANE_Z;
            if(face==FACE_DOWN)     // down
                K = k-1;
            else                    // up                
                K = k+1;            
                    
            x1 = faceCenter[0];
            y1 = faceCenter[1];
            x = center[0];
            y = center[1];
            dx = getCellEdgeLength(0,i,j,k);            
            dy = getCellEdgeLength(1,i,j,k);
            dh = getCellEdgeLength(2,i,j,k);
            break;           
    }
      
    int index, IIndex[4][2];
    double coeff[4], denominator;       
    int corner = getCell(i,j,k)->getCellFaceNearestCorner(face, comp, center);
    
    switch(corner)
    {
    case 2:          //if(x>=x1 && y>=y1)        // north-east
	 x2 = x + dx;
	 y2 = y + dy;        
	 IIndex[0][0] = 0;   IIndex[0][1] = 0;
	 IIndex[1][0] = 1;   IIndex[1][1] = 0;
	 IIndex[2][0] = 1;   IIndex[2][1] = 1;
	 IIndex[3][0] = 0;   IIndex[3][1] = 1;          
	 break;
    case 3:          //else if(x<=x1 && y>=y1)   // north-west
	 x2 = x - dx;
	 y2 = y + dy;
	 IIndex[0][0] = 0;   IIndex[0][1] = 0;
	 IIndex[1][0] =-1;   IIndex[1][1] = 0;
	 IIndex[2][0] =-1;   IIndex[2][1] = 1;
	 IIndex[3][0] = 0;   IIndex[3][1] = 1;          
	 break;
    case 0:          //else if(x<=x1 && y<=y1)   // south-west
	 x2 = x - dx;
	 y2 = y - dy;
	 IIndex[0][0] = 0;   IIndex[0][1] = 0;
	 IIndex[1][0] =-1;   IIndex[1][1] = 0;
	 IIndex[2][0] =-1;   IIndex[2][1] =-1;
	 IIndex[3][0] = 0;   IIndex[3][1] =-1;          
	 break;
    case 1:    //else //if(x>=x1 && y<=y1) // south-east
	 x2 = x + dx;
	 y2 = y - dy;
	 IIndex[0][0] = 0;   IIndex[0][1] = 0;
	 IIndex[1][0] = 1;   IIndex[1][1] = 0;
	 IIndex[2][0] = 1;   IIndex[2][1] =-1;
	 IIndex[3][0] = 0;   IIndex[3][1] =-1;          
	 break;
    default:
	 printf("EBM3D_CARTESIAN::getFluxStencilForPartialFace2: corner = %d\n", corner);
	 exit(0);
    }
    getBilinearInterpolationCoeff(x1,y1,x2,y2,x,y,coeff);
    
    stencil.push_back(std::pair<int,double>(getNeighborByPlane(plane, IIndex[0][0],IIndex[0][1],comp,i,j,k), -1/dh*coeff[0]));        
    stencil.push_back(std::pair<int,double>(getNeighborByPlane(plane, IIndex[1][0],IIndex[1][1],comp,i,j,k), -1/dh*coeff[1]));
    stencil.push_back(std::pair<int,double>(getNeighborByPlane(plane, IIndex[2][0],IIndex[2][1],comp,i,j,k), -1/dh*coeff[2]));
    stencil.push_back(std::pair<int,double>(getNeighborByPlane(plane, IIndex[3][0],IIndex[3][1],comp,i,j,k), -1/dh*coeff[3]));
    stencil.push_back(std::pair<int,double>(getNeighborByPlane(plane, IIndex[0][0],IIndex[0][1],comp,I,J,K),  1/dh*coeff[0]));        
    stencil.push_back(std::pair<int,double>(getNeighborByPlane(plane, IIndex[1][0],IIndex[1][1],comp,I,J,K),  1/dh*coeff[1]));
    stencil.push_back(std::pair<int,double>(getNeighborByPlane(plane, IIndex[2][0],IIndex[2][1],comp,I,J,K),  1/dh*coeff[2]));
    stencil.push_back(std::pair<int,double>(getNeighborByPlane(plane, IIndex[3][0],IIndex[3][1],comp,I,J,K),  1/dh*coeff[3]));
}


double EBM3D_CARTESIAN::getVolumeFraction(int index, std::map<int, double> &volumes)
{
     return m_cellVolume;

     std::map<int, double>::iterator iter;
     iter = volumes.begin();
     while(iter!=volumes.end())
     {
	  if(index==(*iter).first)
	       return (*iter).second/m_cellVolume;
	  iter++;
     }
     return 1;
}


/*
 * In the positive intfcNormal[] direction, find crossing of the coordinates lines; 
 * then find the 2 points stencil to interpolate the unknowns at the crossing.
 */
void EBM3D_CARTESIAN::getSolutionInterpStencil(double intfcCenter[3], double normal[3], int comp,
                std::vector<int> &indices, std::vector<double> &coefficients, std::vector<double> &distances,
                int i, int j, int k)
{    
     // 1) get the crossing along the normal direction with the target plane.
     double crx[3], dist;
     getStencilPointsAndDistance(intfcCenter, normal, comp, crx, dist, i,j,k);
     distances.push_back(dist);
    // 2) get the interpolation coeffs of the solution at the crossing.
    get4PointsStencil(crx, comp, indices, coefficients);
    
}

 /*
  * In the positive normal[] direction, find crossing of the coordinates lines; 
  * then find the 3 points stencil to interpolate the unknowns at the crossing.
  */
void EBM3D_CARTESIAN::getSolutionInterpStencil2(double intfcCenter[3], 
						double normal[3], 
						int comp,
						std::vector<int> &indices, 
						std::vector<double> &coefficients, 
						std::vector<double> &distances, 
						int i, int j, int k)
{    
     EBM_PLANE plane;
     int sign;

     // set the intersecting plane using the maximum of the normal
     double maximum = fabs(normal[0]);
     plane = PLANE_X;
     if(maximum<fabs(normal[1]))
     {
	  plane = PLANE_Y;
	  maximum = fabs(normal[1]);
     }
     if(maximum<fabs(normal[2]))
	  plane = PLANE_Z;
	  
     if(normal[plane/2]>0)
	 sign = 1;
     else
	 sign = -1;

     int targetI,targetJ, targetK;
     targetI = targetJ = targetK = -1;
     double t, crossing[3];
     int index[6];
     double coeffs[6];

     indices.clear();
     coefficients.clear();
     distances.clear();    
   
     for(int m=1; m<=2; m++)
     {
	  if(plane==PLANE_X)
	  {
	       targetI = i + sign * m;
	       crossing[0] = getCellEdgeCenter(0, targetI);
	       t = (crossing[0]-intfcCenter[0])/normal[0];
	       crossing[1] = t * normal[1] + intfcCenter[1];
	       crossing[2] = t * normal[2] + intfcCenter[2];
	  }
	  else if(plane==PLANE_Y)
	  {
	       targetJ = j + sign * m;
	       crossing[1] = getCellEdgeCenter(1, targetJ);
	       t = (crossing[1]-intfcCenter[1])/normal[1];
	       crossing[0] = t * normal[0] + intfcCenter[0];
	       crossing[2] = t * normal[2] + intfcCenter[2];
	  }
	  else
	  {
	       targetK = k + sign * m;
	       crossing[2] = getCellEdgeCenter(2, targetK);
	       t = (crossing[2]-intfcCenter[2])/normal[2];
	       crossing[0] = t * normal[0] + intfcCenter[0];
	       crossing[1] = t * normal[1] + intfcCenter[1];
	  }
	  locateCell(crossing, targetI, targetJ, targetK);
	  
	  get6PointsSolutionInterpCoeffs_Plane(plane,comp,crossing, index, coeffs, targetI,targetJ,targetK);
	  for(int n=0; n<6; n++)
	  {
	       indices.push_back(index[n]);
	       coefficients.push_back(coeffs[n]);
	  }
	  
	  distances.push_back(EBM3D_CELL::getDistance(crossing,intfcCenter));
     }
 }


/*
 * find a point along the normal[3] direction. the component of the point must be the same as comp.
 * 
 * OUTPUT:
 *    coords[][];
 *    distance;
 */
void EBM3D_CARTESIAN::getStencilPointsAndDistance(double crx[3], double normal[3], int comp, double coords[3], double &distance, int i, int j, int k)
{
     int stepsize = 10;                // to be optimized ???
     double min_dist, max_dist;
     min_dist = m_mindd;                     // potential error
     max_dist = sqrt(m_dx*m_dx + m_dy*m_dy + m_dz*m_dz);
     
     double dd = (max_dist - min_dist)/stepsize;

     int I,J,K, dI, dJ, dK;
     boolean bValid = YES;
     distance = min_dist;
     while(distance < 2*max_dist)      // to be optimized ???
     {
	  bValid = YES;
	       
	  coords[0] = crx[0] + normal[0] * distance;
	  coords[1] = crx[1] + normal[1] * distance;
	  coords[2] = crx[2] + normal[2] * distance;
	       
	  locateCell(coords, I,J,K);
	  dI = I-i;
	  dJ = J-j;
	  dK = K-k;
	  if(dI==0&&dJ==0&&dK==0 )    // try not to pick the same cell as (i,j,k).
	       bValid = NO;
	
	  if(bValid)
	       return;
	  distance += dd;
     }
}

/*
 * find 4 unknowns which are nearest to crx[3].
 *
 */
void EBM3D_CARTESIAN::get4PointsStencil(double crx[3], int comp, 
					std::vector<int> &indices, std::vector<double> &coefficients)
{
     // 27 neighbors in total
     static const int neighbors[][3] = {{0,0,0}, 
					{-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1},     // x,y,z
					{-1,-1,0},{1,-1,0},{-1,1,0},{1,1,0},                    // xy
					{-1,0,-1},{1,0,-1},{-1,0,1},{1,0,1},                    // xz
					{0,-1,-1},{0,1,-1},{0,-1,1},{0,1,1},                    // yz
					{-1,-1,-1},{1,-1,-1},{-1,1,-1},{-1,-1,1},{1,1,-1},{1,-1,1},{-1,1,1},{1,1,1}  //xyz
                                       };

     int I, J, K;
     // 1) find the cell (I,J,K) containing crx.
     locateCell(crx, I,J,K);

     int size = 0, index[4], pi, pj, pk;
     double coords[4][3], volume, area;
     boolean bValid = NO;

     // 2) start with (I,J,K) to find 4 neighboring unknowns with the same component as comp.
     for(int m=0; m<27; m++)
     {
	  pi = I + neighbors[m][0];
	  pj = J + neighbors[m][1];
	  pk = K + neighbors[m][2];

	  index[size] = getCellUnknownIndexByTrial(comp, pi,pj,pk);
	  if(index[size]<0)                          // not valid unknown.
	       continue;
	 
	  getCellCenter(coords[size], pi,pj,pk);
	  if(size==2)    // make sure that the 3 points are not on a line.
	  {
	       area = fabs(EBM3D_CELL::getTriArea(coords[0],coords[1],coords[2]));
	       if(area<m_mindd*m_mindd/10)
		    continue;
	  }
	  if(size==3)    // make sure that the 4 points are not on a plane.
	  {
	       volume = fabs(EBM3D_CELL::getTetraVolume(coords[0],coords[1],coords[2],coords[3]));
	       if(volume<m_cellVolume/10)          // a little arbitrary
		    continue;                          // the 4 points are located on one plane
	       else
	       {
		    bValid = YES;
		    break;
	       }
	  }
	  size++;		
     }
     if(!bValid)
     {
	  printf("EBM3D_CARTESIAN::get4PointsStencil: can't find 4 unknowns in the neighbor of (%d,%d,%d)! with comp = %d.\n",I,J,K, comp);
	  debug_printNeighborComp(I,J,K);
	  exit(0);
     }
     
     // get the interpolation coefficients
     double coeffs[4];
     get4PointsSolutionInterpCoeffs(coords, crx, coeffs);
     
     for(int m=0; m<4; m++)
     {
	  indices.push_back(index[m]);
	  coefficients.push_back(coeffs[m]);
     }
}
// find the crossing along the normal direction with the target plane
// important concept is:
//     target plane, normal direction.
void EBM3D_CARTESIAN::getCrossingAlongNormal(double coords[3], double normal[3], int dir, int sign, double crx[3], int i, int j, int k)
{
    // the target plane
	int I, J, K;
	I = i;     J = j;     K = k; 
	switch(dir)
	{
	case 0:
		 I += sign;
		 break;
	case 1:
		 J += sign;
		 break;
	case 2:
		 K += sign;
		 break;
	}
	double cellCenter[3], t;
	getCellCenter(cellCenter, I,J,K);
	
	// crx[3]
	t = (cellCenter[dir] - coords[dir]) / normal[dir];  // normal[dir] can't be zero.
	for(int m=0; m<3; m++)
		 crx[m] = coords[m] + t * normal[m];
}


// linear interpolation of the solution.
void EBM3D_CARTESIAN::get2PointsSolutionInterpCoeffs(double coords[2], double x, double coeffs[2])
{    
    double denominator = coords[1] - coords[0];
    coeffs[0] = (coords[1] - x)/denominator;
    coeffs[1] = (x - coords[0])/denominator;
}

// linear interpolation of the derivative
void EBM3D_CARTESIAN::get2PointsDerivativeInterpCoeffs(double coords[2], double x, double coeffs[2])
{
    double denominator = coords[1] - coords[0];
    coeffs[0] = -1/denominator;
    coeffs[1] =  1/denominator;
}

/*
 * Give 3 points on the x coordinates coords[0], coords[1], coords[2], get the
 * interpolation coefficients for the point with coordinate x.
 * 
 * the same as the function in EBM_CARTESIAN.
 */ 
void EBM3D_CARTESIAN::get3PointsSolutionInterpCoeffs(double coords[3], double x, double coeffs[3])
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
 *
 * the same as the function in EBM_CARTESIAN.
 */ 
void EBM3D_CARTESIAN::get3PointsDerivativeInterpCoeffs(double coords[3], double x, double coeffs[3])
{
    double y0,y1,y, denominator;
    y0 = coords[1] - coords[0];
    y1 = coords[2] - coords[0];
    y  = x - coords[0];
    coeffs[0] = -(-2*y+y0+y1)/(y0*y1);    
    coeffs[1] = (2*y-y1)/(y0*(y0-y1));
    coeffs[2] = (y0-2*y)/(y1*(y0-y1));
}    

// bilinear
void EBM3D_CARTESIAN::get3PointsSolutionInterpCoeffs_Plane(double coords[3][2], double x[2], double coeffs[3])
{
    double denominator =  coords[2][0]*(coords[1][1]-coords[0][1])
                        + coords[1][0]*(coords[0][1]-coords[2][1])
                        + coords[0][0]*(coords[2][1]-coords[1][1]);
    coeffs[0] = (  coords[2][0]*(coords[1][1]-x[1])
                 + coords[1][0]*(x[1]        -coords[2][1])
                 + x[0]        *(coords[2][1]-coords[1][1]) ) / denominator;
    coeffs[1] =-(  coords[2][0]*(coords[0][1]-x[1])
                 + coords[0][0]*(x[1]        -coords[2][1])
                 + x[0]        *(coords[2][1]-coords[0][1]) ) / denominator;
    coeffs[2] = (  coords[1][0]*(coords[0][1]-x[1])
                 + coords[0][0]*(x[1]        -coords[1][1])
                 + x[0]        *(coords[1][1]-coords[0][1]) ) / denominator;
}


/*
 *
 */
void EBM3D_CARTESIAN::get6PointsSolutionInterpCoeffs_Plane(EBM_PLANE plane, 
							   int comp, 
							   double crx[3],
							   int index[6], 
							   double coeffs[6], 
							   int i, int j, int k)
{
     boolean selected[3][3] = {{NO,NO,NO},{NO,NO,NO},{NO,NO,NO}};
     boolean row_selected[3] = {NO,NO,NO};
     double coords[6][3], plane_crx[3], plane_coords[6][2];
     
     int x, y, total = 0;
     EBM3D_INDEX tri;

     // select 6 candiate neighors
     // first, by column
     for(x=-1; x<=1; x++)
	  for(y=-1; y<=1; y++)
	       if((index[total]=getNeighborByPlane(plane,x,y,comp,tri.m_index,i,j,k))>=0)
	       {
		    selected[x+1][y+1] = YES;
		    row_selected[y+1] = YES;
		    getCellCenter(coords[total],tri.m_index[0],tri.m_index[1],tri.m_index[2]);
		    total++;
		    break;
	       }
     // next, by row
     for(y=-1; y<=1; y++)
     {
	  if(row_selected[y+1])
	       continue;
	  for(x=-1; x<=1; x++)
	  {
	       //if(selected[x+1][y+1])          // not necessary.
	       //	    continue;
	       if((index[total]=getNeighborByPlane(plane,x,y,comp,tri.m_index,i,j,k))>=0)
	       {
		    selected[x+1][y+1] = YES;
		    getCellCenter(coords[total],tri.m_index[0],tri.m_index[1],tri.m_index[2]);
		    total++;
		    break;
	       }
	  }
     }

     // total must be less than 6 for now.

     for(x=-1; x<=1; x++)
	  for(y=-1; y<=1; y++)
	  {
	       if(selected[x+1][y+1])
		    continue;
	       if((index[total]=getNeighborByPlane(plane,x,y,comp,tri.m_index,i,j,k))>=0)
	       {
		    getCellCenter(coords[total],tri.m_index[0],tri.m_index[1],tri.m_index[2]);
		    total++;
	       }

	       if(total==6)
		    goto LABEL_6_CANDIDATES;
	  }
     
LABEL_6_CANDIDATES:
     if(total==6)
     {
	  getCoordsByPlane(plane, crx, plane_crx);
	  for(x=0; x<6; x++)
	       getCoordsByPlane(plane, coords[x], plane_coords[x]);
	  get6PointsSolutionInterpCoeffs_Plane(plane_coords,plane_crx,coeffs);
     }
     else
     {
	  printf("EBM3D_CARTESIAN::get6PointsSolutionInterpCoeffs_Plane:"
		 " can't find enough candidate cell.\n");
	  printf("\t current cell: (%d,%d,%d), comp=%d, plane=%d, total=%d\n",
		 i,j,k,comp, plane, total);
	  printf("\t (m_i,m_j,m_k): (%d,%d,%d) \n", m_i,m_j,m_k);
	  debug_printNeighborComp(plane,i,j,k);
	  exit(0);
     }
}

// biquadratic
void EBM3D_CARTESIAN::get6PointsSolutionInterpCoeffs_Plane(double coords[6][2], 
							   double x[2], 
							   double coeffs[6])
{
    double X, Y;
    //PETSc solver;
    //solver.SetTol(10e-9);
    LAPACK solver;
    
    solver.Create(0, 5, 10, 0);
    for(int i=0; i<6; i++)
    {
        X = coords[i][0];
        Y = coords[i][1];
        
        solver.Set_A(0, i, 1);
        solver.Set_A(1, i, X);
        solver.Set_A(2, i, Y);        
        solver.Set_A(3, i, X*X);
        solver.Set_A(4, i, X*Y);
        solver.Set_A(5, i, Y*Y);
        
        //printf("%f %f %f %f %f %f \n",1.0,X,Y,X*X,X*Y,Y*Y);        
    }
    X = x[0];
    Y = x[1];
    solver.Set_b(0, 1);
    solver.Set_b(1, X);
    solver.Set_b(2, Y);
    solver.Set_b(3, X*X);
    solver.Set_b(4, X*Y);
    solver.Set_b(5, Y*Y); 
    
    solver.Solve();
    solver.Get_x(coeffs); 
    //printf("%f %f %f %f %f %f\n", 
//	   coeffs[0],coeffs[1],coeffs[2],coeffs[3],coeffs[4],coeffs[5]);
}

// trilinear
// Reference: mathematics reference (in Chinese).
void EBM3D_CARTESIAN::get4PointsSolutionInterpCoeffs(double coords[4][3], double x[3], double coeffs[4])
{
    double x1,x2,x3, y1,y2,y3,z1,z2,z3;
    double X1, X2, X3, Y1, Y2, Y3, Z1, Z2, Z3;
    x1 = coords[0][0] - coords[3][0];   x2 = coords[1][0] - coords[3][0];   x3 = coords[2][0] - coords[3][0];
    y1 = coords[0][1] - coords[3][1];   y2 = coords[1][1] - coords[3][1];   y3 = coords[2][1] - coords[3][1];
    z1 = coords[0][2] - coords[3][2];   z2 = coords[1][2] - coords[3][2];   z3 = coords[2][2] - coords[3][2];
    
    double V =(  x1*y2*z3 + x2*y3*z1 + x3*y1*z2 
               -(x1*y3*z2 + x2*y1*z3 + x3*y2*z1) );
    X1 = y2*z3 - y3*z2;                 X2 = y3*z1 - y1*z3;                 X3 = y1*z2 - y2*z1;
    Y1 = z2*x3 - z3*x2;                 Y2 = z3*x1 - z1*x3;                 Y3 = z1*x2 - z2*x1;
    Z1 = x2*y3 - x3*y2;                 Z2 = x3*y1 - x1*y3;                 Z3 = x1*y2 - x2*y1;
    
    coeffs[0] = 1.0/V * ( X1*(x[0]-coords[3][0]) + Y1*(x[1]-coords[3][1]) + Z1*(x[2]-coords[3][2]) );
    coeffs[1] = 1.0/V * ( X2*(x[0]-coords[3][0]) + Y2*(x[1]-coords[3][1]) + Z2*(x[2]-coords[3][2]) );
    coeffs[2] = 1.0/V * ( X3*(x[0]-coords[3][0]) + Y3*(x[1]-coords[3][1]) + Z3*(x[2]-coords[3][2]) );
    coeffs[3] = 1 - coeffs[0] - coeffs[1] - coeffs[2];
}

// quadratic
void EBM3D_CARTESIAN::get10PointsSolutionInterpCoeffs(double coords[10][3], double x[3], double coeffs[10])
{
    LAPACK solver;
    //solver.SetTol(10e-9);
    
    solver.Create(0, 9, 10, 0);
    for(int i=0; i<10; i++)
    {
        solver.Set_A(0, i, 1);
        solver.Set_A(1, i, coords[i][0]);
        solver.Set_A(2, i, coords[i][1]);
        solver.Set_A(3, i, coords[i][2]);
        solver.Set_A(4, i, coords[i][0]*coords[i][1]);
        solver.Set_A(5, i, coords[i][1]*coords[i][2]);
        solver.Set_A(6, i, coords[i][0]*coords[i][2]);
        solver.Set_A(7, i, coords[i][0]*coords[i][0]);
        solver.Set_A(8, i, coords[i][1]*coords[i][1]);
        solver.Set_A(9, i, coords[i][2]*coords[i][2]);
        printf("%f %f %f %f %f %f %f %f %f %f\n",
                1.0, coords[i][0], coords[i][1], coords[i][2],
                coords[i][0]*coords[i][1], coords[i][1]*coords[i][2], coords[i][0]*coords[i][2],
                coords[i][0]*coords[i][0], coords[i][1]*coords[i][1], coords[i][2]*coords[i][2]);        
    }
    
    solver.Solve();
    solver.Get_x(coeffs);    
}

void EBM3D_CARTESIAN::getBilinearInterpolationCoeff(double x1, double y1, double x2, double y2, 
                double x, double y,
                double coeff[4])
{
    double denominator = (x2-x1)*(y2-y1);
    coeff[0] = (x2-x)*(y2-y)/denominator;
    coeff[1] = (x-x1)*(y2-y)/denominator;
    coeff[2] = (x-x1)*(y-y1)/denominator;
    coeff[3] = (x2-x)*(y-y1)/denominator;
}

/*
 * x, y can have the same address.
 */
void EBM3D_CARTESIAN::getCoordsByPlane(int plane, double x[3], double y[2])
{
     switch(plane)
     {
     case PLANE_X:
	  y[0] = x[1];	       y[1] = x[2];
	  break;
     case PLANE_Y:
	  y[0] = x[0];	       y[1] = x[2];
	  break;
     case PLANE_Z:
	  y[0] = x[0];	       y[1] = x[1];
	  break;
     }
}

/*
 * set m_lbuf[], m_ubuf[].
 * It is implicitly assumed that 
 *          m_lbuf[0]==m_ubuf[0] && m_lbuf[1]==m_ubuf[1] && m_lbuf[2]==m_ubuf[2]
 */
void EBM3D_CARTESIAN::pp_setup(void)
{
     m_lbuf[0] = m_ubuf[0] = 2;
     m_lbuf[1] = m_ubuf[1] = 2;
     m_lbuf[2] = m_ubuf[2] = 2;
}
/*
 * reset the unknown index at each cell to global index, excluding the buffered cell.
 *  1) calculate the range of unknown indices: [m_ilower, m_iupper].
 *  2) reset each cell's unknown index, excluding buffered cells.
 *
 * this function should be called after m_nLocalUnknowns is set.
 */                          
void EBM3D_CARTESIAN::pp_resetCellUnknownIndex(void)                              
{
     int i,j,k, *n_dist, num_nodes;

     // 1) calculate the range of unknown indices: [m_ilower, m_iupper];
     num_nodes = pp_numnodes();
     FT_VectorMemoryAlloc((POINTER*)&n_dist, num_nodes, sizeof(int));
     n_dist[m_myid] = m_nLocalUnknowns;
     pp_global_imax(n_dist, num_nodes);
     
     for(i=0; i<num_nodes; i++)
	  printf("ID %d: n_dist[%d]=%d \n", m_myid, i, n_dist[i]);

     m_ilower = 0;
     for(i=1; i<=m_myid; i++)
	  m_ilower += n_dist[i-1];
     m_iupper = m_ilower + m_nLocalUnknowns - 1;
     printf("ID %d: m_ilower=%d, m_iupper=%d\n", m_myid, m_ilower, m_iupper); 

     // 2) reset each cell's unknown index.
     for(k=m_lbuf[2]; k<m_lbuf[2]+m_gmax[2]; k++)
     for(j=m_lbuf[1]; j<m_lbuf[1]+m_gmax[1]; j++)      
     for(i=m_lbuf[0]; i<m_lbuf[0]+m_gmax[0]; i++)
	  getCell(i,j,k)->resetCellUnknownIndex(m_ilower);
     
     FT_FreeThese(1, n_dist);
}
/*
 * Reset the unknown index at each buffered cell to global index using 
 * information from the node's neighbors.
 */
void EBM3D_CARTESIAN::pp_resetBufferCellUnknownIndex(void)      
{
     int myid = pp_mynode();
     int **sendbuff, **recvbuff, size;

     // in each cell, there are 8 index.

     int length[] = {m_lbuf[0]+m_gmax[0]+m_ubuf[0],
		     m_lbuf[1]+m_gmax[1]+m_ubuf[1],
		     m_lbuf[2]+m_gmax[2]+m_ubuf[2]};
     int maxBuffSize = std::max( length[1]*length[2]*m_lbuf[0],
				 std::max(length[0]*length[2]*m_lbuf[1],
					  length[0]*length[1]*m_lbuf[2]) );
     size = maxBuffSize * 8;
          
     FT_MatrixMemoryAlloc((POINTER*)&sendbuff, 2, size, sizeof(int));
     FT_MatrixMemoryAlloc((POINTER*)&recvbuff, 2, size, sizeof(int));
     
     //printf("ID %d: size of sendbuff/recvbuff is (%d,%d). \n", pp_mynode(),2,size);

     for(int dir=0; dir<3; dir++)        // for X direction, then Y direction.
     {
	  switch(dir)
	  {
	  case 0:
	       size = length[1]*length[2]*m_lbuf[0] * 8;
	       break;
	  case 1:
	       size = length[0]*length[2]*m_lbuf[1] * 8;
	       break;
	  case 2:
	       size = length[0]*length[1]*m_lbuf[2] * 8;
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
void EBM3D_CARTESIAN::pp_packCellUnknownIndex(int dir, int size, int **data)
{
     int i, j, k, p;
     if(dir==0)        // WEST, EAST    
     {
	  p = 0;
	  for(k=0; k<m_lbuf[2]+m_gmax[2]+m_ubuf[2]; k++)
	  for(j=0; j<m_lbuf[1]+m_gmax[1]+m_ubuf[1]; j++)
	  for(i=m_lbuf[0]; i<m_lbuf[0]+m_lbuf[0]; i++)
	  {
	       getCell(i,j,k)->getCellUnknownIndex(&(data[0][p]));
	       p += 8;
	  }
	  //printf("ID %d: size=%d, p=%d\n", pp_mynode(),size,p);
	  p = 0;
	  for(k=0; k<m_lbuf[2]+m_gmax[2]+m_ubuf[2]; k++)
	  for(j=0; j<m_lbuf[1]+m_gmax[1]+m_ubuf[1]; j++)
	  for(i=m_lbuf[0]+m_gmax[0]-m_ubuf[0]; i<m_lbuf[0]+m_gmax[0]; i++)
	  {
	       getCell(i,j,k)->getCellUnknownIndex(&(data[1][p]));
	       p += 8;
	  }
	  //printf("ID %d: size=%d, p=%d\n", pp_mynode(),size,p);
     }
     else if(dir==1)     // SOUTH, NORTH
     {
	  p = 0;
	  for(k=0; k<m_lbuf[2]+m_gmax[2]+m_ubuf[2]; k++)
	  for(j=m_lbuf[1]; j<m_lbuf[1]+m_lbuf[1]; j++)
	  for(i=0; i<m_lbuf[0]+m_gmax[0]+m_ubuf[0]; i++)
	  {
	       getCell(i,j,k)->getCellUnknownIndex(&(data[0][p]));
	       p += 8;
	  }
	  //printf("ID %d: size=%d, p=%d\n", pp_mynode(),size,p);
	  p = 0;
	  for(k=0; k<m_lbuf[2]+m_gmax[2]+m_ubuf[2]; k++)
	  for(j=m_lbuf[1]+m_gmax[1]-m_ubuf[1]; j<m_lbuf[1]+m_gmax[1]; j++)
	  for(i=0; i<m_lbuf[0]+m_gmax[0]+m_ubuf[0]; i++)
	  {
	       getCell(i,j,k)->getCellUnknownIndex(&(data[1][p]));
	       p += 8;
	  }
	  //printf("ID %d: size=%d, p=%d\n", pp_mynode(),size,p);
     }
     else                 // DOWN, UP
     {
	  p = 0;
	  for(k=m_lbuf[2]; k<m_lbuf[2]+m_lbuf[2]; k++)
	  for(j=0; j<m_lbuf[1]+m_gmax[1]+m_ubuf[1]; j++)
	  for(i=0; i<m_lbuf[0]+m_gmax[0]+m_ubuf[0]; i++)
	  {
	       getCell(i,j,k)->getCellUnknownIndex(&(data[0][p]));
	       p += 8;
	  }
	  //printf("ID %d: size=%d, p=%d\n", pp_mynode(),size,p);
	  p = 0;
	  for(k=m_lbuf[2]+m_gmax[2]-m_ubuf[2]; k<m_lbuf[2]+m_gmax[2]; k++)
	  for(j=0; j<m_lbuf[1]+m_gmax[1]+m_ubuf[1]; j++)
	  for(i=0; i<m_lbuf[0]+m_gmax[0]+m_ubuf[0]; i++)
	  {
	       getCell(i,j,k)->getCellUnknownIndex(&(data[1][p]));
	       p += 8;
	  }
	  //printf("ID %d: size=%d, p=%d\n", pp_mynode(),size,p);
     }
}

void EBM3D_CARTESIAN::pp_unpackCellUnknownIndex(int dir, int size, int **data)
{
          int i, j, k, p;
     if(dir==0)        // WEST, EAST    
     {
	  p = 0;
	  for(k=0; k<m_lbuf[2]+m_gmax[2]+m_ubuf[2]; k++)
	  for(j=0; j<m_lbuf[1]+m_gmax[1]+m_ubuf[1]; j++)
	  for(i=0; i<m_lbuf[0]; i++)
	  {
	       getCell(i,j,k)->setCellUnknownIndex(&(data[0][p]));
	       p += 8;
	  }
	  //printf("ID %d: size=%d, p=%d\n", pp_mynode(),size,p);
	  p = 0;
	  for(k=0; k<m_lbuf[2]+m_gmax[2]+m_ubuf[2]; k++)
	  for(j=0; j<m_lbuf[1]+m_gmax[1]+m_ubuf[1]; j++)
	  for(i=m_lbuf[0]+m_gmax[0]; i<m_lbuf[0]+m_gmax[0]+m_ubuf[0]; i++)
	  {
	       getCell(i,j,k)->setCellUnknownIndex(&(data[1][p]));
	       p += 8;
	  }
	  //printf("ID %d: size=%d, p=%d\n", pp_mynode(),size,p);
     }
     else if(dir==1)     // SOUTH, NORTH
     {
	  p = 0;
	  for(k=0; k<m_lbuf[2]+m_gmax[2]+m_ubuf[2]; k++)
	  for(j=0; j<m_lbuf[1]; j++)
	  for(i=0; i<m_lbuf[0]+m_gmax[0]+m_ubuf[0]; i++)
	  {
	       getCell(i,j,k)->setCellUnknownIndex(&(data[0][p]));
	       p += 8;
	  }
	  //printf("ID %d: size=%d, p=%d\n", pp_mynode(),size,p);
	  p = 0;
	  for(k=0; k<m_lbuf[2]+m_gmax[2]+m_ubuf[2]; k++)
	  for(j=m_lbuf[1]+m_gmax[1]; j<m_lbuf[1]+m_gmax[1]+m_ubuf[1]; j++)
	  for(i=0; i<m_lbuf[0]+m_gmax[0]+m_ubuf[0]; i++)
	  {
	       getCell(i,j,k)->setCellUnknownIndex(&(data[1][p]));
	       p += 8;
	  }
	  //printf("ID %d: size=%d, p=%d\n", pp_mynode(),size,p);
     }
     else                 // DOWN, UP
     {
	  p = 0;
	  for(k=0; k<m_lbuf[2]; k++)
	  for(j=0; j<m_lbuf[1]+m_gmax[1]+m_ubuf[1]; j++)
	  for(i=0; i<m_lbuf[0]+m_gmax[0]+m_ubuf[0]; i++)
	  {
	       getCell(i,j,k)->setCellUnknownIndex(&(data[0][p]));
	       p += 8;
	  }
	  //printf("ID %d: size=%d, p=%d\n", pp_mynode(),size,p);
	  p = 0;
	  for(k=m_lbuf[2]+m_gmax[2]; k<m_lbuf[2]+m_gmax[2]+m_ubuf[2]; k++)
	  for(j=0; j<m_lbuf[1]+m_gmax[1]+m_ubuf[1]; j++)
	  for(i=0; i<m_lbuf[0]+m_gmax[0]+m_ubuf[0]; i++)
	  {
	       getCell(i,j,k)->setCellUnknownIndex(&(data[1][p]));
	       p += 8;
	  }
	  //printf("ID %d: size=%d, p=%d\n", pp_mynode(),size,p);
     }
}					  

void EBM3D_CARTESIAN::pp_resetBufferCellCompList(void)
{
     int **sendbuff, **recvbuff, size;

     // in each cell, there are 2 unknowns.
     int length[] = {m_lbuf[0]+m_gmax[0]+m_ubuf[0],
		     m_lbuf[1]+m_gmax[1]+m_ubuf[1],
		     m_lbuf[2]+m_gmax[2]+m_ubuf[2]};
     int maxBuffSize = std::max( length[1]*length[2]*m_lbuf[0],
				 std::max(length[0]*length[2]*m_lbuf[1],
					  length[0]*length[1]*m_lbuf[2]) );
     size = maxBuffSize * 2;
     FT_MatrixMemoryAlloc((POINTER*)&sendbuff, 2, size, sizeof(int));
     FT_MatrixMemoryAlloc((POINTER*)&recvbuff, 2, size, sizeof(int));
     
     for(int dir=0; dir<2; dir++)
     {
	  switch(dir)
	  {
	  case 0:
	       size = length[1]*length[2]*m_lbuf[0] * 2;
	       break;
	  case 1:
	       size = length[0]*length[2]*m_lbuf[1] * 2;
	       break;
	  case 2:
	       size = length[0]*length[1]*m_lbuf[2] * 2;
	       break;
	  }
	  	  
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
void EBM3D_CARTESIAN::pp_packCellCompList(int dir, int size, int **data)
{
     int i, j, k, p;
     if(dir==0)        // WEST, EAST    
     {
	  p = 0;
	  for(k=0; k<m_lbuf[2]+m_gmax[2]+m_ubuf[2]; k++)
	  for(j=0; j<m_lbuf[1]+m_gmax[1]+m_ubuf[1]; j++)
	  for(i=m_lbuf[0]; i<m_lbuf[0]+m_lbuf[0]; i++)
	  {
	       getCellCompList(&data[0][p],i,j,k);
	       p += 2;
	  }
	  //printf("ID %d: size=%d, p=%d\n", pp_mynode(),size,p);
	  p = 0;
	  for(k=0; k<m_lbuf[2]+m_gmax[2]+m_ubuf[2]; k++)
	  for(j=0; j<m_lbuf[1]+m_gmax[1]+m_ubuf[1]; j++)
	  for(i=m_lbuf[0]+m_gmax[0]-m_ubuf[0]; i<m_lbuf[0]+m_gmax[0]; i++)
	  {
	       getCellCompList(&data[1][p],i,j,k);
	       p += 2;
	  }
	  //printf("ID %d: size=%d, p=%d\n", pp_mynode(),size,p);
     }
     else if(dir==1)     // SOUTH, NORTH
     {
	  p = 0;
	  for(k=0; k<m_lbuf[2]+m_gmax[2]+m_ubuf[2]; k++)
	  for(j=m_lbuf[1]; j<m_lbuf[1]+m_lbuf[1]; j++)
	  for(i=0; i<m_lbuf[0]+m_gmax[0]+m_ubuf[0]; i++)
	  {
	       getCellCompList(&data[0][p],i,j,k);
	       p += 2;
	  }
	  //printf("ID %d: size=%d, p=%d\n", pp_mynode(),size,p);
	  p = 0;
	  for(k=0; k<m_lbuf[2]+m_gmax[2]+m_ubuf[2]; k++)
	  for(j=m_lbuf[1]+m_gmax[1]-m_ubuf[1]; j<m_lbuf[1]+m_gmax[1]; j++)
	  for(i=0; i<m_lbuf[0]+m_gmax[0]+m_ubuf[0]; i++)
	  {
	       getCellCompList(&data[1][p],i,j,k);
	       p += 2;
	  }
	  //printf("ID %d: size=%d, p=%d\n", pp_mynode(),size,p);
     }
     else                 // DOWN, UP
     {
	  p = 0;
	  for(k=m_lbuf[2]; k<m_lbuf[2]+m_lbuf[2]; k++)
	  for(j=0; j<m_lbuf[1]+m_gmax[1]+m_ubuf[1]; j++)
	  for(i=0; i<m_lbuf[0]+m_gmax[0]+m_ubuf[0]; i++)
	  {
	       getCellCompList(&data[0][p],i,j,k);
	       p += 2;
	  }
	  //printf("ID %d: size=%d, p=%d\n", pp_mynode(),size,p);
	  p = 0;
	  for(k=m_lbuf[2]+m_gmax[2]-m_ubuf[2]; k<m_lbuf[2]+m_gmax[2]; k++)
	  for(j=0; j<m_lbuf[1]+m_gmax[1]+m_ubuf[1]; j++)
	  for(i=0; i<m_lbuf[0]+m_gmax[0]+m_ubuf[0]; i++)
	  {
	       getCellCompList(&data[1][p],i,j,k);
	       p += 2;
	  }
	  //printf("ID %d: size=%d, p=%d\n", pp_mynode(),size,p);
     }
}
void EBM3D_CARTESIAN::pp_unpackCellCompList(int dir, int size, int **data)
{
         int i, j, k, p;
     if(dir==0)        // WEST, EAST    
     {
	  p = 0;
	  for(k=0; k<m_lbuf[2]+m_gmax[2]+m_ubuf[2]; k++)
	  for(j=0; j<m_lbuf[1]+m_gmax[1]+m_ubuf[1]; j++)
	  for(i=0; i<m_lbuf[0]; i++)
	  {
	       setCellCompList(&data[0][p],i,j,k);
	       p += 2;
	  }
	  //printf("ID %d: size=%d, p=%d\n", pp_mynode(),size,p);
	  p = 0;
	  for(k=0; k<m_lbuf[2]+m_gmax[2]+m_ubuf[2]; k++)
	  for(j=0; j<m_lbuf[1]+m_gmax[1]+m_ubuf[1]; j++)
	  for(i=m_lbuf[0]+m_gmax[0]; i<m_lbuf[0]+m_gmax[0]+m_ubuf[0]; i++)
	  {
	       setCellCompList(&data[1][p],i,j,k);
	       p += 2;
	  }
	  //printf("ID %d: size=%d, p=%d\n", pp_mynode(),size,p);
     }
     else if(dir==1)     // SOUTH, NORTH
     {
	  p = 0;
	  for(k=0; k<m_lbuf[2]+m_gmax[2]+m_ubuf[2]; k++)
	  for(j=0; j<m_lbuf[1]; j++)
	  for(i=0; i<m_lbuf[0]+m_gmax[0]+m_ubuf[0]; i++)
	  {
	       setCellCompList(&data[0][p],i,j,k);
	       p += 2;
	  }
	  //printf("ID %d: size=%d, p=%d\n", pp_mynode(),size,p);
	  p = 0;
	  for(k=0; k<m_lbuf[2]+m_gmax[2]+m_ubuf[2]; k++)
	  for(j=m_lbuf[1]+m_gmax[1]; j<m_lbuf[1]+m_gmax[1]+m_ubuf[1]; j++)
	  for(i=0; i<m_lbuf[0]+m_gmax[0]+m_ubuf[0]; i++)
	  {
	       setCellCompList(&data[1][p],i,j,k);
	       p += 2;
	  }
	  //printf("ID %d: size=%d, p=%d\n", pp_mynode(),size,p);
     }
     else                 // DOWN, UP
     {
	  p = 0;
	  for(k=0; k<m_lbuf[2]; k++)
	  for(j=0; j<m_lbuf[1]+m_gmax[1]+m_ubuf[1]; j++)
	  for(i=0; i<m_lbuf[0]+m_gmax[0]+m_ubuf[0]; i++)
	  {
	       setCellCompList(&data[0][p],i,j,k);
	       p += 2;
	  }
	  //printf("ID %d: size=%d, p=%d\n", pp_mynode(),size,p);
	  p = 0;
	  for(k=m_lbuf[2]+m_gmax[2]; k<m_lbuf[2]+m_gmax[2]+m_ubuf[2]; k++)
	  for(j=0; j<m_lbuf[1]+m_gmax[1]+m_ubuf[1]; j++)
	  for(i=0; i<m_lbuf[0]+m_gmax[0]+m_ubuf[0]; i++)
	  {
	       setCellCompList(&data[1][p],i,j,k);
	       p += 2;
	  }
	  //printf("ID %d: size=%d, p=%d\n", pp_mynode(),size,p);
     }
}


/*
 * reset the unknown for the buffered cells.
 */
void EBM3D_CARTESIAN::pp_resetBufferCellUnknown(void)
{
     double **sendbuff, **recvbuff, size;

     // in each cell, there are 2 unknowns.
     int length[] = {m_lbuf[0]+m_gmax[0]+m_ubuf[0],
		     m_lbuf[1]+m_gmax[1]+m_ubuf[1],
		     m_lbuf[2]+m_gmax[2]+m_ubuf[2]};
     int maxBuffSize = std::max( length[1]*length[2]*m_lbuf[0],
				 std::max(length[0]*length[2]*m_lbuf[1],
					  length[0]*length[1]*m_lbuf[2]) );

     size = maxBuffSize * 2;

     FT_MatrixMemoryAlloc((POINTER*)&sendbuff, 2, size, sizeof(double));
     FT_MatrixMemoryAlloc((POINTER*)&recvbuff, 2, size, sizeof(double));
     
     for(int dir=0; dir<2; dir++)
     {
	  switch(dir)
	  {
	  case 0:
	       size = length[1]*length[2]*m_lbuf[0] * 2;
	       break;
	  case 1:
	       size = length[0]*length[2]*m_lbuf[1] * 2;
	       break;
	  case 2:
	       size = length[0]*length[1]*m_lbuf[2] * 2;
	       break;
	  }

	  pp_packCellUnknown(dir, size, sendbuff);
	  FronTierPPGrid::sendrecv(dir, m_pFront, (void**)sendbuff, (void**)recvbuff, size, MPI_DOUBLE);
	  pp_unpackCellUnknown(dir, size, recvbuff);
     }
     FT_FreeThese(2, sendbuff, recvbuff);
}
void EBM3D_CARTESIAN::pp_packCellUnknown(int dir, int size, double **data)
{
     int i, j, k, p;
     if(dir==0)        // WEST, EAST    
     {
	  p = 0;
	  for(k=0; k<m_lbuf[2]+m_gmax[2]+m_ubuf[2]; k++)
	  for(j=0; j<m_lbuf[1]+m_gmax[1]+m_ubuf[1]; j++)
	  for(i=m_lbuf[0]; i<m_lbuf[0]+m_lbuf[0]; i++)
	  {
	       getCellUnknown(&data[0][p],i,j,k);
	       p += 2;
	  }
	  //printf("ID %d: size=%d, p=%d\n", pp_mynode(),size,p);
	  p = 0;
	  for(k=0; k<m_lbuf[2]+m_gmax[2]+m_ubuf[2]; k++)
	  for(j=0; j<m_lbuf[1]+m_gmax[1]+m_ubuf[1]; j++)
	  for(i=m_lbuf[0]+m_gmax[0]-m_ubuf[0]; i<m_lbuf[0]+m_gmax[0]; i++)
	  {
	       getCellUnknown(&data[1][p],i,j,k);
	       p += 2;
	  }
	  //printf("ID %d: size=%d, p=%d\n", pp_mynode(),size,p);
     }
     else if(dir==1)     // SOUTH, NORTH
     {
	  p = 0;
	  for(k=0; k<m_lbuf[2]+m_gmax[2]+m_ubuf[2]; k++)
	  for(j=m_lbuf[1]; j<m_lbuf[1]+m_lbuf[1]; j++)
	  for(i=0; i<m_lbuf[0]+m_gmax[0]+m_ubuf[0]; i++)
	  {
	       getCellUnknown(&data[0][p],i,j,k);
	       p += 2;
	  }
	  //printf("ID %d: size=%d, p=%d\n", pp_mynode(),size,p);
	  p = 0;
	  for(k=0; k<m_lbuf[2]+m_gmax[2]+m_ubuf[2]; k++)
	  for(j=m_lbuf[1]+m_gmax[1]-m_ubuf[1]; j<m_lbuf[1]+m_gmax[1]; j++)
	  for(i=0; i<m_lbuf[0]+m_gmax[0]+m_ubuf[0]; i++)
	  {
	       getCellUnknown(&data[1][p],i,j,k);
	       p += 2;
	  }
	  //printf("ID %d: size=%d, p=%d\n", pp_mynode(),size,p);
     }
     else                 // DOWN, UP
     {
	  p = 0;
	  for(k=m_lbuf[2]; k<m_lbuf[2]+m_lbuf[2]; k++)
	  for(j=0; j<m_lbuf[1]+m_gmax[1]+m_ubuf[1]; j++)
	  for(i=0; i<m_lbuf[0]+m_gmax[0]+m_ubuf[0]; i++)
	  {
	       getCellUnknown(&data[0][p],i,j,k);
	       p += 2;
	  }
	  //printf("ID %d: size=%d, p=%d\n", pp_mynode(),size,p);
	  p = 0;
	  for(k=m_lbuf[2]+m_gmax[2]-m_ubuf[2]; k<m_lbuf[2]+m_gmax[2]; k++)
	  for(j=0; j<m_lbuf[1]+m_gmax[1]+m_ubuf[1]; j++)
	  for(i=0; i<m_lbuf[0]+m_gmax[0]+m_ubuf[0]; i++)
	  {
	       getCellUnknown(&data[1][p],i,j,k);
	       p += 2;
	  }
	  //printf("ID %d: size=%d, p=%d\n", pp_mynode(),size,p);
     }
}
void EBM3D_CARTESIAN::pp_unpackCellUnknown(int dir, int size, double **data)
{
     int i, j, k, p;
     if(dir==0)        // WEST, EAST    
     {
	  p = 0;
	  for(k=0; k<m_lbuf[2]+m_gmax[2]+m_ubuf[2]; k++)
	  for(j=0; j<m_lbuf[1]+m_gmax[1]+m_ubuf[1]; j++)
	  for(i=0; i<m_lbuf[0]; i++)
	  {
	       setCellUnknown(&data[0][p],i,j,k);
	       p += 2;
	  }
	  //printf("ID %d: size=%d, p=%d\n", pp_mynode(),size,p);
	  p = 0;
	  for(k=0; k<m_lbuf[2]+m_gmax[2]+m_ubuf[2]; k++)
	  for(j=0; j<m_lbuf[1]+m_gmax[1]+m_ubuf[1]; j++)
	  for(i=m_lbuf[0]+m_gmax[0]; i<m_lbuf[0]+m_gmax[0]+m_ubuf[0]; i++)
	  {
	       setCellUnknown(&data[1][p],i,j,k);
	       p += 2;
	  }
	  //printf("ID %d: size=%d, p=%d\n", pp_mynode(),size,p);
     }
     else if(dir==1)     // SOUTH, NORTH
     {
	  p = 0;
	  for(k=0; k<m_lbuf[2]+m_gmax[2]+m_ubuf[2]; k++)
	  for(j=0; j<m_lbuf[1]; j++)
	  for(i=0; i<m_lbuf[0]+m_gmax[0]+m_ubuf[0]; i++)
	  {
	       setCellUnknown(&data[0][p],i,j,k);
	       p += 2;
	  }
	  //printf("ID %d: size=%d, p=%d\n", pp_mynode(),size,p);
	  p = 0;
	  for(k=0; k<m_lbuf[2]+m_gmax[2]+m_ubuf[2]; k++)
	  for(j=m_lbuf[1]+m_gmax[1]; j<m_lbuf[1]+m_gmax[1]+m_ubuf[1]; j++)
	  for(i=0; i<m_lbuf[0]+m_gmax[0]+m_ubuf[0]; i++)
	  {
	       setCellUnknown(&data[1][p],i,j,k);
	       p += 2;
	  }
	  //printf("ID %d: size=%d, p=%d\n", pp_mynode(),size,p);
     }
     else                 // DOWN, UP
     {
	  p = 0;
	  for(k=0; k<m_lbuf[2]; k++)
	  for(j=0; j<m_lbuf[1]+m_gmax[1]+m_ubuf[1]; j++)
	  for(i=0; i<m_lbuf[0]+m_gmax[0]+m_ubuf[0]; i++)
	  {
	       setCellUnknown(&data[0][p],i,j,k);
	       p += 2;
	  }
	  //printf("ID %d: size=%d, p=%d\n", pp_mynode(),size,p);
	  p = 0;
	  for(k=m_lbuf[2]+m_gmax[2]; k<m_lbuf[2]+m_gmax[2]+m_ubuf[2]; k++)
	  for(j=0; j<m_lbuf[1]+m_gmax[1]+m_ubuf[1]; j++)
	  for(i=0; i<m_lbuf[0]+m_gmax[0]+m_ubuf[0]; i++)
	  {
	       setCellUnknown(&data[1][p],i,j,k);
	       p += 2;
	  }
	  //printf("ID %d: size=%d, p=%d\n", pp_mynode(),size,p);
     }
}

/*
 * Setup the data structure.        
 * 
 */
void EBM3D_CARTESIAN::setup(void)
{
     m_myid = pp_mynode();
     
     int i,j,k;
     

    // the following init is moved here because EBM3D_CARTESIAN::getComp() need m_mindd.
    m_i = m_j = m_k = -1;
    m_dx = getCellEdgeLength(0,0,0,0);
    m_dy = getCellEdgeLength(1,0,0,0);
    m_dz = getCellEdgeLength(2,0,0,0);
    m_mindd = std::min(std::min(m_dx,m_dy), m_dz);
    m_maxdd = std::max(std::max(m_dx,m_dy), m_dz);
    m_cellVolume = m_dx*m_dy*m_dz;
       
    // m_gmax[]
    RECT_GRID *rect_grid = m_pFront->rect_grid;
    m_gmax[0] = rect_grid->gmax[0];
    m_gmax[1] = rect_grid->gmax[1];
    m_gmax[2] = rect_grid->gmax[2];
    
    pp_setup();

    int length[] = {m_lbuf[0]+m_gmax[0]+m_ubuf[0], 
		    m_lbuf[1]+m_gmax[1]+m_ubuf[1],
		    m_lbuf[2]+m_gmax[2]+m_ubuf[2]};

    // m_cells
    FT_TriArrayMemoryAlloc((POINTER*)&m_cells, 
	      length[0], length[1], length[2],
	      sizeof(EBM3D_CELL));
    for(k=0; k<length[2]; k++)
	 for(j=0; j<length[1]; j++)
	      for(i=0; i<length[0]; i++)
		   m_cells[i][j][k].init();    // call the constructor

    
    // m_comps
    double p[3];
    m_nComps = (length[0]+1)*(length[1]+1)*(length[2]+1);
    FT_VectorMemoryAlloc((POINTER*)&m_comps, m_nComps,sizeof(int));
    for(k=0; k<=length[2]; k++)
    for(j=0; j<=length[1]; j++)
    for(i=0; i<=length[0]; i++)   
    {
	 getCellCornerCoords(0,p,i,j,k);
	 m_comps[getCellCornerIndex(0,i,j,k)] = getComp(p);
    }
    
    // the following is used to modify the components from 0 to 2
    // so that we only have components 1, 2 left.
    setEquivalentComponents();  
    
    // m_values
    m_nLocalUnknowns = 0;
    int numCellUnknowns;
    int comp[9];        // comp[8] for the comp of the center;
    
    CELL_TYPE type;
    
    for(k=m_lbuf[2]; k<m_lbuf[2]+m_gmax[2]; k++)
    for(j=m_lbuf[1]; j<m_lbuf[1]+m_gmax[1]; j++)
    for(i=m_lbuf[0]; i<m_lbuf[0]+m_gmax[0]; i++)
    {            
        //numCellUnknowns = getCellUnknownNumber(comp,i,j,k);
        //getCell(i,j,k)->setCellUnknownIndex(m_nUnknowns, comp, numCellUnknowns);        
        getCellComp(comp,i,j,k);
        getCell(i,j,k)->setComp(comp);   
	getCell(i,j,k)->setCellType(comp);
	type = getCell(i,j,k)->getCellType();
		
	if(type==EXTERNAL)       // outside of the domain
	     continue;
	else              // INTERNAL, PARTIAL, BOUNDARY.
	{
	     numCellUnknowns = getCell(i,j,k)->setCellUnknownIndex(m_nLocalUnknowns);
	     numCellUnknowns = numCellUnknowns*2 - 1;
	     	   
	     m_nLocalUnknowns += numCellUnknowns;
	}
    }    
    
    for(int d=0; d<2; d++)
     {
	  FT_TriArrayMemoryAlloc((POINTER*)&m_unknowns[d],length[0],length[1],length[2],sizeof(double));
	  FT_TriArrayMemoryAlloc((POINTER*)&m_compList[d],length[0],length[1],length[2],sizeof(int)); 

	  for(k=0; k<length[2]; k++)
	       for(j=0; j<length[1]; j++)
		    for(i=0; i<length[0]; i++)
		    {
			 m_unknowns[d][i][j][k] = 0;
			 m_compList[d][i][j][k] = -1;
		    }
     }
 
}

// The following functions are used for the EBM methods.     
void EBM3D_CARTESIAN::setMatrix(void)
{    
    // this is very important to speed up the matrix setup!
    m_pSolver->Create(m_ilower, m_iupper, 36, 0);     // is this enough?
    
    for(int k=m_lbuf[2]; k<m_lbuf[2]+m_gmax[2]; k++)    
    for(int j=m_lbuf[1]; j<m_lbuf[1]+m_gmax[1]; j++)
    for(int i=m_lbuf[0]; i<m_lbuf[0]+m_gmax[0]; i++)   
    {   
	 setCurrentCell(i,j,k);      // this is good for debugging too.
	//int numCellUnknowns;
	//numCellUnknowns = getCell(i,j,k)->getCellUnknownNumber();
        //switch(getCell(i,j,k)->getCellUnknownNumber())

	switch(getCell(i,j,k)->getCellType())
        {
	case INTERNAL:
	     setMatrixForInternalCell(i,j,k);
	     break;
	case PARTIAL:  
	     setMatrixForPartialCell_2Unknowns(i,j,k);                
	     break;
	case BOUNDARY:
	     setMatrixForBoundaryCell(i,j,k);
	     break;
	//default:           // EXTERNAL
        }               
    }    
}

EBM3D_CARTESIAN::EBM3D_CARTESIAN(Front *front, SOLVER *solver, EBM3D_LAPLACE *laplace)
{
     m_pFront = front;
     m_pSolver = solver;
     m_pLaplace = laplace;
}
EBM3D_CARTESIAN::~EBM3D_CARTESIAN()
{
    FT_FreeThese(2, m_cells, m_comps);
    FT_FreeThese(4, m_unknowns[0], m_unknowns[1], m_compList[0], m_compList[1]);
}


void EBM3D_CARTESIAN::solve(void)
{
     setup();
     printf("ID %d: EBM3D_CARTESIAN::solve:\n", m_myid);
     printf("ID %d: \tm_gmax = {%d,%d,%d}\n", m_myid, m_gmax[0],m_gmax[1],m_gmax[2]);
     printf("ID %d: \tm_dx = %e, m_dy = %e, m_dz = %e\n", m_myid, m_dx, m_dy, m_dz);
     printf("ID %d: \tm_nLocalUnknowns = %d\n", m_myid, m_nLocalUnknowns);
     
     printf("ID %d: pp_resetCellUnknownIndex().\n", m_myid);
     pp_resetCellUnknownIndex();
     printf("ID %d: pp_resetBufferCellUnknownIndex().\n", m_myid);
     pp_resetBufferCellUnknownIndex();
     
     //debug_saveUnknownIndex();
     //return;
     //exit(0);

     printf("ID %d: \tsetMatrix()\n", m_myid);
     setMatrix();

     printf("ID %d: \tsolve()\n", m_myid);
     m_pSolver->Solve();    
     //m_pSolver->Solve_withPureNeumann();

     int iter;
     double residual;
     m_pSolver->GetNumIterations(&iter);
     m_pSolver->GetFinalRelativeResidualNorm(&residual);            
     printf("ID %d: \tnumber of iterations = %d\n", m_myid, iter);
     printf("ID %d: \trel residual norm = %e\n", m_myid, residual);
     
     double *x = new double[m_nLocalUnknowns];
     m_pSolver->Get_x(x);
     
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
 * see also tecplot_plot_surfaces().
 */
void EBM3D_CARTESIAN::saveInterface_Tecplot(const char *filename)
{    
    INTERFACE *intfc    = m_pFront->interf;
        
    FILE    *hfile;
    POINT   *p;
    SURFACE **s;
    TRI	    *tri;
   
    double  *crds;
    int	    i, j, k, l;
    int     npts, ntris, nsurfs;
    double  *pts = NULL;
    int     *verts = NULL;

    
    for (ntris = 0, s = intfc->surfaces; s && *s; ++s)
    {
        ntris += (*s)->num_tri;
        for (tri=first_tri(*s); !at_end_of_tri_list(tri,*s); tri=tri->next)
        {
            for (k = 0; k < 3; ++k)
                Index_of_point(Point_of_tri(tri)[k]) = -1;
        }
    }

    pts = new double[3*intfc->num_points];
    verts = new int[4*ntris];
    
    for (npts=0, ntris=0, nsurfs=0, s = intfc->surfaces; s && *s; ++s)
    {      
        for (tri = first_tri(*s); !at_end_of_tri_list(tri,*s); 
             tri = tri->next)
        {
            for (k = 0; k < 3; ++k)
            {
                p = Point_of_tri(tri)[k];
                if (Index_of_point(p) == -1)
                {
                    crds = Coords(p);
                    for (l = 0; l < 3; ++l)
                        pts[3*npts+l] = crds[l];
                    Index_of_point(p) = npts++;
                }
                verts[4*ntris+k] = Index_of_point(p);
            }
            verts[4*ntris+3] = nsurfs;
            ++ntris;
        }
        ++nsurfs;        
    }
    if (nsurfs == 0)
        return;

    if ((hfile = fopen(filename,"w")) == NULL)
    {
        (void) printf("EBM3D_CARTESIAN::saveInterface_Tecplot: "
                      "can't open %s\n",filename);
        return;
    }

    fprintf(hfile, "VARIABLES = \"x\", \"y\", \"z\"\n");
    
    for (nsurfs=0, s = intfc->surfaces; s && *s; ++s)
    {      
        if(nsurfs==0)   // first zone
        {
            fprintf(hfile, "ZONE N=%d, E=%d, DATAPACKING=POINT, ZONETYPE=FETRIANGLE\n",
                    npts, (*s)->num_tri);
            for (i = 0; i < npts; ++i)
                fprintf(hfile,"%-9g %-9g %-9g\n",
                           pts[3*i],pts[3*i+1],pts[3*i+2]);            
        }
        else
            fprintf(hfile, "ZONE N=%d, E=%d, DATAPACKING=POINT, ZONETYPE=FETRIANGLE VARSHARELIST=([1 2 3]=1)\n",
                    npts, (*s)->num_tri);
        
        for (j = 0; j < ntris; ++j)
            fprintf(hfile," %-4d %-4d %-4d \n", 
                       verts[4*j]+1,verts[4*j+1]+1,verts[4*j+2]+1);                 
        ++nsurfs;        
    }

    (void) fclose(hfile);
    
    delete [] pts;
    delete [] verts;
    
}

void EBM3D_CARTESIAN::saveReconstructedInterface_Tecplot(const char *filename, int i, int j, int k)
{
    int comp[8], index, size;
    double coords[8][3], crx[19][3];
    
    std::map<int, double> volumes;
    std::map<int, std::vector<EBM3D_TRIANGLE> > intfc;    
    std::map<int, std::vector<EBM3D_TRIANGLE> >::iterator iterMap;
    
    EBM3D_TRIANGLE triangle;
    
    FILE *hfile;
    if((hfile = fopen(filename,"w")) == NULL )
    {
	 (void) printf("EBM3D_CELL::saveIntfc_tecplot: "
		       "can't open %s\n",filename);
	 return;
    }
    
    int range[3][2];
    setPlotRange(range, i,j,k);
    fprintf(hfile, "VARIABLES = X Y Z U V W\n");
    
    for(int pk=range[2][0]; pk<range[2][1]; pk++)
	 for(int pj=range[1][0]; pj<range[1][1]; pj++)
	     for(int pi=range[0][0]; pi<range[0][1]; pi++)
            {
                getCellCompCoordsCrx(comp,coords,crx,pi,pj,pk);
                //getCell(pi,pj,pk)->init();
                getCell(pi,pj,pk)->setCompCoordsCrx(comp,coords,crx);
                                
                volumes.clear();
                intfc.clear();
                getCell(pi,pj,pk)->getVolumeAndIntfc(volumes, intfc);
                
                iterMap = intfc.begin();                    
                while(iterMap!=intfc.end())
                {
                    size = (*iterMap).second.size();
                    fprintf(hfile, "ZONE N=%d, E=%d, DATAPACKING=POINT, ZONETYPE=FETRIANGLE\n",
                                            3*size, size);
                    for(int m=0; m<(*iterMap).second.size(); m++)
                    {            
                        triangle = (*iterMap).second[m];
                        for(int n=0; n<3; n++)
                            fprintf(hfile,"%-9g %-9g %-9g 0 0 0\n", triangle.m_vertex[n][0], triangle.m_vertex[n][1], triangle.m_vertex[n][2]);                     
                    }        
                    for(int m=0; m<(*iterMap).second.size(); m++)
                        fprintf(hfile, "%d %d %d \n", 3*m+1, 3*m+2, 3*m+3);

                    iterMap++;
                }   
            }
        
    (void) fclose(hfile);    
}

/*
 * Save the edges of the emb_grid_intfc.            
 */
void EBM3D_CARTESIAN::saveInterfaceEdges_Tecplot(const char *filename)
{   
/*
    Table *T = table_of_interface(m_pFront->emb_grid_intfc);
    BLK_EDGE	**blk_edge = T->blk_edge;
    
    RECT_GRID *rect_grid = m_pFront->rect_grid;
    int xmax = rect_grid->gmax[0];
    int ymax = rect_grid->gmax[1];
    
    
    FILE *hfile=fopen(filename, "w");
    if(hfile==NULL)
    {
        printf("EBM3D_CARTESIAN::saveInterfaceEdges_Tecplot: can't open %s\n",
                filename);
        exit(0);
    }
    
    fprintf(hfile, "VARIABLES = X, Y, VALUES \n");
    fprintf(hfile, "ZONE I=%d, J=%d, F=POINT\n", xmax, ymax);
    
    
    for(int j=0; j<ymax; j++)
        for(int i=0; i<xmax; i++)
        {
            fprintf(hfile, "%f %f %d \n",
                    cell_center(i,0,rect_grid), cell_center(j,1,rect_grid),
                    blk_edge[i][j].ctype);
        }
    
    fclose(hfile);    
*/
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
void EBM3D_CARTESIAN::savePartialCells_Tecplot(const char *filename,
					       std::map<int, std::vector<EBM3D_PARTIALFACE> > &mapPartialFaces,
					       std::map<int, std::vector<EBM3D_TRIANGLE> > &mapIntfc,
					       int i, int j, int k)
{
     int ntri;
     
     FILE *hfile = fopen(filename, "w");
     if(hfile==NULL)
     {
	  printf("EBM3D_CARTESIAN::savePartialCells_Tecplot: can't open %s for writing!\n", filename);
	  exit(0);
     }

     fprintf(hfile, "VARIABLES = X Y Z \n");
     fprintf(hfile, "# (%d,%d,%d) \n", i,j,k);
     std::map<int, std::vector<EBM3D_PARTIALFACE> >::iterator iterPartialFaces;
     iterPartialFaces = mapPartialFaces.begin();
     while(iterPartialFaces!=mapPartialFaces.end())
     {
	  fprintf(hfile, "# face = %d \n", (*iterPartialFaces).first);
	  
	  std::vector<EBM3D_PARTIALFACE> &partialfaces = (*iterPartialFaces).second;
	  for(int m = 0; m<partialfaces.size(); m++)
	  {
	       std::vector<EBM3D_TRIANGLE> &triangles = partialfaces[m].m_triangles;
	       ntri = triangles.size();
	       if(ntri>0)
	       {
		    fprintf(hfile, "# comp = %d\n", partialfaces[m].m_comp);
		    fprintf(hfile, "ZONE N=%d, E=%d, DATAPACKING=POINT, ZONETYPE=FETRIANGLE \n", ntri*3, ntri);
		    for(int n=0; n<ntri; n++)
		    {
			 for(int nn=0; nn<3; nn++)
			      fprintf(hfile, "%e %e %e\n", 
				      triangles[n].m_vertex[nn][0], 
				      triangles[n].m_vertex[nn][1], 
				      triangles[n].m_vertex[nn][2]);
		    }
		    for(int n=0; n<ntri; n++)
			 fprintf(hfile, "%d %d %d\n", 3*n+1, 3*n+2, 3*n+3);
	       }
	       
	  }
	  iterPartialFaces++;
     }

     fprintf(hfile, "# intfc\n");
     std::map<int, std::vector<EBM3D_TRIANGLE> >::iterator iterIntfc;
     iterIntfc = mapIntfc.begin();
     while(iterIntfc!=mapIntfc.end())
     {
	  std::vector<EBM3D_TRIANGLE> &triangles = (*iterIntfc).second;
	  ntri = triangles.size();
	  fprintf(hfile, "ZONE N=%d, E=%d, DATAPACKING=POINT, ZONETYPE=FETRIANGLE \n", ntri*3, ntri);
	  for(int n=0; n<ntri; n++)
	  {
	       for(int nn=0; nn<3; nn++)
		    fprintf(hfile, "%e %e %e\n", 
			    triangles[n].m_vertex[nn][0], 
			    triangles[n].m_vertex[nn][1], 
			    triangles[n].m_vertex[nn][2]);
	  }
	  for(int n=0; n<ntri; n++)
	       fprintf(hfile, "%d %d %d\n", 3*n+1, 3*n+2, 3*n+3);
	       iterIntfc++;
     }

     fclose(hfile);
}
/*
 * save errors inside the computational domain.
 * the range of the i is the following (j or k is similar):
 *       i belongs to [0,m_gmax[0]) if i<0;
 *                    [i, i]        if i>=0;
 */
void EBM3D_CARTESIAN::saveStates_Tecplot(int i, int j, int k)
{
     char filename[100];
     sprintf(filename, "states_%d.plt", m_myid);

    RECT_GRID *rect_grid = m_pFront->rect_grid;
    INTERFACE *intfc    = m_pFront->interf;

    FILE *hfile=fopen(filename, "w");
    if(hfile==NULL)
    {
        printf("EBM3D_CARTESIAN::saveStates_Tecplot: can't open %s \n", filename);
        exit(0);
    }
    
    fprintf(hfile, "VARIABLES = X, Y, Z, Approximate, Exact, MaxError \n");
    
    
    int comp[8], max_i, max_j, max_k, range[3][2];
    int debug_min[3];
    debug_min[0] = m_gmax[0];
    debug_min[1] = m_gmax[1];
    debug_min[2] = m_gmax[2];

    setPlotRange(range, i,j,k);

    fprintf(hfile, "ZONE I=%d, J=%d, K=%d, F=POINT\n", 
	    range[0][1]-range[0][0], range[1][1]-range[1][0], range[2][1]-range[2][0]);

    double max_error = 0, approximate, exact, error, tmp;
    double cellCenter[3];
    max_i = max_j = max_k = 0;
    for(int pk=range[2][0]; pk<range[2][1]; pk++)
    for(int pj=range[1][0]; pj<range[1][1]; pj++)
    for(int pi=range[0][0]; pi<range[0][1]; pi++)
    {    	 
	getCellComp(comp,pi,pj,pk);            
        getCellCenter(cellCenter,pi,pj,pk);
	approximate = exact = error = 0;
        for(int corner=0; corner<8; corner++)
        {
	     if(comp[corner]<0)
		  continue;
	     //tmp = m_unknowns[getCell(pi,pj,pk)->getCellUnknownIndex(corner)];                            
	     tmp = getCellUnknown(comp[corner],pi,pj,pk);
	     tmp = fabs( tmp - m_pLaplace->getExactSolution(cellCenter, comp[corner]) );
	     if(tmp>error)
	     {
		  //approximate = m_unknowns[getCell(pi,pj,pk)->getCellUnknownIndex(corner)];
		  approximate = getCellUnknown(comp[corner],pi,pj,pk);
		  exact = m_pLaplace->getExactSolution(cellCenter, comp[corner]);
		  error = tmp;
	     }
	     if(error>max_error)
	     {
		  max_i = pi;
		  max_j = pj;
		  max_k = pk;
		  max_error = error;
		  
	     }
	     
        }            
        fprintf(hfile, "%f %f %f %e %e %e\n", 
                cellCenter[0], cellCenter[1], cellCenter[2], approximate, exact, error);                        
    }
            
    fclose(hfile);    
    printf("ID %d: EBM3D_CARTESIAN::saveStates_Tecplot: \n",m_myid);
    printf("ID %d: \tthe max error is %f in cell (%d,%d,%d)\n",m_myid,max_error,max_i,max_j,max_k);
    printf("ID %d: \tdebug_min: {%d,%d,%d} \n", m_myid,debug_min[0],debug_min[1],debug_min[2]);
}

void EBM3D_CARTESIAN::saveStates_VTK(int i, int j, int k)
{
     char filename[100];
     sprintf(filename, "states_%d.vtk", m_myid);

    RECT_GRID *rect_grid = m_pFront->rect_grid;
    INTERFACE *intfc    = m_pFront->interf;

    FILE *hfile=fopen(filename, "w");
    if(hfile==NULL)
    {
        printf("EBM3D_CARTESIAN::saveStates_VTK: can't open %s \n", filename);
        exit(0);
    }
    
    double cellCenter[3];    
    double origin[3];
    
    int comp[8], max_i, max_j, max_k, range[3][2];
    int debug_min[3];
    double max_error = 0, approximate, exact, error, tmp;

    debug_min[0] = m_gmax[0];
    debug_min[1] = m_gmax[1];
    debug_min[2] = m_gmax[2];

    setPlotRange(range, i,j,k);

    getCellCornerCoords(0,origin,range[0][0],range[1][0],range[2][0]);
    
    fprintf(hfile,"# vtk DataFile Version 3.0\n");
    fprintf(hfile,"Max_Error\n");
    fprintf(hfile,"ASCII\n");
    fprintf(hfile,"DATASET STRUCTURED_POINTS\n");
    fprintf(hfile,"DIMENSIONS %d %d %d\n", 
	    range[0][1]-range[0][0]+1,
	    range[1][1]-range[1][0]+1,
	    range[2][1]-range[2][0]+1);
    fprintf(hfile,"SPACING %f %f %f\n", m_dx, m_dy,m_dz);
    fprintf(hfile,"ORIGIN %f %f %f\n", origin[0], origin[1], origin[2]);
    fprintf(hfile,"CELL_DATA %d\n", 
	    (range[0][1]-range[0][0])*
	    (range[1][1]-range[1][0])*
	    (range[2][1]-range[2][0]));
	    
    fprintf(hfile,"SCALARS Max_Error double 1\n");
    fprintf(hfile,"LOOKUP_TABLE DEFAULT\n");


    max_i = max_j = max_k = 0;
    for(int pk=range[2][0]; pk<range[2][1]; pk++)
    for(int pj=range[1][0]; pj<range[1][1]; pj++)
    for(int pi=range[0][0]; pi<range[0][1]; pi++)
    {    	 
	getCellComp(comp,pi,pj,pk);            
        getCellCenter(cellCenter,pi,pj,pk);
	approximate = exact = error = 0;
        for(int corner=0; corner<8; corner++)
        {
	     if(comp[corner]<0)
		  continue;
	     //tmp = m_unknowns[getCell(pi,pj,pk)->getCellUnknownIndex(corner)];                            
	     tmp = getCellUnknown(comp[corner],pi,pj,pk);
	     tmp = fabs( tmp - m_pLaplace->getExactSolution(cellCenter, comp[corner]) );
	     if(tmp>error)
	     {
		  //approximate = m_unknowns[getCell(pi,pj,pk)->getCellUnknownIndex(corner)];
		  approximate = getCellUnknown(comp[corner],pi,pj,pk);
		  exact = m_pLaplace->getExactSolution(cellCenter, comp[corner]);
		  error = tmp;
	     }
	     if(error>max_error)
	     {
		  max_i = pi;
		  max_j = pj;
		  max_k = pk;
		  max_error = error;
		  
	     }
	     
        }            
        fprintf(hfile, "%f\n", error);                        
    }
            
    fclose(hfile);    
    printf("ID %d: EBM3D_CARTESIAN::saveStates_VTK: \n",m_myid);
    printf("ID %d: \tthe max error is %f in cell (%d,%d,%d)\n",m_myid,max_error,max_i,max_j,max_k);
    printf("ID %d: \tdebug_min: {%d,%d,%d} \n", m_myid,debug_min[0],debug_min[1],debug_min[2]);
}

void EBM3D_CARTESIAN::saveComponent_Tecplot(const char *filename)
{

    RECT_GRID *rect_grid = m_pFront->rect_grid;
    INTERFACE *intfc    = m_pFront->interf;

    FILE *hfile=fopen(filename, "w");
    if(hfile==NULL)
    {
        printf("EBM3D_CARTESIAN::saveComponent_Tecplot: can't open %s \n", filename);
        exit(0);
    }
    
    fprintf(hfile, "VARIABLES = X, Y, Z, VALUES \n");
    fprintf(hfile, "ZONE I=%d, J=%d, K=%d, F=POINT\n", m_gmax[0]+1, m_gmax[1]+1, m_gmax[2]+1);
    
    double comp, coords[3];    
    for(int k=0; k<=m_lbuf[2]+m_gmax[2]+m_ubuf[2]; k++)
    for(int j=0; j<=m_lbuf[1]+m_gmax[1]+m_ubuf[1]; j++)
    for(int i=0; i<=m_lbuf[0]+m_gmax[0]+m_ubuf[0]; i++)
    {            
        comp = getCellCornerComp(getCellCornerIndex(0,i,j,k));
	getCellCornerCoords(0, coords, i,j,k);
        
        fprintf(hfile, "%f %f %f %e\n", 
                coords[0], coords[1], coords[2], comp);                        
    }
            
    fclose(hfile);        
}

void EBM3D_CARTESIAN::debug(std::vector< std::pair<int,double> > &stencil)
{
    printf("stencil:\n");
    for(int i=0; i<stencil.size(); i++)    
        printf("\t(%d,%f)\n",stencil[i].first, stencil[i].second);    
}

void EBM3D_CARTESIAN::debug_printNeighborComp(int i, int j, int k)
{
    int comp[8];    
    for(int pk=k-1; pk<=k+1; pk++)
	 for(int pj=j-1; pj<=j+1; pj++)
	      for(int pi=i-1; pi<=i+1; pi++)
	      {
		   getCellComp(comp,pi,pj,pk); 
		   debug_printCellComp(comp,pi,pj,pk);
	      }    
}

void EBM3D_CARTESIAN::debug_printNeighborComp(EBM_PLANE plane, int i, int j, int k)
{
     int comp[8];    
     int range[3][2];
     range[0][0] = i - 1;   range[0][1] = i + 1;
     range[1][0] = j - 1;   range[1][1] = j + 1;
     range[2][0] = k - 1;   range[2][1] = k + 1;
     
     switch(plane)
     {
     case PLANE_X:
	  range[0][0] = range[0][1] = i;
	  break;
     case PLANE_Y:
	  range[1][0] = range[1][1] = j;
	  break;
     case PLANE_Z:
	  range[2][0] = range[2][1] = k;
	  break;
     }
     for(int pk=range[2][0]; pk<=range[2][1]; pk++)
	  for(int pj=range[1][0]; pj<=range[1][1]; pj++)
	       for(int pi=range[0][0]; pi<=range[0][1]; pi++)
	       {
		    getCellComp(comp,pi,pj,pk); 
		    debug_printCellComp(comp,pi,pj,pk);
	       }    

}


void EBM3D_CARTESIAN::debug_printCellComp(int comp[8], int i, int j, int k)
{
    printf("(%d,%d,%d) has comp = {%d,%d,%d,%d,%d,%d,%d,%d}\n",
            i,j,k, 
            comp[0], comp[1], comp[2], comp[3],
            comp[4], comp[5], comp[6], comp[7]);
}

void EBM3D_CARTESIAN::debug_get4PointsSolutionInterpCoeffs(void)
{
    double coords[4][3], x[3], coeffs[4];
    
    coords[0][0] = -1; coords[0][1] = -1; coords[0][2] = -1;
    coords[1][0] =  1; coords[1][1] = -1; coords[1][2] = -1;
    coords[2][0] =  0; coords[2][1] =  1; coords[2][2] = -1;
    coords[3][0] =  0; coords[3][1] =  0; coords[3][2] =  3;
    
    double volume = EBM3D_CELL::getTetraVolume(coords[0],coords[1],coords[2],coords[3]);    
    printf("volume = %e \n", volume);

    x[0] = 1.0/3; x[1] = 1.0/4; x[2] = 1.0/5;
    get4PointsSolutionInterpCoeffs(coords, x, coeffs);
    printf("x: %f, %f, %f\n", x[0], x[1], x[2]);
    printf("coeffs: %f, %f, %f, %f\n", coeffs[0], coeffs[1], coeffs[2], coeffs[3]);      
    for(int i=0; i<3; i++)
    {
        x[i] = 0;
        for(int j=0; j<4; j++)
            x[i] += coeffs[j]*coords[j][i];
    }
    printf("x': %f, %f, %f\n\n", x[0], x[1], x[2]);            
    
}


void EBM3D_CARTESIAN::debug_get6PointsSolutionInterpCoeffs_Plane(void)
{
    const int N = 6;
    int i, j;
    double sum = 0, coords[N][2], x[2], coeffs[N];
    //srand(10);
    //int r = rand();
    
/*
    for(i=0; i<6; i++)
    {        
        coords[i][0] = (i+1)/11.0;
        coords[i][1] = coords[i][0]*coords[i][0];
        //coords[i][2] = coords[i][0]*coords[i][0]*coords[i][0];
        printf("coords[%d]: %f, %f\n", i, coords[i][0], coords[i][1]);
    }
*/
    coords[0][0] = 0; coords[0][1] = 1; 
    coords[1][0] = 1; coords[1][1] = 1; 
    coords[2][0] = 2; coords[2][1] = 1; 
    coords[3][0] = 3; coords[3][1] = 1; 
    coords[4][0] = 4; coords[4][1] = 1; 
    coords[5][0] = 0; coords[5][1] = 2; 
   
    
    
    for(int i=0; i<N; i++)
    {
        coeffs[i] = i;
        sum += coeffs[i];
    }
    
    x[0] = x[1] = 0;
    for(i=0; i<N; i++)
    {
        coeffs[i] /= sum;
        for(j=0; j<2; j++)
            x[j] += coeffs[i]*coords[i][j];
    }
    
    printf("x: %f, %f \n", x[0], x[1]);        
    printf("coeffs: ");
    for(i=0; i<N; i++)
        printf(" %f,", coeffs[i]);
    printf("\n");
    
    get6PointsSolutionInterpCoeffs_Plane(coords, x, coeffs);
    
    printf("x: %f, %f \n", x[0], x[1]);        
    printf("coeffs: ");
    for(i=0; i<N; i++)
        printf(" %f,", coeffs[i]);
    printf("\n");
}

void EBM3D_CARTESIAN::debug_get10PointsSolutionInterpCoeffs(void)
{
    int i, j;
    double sum = 0, coords[10][3], x[3], coeffs[10];
    //srand(10);
    int r = rand();
    
    for(i=0; i<10; i++)
    {
        //coords[i][0] = double(std::rand())/RAND_MAX;
        //coords[i][1] = double(std::rand())/RAND_MAX;
        //coords[i][2] = double(std::rand())/RAND_MAX;
        coords[i][0] = (i+1)/11.0;
        coords[i][1] = coords[i][0]*coords[i][0];
        coords[i][2] = coords[i][0]*coords[i][0]*coords[i][0];
        printf("coords[%d]: %f, %f, %f\n", i, coords[i][0], coords[i][1], coords[i][2]);
    }
    coords[0][0] = 0; coords[0][1] = 0; coords[0][2] = 0;
    coords[1][0] = 1; coords[1][1] = 0; coords[1][2] = 0;
    coords[2][0] = 1; coords[2][1] = 1; coords[2][2] = 0;
    coords[3][0] = 0; coords[3][1] = 1; coords[3][2] = 0;
    coords[4][0] = 0; coords[4][1] = 0; coords[4][2] = 1;
    coords[5][0] = 1; coords[5][1] = 0; coords[5][2] = 1;
    coords[6][0] = 1; coords[6][1] = 1; coords[6][2] = 1;
    coords[7][0] = 0; coords[7][1] = 1; coords[7][2] = 1;
    coords[8][0] = 2; coords[8][1] = 2; coords[8][2] = 1;
    coords[9][0] = 1; coords[9][1] = 1; coords[9][2] = 2;
    
    
    for(int i=0; i<10; i++)
    {
        coeffs[i] = i;
        sum += coeffs[i];
    }
    
    x[0] = x[1] = x[2] = 0;
    for(i=0; i<10; i++)
    {
        coeffs[i] /= sum;
        for(j=0; j<3; j++)
            x[j] += coeffs[i]*coords[i][j];
    }
    
    printf("x: %f, %f, %f\n", x[0], x[1], x[2]);        
    printf("coeffs: ");
    for(i=0; i<10; i++)
        printf(" %f,", coeffs[i]);
    printf("\n");
    
    get10PointsSolutionInterpCoeffs(coords, x, coeffs);
    
    printf("x: %f, %f, %f\n", x[0], x[1], x[2]);        
    printf("coeffs: ");
    for(i=0; i<10; i++)
        printf(" %f,", coeffs[i]);
    printf("\n");
              
}


void EBM3D_CARTESIAN::debug_printNormal(FILE *hfile, double center[3], double normal[3])
{
     fprintf(hfile, "ZONE I=1 J=1 K=1 DATAPACKING=POINT\n");
     fprintf(hfile, "%f %f %f %f %f %f\n", center[0],center[1],center[2],
	     normal[0],normal[1],normal[2]);
}


void EBM3D_CARTESIAN::setPlotRange(int range[3][2], int i, int j, int k)
{
    if(i<0)
    {
	 range[0][0] = m_lbuf[0];
	 range[0][1] = m_lbuf[0]+m_gmax[0];
    }
    else
    {
	 range[0][0] = i;
	 range[0][1] = i+1;
    }
    if(j<0)
    {
	 range[1][0] = m_lbuf[1];
	 range[1][1] = m_lbuf[1]+m_gmax[1];
    }
    else
    {
	 range[1][0] = j;
	 range[1][1] = j+1;
    }
    if(k<0)
    {
	 range[2][0] = m_lbuf[2];
	 range[2][1] = m_lbuf[2]+m_gmax[2];
    }
    else
    {
	 range[2][0] = k;
	 range[2][1] = k+1;
    }
}

void EBM3D_CARTESIAN::debug_saveUnknownIndex(void)
{
     char filename[100];
     sprintf(filename, "unknownindex_%d.plt", m_myid);
    
    RECT_GRID *rect_grid = m_pFront->rect_grid;
    INTERFACE *intfc    = m_pFront->interf;

    FILE *hfile=fopen(filename, "w");
    if(hfile==NULL)
    {
        printf("EBM3D_CARTESIAN::debug_saveUnknownIndex: can't open %s \n", filename);
        exit(0);
    }
    
    fprintf(hfile, "VARIABLES = I, J, K, Index0, Index1, Index2, Index3"
                                        "Index4, Index5, Index6, Index7\n");

    fprintf(hfile, "ZONE I=%d, J=%d, K=%d, F=POINT\n", 
	    m_lbuf[0]+m_gmax[0]+m_ubuf[0], 
	    m_lbuf[1]+m_gmax[1]+m_ubuf[1],
	    m_lbuf[2]+m_gmax[2]+m_ubuf[2]);

    int index[8];
    double center[3];
    for(int k=0; k<m_lbuf[2]+m_gmax[2]+m_ubuf[2]; k++)
    for(int j=0; j<m_lbuf[1]+m_gmax[1]+m_ubuf[1]; j++)
    for(int i=0; i<m_lbuf[0]+m_gmax[0]+m_ubuf[0]; i++)
    {    
	 //getCellCenter(center,i,j);
	 getCell(i,j,k)->getCellUnknownIndex(index);
	 fprintf(hfile, "%d %d %d %d %d %d %d %d %d %d %d\n", 
		 i,j,k,
		 index[0],index[1],index[2],index[3],
		 index[4],index[5],index[6],index[7]);
    }
            
    fclose(hfile);    
}

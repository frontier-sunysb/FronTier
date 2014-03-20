/***************************************************************
FronTier is a set of libraries that implements differnt types of 
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

/*******************************************************************************
 *                  geometry.c
 ******************************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <geometry.h>
/*
 * find the intersection of two lines :
 *      line 1: p1->p2;
 *      line 2: p3->p4.
 * return 
 *      YES,    if successful;
 *      NO,   otherwise.
 * Reference:
 *   http://local.wasp.uwa.edu.au/~pbourke/geometry/lineline2d/
 */
boolean GEOMETRY::intersectLineLine2D(double p1[2], double p2[2],
				     double p3[2], double p4[2],
				     double crx[2])
{
     boolean ret = YES;
     double denominator = 
	  (p4[1]-p3[1])*(p2[0]-p1[0]) - (p4[0]-p3[0])*(p2[1]-p1[1]);
     if(fabs(denominator)<m_tol)      // two lines are almost in parallel.
     {
	  int sign = denominator>0 ? 1 : -1;
	  denominator = sign * m_tol;
     }
     
     double nominator =
	  (p4[0]-p3[0])*(p1[1]-p3[1]) - (p4[1]-p3[1])*(p1[0]-p3[0]);

     
     
     double ua = nominator/denominator;

     crx[0] = p1[0] + ua * (p2[0]-p1[0]);
     crx[1] = p1[1] + ua * (p2[1]-p1[1]);

     return ret;
}


void GEOMETRY::test(void)
{
     double p[][2] = {{0,0},
		     {1,0},
		     {1,1},
		     {0,1},
                     {0,0.5},
		     {0.5,0},
		     {1,0.5},
		     {0.5,1}};
     
     for(int i=0; i<8; i++)
	  printf("p[%d]={%f,%f}\n", i, p[i][0], p[i][1]);
     
     double crx[2];
     
     intersectLineLine2D(p[0],p[2],p[1],p[3],crx);
     printf("lines 0-2 intersects with 1-3 at {%f,%f}\n",crx[0],crx[1]);

     intersectLineLine2D(p[1],p[3],p[0],p[2],crx);
     printf("lines 1-3 intersects with 0-2 at {%f,%f}\n",crx[0],crx[1]);

     intersectLineLine2D(p[4],p[6],p[5],p[7],crx);
     printf("lines 4-6 intersects with 5-7 at {%f,%f}\n",crx[0],crx[1]);

     intersectLineLine2D(p[5],p[7],p[6],p[4],crx);
     printf("lines 5-7 intersects with 6-4 at {%f,%f}\n",crx[0],crx[1]);

     intersectLineLine2D(p[1],p[2],p[2],p[3],crx);
     printf("lines 1-2 intersects with 2-3 at {%f,%f}\n",crx[0],crx[1]);

     


}

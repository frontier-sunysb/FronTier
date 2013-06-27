/***************************************************************
FronTier is a set of libraries that implements differnt types of 
Front Traking algorithms. Front Tracking is a numerical method for 
the solution of partial differential equations whose solutions 
have discontinuities.  


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
*				imap.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*/

#include <intfc/int.h>		/* includes int.h, table.h */

EXPORT	void I_MoveNodeToPoint(
	POINT *pt,
        CURVE *curve)
{
	BOND *b;
	if (curve->start != curve->end)
	    return;
	curve_bond_loop(curve,b)
	{
	    if (pt == b->start)
	    {
		move_closed_loop_node(curve,b);
		return;
	    }
	}
}	/* end I_MoveNodeToPoint */

EXPORT	CURVE **I_SplitCurve(
	POINT *pt,
        CURVE *curve)
{
	BOND *b;	
	COMPONENT ncomp,pcomp;

	if (Dimension(curve->interface) == 2)
	{
	    ncomp = negative_component(curve);
	    pcomp = positive_component(curve);
	}
	curve_bond_loop(curve,b)
	{
	    if (pt == b->start)
		return split_curve(pt,b,curve,ncomp,pcomp,ncomp,pcomp);
	}
	return NULL;
}	/* end I_SplitCurve */

static double ave_color(TRI *tri)
{
        TRI *nbtri;
        int i,n;
        double color;

        n = 0;
        color = tri->color;
        for (i = 0; i < 3; ++i)
        {
            if (is_side_bdry(tri,i)) continue;
            n++;
            nbtri = Tri_on_side(tri,i);
            color += nbtri->color;
        }
        color /= n;
        return color;
}       /* end ave_color */

EXPORT	void I_SmoothSurfColor(
	SURFACE *surf,
	int num_rounds)
{
	double *color;
	int i,n,num_tri = surf->num_tri;
	TRI *tri;

	uni_array(&color,num_tri,FLOAT);

	for (i = 0; i < num_rounds; ++i)
	{
	    n = 0;
	    surf_tri_loop(surf,tri)
            	color[n++] = ave_color(tri);
	    n = 0;
	    surf_tri_loop(surf,tri)
		tri->color = color[n++];
	}
	free_these(1,color);
}	/* end I_SmoothSurfColor */

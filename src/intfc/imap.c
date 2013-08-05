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

LOCAL boolean curve_of_boundary_hs(CURVE*);

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

EXPORT SURFACE *I_CopySurface(
	SURFACE *surf)
{
	return copy_surface(surf,surf->pos_curves,surf->neg_curves,YES);
}	/* end I_CopySurface */

EXPORT SURFACE *I_AddTwoSurfaces(
	SURFACE *surf1,
	SURFACE *surf2)
{
	last_tri(surf1)->next = first_tri(surf2);
	first_tri(surf2)->prev = last_tri(surf1);
	link_tri_list_to_surface(first_tri(surf1),last_tri(surf2),surf1);
	delete_surface(surf2);
	return surf1;	
}	/* end I_AddTwoSurfaces */

EXPORT void I_ShiftSurface(
	SURFACE *surf,
	double *displacement)
{
	TRI *tri;
	POINT *p;
	int i,j;

	surf_tri_loop(surf,tri)
	{
	    for (i = 0; i < 3; ++i)
	    {
		p = Point_of_tri(tri)[i];
		sorted(p) = NO;
	    }
	}
	surf_tri_loop(surf,tri)
	{
	    for (i = 0; i < 3; ++i)
	    {
		p = Point_of_tri(tri)[i];
		if (sorted(p) == YES) continue;
		for (j = 0; j < 3; ++j)
		    Coords(p)[j] += displacement[j];
		sorted(p) = YES;
	    }
	}
}	/* end I_ShiftSurface */

EXPORT void I_TransInteriorIntfcPoints(
	INTERFACE *intfc,
	double *disp)
{
        POINT *p;
	int i,dim = intfc->dim;
	SURFACE **s;
	CURVE **c;
	TRI *t;
	BOND *b;

	intfc_surface_loop(intfc,s)
	{
	    surf_tri_loop(*s,t)
	    {
		for (i = 0; i < 3; ++i)
		{
		    p = Point_of_tri(t)[i];
		    sorted(p) = NO;
		}
	    }
	}
	intfc_curve_loop(intfc,c)
	{
	    b = (*c)->first;	p = b->start;
	    sorted(p) = NO;
	    curve_bond_loop(*c,b)
	    {
		p = b->end;
	    	sorted(p) = NO;
	    }
	}

	intfc_surface_loop(intfc,s)
	{
	    if (is_bdry_hs(Hyper_surf(*s)))
		continue;
	    surf_tri_loop(*s,t)
	    {
		for (i = 0; i < 3; ++i)
		{
		    p = Point_of_tri(t)[i];
		    if (sorted(p)) continue;
	    	    for (i = 0; i < dim; ++i)
			Coords(p)[i] += disp[i];
		    sorted(p) = YES;
		}
	    }
	}
	intfc_curve_loop(intfc,c)
	{
	    if (curve_of_boundary_hs(*c))
	    {
		printf("Skip boundary curve:\n");
		continue;
	    }
	    b = (*c)->first;	p = b->start;
	    if (!sorted(p))
	    {
	    	for (i = 0; i < dim; ++i)
		    Coords(p)[i] += disp[i];
	    	sorted(p) = YES;
	    }
	    curve_bond_loop(*c,b)
	    {
		p = b->end;
		if (sorted(p)) continue;
	    	for (i = 0; i < dim; ++i)
		    Coords(p)[i] += disp[i];
	    	sorted(p) = YES;
	    }
	}
}	/* end I_TransInteriorIntfcPoints */

LOCAL boolean curve_of_boundary_hs(
	CURVE *c)
{
	SURFACE **s;
	curve_pos_surf_loop(c,s)
	{
	    if (Boundary_hs(Hyper_surf(*s)))
		return YES;
	}
	curve_neg_surf_loop(c,s)
	{
	    if (Boundary_hs(Hyper_surf(*s)))
		return YES;
	}
	return NO;
}	/* end curve_of_boundary_hs */

EXPORT void I_SphericalRotatePoint(
	POINT *p,
        double *center,
        double phi,
        double theta,
        boolean first)
{
	rotate_point_with_spherical_angle(p,center,phi,theta,first);
}	/* end I_RoratePoint */

EXPORT void I_PolarRotatePoint(
	POINT *p,
        double *center,
        double phi,
        boolean first)
{
	rotate_point_with_polar_angle(p,center,phi,first);
}	/* end I_RoratePoint */

EXPORT void I_SphericalRotateInteriorIntfcPoints(
	INTERFACE *intfc,
        double *center,
        double phi,
        double theta)
{
        POINT *p;
	int i,dim = intfc->dim;
	SURFACE **s;
	CURVE **c;
	TRI *t;
	BOND *b;
	boolean first = YES;

	intfc_surface_loop(intfc,s)
	{
	    surf_tri_loop(*s,t)
	    {
		for (i = 0; i < 3; ++i)
		{
		    p = Point_of_tri(t)[i];
		    sorted(p) = NO;
		}
	    }
	}
	intfc_curve_loop(intfc,c)
	{
	    b = (*c)->first;	p = b->start;
	    sorted(p) = NO;
	    curve_bond_loop(*c,b)
	    {
		p = b->end;
	    	sorted(p) = NO;
	    }
	}

	intfc_surface_loop(intfc,s)
	{
	    if (is_bdry_hs(Hyper_surf(*s)))
		continue;
	    surf_tri_loop(*s,t)
	    {
		for (i = 0; i < 3; ++i)
		{
		    p = Point_of_tri(t)[i];
		    if (sorted(p)) continue;
		    rotate_point_with_spherical_angle(p,center,phi,theta,first);
		    if (first == YES) first = NO;
		    sorted(p) = YES;
		}
	    }
	}
	intfc_curve_loop(intfc,c)
	{
	    if (curve_of_boundary_hs(*c))
		continue;
	    b = (*c)->first;	p = b->start;
	    if (!sorted(p))
	    {
		rotate_point_with_spherical_angle(p,center,phi,theta,first);
		if (first == YES) first = NO;
	    	sorted(p) = YES;
	    }
	    curve_bond_loop(*c,b)
	    {
		p = b->end;
		if (sorted(p)) continue;
		rotate_point_with_spherical_angle(p,center,phi,theta,first);
		if (first == YES) first = NO;
	    	sorted(p) = YES;
	    }
	}
}	/* end I_SphericalRotateInteriorIntfcPoints */

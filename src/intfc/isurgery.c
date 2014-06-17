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


/*
*				trisurgery.c:
*
*       Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*/


#include <intfc/int.h>

enum {SURF_STITCH,CURV_STITCH};

enum {
	MAX_TRIS_AT_VERTEX =  100
};

struct _STITCH {
  	POINT 		*start_pt;
  	POINT		*end_pt;
	TRI             *detached_tri;
	BOND            *detached_bond;
  	int		detached_side;
};
typedef struct _STITCH STITCH;

struct _SEAM {
  	STITCH 	*stitch;
  	SURFACE *surface;
  	int	number_of_stitches;
  	POINT 	*begin_pt;
  	POINT	*finish_pt;
};
typedef struct _SEAM SEAM;

struct _SURGERY {
  	SEAM 		*seam1;
  	SEAM		*seam2;
  	TRI		**new_tris;
};
typedef struct _SURGERY SURGERY;


LOCAL   void 	align_seam_with_point(SEAM*,POINT*,SEAM*);
LOCAL 	void 	sew_surface_to_curve(SURGERY*);
LOCAL	int	count_null_sides_of_surf(SURFACE*);
LOCAL	int	set_surface_seam(SURFACE*,SEAM*);
LOCAL	void 	set_curve_seam(CURVE*,SEAM*,ORIENTATION);
LOCAL	void 	print_seam(SEAM*);
LOCAL	TRI 	*bdry_adjacent_tri(TRI*,POINT*,int);

EXPORT void identify_detached_surface_curve_pair(
	INTERFACE 	*intfc)
{
  	CURVE 		*detached_curv;
	SURFACE 	**surf,*detached_surf=NULL;
	int		num_sides = 0;
	ORIENTATION	orientation;
	SEAM		seam;
	SEAM		surface_seam;
	POINT		*curv_node_pt;
	SEAM		curve_seam;
	int		num_bonds;
	SURGERY		surgery;
	DEBUG_ENTER(identify_detached_surface_curve_pair)

	/* Assumptions:							*/
	/* The interface contains exactly one surface with tris having 	*/
	/* NULL sides.  This surface is attached to exactly one curve.	*/

	for (surf = intfc->surfaces; surf && *surf; ++surf)
	{  
	    if ((num_sides = count_null_sides_of_surf(*surf)) == 0)
	        continue;

	    detached_surf = *surf;
	    if (detached_surf->pos_curves != NULL)
            {
		detached_curv = detached_surf->pos_curves[0];
		orientation = POSITIVE_ORIENTATION;
	    }
	    else if (detached_surf->neg_curves != NULL)
	    {
		detached_curv = detached_surf->neg_curves[0];
		orientation = NEGATIVE_ORIENTATION;
	    }

	    if (debugging("surgery"))
	    {
	        (void) printf("detached surf = %p  detached curve = %p  ",
			      (void*)detached_surf,(void*)detached_curv);
		print_orientation("orientation = ",orientation,"\n");
	    }
	    break;
	}

	if (detached_surf == NULL)
	{
	    DEBUG_LEAVE(identify_detached_surface_curve_pair)
	    return;
	}

	uni_array(&surface_seam.stitch,num_sides+100,sizeof(STITCH));
	if (set_surface_seam(detached_surf,&surface_seam) < num_sides) 
        {
	    (void) printf("ERROR in identify_detached_surface_curve_pair(),\n"
			  "the seam does not contain all null sides of "
			  "surface.\n");
	    clean_up(ERROR);
	}

	curv_node_pt = detached_curv->start->posn; 
	uni_array(&seam.stitch,surface_seam.number_of_stitches,sizeof(STITCH));
	align_seam_with_point(&surface_seam,curv_node_pt,&seam);
	free(surface_seam.stitch);

	num_bonds = detached_curv->num_points - 1;
	uni_array(&curve_seam.stitch,num_bonds,sizeof(STITCH));
	set_curve_seam(detached_curv,&curve_seam,orientation);

	if (debugging("surgery"))
	{
	    (void) printf("The detached surface SEAM, after alignment.\n");
	    print_seam(&seam);
	    (void) printf("The detached curve SEAM.\n");
	    print_seam(&curve_seam);
	}

	surgery.seam1 = &seam;
	surgery.seam2 = &curve_seam;

	sew_surface_to_curve(&surgery);

	free(seam.stitch);
	free(curve_seam.stitch);
	DEBUG_LEAVE(identify_detached_surface_curve_pair)
}		/*end identify_detached_surface_curve_pair*/

LOCAL 	void	sew_surface_to_curve(
	SURGERY *surgery)
{
  	SEAM 		*surf_seam = surgery->seam1;
	SEAM 		*curv_seam = surgery->seam2;
	POINT		*surf_pt,*curv_pt,*other_pt;
	TRI		*new_tri,*prev_tri = NULL,*detached_tri;
	SURFACE		*detached_surf = surf_seam->surface;
	BOND		*detached_bond;
	BOND_TRI	**btris;
	int		s_index = 0, c_index = 0,i;
	int		work_on;
	double		surf_tri_measure,curv_tri_measure;
	int		s_stitches = surf_seam->number_of_stitches;
	int		c_stitches = curv_seam->number_of_stitches;
	int		num_tris = s_stitches + c_stitches;
	int		detached_side;
	DEBUG_ENTER(sew_surface_to_curve)

	uni_array(&surgery->new_tris,num_tris,sizeof(TRI*));

	for (i = 0; i < num_tris; ++i)
	{
	    surf_pt = surf_seam->stitch[s_index % s_stitches].start_pt;
	    curv_pt = curv_seam->stitch[c_index % c_stitches].start_pt;

	    surf_tri_measure = ((s_index < s_stitches) ?
	        distance_between_positions(Coords(curv_pt),
			Coords(surf_seam->stitch[s_index].end_pt),3) :
		HUGE_VAL);

	    curv_tri_measure = ((c_index < c_stitches) ?
		distance_between_positions(Coords(surf_pt),
			Coords(curv_seam->stitch[c_index].end_pt),3) :
	        HUGE_VAL);

	    if (surf_tri_measure < curv_tri_measure)
	        work_on = SURF_STITCH;
	    else
	        work_on = CURV_STITCH;

	    switch (work_on)
	    {
	    case SURF_STITCH:
	        other_pt = surf_seam->stitch[s_index].end_pt;
		detached_tri = surf_seam->stitch[s_index].detached_tri;
		detached_side = surf_seam->stitch[s_index++].detached_side;

		new_tri = make_tri(other_pt,surf_pt,curv_pt,
				   detached_tri,prev_tri,NULL,0);
		insert_tri_at_tail_of_list(new_tri,detached_surf);
		Tri_on_side(detached_tri,detached_side) = new_tri;
		surgery->new_tris[i] = new_tri;
		break;

	    case CURV_STITCH:
	        other_pt = curv_seam->stitch[c_index].end_pt;
		detached_bond = curv_seam->stitch[c_index++].detached_bond;

		new_tri = make_tri(surf_pt,curv_pt,other_pt,
				   prev_tri,detached_bond,NULL,2);
		insert_tri_at_tail_of_list(new_tri,detached_surf);
		surgery->new_tris[i] = new_tri;
		btris = Btris(detached_bond);
		while (*btris && ((*btris)->surface != detached_surf)) 
		    ++btris;

		(*btris)->tri = new_tri;
		break;
	    }

	    if (prev_tri != NULL)
	        Tri_on_side20(prev_tri) = new_tri;
	    prev_tri = new_tri;
	}	   

	Tri_on_side20(prev_tri) = surgery->new_tris[0];

	if (Boundary_tri(surgery->new_tris[0]))/*first new_tri was CURV_STITCH*/ 
	    Tri_on_side01(surgery->new_tris[0]) = prev_tri;

	else				       /*first new_tri was SURF_STITCH*/
	    Tri_on_side12(surgery->new_tris[0]) = prev_tri;

	free(surgery->new_tris);
	DEBUG_LEAVE(sew_surface_to_curve)
}		/*end sew_surface_to_curve*/

LOCAL   void 	align_seam_with_point(
	SEAM 		*unaligned,
	POINT 		*point,
	SEAM 		*aligned)
{
  	int 	num_stitches,i,index_of_nearest_pt;
	double  	least_dist = 10000, dist;
	POINT	*stitch_pt;

  	num_stitches = unaligned->number_of_stitches;
	for (i = 0; i < num_stitches; ++i)
	{
	    stitch_pt = unaligned->stitch[i].start_pt;
	    if ((dist = distance_between_positions(Coords(point),
	    	 Coords(stitch_pt),3)) < least_dist)
	    {
		least_dist = dist;
		index_of_nearest_pt = i;
	    }
	}

	for (i = 0; i < num_stitches; ++i)
	    ft_assign(&(aligned->stitch[i]),
		   &(unaligned->stitch[(i+index_of_nearest_pt)%num_stitches]),
		   sizeof(STITCH));
	    
	aligned->begin_pt = aligned->stitch[0].start_pt;
	aligned->finish_pt = aligned->stitch[num_stitches - 1].end_pt;
	aligned->number_of_stitches = num_stitches;
	aligned->surface = unaligned->surface;
}		/*end align_seam_with_point*/

EXPORT int next_null_sided_tri(
	TRI 		*start_tri,
	POINT 		*vertex,
	TRI 		**next_tri)
{
	int 		next_side;
	TRI		*new_tri;
	int		count = 0;

	new_tri = start_tri;
	while (new_tri != NULL)
	{
	    next_side = Next_side_at_vertex(new_tri,vertex);
	    if (next_side == ERROR)
	    {
	    	screen("In next_null_sided_tri(), next_side is ERROR\n");
		printf("new_tri vertices: %p %p %p\n",
				(void*)Point_of_tri(new_tri)[0],
				(void*)Point_of_tri(new_tri)[1],
				(void*)Point_of_tri(new_tri)[2]);
		printf("vertex = %p\n",(void*)vertex);
		clean_up(ERROR);
	    }
	    if (Neighbor_on_side(new_tri,next_side) == NULL)
	    {
	    	*next_tri = new_tri;
		return next_side;
	    }
	    if (is_side_bdry(new_tri,next_side))
	    	new_tri = bdry_adjacent_tri(new_tri,vertex,next_side);
	    else
	    	new_tri = Tri_on_side(new_tri,next_side);
	    if (new_tri == start_tri)
	    	break;
	    if (count++ > MAX_TRIS_AT_VERTEX)
	    {
	    	screen("ERROR: dead loop in next_null_sided_tri()\n");
		clean_up(ERROR);
	    }
	}
	*next_tri = NULL;
	return ERROR;
}		/*end next_null_sided_tri*/

EXPORT int prev_null_sided_tri(
	TRI 		*start_tri,
	POINT 		*vertex,
	TRI 		**prev_tri)
{
	int 		prev_side;
	TRI		*new_tri;
	int		count = 0;

	new_tri = start_tri;
	while (new_tri != NULL)
	{
	    prev_side = Prev_side_at_vertex(new_tri,vertex);
	    if (prev_side == ERROR)
	    {
	    	screen("In prev_null_sided_tri(), prev_side is ERROR\n");
		printf("new_tri vertices: %p %p %p\n",
				(void*)Point_of_tri(new_tri)[0],
				(void*)Point_of_tri(new_tri)[1],
				(void*)Point_of_tri(new_tri)[2]);
		printf("vertex = %p\n",(void*)vertex);
		clean_up(ERROR);
	    }
	    if (Neighbor_on_side(new_tri,prev_side) == NULL)
	    {
	    	*prev_tri = new_tri;
		return prev_side;
	    }
	    if (is_side_bdry(new_tri,prev_side))
	    	new_tri = bdry_adjacent_tri(new_tri,vertex,prev_side);
	    else
	    	new_tri = Tri_on_side(new_tri,prev_side);
	    if (new_tri == start_tri)
	    	break;
	    if (count++ > MAX_TRIS_AT_VERTEX)
	    {
	    	screen("ERROR: dead loop in prev_null_sided_tri()\n");
		clean_up(ERROR);
	    }
	}
	*prev_tri = NULL;
	return ERROR;
}		/*end prev_null_sided_tri*/

LOCAL	TRI *bdry_adjacent_tri(
	TRI *tri,
	POINT *p,
	int side)
{
	BOND *b,*adj_b;
	BOND_TRI **btris;
	CURVE *c;

	if (!is_side_bdry(tri,side)) return NULL;

	b = Bond_on_side(tri,side);

	c = (*Btris(b))->curve;

	if (is_closed_curve(c))
	{
	    if (b == c->first && p == b->start)
	    	adj_b = c->last;
	    else if (b == c->last && p == b->end)
	    	adj_b = c->first;
	    else
	    	adj_b = (p == b->start) ? b->prev : b->next;
	}
	else
	    adj_b = (p == b->start) ? b->prev : b->next;

	if (adj_b == NULL) return NULL;

	for (btris = Btris(adj_b); btris && *btris; ++btris)
	{
	    if ((*btris)->surface == tri->surf)
	    	return (*btris)->tri;
	}
}	/* end bdry_adjacent_tri */

LOCAL	int	count_null_sides_of_surf(
	SURFACE *surf)
{
  	int 	i, num_null_sides = 0;
	TRI	*t;

	for (t=first_tri(surf); !at_end_of_tri_list(t,surf); t=t->next)
	    for (i = 0; i < 3; ++i)
		if (Tri_on_side(t,i)==NULL) ++num_null_sides;
	return num_null_sides;
}		/*end count_null_sides_of_surf*/

LOCAL	int	set_surface_seam(
	SURFACE	*surf,
	SEAM	*seam)
{
	TRI	*tri, *next_tri;
	int	side, next_side, i=0;
	STITCH	*stitch;

        for (tri=first_tri(surf); !at_end_of_tri_list(tri,surf); tri=tri->next)
	{
	    for (side = 0; side < 3; ++side)
		if (Tri_on_side(tri,side) == NULL)
		    break;
	    if (side < 3)
	        break;
	}

	seam->number_of_stitches = 1;
	seam->surface = surf;
	stitch = seam->stitch;

	stitch[0].start_pt = Point_of_tri(tri)[side];
	stitch[0].end_pt = Point_of_tri(tri)[Next_m3(side)];
	stitch[0].detached_tri = tri;
	stitch[0].detached_side = side;

	seam->begin_pt = seam->finish_pt = stitch[0].start_pt;

	while (stitch[i].end_pt != seam->finish_pt)
	{
	    POINT *vertex;
	    ++i;
	    vertex = stitch[i].start_pt = stitch[i-1].end_pt;
	    next_side = next_null_sided_tri(tri,vertex,&next_tri);
	    stitch[i].end_pt = Point_of_tri(next_tri)[Next_m3(next_side)];
	    stitch[i].detached_tri = next_tri;
	    stitch[i].detached_side = next_side;
	    ++seam->number_of_stitches;
	    tri = next_tri;
	}
	return seam->number_of_stitches;
}		/*end set_surface_seam*/

LOCAL	void 	set_curve_seam(
	CURVE 	*curve,
	SEAM  	*seam,
	ORIENTATION orientation)
{
	int 	i,num_bonds;
	BOND 	*bond;
	STITCH	*stitch;

	num_bonds = curve->num_points - 1;
	seam->number_of_stitches = num_bonds;
	seam->surface = NULL;
	stitch = seam->stitch;

	switch (orientation)
	{
	case POSITIVE_ORIENTATION:
	    seam->begin_pt = curve->start->posn;
	    seam->finish_pt = curve->end->posn;

	    for (i = 0, bond = curve->first; i < num_bonds;
		 ++i, bond = bond->next)
	    {
	        stitch[i].start_pt = bond->start;
		stitch[i].end_pt = bond->end;
		stitch[i].detached_bond = bond;
		stitch[i].detached_side = -1; /*BOND*/
	    }
	    break;

	case NEGATIVE_ORIENTATION:
	    seam->begin_pt = curve->end->posn;
	    seam->finish_pt = curve->start->posn;

	    for (i = 0, bond = curve->last; i < num_bonds;
		 ++i, bond = bond->prev)
	    {
	        stitch[i].start_pt = bond->end;
		stitch[i].end_pt = bond->start;
		stitch[i].detached_bond = bond;
		stitch[i].detached_side = -1; /*BOND*/
	    }
	    break;

	case ORIENTATION_NOT_SET:
	default:
	    (void) printf("ERROR in set_curve_seam(), ORIENTATION_NOT_SET\n");
	    clean_up(ERROR);
	}
}		/*end set_curve_seam*/

LOCAL	void print_seam(
	SEAM 	*seam)
{
  	int i;
	double *start,*end;

	(void) printf("Start of SEAM %p\n",(void*)seam);
	(void) printf("  surface = %p    "
		      "number_of_stitches = %d\n",
		      (void*)seam->surface,seam->number_of_stitches);
	(void) printf("  begin_pt = ( %g %g %g )  "
		      "finish_pt = ( %g %g %g )\n",
		      Coords(seam->begin_pt)[0],Coords(seam->begin_pt)[1],
		      Coords(seam->begin_pt)[2],
		      Coords(seam->finish_pt)[0],Coords(seam->finish_pt)[1],
		      Coords(seam->finish_pt)[2]);
	(void) printf("\n");
	(void) printf("%3s  %30s  %30s  %8s\n",
		      "num","   start of stitch   ","    end of stitch    ",
		      "    nhbr  ");
	(void) printf("=");
	for (i = 0; i < 100; ++i) printf("=");
	printf("\n");
	for (i = 0; i < seam->number_of_stitches; ++i)
	{
	    start = Coords(seam->stitch[i].start_pt);
	    end = Coords(seam->stitch[i].end_pt);
	    (void) printf("%3d  ( %g %g %g )  "
			  "( %g %g %g )  ",
			  i,start[0],start[1],start[2],end[0],end[1],end[2]);
	    if (seam->stitch[i].detached_side == -1)
	        (void) printf(" BOND\n");
	    else if (seam->stitch[i].detached_side > 2)
	        (void) printf(" unknown\n");
	    else
		(void) printf(" TRI SIDE%d%d\n",
	                      seam->stitch[i].detached_side, 
	                      Next_m3(seam->stitch[i].detached_side)); 
	}
	(void) printf("End of SEAM %p.\n",(void*)seam);
	(void) printf("\n");
}		/*end print_seam*/

EXPORT POINT *insert_point_in_surface(
	int idir,
	double *coords,
	SURFACE *surf)
{
	INTERFACE *intfc = surf->interface;
	RECT_GRID *gr = &topological_grid(intfc);
	TRI *tri;
	double crds_start[MAXD],crds_crx[MAXD];
	double h[MAXD],crds_min,crds_max,dh;
	int i,iv,ie;
	POINT *pt = NULL;

	if (debugging("add_gore"))
	{
	    (void) printf("Entering insert_point_in_surface()\n");
	    gview_plot_surf_within_range("gvertex-0",surf,coords,3*h[0]);
	}

	for (i = 0; i < 3; ++i)
	{
	    crds_start[i] = crds_crx[i] = coords[i];
	    h[i] = gr->h[i];
	}
	dh = h[idir];

	for (tri = first_tri(surf); !at_end_of_tri_list(tri,surf);
				tri = tri->next)
	{
	    crds_min = HUGE;	crds_max = -HUGE;
	    for (i = 0; i < 3; ++i)
	    {
		if (crds_min > Coords(Point_of_tri(tri)[i])[idir])
		    crds_min = Coords(Point_of_tri(tri)[i])[idir];
		if (crds_max < Coords(Point_of_tri(tri)[i])[idir])
		    crds_max = Coords(Point_of_tri(tri)[i])[idir];
	    }
	    h[idir] = 2*dh + crds_max - crds_min;
	    crds_start[idir] = crds_crx[idir] = crds_min - dh;
	    if (tri_edge_crossing(tri,crds_start,crds_crx,idir,&iv,&ie,h))
	    {
		if (iv != ERROR)
		{
		    if (debugging("add_gore"))
			(void) printf("Verter case\n");
		    pt = Point_of_tri(tri)[iv];
		    break;
		}
		else if (ie != ERROR)
		{
		    if (debugging("add_gore"))
			(void) printf("Edge case\n");
		    pt = Point(crds_crx);
		    insert_point_in_tri_side(pt,ie,tri,surf);
		    break;
		}
		else
		{
		    if (debugging("add_gore"))
			(void) printf("Tri case\n");
		    pt = Point(crds_crx);
		    insert_point_in_tri(pt,tri,surf);
		    break;
		}
	    }
	}
	if (debugging("add_gore"))
	{
	    (void) printf("Entering insert_point_in_surface()\n");
	    gview_plot_surf_within_range("gvertex-1",surf,coords,3*h[0]);
	}
	return pt;
}	/* end insert_point_in_surface */


/*	This function insert a curve in a surface, the inputs are two
*	nodes, the start and end nodes, and the surface. This function
*	is rather expensive and is only used in initialization.
*/

EXPORT CURVE *insert_curve_in_surface(
	double *nor,
	NODE *ns,
	NODE *ne,
	SURFACE *surf)
{
	CURVE *curve;
	INTERFACE *intfc = surf->interface;
	TRI *tri,**tris,**pos_tris,**neg_tris;
	double *ps = Coords(ns->posn);
	double *pe = Coords(ne->posn);
	double v1[MAXD],v2[MAXD],d1,d2;
	double *p1,*p2;
	int i,j,k,iv,num_pts,num_tris;
	double plane[4],pc[MAXD];
	POINT **cpts,*pnew;
	boolean within_nodes;
	RECT_GRID *gr = &topological_grid(intfc);
	double *h = gr->h;
	double xdiff;
	int ix;
	BOND_TRI *btri;
	BOND *b;

	if (debugging("insert_curve_in_surf"))
	{
	    (void) printf("Entering insert_curve_in_surface()\n");
	    (void) printf("nor = %f %f %f\n",nor[0],nor[1],nor[2]);
	    (void) printf("Start and end nodes:\n");
	    (void) print_node(ns);
	    (void) print_node(ne);
	    gview_plot_interface("ginsert_curve-0",intfc);
	}
	curve = make_curve(0,0,ns,ne);
	install_curve_in_surface_bdry(surf,curve,POSITIVE_ORIENTATION);

	num_pts = 2;
	plane[3] = 0.0;
	xdiff = 0.0;
	for (i = 0; i < 3; ++i)
	{
	    plane[i] = nor[i];
	    plane[3] += nor[i]*ps[i];
	    if (xdiff < fabs(ps[i] - pe[i]))
	    {
		xdiff = fabs(ps[i] - pe[i]);
		ix = i;
	    }
	}
	reset_surface_points(surf);
	for (tri = first_tri(surf); !at_end_of_tri_list(tri,surf); 
				tri = tri->next)
	{
	    for (i = 0; i < 3; ++i)
	    {
		if (plane_side_intersection(plane,tri,i,pc,&iv))
		{
		    if (within_interval(ps[ix],pe[ix],pc[ix]))
			within_nodes = YES;
		    else
			within_nodes = NO;
		    if (!within_nodes) continue;
		    if (iv != ERROR)
		    {
			pnew = Point_of_tri(tri)[iv];
			if (sorted(pnew)) continue;
			else
			{
			    sorted(pnew) = YES;
			    num_pts++;
			    continue;
			}
		    }
		    num_pts++;
		}
	    }
	}
	FT_VectorMemoryAlloc((POINTER*)&cpts,num_pts,sizeof(POINT*));
	num_pts = 0;
	cpts[num_pts++] = ns->posn;
	for (tri = first_tri(surf); !at_end_of_tri_list(tri,surf); 
				tri = tri->next)
	{
repeat:
	    for (i = 0; i < 3; ++i)
	    {
		if (plane_side_intersection(plane,tri,i,pc,&iv))
		{
		    if (within_interval(ps[ix],pe[ix],pc[ix]))
			within_nodes = YES;
		    else
			within_nodes = NO;
		    if (!within_nodes) continue;
		    if (iv != ERROR)
		    {
			pnew = Point_of_tri(tri)[iv];
			if (pnew == ns->posn || pnew == ne->posn) 
			    continue;
			if (pointer_in_list((POINTER)pnew,num_pts,
				(POINTER*)cpts))
			    continue;
			cpts[num_pts++] = pnew;
		    }
		    else
		    {
			pnew = Point(pc);
			insert_point_in_tri_side(pnew,i,tri,surf);
			cpts[num_pts++] = pnew;
			goto repeat;
		    }
		}
	    }
	}
	cpts[num_pts++] = ne->posn;
	if (debugging("insert_curve_in_surf"))
	    (void) printf("Total number of points in curve: %d\n",num_pts);

	/* Make sure all new points are attached to a tri */
	for (i = 0; i < num_pts; ++i)
	{
	    cpts[i]->hse = NULL;
	    for (tri = first_tri(surf); !at_end_of_tri_list(tri,surf); 
				tri = tri->next)
	    {
	        for (iv = 0; iv < 3; ++iv)
		{
		    pnew = Point_of_tri(tri)[iv];
		    if (pnew == cpts[i])
		    {
			cpts[i]->hse = Hyper_surf_element(tri);	
			break;
		    }
		}
		if (iv < 3) break;
	    }
	    if (cpts[i]->hse == NULL)
	    {
		(void) printf("cpts[%d] cannot find attached tri!\n",i);
		clean_up(ERROR);
	    }
	}
	/* Sort new points along ns-->ne direction */
	FT_VectorMemoryAlloc((POINTER*)&pos_tris,num_pts-1,sizeof(TRI*));
	FT_VectorMemoryAlloc((POINTER*)&neg_tris,num_pts-1,sizeof(TRI*));
	for (i = 0; i < num_pts-1; ++i)
	{
	    boolean next_found = NO;

	    tri = Tri_of_hse(cpts[i]->hse);
	    num_tris = FT_FirstRingTrisAroundPoint(cpts[i],tri,&tris);
	    for (j = 0; j < num_tris; ++j)
	    {
	    	for (iv = 0; iv < 3; ++iv)
		{
		    pnew = Point_of_tri(tris[j])[iv];
		    if (pnew == cpts[i]) continue;	/* skip self */
		    for (k = i+1; k < num_pts; ++k)
		    {
			if (pnew == cpts[k])
			{
			    if (k != i+1)
			    {
				pnew = cpts[k];
				cpts[k] = cpts[i+1];
				cpts[i+1] = pnew;
			    }
			    next_found = YES;
			    break;
			}
		    }
		    if (next_found) break;
		}
		if (next_found) break;
	    }
	    if (!next_found)
	    {
		(void) printf("Consecutive ordering failed!\n");
		clean_up(ERROR);
	    }
	    pos_tris[i] = neg_tris[i] = NULL;
	    for (j = 0; j < num_tris; ++j)
	    {
		for (k = 0; k < 3; ++k)
		{
		    if (cpts[i] == Point_of_tri(tris[j])[k] &&
			cpts[i+1] == Point_of_tri(tris[j])[(k+1)%3])
			pos_tris[i] = tris[j];
		    if (cpts[i+1] == Point_of_tri(tris[j])[k] &&
			cpts[i] == Point_of_tri(tris[j])[(k+1)%3])
			neg_tris[i] = tris[j];
		}
	    }
	    if (pos_tris[i] == NULL)
	    {
		printf("pos_tris[%d] not found\n",i);
	    }
	    if (neg_tris[i] == NULL)
	    {
		printf("neg_tris[%d] not found\n",i);
	    }
	}
	if (debugging("insert_curve_in_surf"))
	{
	    (void) printf("Sorted points:\n");
	    (void) printf("ps = %f %f\n",ps[0],ps[1]);
	    for (i = 1; i < num_pts-1; ++i)
	    	(void) printf("pc[%d] = %f %f\n",i,Coords(cpts[i])[0],
				Coords(cpts[i])[1]);
	    (void) printf("pe = %f %f\n",pe[0],pe[1]);
	}
	for (i = 1; i < num_pts-1; ++i)
	{
	    insert_point_in_bond(cpts[i],curve->last,curve);
	}
	for (i = 0, b = curve->first; b!= NULL; b = b->next, i++)
	{
	    btri = link_tri_to_bond(NULL,pos_tris[i],surf,b,curve);
	    btri = link_tri_to_bond(NULL,neg_tris[i],surf,b,curve);
	}

	FT_FreeThese(3,cpts,pos_tris,neg_tris);
	if (debugging("insert_curve_in_surf"))
	{
	    (void) printf("Leaving insert_curve_in_surf()\n");
	    gview_plot_interface("ginsert_curve-1",intfc);
	}
	reset_intfc_num_points(surf->interface);
	return curve;
}	/* end insert_curve_in_surf */

/*	
 *	This is a two dimensional rotation with polar angles 
 */

EXPORT void rotate_point_with_polar_angle(
	POINT *p,
	double *center,			/* Rotation center */
	double phi,			/* Angle in counter-clock wise dir */
	boolean first)
{
	static double rot_matrix[2][2];
	double x0,y0,x1,y1;

	if (first)
	{
	    rot_matrix[0][0] = cos(phi);	 
	    rot_matrix[0][1] = -sin(phi);	 
	    rot_matrix[1][0] = sin(phi);	 
	    rot_matrix[1][1] = cos(phi);	 
	}
	x0 = Coords(p)[0] - center[0];
	y0 = Coords(p)[1] - center[1];
	x1 = rot_matrix[0][0]*x0 + rot_matrix[0][1]*y0;
	y1 = rot_matrix[1][0]*x0 + rot_matrix[1][1]*y0;
	Coords(p)[0] = x1 + center[0];
	Coords(p)[1] = y1 + center[1];
}	/* end rotate_point_with_polar_angle */

/*	
 *	This is a three dimensional rotation with polar angles 
 */

EXPORT void rotate_point_with_spherical_angles(
	POINT *p,
	double *center,			/* Rotation center */
	double phi,			/* Azimuthal angle */
	double theta,			/* Polar angle */
	boolean first)
{
	static double cc,ss,cs,sc,c,s;
	double pt[3],rho,r;
	int i;

	if (first)
	{
	    cc = cos(theta)*cos(phi);
	    ss = sin(theta)*sin(phi);
	    cs = cos(theta)*sin(phi);
	    sc = sin(theta)*cos(phi);
	    c = cos(theta);
	    s = sin(theta);
	}
	r = 0;
	for (i = 0; i < 3; i++)
	{
	    pt[i] = Coords(p)[i] - center[i];
	    r += sqr(pt[i]);
	}
	if (r == 0.0) return;
	rho = sqrt(sqr(pt[0]) + sqr(pt[1]));
	Coords(p)[0] = cc*pt[0] - cs*pt[1] + 
		sc*pt[2]*pt[0]/rho - ss*pt[2]*pt[1]/rho;
	Coords(p)[1] = cc*pt[1] - cs*pt[0] + 
		sc*pt[2]*pt[1]/rho - ss*pt[2]*pt[0]/rho;
	Coords(p)[2] = c*pt[2] - rho*s;
	for (i = 0; i < 3; i++)
	    Coords(p)[i] += center[i];
}	/* end rotate_point_with_spherical_angle */

EXPORT void rotate_point_with_spherical_angle(
	POINT *p,
	double *center,			/* Rotation center */
	double phi,			/* Azimuthal angle */
	double theta,			/* Polar angle */
	boolean first)
{
	static double roty_matrix[3][3];
	static double rotz_matrix[3][3];
	double pt[3];
	int i,j;

	if (first == YES)
	{
	    /* Set up rotation matrix */
	    roty_matrix[0][0] = cos(theta);
	    roty_matrix[0][1] = 0.0;
	    roty_matrix[0][2] = sin(theta);
	    roty_matrix[1][0] = 0.0;
	    roty_matrix[1][1] = 1.0;
	    roty_matrix[1][2] = 0.0;
	    roty_matrix[2][0] = -sin(theta);
	    roty_matrix[2][1] = 0.0;
	    roty_matrix[2][2] = cos(theta);
	    
	    rotz_matrix[0][0] = cos(phi);
	    rotz_matrix[0][1] = -sin(phi);
	    rotz_matrix[0][2] = 0.0;
	    rotz_matrix[1][0] = sin(phi);
	    rotz_matrix[1][1] = cos(phi);
	    rotz_matrix[1][2] = 0.0;
	    rotz_matrix[2][0] = 0.0;
	    rotz_matrix[2][1] = 0.0;
	    rotz_matrix[2][2] = 1.0;
	}
	for (i = 0; i < 3; i++)
	    Coords(p)[i] -= center[i];
	for (i = 0; i < 3; i++)
	{
	    pt[i] = 0.0; 
	    for (j = 0; j < 3; j++)
	    {
		pt[i] += roty_matrix[i][j]*Coords(p)[j];
	    }
	}
	for (i = 0; i < 3; i++)
	{
	    Coords(p)[i] = 0.0;
	    for (j = 0; j < 3; j++)
	    {
		Coords(p)[i] += rotz_matrix[i][j]*pt[j];
	    }
	}
	for (i = 0; i < 3; i++)
	    Coords(p)[i] += center[i];
}	/* end rotate_point_with_spherical_angle */

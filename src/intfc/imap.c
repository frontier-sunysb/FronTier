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
LOCAL void change_vertex_of_tris(POINT*,TRI*,POINT*);
LOCAL POINT* crdsToPoint(double*,BOND*,CURVE*);
LOCAL BOND* PointOnCurve(double*,CURVE*);
LOCAL boolean PointOnBond(double*,BOND*);

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

EXPORT void I_RotatePointAboutAxis(
	POINT *p,
	double *dir,		/* dir vector is a unit vector */
	double *axis,		/* a point on the axis */
	double phi)
{
	double v[MAXD],vt[MAXD],vn[MAXD],nor[MAXD],cnor[MAXD];
	double dot_prod = 0;
	double mag_vn;
	int i;

	for (i = 0; i < 3; ++i)
	{
	    /* from position vector relative vector */
	    v[i] = Coords(p)[i] - axis[i];
	    dot_prod += v[i]*dir[i];
	}
	/* vt and vn are tangential and normal components of v */
	for (i = 0; i < 3; ++i)
	{
	    vt[i] = dot_prod*dir[i];
	    vn[i] = v[i] - vt[i];
	}
	mag_vn = Mag3d(vn);
	if (mag_vn == 0.0) return; /* no normal component to rotate */
	/* calculating unit normal vector */
	for (i = 0; i < 3; ++i)
	{
	    nor[i] = vn[i]/mag_vn;
	    vn[i] = 0.0;
	}
	/* calculating unit co-normal vector */
	Cross3d(dir,nor,cnor);
	for (i = 0; i < 3; ++i)
	{
	    /* rotate the normal vector */
	    vn[i] += mag_vn*cos(phi)*nor[i] + mag_vn*sin(phi)*cnor[i];
	    /* add the tangential vector */
	    v[i] = vn[i] + vt[i];
	    /* recover the position vector */
	    Coords(p)[i] = v[i] + axis[i];
	}
}	/* end I_RotatePointAboutAxis */

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

	reset_sort_status(intfc);

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

EXPORT void I_SphericalRotateInteriorSurfPoints(
	SURFACE *surf,
        double *center,
        double phi,
        double theta)
{
	TRI *t;
	POINT *p;
	int i;
	boolean first = YES;

	surf_tri_loop(surf,t)
	{
	    for (i = 0; i < 3; ++i)
	    {
		p = Point_of_tri(t)[i];
		sorted(p) = NO;
	    }
	}
	surf_tri_loop(surf,t)
	{
	    for (i = 0; i < 3; ++i)
	    {
		p = Point_of_tri(t)[i];
		if (sorted(p) || Boundary_point(p)) continue;
		rotate_point_with_spherical_angle(p,center,phi,theta,first);
		if (first == YES) first = NO;
		sorted(p) = YES;
	    }
	}
}	/* end I_SphericalRotateInteriorSurfPoints */

EXPORT void I_SphericalRotateInteriorCurvePoints(
	CURVE *curve,
        double *center,
        double phi,
        double theta)
{
	boolean first = YES;
	POINT *p;
	BOND *b;
	for (b = curve->first; b != curve->last; b = b->next)
	{
	    p = b->end;
	    rotate_point_with_spherical_angle(p,center,phi,theta,first);
	    if (first == YES)	first = NO;
	}
}	/* end I_SphericalRotateInteriorCurvePoints */

EXPORT void I_SphericalRotateInteriorNodePoints(
	NODE *node,
        double *center,
        double phi,
        double theta)
{
	rotate_point_with_spherical_angle(node->posn,center,phi,theta,YES);
}	/* end I_SphericalRotateInteriorNodePoints */

EXPORT int I_NumOfSurfInteriorPoints(SURFACE *surf)
{
	TRI *tri;
        POINT *p;
        int i,n;

	surf_tri_loop(surf,tri)
	{
	    for (i = 0; i < 3; ++i)
		sorted(Point_of_tri(tri)[i]) = NO;
	}
	n = 0;
	surf_tri_loop(surf,tri)
	{
	    for (i = 0; i < 3; ++i)
	    {
		p = Point_of_tri(tri)[i];
		if (sorted(p) || Boundary_point(p)) continue;
		++n;
		sorted(p) = YES;
	    }
	}
	return n;
}	/* end I_NumOfSurfInteriorPoints */

EXPORT int I_NumOfCurveInteriorPoints(CURVE *curve)
{
	return curve->num_points - 2;
}	/* end I_NumOfCurveInteriorPoints */

EXPORT void I_FoldSurface(
	SURFACE *surf,
	double *dir,
	double *axis,
	double angle,
	SIDE side,
	boolean first)
{
	TRI *tri;
	POINT *p;
        int i,j;
	double pv[MAXD],cx[MAXD];
	const double *nor;
	double prod;
	double tol = 0.1;

	if (first == YES)
	{
	    surf_tri_loop(surf,tri)
	    {
	    	for (i = 0; i < 3; ++i)
		    sorted(Point_of_tri(tri)[i]) = NO;
	    }
	}
	surf_tri_loop(surf,tri)
	{
	    nor = Tri_normal(tri);
	    for (i = 0; i < 3; ++i)
	    {
		p = Point_of_tri(tri)[i];
		if (sorted(p)) continue;
		for (j = 0; j < 3; ++j)
		    pv[j] = Coords(p)[j] - axis[j];
		Cross3d(pv,dir,cx);
		prod = Mag3d(cx);
		if (prod == 0.0) continue;
		for (j = 0; j < 3; ++j)
		    cx[j] /= prod;
		prod = Dot3d(cx,nor);
		if ((prod < 0  && side == POSITIVE_SIDE) ||
		    (prod >= 0 && side == NEGATIVE_SIDE))
		    continue;
		sorted(p) = YES;
		I_RotatePointAboutAxis(p,dir,axis,angle);
		
	    }
	}
}	/* end I_FoldSurface */

/***********************************************************************
*       Sew surface from crds_start to crds_end. The sewing line must  *
*       be along existing curves with tolerance, else return NO.       *
***********************************************************************/

static boolean find_sewing_segments(SURFACE*,double*,double*,BOND**,
				BOND**,BOND**,BOND**);
static boolean find_seg_start_and_end(CURVE*,ORIENTATION,double*,double*,
				double*,BOND**,BOND**);
static boolean merge_tris_on_bonds(BOND*,BOND*,BOND*,BOND*,SURFACE*);

EXPORT boolean I_SewSurface(
	SURFACE *surf,
	double *crds_start,
	double *crds_end)
{
	BOND *bps,*bpe;
	BOND *bns,*bne;
	boolean status;
	BOND *b;

	if (debugging("sewing"))
	    (void) printf("Entering I_SewSurface()\n");
	status = find_sewing_segments(surf,crds_start,crds_end,&bps,&bpe,
				&bns,&bne);
	if (!status)
	{
	    (void) printf("find_sewing_segments status = %d\n",status);
	    clean_up(0);
	}
	if (debugging("sewing"))
	{
	    (void) printf("find_sewing_segments status = %d\n",status);
	    (void) printf("bps = %p  bpe = %p\n",(void*)bps,(void*)bpe);
	    (void) printf("bns = %p  bne = %p\n",(void*)bns,(void*)bne);
	    (void) printf("Positive segment:\n");
	    (void) printf("%f %f %f\n",Coords(bps->start)[0],
				Coords(bps->start)[1],Coords(bps->start)[2]);
	    (void) printf("%f %f %f\n",Coords(bpe->end)[0],Coords(bpe->end)[1],
				Coords(bpe->end)[2]);
	    (void) printf("Negative segment:\n");
	    (void) printf("%f %f %f\n",Coords(bns->end)[0],Coords(bns->end)[1],
				Coords(bns->end)[2]);
	    (void) printf("%f %f %f\n",Coords(bne->start)[0],
				Coords(bne->start)[1],Coords(bne->start)[2]);
	}
	status = merge_tris_on_bonds(bps,bpe,bns,bne,surf);
	if (debugging("sewing"))
	{
	    (void) printf("merge_tris_on_bonds status = %d\n",status);
	}
	if (debugging("sewing"))
	{
	    (void) printf("Leaving I_SewSurface()\n");
	    (void) printf("Check consistency of interface:\n");
	    if (consistent_interface(surf->interface))
		(void) printf("Interface is consistent!\n");
	}
}	/* end I_SewSurface */

static boolean merge_tris_on_bonds(
	BOND *bps,
	BOND *bpe,
	BOND *bns,
	BOND *bne,
	SURFACE *surf)
{
	INTERFACE *intfc = surf->interface;
	POINT *p,*ps,*pe;
	CURVE **curves,*curve;
	BOND *bp,*bn, *bond;
	BOND *bps_prev,*bpe_next,*bns_next,*bne_prev;
	CURVE *cp,*cps_prev,*cpe_next,*cn,*cns_next,*cne_prev;
	NODE *nps,*npe,*nns,*nne;
	int i,j;
	POINT *pt_ps,*pt_pe,*pt_ns,*pt_ne;
	TRI *trip,*trin;

	bps_prev = bps->prev;
	bpe_next = bpe->next;
	bns_next = bns->next;
	bne_prev = bne->prev;

	p = bps->start;
	curve = curve_of_bond(bps,intfc);
	if (is_closed_curve(curve))
	{
	    move_closed_loop_node(curve,bps);
	}
	else
	{
	    curves = split_curve(p,bps,curve,NO_COMP,NO_COMP,NO_COMP,NO_COMP);
	    curve = curves[1];
	}
	p = bpe->end;
	curves = split_curve(p,bpe,curve,NO_COMP,NO_COMP,NO_COMP,NO_COMP);
	cp = curves[0];

	p = bns->end;
	curve = curve_of_bond(bns,intfc);
	curves = split_curve(p,bns,curve,NO_COMP,NO_COMP,NO_COMP,NO_COMP);
	curve = curves[0];
	p = bne->start;
	curves = split_curve(p,bne,curve,NO_COMP,NO_COMP,NO_COMP,NO_COMP);
	cn = curves[1];

        curve_bond_loop(cp,bond)
        {
            p = bond->end;
            crdsToPoint(Coords(p), PointOnCurve(Coords(p), cn) ,cn);
        }
        curve_bond_loop(cn,bond)
        {
            p = bond->end;
            crdsToPoint(Coords(p), PointOnCurve(Coords(p), cp) ,cp);
        }

	if (debugging("sewing"))
	{
	    print_curve(cp);
	    print_curve(cn);
	}

	bps = cp->first;
	bns = cn->last;
	for (bp = bps, bn = bns; bp != NULL; bp = bp->next, bn = bn->prev)
	{
	    trip = (*Btris(bp))->tri;
	    trin = (*Btris(bn))->tri;
	    pt_ps = bp->start;
	    pt_ns = bn->end;
	    pt_pe = bp->end;
	    pt_ne = bn->start;
	    change_vertex_of_tris(pt_ps,trin,pt_ns);
	    if (bp->next == NULL)
	    {
	    	change_vertex_of_tris(pt_pe,trin,pt_ne);
	    }
	    else
		Boundary_point(pt_pe) = NO;
	}
	for (bp = bps, bn = bns; bp != NULL; bp = bp->next, bn = bn->prev)
	{
	    trip = (*Btris(bp))->tri;
	    trin = (*Btris(bn))->tri;
	    for (i = 0; i < 3; ++i)
	    for (j = 0; j < 3; ++j)
	    {
		if (Point_of_tri(trip)[i] == Point_of_tri(trin)[(j+1)%3] &&
		    Point_of_tri(trip)[(i+1)%3] == Point_of_tri(trin)[j])
		{
		    Tri_on_side(trip,i) = trin;
		    Tri_on_side(trin,j) = trip;
		    set_side_bdry(Boundary_tri(trip),i,NO);
		    set_side_bdry(Boundary_tri(trin),j,NO);
		}
	    }
	}
	cps_prev = curve_of_bond(bps_prev,intfc);
	cpe_next = curve_of_bond(bpe_next,intfc);
	cns_next = curve_of_bond(bns_next,intfc);
	cne_prev = curve_of_bond(bne_prev,intfc);
	nps = cps_prev->end;
	npe = cpe_next->start;
	nns = cns_next->start;
	nne = cne_prev->end;

	change_node_of_curve(cns_next,POSITIVE_ORIENTATION,nps);
	change_node_of_curve(cne_prev,NEGATIVE_ORIENTATION,npe);

	trin = (*Btris(bns_next))->tri;
	delete_from_pointers(*Btris(bns_next),&Btris(bns_next));
	link_tri_to_bond(NULL,trin,surf,bns_next,cns_next);
	trin = (*Btris(bne_prev))->tri;
	delete_from_pointers(*Btris(bne_prev),&Btris(bne_prev));
	link_tri_to_bond(NULL,trin,surf,bne_prev,cne_prev);

	delete_curve(cp);
	delete_curve(cn);
	delete_node(nns);
	delete_node(nne);

	cps_prev = join_curves(cps_prev,cns_next,NO_COMP,NO_COMP,NULL);
	cne_prev = join_curves(cne_prev,cpe_next,NO_COMP,NO_COMP,NULL);
//	delete_node(nps);
//	delete_node(npe);

	return YES;
}	/* end merge_tris_on_bonds */
	

static boolean find_sewing_segments(
	SURFACE *surf,
	double *crds_start,
	double *crds_end,
	BOND **bps,
	BOND **bpe,
	BOND **bns,
	BOND **bne)
{
	int i;
	CURVE **c;
	double dir[MAXD];
	double len;
	boolean status,pos_status,neg_status;
	BOND *bs,*be;

	for (i = 0; i < 3; ++i)
	    dir[i] = crds_end[i] - crds_start[i];
	len = Mag3d(dir);
	for (i = 0; i < 3; ++i) dir[i] /= len;

	*bps = *bpe = *bns = *bne = NULL;
	pos_status = neg_status = NO;
	surf_neg_curve_loop(surf,c)
	{
	    status = find_seg_start_and_end(*c,POSITIVE_ORIENTATION,						dir,crds_start,crds_end,&bs,&be);
	    if (status == YES)
	    {
		*bps = bs;
		*bpe = be;
		pos_status = YES;
	    }
	    status = find_seg_start_and_end(*c,NEGATIVE_ORIENTATION,						dir,crds_start,crds_end,&bs,&be);
	    if (status == YES)
	    {
		*bns = bs;
		*bne = be;
		neg_status = YES;
	    }
	}
	if (pos_status == YES && neg_status == YES)
	    return YES;
	surf_pos_curve_loop(surf,c)
	{
	    status = find_seg_start_and_end(*c,POSITIVE_ORIENTATION,						dir,crds_start,crds_end,&bs,&be);
	    if (status == YES)
	    {
		*bps = bs;
		*bpe = be;
		pos_status = YES;
	    }
	    status = find_seg_start_and_end(*c,NEGATIVE_ORIENTATION,						dir,crds_start,crds_end,&bs,&be);
	    if (status == YES)
	    {
		*bns = bs;
		*bne = be;
		neg_status = YES;
	    }
	}
	if (pos_status == YES && neg_status == YES)
	    return YES;
	return NO;
}	/* end find_sewing_segments */

static boolean find_seg_start_and_end(
	CURVE *curve,
	ORIENTATION orient,
	double *dir,
	double *crds_start,
	double *crds_end,
	BOND **bs,
	BOND **be)
{
	BOND *b,*start,*end;
	int i,j,i_max;
	double pv[MAXD],bv[MAXD],dir_max;
	double *p;
	double lambda,tol = 1.0e-8;
	start = end = NULL;
	boolean positive;

	dir_max = 0.0;
	for (i = 0; i < 3; ++i)
	{
	    if (fabs(dir_max) < fabs(dir[i]))
	    {
		dir_max = dir[i];
		i_max = i;
	    }
	}

	p = crds_start;
	for (b = curve->first; b != NULL; b = b->next)
	{
	    dir_max = (Coords(b->end)[i_max] - Coords(b->start)[i_max])/
			bond_length(b);
	    if (orient == POSITIVE_ORIENTATION && 
		fabs(dir_max - dir[i_max]) > tol) 
		continue;
	    if (orient == NEGATIVE_ORIENTATION && 
		fabs(dir_max + dir[i_max]) > tol) 
		continue;

	    lambda = (p[i_max] - Coords(b->start)[i_max])/
		(Coords(b->end)[i_max] - Coords(b->start)[i_max]);
	    if (lambda < 0.0 || lambda > 1.0) continue;
	    for (i = 0; i < 3; ++i)
	    {
		pv[i] = lambda*Coords(b->end)[i] + 
			(1.0 - lambda)*Coords(b->start)[i];
		if (fabs(pv[i] - p[i]) > tol) break;
	    }
            if (i < 3) continue;
	    start = b;
	    break;
	}

	p = crds_end;
	for (b = curve->first; b != NULL; b = b->next)
	{
	    dir_max = (Coords(b->end)[i_max] - Coords(b->start)[i_max])/
			bond_length(b);
	    if (orient == POSITIVE_ORIENTATION && 
		fabs(dir_max - dir[i_max]) > tol) 
		continue;
	    if (orient == NEGATIVE_ORIENTATION && 
		fabs(dir_max + dir[i_max]) > tol) 
		continue;

	    lambda = (p[i_max] - Coords(b->start)[i_max])/
		(Coords(b->end)[i_max] - Coords(b->start)[i_max]);
	    if (lambda < 0.0 || lambda > 1.0) continue;
	    for (i = 0; i < 3; ++i)
	    {
		pv[i] = lambda*Coords(b->end)[i] + 
			(1.0 - lambda)*Coords(b->start)[i];
		if (fabs(pv[i] - p[i]) > tol) break;
	    }
	    if (i < 3) continue;
	    if (orient == NEGATIVE_ORIENTATION && fabs(1.0-lambda) < tol)
		end = b->next;
	    else
		end = b;
	    break;
	}
	*bs = start;	*be = end;
	if (start == NULL || end == NULL) return NO;
	return YES;
}	/* end find_seg_start_and_end */

LOCAL void change_vertex_of_tris(
	POINT *pnew,
	TRI *tri,
	POINT *pold)
{
	int i,j,nt;
	TRI **tris;

	nt = set_tri_list_around_point(pold,tri,&tris,pold->interface);
	for (i = 0; i < nt; ++i)
	{
	    for (j = 0; j < 3; ++j)
	    {
		if (Point_of_tri(tris[i])[j] == pold)
		    Point_of_tri(tris[i])[j] = pnew;
	    }
	}
}	/* end change_vertex_of_tris */

LOCAL POINT* crdsToPoint(
        double *crds,
        BOND *b,
        CURVE *c)
{
        int i;
        POINT *s = b->start;
        POINT *e = b->end;
        POINT *p;
        double A[3], B[3], C[3];

        for (i = 0; i < 3; i++)
        {
            A[i] = crds[i] - Coords(s)[i];
            B[i] = crds[i] - Coords(e)[i];
        }
        if (Mag3d(A) < 1e-3)
            return s;
        if (Mag3d(B) < 1e-3)
            return e;
        p = Point(crds);
        insert_point_in_bond(p,b,c);
        return p;
}       /* end crdsToPoint */

LOCAL BOND* PointOnCurve(
        double *crds,
        CURVE *c)
{
        BOND *b;

        curve_bond_loop(c,b)
            if (PointOnBond(crds,b))
                return b;
        return NULL;
}       /* end PointOnCurve */

LOCAL boolean PointOnBond(
        double *crds,
        BOND *b)
{
        int i;
        POINT *s = b->start;
        POINT *e = b->end;
        double A[3], B[3], C[3];

        for (i = 0; i < 3; i++)
        {
            A[i] = crds[i] - Coords(s)[i];
            B[i] = crds[i] - Coords(e)[i];
            C[i] = Coords(s)[i] - Coords(e)[i];
        }
        if (fabs(Mag3d(A) + Mag3d(B) - Mag3d(C)) < 1e-8)
            return YES;
        return NO;
}       /* end PointOnBond */

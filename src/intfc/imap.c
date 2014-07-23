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
				BOND**,BOND**,BOND**,CURVE**,CURVE**);
static boolean find_seg_start_and_end(CURVE*,ORIENTATION,double*,double*,
				double*,BOND**,BOND**);
EXPORT boolean I_SewSurface(
	SURFACE *surf,
	double *crds_start,
	double *crds_end)
{
	BOND *bps,*bpe;
	BOND *bns,*bne;
	CURVE *pos_curve,*neg_curve;
	boolean status;
	BOND *b;

	printf("Entering I_SewSurface()\n");
	status = find_sewing_segments(surf,crds_start,crds_end,
			&bps,&bpe,&bns,&bne,&pos_curve,&neg_curve);
	printf("status = %d\n",status);
	printf("bps = %d  bpe = %d  pos_curve = %d\n",bps,bpe,pos_curve);
	printf("bns = %d  bne = %d  neg_curve = %d\n",bns,bne,neg_curve);
	printf("Positive segment:\n");
	if (bps == NULL)
	{
	    printf("bps = %d bpe = %d  bns = %d bne = %d\n",bps,bpe,bns,bns);
	    clean_up(0);
	}
	printf("%f %f %f\n",Coords(bps->start)[0],Coords(bps->start)[1],
					Coords(bps->start)[2]);
	for (b = bps; b != bpe->next; b = b->next);
	    printf("%f %f %f\n",Coords(b->end)[0],Coords(b->end)[1],
					Coords(b->end)[2]);
	printf("Negative segment:\n");
	printf("%f %f %f\n",Coords(bns->end)[0],Coords(bns->end)[1],
					Coords(bns->end)[2]);
	for (b = bns; b != bne->prev; b = b->prev);
	    printf("%f %f %f\n",Coords(b->start)[0],Coords(b->start)[1],
					Coords(b->start)[2]);
	clean_up(0);
}	/* end I_SewSurface */

static boolean find_sewing_segments(
	SURFACE *surf,
	double *crds_start,
	double *crds_end,
	BOND **bps,
	BOND **bpe,
	BOND **bns,
	BOND **bne,
	CURVE **pos_curve,
	CURVE **neg_curve)
{
	int i;
	CURVE **c;
	double dir[MAXD];
	double len;
	boolean status;
	BOND *bs,*be;

	for (i = 0; i < 3; ++i)
	    dir[i] = crds_end[i] - crds_start[i];
	len = Mag3d(dir);
	for (i = 0; i < 3; ++i) dir[i] /= len;

	*pos_curve = *neg_curve = NULL;
	*bps = *bpe = *bns = *bne = NULL;
	surf_neg_curve_loop(surf,c)
	{
	    status = find_seg_start_and_end(*c,POSITIVE_ORIENTATION,						dir,crds_start,crds_end,&bs,&be);
	    if (status == YES)
	    {
		*bps = bs;
		*bpe = be;
		*pos_curve = *c;
	    }
	    status = find_seg_start_and_end(*c,NEGATIVE_ORIENTATION,						dir,crds_start,crds_end,&bs,&be);
	    if (status == YES)
	    {
		*bns = bs;
		*bne = be;
		*neg_curve = *c;
	    }
	}
	if (*pos_curve != NULL && *neg_curve != NULL)
	    return YES;
	surf_pos_curve_loop(surf,c)
	{
	    status = find_seg_start_and_end(*c,POSITIVE_ORIENTATION,						dir,crds_start,crds_end,&bs,&be);
	    if (status == YES)
	    {
		*bps = bs;
		*bpe = be;
		*pos_curve = *c;
	    }
	    status = find_seg_start_and_end(*c,NEGATIVE_ORIENTATION,						dir,crds_start,crds_end,&bs,&be);
	    if (status == YES)
	    {
		*bns = bs;
		*bne = be;
		*neg_curve = *c;
	    }
	}
	if (*pos_curve != NULL && *neg_curve != NULL)
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
	double lambda,tol = 1.0e-5;
	start = end = NULL;
	p = crds_start;

	dir_max = -1.0;
	for (i = 0; i < 3; ++i)
	{
	    if (dir_max < dir[i])
	    {
		dir_max = dir[i];
		i_max = i;
	    }
	}

	switch (orient)
	{
	case POSITIVE_ORIENTATION:
	    printf("positive loop:\n");
	    for (b = curve->first; b != NULL; b = b->next)
	    {
		
		dir_max = (Coords(b->end)[i_max] - Coords(b->start)[i_max])/
				bond_length(b);
		if (fabs(dir_max - dir[i_max]) > tol) 
		{
		    p = crds_start;
		    continue;
		}
		lambda = (p[i_max] - Coords(b->start)[i_max])/
			(Coords(b->end)[i_max] - Coords(b->start)[i_max]);
		for (i = 0; i < 3; ++i)
		{
		    pv[i] = Coords(b->start)[i] + lambda*dir[i];
		    if (fabs(pv[i] - p[i]) > tol) break;
		}
		if (i < 3)
		{
		    p = crds_start;
		    continue;
		}
		if (start == NULL)
		{
		    start = b;
		    p = crds_end;
		}
		else
		{
		    end = b;
		    break;
		}
	    }
	    print_bond(start);
	    print_bond(end);
	    clean_up(0);
	    break;
	case NEGATIVE_ORIENTATION:
	    printf("negative loop:\n");
	    for (b = curve->last; b != NULL; b = b->prev)
	    {
		for (i = 0; i < 3; ++i)
		{
		    if (fabs(dir[i]) < tol) continue;
		    bv[i] = (Coords(b->start)[i] - Coords(b->end)[i])/
				bond_length(b);
		    pv[i] = (p[i] - Coords(b->end)[i])/bond_length(b);
		    if (fabs(bv[i] - dir[i]) > tol ||
		        pv[i]/bv[i] < 0.0 || pv[i]/bv[i] > 1.0)
		    {
		    	start = NULL;		/* restart */
			p = crds_start;
			break;
		    }
		}
		if (i < 3) continue;
		if (start == NULL)
		{
		    start = b;
		    p = crds_end;
		}
		else
		    end = b;
		if (end != NULL) break;
	    }
	    break;
	}
	*bs = start;	*be = end;
	if (start == NULL || end == NULL) return NO;
	return YES;
}	/* end find_seg_start_and_end */

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
*				fmap2.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*/

#include <front/fdecs.h>		/* includes int.h, table.h */

LOCAL double cell_frac_comp2(double***,double***,COMPONENT**,COMPONENT);
LOCAL double corner_cell_frac(double***,double***,COMPONENT**,
			COMPONENT,boolean);
LOCAL void rot_clockwise_90(double***,double***,int**);
LOCAL void FT_ComputeGridVolumeFraction2d(Front*,COMPONENT,POINTER*);
LOCAL double FT_ComputeTotalVolumeFraction2d(Front*,COMPONENT);
LOCAL void forward_curve_seg_len_constr(CURVE*,BOND*,BOND*,int,double);
LOCAL void backward_curve_seg_len_constr(CURVE*,BOND*,BOND*,int,double);
static boolean find_nearest_ring2_cell_with_comp(int*,int*,int*,int*,int);

EXPORT double FT_ComputeTotalVolumeFraction(
	Front *front,
	COMPONENT comp_of_vol)
{
	switch (front->rect_grid->dim)
	{
	case 2:
	    return FT_ComputeTotalVolumeFraction2d(front,comp_of_vol);
	case 3:
	    screen("FT_ComputeTotalVolumeFraction3d() not implemented!\n");
	    return 0.0;
	}
}	/* end FT_ComputeTotalVolumeFraction */

LOCAL double FT_ComputeTotalVolumeFraction2d(
	Front *front,
	COMPONENT comp_of_vol)
{
	double **grid_vol_frac;
	INTERFACE *grid_intfc = front->grid_intfc;
	RECT_GRID *top_grid = &topological_grid(grid_intfc);
	int *lbuf = front->rect_grid->lbuf;
	int *ubuf = front->rect_grid->ubuf;
	int *top_gmax = top_grid->gmax;
	double total_vol = 0.0;
	int i,j,imin,jmin,imax,jmax;

	FT_ComputeGridVolumeFraction2d(front,comp_of_vol,
			(POINTER*)&grid_vol_frac);
	imin = (lbuf[0] == 0) ? 1 : lbuf[0];
	jmin = (lbuf[1] == 0) ? 1 : lbuf[1];
        imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
        jmax = (ubuf[1] == 0) ? top_gmax[1] - 1 : top_gmax[1] - ubuf[1];
	for (i = imin; i <= imax; ++i)
	for (j = jmin; j <= jmax; ++j)
	    total_vol += grid_vol_frac[i][j];

	pp_global_sum(&total_vol,1);
	FT_FreeThese(1,grid_vol_frac);
	return total_vol;
}	/* end FT_ComputeTotalVolumeFraction2d */

EXPORT void FT_ComputeGridVolumeFraction(
	Front *front,
	COMPONENT comp_of_vol,
	POINTER *grid_vol_frac)
{
	switch (front->rect_grid->dim)
	{
	case 2:
	    FT_ComputeGridVolumeFraction2d(front,comp_of_vol,grid_vol_frac);
	    return;
	case 3:
	    screen("FT_ComputeVolumeFraction3d() not implemented!\n");
	    return;
	}
}	/* end FT_ComputeGridVolumeFraction */

LOCAL void FT_ComputeGridVolumeFraction2d(
	Front *front,
	COMPONENT comp_of_vol,
	POINTER *grid_vol_frac)
{
	int i,j,i1,j1,ii,jj;
	INTERFACE *grid_intfc = front->grid_intfc;
	RECT_GRID *top_grid = &topological_grid(grid_intfc);
	int *lbuf = front->rect_grid->lbuf;
	int *ubuf = front->rect_grid->ubuf;
	int *gmax = top_grid->gmax;
	double *L = top_grid->L;
	double *U = top_grid->U;
	double *h = top_grid->h;
	int dim = grid_intfc->dim;
	struct Table *T = table_of_interface(grid_intfc);
	COMPONENT *top_comp = T->components;
	double **vol_frac_2d;
	static double ***corner,***edge_crx;
	static COMPONENT **corner_comp;
	int ic,num_comp,icoords[MAXD];
	POINTER state;
	HYPER_SURF *hs;

	if (corner == NULL)
	{
	    FT_TriArrayMemoryAlloc((POINTER*)&corner,2,2,MAXD,sizeof(double));
	    FT_TriArrayMemoryAlloc((POINTER*)&edge_crx,2,2,MAXD,sizeof(double));
	    FT_MatrixMemoryAlloc((POINTER*)&corner_comp,2,2,sizeof(int));
	}
	FT_MatrixMemoryAlloc((POINTER*)&vol_frac_2d,gmax[0],gmax[1],
				sizeof(double));
	for (i = 0; i < gmax[0]; ++i)
	for (j = 0; j < gmax[1]; ++j)
	{
	    num_comp = 0;
	    for (i1 = 0; i1 < 2; ++i1)
	    for (j1 = 0; j1 < 2; ++j1)
	    {
		ii = i+i1;
		jj = j+j1;
	    	corner[i1][j1][0] = L[0] + ii*h[0];
	    	corner[i1][j1][1] = L[0] + jj*h[0];
	    	ic = d_index2d(ii,jj,gmax);
	    	corner_comp[i1][j1] = top_comp[ic];
		if (top_comp[ic] == comp_of_vol) num_comp++;
	    }
	
	    for (i1 = 0; i1 < 2; ++i1)
	    {
	    	if (corner_comp[0][i1] != corner_comp[1][i1])
		{
		    icoords[0] = i;
		    icoords[1] = j + i1;
		    FT_StateStructAtGridCrossing(front,grid_intfc,icoords,EAST,
				corner_comp[0][i1],&state,&hs,edge_crx[0][i1]);
		}
	    	if (corner_comp[i1][0] != corner_comp[i1][1])
		{
		    icoords[0] = i + i1;
		    icoords[1] = j;
		    FT_StateStructAtGridCrossing(front,grid_intfc,icoords,NORTH,
				corner_comp[i1][0],&state,&hs,edge_crx[1][i1]);
		}
	    }
	    switch (num_comp)
	    {
	    case 0:
		vol_frac_2d[i][j] = 0.0;
		break;
	    case 1:
		vol_frac_2d[i][j] = corner_cell_frac(corner,edge_crx,
				corner_comp,comp_of_vol,NO);
		break;
	    case 2:
		vol_frac_2d[i][j] = cell_frac_comp2(corner,edge_crx,
				corner_comp,comp_of_vol);
		break;
	    case 3:
		vol_frac_2d[i][j] = h[0]*h[1] - corner_cell_frac(corner,
				edge_crx,corner_comp,comp_of_vol,YES);
		break;
	    case 4:
		vol_frac_2d[i][j] = h[0]*h[1];
		break;
	    }
	}
	*grid_vol_frac = (POINTER)vol_frac_2d;
}	/* end FT_ComputeGridVolumeFraction2d */

#define	the_corner(comp,comp_of_vol,use_complement)  \
	(((comp) == (comp_of_vol) && !(use_complement)) || \
	 ((comp) != (comp_of_vol) && (use_complement)))

LOCAL double corner_cell_frac(
	double ***corner,
	double ***edge,
	COMPONENT **comp,
	COMPONENT comp_of_vol,
	boolean use_complement)
{
	int i;
	double a,b;
	/* All triangle cases */
	for (i = 0; i < 4; ++i)
	{
	    if (the_corner(comp[0][0],comp_of_vol,use_complement))
	    {
	    	a = edge[0][0][0] - corner[0][0][0];
	    	b = edge[1][0][1] - corner[0][0][1];
		return 0.5*a*b;
	    }
	    rot_clockwise_90(corner,edge,comp);
	}
}	/* end corner_cell_frac */

LOCAL double cell_frac_comp2(
	double ***corner,
	double ***edge,
	COMPONENT **comp,
	COMPONENT comp_of_vol)
{
	int i;
	double a,b,h;
	double volume = 0.0;
	/* Traingle cases */
	if ((the_corner(comp[0][0],comp_of_vol,NO) &&
	     the_corner(comp[1][1],comp_of_vol,NO)) ||
	    (the_corner(comp[1][0],comp_of_vol,NO) &&
	     the_corner(comp[0][1],comp_of_vol,NO)))
	{
	    for (i = 0; i < 4; ++i)
	    {
	    	if (the_corner(comp[0][0],comp_of_vol,NO))
	    	{
	    	    a = edge[0][0][0] - corner[0][0][0];
	    	    b = edge[1][0][1] - corner[0][0][1];
		    volume += 0.5*a*b;
	    	}
	    	rot_clockwise_90(corner,edge,comp);
	    }
	    return volume;
	}
	/* Trapezoidal cases */
	else 
	{
	    for (i = 0; i < 4; ++i)
	    {
	    	if (the_corner(comp[0][0],comp_of_vol,NO) &&
	    	    the_corner(comp[1][0],comp_of_vol,NO))
		{
		    a = edge[1][0][1] - corner[0][0][1];
		    b = edge[1][1][1] - corner[1][0][1];
		    h = corner[1][0][0] - corner[0][0][0];
		    return 0.5*(a + b)*h;
		}
	    	rot_clockwise_90(corner,edge,comp);
	    }
	}
}	/* end cell_frac_comp2 */

LOCAL void rot_clockwise_90(
	double ***corner,
	double ***edge,
	int **comp)
{
	COMPONENT tmp_comp;
	double tmp_corner[MAXD];
	double tmp_edge[MAXD];
	int i;

	tmp_comp = comp[0][0];
	comp[0][0] = comp[1][0];
	comp[1][0] = comp[1][1];
	comp[1][1] = comp[0][1];
	comp[0][1] = tmp_comp;
	for (i = 0; i < 2; ++i)
	{
	    tmp_corner[i] = corner[0][0][i];
	    corner[0][0][i] = corner[1][0][i];
	    corner[1][0][i] = corner[1][1][i];
	    corner[1][1][i] = corner[0][1][i];
	    corner[0][1][i] = tmp_corner[i];
	    tmp_edge[i] = edge[0][0][i];
	    edge[0][0][i] = edge[1][1][i];
	    edge[1][1][i] = edge[0][1][i];
	    edge[0][1][i] = edge[1][0][i];
	    edge[1][0][i] = tmp_edge[i];
	}
}	/* end rot_clockwise_90 */

EXPORT double FT_CurveLength(
        CURVE *c)
{
	return curve_length(c);
}	/* end FT_CurveLength */

EXPORT void FT_CurveSegLengthConstr(
        CURVE *c,
        BOND *bs,
        BOND *be,
        int nbds,
        double seg_length,
        REDISTRIBUTION_DIRECTION dir)
{
        switch (dir)
        {
        case FORWARD_REDISTRIBUTION:
	    forward_curve_seg_len_constr(c,bs,be,nbds,seg_length);
            return;
        case BACKWARD_REDISTRIBUTION:
	    backward_curve_seg_len_constr(c,bs,be,nbds,seg_length);
            return;
        }
}       /* end FT_CurveSegLengthConstr */


static	void forward_curve_seg_len_constr(
	CURVE		*c,
	BOND		*bs,
	BOND		*be,
	int		nbds,
	double		seg_len)
{
	BOND		*b, *bstart, *bend;
	double		b_len, sc_len, offset, total_offset, s, oms;
	double		coords[MAXD];
	double		space_tol;
	boolean		reset_bend;

	b_len = seg_len/(double)nbds;
	space_tol = b_len*1.0e-8;


	offset = b_len;		total_offset = seg_len;
	bstart = bs;		bend = be;
	while (bstart != bend)
	{
	    b = bstart;
	    while ((sc_len = bond_length(b)) < offset - space_tol)
	    {
	    	offset -= sc_len;
		if (b == bend) break;
	    	b = b->next;
	    }
	    if (sc_len > offset + space_tol)
	    {
	    	s = offset/sc_len;
		interpolate_coords(coords,Coords(b->start),Coords(b->end),s,2);
		if (b == bend) reset_bend = YES;
		else reset_bend = NO;
	    	if (insert_point_in_bond(Point(coords),b,c) !=
		    FUNCTION_SUCCEEDED)
	    	{
	    	    screen("ERROR in forward_curve_seg_len_constr(), "
	    	           "insert_point_in_bond failed\n");
	    	    clean_up(ERROR);
	    	}
		if (reset_bend) bend = b->next;
	    }
	    replace_curve_seg_by_bond(c,bstart,b);
	    total_offset -= b_len;
	    bstart = bstart->next;
	    offset = b_len;
	    if (total_offset < b_len - space_tol) break;
	}
	if (bstart == bend)
	{
	    if (total_offset < space_tol)
	    {
		Coords(bend->end)[0] = Coords(bstart->start)[0];
		Coords(bend->end)[1] = Coords(bstart->start)[1];
		bstart = bstart->prev;
	    	replace_curve_seg_by_bond(c,bstart,bend);
	    }
	    else
	    {
	    	sc_len = bond_length(bend);
	    	s = total_offset/sc_len;      oms = 1.0 - s;
	    	Coords(bend->end)[0] = Coords(bend->start)[0] + s*
			(Coords(bend->end)[0] - Coords(bend->start)[0]);
	    	Coords(bend->end)[1] = Coords(bend->start)[1] + s*
			(Coords(bend->end)[1] - Coords(bend->start)[1]);
	    	set_bond_length(bend,2);
		b = bend;
		while (total_offset/b_len > 2.0 - space_tol)
		{
	    	    s = b_len/total_offset;
		    interpolate_coords(coords,Coords(b->start),
					Coords(b->end),s,2);
	    	    if (insert_point_in_bond(Point(coords),b,c) !=
		    	FUNCTION_SUCCEEDED)
	    	    {
	    	    	screen("ERROR in forward_curve_seg_len_constr(), "
	    	           	"insert_point_in_bond failed\n");
	    	    	clean_up(ERROR);
	    	    }
		    total_offset -= b_len;
		    b = b->next;
		}
	    }
	}
	else
	{
	    if (total_offset < space_tol && bstart != NULL)
	    {
		Coords(bend->end)[0] = Coords(bstart->start)[0];
		Coords(bend->end)[1] = Coords(bstart->start)[1];
		bstart = bstart->prev;
	    	replace_curve_seg_by_bond(c,bstart,bend);
	    }
	}
	return;
}		/*end forward_curve_seg_len_constr*/

static	void backward_curve_seg_len_constr(
	CURVE		*c,
	BOND		*bs,
	BOND		*be,
	int		nbds,
	double		seg_len)
{
	BOND		*b, *bstart, *bend;
	double		b_len, sc_len, offset, total_offset, s, oms;
	double		coords[MAXD];
	double		space_tol;
	boolean		reset_bend;

	b_len = seg_len/(double)nbds;
	space_tol = b_len*1.0e-8;


	offset = b_len;		total_offset = seg_len;
	bstart = bs;		bend = be;
	while (bstart != bend)
	{
	    b = bend;
	    while ((sc_len = bond_length(b)) < offset - space_tol)
	    {
	    	offset -= sc_len;
		if (b == bstart) break;
	    	b = b->prev;
	    }
	    if (sc_len > offset + space_tol)
	    {
	    	s = offset/sc_len;
		interpolate_coords(coords,Coords(b->end),Coords(b->start),s,2);
		if (b == bend) reset_bend = YES;
		else reset_bend = NO;
	    	if (insert_point_in_bond(Point(coords),b,c) !=
		    FUNCTION_SUCCEEDED)
	    	{
	    	    screen("ERROR in backward_curve_seg_len_constr(), "
	    	           "insert_point_in_bond failed\n");
	    	    clean_up(ERROR);
	    	}
		b = b->next;
		if (reset_bend) bend = b;
	    }
	    replace_curve_seg_by_bond(c,b,bend);
	    total_offset -= b_len;
	    bend = b->prev;
	    offset = b_len;
	    if (total_offset < b_len - space_tol) break;
	}
	if (bstart == bend)
	{
	    if (total_offset < space_tol)
	    {
		Coords(bstart->start)[0] = Coords(bstart->end)[0];
		Coords(bstart->start)[1] = Coords(bstart->end)[1];
	    	replace_curve_seg_by_bond(c,bstart,bstart->next);
	    }
	    else
	    {
	    	sc_len = bond_length(bstart);
	    	s = total_offset/sc_len;      oms = 1.0 - s;
	    	Coords(bstart->start)[0] = Coords(bend->end)[0] + s*
			(Coords(bend->start)[0] - Coords(bend->end)[0]);
	    	Coords(bstart->start)[1] = Coords(bend->end)[1] + s*
			(Coords(bend->start)[1] - Coords(bend->end)[1]);
	    	set_bond_length(bstart,2);
		b = bstart;
		while (total_offset/b_len > 2.0 - space_tol)
		{
	    	    s = b_len/total_offset;
		    interpolate_coords(coords,Coords(b->start),Coords(b->end),
					s,2);
	    	    if (insert_point_in_bond(Point(coords),b,c) !=
		    	FUNCTION_SUCCEEDED)
	    	    {
	    	    	screen("ERROR in forward_curve_seg_len_constr(), "
	    	           	"insert_point_in_bond failed\n");
	    	    	clean_up(ERROR);
	    	    }
		    total_offset -= b_len;
		}
	    }
	}
	else
	{
	    if (total_offset < space_tol && bend != NULL)
	    {
		Coords(bstart->start)[0] = Coords(bend->end)[0];
		Coords(bstart->start)[1] = Coords(bend->end)[1];
		bend = bend->next;
	    	replace_curve_seg_by_bond(c,bstart,bend);
	    }
	}
	return;
}		/*end backward_curve_seg_len_constr*/

EXPORT	void FT_SetGlobalIndex(
	Front *front)
{
	set_point_gindex(front);
	set_surface_gindex(front);
	set_curve_gindex(front);
}	/* end FT_SetGlobalIndex */

EXPORT void FT_SetSurfGlobalIndex(
        Front *front)   
{
        set_surface_gindex(front);
}       /* end FT_SetSurfGlobalIndex */

EXPORT void FT_SetCurveGlobalIndex(
        Front *front)   
{
        set_curve_gindex(front);
}       /* end FT_SetCurveGlobalIndex */

EXPORT  CURVE *FT_MakeNodeArrayCurve(
        Front *front,   
        int num_nodes,
        double **node_array,
        COMPONENT   neg_comp,
        COMPONENT   pos_comp,
        boolean is_closed_curve,
	double scale_factor,
        int w_type)
{
	double **dir;
	double *len,*dh,**point_array;
	int i,j,n,*np,num_points;
	CURVE *curve;

	uni_array(&len,num_nodes-1,sizeof(double));
	uni_array(&dh,num_nodes-1,sizeof(double));
	uni_array(&np,num_nodes-1,sizeof(int));
	bi_array(&dir,num_nodes-1,MAXD,sizeof(double));
	num_points = 1;
	for (i = 0; i < num_nodes-1; ++i)
	{
	    dir[i][0] = node_array[i+1][0] - node_array[i][0];
	    dir[i][1] = node_array[i+1][1] - node_array[i][1];
	    len[i] = sqrt(sqr(dir[i][0]) + sqr(dir[i][1]));
	    dir[i][0] /= len[i];
	    dir[i][1] /= len[i];
	    dh[i] = scale_factor*FT_GridSizeInDir(dir[i],front);
	    if (len[i] <= dh[i])
	    {
		np[i] = 1;
		num_points++;
	    }
	    else
	    {
	    	np[i] = irint(len[i]/dh[i]);
	    	dh[i] = len[i]/np[i];
	    	num_points += np[i];
	    }
	}
	bi_array(&point_array,num_points,MAXD,sizeof(double));
	n = 0;
	for (i = 0; i < num_nodes-1; ++i)
	{
	    for (j = 0; j < np[i]; ++j)
	    {
		point_array[n][0] = node_array[i][0] + j*dh[i]*dir[i][0];
		point_array[n][1] = node_array[i][1] + j*dh[i]*dir[i][1];
		n++;
	    }
	}
	point_array[n][0] = node_array[num_nodes-1][0];
	point_array[n][1] = node_array[num_nodes-1][1];
	curve = FT_MakePointArrayCurve(front,num_points,point_array,neg_comp,
			pos_comp,is_closed_curve,w_type);
	free_these(5,point_array,np,len,dir,dh);
	return curve;
}	/* end FT_MakeNodeArrayCurve */

EXPORT	CURVE *FT_MakePointArrayCurve(
	Front *front,
	int num_points,
	double **point_array,
	COMPONENT   neg_comp,
        COMPONENT   pos_comp,
	boolean is_closed_curve,
	int w_type)
{
	CURVE *curve;
	curve = make_array_curve(front->interf,neg_comp,pos_comp,num_points,
			point_array,is_closed_curve);
	wave_type(curve) = w_type;
	if (is_closed_curve)
	    node_type(curve->start) = CLOSED_NODE;
	return curve;
}	/* end FT_MakePointArrayCurve */
		

EXPORT	void FT_MakeEllipticSurf(
	Front *front,
	double *center,
	double *radius,
	COMPONENT   neg_comp,
        COMPONENT   pos_comp,
	int w_type,
	int refinement_level,
	SURFACE **surf)
{
	RECT_GRID *rgr = front->rect_grid;
	RECT_GRID box_rg;
	ELLIP_PARAMS ellip_params;
	int i,dim = rgr->dim;
	double L[MAXD],U[MAXD],dh;
	int gmax[MAXD];
	double *h = rgr->h;

	for (i = 0; i < dim; ++i)
	{
	    ellip_params.cen[i] = center[i];
	    ellip_params.rad[i] = radius[i];
	    L[i] = center[i] - radius[i];
	    U[i] = center[i] + radius[i];
	    gmax[i] = rint((U[i] - L[i])/h[i]);
	    dh = gmax[i]*h[i] - (U[i] - L[i]);
	    L[i] -= sqrt(0.5)*dh;
	    U[i] += sqrt(0.5)*dh;
	    gmax[i] = refinement_level*rint((U[i] - L[i])/h[i]);
	}
	set_box_rect_grid(L,U,gmax,NULL,NULL,dim,&box_rg);
	make_level_surface(&box_rg,front->interf,neg_comp,pos_comp,
			ellipsoid_func,(POINTER)&ellip_params,surf);
	wave_type(*surf) = w_type;
	front->interf->modified = YES;
	interface_reconstructed(front->interf) = NO;
}	/* end FT_MakeEllipticSurf */

EXPORT	void FT_MakeSphericalSurf(
	Front *front,
	double *center,
	double radius,
	COMPONENT   neg_comp,
        COMPONENT   pos_comp,
	int w_type,
	int refinement_level,
	SURFACE **surf)
{
	double radii[MAXD];
	RECT_GRID *rgr = front->rect_grid;
	int i,dim = rgr->dim;

	for (i = 0; i < dim; ++i)
	{
	    radii[i] = radius;
	}
	FT_MakeEllipticSurf(front,center,radii,neg_comp,pos_comp,
			w_type,refinement_level,surf);
}	/* end FT_MakeSphericalSurf */

EXPORT	void FT_MakeDumbBellSurf(
	Front *front,
	double x0,
	double x1,
	double y0,
	double z0,
	double R,
	double r,
	COMPONENT   neg_comp,
        COMPONENT   pos_comp,
	int w_type,
	SURFACE **surf)
{
	RECT_GRID *rgr = front->rect_grid;
	DUMBBELL_PARAMS d_params;
	int i,dim = rgr->dim;
	double coords[MAXD];

	d_params.x0 = x0;
	d_params.x1 = x1;
	d_params.y = y0;
	d_params.z = z0;
	d_params.R = R;
	d_params.rr = r;
	make_level_surface(rgr,front->interf,neg_comp,pos_comp,
			dumbbell_func,(POINTER)&d_params,surf);
	wave_type(*surf) = w_type;
	front->interf->modified = YES;
}	/* end FT_MakeDumbBellSurf */

EXPORT	void FT_MakeProjectileSurf(
	Front *front,
	double *center,
	double R,
	double r,
	double h,
	COMPONENT   neg_comp,
        COMPONENT   pos_comp,
	int w_type,
	SURFACE **surf)
{
	RECT_GRID *rgr = front->rect_grid;
	PROJECTILE_PARAMS proj_params;
	int i,dim = rgr->dim;
	double coords[MAXD];

	proj_params.dim = dim;
	for (i = 0; i < dim; ++i)
	    proj_params.cen[i] = center[i];
	proj_params.R = R;
	proj_params.r = r;
	proj_params.h = h;
	make_level_surface(rgr,front->interf,neg_comp,pos_comp,
			projectile_func,(POINTER)&proj_params,surf);
	wave_type(*surf) = w_type;
	front->interf->modified = YES;
}	/* end FT_MakeProjectileSurf */

EXPORT	void FT_RotateSurface(
	SURFACE *surf,
	double *center,
	double phi,
	double theta)
{
	POINT *p;
	TRI *tri;
	int i;
	boolean first = YES;

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
		if (sorted(p)) continue;
		sorted(p) = YES;
		rotate_point_with_spherical_angle(p,center,phi,theta,first);
		first = NO;
	    }
	}
}	/* end FT_RotateSurface */


EXPORT  void FT_MakeCuboidSurf(
	Front *front,
	double *center,
	double *edge,
	COMPONENT neg_comp,
	COMPONENT pos_comp,
	int w_type,
	int refinement_level,
	SURFACE **surf)
{
	RECT_GRID *rgr = front->rect_grid;
	RECT_GRID box_rg;
	CUBOID_PARAMS cuboid_params;
	int i,dim = rgr->dim;
	double L[MAXD],U[MAXD],dh;
        int gmax[MAXD];
        double *h = rgr->h;
	
	for (i = 0; i < dim; ++i)
	{
	    cuboid_params.center[i] = center[i];
	    cuboid_params.edge[i] = edge[i];
	    L[i] = center[i] - edge[i];
            U[i] = center[i] + edge[i];
            gmax[i] = rint((U[i] - L[i])/h[i]);
            dh = gmax[i]*h[i] - (U[i] - L[i]);
            L[i] -= sqrt(0.5)*dh;
            U[i] += sqrt(0.5)*dh;
            gmax[i] = refinement_level*rint((U[i] - L[i])/h[i]);
	}
	set_box_rect_grid(L,U,gmax,NULL,NULL,dim,&box_rg);
	make_level_surface(&box_rg,front->interf,neg_comp,pos_comp,
			cuboid_func,(POINTER)&cuboid_params,surf);
        wave_type(*surf) = w_type;
        front->interf->modified = YES;
	interface_reconstructed(front->interf) = NO;
}	 /*end FT_MakeCuboidSurf*/

EXPORT  void FT_MakeCylinderSurf(
        Front *front,
        double *center,
        double radius,
	double height,
        COMPONENT neg_comp,
        COMPONENT pos_comp,
        int w_type,
        SURFACE **surf)
{
        RECT_GRID *rgr = front->rect_grid;
        CYLINDER_PARAMS cylinder_params;
        int i,dim = rgr->dim;

        for (i = 0; i < dim; ++i)
        {
            cylinder_params.center[i] = center[i];
            cylinder_params.radius = radius;
	    cylinder_params.height = height;	   
        }
        make_level_surface(rgr,front->interf,neg_comp,pos_comp,
                        cylinder_func,(POINTER)&cylinder_params,surf);
        wave_type(*surf) = w_type;
        front->interf->modified = YES;

}        /*end FT_MakeCylinderSurf*/

EXPORT  void FT_MakeConeSurf(
        Front *front,
        double *center,
        double slope,
	double height,
        COMPONENT neg_comp,
        COMPONENT pos_comp,
        int w_type,
        SURFACE **surf)
{
        RECT_GRID *rgr = front->rect_grid;
        CONE_PARAMS cone_params;
        int i,dim = rgr->dim;

        for (i = 0; i < dim; ++i)
        {
            cone_params.center[i] = center[i];
            cone_params.slope = slope;
	    cone_params.height = height;
        }
        make_level_surface(rgr,front->interf,neg_comp,pos_comp,
                        cone_func,(POINTER)&cone_params,surf);
        wave_type(*surf) = w_type;
        front->interf->modified = YES;

}        /*end FT_MakeConeSurf*/

EXPORT  void FT_MakeTetrahedronSurf(
        Front *front,
        double *center,
        double edge,
        COMPONENT neg_comp,
        COMPONENT pos_comp,
        int w_type,
        SURFACE **surf)
{
        RECT_GRID *rgr = front->rect_grid;
        TETRAHEDRON_PARAMS tetrahedron_params;
        int i,dim = rgr->dim;

        for (i = 0; i < dim; ++i)
        {
            tetrahedron_params.center[i] = center[i];
            tetrahedron_params.edge = edge;
        }
        make_level_surface(rgr,front->interf,neg_comp,pos_comp,
                        tetrahedron_func,(POINTER)&tetrahedron_params,surf);
        wave_type(*surf) = w_type;
        front->interf->modified = YES;

}        /*end FT_MakeTetrahedronSurf*/

EXPORT	void FT_MakePlaneSurf(
	Front *front,
	double *plane_nor,
	double *plane_pt,
	boolean reset_bdry_comp,
	COMPONENT   neg_comp,
        COMPONENT   pos_comp,
	int w_type,
	SURFACE **surf)
{
	RECT_GRID *rgr = front->rect_grid;
	PLANE_PARAMS plane_params;
	int i,dim = rgr->dim;

	for (i = 0; i < dim; ++i)
	{
	    plane_params.N[i] = plane_nor[i];
	    plane_params.P[i] = plane_pt[i];
	}
	make_level_surface(rgr,front->interf,neg_comp,pos_comp,
			plane_func,(POINTER)&plane_params,surf);
	wave_type(*surf) = w_type;
	if (reset_bdry_comp)
	{
	    printf("FT_MakePlaneSurf(): reset_bdry_comp code needed\n");
	    clean_up(ERROR);
	}
	front->interf->modified = YES;
}	/* end FT_MakePlaneSurf */

EXPORT	void FT_InstallSurfEdge(
	SURFACE *surf,
	int hsbdry_type)
{
	install_hsbdry_on_surface(surf,hsbdry_type);
}	/* end FT_InstallSurfEdge */

EXPORT  void FT_CutSurfBdry(
        SURFACE *surf,
        boolean constr_func(POINTER,double*),
        POINTER func_params,
        double **insert_coords,
        int num_pts,
        int insert_idir)
{
        int i;
        for (i = 0; i < num_pts; ++i)
        {
            insert_point_in_surface(insert_idir,insert_coords[i],surf);
        }
        cut_surface(surf,constr_func,func_params,YES);
}       /* end FT_CutSurfBdry */

EXPORT	void FT_MakeEllipticCurve(
	Front *front,
	double *center,
	double *radius,
	COMPONENT   neg_comp,
        COMPONENT   pos_comp,
	int w_type,
	int refinement_level,
	CURVE **curve)
{
	ELLIP2D_PARAMS ellip_params;
	int i,num_seg;
	CURVE **curves;
	RECT_GRID box_rg;
	RECT_GRID *rgr = front->rect_grid;
	int dim = rgr->dim;
        double L[MAXD],U[MAXD],dh;
        int gmax[MAXD];
        double *h = rgr->h;

	ellip_params.x0 = center[0];
	ellip_params.y0 = center[1];
	ellip_params.a = radius[0];
	ellip_params.b = radius[1];
	for (i = 0; i < dim; ++i)
        {
            L[i] = center[i] - radius[i];
            U[i] = center[i] + radius[i];
            gmax[i] = rint((U[i] - L[i])/h[i]);
            dh = gmax[i]*h[i] - (U[i] - L[i]);
            L[i] -= 0.5*dh;
            U[i] += 0.5*dh;
            gmax[i] = refinement_level*rint((U[i] - L[i])/h[i]);
        }
	set_box_rect_grid(L,U,gmax,NULL,NULL,dim,&box_rg);
	curves = make_level_curves(&box_rg,front->interf,neg_comp,
			pos_comp,ellipse_func,(POINTER)&ellip_params,
			YES,&num_seg);
	for (i = 0; i < num_seg; ++i)
	    wave_type(curves[i]) = w_type;
	*curve = curves[0];
}	/* end FT_MakeEllipticCurve */

EXPORT	int FT_Dimension()
{
	return Dimension(current_interface());
}	/* end FT_Dimension */


EXPORT	void FT_XgraphSampleLine(
	char *dirname,
	char *varname,
	boolean data_in_domain,
	int size,
	double *x,
	double *var)
{
	static int alloc_total_size = 0;
	static int alloc_nb_size = 0;
	static double *global_x;
	static double *global_var;
	static double *nb_x;
	static double *nb_var;
	static int *nb_size;
	boolean *data_in_nbdomain;
	int tmp_size,max_nb_size;
	int i,j,n,total_size;
	int xtag = 999;
	char fname[200];
	FILE *xfile;

	if (pp_mynode() != 0)
	{
	    pp_send(xtag,(POINTER)&data_in_domain,sizeof(boolean),0);
	    if (data_in_domain == YES)
	    {
	    	pp_send(xtag,(POINTER)&size,sizeof(int),0);
	    	pp_send(xtag,(POINTER)x,size*sizeof(double),0);
	    	pp_send(xtag,(POINTER)var,size*sizeof(double),0);
	    }
	}
	else
	{
	    uni_array(&nb_size,pp_numnodes(),sizeof(int));
	    uni_array(&data_in_nbdomain,pp_numnodes(),sizeof(boolean));
	    total_size = (data_in_domain == YES) ? size : 0;
	    max_nb_size = 0;
	    for (i = 1; i < pp_numnodes(); ++i)
	    {
		pp_recv(xtag,i,(POINTER)(data_in_nbdomain+i),sizeof(boolean));
		if (data_in_nbdomain[i] == YES)
		{
		    pp_recv(xtag,i,(POINTER)&tmp_size,sizeof(int));
		    total_size += tmp_size;
		    nb_size[i] = tmp_size;
		    if (max_nb_size < tmp_size)
			max_nb_size = tmp_size;
		}
		else
		    nb_size[i] = 0;
	    }
	}
	if (pp_mynode() == 0)
	{
	    if (create_directory(dirname,YES) == FUNCTION_FAILED)
            {
                screen("ERROR directory %s doesn't exist ",
                           "and can't be made\n",dirname);
                clean_up(ERROR);
            }
	    sprintf(fname,"%s/%s",dirname,varname);
	    xfile = fopen(fname,"w");
	    fprintf(xfile,"\n");
            fprintf(xfile,"Next\n");
            fprintf(xfile,"color=red\n");
            fprintf(xfile,"thickness=1.0\n");
	    if (alloc_nb_size < max_nb_size)
	    {
	    	if (nb_x != NULL)
		    free_these(2,nb_x,nb_var);
		uni_array(&nb_x,max_nb_size,sizeof(double));
		uni_array(&nb_var,max_nb_size,sizeof(double));
		alloc_nb_size = max_nb_size;
	    }
	    if (alloc_total_size < total_size)
	    {
	    	if (global_x != NULL)
		    free_these(2,global_x,global_var);
		uni_array(&global_x,total_size,sizeof(double));
		uni_array(&global_var,total_size,sizeof(double));
		alloc_total_size = total_size;
	    }
	    n = 0;
	    if (data_in_domain == YES)
	    {
	    	for (i = 0; i < size; ++i)
	    	{
		    global_x[i]   = x[i];
		    global_var[i] = var[i];
	    	}
		n = size;
	    }
	    for (i = 1; i < pp_numnodes(); ++i)
	    {
		if (data_in_nbdomain[i] == NO) continue;
		pp_recv(xtag,i,(POINTER)nb_x,nb_size[i]*sizeof(double));
		pp_recv(xtag,i,(POINTER)nb_var,nb_size[i]*sizeof(double));
		for (j = 0; j < nb_size[i]; ++j)
		{
		    global_x[n+j] = nb_x[j];
		    global_var[n+j] = nb_var[j];
		}
		n += nb_size[i];
	    }
	    for (i = 0; i < total_size-1; ++i)
	    for (j = i+1; j < total_size; ++j)
	    {
		if (global_x[i] > global_x[j])
		{
		    double tmpx,tmpv;
		    tmpx = global_x[i];
		    global_x[i] = global_x[j];
		    global_x[j] = tmpx;
		    tmpv = global_var[i];
		    global_var[i] = global_var[j];
		    global_var[j] = tmpv;
		}
	    }
	    for (i = 0; i < total_size; ++i)
		fprintf(xfile,"%f %f\n",global_x[i],global_var[i]);
	    fclose(xfile);
	}
	pp_gsync();
}	/* end FT_XgraphSampleLine */

EXPORT void FT_PrintWaveType(
	int w_type)
{
	INTERFACE *intfc = current_interface();
	printf("Wave type = %s\n",wave_type_as_string(w_type,intfc));
}	/* end FT_PrintWaveType */

EXPORT void FT_PrintBoundaryType(
	int dir,
	int side)
{
	const char *bdry_type;
	
	switch (FT_BoundaryType(dir,side))
	{
	case SUBDOMAIN_BOUNDARY:
	    bdry_type = "PERIODIC_BOUNDARY/SUBDOMAIN_BOUNDARY";
	    break;
	case REFLECTION_BOUNDARY:
	    bdry_type = "REFLECTION_BOUNDARY";
	    break;
	case MIXED_TYPE_BOUNDARY:
	    bdry_type = "MIXED_TYPE_BOUNDARY";
	    break;
	case OPEN_BOUNDARY:
	    bdry_type = "OPEN_BOUNDARY";
	    break;
	case FIRST_USER_BOUNDARY_TYPE:
	    bdry_type = "FIRST_USER_BOUNDARY_TYPE";
	    break;
	default:
	    bdry_type = "UNKNOWN_BOUNDARY_TYPE";
	}
	(void) printf("Direction %d side %d: Boundary type = %s\n",
				dir,side,bdry_type);
}	/* end FT_PrintBoundaryType */

EXPORT int FT_BoundaryType(
	int dir,
	int side)
{
	INTERFACE *intfc = current_interface();
	return rect_boundary_type(intfc,dir,side);
}	/* end FT_BoundaryType */

EXPORT	double *FT_GridIntfcTopL(
	Front *front)
{
	if (front->grid_intfc == NULL)
	{
	    (void) printf("In FT_GridIntfcTopL() grid_intfc is NULL!\n");
	    clean_up(ERROR);
	}
	return topological_grid(front->grid_intfc).L;
}	/* end FT_GridIntfcTopL */

EXPORT	double *FT_GridIntfcTopU(
	Front *front)
{
	if (front->grid_intfc == NULL)
	{
	    (void) printf("In FT_GridIntfcTopU() grid_intfc is NULL!\n");
	    clean_up(ERROR);
	}
	return topological_grid(front->grid_intfc).U;
}	/* end FT_GridIntfcTopU */

EXPORT	double *FT_GridIntfcToph(
	Front *front)
{
	if (front->grid_intfc == NULL)
	{
	    (void) printf("In FT_GridIntfcToph() grid_intfc is NULL!\n");
	    clean_up(ERROR);
	}
	return topological_grid(front->grid_intfc).h;
}	/* end FT_GridIntfcToph */

EXPORT	COMPONENT *FT_GridIntfcTopComp(
	Front *front)
{
	Table *T;	
	if (front->grid_intfc == NULL)
	{
	    (void) printf("In FT_GridIntfcTopComp() grid_intfc is NULL!\n");
	    clean_up(ERROR);
	}
	T = table_of_interface(front->grid_intfc);
	return T->components;
}	/* end FT_GridIntfcTopComp */

EXPORT	int *FT_GridIntfcTopGmax(
	Front *front)
{
	if (front->grid_intfc == NULL)
	{
	    (void) printf("In FT_GridIntfcTopGmax() grid_intfc is NULL!\n");
	    clean_up(ERROR);
	}
	return topological_grid(front->grid_intfc).gmax;
}	/* end FT_GridIntfcTopGmax */

EXPORT	RECT_GRID *FT_GridIntfcTopGrid(
	Front *front)
{
	if (front->grid_intfc == NULL)
	{
	    (void) printf("In FT_GridIntfcTopGmax() grid_intfc is NULL!\n");
	    clean_up(ERROR);
	}
	return &topological_grid(front->grid_intfc);
}	/* end FT_GridIntfcTopGrid */


#define		MAX_MOVIE_VARIABLES	20

EXPORT void FT_AddHdfMovieVariable(
	Front *front,
	boolean preset_bound,
	boolean untracked,
	COMPONENT obst_comp,
	const char *var_name,
	int idir,
	double *var_field,
	double (*getStateFunc)(POINTER),
	double max_var,
	double min_var)
{
	static HDF_MOVIE_VAR *hdf_movie_var;
	int i;

        if (debugging("trace"))
            (void) printf("Entering FT_AddHdfMovieVariable()\n");

	if (hdf_movie_var == NULL)
	{
	    FT_ScalarMemoryAlloc((POINTER*)&hdf_movie_var,
				sizeof(HDF_MOVIE_VAR));
	    FT_MatrixMemoryAlloc((POINTER*)&hdf_movie_var->var_name,
				MAX_MOVIE_VARIABLES,100,sizeof(char));
	    FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->top_var,
				MAX_MOVIE_VARIABLES,sizeof(double*));
	    FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->obstacle_comp,
				MAX_MOVIE_VARIABLES,sizeof(COMPONENT));
	    FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->preset_bound,
				MAX_MOVIE_VARIABLES,sizeof(boolean));
	    FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->var_min,
				MAX_MOVIE_VARIABLES,sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->var_max,
				MAX_MOVIE_VARIABLES,sizeof(double));
	    if (FT_Dimension() == 3)
	    	FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->idir,
				MAX_MOVIE_VARIABLES,sizeof(int));
	    hdf_movie_var->num_var = 0;
	    hdf_movie_var->plot_comp = YES;	/* default */
	    hdf_movie_var->plot_bullet = YES;	/* default */
	    front->hdf_movie_var = hdf_movie_var;
	}
	i = front->hdf_movie_var->num_var;
	hdf_movie_var->untracked = untracked;
	hdf_movie_var->get_state_var[i] = getStateFunc;
	hdf_movie_var->top_var[i] = var_field;
	sprintf(hdf_movie_var->var_name[i],"%s",var_name);
	hdf_movie_var->preset_bound[i] = preset_bound;
	hdf_movie_var->obstacle_comp[i] = obst_comp;
	if (preset_bound == YES)
	{
	    hdf_movie_var->var_min[i] = min_var;
	    hdf_movie_var->var_max[i] = max_var;
	}
	if (FT_Dimension() == 3)
	    hdf_movie_var->idir[i] = idir;
	front->hdf_movie_var->num_var += 1;

        if (debugging("trace"))
            (void) printf("Leaving FT_AddHdfMovieVariable()\n");
}	/* end FT_AddHdfMovieVariable */

EXPORT void FT_AddVtkVectorMovieVariable(
	Front *front,
	const char *var_name,
	double **vector_var_field)
{
	static VTK_MOVIE_VAR *vtk_movie_var;
	int i;

        if (debugging("trace"))
            (void) printf("Entering FT_AddVtkVectorMovieVariable()\n");

	if (front->vtk_movie_var == NULL)
	{
	    FT_ScalarMemoryAlloc((POINTER*)&vtk_movie_var,
				sizeof(VTK_MOVIE_VAR));
	    FT_VectorMemoryAlloc((POINTER*)&vtk_movie_var->scalar_var,
				MAX_MOVIE_VARIABLES,sizeof(double*));
	    FT_VectorMemoryAlloc((POINTER*)&vtk_movie_var->vector_var,
				MAX_MOVIE_VARIABLES,sizeof(double**));
	    FT_MatrixMemoryAlloc((POINTER*)&vtk_movie_var->scalar_var_name,
                                MAX_MOVIE_VARIABLES,100,sizeof(char));
	    FT_MatrixMemoryAlloc((POINTER*)&vtk_movie_var->vector_var_name,
                                MAX_MOVIE_VARIABLES,100,sizeof(char));
	    FT_VectorMemoryAlloc((POINTER*)&vtk_movie_var->intfc_var_name,
                                100,sizeof(char));
	    vtk_movie_var->num_scalar_var = 0;
	    vtk_movie_var->num_vector_var = 0;
	    vtk_movie_var->plot_intfc_var = NO;
	    front->vtk_movie_var = vtk_movie_var;
	}
	else
	    vtk_movie_var = front->vtk_movie_var;
	i = front->vtk_movie_var->num_vector_var;
	vtk_movie_var->vector_var[i] = vector_var_field;
	sprintf(vtk_movie_var->vector_var_name[i],"%s",var_name);
	front->vtk_movie_var->num_vector_var += 1;

        if (debugging("trace"))
            (void) printf("Leaving FT_AddVtkVectorMovieVariable()\n");
}	/* end FT_AddVtkVectorMovieVariable */

EXPORT void FT_AddVtkScalarMovieVariable(
	Front *front,
	const char *var_name,
	double *scalar_var_field)
{
	static VTK_MOVIE_VAR *vtk_movie_var;
	int i;

        if (debugging("trace"))
            (void) printf("Entering FT_AddVtkScalarMovieVariable()\n");

	if (front->vtk_movie_var == NULL)
	{
	    FT_ScalarMemoryAlloc((POINTER*)&vtk_movie_var,
				sizeof(VTK_MOVIE_VAR));
	    FT_VectorMemoryAlloc((POINTER*)&vtk_movie_var->scalar_var,
				MAX_MOVIE_VARIABLES,sizeof(double*));
	    FT_VectorMemoryAlloc((POINTER*)&vtk_movie_var->vector_var,
				MAX_MOVIE_VARIABLES,sizeof(double**));
	    FT_MatrixMemoryAlloc((POINTER*)&vtk_movie_var->scalar_var_name,
                                MAX_MOVIE_VARIABLES,100,sizeof(char));
	    FT_MatrixMemoryAlloc((POINTER*)&vtk_movie_var->vector_var_name,
                                MAX_MOVIE_VARIABLES,100,sizeof(char));
	    FT_VectorMemoryAlloc((POINTER*)&vtk_movie_var->intfc_var_name,
                                100,sizeof(char));
	    vtk_movie_var->num_scalar_var = 0;
	    vtk_movie_var->num_vector_var = 0;
	    vtk_movie_var->plot_intfc_var = NO;
	    front->vtk_movie_var = vtk_movie_var;
	}
	else
	    vtk_movie_var = front->vtk_movie_var;
	i = front->vtk_movie_var->num_scalar_var;
	vtk_movie_var->scalar_var[i] = scalar_var_field;
	sprintf(vtk_movie_var->scalar_var_name[i],"%s",var_name);
	front->vtk_movie_var->num_scalar_var += 1;

        if (debugging("trace"))
            (void) printf("Leaving FT_AddVtkScalarMovieVariable()\n");
}	/* end FT_AddVtkScalarMovieVariable */

EXPORT void FT_AddVtkIntfcMovieVariable(
	Front *front,
	const char *var_name)
{
	static VTK_MOVIE_VAR *vtk_movie_var;
	int i;

        if (debugging("trace"))
            (void) printf("Entering FT_AddVtkIntfcMovieVariable()\n");

	if (front->vtk_movie_var == NULL)
	{
	    FT_ScalarMemoryAlloc((POINTER*)&vtk_movie_var,
				sizeof(VTK_MOVIE_VAR));
	    FT_VectorMemoryAlloc((POINTER*)&vtk_movie_var->scalar_var,
				MAX_MOVIE_VARIABLES,sizeof(double*));
	    FT_VectorMemoryAlloc((POINTER*)&vtk_movie_var->vector_var,
				MAX_MOVIE_VARIABLES,sizeof(double**));
	    FT_MatrixMemoryAlloc((POINTER*)&vtk_movie_var->scalar_var_name,
                                MAX_MOVIE_VARIABLES,100,sizeof(char));
	    FT_MatrixMemoryAlloc((POINTER*)&vtk_movie_var->vector_var_name,
                                MAX_MOVIE_VARIABLES,100,sizeof(char));
	    FT_VectorMemoryAlloc((POINTER*)&vtk_movie_var->intfc_var_name,
                                100,sizeof(char));
	    vtk_movie_var->num_scalar_var = 0;
	    vtk_movie_var->num_vector_var = 0;
	    vtk_movie_var->plot_intfc_var = NO;
	    front->vtk_movie_var = vtk_movie_var;
	}
	else
	    vtk_movie_var = front->vtk_movie_var;
	vtk_movie_var->plot_intfc_var = YES;
	sprintf(vtk_movie_var->intfc_var_name,"%s",var_name);

        if (debugging("trace"))
            (void) printf("Leaving FT_AddVtkIntfcMovieVariable()\n");
}	/* end FT_AddVtkIntfcMovieVariable */

EXPORT void FT_ResetDomainAndGrid(
	Front *front,
	double *L,
	double *U,
	int *gmax)
{
	RECT_GRID *gr;
	INTERFACE *intfc = front->interf;
	SURFACE **s;
	CURVE **c;
	const double eps = 10.0*MACH_EPS;

	/* Reset front->rect_grid */
	gr = front->rect_grid;
	set_rect_grid(L,U,L,U,gr->lbuf,gr->ubuf,gmax,gr->dim,&gr->Remap,gr);
	/* Reset interface computational grid */
	gr = computational_grid(front->interf);
	set_rect_grid(L,U,L,U,gr->lbuf,gr->ubuf,gmax,gr->dim,&gr->Remap,gr);
	/* Reset interface topological grid */
	set_topological_grid(front->interf,gr);
	if (Dimension(intfc) == 2)
	{
	    (void) printf("In FT_ResetDomainAndGrid(), code needed!\n");
	    clean_up(ERROR);
	}
	for (s = intfc->surfaces; s && *s; ++s)
	{
	    if (!Boundary(*s))
		continue;
	    for (c = (*s)->pos_curves; c && *c; ++c)
	    {
		if ((*c)->start != (*c)->end)
		{
		    delete_node((*c)->start);
		    delete_node((*c)->end);
		}
		else
		    delete_node((*c)->start);
		delete_curve(*c);
	    }
	    for (c = (*s)->neg_curves; c && *c; ++c)
	    {
		if ((*c)->start != (*c)->end)
		{
		    delete_node((*c)->start);
		    delete_node((*c)->end);
		}
		else
		    delete_node((*c)->start);
		delete_curve(*c);
	    }
	    delete_surface(*s);
	    s = intfc->surfaces;
	}
        (void) set_boundary(front->interf,front->rect_grid,
                        front->interf->default_comp,eps);
        scatter_front(front);
}	/* end FT_ResetDomainAndGrid */

EXPORT	INTERFACE *FT_CollectHypersurfFromSubdomains(
	Front *front,
	int *owner,			/* Destination of collection */
	int w_type)			/* Wave type of surface to collect */
{
	return collect_hyper_surface(front,owner,w_type);
}	/* end FT_CollectHypersurfFromSubdomains */

EXPORT	boolean FT_CoordsInSubdomain(
	Front *front,
	double *coords)			/* Testing coordinate */
{
	RECT_GRID *gr = front->rect_grid;
	int i,dim = gr->dim;
	for (i = 0; i < dim; ++i)
	{
	    if (coords[i] < gr->L[i]) return NO;
	    if (coords[i] >= gr->U[i]) return NO;
	}
	return YES;
}	/* end FT_CoordsInSubdomain */

LOCAL	boolean increment_parametric_function(
	boolean (*func)(POINTER,double,double*),
        POINTER func_params,
	double	ds0,
	double	*dt_next,
	double	t,
	double	*p0,
	double	*p_next,
	int 	dim)
{
	double ds;
	double dt_min,dt_max,dt;
	int    i,imax = 16;
	double ds_min = 0.9*ds0;
	double ds_max = 1.1*ds0;

	dt = *dt_next;
	if (dt <= 0.0)
	{
	    (void) printf("In increment_parametric_function():\n");
	    (void) printf("Initial dt cannot be less or equal to zero\n");
	    clean_up(ERROR);
	}
	(*func)(func_params,t+dt,p_next);
	ds = distance_between_positions(p0,p_next,dim);
	if (ds >= ds_min && ds <= ds_max) return YES;

	if (ds < ds_min)
	{
	    dt_min = t + dt;
	    for (i = 0; i < imax; ++i)
	    {
		dt *= 2.0;
	    	(*func)(func_params,t+dt,p_next);
	    	ds = distance_between_positions(p0,p_next,dim);
		if (ds > ds_min) break;
	    }
	    if (i == imax) return NO;
	    dt_max = t + dt;
	}
	else if (ds > ds_max)
	{
	    dt_max = t + dt;
	    for (i = 0; i < imax; ++i)
	    {
		dt /= 2.0;
	    	(*func)(func_params,t+dt,p_next);
	    	ds = distance_between_positions(p0,p_next,dim);
		if (ds < ds_max) break;
	    }
	    if (i == imax) return NO;
	    dt_min = t + dt;
	}
	if (ds >= ds_min && ds <= ds_max) 
	{
	    *dt_next = dt;
	    return YES;
	}
	for (i = 0; i < imax; ++i)
	{
	    dt = 0.5*(dt_max + dt_min);
	    (*func)(func_params,t+dt,p_next);
	    ds = distance_between_positions(p0,p_next,dim);
	    if (ds > ds_max)
		dt_max = dt;
	    else if (ds < ds_min)
		dt_min = dt;
	    else
		break;
	}
	if (i == imax) return NO;
	else
	{
	    *dt_next = dt;
	    return YES;
	}
}	/* end increment_parametric_function */

EXPORT CURVE *FT_MakeParametricCurve(
	Front	    *front,
	COMPONENT   left_comp,
	COMPONENT   right_comp,
	int	    w_type,
	boolean     (*func)(POINTER,double,double*),
        POINTER     func_params,
	int	    refinement_level,
	boolean	    is_closed)
{
	INTERFACE *intfc = front->interf;
	int     i,dim = FT_Dimension();
	CURVE   *c;
	int	num_points = 0;
	double	**point_array = NULL;
	double	dt,dt_min,dt_max;
	double	t = 0.0;		/* parametric variable */
	int 	max_num_pts = 400;
	double	hmin,ds;
	double	*p0,*p_next;

	hmin = front->rect_grid->h[0];
	for (i = 1; i < dim; ++i)
	    if (hmin > front->rect_grid->h[i])
		hmin = front->rect_grid->h[i];
	ds = hmin/refinement_level; 	/* Average bond length */
	dt = 1.0/400.0;	/* initial increment of parametric variable */

	bi_array(&point_array,max_num_pts,MAXD,sizeof(double));
	/* For all parametric functions, parametric variable must    */
	/* start from 0 and end at 1. If not, needs to be normalized */
	p0 = point_array[0];
	(*func)(func_params,t,point_array[0]);
	num_points++;
	while (t + dt < 1.0)
	{
	    if (num_points >= max_num_pts) /* expand memory */
	    {
		free(point_array);
		max_num_pts += 400;
		bi_array(&point_array,max_num_pts,MAXD,sizeof(double));
	    }
	    p_next = point_array[num_points];
	    if (!increment_parametric_function(func,func_params,ds,
			&dt,t,p0,p_next,dim))
	    {
		(void) printf("ERROR in increment_parametric_function()\n");
		(void) printf("Re-run and turn on parametric_curve\n");
	    }
	    t += dt;
	    p0 = point_array[num_points];
	    num_points++;
	}

	c = make_array_curve(intfc,left_comp,right_comp,num_points,
			point_array,is_closed);

	wave_type(c) = w_type;
	if (point_array != NULL)
	    free(point_array);
	return c;
}	/* end FT_MakeParametricCurve */

EXPORT void FT_PrintTimeStamp(Front *front)
{
	(void) printf("\ntime = %20.14f   step = %5d   next dt = %20.14f\n",
                        front->time,front->step,front->dt);
	fflush(stdout);
}	/* end FT_PrintTimeStamp */

EXPORT  void FT_MakeCrossCylinderSurf(
        Front *front,
        double *center1,
	double *center2,
        double radius1,
	double radius2,
        double height1,
	double height2,
        COMPONENT neg_comp,
        COMPONENT pos_comp,
        int w_type,
        SURFACE **surf)
{
        RECT_GRID *rgr = front->rect_grid;
        CROSSCYLINDER_PARAMS crosscylinder_params;
        int i,dim = rgr->dim;

        for (i = 0; i < dim; ++i)
        {
            crosscylinder_params.center1[i] = center1[i];
	    crosscylinder_params.center2[i] = center2[i];
        }
	
	crosscylinder_params.radius1 = radius1;
	crosscylinder_params.radius2 = radius2;
        crosscylinder_params.height1 = height1;
	crosscylinder_params.height2 = height2;

        make_level_surface(rgr,front->interf,neg_comp,pos_comp,
			crosscylinder_func,
			(POINTER)&crosscylinder_params,surf);
        wave_type(*surf) = w_type;
        front->interf->modified = YES;

}        /*end FT_MakeCrossCylinderSurf*/

EXPORT  void FT_MakeBowlSurf(
        Front *front,
        double *center,
        double radius1,
        double radius2,
	double radius3,
        double height1,
	double height2,
        COMPONENT neg_comp,
        COMPONENT pos_comp,
        int w_type,
        SURFACE **surf)
{
        RECT_GRID *rgr = front->rect_grid;
        BOWL_PARAMS bowl_params;
        int i,dim = rgr->dim;

        for (i = 0; i < dim; ++i)
        {
            bowl_params.center[i] = center[i];
        }
        bowl_params.radius1 = radius1;
        bowl_params.radius2 = radius2;
	bowl_params.radius3 = radius3;
        bowl_params.height1 = height1;
	bowl_params.height2 = height2;

        make_level_surface(rgr,front->interf,neg_comp,pos_comp,
                        bowl_func,(POINTER)&bowl_params,surf);
        wave_type(*surf) = w_type;
        front->interf->modified = YES;

}        /*end FT_MakeBowlSurf*/

EXPORT  void FT_MakeStellatedOctahedronSurf(
        Front *front,
        double *center,
        double edge,
        COMPONENT neg_comp,
        COMPONENT pos_comp,
        int w_type,
        SURFACE **surf)
{
        RECT_GRID *rgr = front->rect_grid;
        STELLATEDOCTAHEDRON_PARAMS stellatedoctahedron_params;
        int i,dim = rgr->dim;

        for (i = 0; i < dim; ++i)
        {
            stellatedoctahedron_params.center[i] = center[i];
        }
	
	stellatedoctahedron_params.edge = edge;

        make_level_surface(rgr,front->interf,neg_comp,pos_comp,
                        stellatedoctahedron_func,
			(POINTER)&stellatedoctahedron_params,
			surf);
        wave_type(*surf) = w_type;
        front->interf->modified = YES;

}        /*end FT_MakeStellatedOctahedronSurf */

EXPORT  void FT_MakePlatformSurf(
	Front *front,
	double *center,
	double radius,
	double height, 
	double slope, 
	COMPONENT neg_comp,
        COMPONENT pos_comp,
        int w_type,
        SURFACE **surf)
{
	RECT_GRID *rgr = front->rect_grid;
        PLATFORM_PARAMS platform_params;
        int i,dim = rgr->dim;

        for (i = 0; i < dim; ++i)
        {
            platform_params.center[i] = center[i];
        }

	platform_params.slope = slope;
        platform_params.height = height;
	platform_params.radius = radius; 

        make_level_surface(rgr,front->interf,neg_comp,pos_comp,
                        platform_func,(POINTER)&platform_params,surf);
        wave_type(*surf) = w_type;
        front->interf->modified = YES;	
}	/* end FT_MakePlatformSurf */

EXPORT void FT_InitSurfVeloFunc(
	SURFACE *surf,
	const char *vfunc_name,
	POINTER vparams,
	int (*vfunc)(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,
			HYPER_SURF*,double*))
{
	VELO_FUNC_PACK *vel_func_pack;
	scalar(&vel_func_pack,sizeof(VELO_FUNC_PACK));
	strcpy(surf->vfunc_name,vfunc_name);
	vel_func_pack->func_params = vparams;
	vel_func_pack->func = vfunc;
	surf->vel_pack = (POINTER)vel_func_pack;
}	/* end FT_InitSurfVeloFunc */

EXPORT void FT_InitCurveVeloFunc(
	CURVE *curve,
	const char *vfunc_name,
	POINTER vparams,
	int (*vfunc)(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,
			HYPER_SURF*,double*))
{
	VELO_FUNC_PACK *vel_func_pack;
	strcpy(curve->vfunc_name,vfunc_name);
	scalar(&vel_func_pack,sizeof(VELO_FUNC_PACK));
	vel_func_pack->func_params = vparams;
	vel_func_pack->func = vfunc;
	curve->vel_pack = (POINTER)vel_func_pack;
}	/* end FT_InitCurveVeloFunc */

EXPORT void FT_InitNodeVeloFunc(
	NODE *node,
	const char *vfunc_name,
	POINTER vparams,
	int (*vfunc)(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,
			HYPER_SURF*,double*))
{
	VELO_FUNC_PACK *vel_func_pack;
	scalar(&vel_func_pack,sizeof(VELO_FUNC_PACK));
	strcpy(node->vfunc_name,vfunc_name);
	vel_func_pack->func_params = vparams;
	vel_func_pack->func = vfunc;
	node->vel_pack = (POINTER)vel_func_pack;
}	/* end FT_InitNodeVeloFunc */


EXPORT void FT_InitFrontVeloFunc(
	Front *front,
	VELO_FUNC_PACK *velo_func_pack)
{
	char s[100];
	int dim = front->rect_grid->dim;

	/* Initialize front velocity field */

	if (debugging("trace")) 
	    (void) printf("Entering FT_InitFrontVeloFunc()\n");
	if (velo_func_pack == NULL)
	{
	    screen("\n\t\tSpecifying Velocity Field\n\n");
	    screen("Supported velocity fields are \n"
	       "\tTranslation (t)\n"
	       "\tRadial motion (r)\n"
	       "\tShear motion (s)\n"
	       "\tSinusiodal motion (w)\n"
	       "\tCircular rotation (c)\n"
	       "\tNormal motion (n)\n"
	       "\tCurvature dependent motion (k)\n"
	       "\tFlame motion (f)\n"
	       "\tBurgers equation solver (b)\n"
	       "\tBi-Polar velocity (h)\n"
	       "\tVortex velocity (v)\n"
	       "\tDouble vortex velocity (d)\n"
	       "Enter choice: ");
	    (void) Scanf("%s\n",s);
	    switch (s[0])
	    {
	    case 'T':
	    case 't':
	    	front->vfunc = translation_vel;
	    	init_translation_params(&front->vparams,dim);
	    	front->_point_propagate = fourth_order_point_propagate;
	    	break;
	    case 'R':
	    case 'r':
	    	front->vfunc = radial_motion_vel;
	    	init_radial_motion_params(&front->vparams,dim);
	    	front->_point_propagate = fourth_order_point_propagate;
	    	break;
	    case 'S':
	    case 's':
	    	front->vfunc = shear_motion_vel;
	    	init_shear_motion_params(&front->vparams,dim);
	    	front->_point_propagate = fourth_order_point_propagate;
	    	break;
	    case 'W':
	    case 'w':
	    	front->vfunc = sine_motion_vel;
	    	init_sine_motion_params(&front->vparams,dim);
	    	front->_point_propagate = fourth_order_point_propagate;
	    break;
	    case 'C':
	    case 'c':
	    	front->vfunc = circular_rotation_vel;
	    	init_circular_rotation_params(&front->vparams,dim);
	    	front->_point_propagate = fourth_order_point_propagate;
	    	break;
	    case 'B':
	    case 'b':
	    	front->vfunc = burgers_vel;
	    	init_burgers_params(&front->vparams,dim);
	    	front->_point_propagate = fourth_order_point_propagate;
	    break;
	    case 'H':
	    case 'h':
	    	front->vfunc = bipolar_vel;
	    	init_bipolar_params(&front->vparams,dim);
	    	front->_point_propagate = fourth_order_point_propagate;
	    	break;
	    case 'N':
	    case 'n':
	    	front->vfunc = normal_vel;
	    	init_norv_params(&front->vparams,dim);
	    	front->_point_propagate = first_order_point_propagate;
	    	break;
	    case 'K':
	    case 'k':
	    	front->vfunc = curvature_vel;
	    	init_curvature_params(&front->vparams,dim);
	    	front->_point_propagate = first_order_point_propagate;
	    	break;
	    case 'F':
	    case 'f':
	    	front->vfunc = flame_vel;
	    	init_flame_params(&front->vparams,dim);
	    	front->_point_propagate = first_order_point_propagate;
	    	break;
	    case 'V':
	    case 'v':
	    	front->vfunc = vortex_vel;
	    	init_vortex_params(&front->vparams,dim);
	    	front->_point_propagate = fourth_order_point_propagate;
	    	break;
	    case 'D':
	    case 'd':
	    	front->vfunc = double_vortex_vel;
	    	init_bipolar_params(&front->vparams,dim);
	    	front->_point_propagate = fourth_order_point_propagate;
	    	break;
	    default:
	    	screen("Unsupported velocity field for first "	
		       "order point propagation\n");
	    	clean_up(ERROR);
	    }
	}
	else
	{
	    front->vfunc = velo_func_pack->func;
	    front->vparams = velo_func_pack->func_params;
	    if (velo_func_pack->point_propagate != NULL)
	    	front->_point_propagate = velo_func_pack->point_propagate;
	    else
	    	front->_point_propagate = first_order_point_propagate;
	}
	if (debugging("trace")) 
	    (void) printf("Leaving FT_InitFrontVeloFunc()\n");
}	/* end FT_InitFrontVeloFunc_function */

EXPORT boolean FT_CheckSurfCompConsistency(
	Front *front,
	SURFACE *surf)
{
	INTERFACE *intfc = front->interf;
	INTERFACE *test_intfc = copy_interface(intfc);
	SURFACE *test_surf;
	SURFACE **ss,**ts;
	TRI *tri;
	POINT *p;
	int outer_comp = positive_component(surf);
	int test_comp;
	int i;

	for (ss = intfc->surfaces, ts = test_intfc->surfaces; 
	     ss && *ss; ++ss, ++ts)
	{
	    if (*ss == surf)
	    {
		test_surf = *ts;
		break;
	    }
	}
	delete_surface(test_surf);
	reset_surface_points(surf);
	surf_tri_loop(surf,tri)
	{
	    for (i = 0; i < 3; ++i)
	    {
		p = Point_of_tri(tri)[i];
		if (sorted(p)) continue;
		test_comp = component(Coords(p),test_intfc);
		if (outer_comp != test_comp)
		{
		    (void) printf("Inconsistent component at (%f %f %f)\n",
				Coords(p)[0],Coords(p)[1],Coords(p)[2]);
		    (void) printf("Surface outer component = %d\n",outer_comp);
		    (void) printf("Ambient component = %d\n",test_comp);
		    delete_interface(test_intfc);
		    return NO;
		}
		sorted(p) = YES;
	    }
	}
	delete_interface(test_intfc);
	return YES;
}	/* end FT_CheckSurfCompConsistency */

EXPORT boolean FT_NextTopGridIcoordsInDir(
	Front *front,
	int *icoords,
	GRID_DIRECTION dir,
	int *icoords_next)
{
	RECT_GRID *rgr = &topological_grid(front->grid_intfc);
	int gmin[3] = {0,0,0};
	int *gmax = rgr->gmax;
	boolean status;
	status = next_ip_in_dir(icoords,dir,icoords_next,gmin,gmax);
	return status;
}	/* end FT_NextTopGridIcoordsInDir */

EXPORT boolean FT_AdjTopGridIcoords(
	Front *front,
	int *icoords,
	GRID_DIRECTION *dir,
	int *icoords_next)
{
	RECT_GRID *rgr = &topological_grid(front->grid_intfc);
	int gmin[3] = {0,0,0};
	int *gmax = rgr->gmax;
	boolean status;
	int ic_tmp[MAXD];
	int i,dim = rgr->dim;
	for (i = 0; i < dim; ++i)
	    ic_tmp[i] = icoords[i];
	status = YES;
	for (i = 0; i < dim; ++i)
	{
	    if (dir[i] != CENTER)
	    	status = next_ip_in_dir(ic_tmp,dir[i],ic_tmp,gmin,gmax);
	    if (!status) return NO;
	}
	for (i = 0; i < dim; ++i)
	    icoords_next[i] = ic_tmp[i];
	return YES;
}	/* end FT_AdjTopGridIcoords */

EXPORT boolean FT_NearestGridIcoordsWithComp(
	Front *front,
	int *icoords,
	int comp,
	int *comp_map,
	int *gmax,
	int *icoords_next)
{
	boolean status;
	status = find_nearest_ring2_cell_with_comp(icoords,icoords_next,
                      gmax,comp_map,comp);
	return status;
}	/* end FT_NearestGridIcoordsWithComp */

struct _CELL_INFO_2D {
	/* Input vatiables */
	boolean is_corner;
	int comp[2][2];
	int icrds[2][2][MAXD];
	double crds[2][2][MAXD];
	int soln_comp;
	int cell_comp;
	double crx_coords[2][2][MAXD];
	int nv;
	/* output variables */
	int side_flag;		/* used in corner case association */
	boolean nb_flag[3][3];	/* Yes for side of legal transfer */
	double nb_frac[3][3];
	double orphan;
	double full_cell_vol;
};
typedef struct _CELL_INFO_2D CELL_INFO_2D;

struct _CELL_INFO_1D {
	/* Input vatiables */
	int comp[2];
	int icrds[2][MAXD];
	double crds[2][MAXD];
	int soln_comp;
	int cell_comp;
	double crx_coords[MAXD];
	int nv;
	/* output variables */
	double nb_frac[3];
	boolean nb_flag[3];
	double full_cell_vol;
};
typedef struct _CELL_INFO_1D CELL_INFO_1D;

static void cell_area(CELL_INFO_2D*);
static void cell_length(CELL_INFO_1D*);
static void rotate_cell(CELL_INFO_2D*,boolean);
static void reverse_cell(CELL_INFO_2D*);

EXPORT void FT_ComputeVolumeFraction(
	Front *front,
	int num_phases,
	COMPONENT *comps,
	CELL_PART *cell_part)
{
	INTERFACE *grid_intfc,*comp_grid_intfc;
	RECT_GRID *top_grid,*comp_grid;
	int *top_comp,*comp_comp;
	int *top_gmax,*comp_gmax;
	double *top_L,*top_h,*comp_L,*comp_U,*comp_h;
	Table *T,*compT;
	int *lbuf,*ubuf;
	int imin[MAXD],imax[MAXD];
	int dim;

	int i,j,k,ii,jj,kk,l,m,n,nr,ns;
	int ic,icn,icc,icc1,icc2;
	HYPER_SURF *hs;
	static CELL_INFO_2D cell2d;
	static CELL_INFO_1D cell1d;
	int lm[MAXD],um[MAXD]; 	/* lower and upper margin */
	int icoords[MAXD],ipn[MAXD];
	static boolean first = YES;
	const GRID_DIRECTION dir[3][3] = {{WEST,CENTER,EAST},
					  {SOUTH,CENTER,NORTH},
					  {LOWER,CENTER,UPPER}};
	double S,Sn,new_S,v_old,v_oldn;
	double full_cell_vol;
	GRID_DIRECTION dirs[MAXD];
	int size;

	if (debugging("trace"))
	    (void) printf("Entering computeVolumeFraction()\n");

	if (num_phases > 2)
	{
	    (void) printf("Too many phases %d, code needed\n",num_phases);
	    (void) printf("Current code can only handle two phases!\n");
	    clean_up(ERROR);
	}
	/* Set dual grid parameters */
	grid_intfc = front->grid_intfc;
        top_grid = &topological_grid(grid_intfc);
	top_gmax = top_grid->gmax;
        lbuf = front->rect_grid->lbuf;
        ubuf = front->rect_grid->ubuf;
        top_L = top_grid->L;
        top_h = top_grid->h;
        dim = grid_intfc->dim;
        T = table_of_interface(grid_intfc);
        top_comp = T->components;

	/* Set computational grid parameters */
	comp_grid_intfc = front->comp_grid_intfc;
	comp_comp = table_of_interface(comp_grid_intfc)->components;
	comp_grid = &topological_grid(comp_grid_intfc);
	comp_L = comp_grid->L;
	comp_U = comp_grid->U;
	comp_h = comp_grid->h;
	comp_gmax = comp_grid->gmax;
	size = 1;
	for (i = 0; i < dim; ++i)
	{
	    lm[i] = (lbuf[i] == 0) ? 0 : 1;
	    um[i] = (ubuf[i] == 0) ? 0 : 1;
            imin[i] = (lbuf[i] == 0) ? 1 : lbuf[i];
            imax[i] = (ubuf[i] == 0) ? top_gmax[i] - 1 : top_gmax[i] - ubuf[i];
	    size *= (top_gmax[i] + 1);
	}
	for (ic = 0; ic < size; ++ic)
	{
	    cell_part[ic].nr_new = 0;
	    cell_part[ic].ns_new = 0;
	    for (i = 0; i < 2; ++i)
	    {
	    	cell_part[ic].vol_old[i] = cell_part[ic].vol_new[i];
	    	cell_part[ic].vol_new[i] = 0.0;
	    }
	}

      for (n = 0; n < num_phases; ++n)
      {
	switch (dim)
	{
	case 1:
	    cell1d.soln_comp = comps[n];
	    cell1d.full_cell_vol = full_cell_vol = comp_h[0];
	    for (i = imin[0] - lm[0]; i <= imax[0] + um[0]; ++i)
	    {
		icoords[0] = i;
		ic = d_index(icoords,top_gmax,dim);
		cell1d.cell_comp = top_comp[ic];
		cell1d.nv = 0;
	    	for (l = 0; l < 3; ++l)
		{
		    cell1d.nb_flag[l] = NO;
		    cell1d.nb_frac[l] = 0.0;
		}
	    	for (ii = 0; ii < 2; ++ii)
		{
		    cell1d.icrds[ii][0] = i + ii - 1;
		    cell1d.crds[ii][0] = comp_L[0] + (i + ii - 1)*comp_h[0];
		    icc = d_index(cell1d.icrds[ii],comp_gmax,dim);
		    cell1d.comp[ii] = comp_comp[icc];
		    if (comp_comp[icc] == cell1d.soln_comp) 
		    {
		    	cell1d.nv++;
			for (l = ii; l < ii+2; ++l)
			    cell1d.nb_flag[l] = YES;
		    }
		}
		if (cell1d.nv == 0 || cell1d.nv == 2)
		{
		    if (top_comp[ic] == comps[n])
		    {
			cell_part[ic].vol_new[n] += full_cell_vol;
		    }
		    continue;
		}
		icc1 = d_index(cell1d.icrds[0],comp_gmax,dim);
		icc2 = d_index(cell1d.icrds[1],comp_gmax,dim);
		cell1d.crx_coords[0] = HUGE;
		if (comp_comp[icc1] != comp_comp[icc2])
		{
		    boolean status;
		    icc1 = d_index(cell1d.icrds[0],comp_gmax,dim);
		    status = FT_CoordsAtGridCrossing(front,comp_grid_intfc,
		    		cell1d.icrds[0],EAST,comp_comp[icc1],
				&hs,cell1d.crx_coords);
		}
	    	for (l = 0; l < 3; ++l)
		{
		    if (cell1d.nb_flag[l] == NO) continue;
		    dirs[0] = dir[0][l];
		    if (!FT_AdjTopGridIcoords(front,icoords,dirs,ipn))
		    	cell1d.nb_flag[l] = NO;
		    icn = d_index(ipn,top_gmax,dim);
		    if (top_comp[icn] != comps[n])
		    	cell1d.nb_flag[l] = NO;
		}
		cell_length(&cell1d);
		for (l = 0; l < 3; ++l)
		{
		    if (cell1d.nb_frac[l] == 0.0) continue;
		    dirs[0] = dir[0][l];
		    FT_AdjTopGridIcoords(front,icoords,dirs,ipn);
		    icn = d_index(ipn,top_gmax,dim);
                    cell_part[icn].vol_new[n] += cell1d.nb_frac[l];
		    if (ic != icn)
		    {
		    	nr = cell_part[icn].nr_new;
                    	cell_part[icn].icn_new[nr] = ic;
		    	cell_part[icn].nr_new++;
		    	ns = cell_part[ic].ns_new;
                    	cell_part[ic].icn_new_send[ns] = icn;
		    	cell_part[ic].ns_new++;
		    }
		}
	    }
	    break;
	case 2:
	    cell2d.soln_comp = comps[n];
	    cell2d.full_cell_vol = full_cell_vol = comp_h[0]*comp_h[1];
	    for (i = imin[0] - lm[0]; i <= imax[0] + um[0]; ++i)
	    for (j = imin[1] - lm[1]; j <= imax[1] + um[1]; ++j)
	    {
		icoords[0] = i;
		icoords[1] = j;
		ic = d_index(icoords,top_gmax,dim);
		cell2d.cell_comp = top_comp[ic];
		cell2d.nv = 0;
		cell2d.orphan = 0.0;
		cell2d.is_corner = NO;
	    	for (l = 0; l < 3; ++l)
	    	for (m = 0; m < 3; ++m)
		{
		    cell2d.nb_flag[l][m] = NO;
		    cell2d.nb_frac[l][m] = 0.0;
		}
	    	for (ii = 0; ii < 2; ++ii)
	    	for (jj = 0; jj < 2; ++jj)
		{
		    cell2d.icrds[ii][jj][0] = i + ii - 1;
		    cell2d.icrds[ii][jj][1] = j + jj - 1;
		    cell2d.crds[ii][jj][0] = comp_L[0] + (i + ii - 1)*comp_h[0];
		    cell2d.crds[ii][jj][1] = comp_L[1] + (j + jj - 1)*comp_h[1];
		    icc = d_index(cell2d.icrds[ii][jj],comp_gmax,dim);
		    cell2d.comp[ii][jj] = comp_comp[icc];
		    if (comp_comp[icc] == cell2d.soln_comp) 
		    {
		    	cell2d.nv++;
			for (l = ii; l < ii+2; ++l)
			for (m = jj; m < jj+2; ++m)
			    cell2d.nb_flag[l][m] = YES;
		    }
		}
		if (cell2d.nv == 0 || cell2d.nv == 4)
		{
		    if (top_comp[ic] == comps[n])
		    {
			cell_part[ic].vol_new[n]  += full_cell_vol;
		    }
		    continue;
		}
	    	for (ii = 0; ii < 2; ++ii)
		{
		    icc1 = d_index(cell2d.icrds[0][ii],comp_gmax,dim);
		    icc2 = d_index(cell2d.icrds[1][ii],comp_gmax,dim);
		    cell2d.crx_coords[0][ii][0] = cell2d.crx_coords[0][ii][1] 
		    			= HUGE;
		    if (comp_comp[icc1] != comp_comp[icc2])
		    {
		    	FT_CoordsAtGridCrossing(front,comp_grid_intfc,
		    		cell2d.icrds[0][ii],EAST,comp_comp[icc1],
				&hs,cell2d.crx_coords[0][ii]);
		    }
		    icc1 = d_index(cell2d.icrds[ii][0],comp_gmax,dim);
		    icc2 = d_index(cell2d.icrds[ii][1],comp_gmax,dim);
		    cell2d.crx_coords[1][ii][0] = cell2d.crx_coords[1][ii][1] 
		    			= HUGE;
		    if (comp_comp[icc1] != comp_comp[icc2])
		    {
		    	FT_CoordsAtGridCrossing(front,comp_grid_intfc,
		    		cell2d.icrds[ii][0],NORTH,comp_comp[icc1],
				&hs,cell2d.crx_coords[1][ii]);
		    }
		}
	    	for (l = 0; l < 3; ++l)
	    	for (m = 0; m < 3; ++m)
		{
		    if (cell2d.nb_flag[l][m] == NO) continue;
		    dirs[0] = dir[0][l];
		    dirs[1] = dir[1][m];
		    if (!FT_AdjTopGridIcoords(front,icoords,dirs,ipn))
		    	cell2d.nb_flag[l][m] = NO;
		    icn = d_index(ipn,top_gmax,dim);
		    if (top_comp[icn] != comps[n])
		    	cell2d.nb_flag[l][m] = NO;
		}
		if ((i == imin[0] || i == imax[0]) && 
		    (j == imin[1] || j == imax[1]))
		    cell2d.is_corner = YES;
		cell_area(&cell2d);
		for (l = 0; l < 3; ++l)
		for (m = 0; m < 3; ++m)
		{
		    if (cell2d.nb_frac[l][m] == 0.0) continue;
		    dirs[0] = dir[0][l];
		    dirs[1] = dir[1][m];
		    FT_AdjTopGridIcoords(front,icoords,dirs,ipn);
		    icn = d_index(ipn,top_gmax,dim);
		    cell_part[icn].vol_new[n] += cell2d.nb_frac[l][m];
		    if (icn != ic)
		    {
		    	nr = cell_part[icn].nr_new;
		    	cell_part[icn].icn_new[nr] = ic;
		    	cell_part[icn].nr_new++;
		    	ns = cell_part[ic].ns_new;
                    	cell_part[ic].icn_new_send[ns] = icn;
		    	cell_part[ic].ns_new++;
		    }
		}
		if (cell2d.orphan != 0.0)
		{
		    if (!find_nearest_ring2_cell_with_comp(icoords,ipn,
		    		top_gmax,top_comp,comps[n]))
		    {
		    	printf("Cannot find parent cell in ring 2!\n");
			clean_up(ERROR);
		    }
		    icn = d_index(ipn,top_gmax,dim);
		    cell_part[icn].vol_new[n] += cell2d.orphan;
		    nr = cell_part[icn].nr_new;
		    cell_part[ic].icn_new[nr] = icn;
		    cell_part[icn].nr_new++;
		    ns = cell_part[ic].ns_new;
                    cell_part[ic].icn_new_send[ns] = icn;
		    cell_part[ic].ns_new++;
		}
	    }
	    break;
	case 3:
	    for (i = imin[0] - lm[0]; i <= imax[0] + um[0]; ++i)
	    for (j = imin[1] - lm[1]; j <= imax[1] + um[1]; ++j)
	    for (k = imin[2] - lm[2]; k <= imax[2] + um[2]; ++k)
	    {
	    	(void) printf("In computeVolumeFraction(), 3D code needed!\n");
		clean_up(ERROR);
	    }
	    break;
	default:
	    (void) printf("Unknown dimension\n");
	    clean_up(ERROR);
	}
      }

	/*
	FT_ParallelExchGridStructArrayBuffer((POINTER)cell_part,front,
				sizeof(CELL_PART));
	*/
	if (first)
	{
	    for (i = 0; i < size; ++i)
	    {
		cell_part[i].vol_old[0] = cell_part[i].vol_new[0];
		cell_part[i].vol_old[1] = cell_part[i].vol_new[1];
	    }
	}
	first = NO;
	if (debugging("trace"))
	    (void) printf("Leaving computeVolumeFraction()\n");
}	/* end computeVolumeFraction */

static boolean find_nearest_ring2_cell_with_comp(
	int *icoords,
	int *ipn,
	int *top_gmax,
	int *top_comp,
	int soln_comp)
{
	int i,j,ic;
	int dist,min_dist = 100;
	int ip[MAXD];
	boolean status = NO;

	for (i = -8; i <= 8; ++i)
	for (j = -8; j <= 8; ++j)
	{
	    ip[0] = icoords[0] + i;
	    ip[1] = icoords[1] + j;
	    if (ip[0] < 0 || ip[0] > top_gmax[0]) continue;
	    if (ip[1] < 0 || ip[1] > top_gmax[1]) continue;
	    ic = d_index(ip,top_gmax,2);
	    if (top_comp[ic] != soln_comp) continue;
	    status = YES;
	    dist = i*i + j*j;
	    if (dist < min_dist)
	    {
	    	min_dist = dist;
		ipn[0] = ip[0];
		ipn[1] = ip[1];
	    }
	}
	return status;
}	/* end find_nearest_ring2_cell_with_comp */

static double nv1_area(CELL_INFO_2D *cell)
{
	double a,h;
	int i;
	a = h = 0.0;
	for (i = 0; i < 2; ++i)
	{
	    a += fabs(cell->crx_coords[0][0][i] - cell->crds[0][0][i]);
	    h += fabs(cell->crx_coords[1][0][i] - cell->crds[0][0][i]);
	    if (cell->crx_coords[0][0][i] == HUGE ||
	    	cell->crx_coords[1][0][i] == HUGE)
	    {
	    	printf("ERROR in nv1_area(), crx_coords unassigned!\n");
		clean_up(ERROR);
	    }
	}
	if (a > h) cell->side_flag = 0;
	else cell->side_flag = 1;
	return 0.5*a*h;
}	/* end nv1_area */

static double nv2_area(CELL_INFO_2D *cell)
{
	double a,b,h;
	int i;
	a = b = h = 0.0;
	for (i = 0; i < 2; ++i)
	{
	    a += fabs(cell->crx_coords[1][0][i] - cell->crds[0][0][i]);
	    b += fabs(cell->crx_coords[1][1][i] - cell->crds[1][0][i]);
	    h += fabs(cell->crds[1][0][i] - cell->crds[0][0][i]);
	    if (cell->crx_coords[1][0][i] == HUGE ||
	    	cell->crx_coords[1][1][i] == HUGE)
	    {
	    	printf("ERROR in nv2_area(), crx_coords unassigned!\n");
		clean_up(ERROR);
	    }
	}
	if (a > b) cell->side_flag = 0;
	else cell->side_flag = 1;
	return 0.5*(a + b)*h;
}	/* end nv2_area */

static double nv3_area(CELL_INFO_2D *cell)
{
	double a,h,l1,l2;
	int i;
	a = h = l1 = l2 = 0.0;
	for (i = 0; i < 2; ++i)
	{
	    a += fabs(cell->crx_coords[0][0][i] - cell->crds[0][0][i]);
	    h += fabs(cell->crx_coords[1][0][i] - cell->crds[0][0][i]);
	    l1 += fabs(cell->crds[1][0][i] - cell->crds[0][0][i]);
	    l2 += fabs(cell->crds[0][1][i] - cell->crds[0][0][i]);
	    if (cell->crx_coords[0][0][i] == HUGE ||
	    	cell->crx_coords[1][0][i] == HUGE)
	    {
	    	printf("ERROR in nv3_area(), crx_coords unassigned!\n");
		clean_up(ERROR);
	    }
	}
	if (a > h) cell->side_flag = 0;
	else cell->side_flag = 1;
	return l1*l2 - 0.5*a*h;
}	/* end nv3_area */

static void cell_area(
	CELL_INFO_2D *cell)
{
	double area;
	double a,b,h;
	int i,j;
	boolean flag[3][3];
	double frac[3][3];
	boolean first = YES;
	boolean rotate_vfrac = NO;

	for (i = 0; i < 2; ++i)
	for (j = 0; j < 2; ++j)
	    cell->nb_frac[i][j] = 0.0;
	cell->orphan = 0.0;

	if (cell->nv == 1)
	{
	    if (cell->is_corner)
	    {
	    	cell->nb_frac[1][1] = (cell->cell_comp == cell->soln_comp) ?
			cell->full_cell_vol : 0.0;
	    	return;
	    }
	    for (i = 0; i < 4; ++i)
	    {
	    	if (cell->comp[0][0] == cell->soln_comp)
		{
	    	    area = nv1_area(cell);
		    if (cell->cell_comp == cell->soln_comp)
		    	area = 0.5*cell->full_cell_vol;
		    if (cell->cell_comp == cell->soln_comp)
		    	cell->nb_frac[1][1] = area;
		    else if ((cell->side_flag == 0 || !cell->nb_flag[0][1])
			&& cell->nb_flag[1][0])
		    	cell->nb_frac[1][0] = area;
		    else if (cell->nb_flag[0][1])
		    	cell->nb_frac[0][1] = area;
		    else if (cell->nb_flag[0][0])
		    	cell->nb_frac[0][0] = area;
		    else
		    	cell->orphan = area;
		    if (cell->nb_frac[1][1] != 0.0 || cell->orphan != 0.0)
		    	break;
		    for (j = 0; j < i; ++j)
		    	reverse_cell(cell);
		    break;
		}
		rotate_cell(cell,NO);
	    }
	}
	else if (cell->nv == 2)
	{
	    double area_sum = 0.0;
	    for (i = 0; i < 4; ++i)
	    {
	    	if (cell->comp[0][0] == cell->soln_comp &&
	    	    cell->comp[1][0] == cell->soln_comp)
	    	{
	    	    area = nv2_area(cell);
		    if (cell->cell_comp == cell->soln_comp &&
		    	area < 0.5*cell->full_cell_vol)
		    	area = 0.5*cell->full_cell_vol;
		    if (cell->cell_comp != cell->soln_comp &&
		    	area > 0.5*cell->full_cell_vol)
		    	area = 0.5*cell->full_cell_vol;

		    if (cell->cell_comp == cell->soln_comp)
		    	cell->nb_frac[1][1] = area;
		    else if (cell->nb_flag[1][0])
		    	cell->nb_frac[1][0] = area;
		    else if ((cell->side_flag == 0 || !cell->nb_flag[2][1]) &&
		    	cell->nb_flag[0][1])
		    	cell->nb_frac[0][1] = area;
		    else if (cell->nb_flag[2][1])
		    	cell->nb_frac[2][1] = area;
		    else if (cell->nb_flag[0][0])
		    	cell->nb_frac[0][0] = area;
		    else if (cell->nb_flag[2][0])
		    	cell->nb_frac[2][0] = area;
		    else
    		    	cell->orphan = area;
		    if (cell->nb_frac[1][1] != 0.0 || cell->orphan != 0.0)
		    	break;
		    for (j = 0; j < i; ++j)
		    	reverse_cell(cell);
		    break;
	    	}
	    	else if (cell->comp[0][0] == cell->soln_comp &&
	    	         cell->comp[1][1] == cell->soln_comp && 
		        !cell->nb_flag[1][1])
	    	{
		    rotate_vfrac = YES;
	    	    area = nv1_area(cell);
		    area_sum += area;
		    if ((cell->side_flag == 0 || !cell->nb_flag[0][1])
			&& cell->nb_flag[1][0])
		    	cell->nb_frac[1][0] = area;
		    else if (cell->nb_flag[0][1])
		    	cell->nb_frac[0][1] = area;
		    else if (cell->nb_flag[0][0])
		    	cell->nb_frac[0][0] = area;
		    else
		    	cell->orphan = area;
		    if (first) 
		    	first = NO;
		    else 
		    {
		    	for (j = 0; j < i; ++j)
		    	    reverse_cell(cell);
		    	break;
		    }
		}
	    	else if (cell->comp[0][0] != cell->soln_comp &&
	    	         cell->comp[1][1] != cell->soln_comp && 
		         cell->nb_flag[1][1])
	    	{
		    rotate_vfrac = YES;
		    if (first) 
		    {
		    	first = NO;
			cell->nb_frac[1][1] = nv3_area(cell);
		    }
		    else 
		    {
			cell->nb_frac[1][1] -= nv1_area(cell);
		    	area_sum = cell->nb_frac[1][1];
		    	break;
		    }
		}
	    	rotate_cell(cell,rotate_vfrac);
	    }
	    if (rotate_vfrac)
	    {
		if (debugging("DD"))
		{
		    printf("cell->cell_comp = %d cell->soln_comp = %d\n",
			cell->cell_comp,cell->soln_comp);
		    printf("area_sum = %20.14f\n",area_sum);
		    printf("half_vol = %20.14f\n",0.5*cell->full_cell_vol);
		}
	    	if ((cell->cell_comp != cell->soln_comp &&
		     area_sum > 0.5*cell->full_cell_vol) ||
		    (cell->cell_comp == cell->soln_comp &&
		     area_sum < 0.5*cell->full_cell_vol))
		{
		    double lambda = 0.5*cell->full_cell_vol/area_sum;
		    if (debugging("DD"))
		    	printf("Case caught: half cell vol = %20.14f\n",
		    		0.5*cell->full_cell_vol);
		    for (i = 0; i < 3; ++i)
		    for (j = 0; j < 3; ++j)
		    {
		    	cell->nb_frac[i][j] *= lambda;
		    }
		    cell->orphan *= lambda;
		}
	    }
	}
	else if (cell->nv == 3)
	{
	    for (i = 0; i < 4; ++i)
	    {
	    	if (cell->comp[0][0] != cell->soln_comp)
	    	{
		    area = nv3_area(cell);
		    if (cell->cell_comp != cell->soln_comp)
		    	area = 0.5*cell->full_cell_vol;
		    if (cell->cell_comp == cell->soln_comp)
	    	    	cell->nb_frac[1][1] = area;
		    else
		    {
		    	if ((cell->side_flag == 0 || !cell->nb_flag[2][1])
		    	    && cell->nb_flag[1][2])
		    	    cell->nb_frac[1][2] = area;
		    	else if (cell->nb_flag[2][1])
		    	    cell->nb_frac[2][1] = area;
		    	else if ((cell->side_flag == 0 || !cell->nb_flag[0][1])
		    	    && cell->nb_flag[1][0])
		    	    cell->nb_frac[1][0] = area;
		    	else if (cell->nb_flag[0][1])
		    	    cell->nb_frac[0][1] = area;
		    	else
			     cell->orphan = area;
		    }
		    if (cell->nb_frac[1][1] != 0.0 || cell->orphan != 0.0)
		    	break;
		    for (j = 0; j < i; ++j)
			reverse_cell(cell);
		    break;
	    	}
	    	rotate_cell(cell,NO);
	    }
	}
}	/* end cell_area */

static void reverse_cell(
	CELL_INFO_2D *cell)
{
	boolean tmp_flag;
	double tmp_frac;
	/* Rotate corners */
	tmp_flag = cell->nb_flag[0][0];
	tmp_frac = cell->nb_frac[0][0];
	cell->nb_flag[0][0] = cell->nb_flag[2][0];
	cell->nb_frac[0][0] = cell->nb_frac[2][0];
	cell->nb_flag[2][0] = cell->nb_flag[2][2];
	cell->nb_frac[2][0] = cell->nb_frac[2][2];
	cell->nb_flag[2][2] = cell->nb_flag[0][2];
	cell->nb_frac[2][2] = cell->nb_frac[0][2];
	cell->nb_flag[0][2] = tmp_flag;
	cell->nb_frac[0][2] = tmp_frac;
	/* Rotate sides */
	tmp_flag = cell->nb_flag[1][0];
	tmp_frac = cell->nb_frac[1][0];
	cell->nb_flag[1][0] = cell->nb_flag[2][1];
	cell->nb_frac[1][0] = cell->nb_frac[2][1];
	cell->nb_flag[2][1] = cell->nb_flag[1][2];
	cell->nb_frac[2][1] = cell->nb_frac[1][2];
	cell->nb_flag[1][2] = cell->nb_flag[0][1];
	cell->nb_frac[1][2] = cell->nb_frac[0][1];
	cell->nb_flag[0][1] = tmp_flag;
	cell->nb_frac[0][1] = tmp_frac;
}	/* end reverse_cell */

static void rotate_cell(	/* Rotate cell counter-clock wise */
	CELL_INFO_2D *cell,
	boolean rotate_vfrac)
{
	int tmp_comp;
	double tmp_crds[MAXD];
	boolean tmp_flag;
	double tmp_frac;
	int i;

	/* Rotate cornres */
	tmp_comp = cell->comp[0][0];
	cell->comp[0][0] = cell->comp[0][1];
	cell->comp[0][1] = cell->comp[1][1];
	cell->comp[1][1] = cell->comp[1][0];
	cell->comp[1][0] = tmp_comp;

	for (i = 0; i < 2; ++i)
	{
	    tmp_crds[i] = cell->crds[0][0][i];
	    cell->crds[0][0][i] = cell->crds[0][1][i];
	    cell->crds[0][1][i] = cell->crds[1][1][i];
	    cell->crds[1][1][i] = cell->crds[1][0][i];
	    cell->crds[1][0][i] = tmp_crds[i];

	/* Rotate edges */
	    tmp_crds[i] = cell->crx_coords[0][0][i];
	    cell->crx_coords[0][0][i] = cell->crx_coords[1][0][i];
	    cell->crx_coords[1][0][i] = cell->crx_coords[0][1][i];
	    cell->crx_coords[0][1][i] = cell->crx_coords[1][1][i];
	    cell->crx_coords[1][1][i] = tmp_crds[i];
	}

	/* Rotate neighbors */

	/* Rotate corners */
	tmp_flag = cell->nb_flag[0][0];
	cell->nb_flag[0][0] = cell->nb_flag[0][2];
	cell->nb_flag[0][2] = cell->nb_flag[2][2];
	cell->nb_flag[2][2] = cell->nb_flag[2][0];
	cell->nb_flag[2][0] = tmp_flag;
	/* Rotate sides */
	tmp_flag = cell->nb_flag[1][0];
	cell->nb_flag[1][0] = cell->nb_flag[0][1];
	cell->nb_flag[0][1] = cell->nb_flag[1][2];
	cell->nb_flag[1][2] = cell->nb_flag[2][1];
	cell->nb_flag[2][1] = tmp_flag;

	if (rotate_vfrac)	/* Needed for diagonal case */
	{
	    /* Rotate corners */
	    tmp_frac = cell->nb_frac[0][0];
	    cell->nb_frac[0][0] = cell->nb_frac[0][2];
	    cell->nb_frac[0][2] = cell->nb_frac[2][2];
	    cell->nb_frac[2][2] = cell->nb_frac[2][0];
	    cell->nb_frac[2][0] = tmp_frac;
	    /* Rotate sides */
	    tmp_frac = cell->nb_frac[1][0];
	    cell->nb_frac[1][0] = cell->nb_frac[0][1];
	    cell->nb_frac[0][1] = cell->nb_frac[1][2];
	    cell->nb_frac[1][2] = cell->nb_frac[2][1];
	    cell->nb_frac[2][1] = tmp_frac;
	}

}	/* end rotate_cell */

static void cell_length(
	CELL_INFO_1D *cell)
{
	double length;

	if (cell->crx_coords[0] == HUGE)
	{
	    printf("ERROR in cell_length(), crx_coords unassigned!\n");
	    clean_up(ERROR);
	}
	if (cell->comp[0] == cell->soln_comp)
	{
	    length = fabs(cell->crx_coords[0] - cell->crds[0][0]);
	    if (cell->nb_flag[1])
	    	cell->nb_frac[1] = length;
	    else if (cell->nb_flag[0])
	    	cell->nb_frac[0] = length;
	    else
	    {
	    	printf("In cell_length(): case 11\n");
		clean_up(ERROR);
	    }
	}
	else if (cell->comp[1] == cell->soln_comp)
	{
	    length = fabs(cell->crx_coords[0] - cell->crds[1][0]);
	    if (cell->nb_flag[1])
	    	cell->nb_frac[1] = length;
	    else if (cell->nb_flag[2])
	    	cell->nb_frac[2] = length;
	    else
	    {
	    	printf("In cell_length(): case 12\n");
		clean_up(ERROR);
	    }
	}
}	/* end cell_length */

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

#include <iFluid.h>
#include <airfoil.h>
#include "solver.h"

#define 	MAX_NUM_RING1		30

extern void count_node_neighbors(
	NODE *node,
	SPRING_VERTEX *sv,
	int *n)
{
	CURVE **c;
        BOND *b;
	POINT *p,*p_nb;
	int num_nb;
	int i,j,dim;

	dim = Dimension(node->interface);
	num_nb = 0;
	for (c = node->out_curves; c && *c; ++c)
	    num_nb++;
	for (c = node->in_curves; c && *c; ++c)
	    num_nb++;
	node->posn->indx = *n;
	if (dim == 3)
	{
	    BOND_TRI **btris;
	    TRI **tris,*tri_list[500];
	    int k,side,nt,num_tris;
	    TRI *tri;

	    num_tris = 0;
	    p = node->posn;
	    for (c = node->out_curves; c && *c; ++c)
	    {
		b = (*c)->first;
		for (btris = Btris(b); btris && *btris; ++btris)
		{
		    nt = FT_FirstRingTrisAroundPoint(p,(*btris)->tri,&tris);
		    for (j = 0; j < nt; ++j)
		    {
			if (!pointer_in_list((POINTER)tris[j],num_tris,
					(POINTER*)tri_list))
			    tri_list[num_tris++] = tris[j];
		    }
		}
	    }
	    for (c = node->in_curves; c && *c; ++c)
	    {
		b = (*c)->last;
		for (btris = Btris(b); btris && *btris; ++btris)
		{
		    nt = FT_FirstRingTrisAroundPoint(p,(*btris)->tri,&tris);
		    for (j = 0; j < nt; ++j)
		    {
			if (!pointer_in_list((POINTER)tris[j],num_tris,
					(POINTER*)tri_list))
			    tri_list[num_tris++] = tris[j];
		    }
		}
	    }
	    for (i = 0; i < num_tris; ++i)
	    {
		tri = tri_list[i];
		for (side = 0; side < 3; ++side)
		{
		    if (p == Point_of_tri(tri)[side])
		    {
			if (is_side_bdry(tri,side))
			    continue;
			num_nb++;
		    }
		}
	    }
	}
	sv[*n].num_nb = num_nb;
	(*n)++;
}	/* end count_node_neighbors */

extern void count_curve_neighbors(
	CURVE *curve,
	SPRING_VERTEX *sv,
	int *n)
{
	int i,j;
	BOND *b;
	int dim = Dimension(curve->interface);

	i = *n;
	for (b = curve->first; b != curve->last; b = b->next)
	{
	    sv[i].num_nb = 2;
	    i++;
	}

	if (dim == 3)
	{
	    POINT *p,*p_nb;
	    BOND_TRI **btris;
	    TRI **tris;
	    int j,side,nt;
	    i = *n;
	    for (b = curve->first; b != curve->last; b = b->next)
	    {
		p = b->end;
		for (btris = Btris(b); btris && *btris; ++btris)
		{
		    nt = FT_FirstRingTrisAroundPoint(p,(*btris)->tri,&tris);
		    for (j = 0; j < nt; ++j)
		    {
			for (side = 0; side < 3; ++side)
			{
			    if (p == Point_of_tri(tris[j])[side])
			    {
				if (is_side_bdry(tris[j],side))
				    continue;
				sv[i].num_nb++;
			    }
			}
		    }
		}
		p->indx = i;
		i++;
	    }
	}
	*n = i;
}	/* end count_curve_neighbors */

extern void count_surf_neighbors(
	SURFACE *surf,
	SPRING_VERTEX *sv,
	int *n)
{
	int i,j,k,nt;
	TRI *tri;
	POINT *p;
	int dim = 3;
	TRI *tris[MAX_NUM_RING1];

	unsort_surf_point(surf);
	i = *n;
	for (tri = first_tri(surf); !at_end_of_tri_list(tri,surf); 
			tri = tri->next)
	{
	    for (j = 0; j < 3; ++j)
	    {
		p = Point_of_tri(tri)[j];
		if (sorted(p) || Boundary_point(p)) continue;
		PointAndFirstRingTris(p,Hyper_surf_element(tri),
			Hyper_surf(surf),&nt,tris);
		sv[i].num_nb = nt;
		sorted(p) = YES;
		p->indx = i;
	    	++i;
	    }
	}
	*n = i;
}	/* end count_surf_neighbors */

extern void set_spring_vertex_memory(
	SPRING_VERTEX *sv,
	int size)
{
	int i,j,num_nb;
	for (i = 0; i < size; ++i)
	{
	    num_nb = sv[i].num_nb;
	    FT_VectorMemoryAlloc((POINTER*)&sv[i].x_nb,num_nb,sizeof(double*));
	    FT_VectorMemoryAlloc((POINTER*)&sv[i].k,num_nb,sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&sv[i].len0,num_nb,sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&sv[i].ix_nb,num_nb,sizeof(int));
	    for (j = 0; j < MAXD; ++j)	// reset external acceleration
		sv[i].ext_accel[j] = 0.0;
	}
}	/* end set_spring_vertex_memory */

extern void set_node_spring_vertex(
	PARACHUTE_SET *geom_set,
	NODE *node,
	double **x,
	double **v,
	SPRING_VERTEX *sv,
	int *n)
{
	CURVE **c;
	BOND *b;
	POINT *p,*p_nb;
	INTERFACE *intfc = geom_set->front->interf;
	int i,j,nn,dim = Dimension(intfc);
	double ks = geom_set->ks;
	double kl = geom_set->kl;
	double kg = geom_set->kg;
	double mass;
	double lambda_s = geom_set->lambda_s;
	double lambda_l = geom_set->lambda_l;
	double lambda_g = geom_set->lambda_g;
	boolean is_fixed = NO;

	if (dim == 3)
	{
	    AF_NODE_EXTRA *extra = (AF_NODE_EXTRA*)node->extra;
	    if (extra != NULL)
	    {
		if (extra->af_node_type == PRESET_NODE)
		    is_fixed = YES;
		else if (extra->af_node_type == LOAD_NODE)
		{
	    	    Front *front = geom_set->front;
	    	    AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
		    mass = af_params->payload;
		}
		else if (extra->af_node_type == GORE_NODE)
                    mass = geom_set->m_g;
		else if (extra->af_node_type == STRING_NODE)
                    mass = geom_set->m_s;
	    }
	    else
                mass = geom_set->m_s;
	}
	else
            mass = geom_set->m_l;

	nn = 0;
	x[*n] = sv[*n].x = Coords(node->posn);
	v[*n] = sv[*n].v = node->posn->vel;
	sv[*n].f = node->posn->force;
	for (c = node->out_curves; c && *c; ++c)
	{
	    b = (*c)->first;
	    sv[*n].x_nb[nn] = Coords(b->end);
	    sv[*n].ix_nb[nn] = b->end->indx;
	    sv[*n].len0[nn] = bond_length0(b);
	    sv[*n].m = mass;
	    if (dim == 3)
	    {
		if (is_fixed)
		    sv[*n].k[nn] = 0.0;
		else if (is_load_node(node) == YES)
		    sv[*n].k[nn] = kl;
		else if (hsbdry_type(*c) == STRING_HSBDRY)
		    sv[*n].k[nn] = kl;
		else if (hsbdry_type(*c) == MONO_COMP_HSBDRY)
		    sv[*n].k[nn] = ks;
		else if (hsbdry_type(*c) == GORE_HSBDRY)
		    sv[*n].k[nn] = kg;
		else if (hsbdry_type(*c) == FIXED_HSBDRY)
		    is_fixed = YES;
	    }
	    else
		sv[*n].k[nn] = kl;
	    ++nn;
	}
	for (c = node->in_curves; c && *c; ++c)
	{
	    b = (*c)->last;
	    sv[*n].x_nb[nn] = Coords(b->start);
	    sv[*n].ix_nb[nn] = b->start->indx;
	    sv[*n].len0[nn] = bond_length0(b);
	    sv[*n].m = mass;
	    if (dim == 3)
	    {
		if (is_fixed)
		    sv[*n].k[nn] = 0.0;
		else if (is_load_node(node) == YES)
		    sv[*n].k[nn] = kl;
		else if (hsbdry_type(*c) == STRING_HSBDRY)
		    sv[*n].k[nn] = kl;
		else if (hsbdry_type(*c) == MONO_COMP_HSBDRY)
		    sv[*n].k[nn] = ks;
		else if (hsbdry_type(*c) == GORE_HSBDRY)
		    sv[*n].k[nn] = kg;
		else if (hsbdry_type(*c) == FIXED_HSBDRY)
		    is_fixed = YES;
	    }
	    else
		sv[*n].k[nn] = kl;
	    ++nn;
	}
	if (dim == 3)
	{
	    BOND_TRI **btris;
	    TRI **tris,*tri_list[500];
	    int k,side,nt,num_tris;
	    TRI *tri;

	    num_tris = 0;
	    p = node->posn;
	    for (c = node->out_curves; c && *c; ++c)
	    {
		b = (*c)->first;
		for (btris = Btris(b); btris && *btris; ++btris)
		{
		    nt = FT_FirstRingTrisAroundPoint(p,(*btris)->tri,&tris);
		    for (j = 0; j < nt; ++j)
		    {
			if (!pointer_in_list((POINTER)tris[j],num_tris,
					(POINTER*)tri_list))
			    tri_list[num_tris++] = tris[j];
		    }
		}
	    }
	    for (c = node->in_curves; c && *c; ++c)
	    {
		b = (*c)->last;
		for (btris = Btris(b); btris && *btris; ++btris)
		{
		    nt = FT_FirstRingTrisAroundPoint(p,(*btris)->tri,&tris);
		    for (j = 0; j < nt; ++j)
		    {
			if (!pointer_in_list((POINTER)tris[j],num_tris,
					(POINTER*)tri_list))
			    tri_list[num_tris++] = tris[j];
		    }
		}
	    }
	    for (i = 0; i < num_tris; ++i)
	    {
		tri = tri_list[i];
		for (side = 0; side < 3; ++side)
		{
		    if (p == Point_of_tri(tri)[side])
		    {
			if (is_side_bdry(tri,side))
			    continue;
			p_nb = Point_of_tri(tri)[(side+1)%3];
			sv[*n].x_nb[nn] = Coords(p_nb);
			sv[*n].ix_nb[nn] = p_nb->indx;
			sv[*n].k[nn] = ks;
			if (is_fixed) sv[*n].k[nn] = 0.0;
			sv[*n].len0[nn] = tri->side_length0[side];
			++nn;
		    }
		}
	    }
	    if (!is_load_node(node))
	    {
		sv[*n].lambda = lambda_s;
	    }
	    if (is_fixed) 
	    {
		sv[*n].lambda = 0.0;
	    	for (i = 0; i < sv[*n].num_nb; ++i)
		    sv[*n].k[nn] = 0.0;
	    }
	}
	else
	{
	    sv[*n].lambda = lambda_l;
	}
	(*n)++;
}	/* end set_node_spring_vertex */

extern void set_curve_spring_vertex(
	PARACHUTE_SET *geom_set,
	CURVE *curve,
	double **x,
	double **v,
	SPRING_VERTEX *sv,
	int *n)
{
	int i,j,nn;
	BOND *b;
	int dim = Dimension(curve->interface);
	double kl,m_l,lambda_l;

	if (dim == 3)
	{
	    if (hsbdry_type(curve) == STRING_HSBDRY)
	    {
	    	kl = geom_set->kl;
	    	m_l = geom_set->m_l;
	    	lambda_l = geom_set->lambda_l;
	    }
	    else if (hsbdry_type(curve) == GORE_HSBDRY)
	    {
	    	kl = geom_set->kg;
	    	m_l = geom_set->m_g;
	    	lambda_l = geom_set->lambda_g;
	    }
	    else if (hsbdry_type(curve) == FIXED_HSBDRY)
	    {
	    	kl = 0.0;
	    	m_l = geom_set->m_l;
	    	lambda_l = 0.0;
	    }
	    else
	    {
	    	kl = geom_set->ks;
	    	m_l = geom_set->m_s;
	    	lambda_l = geom_set->lambda_s;
	    }
	}
	else
	{
	    kl = geom_set->kl;
	    m_l = geom_set->m_l;
	    lambda_l = geom_set->lambda_l;
	}

	i = *n;
	for (b = curve->first; b != curve->last; b = b->next)
	{
	    x[i] = sv[i].x = Coords(b->end);
	    v[i] = sv[i].v = b->end->vel;
	    sv[i].f = b->end->force;
	    sv[i].x_nb[0] = Coords(b->start);
	    sv[i].x_nb[1] = Coords(b->next->end);
	    sv[i].ix_nb[0] = b->start->indx;
	    sv[i].ix_nb[1] = b->next->end->indx;
	    sv[i].len0[0] = bond_length0(b);
	    sv[i].len0[1] = bond_length0(b->next);
	    sv[i].k[0] = sv[i].k[1] = kl;
	    sv[i].m = m_l;
	    sv[i].num_nb = 2;
	    sv[i].lambda = lambda_l;
	    ++i;
	}

	if (dim == 3)
	{
	    POINT *p,*p_nb;
	    BOND_TRI **btris;
	    TRI **tris;
	    int j,k,side,nt;
	    double ks;

	    if (hsbdry_type(curve) == FIXED_HSBDRY)
		ks = 0.0;
	    else
		ks = geom_set->ks;
	    i = *n;
	    for (b = curve->first; b != curve->last; b = b->next)
	    {
		p = b->end;
		nn = sv[i].num_nb;
		sv[i].m = m_l;
		for (btris = Btris(b); btris && *btris; ++btris)
		{
		    nt = FT_FirstRingTrisAroundPoint(p,(*btris)->tri,&tris);
		    for (j = 0; j < nt; ++j)
		    {
			for (side = 0; side < 3; ++side)
			{
			    if (p == Point_of_tri(tris[j])[side])
			    {
				if (is_side_bdry(tris[j],side))
				    continue;
				p_nb = Point_of_tri(tris[j])[(side+1)%3];
				sv[i].x_nb[nn] = Coords(p_nb);
				sv[i].ix_nb[nn] = p_nb->indx;
				sv[i].k[nn] = ks;
				sv[i].len0[nn] = tris[j]->side_length0[side];
				++nn;
			    }
			}
		    }
		}
		sv[i].num_nb = nn;
		i++;
	    }
	}
	*n = i;
}	/* end set_curve_spring_vertex */

extern void set_surf_spring_vertex(
	PARACHUTE_SET *geom_set,
	SURFACE *surf,
	double **x,
	double **v,
	SPRING_VERTEX *sv,
	int *n)
{
	int i,j,k,l,nt;
	TRI *tri;
	TRI *tris[MAX_NUM_RING1];
	POINT *p,*p_nb;
	double ks = geom_set->ks;
	double m_s = geom_set->m_s;
	double lambda_s = geom_set->lambda_s;

	unsort_surf_point(surf);
	i = *n;
	for (tri = first_tri(surf); !at_end_of_tri_list(tri,surf); 
			tri = tri->next)
	{
	    for (j = 0; j < 3; ++j)
	    {
		p = Point_of_tri(tri)[j];
		if (sorted(p) || Boundary_point(p)) continue;
		sv[i].m = m_s;
		sv[i].lambda = lambda_s;
		x[i] = sv[i].x = Coords(p);
		v[i] = sv[i].v = p->vel;
		sv[i].f = p->force;
		PointAndFirstRingTris(p,Hyper_surf_element(tri),
				Hyper_surf(surf),&nt,tris);
		sv[i].num_nb = nt;
		for (k = 0; k < nt; ++k)
		for (l = 0; l < 3; ++l)
		if (Point_of_tri(tris[k])[l] == p)
		{
		    p_nb = Point_of_tri(tris[k])[(l+1)%3];
		    sv[i].x_nb[k] = Coords(p_nb);
		    sv[i].ix_nb[k] = p_nb->indx;
		    sv[i].k[k] = ks;
		    sv[i].len0[k] = tris[k]->side_length0[l];;
		}
		sorted(p) = YES;
	    	++i;
	    }
	}
	*n = i;
}	/* end set_surf_spring_vertex */

extern void compute_spring_accel1(
	SPRING_VERTEX sv,
	double *f,
	int dim)
{
	int i,k;
	double len,vec[MAXD];

	for (k = 0; k < dim; ++k)
	    f[k] = 0.0;
	for (i = 0; i < sv.num_nb; ++i)
	{
	    len = 0.0;
	    for (k = 0; k < dim; ++k)
	    {
		vec[k] = sv.x_nb[i][k] - sv.x[k];
		len += sqr(vec[k]);
	    }
	    len = sqrt(len);
	    for (k = 0; k < dim; ++k)
	    {
		vec[k] /= len;
		f[k] += sv.k[i]*((len - sv.len0[i])*vec[k])/sv.m;
	    }
	}
	for (k = 0; k < dim; ++k)
	    sv.f[k] = f[k]*sv.m;
	for (k = 0; k < dim; ++k)
	{
	    f[k] += -sv.lambda*sv.v[k]/sv.m + sv.ext_accel[k];
	}
}	/* end compute_spring_accel */

extern void count_canopy_spring_neighbors(
	PARACHUTE_SET *geom_set,
	SPRING_VERTEX *sv)
{
	int i,n = 0;
	int ns,nbc,ngc,ng;
	int n_start,n_end;

	if (debugging("canopy"))
	    (void) printf("Entering count_canopy_spring_neighbors()\n");
	ns = geom_set->num_strings;
	nbc = geom_set->num_mono_hsbdry;
	ngc = geom_set->num_gore_hsbdry;
	ng = geom_set->num_gore_nodes;

	count_surf_neighbors(geom_set->canopy,sv,&n);
	for (i = 0; i < ng; ++i)
            count_node_neighbors(geom_set->gore_nodes[i],sv,&n);
	for (i = 0; i < ns; ++i)
	    count_node_neighbors(geom_set->string_node[i],sv,&n);
	for (i = 0; i < ngc; ++i)
	    count_curve_neighbors(geom_set->gore_hsbdry[i],sv,&n);
	for (i = 0; i < nbc; ++i)
	{
	    count_curve_neighbors(geom_set->mono_hsbdry[i],sv,&n);
	    if (is_closed_curve(geom_set->mono_hsbdry[i]))
		count_node_neighbors(geom_set->mono_hsbdry[i]->start,sv,&n);	
	}
	if (debugging("canopy"))
	    (void) printf("Leaving count_canopy_spring_neighbors()\n");
}	/* end count_canopy_spring_neighbors */

extern void count_string_spring_neighbors(
	PARACHUTE_SET *geom_set,
	SPRING_VERTEX *sv)
{
	int i,n;
	int ns = geom_set->num_strings;

	if (debugging("canopy"))
	    (void) printf("Entering count_string_spring_neighbors()\n");

	n = geom_set->n_cps;
	count_node_neighbors(geom_set->load_node,sv,&n);
	for (i = 0; i < ns; ++i)
	{
	    count_curve_neighbors(geom_set->string_curves[i],sv,&n);
	}

	if (debugging("canopy"))
	    (void) printf("Leaving count_string_spring_neighbors()\n");
}	/* end  count_string_spring_neighbors */

extern void set_canopy_spring_vertex(
	PARACHUTE_SET *geom_set,
	double **x,
	double **v,
	SPRING_VERTEX *sv)
{
	int i,n = 0;
	int ns,nbc,ngc,ng;
	int n_start,n_end;
	Front *fr = geom_set->front;

	if (debugging("canopy"))
	    (void) printf("Entering set_canopy_spring_vertex()\n");

	ns = geom_set->num_strings;
	nbc = geom_set->num_mono_hsbdry;
	ngc = geom_set->num_gore_hsbdry;
	ng = geom_set->num_gore_nodes;

	set_surf_spring_vertex(geom_set,geom_set->canopy,x,v,sv,&n);
	for (i = 0; i < ng; ++i)
            set_node_spring_vertex(geom_set,geom_set->gore_nodes[i],x,v,sv,&n);
	for (i = 0; i < ns; ++i)
	    set_node_spring_vertex(geom_set,geom_set->string_node[i],x,v,sv,&n);
	for (i = 0; i < ngc; ++i)
	    set_curve_spring_vertex(geom_set,geom_set->gore_hsbdry[i],
					x,v,sv,&n);
	for (i = 0; i < nbc; ++i)
	{
	    set_curve_spring_vertex(geom_set,geom_set->mono_hsbdry[i],x,
				v,sv,&n);
	    if (is_closed_curve(geom_set->mono_hsbdry[i]))
		set_node_spring_vertex(geom_set,
				geom_set->mono_hsbdry[i]->start,x,v,sv,&n);	
	}
	if (debugging("canopy"))
	    (void) printf("Leaving set_canopy_spring_vertex()\n");
}	/* end set_canopy_spring_vertex */

extern void set_string_spring_vertex(
	PARACHUTE_SET *geom_set,
	double **x,
	double **v,
	SPRING_VERTEX *sv)
{
	int i,n;
	int ns = geom_set->num_strings;

	if (debugging("canopy"))
	    (void) printf("Entering set_string_spring_vertex()\n");

	n = geom_set->n_cps;
	set_node_spring_vertex(geom_set,geom_set->load_node,x,v,sv,&n);
	for (i = 0; i < ns; ++i)
	{
	    set_curve_spring_vertex(geom_set,geom_set->string_curves[i],x,v,
					sv,&n);
	}

	if (debugging("canopy"))
	    (void) printf("Leaving set_string_spring_vertex()\n");
}	/* end  set_string_spring_vertex */

extern void generic_spring_solver(
	SPRING_VERTEX *sv,
	double **x_pos,
	double **v_pos,
	int dim,
	int size,
	int n_loop,
	double dt)
{
	static double **x_old,**x_new,**v_old,**v_new,**accel;
	int i,j,n;
	
	if (debugging("trace"))
	    (void) printf("Entering generic_spring_solver()\n");
	if (x_old == NULL)
	{
	    FT_MatrixMemoryAlloc((POINTER*)&x_old,size,3,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&v_old,size,3,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&x_new,size,3,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&v_new,size,3,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&accel,size,3,sizeof(double));
	}

	for (n = 0; n < n_loop; ++n)
	{
	    for (i = 0; i < size; ++i)
	    for (j = 0; j < dim; ++j)
	    {
		x_old[i][j] = x_pos[i][j];
		v_old[i][j] = v_pos[i][j];
	    }

	    for (i = 0; i < size; ++i)
	    {
		compute_spring_accel1(sv[i],accel[i],dim);
	    }
	    for (i = 0; i < size; ++i)
	    for (j = 0; j < dim; ++j)
	    {
		x_new[i][j] = x_old[i][j] + dt*v_old[i][j]/6.0;
                v_new[i][j] = v_old[i][j] + dt*accel[i][j]/6.0;
	    	x_pos[i][j] = x_old[i][j] + 0.5*v_old[i][j]*dt;
	    	v_pos[i][j] = v_old[i][j] + 0.5*accel[i][j]*dt;
	    }

	    for (i = 0; i < size; ++i)
		compute_spring_accel1(sv[i],accel[i],dim);
	    for (i = 0; i < size; ++i)
	    for (j = 0; j < dim; ++j)
	    {
		x_new[i][j] += dt*v_pos[i][j]/3.0;
                v_new[i][j] += dt*accel[i][j]/3.0;
	    	x_pos[i][j] = x_old[i][j] + 0.5*v_pos[i][j]*dt;
	    	v_pos[i][j] = v_old[i][j] + 0.5*accel[i][j]*dt;
	    }
	
	    for (i = 0; i < size; ++i)
		compute_spring_accel1(sv[i],accel[i],dim);
	    for (i = 0; i < size; ++i)
	    for (j = 0; j < dim; ++j)
	    {
		x_new[i][j] += dt*v_pos[i][j]/3.0;
                v_new[i][j] += dt*accel[i][j]/3.0;
	    	x_pos[i][j] = x_old[i][j] + v_pos[i][j]*dt;
	    	v_pos[i][j] = v_old[i][j] + accel[i][j]*dt; 
	    }

	    for (i = 0; i < size; ++i)
		compute_spring_accel1(sv[i],accel[i],dim);
	    for (i = 0; i < size; ++i)
	    for (j = 0; j < dim; ++j)
	    {
		x_new[i][j] += dt*v_pos[i][j]/6.0;
                v_new[i][j] += dt*accel[i][j]/6.0;
	    }
	    for (i = 0; i < size; ++i)
	    for (j = 0; j < dim; ++j)
	    {
		x_pos[i][j] = x_new[i][j];
                v_pos[i][j] = v_new[i][j];
	    }

	    if (n != n_loop-1)
	    {
	    	for (i = 0; i < size; ++i)
		    compute_spring_accel1(sv[i],accel[i],dim);
	    }
	}
	if (debugging("trace"))
	    (void) printf("Leaving generic_spring_solver()\n");
}	/* end generic_spring_solver */

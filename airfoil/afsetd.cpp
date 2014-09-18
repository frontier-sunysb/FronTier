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

static void link_surf_point_set(ELASTIC_SET*,SURFACE*,GLOBAL_POINT**,
				GLOBAL_POINT*,int*);
static void link_curve_point_set(ELASTIC_SET*,CURVE*,GLOBAL_POINT**,
				GLOBAL_POINT*,int*);
static void link_node_point_set(ELASTIC_SET*,NODE*,GLOBAL_POINT**,
				GLOBAL_POINT*,int*);

static void surf_get_point_set_from(SURFACE*,GLOBAL_POINT**);
static void curve_get_point_set_from(CURVE*,GLOBAL_POINT**);
static void node_get_point_set_from(NODE*,GLOBAL_POINT**);

static void surf_put_point_set_to(SURFACE*,GLOBAL_POINT**);
static void curve_put_point_set_to(CURVE*,GLOBAL_POINT**);
static void node_put_point_set_to(NODE*,GLOBAL_POINT**);

static void count_surf_neighbors(SURFACE*,SPRING_VERTEX*,int*);
static void count_curve_neighbors(CURVE*,SPRING_VERTEX*,int*);
static void count_node_neighbors(NODE*,SPRING_VERTEX*,int*);

static void set_surf_impulse(ELASTIC_SET*,SURFACE*,SPRING_VERTEX*,int*);
static void set_curve_impulse(ELASTIC_SET*,CURVE*,SPRING_VERTEX*,int*);
static void set_node_impulse(ELASTIC_SET*,NODE*,SPRING_VERTEX*,int*);

static void get_point_value_from(POINT*,GLOBAL_POINT**);
static void put_point_value_to(POINT*,GLOBAL_POINT**);

#define 	MAX_NUM_RING1		30

static void count_node_neighbors(
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
	sv[*n].ix = node->posn->indx = *n;
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

static void count_curve_neighbors(
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
		sv[i].ix = p->indx = i;
		i++;
	    }
	}
	*n = i;
}	/* end count_curve_neighbors */

static void count_surf_neighbors(
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
		sv[i].ix = p->indx = i;
	    	++i;
		sorted(p) = YES;
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

extern void compute_spring_accel1(
	SPRING_VERTEX *sv,
	double *f,
	int dim)
{
	int i,k;
	double len,vec[MAXD];

	for (k = 0; k < dim; ++k)
	    f[k] = 0.0;
	for (i = 0; i < sv->num_nb; ++i)
	{
	    len = 0.0;
	    for (k = 0; k < dim; ++k)
	    {
		vec[k] = sv->x_nb[i][k] - sv->x[k];
		len += sqr(vec[k]);
	    }
	    len = sqrt(len);
	    for (k = 0; k < dim; ++k)
	    {
		vec[k] /= len;
		f[k] += sv->k[i]*((len - sv->len0[i])*vec[k])/sv->m;
	    }
	}
	for (k = 0; k < dim; ++k)
	    sv->f[k] = f[k]*sv->m;
	for (k = 0; k < dim; ++k)
	{
	    f[k] += -sv->lambda*sv->v[k]/sv->m;
	}
}	/* end compute_spring_accel */

extern void generic_spring_solver(
	SPRING_VERTEX *sv,
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

	for (i = 0; i < size; ++i)
	{
	    compute_spring_accel1(&sv[i],accel[i],dim);
	}
	for (i = 0; i < size; ++i)
	for (j = 0; j < dim; ++j)
	{
	    x_old[i][j] = sv[i].x[j];
	    v_old[i][j] = sv[i].v[j];
	}
	for (n = 0; n < n_loop; ++n)
	{
	    for (i = 0; i < size; ++i)
	    for (j = 0; j < dim; ++j)
	    {
		x_new[i][j] = x_old[i][j] + dt*v_old[i][j]/6.0;
                v_new[i][j] = v_old[i][j] + dt*accel[i][j]/6.0;
	    	sv[i].x[j] = x_old[i][j] + 0.5*v_old[i][j]*dt;
	    	sv[i].v[j] = v_old[i][j] + 0.5*accel[i][j]*dt;
	    }

	    for (i = 0; i < size; ++i)
	    {
		compute_spring_accel1(&sv[i],accel[i],dim);
	    }
	    for (i = 0; i < size; ++i)
	    for (j = 0; j < dim; ++j)
	    {
		x_new[i][j] += dt*sv[i].v[j]/3.0;
                v_new[i][j] += dt*accel[i][j]/3.0;
	    	sv[i].x[j] = x_old[i][j] + 0.5*sv[i].v[j]*dt;
	    	sv[i].v[j] = v_old[i][j] + 0.5*accel[i][j]*dt;
	    }
	
	    for (i = 0; i < size; ++i)
	    {
		compute_spring_accel1(&sv[i],accel[i],dim);
	    }
	    for (i = 0; i < size; ++i)
	    for (j = 0; j < dim; ++j)
	    {
		x_new[i][j] += dt*sv[i].v[j]/3.0;
                v_new[i][j] += dt*accel[i][j]/3.0;
	    	sv[i].x[j] = x_old[i][j] + sv[i].v[j]*dt;
	    	sv[i].v[j] = v_old[i][j] + accel[i][j]*dt; 
	    }

	    for (i = 0; i < size; ++i)
	    {
		compute_spring_accel1(&sv[i],accel[i],dim);
	    }
	    for (i = 0; i < size; ++i)
	    for (j = 0; j < dim; ++j)
	    {
		x_new[i][j] += dt*sv[i].v[j]/6.0;
                v_new[i][j] += dt*accel[i][j]/6.0;
	    }
	    for (i = 0; i < size; ++i)
            for (j = 0; j < 3; ++j)
            {
                x_new[i][j] += (sv[i].ext_impul[j]
                                + 0.5*sv[i].ext_accel[j]*dt)*dt;
                sv[i].ext_impul[j] += sv[i].ext_accel[j]*dt;
            }
	    for (i = 0; i < size; ++i)
	    for (j = 0; j < dim; ++j)
	    {
		sv[i].x[j] = x_new[i][j];
                sv[i].v[j] = v_new[i][j];
	    }

	    if (n != n_loop-1)
	    {
		for (i = 0; i < size; ++i)
                for (j = 0; j < 3; ++j)
                {
                    x_old[i][j] = sv[i].x[j];
                    v_old[i][j] = sv[i].v[j];
                }
	    	for (i = 0; i < size; ++i)
		{
		    compute_spring_accel1(&sv[i],accel[i],dim);
		}
	    }
	}
	if (debugging("trace"))
	    (void) printf("Leaving generic_spring_solver()\n");
}	/* end generic_spring_solver */

extern void set_vertex_impulse(
        ELASTIC_SET *geom_set,
        SPRING_VERTEX *sv)
{
        int i,n,ns,nc,nn;

        ns = geom_set->num_surfs;
        nc = geom_set->num_curves;
        nn = geom_set->num_nodes;
        n = 0;
        for (i = 0; i < ns; ++i)
            set_surf_impulse(geom_set,geom_set->surfs[i],sv,&n);
        for (i = 0; i < nc; ++i)
            set_curve_impulse(geom_set,geom_set->curves[i],sv,&n);
        for (i = 0; i < nn; ++i)
            set_node_impulse(geom_set,geom_set->nodes[i],sv,&n);

}       /* end set_vertex_impulse */

static void set_node_impulse(
	ELASTIC_SET *geom_set,
	NODE *node,
	SPRING_VERTEX *sv,
	int *n)
{
	int i,dim;
	STATE *sl,*sr;

	dim = FT_Dimension();
	sl = (STATE*)left_state(node->posn);
	sr = (STATE*)right_state(node->posn);
	for (i = 0; i < dim; ++i)
	    sl->impulse[i] = sr->impulse[i] = sv[*n].ext_impul[i];
	(*n)++;
}	/* end set_node_impulse */

static void set_curve_impulse(
	ELASTIC_SET *geom_set,
	CURVE *curve,
	SPRING_VERTEX *sv,
	int *n)
{
	int i,j,dim;
	STATE *sl,*sr;
	BOND *b;

	dim = FT_Dimension();

	i = *n;
	for (b = curve->first; b != curve->last; b = b->next)
	{
	    sl = (STATE*)left_state(b->end);
	    sr = (STATE*)right_state(b->end);
            for (j = 0; j < dim; ++j)
            {
	    	sl->impulse[j] = sr->impulse[j] = sv[i].ext_impul[j];
	    }
	    ++i;
	}
	*n = i;
}	/* end set_curve_impulse */

static void set_surf_impulse(
	ELASTIC_SET *geom_set,
	SURFACE *surf,
	SPRING_VERTEX *sv,
	int *n)
{
	int i,j,k;
	TRI *tri;
	POINT *p;
	HYPER_SURF *hs;
	HYPER_SURF_ELEMENT *hse;
	STATE *sl,*sr;

	unsort_surf_point(surf);
	hs = Hyper_surf(surf);
	i = *n;
	for (tri = first_tri(surf); !at_end_of_tri_list(tri,surf); 
			tri = tri->next)
	{
	    hse = Hyper_surf_element(tri);
	    for (j = 0; j < 3; ++j)
	    {
		p = Point_of_tri(tri)[j];
		if (sorted(p) || Boundary_point(p)) continue;
		sorted(p) = YES;
		FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
            	for (k = 0; k < 3; ++k)
            	{
	    	    sl->impulse[k] = sr->impulse[k] = sv[i].ext_impul[k];
	    	}
	    	++i;
	    }
	}
	*n = i;
}	/* end set_surf_impulse */

extern void count_vertex_neighbors(
	ELASTIC_SET *geom_set,
	SPRING_VERTEX *sv)
{
	int i,n,ns,nc,nn;

	if (debugging("canopy"))
	    (void) printf("Entering count_vertex_neighbors()\n");

	ns = geom_set->num_surfs;
	nc = geom_set->num_curves;
	nn = geom_set->num_nodes;
	n = 0;
	for (i = 0; i < ns; ++i)
	    count_surf_neighbors(geom_set->surfs[i],sv,&n);
	for (i = 0; i < nc; ++i)
	    count_curve_neighbors(geom_set->curves[i],sv,&n);
	for (i = 0; i < nn; ++i)
	    count_node_neighbors(geom_set->nodes[i],sv,&n);	

	if (debugging("canopy"))
	    (void) printf("Leaving count_vertex_neighbors()\n");
}	/* end  count_vertex_neighbors */

extern void link_point_set(
	ELASTIC_SET *geom_set,
	GLOBAL_POINT **point_set,
	GLOBAL_POINT *point_set_store)
{
	int i,n,ns,nc,nn;

	if (debugging("canopy"))
	    (void) printf("Entering link_point_set()\n");

	ns = geom_set->num_surfs;
	nc = geom_set->num_curves;
	nn = geom_set->num_nodes;
	n = 0;
	for (i = 0; i < ns; ++i)
	    link_surf_point_set(geom_set,geom_set->surfs[i],point_set,
				point_set_store,&n);
	for (i = 0; i < nc; ++i)
	    link_curve_point_set(geom_set,geom_set->curves[i],point_set,
				point_set_store,&n);
	for (i = 0; i < nn; ++i)
	    link_node_point_set(geom_set,geom_set->nodes[i],point_set,
				point_set_store,&n);

	if (debugging("canopy"))
	{
	    (void) printf("Final n = %d\n",n);
	    (void) printf("Leaving link_point_set()\n");
	}
}	/* end link_point_set */

static void link_surf_point_set(
	ELASTIC_SET *geom_set,
	SURFACE *surf,
	GLOBAL_POINT **point_set,
	GLOBAL_POINT *point_set_store,
	int *n)
{
	TRI *tri;
	POINT *p;
	long gindex;
	int i,j;

	unsort_surf_point(surf);
	i = *n;
	surf_tri_loop(surf,tri)
	{
	    for (j = 0; j < 3; ++j)
	    {
		p = Point_of_tri(tri)[j];
		if (sorted(p) || Boundary_point(p)) continue;
		gindex = Gindex(p);
		point_set[gindex] = point_set_store + i;
		sorted(p) = YES;
		i++;
	    }
	}
	*n = i;
}	/* end link_surf_point_set */

static void link_curve_point_set(
	ELASTIC_SET *geom_set,
	CURVE *curve,
	GLOBAL_POINT **point_set,
	GLOBAL_POINT *point_set_store,
	int *n)
{
	BOND *b;
	POINT *p;
	long gindex;
	int i = *n;

	for (b = curve->first; b != curve->last; b = b->next)
	{
	    p = b->end;
	    gindex = Gindex(p);
	    point_set[gindex] = point_set_store + i;
	    point_set[gindex]->gindex = gindex;
	    i++;
	}
	*n = i;
}	/* end link_curve_point_set */

static void link_node_point_set(
	ELASTIC_SET *geom_set,
	NODE *node,
	GLOBAL_POINT **point_set,
	GLOBAL_POINT *point_set_store,
	int *n)
{
	long gindex;
	POINT *p = node->posn;

	gindex = Gindex(p);
	point_set[gindex] = point_set_store + (*n);
	point_set[gindex]->gindex = gindex;
	(*n)++;
}	/* end link_node_point_set */

extern void set_vertex_neighbors(
	ELASTIC_SET *geom_set,
	SPRING_VERTEX *sv,
	GLOBAL_POINT **point_set)
{
	int i,n,ns,nc,nn;

	if (debugging("canopy"))
	    (void) printf("Entering set_vertex_neighbors()\n");

	ns = geom_set->num_surfs;
	nc = geom_set->num_curves;
	nn = geom_set->num_nodes;
	n = 0;
	for (i = 0; i < ns; ++i)
	    set_surf_spring_vertex(geom_set,geom_set->surfs[i],sv,&n,
					point_set);
	for (i = 0; i < nc; ++i)
	    set_curve_spring_vertex(geom_set,geom_set->curves[i],sv,&n,
					point_set);
	for (i = 0; i < nn; ++i)
	    set_node_spring_vertex(geom_set,geom_set->nodes[i],sv,&n,
					point_set);

	if (debugging("canopy"))
	    (void) printf("Leaving set_vertex_neighbors()\n");
}	/* end  set_vertex_neighbors */

extern void set_node_spring_vertex(
	ELASTIC_SET *geom_set,
	NODE *node,
	SPRING_VERTEX *sv,
	int *n,
	GLOBAL_POINT **point_set)
{
	CURVE **c;
	BOND *b;
	POINT *p,*p_nb;
	Front *front = geom_set->front;
	INTERFACE *intfc = front->interf;
	int i,j,nn,dim = Dimension(intfc);
	double ks = geom_set->ks;
	double kl = geom_set->kl;
	double kg = geom_set->kg;
	double mass;
	double lambda_s = geom_set->lambda_s;
	double lambda_l = geom_set->lambda_l;
	double lambda_g = geom_set->lambda_g;
	boolean is_fixed = NO;
	STATE *sl,*sr;
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	double *g;
	long gindex,gindex_nb;

	if (iFparams != NULL)
 	    g = iFparams->gravity;
	else
	    g == NULL;
	if (dim == 3)
	{
	    AF_NODE_EXTRA *extra = (AF_NODE_EXTRA*)node->extra;
	    if (extra != NULL)
	    {
		if (extra->af_node_type == PRESET_NODE)
		{
                    mass = geom_set->m_s;
		    is_fixed = YES;
		}
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
	gindex = Gindex(node->posn);
	sv[*n].x = point_set[gindex]->x;
	sv[*n].v = point_set[gindex]->v;
	sv[*n].f = point_set[gindex]->f;
	sv[*n].ext_impul = point_set[gindex]->impuls;
	for (c = node->out_curves; c && *c; ++c)
	{
	    b = (*c)->first;
	    gindex_nb = Gindex(b->end);
	    sv[*n].x_nb[nn] = point_set[gindex_nb]->x;
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
	    gindex_nb = Gindex(b->start);
	    sv[*n].x_nb[nn] = point_set[gindex_nb]->x;
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
			gindex_nb = Gindex(p_nb);
			sv[*n].x_nb[nn] = point_set[gindex_nb]->x;
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
		    sv[*n].k[i] = 0.0;
	    }
	}
	else
	{
	    sv[*n].lambda = lambda_l;
	}
	if (dim == 3)
	{
            for (i = 0; i < dim; ++i)
            {
		if (is_fixed || g == NULL)
	    	    sv[*n].ext_accel[i] = 0;
		else
	    	    sv[*n].ext_accel[i] = g[i];
	    }
	}
	(*n)++;
}	/* end set_node_spring_vertex */

extern void set_curve_spring_vertex(
	ELASTIC_SET *geom_set,
	CURVE *curve,
	SPRING_VERTEX *sv,
	int *n,
	GLOBAL_POINT **point_set)
{
	Front *front = geom_set->front;
	int i,j,nn;
	BOND *b;
	double kl,m_l,lambda_l;
	int dim = front->rect_grid->dim;
        IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
        double *g;
	long gindex,gindex_nb;

	if (iFparams != NULL)
	    g = iFparams->gravity;
	else
	    g = NULL;
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
	    gindex = Gindex(b->end);
	    sv[i].x = point_set[gindex]->x;
	    sv[i].v = point_set[gindex]->v;
	    sv[i].f = point_set[gindex]->f;
	    sv[i].ext_impul = point_set[gindex]->impuls;
	    gindex_nb = Gindex(b->start);
	    sv[i].x_nb[0] = point_set[gindex_nb]->x;
	    gindex_nb = Gindex(b->next->end);
	    sv[i].x_nb[1] = point_set[gindex_nb]->x;
	    sv[i].ix_nb[0] = b->start->indx;
	    sv[i].ix_nb[1] = b->next->end->indx;
	    sv[i].len0[0] = bond_length0(b);
	    sv[i].len0[1] = bond_length0(b->next);
	    sv[i].k[0] = sv[i].k[1] = kl;
	    sv[i].m = m_l;
	    sv[i].num_nb = 2;
	    sv[i].lambda = lambda_l;
	    if (dim == 3)
	    {
	    	if (hsbdry_type(curve) == FIXED_HSBDRY || g == NULL)
	    	{
            	    for (j = 0; j < dim; ++j)
	    	    	sv[i].ext_accel[j] = 0;
	    	}
		else
		{
            	    for (j = 0; j < dim; ++j)
	    	    	sv[i].ext_accel[j] = g[j];
		}
	    }
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
				gindex_nb = Gindex(p_nb);
				sv[i].x_nb[nn] = point_set[gindex_nb]->x;
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
	ELASTIC_SET *geom_set,
	SURFACE *surf,
	SPRING_VERTEX *sv,
	int *n,
	GLOBAL_POINT **point_set)
{
	Front *front = geom_set->front;
	int i,j,k,l,nt;
	TRI *tri;
	TRI *tris[MAX_NUM_RING1];
	POINT *p,*p_nb;
	double ks = geom_set->ks;
	double m_s = geom_set->m_s;
	double lambda_s = geom_set->lambda_s;
	boolean is_stationary_point;
	int dim = front->rect_grid->dim;
        IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
        double *g;
	STATE *sl,*sr;
	HYPER_SURF_ELEMENT *hse;
        HYPER_SURF         *hs = Hyper_surf(surf);
	long gindex,gindex_nb;

	if (iFparams != NULL)
 	    g = iFparams->gravity;
	else g = NULL;

	unsort_surf_point(surf);
	i = *n;
	for (tri = first_tri(surf); !at_end_of_tri_list(tri,surf); 
			tri = tri->next)
	{
	    hse = Hyper_surf_element(tri);
	    for (j = 0; j < 3; ++j)
	    {
		p = Point_of_tri(tri)[j];
		if (sorted(p) || Boundary_point(p)) continue;
		is_stationary_point = is_registered_point(surf,p);
		sv[i].m = m_s;
		sv[i].lambda = lambda_s;
		if (is_stationary_point == YES)
		    sv[i].lambda = 0.0;
            	for (k = 0; k < dim; ++k)
            	{
		    if (is_stationary_point == YES || g == NULL)
	    	    	sv[i].ext_accel[k] = 0.0;
		    else
			sv[i].ext_accel[k] = g[k];
	    	}
		gindex = Gindex(p);
		sv[i].x = point_set[gindex]->x;
		sv[i].v = point_set[gindex]->v;
		sv[i].f = point_set[gindex]->f;
	    	sv[i].ext_impul = point_set[gindex]->impuls;
		PointAndFirstRingTris(p,Hyper_surf_element(tri),
				Hyper_surf(surf),&nt,tris);
		sv[i].num_nb = nt;
		for (k = 0; k < nt; ++k)
		for (l = 0; l < 3; ++l)
		if (Point_of_tri(tris[k])[l] == p)
		{
		    p_nb = Point_of_tri(tris[k])[(l+1)%3];
		    gindex_nb = Gindex(p_nb);
		    sv[i].x_nb[k] = point_set[gindex_nb]->x;
		    sv[i].ix_nb[k] = p_nb->indx;
		    if (is_stationary_point == YES)
		    	sv[i].k[k] = 0.0;
		    else
		    	sv[i].k[k] = ks;
		    sv[i].len0[k] = tris[k]->side_length0[l];;
		}
		sorted(p) = YES;
	    	++i;
	    }
	}
	*n = i;
}	/* end set_surf_spring_vertex */

static void get_point_value_from(
	POINT *p,
	GLOBAL_POINT **point_set)
{
	int i;
	STATE *state = (STATE*)left_state(p);
	long gindex = Gindex(p);

	point_set[gindex]->gindex = gindex;
	for (i = 0; i < 3; ++i)
	{
	    point_set[gindex]->x[i] = Coords(p)[i];
	    point_set[gindex]->v[i] = p->vel[i];
	    point_set[gindex]->f[i] = p->force[i];
	    point_set[gindex]->impuls[i] = state->impulse[i];
	}
}	/* end get_point_value_from */
	
static void put_point_value_to(
	POINT *p,
	GLOBAL_POINT **point_set)
{
	int i;
	long gindex = Gindex(p);
	STATE *state = (STATE*)left_state(p);

	for (i = 0; i < 3; ++i)
	{
	    Coords(p)[i] = point_set[gindex]->x[i];
	    p->vel[i] = point_set[gindex]->v[i];
	    p->force[i] = point_set[gindex]->f[i];
	    state->impulse[i] = point_set[gindex]->impuls[i];
	}
}	/* end put_point_value_to */
	
extern void get_point_set_from(
	ELASTIC_SET *geom_set,
	GLOBAL_POINT **point_set)
{
	int i,ns,nc,nn;

	if (debugging("canopy"))
	    (void) printf("Entering get_point_set_from()\n");

	ns = geom_set->num_surfs;
	nc = geom_set->num_curves;
	nn = geom_set->num_nodes;
	for (i = 0; i < ns; ++i)
	    surf_get_point_set_from(geom_set->surfs[i],point_set);
	for (i = 0; i < nc; ++i)
	    curve_get_point_set_from(geom_set->curves[i],point_set);
	for (i = 0; i < nn; ++i)
	    node_get_point_set_from(geom_set->nodes[i],point_set);

	if (debugging("canopy"))
	    (void) printf("Leaving get_point_set_from()\n");
}	/* end  get_point_set_from */
	
extern void put_point_set_to(
	ELASTIC_SET *geom_set,
	GLOBAL_POINT **point_set)
{
	int i,ns,nc,nn;

	if (debugging("canopy"))
	    (void) printf("Entering put_point_set_to()\n");

	ns = geom_set->num_surfs;
	nc = geom_set->num_curves;
	nn = geom_set->num_nodes;
	for (i = 0; i < ns; ++i)
	    surf_put_point_set_to(geom_set->surfs[i],point_set);
	for (i = 0; i < nc; ++i)
	    curve_put_point_set_to(geom_set->curves[i],point_set);
	for (i = 0; i < nn; ++i)
	    node_put_point_set_to(geom_set->nodes[i],point_set);

	if (debugging("canopy"))
	    (void) printf("Leaving put_point_set_to()\n");
}	/* end  put_point_set_to */

static void surf_get_point_set_from(
	SURFACE *surf,
	GLOBAL_POINT **point_set)
{
	TRI *tri;
	POINT *p;
	int j;

	unsort_surf_point(surf);
	surf_tri_loop(surf,tri)
	{
	    for (j = 0; j < 3; ++j)
	    {
		p = Point_of_tri(tri)[j];
		if (sorted(p) || Boundary_point(p)) continue;
		get_point_value_from(p,point_set);
		sorted(p) = YES;
	    }
	}
}	/* end surf_get_point_set_from */

static void curve_get_point_set_from(
	CURVE *curve,
	GLOBAL_POINT **point_set)
{
	BOND *b;
	POINT *p;

	for (b = curve->first; b != curve->last; b = b->next)
	{
	    p = b->end;
	    get_point_value_from(p,point_set);
	}
}	/* end curve_get_point_set_from */

static void node_get_point_set_from(
	NODE *node,
	GLOBAL_POINT **point_set)
{
	POINT *p = node->posn;
	get_point_value_from(p,point_set);
}	/* end node_get_point_set_from */

static void surf_put_point_set_to(
	SURFACE *surf,
	GLOBAL_POINT **point_set)
{
	TRI *tri;
	POINT *p;
	int j;

	unsort_surf_point(surf);
	surf_tri_loop(surf,tri)
	{
	    for (j = 0; j < 3; ++j)
	    {
		p = Point_of_tri(tri)[j];
		if (sorted(p) || Boundary_point(p)) continue;
		put_point_value_to(p,point_set);
		sorted(p) = YES;
	    }
	}
}	/* end surf_put_point_set_to */

static void curve_put_point_set_to(
	CURVE *curve,
	GLOBAL_POINT **point_set)
{
	BOND *b;
	POINT *p;

	for (b = curve->first; b != curve->last; b = b->next)
	{
	    p = b->end;
	    put_point_value_to(p,point_set);
	}
}	/* end curve_put_point_set_to */

static void node_put_point_set_to(
	NODE *node,
	GLOBAL_POINT **point_set)
{
	POINT *p = node->posn;
	put_point_value_to(p,point_set);
}	/* end node_put_point_set_to */

extern void set_elastic_params(
	ELASTIC_SET *geom_set,
	double fr_dt)
{
	Front *front = geom_set->front;
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	double dt_tol;

	/* Set elastic set kinetic parameters */
        geom_set->ks = af_params->ks;
        geom_set->lambda_s = af_params->lambda_s;
        geom_set->m_s = af_params->m_s;
        geom_set->kl = af_params->kl;
        geom_set->lambda_l = af_params->lambda_l;
        geom_set->m_l = af_params->m_l;
        geom_set->kg = af_params->kg;
        geom_set->lambda_g = af_params->lambda_g;
        geom_set->m_g = af_params->m_g;

	/* Set elastic set time step */
        dt_tol = sqrt((af_params->m_s)/(af_params->ks))/10.0;
        if (af_params->m_l != 0.0 &&
            dt_tol > sqrt((af_params->m_l)/(af_params->kl))/10.0)
            dt_tol = sqrt((af_params->m_l)/(af_params->kl))/10.0;
        if (af_params->m_g != 0.0 &&
            dt_tol > sqrt((af_params->m_g)/(af_params->kg))/10.0)
            dt_tol = sqrt((af_params->m_g)/(af_params->kg))/10.0;
	geom_set->dt_tol = dt_tol;
}	/* end set_elastic_params */

extern void merge_global_point_set(
	GLOBAL_POINT **point_set,
	GLOBAL_POINT *gpoint_store,
	int num_gpoint)
{
	int i,k,gindex;
	for (i = 0; i < num_gpoint; ++i)
	{
	    gindex = gpoint_store[i].gindex;
	    *point_set[gindex] = gpoint_store[i];
	}
}	/* end merge_global_point_set */

extern void assembleParachuteSet(
	INTERFACE *intfc,
	ELASTIC_SET *geom_set)
{
	SURFACE **s;
	CURVE **c;
	NODE **n;
	int i,l,ns,nc,nn;
	SURFACE **surfs = geom_set->surfs;
	CURVE **curves = geom_set->curves;
	NODE **nodes = geom_set->nodes;
	int num_layers = 3;

	ns = nc = nn = 0;
	/* Assemble canopy surfaces */
	intfc_surface_loop(intfc,s)
	{
	    if (wave_type(*s) != ELASTIC_BOUNDARY)
		continue;
	    surfs[ns++] = *s;
	    surf_pos_curve_loop(*s,c)
	    {
	    	if (!pointer_in_list(*c,nc,(POINTER*)curves))
	    	{
		    curves[nc++] = *c;
		    if (!pointer_in_list((*c)->start,nn,(POINTER*)nodes))
		    	nodes[nn++] = (*c)->start;
		    if (!pointer_in_list((*c)->end,nn,(POINTER*)nodes))
		    	nodes[nn++] = (*c)->end;
	    	}
	    }
	    surf_neg_curve_loop(*s,c)
	    {
	    	if (!pointer_in_list(*c,nc,(POINTER*)curves))
	    	{
		    curves[nc++] = *c;
		    if (!pointer_in_list((*c)->start,nn,(POINTER*)nodes))
		    	nodes[nn++] = (*c)->start;
		    if (!pointer_in_list((*c)->end,nn,(POINTER*)nodes))
		    	nodes[nn++] = (*c)->end;
	    	}
	    }
	}

	/* Assemble curves and nodes */
	for (l = 0; l < num_layers; ++l)
	{
	    for (i = 0; i < nn; ++i)
	    {
	    	node_in_curve_loop(nodes[i],c)
	    	{
		    if (!pointer_in_list(*c,nc,(POINTER*)curves))
		    {
		    	curves[nc++] = *c;
		    	if (!pointer_in_list((*c)->start,nn,(POINTER*)nodes))
		    	    nodes[nn++] = (*c)->start;
		    	if (!pointer_in_list((*c)->end,nn,(POINTER*)nodes))
		    	    nodes[nn++] = (*c)->end;
		    }
	    	}
	    	node_out_curve_loop(nodes[i],c)
	    	{
		    if (!pointer_in_list(*c,nc,(POINTER*)curves))
		    {
		    	curves[nc++] = *c;
		    	if (!pointer_in_list((*c)->start,nn,(POINTER*)nodes))
		    	    nodes[nn++] = (*c)->start;
		    	if (!pointer_in_list((*c)->end,nn,(POINTER*)nodes))
		    	    nodes[nn++] = (*c)->end;
		    }
	    	}
	    }
	}
	geom_set->num_surfs = ns;
	geom_set->num_curves = nc;
	geom_set->num_nodes = nn;
	geom_set->num_verts = 0;
	for (i = 0; i < ns; ++i)
	    geom_set->num_verts += I_NumOfSurfInteriorPoints(surfs[i]);
	for (i = 0; i < nc; ++i)
	    geom_set->num_verts += I_NumOfCurveInteriorPoints(curves[i]);
	geom_set->num_verts += nn;
	for (i = 0; i < nn; ++i)
	    if (is_load_node(nodes[i]))
		geom_set->load_node = nodes[i];
}	/* end assembleParachuteSet */

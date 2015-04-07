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
*				fglobal.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	This file contains routines related to global index operations.
*/


#include <front/fdecs.h>

struct _GTRI
{
	long gtri[3];
};
typedef struct _GTRI GTRI;

struct _GBOND
{
	long gbond[2];
};
typedef struct _GBOND GBOND;


	/* LOCAL Function Declarations */
LOCAL boolean point_out_domain(POINT*,INTERFACE*,double*,double*);
LOCAL void merge_with_nb_surf_gindex(INTERFACE*,int,long*,int*,GTRI**);
LOCAL void merge_with_nb_curve_gindex(INTERFACE*,int,long*,int*,GBOND**);
LOCAL void exchange_surf_gindex(Front*);
LOCAL void exchange_curve_gindex(Front*);

LOCAL boolean point_out_domain(
	POINT *p,
	INTERFACE *intfc,
	double *L,
	double *U)
{
	int i,dim;
	dim = Dimension(intfc);
	for (i = 0; i < dim; ++i)
	{
	    if (rect_boundary_type(intfc,i,0) == SUBDOMAIN_BOUNDARY &&
		Coords(p)[i] < L[i])
		return YES;
	    if (rect_boundary_type(intfc,i,1) == SUBDOMAIN_BOUNDARY &&
		Coords(p)[i] >= U[i])
		return YES;
	}
	return NO;
}	/* end point_out_domain */

EXPORT	void set_point_gindex(
	Front *front)
{
	INTERFACE *intfc = front->interf;
	RECT_GRID *gr = front->rect_grid;
	SURFACE **s;
	CURVE **c;
	NODE **n;
	TRI *t;
	BOND *b;
	POINT *p;
	double *L = gr->L;
	double *U = gr->U;
	int i,j,iv,dim = gr->dim;
	long gindex = 0;
	long ilower,iupper,*n_dist;
	int num_nodes = pp_numnodes();
	int myid = pp_mynode();

	if (debugging("global_index"))
	    (void) printf("Entering set_point_gindex()\n");

        FT_VectorMemoryAlloc((POINTER*)&n_dist,num_nodes,sizeof(long));

	reset_sort_status(intfc);
	intfc_node_loop(intfc,n)
	{
	    p = (*n)->posn;
	    if (sorted(p)) continue;
	    sorted(p) = YES;
	    if (point_out_domain(p,intfc,L,U))
	    	continue;
	    Gindex(p) = gindex++;
	}
	intfc_curve_loop(intfc,c)
	{
	    for (b = (*c)->first; b != (*c)->last; b = b->next)
	    {
		p = b->end;
	    	if (sorted(p)) continue;
	    	sorted(p) = YES;
		if (point_out_domain(p,intfc,L,U))
		    continue;
	    	Gindex(p) = gindex++;
	    }
	}
	intfc_surface_loop(intfc,s)
	{
	    surf_tri_loop(*s,t)
	    {
		for (iv = 0; iv < 3; ++iv)
		{
		    p = Point_of_tri(t)[iv];
	    	    if (sorted(p)) continue;
	    	    sorted(p) = YES;
		    if (point_out_domain(p,intfc,L,U))
		    	continue;
		    Gindex(p) = gindex++;
		}
	    }
	}

	for (i = 0; i < num_nodes; ++i) n_dist[i] = 0;
	n_dist[myid] = gindex;
	pp_global_lmax(n_dist,num_nodes);
	ilower = 0;
        iupper = n_dist[0];
	for (i = 1; i <= myid; i++)
        {
            ilower += n_dist[i-1];
            iupper += n_dist[i];
        }
	intfc->max_point_gindex = iupper;
	pp_global_lmax(&intfc->max_point_gindex,1);
	reset_sort_status(intfc);
	intfc_node_loop(intfc,n)
	{
	    p = (*n)->posn;
	    if (sorted(p)) continue;
	    sorted(p) = YES;
	    if (point_out_domain(p,intfc,L,U))
	    	continue;
	    Gindex(p) += ilower;
	}
	intfc_curve_loop(intfc,c)
	{
	    for (b = (*c)->first; b != (*c)->last; b = b->next)
	    {
		p = b->end;
	    	if (sorted(p)) continue;
	    	sorted(p) = YES;
	    	if (point_out_domain(p,intfc,L,U))
	    	    continue;
	    	Gindex(p) += ilower;
	    }
	}
	intfc_surface_loop(intfc,s)
	{
	    surf_tri_loop(*s,t)
	    {
		for (iv = 0; iv < 3; ++iv)
		{
		    p = Point_of_tri(t)[iv];
	    	    if (sorted(p)) continue;
	    	    sorted(p) = YES;
	    	    if (point_out_domain(p,intfc,L,U))
	    	    	continue;
	    	    Gindex(p) += ilower;
		}
	    }
	}
	free_these(1,n_dist);

	static_mesh(intfc) = NO;
	scatter_front(front);
	if (debugging("global_index"))
	{
	    (void) printf("Call check_global_index()\n");
	    check_global_index(intfc);
	    (void) printf("Leaving set_point_gindex()\n");
	}
}	/* end set_point_gindex */

EXPORT	void set_surface_gindex(
	Front *front)
{
	RECT_GRID *gr = front->rect_grid;
	INTERFACE *intfc = front->interf;
	SURFACE **s;
	int num_nodes = pp_numnodes();
	int myid = pp_mynode();
	long ilower,iupper,*n_dist;
	int i,num_surfs;
	long max_surf_gindex;

        FT_VectorMemoryAlloc((POINTER*)&n_dist,num_nodes,sizeof(long));
	num_surfs = I_NumOfIntfcSurfaces(intfc);	
	for (i = 0; i < num_nodes; ++i) n_dist[i] = 0;
	n_dist[myid] = num_surfs;
	pp_global_lmax(n_dist,num_nodes);
	ilower = 0;
        iupper = n_dist[0];
	for (i = 1; i <= myid; i++)
        {
            ilower += n_dist[i-1];
            iupper += n_dist[i];
        }

	num_surfs = 0;
	intfc_surface_loop(intfc,s)
	{
	    Gindex(*s) = num_surfs + ilower;
	    if (intfc->max_surf_gindex < Gindex(*s))
		intfc->max_surf_gindex = Gindex(*s);
	    num_surfs++;
	}
	free_these(1,n_dist);

	pp_global_lmax(&intfc->max_surf_gindex,1);

	max_surf_gindex = intfc->max_surf_gindex + 1;
	while (max_surf_gindex > intfc->max_surf_gindex)
	{
	    max_surf_gindex = intfc->max_surf_gindex;
	    exchange_surf_gindex(front);
	}
}	/* end set_surface_gindex */

EXPORT	void set_curve_gindex(
	Front *front)
{
	RECT_GRID *gr = front->rect_grid;
	INTERFACE *intfc = front->interf;
	CURVE **c;
	int num_nodes = pp_numnodes();
	int myid = pp_mynode();
	long ilower,iupper,*n_dist;
	int i,num_curves;
	int dim = gr->dim;	
	long max_curve_gindex;

	if (debugging("curve_gindex"))
	    (void) printf("Entering set_curve_gindex()\n");
	num_curves = 0;
	intfc_curve_loop(intfc,c)
	{
	    if (dim == 3 && hsbdry_type(*c) == SUBDOMAIN_HSBDRY)
		continue;
	    num_curves++;
	}
        FT_VectorMemoryAlloc((POINTER*)&n_dist,num_nodes,sizeof(long));
	for (i = 0; i < num_nodes; ++i) n_dist[i] = 0;
	n_dist[myid] = num_curves;
	pp_global_lmax(n_dist,num_nodes);
	ilower = 0;
        iupper = n_dist[0];
	for (i = 1; i <= myid; i++)
        {
            ilower += n_dist[i-1];
            iupper += n_dist[i];
        }

	num_curves = 0;
	intfc_curve_loop(intfc,c)
	{
	    if (dim == 3 && hsbdry_type(*c) == SUBDOMAIN_HSBDRY)
		continue;
	    Gindex(*c) = num_curves + ilower;
	    if (intfc->max_curve_gindex < Gindex(*c))
		intfc->max_curve_gindex = Gindex(*c);
	    num_curves++;
	}
	free_these(1,n_dist);

	pp_global_lmax(&intfc->max_curve_gindex,1);

	max_curve_gindex = intfc->max_curve_gindex + 1;
	while (max_curve_gindex > intfc->max_curve_gindex)
	{
	    max_curve_gindex = intfc->max_curve_gindex;
	    exchange_curve_gindex(front);
	}
	if (debugging("curve_gindex"))
	    (void) printf("Leaving set_curve_gindex()\n");
}	/* end set_curve_gindex */

LOCAL void merge_with_nb_surf_gindex(
	INTERFACE *intfc,
	int num_surfs,
	long *surf_gindex,
	int *num_out_tris,
	GTRI **gtris)
{
	SURFACE **s;
	TRI *t;
	int i,j;
	boolean match_found;

	intfc_surface_loop(intfc,s)
	{
	    match_found = NO;
	    surf_tri_loop(*s,t)
	    {
		for (i = 0; i < num_surfs; ++i)
		{
		    for (j = 0; j < num_out_tris[i]; ++j)
		    {
			if (Gindex(Point_of_tri(t)[0]) == gtris[i][j].gtri[0] &&
			    Gindex(Point_of_tri(t)[1]) == gtris[i][j].gtri[1] &&
			    Gindex(Point_of_tri(t)[2]) == gtris[i][j].gtri[2])
			{
			    match_found = YES;
			    Gindex(*s) = min(Gindex(*s),surf_gindex[i]);
			    break;
			}
		    }
		    if (match_found) break;
		}
		if (match_found) break;
	    }
	}
}	/* end merge_with_nb_surf_gindex */
	
LOCAL void merge_with_nb_curve_gindex(
	INTERFACE *intfc,
	int num_curves,
	long *curve_gindex,
	int *num_out_bonds,
	GBOND **gbonds)
{
	CURVE **c;
	BOND *b;
	int i,j,dim = intfc->dim;
	boolean match_found;

	intfc_curve_loop(intfc,c)
	{
	    match_found = NO;
	    if (dim == 3 && hsbdry_type(*c) == SUBDOMAIN_HSBDRY)
		continue;
	    curve_bond_loop(*c,b)
	    {
		for (i = 0; i < num_curves; ++i)
		{
		    for (j = 0; j < num_out_bonds[i]; ++j)
		    {
			if (Gindex(b->start) == gbonds[i][j].gbond[0] &&
			    Gindex(b->end)   == gbonds[i][j].gbond[1]) 
			{
			    match_found = YES;
			    Gindex(*c) = min(Gindex(*c),curve_gindex[i]);
			    break;
			}
		    }
		    if (match_found) break;
		}
		if (match_found) break;
	    }
	}
}	/* end merge_with_nb_curve_gindex */

LOCAL void exchange_surf_gindex(
	Front *front)
{
	INTERFACE *intfc = front->interf;
	RECT_GRID *gr = front->rect_grid;
	double L[MAXD],U[MAXD];
	SURFACE **s;
	TRI *t;
	POINT *p;
	PP_GRID *pp_grid = front->pp_grid;
	int *G = pp_grid->gmax;
	int me[MAXD],him[MAXD];
	int dst_id,myid = pp_mynode();
	int num_nodes = pp_numnodes();
	int *num_out_tris;
	int num_surfs;
	GTRI **gtris;
	long *surf_gindex,*g_list;
	int i,j,k,imax,dim = gr->dim;

	num_surfs = I_NumOfIntfcSurfaces(intfc);	
	FT_VectorMemoryAlloc((POINTER*)&num_out_tris,num_surfs,sizeof(int));
	FT_VectorMemoryAlloc((POINTER*)&surf_gindex,num_surfs,sizeof(long));
	FT_VectorMemoryAlloc((POINTER*)&gtris,num_surfs,sizeof(GTRI*));

	for (i = 0; i < dim; i++)
	{
	    L[i] = gr->L[i]+(gr->lbuf[i]*gr->h[i]);
	    U[i] = gr->U[i]-(gr->ubuf[i]*gr->h[i]);
	}
	num_surfs = 0;
	intfc_surface_loop(intfc,s)
	{
	    num_out_tris[num_surfs] = 0;
	    surf_tri_loop(*s,t)
	    {
		for (i = 0; i < 3; ++i)
                {
                    p = Point_of_tri(t)[i];
                    if (point_out_domain(p,intfc,L,U))
                    {
                        num_out_tris[num_surfs]++;
                        break;
                    }
                }
	    }
	    num_surfs++;
	}
	for (i = 0; i < num_surfs; ++i)
	    FT_VectorMemoryAlloc((POINTER*)&gtris[i],num_out_tris[i],
				sizeof(GTRI));

	num_surfs = 0;
	intfc_surface_loop(intfc,s)
	{
	    surf_gindex[num_surfs] = Gindex(*s);
	    i = 0;
	    surf_tri_loop(*s,t)
	    {
		for (j = 0; j < 3; j++)
                {
                    p = Point_of_tri(t)[j];
                    if (point_out_domain(p,intfc,L,U))
                        break;
                }
                if (j < 3)
                {
                    gtris[num_surfs][i].gtri[0] = Gindex(Point_of_tri(t)[0]);
                    gtris[num_surfs][i].gtri[1] = Gindex(Point_of_tri(t)[1]);
                    gtris[num_surfs][i].gtri[2] = Gindex(Point_of_tri(t)[2]);
                    i++;
                }
	    }
	    num_out_tris[num_surfs] = i;
	    num_surfs++;
	}

	find_Cartesian_coordinates(myid,pp_grid,me);
	for (i = 0; i < dim; ++i)
	{
	    if (G[i] == 1) continue;
	    for (j = 0; j < 2; ++j)
	    {
	    	if (rect_boundary_type(intfc,i,j) != SUBDOMAIN_BOUNDARY)
		    continue;
		for (k = 0; k < dim; ++k)  him[k] = me[k];
		him[i] = me[i] + 2*j - 1;
                him[i] = (him[i]+G[i])%G[i];
		dst_id = domain_id(him,G,dim);
		pp_send(50,&num_surfs,sizeof(int),dst_id);
		for (k = 0; k < num_surfs; ++k)
		{
		    pp_send(51,&surf_gindex[k],sizeof(long),dst_id);
		    pp_send(52,&num_out_tris[k],sizeof(int),dst_id);
		    pp_send(53,gtris[k],num_out_tris[k]*sizeof(GTRI),dst_id);
		}
	    }
	}
	for (i = 0; i < num_surfs; ++i)
	    free_these(1,gtris[i]);
	free_these(3,surf_gindex,num_out_tris,gtris);

	for (i = 0; i < dim; ++i)
	{
	    if (G[i] == 1) continue;
	    for (j = 1; j >= 0; --j)
	    {
	    	if (rect_boundary_type(intfc,i,j) != SUBDOMAIN_BOUNDARY)
		    continue;
		for (k = 0; k < dim; ++k)  him[k] = me[k];
		him[i] = me[i] + 2*j - 1;
                him[i] = (him[i]+G[i])%G[i];
		dst_id = domain_id(him,G,dim);
		pp_recv(50,dst_id,&num_surfs,sizeof(int));
		FT_VectorMemoryAlloc((POINTER*)&num_out_tris,num_surfs,
						sizeof(int));
		FT_VectorMemoryAlloc((POINTER*)&surf_gindex,num_surfs,
						sizeof(long));
		FT_VectorMemoryAlloc((POINTER*)&gtris,num_surfs,
						sizeof(GTRI*));
		for (k = 0; k < num_surfs; ++k)
		{
		    pp_recv(51,dst_id,&surf_gindex[k],sizeof(long));
		    pp_recv(52,dst_id,&num_out_tris[k],sizeof(int));
	    	    FT_VectorMemoryAlloc((POINTER*)&gtris[k],num_out_tris[k],
				sizeof(GTRI));
		    pp_recv(53,dst_id,gtris[k],num_out_tris[k]*sizeof(GTRI));
		}
		merge_with_nb_surf_gindex(intfc,num_surfs,surf_gindex,
				num_out_tris,gtris);
		for (k = 0; k < num_surfs; ++k)
		    free_these(1,gtris[k]);
		free_these(3,surf_gindex,num_out_tris,gtris);
	    }
	}
	intfc->max_surf_gindex = 0;
	intfc_surface_loop(intfc,s)
	{
	    if (intfc->max_surf_gindex < Gindex(*s))
		intfc->max_surf_gindex = Gindex(*s);
	}
	pp_global_lmax(&intfc->max_surf_gindex,1);
	FT_VectorMemoryAlloc((POINTER*)&g_list,intfc->max_surf_gindex+1,
				sizeof(long));
	for (i = 0; i <= intfc->max_surf_gindex; ++i) g_list[i] = -1;
	intfc_surface_loop(intfc,s)
	{
	    i = (int)Gindex(*s);
	    g_list[i] = Gindex(*s);
	}
	pp_global_lmax(g_list,intfc->max_surf_gindex+1);
	if (debugging("surf_gindex"))
	{
	    printf("Global index after merging:\n");
	    for (i = 0; i <= intfc->max_surf_gindex; ++i) 
	    	(void) printf("g_list[%d] = %ld\n",i,g_list[i]);
	}
	/* Shift global index to the lower end */
	for (i = 0; i <= intfc->max_surf_gindex; ++i) 
	{
	    if (g_list[i] == -1)
	    {
		for (j = 0; j <= intfc->max_surf_gindex-i; ++j)
		    if (g_list[i+j] != -1) break;
		if (j > intfc->max_surf_gindex) 
		    break;
		else
		{
		    for (k = i; k <= intfc->max_surf_gindex; ++k)
		    {
			if (k+j > intfc->max_surf_gindex)
			    g_list[k] = -1;
			else
			    g_list[k] = g_list[k+j];
		    }
		}
	    }
	}
	if (debugging("surf_gindex"))
	{
	    printf("Global index after shifting:\n");
	    for (i = 0; i <= intfc->max_surf_gindex; ++i) 
	    	(void) printf("g_list[%d] = %ld\n",i,g_list[i]);
	}
	intfc_surface_loop(intfc,s)
	{
	    for (i = 0; i <= intfc->max_surf_gindex; ++i)
		if (Gindex(*s) == g_list[i])
		    Gindex(*s) = i;
	}
	for (i = 0; i <= intfc->max_surf_gindex; ++i) 
	    if (g_list[i] != -1) 
	{
	    g_list[i] = i;
	    imax = i;
	}
	if (debugging("surf_gindex"))
	{
	    printf("Global index after resetting:\n");
	    for (i = 0; i <= intfc->max_surf_gindex; ++i) 
	    	(void) printf("g_list[%d] = %ld\n",i,g_list[i]);
	}
	intfc->max_surf_gindex = imax;
	if (debugging("surf_gindex"))
	{
	    (void) printf("Leaving exchange_surf_gindex():\n");
	    (void) printf("max_surf_gindex = %ld\n",intfc->max_surf_gindex);
	}
}	/* end exchange_surf_gindex */

LOCAL void exchange_curve_gindex(
	Front *front)
{
	INTERFACE *intfc = front->interf;
	RECT_GRID *gr = front->rect_grid;
	double L[MAXD],U[MAXD];
	CURVE **c;
	BOND *b;
	POINT *p;
	PP_GRID *pp_grid = front->pp_grid;
	int *G = pp_grid->gmax;
	int me[MAXD],him[MAXD];
	int dst_id,myid = pp_mynode();
	int num_nodes = pp_numnodes();
	int *num_out_bonds;
	int num_curves;
	GBOND **gbonds;
	long *curve_gindex,*g_list;
	int i,j,k,imax,dim = gr->dim;

	num_curves = I_NumOfIntfcCurves(intfc);	
	FT_VectorMemoryAlloc((POINTER*)&num_out_bonds,num_curves,sizeof(int));
	FT_VectorMemoryAlloc((POINTER*)&curve_gindex,num_curves,sizeof(long));
	FT_VectorMemoryAlloc((POINTER*)&gbonds,num_curves,sizeof(GBOND*));

	for (i = 0; i < dim; i++)
	{
	    L[i] = gr->L[i]+(gr->lbuf[i]*gr->h[i]);
	    U[i] = gr->U[i]-(gr->ubuf[i]*gr->h[i]);
	}
	num_curves = 0;
	intfc_curve_loop(intfc,c)
	{
	    if (dim == 3 && hsbdry_type(*c) == SUBDOMAIN_HSBDRY)
		continue;
	    num_out_bonds[num_curves] = 0;
	    curve_bond_loop(*c,b)
	    {
		if (point_out_domain(b->start,intfc,L,U) ||
                    point_out_domain(b->end,intfc,L,U))
                {
                    num_out_bonds[num_curves]++;
                }
	    }
	    num_curves++;
	}
	for (i = 0; i < num_curves; ++i)
	    FT_VectorMemoryAlloc((POINTER*)&gbonds[i],num_out_bonds[i],
				sizeof(GBOND));

	num_curves = 0;
	intfc_curve_loop(intfc,c)
	{
	    if (dim == 3 && hsbdry_type(*c) == SUBDOMAIN_HSBDRY)
		continue;
	    curve_gindex[num_curves] = Gindex(*c);
	    i = 0;
	    curve_bond_loop(*c,b)
	    {
		if (point_out_domain(b->start,intfc,L,U) ||
                    point_out_domain(b->end,intfc,L,U))
                {
                    num_out_bonds[num_curves]++;
		    gbonds[num_curves][i].gbond[0] = Gindex(b->start);
		    gbonds[num_curves][i].gbond[1] = Gindex(b->end);
		    i++;
                }
	    }
	    num_out_bonds[num_curves] = i;
	    num_curves++;
	}

	find_Cartesian_coordinates(myid,pp_grid,me);
	for (i = 0; i < dim; ++i)
	{
	    if (G[i] == 1) continue;
	    for (j = 0; j < 2; ++j)
	    {
	    	if (rect_boundary_type(intfc,i,j) != SUBDOMAIN_BOUNDARY)
		    continue;
		for (k = 0; k < dim; ++k)  him[k] = me[k];
		him[i] = me[i] + 2*j - 1;
                him[i] = (him[i]+G[i])%G[i];
		dst_id = domain_id(him,G,dim);
		pp_send(50,&num_curves,sizeof(int),dst_id);
		for (k = 0; k < num_curves; ++k)
		{
		    pp_send(51,&curve_gindex[k],sizeof(long),dst_id);
		    pp_send(52,&num_out_bonds[k],sizeof(int),dst_id);
		    pp_send(53,gbonds[k],num_out_bonds[k]*sizeof(GBOND),dst_id);
		}
	    }
	}
	for (i = 0; i < num_curves; ++i)
	    free_these(1,gbonds[i]);
	free_these(3,curve_gindex,num_out_bonds,gbonds);

	for (i = 0; i < dim; ++i)
	{
	    if (G[i] == 1) continue;
	    for (j = 1; j >= 0; --j)
	    {
	    	if (rect_boundary_type(intfc,i,j) != SUBDOMAIN_BOUNDARY)
		    continue;
		for (k = 0; k < dim; ++k)  him[k] = me[k];
		him[i] = me[i] + 2*j - 1;
                him[i] = (him[i]+G[i])%G[i];
		dst_id = domain_id(him,G,dim);
		pp_recv(50,dst_id,&num_curves,sizeof(int));
		FT_VectorMemoryAlloc((POINTER*)&num_out_bonds,num_curves,
						sizeof(int));
		FT_VectorMemoryAlloc((POINTER*)&curve_gindex,num_curves,
						sizeof(long));
		FT_VectorMemoryAlloc((POINTER*)&gbonds,num_curves,
						sizeof(GBOND*));
		for (k = 0; k < num_curves; ++k)
		{
		    pp_recv(51,dst_id,&curve_gindex[k],sizeof(long));
		    pp_recv(52,dst_id,&num_out_bonds[k],sizeof(int));
	    	    FT_VectorMemoryAlloc((POINTER*)&gbonds[k],num_out_bonds[k],
				sizeof(GBOND));
		    pp_recv(53,dst_id,gbonds[k],num_out_bonds[k]*sizeof(GBOND));
		}
		merge_with_nb_curve_gindex(intfc,num_curves,curve_gindex,
				num_out_bonds,gbonds);
		for (k = 0; k < num_curves; ++k)
		    free_these(1,gbonds[k]);
		free_these(3,curve_gindex,num_out_bonds,gbonds);
	    }
	}
	intfc->max_curve_gindex = 0;
	intfc_curve_loop(intfc,c)
	{
	    if (dim == 3 && hsbdry_type(*c) == SUBDOMAIN_HSBDRY)
		continue;
	    if (intfc->max_curve_gindex < Gindex(*c))
		intfc->max_curve_gindex = Gindex(*c);
	}
	pp_global_lmax(&intfc->max_curve_gindex,1);
	FT_VectorMemoryAlloc((POINTER*)&g_list,intfc->max_curve_gindex+1,
				sizeof(long));
	for (i = 0; i <= intfc->max_curve_gindex; ++i) g_list[i] = -1;
	intfc_curve_loop(intfc,c)
	{
	    i = (int)Gindex(*c);
	    g_list[i] = Gindex(*c);
	}
	pp_global_lmax(g_list,intfc->max_curve_gindex+1);
	if (debugging("curve_gindex"))
	{
	    printf("Global index after merging:\n");
	    for (i = 0; i <= intfc->max_curve_gindex; ++i) 
	    	(void) printf("g_list[%d] = %ld\n",i,g_list[i]);
	}
	/* Shift global index to the lower end */
	for (i = 0; i <= intfc->max_curve_gindex; ++i) 
	{
	    if (g_list[i] == -1)
	    {
		for (j = 0; j <= intfc->max_curve_gindex-i; ++j)
		    if (g_list[i+j] != -1) break;
		if (j > intfc->max_curve_gindex) 
		    break;
		else
		{
		    for (k = i; k <= intfc->max_curve_gindex; ++k)
		    {
			if (k+j > intfc->max_curve_gindex)
			    g_list[k] = -1;
			else
			    g_list[k] = g_list[k+j];
		    }
		}
	    }
	}
	if (debugging("curve_gindex"))
	{
	    printf("Global index after shifting:\n");
	    for (i = 0; i <= intfc->max_curve_gindex; ++i) 
	    	(void) printf("g_list[%d] = %ld\n",i,g_list[i]);
	}
	intfc_curve_loop(intfc,c)
	{
	    if (dim == 3 && hsbdry_type(*c) == SUBDOMAIN_HSBDRY)
		continue;
	    for (i = 0; i <= intfc->max_curve_gindex; ++i)
		if (Gindex(*c) == g_list[i])
		    Gindex(*c) = i;
	}
	for (i = 0; i <= intfc->max_curve_gindex; ++i) 
	    if (g_list[i] != -1) 
	{
	    g_list[i] = i;
	    imax = i;
	}
	if (debugging("curve_gindex"))
	{
	    printf("Global index after resetting:\n");
	    for (i = 0; i <= intfc->max_curve_gindex; ++i) 
	    	(void) printf("g_list[%d] = %ld\n",i,g_list[i]);
	}
	intfc->max_curve_gindex = imax;
	if (debugging("curve_gindex"))
	{
	    (void) printf("Leaving exchange_curve_gindex():\n");
	    (void) printf("max_curve_gindex = %ld\n",intfc->max_curve_gindex);
	}
}	/* end exchange_curve_gindex */

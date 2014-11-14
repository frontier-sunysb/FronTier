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
*			fscat3d3.c
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*/


#define DEBUG_STRING    "fscatter"
#include <front/fdecs.h>

struct _POINT_LIST {
	POINT              *p;
	HYPER_SURF         *hs;
	HYPER_SURF_ELEMENT *hse;
	struct _POINT_LIST *prev, *next;
};
typedef struct _POINT_LIST POINT_LIST;

struct _POINT_LIST_STORE {
	POINT_LIST *pl;
	int        len;
};
typedef struct _POINT_LIST_STORE POINT_LIST_STORE;

	/* LOCAL Function Declarations */
LOCAL	POINT_LIST *set_point_list(TRI**,int,HYPER_SURF*,
					POINT_LIST_STORE*);
LOCAL	boolean	add_matching_pt_to_hash_table(TRI**,TRI**,int,int,SURFACE*,
					      SURFACE*,P_LINK*,int);
LOCAL	boolean	buffer_extension3d3(INTERFACE*,INTERFACE*,int,int,boolean);
LOCAL	boolean	match_tris_at_subdomain_bdry(SURFACE*,SURFACE*,TRI**,TRI**,
	                                     int,int);
LOCAL	boolean	match_two_tris(TRI*,TRI*);
LOCAL	boolean	tri_cross_line(TRI*,double,int);
LOCAL	int	append_adj_intfc_to_buffer3(INTERFACE*,INTERFACE*,
					   RECT_GRID*,int,int);
LOCAL	int	append_buffer_surface3(SURFACE*,SURFACE*,RECT_GRID*,int,int,
				      P_LINK*,int);
LOCAL	void	clip_intfc_at_grid_bdry1(INTERFACE*);
LOCAL	boolean	tri_bond_cross_line(TRI*,double,int);
LOCAL	boolean	tri_bond_cross_test(TRI*,double,int);
LOCAL	void	synchronize_tris_at_subdomain_bdry(TRI**,TRI**,int,P_LINK*,int);
LOCAL 	INTERFACE *cut_intfc_to_wave_type(INTERFACE*,int);
LOCAL 	void 	delete_surface_set(SURFACE*);
LOCAL   boolean append_other_curves3(INTERFACE*,INTERFACE*,P_LINK*,int);
LOCAL   boolean bond_match3(BOND*,BOND*);
LOCAL   void merge_overlap_nodes(INTERFACE*);
LOCAL	void print_unmatched_tris(TRI**,TRI**,int,int);
LOCAL 	void exchange_intfc_extra(INTERFACE*);

LOCAL	double	tol1[MAXD]; /*TOLERANCE*/

#define MAX_SUBDOMAIN_TRIS      3000

/*	New LOCAL functions	*/
LOCAL  void    set_floating_point_tolerance3(double*);

/*ARGSUSED*/

EXPORT	boolean f_intfc_communication3d3(
	Front		*fr)
{
	INTERFACE    *intfc = fr->interf;
	INTERFACE    *adj_intfc[2], *sav_intfc, *buf_intfc;
	PP_GRID	     *pp_grid = fr->pp_grid;
	RECT_GRID    *gr = fr->rect_grid;
	boolean	     sav_copy;
	boolean      status = FUNCTION_SUCCEEDED;
	double        *U = gr->U, *L = gr->L;
	double        *nor, p[3];
	double        T;
	int	     me[MAXD], him[MAXD];
	int	     myid, dst_id;
	int	     *G;
	int	     i, j, k;
	int	     dim = intfc->dim;
	static double nors[] = {  1.0,  0.0,  0.0,
			         0.0,  1.0,  0.0,
			         0.0,  0.0,  1.0,
			        -1.0,  0.0,  0.0,
			         0.0, -1.0,  0.0,
			         0.0,  0.0, -1.0};

	DEBUG_ENTER(f_intfc_communication3d3)

	set_floating_point_tolerance3(fr->rect_grid->h);
	sav_copy = copy_intfc_states();
	sav_intfc = current_interface();
	set_copy_intfc_states(YES);

	myid = pp_mynode();
	G = pp_grid->gmax;
	find_Cartesian_coordinates(myid,pp_grid,me);

	if (DEBUG)
	{
	    (void) printf("myid = %d, ",myid);
	    print_int_vector("me = ",me,dim,"\n");
	    print_PP_GRID_structure(pp_grid);
	    (void) printf("Input interface:\n");
	    print_interface(intfc);
	}

	/* Extend interface in three directions */

        construct_reflect_bdry(fr);
	clip_intfc_at_grid_bdry1(intfc);

	for (i = 0; i < dim; ++i)
	{
	    adj_intfc[0] = adj_intfc[1] = NULL;
	    for (j = 0; j < 2; ++j)
	    {
	    	pp_gsync();

		for (k = 0; k < dim; ++k)
		    him[k] = me[k];

		if (rect_boundary_type(intfc,i,j) == SUBDOMAIN_BOUNDARY)
		{
		    him[i] = (me[i] + 2*j - 1 + G[i])%G[i];
		    dst_id = domain_id(him,G,dim);
		    buf_intfc = cut_buf_interface1(intfc,i,j,me,him);
		    if ((j == 0) && (me[i] == 0))
		    {
			T = gr->GU[i] - gr->GL[i];
	                shift_interface(buf_intfc,T,i);
		    }
		    else if ((j == 1) && (me[i] == (G[i]-1)))
		    {
			T = gr->GL[i] - gr->GU[i];
	                shift_interface(buf_intfc,T,i);
		    }
		    if (dst_id != myid)
		    {
		    	send_interface(buf_intfc,dst_id);
		        (void) delete_interface(buf_intfc);
		    }
		    else
		        adj_intfc[(j+1)%2] = buf_intfc;
		}
                else if (rect_boundary_type(intfc,i,j) == REFLECTION_BOUNDARY)
                {
                    adj_intfc[j] = NULL;
                    set_current_interface(intfc);
                }
		if (rect_boundary_type(intfc,i,(j+1)%2) == SUBDOMAIN_BOUNDARY)
		{
		    him[i] = (me[i] - 2*j + 1 + G[i])%G[i];
		    dst_id = domain_id(him,G,dim);
		    if (dst_id != myid)
		    {
			adj_intfc[(j+1)%2] = receive_interface(dst_id);
		    }
		}

	    }
	    for (j = 0; j < 2; ++j)
	    {
		
		status = FUNCTION_SUCCEEDED;
		if (adj_intfc[j] != NULL)
		{
		    status = buffer_extension3d3(intfc,adj_intfc[j],i,j,status);
                   
		    if (!status)
                    {
		        (void) printf("WARNING in f_intfc_communication3d3 "
				      "buffer_extension3d3 failed for "
                                        "i = %d, j = %d \n", i, j);
                    }

		    (void) delete_interface(adj_intfc[j]);
		    set_current_interface(intfc);
		}    /*if (adj_intfc[j] != NULL) */
		
		/*surface in this direction can not mathch, return */
		status = pp_min_status(status);
		if(!status)
		    goto comm_exit;

	    }	
	}
	exchange_intfc_extra(intfc);

	if (status == FUNCTION_SUCCEEDED)
	{
	    sep_common_pt_for_open_bdry(intfc);

	    install_subdomain_bdry_curves(intfc);
	    if (DEBUG)
	    {
	    	(void) printf("Final intfc:\n");
	    	print_interface(intfc);
	    }
	}

comm_exit:
	set_copy_intfc_states(sav_copy);
	set_current_interface(sav_intfc);
	DEBUG_LEAVE(f_intfc_communication3d3)
	return status;
}	/*end f_intfc_communication3d3 */

LOCAL boolean buffer_extension3d3(
	INTERFACE	*intfc,
	INTERFACE	*adj_intfc,
	int		dir,
	int		nb,
	boolean		status)
{
	RECT_GRID	*gr = computational_grid(intfc);

	DEBUG_ENTER(buffer_extension3d3)

	set_current_interface(intfc);

		/* Patch tris from adj_intfc to intfc */

	if (!append_adj_intfc_to_buffer3(intfc,adj_intfc,gr,dir,nb))
	{
	    status = FUNCTION_FAILED;
	    (void) printf("WARNING in buffer_extension3d3(), "
	                  "append_adj_intfc_to_buffer3() failed\n");
	}

	DEBUG_LEAVE(buffer_extension3d3)
	return status;
}	/*end buffer_extension3d3*/

LOCAL 	int append_adj_intfc_to_buffer3(
	INTERFACE	*intfc,		/* my interface 	       */
	INTERFACE	*adj_intfc,	/* received interface	       */
	RECT_GRID	*grid,		/* Rectangular grid for region */
	int		dir,
	int		nb)
{
	INTERFACE	*cur_intfc;
	SURFACE		**s, **as;
	CURVE		**ac;
	int		p_size;		/*Size of space allocated for p_table*/
	static P_LINK	*p_table = NULL;/* Table of matching points on intfc
					 * and adj_intfc*/
	static int      len_p_table = 0;
	boolean      	corr_surf_found;

	DEBUG_ENTER(append_adj_intfc_to_buffer3)
	if (debugging("append_curves"))
        {
            (void) printf("Entering append_adj_intfc_to_buffer1()\n");
            (void) printf("Number of curves = %d\n",I_NumOfIntfcCurves(intfc));
        }

	cur_intfc = current_interface();
	set_current_interface(intfc);

	p_size = 4*(adj_intfc->num_points) + 1;
	if (len_p_table < p_size)
	{
	    len_p_table = p_size;
	    if (p_table != NULL)
		free(p_table);
	    uni_array(&p_table,len_p_table,sizeof(P_LINK));
	}
	
	reset_hash_table(p_table,p_size);
	
	/* Begin patching adj_intfc to current interface */

	/* append curves not on surfaces */
	append_other_curves3(intfc,adj_intfc,p_table,p_size);

	/* append surfaces and attached curves */
	for (as = adj_intfc->surfaces; as && *as; ++as)
	{
	    corr_surf_found = NO;
	    for (s = intfc->surfaces; s && *s; ++s)
	    {
		/*  assume one neg and pos comp gives ONLY ONE surf */
		if (surfaces_matched(*as,*s))
		{
		    corr_surf_found = YES;
		    if (!append_buffer_surface3(*s,*as,grid,dir,nb,p_table,
						   p_size))
		    {
			set_current_interface(cur_intfc);
			
			(void) printf("WARNING in "
			              "append_adj_intfc_to_buffer3(), "
			              "append surface failed\n");
			
			DEBUG_LEAVE(append_adj_intfc_to_buffer3)
			return NO;
		    }
		}
	    }
	    if (!corr_surf_found && as && *as)
	    {
	    	SURFACE *surf;
		if (as==NULL)
		    continue;
		surf = copy_buffer_surface(*as,p_table,p_size);
		Hyper_surf_index(surf) = Hyper_surf_index((*as));
	    }
	}
	for (ac = adj_intfc->curves; ac && *ac; ++ac)
	{
	    matching_curve(*ac,p_table,p_size); 
	}
	merge_curves(intfc,adj_intfc);
	reset_intfc_num_points(intfc);
	
	set_current_interface(cur_intfc);
	
	if (debugging("append_curves"))
        {
            (void) printf("Leaving append_adj_intfc_to_buffer3()\n");
            (void) printf("Number of curves = %d\n",I_NumOfIntfcCurves(intfc));
        }
	DEBUG_LEAVE(append_adj_intfc_to_buffer3)
	return YES;
}		/*end append_adj_intfc_to_buffer3*/

LOCAL void set_floating_point_tolerance3(
	double		*h)
{
	double eps;
	static const double mfac = 10.0;/*TOLERANCE*/
	static const double gfac = 0.00004;/*TOLERANCE*/
	int   i;

	eps = mfac*MACH_EPS;
	for (i = 0; i < 3; ++i)
	{
	    tol1[i] = gfac*h[i];/*TOLERANCE*/
	    tol1[i] = max(tol1[i],eps);
	}
}		/*end set_floating_point_tolerance3*/

LOCAL int append_buffer_surface3(
	SURFACE		*surf,
	SURFACE		*adj_surf,
	RECT_GRID	*grid,
	int		dir,
	int		nb,
	P_LINK		*p_table,
	int		p_size)
{
	SURFACE    *new_adj_surf;
	TRI	   *tri;
	CURVE	   **c;
	int	   ns, na, new_na;
	double	   crx_coord;
	static TRI **tris_s = NULL, **tris_a = NULL;
	static int len_tris_s = 0, len_tris_a = 0;

	DEBUG_ENTER(append_buffer_surface3)
	
	if (DEBUG)
	{
	    (void) printf("dir = %d,nb = %d\n",dir,nb);
	    (void) printf("My surface\n");
	    print_surface(surf);
	    (void) printf("Buffer surface\n");
	    print_surface(adj_surf);
	}
	
	if (len_tris_s < surf->num_tri)
	{
	    len_tris_s = surf->num_tri;
	    if (tris_s != NULL)
		free(tris_s);
	    uni_array(&tris_s,len_tris_s,sizeof(TRI *));
	}
	if (len_tris_a < adj_surf->num_tri)
	{
	    len_tris_a = adj_surf->num_tri;
	    if (tris_a != NULL)
		free(tris_a);
	    uni_array(&tris_a,len_tris_a,sizeof(TRI *));
	}

	crx_coord = (nb == 0) ? grid->L[dir] : grid->U[dir];

	ns = na = 0;
	for (tri=first_tri(surf); !at_end_of_tri_list(tri,surf); tri=tri->next)
	{
	    if (tri_cross_line(tri,crx_coord,dir) == YES || 
		tri_bond_cross_test(tri,crx_coord,dir) == YES)
	    {
		tris_s[ns++] = tri;
	    }
	}

	for (tri = first_tri(adj_surf); !at_end_of_tri_list(tri,adj_surf); 
	     tri = tri->next)
	{
	    if (tri_cross_line(tri,crx_coord,dir) == YES ||
		tri_bond_cross_test(tri,crx_coord,dir) == YES)
	    {
		tris_a[na++] = tri;
	    }
	}
	if (ns != na)
	{
	    print_unmatched_tris(tris_s,tris_a,ns,na);
	    clean_up(ERROR);
	}

	/* Add matching points to the hashing table p_table */

	if (add_matching_pt_to_hash_table(tris_s,tris_a,ns,na,surf,
		      adj_surf,p_table,p_size))
	{
	    new_adj_surf = copy_buffer_surface(adj_surf,p_table,p_size);
	    adj_surf = new_adj_surf;
	}
	else 
	{
	    printf("dir = %d  nb = %d\n",dir,nb);
	    gview_plot_surface("surf_s",surf);
	    gview_plot_surface("surf_a",adj_surf);
	    gview_plot_tri_list("tris_s",tris_s,ns);
	    gview_plot_tri_list("tris_a",tris_a,na);
	    return NO;
	}

	synchronize_tris_at_subdomain_bdry(tris_a,tris_s,ns,p_table,p_size);
	
	na = 0;
	for (tri = first_tri(adj_surf); !at_end_of_tri_list(tri,adj_surf);
	     tri = tri->next)
	{
	    if (tri_cross_line(tri,crx_coord,dir) == YES || 
		tri_bond_cross_test(tri,crx_coord,dir) == YES)
	    {
	        tris_a[na++] = tri;
	    }
	}

	/*adjoin adj_surf tri list to surf tri list */
	last_tri(surf)->next = first_tri(adj_surf);
	first_tri(adj_surf)->prev = last_tri(surf);
	link_tri_list_to_surface(first_tri(surf),last_tri(adj_surf),surf);
	
	/* BEGIN curves in adj_surf should be added to surf. */
	for (c = adj_surf->pos_curves; c && *c; c++)
	    if (!delete_from_pointers(adj_surf, &(*c)->pos_surfaces))
	    {
	        printf("ERROR: in append_buffer_surface3, "
		       "adj_surf and pos_curves are not paired.\n");
		clean_up(ERROR);
	    }
	    else
	        install_curve_in_surface_bdry(surf, *c, POSITIVE_ORIENTATION);

	for (c = adj_surf->neg_curves; c && *c; c++)
	    if (!delete_from_pointers(adj_surf, &(*c)->neg_surfaces))
	    {
	        printf("ERROR: in append_buffer_surface3, "
		       "adj_surf and neg_curves are not paired.\n");
		clean_up(ERROR);
	    }
	    else
	        install_curve_in_surface_bdry(surf, *c, NEGATIVE_ORIENTATION);

	adj_surf->pos_curves = adj_surf->neg_curves = NULL;
	(void) delete_surface(adj_surf);
	adj_surf = NULL;

	/* average_btris*/
	if (!match_tris_at_subdomain_bdry(surf,adj_surf,tris_s,tris_a,ns,na))
	{
	    (void) printf("WARNING in append_buffer_surface3(), "
	                  "no match of tris at subdomain\n");
	    (void) printf("dir = %d, nb = %d\n",dir,nb);
	    DEBUG_LEAVE(append_buffer_surface3)
	    return NO;
	}

	DEBUG_LEAVE(append_buffer_surface3)
	return YES;
}		/*end append_buffer_surface3*/


LOCAL	void	synchronize_tris_at_subdomain_bdry(
	TRI    **tris_a,
	TRI    **tris_s,
	int    nt,
	P_LINK *p_table,
	int    p_size)
{
	POINT **ps, **pa, *p0, *p1, *p2;
	TRI   *ts, *ta;
	int   i, j, id, idp, idn;

	for (i = 0; i < nt; ++i)
	{
	    ts = tris_s[i];
	    ps = Point_of_tri(ts);
	    for (j = 0; j < nt; ++j)
	    {
		ta = tris_a[j];
		pa = Point_of_tri(ta);
		for (id = 0; id < 3; ++id)
		{
		    p0 = (POINT*) find_from_hash_table((POINTER)pa[id],
						       p_table,p_size);
		    if (p0 == ps[0])
		    {
		        idn = Next_m3(id);
		        p1 = (POINT*) find_from_hash_table((POINTER)pa[idn],
						           p_table,p_size);
		        if (p1 == ps[1])
		        {
		            idp = Prev_m3(id);
		            p2 = (POINT*) find_from_hash_table((POINTER)pa[idp],
						               p_table,p_size);
			    if (p2 == ps[2])
			    {
			        rotate_triangle(ts,(3-id)%3);
				
				if(debugging("ts_tst"))
				{
				    printf("#sync tris af\n");
				    print_tri(ts, ts->surf->interface);
				    remove_from_debug("ts_tst");
				}
			        set_normal_of_tri(ts);
			        set_normal_of_tri(ta);
			        break;
			    }
		        }
		    }
		}
		if (id < 3)
		    break;
	    }
	    if(j == nt)
	    {
	        printf("WARNING, synchronize_tris_at_subdomain_bdry, "
		       "suitable triangle is not found.\n");
	        print_tri(ts, ts->surf->interface);
	    }
	}
}		/*end synchronize_tris_at_subdomain_bdry*/

LOCAL	boolean	tri_cross_line(
	TRI		*tri,
	double		crx_coord,
	int		dir)
{
	double	min_coord, max_coord;
	double	crx_tol = tol1[dir];

	min_coord = max_coord = Coords(Point_of_tri(tri)[0])[dir];

	if (min_coord > Coords(Point_of_tri(tri)[1])[dir])
	    min_coord = Coords(Point_of_tri(tri)[1])[dir];
	if (min_coord > Coords(Point_of_tri(tri)[2])[dir])
	    min_coord = Coords(Point_of_tri(tri)[2])[dir];

	if (max_coord < Coords(Point_of_tri(tri)[1])[dir])
	    max_coord = Coords(Point_of_tri(tri)[1])[dir];
	if (max_coord < Coords(Point_of_tri(tri)[2])[dir])
	    max_coord = Coords(Point_of_tri(tri)[2])[dir];
	
	return (((min_coord - crx_coord) <= crx_tol) && 
	        ((crx_coord - max_coord) <= crx_tol)) ? YES : NO;
}		/*end tri_cross_line*/

/*NO    no bond side or no tri crx line
  YES   one tri crx line
*/
LOCAL	boolean	tri_bond_cross_line(
	TRI		*tri,
	double		crx_coord,
	int		dir)
{
	int		i;
	BOND		*b;
	BOND_TRI	**btris;

	for(i=0; i<3; i++)
	{
	    if(!is_side_bdry(tri, i))
	        continue;
	    b = Bond_on_side(tri, i);
	    for(btris = Btris(b); btris && *btris; btris++)
	    	if(tri_cross_line((*btris)->tri,crx_coord,dir))
	            return YES;
	}
	return NO;
}

LOCAL boolean tri_bond_cross_test(
	TRI		*tri,
	double		crx_coord,
	int		dir)
{
	int	i, j, n;
	TRI	**tris;
	POINT	*p;

	for(i = 0; i < 3; i++)
	{
	    if(Boundary_point(Point_of_tri(tri)[i]))
	    {
	        p = Point_of_tri(tri)[i];
		n = set_tri_list_around_point(p,tri,&tris,tri->surf->interface);
		if (tri_bond_cross_line(tris[0],crx_coord,dir) && 
	   	    tri_bond_cross_line(tris[n-1],crx_coord,dir))
	    	    return YES;
	    }
	}
	return NO;
}


LOCAL boolean match_tris_at_subdomain_bdry(
	SURFACE		*surf,
	SURFACE		*adj_surf,
	TRI		**tri_s,
	TRI		**tri_a,
	int		ns,
	int		na)
{
	TRI		*ta, *ts;
	int		i, j;
	int		ums,uma;
	static boolean	*ms = NULL, *ma = NULL;
	static int      ms_len = 0, ma_len = 0;

	DEBUG_ENTER(match_tris_at_subdomain_bdry)
	if (DEBUG)
	{
	    (void) printf("Entered match_tris_at_subdomain_bdry()\n");
	    (void) printf("tri_a: na=%d\n",na);
	    for (i = 0; i < na; ++i)
		print_tri(tri_a[i],adj_surf->interface);
	    (void) printf("tri_s: ns=%d\n",ns);
	    for (i = 0; i < ns; ++i)
		print_tri(tri_s[i],surf->interface);
	}

	if (ms_len < ns)
	{
	    ms_len = ns;
	    if (ms != NULL)
		free(ms);
	    uni_array(&ms,ms_len,sizeof(boolean));
	}
	if (ma_len < na)
	{
	    ma_len = na;
	    if (ma != NULL)
		free(ma);
	    uni_array(&ma,ma_len,sizeof(boolean));
	}

	for (i = 0; i < ns; ++i)
	    ms[i] = NO;
	for (i = 0; i < na; ++i)
	    ma[i] = NO;

	for (i = 0; i < na; ++i)
	{
	    ta = tri_a[i];
	    
	    for (j = 0; j < ns; ++j)
	    {
	        ts = tri_s[j];
	        if (match_two_tris(ts,ta))
		{
		    ma[i] = ms[j] = YES;
		    
		    if(debugging("app_tri"))
		        printf("ns  %d \n", j);

	    	    average_btris(ts, surf, ta, adj_surf);
		    merge_two_tris(ts,ta,surf,adj_surf);
		    break;
		}
	    }

	    remove_from_debug("app_tri");
	}
	ums = uma = 0;
	for (i = 0; i < ns; ++i)
	    if (ms[i] == NO)
		++ums;
	for (i = 0; i < na; ++i)
	    if (ma[i] == NO)
		++uma;
	if (ums != 0 || uma != 0)
	{
	    return NO;
	    (void) printf("WARNING in match_tris_at_subdomain_bdry(), "
	                  "unmatched local tris\n"
	                  "na = %d ns = %d\n"
	                  "ums = %d uma = %d\n",na,ns,ums,uma);
	    for (i = 0; i < ns; ++i)
	    	if (ms[i] == NO)
		    print_tri(tri_s[i],surf->interface);
	    (void) printf("unmatched adj tris:\n");
	    for (i = 0; i < na; ++i)
	    	if (ma[i] == NO)
		    print_tri(tri_a[i],adj_surf->interface);
	    
	    DEBUG_LEAVE(match_tris_at_subdomain_bdry)
	    clean_up(ERROR);
	    return NO;
	}

	DEBUG_LEAVE(match_tris_at_subdomain_bdry)
	return YES;
}		/*end match_tris_at_subdomain_bdry*/

/*
*			match_two_tris():
*
*	Determines whether two triangles are the same triangle up to a rotation
*	of the vertex indices.  Returns YES if the tris are the same,
*	otherwise returns NO.
*/

LOCAL boolean match_two_tris(
	TRI		*tri,
	TRI		*atri)
{
	int		j;
	POINT		**p, **ap;

	p = Point_of_tri(tri);
	ap = Point_of_tri(atri);
	for (j = 0; j < 3; ++j)
	    if (p[0] == ap[j])
		break;
	if (j == 3)
	    return NO;
	return ((p[1]==ap[Next_m3(j)]) && (p[2]==ap[Prev_m3(j)])) ? YES : NO;
}		/*end match_two_tris*/

/*
*			merge_two_tris3():
*
*	Merges two triangles ts and ta by replacing all linking
*	information to and from ta by linking with ts  while preserving
*	any nontrivial linking information in ts.
*/

LOCAL	void merge_two_tris3(
	TRI	*ts,
	TRI	*ta,
	SURFACE	*s,
	SURFACE *as)
{
	TRI		*tri;
	BOND_TRI	*btri;
	int		i, j;

	DEBUG_ENTER(merge_two_tris3)
	if (DEBUG)
	{
	    print_tri(ts,s->interface);
	    print_tri(ta,as->interface);
	}

	/* Set the links on the null side of ts by corresponding links to ta*/
	for (i = 0; i < 3; ++i)
	{
	    /* for grid based, this part is not necessary*/
            if(is_side_bdry(ta,i))
            {
	       if(Bond_on_side(ts,i) == NULL)
               {
                       if(Bond_on_side(ta,i) == NULL)
                           continue;
		       
	       	       printf("#merge_two_tris3: bond on block face.\n");
		       /*set tri_neighbor btri*/
                       Bond_tri_on_side(ts,i) = btri = Bond_tri_on_side(ta,i);
                       btri->tri = ts;
                       btri->surface = s;
               }
	       continue;
            }

	    if (Tri_on_side(ts,i) == NULL)
	    {
	    	if (Tri_on_side(ta,i) == NULL)
		{
		    continue;
		}
		/* set tri_neighbor tri*/
	    	Tri_on_side(ts,i) = tri = Tri_on_side(ta,i);
		
		for (j = 0; j < 3; ++j)
	    	{
	    	    if (Tri_on_side(tri,j) == ta)
	    	    	Tri_on_side(tri,j) = ts;
	    	}
	    }
	}
	
	/*remove ta from s tri list */
	remove_tri_from_surface(ta,s,NO);
	if (DEBUG)
	{
	    (void) printf("after merge_two_tris3\n");
	    print_tri(ts,s->interface);
	}
	DEBUG_LEAVE(merge_two_tris3)
}		/*end merge_two_tris3*/



/*
*			add_matching_pt_to_hash_table():
*
*	Creates a hashed list of matching points on intfc and adj_intfc,
*	where two points match if their positions are within a prescribed
*	floating point tolerance.  This is an extremely expensive function
*	and should be a candidate for an improved version.
*/

LOCAL boolean add_matching_pt_to_hash_table(
	TRI 	**tris_s,
	TRI	**tris_a,
	int	ns,
	int 	na,
	SURFACE	*ss,
	SURFACE	*sa,
	P_LINK	*p_table,
	int	p_size)
{
	boolean                  status;
	int 	                 i, tstn;
	POINT_LIST               *plists, *plista, *pls, *pla;
	static POINT_LIST_STORE  Pslist_store, Palist_store;

	plists = set_point_list(tris_s,ns,Hyper_surf(ss),&Pslist_store);
	plista = set_point_list(tris_a,na,Hyper_surf(sa),&Palist_store);

	tstn = 0;
	status = YES;
	for (pls = plists; pls != NULL; pls = pls->next)
	{
	    for (pla = plista; pla != NULL; pla = pla->next)
	    {
		/* Using global index matching */
		if (Gindex(pla->p) == Gindex(pls->p))
		{
		    tstn++;
	            (void) add_to_hash_table((POINTER)pla->p,(POINTER)pls->p,
				             p_table,p_size);
		    set_full_average(NO);
		    (void) average_points(NO,pla->p,pla->hse,pla->hs,
				          pls->p,pls->hse,pls->hs);
		    set_full_average(YES);
		    
		    if (pla->prev == NULL)
			plista = pla->next;
		    else
			pla->prev->next = pla->next;
		    if (pla->next != NULL)
			pla->next->prev = pla->prev;
		    break;
		}
	    }
	    /*can not find a corresponding point in plista, in this case, 
	      we continue search just set status = NO.  */
	    if (pla == NULL) 
	        status = NO;
	}
	if (plista != NULL)
	    status = NO;

	if(status == NO)
	{
	    double	*pa, *ps;
	    double	len, min_len;

	    printf("add_matching_pt_to_hash_table: check nearest point pairs\n");

	    for (pla = plista; pla != NULL; pla = pla->next)
	    {
	        min_len = HUGE_VAL;
		if (plists == NULL)
		{
		    (void) printf("plists is NULL\n");
		    continue;
		}
		for (pls = plists; pls != NULL; pls = pls->next)
	        {
		    len = distance_between_positions(Coords(pla->p),
					Coords(pls->p),3);
		    if(len < min_len)
		    {
		        min_len = len;
			pa = Coords(pla->p);
			ps = Coords(pls->p);
		    }
		}
		(void) printf("Gindex(pla->p) = %ld\n",Gindex(pla->p));
		(void) printf("Gindex(pls->p) = %ld\n",Gindex(pls->p));
		print_general_vector("pa", pa, 3, "\n");
		print_general_vector("ps", ps, 3, "\n");
	    }
	}

	return (ns == na) ? status : NO;
}		/*end add_matching_pt_to_hash_table*/

LOCAL	POINT_LIST	*set_point_list(
	TRI	         **tris,
	int              nt,
	HYPER_SURF       *hs,
	POINT_LIST_STORE *plist_store)
{
	POINT      **p;
	POINT_LIST *plist, *pl, Plist;
	TRI        *tri;
	int        i, j, max_np, tstnum;

	tstnum = 0;
	if (nt == 0)
	    return NULL;
	max_np = 3*nt + 1;
	if (max_np > plist_store->len)
	{
	    if (plist_store->pl != NULL)
		free(plist_store->pl);
	    plist_store->len = max_np;
	    uni_array(&plist_store->pl,sizeof(POINT_LIST),plist_store->len);
	}
	zero_scalar(plist_store->pl,sizeof(POINT_LIST)*plist_store->len);
	for (i = 0; i < nt; ++i)
	{
	    p = Point_of_tri(tris[i]);
	    for (j = 0; j < 3; ++j)
	        sorted(p[j]) = NO;
	}

	plist = plist_store->pl;
	pl = &Plist;
	for (i = 0; i < nt; ++i)
	{
	    tri = tris[i];
	    
	    p = Point_of_tri(tri);
	    for (j = 0; j < 3; ++j)
	    {
		if (!sorted(p[j]))
		{
		    pl->next = plist++;
		    pl->next->prev = pl;
		    pl = pl->next;
	            pl->p = p[j];
	            pl->hs = hs;
	            pl->hse = Hyper_surf_element(tri);
	            sorted(pl->p) = YES;
		    tstnum++; 
		}
	    }
	}
	Plist.next->prev = NULL;

	return Plist.next;
}		/*end set_point_list*/

LOCAL void clip_intfc_at_grid_bdry1(
	INTERFACE	*intfc)
{
	INTERFACE	*cur_intfc = current_interface();
	RECT_GRID	*gr = computational_grid(intfc);
	int		dim = intfc->dim;
	int		dir, nb;
	double		L[MAXD],U[MAXD];

	DEBUG_ENTER(clip_intfc_at_grid_bdry1)
	strip_subdomain_bdry_curves(intfc);
	set_current_interface(intfc);
	for (dir = 0; dir < dim; ++dir)
	{
	    L[dir] = gr->L[dir];
	    U[dir] = gr->U[dir];
	    if (gr->lbuf[dir] == 0) L[dir] -= 0.5*gr->h[dir];
	    if (gr->ubuf[dir] == 0) U[dir] += 0.5*gr->h[dir];
	    
	    if(rect_boundary_type(intfc,dir,0) == OPEN_BOUNDARY)
		L[dir] = gr->VL[dir];
	    if(rect_boundary_type(intfc,dir,1) == OPEN_BOUNDARY)
		U[dir] = gr->VU[dir];

            /* do not cut the reflect part */
            if(rect_boundary_type(intfc,dir,0) == REFLECTION_BOUNDARY)
                L[dir] = gr->VL[dir];
            if(rect_boundary_type(intfc,dir,1) == REFLECTION_BOUNDARY)
                U[dir] = gr->VU[dir];
	}
	
	for (dir = 0; dir < dim; ++dir)
	{
	    for (nb = 0; nb < 2; ++nb)
	    	open_null_sides1(intfc,L,U,dir,nb);
	}
	for (dir = 0; dir < dim; ++dir)
        for (nb = 0; nb < 2; ++nb)
            open_null_bonds(intfc,L,U,dir,nb);
	
	reset_intfc_num_points(intfc);
	set_current_interface(cur_intfc);
	DEBUG_LEAVE(clip_intfc_at_grid_bdry1)
}		/*end clip_intfc_at_grid_bdry1*/

EXPORT INTERFACE *collect_hyper_surface(
	Front *fr,
        int *owner,                     /* Destination of collection */
        int w_type)
{
	PP_GRID *pp_grid = fr->pp_grid;
	int 	*G = pp_grid->gmax;	
	int	dim = FT_Dimension();
	int 	dst_id,owner_id = domain_id(owner,G,dim);
	int	myid = pp_mynode();
	int	i,num_ids = pp_numnodes();
	boolean      sav_copy;
	INTERFACE *intfc = fr->interf;
	INTERFACE *cut_intfc,*recv_intfc;
	int myic[MAXD],ip[MAXD];
	int imax = 0;
	char fname[100];
	RECT_GRID *gr,*recv_gr;
	boolean status;

	for (i = 0; i < dim; ++i)
	    if (imax < G[i])  imax = G[i];
	find_Cartesian_coordinates(myid,pp_grid,myic);

	if (debugging("collect_intfc"))
	{
	    (void) printf("Entering collect_hyper_surface()\n");
	    (void) printf("myid = %d  owner_id = %d\n",myid,owner_id);
	    (void) printf("myic = (%d %d %d) owner = (%d %d %d)\n",myic[0],
				myic[1],myic[2],owner[0],owner[1],owner[2]);
	}

	    /* prepare interface to send */
	sav_copy = copy_intfc_states();
        set_copy_intfc_states(YES);

	cut_intfc = cut_intfc_to_wave_type(intfc,w_type);
	gr = computational_grid(cut_intfc);

	/* Patch in x-direction */
	if (myic[0] != owner[0]) 
	{
	    ip[0] = owner[0];
	    ip[1] = myic[1];
	    ip[2] = myic[2];
	    dst_id = domain_id(ip,G,dim);
	    send_interface(cut_intfc,dst_id);
	}
	else if (myic[0] == owner[0])
	{
	    ip[1] = myic[1];
	    ip[2] = myic[2];
	    for (i = 1; i < G[0]; ++i)
	    {
		ip[0] = myic[0] + i;
		if (ip[0] < G[0])
		{
	    	    dst_id = domain_id(ip,G,dim);
		    recv_intfc = receive_interface(dst_id);
		    status = buffer_extension3d3(cut_intfc,recv_intfc,
					0,1,status);
		    recv_gr = computational_grid(recv_intfc);
		    merge_rect_grids(gr,gr,recv_gr);
		    delete_interface(recv_intfc);
		}
		ip[0] = myic[0] - i;
		if (ip[0] >= 0)
		{
	    	    dst_id = domain_id(ip,G,dim);
		    recv_intfc = receive_interface(dst_id);
		    status = buffer_extension3d3(cut_intfc,recv_intfc,
					0,0,status);
		    recv_gr = computational_grid(recv_intfc);
		    merge_rect_grids(gr,gr,recv_gr);
		    delete_interface(recv_intfc);
		}
	    }
	}
	pp_gsync();

	/* Patch in y-direction */
	if (myic[0] == owner[0] && myic[1] != owner[1]) 
	{
	    ip[0] = owner[0];
	    ip[1] = owner[1];
	    ip[2] = myic[2];
	    dst_id = domain_id(ip,G,dim);
	    send_interface(cut_intfc,dst_id);
	}
	else if (myic[0] == owner[0] && myic[1] == owner[1])
	{
	    ip[0] = owner[0];
	    ip[2] = myic[2];
	    for (i = 1; i < G[1]; ++i)
	    {
		ip[1] = myic[1] + i;
		if (ip[1] < G[1])
		{
	    	    dst_id = domain_id(ip,G,dim);
		    recv_intfc = receive_interface(dst_id);
		    status = buffer_extension3d3(cut_intfc,recv_intfc,
					1,1,status);
		    recv_gr = computational_grid(recv_intfc);
		    merge_rect_grids(gr,gr,recv_gr);
		    delete_interface(recv_intfc);
		}
		ip[1] = myic[1] - i;
		if (ip[1] >= 0)
		{
	    	    dst_id = domain_id(ip,G,dim);
		    recv_intfc = receive_interface(dst_id);
		    status = buffer_extension3d3(cut_intfc,recv_intfc,
					1,0,status);
		    recv_gr = computational_grid(recv_intfc);
		    merge_rect_grids(gr,gr,recv_gr);
		    delete_interface(recv_intfc);
		}
	    }
	}
	pp_gsync();

	/* Patch in z-direction */
	if (myic[0] == owner[0] && myic[1] == owner[1] && myic[2] != owner[2]) 
	{
	    ip[0] = owner[0];
	    ip[1] = owner[1];
	    ip[2] = owner[2];
	    dst_id = domain_id(ip,G,dim);
	    send_interface(cut_intfc,dst_id);
	}
	else if (myic[0] == owner[0] && myic[1] == owner[1] && 
		 myic[2] == owner[2])
	{
	    ip[0] = owner[0];
	    ip[1] = owner[1];
	    for (i = 1; i < G[2]; ++i)
	    {
		ip[2] = myic[2] + i;
		if (ip[2] < G[2])
		{
	    	    dst_id = domain_id(ip,G,dim);
		    recv_intfc = receive_interface(dst_id);
		    status = buffer_extension3d3(cut_intfc,recv_intfc,
					2,1,status);
		    recv_gr = computational_grid(recv_intfc);
		    merge_rect_grids(gr,gr,recv_gr);
		    delete_interface(recv_intfc);
		}
		ip[2] = myic[2] - i;
		if (ip[2] >= 0)
		{
	    	    dst_id = domain_id(ip,G,dim);
		    recv_intfc = receive_interface(dst_id);
		    status = buffer_extension3d3(cut_intfc,recv_intfc,
					2,0,status);
		    recv_gr = computational_grid(recv_intfc);
		    merge_rect_grids(gr,gr,recv_gr);
		    delete_interface(recv_intfc);
		}
	    }
	}
	pp_gsync();
	if (debugging("collect_intfc"))
	{
	    (void) printf("Leaving collect_hyper_surface()\n");
	    sprintf(fname,"final-intfc.%d",myid);
	    gview_plot_interface(fname,cut_intfc);
	    (void) printf("Checking consistency:\n");
	    null_sides_are_consistent();
	    consistent_interface(cut_intfc);
	    (void) printf("Passed consistent_interface()\n");
	}
	if (myid == owner_id)
	{
	    install_subdomain_bdry_curves(cut_intfc);    
	    return cut_intfc;
	}
	else
	{
	    delete_interface(cut_intfc);
	    return NULL;
	}
}	/* end collect_hyper_surface */

#define		MAX_DELETE	20

LOCAL INTERFACE *cut_intfc_to_wave_type(
	INTERFACE *intfc,
	int w_type)
{
	INTERFACE *tmp_intfc = copy_interface(intfc);
	SURFACE **s,*surfs_del[MAX_DELETE];
	CURVE **c;
	int i,dir,nb,num_delete = 0;
	INTERFACE *cut_intfc;
	RECT_GRID *gr = computational_grid(intfc);
	int dim = gr->dim;
	char fname[100];

	set_floating_point_tolerance3(computational_grid(intfc)->h);
	intfc_surface_loop(tmp_intfc,s)
	{
	    if (wave_type(*s) != w_type && w_type != ANY_WAVE_TYPE)
		surfs_del[num_delete++] = *s;
	    else
	    {
		for (dir = 0; dir < dim; ++dir)
		for (nb = 0; nb < 2; ++nb)
		{
		    open_surf_null_sides(*s,gr->L,gr->U,dir,nb);
		}
	    }
	}
	for (i = 0; i < num_delete; ++i)
	    delete_surface_set(surfs_del[i]);

	for (dir = 0; dir < dim; ++dir)
	for (nb = 0; nb < 2; ++nb)
	    open_null_bonds(tmp_intfc,gr->L,gr->U,dir,nb);

	reset_intfc_num_points(tmp_intfc);
	cut_intfc = copy_interface(tmp_intfc);
	delete_interface(tmp_intfc);
	return cut_intfc;
}	/* end cut_intfc_to_wave_type */

LOCAL void delete_surface_set(
	SURFACE *surf)
{
	CURVE **c,**curves_del = NULL;
	NODE **n,**nodes_del = NULL;

	surf_pos_curve_loop(surf,c)
	{
	    unique_add_to_pointers((POINTER)*c,(POINTER**)&curves_del);
	    unique_add_to_pointers((POINTER)(*c)->start,(POINTER**)&nodes_del);
	    unique_add_to_pointers((POINTER)(*c)->end,(POINTER**)&nodes_del);
	}
	surf_neg_curve_loop(surf,c)
	{
	    unique_add_to_pointers((POINTER)*c,(POINTER**)&curves_del);
	    unique_add_to_pointers((POINTER)(*c)->start,(POINTER**)&nodes_del);
	    unique_add_to_pointers((POINTER)(*c)->end,(POINTER**)&nodes_del);
	}
	delete_surface(surf);
	for (c = curves_del; c && *c; ++c)
	{
	    if ((*c)->pos_surfaces == NULL && (*c)->neg_surfaces == NULL)
	    {
		delete_curve(*c);
		delete_from_pointers((POINTER)(*c),(POINTER**)&curves_del);
		c--;
	    }
	}
	for (n = nodes_del; n && *n; ++n)
	{
	    if ((*n)->in_curves == NULL && (*n)->out_curves == NULL)
	    {
		delete_node(*n);
		delete_from_pointers((POINTER)(*n),(POINTER**)&nodes_del);
		n--;
	    }
	}
}	/* end delete_surface_set */

#define		MAX_NUM_OTHER_CURVES	500

LOCAL boolean append_other_curves3(
	INTERFACE *intfc,
	INTERFACE *adj_intfc,
	P_LINK *p_table,
	int p_size)
{
	CURVE **cc;
	CURVE *c[MAX_NUM_OTHER_CURVES];
	CURVE *ac[MAX_NUM_OTHER_CURVES];
	int i,j,num_c,num_ac;
	BOND *b,*ba;
	boolean bond_matched;
	POINT *p;

	num_c = num_ac = 0;
	for (cc = intfc->curves; cc && *cc; ++cc)
	{
	    if (I_NumOfCurveSurfaces(*cc) != 0)
	    	continue;
	    if (num_c >= MAX_NUM_OTHER_CURVES)
	    {
		printf("In append_other_curves3(): num_c = %d\n",num_c);
		printf("MAX_NUM_OTHER_CURVES too small!\n");
		clean_up(ERROR);
	    }
	    c[num_c++] = *cc;
	}
	for (cc = adj_intfc->curves; cc && *cc; ++cc)
	{
	    if (I_NumOfCurveSurfaces(*cc) != 0)
	    	continue;
	    if (num_ac >= MAX_NUM_OTHER_CURVES)
	    {
		printf("In append_other_curves3(): num_ac = %d\n",num_ac);
		printf("MAX_NUM_OTHER_CURVES too small!\n");
		clean_up(ERROR);
	    }
	    ac[num_ac++] = *cc;
	}
	for (i = 0; i < num_ac; ++i)
	{
	    for (j = 0; j < num_c; ++j)
	    {
		if (c[j] == NULL)	/* already matched */
		    continue;
		if (Gindex(ac[i]) != Gindex(c[j]))
		    continue;
	    	for (ba = ac[i]->first; ba != NULL; ba = ba->next)
	    	{
		    for (b = c[j]->first; b != NULL; b = b->next)
		    {
			bond_matched = NO;
		    	if (bond_match3(b,ba))
		    	{
			    bond_matched = YES;
			    p = (POINT*)find_from_hash_table((POINTER)ba->start,
                                                p_table,p_size);
			    if (p == NULL)
			    {
				(void) add_to_hash_table((POINTER)ba->start,
						(POINTER)b->start,
                                                p_table,p_size);
			    }
			    else if (p != b->start)
			    {
				screen("Bond start not from hashing table!\n");
				clean_up(ERROR);
			    }
			    p = (POINT*)find_from_hash_table((POINTER)ba->end,
                                                p_table,p_size);
			    if (p == NULL)
			    {
				(void) add_to_hash_table((POINTER)ba->end,
						(POINTER)b->end,
                                                p_table,p_size);
			    }
			    else if (p != b->end)
			    {
				screen("Bond end not from hashing table!\n");
				clean_up(ERROR);
			    }
		    	}
			if (bond_matched) break;
		    }
	    	}	    
	    }
	}
}	/* end append_other_curves1 */

LOCAL	boolean bond_match3(
	BOND *b,
	BOND *ba)
{
	int i;
	if (Gindex(b->start) != Gindex(ba->start)) return NO;
	if (Gindex(b->end) != Gindex(ba->end)) return NO;
	return YES;
}	/* end bond_match3 */

LOCAL  void  merge_overlap_nodes(	
	INTERFACE *intfc)
{
	int i,j,num_nodes = I_NumOfIntfcNodes(intfc);
	boolean *node_merged;
	NODE **n,**node_list,*n1,*n2;
	CURVE **c;

	uni_array(&node_merged,num_nodes,sizeof(boolean));
	uni_array(&node_list,num_nodes,sizeof(NODE*));
	for (i = 0; i < num_nodes; ++i)
	    node_merged[i] = NO;
	i = 0;
	intfc_node_loop(intfc,n)
	{
	    node_list[i++] = *n;
	}
	for (i = 0; i < num_nodes-1; ++i)
	{
	    if (node_merged[i]) continue;
	    n1 = node_list[i];
	    node_merged[i] = YES;
	    for (j = i+1; j < num_nodes; ++j)
	    {
	    	n2 = node_list[j];
		if (n1->posn != n2->posn) continue;
	    	if (node_merged[j]) continue;
		/*
		printf("merging i = %d  j = %d\n",i,j);
		printf("p1 = %d  p2 = %d\n",n1->posn,n2->posn);
		printf("p1 = %f %f %f\n",Coords(n1->posn)[0],
				Coords(n1->posn)[1],Coords(n1->posn)[2]);
		printf("Node n1:\n");
		print_node(n1);
		printf("Node n2:\n");
		print_node(n2);
		*/
	    	node_merged[j] = YES;
		node_in_curve_loop(n2,c)
		{
		    change_node_of_curve(*c,NEGATIVE_ORIENTATION,n1);
		}
		node_out_curve_loop(n2,c)
		{
		    change_node_of_curve(*c,POSITIVE_ORIENTATION,n1);
		}
	    }
	}
}	/* end merge_overlap_nodes */

LOCAL	void print_unmatched_tris(
	TRI **tris_s,
	TRI **tris_a,
	int ns,
	int na)
{
	int i,j,k;
	long gindex_s[3],gindex_a[3];
	boolean match_found;

	(void) printf("Number of local tris = %d\n",ns);
	(void) printf("Number of adj   tris = %d\n",na);
	if (ns > na)
	{
	    for (i = 0; i < ns; ++i)
	    {
		match_found = NO;
		for (k = 0; k < 3; ++k)	
		    gindex_s[k] = Gindex(Point_of_tri(tris_s[i])[k]);
		for (j = 0; j < na; ++j)
		{
		    for (k = 0; k < 3; ++k)	
		    	gindex_a[k] = Gindex(Point_of_tri(tris_a[j])[k]);
		    for (k = 0; k < 3; ++k)	
		    {
			if (gindex_s[0] == gindex_a[k%3] &&
			    gindex_s[1] == gindex_a[(k+1)%3] &&
			    gindex_s[2] == gindex_a[(k+2)%3])
			{
			    match_found = YES;
			    break;
			}
		    }
		    if (match_found) break;
		}
		if (match_found == NO)
		{
		    (void) printf("Unmatched local tri:\n");
		    (void) printf("Global indices: %ld %ld %ld\n",
				Gindex(Point_of_tri(tris_s[i])[0]),
				Gindex(Point_of_tri(tris_s[i])[1]),
				Gindex(Point_of_tri(tris_s[i])[2]));
		    (void) printf("Coordinates:\n");
		    print_tri_coords(tris_s[i]);
		}
	    }
	}
	else
	{
	    for (i = 0; i < na; ++i)
	    {
		match_found = NO;
		for (k = 0; k < 3; ++k)	
		    gindex_a[k] = Gindex(Point_of_tri(tris_a[i])[k]);
		for (j = 0; j < ns; ++j)
		{
		    for (k = 0; k < 3; ++k)	
		    	gindex_s[k] = Gindex(Point_of_tri(tris_s[j])[k]);
		    for (k = 0; k < 3; ++k)	
		    {
			if (gindex_s[0] == gindex_a[k%3] &&
			    gindex_s[1] == gindex_a[(k+1)%3] &&
			    gindex_s[2] == gindex_a[(k+2)%3])
			{
			    match_found = YES;
			    break;
			}
		    }
		    if (match_found) break;
		}
		if (match_found == NO)
		{
		    (void) printf("Unmatched adj tri:\n");
		    (void) printf("Global indices: %ld %ld %ld\n",
				Gindex(Point_of_tri(tris_a[i])[0]),
				Gindex(Point_of_tri(tris_a[i])[1]),
				Gindex(Point_of_tri(tris_a[i])[2]));
		    (void) printf("Coordinates:\n");
		    print_tri_coords(tris_a[i]);
		}
	    }
	}
}	/* end print_unmatched_tris */

LOCAL void exchange_intfc_extra(INTERFACE *intfc)
{
	int i,j,num_nodes;
	int global_index;
	int size_of_extra;
	NODE **n;
	POINTER extra;
	boolean extra_assigned;

	num_nodes = 0;
	intfc_node_loop(intfc,n)
	    if ((*n)->extra != NULL) num_nodes++;

	for (i = 0; i < pp_numnodes(); ++i)
	{
	    if (i == pp_mynode()) continue;
	    pp_send(30,&num_nodes,sizeof(int),i);
	    intfc_node_loop(intfc,n)
	    {
	    	if ((*n)->extra == NULL) continue;
	    	pp_send(31,&(Gindex((*n)->posn)),sizeof(long),i);
	    	pp_send(32,&((*n)->size_of_extra),sizeof(int),i);
	    	pp_send(33,(*n)->extra,(*n)->size_of_extra,i);
	    }
	}
	pp_gsync();
	for (i = 0; i < pp_numnodes(); ++i)
	{
	    if (i == pp_mynode()) continue;
	    pp_recv(30,i,&num_nodes,sizeof(int));
	    for (j = 0; j < num_nodes; ++j)
	    {
	    	pp_recv(31,i,&global_index,sizeof(long));
	    	pp_recv(32,i,&size_of_extra,sizeof(int));
		FT_ScalarMemoryAlloc((POINTER*)&extra,size_of_extra);
	    	pp_recv(33,i,extra,size_of_extra);
		extra_assigned = NO;
		intfc_node_loop(intfc,n)
		{
		    if ((*n)->extra == NULL && 
			Gindex((*n)->posn) == global_index)
		    {
			(*n)->extra = extra;
			(*n)->size_of_extra = size_of_extra;
			extra_assigned = YES;
		    }
		}
		if (extra_assigned == NO)
		    free_these(1,extra);
	    }
	}
}	/* end exchange_intfc_extra */

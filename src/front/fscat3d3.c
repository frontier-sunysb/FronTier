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
LOCAL	POINT_LIST	*set_point_list(TRI**,int,HYPER_SURF*,
					POINT_LIST_STORE*);
LOCAL	boolean	add_matching_pt_to_hash_table(TRI**,TRI**,int,int,SURFACE*,
					      SURFACE*,P_LINK*,int);
LOCAL	boolean	buffer_extension3d3(INTERFACE*,INTERFACE*,int,int,boolean);
LOCAL	boolean	match_tris_at_subdomain_bdry(SURFACE*,SURFACE*,TRI**,TRI**,
	                                     int,int);
LOCAL	boolean	match_two_tris(TRI*,TRI*);
LOCAL	boolean	tri_cross_line(TRI*,double,int);
LOCAL	boolean	tri_out_domain1(TRI*,double*,double*,int,int);
LOCAL	int	append_adj_intfc_to_buffer3(INTERFACE*,INTERFACE*,
					   RECT_GRID*,int,int);
LOCAL	int	append_buffer_surface3(SURFACE*,SURFACE*,RECT_GRID*,int,int,
				      P_LINK*,int);
LOCAL	void	clip_intfc_at_grid_bdry1(INTERFACE*);
LOCAL	boolean	tri_bond_out_domain(TRI*,double*,double*,int,int);
LOCAL	boolean	tri_bond_cross_line(TRI*,double,int);
LOCAL	boolean	tri_bond_cross_test(TRI*,double,int);
LOCAL	void	copy_tri_state_to_btri(BOND_TRI*,BOND*,ORIENTATION,INTERFACE *);
LOCAL	void	merge_point_pointers_at_subdomain_bdry(TRI**,TRI**,
						       int,P_LINK*,int);
LOCAL	void	strip_curve_from_surf(CURVE*,SURFACE*,ORIENTATION);
LOCAL	void	synchronize_tris_at_subdomain_bdry(TRI**,TRI**,int,P_LINK*,int);
	void 	sep_common_pt_for_open_bdry(INTERFACE*);
LOCAL	INTERFACE	*cut_buf_interface1(INTERFACE*,int,int,int*,int*);
LOCAL	COMPONENT buffer_component(INTERFACE*,int,int);
LOCAL   double   *constr_position(double*,double*,
                        boolean (*constr_func)(POINTER,double*),POINTER);

LOCAL	double	tol1[MAXD]; /*TOLERANCE*/
#define MAX_NULL_SIDE      8000

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

	set_floating_point_tolerance1(fr->rect_grid->h);
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

	if(fr->step == -1)
	    add_to_debug("out_surf");
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
			adj_intfc[(j+1)%2] = receive_interface(dst_id);
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
	    reset_intfc_num_points(intfc);
	}

	if (status == FUNCTION_SUCCEEDED)
	{
	    sep_common_pt_for_open_bdry(intfc);

	    install_subdomain_bdry_curves(intfc);

	    reset_intfc_num_points(intfc);
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
	int		p_size;		/*Size of space allocated for p_table*/
	static P_LINK	*p_table = NULL;/* Table of matching points on intfc
					 * and adj_intfc*/
	static int      len_p_table = 0;
	boolean      	corr_surf_found;

	DEBUG_ENTER(append_adj_intfc_to_buffer3)

	if (DEBUG)
	{
	    char dname[256];
	    static int ntimes[3][2];
	    static const char *strdir[3] = { "X", "Y", "Z" };
	    static const char *strnb[2] = { "LOWER", "UPPER" };

	    (void) sprintf(dname,"fscatter/"
				 "append_adj_intfc_to_buffer3/Data%d-%s_%s/%s",
			         ntimes[dir][nb],strdir[dir],strnb[nb],
				 "intfc_gv");
	    gview_plot_interface(dname,intfc);
	    (void) sprintf(dname,"fscatter/"
				 "append_adj_intfc_to_buffer3/Data%d-%s_%s/%s",
			         ntimes[dir][nb],strdir[dir],strnb[nb],
				 "adj_intfc_gv");
	    gview_plot_interface(dname,adj_intfc);

	    ++ntimes[dir][nb];
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
	for (as = adj_intfc->surfaces; as && *as; ++as)
	{
	    corr_surf_found = NO;
	    if (wave_type(*as) == FIRST_SCALAR_PHYSICS_WAVE_TYPE && 
		debugging("step_out"))
	    {    
		add_to_debug("out_matching");
	        printf("#out_mathcing\n");
	    }
	    for (s = intfc->surfaces; s && *s; ++s)
	    {
		/*
		*  COMMENT -
		*  The Hyper_surf_index() function is not
		*  fully supported.  This will fail in the
		*  presence of interface changes in topology
		*  TODO: FULLY SUPPORT THIS OBJECT
		*/
                
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
	
	merge_curves(intfc,adj_intfc);
	
	set_current_interface(cur_intfc);
	
	DEBUG_LEAVE(append_adj_intfc_to_buffer3)
	return YES;
}		/*end append_adj_intfc_to_buffer3*/

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
	
	if(debugging("out_surf"))
	{
	    char   s[40];
	    
	    printf("#check curve_surf\n");
	    check_surface_curve(surf);
	    printf("#check_curve adj_surf\n");
	    check_surface_curve(adj_surf);

	    sprintf(s, "surf_%d", pp_mynode());
	    tecplot_surface(s, NULL, surf);
	    sprintf(s, "adj_surf_%d", pp_mynode());
	    tecplot_surface(s, NULL, adj_surf);
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
		tris_a[na++] = tri;
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
	    screen("WARNING in append_buffer_surface3(), can't match "
	           "tris on with buffer surface\n");
	    FILE  *fp;
	    char  s[50];
	
            sprintf(s, "%sapp_nsna_%d.plt", get_directory(), pp_mynode());
	    printf("#error surface file %s\n", s);

	    fp = fopen(s, "w");
	    tecplot_show_tris("ns", tris_s, ns, fp);
	    tecplot_show_tris("na", tris_a, na, fp);
	    
	    tecplot_surface(NULL, fp, surf);
	    tecplot_surface(NULL, fp, adj_surf);
	    fclose(fp);
	    
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
		/*
	        if(the_tri(tri))
		{
		    printf("#copied tri\n");
		    print_tri(tri, adj_surf->interface);
		}
		*/
	    }
	}

	/*adjoin adj_surf tri list to surf tri list */
	last_tri(surf)->next = first_tri(adj_surf);
	first_tri(adj_surf)->prev = last_tri(surf);
	link_tri_list_to_surface(first_tri(surf),last_tri(adj_surf),surf);
	
	/* BEGIN curves in adj_surf should be added to surf. */
	for(c=adj_surf->pos_curves; c && *c; c++)
	    if(! delete_from_pointers(adj_surf, &(*c)->pos_surfaces))
	    {
	        printf("ERROR: in append_buffer_surface3, "
		       "adj_surf and pos_curves are not paired.\n");
		clean_up(ERROR);
	    }
	    else
	        install_curve_in_surface_bdry(surf, *c, POSITIVE_ORIENTATION);

	for(c=adj_surf->neg_curves; c && *c; c++)
	    if(! delete_from_pointers(adj_surf, &(*c)->neg_surfaces))
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
				/*
			        if(the_tri(ta))
				{
				    printf("#sync tris  \n");
				    print_tri(ts, ts->surf->interface);
				    print_tri(ta, ta->surf->interface);
				    add_to_debug("ts_tst");
				}
				*/

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

//NO    no bond side or no tri crx line
//YES   one tri crx line
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

	if(tri_bond_cross_line(tri,crx_coord,dir))
	    return YES;
	j = 0;
	for(i=0; i<3; i++)
	{
	    if(Boundary_point(Point_of_tri(tri)[i]))
	    {
	        p = Point_of_tri(tri)[i];
	        j++;
	    }
	}

	if(j != 1)
	    return NO;
	n = set_tri_list_around_point(p, tri, &tris, tri->surf->interface);
	
	if(tri_bond_cross_line(tris[0],crx_coord,dir) && 
	   tri_bond_cross_line(tris[n-1],crx_coord,dir) )
	    return YES;
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
	    /*
	    if(the_tri(ta))
	    {
	        printf("#append tri\n");
		print_tri(ta, surf->interface);
		add_to_debug("app_tri");
	    }
	    */
	    
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
		for (i = 0; i < 3; ++i) /*Floating point TOLERANCE test*/
	            if (fabs(Coords(pla->p)[i]-Coords(pls->p)[i]) > tol1[i])
			break;
		if (i == 3)
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
		for (pls = plists; pls != NULL; pls = pls->next)
	        {
		    len = distance_between_positions(Coords(pla->p), Coords(pls->p), 3);
		    if(len < min_len)
		    {
		        min_len = len;
			pa = Coords(pla->p);
			ps = Coords(pls->p);
		    }
		}
		
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

LOCAL boolean tri_out_domain1(
	TRI		*tri,
	double		*L,
	double		*U,
	int		dir,
	int		nb)
{
	POINT **p;
	int   i;

	p = Point_of_tri(tri);
	if (nb == 0)
	{
	    for (i = 0; i < 3; ++i)
	    {
		if ((L[dir] - Coords(p[i])[dir]) <= tol1[dir])
	    	    return NO;
	    }
	}
	else
	{
	    for (i = 0; i < 3; ++i)
	    {
		if ((Coords(p[i])[dir] - U[dir]) <= tol1[dir])
	    	    return NO;
	    }
	}
	return YES;
}	/* end tri_out_domain1 */

LOCAL boolean tri_bond_out_domain(
	TRI		*tri,
	double		*L,
	double		*U,
	int		dir,
	int		nb)
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
	        if(!tri_out_domain1((*btris)->tri, L, U, dir, nb))
		    return NO;
	}
	return YES;
}

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

            //do not cut the reflect part
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
	cut_out_curves_in_buffer(intfc);
	reset_intfc_num_points(intfc);
	set_current_interface(cur_intfc);
	DEBUG_LEAVE(clip_intfc_at_grid_bdry1)
}		/*end clip_intfc_at_grid_bdry1*/

LOCAL void merge_point_pointers_at_subdomain_bdry(
	TRI             **tris_a,
	TRI		**tris_s,
	int		nt,
	P_LINK		*p_table,
	int		p_size)
{
	POINT		*p, *ap;
	int		i, j;

        for (i = 0; i < nt; ++i)
	{
	    sorted(Point_of_tri(tris_s[i])[0]) = NO; 
	    sorted(Point_of_tri(tris_s[i])[1]) = NO;
	    sorted(Point_of_tri(tris_s[i])[2]) = NO;
	    sorted(Point_of_tri(tris_a[i])[0]) = NO; 
	    sorted(Point_of_tri(tris_a[i])[1]) = NO;
	    sorted(Point_of_tri(tris_a[i])[2]) = NO;
	}
	for (i = 0; i < nt; ++i)
	{
	    for (j = 0; j < 3; ++j)
	    {
		ap = Point_of_tri(tris_a[i])[j];
		if (sorted(ap) == NO)
		{
	            p = (POINT*)find_from_hash_table((POINTER)ap,
						     p_table,p_size);
		    if(p != Point_of_tri(tris_a[i])[j])
		        printf("#merge_point  \n");
		    
		    Point_of_tri(tris_a[i])[j] = p;
		    sorted(p) = YES;
		}
	    }
	}
}	/*end merge_point_pointers_at_subdomain_bdry*/


LOCAL	void	copy_tri_state_to_btri(
	BOND_TRI	*btri,
	BOND		*b,
	ORIENTATION	orient,
	INTERFACE	*intfc)
{
	POINT		*p;
	Locstate	sl, sr;
	NODE		*n;
	CURVE		**c;
	BOND_TRI	**bt;
	boolean		found;
	size_t	 	sizest = size_of_state(intfc);

	if (orient == POSITIVE_ORIENTATION)
	{
	    p = b->start;
	    sl = left_start_btri_state(btri);
	    sr = right_start_btri_state(btri);
	}
	else
	{
	    p = b->end;
	    sl = left_end_btri_state(btri);
	    sr = right_end_btri_state(btri);
	}

	//in install_subdomain_bdry_curve the new inserted point is always 
	//a node point.
	if ((n = node_of_point(p,intfc)) == NULL)
	{
	    ft_assign(sl,left_state(p),sizest);
	    ft_assign(sr,right_state(p),sizest);
	    return;
	}

	found = NO;
	// deal with 3 comp curve
	for(c = n->out_curves; c && *c; c++)
	{
	    if(*c == btri->curve)
	        continue;
	    found = YES;
	    for(bt=Btris((*c)->first); bt && *bt; bt++)
	        if((*bt)->surface == btri->surface)
		{
		    ft_assign(sl, left_start_btri_state(*bt), sizest);
		    ft_assign(sr, right_start_btri_state(*bt), sizest);
		    return;
		}
	}
	
	for(c = n->in_curves; c && *c; c++)
	{
	    if(*c == btri->curve)
	        continue;
	    found = YES;
	    for(bt=Btris((*c)->last); bt && *bt; bt++)
	        if((*bt)->surface == btri->surface)
		{
		    ft_assign(sl, left_end_btri_state(*bt), sizest);
		    ft_assign(sr, right_end_btri_state(*bt), sizest);
		    return;
		}
	}

	//Now curve is the only curve in node n
	if(!found)
	{
	    ft_assign(sl,left_state(p),sizest);
	    ft_assign(sr,right_state(p),sizest);
	    return;
	}

	printf("ERROR copy_tri_state_to_btri, no suitable btri found.\n");
	clean_up(ERROR);

}		/*end copy_tri_state_to_btri*/

LOCAL void strip_curve_from_surf(
	CURVE		*curve,
	SURFACE		*surf,
	ORIENTATION	orient)
{
	BOND		*b;
	BOND_TRI	**btris;
	POINT		*p;
	TRI		*tri;
	int		i;
	size_t	 	sizest = size_of_state(surf->interface);
	
	DEBUG_ENTER(strip_curve_from_surf)
	for (b = curve->first; b != NULL; b = b->next)
	{
	    Boundary_point(b->start) = Boundary_point(b->end) = 0;
	    for (btris = Btris(b); btris && *btris; ++btris)
	    {
		tri = (*btris)->tri;
		for (i = 0; i < 3; ++i)
		{
		    if (Bond_tri_on_side(tri,i) == *btris)
		    {
			if(sizest > 0)
			{
			    p = (*btris)->bond->start;
			    ft_assign(left_state(p),
				left_start_btri_state(*btris),sizest);
			    ft_assign(right_state(p),
				right_start_btri_state(*btris),sizest);
			    
			    p = (*btris)->bond->end;
			    ft_assign(left_state(p),
				left_end_btri_state(*btris),sizest);
			    ft_assign(right_state(p),
				right_end_btri_state(*btris),sizest);
			}

		    	Bond_tri_on_side(tri,i) = NULL;
		    	set_side_bdry(Boundary_tri(tri),i,NO);
		    	break;
		    }
		}
		if (!delete_from_pointers(*btris,&Btris(b)))
		{
		    screen("ERROR in strip_curve_from_surf(), "
		           "delete_from_pointers() failed\n");
		    clean_up(ERROR);
		}
	    }
	}
	if (orient == POSITIVE_ORIENTATION)
	{
	    if (!(delete_from_pointers(curve,&surf->pos_curves) &&
		     delete_from_pointers(surf,&curve->pos_surfaces)))
	    {
		screen("ERROR in strip_curve_from_surf(), "
		       "delete_from_pointers() failed\n");
		clean_up(ERROR);
	    }
	}
	else
	{
	    if (!(delete_from_pointers(curve,&surf->neg_curves) &&
	    	     delete_from_pointers(surf,&curve->neg_surfaces)))
	    {
	    	screen("ERROR in strip_curve_from_surf(), "
		       "delete_from_pointers() failed\n");
		clean_up(ERROR);
	    }
	}
	DEBUG_LEAVE(strip_curve_from_surf)
}		/*end strip_curve_from_surf*/


LOCAL	INTERFACE  *cut_buf_interface1(
	INTERFACE	*intfc,
	int		dir,
	int		nb,
	int		*me,
	int		*him)
{
	INTERFACE	  *sav_intfc, *tmp_intfc, *buf_intfc;
	RECT_GRID	  *gr = computational_grid(intfc);
	RECT_GRID	  dual_gr;
	double		  L[3], U[3];
	boolean		  sav_copy = copy_intfc_states();
	char              dname[1024]; 
	static const char *strdir[3] = { "X", "Y", "Z" };
	static const char *strnb[2] = { "LOWER", "UPPER" };

	DEBUG_ENTER(cut_buf_interface1)

	if (DEBUG)
	{
	    static int into[3][2];
	    (void) sprintf(dname,"fscatter/cut_buf_interface1/Into%d-%s_%s",
			   into[dir][nb],strdir[dir],strnb[nb]);
	    (void) printf("cut_buf_interface1() at ENTRY\n"
			  "  dir = %s  nb = %s  me = ( %d, %d, %d )  "
			  "him = ( %d, %d, %d )\n",
			  strdir[dir],strnb[nb],me[0],me[1],me[2],
			  him[0],him[1],him[2]);
	    summarize_interface(dname,"intfc",intfc,XY_PLANE,
				"cut_buf_interface1","intfc at ENTER");
	    ++into[dir][nb];
	}

	set_copy_intfc_states(YES);
	sav_intfc = current_interface();

	set_size_of_intfc_state(size_of_state(intfc));
	tmp_intfc = copy_interface(intfc);
	set_dual_grid(&dual_gr,gr);

	if (nb == 0)
	{
	    L[dir] = gr->L[dir];
	    U[dir] = gr->L[dir]+(gr->L[dir]-dual_gr.VL[dir])+0.5*gr->h[dir];
	}
	else
	{
	    L[dir] = gr->U[dir]-(dual_gr.VU[dir]-gr->U[dir])-0.5*gr->h[dir];
	    U[dir] = gr->U[dir];
	}

	open_null_sides1(tmp_intfc,L,U,dir,(nb+1)%2);

	/*
	 * Split curves that have been disconnected from surfaces
	 * and delete those sections that lie entirely within the
	 * subdomain boundary.
	 */

	cut_out_curves_in_buffer(tmp_intfc);

	reset_intfc_num_points(tmp_intfc);

	if (me[dir] == him[dir])
	{
	    buf_intfc = tmp_intfc;
	}
	else
	{
	    set_size_of_intfc_state(size_of_state(intfc));
	    buf_intfc = copy_interface(tmp_intfc);
	    (void) delete_interface(tmp_intfc);
	}
	set_copy_intfc_states(sav_copy);
	set_current_interface(sav_intfc);

	if (DEBUG)
	{
	    static int outof[3][2];
	    (void) sprintf(dname,"fscatter/cut_buf_interface1/Outof%d-%s_%s",
			   outof[dir][nb],strdir[dir],strnb[nb]);
	    summarize_interface(dname,"intfc",intfc,XY_PLANE,
				"cut_buf_interface1","intfc at EXIT");
	    summarize_interface(dname,"buf_intfc",buf_intfc,XY_PLANE,
				"cut_buf_interface1","buf_intfc at EXIT");
	    ++outof[dir][nb];
	}
	DEBUG_LEAVE(cut_buf_interface1)
	return buf_intfc;
}		/*end cut_buf_interface1*/

LOCAL COMPONENT buffer_component(
	INTERFACE *intfc,
	int dir,
	int nb)
{
	RECT_GRID *gr = computational_grid(intfc);
	double *L = gr->L;
	double *U = gr->U;
	double *h = gr->h;
	double coords[MAXD];
	int i, comp;

	if (intfc->surfaces == NULL) 
	    return intfc->default_comp;
	for (i = 0; i < 3; ++i)
	    coords[i] = 0.5*(L[i] + U[i]);
	coords[dir] = (nb == 0) ? (L[dir] - 0.5*h[dir]) : 
			(U[dir] + 0.5*h[dir]);

	comp = component(coords,intfc);
	if(comp == NO_COMP)
	    return intfc->default_comp;
	else
	    return comp;

}	/* end buffer_component */

LOCAL 	double *constr_position(
	double *p1,
	double *p2,
	boolean (*constr_func)(POINTER,double*), // Constraint function
        POINTER func_params)		     // Constraint function params
{
	static double p_in[MAXD];
	double p_out[MAXD],pc[MAXD];
	int i,j,N = 20;		// accurate to 2^-20*(p1-p2)
	
	if (constr_func(func_params,p1) && !constr_func(func_params,p2))
	{
	    for (i = 0; i < MAXD; ++i)
	    {
		p_in[i]  = p1[i];
		p_out[i] = p2[i];
	    }
	}
	else if (constr_func(func_params,p2) && !constr_func(func_params,p1))
	{
	    for (i = 0; i < MAXD; ++i)
	    {
		p_in[i]  = p2[i];
		p_out[i] = p1[i];
	    }
	}
	else
	{
	    screen("ERROR: In constr_position(), same side!\n");
	    clean_up(ERROR);
	}
	for (i = 0; i < N; ++i)
	{
	    for (j = 0; j < MAXD; ++j)
		pc[j] = 0.5*(p_in[j] + p_out[j]);
	    if (constr_func(func_params,pc))
		for (j = 0; j < MAXD; ++j)
		    p_in[j] = pc[j];
	    else
		for (j = 0; j < MAXD; ++j)
		    p_out[j] = pc[j];
	}
	return p_in;
}	/* end constr_position */

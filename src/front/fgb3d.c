/*********************************************************************
FronTier is a set of libraries that implements differnt types of 
Front Traking algorithms. Front Tracking is a numerical method 
for the solution of partial differential equations whose solutions 
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
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

************************************************************************/


/*
*                               tricrx.c:
*
*       Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*/


#define DEBUG_STRING "crx_intfc"
#define DB_TEC  true

#include <front/fdecs.h>
#include <front/fpatrecon.h>

LOCAL	boolean track_comp_and_repair3d(int*,int*,int*,INTERFACE*,Front*);
LOCAL	boolean track_comp_and_set_boxes(int*,int*,int*,INTERFACE*,Front*);
	boolean	communicate_boxes(RECT_BOX*,Front*);
LOCAL	boolean	repair_intfc3d_in_box(INTERFACE*,int*,int*,Front*);
LOCAL	boolean set_reconstruction_boxes(int*,int*,int**,int,RECT_BOX**,
                                        INTERFACE*);
LOCAL	boolean seal_strip_with_tris(RECT_BOX*,TRI**,int*,int,POINT**,
				double*,char*,boolean,INTERFACE*);
LOCAL   boolean seal_strip_with_tris_once(RECT_BOX*,TRI**,int*,int,POINT**,
                          double*,char*,boolean,INTERFACE*,TRI**,int,TRI**,int);

LOCAL   boolean null_sides_sharing_same_vertex(TRI*,int,TRI**,int);
	boolean null_side_loop(TRI*,int,ORIENTATION,TRI***,int**,POINT***,
                                int*,double**);
LOCAL	boolean null_side_tri_in_list(TRI**,int,TRI**,int*);
LOCAL	boolean grid_based_box_untangle(INTERFACE*,RECT_BOX*,boolean);
LOCAL   boolean tri_in_box(TRI*,INTERFACE*,RECT_BOX*);
LOCAL	boolean check_extension_of_surface(TRI**, int, SURFACE*);
LOCAL	void 	gview_show_box_tri(RECT_BOX*,TRI**,int,FILE*);
LOCAL   int     find_nearest_pts(POINT**,int,POINT**,int,int,TRI**,TRI**,
				RECT_GRID,RECT_GRID);
LOCAL   int     find_nearest_pts_pair(POINT**,int,POINT**,int,int,RECT_GRID,
				RECT_GRID);
LOCAL	void 	find_nearest_tri_pair(TRI**,int,TRI**,int,TRI**,
			int*,TRI**,int*,RECT_BOX*);
LOCAL   void    find_nearest_tri_pair_new(TRI**,int*,TRI**,int*,TRI**,
                int*,TRI**,int*,TRI**,int*,TRI**,RECT_BOX*,INTERFACE*);
LOCAL   void    find_nearest_tri_pair_crx(TRI**,int*,TRI**,int*,TRI**,
                int*,TRI**,int*,TRI**,int*,TRI**,RECT_BOX*,INTERFACE*);
LOCAL	boolean	null_sides_with_suitable_angle(TRI*, int, TRI*, int);
LOCAL	boolean	bifurcation_detected(TRI**, int, TRI***, int*, TRI***, int*);
LOCAL   boolean pts_in_cross_zone(double*,double,double,int);
LOCAL   boolean pts_cross_line(double*,double,double,int);
LOCAL   boolean tri_cross_zone(TRI*,double,double,int);
LOCAL   boolean tri_cross_line1(TRI*,double,double,int);
LOCAL   void    find_comm_box(INTERFACE*,RECT_BOX*,int*,COMM_BOX*);
LOCAL   void    initialize_comm_box(COMM_BOX*,boolean);
LOCAL   void    pp_send_box(PP_GRID*,int*,int,int,COMM_BOX*,INTERFACE*);
LOCAL   void    pp_receive_box(PP_GRID*,int*,int,int,COMM_BOX*,INTERFACE*);
LOCAL   boolean compare_comm_box(INTERFACE*,RECT_BOX*,
				COMM_BOX*,COMM_BOX*,int*,int);
LOCAL   boolean boxes_sharing_tris(RECT_BOX*,RECT_BOX*, int*, INTERFACE*);
LOCAL   boolean boxes_crxing_same_boundary(RECT_BOX*,RECT_BOX*,INTERFACE*);
LOCAL   void    extend_boundary_side(INTERFACE*,int*,int*,RECT_GRID*);
LOCAL   void    check_and_extend_point(POINT*,TRI*,SURFACE*,double*,
                                double*,double*);
LOCAL   void    check_and_extend_side(TRI*,int,SURFACE*,double*,
                                double*,double*);
LOCAL   boolean crx_bnd_out_tris(TRI**,int,TRI**,int,TRI*,INTERFACE*,RECT_BOX*);
LOCAL   boolean crx_out_tris_twice(TRI**,int,TRI**,int,TRI*,INTERFACE*,
				RECT_BOX*);
LOCAL	boolean ip_connected(int **, int, int*);
LOCAL	boolean overlapping_boxes(RECT_BOX*, RECT_BOX*);
LOCAL	void 	tecplot_show_box(char*,RECT_BOX*,FILE*);
LOCAL	void	tecplot_show_null_tris(TRI**, POINT**, int, FILE*);
LOCAL	void 	remove_crx_tri_on_edge(TRI*,TRI*,int,TRI**,int*);
LOCAL	boolean	check_adjecant_constrain0(POINT*,POINT*,POINT*,POINT*,
				TRI**,int);
LOCAL	void	seal_null_loop(TRI**,int*,POINT**,int);
LOCAL	boolean rm_bad_crxs_in_box(int*,int*,int**,int,RECT_BOX*,INTERFACE*);
LOCAL	void  	smooth_null_loop(TRI**,int*,POINT**,int);
LOCAL	boolean set_tst_recon_boxes(int*,int*,RECT_BOX**,INTERFACE*);
	void 	smooth_tris(TRI **,int);

boolean    set_use_rect_tris(boolean);

int	neighbor_tri_side(TRI*, TRI*);
void    merge_near_tri_pairs(int*,int*,INTERFACE*);
boolean	connect_tris_holes(TRI*,TRI*,INTERFACE*,char*);
boolean	merge_near_interface(Front*);

EXPORT	boolean  skip_tag_tri(TRI*);

#define		MAX_NULL_SIDE_LOOP		1000


struct _FBOX {
	double		f[2][3];
	struct _FBOX	*prev, *next;
};

typedef struct _FBOX FBOX;
	
	boolean	box_sect(
	double	fbox1[2][3],
	double	fbox2[2][3],
	double	tol)
{
	int	i;

	for(i=0; i<3; i++)
	    if(fbox1[1][i] + tol < fbox2[0][i] || 
	       fbox2[1][i] + tol < fbox1[0][i])
		return NO;

	return YES;
}

	void	box_merge(
	double	fbox1[2][3],
	double	fbox2[2][3])
{
	int	i;

	for(i=0; i<3; i++)
	{
	    fbox1[0][i] = min(fbox1[0][i], fbox2[0][i]);
	    fbox1[1][i] = max(fbox1[1][i], fbox2[1][i]);
	}
}

	void	box_of_tri(
	double	fbox[2][3],
	TRI	*t)
{
	int	i, j;
	double	*p;

	for(j=0; j<3; j++)
	{
	    fbox[0][j] = HUGE_VAL;
	    fbox[1][j] = -HUGE_VAL;
	}

	for(i=0; i<3; i++)
	{
	    p = Coords(Point_of_tri(t)[i]);
	    for(j=0; j<3; j++)
	    {
		fbox[0][j] = min(fbox[0][j], p[j]);
		fbox[1][j] = max(fbox[1][j], p[j]);
	    }
	}
}

void	show_tris_bound_tris(char*,char*,TRI**,int,INTERFACE*);

#define	MAX_SMOOTH_TRIS	2000

	void	smooth_tris_in_boxes(
	int		*smin,
	int		*smax,
	int		**ips,
	int		num_ip,
	INTERFACE	*intfc)
{
	RECT_BOX	*boxes, *pb;
	FBOX		Fbox, Box, *box, *nbox, *fboxes;
	int		i,j,k,num_box, nproc, num, nt;
	static	FBOX	*fbox = NULL;
	double		tol;
	boolean		box_merged;
	TRI		*tris[MAX_SMOOTH_TRIS], *tri;
	SURFACE		**s;
	RECT_GRID	*gr = &topological_grid(intfc);

	DEBUG_ENTER(smooth_tris_in_boxes)

	nproc = pp_numnodes();
	if(fbox == NULL)
	    uni_array(&fbox, nproc, sizeof(FBOX));

	/* tst lgb */
	if(num_ip != 0)
	{
	    set_reconstruction_boxes(smin,smax,ips,num_ip,&boxes,intfc);
	    num_box = 0;
	    for(pb=boxes; pb!=NULL; pb=pb->next)
		num_box++;
	    pb = boxes;
	}
	else
	{
	    num_box = 0;
	    pb = NULL;
	}

	pp_global_imax(&num_box, 1L);
	printf("#max num_box = %d\n", num_box);
	if(num_box == 0)
	{
	    DEBUG_LEAVE(smooth_tris_in_boxes)
	    return;
	}

	Box.prev = Box.next = NULL;
	box = &Box;
	num = 0;
	for(i=0; i<num_box; i++)
	{
	    if(pb != NULL)
	    {
		for(j=0; j<3; j++)
		{
		    Fbox.f[0][j] = gr->L[j] + gr->h[j]*pb->bmin[j];
		    Fbox.f[1][j] = gr->L[j] + gr->h[j]*pb->bmax[j];
		}
		pb = pb->next;
	    }
	    else
	    {
		for(j=0; j<3; j++)
		{
		    Fbox.f[0][j] = 1.0;
		    Fbox.f[1][j] = -1.0;
		}
	    }
	    
	    pp_all_gather((POINTER)&Fbox,sizeof(FBOX),(POINTER)fbox,sizeof(FBOX));
	    
	    for(j=0; j<nproc; j++)
	    {
		if(fbox[j].f[0][0] >= fbox[j].f[1][0])
		    continue;
		
		num++;
		box->next = (FBOX *)store(sizeof(FBOX));
		box->next->prev = box;
		box->next->next = NULL;
		box = box->next;

		for(k=0; k<3; k++)
		{
		    box->f[0][k] = fbox[j].f[0][k];
		    box->f[1][k] = fbox[j].f[1][k];
		}
	    }
	}

	printf("#total num box = %d\n", num);

	tol = 1.0e-6*min3(gr->h[0], gr->h[1], gr->h[2]);
	fboxes = Box.next;
	box_merged = YES;
	
	while(box_merged)
	{
	    box_merged = NO;
	    for(box=fboxes; box!=NULL; box=box->next)
	    {
		for(nbox=box->next; nbox!=NULL; nbox=nbox->next)
		{
		    if(box_sect(box->f, nbox->f, tol))
		    {
			/* merge nbox with box and delete nbox. */
			box_merge(box->f, nbox->f);
			nbox->prev->next = nbox->next;
			if (nbox->next != NULL)
			    nbox->next->prev = nbox->prev;
			
			box_merged = YES;
			break;
		    }
		}   /* for(nbox= ) */
		if(box_merged)
		    break;
	    }   /* for(box=  ) */
	}

	for(i=0; i<3; i++)
	{
	    Fbox.f[0][i] = gr->L[i];
	    Fbox.f[1][i] = gr->U[i];
	}

	fboxes = Box.next;
	/* remove outside domain boxes */
	for(box=fboxes; box!=NULL; box=box->next)
	{
	    if(!box_sect(box->f, Fbox.f, tol))
	    {
		box->prev->next = box->next;
		if(box->next != NULL)
		    box->next->prev = box->prev;
		box = box->prev;
	    }
	}

	fboxes = Box.next;
	num = 0;
	for(box=fboxes; box!=NULL; box=box->next)
	{
	    print_general_vector("l", box->f[0], 3, "\n"); 
	    print_general_vector("r", box->f[1], 3, "\n");

	    printf("#size %d: ", num);
	    num++;
	    for(i=0; i<3; i++)
		printf("%3d  ", (int)((box->f[1][i]-box->f[0][i] + tol)/gr->h[i]));
	    printf("\n");
	}

	printf("#num inside %d\n", num);

	print_rectangular_grid(gr);

	nt = 0;
	for(s=intfc->surfaces; s && *s; s++)
	{
	    for(tri=first_tri(*s); !at_end_of_tri_list(tri,*s); tri=tri->next)
	    {
		box_of_tri(Fbox.f, tri);
		for(box=fboxes; box!=NULL; box=box->next)
		{
		    if(box_sect(box->f, Fbox.f, tol))
		    {
			tris[nt] = tri;
			nt++;
			if(nt > MAX_SMOOTH_TRIS)
			{
			    printf("ERROR smooth_tris_in_boxes, too many tris.\n");
			    clean_up(ERROR);
			}
		    }
		}
	    }
	}

if(debugging("smotribox"))
{
	if(nt > 0)
	{
	    char	fname[500];
	    FILE	*fp;
	    
	    add_to_debug("box_intfc");
	    sprintf(fname, "box_tris_%d.plt", pp_mynode());
	    
	    printf("#show tris in box %s\n", fname);
	    
	    fp = fopen(fname, "w");
	    fprintf(fp,"TITLE = \"tecplot tris\"\n"
		     "VARIABLES = \"x\", \"y\", \"z\"\n");

	    tecplot_box(NULL, fp, gr->L, gr->U);
	    for(box=fboxes; box!=NULL; box=box->next)
	    {
		tecplot_box(NULL, fp, box->f[0], box->f[1]);
	    }
	    fclose(fp);

	    show_tris_bound_tris("con", fname, tris, nt, intfc);

	    smooth_tris(tris, nt);
	    
	    show_tris_bound_tris("smo_con", fname, tris, nt, intfc);
	    
	    remove_from_debug("box_intfc");
	}
}
else
{
	if(nt > 0)
	    smooth_tris(tris, nt);
}

	DEBUG_LEAVE(smooth_tris_in_boxes)
}

EXPORT	boolean	merge_near_interface(
	Front     *front)
{
	INTERFACE *intfc = front->interf;
	boolean      sav_intrp = interpolate_intfc_states(intfc);
	int 	  i, smin[3], smax[3], *gmax;
	int	  n_reg_node, num_ip;
	size_t	  sizest = front->sizest;
	RECT_GRID Dual_grid, *comp_grid = front->rect_grid;
	RECT_GRID gr_save = topological_grid(intfc), *gr;
	Table 	  *T;
	COMPONENT *comp;
	static int 	**ips = NULL;

	DEBUG_ENTER(merge_near_interface)
	
	if(ips == NULL)
	    stat_matrix(&ips,MAX_NUM_UNPHY_IP,3,INT);

	set_dual_grid(&Dual_grid,comp_grid);
	set_size_of_intfc_state(sizest);
	set_copy_intfc_states(YES);

	gr = &topological_grid(intfc);
	set_expanded_grid(&Dual_grid, gr);

	make_interface_topology_lists(intfc);
	set_crx_storage_for_reconstruction(intfc,NULL);

	/* set reconstruction boundary and tolerance */

	gmax = gr->gmax;
	for (i = 0; i < 3; ++i)
	{
	    smin[i] = 0;
            smax[i] = gmax[i];
	}
	
	start_clock("insert_grid_crossings3d");
	
	interpolate_intfc_states(intfc) = YES;
	insert_grid_intfc_crossings(intfc);

	stop_clock("insert_grid_crossings3d");
	
	T = table_of_interface(intfc);
	comp = T->components;

	n_reg_node = (gmax[0]+1)*(gmax[1]+1)*(gmax[2]+1);
	for (i = 0; i < n_reg_node; ++i)
	    comp[i] = NO_COMP;
	
	adjust_crossings(smin,smax,intfc);
	fill_comp_from_prev_intfc(intfc, smin, smax);
	
	fill_physical_comps(smin, smax, gmax, intfc);

	fill_comp_with_component3d(smin,smax,gmax,intfc);
	
	if(front->step == 3780)
	    add_to_debug("dbmgtris");

	num_ip = record_unphysical_ips(smin,smax,intfc,ips);
	
	set_use_rect_tris(NO);
	smooth_tris_in_boxes(smin, smax, ips, num_ip, intfc);

	set_use_rect_tris(YES);

	interpolate_intfc_states(intfc) = sav_intrp;

	reset_intfc_num_points(intfc);
	topological_grid(intfc) = gr_save;
        free_crx_storage(intfc);
	
	front->interf = copy_interface(intfc);
	
	delete_interface(intfc);
	
	DEBUG_LEAVE(merge_near_interface)
	
	return YES;
}

void	reset_tri_order(
	INTERFACE	*intfc)
{
	SURFACE		**s;
	TRI		*tri;

	for(s=intfc->surfaces; s && *s; ++s)
	    for(tri=first_tri(*s); !at_end_of_tri_list(tri,*s); tri=tri->next)
	    {
		Tri_order(tri) = 0;
	    }
}

#define	MAX_MERGE_CRX	5000

typedef  struct {
	CRXING		*crx1, *crx2;
	double		dist;
}    CRX_PAIR;

/*dist decrease order.*/
int compare_crx_pair(const void *a, const void *b)
{
	CRX_PAIR	*c1=(CRX_PAIR*)a, *c2=(CRX_PAIR*)b;
	double	 	dist1 = c1->dist, dist2 = c2->dist;

	return dist1 < dist2 ? -1 : 1;
}

boolean	check_two_tris_cond(TRI*,TRI*,INTERFACE*);
	
	void merge_near_tri_pairs(
	int      *smin,
        int      *smax,
	INTERFACE *intfc)
{
	RECT_GRID      	*gr = &topological_grid(intfc);
	GRID_DIRECTION	dir[3] = {EAST,NORTH,UPPER};
	CRXING 	       	*crx1, *crx2;
	int 	       	ix,iy,iz, ip[3], ipn[3];
	int 	       	i,k,j,num,nc,list1;
	double		minh, merge_tol, dist;
	double	       	*L = gr->L, *h = gr->h;
	int		*gmax = gr->gmax;
	Table		*T = table_of_interface(intfc);
	COMPONENT 	*comp = T->components;
	CRX_PAIR	crx_pairs[MAX_MERGE_CRX];
	boolean		status;
	char		fname[500];

	minh = min3(h[0], h[1], h[2]);
	merge_tol = 0.4*minh;
	num = 0;

	for (iz = smin[2]; iz <= smax[2]; ++iz)
	{
	    for (iy = smin[1]; iy <= smax[1]; ++iy)
	    {
	        for (ix = smin[0]; ix <= smax[0]; ++ix)
	        {
	            for (i = 0; i < 3; ++i)
	            {
	                if (ix == smax[0] && dir[i] == EAST)
	                    continue;
	                if (iy == smax[1] && dir[i] == NORTH)
	                    continue;
	                if (iz == smax[2] && dir[i] == UPPER)
	                    continue;

	                ip[0] = ix;
	                ip[1] = iy;
	                ip[2] = iz;
			next_ip_in_dir(ip,dir[i],ipn,smin,smax);

			/*no unphysical crxings on the edge.*/
			if(comp[d_index3d(ip[0],ip[1],ip[2],gmax)] == NO_COMP || 
			   comp[d_index3d(ipn[0],ipn[1],ipn[2],gmax)] == NO_COMP)
			    continue;
			
			k = seg_index3d(ix,iy,iz,dir[i],gmax);
	                nc = T->seg_crx_count[k];

			if(nc <= 1)
			    continue;

			for(j=0; j<nc-1; j++)
			{
			    list1 = T->seg_crx_lists[k][j];
	                    crx1 = T->crx_store+list1;
			    
			    list1 = T->seg_crx_lists[k][j+1];
	                    crx2 = T->crx_store+list1;
			    
			    dist = distance_between_positions(Coords(crx1->pt), 
			    	Coords(crx2->pt), 3);
			    
	    		    if(!check_two_tris_cond(crx1->tri, crx2->tri, intfc))
				continue;

			    if(dist < merge_tol)
			    {
				crx_pairs[num].crx1 = crx1;
				crx_pairs[num].crx2 = crx2;
				crx_pairs[num].dist = dist;
				num++;
				if(num >= MAX_MERGE_CRX)
				{
				    printf("ERROR merge_near_interface"
				    	   "num is too large.\n");
				    clean_up(ERROR);
				}
			    }
			} /* for j: crx index */
	            } /* for i: directions */
	        }
	    }
	}

	printf("#merge num = %d\n", num);

	if(num == 0)
	    return;

	qsort((POINTER)crx_pairs, num, sizeof(CRX_PAIR), compare_crx_pair);
	
	if(NO)
	{
	    char	fname[500];
	    FILE	*fp;
	    TRI		**tris;

	    uni_array(&tris, num*2, sizeof(TRI*));
	    sprintf(fname, "merge_tris_%d.plt", pp_mynode());
	   
	    for(i=0; i<num; i++)
	    {
		crx1 = crx_pairs[i].crx1;
		crx2 = crx_pairs[i].crx2;

		tris[i*2] = crx1->tri;
		tris[i*2+1] = crx2->tri;
	    }

	    fp = fopen(fname, "w");
	    tecplot_show_tris("merge", tris, num*2, fp);
	    fclose(fp);

	    free(tris);
	}

	reset_tri_order(intfc);

	sprintf(fname, "conn_tris_%d.plt", pp_mynode());

	j = 0;
	for(i=0; i<num; i++)
	{
	    crx1 = crx_pairs[i].crx1;
	    crx2 = crx_pairs[i].crx2;
	    printf("#crx_pairs  %5d  %15.8e\n", i, crx_pairs[i].dist/minh);
	    status = connect_tris_holes(crx1->tri, crx2->tri, intfc, fname);
	    if(status)
		j++;
	}

	printf("#connect pairs %d\n", j);
}


/*
*		repair_intfc_at_crossings3d():
*/

EXPORT	boolean	repair_intfc_at_crossings3d(
	Front     *front)
{
	INTERFACE *intfc = front->interf;
	boolean      sav_intrp = interpolate_intfc_states(intfc);
	int 	  i, smin[3], smax[3];
	size_t	  sizest = front->sizest;
	RECT_GRID Dual_grid, *comp_grid = front->rect_grid;
	RECT_GRID gr_save = topological_grid(intfc);

	DEBUG_ENTER(repair_intfc_at_crossings3d)

	set_dual_grid(&Dual_grid,comp_grid);
	set_size_of_intfc_state(sizest);
	set_copy_intfc_states(YES);

	print_storage("Entering repair_intfc_at_crossings3d","crx_store");

	set_expanded_grid(&Dual_grid,&topological_grid(intfc));

	make_interface_topology_lists(intfc);
	set_crx_storage_for_reconstruction(intfc,NULL);

	/* set reconstruction boundary and tolerance */

	for (i = 0; i < 3; ++i)
	{
	    smin[i] = 0;
            smax[i] = topological_grid(intfc).gmax[i];
	}
	
	start_clock("insert_grid_crossings3d");
	interpolate_intfc_states(intfc) = YES;
	insert_grid_intfc_crossings(intfc);

	stop_clock("insert_grid_crossings3d");
	print_storage("After insert_grid_crossings3d","crx_store");

	start_clock("reconstruct_intfc3d_in_box");
	if (!repair_intfc3d_in_box(intfc,smin,smax,front))
	{
	    interpolate_intfc_states(intfc) = sav_intrp;
	    DEBUG_LEAVE(repair_intfc_at_crossings3d)
	    return NO;
	}
	stop_clock("reconstruct_intfc3d_in_box");
	print_storage("After reconstruct_intfc3d_in_box","crx_store");

	interpolate_intfc_states(intfc) = sav_intrp;

	reset_intfc_num_points(intfc);
	topological_grid(intfc) = gr_save;
        free_crx_storage(intfc);
	
	set_size_of_intfc_state(sizest);
	set_copy_intfc_states(YES);
	front->interf = copy_interface(intfc);
	
	delete_interface(intfc);
	print_storage("After copy_interface","crx_store");
	DEBUG_LEAVE(repair_intfc_at_crossings3d)
	return YES;
}	/*end repair_intfc_at_crossings3d*/

/*
*		rebuild_intfc_at_crossings3d():
*/
EXPORT	boolean	rebuild_intfc_at_crossings3d(
	Front     *front)
{
	if(Tracking_algorithm(front) == THREE_COMP_GRID_BASED_TRACKING)
	    return rebuild_intfc_at_crossings3d3(front);
	else if(Tracking_algorithm(front) == GRID_BASED_TRACKING)
	    return rebuild_intfc_at_crossings3d2(front);
	else
	{
	    /* LGB and GB mixed prop */
	    printf("WARNING: rebuild_intfc_at_crossings3d, call 3 "
		   "comp recon alg.\n");
	    return rebuild_intfc_at_crossings3d3(front);
	}
}

/* 2 comp grid based track */

EXPORT	boolean	rebuild_intfc_at_crossings3d2(
	Front     *front)
{
	INTERFACE *intfc = front->interf;
	boolean   	  sav_intrp = interpolate_intfc_states(intfc);
	int 	  i, smin[3], smax[3];
	size_t	  sizest = front->sizest;
	RECT_GRID Dual_grid, *comp_grid = front->rect_grid;
	RECT_GRID gr_save = topological_grid(intfc);
	RECT_GRID *gr;

	DEBUG_ENTER(rebuild_intfc_at_crossings3d)

	set_dual_grid(&Dual_grid,comp_grid);
	set_size_of_intfc_state(sizest);
	set_copy_intfc_states(YES);

	print_storage("Entering rebuild_intfc_at_crossings3d","crx_store");

	set_expanded_grid(&Dual_grid,&topological_grid(intfc));
	gr = &topological_grid(intfc);

	/* set reconstruction boundary and tolerance */

	for (i = 0; i < 3; ++i)
	{
	    smin[i] = 0;	
	    smax[i] = topological_grid(intfc).gmax[i]; /*expanded_dual_grid*/
            if (Dual_grid.lbuf[i] != 0)
		smin[i] += Dual_grid.lbuf[i] - 1;
	    else if (curves_on_bdry_side(i,0,intfc) == YES)
	        smin[i] += 1;
            if (Dual_grid.ubuf[i] != 0)
		smax[i] -= Dual_grid.ubuf[i] - 1;
	    else if (curves_on_bdry_side(i,1,intfc) == YES)
	        smax[i] -= 1;
	}
	
	extend_boundary_side(intfc,smin,smax,gr);

	make_interface_topology_lists(intfc);
        set_crx_storage_for_reconstruction(intfc,NULL);

	start_clock("insert_grid_crossings3d");
	interpolate_intfc_states(intfc) = YES;
	
	insert_grid_intfc_crossings(intfc);
	stop_clock("insert_grid_crossings3d");
	print_storage("After insert_grid_crossings3d","crx_store");

	start_clock("remove_unphysical_crossings3d");
	if (!remove_unphysical_crossings3d(intfc,smin,smax))
	{
	    interpolate_intfc_states(intfc) = sav_intrp;
	    DEBUG_LEAVE(rebuild_intfc_at_crossings3d)
	    return NO;
	}
	stop_clock("remove_unphysical_crossings3d");
	strip_subdomain_bdry_curves(intfc);

	start_clock("reconstruct_intfc3d_in_box");
	if (!reconstruct_intfc3d_in_box(intfc,smin,smax,YES,NULL))
	{
	    interpolate_intfc_states(intfc) = sav_intrp;
	    DEBUG_LEAVE(rebuild_intfc_at_crossings3d)
	    return NO;
	}
	stop_clock("reconstruct_intfc3d_in_box");
	print_storage("After reconstruct_intfc3d_in_box","crx_store");

	interpolate_intfc_states(intfc) = NO;
	install_subdomain_bdry_curves(intfc);
	interpolate_intfc_states(intfc) = sav_intrp;

	identify_detached_surface_curve_pair(intfc);

	reset_intfc_num_points(intfc);
	topological_grid(intfc) = gr_save;
        free_crx_storage(intfc);
	front->interf = copy_interface(intfc);
	delete_interface(intfc);

	interface_reconstructed(front->interf) = YES;
	print_storage("After copy_interface","crx_store");
	DEBUG_LEAVE(rebuild_intfc_at_crossings3d)
	return YES;
}	/*end rebuild_intfc_at_crossings3d*/

LOCAL	boolean repair_intfc3d_in_box(
	INTERFACE  *intfc,
	int *smin,
	int *smax,
	Front *front)
{
	Table *T = table_of_interface(intfc);
	RECT_GRID gr = topological_grid(intfc);
	COMPONENT *comp = T->components;
	int i,n_reg_node;
	int *gmax = gr.gmax;
	boolean status;

	if (debugging("set_component"))
	    (void) printf("Entering repair_intfc3d_in_box()\n");

	n_reg_node = (gmax[0]+1)*(gmax[1]+1)*(gmax[2]+1);
	for (i = 0; i < n_reg_node; ++i)
	    comp[i] = NO_COMP;

	/* new func */
	status = track_comp_and_repair3d(smin,smax,gmax,intfc,front);

	if (status == FUNCTION_FAILED)
	{
	    screen("track_comp_and_repair3d() failed!\n");
	    clean_up(ERROR);
	    return status;
	}
	if (debugging("set_component"))
	    (void) printf("Leaving repair_intfc3d_in_box()\n");
	return status;
}	/* end repair_intfc3d_in_box */



LOCAL	int box_index;
LOCAL   int test_tri = -10;
LOCAL   double cosu = HUGE;
LOCAL   int in_tri_recorded[MAX_NUM_UNPHY_IP];
LOCAL   int out_tri_recorded[MAX_NUM_UNPHY_IP];
LOCAL   int num_tri_recorded;

void   print_edge_crossings(int *, int*, INTERFACE *);
LOCAL	boolean set_tst_recon_boxes(int*,int*,RECT_BOX**,INTERFACE*);
LOCAL	void communicate_box_recon(RECT_BOX*,INTERFACE*, Front*);

LOCAL	boolean track_comp_and_repair3d(
	int		*smin,
	int		*smax,
        int		*gmax,
        INTERFACE	*intfc,
	Front		*fr)
{
	int		num_ip, i;
	RECT_BOX 	*boxes, *pb;
	static int 	**ips = NULL;
	RECT_GRID	*gr = &topological_grid(intfc);

	if (debugging("trace"))
	    (void) printf("Entering track_comp_and_repair3d()\n");
	if(ips == NULL)
	    stat_matrix(&ips,MAX_NUM_UNPHY_IP,3,INT);

	adjust_crossings(smin,smax,intfc);
	fill_comp_from_prev_intfc(intfc, smin, smax);
	
	fill_physical_comps(smin,smax,gmax,intfc);

	remove_unphysical_crxings(smin,smax,gmax,intfc,SINGLE,&num_ip,ips);

	if (debugging("compcrx"))
	{
	    printf("After remove_up\n");
	    printf("#num_ip  %d\n", num_ip);
	    for (i = 0; i < num_ip; ++i)
		printf("ip[%d] = %d %d %d\n",i,ips[i][0],ips[i][1],ips[i][2]);
	}
	
	if(!set_reconstruction_boxes(smin,smax,ips,num_ip,&boxes,intfc))
	{
	    if (debugging("trace"))
	    {
		(void) printf("Function failed!\n");
	    	(void) printf("Leaving track_comp_and_repair3d()\n");
	    }
	    return FUNCTION_FAILED;
	}

	communicate_box_recon(boxes, intfc, fr);
	
	if (debugging("trace"))
	    (void) printf("Leaving track_comp_and_repair3d()\n");
	return FUNCTION_SUCCEEDED;

}	/* end track_comp_and_repair3d */


boolean	rbox_communication_boxes(RECT_BOX*,Front*);

LOCAL void communicate_box_recon(
	RECT_BOX	*boxes,
	INTERFACE	*intfc,
	Front		*fr)
{
	RECT_BOX	*pb;

	set_current_interface(intfc);
	
	if(debugging("newrecon"))
	{
	  rbox_communication_boxes(boxes, fr);
	  return;
	}

	box_index = 0;
	
	for (pb = boxes; pb != NULL; pb = pb->next)
	{
	    RECT_GRID	*gr = &topological_grid(intfc);
	    double	*L = gr->L, *h = gr->h;

	    if (debugging("box_intfc"))
	    {
	    	printf("#untangle box_index = %d\n", box_index);

	    	printf("bmin = %-2d  %-2d  %-2d\n",pb->bmin[0],
		   pb->bmin[1],pb->bmin[2]);

	    	printf("bmax = %-2d  %-2d  %-2d\n",pb->bmax[0],
		   pb->bmax[1],pb->bmax[2]);

	    	printf("bmin = %15.8e, %15.8e, %15.8e\n", 
		    L[0] + h[0]*pb->bmin[0], L[1] + h[1]*pb->bmin[1], 
		    L[2] + h[2]*pb->bmin[2]);
	    
	    	printf("bmax = %15.8e, %15.8e, %15.8e\n", 
		    L[0] + h[0]*pb->bmax[0], L[1] + h[1]*pb->bmax[1], 
		    L[2] + h[2]*pb->bmax[2]);
	    }

	    if (!grid_based_box_untangle(intfc,pb, NO))
	    {
	    	printf("WARNING, track_comp_and_repair3d, "
		       "grid_based_box_untangle fails. \n");
		DEBUG_LEAVE(track_comp_and_repair3d)
	        return;
	    }
	    
	    box_index++;
	}
	
	if (debugging("box_intfc"))
	{
	    remove_from_debug("box_intfc");
	    if (consistent_interface(intfc))
	    {
	    	(void) printf("After untangle, interface is consistent!\n");
	    }
	}
}

LOCAL	boolean track_comp_and_set_boxes(
	int		*smin,
	int		*smax,
        int		*gmax,
        INTERFACE	*intfc,
	Front		*fr)
{
	int		num_ip, i;
	RECT_BOX 	*boxes, *pb;
	static int 	**ips = NULL;
	RECT_GRID	*gr = &topological_grid(intfc);
	double		*L = gr->L, *h = gr->h;

	DEBUG_ENTER(track_comp_and_set_boxes)

	if(ips == NULL)
	    stat_matrix(&ips,MAX_NUM_UNPHY_IP,3,INT);

	print_rectangular_grid(gr);

	/* see function track_comp_through_crxings3d */
	
	adjust_crossings(smin,smax,intfc);
	fill_comp_from_prev_intfc(intfc, smin, smax);
	
	fill_physical_comps(smin,smax,gmax,intfc);

	fill_comp_with_component3d(smin,smax,gmax,intfc);

	if(debugging("compcrx"))
	{
	    int  tmin[3] = {0, 12, 20}, tmax[3] = {2, 14, 22};
	    
	    printf("#check comp physical\n");

	    show_grid_components(tmin,tmax,0,intfc);
	    show_grid_components(tmin,tmax,1,intfc);
	    show_grid_components(tmin,tmax,2,intfc);
	}

	/* record all unphysical ip's */
	num_ip = record_unphysical_ips(smin,smax,intfc,ips);

	if (debugging("compcrx"))
	    printf("#num_ip  %d\n", num_ip);

	/* make the boxes for lgb reconstruction. */
	if(!set_reconstruction_boxes(smin,smax,ips,num_ip,&boxes,intfc))
	{
	    DEBUG_LEAVE(track_comp_and_set_boxes)
	    return FUNCTION_FAILED;
	}
	
	if(!communicate_boxes(boxes, fr))
	{
	    DEBUG_LEAVE(track_comp_and_set_boxes)
	    return FUNCTION_FAILED;
	}

	clean_up(0);

	DEBUG_LEAVE(track_comp_and_set_boxes)
	return FUNCTION_SUCCEEDED;

}	/* end track_comp_and_repair3d */

LOCAL   int  tmp_n;
LOCAL   TRI  *tmp_tris[200];

EXPORT	boolean	check_normal_on_intfc(
	INTERFACE	*intfc)
{
	const double	*nor;
	SURFACE		**s;
	TRI		*t;

	for (s = intfc->surfaces; s && *s; ++s)
	    for (t = first_tri(*s); !at_end_of_tri_list(t,*s); t = t->next)
	    {
		nor = Tri_normal(t);
		if(isnan(Mag3d(nor)) || 
		   (nor[0] == 0.0 && nor[1] == 0.0 && nor[2] == 0.0))
		{
		    printf("ERROR check_normal_on_intfc,  "
		    	   "bad normal appears.\n");
		    print_tri(t, intfc);
		    clean_up(ERROR);
		}
	    }

	return YES;
}

void	tecplot_debug_tris(const char*,char*,TRI**,int,INTERFACE*);

	void	tecplot_debug_tris(
	const char	*tname,
	char		*fname,
	TRI		**tris,
	int		num_tris,
	INTERFACE	*intfc)
{
	FILE	*file;

	if(!debugging("box_intfc") && !debugging("dbmgtris"))
	    return;
	/*
	if(!debugging("no_valid_test"))
	{
	    printf("#check valid tris %s \n", tname);
	    if(!check_valid_tris(tris, num_tris, intfc))
	    {
		printf("#check valid tris, invalid tris %s\n",
		       tname);
	    }
	}
	*/

	file = fopen(fname,"a");
	tecplot_show_tris(tname, tris, num_tris, file);
	fclose(file);
}

	boolean	set_tri_area_tol(double);

LOCAL	boolean grid_based_box_untangle(
	INTERFACE	*intfc,
	RECT_BOX 	*box,
	boolean		tri_tag)
{
	SURFACE	**surfs,*ref_surfs[100];
	RECT_GRID *gr;
	int 	*gmax,*smin,*smax;
	TRI	*tri, *last_tris[100];
	TRI 	**test_tris,**in_tris,**out_tris;
	TRI	**ref_tris, **new_tris, **deg_tris, **sep_new_tris;
	int	num_test_tris, num_ref_tris, num_in_tris, num_seal_tris;
	int	num_out_tris, num_deg_tris, num_new_tris, max_n_new;
	int 	i,j,k,num_surfs,total_nt, kmin[3], kmax[3];
	double	hmin;
	FILE 	*file;
	char 	dname[100],fname[100];
	
	DEBUG_ENTER(grid_based_box_untangle)
	
	gr = box->grid;
	smin = box->bmin;
	smax = box->bmax;
	gmax = gr->gmax;
	hmin = min3(gr->h[0], gr->h[1], gr->h[2]);
	set_tri_area_tol(hmin*hmin*1.0e-12);
	
	for(i=0; i<3; i++)
	{
	    kmin[i] = max(smin[i]-1, 0);
	    kmax[i] = min(smax[i]+1, gmax[i]);
	}

	total_nt = count_tris_in_top_box(kmin, kmax, intfc);

	if (debugging("gbd3d"))
	    printf("#total_nt = %d\n", total_nt);

	if(total_nt == 0)
	    return FUNCTION_SUCCEEDED;
	
	max_n_new = total_nt*6;
	uni_array(&test_tris,total_nt,sizeof(TRI*));
	uni_array(&ref_tris,total_nt,sizeof(TRI*));
	uni_array(&out_tris,total_nt*3,sizeof(TRI*));
	uni_array(&in_tris,total_nt*3,sizeof(TRI*));
	uni_array(&deg_tris,total_nt,sizeof(TRI*));
	uni_array(&new_tris, max_n_new, sizeof(TRI*));
	
	/*getting tris directly from the intfc. When there are two rect boxes
	  the new generated tris are not in the intfc table, getting tris 
	  from intfc table may skip these tris.
	*/
	num_test_tris = tris_set_in_top_box(test_tris, total_nt, 
					    kmin, kmax, intfc);
	
	if(debugging("box_intfc"))
	{
	    TRI	*sect_tris[500];
	    int	nstris;
	    static int	cnt = 0;
	    
	    sprintf(fname,"lgb_debug%d_%s_%d.plt", cnt,
	    	right_flush(pp_mynode(),PP_NODE_FIELD_WIDTH), box_index);
	    
	    cnt++;
	    printf("\n Enter grid_based_box_untangle, file name: %s\n", fname);
	    
	    if(pp_mynode() == 6)
		set_shift_for_tecplot(0.0, -1.0, 0.0);
	    else
		set_shift_for_tecplot(0.0, 0.0, 0.0);

	    /* open the tecplot debug file and output test_tris */
	    file = fopen(fname,"w");
	    fclose(file);
	   
	    /* output tris with intersection. */
	    nstris = tris_intersection(sect_tris, test_tris, num_test_tris);
	    
	    tecplot_debug_tris("sect_tris", fname, sect_tris, nstris, intfc);
	}

	/* finding ref_tris */
	num_ref_tris = 0;
	for(i=0; i<num_test_tris; i++)
	{
	    tri = test_tris[i];
	    if(skip_bdry_tri(tri) || !tri_in_box(tri,intfc,box) || 
	       (tri_tag && skip_tag_tri(tri)))
		continue;
	    if(!tri_recorded(tri,ref_tris,num_ref_tris))
		ref_tris[num_ref_tris++] = tri;
	}
	tecplot_debug_tris("ref_tris", fname, ref_tris, num_ref_tris, intfc);
	
	/* finding out_tris */

	num_out_tris = bound_tris_set(out_tris, ref_tris, num_ref_tris);

	/* removing ref_tris */
	for(i=0; i<num_ref_tris; i++)
	    remove_tri_from_surface(ref_tris[i],ref_tris[i]->surf,NO);
	
	tecplot_debug_tris("out_tris", fname, out_tris, num_out_tris,intfc);
	
	if (debugging("box_intfc"))
	{
	    file = fopen(fname,"a");
	    tecplot_interface_in_box("removed_surf",file,kmin,kmax,intfc);
	    fclose(file);
	}

	/* finding in_tris */
	i = 0;
	for(surfs=intfc->surfaces; surfs && *surfs; ++surfs)
	{
	    last_tris[i] = last_tri(*surfs);
	    ref_surfs[i] = *surfs;
	    i++;
	}
	num_surfs = i;
	reconstruct_intfc3d_in_box_lgb(intfc,smin,smax,NO,NULL);
	/* Check if any surface has been deleted */
	for (i = 0; i < num_surfs; ++i)
	{
	    if (!surf_in_interface(ref_surfs[i],intfc))
		last_tris[i] = NULL;
	}
	i = 0;
	num_in_tris = 0;
	for(surfs=intfc->surfaces; surfs && *surfs; ++surfs)
	{
	    if (last_tris[i] == NULL) continue; /* surface already deleted */
	    if (at_end_of_tri_list(last_tris[i],*surfs))
		continue;
	    for (tri=last_tris[i]->next; !at_end_of_tri_list(tri,*surfs);
		    	tri = tri->next)
		in_tris[num_in_tris++] = tri;
	    i++;
	}
	tecplot_debug_tris("in_tris", fname, in_tris, num_in_tris,intfc);
	
	/*NOTE: The following part only has relations with in_tris, out_tris, 
	  ref_tris. It can be extended to a more general form.
	  sep double loops for out_tris
	*/
	if (debugging("box_intfc"))
	    printf("#count in out ref tris  %d %d %d \n", 
		 num_in_tris, num_out_tris, num_ref_tris);
	
	sep_common_point_from_loop(out_tris, num_out_tris, 
				   ref_tris, &num_ref_tris, intfc);
	if (debugging("box_intfc"))
	    tecplot_debug_tris("sep comm pt out_tris", fname, 
			   out_tris, num_out_tris,intfc);
	
	if(debugging("comm_pt_inout"))
	{
	    num_new_tris = 0;
	    num_new_tris = merge_tris_set(new_tris, num_new_tris, 
					  out_tris, num_out_tris);
	    num_new_tris = merge_tris_set(new_tris, num_new_tris, 
					  in_tris, num_in_tris);
	    printf("#check in and out common point.\n");
	    if(!check_valid_tris(new_tris, num_new_tris, intfc))
	    {
		printf("ERROR grid_based_box_untangle "
		       "in_tris and out_tris share points.\n");
		tecplot_tris("in_tris", in_tris, num_in_tris);
		tecplot_tris("out_tris", out_tris, num_out_tris);
		clean_up(ERROR);
	    }
	}

	/* remove degenerated loops */
	num_deg_tris = 0;
	num_in_tris = seal_all_loops_wo_constraint(deg_tris, &num_deg_tris, 
					in_tris, num_in_tris, 1, YES);
	tecplot_debug_tris("in_tris deg", fname, in_tris, num_in_tris,intfc);
	
	num_out_tris = seal_all_loops_wo_constraint(deg_tris, &num_deg_tris, 
					out_tris, num_out_tris, 1, YES);
	tecplot_debug_tris("out_tris deg", fname, 
			     out_tris, num_out_tris, intfc);
	
	tecplot_debug_tris("deg_tris", fname, deg_tris, num_deg_tris, intfc);

	/* linking suitable pairs */
	if(debugging("pairsfix"))
	{
	    num_new_tris = linking_tris_with_pairs_fix(new_tris, max_n_new, 
		out_tris, num_out_tris, in_tris, num_in_tris, 
		ref_tris, num_ref_tris, intfc);

	}
	else
	{
	    num_new_tris = linking_tris_with_pairs(new_tris, max_n_new, 
		out_tris, num_out_tris, in_tris, num_in_tris, 
		ref_tris, num_ref_tris);
	}

	tecplot_debug_tris("new_tris pairs", fname, 
			   new_tris, num_new_tris, intfc);

	/*merge all current tris, do not merge deg_tris because 
	  they do not have null sides */
	num_new_tris = merge_tris_set(new_tris, num_new_tris, 
				      in_tris, num_in_tris);
	num_new_tris = merge_tris_set(new_tris, num_new_tris, 
				      out_tris, num_out_tris);
	tecplot_debug_tris("merge in out new tris", fname, 
			    new_tris, num_new_tris, intfc);

	/*removing all linking non-null edges in loops */
	num_new_tris = sep_common_edge_from_tris(&sep_new_tris, new_tris, 
						 num_new_tris, intfc);
	tecplot_debug_tris("sep comm edge merge", fname,
			    sep_new_tris, num_new_tris, intfc);

	/* sealing all the null loops */
	num_seal_tris = 0;
	num_new_tris = seal_all_loops_wo_constraint(new_tris, &num_seal_tris, 
					sep_new_tris, num_new_tris, 1, NO);
	
	tecplot_debug_tris("seal_all_loops sep_new_tris", fname,
			     sep_new_tris, num_new_tris, intfc);
	tecplot_debug_tris("seal_all_loops new_tris", fname,
			     new_tris, num_seal_tris, intfc);

	/* combining everything */
	num_new_tris = merge_tris_set(new_tris, num_seal_tris, 
				      sep_new_tris, num_new_tris);
	num_new_tris = merge_tris_set(new_tris, num_new_tris, 
				      deg_tris, num_deg_tris);

	/* smoothing in_tris because they are grid based */
	if(num_in_tris > 0)
	    smooth_tris(in_tris, num_in_tris);
	tecplot_debug_tris("smooth in_tris", fname,
			    new_tris, num_new_tris, intfc);

	if(tri_tag)
	{
	    for(i=0; i<num_new_tris; i++)
	      Tri_index(new_tris[i]) = 0;
	}

	if(debugging("box_intfc"))
	{
	    num_test_tris = tris_set_in_top_box(new_tris, max_n_new,
					    kmin, kmax, intfc);
	    file = fopen(fname,"a");
	    tecplot_show_tris("recon tris",new_tris,num_test_tris,file);
	    fclose(file);
	}

delete_s:
	/* no tris */
	for(surfs = intfc->surfaces; surfs && *surfs; surfs++)
	    if(no_tris_on_surface(*surfs) || (*surfs)->num_tri == 0)
	    {
	    	delete_surface(*surfs);
		goto delete_s;
	    }

	free_these(6,test_tris,ref_tris,in_tris,out_tris,new_tris,deg_tris);
	
	DEBUG_LEAVE(grid_based_box_untangle)
	return FUNCTION_SUCCEEDED;
}	/* end grid_based_box_untangle */


#define  LOOP_SMOOTH_PARA	100

LOCAL	void  smooth_null_loop(
	TRI	**null_tris,
	int	*null_sides,
	POINT	**pts,
	int	num_null_sides)
{
	int		i, j, num;
	TRI		*tri;
	POINT		*p;
	SURFACE		*s;
	SMOOTH_PARA	smooth_que[LOOP_SMOOTH_PARA];
	SMOOTH_TOL	stol;

	if(num_null_sides >= LOOP_SMOOTH_PARA)
	{
	    printf("ERROR smooth_null_loop, too many points in the null loop %d.\n", num_null_sides);
	    clean_up(ERROR);
	}

	/* smooth paramaters. */
	stol.cone_ratio = 0.1;
	stol.max_cos = 0.6;
	stol.alpha = sqrt(0.65);
	s = null_tris[0]->surf;

	/* Compute the the parameters in each points */
	num = 0;
	for(i=0; i<num_null_sides; i++)
	{
	    tri = null_tris[i];
	    p = Point_of_tri(tri)[null_sides[i]];

	    if(!compute_smooth_para(&smooth_que[num], p,tri,s,&stol))
		continue;
	    num++;
	}

	for(i=0; i<num; i++)
	    compute_point_smooth(&smooth_que[i], &stol, s->interface);
}

#define MAX_SMOOTH_POINT    5000
boolean	compute_average_point(SMOOTH_PARA*,POINT*,TRI*,SURFACE*,SMOOTH_TOL*);

	void smooth_tris(
	TRI	**tris,
	int	num_tris)
{
	int		i, j, k, num;
	TRI		*tri;
	POINT		*p;
	SURFACE		*s;
	SMOOTH_PARA	smooth_que[MAX_SMOOTH_POINT];
	SMOOTH_TOL	stol;

	if(num_tris == 0)
	    return;

	if(3*num_tris >= MAX_SMOOTH_POINT)
	{
	    printf("ERROR smooth_tris, too many points  %d.\n", num_tris);
	    clean_up(ERROR);
	}

	/* smooth paramaters. */
	stol.cone_ratio = 0.1;
	stol.max_cos = 0.6;
	stol.alpha = sqrt(0.65);

	for(k=0; k<num_tris; k++)
	{
	    for (i = 0; i < 3; i++)
	    {
		Index_of_point(Point_of_tri(tris[k])[i]) = -1;
	    }
	}

	s = tris[0]->surf;
	num = 0;
	
	/* Compute the the parameters in each point */
	for(k=0; k<num_tris; k++)
	{
	    tri = tris[k];
	    for (i = 0; i < 3; i++)
	    {
		p = Point_of_tri(tri)[i];
		if(Boundary_point(p) || Index_of_point(p) != -1)
		    continue;

		Index_of_point(p) = 1;
		if(!compute_average_point(&smooth_que[num], p,tri,s,&stol))
		    continue;
		
		num++;
	    }
	}

	/* Apply Laplacian smooth */
	for(i=0; i<num; i++)
	    compute_point_smooth(&smooth_que[i], &stol, s->interface);
}

	void  seal_null_loop(
	TRI	**null_tris,
	int	*null_sides,
	POINT	**pts,
	int	num_null_sides)
{
	int	i, j, side;
	double	avep[3];
	POINT	*midp, *pt, *ptn;
	TRI	*tri, *new_tri, *new_tri1, *prev_tri;
	SURFACE	*surf = null_tris[0]->surf;
	
	tmp_n = 0;

	/* find the center point */
	for(j=0; j<3; j++)
	    avep[j] = 0;
	for(i=0; i<num_null_sides; i++)
	{
	    for(j=0; j<3; j++)
		avep[j] += Coords(pts[i])[j];
	}
	for(j=0; j<3; j++)
	    avep[j] /= num_null_sides;
	
	tri = null_tris[0];
	side = null_sides[0];
	
	/*TMP copy the states from the first point, 
	  otherwise sl, sr will be NULL */
	midp = copy_point(Point_of_tri(tri)[side]);
	ft_assign(Coords(midp), avep, 3*FLOAT);
	
	printf("#seal_null_loop is called. num_null_sides = %d\n", num_null_sides);
	print_general_vector("avep= ", avep, 3, "\n");
	
	for(i=0; i<num_null_sides; i++)
	{
	    tri = null_tris[i];
	    side = null_sides[i];
	    pt = Point_of_tri(tri)[side];
	    ptn = Point_of_tri(tri)[Next_m3(side)];
	    
	    new_tri = make_tri(midp, ptn, pt, NULL, NULL, NULL, NO);
	    
	    tmp_tris[tmp_n++] = new_tri;
	    
	    /* link with the tri on the loop */
	    link_neighbor_tris(new_tri, tri);
	    
	    /* link with the new tri */
	    if(i == 0)
		new_tri1 = new_tri;
	    else  if(!link_neighbor_tris(new_tri, prev_tri))
	    {
		FILE    *fp;
		char    s[50];

		sprintf(s, "seal_null_%d", pp_mynode());
		printf("ERROR in seal_null_loop, two new tris do not share points, filename %s.\n", s);
	
		fp = fopen(s, "w");
		tecplot_show_null_tris(null_tris, pts, num_null_sides, fp);
		fclose(fp);
	
		clean_up(ERROR);
	    }

	    insert_tri_at_tail_of_list(new_tri,surf);
	    prev_tri = new_tri;
	}

	/* new_tri1 is the first new tri, new_tri is the last new_tri. */
	if(!link_neighbor_tris(new_tri, new_tri1))
	{
	    printf("ERROR in seal_null_loop, two new tris do not share "
		   "points for the first tri.\n");
	    clean_up(ERROR);
	}
}


/*************************************************************************
*                                                                        *
*                         /\                                             *
*                        /  \ in_tri                                     *
*                       /    \                                           *
*                    vi---------------------------------|                *
*                     /| in_side                        |                *
*                    / |                                |                *
*                   /  |                                |                *
*         base_tri  \  |      strip to be triangulated  |                *
*                    \ |                                |                *
*                     \| out_side                       |                *
*                    vo---------------------------------|                *
*                       \    /                                           *
*                        \  / out_tri                                    *
*                         \/                                             *
*                                                                        *
*************************************************************************/


LOCAL   boolean seal_strip_with_tris_once(
        RECT_BOX *box,
	TRI **null_tris,
	int *null_sides,
	int num_null_sides,
	POINT **pts,
	double *tnor,
	char *base_name,
	boolean zero_base,
	INTERFACE *intfc,
	TRI **out_tris,
	int num_out_tris,
	TRI **in_tris,
	int num_in_tris)
{
	RECT_GRID gr = topological_grid(intfc);
	RECT_GRID comp_grid = Computational_grid(intfc);
        POINT *p0,*p1,*p2,*p3;
	TRI *new_tri,**newtris,**debug_tris,*new_tri1,*new_tri2;
	TRI *in_tri,*out_tri,*base_tri,*test_tri1=NULL,*test_tri2=NULL,
	    *test_tri;
	int in_side,out_side;
	int i,j,num_newtris,new_num_null_sides,tmp_num,m;
	int i_base,i_in,i_out,count = 0,tmp,min_num,test_base;
	boolean end_in_stitch = NO,status_in = NO,done = NO,test = NO,
	     status_out = NO, crx_buffer2 = NO;
	double dist1,dist2,min_dist = HUGE,min_dist_buffer = HUGE;
	SURFACE *surf;
	FILE *file;
	char fname[100];
	double *h,*L,crx_coord[2],min_dist_buffer2 = HUGE;
	
	tmp_n = 0;

	L = gr.L;
	h = gr.h;
	new_num_null_sides = num_null_sides;
	for (i = 0; i < new_num_null_sides; ++i)
	{
	    if (!tri_recorded(null_tris[i],in_tris,num_in_tris))
	        status_in = YES;
	    if (!tri_recorded(null_tris[i],out_tris,num_out_tris))
	        status_out = YES;
	}
	while (!end_in_stitch)
	{
	    i_base = -1;
	    min_dist = HUGE;
	    min_dist_buffer = HUGE;
	    test = NO;
	    min_dist_buffer2 = HUGE;
	    crx_buffer2 = NO;
	   
	    /*try to determine from which triangle we begin to 
	      seal the null loop */
	    for (i = 0; i < new_num_null_sides; ++i)
	    {
		i_in  = (i+1)%new_num_null_sides;
		in_tri   = null_tris[i_in];
		if (null_tris[i] == in_tri) 
		    continue;
		if (new_num_null_sides <= 5)
		{
		    i_in  = (i+new_num_null_sides-2)%new_num_null_sides;
		    in_tri   = null_tris[i_in];
		    i_out  = (i+new_num_null_sides-1)%new_num_null_sides;
		    out_tri   = null_tris[i_out];

		    if (in_tri == out_tri)
		        continue;

		    i_in  = (i+2)%new_num_null_sides;
		    in_tri   = null_tris[i_in];
		    i_out  = (i+3)%new_num_null_sides;
		    out_tri   = null_tris[i_out];
		    if (in_tri == out_tri)
		        continue;
		}
		i_in  = (i+1)%new_num_null_sides;
		in_tri   = null_tris[i_in];
		p0 = Point_of_tri(null_tris[i])[null_sides[i]];
		p3 = Point_of_tri(null_tris[i])[Prev_m3(null_sides[i])];
		p1 = Point_of_tri(null_tris[i_in])[Next_m3(null_sides[i_in])];
		p2 = Point_of_tri(null_tris[i_in])[Prev_m3(null_sides[i_in])];
		
		/*never use 3 in_tri points to get a triangle, 
		  otherwise the tri will be on a face of a block. */
		if (status_in)
		    if (tri_recorded(null_tris[i],in_tris,num_in_tris) &&
		        tri_recorded(in_tri,in_tris,num_in_tris))
		        continue;
		if (two_points_share_side(p0, null_tris[i], p1,
		              (null_tris[i]->surf)->interface)==-1)
		    return NO;	/*set_tri_list (p0)  fails  (1) */
		if (two_points_share_side(p0, null_tris[i], p1,
		              (null_tris[i]->surf)->interface) == 1)
		    continue;
		
		if (two_points_share_side(p0, null_tris[i], p2,
		           (null_tris[i]->surf)->interface) == -1 &&
			   new_num_null_sides > 5)
		    return NO;  /*set_tri_list (p0)  fails, same as (1) */
		if (two_points_share_side(p0, null_tris[i], p2,
		           (null_tris[i]->surf)->interface) == 1 &&
			   new_num_null_sides > 5)
	            continue;   /*this is not the bad case  (4) */
		
		if (two_points_share_side(p3, null_tris[i], p1,
		           (null_tris[i]->surf)->interface) == -1 &&
			             new_num_null_sides > 5)
		    return NO;	/*set_tri_list (p3) fails */
		if (two_points_share_side(p3, null_tris[i], p1,
		           (null_tris[i]->surf)->interface) == 1 &&
			              new_num_null_sides > 5)
		    continue;   /*same as (4) */
		
		dist1 = distance_between_positions(Coords(p0),Coords(p1),3);
		/* seal from min dist two points and from the buffer zone */
		for (j = 0; j < 3; j++)
		{
		    crx_coord[0] = comp_grid.L[j];
		    crx_coord[1] = comp_grid.U[j];
                    /*this can never be satisfied, null_tris[i] 
                      is larger than the comp domain */
		    if ((tri_cross_zone(null_tris[i],crx_coord[0],h[j],j)==YES&&
		        tri_cross_zone(null_tris[i],crx_coord[1],h[j],j)==YES)&&
			(tri_cross_zone(in_tri,crx_coord[0],h[j],j) == YES &&
			tri_cross_zone(in_tri,crx_coord[1],h[j],j) == YES))
		    {
		        if (dist1 < min_dist_buffer2)
			{
			    test = YES;
			    crx_buffer2 = YES;
			    i_base = i;
			    min_dist_buffer2 = dist1;
			}
		    }
		    
		    /*always crx_buffer2 = NO; */
		    /*if a tri lies between two procs, we construct it first. */
		    if (!crx_buffer2)
		        for (m = 0; m < 2; m++)
		            if (tri_cross_zone(null_tris[i],
				    crx_coord[m],h[j],j) == YES &&
		                tri_cross_zone(in_tri,
				    crx_coord[m],h[j],j) == YES)
		        {
		            if (dist1 < min_dist_buffer)
		            {
		                test = YES;
		                i_base = i;
		                min_dist_buffer = dist1;
		            }
		        }
		}
		if (!test)
		    if (dist1 < min_dist)
		    {
		        i_base = i;
		        min_dist = dist1;
		    }
	    }

	    if (i_base == -1)
	    {
	        printf("WARNING in seal_strip_with_tris_once, "
		       "can not find the starting tri to seal.\n");
		return FUNCTION_FAILED;
	    }

	    /* form a triangle */
	    i_in  = (i_base+1)%new_num_null_sides;
	    in_tri   = null_tris[i_in];

	    i_out = (i_in+1)%new_num_null_sides;
	    out_tri = null_tris[i_out];

	    /* de part b */
	    if (out_tri == in_tri)
	    {
		test_tri = out_tri;
	    }
	    /* de part e */

	    i_in  = (i_base+1)%new_num_null_sides;
	    in_tri   = null_tris[i_in];
	    base_tri = null_tris[i_base];
	    surf = base_tri->surf;

	    p0 = Point_of_tri(null_tris[i_base])[null_sides[i_base]];
	    p1 = Point_of_tri(null_tris[i_in])[Next_m3(null_sides[i_in])];
	    p2 = Point_of_tri(null_tris[i_base])[Next_m3(null_sides[i_base])];

	    new_tri = make_tri(p0,p1,p2,NULL,(POINTER)in_tri,
	                    (POINTER)base_tri,NO);
	    
	    tmp_tris[tmp_n++] = new_tri;

	    
	    insert_tri_at_tail_of_list(new_tri,surf);
	    link_neighbor_tris(new_tri,in_tri);
	    link_neighbor_tris(new_tri,base_tri);

	    /*from another triangle if the remaining is a null triangle. */
	    i_in = (i_in+1)%new_num_null_sides;
	    i_out = (i_base+new_num_null_sides-1)%new_num_null_sides;
	    in_tri   = null_tris[i_in];
	    out_tri = null_tris[i_out];

	    /*we will never reach this part except the input null loop 
	      is a triangle. */
	    if (out_tri == in_tri)
	    {
	        end_in_stitch = YES;
		break;
	    }

	    if (Point_of_tri(in_tri)[Next_m3(null_sides[i_in])] ==
	        Point_of_tri(out_tri)[null_sides[i_out]])
	    {
		base_tri = new_tri;
		p0 = Point_of_tri(null_tris[i_in])[null_sides[i_in]];
		p1 = Point_of_tri(null_tris[i_out])[Next_m3(null_sides[i_out])];
		p2 = Point_of_tri(null_tris[i_in])[Next_m3(null_sides[i_in])];

		new_tri = make_tri(p0,p1,p2,(POINTER)base_tri,
		                (POINTER)out_tri,(POINTER)in_tri,NO);
	        
	        tmp_tris[tmp_n++] = new_tri;
		
		insert_tri_at_tail_of_list(new_tri,surf);
		link_neighbor_tris(new_tri,in_tri);
		link_neighbor_tris(new_tri,out_tri);
		link_neighbor_tris(new_tri,base_tri);
		end_in_stitch = YES;
		break;
	    }

	    new_num_null_sides--;
	    null_tris[i_base%new_num_null_sides] = new_tri;
	    null_sides[i_base%new_num_null_sides] = 0;
	    for (i = i_base+1; i < new_num_null_sides; i++)
	    {
		null_tris[i] = null_tris[i+1];
		null_sides[i] = null_sides[i+1];
	    }
	}

	return FUNCTION_SUCCEEDED;
	
}  /*end seal_strip_with_tris_once*/


LOCAL	boolean seal_strip_with_tris(
	RECT_BOX *box,
	TRI **null_tris,
	int *null_sides,
	int num_null_sides,
	POINT **pts,
	double *tnor,
	char *base_name,
	boolean zero_base,
	INTERFACE *intfc)
{
	POINT *p0,*p1,*p2,*vi,*vo,*p3,*p4;
	TRI *new_tri,**newtris,**debug_tris;
	TRI *in_tri,*out_tri,*base_tri;
	int in_side,out_side;
	int i,j,num_newtris;
	int i_base,i_in,i_out,count = 0;
	boolean end_in_stitch = NO,status = NO;
	double dist1,dist2,min_dist = HUGE;
	SURFACE *surf;
	FILE *file;
	char fname[100];
	

	if (debugging("box_intfc"))
	{
	    sprintf(fname,"%s-nullsides.list",base_name);
	    file = fopen(fname,"w");
	    gview_show_box_tri(box,null_tris,num_null_sides,file);
	    fclose(file);
	}
	    i_base = -1;
	    if (zero_base) i_base = 0;
	    else
	    {
		/* Determine base_tri:*/
		        for (i = 0; i < num_null_sides; ++i)
			{
			    i_out = (i+num_null_sides-1)%num_null_sides;
			    i_in  = (i+1)%num_null_sides;
			    p2 = Point_of_tri(null_tris[i_in])[Next_m3(
					    null_sides[i_in])];
			    p3 = Point_of_tri(null_tris[i_out])[
				    	    null_sides[i_out]];
			    p1 = Point_of_tri(null_tris[i])[null_sides[i]];
			    p2 = Point_of_tri(null_tris[i])[Next_m3(
					    null_sides[i])];
			    dist1 = distance_between_positions(
					    Coords(p1),Coords(p2),3);
		            if (dist1 < min_dist)
		            {
		                i_base = i;
		                min_dist = dist1;
			    }
			}
	    }
	    if (i_base == -1)
	    {
		printf("All the triangles share the same side\n");
		exit(0);
		return FUNCTION_SUCCEEDED;
	    }
	    i_out = (i_base+num_null_sides-1)%num_null_sides;
	    i_in  = (i_base+1)%num_null_sides;

	    in_tri   = null_tris[i_in];
	    base_tri = null_tris[i_base];
	    out_tri  = null_tris[i_out];
	    in_side  = null_sides[i_in];
	    out_side = null_sides[i_out];
	    
	    check_tri_and_neighbor(out_tri);
	    
      	    vi = Point_of_tri(in_tri)[null_sides[i_in]];
	    vo = Point_of_tri(out_tri)[Next_m3(null_sides[i_out])];
	    surf = base_tri->surf;
	
	if (debugging("box_intfc"))
	{
	    uni_array(&debug_tris,300,sizeof(TRI*));
	    num_newtris = 0;
	    for (i = 0; i < num_null_sides; ++i)
		debug_tris[num_newtris++] = null_tris[i];
	}
	while (!end_in_stitch)
	{
	    p2 = Point_of_tri(null_tris[i_in])[Next_m3(null_sides[i_in])];
	    p3 = Point_of_tri(null_tris[i_out])[null_sides[i_out]];
	    dist1 = distance_between_positions(Coords(vo),
			Coords(Point_of_tri(in_tri)[Next_m3(in_side)]),3);
	    dist2 = distance_between_positions(Coords(vi),
			Coords(Point_of_tri(out_tri)[out_side]),3);

	    if (((null_tris[i_in] == null_tris[(i_in+1)%num_null_sides])||
		(null_tris[i_base] == null_tris[i_out])) && count == 0)
	    {
		p0 = vi;
		p1 = vo;
		p2 = Point_of_tri(in_tri)[Next_m3(in_side)];
		new_tri = make_tri(p0,p1,p2,(POINTER)base_tri,	
		    		NULL,(POINTER)in_tri,NO);
	    	insert_tri_at_tail_of_list(new_tri,surf);
		link_neighbor_tris(new_tri,in_tri);
		link_neighbor_tris(new_tri,base_tri);
		if (debugging("box_intfc"))
		{    
		    debug_tris[num_newtris++] = new_tri;
		    sprintf(fname,"%s-debug-%d.list",base_name,
				    num_newtris-num_null_sides);
            	    file = fopen(fname,"w");
            	    gview_show_box_tri(box,debug_tris,num_newtris,file);
            	    fclose(file);
		}
		count++;
		vi = p2;
	    	i_in  = (i_in+1)%num_null_sides;
		in_tri = null_tris[i_in];
		in_side = null_sides[i_in];
		if (in_tri == out_tri)
		{
		    end_in_stitch = YES;
		    break;
		}
	    }
	    else if (((null_tris[i_out] == 
	    	null_tris[(i_out-1+num_null_sides)%num_null_sides])||
		    (null_tris[i_base] == null_tris[i_in])) && count == 0)
	    {
		p0 = vo;
		p1 = Point_of_tri(out_tri)[out_side];
		p2 = vi;
		new_tri = make_tri(p0,p1,p2,(POINTER)out_tri,
		    		NULL,(POINTER)base_tri,NO);
	    	insert_tri_at_tail_of_list(new_tri,surf);
		link_neighbor_tris(new_tri,out_tri);
		link_neighbor_tris(new_tri,base_tri);

 		if (debugging("box_intfc"))
                {
                    debug_tris[num_newtris++] = new_tri;
                    sprintf(fname,"%s-debug-%d.list",base_name,
				    num_newtris-num_null_sides);
                    file = fopen(fname,"w");
                    gview_show_box_tri(box,debug_tris,num_newtris,file);
                    fclose(file);
                }

		count++;
		vo = p1;
	    	i_out = (i_out+num_null_sides-1)%num_null_sides;
		out_tri = null_tris[i_out];
		out_side = null_sides[i_out];
		if (out_tri == in_tri)
		{
		    end_in_stitch = YES;
		    break;
		}
	    }
	    else if (dist1 < dist2)
	    {
		p0 = vi;
		p1 = vo;
		p2 = Point_of_tri(in_tri)[Next_m3(in_side)];
		new_tri = make_tri(p0,p1,p2,(POINTER)base_tri,	
		    		NULL,(POINTER)in_tri,NO);
	    	insert_tri_at_tail_of_list(new_tri,surf);
		link_neighbor_tris(new_tri,in_tri);
		link_neighbor_tris(new_tri,base_tri);

	
 		if (debugging("box_intfc"))
                {
                    debug_tris[num_newtris++] = new_tri;
                    sprintf(fname,"%s-debug-%d.list",base_name, 
				    num_newtris-num_null_sides);
                    file = fopen(fname,"w");
                    gview_show_box_tri(box,debug_tris,num_newtris,file);
                    fclose(file);
                }
	
		vi = p2;
	    	i_in  = (i_in+1)%num_null_sides;
		in_tri = null_tris[i_in];
		in_side = null_sides[i_in];
		if (in_tri == out_tri)
		{
		    end_in_stitch = YES;
		    break;
		}
	    }
	    else
	    {
		p0 = vo;
		p1 = Point_of_tri(out_tri)[out_side];
		p2 = vi;
		new_tri = make_tri(p0,p1,p2,(POINTER)out_tri,
		    		NULL,(POINTER)base_tri,NO);
	    	insert_tri_at_tail_of_list(new_tri,surf);
		link_neighbor_tris(new_tri,out_tri);
		link_neighbor_tris(new_tri,base_tri);

 		if (debugging("box_intfc"))
                {
                    debug_tris[num_newtris++] = new_tri;
                    sprintf(fname,"%s-debug-%d.list",base_name, 
				    num_newtris-num_null_sides);
                    file = fopen(fname,"w");
                    gview_show_box_tri(box,debug_tris,num_newtris,file);
                    fclose(file);
                }


		vo = p1;
	    	i_out = (i_out+num_null_sides-1)%num_null_sides;
		out_tri = null_tris[i_out];
		out_side = null_sides[i_out];
		if (out_tri == in_tri)
		{
		    end_in_stitch = YES;
		    break;
		}
	    }
	    base_tri = new_tri;
	    if (Point_of_tri(in_tri)[Next_m3(in_side)] ==
		Point_of_tri(out_tri)[out_side])
	    {
		p0 = vi;
		p1 = vo;
		p2 = Point_of_tri(out_tri)[out_side];
		base_tri = new_tri;
		new_tri = make_tri(p0,p1,p2,(POINTER)base_tri,
		    		(POINTER)out_tri,(POINTER)in_tri,NO);
	    	insert_tri_at_tail_of_list(new_tri,surf);
		link_neighbor_tris(new_tri,in_tri);
		link_neighbor_tris(new_tri,out_tri);
		link_neighbor_tris(new_tri,base_tri);
		
 		if (debugging("box_intfc"))
                {
                    debug_tris[num_newtris++] = new_tri;
                    sprintf(fname,"%s-debug-%d.list",base_name,
				    num_newtris-num_null_sides);
                    file = fopen(fname,"w");
                    gview_show_box_tri(box,debug_tris,num_newtris,file);
                    fclose(file);
                }
		end_in_stitch = YES;
	    }
	}
	if (debugging("box_intfc"))
	{
	    sprintf(fname,"%s-newtris.list",base_name);
	    file = fopen(fname,"w");
	    gview_show_box_tri(box,debug_tris,num_newtris,file);
	    fclose(file);
	    free_these(1,debug_tris);
	}
	return FUNCTION_SUCCEEDED;
}	/* end seal_strip_with_tris */

/*      
	Note: Function find_nearest_tri_pair() is gonna replace by 
             find_nearest_tri_pair_new()
*/ 
LOCAL	void find_nearest_tri_pair(
	TRI **in_tris,
	int num_in_tris,
	TRI **out_tris,
	int num_out_tris,
	TRI **in_tri,
	int *in_side,
	TRI **out_tri,
	int *out_side,
	RECT_BOX *box)
{
	double *pi,*po;
	int i,j,k,l;
	double dist;
	double min_dist = HUGE;
	RECT_GRID *rgr = box->grid;
	int ipo[3],*ipi;

	*in_tri = *out_tri = NULL;

	for (i = 0; i < num_out_tris; ++i)
	{
	    for (j = 0; j < 3; ++j)
	    {
	    	if (Tri_on_side(out_tris[i],j) == NULL)
		{
		    po = Coords(Point_of_tri(out_tris[i])[j]);
		    rect_in_which(po,ipo,rgr);
		    for (k = 0; k < num_in_tris; ++k)
		    {
		    	ipi = Tri_icoords(in_tris[k]);
			if (!adjacent_cell(ipo,ipi))
			    continue;
		    	for (l = 0; l < 3; ++l)
			{
			    if (Tri_on_side(in_tris[k],l) != NULL)
			    	continue;
			    pi = Coords(Point_of_tri(in_tris[k])[Next_m3(l)]);
		    	    dist = distance_between_positions(pi,po,3);
		    	    if (dist < min_dist)
		    	    {
				min_dist = dist;
				*out_tri = out_tris[i];
				*out_side = j;
				*in_tri = in_tris[k];
				*in_side = l;
		    	    }
			}
		    }
		}
	    }
	}
}	/* end find_nearest_tri_pair */


boolean the_tri1(TRI *);

LOCAL	void find_nearest_tri_pair_new(
	TRI **in_tris,
	int *num_in_tris,
	TRI **out_tris,
	int *num_out_tris,
	TRI **crx_tris,
	int *num_crx_tris,
	TRI **in_tri,
	int *in_side,
	TRI **out_tri,
	int *out_side,
        TRI **crx_tri,
	RECT_BOX *box,
	INTERFACE *intfc)
{
	RECT_GRID gr = topological_grid(intfc);
	double *pi,*po, *pi1, *po1,*tp;
	int i,j,k,m,l,test1 = -1,test2,n,test;
	double dist=0.0,dist1=0.0,dist2=0.0;
	double min_dist = HUGE,min_dist_buffer = HUGE;
	double min_dist_buffer2 = HUGE;
	RECT_GRID *rgr = box->grid;
	int ipo[3],tipo[3],ipo1[3],*ipi,*tipi;
	double tol = 1.0e-5,*L,*h,crx_coord[2];
	boolean status = NO,crx_buffer2 = NO;

	*in_tri = *out_tri = NULL;
	L = gr.L;
	h = gr.h;

	test = -1;
	for (i = 0; i < *num_out_tris; ++i)
	{
	    for (j = 0; j < 3; ++j)
	    {
	    	if (Tri_on_side(out_tris[i],j) != NULL)
		    continue;
		for (m = 0; m < *num_crx_tris; ++m)
		    if (two_tris_share_side(crx_tris[m],out_tris[i],j))
		    {
		        if(debugging("new_tris"))
			{
			    printf("crx found\n");
			    print_tri(crx_tris[m], intfc);
			}

		        for (k = 0; k < *num_in_tris; ++k)
		        {
			    for (l = 0; l < 3; ++l)
			    {
                                if (Tri_on_side(in_tris[k],l) != NULL)
				    continue;
				pi = Coords(Point_of_tri(in_tris[k])[l]);
				pi1 = Coords(Point_of_tri(
					in_tris[k])[Next_m3(l)]);
			        po = Coords(Point_of_tri(out_tris[i])[j]);
				po1 = Coords(Point_of_tri(
					out_tris[i])[Next_m3(j)]);
				
				if(debugging("new_tris") && the_tri1(in_tris[k]))
				{
				    printf("#in_tri check\n");
				    print_tri(in_tris[k], intfc);

				}

				if (!point_in_crx_tri(pi,crx_tris[m]) &&
				   !point_in_crx_tri(pi1,crx_tris[m]))
				    continue;

				if(debugging("new_tris"))
				{
				    printf("#in_tri found\n");
				    print_tri(in_tris[k], intfc);
				}

				if (point_in_crx_tri(pi,crx_tris[m]))
				{
				    dist = distance_between_positions(pi,po,3) +
				          distance_between_positions(pi,po1,3) +
					  distance_between_positions(pi1,po1,3);
				}
				else
				{
				    dist=distance_between_positions(pi1,po,3) +
				         distance_between_positions(pi1,po1,3) +
					 distance_between_positions(pi,po1,3);
				}
				
				/*Find suitable null sides pair*/
				if (dist < min_dist) 
				{
				    min_dist = dist;
				    *out_tri = out_tris[i];
				    *out_side = j;
				    *in_tri = in_tris[k];
				    *in_side = l;
                                    *crx_tri = crx_tris[m];
				}
			    }
			}
		    }
	    }

	    if(debugging("new_tris") && *out_tri != NULL)
	    {
	        printf("min_dist %24.16e \n", min_dist);
		print_tri(*out_tri, intfc);
		if(the_tri(*out_tri))
		{
	            printf("#out_tri choice\n");
	        }
	    }
	}
	return;
}	/* end find_nearest_tri_pair_new */

LOCAL   void find_nearest_tri_pair_crx(
        TRI **in_tris,
        int *num_in_tris,
        TRI **out_tris,
        int *num_out_tris,
        TRI **tri_list,
        int *num_tris,
        TRI **in_tri,
        int *in_side,                                                                            
        TRI **out_tri,
        int *out_side,
        TRI **crx_tri,
        RECT_BOX *box,
        INTERFACE *intfc)
{       
	RECT_GRID gr = topological_grid(intfc);
        double *pi,*po, *pi1, *po1,*tp;
        int i,j,k,m,l,test1 = -1,test2,n,test;
        double dist=0.0,dist1=0.0,dist2=0.0;
        double min_dist = HUGE,min_dist_buffer = HUGE;
        double min_dist_buffer2 = HUGE;
        RECT_GRID *rgr = box->grid; 
        int ipo[3],tipo[3],ipo1[3],*ipi,*tipi;
        double tol = 1.0e-5,*L,*h,crx_coord[2];
        boolean status = NO,crx_buffer2 = NO;
    
        *in_tri = *out_tri = NULL;
        L = gr.L;
        h = gr.h;

     
        test = -1;
        for (i = 0; i < *num_out_tris; ++i)
        {
            for (j = 0; j < 3; ++j)
            {
                if (Tri_on_side(out_tris[i],j) != NULL)
                    continue;
                for (m = 0; m < *num_tris; ++m)
                    if (two_tris_share_pts(tri_list[m],out_tris[i],j))
                    {
                        for (k = 0; k < *num_in_tris; ++k)
                        {
                            for (l = 0; l < 3; ++l)
                            {
                                if (Tri_on_side(in_tris[k],l) != NULL)
                                    continue;
                                pi = Coords(Point_of_tri(in_tris[k])[l]);
                                pi1 = Coords(Point_of_tri(
					in_tris[k])[Next_m3(l)]);
                                po = Coords(Point_of_tri(out_tris[i])[j]);
                                po1 = Coords(Point_of_tri(
					out_tris[i])[Next_m3(j)]);
                                if (!point_in_crx_tri(pi,tri_list[m]) &&
                                   !point_in_crx_tri(pi1,tri_list[m]))
                                    continue;
                                if (point_in_crx_tri(pi,tri_list[m]))
                                {
                                    dist=distance_between_positions(pi,po,3) +
                                         distance_between_positions(pi,po1,3) +
                                         distance_between_positions(pi1,po1,3);
                                }
                                else
                                {
                                    dist=distance_between_positions(pi1,po,3) +
                                         distance_between_positions(pi1,po1,3) +
                                         distance_between_positions(pi,po1,3);
                                }

                                /*Find suitable null sides pair*/
                                if (dist < min_dist)
                                {
                                    min_dist = dist;
                                    *out_tri = out_tris[i];
                                    *out_side = j;
                                    *in_tri = in_tris[k];
                                    *in_side = l;
                                    *crx_tri = tri_list[m];
                                }
                            }
                        }
                    }
            }
        }
	return;
}       /* end find_nearest_tri_pair_new */


LOCAL	void remove_crx_tri_on_edge(
	TRI *crx_tri,
	TRI *in_tri,
	int in_side,
	TRI **in_tris,
	int *num_in_tris)
{
	double	*p;
	int	i,j,k;

	for(i=0; i<2; i++)
	{
	    p = Coords(Point_of_tri(in_tri)[(in_side + i)%3]);
	    if(!point_in_crx_tri(p, crx_tri))
	        continue;
	    
	    /*find all the tris in in_tris s.t. point p is in the tris.
	      p lies in the edge of in_tris[j] and crx_tri */
	    for(j=0; j<*num_in_tris; j++)
	    {
	        if(!point_in_crx_tri(p, in_tris[j]))
		    continue;

		printf("#remove_crx_tri_on_edge active.\n");
		for(k=j; k<(*num_in_tris)-1; k++)
		    in_tris[k] = in_tris[k+1];
		(*num_in_tris)--;
	    }
	}
}

LOCAL  boolean  check_adjecant_constrain0(
	POINT	*p0,
	POINT	*p1,
	POINT	*p2,
	POINT	*p3, 
	TRI	**new_tris, 
	int	num_new_tris)
{
	int	i, j;
	TRI	*tri;
	POINT	*pt, *ptn;
	boolean	status_in, status_out;

	if(NO && debugging("box_intfc"))
	{
	if(min_null_pair_angle(Coords(p0), Coords(p1), Coords(p2), Coords(p3)) > 0.87 ||
	   min_null_pair_angle(Coords(p1), Coords(p0), Coords(p3), Coords(p2)) > 0.87)
	return NO;
	}

	status_in = NO;
	status_out = NO;
	for(i=0; i<num_new_tris; i++)
	{
	    tri = new_tris[i];
	    for(j=0; j<3; j++)
	    {
		if(Tri_on_side(tri, j) != NULL)
		    continue;
		pt = Point_of_tri(tri)[j];
		ptn = Point_of_tri(tri)[(j+1)%3];

		if( ((p0 == pt) && (p2 == ptn)) || ((p0 == pt) && (p2 == ptn)) ||
		    ((p3 == pt) && (p1 == ptn)) || ((p3 == pt) && (p1 == ptn)) )
		    return YES;

		/*only the out_side or the in_side have common points with 
		  new_tris, this constraint can avoid tangled surface. */
		if(p0 == pt || p0 == ptn || p1 ==pt || p1 == ptn)
		    status_out = YES;
		if(p2 == pt || p2 == ptn || p3 ==pt || p3 == ptn)
		    status_in = YES;
	    }
	}

	if(status_in && status_out)
	    return NO;
	else
	    return YES;
}

LOCAL  boolean tri_in_box_pre(
       TRI *tri,
       INTERFACE *intfc,
       RECT_BOX *box)
{
       double *L,*h,min_coord,max_coord,p_coord;
       int i,j,k,m,count = 0;
       RECT_GRID gr = topological_grid(intfc);

       L = gr.L;
       h = gr.h;

       for (j = 0; j < 3; ++j)
       {
           for (i = 0; i < 3; ++i)
	   {
	       min_coord = L[i]+box->bmin[i]*h[i];
	       max_coord = L[i]+box->bmax[i]*h[i];
	       p_coord = Coords(Point_of_tri(tri)[j])[i];
	       if (p_coord >= min_coord && p_coord <= max_coord)
	           count++;   
	   }
	   if (count == 3)
	       return YES;
	   else
	       count = 0;
           
       }

       return NO;
}


LOCAL  boolean tri_in_box(
       TRI *tri,
       INTERFACE *intfc,
       RECT_BOX *box)
{
	RECT_GRID 	*gr = &topological_grid(intfc);
	
	return tri_in_grid_block(tri, box->bmin, box->bmax, gr);
}


/* remove all triangles connected with boundary tris */
LOCAL   boolean crx_bnd_out_tris(
        TRI **tri_list,
        int num_tris,
        TRI **test_tris,
        int num_test_tris,
        TRI *tri,
        INTERFACE *intfc,
        RECT_BOX *box)
{
        double *L,*h,min_coord,max_coord,p_coord;
        int i,j,k,m,count = 0;
        boolean status; 
	RECT_GRID gr = topological_grid(intfc);
            
        L = gr.L;
        h = gr.h;

        status = NO;

        for (j = 0; j < 3; ++j)
        {
            for (i = 0; i < num_test_tris; ++i)
            {
                for (k = 0; k < 3; ++k)
                    if (Point_of_tri(tri)[j] == Point_of_tri(test_tris[i])[k])
                        break;
                if (k < 3)
                {
                    status = YES;
                    break;
                }
            }
        }

        if (status)
            return YES;
        else
            return NO;

} /* end  crx_bnd_out_tris */

/*make sure crx_tri has only one common point with out_tri, 
  so perform find_nearest_tri_pair_crx in the following is valid.
*/
LOCAL   boolean crx_out_tris_twice(
        TRI **tri_list,		/*out_tris */
	int num_tris,		/*test_out_tris tris with bdry curves */
	TRI **test_tris,
	int num_test_tris,
	TRI *tri,
	INTERFACE *intfc,
	RECT_BOX *box)
{
	RECT_GRID gr = topological_grid(intfc);
        double *L,*h,min_coord,max_coord,p_coord;
	int i,j,k,m,count = 0;
	boolean status;

	L = gr.L;
	h = gr.h;

	status = NO;

	for (j = 0; j < 3; ++j)
	{
	    for (i = 0; i < num_tris; ++i)
	    {
	        for (k = 0; k < 3; ++k)
		    if (Point_of_tri(tri)[j] == Point_of_tri(tri_list[i])[k])
		        break;
	        if (k < 3)    /*k==3 means is the out_tri, can not remove it */
		{    
		    ++count;
		    break;
		}
	    }
	    for (i = 0; i < num_test_tris; ++i)
	    {
	        for (k = 0; k < 3; ++k)
		    if (Point_of_tri(tri)[j] == Point_of_tri(test_tris[i])[k])
		        break;
	        if (k < 3)    /*k==3 means boundary tri, can not remove it */
		{
		    status = YES;
		    break;
		}
	    }
	}

	/* count: number of points of tri lie in the out_tri */
	if (count >= 2 || status)
	    return YES;
	else
	    return NO;

} /* end  crx_out_tris_twice */

LOCAL	void gview_show_box_tri(
	RECT_BOX *box,
	TRI **tris,
	int num_tris,
	FILE *file)
{
	static const char *indent = "    ";
	POINT *p;
	int i,j,k;
	int nls;	/* number of grid lines */
	double *L = box->grid->L;
	double *U = box->grid->U;
	double *h = box->grid->h;

	(void) fprintf(file,"{ LIST\n");

	nls = (box->bmax[0] - box->bmin[0] + 1)*
	      (box->bmax[1] - box->bmin[1] + 1) +
	      (box->bmax[1] - box->bmin[1] + 1)*
	      (box->bmax[2] - box->bmin[2] + 1) +
	      (box->bmax[2] - box->bmin[2] + 1)*
	      (box->bmax[0] - box->bmin[0] + 1);

	(void) fprintf(file,"%s{\n%s%sOFF\n%s%s%6d %6d %6d\n",
	    		indent,indent,indent,indent,indent,
			2*nls,nls,0);
	for (i = box->bmin[0]; i <= box->bmax[0]; ++i)
	{
	    for (j = box->bmin[1]; j <= box->bmax[1]; ++j)
	    {
	    	(void) fprintf(file, "%s%s%-9g %-9g %-9g\n",indent,indent,
			L[0] + i*h[0],L[1] + j*h[1],
			L[2] + box->bmin[2]*h[2]);
	    	(void) fprintf(file, "%s%s%-9g %-9g %-9g\n",indent,indent,
			L[0] + i*h[0],L[1] + j*h[1],
			L[2] + box->bmax[2]*h[2]);
	    }
	}
	for (j = box->bmin[1]; j <= box->bmax[1]; ++j)
	{
	    for (k = box->bmin[2]; k <= box->bmax[2]; ++k)
	    {
	    	(void) fprintf(file, "%s%s%-9g %-9g %-9g\n",indent,indent,
			L[0] + box->bmin[0]*h[0],L[1] + j*h[1],
			L[2] + k*h[2]);
	    	(void) fprintf(file, "%s%s%-9g %-9g %-9g\n",indent,indent,
			L[0] + box->bmax[0]*h[0],L[1] + j*h[1],
			L[2] + k*h[2]);
	    }
	}
	for (k = box->bmin[2]; k <= box->bmax[2]; ++k)
	{
	    for (i = box->bmin[0]; i <= box->bmax[0]; ++i)
	    {
	    	(void) fprintf(file, "%s%s%-9g %-9g %-9g\n",indent,indent,
			L[0] + i*h[0],L[1] + box->bmin[1]*h[1],
			L[2] + k*h[2]);
	    	(void) fprintf(file, "%s%s%-9g %-9g %-9g\n",indent,indent,
			L[0] + i*h[0],L[1] + box->bmax[1]*h[1],
			L[2] + k*h[2]);
	    }
	}
	for (i = 0; i < nls; ++i)
	{
	    (void) fprintf(file,"%s%s%-4d %-4d %-4d\n",indent,indent,
			2,2*i,2*i+1);
	}
	(void) fprintf(file,"%s}\n",indent);

	(void) fprintf(file,"%s{\n%s%sOFF\n%s%s%6d %6d %6d\n",
	    		indent,indent,indent,indent,indent,
			3*num_tris,num_tris,0);
	for (i = 0; i < num_tris; ++i)
	{
	    for (j = 0; j < 3; ++j)
	    {
	    	p = Point_of_tri(tris[i])[j];
	    	(void) fprintf(file, "%s%s%-9g %-9g %-9g\n",indent,indent,
		    	Coords(p)[0],Coords(p)[1],Coords(p)[2]);
	    }
	}
	for (i = 0; i < num_tris; ++i)
	{
	    (void) fprintf(file,"%s%s%-4d %-4d %-4d %-4d\n",indent,indent,
			3,3*i,3*i+1,3*i+2);
	}
	(void) fprintf(file,"%s}\n",indent);
	(void) fprintf(file,"}\n");
}	/* end gview_show_box_tri */


LOCAL	void tecplot_show_box(
	char	 *bname,
	RECT_BOX *box,
	FILE 	 *file)
{
	int 	i;
	double	lc[3], rc[3];
	double 	*L = box->grid->L;
	double 	*U = box->grid->U;
	double 	*h = box->grid->h;

	for(i=0; i<3; i++)
	{
	    lc[i] = L[i] + box->bmin[i]*h[i];
	    rc[i] = L[i] + box->bmax[i]*h[i];
	}

	fprintf(file, "ZONE T=\"%s\" I=16\n", bname);
	
	fprintf(file, "%-9g  %-9g  %-9g\n", lc[0], lc[1], lc[2]);
	fprintf(file, "%-9g  %-9g  %-9g\n", lc[0], rc[1], lc[2]);
	fprintf(file, "%-9g  %-9g  %-9g\n", rc[0], rc[1], lc[2]);
	fprintf(file, "%-9g  %-9g  %-9g\n", rc[0], lc[1], lc[2]);
	fprintf(file, "%-9g  %-9g  %-9g\n", lc[0], lc[1], lc[2]);
	
	fprintf(file, "%-9g  %-9g  %-9g\n", lc[0], lc[1], rc[2]);
	
	fprintf(file, "%-9g  %-9g  %-9g\n", lc[0], rc[1], rc[2]);
	fprintf(file, "%-9g  %-9g  %-9g\n", lc[0], rc[1], lc[2]);
	fprintf(file, "%-9g  %-9g  %-9g\n", lc[0], rc[1], rc[2]);
	
	fprintf(file, "%-9g  %-9g  %-9g\n", rc[0], rc[1], rc[2]);
	fprintf(file, "%-9g  %-9g  %-9g\n", rc[0], rc[1], lc[2]);
	fprintf(file, "%-9g  %-9g  %-9g\n", rc[0], rc[1], rc[2]);


	fprintf(file, "%-9g  %-9g  %-9g\n", rc[0], lc[1], rc[2]);
	fprintf(file, "%-9g  %-9g  %-9g\n", rc[0], lc[1], lc[2]);
	fprintf(file, "%-9g  %-9g  %-9g\n", rc[0], lc[1], rc[2]);
	
	fprintf(file, "%-9g  %-9g  %-9g\n", lc[0], lc[1], rc[2]);

}

/*
	void tecplot_show_box_tri(
	char		*tname,
	RECT_BOX	*box,
	TRI		**tris,
	int		num_tris,
	FILE 		*file)
{
	int 	i,j;
	double	lc[3], rc[3];
	double 	*L = box->grid->L;
	double 	*U = box->grid->U;
	double	*h = box->grid->h;

	for(i=0; i<3; i++)
	{
	    lc[i] = L[i] + box->bmin[i]*h[i];
	    rc[i] = L[i] + box->bmax[i]*h[i];
	}
	fprintf(file,"TITLE = \"tecplot tris\"\n"
		     "VARIABLES = \"x\", \"y\", \"z\"\n");

	tecplot_box(NULL, file, lc, rc);
	tecplot_show_tris(tname, tris, num_tris, file);
}
*/

LOCAL	void tecplot_show_null_tris(
	TRI	**tris,
	POINT	**pts,
	int	num_tris,
	FILE 	*file)
{
int	i;

	tecplot_show_tris("null_tris", tris, num_tris, file);
	
	fprintf(file, "ZONE T=\"pts\" I=%d\n", num_tris);
	for (i = 0; i < num_tris; ++i)
	{
	    fprintf(file, "%-9g %-9g %-9g\n",
		    	Coords(pts[i])[0],Coords(pts[i])[1],Coords(pts[i])[2]);
	}
}

LOCAL   boolean null_sides_sharing_same_vertex(
        TRI *tri,
	int side,
	TRI **tri_list,
	int num_tri)
{
       int i,l;

       for (i = 0; i < num_tri; i++)
           for (l = 0; l < 3; l++)
	   {
	       if (Tri_on_side(tri_list[i],l) == NULL && tri != tri_list[i])
	       {
	           if (Point_of_tri(tri)[side] == 
		       Point_of_tri(tri_list[i])[Next_m3(l)] ||
		       Point_of_tri(tri)[Next_m3(side)] ==
		       Point_of_tri(tri_list[i])[l])
		       return YES;
	       }
	   }
       return NO;

       
}  /*end null_sides_sharing_same_vertex*/

LOCAL	double box_dist(
	RECT_BOX *box)
{
	double	cen, sum;
	int	i;

	sum = 0.0;

	for(i=0; i<3; i++)
	{
	    cen = (box->bmin[i] + box->bmax[i])*0.5;
	    sum += min(box->smax[i]-cen, cen-box->smin[i]);
	}

	return sum;
}

LOCAL	void	make_boxes_from_ips(
	int		**ips_in,
	int		num_ip,
	RECT_BOX	**boxes)
{
	RECT_BOX 	Box,*box;
	int 		**pip,nb, itmp[3], i,j,k;
	static	int	**ips = NULL;

	DEBUG_ENTER(rm_bad_crxs_in_box)
	
	if(ips == NULL)
	    stat_matrix(&ips,MAX_NUM_UNPHY_IP,3,INT);

	for(i=0; i<num_ip; i++)
	    ft_assign(ips[i], ips_in[i], 3*INT);

	Box.prev = Box.next = NULL;
	box = &Box;
	
	i = 0;
	while (i < num_ip)
	{
	    box->next = (RECT_BOX *)store(sizeof(RECT_BOX));
	    box->next->prev = box;
	    box->next->next = NULL;
	    box = box->next;

	    /* move connected ips together, begin with pip, number is nb. */
	    pip = ips + i;
	    nb = 1;
	    for (j = i+1; j < num_ip; ++j)
	    {
		if (ip_connected(pip,nb,ips[j]))
		{
		    if (j > i+nb)
		    {
			for (k = 0; k < 3; ++k)
			{
			    itmp[k] = ips[j][k];
			    ips[j][k] = pip[nb][k];
			    pip[nb][k] = itmp[k];
			}
		    }
		    j = i + nb;
		    nb++;
		}
	    }

	    /* makeing the box from min and max. */
	    for (k = 0; k < 3; k++)
	    	box->bmin[k] = box->bmax[k] = pip[0][k];
	    for (j = 1; j < nb; j++)
	    {
	    	for (k = 0; k < 3; k++)
		{
		    if (box->bmin[k] > pip[j][k])
		    	box->bmin[k] = pip[j][k];
		    if (box->bmax[k] < pip[j][k])
		    	box->bmax[k] = pip[j][k];
		}
	    }
	    i += nb;
	}
	
	*boxes = Box.next;
}

LOCAL	void	copy_rect_boxes(
	RECT_BOX	**boxes,
	RECT_BOX	*boxes_in)
{
	RECT_BOX 	Box,*box,*box_in;
	
	Box.prev = Box.next = NULL;
	box = &Box;
	
	for(box_in = boxes_in; box_in != NULL; box_in = box_in->next)
	{
	    box->next = (RECT_BOX *)store(sizeof(RECT_BOX));
	    box->next->prev = box;
	    box->next->next = NULL;
	    box = box->next;
	    
	    ft_assign(box->bmin, box_in->bmin, 3*INT);
	    ft_assign(box->bmax, box_in->bmax, 3*INT);
	    ft_assign(box->smin, box_in->smin, 3*INT);
	    ft_assign(box->smax, box_in->smax, 3*INT);
	    box->grid = box_in->grid;
	}

	*boxes = Box.next;
}

LOCAL	boolean	merge_adjacent_boxes(
	RECT_BOX	*boxes)
{
	RECT_BOX 	*box, *nbox;
	boolean		boxes_merged;
	int		i;

	/* Merge overlapping boxes */
	boxes_merged = YES;
	while(boxes_merged)
	{
	    boxes_merged = NO;
	    for(box=boxes; box!=NULL; box=box->next)
	    {
		if(max3(box->bmax[2]-box->bmin[2],
			box->bmax[1]-box->bmin[1],
			box->bmax[0]-box->bmin[0]) > 10)
			{
			    printf("WARNING merge_adjacent_boxes "
			    	    "box is too large when merging.\n");
			    return NO;
			}

		for(nbox=box->next; nbox!=NULL; nbox=nbox->next)
		{
		    if(overlapping_boxes(box, nbox))
		    {
			for (i = 0; i < 3; ++i)
			{
			    box->bmin[i] = min(box->bmin[i],nbox->bmin[i]);
			    box->bmax[i] = max(box->bmax[i],nbox->bmax[i]);
			}
			nbox->prev->next = nbox->next;
			if (nbox->next != NULL)
			    nbox->next->prev = nbox->prev;
			boxes_merged = YES;
			break;
		    }
		}   /* for(nbox= ) */
		if(boxes_merged)
		    break;
	    }   /* for(box=  ) */
	}
	
	return YES;
}

LOCAL   boolean    use_rect_tris = YES;
boolean    set_use_rect_tris(boolean);
boolean    set_use_rect_tris(boolean fg)
{           
	use_rect_tris = fg;
}


/* orig version set_reconstruction_boxes_prev */
LOCAL	boolean  set_reconstruction_boxes(
	int *smin,
	int *smax,
	int **ips,
	int num_ip,
	RECT_BOX **boxes,
	INTERFACE *intfc)
{
	RECT_BOX 	*box,*nbox, *box0;
	boolean		boxes_merged;
	RECT_GRID	*gr = &topological_grid(intfc);
	int		i;

	DEBUG_ENTER(set_reconstruction_boxes)

	make_boxes_from_ips(ips, num_ip, &box0);

	/* check box size. */
	for(box=box0; box!=NULL; box=box->next)
	{
	    if(max3(box->bmax[2]-box->bmin[2],
	    	    box->bmax[1]-box->bmin[1], 
	    	    box->bmax[0]-box->bmin[0]) > 10)
	    {
		printf("WARNING set_reconstruction_boxes, too many bad ips");
                printf("box size = %d %d %d\n",
                            box->bmax[0]-box->bmin[0],
                            box->bmax[1]-box->bmin[1],
                            box->bmax[2]-box->bmin[2]);
	    }
	}

	if (debugging("tri_set") && use_rect_tris &&
	    rect_boxes_from_tangled_tris(box0, intfc) == 0)
	{
	    printf("#no tangled tris in box.\n");
	}

	for(box=box0; box!=NULL; box=box->next)
	{
	    for (i = 0; i < 3; i++)
	    {
	        box->smin[i] = smin[i];
		box->smax[i] = smax[i];
	    }

	    for (i = 0; i < 3; i++)
	    {
		if(debugging("lgbox"))
		{
		    box->bmin[i]--;
		    box->bmax[i]++;
		}
		
		if(box->bmax[i] - box->bmin[i] == 0)
		{
		    box->bmin[i]--;
		    box->bmax[i]++;
		}
		if (box->bmin[i] < smin[i]) box->bmin[i] = smin[i];
		if (box->bmax[i] > smax[i]) box->bmax[i] = smax[i];
	    }
	    box->grid = gr;
	}

	if(!merge_adjacent_boxes(box0))
	{
	    DEBUG_LEAVE(set_reconstruction_boxes)
	    return NO;
	}

	i = 0;
	for (box=box0; box!=NULL; box=box->next)
	{
	    if (debugging("lgbox"))
	    {
	    	print_int_vector("bminn", box->bmin, 3, "\n");
	    	print_int_vector("bmaxn", box->bmax, 3, "\n");
	    	printf("#blen %d    %3d %3d %3d\n", i,
	    	   	box->bmax[0] - box->bmin[0],
	    	   	box->bmax[1] - box->bmin[1],
	    	   	box->bmax[2] - box->bmin[2]);
	    }
            for(i=0; i<3; i++)
                if(box->bmax[i] - box->bmin[i] > 7)
                    box->bmax[i] = box->bmin[i] + 7;
	}
	*boxes = box0;
	
	DEBUG_LEAVE(set_reconstruction_boxes)
	return YES;
}	/* end set_reconstruction_boxes */

/* ref: rm_bad_crxs_in_box_prev */
LOCAL	boolean  rm_bad_crxs_in_box(
	int 		*smin,
	int 		*smax,
	int 		**ips,
	int		num_ip,
	RECT_BOX	*box0,
	INTERFACE	*intfc)
{
	RECT_BOX 	*box1, *box2, *box, *last_box;
	int 		i;
	RECT_GRID	*gr = &topological_grid(intfc);

	DEBUG_ENTER(rm_bad_crxs_in_box)
	
	make_boxes_from_ips(ips, num_ip, &box1);
	for(box = box1; box != NULL; box = box->next)
	{
	    for (i = 0; i < 3; i++)
	    {
	        box->smin[i] = smin[i];
		box->smax[i] = smax[i];
	    }
	
	    /*bad ips are in the faces of [bmin, bmax], 
	      should increase by one to make consistent surf
	      in the face of the box.
	    */
	    for (i = 0; i < 3; i++)
	    {
		box->bmin[i]--;
		box->bmax[i]++;
		if (box->bmin[i] < smin[i]) box->bmin[i] = smin[i];
		if (box->bmax[i] > smax[i]) box->bmax[i] = smax[i];
	    }
	    box->grid = gr;
	    if(box->next == NULL)
		last_box = box;
	}
	copy_rect_boxes(&box2, box0);
	
	/* since num_ip != 0 box1 != NULL, but box2 can be NULL. */
	
	last_box->next = box2;
	if(box2 != NULL)
	    box2->prev = last_box;

	merge_adjacent_boxes(box1);
	
	i = 0;
	for(box = box1; box != NULL; box = box->next)
	{
	    i++;

	    /* fix crossings in a box */
	    remove_unphysical_crxings(box->bmin,box->bmax,gr->gmax,
			intfc,SINGLE,&num_ip,ips);
	    check_and_repair_crx(intfc,box->bmin,box->bmax);
	}

	DEBUG_LEAVE(rm_bad_crxs_in_box)
	return YES;
}


LOCAL	boolean  set_tst_recon_boxes(
	int *smin,
	int *smax,
	RECT_BOX **boxes,
	INTERFACE *intfc)
{
	RECT_BOX 	Box,*box,*nbox;
	int 		**pip,nb,itmp[3],i,j,k;
	boolean		boxes_merged, status;
	RECT_GRID	*gr = &topological_grid(intfc);
	int		tmin[3], tmax[3];

	DEBUG_ENTER(set_reconstruction_boxes)
	
	Box.prev = Box.next = NULL;
	box = &Box;

	printf("#random boxes\n");
	for(i=0; i<20; i++)
	{
	    box->next = (RECT_BOX *)store(sizeof(RECT_BOX));
	    box->next->prev = box;
	    box->next->next = NULL;
	    box = box->next;
	    
	    for (k = 0; k < 3; k++)
	    	box->bmin[k] = box->bmax[k] = smin[k] + 
			    (int)((smax[k] - smin[k])*drand48());
	    printf("#randbx  %d %d %d\n ", 
	    	box->bmin[0], box->bmin[1], box->bmin[2]);
	    fflush(NULL);
	}
	printf("\n");
	
	/*bad ips are in the faces of [bmin, bmax], 
	  should increase by one to make consistent surf
	  in the face of the box.
	*/
	for (box = Box.next; box != NULL; box = box->next)
	{
	    for (i = 0; i < 3; i++)
	    {
		box->smin[i] = smin[i];
		box->smax[i] = smax[i];
	
		if(i == 0)
		{
		    box->bmin[i] -= 2;
		    box->bmax[i] += 1;
		}
		else if(i == 1)
		{
		    box->bmin[i] -= 2;
		    box->bmax[i] += 1;
		}
		else
		{
		    box->bmin[i] -= 2;
		    box->bmax[i] += 2;
		}

		if (box->bmin[i] < smin[i]) box->bmin[i] = smin[i];
		if (box->bmax[i] > smax[i]) box->bmax[i] = smax[i];
		
		/* avoid boundary surfaces, assume box size is at least 2. */
		if(!buffered_boundary_type(rect_boundary_type(intfc,i,0)))
		    if (box->bmin[i] == smin[i])
			box->bmin[i] = smin[i]+1;
		if(!buffered_boundary_type(rect_boundary_type(intfc,i,1)))
		    if (box->bmax[i] == smax[i])
			box->bmax[i] = smax[i]-1;
	    }
	    box->grid = gr;
	}
	
	/* Merge overlapping boxes */
	boxes_merged = YES;
	while (boxes_merged)
	{
	    boxes_merged = NO;
	    for (box = Box.next; box != NULL; box = box->next)
	    {
		if(max3(box->bmax[2]-box->bmin[2],
			box->bmax[1]-box->bmin[1],
			box->bmax[0]-box->bmin[0]) > 9)
			{
			    printf("WARN set_tst_recon_boxes "
			    	    "box is too large.\n");
			    DEBUG_LEAVE(set_tst_recon_boxes)
			    return NO;
			}

		for (nbox = box->next; nbox != NULL; nbox = nbox->next)
		{
		    if (overlapping_boxes(box, nbox) || 
			boxes_sharing_tris(box, nbox, smax, intfc)) 
		    {
			/* merge nbox with box and delete nbox. */
			for (i = 0; i < 3; ++i)
			{
			    box->bmin[i] = min(box->bmin[i],nbox->bmin[i]);
			    box->bmax[i] = max(box->bmax[i],nbox->bmax[i]);
			}
			nbox->prev->next = nbox->next;
			if (nbox->next != NULL)
			    nbox->next->prev = nbox->prev;
			boxes_merged = YES;
			break;
		    }
		}   /* for(nbox= ) */
		if(boxes_merged)
		    break;
	    }   /* for(box=  ) */
	}
	*boxes = Box.next;

	i = 0;
	for (box = Box.next; box != NULL; box = box->next)
	{
	    i++;
	}

	DEBUG_LEAVE(set_tst_recon_boxes)
	return YES;

}	/* end set_tst_recon_boxes */



LOCAL	boolean overlapping_boxes(
	RECT_BOX *box1,
	RECT_BOX *box2)
{
	int i,j,k,ip[3];
	for (i = 0; i < 2; ++i)
	{
	    ip[0] = (i == 0) ? box1->bmin[0] : box1->bmax[0];
	    if (ip[0] > box2->bmax[0] || ip[0] < box2->bmin[0])
	    	continue;
	    for (j = 0; j < 2; ++j)
	    {
	    	ip[1] = (j == 0) ? box1->bmin[1] : box1->bmax[1];
	    	if (ip[1] > box2->bmax[1] || ip[1] < box2->bmin[1])
	    	    continue;
	    	for (k = 0; k < 2; ++k)
		{
	    	    ip[2] = (k == 0) ? box1->bmin[2] : box1->bmax[2];
	    	    if (ip[2] > box2->bmax[2] || ip[2] < box2->bmin[2])
	    	    	continue;
		    return YES;
		}
	    }
	}
	for (i = 0; i < 2; ++i)
	{
	    ip[0] = (i == 0) ? box2->bmin[0] : box2->bmax[0];
	    if (ip[0] > box1->bmax[0] || ip[0] < box1->bmin[0])
	    	continue;
	    for (j = 0; j < 2; ++j)
	    {
	    	ip[1] = (j == 0) ? box2->bmin[1] : box2->bmax[1];
	    	if (ip[1] > box1->bmax[1] || ip[1] < box1->bmin[1])
	    	    continue;
	    	for (k = 0; k < 2; ++k)
		{
	    	    ip[2] = (k == 0) ? box2->bmin[2] : box2->bmax[2];
	    	    if (ip[2] > box1->bmax[2] || ip[2] < box1->bmin[2])
	    	    	continue;
		    return YES;
		}
	    }
	}
	return NO;
}	/* end end overlapping_box */

#define	MAX_BOUND_TRIS	10000
LOCAL   boolean boxes_sharing_tris(
        RECT_BOX  *box1,
	RECT_BOX  *box2,
	int       *smax,
	INTERFACE *intfc)
{
        TRI      **tris;
	TRI	 *tri_list1[MAX_BOUND_TRIS], *tri_list2[MAX_BOUND_TRIS];
	int      i, j, k, l, m,  nt;
	int      num_tris1, num_tris2;
	struct   Table *T;
	boolean     flag = NO;
	int      upbound[3];

	num_tris1 = 0;
	num_tris2 = 0;
	T = table_of_interface(intfc);

	for (i = 0; i < 3; ++i)
	    upbound[i] = (box1->bmax[i] == smax[i])? box1->bmax[i]-1 :
	                        box1->bmax[i];

	/* Forming the tri_list on the boundary of box1 */
	for (m = 0; m < 2; ++m)
	{
	    k = (m == 0) ? box1->bmin[2] : upbound[2];
	    for (j = box1->bmin[1]; j <= upbound[1]; ++j)
	    {
	        for (i = box1->bmin[0]; i <= upbound[0]; ++i)
		{
		    tris = T->tris[k][j][i];
		    nt = T->num_of_tris[k][j][i];
		    for (l = 0; l < nt; ++l)
		    {
		        if (!tri_recorded(tris[l],tri_list1,num_tris1))
			        tri_list1[num_tris1++] = tris[l];
				
		    }
		}
	    }

        j = (m == 0) ? box1->bmin[1] : upbound[1];
	for (k = box1->bmin[2]; k <= upbound[2]; ++k)
	{
	    for (i = box1->bmin[0]; i <= upbound[0]; ++i)
	    {
	        tris = T->tris[k][j][i];
		nt = T->num_of_tris[k][j][i];
		for (l = 0; l < nt; ++l)
		{
		    if (!tri_recorded(tris[l],tri_list1,num_tris1))
		            tri_list1[num_tris1++] = tris[l];

		}
	    }
	}

	i = (m == 0) ? box1->bmin[0] : upbound[0];
	for (k = box1->bmin[2]; k <= upbound[2]; ++k)
	{
	    for (j = box1->bmin[1]; j <= upbound[1]; ++j)
	    {
	        tris = T->tris[k][j][i];
		nt = T->num_of_tris[k][j][i];
		for (l = 0; l < nt; ++l)
		{
		    if (!tri_recorded(tris[l],tri_list1,num_tris1))
		            tri_list1[num_tris1++] = tris[l];
		}
	    }
	}
	}

	/* Form the tri_list on the boundary of box2 */

	for (i = 0; i < 3; ++i)
	    upbound[i] = (box2->bmax[i] == smax[i])? box2->bmax[i]-1 :
	                        box2->bmax[i];

	for (m = 0; m < 2; ++m)
	{
	    k = (m == 0) ? box2->bmin[2] : upbound[2];
	    for (j = box2->bmin[1]; j <= upbound[1]; ++j)
	    {
	        for (i = box2->bmin[0]; i <= upbound[0]; ++i)
		{
		    tris = T->tris[k][j][i];
		    nt = T->num_of_tris[k][j][i];
		    for (l = 0; l < nt; ++l)
		    {
		        if (!tri_recorded(tris[l],tri_list2,num_tris2))
			        tri_list2[num_tris2++] = tris[l];
		    }
		}
	    }


	j = (m == 0) ? box2->bmin[1] : upbound[1];
	for (k = box2->bmin[2]; k <= upbound[2]; ++k)
	{
	    for (i = box2->bmin[0]; i <= upbound[0]; ++i)
	    {
	        tris = T->tris[k][j][i];
		nt = T->num_of_tris[k][j][i];
		for (l = 0; l < nt; ++l)
		{
		    if (!tri_recorded(tris[l],tri_list2,num_tris2))
		            tri_list2[num_tris2++] = tris[l];

		}
	    }
	}

	i = (m == 0) ? box2->bmin[0] : upbound[0];
	for (k = box2->bmin[2]; k <= upbound[2]; ++k)
	{
	    for (j = box2->bmin[1]; j <= upbound[1]; ++j)
	    {
	        tris = T->tris[k][j][i];
		nt = T->num_of_tris[k][j][i];
		for (l = 0; l < nt; ++l)
		{
		    if (!tri_recorded(tris[l],tri_list2,num_tris2))
		            tri_list2[num_tris2++] = tris[l];

		}
	    }
	}
        }

	/* Comparing tri_list */
	if ((num_tris1 > 0) && (num_tris2 > 0))
	{
	    for (i = 0; i < num_tris1; ++i)
	        for (j = 0; j < num_tris2; ++j)
		    if ( tri_list1[i] == tri_list2[j])
		        flag = YES;
	}

	if(num_tris1 > MAX_BOUND_TRIS || num_tris2 > MAX_BOUND_TRIS)
	{
	    printf("ERROR boxes_sharing_tris, num %d  %d\n", 
	    	   num_tris1, num_tris2);
	    clean_up(ERROR);
	}
	
	if (flag)
	    printf("Two boxes are sharing tris, will be merged!\n");
	return flag;
}               /*end boxes_sharing_tris*/

LOCAL   boolean boxes_crxing_same_boundary(
        RECT_BOX  *box,
	RECT_BOX  *nbox,
	INTERFACE *intfc)
{
        int i;
	double min_coord,max_coord,nmin_coord,nmax_coord;
	double *L,crx_tol,*h;
	RECT_GRID gr = topological_grid(intfc);
	RECT_GRID comp_grid = Computational_grid(intfc);

	L = gr.L;
	h = gr.h;

	for (i = 0; i < 2; i++)
	{
	    min_coord = L[i] + box->bmin[i]*h[i];
	    max_coord = L[i] + box->bmax[i]*h[i];
	    nmin_coord = L[i] + nbox->bmin[i]*h[i];
	    nmax_coord = L[i] + nbox->bmax[i]*h[i];
            crx_tol = 0.004*h[i];

	    if ((min_coord <  comp_grid.L[i] + h[i]/2 + crx_tol) &&
	        (max_coord >  comp_grid.L[i] - h[i]/2 - crx_tol) &&
		(nmin_coord <  comp_grid.L[i] + h[i]/2 + crx_tol) &&
		(nmax_coord >  comp_grid.L[i] - h[i]/2 - crx_tol))
	    {
	        return YES;
	    }

	    if ((min_coord <  comp_grid.U[i] + h[i]/2 + crx_tol) &&
	        (max_coord >  comp_grid.U[i] - h[i]/2 - crx_tol) &&
		(nmin_coord <  comp_grid.U[i] + h[i]/2 + crx_tol) &&
		(nmax_coord >  comp_grid.U[i] - h[i]/2 - crx_tol))
	    {
	        return YES;
	    }
	    
	}

	return NO;
}

LOCAL	boolean ip_connected(
	int **ips,
	int nb,
	int *ip)
{
	int i,k,ni,na;
	for (i = 0; i < nb; ++i)
	{
	    ni = na = 0;
	    for (k = 0; k < 3; ++k)
	    {
	    	if (ips[i][k] == ip[k]) 
		    ni++;
		else if (ips[i][k] - ip[k] == 1 ||
		         ips[i][k] - ip[k] == -1) 
		    na++;
	    }
	    if (ni == 2 && na == 1) return YES;
	}
	return NO;
}	/* end ip_connected */

LOCAL   int  find_nearest_pts(
        POINT **sub_in_pts,
        int   end_in_pts,
        POINT **sub_out_pts,
        int start_out_pts,
        int end_out_pts,
        TRI **sub_in_tris,
        TRI **sub_out_tris,
        RECT_GRID rgr,
	RECT_GRID comp_gr)
{
        int i,j,m,n;
        int start_pts = 1000;
        POINT *p0,*p1,*p2,*p3;
        double dist,min_dist = HUGE;
        double *L,*h,crx_coord[2];
        boolean pts_in_match = YES;
	double r0,r1,u0[3],u1[3];
	double lgbtol = 10*MACH_EPS;

        L = rgr.L;
        h = rgr.h;

        p0 = sub_in_pts[end_in_pts];
        p1 = sub_in_pts[end_in_pts-1];
	r0 = distance_between_positions(Coords(p0), Coords(p1), 3);
	for (i = 0; i < 3; ++i)
	    u0[i] = (Coords(p1)[i] - Coords(p0)[i])/r0;
        
	for (i = start_out_pts; i < end_out_pts; i++)
        {
            
	    p2 = sub_out_pts[i];
	    p3 = sub_out_pts[i+1];
	    r1 = distance_between_positions(Coords(p2), Coords(p3), 3);
	    for (j = 0; j < 3; ++j)
	        u1[j] = (Coords(p3)[j] - Coords(p2)[j])/r1;
	    cosu = Dot3d(u0, u1);
	    if (cosu >= -lgbtol)
	        continue;
	    pts_in_match = YES;
            for (n = 0; n < 3; n++)
            {
                crx_coord[0] = comp_gr.L[n];
                crx_coord[1] = comp_gr.U[n];
                if (tri_cross_zone(sub_in_tris[end_in_pts-1],crx_coord[0],h[n],n) == YES &&
                    tri_cross_zone(sub_in_tris[end_in_pts-1],crx_coord[1],h[n],n) == YES)
                    if (tri_cross_zone(sub_out_tris[start_out_pts],crx_coord[0],h[n],n) == NO ||
                       tri_cross_zone(sub_out_tris[start_out_pts],crx_coord[1],h[n],n) == NO)
                    {
                        pts_in_match = NO;
                        break;
                    }
                for (m = 0; m < 2; m++)
                {
                    if (tri_cross_zone(sub_in_tris[end_in_pts-1],crx_coord[m],h[n],n) == YES)
                        if (tri_cross_zone(sub_out_tris[start_out_pts],crx_coord[m],h[n],n) == NO)
                        {
                            pts_in_match = NO;
                            break;
                        }
                    if (tri_cross_line1(sub_in_tris[end_in_pts-1],crx_coord[m],h[n],n) == YES)
                        if (tri_cross_line1(sub_out_tris[start_out_pts],crx_coord[m],h[n],n) == NO)
                        {
                            pts_in_match = NO;
                            break;
                        }
                    if (tri_cross_line1(sub_in_tris[end_in_pts-1],crx_coord[m],-h[n],n) == YES)
                        if (tri_cross_line1(sub_out_tris[start_out_pts],crx_coord[m],-h[n],n) == NO)
                        {
                            pts_in_match = NO;
                            break;
                        }
                }
            }
            if (!pts_in_match) continue;
            dist = distance_between_positions(Coords(p0),Coords(p2),3) + 
	           distance_between_positions(Coords(p0),Coords(p3),3) +
		   distance_between_positions(Coords(p1),Coords(p2),3) +
		   distance_between_positions(Coords(p1),Coords(p3),3);
            if (dist < min_dist)
            {
                start_pts = i;
                min_dist = dist;
            }
        }

        return start_pts;
}  /*end find_nearest_pts*/

LOCAL   int  find_nearest_pts_pair(
        POINT **sub_in_pts,
	int   end_in_pts,
	POINT **sub_out_pts,
	int num_sub_out_pts,
	int start_out_pts,
	RECT_GRID rgr,
	RECT_GRID comp_gr)
{
        int i,j,m,n;
	int end_out_pts,start_pts = start_out_pts;
	POINT *p0,*p1;
	double dist,min_dist = HUGE;
	double *L,*h,crx_coord[2];
	boolean pts_in_match = YES;

	L = rgr.L;
	h = rgr.h;
	
	end_out_pts = 10000;

	p0 = sub_in_pts[end_in_pts];
	/*
	if (sub_in_pts[end_in_pts+2] == sub_in_pts[0] &&
	    start_out_pts <= 1)
	{
	    printf("end_in_pts = %d\n",end_in_pts);
	    start_pts = start_out_pts+1;
	    printf("start_pts = %d\n",start_pts);
	}
	*/
	for (i = start_pts; i < num_sub_out_pts; i++)
	{
	    pts_in_match = YES;
	    p1 = sub_out_pts[i];
	    for (n = 0; n < 3; n++)
	    {
	        crx_coord[0] = comp_gr.L[n];
		crx_coord[1] = comp_gr.U[n];
		if (pts_in_cross_zone(Coords(p0),crx_coord[0],h[n],n) == YES &&
		    pts_in_cross_zone(Coords(p0),crx_coord[1],h[n],n) == YES)
		    if (pts_in_cross_zone(Coords(p1),crx_coord[0],h[n],n) == NO ||
		       pts_in_cross_zone(Coords(p0),crx_coord[1],h[n],n) == NO)
		    {
		        pts_in_match = NO;
		        break;
		    }
		for (m = 0; m < 2; m++)
		{
		    if (pts_in_cross_zone(Coords(p0),crx_coord[m],h[n],n) == YES)
		        if (pts_in_cross_zone(Coords(p1),crx_coord[m],h[n],n) == NO)
	                {
			    pts_in_match = NO;
			    break;
			}
		    if (pts_cross_line(Coords(p0),crx_coord[m],h[n],n) == YES)
		        if (pts_cross_line(Coords(p1),crx_coord[m],h[n],n) == NO)
			{
			    pts_in_match = NO;
			    break;
			}
	            if (pts_cross_line(Coords(p0),crx_coord[m],-h[n],n) == YES)
		        if (pts_cross_line(Coords(p1),crx_coord[m],-h[n],n) == NO)
			{
			    pts_in_match = NO;
			    break;
			}
		}
	    }
	    if (!pts_in_match) continue;
	    dist = distance_between_positions(Coords(p0),Coords(p1),3);
	    if (dist < min_dist)
	    {
	        end_out_pts = i;
		min_dist = dist;
	    }
	}
	return end_out_pts;
}  /*end find_nearest_pts_pair*/


EXPORT	boolean	check_degenerated_loop(
	TRI	**tris,
	int	*sides,
	POINT	**pts,
	int	num_tris)
{
	TRI	*new_tri;
	SURFACE	*surf;

	if(num_tris > 3)
	    return NO;
	
	if(num_tris < 3)
	{
	    printf("ERROR in check_degenerated_loops, num_tris = %d.\n", num_tris);
	    clean_up(ERROR);
	}

	surf = tris[0]->surf;
	
	if (tris[0] == tris[1] && tris[1] == tris[2])
	{
	    remove_tri_from_surface(tris[0],surf,NO);
	    tmp_n = 0;
        }
	else
	{
	    /* assume the null loop has POSITIVE_ORIENTATION */
	    
	    new_tri = make_tri(pts[0],pts[2],pts[1], NULL,NULL,NULL, NO);
	    link_neighbor_tris(new_tri, tris[0]);
	    link_neighbor_tris(new_tri, tris[1]);
	    link_neighbor_tris(new_tri, tris[2]);

	    tmp_tris[0] = new_tri;
	    tmp_n = 1;
	    
	    insert_tri_at_tail_of_list(new_tri,surf);
	}

	return YES;
}


	boolean null_side_loop(
	TRI *start,
	int start_side,
	ORIENTATION orient,
	TRI ***tris,
	int **sides,
	POINT ***side_points,
	int *num_sides,
	double **tnor)
{
	POINT *p;
	int i, j, k, i_diff, side = start_side;
	TRI *next_tri = start;
	static TRI **null_tris = NULL;
	static int *null_sides =  NULL;
	static POINT **pts;
	static double normal[MAXD];
	double v1[MAXD],v2[MAXD],cprod[MAXD];

	DEBUG_ENTER(null_side_loop)

	if (null_tris == NULL)
	{
	    
	    uni_array(&null_tris,MAX_NULL_SIDE_LOOP,sizeof(TRI*));
	    uni_array(&null_sides,MAX_NULL_SIDE_LOOP,INT);
	    uni_array(&pts,MAX_NULL_SIDE_LOOP,sizeof(POINT*));
	}
	for (i = 0; i < 3; ++i) normal[i] = 0.0;
	*num_sides = 0;
	pts[0] = (orient == POSITIVE_ORIENTATION) ? Point_of_tri(start)[side] :
			Point_of_tri(start)[Next_m3(side)];
	
	do
	{
	    null_tris[*num_sides] = next_tri;
	    null_sides[*num_sides] = side;
	    ++(*num_sides);
	    if (orient == POSITIVE_ORIENTATION)
	    {
		if (debugging("null_loop"))
		{
	    	    printf("side:  %d\n",side);
		    print_point(Point_of_tri(next_tri)[side]);
		    print_point(Point_of_tri(next_tri)[Next_m3(side)]);
		}
		pts[*num_sides] = p = Point_of_tri(next_tri)[Next_m3(side)];
		if (*num_sides > 300) 
		{
	    	    int i;
	    	    for (i = 0; i < *num_sides; ++i)
	    	        printf("%d %d %d\n",2,i,i+1);
		    printf("%d %d %d\n",*num_sides,*num_sides,0);
		    clean_up(ERROR);
		}
		if (p == pts[0]) break;  
			/*the only way to exit do while except 
			  return NO, pts must loop
			*/
	    	side = next_null_sided_tri(next_tri,p,&next_tri);
		if (side == -1)
		{
		    DEBUG_LEAVE(null_side_loop)
		    return NO;
		}
	    }
	    else
	    {
		if (debugging("null_loop"))
		{
	    	    printf("side:  %d\n",side);
		    print_point(Point_of_tri(next_tri)[Next_m3(side)]);
		    print_point(Point_of_tri(next_tri)[side]);
		}
		pts[*num_sides] = p = Point_of_tri(next_tri)[side];
		if (p == pts[0]) break;
	    	side = prev_null_sided_tri(next_tri,p,&next_tri);
	    }
	    if (next_tri == NULL)
	    {
		printf("Cannot continue the loop\n");
		clean_up(ERROR);
	    }

	}
	while (next_tri != start || side != start_side);

	/* Check for double loop */
	for (i = 1; i < *num_sides-1; ++i)
	{
	    for (j = i+1; j < *num_sides; ++j)
	    {
	    	if (pts[i] == pts[j])
		{
		    i_diff = j - i;
		    for (k = i; k < *num_sides-i_diff; ++k)
		    {
	    		null_tris[k] = null_tris[k+i_diff];
	        	null_sides[k] = null_sides[k+i_diff];
			pts[k] = pts[k+i_diff];
		    }
	    	    pts[*num_sides-i_diff] = pts[*num_sides];
		    (*num_sides) -= i_diff;
		    --i;
		    break;
		}
	    }
	}
	

	for (j = 0; j < 3; ++j) normal[j] = 0.0;
	for (i = 0; i < *num_sides; ++i)
	{
	    difference(Coords(pts[i]),Coords(pts[i+1]),v1,3);
	    difference(Coords(pts[(i+2)%(*num_sides+1)]),
	    		Coords(pts[i+1]),v2,3);
	    Cross3d(v2,v1,cprod);
	    for (j = 0; j < 3; ++j) normal[j] += cprod[j];
	}

	if (orient == NEGATIVE_ORIENTATION)
	    for (j = 0; j < 3; ++j) normal[j] *= -1.0;
	*tris = null_tris;
	*sides = null_sides;
	*side_points = pts;
	*tnor = normal;
	
	DEBUG_LEAVE(null_side_loop)
	return YES;
}	/* end null_side_loop */


LOCAL	boolean null_side_tri_in_list(
	TRI **tris,
	int num_tris,
	TRI **tri,
	int *side)
{
	int i,j;
	POINT **test;

	
	for (i = 0; i < num_tris; ++i)
	{
	    for (j = 0; j < 3; ++j)
	    {
	    	if (Neighbor_on_side(tris[i],j) == NULL)
		{
		    *tri = tris[i];
		    *side = j;
		    return YES;
		}
	    }
	}
	*tri = NULL;
	return NO;
}	/* end null_side_tri_in_list */

LOCAL	boolean check_extension_of_surface(
	TRI **tri_list,
	int num_tris,
	SURFACE *s)
{
	TRI *tri1,*tri2,*t[2];
	POINT *p1,*p2;
	int i,j,k,count;
	
	printf("Entering check_extension_of_surface!\n");
	
	for (k = 0; k < num_tris; ++k)
	{
	    tri1 = tri_list[k];
	    for (i = 0; i < 3; ++i)
	    {
		p1 = Point_of_tri(tri1)[i];
		p2 = Point_of_tri(tri1)[(i+1)%3];
		count = 0;
	    	for (tri2 = first_tri(s); !at_end_of_tri_list(tri2,s);
	    	    	tri2 = tri2->next)
	    	{
		    if (tri1 == tri2) continue;
		    for (j = 0; j < 3; ++j)
		    {
		    	if ((p1 == Point_of_tri(tri2)[j] &&
			     p2 == Point_of_tri(tri2)[(j+1)%3]) ||
			    (p1 == Point_of_tri(tri2)[(j+1)%3] &&
			     p2 == Point_of_tri(tri2)[j]))
			{
			    t[count++] = tri2;
			    if (count == 2)
			    {
			    	printf("Triple edge found\n");
				printf("Triangle 1:\n");
				print_tri(tri1,s->interface);
				printf("Triangle 2:\n");
				print_tri(t[0],s->interface);
				printf("Triangle 3:\n");
				print_tri(t[1],s->interface);
				clean_up(ERROR);
			    }
			}
		    }
		}
	    }
	}
	printf("Surface contains no triple edge "
	       "and is topologically sound!\n");
	return YES;
}	/* end check_extension_of_surface */


EXPORT  boolean check_extension_of_surface_global(
        SURFACE *s)
{
        TRI *tri1,*tri2,*t[2];
        POINT *p1,*p2;
        int i,j,count;

        printf("Entering check_extension_of_surface_global!\n");

        for (tri1 = first_tri(s); !at_end_of_tri_list(tri1->next, s);
                tri1 = tri1->next)
        {
            for (i = 0; i < 3; ++i)
            {
                p1 = Point_of_tri(tri1)[i];
                p2 = Point_of_tri(tri1)[(i+1)%3];
                count = 0;
                for (tri2 = tri1->next; !at_end_of_tri_list(tri2,s);
                        tri2 = tri2->next)
                {
                    if (tri1 == tri2) continue;
                    for (j = 0; j < 3; ++j)
                    {
                        if ((p1 == Point_of_tri(tri2)[j] &&
                             p2 == Point_of_tri(tri2)[(j+1)%3]) ||
                            (p1 == Point_of_tri(tri2)[(j+1)%3] &&
                             p2 == Point_of_tri(tri2)[j]))
                        {
                            t[count++] = tri2;
                            if (count == 2)
                            {
                                printf("Triple edge found\n");
                                printf("Triangle 1:\n");
                                print_tri(tri1,s->interface);
                                printf("Triangle 2:\n");
                                print_tri(t[0],s->interface);
                                printf("Triangle 3:\n");
                                print_tri(t[1],s->interface);
                                clean_up(ERROR);
                            }
                        }
                    }
                }
            }
        }
        printf("Surface contains no triple edge "
               "and is topologically sound!\n");
        return YES;
}       /* end check_extension_of_surface_global */



LOCAL	boolean null_sides_with_suitable_angle(
	TRI *tri1,
	int side1,
	TRI *tri2,
	int side2)
{
	double *p0, *p1, *q0, *q1;
	double u0[3], u1[3], pm[3], qm[3], pq[3], tm[3];
	double r0, r1, ke, dm, r2,r3;
	double r4,r5,d,d2,d3,d4,d5;
	double u2[3],u3[3],u4[3],u5[3];
	int i,pi[3],po[3];
	double lgbtol = 10*MACH_EPS,thres = -0.5;
	
	p0 = Coords(Point_of_tri(tri1)[side1]);
	p1 = Coords(Point_of_tri(tri1)[Next_m3(side1)]);
	q0 = Coords(Point_of_tri(tri2)[side2]);
	q1 = Coords(Point_of_tri(tri2)[Next_m3(side2)]);
	r0 = distance_between_positions(p0, p1, 3);
	r1 = distance_between_positions(q0, q1, 3);
	r2 = distance_between_positions(p0, q0, 3);
	r3 = distance_between_positions(p1, q0, 3);
	r4 = distance_between_positions(q0, p0, 3);
	r5 = distance_between_positions(q1, p0, 3);
	
	for (i=0; i<3; ++i)
	{
	    pm[i] = (p1[i] + p0[i])/2.0;
	    qm[i] = (q1[i] + q0[i])/2.0;
	    u0[i] = (p1[i] - p0[i])/r0;
	    u1[i] = (q1[i] - q0[i])/r1;
	    u2[i] = (p0[i] - q0[i])/r2;
	    u3[i] = (p1[i] - q0[i])/r3;
	    u4[i] = (q0[i] - p0[i])/r4;
	    u5[i] = (q1[i] - p0[i])/r5;
	    tm[i] = pm[i] - q0[i];
	    pi[i] = (int)(p0[i] - q0[i]);
	    po[i] = (int)(q0[i] - p0[i]);
	}
	ke = Dot3d(tm, u1);
	for (i=0; i<3; ++i)
	    pq[i] = q0[i] + ke*u1[i];
	
	ke = Dot3d(u2,u1);
	d2 = r2*sqrt(1-ke*ke);
	ke = Dot3d(u3,u1);
	d3 = r3*sqrt(1-ke*ke);
	ke = Dot3d(u4,u0);
	d4 = r4*sqrt(1-ke*ke);
	ke = Dot3d(u5,u0);
	d5 = r5*sqrt(1-ke*ke);
	d = min(d2,d3);
	d = min(d,d4);
	d = min(d,d5);
	
	/*dm is the distance between the midpoint of side2 and the projection
	  of the mid point of side1 on side2*/
	dm = distance_between_positions(pq, qm, 3);
	
	cosu = Dot3d(u0, u1);
	
	if (cosu <= -lgbtol)
	    return YES;
	else
	    return NO;

	/*cosu > 0 ? NO: YES;*/
}	/* end null_sides_with_suitable_angle */	

LOCAL   boolean bifurcation_detected(
	TRI **tri_list,
	int num_tris,
	TRI ***sub_tris1,
	int *num_sub_tris1,
	TRI ***sub_tris2,
	int *num_sub_tris2)
{
	TRI **null_tris;
	static TRI **sub_list1 = NULL, **sub_list2 = NULL;
	TRI **tmp_tris,**mom_list;
	TRI *start_tri,**recorded_tris;
	int i, start_side, num_null_sides, num_mom_tris,num_sub_list1, 
	    num_sub_list2, *null_sides=0,  count=0, min_num = 1000,
	    max_num = -1000,num_tmp_tris=0;
	int num_sub_tris[10],j,num_recorded_tris;
	POINT **side_pts;
	double *tnor;
	boolean bifurcation_found = NO;

	
	uni_array(&mom_list, 3*num_tris, sizeof(TRI*));
	if (sub_list1 == NULL)
	{
	    uni_array(&sub_list1, MAX_NULL_SIDE_LOOP, sizeof(TRI*));
	    uni_array(&sub_list2, MAX_NULL_SIDE_LOOP, sizeof(TRI*));
	}
	uni_array(&tmp_tris, 3*num_tris, sizeof(TRI*));
	uni_array(&recorded_tris, 3*num_tris, sizeof(TRI*));



	num_sub_list1 = num_sub_list2 = 0;
	for (i = 0; i < num_tris; ++i)
		mom_list[i] = tri_list[i];
	num_mom_tris = num_tris;

	num_recorded_tris = 0;
	while (null_side_tri_in_list(mom_list,num_mom_tris,
				&start_tri,&start_side) )
	{
	    
	    null_side_loop(start_tri,start_side,POSITIVE_ORIENTATION,
	        &null_tris,&null_sides,&side_pts,&num_null_sides,&tnor);
	    
	    
	    for (i = 0; i < num_null_sides; i++)
	        if (!tri_recorded(null_tris[i],recorded_tris,num_recorded_tris))
	            recorded_tris[num_recorded_tris++] = null_tris[i];
	    if (num_null_sides < min_num)
	    {
	        for (i = 0; i < num_null_sides; i++)
		    sub_list1[i] = null_tris[i];
		num_sub_list1 = num_null_sides;
		min_num = num_null_sides;
	    }
	    if (num_null_sides > max_num)
	    {
	        for (i = 0; i < num_null_sides; i++)
		    sub_list2[i] = null_tris[i];
		num_sub_list2 = num_null_sides;
		max_num = num_null_sides;
	    }
	    count++;
	    
	    num_mom_tris = 0;
	    for (i = 0; i < num_tris; i++)
	    {
		if (!tri_recorded(tri_list[i],recorded_tris,num_recorded_tris))
		{
		    mom_list[num_mom_tris++] = tri_list[i];
		}
	    }
	}

	
	printf("num_sub_tris2 = %d\n",num_sub_list2);
	printf("num_sub_tris1 = %d\n",num_sub_list1);
	
	/*
	if (count > 2)
	{
	    screen("Bifurcation:More than two null_side_loop in tri_list!\n");
	    clean_up(ERROR);
	}
	else if (count == 2)
		bifurcation_found = YES;
	else
		bifurcation_found = NO;	
	if (count >= 2)
	    bifurcation_found = YES;
	else
	*/
	bifurcation_found = YES;
	printf("count = %d\n",count);
	*sub_tris1 = sub_list1;
	*num_sub_tris1 = num_sub_list1;
	*sub_tris2 = sub_list2;
	*num_sub_tris2 = num_sub_list2;

	free_these(3,recorded_tris,mom_list,tmp_tris);
	return bifurcation_found;
}	/* end bifurcation_detected */ 


LOCAL   boolean    pts_in_cross_zone(
        double         *pt,
        double          crx_coord,
        double          h,
        int             dir)
{
        double   min_coord,max_coord;
        double   crx_tol = 0.004*h;
        int     i;
        
	min_coord = max_coord = pt[dir];

        if (((min_coord - crx_coord) <= h+crx_tol) &&
            ((crx_coord - max_coord) <= h+crx_tol))
            return YES;
        else
            return NO;
}    /*end tri_cross_zone*/

/*if one coord of dir direction for a point of a tri lies in the 
  zone [crx_coord - h, crx_coord + h] it will return YES;
*/

LOCAL   boolean    tri_cross_zone(
        TRI             *tri,
	double           crx_coord,
	double           h,
	int             dir)
{
        double   min_coord,max_coord;
        double   crx_tol = 0.004*h;	
        int     i;
	
	min_coord = max_coord = Coords(Point_of_tri(tri)[0])[dir];

	if (min_coord > Coords(Point_of_tri(tri)[1])[dir])
	    min_coord = Coords(Point_of_tri(tri)[1])[dir];
	if (min_coord > Coords(Point_of_tri(tri)[2])[dir])
	    min_coord = Coords(Point_of_tri(tri)[2])[dir];
	
	if (max_coord < Coords(Point_of_tri(tri)[1])[dir])
	    max_coord = Coords(Point_of_tri(tri)[1])[dir];
	if (max_coord < Coords(Point_of_tri(tri)[2])[dir])
	    max_coord = Coords(Point_of_tri(tri)[2])[dir];
	
	if (((min_coord - crx_coord) <= h+crx_tol) &&
	    ((crx_coord - max_coord) <= h+crx_tol))
	    return YES;
	else
	    return NO;
}    /*end tri_cross_zone*/


LOCAL   boolean    pts_cross_line(
        double           *pt,
        double           crx_coord,
        double           h,
        int             dir)
{
        double   min_coord, max_coord;
        double   crx_tol = 0.004*h;

        min_coord = max_coord = pt[dir];

        if (h > 0)
            return ((min_coord - crx_coord) >= h+crx_tol) ? YES : NO;
        else
            return ((crx_coord - max_coord) >= -h+crx_tol) ? YES : NO;
}


LOCAL   boolean    tri_cross_line1(
        TRI             *tri,
	double           crx_coord,
	double           h,
	int             dir)
{
        double   min_coord, max_coord;
	double   crx_tol = 0.004*h;

	min_coord = max_coord = Coords(Point_of_tri(tri)[0])[dir];

	if (min_coord > Coords(Point_of_tri(tri)[1])[dir])
	    min_coord = Coords(Point_of_tri(tri)[1])[dir];
	if (min_coord > Coords(Point_of_tri(tri)[2])[dir])
	    min_coord = Coords(Point_of_tri(tri)[2])[dir];
	
	if (max_coord < Coords(Point_of_tri(tri)[1])[dir])
	    max_coord = Coords(Point_of_tri(tri)[1])[dir];
	if (max_coord < Coords(Point_of_tri(tri)[2])[dir])
	    max_coord = Coords(Point_of_tri(tri)[2])[dir];
        
        if (h > 0)
	    return ((min_coord - crx_coord) >= h+crx_tol) ? YES : NO;
	else
	    return ((crx_coord - max_coord) >= -h+crx_tol) ? YES : NO;
}

LOCAL   void    initialize_comm_box(
        COMM_BOX         *comm_box,
	boolean             set_flag)
{
	int    inter=1000,i,j,k;
	
	for (i = 0; i < 4; i++)
	    for (k = 0; k < 2; k++)
	        for (j = 0; j < 3; j++)
	        {
	            comm_box->lmin[i][k][j] = inter;
		    comm_box->lmax[i][k][j] = -inter;
		    comm_box->umin[i][k][j] = inter;
		    comm_box->umax[i][k][j] = -inter;
	        }
	comm_box->L0 = HUGE;
	comm_box->U0 = -HUGE;
	comm_box->L1 = HUGE;
	comm_box->U1 = -HUGE;
	comm_box->lx = 0;
	comm_box->ly = 0;
	comm_box->ux = 0;
	comm_box->uy = 0;
	if (!set_flag)
	    comm_box->flag = 0;

}    /*end initailize_comm_box*/

LOCAL   void    find_comm_box(
        INTERFACE	 *intfc,
	RECT_BOX         *box,
        int              *gmax,
	COMM_BOX         *comm_box)
{
	RECT_GRID gr = topological_grid(intfc);
	RECT_GRID comp_grid = Computational_grid(intfc);
        double *L,*h,crx_tol,min_coord,max_coord;
	int i,j,check_crx_buffer = 0,check_crx_boundary;

	L = gr.L;
	h = gr.h;

        check_crx_buffer = 0;
	check_crx_boundary = 0;
	for (i = 0; i < 3; i++)
        {
            min_coord = L[i] + box->bmin[i]*h[i];
            max_coord = L[i] + box->bmax[i]*h[i];
	    crx_tol = 0.004*h[i];
            if (((min_coord <  comp_grid.L[i] + 3*h[i]/2+crx_tol) &&
                (max_coord >  comp_grid.L[i] - 3*h[i]/2-crx_tol)) ||
                ((min_coord <  comp_grid.U[i] + 3*h[i]/2+crx_tol) &&
                (max_coord >  comp_grid.U[i] - 3*h[i]/2-crx_tol)))
                check_crx_buffer++;
	    if (box->bmin[i] == 0 || box->bmax[i] == gmax[i])
	        check_crx_boundary++;
        }

        /*
	if (check_crx_buffer >= 2 && check_crx_boundary >= 2)
            comm_box->flag = 1;
        */
	
	for (i = 0; i < 2; i++)
	{
	    min_coord = L[i] + box->bmin[i]*h[i];
	    max_coord = L[i] + box->bmax[i]*h[i];
	    crx_tol = 0.004*h[i];

	    if ((min_coord <  comp_grid.L[i] + 3*h[i]/2+crx_tol) &&
	        (max_coord >  comp_grid.L[i] - 3*h[i]/2-crx_tol))
	    {
	        if (i == 0)
		{
	            comm_box->L0 = comp_grid.L[i];
		    comm_box->lmin[comm_box->lx][0][0] = box->bmin[0];
		    comm_box->lmax[comm_box->lx][0][0] = box->bmax[0];
		    comm_box->lmin[comm_box->lx][0][1] = box->bmin[1];
		    comm_box->lmax[comm_box->lx][0][1] = box->bmax[1];
		    comm_box->lmin[comm_box->lx][0][2] = box->bmin[2];
		    comm_box->lmax[comm_box->lx][0][2] = box->bmax[2];
		    comm_box->lx++;
		}
		if (i == 1)
		{
		    comm_box->L1 = comp_grid.L[i];
		    comm_box->lmin[comm_box->ly][1][1] = box->bmin[1];
		    comm_box->lmax[comm_box->ly][1][1] = box->bmax[1];
		    comm_box->lmin[comm_box->ly][1][0] = box->bmin[0];
		    comm_box->lmax[comm_box->ly][1][0] = box->bmax[0];
		    comm_box->lmin[comm_box->ly][1][2] = box->bmin[2];
		    comm_box->lmax[comm_box->ly][1][2] = box->bmax[2];
		    comm_box->ly++;
		}
	    }
	    if ((min_coord <  comp_grid.U[i] + 3*h[i]/2+crx_tol) &&
	        (max_coord >  comp_grid.U[i] - 3*h[i]/2-crx_tol))
	    {
	        if (i == 0)
		{
		    comm_box->U0 = comp_grid.U[i];
		    comm_box->umin[comm_box->ux][0][0] = box->bmin[0];
		    comm_box->umax[comm_box->ux][0][0] = box->bmax[0];
		    comm_box->umin[comm_box->ux][0][1] = box->bmin[1];
		    comm_box->umax[comm_box->ux][0][1] = box->bmax[1];
		    comm_box->umin[comm_box->ux][0][2] = box->bmin[2];
		    comm_box->umax[comm_box->ux][0][2] = box->bmax[2];
		    comm_box->ux++;
		}
		if (i == 1)
		{
		    comm_box->U1 = comp_grid.U[i];
		    comm_box->umin[comm_box->uy][1][1] = box->bmin[1];
		    comm_box->umax[comm_box->uy][1][1] = box->bmax[1];
		    comm_box->umin[comm_box->uy][1][0] = box->bmin[0];
		    comm_box->umax[comm_box->uy][1][0] = box->bmax[0];
		    comm_box->umin[comm_box->uy][1][2] = box->bmin[2];
		    comm_box->umax[comm_box->uy][1][2] = box->bmax[2];
		    comm_box->uy++;
		}
	    }
	}
}  /*end find_comm_box */

LOCAL    boolean    compare_comm_box(
        INTERFACE	 *intfc,
	RECT_BOX         *box,
	COMM_BOX         *adj_box_l,
	COMM_BOX         *adj_box_u,
	int              *gmax,
	int              dir)
{
	RECT_GRID gr = topological_grid(intfc);
	RECT_GRID comp_grid = Computational_grid(intfc);
        double *L,*h,crx_tol,min_coord, max_coord;
	int i,j,k,test_box,check_crx_buffer = 0,
	     check_crx_boundary = 0;
	double icoord;
	int comm_grid_size,grid_size;

	L = gr.L;
	h = gr.h;


	for (i = 0; i < 3; i++)
	{
	    min_coord = L[i] + box->bmin[i]*h[i];
	    max_coord = L[i] + box->bmax[i]*h[i];
	    crx_tol = 0.004*h[i];
	    comm_grid_size = gmax[i] - comp_grid.lbuf[i] -
	                     comp_grid.ubuf[i] + 1;
	    if (((min_coord <  comp_grid.L[i] + 3*h[i]/2+crx_tol) &&
	        (max_coord >  comp_grid.L[i] - 3*h[i]/2-crx_tol)) ||
		((min_coord <  comp_grid.U[i] + 3*h[i]/2+crx_tol) &&
		(max_coord >  comp_grid.U[i] - 3*h[i]/2-crx_tol)))
	        check_crx_buffer++;
	    if (box->bmin[i] == 0 || box->bmax[i] == gmax[i])
	        check_crx_boundary++;
	}

	if (check_crx_buffer >= 2 && check_crx_boundary >= 2)
	    return YES;
	for (i = 0; i < 2; i++)
	{
	    min_coord = L[i] + box->bmin[i]*h[i];
	    max_coord = L[i] + box->bmax[i]*h[i];
	    crx_tol = 0.004*h[i];

	    if (((min_coord <  comp_grid.L[i] + 3*h[i]/2+crx_tol) &&
	        (max_coord >  comp_grid.L[i] - 3*h[i]/2-crx_tol)) ||
		((min_coord <  comp_grid.U[i] + 3*h[i]/2+crx_tol) &&
		(max_coord >  comp_grid.U[i] - 3*h[i]/2-crx_tol)))
	    {
	        comm_grid_size = gmax[i] - comp_grid.lbuf[i] -
		                 comp_grid.ubuf[i] + 1;
		if ((min_coord <  comp_grid.L[i] + 3*h[i]/2+crx_tol) &&
		    (max_coord >  comp_grid.L[i] - 3*h[i]/2-crx_tol))
		{
		    icoord = comp_grid.L[i];
		    if (i == 0 && dir == 0)
		    {
		        if (fabs(icoord - adj_box_u->U0) < crx_tol ||
		            (fabs(icoord - adj_box_u->U0) > 2-crx_tol &&
			    fabs(icoord - adj_box_u->U0) < 2+crx_tol))
		        {
		            test_box = abs(adj_box_u->umax[0][0][1] - 
				  box->bmax[1]) + abs(adj_box_u->umax[0][0][2] 
				  - box->bmax[2]);
			    k = 0;
		            for (j = 1; j < adj_box_u->ux; j++)
			    {
			        if (test_box > fabs(adj_box_u->umax[j][0][1] - 
				  box->bmax[1]) + fabs(adj_box_u->umax[0][0][2] 
				  - box->bmax[2]))
			        {
			            test_box = abs(adj_box_u->umax[j][0][1] - 
				 	box->bmax[1]) + abs(
					adj_box_u->umax[0][0][2] - box->bmax[2]);
				    k = j;
			        }
			    }
			    grid_size = min(box->bmin[0],adj_box_u->umin[k][0][0] -
			                           comm_grid_size);
			    box->bmin[0] = max(0,grid_size);
			    box->bmax[0] = max(box->bmax[0], adj_box_u->umax[k][0][0] -
			                   comm_grid_size);
			    box->bmin[1] = min(box->bmin[1],adj_box_u->umin[k][0][1]);
		            box->bmax[1] = max(box->bmax[1],adj_box_u->umax[k][0][1]);
		            box->bmin[2] = min(box->bmin[2],adj_box_u->umin[k][0][2]);
		            box->bmax[2] = max(box->bmax[2],adj_box_u->umax[k][0][2]);
		        }

		        if (fabs(icoord - adj_box_l->U0) < crx_tol ||
		            (fabs(icoord - adj_box_l->U0) > 2-crx_tol &&
			    fabs(icoord - adj_box_l->U0) < 2+crx_tol))
	                {
		            test_box = abs(adj_box_l->umax[0][0][1] - box->bmax[1]) +
			               abs(adj_box_l->umax[0][0][2] - box->bmax[2]);
			    k = 0;
			    for (j = 1; j < adj_box_l->ux; j++)
			        if (test_box > fabs(adj_box_l->umax[j][0][1] - box->bmax[1]) +
				              fabs(adj_box_l->umax[j][0][2] - box->bmax[2]))
			        {
			            test_box = abs(adj_box_l->umax[j][0][1] - box->bmax[1]) +
				                abs(adj_box_l->umax[j][0][2] - box->bmax[2]);
				    k = j;
			        }
			    grid_size = min(box->bmin[0],adj_box_l->umin[k][0][0] -
			                           comm_grid_size);
			    box->bmin[0] = max(0,grid_size);
			    box->bmax[0] = max(box->bmax[0], adj_box_l->umax[k][0][0] -
			                   comm_grid_size);
			    box->bmin[1] = min(box->bmin[1],adj_box_l->umin[k][0][1]);
			    box->bmax[1] = max(box->bmax[1],adj_box_l->umax[k][0][1]);
			    box->bmin[2] = min(box->bmin[2],adj_box_l->umin[k][0][2]);
			    box->bmax[2] = max(box->bmax[2],adj_box_l->umax[k][0][2]);
		        }
		    }
		    if (i == 1 && dir == 1)
		    {
		        if (fabs(icoord - adj_box_u->U1) < crx_tol ||
		             (fabs(icoord - adj_box_u->U1) > 2-crx_tol &&
			     fabs(icoord - adj_box_u->U1) < 2+crx_tol))
		        {
		            test_box = abs(adj_box_u->umax[0][1][0] - box->bmax[0]) +
			               abs(adj_box_u->umax[0][1][2] - box->bmax[2]);
			    k = 0;
			    for (j = 1; j < adj_box_u->uy; j++)
			        if (test_box > fabs(adj_box_u->umax[j][1][0] - box->bmax[0]) +
				    fabs(adj_box_u->umax[0][1][2] - box->bmax[2]))
			        {
			            test_box = abs(adj_box_u->umax[j][1][0] - box->bmax[0]) +
				          abs(adj_box_u->umax[0][1][2] - box->bmax[2]);
				    k = j;
			        }
		            grid_size = min(box->bmin[1],adj_box_u->umin[k][1][1] -
                                                   comm_grid_size);
                            box->bmin[1] = max(0,grid_size);
                            box->bmax[1] = max(box->bmax[1], adj_box_u->umax[k][1][1] -
                                           comm_grid_size);
			    box->bmin[0] = min(box->bmin[0],adj_box_u->umin[k][1][0]);
		            box->bmax[0] = max(box->bmax[0],adj_box_u->umax[k][1][0]);
		            box->bmin[2] = min(box->bmin[2],adj_box_u->umin[k][1][2]);
		            box->bmax[2] = max(box->bmax[2],adj_box_u->umax[k][1][2]);
		        }
		        if (fabs(icoord - adj_box_l->U1) < crx_tol ||
		            (fabs(icoord - adj_box_l->U1) > 2-crx_tol &&
			    fabs(icoord - adj_box_l->U1) < 2+crx_tol))
		        {
		            test_box = abs(adj_box_l->umax[0][1][0] - box->bmax[0]) +
			               abs(adj_box_l->umax[0][1][2] - box->bmax[2]);
			    k = 0;
			    for (j = 1; j < adj_box_l->uy; j++)
			        if (test_box > fabs(adj_box_l->umax[j][1][0] - box->bmax[0]) +
				     fabs(adj_box_l->umax[j][1][2] - box->bmax[2]))
			        {
			            test_box = abs(adj_box_l->umax[j][1][0] - box->bmax[0]) +
				         abs(adj_box_l->umax[j][1][2] - box->bmax[2]);
				    k = j;
			        }
                            grid_size = min(box->bmin[1],adj_box_l->umin[k][1][1] -
                                                   comm_grid_size);
                            box->bmin[1] = max(0,grid_size);
                            box->bmax[1] = max(box->bmax[1], adj_box_l->umax[k][1][1] -
                                           comm_grid_size);
			    box->bmin[0] = min(box->bmin[0],adj_box_l->umin[k][1][0]);
			    box->bmax[0] = max(box->bmax[0],adj_box_l->umax[k][1][0]);
			    box->bmin[2] = min(box->bmin[2],adj_box_l->umin[k][1][2]);
			    box->bmax[2] = max(box->bmax[2],adj_box_l->umax[k][1][2]);
		        }
		    }
		}
		else    /* box in buffer left */
		{
		    icoord = comp_grid.U[i];
		    if (i == 0 && dir == 0)
		    {
		        if (fabs(icoord - adj_box_u->L0) < crx_tol ||
			    (fabs(icoord - adj_box_u->L0) > 2-crx_tol &&
			    fabs(icoord - adj_box_u->L0) < 2+crx_tol))
			{
			    test_box = abs(adj_box_u->lmax[0][0][1] - box->bmax[1]) +
			             abs(adj_box_u->lmax[0][0][2] - box->bmax[2]);
			    k = 0;
			    for (j = 1; j < adj_box_u->lx; j++)
			        if (test_box > fabs(adj_box_u->lmax[j][0][1] - box->bmax[1]) +
				     fabs(adj_box_u->lmax[0][0][2] - box->bmax[2]))
				{
				    test_box = abs(adj_box_u->lmax[j][0][1] - box->bmax[1]) +
				           abs(adj_box_u->lmax[0][0][2] - box->bmax[2]);
				    k = j;
				}
                             box->bmin[0] = min(box->bmin[0],adj_box_u->lmin[k][0][0] +
                                                   comm_grid_size);
                             grid_size = max(box->bmax[0], adj_box_u->lmax[k][0][0] +
                                           comm_grid_size);
                             box->bmax[0] = min(gmax[i],grid_size);
                             box->bmin[1] = min(box->bmin[1],adj_box_u->lmin[k][0][1]);
			     box->bmax[1] = max(box->bmax[1],adj_box_u->lmax[k][0][1]);
			     box->bmin[2] = min(box->bmin[2],adj_box_u->lmin[k][0][2]);
			     box->bmax[2] = max(box->bmax[2],adj_box_u->lmax[k][0][2]);
			}
		    if (fabs(icoord - adj_box_l->L0) < crx_tol ||
		        (fabs(icoord - adj_box_l->L0) > 2-crx_tol &&
			fabs(icoord - adj_box_l->L0) < 2+crx_tol))
		    {
		        test_box = abs(adj_box_l->lmax[0][0][1] - box->bmax[1]) +
			         abs(adj_box_l->lmax[0][0][2] - box->bmax[2]);
			k = 0;
			for (j = 1; j < adj_box_l->lx; j++)
			    if (test_box > fabs(adj_box_l->lmax[j][0][1] - box->bmax[1]) +
			         fabs(adj_box_l->lmax[0][0][2] - box->bmax[2]))
			    {
			        test_box = abs(adj_box_l->lmax[j][0][1] - box->bmax[1]) +
				           abs(adj_box_l->lmax[0][0][2] - box->bmax[2]);
				k = j;
			    }
                        box->bmin[0] = min(box->bmin[0],adj_box_l->lmin[k][0][0] +
                                                   comm_grid_size);
                        grid_size = max(box->bmax[0], adj_box_l->lmax[k][0][0] +
                                           comm_grid_size);
                             box->bmax[0] = min(gmax[i],grid_size);
			box->bmin[1] = min(box->bmin[1],adj_box_l->lmin[k][0][1]);
			box->bmax[1] = max(box->bmax[1],adj_box_l->lmax[k][0][1]);
			box->bmin[2] = min(box->bmin[2],adj_box_l->lmin[k][0][2]);
			box->bmax[2] = max(box->bmax[2],adj_box_l->lmax[k][0][2]);
		    }
                    }
		    if (i == 1 && dir == 1)
		    {
		        if (fabs(icoord - adj_box_u->L1) < crx_tol ||
			    (fabs(icoord - adj_box_u->L1) > 2-crx_tol &&
			    fabs(icoord - adj_box_u->L1) < 2+crx_tol))
			{
			    test_box = abs(adj_box_u->lmax[0][1][0] - box->bmax[0]) +
			               abs(adj_box_u->lmax[0][1][2] - box->bmax[2]);
			    k = 0;
			    for (j = 1; j < adj_box_u->ly; j++)
			    {
			        if (test_box > fabs(adj_box_u->lmax[j][1][0] - box->bmax[0]) +
				              fabs(adj_box_u->lmax[j][1][2] - box->bmax[2]))
				{
				    test_box = abs(adj_box_u->lmax[j][1][0] - box->bmax[0]) +
				              abs(adj_box_u->lmax[j][1][2] - box->bmax[2]);
				    k = j;
				}
			    }
			    box->bmin[1] = min(box->bmin[1],adj_box_u->lmin[k][1][1] +
                                                   comm_grid_size);
                            grid_size = max(box->bmax[1], adj_box_u->lmax[k][1][1] +
                                           comm_grid_size);
                            box->bmax[1] = min(gmax[i],grid_size);

			    box->bmin[0] = min(box->bmin[0],adj_box_u->lmin[k][1][0]);
			    box->bmax[0] = max(box->bmax[0],adj_box_u->lmax[k][1][0]);
			    box->bmin[2] = min(box->bmin[2],adj_box_u->lmin[k][1][2]);
			    box->bmax[2] = max(box->bmax[2],adj_box_u->lmax[k][1][2]);
			}

			if (fabs(icoord - adj_box_l->L1) < crx_tol ||
			    (fabs(icoord - adj_box_l->L1) > 2-crx_tol &&
			    fabs(icoord - adj_box_l->L1) < 2+crx_tol))
			{
			    test_box = abs(adj_box_l->lmax[0][1][0] - box->bmax[0]) +
			               abs(adj_box_l->lmax[0][1][2] - box->bmax[2]);
			    k = 0;
			    for (j = 1; j < adj_box_l->ly; j++)
			        if (test_box > fabs(adj_box_l->lmax[j][1][0] - box->bmax[0]) +
				              fabs(adj_box_l->lmax[j][1][2] - box->bmax[2]))
				{
				    test_box = abs(adj_box_l->lmax[j][1][0] - box->bmax[0]) +
				              abs(adj_box_l->lmax[j][1][2] - box->bmax[2]);
				    k = j;
				}
                            box->bmin[1] = min(box->bmin[1],adj_box_l->lmin[k][1][1] +
                                                   comm_grid_size);
                            grid_size = max(box->bmax[1], adj_box_l->lmax[k][1][1] +
                                           comm_grid_size);
                            box->bmax[1] = min(gmax[i],grid_size);

			    box->bmin[0] = min(box->bmin[0],adj_box_l->lmin[k][1][0]);
			    box->bmax[0] = max(box->bmax[0],adj_box_l->lmax[k][1][0]);
			    box->bmin[2] = min(box->bmin[2],adj_box_l->lmin[k][1][2]);
			    box->bmax[2] = max(box->bmax[2],adj_box_l->lmax[k][1][2]);
			}
				
		    }
		}
	    }
        }
	return NO;
}


LOCAL  void  pp_send_box(
	PP_GRID *pp_grid,
        int    *me,
	int    dir,
	int    side,
	COMM_BOX  *comm_box,
	INTERFACE *interfa)
{
        int       him[MAXD],i,j,k;
	int       myid, dst_id;
	size_t    len;
	byte      *buf, *storage = NULL;
	POINTER   info;


	if (rect_boundary_type(interfa,dir,side) != SUBDOMAIN_BOUNDARY)
	    return;
	
	myid = pp_mynode();
	dst_id = neighbor_id(him, me, dir, side, pp_grid);

	if (myid == dst_id) return;

	len = sizeof(COMM_BOX);

	scalar(&storage, len);
	buf = storage;

	for (i = 0; i < 4; i++)
	    for (k = 0; k < 2; k++)
	        for (j = 0; j < 3; j++)
	        {
	            info = (POINTER) buf;
		    ft_assign(info, &(comm_box->lmin[i][k][j]),sizeof(int));
		    buf += sizeof(int);

		    info = (POINTER) buf;
		    ft_assign(info, &(comm_box->lmax[i][k][j]),sizeof(int));
		    buf += sizeof(int);

		    info = (POINTER) buf;
		    ft_assign(info, &(comm_box->umin[i][k][j]),sizeof(int));
		    buf += sizeof(int);

		    info = (POINTER) buf;
		    ft_assign(info, &(comm_box->umax[i][k][j]),sizeof(int));
		    buf += sizeof(int);
	        }

	info = (POINTER) buf;
	ft_assign(info, &(comm_box->L0),sizeof(double));
	buf += sizeof(double);

	info = (POINTER) buf;
	ft_assign(info, &(comm_box->U0),sizeof(double));
	buf += sizeof(double);

	info = (POINTER) buf;
	ft_assign(info, &(comm_box->L1),sizeof(double));
	buf += sizeof(double);

	info = (POINTER) buf;
	ft_assign(info, &(comm_box->U1),sizeof(double));
	buf += sizeof(double);

	info = (POINTER) buf;
	ft_assign(info, &(comm_box->lx),sizeof(int));
	buf += sizeof(int);

	info = (POINTER) buf;
	ft_assign(info, &(comm_box->ly),sizeof(int));
	buf += sizeof(int);

	info = (POINTER) buf;
	ft_assign(info, &(comm_box->ux),sizeof(int));
	buf += sizeof(int);

	info = (POINTER) buf;
	ft_assign(info, &(comm_box->uy),sizeof(int));
	buf += sizeof(int);

	info = (POINTER) buf;
	ft_assign(info, &(comm_box->flag),sizeof(int));
	buf += sizeof(int);
	
	
	
	

        buf = storage + sizeof(COMM_BOX);

	pp_send(0, (POINTER)storage, len, dst_id);

	free(storage);
}  /*end pp_send_box */


LOCAL    void    pp_receive_box(
	PP_GRID    *pp_grid,
        int        *me,
	int        dir,
	int        side,
	COMM_BOX   *comm_box,
	INTERFACE *interfa)
{
        int        him[MAXD],i,j,k;
	int        myid, src_id;
	size_t     len;
	byte       *buf, *storage = NULL;
	POINTER    info;

	if (rect_boundary_type(interfa,dir,side) != SUBDOMAIN_BOUNDARY)
	    return;

	myid = pp_mynode();
	src_id = neighbor_id(him,me,dir,side,pp_grid);

	if (myid == src_id)
	    return;
        
	len = sizeof(COMM_BOX);
	scalar(&storage, len);
	pp_recv(0, src_id, (POINTER)(storage), len);

	buf = storage;
	
        
	for (i = 0; i < 4; i++)
	    for (k = 0; k < 2; k++)
	    for (j = 0; j < 3; j++)
	    {
	        info = (POINTER) buf;
		ft_assign(&(comm_box->lmin[i][k][j]), info, sizeof(int));
		buf += sizeof(int);

		info = (POINTER) buf;
		ft_assign(&(comm_box->lmax[i][k][j]), info, sizeof(int));
		buf += sizeof(int);

		info = (POINTER) buf;
		ft_assign(&(comm_box->umin[i][k][j]), info, sizeof(int));
		buf += sizeof(int);

		info = (POINTER) buf;
		ft_assign(&(comm_box->umax[i][k][j]), info, sizeof(int));
		buf += sizeof(int);
	    }

	info = (POINTER) buf;
	ft_assign(&(comm_box->L0), info, sizeof(double));
	buf += sizeof(double);

	info = (POINTER) buf;
	ft_assign(&(comm_box->U0), info, sizeof(double));
	buf += sizeof(double);

	info = (POINTER) buf;
	ft_assign(&(comm_box->L1), info, sizeof(double));
	buf += sizeof(double);

	info = (POINTER) buf;
	ft_assign(&(comm_box->U1), info, sizeof(double));
	buf += sizeof(double);

	info = (POINTER) buf;
	ft_assign(&(comm_box->lx), info, sizeof(int));
	buf += sizeof(int);

	info = (POINTER) buf;
	ft_assign(&(comm_box->ly), info, sizeof(int));
	buf += sizeof(int);

	info = (POINTER) buf;
	ft_assign(&(comm_box->ux), info, sizeof(int));
	buf += sizeof(int);

	info = (POINTER) buf;
	ft_assign(&(comm_box->uy), info, sizeof(int));
	buf += sizeof(int);

	info = (POINTER) buf;
	ft_assign(&(comm_box->flag), info, sizeof(int));
	buf += sizeof(int);
	
	
	
	

	buf = storage + sizeof(COMM_BOX);

	free(storage);

} /* end pp_receive_box */

LOCAL   void extend_boundary_side(
        INTERFACE *intfc,
        int *smin,
        int *smax,
        RECT_GRID *gr)
{
        SURFACE **s;
        TRI *t;
        POINT *p;
        int i;
        double L[MAXD],U[MAXD];

        for (i = 0; i < 3; ++i)
        {
            L[i] = gr->L[i] + smin[i]*gr->h[i];
            U[i] = gr->L[i] + smax[i]*gr->h[i];
        }
        next_point(intfc,NULL,NULL,NULL);
        for (s = intfc->surfaces; s && *s; ++s)
        {
            for (t = first_tri(*s); !at_end_of_tri_list(t,*s); t = t->next)
            {
                for (i = 0; i < 3; ++i)
                {
                    if (is_side_bdry(t,i))
                    {
                        p = Point_of_tri(t)[i];
                        if (!sorted(p))
                            check_and_extend_point(p,t,*s,L,U,gr->h);
                        p = Point_of_tri(t)[(i+1)%3];
                        if (!sorted(p))
                            check_and_extend_point(p,t,*s,L,U,gr->h);
                        check_and_extend_side(t,i,*s,L,U,gr->h);
                    }
                }
            }
        }
}       /* end extend_boundary_side */

LOCAL void check_and_extend_point(
        POINT *p,
        TRI *t,
        SURFACE *s,
        double *L,
        double *U,
        double *h)
{
        int i,num_tris;
        TRI **tris;
        int dir,nb;
        double dmin = HUGE;
        for (i = 0; i < 3; ++i)
        {
            if (Coords(p)[i] < L[i] || Coords(p)[i] > U[i])
                return;
            else
            {
                if (fabs(Coords(p)[i] - L[i]) < dmin)
                {
                    dmin = fabs(Coords(p)[i] - L[i]);
                    dir = i;    nb = 0;
                }
                if (fabs(Coords(p)[i] - U[i]) < dmin)
                {
                    dmin = fabs(Coords(p)[i] - U[i]);
                    dir = i;    nb = 1;
                }
            }
        }

        num_tris = set_tri_list_around_point(p,t,&tris,s->interface);
        if (nb == 0)  Coords(p)[dir] = L[dir] - 0.001*h[dir];
        else if (nb == 1)  Coords(p)[dir] = U[dir] + 0.001*h[dir];
        s->interface->modified = YES;
        sorted(p) = YES;
        for (i = 0; i < num_tris; ++i)
            set_normal_of_tri(tris[i]);
}       /* end check_and_extend_point */


LOCAL   void check_and_extend_side(
        TRI *tri,
        int side,
        SURFACE *surf,
        double *L,
        double *U,
        double *h)
{
        int i,j,k;
        POINT *p1,*p2,*p;
        double d1,d2,coords[MAXD];
        const double *nor;
        int l,num_tris;
        TRI **tris;

	p = NULL;
        p1 = Point_of_tri(tri)[side];
        p2 = Point_of_tri(tri)[Next_m3(side)];
        nor = Tri_normal(tri);
        for (i = 0; i < 3; ++i)
        {
            j = (i+1)%3;
            if ((Coords(p1)[i] < L[i] && Coords(p1)[j] > L[j] &&
                 Coords(p2)[i] > L[i] && Coords(p2)[j] < L[j]) ||
                (Coords(p2)[i] < L[i] && Coords(p2)[j] > L[j] &&
                 Coords(p1)[i] > L[i] && Coords(p1)[j] < L[j]))
            {
                k = Next_m3(j);
                d1 = sqr(Coords(p1)[i] - L[i]) + sqr(Coords(p1)[j] - L[j]);
                d2 = sqr(Coords(p2)[i] - L[i]) + sqr(Coords(p2)[j] - L[j]);
                if (d1 < d2) p = p1;
                else p = p2;
                coords[i] = L[i] - 0.0001*h[i];
                coords[j] = L[j] - 0.0001*h[j];
                coords[k] = Coords(p)[k] - (nor[i]*(coords[i] - Coords(p)[i])
                        + nor[j]*(coords[j] - Coords(p)[j]))/nor[k];
                for (l = 0; l < 3; ++l)
                    Coords(p)[l] = coords[l];
            }
            else if ((Coords(p1)[i] > U[i] && Coords(p1)[j] > L[j] &&
                 Coords(p2)[i] < U[i] && Coords(p2)[j] < L[j]) ||
                (Coords(p2)[i] > U[i] && Coords(p2)[j] > L[j] &&
                 Coords(p1)[i] < U[i] && Coords(p1)[j] < L[j]))
            {
                k = Next_m3(j);
                d1 = sqr(Coords(p1)[i] - U[i]) + sqr(Coords(p1)[j] - L[j]);
                d2 = sqr(Coords(p2)[i] - U[i]) + sqr(Coords(p2)[j] - L[j]);
                if (d1 < d2) p = p1;
                else p = p2;
                coords[i] = U[i] + 0.0001*h[i];
                coords[j] = L[j] - 0.0001*h[j];
                coords[k] = Coords(p)[k] - (nor[i]*(coords[i] - Coords(p)[i])
                        + nor[j]*(coords[j] - Coords(p)[j]))/nor[k];
                for (l = 0; l < 3; ++l)
                    Coords(p)[l] = coords[l];
            }
            else if ((Coords(p1)[i] < L[i] && Coords(p1)[j] < U[j] &&
                 Coords(p2)[i] > L[i] && Coords(p2)[j] > U[j]) ||
                (Coords(p2)[i] < L[i] && Coords(p2)[j] < U[j] &&
                 Coords(p1)[i] > L[i] && Coords(p1)[j] > U[j]))
            {
                k = Next_m3(j);
                d1 = sqr(Coords(p1)[i] - L[i]) + sqr(Coords(p1)[j] - U[j]);
                d2 = sqr(Coords(p2)[i] - L[i]) + sqr(Coords(p2)[j] - U[j]);
                if (d1 < d2) p = p1;
                else p = p2;
                coords[i] = L[i] - 0.0001*h[i];
                coords[j] = U[j] + 0.0001*h[j];
                coords[k] = Coords(p)[k] - (nor[i]*(coords[i] - Coords(p)[i])
                        + nor[j]*(coords[j] - Coords(p)[j]))/nor[k];
                for (l = 0; l < 3; ++l)
                    Coords(p)[l] = coords[l];
            }
            else if ((Coords(p1)[i] > U[i] && Coords(p1)[j] < U[j] &&
                 Coords(p2)[i] < U[i] && Coords(p2)[j] > U[j]) ||
                (Coords(p2)[i] > U[i] && Coords(p2)[j] < U[j] &&
                 Coords(p1)[i] < U[i] && Coords(p1)[j] > U[j]))
            {
                k = Next_m3(j);
                d1 = sqr(Coords(p1)[i] - U[i]) + sqr(Coords(p1)[j] - U[j]);
                d2 = sqr(Coords(p2)[i] - U[i]) + sqr(Coords(p2)[j] - U[j]);
                if (d1 < d2) p = p1;
                else p = p2;
                coords[i] = U[i] + 0.0001*h[i];
                coords[j] = U[j] + 0.0001*h[j];
                coords[k] = Coords(p)[k] - (nor[i]*(coords[i] - Coords(p)[i])
                        + nor[j]*(coords[j] - Coords(p)[j]))/nor[k];
                for (l = 0; l < 3; ++l)
                    Coords(p)[l] = coords[l];
            }
        }
	if (p != NULL)
	{
            num_tris = set_tri_list_around_point(p,tri,&tris,
                                surf->interface);
            for (l = 0; l < num_tris; ++l)
            {
            	set_normal_of_tri(tris[l]);
            }
	}
}       /* end check_and_extend_side */

LOCAL	boolean  set_reconstruction_boxes_prev(
	int *smin,
	int *smax,
	int **ips,
	int num_ip,
	RECT_BOX **boxes,
	INTERFACE *intfc)
{
	RECT_BOX 	Box,*box,*nbox, *btmp;
	int 		**pip,nb,itmp[3],i,j,k;
	boolean		boxes_merged, status;
	RECT_GRID	*gr = &topological_grid(intfc);
	int		tmin[3], tmax[3];
	double		dist1, dist2;

	DEBUG_ENTER(set_reconstruction_boxes)
	
	Box.prev = Box.next = NULL;
	box = &Box;
	

	/* finding the actural box by checking tangled tris. */
	for (box = Box.next; box != NULL; box = box->next)
	{
	    print_int_vector("bminp", box->bmin, 3, "\n");
	    print_int_vector("bmaxp", box->bmax, 3, "\n");
	    
	    if(max3(box->bmax[2]-box->bmin[2],
	    	    box->bmax[1]-box->bmin[1], 
	    	    box->bmax[0]-box->bmin[0]) > 7)
	    {
		printf("ERROR set_reconstruction_boxes, too many bad ips");
		clean_up(ERROR);
	    }

	    /*
	    status = tangled_tris_bound_box(tmin, tmax, box->bmin, 
	    				    box->bmax, intfc);
	    printf("#tangle status %d\n", status);

	    if(status)
	    {
		ft_assign(box->bmin, tmin, 3*INT);
		ft_assign(box->bmax, tmax, 3*INT);
	    }
	    else
	    {
		printf("#no tris tangle.\n");
		for (i = 0; i < 3; i++)
		{
		    box->bmax[i]++;
		}
	    }
	    */
	}

	if(rect_boxes_from_tangled_tris(Box.next, intfc) == 0)
	{
	    *boxes = NULL;
	    DEBUG_LEAVE(set_reconstruction_boxes)
	    return YES;
	}

	/*bad ips are in the faces of [bmin, bmax], 
	  should increase by one to make consistent surf
	  in the face of the box.
	*/
	for (box = Box.next; box != NULL; box = box->next)
	{
	    for (i = 0; i < 3; i++)
	    {
	        box->smin[i] = smin[i];
		box->smax[i] = smax[i];
	    }
	    for (i = 0; i < 3; i++)
	    {
		if (box->bmin[i] < smin[i]) box->bmin[i] = smin[i];
		if (box->bmax[i] > smax[i]) box->bmax[i] = smax[i];
	    }
	    box->grid = gr;
	}


	/* Merge overlapping boxes */
	boxes_merged = YES;
	while (boxes_merged)
	{
	    boxes_merged = NO;
	    for (box = Box.next; box != NULL; box = box->next)
	    {
		if(max3(box->bmax[2]-box->bmin[2],
			box->bmax[1]-box->bmin[1],
			box->bmax[0]-box->bmin[0]) > 10)
			{
			    printf("WARNING set_reconstruction_boxes"
			    	    "box is too large when merging.\n");
			    DEBUG_LEAVE(set_reconstruction_boxes)
			    return NO;
			}

		for (nbox = box->next; nbox != NULL; nbox = nbox->next)
		{
		    if (overlapping_boxes(box, nbox))
		    {
			/* merge nbox with box and delete nbox. */
			for (i = 0; i < 3; ++i)
			{
			    box->bmin[i] = min(box->bmin[i],nbox->bmin[i]);
			    box->bmax[i] = max(box->bmax[i],nbox->bmax[i]);
			}
			nbox->prev->next = nbox->next;
			if (nbox->next != NULL)
			    nbox->next->prev = nbox->prev;
			boxes_merged = YES;
			break;
		    }
		}   /* for(nbox= ) */
		if(boxes_merged)
		    break;
	    }   /* for(box=  ) */
	}

	*boxes = Box.next;

	i = 0;
	for (box = Box.next; box != NULL; box = box->next)
	{
	    i++;
	}

	DEBUG_LEAVE(set_reconstruction_boxes)
	return YES;
}	/* end set_reconstruction_boxes */

/* rm_bad_crxs_in_box */
LOCAL	boolean  rm_bad_crxs_in_box_prev(
	int *smin,
	int *smax,
	int **ips_in,
	int num_ip,
	INTERFACE *intfc)
{
	RECT_BOX 	*box0,*box;
	int 		i,j,k;
	RECT_GRID	*gr = &topological_grid(intfc);
	static	int	**ips = NULL;

	DEBUG_ENTER(rm_bad_crxs_in_box)
	
	if(ips == NULL)
	    stat_matrix(&ips,MAX_NUM_UNPHY_IP,3,INT);

	for(i=0; i<num_ip; i++)
	    ft_assign(ips[i], ips_in[i], 3*INT);

	make_boxes_from_ips(ips, num_ip, &box0);

	for(box=box0; box!=NULL; box=box->next)
	{
	    print_int_vector("bminp", box->bmin, 3, "\n");
	    print_int_vector("bmaxp", box->bmax, 3, "\n");
	    
	    for (i = 0; i < 3; i++)
	    {
	        box->smin[i] = smin[i];
		box->smax[i] = smax[i];
	    }
	
	    /*bad ips are in the faces of [bmin, bmax], 
	      should increase by one to make consistent surf
	      in the face of the box.
	    */
	    for (i = 0; i < 3; i++)
	    {
		box->bmin[i]--;
		box->bmax[i]++;
		if (box->bmin[i] < smin[i]) box->bmin[i] = smin[i];
		if (box->bmax[i] > smax[i]) box->bmax[i] = smax[i];
	    }
	    box->grid = gr;
	}

	merge_adjacent_boxes(box0);
	
	i = 0;
	for(box=box0; box!=NULL; box=box->next)
	{
	    i++;

	    /* fix crossings in a box */
	    remove_unphysical_crxings(box->bmin,box->bmax,gr->gmax,
			intfc,SINGLE,&num_ip,ips);
	    check_and_repair_crx(intfc,box->bmin,box->bmax);
	}

	DEBUG_LEAVE(rm_bad_crxs_in_box)
	return YES;
}	/* end rm_bad_crxs_in_box_prev */

void	make_ggrid(
	GGRID		*ggr,
	Front		*fr)
{
	RECT_GRID	Dgr;
	RECT_GRID	*pgr = &fr->pp_grid->Global_grid;
	INTERFACE	*intfc = fr->interf;
	RECT_GRID	*tgr = &topological_grid(intfc);
	RECON_BOX	*rbox;
	int		i;

	set_dual_grid(&Dgr, pgr);
	
	ft_assign(ggr->h, Dgr.h, 3*FLOAT);
	ft_assign(ggr->GL, Dgr.L, 3*FLOAT);
	ft_assign(ggr->GU, Dgr.U, 3*FLOAT);
	ft_assign(ggr->gmax, Dgr.gmax, 3*INT);

	rbox = &ggr->rbox;
	for(i=0; i<3; i++)
	{
	    rbox->bmin[i] = floor((tgr->L[i]-ggr->GL[i])/ggr->h[i] + 0.5);
	    rbox->bmax[i] = floor((tgr->U[i]-ggr->GL[i])/ggr->h[i] + 0.5);
	    rbox->fmin[i] = tgr->L[i];
	    rbox->fmax[i] = tgr->U[i];
	}
	
	print_int_vector("rb bmin", rbox->bmin, 3, "\n");
	print_int_vector("rb bmax", rbox->bmax, 3, "\n");
}


int	make_recon_boxes(
	RECON_BOX	*rbox,
	RECT_BOX	*boxes,
	GGRID		*ggr)
{
	int		i, k;
	RECT_BOX	*box;

	i = 0;
	for(box=boxes; box!=NULL; box=box->next)
	{
	  for(k=0; k<3; k++)
	  {
	    rbox[i].bmin[k] = box->bmin[k] + (ggr->rbox).bmin[k];
	    rbox[i].bmax[k] = box->bmax[k] + (ggr->rbox).bmin[k];
	    rbox[i].fmin[k] = ggr->GL[k] + rbox[i].bmin[k]*ggr->h[k];
	    rbox[i].fmax[k] = ggr->GL[k] + rbox[i].bmax[k]*ggr->h[k];
	    rbox[i].shift[0][k] = 0;
	  }
	  rbox[i].flag = 0;
	  rbox[i].proc = rbox[i].procs[0] = pp_mynode();
	  rbox[i].np = 1;
	  
	  i++;
	  if(i >= MAX_RBOX)
	  {
	    printf("ERROR communicate_boxes, too many boxes.\n");
	    clean_up(ERROR);
	  }
	}

	return i;
}

int	rbox_arrange(
	RECON_BOX	*rbox,
	int		nb)
{
	int	i, ni;

	ni = 0;
	for(i=0; i<nb; i++)
	{
	  if(i>ni)
	    rbox[ni] = rbox[i];
	  if(rbox[i].flag == 0)
	    ni++;
	}

	return ni;
}

boolean	rbox_same(
	RECON_BOX	*b1,
	RECON_BOX	*b2)
{
	int	i;

	for(i=0; i<3; i++)
	  if(b1->bmin[i] != b2->bmin[i] ||
	     b1->bmax[i] != b2->bmax[i])
	    return NO;

	return YES;
}

boolean	rbox_sect(
	RECON_BOX	*b1,
	RECON_BOX	*b2,
	boolean		bd)
{
	int	i;

	/* bd == YES mean boundary overlap is treated as sect */
	if(bd)
	{
	  for(i=0; i<3; i++)
	    if(b1->bmin[i] > b2->bmax[i] ||
	       b2->bmin[i] > b1->bmax[i])
	    return NO;
	}
	else
	{
	  for(i=0; i<3; i++)
	    if(b1->bmin[i] >= b2->bmax[i] ||
	       b2->bmin[i] >= b1->bmax[i])
	    return NO;
	}
	return YES;
}

/* box b1 is contained in box b2 */
boolean	rbox_contain(
	RECON_BOX	*b1,
	RECON_BOX	*b2)
{
	int	i;

	for(i=0; i<3; i++)
	  if(b1->bmin[i] < b2->bmin[i] ||
	     b1->bmax[i] > b2->bmax[i])
	    return NO;
  
	return YES;
}


int	rbox_add_one_proc(
	RECON_BOX	*rbox,
	int		proc,
	int		*sv)
{
	int	j, k;

	for(j=0; j<rbox->np; j++)
	  if(proc >= rbox->procs[j])
	    break;

	/* proc exist in b1 */
	if(j<rbox->np && proc == rbox->procs[j])
	  return -1;

	/* now b1->procs[j-1] > b2->procs[i] > b1->procs[j] */
	if(rbox->np >= MAX_RBOX_PROC-1)
	{
	  printf("ERROR rbox_add_one_proc, too many procs.\n");
	  clean_up(ERROR);
	}
	  
	/* proc does not exist in b1 */
	for(k=rbox->np-1; k>=j; k--)
	{
	  rbox->procs[k+1] = rbox->procs[k];
	  ft_assign(rbox->shift[k+1], rbox->shift[k], 3*INT);
	  rbox->bd_sum[k+1] = rbox->bd_sum[k];
	}
	rbox->procs[j] = proc;

	/* rbox belongs to processor proc, do not need shift */
	for(k=0; k<3; k++)
	  rbox->shift[j][k] = sv[k];

	rbox->np++;

	return j;
}

/* merge b2 to b1 */
void	rbox_merge_procs(
	RECON_BOX	*b1,
	RECON_BOX	*b2)
{
	int	i, ind;

	for(i=0; i<b2->np; i++)
	{
	  ind = rbox_add_one_proc(b1, b2->procs[i], b2->shift[i]);
	  if(ind != -1)
	    b1->bd_sum[ind] = b2->bd_sum[i];
	}
}

void	rbox_merge(
	RECON_BOX	*b1,
	RECON_BOX	*b2)
{
	int	i;

	for(i=0; i<3; i++)
	{
	  b1->bmin[i] = min(b1->bmin[i], b2->bmin[i]);
	  b1->bmax[i] = max(b1->bmax[i], b2->bmax[i]);
	  b1->fmin[i] = min(b1->fmin[i], b2->fmin[i]);
	  b1->fmax[i] = max(b1->fmax[i], b2->fmax[i]);
	}
	rbox_merge_procs(b1, b2);
	b2->flag = -1;
}

void	rbox_shift(
	RECON_BOX	*rbox,
	int		st,
	int		ed,
	RECT_GRID	*gr,
	int		dir,
	int		sign)
{
	RECON_BOX	*b;
	double		T;
	int		Tc, i, j;

	T = (gr->U[dir] - gr->L[dir])*sign;
	Tc = gr->gmax[dir]*sign;
	
	for(i=st; i<st+ed; i++)
	{
	  b = &rbox[i];
	  
	  b->bmin[dir] += Tc;
	  b->bmax[dir] += Tc;
	  b->fmin[dir] += T;
	  b->fmax[dir] += T;
	  
	  for(j=0; j<b->np; j++)
	    b->shift[j][dir] += sign;
	}
}

/* send rbox from 0 to nbox and recv from nbox to bst */
int	rbox_communication_in_dir(
	RECON_BOX	*rbox,
	int		nbox,
	int		dir,
	Front		*fr)
{
	RECT_GRID	*ggrid;
	INTERFACE    	*intfc = fr->interf;
	PP_GRID	     	*pp_grid = fr->pp_grid;
	int	     	me[3], him[3];
	int	     	dst_id;
	int		*G, bst, nbr;
	int		i, j, k, jp;

	G = pp_grid->gmax;
	find_Cartesian_coordinates(pp_mynode(),pp_grid,me);
	ggrid = &pp_grid->Global_grid;

	bst = nbox;
	
	pp_gsync();
	
	/* send and recv boxes */
	for(j=0; j<2; j++)
	{
	  if(rect_boundary_type(intfc,dir,j) == SUBDOMAIN_BOUNDARY)
	  {
	    dst_id = neighbor_id(him, me, dir, j, pp_grid);
	      
	    pp_send(2*j, &nbox, INT, dst_id);
	    if(nbox != 0)
	      pp_send(2*j+1, rbox, nbox*sizeof(RECON_BOX), dst_id);
	  }

	  jp = (j+1)%2;
	  if(rect_boundary_type(intfc,dir,jp) == SUBDOMAIN_BOUNDARY)
	  {
	    dst_id = neighbor_id(him, me, dir, jp, pp_grid);
	      
	    pp_recv(2*j, dst_id, &nbr, INT);
	    if(nbr == 0)
	      continue;
	    if(bst + nbr >= MAX_RBOX)
	    {
	      printf("ERROR recon_box_communication, too many boxes\n");
	      clean_up(ERROR);
	    }
	    pp_recv(2*j+1, dst_id, &rbox[bst], nbr*sizeof(RECON_BOX));
	      
	    /* shift for periodic case */
	    if(me[dir] == 0 && jp == 0)
	      rbox_shift(rbox, bst, nbr, ggrid, dir, -1);
	    else if(me[dir] == G[dir]-1 && jp == 1)
	      rbox_shift(rbox, bst, nbr, ggrid, dir, 1);
	      
	    bst += nbr;
	  }
	}

	return bst;
}

/*WARNING proc send to itself is not allowed 
  after calling this function, each proc knows all its intersected boxes
*/
int	recon_box_communication(
	RECON_BOX	*rbox,
	int		nbox,
	GGRID		*ggr,
	Front		*fr)
{
	int		bst;
	int		i, k, n;
	static int	zv[] = {0,0,0};

	DEBUG_ENTER(recon_box_communication)

	for(i=0; i<3; i++)
	{
	  bst = rbox_communication_in_dir(rbox, nbox, i, fr);

	  /* add my proc number to the procs list */
	  for(k=nbox; k<bst; k++)
	    rbox_add_one_proc(&rbox[k], pp_mynode(), zv);
  
	  /* remove outside boxes */
	  for(k=0; k<bst; k++)
	    if(rbox_sect(&rbox[k], &ggr->rbox, YES) == NO)
	      rbox[k].flag = -1;

	  /* merge same boxes */
	  for(k=0; k<bst; k++)
	  {
	    if(rbox[k].flag == -1)
	      continue;
	    for(n=k+1; n<bst; n++)
	    {
	      if(rbox[n].flag == -1)
	        continue;
	      if(rbox_same(&rbox[k],&rbox[n]))
	        rbox_merge(&rbox[k], &rbox[n]);
	    }
	  }
	  nbox = rbox_arrange(rbox, bst);
	}
	
	DEBUG_LEAVE(recon_box_communication)
	return nbox;
}

void	rbox_copy(
	RECON_BOX	*b1,
	RECON_BOX	*b2)
{
	int	i, np;

	ft_assign(b1->bmin, b2->bmin, 3*INT);
	ft_assign(b1->bmax, b2->bmax, 3*INT);
	ft_assign(b1->fmin, b2->fmin, 3*FLOAT);
	ft_assign(b1->fmax, b2->fmax, 3*FLOAT);
	
	np = b2->np;
	b1->np = np;
	ft_assign(b1->procs, b2->procs, np*INT);
	for(i=0; i<np; i++)
	  ft_assign(b1->shift[i], b2->shift[i], 3*INT);
	b1->flag = b2->flag;
	b1->number = b2->flag;
}
/* when calling this function, we assume b1 is in proc pp_mynode() */
void	rbox_buffer_box(
	RECON_BOX	*b,
	RECON_BOX	*b1,
	GGRID		*ggr)
{
	int		k;
	RECON_BOX	*rbox;
	double		tol;

	tol = 1e-2;
	rbox = &ggr->rbox;

	rbox_copy(b, b1);
	for(k=0; k<3; k++)
	{
	    b->bmin[k] -= RECON_BUFFER;
	    b->bd_flag[k][0] = NO;
	    
	    /* proc is in the left of the domain */
	    if(rbox->fmin[k] <= ggr->GL[k] + ggr->h[k]*tol)
	    {
	      if(b->bmin[k] <= rbox->bmin[k])
	      {
	        b->bd_flag[k][0] = YES;
	        b->bmin[k] = rbox->bmin[k];
	      }
	    }
	    
	    b->bmax[k] += RECON_BUFFER;
	    b->bd_flag[k][1] = NO;
	    
	    /* proc is in the right of the domain */
	    if(rbox->fmax[k] >= ggr->GU[k] - ggr->h[k]*tol)
	    {
	      if(b->bmax[k] >= rbox->bmax[k])
	      {
	        b->bd_flag[k][1] = YES;
	        b->bmax[k] = rbox->bmax[k];
	      }
	    }

	    b->fmin[k] = ggr->GL[k] + b->bmin[k]*ggr->h[k];
	    b->fmax[k] = ggr->GL[k] + b->bmax[k]*ggr->h[k];
	    b->shift[0][k] = 0;
	}
	b->flag = 0;
	b->proc = b->procs[0] = pp_mynode();
	b->np = 1;
}

boolean	rbox_communication_interface(RECON_BOX*,int,RECON_BOX*,GGRID*,Front*);
void	tecplot_rboxes(char*, GGRID*,RECON_BOX*,int,boolean);

boolean	rbox_one_box_communication(
	RECON_BOX	*b1,
	GGRID		*ggr,
	Front		*fr,
	RECON_BOX	*rboxin,
	int		nbin)
{
	RECON_BOX	rbox[MAX_RBOX], *b;
	INTERFACE	*intfc;
	int		i,j, nbox;
	boolean		recon;

	printf("#rbox_one_box_communication\n");

	/* copy b1 to rbox[0] if exists */
	if(b1 == NULL)
	  nbox = 0;
	else
	{
	  /* enlarge b1 by 2 cells */
	  b = &rbox[0];
	  rbox_buffer_box(b, b1, ggr);
	  nbox = 1;
	}
	
	for(i=0; i<3; i++)
	  nbox = rbox_communication_in_dir(rbox, nbox, i, fr);

	recon = NO;
	if(b1 != NULL)
	{
	  /*check overlap, if overlapping, activate the one with 
	    small proc number. */
	  for(i=1; i<nbox; i++)
	    if(rbox_sect(b, &rbox[i], YES) && b->proc < rbox[i].proc)
	      break;

	  /* activate b */
	  if(i == nbox)
	  {
	    recon = YES;
	    nbox = 1;
	  }
	  else
	    nbox = 0;
	}
	else
	  nbox = 0;

	/* comm all boxes */
	nbox = recon_box_communication(rbox,nbox,ggr,fr);

	/* must be called TWICE so each box knows which procs it belongs to */
	nbox = recon_box_communication(rbox,nbox,ggr,fr);

	tecplot_rboxes("rbox", ggr, rbox, nbox, YES);
	
	add_time_start(12);
	if(recon)
	  rbox_communication_interface(rbox, nbox, &rbox[0], ggr, fr);
	else
	  rbox_communication_interface(rbox, nbox, NULL, ggr, fr);
	add_time_end(12);

	/* if a box in rboxin is used, deact it */
	for(i=0; i<nbin; i++)
	{
	  b = &rboxin[i];
	  for(j=0; j<nbox; j++)
	  {
	    if(rbox_contain(b, &rbox[j]))
	      b->flag = -1;
	  }
	}

	return YES;
}

void	tecplot_rboxes(
	char		*bname,
	GGRID		*ggr,
	RECON_BOX	*rbox,
	int		nbox,
	boolean		info)
{
	int	i, j;
	char	fname[128];
	FILE	*fp;
	
	printf("nbox = %d\n", nbox);
	for(i=0; i<nbox; i++)
	{
	  printf("%s %2d number %2d proc %5d  procs (%2d):  ", 
	        bname, i, rbox[i].number, rbox[i].proc, rbox[i].np);
	  
	  for(j=0; j<rbox[i].np; j++)
	    printf("%5d(%2d%2d%2d) %3d", rbox[i].procs[j],
	          rbox[i].shift[j][0], rbox[i].shift[j][1], rbox[i].shift[j][2],
		  rbox[i].bd_sum[j]);
	  printf("\n");
	  
	  print_int_vector("bmin", rbox[i].bmin, 3, "\n");
	  print_int_vector("bmax", rbox[i].bmax, 3, "\n");
	}
	
	if(info)
	  return;

	if(nbox == 0)
	  return;

	sprintf(fname, "%s_%07d.plt", bname, pp_mynode());
	fp = fopen(fname, "w");
	if(fp == NULL)
	{
	  printf("ERROR tecplot_rboxes, can not open %s\n", fname);
	  clean_up(ERROR);
	}
	printf("box file name %s nbox = %d\n", fname, nbox);

	fprintf(fp,"TITLE = \"tecplot box\"\n"
	   	"VARIABLES = \"x\", \"y\", \"z\"\n");

	tecplot_box(NULL, fp, ggr->rbox.fmin, ggr->rbox.fmax);
	for(i=0; i<nbox; i++)
	  tecplot_box(NULL, fp, rbox[i].fmin, rbox[i].fmax);
	
	fclose(fp);
}

EXPORT	INTERFACE  *bboxes_make_patch_intfc(RECON_BOX*,INTERFACE*,GGRID*);

boolean	communicate_boxes(RECT_BOX*,Front*);

boolean	communicate_boxes(
	RECT_BOX	*boxes,
	Front		*fr)
{
	GGRID		ggr;
	RECT_BOX	*box;
	RECON_BOX	rbox[MAX_RBOX];
	INTERFACE	*intfc;
	int		nbox, max_nbox;
	int		i;

	add_time_clear(11);
	add_time_clear(12);
	add_time_clear(13);
	add_time_start(11);

	make_ggrid(&ggr, fr);
	nbox = make_recon_boxes(rbox, boxes, &ggr);

	/*recon_box_communication must be called TWICE to make sure each rbox
	  knows which procs it belongs to
	  1. send rbox to other procs, after this call, each rbox knows 
	  which proc it comes from, but does not know if it also belongs 
	  to other procs
	*/
	nbox = recon_box_communication(rbox, nbox, &ggr, fr);
	
	/*2. each rbox will know which procs it belongs to */
	nbox = recon_box_communication(rbox, nbox, &ggr, fr);

	tecplot_rboxes("rbox", &ggr, rbox, nbox, NO);
	
	printf("#rbox one comm bf nbox = %d\n", nbox);
	
	clean_up(0);

	max_nbox = 1;
	while(max_nbox > 0)
	{
	for(i=0; i<nbox; i++)
	{
	  if(rbox_contain(&rbox[i],&ggr.rbox))
	    break;
	}
	
	if(i == nbox)
	  rbox_one_box_communication(NULL, &ggr, fr, rbox, nbox);
	else
	  rbox_one_box_communication(&rbox[i], &ggr, fr, rbox, nbox);

	nbox = rbox_arrange(rbox, nbox);
	
	max_nbox = nbox;
	pp_global_imax(&max_nbox, 1);
	printf("#rbox one comm max_nbox = %d nbox = %d\n", max_nbox, nbox);

	tecplot_rboxes("rbox", &ggr, rbox, nbox, YES);
	}

	add_time_end(11);
	printf("#tol rbox_comm_intfc copy_intfc \n");
	printf("#comm boxes time %15.8e  %15.8e  %15.8e\n", 
	    add_time_get(11), add_time_get(12), add_time_get(13));

	return YES;
}

int	rbox_set_buffer(
	RECON_BOX	*rbox,
	GGRID		*ggr,
	int		buf,
	boolean		bdry)
{
	RECON_BOX	*gbox;
	int		i, tmp;

	gbox = &ggr->rbox;

	for(i=0; i<3; i++)
	{
	  tmp = rbox->bmin[i];
	  rbox->bd_flag[i][0] = bdry ? tmp - max(gbox->bmin[i], tmp - buf) : buf;
	  
	  tmp = rbox->bmax[i];
	  rbox->bd_flag[i][1] = bdry ? min(gbox->bmax[i], tmp + buf) - tmp : buf;
	}
}

int	rbox_global_number(
        RECON_BOX       *rbox,
        int             nbox)
{
        int     *procs;
        int     nproc, myid, i, num;

        nproc = pp_numnodes();
        myid = pp_mynode();

        uni_array(&procs, nproc, INT);
        
	procs[myid] = nbox;
        pp_global_isum(procs, nproc);

        num = 0;
        for(i=0; i<myid; i++)
          num += procs[i];

        for(i=0; i<nbox; i++)
        {
          rbox[i].number = num;
          num++;
        }

        free(procs);
}

/*sign = -1:  decrease rbox
  sign = 1 :  increase rbox
*/
void	rbox_set_buffer_box(
	RECON_BOX	*rbox,
	GGRID		*ggr,
	int		sign)
{
	int	i;
	
	for(i=0; i<3; i++)
	{
	  rbox->bmin[i] -= sign*rbox->bd_flag[i][0];
	  rbox->bmax[i] += sign*rbox->bd_flag[i][1];
	  rbox->fmin[i] = ggr->GL[i] + rbox->bmin[i]*ggr->h[i];
	  rbox->fmax[i] = ggr->GL[i] + rbox->bmax[i]*ggr->h[i];
	}
}

/* after the function rbox are contained in pp_mynode() */
int	rbox_contain_communication(
	RECON_BOX	*rbox,
	int		nbox,
	GGRID		*ggr,
	Front		*fr)
{
	int		bst;
	int		i, k, n;
	static int	zv[] = {0,0,0};

	DEBUG_ENTER(rbox_contain_communication)

	for(i=0; i<3; i++)
	{
	  bst = rbox_communication_in_dir(rbox, nbox, i, fr);

	  /* remove outside boxes */
	  for(k=0; k<bst; k++)
	    if(rbox_contain(&rbox[k], &ggr->rbox) == NO)
	      rbox[k].flag = -1;

	  /* merge same boxes */
	  for(k=0; k<bst; k++)
	  {
	    if(rbox[k].flag == -1)
	      continue;
	    for(n=k+1; n<bst; n++)
	    {
	      if(rbox[n].flag == -1)
	        continue;
	      if(rbox_same(&rbox[k],&rbox[n]))
	        rbox_merge(&rbox[k], &rbox[n]);
	    }
	  }
	  nbox = rbox_arrange(rbox, bst);
	}
	
	DEBUG_LEAVE(rbox_contain_communication)
	return nbox;
}

#define	dboxindex(pd,i,j)	(pd)->dtabstore[(pd)->dtab[i][1]+(j)]

void	make_dtable(
	DTABLE		*pd,
	RECON_BOX	*rbox,
	int		nbox,
	RECON_BOX	*dbox,
	int		ndbox)
{
	int	i, j, ds, ds0;

	pd->nbox = nbox;
	ds = 0;
	for(i=0; i<nbox; i++)
	{
	  ds0 = ds;
	  for(j=0; j<ndbox; j++)
	  {
	    if(rbox == dbox && i == j)
	      continue;

	    if(rbox_sect(&rbox[i], &dbox[j], YES) && 
	      ((dbox[j].proc != rbox[i].proc && dbox[j].proc > rbox[i].proc) ||
	       (dbox[j].proc == rbox[i].proc && dbox[j].number < rbox[i].number)))
	    {
	      pd->dtabstore[ds] = j;
	      ds++;
	      if(ds >= MAX_DSTORE)
	      {
	        printf("ERROR make_dtable, ds is too large.\n");
		clean_up(ERROR);
	      }
	    }
	  }
	  pd->dtab[i][0] = ds-ds0;
	  if(ds-ds0 != 0)
	    pd->dtab[i][1] = ds0;
	}

}


void	print_dtable(
	char		*str,
	DTABLE		*pd,
	RECON_BOX	*rbox,
	RECON_BOX	*dbox)
{
	int	i, j, nd, ind;

	printf("#dtable %s\n", str);
	for(i=0; i<pd->nbox; i++)
	{
	  nd = pd->dtab[i][0];
	  printf("#dtable  %3d %4d %4d | ", i, rbox[i].number, nd);
	  for(j=0; j<nd; j++)
	  {
	    ind = dboxindex(pd,i,j);
	    printf("%3d %3d %3d | ", ind, dbox[ind].number, dbox[ind].proc);
	  }
	  printf("\n");
	}
}

boolean	check_dtable(
	DTABLE		*pd,
	int		num,
	RECON_BOX	*dbox)
{
	int	i, ind;

	for(i=0; i<pd->dtab[num][0]; i++)
	{
	  ind = dboxindex(pd,num,i);
	  if(dbox[ind].flag != -1)
	    return NO;
	}
	return YES;
}

void	set_tag_shift(int);
int	rbox_proc_index(int,RECON_BOX*);
void	rbox_shift_interface(INTERFACE*,int*,RECT_GRID*);

INTERFACE	*rbox_wait_recv(
	int		*bnum,
	DTABLE		*ddtab,
	RECON_BOX	*rbox,
	int		nbox,
	RECT_GRID	*gr)
{
	boolean		found;
	int		i,ind;
	RECON_BOX	*b;
	INTERFACE	*intfc;
	int		cnt;

	cnt = 0;
	found = NO;
	while(!found)
	{
	  for(i=0; i<nbox; i++)
	  {
	    b = &rbox[i];
	    if(b->flag == -1)
	      continue;
	    if(pp_iprobe(b->proc, TABFACTOR*b->number + TABLE_ID) &&
	       check_dtable(ddtab, i, rbox))
	    {
	      add_time_start(328);
	      printf("#recv box %d from proc %d number %d bf\n", i, b->proc, b->number);
	      printf("#cnt = %d\n", cnt);
	      fflush(NULL);

	      set_tag_shift(TABFACTOR*b->number);
	      intfc = receive_interface(b->proc);
	      
	      printf("#recv box %d from proc af\n", i);
	      fflush(NULL);
	      
	      ind = rbox_proc_index(b->proc,b);
	      rbox_shift_interface(intfc,b->shift[ind],gr);
	      found = YES;
	      add_time_end(328);
	      break;
	    }
	  }
	  cnt++;
	}

	*bnum = i;
	return intfc;
}


boolean	rbox_active(
	RECON_BOX	*rbox,
	int		nbox)
{
	int	i;

	for(i=0; i<nbox; i++)
	  if(rbox[i].flag != -1)
	    return YES;
	return NO;
}


void	rbox_set_rect_box(
	RECT_BOX	*box,
	RECON_BOX	*rbox,
	GGRID		*ggr,
	INTERFACE	*intfc)
{

	RECT_GRID	*gr = &topological_grid(intfc);
	int		i;
	RECON_BOX	*gbox;

	gbox = &ggr->rbox;
	
	for(i=0; i<3; i++)
	{
	  box->bmin[i] = rbox->bmin[i] + rbox->bd_flag[i][0] - gbox->bmin[i];
	  box->bmax[i] = rbox->bmax[i] - rbox->bd_flag[i][1] - gbox->bmin[i];
	  box->smin[i] = 0;
	  box->smax[i] = gr->gmax[i];
	  
	  printf("%d %d   %d %d\n", box->bmin[i], box->bmax[i], 
	                            rbox->bd_flag[i][0], rbox->bd_flag[i][1]);
	}
	box->grid = gr;
	box->prev = box->next = NULL;
}

EXPORT	boolean  skip_tag_tri(TRI*);

EXPORT	boolean  skip_tag_tri(
	TRI	*tri)
{
	int	i;

	for(i=0; i<3; i++)
	    if(sorted(Point_of_tri(tri)[i]))
	        return YES;
	return NO;
}

EXPORT	boolean	bboxes_tri_outside(BBOX*,TRI*);

void	rbox_tag_recon_tris(
	INTERFACE	*intfc,
	BBOX		*bbox)
{
	SURFACE		**s;
	TRI		*tri;
	POINT		**p;
	int		i;

	next_point(intfc, NULL, NULL, NULL);
	for(s = intfc->surfaces; s && *s; ++s)
	{
	    for (tri=first_tri(*s); !at_end_of_tri_list(tri,*s); tri=tri->next)
	    {
		if(bboxes_tri_outside(bbox, tri))
		{
		  Tri_index(tri) = -1;
		  p = Point_of_tri(tri);
		  for(i=0; i<3; i++)
		    sorted(p[i]) = YES;
		}
		else
		  Tri_index(tri) = 0;
	    }
	}
}

void	rbox_copy_tris_index(
	INTERFACE	*intfc)
{
	SURFACE		**s;
	TRI		*tri;
	POINT		**p;
	int		i;

	next_point(intfc, NULL, NULL, NULL);
	for(s = intfc->surfaces; s && *s; ++s)
	{
	    for (tri=first_tri(*s); !at_end_of_tri_list(tri,*s); tri=tri->next)
	    {
	      Tri_order(tri) = Tri_index(tri);
	    }
	}
}

int	get_tag_shift();
void	assign_point_index(INTERFACE*,int);
void 	bbox_set(BBOX*,RECON_BOX*,GGRID*,double);
EXPORT	INTERFACE  *bboxes_intfc_sect(BBOX*,INTERFACE*,TRI_POS);
EXPORT void bboxes_intfc_cut(BBOX*,INTERFACE*,TRI_POS);
boolean	merge_buffer_interface(INTERFACE*,INTERFACE*,int);
void	tecplot_merge_surfaces(char*,INTERFACE*,INTERFACE**,int,int);
void	cut_intfc_in_grid(INTERFACE*,RECT_GRID*);

boolean	recon_communication_intfc(
	Front		*fr,
	GGRID		*ggr,
	RECON_BOX	*rbox,
	int		nbox,
	RECON_BOX	*dbox,
	int		ndbox)
{
	INTERFACE	*intfc, *send_intfc, *recv_intfc;
	RECT_BOX	box;
	BBOX		bbox;
	RECON_BOX	*b;
	DTABLE		dtab, rtab, ddtab;
	double		tol;
	boolean		found;
	int		i,j,ind,proc,tg_shf;
	RECT_GRID	*gr = fr->rect_grid;

	DEBUG_ENTER(recon_communication_intfc)

	add_time_clear(320);
	add_time_clear(321);
	add_time_clear(322);
	add_time_clear(323);
	add_time_clear(324);
	add_time_clear(325);
	add_time_clear(326);
	add_time_clear(327);
	add_time_clear(328);
	
	add_time_clear(331);
	add_time_clear(332);
	add_time_clear(333);
	add_time_clear(334);
	add_time_clear(335);
	
	add_time_clear(350);
	add_time_clear(351);

	tg_shf = get_tag_shift();

	tol = 1.0e-4;

	intfc = fr->interf;
	
	make_dtable(&dtab, rbox, nbox, dbox, ndbox);
	make_dtable(&rtab, rbox, nbox, rbox, nbox);
	make_dtable(&ddtab, dbox, ndbox, dbox, ndbox);

	print_dtable("dtab", &dtab, rbox, dbox);
	print_dtable("rtab", &rtab, rbox, rbox);
	print_dtable("ddtab", &ddtab, dbox, dbox);

	assign_point_index(intfc, pp_mynode());
	strip_subdomain_bdry_curves(intfc);
	
	while(rbox_active(rbox,nbox) || rbox_active(dbox,ndbox))
	{
	  add_time_start(320);
	  if(rbox_active(rbox,nbox))
	  {
	    while(1)
	    {
	      /* find a box which can recon */
	      found = NO;
	      for(i=0; i<nbox; i++)
	      {
	        if(rbox[i].flag != -1 && 
		   check_dtable(&dtab, i, dbox) && 
		   check_dtable(&rtab, i, rbox))
		{
		  found = YES;
		  break;
		}
	      }

	      if(!found)
	        break;

	      b = &rbox[i];
	      b->flag = -1;

	      install_subdomain_bdry_curves(intfc);
	        
	      printf("#gb recon %d  number = %d\n", i, b->number);

	      bbox_set(&bbox, b, ggr, tol);
	      rbox_tag_recon_tris(intfc, &bbox);

	      rbox_set_rect_box(&box, b, ggr, intfc);
	      
	      if(fr->step == 721)
	        add_to_debug("box_intfc");
	      grid_based_box_untangle(intfc, &box, YES);
	        
	      strip_subdomain_bdry_curves(intfc);
	      
	      rbox_copy_tris_index(intfc);

	      send_intfc = bboxes_intfc_sect(&bbox, intfc, TRIFLAGA);
	  
	      for(j=b->np-1; j>=0; j--)
	      {
		proc = b->procs[j];
		if(proc == pp_mynode())
		  continue;

		printf("#send intfc to proc %d bf\n", proc);
		fflush(NULL);
	          
		set_tag_shift(TABFACTOR*b->number);
		send_interface(send_intfc, proc);
	          
		printf("#send intfc to proc %d \n", proc);
	        fflush(NULL);
	      }
	        
	      delete_interface(send_intfc);
	    } /* while(1) */
	  } /* rbox active */
	  add_time_end(320);
	  
	  add_time_start(321);
	  if(rbox_active(dbox,ndbox))
	  {
	    add_time_start(322);
	    recv_intfc = rbox_wait_recv(&ind, &ddtab, dbox, ndbox, gr);
	    dbox[ind].flag = -1;
	
	    cut_intfc_in_grid(recv_intfc, computational_grid(intfc));
	    add_time_end(322);
	
    
	    add_time_start(323);
	    bbox_set(&bbox, &dbox[ind], ggr, tol);
	    bboxes_intfc_cut(&bbox, intfc, TRIINSIDE);
	    add_time_end(323);
	    

	    printf("#merge recv_intfc %d number = %d\n", ind, dbox[ind].number);
	    
	    if(dbox[ind].number == 39)
	    {
	      INTERFACE	*intfcs[2];
	      intfcs[0] = recv_intfc;
	      tecplot_merge_surfaces("bpintfc", intfc, intfcs, 1, pp_mynode());
	    }

	    add_time_start(324);
	    merge_buffer_interface(intfc, recv_intfc, dbox[ind].proc);
	    add_time_end(324);

	    if(NO)
	    {
	      null_sides_are_consistent();
	      check_print_intfc("After merge_buffer", "mgbuff", 
	            's', intfc, fr->step, fr->step, NO);
	    }

	    delete_interface(recv_intfc);
	  }
	  add_time_end(321);
	}
	
	install_subdomain_bdry_curves(intfc);

	set_tag_shift(tg_shf);
	
	printf("#                 recon  merge_tot  recv  cut  merge recv_intfc\n");
	printf("#recon comm intfc %15.8e  %15.8e  %15.8e  %15.8e  %15.8e  %15.8e\n",
	      add_time_get(320), add_time_get(321), add_time_get(322),
	      add_time_get(323), add_time_get(324), add_time_get(328));
	
	printf("#                 pp_recv   recon_pts\n");
	printf("#merge intfc %15.8e  %15.8e  \n",
	      add_time_get(350), add_time_get(351));
	
	printf("#                 point  tris  edge\n");
	printf("#merge intfc %15.8e  %15.8e  %15.8e  \n",
	      add_time_get(325), add_time_get(326), add_time_get(327));
	
	printf("#                 count  sort  compare  rotate  replace\n");
	printf("#merge intfc %15.8e  %15.8e  %15.8e  %15.8e %15.8e\n",
	      add_time_get(331), add_time_get(332), add_time_get(333),
	      add_time_get(334), add_time_get(335));

	DEBUG_LEAVE(recon_communication_intfc)
	return YES;
}

int	rbox_contain_proc(
	RECON_BOX	*rbox)
{
	int	bd_max, i, ind;

	bd_max = -1;
	for(i=0; i<rbox->np; i++)
	{
	  if(rbox->bd_sum[i] > bd_max)
	  {
	    ind = i;
	    bd_max = rbox->bd_sum[i];
	  }
	}

	return ind;
}

void	set_rbox_param(
	RECON_BOX	*rbox,
	GGRID		*ggr)
{
	int	j;

	rbox_set_buffer(rbox, ggr, 1, YES);
	rbox->proc = rbox->procs[0] = pp_mynode();
	rbox->np = 1;
	rbox->flag = 0;
	rbox->bd_sum[0] = 0;
	for(j=0; j<3; j++)
	{
	  rbox->bd_sum[0] += rbox->bd_flag[j][0] + rbox->bd_flag[j][1];
	  rbox->shift[0][j] = 0;
	}
}

int	global_set_rbox(
	RECON_BOX	*rbox,
	int		nbox,
	GGRID		*ggr,
	Front		*fr)
{
	int	i, j, ind;

	DEBUG_ENTER(global_set_rbox)

	/*1. communicate boxes and delete contained boxes
	  collect box info from other procs
	*/
	nbox = recon_box_communication(rbox, nbox, ggr, fr);
	nbox = recon_box_communication(rbox, nbox, ggr, fr);
	
	printf("#nbox 1 = %d\n", nbox);
	fflush(NULL);

	/* delete contained boxes */
	for(i=0; i<nbox; i++)
	  for(j=0; j<nbox; j++)
	    if(j != i && rbox_contain(&rbox[i], &rbox[j]))
	      rbox[i].flag = -1;
	nbox = rbox_arrange(rbox, nbox);
	
	/* find all boxes in this proc from positions */
	for(i=0; i<nbox; i++)
	  if(!rbox_contain(&rbox[i], &ggr->rbox))
	    rbox[i].flag = -1;
	nbox = rbox_arrange(rbox, nbox);

	/* 2. determine which box belongs to which proc */
	for(i=0; i<nbox; i++)
	  set_rbox_param(&rbox[i], ggr);

	printf("#nbox 2 = %d\n", nbox);
	fflush(NULL);

	nbox = rbox_contain_communication(rbox, nbox, ggr, fr);
	nbox = rbox_contain_communication(rbox, nbox, ggr, fr);

	/* if a box is contained in many proc. */
	for(i=0; i<nbox; i++)
	{
	  ind = rbox_contain_proc(&rbox[i]);
	  if(rbox[i].procs[ind] == pp_mynode())
	    rbox[i].proc = pp_mynode();
	  else
	    rbox[i].flag = -1;
	}
	nbox = rbox_arrange(rbox, nbox);
	
	printf("#nbox 3 = %d\n", nbox);
	fflush(NULL);

	return nbox;
	DEBUG_LEAVE(global_set_rbox)
}

boolean	merge_rbox(
	RECON_BOX	*rbox,
	int		nbox)
{
	boolean	merged, found;
	int	i, j;

	found = NO;
	while(1)
	{
	  merged = NO;
	  for(i=0; i<nbox; i++)
	  {
	    for(j=i+1; j<nbox; j++)
	    {
	      if(rbox_sect(&rbox[i], &rbox[j], YES))
	      {
	        rbox_merge(&rbox[i], &rbox[j]);
		merged = YES;
	        found = YES;
		rbox[j].flag = -1;
	      }
	    }
	    if(merged)
	      break;
	  }
	  nbox = rbox_arrange(rbox, nbox);
	  if(!merged)
	    break;
	}

	return found;
}

boolean ng_form_patch_subintfc_via_cut3d(Front *);

boolean	rbox_communication_boxes(RECT_BOX*,Front*);

boolean	rbox_communication_boxes(
	RECT_BOX	*boxes,
	Front		*fr)
{
	GGRID		ggr;
	RECT_BOX	*box;
	RECON_BOX	rbox[MAX_RBOX], dbox[MAX_RBOX];
	INTERFACE	*intfc;
	int		nbox, ndbox, max_nbox, max_ndbox;
	int		i,j;
	boolean		status, merged;

	DEBUG_ENTER(rbox_communication_boxes)

	pp_gsync();
	
	add_time_clear(311);
	add_time_clear(312);
	add_time_clear(313);
	add_time_clear(314);
	
	add_time_start(311);
	add_time_start(312);
	
	make_ggrid(&ggr, fr);
	
	nbox = make_recon_boxes(rbox, boxes, &ggr);
	
	/* box comm */
	j = 0;
	status = YES;
	while(status)
	{
	  for(i=0; i<nbox; i++)
	    set_rbox_param(&rbox[i], &ggr);
	  nbox = global_set_rbox(rbox, nbox, &ggr, fr);
	  merged = merge_rbox(rbox, nbox);
	  status = pp_max_status(merged);
	  printf("#merge j = %d status = %d merged = %d\n", j, status, merged);
	  j++;
	}
	
	for(i=0; i<nbox; i++)
	  rbox_set_buffer_box(&rbox[i], &ggr, 1);

	/* global index of boxes */
	rbox_global_number(rbox, nbox);
	
	tecplot_rboxes("rbox", &ggr, rbox, nbox, YES);
	
	/* 3. how many procs a box belongs to and get dbox */
	nbox = recon_box_communication(rbox, nbox, &ggr, fr);
	nbox = recon_box_communication(rbox, nbox, &ggr, fr);

	ndbox = 0;
	for(i=0; i<nbox; i++)
	  if(rbox[i].proc != pp_mynode())
	  {
	    dbox[ndbox] = rbox[i];
	    ndbox++;
	    rbox[i].flag = -1;
	  }
	nbox = rbox_arrange(rbox, nbox);
	
	printf("#nbox 3 = %d ndbox = %d\n", nbox, ndbox);
	
	tecplot_rboxes("rboxa", &ggr, rbox, nbox, NO);
	tecplot_rboxes("dbox", &ggr, dbox, ndbox, NO);
	fflush(NULL);
	
	/* print out box info */
	max_nbox = nbox;
	pp_global_imax(&max_nbox, 1);
	if(max_nbox == nbox)
	  printf("#max nbox = %d is in proc %d\n", nbox, pp_mynode());
	
	max_ndbox = ndbox;
	pp_global_imax(&max_ndbox, 1);
	if(max_ndbox == ndbox)
	  printf("#max ndbox = %d is in proc %d\n", ndbox, pp_mynode());

	add_time_end(312);
	
	add_time_start(313);
	
	/* 4. box reconstruction */
	status = recon_communication_intfc(fr,&ggr,rbox,nbox,dbox,ndbox);

	status = pp_min_status(status);
	
	add_time_end(313);
	
	add_time_start(314);
	
	set_current_interface(fr->interf);
	reset_normal_on_intfc(fr->interf);
	install_subdomain_bdry_curves(fr->interf);
	
	add_time_end(314);
	
	add_time_end(311);
	
	printf("#comm boxes time %d %15.8e  %15.8e  %15.8e  %15.8e\n", pp_mynode(),
	    add_time_get(311), add_time_get(312), add_time_get(313), add_time_get(314));

	DEBUG_LEAVE(rbox_communication_boxes)
}


LOCAL	boolean rbox_repair_intfc_in_box(
	INTERFACE  *intfc,
	int *smin,
	int *smax,
	int  flag[3][2])
{
	Table		*T = table_of_interface(intfc);
	RECT_GRID 	*gr = &topological_grid(intfc);
	COMPONENT 	*comp = T->components;
	RECT_BOX	box;
	int 		i,n_reg_node;
	int 		*gmax = gr->gmax;

	DEBUG_ENTER(rbox_repair_intfc_in_box)

	n_reg_node = (gmax[0]+1)*(gmax[1]+1)*(gmax[2]+1);
	for(i = 0; i < n_reg_node; ++i)
	    comp[i] = NO_COMP;

	adjust_crossings(smin,smax,intfc);
	
	fill_physical_comps(smin,smax,gmax,intfc);

	remove_unphysical_crossings3d(intfc, smin, smax);
	
	for(i=0; i<3; i++)
	{
	  box.smin[i] = smin[i];
	  box.smax[i] = smax[i];
	  
	  box.bmin[i] = smin[i];
	  if(!flag[i][0])
	    box.bmin[i] += RECON_BUFFER;
	  
	  box.bmax[i] = smax[i];
	  if(!flag[i][1])
	    box.bmax[i] -= RECON_BUFFER;
	}
	box.prev = box.next = NULL;
	box.grid = gr;

	print_int_vector("bmin", box.bmin, 3, "\n");
	print_int_vector("bmax", box.bmax, 3, "\n");
	
	set_current_interface(intfc);

	if(NO)
	{
	    add_to_debug("box_intfc");
	}
	
	box_index = 0;
	if (!grid_based_box_untangle(intfc,&box,NO))
	{
	  printf("WARNING, rbox_repair_intfc_in_box, "
	       "grid_based_box_untangle fails. \n");
	  DEBUG_LEAVE(rbox_repair_intfc_in_box)
	  return FUNCTION_FAILED;
	}
	    
	if (debugging("box_intfc"))
	{
	    remove_from_debug("box_intfc");
	    if (consistent_interface(intfc))
	    {
	    	(void) printf("After untangle, interface is consistent!\n");
	    }
	}

	DEBUG_LEAVE(rbox_repair_intfc_in_box)
	
	return FUNCTION_SUCCEEDED;
}

EXPORT	boolean	rbox_repair_intfc(Front*,int flag[3][2]);

/*topological_grid of the interface is already set, the boundary flags are flag
*/
EXPORT	boolean	rbox_repair_intfc(
	Front     *front,
	int	  flag[3][2])
{
	INTERFACE *intfc = front->interf;
	boolean      sav_intrp = interpolate_intfc_states(intfc);
	int 	  i, smin[3], smax[3];
	size_t	  sizest = front->sizest;

	DEBUG_ENTER(rbox_repair_intfc)

	set_size_of_intfc_state(sizest);
	set_copy_intfc_states(YES);

	make_interface_topology_lists(intfc);
	set_crx_storage_for_reconstruction(intfc,NULL);

	/* set reconstruction boundary and tolerance */

	for (i = 0; i < 3; ++i)
	{
	    smin[i] = 0;
            smax[i] = topological_grid(intfc).gmax[i];
	}
	
	start_clock("insert_grid_crossings3d");
	interpolate_intfc_states(intfc) = YES;
	insert_grid_intfc_crossings(intfc);

	stop_clock("insert_grid_crossings3d");
	print_storage("After insert_grid_crossings3d","crx_store");

	start_clock("reconstruct_intfc3d_in_box");
	
	if(!rbox_repair_intfc_in_box(intfc,smin,smax,flag))
	{
	    interpolate_intfc_states(intfc) = sav_intrp;
	    DEBUG_LEAVE(rbox_repair_intfc)
	    return NO;
	}
	
	stop_clock("reconstruct_intfc3d_in_box");

	interpolate_intfc_states(intfc) = sav_intrp;

	reset_intfc_num_points(intfc);
        free_crx_storage(intfc);
	
	front->interf = copy_interface(intfc);
	
	delete_interface(intfc);
	DEBUG_LEAVE(rbox_repair_intfc)
	return YES;
}

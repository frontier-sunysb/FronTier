/*
*                               fpatch3d.c:
*
*       Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*
*       this file includes functions set_patch_front, set_patch_topo_grid,
*                       and set_patch_comput_grid
*/

#define DEBUG_STRING    "fpatch3d"

#include <front/fdecs.h>
#include <front/fpatrecon.h>

//using namespace std;
//#include <map>

//AMR3D
EXPORT boolean ng_form_patch_subintfc_via_cut3d(Front*);
EXPORT void clip_intfc_at_grid_bdry1(INTERFACE*);
INTERFACE *cut_buf_interface1(INTERFACE*,int,int,int*,int*);
void	assign_point_index(INTERFACE*,int);
void	assign_index_of_point(INTERFACE*);
boolean	merge_buffer_interface(INTERFACE*,INTERFACE*,int);
int     compare_points(const void *, const void *);
//EXPORT  void    ft_tecplot_interface(FILE*,const char*,INTERFACE*);

struct _EDGETRI {
	TRI	*tri;
	int	edge;
};
typedef struct _EDGETRI	EDGETRI;

int compare_null_edges(const void *a, const void *b)
{
	EDGETRI	*a1=(EDGETRI*)a, *b1=(EDGETRI*)b;
	TRI	*t1 = a1->tri, *t2 = b1->tri;
	int	e1 = a1->edge, e2 = b1->edge;
	int	p11,p12, p21,p22, m1,m2;

	p11 = Index_of_point(Point_of_tri(t1)[e1]);
	p12 = Index_of_point(Point_of_tri(t1)[Next_m3(e1)]);
	
	p21 = Index_of_point(Point_of_tri(t2)[e2]);
	p22 = Index_of_point(Point_of_tri(t2)[Next_m3(e2)]);

	m1 = min(p11,p12);
	m2 = min(p21,p22);

	if(m1 < m2)
	  return -1;
	if(m1 > m2)
	  return 1;
	
	m1 = max(p11,p12);
	m2 = max(p21,p22);

	if(m1 < m2)
	  return -1;
	if(m1 > m2)
	  return 1;

	return 0;
}

//The compared variables should be rotationally invariant
int compare_tris_pointers(const void *a, const void *b)
{
	TRI	**a1=(TRI**)a, **b1=(TRI**)b;
	TRI	*t1 = *a1, *t2 = *b1;
	POINT	**p1, **p2;
	int	p1i[3], p2i[3], p1m, p2m, i, s1, s2;

	p1 = Point_of_tri(t1);
	p2 = Point_of_tri(t2);
	for(i=0; i<3; i++)
	{
	  p1i[i] = Index_of_point(p1[i]);
	  p2i[i] = Index_of_point(p2[i]);
	}

	//comp min point index
	p1m = min3(p1i[0],p1i[1],p1i[2]);
	p2m = min3(p2i[0],p2i[1],p2i[2]);
	if(p1m < p2m)
	    return -1;
	if(p1m > p2m)
	    return 1;

	//comp max point index
	p1m = max3(p1i[0],p1i[1],p1i[2]);
	p2m = max3(p2i[0],p2i[1],p2i[2]);
	if(p1m < p2m)
	    return -1;
	if(p1m > p2m)
	    return 1;

	s1 = p1i[0] + p1i[1] + p1i[2];
	s2 = p2i[0] + p2i[1] + p2i[2];
	
	if(s1 < s2)
	    return -1;
	if(s1 > s2)
	    return 1;

	return 0;
}


boolean	check_edge_orient(
	TRI	*t1,
	int	e1,
	TRI	*t2,
	int	e2)
{
	POINT	**p1, **p2;

	p1 = Point_of_tri(t1);
	p2 = Point_of_tri(t2);
	
	if(p1[e1] != p2[Next_m3(e2)] || p1[Next_m3(e1)] != p2[e2])
	  return NO;
	return YES;
}

void	check_double_edge_consistent(
	INTERFACE	*intfc)
{
	SURFACE			**s;
	EDGETRI			*edges, *e1, *e2;
	TRI			*tri;
	int			i, j, num, np;
	boolean			pteq;

	DEBUG_LEAVE(check_double_edge_consistent)

	//count number of edges
	num = 0;
	for(s=intfc->surfaces; s && *s; s++)
	  for(tri=first_tri(*s); !at_end_of_tri_list(tri,*s); tri=tri->next)
	    num++;
	
	num *= 3;

	if(num == 0)
	{
	  DEBUG_LEAVE(check_double_edge_consistent)
	  return;
	}

	uni_array(&edges, num, sizeof(EDGETRI));
	
	//put all edges in edges and sort
	i = 0;
	for(s=intfc->surfaces; s && *s; s++)
	  for(tri=first_tri(*s); !at_end_of_tri_list(tri,*s); tri=tri->next)
	    for(j=0; j<3; j++)
	    {
	      edges[i].edge = j;
	      edges[i].tri = tri;
	      i++;
	    }
	
	printf("#check_double_edge_consistent %d edges\n", num);

	assign_index_of_point(intfc);
	qsort(edges, num, sizeof(EDGETRI), compare_null_edges);

	np = 0;
	pteq = NO;
	for(i=0; i<num-1; i++)
	{
	    e1 = &edges[i];
	    e2 = &edges[i+1];
	    
	    if(compare_null_edges((const void*)e1, (const void*)e2) == 0)
	    {
	      if(pteq)
	      {
	        printf("ERROR check_double_edge_consistent, 3 edges are same.\n");
		
		printf("i %d\n", e1->edge);
		print_tri(e1->tri, intfc);
		
		printf("i+1 %d\n", e2->edge);
		print_tri(e2->tri, intfc);
		
		e2 = &edges[i-1];
		printf("i-1 %d\n", e2->edge);
		print_tri(e2->tri, intfc);

		clean_up(ERROR);
	      }
	      else
	      {
	        np++;
		pteq = YES;
		
		if(Tri_on_side(e1->tri,e1->edge) == NULL)
		{
		  printf("ERROR check_double_edge_consistent, e1 is null.\n");
		  clean_up(ERROR);
		}
		if(Tri_on_side(e2->tri,e2->edge) == NULL)
		{
		  printf("ERROR check_double_edge_consistent, e2 is null.\n");
		  clean_up(ERROR);
		}
		if(Tri_on_side(e1->tri,e1->edge) != e2->tri)
		{
		  printf("ERROR check_double_edge_consistent, tri on e1 != e2.\n");
		  clean_up(ERROR);
		}
		if(Tri_on_side(e2->tri,e2->edge) != e1->tri)
		{
		  printf("ERROR check_double_edge_consistent, tri on e2 != e1.\n");
		  clean_up(ERROR);
		}
	      } //if pteq
	    } //if compare_null_edges
	    else
	      pteq = NO;
	}

	printf("#check_double_edge_consistent, "
	       "%d total edges, %d edges pairs\n", num, np);

	free(edges);

	DEBUG_LEAVE(check_double_edge_consistent)
}

void	cut_intfc_in_grid(INTERFACE*,RECT_GRID*);

//gr is the topological grid of the domain (no buffer)
void	cut_intfc_in_grid(
	INTERFACE	*intfc,
	RECT_GRID	*gr)
{
	double	L[3], U[3];
	int	dir, nb;

	for(dir = 0; dir < 3; ++dir)
	{
	    L[dir] = gr->VL[dir];
	    U[dir] = gr->VU[dir];
	}
	
	for(dir = 0; dir < 3; ++dir)
	{
	    for (nb = 0; nb < 2; ++nb)
	    	open_null_sides1(intfc,L,U,dir,nb);
	}
	
	cut_out_curves_in_buffer(intfc);
	reset_intfc_num_points(intfc);
}

EXPORT boolean ng_form_patch_subintfc_via_cut3d(Front *fr)
{
	INTERFACE	*intfc = fr->interf;
	INTERFACE	*cur_intfc = current_interface();
	RECT_GRID	*gr = computational_grid(intfc);
	int		dim = intfc->dim;
	int		dir, nb;
	double		L[MAXD],U[MAXD];

	DEBUG_ENTER(ng_form_patch_subintfc_via_cut3d)
	
	strip_subdomain_bdry_curves(intfc);
	set_current_interface(intfc);
	for (dir = 0; dir < dim; ++dir)
	{
	    L[dir] = gr->VL[dir];
	    U[dir] = gr->VU[dir];
	    if (gr->lbuf[dir] == 0) L[dir] -= 0.5*gr->h[dir];
	    if (gr->ubuf[dir] == 0) U[dir] += 0.5*gr->h[dir];
	}
	
	for (dir = 0; dir < dim; ++dir)
	{
	    for (nb = 0; nb < 2; ++nb)
	    	open_null_sides1(intfc,L,U,dir,nb);
	}
	cut_out_curves_in_buffer(intfc);
	reset_intfc_num_points(intfc);
	printf("#ng_form_patch_subintfc_via_cut3d, recon %d \n", 
	    interface_reconstructed(intfc));
	interface_reconstructed(intfc) = NO;
	set_current_interface(cur_intfc);
	
	DEBUG_LEAVE(ng_form_patch_subintfc_via_cut3d)

	return YES;
}

#if defined(USE_AMR)

EXPORT	boolean f_patch_intfc_communication3d(Front**,FT_NeighborPatchInfo*,int,int);

#define	MAX_NEIGHBOR_PROCS	20
#define	MAX_NEIGHBOR_INTFC	20

//remove repeated proc numbers.
int	get_procs(
	int			*procs,
	FT_NeighborFaceInfo	*dnp,
	int			dnn)
{
	int	i, j, n;

	n = 0;
	for(i=0; i<dnn; i++)
	{
	  for(j=0; j<n; j++)
	    if(procs[j] == dnp[i].proc)
	      break;
	  if(j == n)
	  {
	    procs[j] = dnp[i].proc;
	    n++;
	  }
	  
	  if(n >= MAX_NEIGHBOR_PROCS)
	  {
	    printf("ERROR get_procs, too many procs.\n");
	    clean_up(ERROR);
	  }
	}
	
	return n;
}


EXPORT	boolean f_patch_intfc_communication3d(
	Front			**frs,
	FT_NeighborPatchInfo	*npi,
	int			nfr,
	int			myid)
{
	INTERFACE	*sav_intfc, *intfc, *buf_intfc;
	INTERFACE	*rintfcs[MAX_NEIGHBOR_INTFC];
	RECT_GRID	*gr;
	boolean		sav_copy;
	double		T;
	int		me[3], him[3];
	int		i, j, k, jp, n;
	int  	        dnn, patch_num;
	int		nproc, proc, procs[MAX_NEIGHBOR_PROCS];
	FT_NeighborFaceInfo	*dnp;
	map<int,int>		procmap;  	   //map patch number to proc
	map<int,int>::iterator	it;
	map<int,INTERFACE*>	recvmap[2];	   //map patch number to intfc
	map<int,INTERFACE*>::iterator	itr;
	char  fname[128];
	static	int	cnt = 0;

	DEBUG_ENTER(f_patch_intfc_communication3d)

	cnt++;

	sav_copy = copy_intfc_states();
	sav_intfc = current_interface();

	//only h and GL GU are used in gr
	gr = frs[0]->rect_grid;
	set_floating_point_tolerance1(gr->h);
	set_copy_intfc_states(YES);
	
	//just use me and him to call same functions.
	for(i=0; i<3; i++)
	{
	  me[i] = myid;
	  him[i] = myid+1;
	}
	
	//clip patch intfcs to patch interior
	for(k=0; k<nfr; k++)
	{
	  clip_intfc_at_grid_bdry1(frs[k]->interf);
	  assign_point_index(frs[k]->interf, npi[k].patch_number);
	}

	for (i = 0; i < 3; ++i)
	{
	  for (j = 0; j < 2; ++j)
	  {
	    //should have sync for multiple proc run
	    //pp_gsync();

	    printf("#before send %d %d\n", i, j);
	    
	    jp = (j+1)%2;
	    
	    recvmap[jp].clear();
	    
	    //send all patch intfcs in this proc to the other procs on side j
	    for(k=0; k<nfr; k++)
	    {
	      dnp = npi[k].d_neighbor_patches[i][j];
	      dnn = npi[k].d_num_neighbors[i][j];
	      patch_num = npi[k].patch_number;
	      intfc = frs[k]->interf;

	      //has neighbor patch in direction i, side j
	      if(dnn != 0)
	      {
	        buf_intfc = cut_buf_interface1(intfc,i,j,me,him);
		
		//patch boundary is periodic
		if(rect_boundary_type(intfc,i,j) == SUBDOMAIN_BOUNDARY)
		{
		  T = gr->GU[i] - gr->GL[i];
		  if(j == 1) 
		    T = -T;
	          shift_interface(buf_intfc,T,i);
		}
		
		//mapping from patch_number to interface
		recvmap[jp][patch_num] = buf_intfc;
		
		//send buf_intfc to all the neighbor patches
		nproc = get_procs(procs, dnp, dnn);
		for(n=0; n<nproc; n++)
		{
		  if(procs[n] != myid)
		  {
		    //amr send interface to procs[n]
		    //send_interface(buf_intfc,dst_id);
		    //(void) delete_interface(buf_intfc);
		  }
		} //for n=0 to dnn
	      } //if dnn != 0
	    } //for k=0 to nfr

	    //receive all patch intfcs from other procs on side (j+1)%2
	    printf("#before recv %d %d\n", i, j);
	    
	    //construct patch_number - processor mapping for all
	    //the patch interfaces
	    procmap.clear();
	    
	    for(k=0; k<nfr; k++)
	    {
	      dnp = npi[k].d_neighbor_patches[i][jp];
	      dnn = npi[k].d_num_neighbors[i][jp];

	      for(n=0; n<dnn; n++)
	        procmap[dnp[n].patch_number] = dnp[n].proc;
	    }

	    //receive intfcs
	    for(it=procmap.begin(); it != procmap.end(); it++)
	    {
	      patch_num = (*it).first;
	      proc = (*it).second;
	      
	      //if patch_num is in this proc, recvmap is assigned in the send part
	      if(proc != myid)
	      {
	        //amr receive intfc with patch number patch_num from proc
	        //recvmap[jp][patch_num] = receive_interface(proc, patch_num);
	      }
	    }

	    printf("#after recv\n");
          } // for j=0 to 2 side

	  printf("\n");

          for(j=0; j<2; j++)
	  {
	    printf("#before merge %d %d\n", i, j);
	    
	    //merge interfaces
	    for(k=0; k<nfr; k++)
	    {
	      dnp = npi[k].d_neighbor_patches[i][j];
	      dnn = npi[k].d_num_neighbors[i][j];

	      if(dnn == 0)
	        continue;
	      if(dnn > MAX_NEIGHBOR_INTFC)
	      {
	        printf("ERROR f_patch_intfc_communication3d, dnn=%d\n", dnn);
		clean_up(ERROR);
	      }

	      //get buffer interfaces for frs[k] from recvmap
	      for(n=0; n<dnn; n++)
	        rintfcs[n] = recvmap[j][dnp[n].patch_number];

	      printf("#patch intfc %d dnn=%d\n", k, dnn);

	      if(k == 4 && debugging("surfchk"))
	      {
	        sprintf(fname, "merge_intfc_%02d_%02d_%02d", k,i,j);
	        tecplot_merge_surfaces(fname,frs[k]->interf,rintfcs,dnn,myid);
		add_to_debug("nulledge");
	      }
	      
	      //merge frs[k]->interf with the received intfcs
	      for(n=0; n<dnn; n++)
	        merge_buffer_interface(frs[k]->interf,
		         rintfcs[n],dnp[n].patch_number);
	      
	      //buffer_extension3d1(frs[k]->interf,rintfcs[0],i,j,0);
	      
	      if(k == 4 && debugging("surfchk"))
	      {
		remove_from_debug("nulledge");
	        sprintf(fname, "merged_intfc_%02d_%02d_%02d_", k,i,j);
	        null_sides_are_consistent();
	        check_print_intfc("After merge_surfaces", fname, 
	                      's', frs[k]->interf, 0,0, NO);
	      }
	      
	      printf("#patch intfc end %d \n\n", k);
	    }

/*
	    for(k=0; k<nfr; k++)
	    {
	      sprintf(fname, "merged_intfc_%02d_%02d_%02d_", k,i,j);
	      null_sides_are_consistent();
	      check_print_intfc("After merge_surfaces", fname, 
	                   's', frs[k]->interf, 0,0, NO);
	    }
*/
	    printf("#after merge %d %d\n\n", i, j);

	    for(itr=recvmap[j].begin(); itr != recvmap[j].end(); itr++)
	      (void) delete_interface((*itr).second);
	  } // for j = 0 to 2 sides
	
	} //for i = 0 to 3 directions

	for(k=0; k<nfr; k++)
	{
	  ng_form_patch_subintfc_via_cut3d(frs[k]);
	  reset_intfc_num_points(frs[k]->interf);
	  reset_normal_on_intfc(frs[k]->interf);
	}

comm_exit:
	set_copy_intfc_states(sav_copy);
	set_current_interface(sav_intfc);
	
	DEBUG_LEAVE(f_patch_intfc_communication3d)
	
	return YES;
}	//f_patch_intfc_communication3d


EXPORT boolean InitNeighborPatchInfo(FT_NeighborPatchInfo *d_neighbor_patch_info)
{
	FT_NeighborFaceInfo *d_neighbor_store = d_neighbor_patch_info->d_neighbor_store;
	int		    n_patches = d_neighbor_patch_info->num_patches;
	int		i,j,ind,side,nb;
	static	int	sam_side[] = {0,0,1,1,2,2};  //SAMRAI conversion
	static	int	sam_nb[] = {0,1,0,1,0,1};
	

	for(side=0; side<3; side++)
	    for(nb=0; nb<2; nb++)
	    {
		d_neighbor_patch_info->d_num_neighbors[side][nb] = 0;
		d_neighbor_patch_info->d_neighbor_patches[side][nb] = 0;
	    }

	if(n_patches == 0)
	    return YES;

	//sort neighbor face patches according to their index;
	qsort((void*)d_neighbor_store, n_patches,
	      sizeof(FT_NeighborFaceInfo), compare_face_patches);
	
	i = 0;
	while(i<n_patches)
	{
	    ind = d_neighbor_store[i].bdry_index;
	    //find how many patches in the side with index ind;
	    for(j=i+1; j<n_patches && d_neighbor_store[j].bdry_index==ind; ++j);
	    
	    side = sam_side[ind];
	    nb = sam_nb[ind];
	    d_neighbor_patch_info->d_neighbor_patches[side][nb] = &d_neighbor_store[i];
	    d_neighbor_patch_info->d_num_neighbors[side][nb] = j - i;
	    i = j;
	}

	return YES;
}

int	compare_face_patches(const void *a, const void *b)
{
	FT_NeighborFaceInfo *c1=(FT_NeighborFaceInfo*)a, *c2=(FT_NeighborFaceInfo*)b;
	return c1->bdry_index - c2->bdry_index;
}

EXPORT	void PrintNeighborFaceInfo(FT_NeighborFaceInfo *nfi)
{
	printf("nb face info %d %d %d\n", nfi->bdry_index, nfi->proc,
	       nfi->patch_number);
}


EXPORT void PrintNeighborPatchInfo(FT_NeighborPatchInfo *npi)
{
	int	dir,nb,num,n;

	printf("Total number of neighbor patches %d\n", npi->num_patches);
	
	//no neighbor patches
	if(npi->num_patches == 0)
	    return;

	for(dir=0; dir<3; dir++)
	    for(nb=0; nb<2; nb++)
	    {
		num = npi->d_num_neighbors[dir][nb];
		printf("%d %d  num=%d\n", dir, nb, num);
		for(n=0; n<num; n++)
		{
		    PrintNeighborFaceInfo(&npi->d_neighbor_patches[dir][nb][n]);
		}
	    }
	printf("\n");
}

#endif  //#if defined(USE_AMR)

void	tecplot_one_interface(
	FILE		*file,
	INTERFACE	*intfc)
{
	SURFACE	**s;
	CURVE	**cc;

	//print curves
	for(cc = intfc->curves; cc && *cc; ++cc)
	{
	    //printf("Curve %d: cc %p *cc %p\n",i,cc,*cc);
	    tecplot_curve(NULL,file,*cc);
	}

	//print surfaces
	for (s = intfc->surfaces; s && *s; ++s)
	    tecplot_surface(NULL,file,*s);
}

void	tecplot_merge_surfaces(char*,INTERFACE*,INTERFACE**,int,int);

void	tecplot_merge_surfaces(
	char		*bname,
	INTERFACE	*intfc,
	INTERFACE	**rintfc,
	int		n,
	int		myid)
{
	FILE	*file;
	char	fname[128];
	int	i;

	sprintf(fname, "%s_%05d.plt", bname, myid);
	printf("#tecplot_merge_surfaces %s\n", fname);

	if ((file = fopen(fname,"w")) == NULL)
	{
	    screen("WARNING in tecplot_merge_surfaces(), "
	           "can't open %s\n",fname);
	    return;
	}
	
	(void) fprintf(file,"TITLE = \"tecplot merge surfaces\"\n"
		   	    "VARIABLES = \"x\", \"y\", \"z\"\n");

	tecplot_one_interface(file,intfc);
	for(i=0; i<n; i++)
	    tecplot_one_interface(file,rintfc[i]);

	fclose(file);
}

void	tecplot_box_interface(char*,INTERFACE*,RECON_BOX*,int);

void	tecplot_box_interface(
	char		*bname,
	INTERFACE	*intfc,
	RECON_BOX	*rbox,
	int		myid)
{
	FILE	*file;
	char	fname[128];
	int	i;

	sprintf(fname, "%s_%05d.plt", bname, myid);
	printf("#tecplot_merge_surfaces %s\n", fname);

	if ((file = fopen(fname,"w")) == NULL)
	{
	    screen("WARNING in tecplot_merge_surfaces(), "
	           "can't open %s\n",fname);
	    return;
	}
	
	(void) fprintf(file,"TITLE = \"tecplot merge surfaces\"\n"
		   	    "VARIABLES = \"x\", \"y\", \"z\"\n");

	tecplot_box(NULL, file, rbox->fmin, rbox->fmax);

	tecplot_one_interface(file,intfc);

	fclose(file);
}

void	bucket_set_grid(
	RECT_GRID	*ggr,
	RECT_GRID	*gr,
	POINT		**pts,
	int		num)
{
	int	i,j;
	double	L[3], U[3];

	for(i = 0; i < 3; ++i)
	{
	  L[i] = HUGE_VAL;
	  U[i] = -HUGE_VAL;
	}

	for(j=0; j<num; j++)
	{
	  for (i = 0; i < 3; ++i)
	  {
	    L[i] = min(L[i], Coords(pts[j])[i]);
	    U[i] = max(U[i], Coords(pts[j])[i]);
	  }
	}

	for(i=0; i<3; i++)
	{
	  L[i] -= gr->h[i];
	  U[i] += gr->h[i];
	  ggr->gmax[i] = ceil((U[i]-L[i])/gr->h[i]);
	  if(ggr->gmax[i] > 2*gr->gmax[i])
	  {
	    printf("ERROR set_bucket_grid, grid is too large.\n");
	    clean_up(ERROR);
	  }
	  U[i] = L[i] + ggr->gmax[i]*gr->h[i];
	}
	
	set_rect_grid(L,U,gr->GL,gr->GU,NULL,NULL,ggr->gmax,3,&gr->Remap,ggr);
}

void	bucket_init(
	POINT_BUCKET	*bpts,
	RECT_GRID	*gr,
	POINT		**pts,
	int		num)
{
	RECT_GRID	*ggr;
	POINT		*p;
	int		*gmax;
	int		i, j, k, n, np, pst;
	int		**ptcrds;
	double		tol;

	DEBUG_ENTER(bucket_init)
	
	//tol is set here
	tol = 4.0e-3;
	for(i=0; i<3; i++)
	  bpts->btol[i] = gr->h[i]*tol;

	ggr = &bpts->ggr;
	
	bucket_set_grid(ggr, gr, pts, num);
	
	gmax = ggr->gmax;
	
	tri_array(&bpts->np, gmax[0], gmax[1], gmax[2], INT);
	tri_array(&bpts->pt, gmax[0], gmax[1], gmax[2], sizeof(POINT**));
	uni_array(&bpts->bpt, num, sizeof(POINT*));
	bi_array(&ptcrds, num, 3, INT);

	//pts[i] belongs to which bucket
	//how many points in each bucket
	//ptcrds[i] store the block info for pts[i]
	for(n=0; n<num; n++)
	{
	  p = pts[n];

	  for(j=0; j<3; j++)
	    ptcrds[n][j] = cell_index(Coords(p)[j],j,ggr);
	  bpts->np[ptcrds[n][0]][ptcrds[n][1]][ptcrds[n][2]]++;
	}

	//set pointers for each bucket
	//set np to zero
	pst = 0;
	for(i=0; i<gmax[0]; i++)
	  for(j=0; j<gmax[1]; j++)
	    for(k=0; k<gmax[2]; k++)
	    {
	      if(bpts->np[i][j][k] == 0)
	        continue;
	      bpts->pt[i][j][k] = &bpts->bpt[pst];
	      pst += bpts->np[i][j][k];
	      bpts->np[i][j][k] = 0;
	    }
	
	//save pts to bpts->pt
	//compute np
	for(n=0; n<num; n++)
	{
	  i = ptcrds[n][0];
	  j = ptcrds[n][1];
	  k = ptcrds[n][2];
	  
	  np = bpts->np[i][j][k];
	  bpts->pt[i][j][k][np] = pts[n];
	  bpts->np[i][j][k]++;
	}

	if(NO)
	{
	printf("#bucket print\n");
	for(i=0; i<gmax[0]; i++)
	  for(j=0; j<gmax[1]; j++)
	    for(k=0; k<gmax[2]; k++)
	    {
	      if(bpts->np[i][j][k] == 0)
	        continue;
	      printf("%3d %3d %3d  %3d  | ", i, j, k, bpts->np[i][j][k]);
	      for(n=0; n<bpts->np[i][j][k]; n++)
	        printf("%d  ", bpts->pt[i][j][k][n]);
	      printf("\n");
	    }
	}

	free(ptcrds);
	
	DEBUG_LEAVE(bucket_init)
}

boolean	bucket_match_pts(
	POINT_BUCKET	*bpts,
	POINT		*p,
	int		ix,
	int		iy,
	int		iz)
{
	int	i,j,k,n,d;
	int	imin,imax,jmin,jmax,kmin,kmax,*gmax;
	boolean	found;
	POINT	*npt, *p1;
	double	dist, dist1;

	DEBUG_ENTER(bucket_match_pts)
	
	gmax = bpts->ggr.gmax;
	imin = max(ix-1,0);
	jmin = max(iy-1,0);
	kmin = max(iz-1,0);
	imax = min(ix+1,gmax[0]-1);
	jmax = min(iy+1,gmax[1]-1);
	kmax = min(iz+1,gmax[2]-1);

	//printf("#match %d %d %d  %d\n", ix, iy, iz, p);
	
	dist = HUGE_VAL;
	found = NO;

	for(i=imin; i<=imax; i++)
	  for(j=jmin; j<=jmax; j++)
	    for(k=kmin; k<=kmax; k++)
	      for(n=0; n<bpts->np[i][j][k]; n++)
	      {
	        //printf("%d %d %d %d  %d\n", i,j,k,n,bpts->pt[i][j][k]);

	        npt = bpts->pt[i][j][k][n];
		
		//printf("npt = %d\n", npt);

		//npt already matched with one point
		if(Cor_point(npt) != NULL)
		  continue;
		
		for(d=0; d<3; d++)
		{
		  if(fabs(Coords(npt)[d]-Coords(p)[d]) > (bpts->btol)[d])
		    break;
		}

		//If there are many points in the range of tol, we pick the nearest
		if(d == 3)
		{
		  dist1 = distance_between_positions(Coords(p),Coords(npt),3);
		  if(dist1 < dist)
		  {
		    p1 = npt;
		    dist = dist1;
		  }
		  found = YES;
		}
	      }

	if(found)
	{
	  if(the_point(p))
	  {
	    printf("#pt fd %d %d\n", p, p1);
	    print_general_vector("p=", Coords(p), 3, "\n");
	    print_general_vector("p1=", Coords(p1), 3, "\n");
	  }

	  if(p->indx == p1->indx)
	  {
	    average_points(NO,p,p->hse,p->hs,p1,p1->hse,p1->hs);
	    Cor_point(p) = p1;
	    Cor_point(p1) = p1;
	  }
	  else if(p->indx < p1->indx)
	  {
	    Cor_point(p) = p1;
	    Cor_point(p1) = p1;
	  }
	  else
	  {
	    Cor_point(p) = p;
	    Cor_point(p1) = p;
	  }
	}

	DEBUG_LEAVE(bucket_match_pts)
	
	return found;
}

void	bucket_release(
	POINT_BUCKET    *bpts)
{
	free_these(3, bpts->np, bpts->pt, bpts->bpt);
}

int	count_surface_points(
	SURFACE		*s)
{
	int	npt, i;
	TRI	*tri;

	for (tri=first_tri(s); !at_end_of_tri_list(tri,s); tri=tri->next)
	{
	    for (i = 0; i < 3; i++)
		Index_of_point(Point_of_tri(tri)[i]) = -1;
	}
	
	npt = 0;
	for (tri=first_tri(s); !at_end_of_tri_list(tri,s); tri=tri->next)
	{
	    for (i = 0; i < 3; i++)
	    {
		if(Index_of_point(Point_of_tri(tri)[i]) == -1)
		{
		  Index_of_point(Point_of_tri(tri)[i]) = 0;
		  npt++;
	        }
	    }
	}

	return npt;
}

int	get_surface_points(
	POINT		**pts,
	SURFACE		*s)
{
	int	npt, i;
	TRI	*tri;

	for (tri=first_tri(s); !at_end_of_tri_list(tri,s); tri=tri->next)
	{
	    for (i = 0; i < 3; i++)
		Index_of_point(Point_of_tri(tri)[i]) = -1;
	}
	
	npt = 0;
	for (tri=first_tri(s); !at_end_of_tri_list(tri,s); tri=tri->next)
	{
	    for (i = 0; i < 3; i++)
	    {
		if(Index_of_point(Point_of_tri(tri)[i]) == -1)
		{
		  Index_of_point(Point_of_tri(tri)[i]) = 0;
		  pts[npt] = Point_of_tri(tri)[i];
		  pts[npt]->hs = Hyper_surf(s);
		  pts[npt]->hse = Hyper_surf_element(tri);
		  npt++;
	        }
	    }
	}

	return npt;
}

void	replace_tri_points(
	SURFACE		*s)
{
	int	i;
	TRI	*tri;
	POINT	**p;

	for(tri=first_tri(s); !at_end_of_tri_list(tri,s); tri=tri->next)
	{
	  for(i = 0; i < 3; i++)
	  {
	    p = Point_of_tri(tri);
	    p[i] = (POINT*)Cor_point(p[i]);
	  }
	}
}

void	merge_surface_points(
	INTERFACE	*intfc,
	SURFACE		*s,
	SURFACE		*sa)
{
	RECT_GRID		*gr, *ggr;
	HYPER_SURF              *hs;
        HYPER_SURF_ELEMENT      *hse;
        POINT                   *p, **pts, **ptsa;
	POINT_BUCKET		bpts;
	int			i, j, k, n, num, numa, np;

	DEBUG_ENTER(merge_surface_points)
	
	num = count_surface_points(s);
	if(num == 0)
	{
	    DEBUG_LEAVE(merge_surface_points)
	    return;
	}
	
	numa = count_surface_points(sa);
	if(numa == 0)
	{
	    DEBUG_LEAVE(merge_surface_points)
	    return;
	}

	uni_array(&pts, num, sizeof(POINT*));
	uni_array(&ptsa, numa, sizeof(POINT*));

	get_surface_points(pts, s);
	get_surface_points(ptsa, sa);
	
	printf("#total pts num = %d  numa = %d\n", num, numa);

	gr = &topological_grid(intfc);
	//make bucket for pts
	bucket_init(&bpts, gr, pts, num);
	ggr = &bpts.ggr;
	
	for(i=0; i<num; i++)
	  Cor_point(pts[i]) = NULL;
	
	for(i=0; i<numa; i++)
	  Cor_point(ptsa[i]) = NULL;
	
	np = 0;
	for(n=0; n<numa; n++)
	{
	  p = ptsa[n];

	  //which bucket p belongs to
	  i = cell_index(Coords(p)[0],0,ggr);
	  j = cell_index(Coords(p)[1],1,ggr);
	  k = cell_index(Coords(p)[2],2,ggr);
	  if(bucket_match_pts(&bpts,p,i,j,k))
	    np++;
	}

	for(i=0; i<num; i++)
	  if(Cor_point(pts[i]) == NULL)
	    Cor_point(pts[i]) = pts[i];
	
	for(i=0; i<numa; i++)
	  if(Cor_point(ptsa[i]) == NULL)
	    Cor_point(ptsa[i]) = ptsa[i];

	printf("#merged points %d\n", np);
	
	if(NO)
	{
	  printf("#ptsout af\n");
	  
	  for(i=0;i<num;i++)
	  {
	    p = pts[i];
	    if(Coords(p)[0] > 2.45 && p != Cor_point(p))
	    {
	      printf("%d  %d  %2d   % 15.8e, % 15.8e, % 15.8e\n", 
	           p, Cor_point(p), p->indx,
	           Coords(p)[0], Coords(p)[1], Coords(p)[2]);
	      
	      p = (POINT*)Cor_point(p);
	      printf("%d  %d  %2d   % 15.8e, % 15.8e, % 15.8e\n", 
	           p, Cor_point(p), p->indx,
	           Coords(p)[0], Coords(p)[1], Coords(p)[2]);
	    }
	  }
	}

	bucket_release(&bpts);
	free_these(2, pts, ptsa);
	
	DEBUG_LEAVE(merge_surface_points)
}

void	assign_surface_point_index(
	SURFACE		*s,
	int		index)
{
	TRI		*tri;
	int		i;

	for(tri=first_tri(s); !at_end_of_tri_list(tri,s); tri=tri->next)
	  for(i=0; i<3; i++)
	  {
	    Point_of_tri(tri)[i]->indx = index;
	  }
}

void	assign_point_index(INTERFACE*,int);

void	assign_point_index(
	INTERFACE	*intfc,
	int		index)
{
	SURFACE		**s;

	for(s = intfc->surfaces; s && *s; ++s)
	  assign_surface_point_index(*s, index);
}

void	merge_surface_points_prev(
	INTERFACE	*intfc)
{
	HYPER_SURF              *hs;
        HYPER_SURF_ELEMENT      *hse;
        POINT                   *p, **pts, *p1, *p2;
	int			i, num, np;
	boolean			pteq;

	DEBUG_ENTER(merge_surface_points)

	reset_intfc_num_points(intfc);

	num = intfc->num_points;
	if(num == 0)
	{
	    DEBUG_LEAVE(merge_surface_points)
	    return;
	}

	uni_array(&pts,num,sizeof(POINT*));

	i = 0;
	next_point(intfc,NULL,NULL,NULL);
	while(next_point(intfc,&p,&hse,&hs))
	{
	    p->hse = hse;
	    p->hs = hs;
	    pts[i] = p;
	    i++;
	}

	{
	  printf("#ptsout bf\n");
	  
	  for(i=0;i<num;i++)
	  {
	    p = pts[i];
	    printf("%d  %d  %2d   % 15.8e, % 15.8e, % 15.8e\n", 
	           p, Cor_point(p), p->indx,
	           Coords(p)[0], Coords(p)[1], Coords(p)[2]);
	  }
	}

	qsort(pts, num, sizeof(POINT*), compare_points);
	
	//if(debugging("ptsout"))
	{
	  printf("#ptsout\n");
	  
	  for(i=0;i<num;i++)
	  {
	    p = pts[i];
	    printf("%d  %d  %2d   % 15.8e, % 15.8e, % 15.8e\n", 
	           p, Cor_point(p), p->indx,
	           Coords(p)[0], Coords(p)[1], Coords(p)[2]);
	    if(i == 7278)
	    {
	      printf("#compare %d\n", compare_points((const void*)&pts[i], (const void*)&pts[i+1]));
	      printf("#compare %d\n", compare_points((const void*)&pts[i], (const void*)&pts[i+2]));
	      printf("#compare %d\n", compare_points((const void*)&pts[i], (const void*)&pts[i+3]));
	      printf("#compare %d\n", compare_points((const void*)&pts[i], (const void*)&pts[i+4]));
	      printf("#compare %d\n", compare_points((const void*)&pts[i+1], (const void*)&pts[i+2]));
	    }
	  }
	}

	printf("#total pts %d\n", num);

	for(i=0; i<num; i++)
	  Cor_point(pts[i]) = pts[i];
	
	np = 0;
	pteq = NO;
	for(i=0; i<num-1; i++)
	{
	    p1 = pts[i];
	    p2 = pts[i+1];
	    
	    if(compare_points((const void*)&p1, (const void*)&p2) == 0)
	    {
	      if(pteq)
	      {
	        printf("ERROR merge_surface_points, 3 points are same.\n");
		clean_up(ERROR);
	      }
	      else
	      {
	        np++;
		pteq = YES;
	        //p1 and p2 are from the same patch, periodic bdry case.
		if(p1->indx == p2->indx)
		{
		  average_points(NO,p1,p1->hse,p1->hs,p2,p2->hse,p2->hs);
		  Cor_point(p1) = p2;
		  Cor_point(p2) = p2;
		}
		else if(p1->indx < p2->indx)
		{
		  Cor_point(p1) = p2;
		  Cor_point(p2) = p2;
		}
		else
		{
		  Cor_point(p1) = p1;
		  Cor_point(p2) = p1;
		}
	      } //if pteq
	    } //if compare_points
	    else
	      pteq = NO;
	}

	printf("#merged points %d\n", np);

	free(pts);
	
	DEBUG_LEAVE(merge_surface_points)
}


void	assign_index_of_point(
	INTERFACE	*intfc)
{
	SURFACE		**s;
	POINT		**p;
	TRI		*tri;
	int		i,j;

	//init point index
	for(s=intfc->surfaces; s && *s; s++)
	  for(tri=first_tri(*s); !at_end_of_tri_list(tri,*s); tri=tri->next)
	    for(j=0; j<3; j++)
	      Index_of_point(Point_of_tri(tri)[j]) = -1;
	
	//assign point index, will be used by compare_tris_pointers
	i = 0;
	for(s=intfc->surfaces; s && *s; s++)
	  for(tri=first_tri(*s); !at_end_of_tri_list(tri,*s); tri=tri->next)
	    for(j=0; j<3; j++)
	    {
	      p = Point_of_tri(tri);
	      if(Index_of_point(p[j]) == -1)
	      {
		Index_of_point(p[j]) = i;
	        i++;
	      }
	    }
}

boolean	merge_tris(
	INTERFACE	*intfc)
{
	SURFACE	**s;
	TRI	**tris,**trisa,**trisb;
	TRI	*ta,*tb,*tri;
	POINT	**p;
	int	i,j,k,num,nt;
	boolean	pteq;

	DEBUG_ENTER(merge_tris)
	
	add_time_start(331);
	num = 0;
	for(s=intfc->surfaces; s && *s; s++)
	  for(tri=first_tri(*s); !at_end_of_tri_list(tri,*s); tri=tri->next)
	    num++;
	
	if(num == 0)
	{
	  DEBUG_LEAVE(merge_tris)
	  return YES;
	}

	uni_array(&tris,num,sizeof(TRI*));
	uni_array(&trisa,num,sizeof(TRI*));
	uni_array(&trisb,num,sizeof(TRI*));

	//put all triangles in tris[]
	i = 0;
	for(s=intfc->surfaces; s && *s; s++)
	  for(tri=first_tri(*s); !at_end_of_tri_list(tri,*s); tri=tri->next)
	  {
	    tris[i] = tri;
	    i++;
	  }
	add_time_end(331);

	add_time_start(332);
	//sort tris
	assign_index_of_point(intfc);
	qsort(tris, num, sizeof(TRI*), compare_tris_pointers);
	add_time_end(332);
	
	printf("#num tris %d\n", num);

	add_time_start(333);
	//put matched tris in trisa[] and trisb[] (nt)
	nt = 0;
	pteq = NO;
	for(i=0; i<num-1; i++)
	{
	  ta = tris[i];
	  tb = tris[i+1];
	  if(compare_tris_pointers((const void*)&ta, (const void*)&tb) == 0)
	  {
	    if(pteq)
	    {
	      printf("ERROR merge_tris, 3 tris are same.\n");
	      clean_up(ERROR);
	    }
	    else
	    {
	      pteq = YES;
	      trisa[nt] = ta;
	      trisb[nt] = tb;
	      nt++;
	    } //if pteq
	  } //if compare_points
	  else
	    pteq = NO;
	}
	
	printf("#nt tris %d\n", nt);
	add_time_end(333);
	
	add_time_start(334);
	//assign Tri_index and sync tris
	for(i=0; i<num; i++)
	  Tri_index(tris[i]) = -1;

	for(i=0; i<nt; i++)
	{
	  ta = trisa[i];
	  tb = trisb[i];
	  Tri_index(ta) = i;
	  Tri_index(tb) = i;

	  for(j=0; j<3; j++)
	    if(Point_of_tri(ta)[0] == Point_of_tri(tb)[j])
	      break;

	  rotate_triangle(ta,(3-j)%3);
	}
	add_time_end(334);

// 4 cases while merging
//	Tri_on_side(ta,j)	Tri_on_side(tb,j)	operations
// 1		NULL		      NULL		  NO
// 2		NULL		     NOT NULL 		connect to ta
// 3	       NOT NULL		      NULL		  NO
// 4	       NOT NULL		     NOT NULL	neighbor j must have matched tris

	add_time_start(335);
	for(i=0; i<nt; i++)
	{
	  ta = trisa[i];
	  tb = trisb[i];
	  for(j=0; j<3; j++)
	  {
	    tri = Tri_on_side(tb,j);
	    if(tri == NULL)
	      continue;
	    
	    if(Tri_on_side(ta,j) == NULL)
	    {
	      Tri_on_side(ta,j) = tri;
	      for(k=0; k<3; k++)
	      {
	        if(Tri_on_side(tri,k) == tb)
	    	  Tri_on_side(tri,k) = ta;
	      }
	    }
	    else
	    {
	      if(Tri_index(Tri_on_side(ta,j)) != Tri_index(tri))
	      {
	        printf("ERROR merge_tris, "
		      "can not match tris with buffer surface.\n");
		print_tri(Tri_on_side(ta,j),intfc);
		print_tri(tri,intfc);
		clean_up(ERROR);
	      }
	    } //if Tri_on_side(ta,j) == NULL
	  } //for j=0; j<3
	  
	  remove_tri_from_surface(tb,tb->surf,NO);
	} //for i=0; i<nt
	add_time_end(335);
	
	free_these(3, tris, trisa, trisb);
	
	DEBUG_LEAVE(merge_tris)
	return YES;
}

void	merge_null_edges(
	INTERFACE	*intfc)
{
	SURFACE			**s;
	EDGETRI			*edges, *e1, *e2;
	TRI			*tri;
	int			i, j, num, np;
	boolean			pteq;

	DEBUG_ENTER(merge_null_edges)

	//count and put all null edges in edges[] (num)
	num = 0;
	for(s=intfc->surfaces; s && *s; s++)
	  for(tri=first_tri(*s); !at_end_of_tri_list(tri,*s); tri=tri->next)
	    for(j=0; j<3; j++)
	      if(Tri_on_side(tri,j) == NULL)
	        num++;
	
	if(num == 0)
	{
	  DEBUG_LEAVE(merge_null_edges)
	  return;
	}

	uni_array(&edges, num, sizeof(EDGETRI));

	//assign EDGETRI structure
	i = 0;
	for(s=intfc->surfaces; s && *s; s++)
	  for(tri=first_tri(*s); !at_end_of_tri_list(tri,*s); tri=tri->next)
	    for(j=0; j<3; j++)
	      if(Tri_on_side(tri,j) == NULL)
	      {
	        edges[i].edge = j;
	        edges[i].tri = tri;
		i++;
	      }
	
	//sort edges
	assign_index_of_point(intfc);
	qsort(edges, num, sizeof(EDGETRI), compare_null_edges);

	printf("#num null edges %d\n", num);
	
	np = 0;
	pteq = NO;
	for(i=0; i<num-1; i++)
	{
	    e1 = &edges[i];
	    e2 = &edges[i+1];
	    
	    if(compare_null_edges((const void*)e1, (const void*)e2) == 0)
	    {
	      if(pteq)
	      {
	        printf("ERROR merge_null_edges, 3 edges are same.\n");
		clean_up(ERROR);
	      }
	      else
	      {
	        np++;
		pteq = YES;
		
		if(!check_edge_orient(e1->tri,e1->edge,e2->tri,e2->edge))
		{
		  printf("ERROR merge_null_edges, wrong orientation.\n");
		  clean_up(ERROR);
		}
		
		Tri_on_side(e1->tri,e1->edge) = e2->tri;
		Tri_on_side(e2->tri,e2->edge) = e1->tri;
	      } //if pteq
	    } //if compare_null_edges
	    else
	      pteq = NO;
	}

	printf("#merged null edges %d\n", np);

	if(debugging("nulledge"))
	{
	  POINT	*p1, *p2;

	  printf("#null edge debug\n");
	  for(i=0; i<num; i++)
	  {
	    tri = edges[i].tri;
	    j = edges[i].edge;
	    p1 = Point_of_tri(tri)[j];
	    p2 = Point_of_tri(tri)[Next_m3(j)];
	    
	    if(the_point(p1) || the_point(p2))
	    {
	      print_tri(tri, intfc);
	    }
	    
	    printf("%d %2d  %d(%6d) %d(%6d)\n", tri, j, 
	        p1, Index_of_point(p1), p2, Index_of_point(p2));
	  }
	}

	free(edges);
	
	DEBUG_LEAVE(merge_null_edges)
}

boolean	merge_buffer_interface(INTERFACE*,INTERFACE*,int);

//WARNING before calling this function, all points on intfc must 
//have point indx assigned
boolean	merge_buffer_interface(
	INTERFACE	*intfc,
	INTERFACE	*buf_intfc,
	int		index)
{
	INTERFACE	*sav_intfc;
	SURFACE		**s1,**s,*sa;
	boolean		found;

	DEBUG_ENTER(merge_buffer_interface)
	
	sav_intfc = current_interface();
	set_current_interface(intfc);

	printf("#merge buffer interface patch number %d\n",index);
	printf("#before copy surfaces\n");

	add_time_start(325);
	for(s1 = buf_intfc->surfaces; s1 && *s1; ++s1)
	{
	  //check if s has a surface with the same pos-neg comp
	  found = NO;
	  for(s = intfc->surfaces; s && *s; ++s)
	  {
	    if(negative_component(*s1) == negative_component(*s) &&
               positive_component(*s1) == positive_component(*s))
	    {
	      found = YES;
	      break;
	    }
	  }

	  //append s1 from buf_intfc to intfc
	  sa = copy_surface(*s1, NULL, NULL, YES);
	  assign_surface_point_index(sa, index);
	  Hyper_surf_index(sa) = Hyper_surf_index((*s1));
	  
	  //link tri lists for *s and s1
	  if(found)
	  {
	    //merge points
	    printf("#before merge points, surface %d %d\n", 
	          negative_component(*s), positive_component(*s));
	    
	    merge_surface_points(intfc, *s, sa);

	    last_tri(*s)->next = first_tri(sa);
	    first_tri(sa)->prev = last_tri(*s);
	    link_tri_list_to_surface(first_tri(*s),last_tri(sa),*s);

	    delete_surface(sa);

	    //replace merged points for tris
	    replace_tri_points(*s);
	  }
	}
	add_time_end(325);

	//merge tris
	//printf("#before merge tris\n");
	//add_time_start(326);
	//merge_tris(intfc);
	//add_time_end(326);

	//merge edges
	printf("#before merge null edges\n");
	add_time_start(327);
	merge_null_edges(intfc);
	add_time_end(327);

	printf("#merge finish\n");
	
	if(NO)
	{
	  char fname[128];
	  
	  //tecplot_interface_in_ball("tkmintfc", intfc);
	  
	  sprintf(fname, "merged_buff_");
	  null_sides_are_consistent();
	  check_print_intfc("After merge_buffer_interface", fname, 
	            's', intfc, 0, -1, NO);
	}

	//check_double_edge_consistent(intfc);
	
	reset_intfc_num_points(intfc);
	//reset_normal_on_intfc(intfc);
	
	set_current_interface(sav_intfc);

	DEBUG_LEAVE(merge_buffer_interface)
	
	return YES;
}

//TMP function for regriding

void	tmp_merge_patch_fronts(Front*,Front**,int*,int);

void	tmp_merge_patch_fronts(
	Front	*mfr,
	Front	**frs,
	int	*indices,
	int	nfr)
{
	INTERFACE	*mintfc, *sav_intfc;
	SURFACE		**s;
	RECT_GRID	*gr;
	int		i;
	boolean		sav_copy;

	DEBUG_ENTER(tmp_merge_patch_fronts)
	
	mintfc = mfr->interf;
	
	sav_copy = copy_intfc_states();
	sav_intfc = current_interface();

	//use frs[0] because gr should be the finest level grid
	gr = frs[0]->rect_grid;
	set_floating_point_tolerance1(gr->h);
	set_copy_intfc_states(YES);

	null_sides_are_consistent();
	check_print_intfc("Before merge","mintfc",'s',mintfc,0,-1, NO);

	strip_subdomain_bdry_curves(mintfc);
	for(s=mintfc->surfaces; s && *s; s=mintfc->surfaces)
	  delete_surface(*s);

	//only need to cut AMR_SUBDOMAIN_BOUNDARY side intfc
	printf("#before merge patch frs %d\n", nfr);

	for(i=0; i<nfr; i++)
	{
	  strip_subdomain_bdry_curves(frs[i]->interf);
	  merge_buffer_interface(mintfc,frs[i]->interf,indices[i]);
	}

	printf("#after merge patch frs.\n");
	
	ng_form_patch_subintfc_via_cut3d(mfr);
	
	{
	  char fname[128];
	  
	  sprintf(fname, "merged_all_");
	  null_sides_are_consistent();
	  check_print_intfc("After tmp_merge_patch_fronts", fname, 
	            's', mintfc, 0,-1, NO);
	}

	reset_intfc_num_points(mintfc);
	reset_normal_on_intfc(mintfc);

	set_copy_intfc_states(sav_copy);
	set_current_interface(sav_intfc);

	DEBUG_LEAVE(tmp_merge_patch_fronts)
}

char	*getIntfcRestartName(const char*,int,int);

char	*getIntfcRestartName(
	const char	*dirname,
	int		step,
	int		node)
{
	static char  intfc_name[256];

        sprintf(intfc_name,"%s/intfc.ts%s",dirname,right_flush(step,7));
        sprintf(intfc_name,"%s-nd%s",intfc_name,right_flush(node,4));
	
	return intfc_name;
}

void    TecplotFronts(char *,Front**,int);

void    TecplotFronts(
	char	*fname,
	Front	**frs,
	int	nfrs)
{

	FILE	*fp;
	int	i;
	char	bname[128];

        printf("Tec Level Fronts name %s\n", fname);
	
	fp = fopen(fname, "w");
	if(fp == NULL)
	{
	  printf("ERROR TecplotFronts, can not open %s\n", fname);
	  clean_up(ERROR);
	}
	
	fprintf(fp, "TITLE = \"TecplotFronts\"\n"
		   	    "VARIABLES = \"x\", \"y\", \"z\"\n");
	fprintf(fp, "#nfrs = %d\n", nfrs);

	for(i=0; i<nfrs; i++)
	{
	  sprintf(bname, "intfc_%02d", i);
	  //ft_tecplot_interface(fp, bname, frs[i]->interf);
	  
	  null_sides_are_consistent();
          check_print_intfc("Check consistent ", bname, 
                   's', frs[i]->interf, 0, -1, NO);
	}

	fclose(fp);
}


//all three points of tri is outside one side of bbox
EXPORT	boolean	bbox_tri_outside(
	BBOX	*bbox,
	TRI	*tri)
{
	int	dir, cnt, i;
	POINT	**p;

	p = Point_of_tri(tri);

	for(dir=0; dir<3; dir++)
	{
	    cnt = 0;
	    if(bbox->side[0][dir])
	    {
		for(i = 0; i < 3; ++i)
		{
		    if(Coords(p[i])[dir] < bbox->f[0][dir])
			cnt++;
		}
	    }
	    if(cnt == 3)
		return YES;
	    
	    cnt = 0;
	    if(bbox->side[1][dir])
	    {
		for(i = 0; i < 3; ++i)
		{
		    if(Coords(p[i])[dir] > bbox->f[1][dir])
			cnt++;
		}
	    }
	    if(cnt == 3)
		return YES;
	}

	return NO;
}

//all three points of tri are inside the box
EXPORT	boolean	bbox_tri_inside(
	BBOX	*bbox,
	TRI	*tri)
{
	int	dir, cnt, i;
	POINT	**p;

	p = Point_of_tri(tri);

	for(dir=0; dir<3; dir++)
	{
	    cnt = 0;
	    if(bbox->side[0][dir])
	    {
		for(i = 0; i < 3; ++i)
		{
		    if(Coords(p[i])[dir] < bbox->f[0][dir])
			cnt++;
		}
	    }
	    //more than one point on the left
	    if(cnt >= 1)
		return NO;
	    
	    cnt = 0;
	    if(bbox->side[1][dir])
	    {
		for(i = 0; i < 3; ++i)
		{
		    if(Coords(p[i])[dir] > bbox->f[1][dir])
			cnt++;
		}
	    }
	    //more than one point on the right
	    if(cnt >= 1)
		return NO;
	}

	return YES;
}

//all three points are inside the box or tri has intersection with box faces
EXPORT	boolean	bbox_tri_inside1(
	BBOX	*bbox,
	TRI	*tri)
{
	int	dir, cnt, i;
	POINT	**p;

	p = Point_of_tri(tri);

	for(dir=0; dir<3; dir++)
	{
	    cnt = 0;
	    if(bbox->side[0][dir])
	    {
		for(i = 0; i < 3; ++i)
		{
		    if(Coords(p[i])[dir] < bbox->f[0][dir])
			cnt++;
		}
	    }
	    //more than one point on the left
	    if(cnt == 3)
		return NO;
	    
	    cnt = 0;
	    if(bbox->side[1][dir])
	    {
		for(i = 0; i < 3; ++i)
		{
		    if(Coords(p[i])[dir] > bbox->f[1][dir])
			cnt++;
		}
	    }
	    //more than one point on the right
	    if(cnt == 3)
		return NO;
	}

	return YES;
}


EXPORT	boolean	bboxes_tri_outside(BBOX*,TRI*);

// outside means outside ALL the bboxes.
EXPORT	boolean	bboxes_tri_outside(
	BBOX	*bbox,
	TRI	*tri)
{
	while(bbox != NULL)
	{
	    if(!bbox_tri_outside(bbox, tri))
		return NO;
	    bbox = bbox->next;
	}
	return YES;
}

// inside means inside ANY ONE bboxes.
EXPORT	boolean	bboxes_tri_inside(
	BBOX	*bbox,
	TRI	*tri)
{
	while(bbox != NULL)
	{
	    if(bbox_tri_inside(bbox, tri))
		return YES;
	    bbox = bbox->next;
	}
	return NO;
}

EXPORT	boolean	bboxes_tri_inside1(
	BBOX	*bbox,
	TRI	*tri)
{
	while(bbox != NULL)
	{
	    if(bbox_tri_inside1(bbox, tri))
		return YES;
	    bbox = bbox->next;
	}
	return NO;
}

EXPORT void bboxes_intfc_cut(BBOX*,INTERFACE*,TRI_POS);
EXPORT void bboxes_intfc_cut(
	BBOX		*bbox,
	INTERFACE	*intfc,
	TRI_POS		tri_pos)
{
	TRI		*tri, *ntri;
	SURFACE 	**s;
	
	DEBUG_ENTER(bboxes_intfc_cut)

	for(s = intfc->surfaces; s && *s; ++s)
	{
	    ntri = NULL;
	    for (tri=first_tri(*s); !at_end_of_tri_list(tri,*s); tri=ntri)
	    {
		ntri = tri->next;
		if(tri_pos == TRIOUT && bboxes_tri_outside(bbox, tri))
	    	    remove_out_domain_tri(tri,*s);
		else if(tri_pos == TRIIN && bboxes_tri_inside(bbox, tri))
	    	    remove_out_domain_tri(tri,*s);
		else if(tri_pos == TRIINSIDE && bboxes_tri_inside1(bbox, tri))
	    	    remove_out_domain_tri(tri,*s);
		else if(tri_pos == TRIFLAGA && Tri_order(tri) == -1)
	    	    remove_out_domain_tri(tri,*s);
	    }
	}
	
	for(s = intfc->surfaces; s && *s; ++s)
	{
	    if (no_tris_on_surface(*s))
	    {
	    	(void) delete_surface(*s);
		--s;
	    }
	}
	
	reset_intfc_num_points(intfc);
	
	DEBUG_LEAVE(bboxes_intfc_cut)
}

EXPORT	INTERFACE  *bboxes_intfc_sect(BBOX*,INTERFACE*,TRI_POS);

// cut all tri_pos tris 
EXPORT	INTERFACE  *bboxes_intfc_sect(
	BBOX		*bbox,
	INTERFACE	*intfc,
	TRI_POS		tri_pos)
{
	INTERFACE	*sav_intfc, *tmp_intfc, *buf_intfc;
	boolean		sav_copy;

	DEBUG_ENTER(bboxes_intfc_intersection)

	sav_copy = copy_intfc_states();
	sav_intfc = current_interface();

	set_current_interface(intfc);
	
	set_size_of_intfc_state(size_of_state(intfc));
	set_copy_intfc_states(YES);
	tmp_intfc = copy_interface(intfc);

	bboxes_intfc_cut(bbox, tmp_intfc, tri_pos);

	set_size_of_intfc_state(size_of_state(intfc));
	buf_intfc = copy_interface(tmp_intfc);
	delete_interface(tmp_intfc);

	set_copy_intfc_states(sav_copy);
	set_current_interface(sav_intfc);

	DEBUG_LEAVE(bboxes_intfc_intersection)
	
	return buf_intfc;
}

void 	bbox_set(BBOX*,RECON_BOX*,GGRID*,double);
void 	bbox_set(
	BBOX		*bbox,
	RECON_BOX	*rbox,
	GGRID		*ggr,
	double		tol)
{
	int	i, buf;
        
	buf = 0;

	for(i=0; i<3; i++)
	{
          bbox->side[0][i] = YES;
          bbox->side[1][i] = YES;
          bbox->f[0][i] = rbox->fmin[i] - (tol+buf)*ggr->h[i];
          bbox->f[1][i] = rbox->fmax[i] + (tol+buf)*ggr->h[i];
	}

	bbox->prev = NULL;
	bbox->next = NULL;
}

EXPORT	INTERFACE  *bboxes_make_patch_intfc(RECON_BOX*,INTERFACE*,GGRID*);

EXPORT	INTERFACE  *bboxes_make_patch_intfc(
	RECON_BOX	*rbox,
	INTERFACE	*intfc,
	GGRID		*ggr)
{
	double		tol;
	BBOX		bbox;
	INTERFACE	*cut_intfc, *intfcs[10];

	tol = 1e-4;
	
	strip_subdomain_bdry_curves(intfc);
	
	null_sides_are_consistent();
	check_print_intfc("Before cut intfc", "intfc", 
	            's', intfc, 0, 0, NO);
	
	bbox_set(&bbox, rbox, ggr, tol);
	cut_intfc = bboxes_intfc_sect(&bbox, intfc, TRIOUT);
	
	bbox_set(&bbox, rbox, ggr, -tol);
	bboxes_intfc_cut(&bbox, intfc, TRIIN);
	//bboxes_intfc_cut(&bbox, intfc, TRIINSIDE);

	//check
	null_sides_are_consistent();
	check_print_intfc("Before merge_buffer_interface intfc", "bintfc", 
	            's', intfc, 0,-1, NO);
	
	null_sides_are_consistent();
	check_print_intfc("Before merge_buffer_interface cut_intfc", "cintfc", 
	            's', cut_intfc, 0,-1, NO);

	tecplot_interface_in_ball("tkintfc", intfc);
	tecplot_interface_in_ball("tkcintfc", cut_intfc);
	
	intfcs[0] = cut_intfc;
	tecplot_merge_surfaces("bpintfc", intfc, intfcs, 1, pp_mynode());

	assign_point_index(intfc, 0);
	merge_buffer_interface(intfc, cut_intfc, 0);
	
	null_sides_are_consistent();
	check_print_intfc("After merge_buffer_interface", "mintfc", 
	            's', intfc, 0,0, NO);
}

int	rbox_proc_index(int,RECON_BOX*);

int	rbox_proc_index(
	int		n,
	RECON_BOX	*rbox)
{
	int	i;

	for(i=0; i<rbox->np; i++)
	  if(rbox->procs[i] == n)
	    return i;
	return -1;
}

void	rbox_set_recon_grid(
	RECT_GRID	*rgr,
	RECON_BOX	*rbox,
	INTERFACE	*intfc)
{
	RECT_GRID	*gr = &topological_grid(intfc);
	int		i,gmax[3];

	for(i=0; i<3; i++)
	  gmax[i] = rbox->bmax[i] - rbox->bmin[i];
	set_rect_grid(rbox->fmin, rbox->fmax, gr->GL, gr->GU, NULL, NULL,
		      gmax, 3, &gr->Remap, rgr);
}

EXPORT void rbox_set_recon_intfc(
        INTERFACE       *intfc,
        Front           *patch_front)
{
        INTERFACE               *sav_intfc;
        boolean                    sav_copy; 

        DEBUG_ENTER(rbox_set_recon_intfc)

        printf("#set_patch_intfc bf\n");

        sav_intfc = current_interface();
        sav_copy = copy_intfc_states();
        
	set_size_of_intfc_state(patch_front->sizest);
        set_copy_intfc_states(YES);

        patch_front->interf = copy_interface(intfc);
	patch_front->interf->modified = YES;

        set_current_interface(sav_intfc);
        set_copy_intfc_states(sav_copy);

        printf("#rbox_set_recon_intfc af\n");
        
	DEBUG_LEAVE(rbox_set_recon_intfc)
}

//assume patch_front->rect_grid is initialized and bdry_side is initialized
EXPORT void rbox_set_recon_intfc_grid(
	INTERFACE	*intfc,
	RECT_GRID	*rgr)
{
        int                     i, dir;
        RECT_GRID               *tgr;

        DEBUG_ENTER(rbox_set_recon_intfc_grid)
        
        tgr = &topological_grid(intfc);
        copy_rect_grid(tgr,rgr);
        
	for(dir = 0; dir < 3; dir++)
        {
          for(i = 0; i < 2; i++)
	    rect_boundary_type(intfc,dir,i) = SUBDOMAIN_BOUNDARY;
        }

        DEBUG_LEAVE(rbox_set_recon_intfc_grid)
}

void	tecplot_box_interface(char*,INTERFACE*,RECON_BOX*,int);
EXPORT	boolean	rbox_repair_intfc(Front*,int flag[3][2]);

INTERFACE	*rbox_recon_interface(
	INTERFACE	*intfc,
	RECON_BOX	*rbox,
	Front		*fr)
{
	Front			*rfr;
        INTERFACE               *sav_intfc;
	boolean			sav_copy;

	DEBUG_ENTER(rbox_recon_interface)

	rfr = copy_front(fr);
	scalar(&(rfr->rect_grid),sizeof(RECT_GRID));
	
	rbox_set_recon_grid(rfr->rect_grid, rbox, fr->interf);
	
	rbox_set_recon_intfc(intfc, rfr);
	delete_interface(intfc);
	
	rbox_set_recon_intfc_grid(rfr->interf, rfr->rect_grid);
	
	sav_intfc = current_interface();
        sav_copy = copy_intfc_states();

	set_current_interface(rfr->interf);
	//grid_based_box_untangle check Boundary_point to determine point in boundary
	//it can be proved that the subdomain curves are outside the recon box
	install_subdomain_bdry_curves(rfr->interf);
	
	//tecplot_box_interface("rfr", rfr->interf, rbox, pp_mynode());
	
	check_print_intfc("Before rbox_recon_interface", "rbrbf", 
	            's', rfr->interf, fr->step, -1, NO);

	rbox_repair_intfc(rfr, rbox->bd_flag);
	
	strip_subdomain_bdry_curves(rfr->interf);
	intfc = rfr->interf;
	rfr->interf = NULL;

	free(rfr->rect_grid);
	free_front(rfr);
	
	set_current_interface(sav_intfc);
        set_copy_intfc_states(sav_copy);

	DEBUG_LEAVE(rbox_recon_interface)
	
	return intfc;
}

void	rbox_shift_interface(INTERFACE*,int*,RECT_GRID*);

void	rbox_shift_interface(
	INTERFACE	*intfc,
	int		*shift,
	RECT_GRID	*gr)
{
	int	j;
	double	T;

	for(j=0; j<3; j++)
	  if(shift[j] != 0)
	  {
	    T = (gr->GU[j] - gr->GL[j])*shift[j];
	    printf("#shift in dir %d\n", j);
	    shift_interface(intfc,T,j);
	  }
}

int	rbox_recv_ready(
	RECON_BOX	*rbox)
{
	int	i, proc;

	i = rbox->np;
	while(i == rbox->np)
	{
	  for(i=0; i<rbox->np; i++)
	  {
	    proc = rbox->procs[i];
	    if(proc != pp_mynode() && pp_iprobe(proc, TABLE_ID))
	      break;
	  }
	}

	return i;
}

boolean	rbox_communication_interface(RECON_BOX*,int,RECON_BOX*,GGRID*,Front*);

boolean	rbox_communication_interface(
	RECON_BOX	*sbox,
	int		nbox,
	RECON_BOX	*rbox,
	GGRID		*ggr,
	Front		*fr)
{
	RECON_BOX	*b;
	BBOX		bbox;
	INTERFACE	*intfc, *send_intfc, *recv_intfc[MAX_RBOX_PROC], *mintfc;
	RECT_GRID	*gr;
	int		i,j,ind,proc,mbox,mb;
	int		recv_procs[MAX_RBOX_PROC];
	double		tol, T;

	DEBUG_ENTER(rbox_communication_interface)

	printf("\n#rbox comm intfc nbox = %d rbox=%d\n", nbox, rbox);
	tol = 1.0e-4;

	gr = fr->rect_grid;
	intfc = fr->interf;
	strip_subdomain_bdry_curves(intfc);
	
	//cut sbox and send
	for(i=0; i<nbox; i++)
	{
	  //ex: bboxes_make_patch_intfc
	  b = &sbox[i];
	  bbox_set(&bbox, b, ggr, tol);
	  send_intfc = bboxes_intfc_sect(&bbox, intfc, TRIOUT);
	  
	  proc = b->proc;
	  if(proc != pp_mynode())
	  {
	    send_interface(send_intfc, proc);
	    printf("#send box %d to proc %d\n", i, proc);
	    fflush(NULL);
	    delete_interface(send_intfc);
	  }
	  else
	  {
	    ind = rbox_proc_index(proc, rbox);
	    recv_intfc[ind] = send_intfc;
	    mintfc = send_intfc;
	    printf("#send box %d to itself %d\n", i, proc);
	  }

	  bbox_set(&bbox, b, ggr, -tol);
	  bboxes_intfc_cut(&bbox, intfc, TRIINSIDE);
	}

	printf("#send finish\n");
	fflush(NULL);

	printf("\n#recv intfc\n");
	//recv intfc
	if(rbox != NULL)
	{
	  printf("#rbox->np=%d\n", rbox->np);
	  for(j=0; j<rbox->np-1; j++)
	  {
	    i = rbox_recv_ready(rbox);
	    proc = rbox->procs[i];
	    
	    if(proc != pp_mynode())
	    {
	      ind = rbox_proc_index(proc,rbox);
	      
	      printf("#recv box %d from proc %d bf\n", ind, proc);
	      fflush(NULL);
	      recv_intfc[ind] = receive_interface(proc);
	      printf("#recv box %d from proc %d\n", ind, proc);
	      fflush(NULL);
	      rbox_shift_interface(recv_intfc[ind],rbox->shift[i],gr);
	    }
	  }  //for i<rbox->np
	}
	
	null_sides_are_consistent();
	check_print_intfc("After rbox comm intfc", "intfc", 
	          's', intfc, fr->step, -1, NO);

	//merge recv_intfc and call recon function to recon
	if(rbox != NULL)
	{
	  int np = rbox->np;

	  for(i=0; i<np; i++)
	  {
	    null_sides_are_consistent();
	    check_print_intfc("After rbox comm recv_intfc", "intfc", 
	            's', recv_intfc[i], fr->step, -1, NO);
	  }
	  
	  //tecplot_merge_surfaces("mintfc",intfc, recv_intfc, np, pp_mynode());
	  
	  assign_point_index(mintfc, pp_mynode());

	  for(i=0; i<np; i++)
	  {
	    proc = rbox->procs[i];
	    
	    printf("#merge intfc from proc %d\n", proc);
	    
	    if(proc != rbox->proc)
	    {
	      merge_buffer_interface(mintfc, recv_intfc[i], proc);
	      delete_interface(recv_intfc[i]);
	    }
	  }
	  
	  null_sides_are_consistent();
	  check_print_intfc("After merge_buffer_interface", "mintfcd", 
	            's', mintfc, fr->step, -1, NO);

	  mintfc = rbox_recon_interface(mintfc, rbox, fr);
	  
	  null_sides_are_consistent();
	  check_print_intfc("After rbox_recon_interface", "mintfce", 
	            's', mintfc, fr->step, fr->step-1, NO);
	}
	//else
	//  tecplot_merge_surfaces("mintfc",intfc, NULL, 0, pp_mynode());

	pp_gsync();

	mbox = 0;
	printf("\n#send recon intfc\n");
	//send reconstructed intfc
	if(rbox != NULL)
	{
	  for(i=0; i<rbox->np; i++)
	  {
	    proc = rbox->procs[i];
	    if(proc == pp_mynode())
	    {
	      printf("#send intfc to itself %d\n", proc);
	      fflush(NULL);
	      recv_intfc[mbox] = mintfc;
	      recv_procs[mbox] = proc;
	      mbox++;
	    }
	    else
	    {
	      printf("#send intfc to proc %d bf\n", proc);
	      fflush(NULL);
	      send_interface(mintfc, proc);
	      printf("#send intfc to proc %d \n", proc);
	      fflush(NULL);
	    }
	  }  //for i<rbox->np
	}
	
	{
	int tmbox = mbox;
	
	pp_global_isum(&tmbox, 1);
	printf("#tmbox = %d\n", tmbox);
	}

	printf("\n#recv recon intfc\n");
	
	mb = mbox;
	//recv from sbox
	for(j=0; j<nbox-mb; j++)
	{
	  i = nbox;
	  while(i == nbox)
	  {
	    for(i=0; i<nbox; i++)
	    {
	      proc = sbox[i].proc;
	      if(proc != pp_mynode() && pp_iprobe(proc, TABLE_ID))
	        break;
	    }
	  }

	  b = &sbox[i];
	  proc = b->proc;
	  if(proc != pp_mynode())
	  {
	    printf("#recv box %d from proc %d bf\n", i, proc);
	    fflush(NULL);
	    recv_intfc[mbox] = receive_interface(proc);
	    recv_procs[mbox] = proc;
	    printf("#recv box %d from proc %d \n", i, proc);
	    fflush(NULL);
	      
	    ind = rbox_proc_index(proc,b);
	    rbox_shift_interface(recv_intfc[mbox],b->shift[ind],gr);
	    mbox++;
	  }
	}

	//tecplot_merge_surfaces("recmintfc", intfc, recv_intfc, mbox, pp_mynode());

	printf("\n#merge recon recv_intfc\n");
	assign_point_index(intfc, pp_mynode());

	for(i=0; i<mbox; i++)
	{
	  printf("#merge recv_intfc %d\n", i);
	  merge_buffer_interface(intfc, recv_intfc[i], recv_procs[i]);
	  delete_interface(recv_intfc[i]);
	}
	
	null_sides_are_consistent();
	check_print_intfc("After merge recv_intfc", "mrintfc", 
	          's', intfc, fr->step, fr->step-1, NO);

	printf("#rbox comm intfc finish\n\n");
	
	DEBUG_LEAVE(rbox_communication_interface)
	
	return YES;
}



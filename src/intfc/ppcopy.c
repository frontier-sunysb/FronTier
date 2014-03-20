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
*				ppcopy.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
************************************************************************
*
*NAME
*	ppcopy.c - 	broadcast an INTERFACE onto all the participating 
*			processors from one master processor in 
*			a daisy chain fashion.
*
*SYNOPSIS
*	This module contains the major functions:
*	(1)	i_send_interface(intfc,dst_id)
*		send an interface from processor src_id to processor dst_id;
*
*	(2)	i_receive_interface(src_id)
*		 receive an interface from processor src_id to processor dst_id;
*
*
*DESCRIPTION
* 
*	Broadcast is a log(n) operation.
*
* 			I.	GENERAL IDEAS:
*
*	The following is a tree, log(n) depth, operation: 
*	It broadcast a copy of interface from processor 0 to the rest.
*
*			II. 	CONVENTIONS:  
*
*	Any variable explicitly marked with prefix src_ must be on 
*	the source processor, eg, src_intfc is a source variable.
*	While the counterpart will have conventional name. The following 
*	list shows some more examples:
*		
*		on source processor:		on destination processor:
*			src_intfc     <------>	intfc
*			src_chunk     <------>	chunk
*			src_id        <------>	id
*
*				... ...
*
*		nn  = nearest neighbors
*		nnn = nearest neighboring node
*		nns = nearest neighboring subdomains
*
*
* 			III. CAUTION
*
*	Under development.
*
*
*OPTIONS
*	
*	There are two possible schemes to do this initial copy and the
*	timing depends on the processor numbers.
*	(1) First, using one processor to do the decomposition,
*	    then, copying the decomposed subdomains onto all processors
*	    the total timing for this is (n = total processors)
*		t = t0 * n + (t1/n)*n = t0 * n + t1;
*	    where t0 = decomposition time (for each subdomain)
*	    where t1 = copy time (for each subdomain)
*
*	(2) First, copy the original interface to all processors, 
*	    then, decompose the original interface on all processors
*	    in PARALLEL.
*	    the total timing for this is 
*		t = t0 + t1* ln(n);
*
*	The belief is t0 >> t1 and n is not very small, thus
*	the second scheme is adopted.
*
* 
*ENVIRONMENT
*	Parallel computing with more than one processor. This code
*	ought to work in an unprocessor environment, but in a very
*	subtle way. Cautions should be taken if taking the orginal
*	code to unprocessor environment since (1) it may be slower 
*	than regular sequential code; (2) it may break in very subtle ways.
*	
*SEE ALSO
*	files: 
*		iscatter.c, int.h
* 
*HISTORY
*
*************************************************************************/



#define	DEBUG_STRING	"iscatter"
#include <intfc/iloc.h>

int	tag_shf = 100;

void	set_tag_shift(int);

void	set_tag_shift(
	int	shf)
{
	tag_shf = shf;
}

int	get_tag_shift()
{
	return tag_shf;
}

EXPORT	void i_send_interface(
	INTERFACE	*intfc,
	int		dst_id)
{
	POINTER		*top_addr;
	struct Chunk 	*chunk;
	ALIGN		*top;
	int		i;
	int		nchunks;
	CURVE		**c;

	DEBUG_ENTER(i_send_interface)
	if (exists_interface(intfc) != YES)
	{
	    (void) printf("WARNING in i_send_interface(), "
	                  "NO intfc to be sent.\n");
	    DEBUG_LEAVE(i_send_interface)
	    return;
	}

	/*
	 ensure that sort status of points is set to NO prior to send
	 upon receive it is assumed this flag is already set. One proc
	 set sorted(p), another proc will check it in 
	 	g_reconstruct_point_pointers 
	*/
	if (intfc->dim == 3)
	{
	    (void) next_point(intfc,NULL,NULL,NULL);
	    /*  make sure the curves are consistent. 
	        g_reconstruct_bond_pointers
	        f_reconstruct_bond_pointers */
	    for(c=intfc->curves; c && *c; c++)
	        reorder_curve_link_list(*c);
	}

	/* send out the interface table */

	pp_send(TABLE_ID+tag_shf,(POINTER) intfc->table,
			sizeof(struct Table),dst_id);

	/* count the number of chunks and send this number out */

	nchunks = intfc->table->num_chunks;

	/* send an array with the old chunk addresses */

	uni_array(&top_addr,nchunks,sizeof(POINTER));

	for (i = 0, chunk = intfc->table->first_chunk; chunk; 
			chunk = chunk->next, i++)
	    top_addr[i] = (POINTER) ChunkTop(chunk);

	pp_send(CHUNK_ADDR_ID+tag_shf,(POINTER)top_addr,
			nchunks*sizeof(POINTER),dst_id);
	free(top_addr);
		
	/* send out the chunks */

	for (i = 0, chunk = intfc->table->first_chunk; chunk;
						chunk = chunk->next, i++)
	{
	    top = ChunkTop(chunk);
	    pp_send(chunk_id(i)+tag_shf,(POINTER)top,ChunkSize(intfc),dst_id);
	}
	DEBUG_LEAVE(i_send_interface)
}		/*end i_send_interface*/

EXPORT	INTERFACE *i_receive_interface(
	int		src_id)
{
	INTERFACE	*new_intfc;
	struct Table 	*new_table, *tmp_table;
	ALIGN		*top;
	POINTER		*ncaddr, *ocaddr;
	int		nchunks;
	int		i, dim;

	DEBUG_ENTER(i_receive_interface)

	scalar(&tmp_table,sizeof(struct Table));

	add_time_start(350);
	pp_recv(TABLE_ID+tag_shf,src_id,(POINTER)tmp_table,sizeof(struct Table)) ;

	dim = tmp_table->rect_grid.dim;
	new_intfc = make_interface(dim);
	new_table = new_intfc->table;
	nchunks = tmp_table->num_chunks;
	uni_array(&ocaddr,nchunks,sizeof(POINTER));
	uni_array(&ncaddr,nchunks,sizeof(POINTER));
	pp_recv(CHUNK_ADDR_ID+tag_shf,src_id,(POINTER)ocaddr,
			nchunks*sizeof(POINTER));

	for (i = 0; i < nchunks; i++)
	{
	    top = ChunkTop(new_chunk(new_table));
	    pp_recv(chunk_id(i)+tag_shf,src_id,(POINTER)top,
			ChunkSize(new_intfc));
	    ncaddr[i] = (POINTER) top;
	}
	
	add_time_end(350);

	add_time_start(351);

	reconstruct_interface_pointers(new_intfc,tmp_table,ocaddr,ncaddr);

	add_time_end(351);
	
	free_these(3,tmp_table,ocaddr,ncaddr);

	DEBUG_LEAVE(i_receive_interface)
	return new_intfc;
}		/*end i_receive_interface*/

/*
*		i_reconstruct_interface_pointers():
*
*	Given an interface obtained by a direct memory copy of
*	the interface associated with otbl,  the functions
*	reconstructs  all of the addresses associated with the
*	interface data structures.  The input arrays ocad
*	and ncad contain the addresses of the starting positions
*	of the chunks associated with the two interface tables.
*	It is assumed that the chunks for nintfc have already been
*	allocated and are direct binary copies of the chunks from
*	the old table.
*/

void	check_curve_states(INTERFACE *);

EXPORT	void i_reconstruct_interface_pointers(
	INTERFACE	*nintfc,
	struct Table	*otbl,
	POINTER		*ocad,
	POINTER		*ncad)
{
	INTERFACE	*ointfc;
	struct Table	*ntbl = nintfc->table;
	int		nchks = otbl->num_chunks;
	POINT		**pp;
	CURVE		**c;
	NODE		**n;
	SURFACE		**s;

	nintfc->modified = YES;

	ointfc = (INTERFACE *)new_address(nintfc,
					  otbl->interface,ocad,ncad,nchks);
	otbl->interface = ointfc;

	/* Copy top level data from tbl to ntbl */

	copy_rect_grid(&ntbl->rect_grid,&otbl->rect_grid);
	ntbl->fixed_grid = otbl->fixed_grid;
	ntbl->new_grid = otbl->new_grid;
	ntbl->max_comp = otbl->max_comp;
	ntbl->min_comp = otbl->min_comp;
	ntbl->ext_comp = otbl->ext_comp;
	ntbl->remainder = otbl->remainder;
	ntbl->num_chunks = nchks;
	ntbl->top = (ALIGN *)new_address(nintfc,otbl->top,ocad,ncad,nchks);

	/* Reset top level in the interface pointers */

	nintfc->num_points = ointfc->num_points;
	nintfc->points = (POINT **)new_address(nintfc,
					       ointfc->points,ocad,ncad,nchks);
	nintfc->nodes = (NODE **)new_address(nintfc,
					     ointfc->nodes,ocad,ncad,nchks);
	nintfc->curves = (CURVE **)new_address(nintfc,
					       ointfc->curves,ocad,ncad,nchks);
	nintfc->surfaces =
		(SURFACE **)new_address(nintfc,
					ointfc->surfaces,ocad,ncad,nchks);

	/* Reset addresses in top level arrays */

	/*for 3d code nintfc->points = NULL */
	for (pp = nintfc->points; pp && *pp; pp++)
	{
		*pp = (POINT *)new_address(nintfc,
					   *pp,ocad,ncad,nchks);
		reconstruct_point_pointers(*pp,nintfc,ointfc,ocad,ncad,nchks);
	}

	for (n = nintfc->nodes; n && *n; n++)
	{
		*n = (NODE *)new_address(nintfc,*n,ocad,ncad,nchks);
		reconstruct_node_pointers(*n,nintfc,ointfc,ocad,ncad,nchks);
	}
	for (c = nintfc->curves; c && *c; c++)
	{
		*c = (CURVE *)new_address(nintfc,*c,ocad,ncad,nchks);
		reconstruct_curve_pointers(*c,nintfc,ointfc,ocad,ncad,nchks);
	}

	hyper_surf_list(nintfc) = NULL;
	for (s = nintfc->surfaces; s && *s; s++)
	{
		*s = (SURFACE *)new_address(nintfc,*s,ocad,ncad,nchks);
		reconstruct_surface_pointers(*s,nintfc,ointfc,ocad,ncad,nchks);
	}
	
}		/*end i_reconstruct_interface_pointers*/

/*ARGSUSED*/
EXPORT	void i_reconstruct_point_pointers(
	POINT		*p,
	INTERFACE	*nintfc,
	INTERFACE	*ointfc,
	POINTER		*ocad,
	POINTER		*ncad,
	int		nchks)
{
	if (nintfc->dim == 1)
	{
		p->interface = nintfc;
		p->obj = p;
		Hyper_surf(p) = (HYPER_SURF *)new_address(nintfc,Hyper_surf(p),
						          ocad,ncad,nchks);
	}
}		/*end i_reconstruct_point_pointers*/

/*ARGSUSED*/
EXPORT	void i_reconstruct_node_pointers(
	NODE		*n,
	INTERFACE	*nintfc,
	INTERFACE	*ointfc,
	POINTER		*ocad,
	POINTER		*ncad,
	int		nchks)
{
	CURVE		**c;

	n->interface = nintfc;
	n->obj = n;
	n->posn = (POINT *)new_address(nintfc,n->posn,ocad,ncad,nchks);
	reconstruct_point_pointers(n->posn,nintfc,ointfc,ocad,ncad,nchks);
	n->in_curves = (CURVE **)new_address(nintfc,n->in_curves,
					     ocad,ncad,nchks);
	for (c = n->in_curves; c && *c; c++)
	    *c = (CURVE *)new_address(nintfc,*c,ocad,ncad,nchks);
	n->out_curves = (CURVE **)new_address(nintfc,n->out_curves,
					      ocad,ncad,nchks);
	for (c = n->out_curves; c && *c; c++)
	    *c = (CURVE *)new_address(nintfc,*c,ocad,ncad,nchks);
	if (nintfc->dim == 2)
	{
	    HYPER_SURF_BDRY *hsb;
	    hsb = Hyper_surf_bdry(n);
	    Hyper_surf_bdry(n) = hsb = 
		(HYPER_SURF_BDRY *)new_address(nintfc,hsb,ocad,ncad,nchks);
	    hsb->interface = nintfc;
	    Node_of_hsb(hsb) = n;
	    Hyper_surf_bdry(hsb) = hsb;
	}
}		/*end i_reconstruct_node_pointers*/

/*ARGSUSED*/
EXPORT	void i_reconstruct_bond_pointers(
	BOND		*b,
	INTERFACE	*nintfc,
	INTERFACE	*ointfc,
	POINTER		*ocad,
	POINTER		*ncad,
	int		nchks)
{
	BOND_TRI **btris;

	b->next = (BOND *)new_address(nintfc,b->next,ocad,ncad,nchks);
	b->prev = (BOND *)new_address(nintfc,b->prev,ocad,ncad,nchks);
	b->start = (POINT *)new_address(nintfc,b->start,ocad,ncad,nchks);
	b->end = (POINT *)new_address(nintfc,b->end,ocad,ncad,nchks);

	
	Btris(b) = (BOND_TRI **)new_address(nintfc,Btris(b),ocad,ncad,nchks);
	for (btris = Btris(b); btris && *btris; btris++)
	{
		*btris = (BOND_TRI *)new_address(nintfc,*btris,
						 ocad,ncad,nchks);
		(*btris)->tri =
			(TRI *)new_address(nintfc,(*btris)->tri,
	     				   ocad,ncad,nchks);
	        (*btris)->bond =
                        (BOND *)new_address(nintfc,(*btris)->bond,
                                           ocad,ncad,nchks); 	
                  
                (*btris)->curve =
                        (CURVE *)new_address(nintfc,(*btris)->curve,
                                           ocad,ncad,nchks);
                 
                (*btris)->surface =
                        (SURFACE *)new_address(nintfc,(*btris)->surface,
                                           ocad,ncad,nchks);
	}
}		/*end i_reconstruct_bond_pointers*/

EXPORT	void i_reconstruct_curve_pointers(
	CURVE		*c,
	INTERFACE	*nintfc,
	INTERFACE	*ointfc,
	POINTER		*ocad,
	POINTER		*ncad,
	int		nchks)
{
	HYPER_SURF      *hs;
	HYPER_SURF_BDRY *hsb;
	BOND		*b;

	c->interface = nintfc;
	c->obj = c;
	c->start = (NODE *)new_address(nintfc,c->start,ocad,ncad,nchks);
	c->end = (NODE *)new_address(nintfc,c->end,ocad,ncad,nchks);
	hs = Hyper_surf(c);
	Hyper_surf(c) = hs =
	    (HYPER_SURF *)new_address(nintfc,hs,ocad,ncad,nchks);
	hsb = Hyper_surf_bdry(c);
	Hyper_surf_bdry(c) = hsb =
	    (HYPER_SURF_BDRY *)new_address(nintfc,hsb,ocad,ncad,nchks);
	switch (nintfc->dim)
	{
	case 1:
	    break;
	case 2:
	    hs->interface = nintfc;
	    Curve_of_hs(hs) = c;
	    Hyper_surf(hs) = hs;
	    c->pos_surfaces = NULL;
	    c->neg_surfaces = NULL;
	    break;
	case 3:
	{
	    SURFACE		**s;

	    hsb->interface = nintfc;
	    Curve_of_hsb(hsb) = c;
	    Hyper_surf_bdry(hsb) = hsb;
	    c->pos_surfaces = (SURFACE **)new_address(nintfc,c->pos_surfaces,
						      ocad,ncad,nchks);
	    for (s = c->pos_surfaces; s && *s; s++)
	    	*s = (SURFACE *)new_address(nintfc,*s,ocad,ncad,nchks);
	    c->neg_surfaces = (SURFACE **)new_address(nintfc,c->neg_surfaces,
						      ocad,ncad,nchks);
	    for (s = c->neg_surfaces; s && *s; s++)
	    	*s = (SURFACE *)new_address(nintfc,*s,ocad,ncad,nchks);
	    break;
	}
	}
	
	c->first = (BOND *)new_address(nintfc,c->first,ocad,ncad,nchks);
	c->last = (BOND *)new_address(nintfc,c->last,ocad,ncad,nchks);
	for (b = c->first; b != NULL; b = b->next)
	{
		reconstruct_bond_pointers(b,nintfc,ointfc,ocad,ncad,nchks);
		if (b != c->last)
			reconstruct_point_pointers(b->end,nintfc,ointfc,ocad,
				ncad,nchks);
	}
}		/*end i_reconstruct_curve_pointers*/

/* ARGSUSED */
LIB_LOCAL	void i_reconstruct_surface_pointers(
	SURFACE		*s,
	INTERFACE	*nintfc,
	INTERFACE	*ointfc,
	POINTER		*ocad,
	POINTER		*ncad,
	int		nchks)
{
	TRI	   *tri;
	CURVE	   **c;
	HYPER_SURF *hs;

	s->interface = nintfc;
	s->obj = s;
	hs = Hyper_surf(s);
	Hyper_surf(s) = hs =
	    (HYPER_SURF *)new_address(nintfc,hs,ocad,ncad,nchks);
	if (!add_to_pointers(hs,&hyper_surf_list(nintfc)))
	{
	    printf("ERROR in i_reconstruct_surface_pointers, add_to_pointers fails.\n");
	    clean_up(ERROR);
	}
	
	Hyper_surf(hs) = hs;
	hs->interface = nintfc;
	Surface_of_hs(hs) = s;
	s->pos_curves = (CURVE **)new_address(nintfc,s->pos_curves,
					      ocad,ncad,nchks);
	for (c = s->pos_curves; c && *c; c++)
	    *c = (CURVE *)new_address(nintfc,*c,ocad,ncad,nchks);
	s->neg_curves = (CURVE **)new_address(nintfc,s->neg_curves,
					      ocad,ncad,nchks);
	for (c = s->neg_curves; c && *c; c++)
	    *c = (CURVE *)new_address(nintfc,*c,ocad,ncad,nchks);
	first_tri(s) = (TRI *)new_address(nintfc,first_tri(s),ocad,ncad,nchks);
	last_tri(s) = (TRI *)new_address(nintfc,last_tri(s),ocad,ncad,nchks);
	for (tri = first_tri(s); !at_end_of_tri_list(tri,s); tri = tri->next)
		reconstruct_tri_pointers(tri,nintfc,ointfc,ocad,ncad,nchks);
}		/*end i_reconstruct_surface_pointers*/

/*ARGSUSED*/
LIB_LOCAL	void i_reconstruct_tri_pointers(
	TRI		*t,
	INTERFACE	*nintfc,
	INTERFACE	*ointfc,
	POINTER		*ocad,
	POINTER		*ncad,
	int		nchks)
{
	Point_of_tri(t)[0] =
	    (POINT *)new_address(nintfc,Point_of_tri(t)[0],ocad,ncad,nchks);
	Point_of_tri(t)[1] =
	    (POINT *)new_address(nintfc,Point_of_tri(t)[1],ocad,ncad,nchks);
	Point_of_tri(t)[2] =
	    (POINT *)new_address(nintfc,Point_of_tri(t)[2],ocad,ncad,nchks);
	reconstruct_point_pointers(Point_of_tri(t)[0],
				   nintfc,ointfc,ocad,ncad,nchks);
	reconstruct_point_pointers(Point_of_tri(t)[1],
				   nintfc,ointfc,ocad,ncad,nchks);
	reconstruct_point_pointers(Point_of_tri(t)[2],
				   nintfc,ointfc,ocad,ncad,nchks);
	Neighbor_on_side01(t) = (POINTER)new_address(nintfc,
						     Neighbor_on_side01(t),
						     ocad,ncad,nchks);
	Neighbor_on_side12(t) = (POINTER)new_address(nintfc,
						     Neighbor_on_side12(t),
						     ocad,ncad,nchks);
	Neighbor_on_side20(t) = (POINTER)new_address(nintfc,
						     Neighbor_on_side20(t),
						     ocad,ncad,nchks);
	t->prev = (TRI *)new_address(nintfc,t->prev,ocad,ncad,nchks);
	t->next = (TRI *)new_address(nintfc,t->next,ocad,ncad,nchks);
	t->surf = (SURFACE *)new_address(nintfc,t->surf,ocad,ncad,nchks);
}		/*end i_reconstruct_tri_pointers*/

/*	given the address in other processor oaddr and the old trunk 
 *	starting address ocad, the new trunk starting address ncad. 
 *	_new_address will return the address of ocad in the current processor.
*/

EXPORT	POINTER _new_address(
	INTERFACE	*nintfc,
	POINTER		oaddr,
	POINTER		*ocad,
	POINTER		*ncad,
	int		nchks)
{
	byte		*naddr;
	long		rel_addr;
	int		i;

	if (oaddr == NULL) return NULL;

	for (i = 0; i < nchks; i++)
	{
		rel_addr = ((byte *) oaddr) - ((byte *) ocad[i]);
		if ((0 <= rel_addr) && (rel_addr < ChunkSize(nintfc)))
		{
			naddr = ((byte *) ncad[i]) + rel_addr;
			return (POINTER) naddr;
		}
	}
	return NULL;
}		/*end _new_address*/

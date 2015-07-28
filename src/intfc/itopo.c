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

#include <intfc/iloc.h>

LOCAL void mark_topolinked_tris(TRI*,TRI***,int*);
LOCAL SURFACE *merge_surfs_exist(INTERFACE*);

/**********************************************************************
*	This function start from a seed tri and traversal through     *
*	topological links to mark all connected tris on a surface.    *
*	Prior to this call, relevant tris' sort status must be NO.    *
**********************************************************************/

EXPORT void merge_surfaces(
	INTERFACE *intfc)
{
	SURFACE **s,*surf;
	TRI *tri,**tri_linked_list;
	int num_tri;
	boolean to_merge = NO;
	boolean preserve_neighbors = YES;

	printf("Entering merge_surfaces()\n");
	intfc_surface_loop(intfc,s)
	{
	    surf_tri_loop(*s,tri)
		sorted(tri) = NO;
	}

	while ((surf = merge_surfs_exist(intfc)) != NULL)
	{
	    to_merge = YES;
	    surf_tri_loop(surf,tri)
	    {
		if (sorted(tri) == NO) break;
	    }
	    num_tri = surf->num_tri;
	    mark_topolinked_tris(tri,&tri_linked_list,&num_tri);
	    printf("surf->num_tri = %d linked num_tri = %d\n",
			surf->num_tri,num_tri);
	    surf_tri_loop(surf,tri)
	    {
		if (sorted(tri) == NO)
		    remove_tri_from_surface(tri,surf,preserve_neighbors);
	    }
	    for (s = surf->merge_surfs; s && *s; ++s)
		delete_from_pointers(surf,&(*s)->merge_surfs);
	    surf->merge_surfs = NULL;
	    printf("Step 1: surf->num_tri = %d\n",surf->num_tri);
	    
	    break;
	}
	if (to_merge)
	{
	    printf("Merge surface needed!\n");
	    clean_up(0);
	}
}	/* end merge_surfaces */

LOCAL void mark_topolinked_tris(
	TRI *tri,
	TRI ***tri_linked_list,
	int *num_tri)
{
	TRI *t,*nbt;
	POINTER_Q *seed_queue,*pq;
	int i,nt = 0;
	TRI **tri_list;

	set_pointer_queue_opts(PQ_BLOCK_SIZE,200,PQ_ALLOC_TYPE,"vmalloc",
                                PQ_ALLOC_SIZE_FOR_POINTERS,sizeof(TRI*),0);

	if ((*num_tri) > 0)	/* initial number given */
	    uni_array(&tri_list,(*num_tri),sizeof(TRI*));
	t = tri;
	seed_queue = add_to_pointer_queue(NULL,NULL);
	seed_queue->pointer = (POINTER)t;
	sorted(t) = YES;
	tri_list[nt] = t;
	nt++;
	while (seed_queue)
	{
	    pq = head_of_pointer_queue(seed_queue);
	    t = (TRI*)(pq->pointer);
	    for (i = 0; i < 3; ++i)
	    {
		if (is_side_bdry(t,i) || Tri_on_side(t,i) == NULL)
		    continue;
		nbt = Tri_on_side(t,i);
		if (sorted(nbt)) continue;
		seed_queue = add_to_pointer_queue(NULL,seed_queue);
		seed_queue->pointer = (POINTER)nbt;
		sorted(nbt) = YES;
		if (nt >= (*num_tri))
		{
		    free_these(1,tri_list);
		    (*num_tri) += 2;
	    	    uni_array(tri_list,(*num_tri),sizeof(TRI*));
		}
		tri_list[nt] = t;
		nt++;
	    }
	    seed_queue = delete_from_pointer_queue(pq);
	}
	*num_tri = nt;
	*tri_linked_list = tri_list;
}	/* end mark_topolinked_tris */

LOCAL SURFACE *merge_surfs_exist(
	INTERFACE *intfc)
{
	SURFACE **s;
	intfc_surface_loop(intfc,s)
	{
	    if ((*s)->merge_surfs != NULL)
		return *s;
	}
	return NULL;
}	/* end merge_surfs_exist */

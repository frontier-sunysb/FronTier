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
*				fadv.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*/

#define DEBUG_STRING    "fadv"
#include <front/fdecs.h>		/* includes int.h, table.h */

	/* LOCAL Function Declarations */
LOCAL	boolean	BothSidesActive(HYPER_SURF*,Front*);

LOCAL	int	advance_front1d(double,double*,Front*,Front**,POINTER);
LOCAL	int	advance_front2d(double,double*,Front*,Front**,POINTER);
LOCAL	void	EnforceFlowSpecifedStates1d(Front*);
LOCAL	void	EnforceFlowSpecifedStates2d(Front*);
LOCAL	int	delete_phys_curves_at_old_node(NODE*,int);
LOCAL	int	node_modify_time_step(NODE*,Front*,double*,int);
LOCAL	void	check_bond_lengths(INTERFACE*);
LOCAL	void 	check_2d_orientation(INTERFACE*);

EXPORT	void	set_advance_front(
	INIT_DATA	*init,
	Front		*front)
{
	int dim = front->rect_grid->dim;
	boolean efss,mvfm;

	front->_EnforceFlowSpecifedStates = NULL;
	efss = ((init==NULL) || (enforce_flow_specified_states(init)==YES)) ?
		YES : NO;
	front->movingframe = ((init==NULL) || (movingframe(init)==YES)) 
				? YES : NO;
	switch (dim)
	{
	case 1:
	    front->_pre_advance_front = NULL;
	    front->_advance_front = advance_front1d;
	    if (efss == YES)
	        front->_EnforceFlowSpecifedStates = EnforceFlowSpecifedStates1d;
	    break;
	case 2:
	    front->_advance_front = advance_front2d;
	    front->_pre_advance_front = FrontPreAdvance;
	    if (efss == YES)
	        front->_EnforceFlowSpecifedStates = EnforceFlowSpecifedStates2d;
	    break;
	case 3:
	    front->_pre_advance_front = NULL;
	    front->_advance_front = advance_front3d_tracking_control;
	    if (efss == YES)
	        front->_EnforceFlowSpecifedStates = EnforceFlowSpecifedStates3d;
	    break;
	}
}		/*end set_advance_front*/


EXPORT	int	return_advance_front(
	Front		*front,
	Front		**newfront,
	int		status,
	const char	*fname)
{
	if (front->pp_grid && front->pp_grid->nn > 1)
 	    status = syncronize_time_step_status(status,front->pp_grid);

	if (debugging("final_front"))
	{
	    print_Front_structure(front);
	    print_Front_structure(*newfront);
	}
	if (status != GOOD_STEP)
	{
	    free_front(*newfront);
	    *newfront = NULL;
	}
	else if (front->_EnforceFlowSpecifedStates != NULL)
	{
	    (*front->_EnforceFlowSpecifedStates)(*newfront);
	}
	if (debugging("trace"))
        {
            int dim = front->rect_grid->dim;
            (void) printf("Maximum propagated scaled distance = %f\n",
                        *(front->max_scaled_propagation));
            print_general_vector("Max propagated point: ",
                        front->max_prop_point,dim,"\n");
        }
	debug_front("final_front","after EnforceFlowSpecifedStates():",
	           *newfront);
	debug_print("front","Left %s\n",fname);
	return status;
}		/*end return_advance_front*/

LOCAL	void	EnforceFlowSpecifedStates1d(
	Front	*fr)
{
	INTERFACE	*intfc;
	HYPER_SURF	*hs;
	POINT		**p;
	Locstate	sl, sr;

	if ((fr==NULL) || (Fsr_list(fr)==NULL) || (fr->rect_grid->dim!=1))
	    return;

	intfc = fr->interf;
	for (p = intfc->points; p && *p; ++p)
	{
	    hs = Hyper_surf(*p);
	    if (BothSidesActive(hs,fr) == NO)
	    {
	       	slsr(*p,NULL,Hyper_surf(*p),&sl,&sr);
	    	(void) RegionIsFlowSpecified(sl,NULL,Coords(*p),
	                                 negative_component(hs),NO_COMP,fr);
	    	(void) RegionIsFlowSpecified(sr,NULL,Coords(*p),
	                                 positive_component(hs),NO_COMP,fr);
	    }
	}
}		/*end EnforceFlowSpecifedStates1d*/

LOCAL	void	EnforceFlowSpecifedStates2d(
	Front	*fr)
{
	INTERFACE	*intfc;
	CURVE		**c;
	POINT		*p;
	HYPER_SURF	*hs;
	BOND		*b;
	Locstate	sl, sr;

	if ((fr==NULL) || (Fsr_list(fr)==NULL) || (fr->rect_grid->dim!=2))
	    return;

	intfc = fr->interf;
	for (c = intfc->curves; c && *c; ++c)
	{
	    hs = Hyper_surf(*c);
	    if (is_subdomain_boundary(hs) || is_passive_boundary(hs))
	        continue;
	    if (BothSidesActive(hs,fr) == NO)
	    {
	    	b = (*c)->first;
	    	p = b->start;
	       	slsr(p,Hyper_surf_element(b),hs,&sl,&sr);
	        (void) RegionIsFlowSpecified(sl,NULL,Coords(p),
	                                 negative_component(hs),NO_COMP,fr);
	    	(void) RegionIsFlowSpecified(sr,NULL,Coords(p),
	                                 positive_component(hs),NO_COMP,fr);
		for (; b != NULL; b = b->next)
		{
		    p = b->end;
	    	    slsr(p,Hyper_surf_element(b),hs,&sl,&sr);
	    	    (void) RegionIsFlowSpecified(sl,NULL,Coords(p),
	                                 negative_component(hs),NO_COMP,fr);
	    	    (void) RegionIsFlowSpecified(sr,NULL,Coords(p),
	                                 positive_component(hs),NO_COMP,fr);
		}
	    }
	}
}		/*end EnforceFlowSpecifedStates2d*/

LOCAL	boolean	BothSidesActive(
	HYPER_SURF	*hs,
	Front		*front)
{
	FlowSpecifiedRegion     *fsr;
	INTERFACE	*intfc = front->interf;
	COMPONENT	pcomp = positive_component(hs);
	COMPONENT	ncomp = negative_component(hs);

	if (hs == NULL)
	    return NO;

	for (fsr = Fsr_list(front); fsr != NULL; fsr = fsr->next)
	{
	    if (equivalent_comps(fsr->comp,pcomp,intfc))
	    	return NO;
	    if (equivalent_comps(fsr->comp,ncomp,intfc))
	    	return NO;
	}
	return YES;
}		/*end BothSidesActive*/


EXPORT	INTERFACE	*pp_copy_interface(
	INTERFACE	*intfc)
{
	INTERFACE	*new_intfc;
	boolean		stat;
	boolean		delete_status;

	DEBUG_ENTER(pp_copy_interface)

	new_intfc = copy_interface(intfc);
	stat = (new_intfc != NULL) ? YES : NO;

	if (stat == NO)
	{
	    (void) printf("WARNING in pp_copy_interface(), "
		          "unable to copy interface");
	    if (pp_numnodes() > 1)
		(void) printf(" on processor %d\n",pp_mynode());
	    else
		(void) printf("\n");
	}

	delete_status = YES;
	if (pp_min_status(stat) == NO)
	{
	    if (stat == YES)
	    {
	    	(void) printf("WARNING in pp_copy_interface(), "
		              "unable to copy interface "
		              "on a remote processor\n");
	        delete_status = (delete_interface(new_intfc)) ? YES : NO;
		if (delete_status == NO)
		{
		    screen("ERROR in pp_copy_interface() "
			   "unable to delete interface ");
	            if (pp_numnodes() > 1)
		        screen(" on processor %d\n",pp_mynode());
	            else
		        screen("\n");
		}
	        new_intfc = NULL;
	    }
	}
	if (pp_min_status(delete_status) == NO)
	{
	    if (delete_status == YES)
	    {
	        screen("ERROR in pp_copy_interface(), unable to delete "
		       "interface on a remote processor\n");
	    }
	    clean_up(ERROR);
	}
	DEBUG_LEAVE(pp_copy_interface)
	return new_intfc;
}	/*end pp_copy_interface*/


/*ARGSUSED**/
EXPORT	void set_propagation_limits(
	Front		*front,
	Front		*newfront)
{
	HYPER_SURF_ELEMENT *hse;
	HYPER_SURF         *hs;
	POINT	           *p;

	(void) next_point(newfront->interf,NULL,NULL,NULL);
	while (next_point(newfront->interf,&p,&hse,&hs))
	{
	    n_pt_propagated(p) = NO;
	    t_pt_propagated(p) = NO;
	}
	if (front->interf->dim == 2)
	{
	    NODE **n;
	    CURVE **c;
	    int  i;
	    for (n = front->interf->nodes; *n; ++n) 
	        propagation_status(*n) = PROPAGATED_NODE;
	    for (n = newfront->interf->nodes; *n; ++n) 
	    {
	        for (i = 0; i < 2; ++i)
	            Node_vel(*n)[i] = 0.0;
	        propagation_status(*n) = UNPROPAGATED_NODE;
	    }
	    for (c = newfront->interf->curves; c && *c; ++c)
	    {
		if (!is_closed_curve(*c)) 
		    (*c)->orientation = 0;
		else
		    (*c)->orientation = (area_of_closed_curve(*c) > 0.0) ?
						1 : -1;
	    }
	}
	initialize_max_front_speed(front);
}		/*end set_propagation_limits*/

EXPORT	void copy_hypersurface_flags(
	INTERFACE	*intfc)
{
	HYPER_SURF	**hs;

	for (hs = hyper_surf_list(intfc); hs && *hs; ++hs)
	{
	    if (correspond_hyper_surf(*hs) != NULL)
	    	Hs_flag(*hs) = Hs_flag(correspond_hyper_surf(*hs));
	}
}		/*end copy_hypersurface_flags*/



/*ARGSUSED*/
LOCAL int advance_front1d(
	double		dt,
	double		*dt_frac,
	Front		*front,
	Front		**newfront,
	POINTER		wave)
{
	POINT              *oldp, *newp;
	HYPER_SURF_ELEMENT *oldhse, *newhse;
	HYPER_SURF         *oldhs, *newhs;
	INTERFACE          *intfc_old, *intfc_new;
	int		   status;
	double              V[MAXD];
	boolean		   has_tracked_points;
	static const char	   *fname = "advance_front1d";

	debug_print("front","Entered %s(step %d time %g dt %g)\n",fname,
				front->step,front->time,dt);
	debug_front("old_front","into advance front",front);

	*newfront = copy_front(front);
	Interface_redistributed(*newfront) = NO;

	has_tracked_points = (front->interf->points != NULL) ? YES : NO;
	if (pp_max_status(has_tracked_points) == NO)
	{
	    set_size_of_intfc_state(size_of_state(front->interf));
	    set_copy_intfc_states(YES);
	    (*newfront)->interf = pp_copy_interface(front->interf);
	    status = ((*newfront)->interf != NULL) ? GOOD_STEP : ERROR_IN_STEP;
	    return return_advance_front(front,newfront,status,fname);
	}

	start_clock("propagate");

		/* Initialize Newfront */

	start_clock("init_new_front");
	set_size_of_intfc_state(size_of_state(front->interf));
	set_copy_intfc_states(NO);
	set_add_to_correspond_list(YES);
	(*newfront)->interf = pp_copy_interface(front->interf);
	if ((*newfront)->interf == NULL)
	{
	    (void) printf("ERROR in advance_front1d(), "
			  "unable to copy interface\n");
	    return return_advance_front(front,newfront,ERROR_IN_STEP,fname);
	}
	stop_clock("init_new_front");

		/* Propagate the points */

	set_propagation_limits(front,*newfront);
	set_copy_intfc_states(YES);
	intfc_old = front->interf;
	intfc_new = (*newfront)->interf;

	(void) next_point(intfc_old,NULL,NULL,NULL);
	(void) next_point(intfc_new,NULL,NULL,NULL);
	while (next_point(intfc_old,&oldp,&oldhse,&oldhs) && 
	       next_point(intfc_new,&newp,&newhse,&newhs))
	{
	    point_propagate(front,wave,oldp,newp,oldhse,oldhs,dt,V);
	}
	copy_hypersurface_flags(intfc_new);
	debug_front("pt_front","after point propagate",*newfront);

	switch (redistribute(*newfront,YES,NO)) 
	{
	case	GOOD_REDISTRIBUTION:
	    status = GOOD_STEP;
	    break;
	
	case	MODIFY_TIME_STEP_REDISTRIBUTE:
	    (void) printf("WARNING in advance_front1d(), redistribution "
			  "of front failed, reducing time step\n");
	    debug_front("ERROR_front","after error",*newfront);
	    *dt_frac = max(Min_time_step_modification_factor(front),*dt_frac);
	    status = MODIFY_TIME_STEP;
	    break;

	case	UNABLE_TO_UNTANGLE:
	case	BAD_REDISTRIBUTION:
	default:
	    (void) printf("WARNING in advance_front1d(), "
	                  "redistribution of front failed\n");
	    debug_front("ERROR_front","after error",*newfront);
	    *dt_frac = Min_time_step_modification_factor(front);
	    status = ERROR_IN_STEP;
	    break;
	}
	debug_front("redist_front","after redistribute",*newfront);
	if (status != GOOD_STEP)
	    return return_advance_front(front,newfront,status,fname);
	if (!scatter_front(*newfront))
	{
	    (void) printf("WARNING in advance_front1d(), "
	    	          "scatter_front() failed for "
	    	          "normally propagated front\n");
	    return return_advance_front(front,newfront,ERROR_IN_STEP,fname);
	}

	(*newfront)->step = front->step + 1;
	(*newfront)->time = front->time + dt;
	interpolate_intfc_states(intfc_new) = YES;
	set_size_of_intfc_state(size_of_state(intfc_new));

	if (intfc_new->modified)
	    (void) make_point_comp_lists(intfc_new);

	stop_clock("propagate");
	debug_front("new_front","from advance front",*newfront);
	return return_advance_front(front,newfront,status,fname);
}		/*end advance_front1d*/

/*
*			advance_front2d():
*
*	Advances the Front structure front by copying its interface
*	to the interface for newfront and calling
*		
*		(*front->node_propagate)() for each node
*		(*front->curve_propagate)() for each curve
*		point_propagate() for (certain) points on each curve
*		(*front->bond_propagate)() for (certain) bonds on each curve
*		(*front->snd_node_propagate)() for each node
*		(*front->tan_curve_propagate)() for each curve
*	and
*		redistribute()
*
*	in this order to modify the new interface.
*	
*	To define the points and bonds to be propagated in the normal or
*	tangential propagation loop, the flags n_pt_propagated(p)
*	and t_pt_propagated(p) respectively are set to NO.  As each point
*	is processed by the respective propagate step the corresponding
*	flag is set to YES.  To supress a given propagation step,  the
*	the corresponding flag can be turned on prior to propagation.
*
*	Returns 1 if successful, or 0 if copying the interface fails or
*	one of the propagation, or redistribution routines returns 0.
*
*	For documentation of debug_front(), see the documentation for that
*	function in another file (currently fprint.c).
*/


#define is_bad_status(status)   (					\
				((status) == BIFURCATION_NODE) ||	\
				((status) == CROSS_PAST_CURVE_NODE) ||	\
				((status) == MODIFY_TIME_STEP_NODE) ||	\
				((status) == REPEAT_TIME_STEP_NODE) ||	\
				((status) == NO_CROSS_NODE) ||		\
				((status) == ERROR_NODE)		)

LOCAL int advance_front2d(
	double    dt,
	double    *dt_frac,
	Front    *front,
	Front    **newfront,
	POINTER  wave)
{
	CURVE      *oldc,*tempc,*newc;
	CURVE	   **c;
	INTERFACE  *tempintfc;
	NODE 	   *oldn,*tempn,*newn;
	NODE_FLAG  flag;
	RPROBLEM   *rp;
	RPROBLEM   *rp1;
	boolean	   scatter_normally_propagated_front = YES;
	boolean	   scatter_tangentially_propagated_front = YES;
	boolean	   stat;
	boolean       do_redist;
	int        status;
	long       intfc_modified;
	long       redo_advance_front;
	static const char *fname = "advance_front2d()";
	int	   debug_flag = NO;

	debug_print("front","Entered %s(step %d time %g dt %g)\n",fname,
	        	        front->step,front->time,dt);
	debug_front("old_front","into advance front",front);

	*newfront = copy_front(front);
	Interface_redistributed(*newfront) = NO;
	do_redist = (front->num_mts == 0) ? YES : NO;

begin_advance_front2d:
	redo_advance_front = 0;
	tempintfc = NULL;
	rp = NULL;
	set_to_next_node_only(flag);

	set_node_doubly_linked_list(front->interf);

	        /* Initialize Newfront */

	start_clock("init_new_front");
	capture_waves(front);
	print_storage("before init_new_front","ADV_storage");
	        /* TODO: Remove this option!!!!! */
	if (front->init_topology_of_new_interface)
	    status = (*front->init_topology_of_new_interface)(front,*newfront);
	else 
	{
	    set_size_of_intfc_state(size_of_state(front->interf));
	    set_copy_intfc_states(NO);
	    set_add_to_correspond_list(YES);
	    (*newfront)->interf = pp_copy_interface(front->interf);
	    reset_hs_flags_on_intfc((*newfront)->interf);
	    status = ((*newfront)->interf != NULL) ? GOOD_STEP : ERROR_IN_STEP;
	    set_copy_intfc_states(YES);
	}
	if (front->pp_grid)
	    status = syncronize_time_step_status(status,front->pp_grid);
	if (status != GOOD_STEP) 
	{
	    (void) printf("WARNING in advance_front2d(), "
	                  "unable to copy interface\n");
	    status = ERROR_IN_STEP;
	    stop_clock("init_new_front");
	    return return_advance_front(front,newfront,status,fname);
	}
	print_storage("after init_new_front","ADV_storage");
	stop_clock("init_new_front");

	        /* Set Default Propagation Limits */

	set_propagation_limits(front,*newfront);

	        /* Propagate the Curves */

	if (front->intfc_propagate != NULL)
	{
	    start_clock("intfc_propagate");
	    intfc_propagate(front,wave,front->interf,(*newfront)->interf,dt);
	    debug_front("cp_front","after intfc prop",*newfront);
	    stop_clock("curve_propagate");
	}
	else if (front->curve_propagate != NULL) 
	{
	    start_clock("curve_propagate");
	    if (debugging("front"))
	    	(void) printf("Loop over Curves\n");
	    for (c = front->interf->curves; c && *c; ++c)
	    {
	        oldc = *c;
	        if (((newc = correspond_curve(oldc)) != NULL) &&
	    	     (correspond_curve(newc) != NULL))
	        {
	    	    if (debugging("propagate"))
	                (void) printf("\t\tpropagating curve %llu\n",
		                      (long long unsigned int)curve_number(oldc));
		    curve_propagate(front,wave,oldc,newc,dt);
		    /*f_curve_propagate2d */
		}
	    }
	    debug_front("cp_front","after curve prop",*newfront);
	    stop_clock("curve_propagate");
	}

		/* Propagate the Nodes */

	if (debugging("front"))
	{
	    print_correspond_hyper_surf_list(front->interf);
	    print_correspond_hyper_surf_list((*newfront)->interf);
	}
	if (front->node_propagate != NULL) 
	{
	    start_clock("node_propagate");
	    set_corresponds_for_node_prop(front->interf,(*newfront)->interf);
	    oldn = first_node(front->interf);
	    while (oldn != NULL) 
	    {
	        newn = correspond_node(oldn);
	        if (debugging("crx_status"))
	            print_linked_node_list((*newfront)->interf);
	        status = (newn != NULL) ?
	            (*front->node_propagate)(front,wave,oldn,newn,&rp,
	        			     dt,dt_frac,flag,NULL) : GOOD_NODE;

	        if (debugging("crx_status"))
	        if (is_bad_status(status) &&
	          (point_in_buffer(Coords(oldn->posn),front->rect_grid) == YES))
	        {
	            print_node_status("WARNING in advance_front2d(), "
	                              "node_propagation returns ",status,"\n");
	            (void) printf("Problem occurs in buffer zone - ignoring\n");
	            if (set_node_states_and_continue(oldn,newn,front))
	                status = GOOD_NODE;
	        }

	        switch (status) 
	        {
	        case GOOD_NODE:
	            oldn = adv_node_loop_after_good_prop(oldn,newn,&rp);
	            break;
	        case PSEUDOCROSS_NODE_NODE:
	            debug_print("PSEUDOCROSS","PSEUDOCROSS case\n");
	            oldn = reorder_node_loop(oldn,newn);
	            break;
	        case CROSS_NODE_NODE:
	        case BIFURCATION_NODE:
	            debug_print("CROSS","CROSS case\n");
	            oldn = next_node(oldn);
	            break;
	        case CROSS_PAST_CURVE_NODE:
	            print_node_status("WARNING in advance_front2d(), "
	                              "node_propagate failed with status ",
				      status,"\n");
	            print_node(oldn);
	            if (debugging("CROSS_PAST"))
	            {
	                (void) printf("Cross past curve case\n"
	                              "dt_frac = %g\n",*dt_frac);
	                (void) printf("Reducing time step\n");
	            }
		    status = node_modify_time_step(oldn,front,dt_frac,
					           MODIFY_TIME_STEP);
	            free_rp_list(&rp);
	            goto sync_prop_stat1;
	        case MODIFY_TIME_STEP_NODE:
	            (void) printf("WARNING in advance_front2d(), "
	                          "node_propagate returns "
	                          "MODIFY_TIME_STEP_NODE\n");
	            free_rp_list(&rp);
		    status = node_modify_time_step(oldn,front,NULL,
						   MODIFY_TIME_STEP);
	            goto sync_prop_stat1;
	        case REPEAT_TIME_STEP_NODE:
	            (void) printf("WARNING in advance_front2d(), "
	                          "node_propagate returns "
	                          "REPEAT_TIME_STEP_NODE\n");
	            free_rp_list(&rp);
		    status = node_modify_time_step(oldn,front,NULL,
					           REPEAT_TIME_STEP);
	            goto sync_prop_stat1;
	        case NO_CROSS_NODE:
	            print_node_status("WARNING in advance_front2d(), "
	                              "node_propagate failed with status ",
	                              status,"\n");
	            print_node(oldn);
	            if (debugging("NO_CROSS"))
	            {
	                (void) printf("No cross case\n");
	                (void) printf("dt_frac = %g\n",*dt_frac);
	                (void) printf("Reducing time step\n");
	            }
	            free_rp_list(&rp);
		    status = node_modify_time_step(oldn,front,dt_frac,
					           MODIFY_TIME_STEP);
	            goto sync_prop_stat1;
	        case ERROR_NODE:
	        default:
	            print_node_status("WARNING in advance_front2d(), "
	                              "node_propagate failed with status ",
	                              status,"\n");
	            print_node(oldn);
	            if (debugging("ERROR_NODE"))
	            {
	                (void) printf("Old interface:\n");
	                print_interface(front->interf);
	                print_correspond_hyper_surf_list(front->interf);
	                (void) printf("New interface:\n");
	                print_interface((*newfront)->interf);
	                print_correspond_hyper_surf_list((*newfront)->interf);
	            }
		    status = node_modify_time_step(oldn,front,dt_frac,
					           ERROR_IN_STEP);
	            free_rp_list(&rp);
	            goto sync_prop_stat1;
	        }
	    } /* end of while (oldn != NULL) */
	    set_correspond_hyper_surf_bdrys_to_NULL(front->interf);
	    set_correspond_hyper_surf_bdrys_to_NULL((*newfront)->interf);
	    if (rp && (front->twodrproblem != NULL)) 
	    {
	        for (rp1 = rp; rp1; rp1 = rp1->prev) 
	        {
	            debug_front("2drp_front",
	                "new between node loop and rp loop",*newfront);
	                    
	            status = (*front->twodrproblem)(front,*newfront,wave,&rp1);

	            /* At this point, rp is nothing more than a valid element
	             * of the list which provides a starting point
	             * for deleting the list.  If we delete an element of
	             * the list in front->twodrproblem (presumably due to
	             * merging two RPROBLEM's), then rp may point to freed
	             * storage and will need to be updated.  rp1 should still
	             * be a valid element of the list.
	             */
	            rp = rp1;

		    if (status != GOOD_STEP)
		    {
	                print_time_step_status("WARNING in advance_front2d(), "
					       "rp failed with status = ",
					       status,"\n");
	                switch (status) 
	                {
	                case GOOD_STEP:
			    break;

		        case REPEAT_TIME_STEP:
	                    break;

	                case MODIFY_TIME_STEP:
		            status = rp_modify_time_step(rp1,front,status);
			    if (status == MODIFY_TIME_STEP)
			    {
	                        *dt_frac = rp1->dt_frac;
	                        if (debugging("2drp"))
	                        {
	                            print_rproblem(rp1);
	                            (void) printf("dt_frac %g\n",*dt_frac);
	                            (void) printf("Reducing time step\n");
	                        }
			        *dt_frac = limit_dt_frac(*dt_frac,front);
			    }
	                    break;

	                case ERROR_IN_STEP:
	                default:
	                    print_rproblem(rp1);
	                    /* Try reducing the time step */
		            status = rp_modify_time_step(rp1,front,status);
	                    if (status == MODIFY_TIME_STEP)
	                        *dt_frac *=
				    TIME_STEP_REDUCTION_FACTOR(front->interf);
	                    break;
	                }
		    }
	            if (status != GOOD_STEP)
			break;
	        }
	        free_rp_list(&rp);
	        debug_front("2drp_front","after 2drp loop",*newfront);
	    }
	    else if (rp) 
	    {
	        for (rp1 = rp; rp1; rp1 = rp1->prev) 
	            print_rproblem(rp1);
	        free_rp_list(&rp);
	        (void) printf("WARNING in advance_front2d(), "
	                      "CROSS code needed\n");
	        status = ERROR_IN_STEP;
	    }

sync_prop_stat1:
	    stop_clock("node_propagate");
	    if (front->pp_grid)
	    	status = syncronize_time_step_status(status,front->pp_grid);
	    if (status != GOOD_STEP)
	        return return_advance_front(front,newfront,status,fname);
	}
	if (*front->max_scaled_propagation > 0.5)
	{
	    (void) printf("WARNING in advance_front2d(), "
	                  "front->max_scaled_propagation = %f\n",
			  *(front->max_scaled_propagation));
	    *dt_frac = 0.4/(*front->max_scaled_propagation);
	    status = MODIFY_TIME_STEP;
	    goto sync_prop_stat2;
	}

	stat = consistent_propagated_loop_orientations(dt,dt_frac,front,wave);
	if (stat == NO)
	{
	    (void) printf("WARNING in advance_front2d(), "
		          "Inconsistent orientation of propagated loop "
	                  "detected after point and node propagations");
	    if (pp_numnodes() > 1)
		(void) printf(" on processor %d\n",pp_mynode());
	    else
		(void) printf("\n");
	}
	if (pp_min_status(stat) == NO)
	{
	    if (stat == YES)
	    {
	        (void) printf("WARNING in advance_front2d(), "
		              "Inconsistent orientation of propagated loop "
	                      "detected on a remote processor "
			      "after point and node propagations ");
	    }
	    status = MODIFY_TIME_STEP;
	    goto sync_prop_stat2;
	}

	/* Make Temp Interface for Tangential Propagation */

	set_node_doubly_linked_list((*newfront)->interf);
	if (front->snd_node_propagate) 
	{
	    start_clock("snd_copy_interface");
	    print_storage("before snd_copy_interface","ADV_storage");
	    tempintfc = (*newfront)->interf;
	    set_size_of_intfc_state(size_of_state(tempintfc));
	    set_add_to_correspond_list(YES);
	    if (((*newfront)->interf = pp_copy_interface(tempintfc)) == NULL)
	    {
	    	(void) printf("WARNING in advance_front2d(), "
		              "unable to copy interface\n");
		status = ERROR_IN_STEP;
		goto sync_prop_stat2;
	    }
	    copy_hypersurface_flags((*newfront)->interf);
	    print_storage("after snd_copy_interface","ADV_storage");
	    stop_clock("snd_copy_interface");
	}
	interpolate_intfc_states((*newfront)->interf) = YES;

	/* Second Propagation for the States Around the Nodes */

	if (front->snd_node_propagate) 
	{
	    start_clock("snd_node_propagate");
	    if (debugging("front"))
	    	(void) printf("Second Loop over Nodes\n");

	    tempn = first_node(tempintfc);
	    newn = first_node((*newfront)->interf);
	    while (newn != NULL)
	    {
	    	(*front->snd_node_propagate)(front,*newfront,wave,
	    				     tempintfc,tempn,newn,dt);
	    	tempn = next_node(tempn);
	    	newn = next_node(newn);
	    }

	    debug_front("snd_front","after snd_node prop",*newfront);
	    stop_clock("snd_node_propagate");
	}

	if (tempintfc)
	    (void) delete_interface(tempintfc);
	print_storage("after delete tempintfc","ADV_storage");

		/* Redistribute the New Front */

	switch (redistribute(*newfront,do_redist,NO)) 
	{
	case	GOOD_REDISTRIBUTION:
	    status = GOOD_STEP;
	    break;
	
	case	UNABLE_TO_UNTANGLE:
	    (void) printf("WARNING in advance_front2d(), "
	                  "redistribution of front failed\n"
	                  "Restarting advance_front2d()\n");
	    *dt_frac = Min_time_step_modification_factor(front);
	    status = MODIFY_TIME_STEP;
	    break;

	case	MODIFY_TIME_STEP_REDISTRIBUTE:
	    (void) printf("WARNING in advance_front2d(), "
	                  "redistribute returns\n"
	                  "\t\tMODIFY_TIME_STEP_REDISTRIBUTE, dt_frac = %g\n",
			  *dt_frac);
	    *dt_frac = Min_time_step_modification_factor(front);
	    status = MODIFY_TIME_STEP;
	    break;
		
	case	BAD_REDISTRIBUTION:
	default:
	    (void) printf("WARNING in advance_front2d(), "
	                  "redistribution of front failed\n");
	    debug_front("ERROR_front","after error",*newfront);
	    *dt_frac = Min_time_step_modification_factor(front);
	    status = MODIFY_TIME_STEP;
	    break;
	}
	if (front->pp_grid)
	    status = syncronize_time_step_status(status,front->pp_grid);
	if (status != GOOD_STEP)
	    return return_advance_front(front,newfront,status,fname);

	Redistribution_count(front) = Redistribution_count(*newfront);
	(*newfront)->step = front->step + 1;
	(*newfront)->time = front->time + dt;
	debug_front("redist_front","after redistribution",*newfront);

	/* Communicate topologically propagated front */
	if (scatter_normally_propagated_front == YES)
	{
	    start_clock("scatter_front");
	    if (!scatter_front(*newfront))
	    {
	    	(void) printf("WARNING in advance_front2d(), "
	    	              "scatter_front() failed for "
	    	              "normally propagated front\n");
	    	scatter_normally_propagated_front = NO;
	    	scatter_tangentially_propagated_front = NO;
	    	(void) delete_interface((*newfront)->interf);
	    	(*newfront)->interf = NULL;
	    	goto begin_advance_front2d;
	    }
	    stop_clock("scatter_front");
	}

	debug_front("node_front","after node loop",*newfront);
	if (debugging("front"))
	{
	    print_correspond_hyper_surf_list(front->interf);
	    print_correspond_hyper_surf_list((*newfront)->interf);
	}

	if (front->mass_consv_diagn_driver)
	    (*front->mass_consv_diagn_driver)(front,wave,dt);

	if (debugging("bond_lengths"))
	    check_bond_lengths((*newfront)->interf);

	/* Check for the geometric orientation of loops */

	/* ONLY check loops that will not be deleted !!!! */
	delete_small_loops(*newfront);

		/* Delete non-boundary curves that lie  */
		/* fully on or exterior to the boundary */

	delete_exterior_curves(*newfront,front->interf);
	intfc_delete_fold_back_bonds(*newfront);
	debug_front("dec_front","after delete_exterior_curves:",*newfront);

	interpolate_intfc_states((*newfront)->interf) = YES;

	/* Make Temp Interface for Tangential Propagation */

	if (front->tan_curve_propagate) 
	{
	    start_clock("snd_copy_interface");
	    print_storage("before snd_copy_interface","ADV_storage");
	    tempintfc = (*newfront)->interf;
	    set_size_of_intfc_state(size_of_state(tempintfc));
	    set_add_to_correspond_list(YES);
	    if (((*newfront)->interf = pp_copy_interface(tempintfc)) == NULL)
	    {
	    	(void) printf("WARNING in advance_front2d(), "
		              "unable to copy interface\n");
		status = ERROR_IN_STEP;
		goto sync_prop_stat2;
	    }
	    copy_hypersurface_flags((*newfront)->interf);
	    interpolate_intfc_states((*newfront)->interf) = YES;
	    print_storage("after snd_copy_interface","ADV_storage");
	    stop_clock("snd_copy_interface");
	}

	/* Tangential Sweep for States on the Curves */

	if (front->tan_curve_propagate) 
	{
	    start_clock("tan_curve_propagate");
	    if (debugging("front"))
	    	(void) printf("Second Loop over Curves\n");
	    for (c = tempintfc->curves; c && *c; ++c)
	    {
	    	tempc = *c;
	    	newc = correspond_curve(tempc);
	    	(*front->tan_curve_propagate)(front,*newfront,
	    				      tempintfc,tempc,newc,dt);
	    }
	    debug_front("tcp_front","after tan_curve_propagate:",*newfront);
	    stop_clock("tan_curve_propagate");
	}
	if (tempintfc)
	    (void) delete_interface(tempintfc);
	print_storage("after delete tempintfc","ADV_storage");


		/* Provide robustness for untangle algorithms */

		/*   delete remnants of scalar physical   */
		/*  curves sticking to NEUMANN boundaries */
		/* Add to delete_exterior_curves()? */

	if (pp_min_status(delete_phys_remn_on_bdry(*newfront)) == NO)
	{
	    (void) printf("WARNING in advance_front2d(), "
	                  "delete_phys_remn_on_bdry() detected error\n");
	    debug_front("ERROR_front","after error",*newfront);
	    *dt_frac = Min_time_step_modification_factor(front);
	    status = MODIFY_TIME_STEP;
	    goto sync_prop_stat2;
	}
	debug_front("dspr_front",
		    "after 1st delete_phys_remn_on_bdry():",*newfront);

sync_prop_stat2:
	if (front->pp_grid)
	    status = syncronize_time_step_status(status,front->pp_grid);
	if (status != GOOD_STEP)
	    return return_advance_front(front,newfront,status,fname);

	/* Communicate tangentially propagated front */
	if (scatter_tangentially_propagated_front == YES)
	{
	    start_clock("scatter_front");
	    if (!scatter_front(*newfront))
	    {
	    	(void) printf("WARNING in advance_front2d(), "
	    	              "scatter_front() failed for "
	    	              "tangentially propagated front\n");
	    	scatter_normally_propagated_front = NO;
	    	scatter_tangentially_propagated_front = NO;
	    	(void) delete_interface((*newfront)->interf);
	    	(*newfront)->interf = NULL;
	    	goto begin_advance_front2d;
	    }
	    stop_clock("scatter_front");
	}

	if (status != GOOD_STEP)
	    return return_advance_front(front,newfront,status,fname);


		/* 	Post-process newfront->interf	   */
		/* Provide robustness after redistribution */
		/*   for node propagate on next time step  */

		/* Delete non-boundary curves that lie  */
		/* fully on or exterior to the boundary */

	delete_exterior_curves(*newfront,front->interf);
	debug_front("dec_front","after delete_exterior_curves:",*newfront);

		/*  delete remnants of scalar physical    */
		/* curves sticking to NEUMANN boundaries  */
		/* Add to delete_exterior_curves()? */

	if (pp_min_status(delete_phys_remn_on_bdry(*newfront)) == NO)
	{
	    (void) printf("WARNING in advance_front2d(), "
	                  "delete_phys_remn_on_bdry() detected error\n");
	    debug_front("ERROR_front","after error",*newfront);
	    *dt_frac = Min_time_step_modification_factor(front);
	    status = MODIFY_TIME_STEP;
	    return return_advance_front(front,newfront,status,fname);
	}
	debug_front("dspr_front",
		    "after 2nd delete_phys_remn_on_bdry():",*newfront);


	/* These guys keep sneaking through !! */
	/* This should be the most effective place for this call */
	/* Brent - I believe it is better to have the function at
	*  the end of advance_front2d() applied to the newfront
	*  instead of at the beginning applied to front.
	*  In general our policy should be never to modify the
	*  old interface data.
	*/
	delete_small_loops(*newfront);
	debug_front("dsloop_front","after delete_small_loops():",*newfront);

	test_for_mono_comp_curves((*newfront)->interf);

	/* Check if post processing has changed topology */

	intfc_modified = (*newfront)->interf->modified;
	pp_global_lmax(&intfc_modified,1L);
	if (intfc_modified)
	{
	    if (!scatter_front(*newfront))
	    {
	    	(void) printf("WARNING in advance_front2d(), "
	    	              "final scatter_front() failed\n");
	    	*dt_frac = Max_time_step_modification_factor(front);
	        return return_advance_front(front,newfront,
		                            MODIFY_TIME_STEP,fname);
	    }
	    stat = make_bond_comp_lists((*newfront)->interf);
	    if (pp_min_status(stat) == FUNCTION_FAILED)
	    {
	    	screen("ERROR in advance_front2d(), "
	    	       "make_bond_comp_lists() failed\n");
	    	clean_up(ERROR);
	    }
	}

	return return_advance_front(front,newfront,GOOD_STEP,fname);
}		/*end advance_front2d*/

LOCAL	int node_modify_time_step(
	NODE  *oldn,
	Front *front,
	double *dt_frac,
	int   status)
{
	if (status == GOOD_STEP)
	    return status;
	if (last_time_step_modification(front) == YES)
	{
	    return delete_phys_curves_at_old_node(oldn,status);
	}
	if (dt_frac != NULL)
	{
	    if (status == ERROR_IN_STEP)
	    {
	        *dt_frac = Min_time_step_modification_factor(front);
		status = MODIFY_TIME_STEP;
	    }
	    else
	        *dt_frac *= TIME_STEP_REDUCTION_FACTOR(front->interf);
	}
	return status;
}		/*end node_modify_time_step*/

EXPORT	int rp_modify_time_step(
	RPROBLEM *rp,
	Front    *front,
	int      status)
{
	RP_NODE   *rpn;

	if (status == ERROR_IN_STEP)
	    status = MODIFY_TIME_STEP;
	if ((status==GOOD_STEP) || (last_time_step_modification(front)!=YES))
	    return status;
	for (rpn = rp->first_rp_node; rpn != NULL; rpn = rpn->next)
	    status = delete_phys_curves_at_old_node(rpn->old_node,status);
	return status;
}		/*end rp_modify_time_step*/

EXPORT	double	limit_dt_frac(
	double	dt_frac,
	Front	*front)
{
	if (dt_frac < Min_time_step_modification_factor(front))
	    dt_frac = Min_time_step_modification_factor(front);
	if (dt_frac > Max_time_step_modification_factor(front))
	    dt_frac = Max_time_step_modification_factor(front);
	return dt_frac;
}		/*end limit_dt_frac*/

LOCAL	int delete_phys_curves_at_old_node(
	NODE  *oldn,
	int   status)
{
	boolean untrack = NO;
	CURVE   **c;
	int     wtype = FIRST_VECTOR_PHYSICS_WAVE_TYPE;

	for (c = oldn->in_curves; c && *c; ++c)
	{
	    if (wave_type(*c) < wtype)
		continue;
	    untracked_hyper_surf(*c) = YES;
	    untrack = YES;
	    status = REPEAT_TIME_STEP;
	}
	for (c = oldn->out_curves; c && *c; ++c)
	{
	    if (wave_type(*c) < wtype)
		continue;
	    untracked_hyper_surf(*c) = YES;
	    untrack = YES;
	    status = REPEAT_TIME_STEP;
	}
	if (untrack == YES)
	{
	    (void) printf("WARNING in delete_phys_curves_at_old_node(), "
			  "turning off tracking of vector curves at node\n");
	    print_node(oldn);
	}
	return status;
}		/*end delete_phys_curves_at_old_node*/

EXPORT void reset_hs_flags_on_intfc(
	INTERFACE	*intfc)
{
	HYPER_SURF	**hs;

	if (intfc == NULL)
	    return;

	for (hs = hyper_surf_list(intfc); hs && *hs; ++hs)
	{
	    /*Turn off flags force || preventing redistribution*/
	    /*never_redistribute(*hs) = NO;*/
	    do_not_redistribute(*hs) = NO;
	    redistribute_hyper_surface(*hs) = NO;
	    redistributed(*hs) = NO;
	    perform_redistribution_function(*hs) =
	        default_perform_redistribution_function(intfc);
	}
}		/*end reset_hs_flags_on_intfc*/

EXPORT	NODE *adv_node_loop_after_good_prop(
	NODE		*oldn,
	NODE		*newn,
	RPROBLEM	**rp)
{
	NODE		*n;
	NODE		*next_oldn;
	RPROBLEM	*rp1;
	RP_NODE		*rpn;

	if ((rp != NULL) && (*rp != NULL) && 
	    (newn != NULL) && (node_type(newn) != SUBDOMAIN_NODE) &&
	    (!is_virtual_fixed_node(newn)))
	{
	    for (rp1 = *rp; rp1; rp1 = rp1->prev) 
	    {
	    	for (rpn = rp1->first_rp_node; rpn; rpn = rpn->next)
	    	{
	    	    if ((rpn->node == newn)) 
	    	    	break;
	    	}
	    	if (rpn == NULL)
		    continue;

	    		/* Delete rp1 from rp list */

	    	if (rp1->prev)
		    rp1->prev->next = rp1->next;
	    	if (rp1->next)
		    rp1->next->prev = rp1->prev;
	    	else
		    *rp = rp1->prev;
				
	    	    /* Push other nodes in rp1 to end of node list */

	    	for (rpn = rp1->first_rp_node; rpn; rpn = rpn->next)
	    	{
	    	    if ((rpn->node == newn)) 
	    	    	continue;

	    	    if (propagation_status(rpn->node) == PROPAGATED_NODE)
	    		continue;

	    	    /* Check if node already follows newn */

	       	    for (n = next_node(newn); n != NULL; n = next_node(n))
	    	    {
	    	    	if (n == rpn->node)
			    break;
	    	    }
	    	    if (n != NULL)
			continue;

	    	    n = rpn->node;
	    	    (void) reorder_node_loop(rpn->old_node,rpn->node);
	    	}
	    	free_rp(rp1);
	    }
	}

	/*
	*  Reset propagation status to UNPROPAGATED_NODE
	*  after successful node propagation.
	*/
	next_oldn = next_node(oldn);
	for (oldn = next_oldn; oldn != NULL; oldn = next_node(oldn))
	{
	    if ((newn = correspond_node(oldn)) == NULL)
		continue;
	    if (propagation_status(newn) != PROPAGATED_NODE)
	    	propagation_status(newn) = UNPROPAGATED_NODE;
	}
	return next_oldn;
}		/*end adv_node_loop_after_good_prop*/



/*
*			reorder_node_loop():
*
*	Rearranges the node loop so that *n is placed at the end of the list,
*	and then reset to what was the value of (*n)->next when the function
*	was entered.  Returns the next node in the list after n.
*/

EXPORT NODE *reorder_node_loop(
	NODE *oldn,
	NODE *newn)
{
	INTERFACE *intfc;
	NODE	  *next, *prev;
	NODE      *next_old_node;
	
	if (oldn != NULL)
	{
	    intfc = oldn->interface;
	    next = next_node(oldn);
	    if (next == NULL) 
	    {
	        screen("ERROR in reorder_node_loop(), "
		       "next_node(oldn) == NULL\n");
	        print_node(oldn);
		(void) printf("Old interface information\n");
	        print_linked_node_list(intfc);
	        print_interface(intfc);
		if (newn != NULL)
		{
	            intfc = newn->interface;
		    (void) printf("New interface information\n");
	            print_linked_node_list(intfc);
	            print_interface(intfc);
		}
	        clean_up(ERROR);
	    }
	    prev = prev_node(oldn);
	    prev_node(next) = prev;
	    if (prev != NULL)
	        next_node(prev) = next;
	    else
	        first_node(intfc) = next;
	    next_node(last_node(intfc)) = oldn;
	    next_node(oldn) = NULL;
	    prev_node(oldn) = last_node(intfc);
	    last_node(intfc) = oldn;
	    next_old_node = next;
	}
	else
	    next_old_node = NULL;

	if (newn != NULL)
	{
	    intfc = newn->interface;
	    next = next_node(newn);
	    if (next == NULL) 
	    {
	        screen("ERROR in reorder_node_loop(), "
		       "next_node(newn) == NULL\n");
	        print_node(newn);
		if (oldn != NULL)
		{
		    (void) printf("Old interface information\n");
	            print_linked_node_list(intfc);
	            print_interface(intfc);
		}
	        intfc = newn->interface;
		(void) printf("New interface information\n");
	        print_linked_node_list(intfc);
	        print_interface(intfc);
	        clean_up(ERROR);
	    }
	    prev = prev_node(newn);
	    prev_node(next) = prev;
	    if (prev != NULL)
	        next_node(prev) = next;
	    else
	        first_node(intfc) = next;
	    next_node(last_node(intfc)) = newn;
	    next_node(newn) = NULL;
	    prev_node(newn) = last_node(intfc);
	    last_node(intfc) = newn;
	}
	return next_old_node;
}		/*end reorder_node_loop*/


EXPORT	void set_corresponds_for_node_prop(
	INTERFACE	*ointfc,
	INTERFACE	*nintfc)
{
	NODE		**on, **nn;

	on = ointfc->nodes;	nn = nintfc->nodes;
	for (;on && *on; ++on, ++nn)
	{
	    correspond_hyper_surf_bdry(*on) = Hyper_surf_bdry(*nn);
	    correspond_hyper_surf_bdry(*nn) = Hyper_surf_bdry(*on);
	}
}		/*end set_corresponds_for_node_prop*/

EXPORT	void print_linked_node_list(
	INTERFACE	*intfc)
{
	NODE		**n, *m;
	int		i, node_count;

	n = intfc->nodes;
	if (! n)
	{
	    (void) printf("NULL node list on intfc\n");
	    return;
	}

	(void) printf("\tnode list - intfc %llu\n",(long long unsigned int)interface_number(intfc));
	for (node_count = 0;  n && *n;  ++n, ++node_count)
		;

	m = first_node(intfc);
	for (i = 0;  i <= node_count + 2;  ++i) 
	{
	    if (m != NULL)
	    {
	    	(void) printf("prev %llu  m %llu  next %llu  ",
	    		      (long long unsigned int)node_number(prev_node(m)),
			      (long long unsigned int)node_number(m),
	    		      (long long unsigned int)node_number(next_node(m)));
	    	print_propagation_status(m);
	    }
	    else
	    	break;
	    m = next_node(m);
	}
}		/*end print_linked_node_list*/


EXPORT void set_node_doubly_linked_list(
	INTERFACE	*intfc)
{
	NODE		**n;

	n = intfc->nodes;
	first_node(intfc) = *n;
	prev_node(*n) = NULL;
	for (; *(n+1); ++n) 
	{
	    next_node(*n) = *(n+1);
	    prev_node(*(n+1)) = (*n);
	}
	last_node(intfc) = *n;
	next_node(*n) = NULL;
}		/*end set_node_doubly_linked_list*/



/*ARGSUSED*/
EXPORT	boolean consistent_propagated_loop_orientations(
	double		dt,
	double		*dt_frac,
	Front		*fr,
	POINTER		wave)
{
	static const double AREA_FAC = 0.01;	/*TOLERANCE*/
	static const int   num_attempts = 5;	/*TOLERANCE*/
	CURVE	**c,*oldc,*newc;
	NODE	*oldn, *newn;
	double	dt_f, dT;
	double	A_old, A_new;
	double	min_area = AREA_FAC*fr->rect_grid->h[0]*fr->rect_grid->h[1];
	boolean	status = YES;
	int	i;

	debug_print("orient_consistency",
	      "Entered consistent_propagated_loop_orientations()\n");

	for (c = fr->interf->curves; c && *c; ++c)
	{
	    oldc = *c;
	    if (!is_closed_curve(oldc))
		continue;
	    newc = correspond_curve(oldc);

	    /*TODO: this needs to be solved correctly.*/
	    if ((newc == NULL) || (newc->interface == NULL) ||
	        (!is_closed_curve(newc)))
	    	continue;
	    if (newc->num_points > 20)
	    	continue;	/* sufficiently long curve cannot fully flip */

	    A_old = area_of_closed_curve(oldc);
	    A_new = area_of_closed_curve(newc);

	    if (debugging("orient_consistency"))
	    	(void) printf("A_old = %g, A_new = %g\n",A_old,A_new);

	    if (A_old*A_new > 0.0)
		continue;

	    status = NO;
	    if (debugging("orient_consistency"))
	    {
	    	(void) printf("Loop %llu inverted after propagation\n",
	    		      (long long unsigned int)curve_number(oldc));
	    	(void) printf("Old loop\n");
	    	print_curve(oldc);
	    	(void) printf("Propagated loop\n");
	    	print_curve(newc);
	    }

	    /* Determine dt_frac */

	    if (A_old < min_area)
	    {
	    	*dt_frac *= 0.9;
	    	continue;
	    }

	    dt_f = 1.0;
	    oldn = oldc->start;	newn = newc->start;
	    for (i = 0; i < num_attempts; ++i)
	    {
	    	dt_f *= 0.5;
	    	dT = dt_f * dt;
	    	if (fr->curve_propagate != NULL) 
	    	    curve_propagate(fr,wave,oldc,newc,dT);
		(void) closed_node_propagate(fr,wave,oldn,newn,dT);
		A_new = area_of_closed_curve(newc);

		if (debugging("orient_consistency"))
		{
	            (void) printf("In loop, dt_f = %g, dT = %g, ",dt_f,dT);
		    (void) printf("A_old = %g, A_new = %g\n",A_old,A_new);
		}
		if (A_old*A_new >= 0.0)
		    break;
	    }
	    if (i == num_attempts)
	    {
	        (void) printf("WARNING in "
			      "consistent_propagated_loop_orientations(), "
		              "cannot produce consistent loop "
		              "orientations by reduction of time step\n");
		(void) printf("modifying interface\n");
		(void) delete_curve(newc);
		status = YES;
		continue;
	    }
	    *dt_frac = min(*dt_frac,dt_f);

		/* Restore newc to original state ? 			*/
		/* This is unnecessary if reduce time step is called,   */
		/* as the newc->interface will be deleted, but this	*/
		/* may be desirable for consistency reasons if future	*/
		/* resolutions solve this problem without having to	*/
		/* reduce the time step.  For now this section is	*/
		/* commented out.					*/

	}
	debug_print("orient_consistency",
		"Left consistent_propagated_loop_orientations()\n");
	return status;
}		/*end consistent_propagated_loop_orientation*/



LOCAL	void check_bond_lengths(
	INTERFACE	*intfc)
{
	CURVE		*c;
	BOND		*b;
	double		oldbl, newbl;

	(void) next_bond(intfc,NULL,NULL);
	while (next_bond(intfc,&b,&c))
	{
	    oldbl = bond_length(b);
	    newbl = separation(b->end,b->start,intfc->dim);

	    if (oldbl != newbl)
	    {
	    	(void) printf("WARNING in check_bond_lengths(), "
	    	              "bond lengths not correct at end of "
	    	              "advance\nold %g new %g dlen %g\n",
	    		      oldbl,newbl,oldbl-newbl);
	    	(void) printf("bond - ");
		print_bond(b);
	    	(void) printf("\n");
	    }
	    bond_length(b) = newbl;
	}
}		/*end check_bond_lengths*/

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
*				fadv3d.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*/

#include <front/fdecs.h>		/* includes int.h, table.h */

	/* LOCAL Function Declarations */
LOCAL	boolean	BothSidesActive(HYPER_SURF*,Front*);
LOCAL	int	propagate_3d_front(POINTER,Front*,Front*,double,double*,boolean);
LOCAL   int	propagate_node_points(Front*,Front*,POINTER,double,double*);
LOCAL	void	propagate_curve_points(Front*,Front*,POINTER,double);
LOCAL   void    propagate_surface_points(Front*,Front*,POINTER,double,double*);
LOCAL 	int 	propagate_points_tangentially(Front*,Front*,int,double,double*,
						  int);
LOCAL	int	simple_advance_front3d(double,double*,Front*,Front**,POINTER);
LOCAL	int	struct_advance_front3d(double,double*,Front*,Front**,POINTER);
LOCAL 	int 	hybrid_advance_front3d(double,double*,Front*,Front**,POINTER);
LOCAL 	int 	mixed_advance_front3d(double,double*,Front*,Front**,POINTER);
LOCAL	int	spherical_advance_front3d(double,double*,Front*,
					 	Front**,POINTER);
LOCAL	int	preserve_advance_front3d(double,double*,Front*,
					 	Front**,POINTER);
LOCAL	int	reconstruct_advance_front3d(double,double*,Front*,
						Front**,POINTER);

LOCAL	void	debug_propagate_3d_front(Front*);
LOCAL	void	unfold_interface_section(POINT*,POINT*,TRI*,TRI*,
	                                 SURFACE*,SURFACE*);
LOCAL	void	detach_and_propagate_curves(Front*,Front*,POINTER,double);
LOCAL 	int 	resolve_front_collision(Front*);

EXPORT	void	EnforceFlowSpecifedStates3d(
	Front	*fr)
{
	HYPER_SURF		*hs;
	INTERFACE		*intfc;
	SURFACE			**s;
	POINT			*p;
	TRI			*t;
	Locstate		sl, sr;
	int			i;

	if ((fr==NULL) || (Fsr_list(fr)==NULL) || (fr->rect_grid->dim!=3) )
	    return;

	intfc = fr->interf;
	for (s = intfc->surfaces; s && *s; ++s)
	{
	    hs = Hyper_surf(*s);
	    if (is_subdomain_boundary(hs) || is_passive_boundary(hs))
	        continue;
	    if (BothSidesActive(hs,fr) == YES)
	    	continue;
	    for (t = first_tri(*s); !at_end_of_tri_list(t,*s); t = t->next)
	    {
		for (i = 0; i < 3; ++i)
		{
		    p = Point_of_tri(t)[i];
	    	    slsr(p,Hyper_surf_element(t),hs,&sl,&sr);
	    	    (void) RegionIsFlowSpecified(sl,NULL,Coords(p),
	                                         negative_component(hs),
						 NO_COMP,fr);
	    	    (void) RegionIsFlowSpecified(sr,NULL,Coords(p),
	                                         positive_component(hs),
						 NO_COMP,fr);
		}
	    }
	}
}		/*end EnforceFlowSpecifedStates3d*/

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


EXPORT 	int 	advance_front3d_tracking_control(
	double	dt,
	double	*dt_frac,
	Front	*front,
	Front	**newfront,
	POINTER	wave)
{
	int step_status;

	DEBUG_ENTER(advance_front3d_tracking_control)

  	switch(Tracking_algorithm(front))
	{
	case NO_DYNAMIC_TRACKING:
	    step_status = GOOD_STEP;
	    break;

	case SIMPLE_TRACKING:
	    step_status = simple_advance_front3d(dt,dt_frac,front,
                                                         newfront,wave);
	    break;

	case STRUCTURE_TRACKING:
	    step_status = struct_advance_front3d(dt,dt_frac,front,
							 newfront,wave);
	    break;

	case SPHERICAL_TRACKING:
	    step_status = spherical_advance_front3d(dt,dt_frac,front,
							 newfront,wave);
	    break;

	case GRID_FREE_TRACKING:
	    step_status = preserve_advance_front3d(dt,dt_frac,front,
							 newfront,wave);
	    break;

	case GRID_BASED_TRACKING:
	    step_status = reconstruct_advance_front3d(dt,dt_frac,front,
						            newfront,wave);
	    break;
	
	case THREE_COMP_GRID_BASED_TRACKING:
	    step_status = reconstruct_advance_front3d(dt,dt_frac,front,
						            newfront,wave);
	    break;

	case MIXED_TRACKING:
	case LOCALLY_GRID_BASED_TRACKING:
	    step_status = mixed_advance_front3d(dt,dt_frac,front,
						newfront,wave);
	    break;

	case HYBRID_TRACKING:
	    step_status = hybrid_advance_front3d(dt,dt_frac,front,
						 newfront,wave);
	    break;

	default:
	    screen("ERROR in advance_front3d_tracking_control()\n"
		   "Unknown Tracking algorithm = %d\n",
		   Tracking_algorithm(front));
	    clean_up(ERROR);
	}
	if (step_status == GOOD_STEP)
            front->redis_flag = (*newfront)->redis_flag;
	else
	    front->redis_flag = NO;
	DEBUG_LEAVE(advance_front3d_tracking_control)
	return step_status;
}		/*end advance_front3d_tracking_control*/

/*
*		reconstruct_advance_front3d():
*
*	   The following program is the routine for three-dimensional
*	front advance.  It is called from each of the hyp_drivers in 
*	hyp/hdriver.c, since advancing the front is the first step 
*	of the hyperbolic solution.
*	   The advance front involves three steps; copy the front and 
*	interface,  perform a propagate 3d front, and redistribute the
*	front. 
*
*	   Note - The function appears here to be LOCAL, but access is 
*	provided by a pointer set in the function set_advance_front() 
*	above.	   
*
*
*/

/*ARGSUSED*/
LOCAL int reconstruct_advance_front3d(
	double		dt,
	double		*dt_frac,
	Front		*front,
	Front		**newfront,
	POINTER		wave)
{
	static const char *fname = "reconstruct_advance_front3d";
	boolean	   has_tracked_surfaces;
	int	   status;

	DEBUG_ENTER(reconstruct_advance_front3d)

	debug_print("front","Entered %s(step %d time %g dt %g)\n",fname,
				front->step,front->time,dt);
	debug_front("old_front","into advance front",front);

	*newfront = copy_front(front);
	(*newfront)->redis_flag = NO;
	Interface_redistributed(*newfront) = NO;
	has_tracked_surfaces = (front->interf->surfaces != NULL) ? YES : NO;
	if (pp_max_status(has_tracked_surfaces) == NO)
	{
	    set_size_of_intfc_state(size_of_state(front->interf));
	    set_copy_intfc_states(YES);
	    (*newfront)->interf = pp_copy_interface(front->interf);
	    status = ((*newfront)->interf != NULL) ? GOOD_STEP : ERROR_IN_STEP;
	    return return_advance_front(front,newfront,status,fname);
	}

	start_clock("copy_interface");

		/* Initialize Newfront */

	set_size_of_intfc_state(size_of_state(front->interf));
	set_copy_intfc_states(NO);
	(*newfront)->interf = pp_copy_interface(front->interf);
	if ((*newfront)->interf == NULL)
	{
	    status = ERROR_IN_STEP;
	    (void) printf("ERROR in advance_3d_front(), "
		          "unable to copy interface\n");
	    return return_advance_front(front,newfront,status,fname);
	}
	stop_clock("copy_interface");

		/* Propagate points on surfaces */

	start_clock("propagate");
	set_copy_intfc_states(YES);
	status = propagate_3d_front(wave,front,*newfront,dt,dt_frac,YES);
	debug_front("final_front","after propagate_3d_front front:",*newfront);
	stop_clock("propagate");

	if (status != GOOD_STEP)
	{
	    (void) printf("WARNING in %s(), "
	                  "propagate_3d_front() failed\n",fname);
	    return return_advance_front(front,newfront,status,fname);
	}
	interface_reconstructed((*newfront)->interf) = YES;
	return return_advance_front(front,newfront,status,fname);
}		/*end reconstruct_advance_front3d*/

LOCAL 	int 	hybrid_advance_front3d(
	double	dt,
	double	*dt_frac,
	Front	*front,
	Front	**newfront,
	POINTER	wave)
{
	int	step_status;
	double	dt_frac_restore = *dt_frac;
	DEBUG_ENTER(hybrid_advance_front3d)

	step_status = preserve_advance_front3d(dt,dt_frac,front,
						     newfront,wave);
	if (debugging("hybrid"))
	{
	    (void) printf("HYBRID  preserve_advance_front(): ");
	    print_time_step_status("time step status = ",step_status,"\n");
	}

	if (step_status != GOOD_STEP)
	{
	    *dt_frac = dt_frac_restore;
	    step_status = reconstruct_advance_front3d(dt,dt_frac,front,
							    newfront,wave);
	}

	DEBUG_LEAVE(hybrid_advance_front3d)
	return step_status;
}		/*end hybrid_advance_front3d*/

LOCAL int mixed_advance_front3d(
	double		dt,
	double		*dt_frac,
	Front		*front,
	Front		**newfront,
	POINTER		wave)
{
	int status;
	double   dt_frac_restore = *dt_frac;
	DEBUG_ENTER(mixed_advance_front3d)

	if (debugging("mixed"))
	{
	    (void) printf("TIMESTEP = %d\n",front->step);
	    (void) printf("MIXED  Redistribution_count(front) = %d\n",
			  Redistribution_count(front));
	    (void) printf("MIXED  Frequency_of_reconstruction(front) = %d\n",
			  Frequency_of_reconstruction(front));
	}
	status = preserve_advance_front3d(dt,dt_frac,
				    front,newfront,wave);

	/*do not use grid based */

	DEBUG_LEAVE(mixed_advance_front3d)
	return status;
}		/*end mixed_advance_front3d*/

LOCAL int preserve_advance_front3d(
	double		dt,
	double		*dt_frac,
	Front		*front,
	Front		**newfront,
	POINTER		wave)
{
	static const char *fname = "preserve_advance_front3d";
	boolean	   has_tracked_surfaces;
	int	   status;
	boolean	   do_redist = YES;
	int	   i;

	DEBUG_ENTER(preserve_advance_front3d)

	debug_print("front","Entered %s(step %d time %g dt %g)\n",fname,
				front->step,front->time,dt);
	debug_front("old_front","into advance front",front);

	if (debugging("trace"))
	    (void) printf("Entering preserve_advance_front3d()\n");
	do_redist = NO;
	if (front->step % Frequency_of_redistribution(front,GENERAL_WAVE) == 0)
	    do_redist = YES;

	*newfront = copy_front(front);
	(*newfront)->redis_flag = front->redis_flag;
	Interface_redistributed(*newfront) = NO;
	has_tracked_surfaces = (front->interf->surfaces != NULL) ? YES : NO;
	if (pp_max_status(has_tracked_surfaces) == NO)
	{
	    set_size_of_intfc_state(size_of_state(front->interf));
	    set_copy_intfc_states(YES);
	    (*newfront)->interf = pp_copy_interface(front->interf);
	    status = ((*newfront)->interf != NULL) ? GOOD_STEP : ERROR_IN_STEP;
	    return return_advance_front(front,newfront,status,fname);
	}

	start_clock("copy_interface");

		/* Initialize Newfront */

	set_size_of_intfc_state(size_of_state(front->interf));
	set_copy_intfc_states(NO);
	(*newfront)->interf = pp_copy_interface(front->interf);
	if ((*newfront)->interf == NULL)
	{
	    (void) printf("WARNING in advance_3d_front(), "
	                  "unable to copy interface\n");
	    return return_advance_front(front,newfront,ERROR_IN_STEP,fname);
	}
	stop_clock("copy_interface");

		/* Propagate points on surfaces */

	start_clock("propagate");
	set_copy_intfc_states(YES);

	status = propagate_3d_front(wave,front,*newfront,dt,dt_frac,NO);

	if (status != GOOD_STEP)
	{
	    (void) printf("WARNING in preserve_advance_front3d(), "
			  "status from propagate_3d_front() != GOOD_STEP\n");
	    print_time_step_status("time step status = ",status,"\n");
	    return return_advance_front(front,newfront,status,fname);
	}

	debug_front("final_front","after scatter front:",*newfront);
	stop_clock("propagate");

		/* Redistribute the New Front */

	interpolate_intfc_states((*newfront)->interf) = YES;
	switch (redistribute(*newfront,do_redist,NO)) 
	{
	case GOOD_REDISTRIBUTION:
	    Redistribution_count(front) = Redistribution_count(*newfront);
	    status = GOOD_STEP;
	    break;
	
	case UNABLE_TO_UNTANGLE:
	    *dt_frac *= TIME_STEP_REDUCTION_FACTOR(front->interf);
	    status = MODIFY_TIME_STEP;
	    (void) printf("WARNING in %s(), UNABLE_TO_UNTANGLE, "
	                  "redistribution of front failed\n",fname);
	    (void) printf("dt_frac set to %g\n",*dt_frac);
	    break;

	case MODIFY_TIME_STEP_REDISTRIBUTE: /* Not supported in 3d */
	    *dt_frac *= TIME_STEP_REDUCTION_FACTOR(front->interf);
	    status = MODIFY_TIME_STEP;
	    (void) printf("WARNING in %s(), "
			  "MODIFY_TIME_STEP_REDISTRIBUTE, ",fname);
	    (void) printf("dt_frac set to %g\n",*dt_frac);
	    break;
		
	case BAD_REDISTRIBUTION:
	    *dt_frac = Min_time_step_modification_factor(front);
	    status = MODIFY_TIME_STEP;
	    (void) printf("WARNING in %s(), BAD_REDISTRIBUTION, ",fname);
	    (void) printf("dt_frac set to %g\n",*dt_frac);
	    break;
	
	case INCONSISTENT_RECONSTRUCTION:
	    status = REPEAT_TIME_STEP;
	    (void) printf("WARNING in %s(), INCONSISTENT_RECONSTRUCTION\n",
	    		   fname);
	    break;

	default:
	    (void) printf("WARNING in %s(), "
	                  "unknown redistribution status\n",fname);
	    debug_front("ERROR_front","after error",*newfront);
	    clean_up(ERROR);
	    break;
	}
	if (status == GOOD_STEP)
	    init_intfc_curvature3d(*newfront,(*newfront)->interf);
	if (debugging("final_front"))
	    print_Front_structure(front);
	
	if(status == REPEAT_TIME_STEP)
	{
	    if(recon_repeat() >= 1)
	    {
		printf("ERROR preserve_advance_front3d, "
		       "repeat time step twice.\n");
		clean_up(ERROR);
	    }
	    else
	    {
		set_repeat(1);
	    }
	}
	else
	    set_repeat(0);

	return return_advance_front(front,newfront,status,fname);
}		/*end preserve_advance_front3d*/

/*
*
*
*			propagate_3d_front():
*
*	   This function performs the actual propagation of the front
*	to the new front location.  It consists of both a normal 
*	propagate and a tangential propagate, both of which are followed
*	by a scatter front step; for a total of four steps altogether.
*
*	Returns one of the following:
*		GOOD_STEP
*		ERROR_IN_STEP  (only from failure in node prop)
*		MODIFY_TIME_STEP
*/

void	tecplot_interface_states(const char *, INTERFACE *);

void    g_check_params_comp(POINTER, INTERFACE*);

LOCAL int propagate_3d_front(
	POINTER		wave,
	Front		*front,
	Front		*newfront,
	double		dt,
	double		*dt_frac,
	boolean		reconstruct_front)
{
	INTERFACE	*intfc_old;
	double		V[MAXD];
	int             step_status, tangent_status;
	DEBUG_ENTER(propagate_3d_front)

	debug_print("front","Entered propagate_3d_front().\n");

	intfc_old = front->interf;

	/*after redistribute or restart, the curvature is not calculated. */
        init_intfc_curvature3d(front, front->interf);
	
	if (front->_point_propagate != NULL)
	{
	    boolean tri_list_status;
	    start_clock("normal_propagate");

	    tri_list_status = 
	        make_tri_comp_lists(intfc_old) == FUNCTION_FAILED ? NO : YES;

	    if (pp_min_status(tri_list_status) == NO)
	    {
	        stop_clock("normal_propagate");
	    	screen("ERROR in propagate_3d_front(), "
	    	       "make_tri_comp_lists() failed\n");
	    	clean_up(ERROR);
	    }
	    
	    /* Set Default Propagation Limits */

	    set_propagation_limits(front,newfront);

	    start_clock("init_propagate");
	    init_propagate(front);
	    stop_clock("init_propagate");

	    /*set_wall_flag_for_surface(front->interf); */

	    if (front->_point_propagate != NULL)
                propagate_surface_points(front,newfront,wave,dt,V);
	    else if (front->surface_propagate != NULL)
	        surface_propagate(front,newfront,wave,dt,V);
	    
	    if (front->curve_propagate != NULL)
	        propagate_curve_points(front,newfront,wave,dt);
	
	    if (front->node_propagate != NULL)
	    {
	        step_status = propagate_node_points(front,newfront,
						    wave,dt,dt_frac);
		if (step_status != GOOD_STEP)
		{
		    DEBUG_LEAVE(propagate_3d_front)
		    return step_status;
		}
	    }
	    init_intfc_curvature3d(newfront,newfront->interf);
	    debug_front("np_front","after normal propagation",newfront);
	    stop_clock("normal_propagate");
	}
	if (*front->max_scaled_propagation > 0.5)
        {
            (void) printf("WARNING in advance_front3d(), "
                          "front->max_scaled_propagation = %f\n",
                          *front->max_scaled_propagation);
            *dt_frac = 0.4/(*front->max_scaled_propagation);
            step_status = MODIFY_TIME_STEP;
	    return step_status;
        }
	
	debug_propagate_3d_front(newfront);

	interface_reconstructed(newfront->interf) = NO;
	/*prev_interface(newfront->interf) = front->interf; */

	if (reconstruct_front == YES)
	{
	    boolean reconstruct_status;
	    start_clock("reconstruct_front");
	    
	    reconstruct_status = reconstruct_front_at_grid_crossing(newfront);
	    
	    reconstruct_status = pp_min_status(reconstruct_status);
	    if (!reconstruct_status)
	    {
	    	stop_clock("reconstruct_front");
	    	(void) printf("WARNING in propagate_3d_front(), "
	    	              "reconstruct_front failed, "
	    	              "MODIFY_TIME_STEP\n");
	    	*dt_frac *= TIME_STEP_REDUCTION_FACTOR(front->interf);
		(void) printf("dt_frac = %g after scaling by "
			      "TIME_STEP_REDUCTION_FACTOR(front->interf)"
			      " = %g\n",*dt_frac,
			      TIME_STEP_REDUCTION_FACTOR(front->interf));
		DEBUG_LEAVE(propagate_3d_front)
		return MODIFY_TIME_STEP; 
	    }
	    stop_clock("reconstruct_front");
	}

	start_clock("scatter_front");
	if (!scatter_front(newfront))
	{
	    stop_clock("scatter_front");
	    (void) printf("WARNING in propagate_3d_front(), "
	                  "1st scatter_front() failed\n"
	                  "MODIFY_TIME_STEP\n");  
	    *dt_frac *= TIME_STEP_REDUCTION_FACTOR(front->interf);
	    DEBUG_LEAVE(propagate_3d_front)
	    return MODIFY_TIME_STEP; 
	}
	stop_clock("scatter_front");

        init_intfc_curvature3d(newfront,newfront->interf);
	if (front->_tan_point_propagate != NULL)
	{
	    start_clock("tangentiall");
	    
            front->tan_sec = NO;
	    tangent_status = propagate_points_tangentially(front,newfront,
							   reconstruct_front,
							   dt,dt_frac,
							   (front->step)%2);
	    if (tangent_status != GOOD_STEP)
	    {
		DEBUG_LEAVE(propagate_3d_front)
		return tangent_status;
	    }
	    front->tan_sec = YES;
	    tangent_status = propagate_points_tangentially(front,newfront,
							   reconstruct_front,
							   dt,dt_frac,
							   (front->step+1)%2);
	    if (tangent_status != GOOD_STEP)
	    {
		DEBUG_LEAVE(propagate_3d_front)
		return tangent_status;
	    }

	    stop_clock("tangentiall");
	}
	else if (front->tan_surface_propagate != NULL)
	{
	    INTERFACE *intfc = front->interf;
	    INTERFACE *new_intfc = newfront->interf;
	    SURFACE **olds,**news;

	    start_clock("tangentiall");
	    for (olds = intfc->surfaces, news = new_intfc->surfaces; 
			olds && *olds; ++olds, ++news)
            {
                (*front->tan_surface_propagate)(front,newfront,
                                              intfc,*olds,*news,dt);
            }
	    stop_clock("tangentiall");
	}

	debug_print("front","Left propagate_3d_front()\n");
	DEBUG_LEAVE(propagate_3d_front)
	return GOOD_STEP; 
}		/*end propagate_3d_front*/

EXPORT void   init_intfc_curvature3d(
       Front               *front,
       INTERFACE           *intfc)
{
       HYPER_SURF              *hs;
       HYPER_SURF_ELEMENT      *hse;
       POINT                   *p;
     
      
	intfc->normal_unset = YES;
	reset_normal_on_intfc(intfc);
       	(void) next_point(intfc,NULL,NULL,NULL);
       	while (next_point(intfc,&p,&hse,&hs))
           p->curvature = mean_curvature_at_point(p,hse,hs,front);
       
       	/*fix 3 comp curve case and also a nan bug */
       	(void) next_point(intfc,NULL,NULL,NULL);
       	while (next_point(intfc,&p,&hse,&hs))
       	{
	   if(Boundary_point(p))
	   {
               p->curvature = 0.0;
	       continue;
	   }
	   
	   if(isnan(p->curvature))
	   {
		int	i, nt;
		TRI	**ptris;

		/*in this case, triangle with 0 area may appear. */
		printf("WARNING init_intfc_curvature3d, "
	       	       "p->curvature=%24.16e, set to 0.0\n", p->curvature);
		print_general_vector("p=", Coords(p), 3, "\n");
		nt = set_tri_list_around_point(p,Tri_of_hse(hse),&ptris,intfc);
		p->curvature = 0.0;
		for(i=0; i<nt; i++)
		    print_tri(ptris[i], intfc);
		tecplot_tris("cur_tst", ptris, nt);
	   }
        }
	intfc->normal_unset = NO;
        intfc->curvature_unset = NO;
}      /* end init_intfc_curvature3d */

	  
/*
*			propagate_points_tangentially():
*
*	This is the main driver program for the tangential sweep in 3D.
*	The procedure is to copy the incoming interface of newfront
*	into intfc_new, perform a loop over the two (geometrically identical)
*	interfaces intfc_tmp = newfront->interf,  and intfc_new and install
*	the output from the pointwise tan_point_propagate function into the
*	states on intfc_new.  At the end of the loop the original
*	newfront->interf = intfc_tmp is deleted and intfc_new is installed
*	as newfront->interf.  Finally the updated interface is marked as
*	reconstructed (since it is a copy of a reconstructed interface), and
*	the new front states are scattered across processors.
*
*	NOTE: The current tangential sweep function, f_tan_point_propagate,
*	implicitly assumes that the triangles about a given point project onto
*	the tangent plane defined by the discrete normal vector in such a
*	way as to produce a regular non-intersecting polygon.  This assumption
*	can be violated when a point is the vertex of a cone-like object.
*	Whether this occurs or not depends on the choice of the normal uni_array.
*	Internally, f_tan_point_propagate will attempt to compute an alternative
*	normal vector that does not produce a tangled projection.  If this is
*	not possible,  then f_tan_point_propagate will return a failure status.
*	If this occurs propagate_points_tangentially will remove the
*	offending point,  rescatter the modified interfaces,  recopy intfc_new,
*	and restart the tangential sweep loop.  Needless to say,  this is a
*	rather inefficient solution to the problem.  The hope is that such
*	irregular cases will be rare and this problem will only occur
*	infrequently in a given run.  Thus the goal is to have a robust,  if
*	costly,  way of treating such geometric degeneracies.
*/

LOCAL int propagate_points_tangentially(
	Front		*front,
	Front		*newfront,
	int		reconstruct_front,
	double		dt,
	double		*dt_frac,
	int		dir)
{
	INTERFACE               *intfc_tmp, *intfc_new;
	HYPER_SURF		*tmphs;
	HYPER_SURF_ELEMENT 	*tmphse;
	SURFACE                 **s_tmp, **s_new;
	POINT			*tmpp, *newp;
	TRI                     *tri_tmp, *tri_new;
	int                     i;
	boolean                    modified;
	DEBUG_ENTER(propagate_points_tangentially)

	start_clock("copy_interface");
	intfc_tmp = newfront->interf;
	set_size_of_intfc_state(size_of_state(front->interf));
	if ((intfc_new = pp_copy_interface(intfc_tmp)) == NULL)
	{
	    screen("ERROR in propagate_points_tangentially(), "
		   "unable to copy interface\n");
	    clean_up(ERROR);
	}
	stop_clock("copy_interface");
	if (debugging("consistency") && !consistent_interface(intfc_new))
	{
	    screen("ERROR in propagate_points_tangentially(), "
		   "intfc_new is inconsistent\n");
	    clean_up(ERROR);
	}
	
	start_clock("tan_propagate");

	modified = NO;
	/* Reset sort status for points on intfc_tmp and intfc_new */
	(void) next_point(intfc_tmp,NULL,NULL,NULL);
	(void) next_point(intfc_new,NULL,NULL,NULL);
	for (s_tmp = intfc_tmp->surfaces, s_new = intfc_new->surfaces;
		s_tmp && *s_tmp && s_new && *s_new; ++s_tmp, ++s_new)
	{
	  tmphs = Hyper_surf(*s_tmp);
	  for (tri_tmp = first_tri(*s_tmp), tri_new = first_tri(*s_new);
	       !at_end_of_tri_list(tri_tmp,*s_tmp) &&
	       !at_end_of_tri_list(tri_new,*s_new);
	       tri_tmp = tri_tmp->next, tri_new = tri_new->next)
	  {
	    tmphse = Hyper_surf_element(tri_tmp);
	    for (i = 0; i < 3; ++i)
	    {
	      tmpp = Point_of_tri(tri_tmp)[i];
	      newp = Point_of_tri(tri_new)[i];
	      if (!sorted(tmpp) && !sorted(newp))
	      {
	        if (tan_point_propagate(front,tmpp,newp,tmphse,tmphs,dt,dir))
		  sorted(tmpp) = sorted(newp) = YES;
		else
		{
		  sorted(tmpp) = sorted(newp) = YES;
		}
	      }
	      else if (!sorted(tmpp) || !sorted(newp))
	      {
	        screen("ERROR in propagate_points_tangentially(), "
		       "point lists in intfc_tmp and intfc_new "
		       "inconsistent\n");
	        clean_up(ERROR);
	      }
	    }
	  }
	}
	stop_clock("tan_propagate");

	(void) delete_interface(intfc_tmp);
	newfront->interf = intfc_new;
	debug_front("tp_front","after tangential propagation",newfront);

	start_clock("scatter_front");
	start_clock("copy_interface");
	set_size_of_intfc_state(size_of_state(front->interf));
	if ((intfc_new = pp_copy_interface(newfront->interf)) == NULL)
	{
	    screen("ERROR in propagate_points_tangentially(), "
		   "unable to copy interface\n");
	    clean_up(ERROR);
	}
	stop_clock("copy_interface");

        init_intfc_curvature3d(newfront,newfront->interf);
	
	if(NO)
	{
	    char    s[50], sn[50];

	    sprintf(sn, "%s", right_flush(pp_mynode(),PP_NODE_FIELD_WIDTH));
	    sprintf(s, "curvature%s_%s", right_flush(newfront->step, 
	    			TSTEP_FIELD_WIDTH), sn);
	    /*tecplot_interface_states(s, newfront->interf);*/
	}
	
	for (i = 0; i < 2; ++i)
	{
	  if (scatter_front(newfront))
	      break;
	  else
	  {
	    printf("\n entered else part!");
	    delete_interface(newfront->interf);
	    newfront->interf = intfc_new;
	    start_clock("copy_interface");
	    set_size_of_intfc_state(size_of_state(front->interf));
	    if ((intfc_new = pp_copy_interface(newfront->interf)) == NULL)
	    {
	      screen("ERROR in propagate_points_tangentially(), "
		     "unable to copy interface\n");
	      clean_up(ERROR);
	    }
	    stop_clock("copy_interface");
	    interface_reconstructed(newfront->interf) = NO;
	    if (i == 1)
	    {
	      if (redistribute(newfront,YES,NO) != GOOD_REDISTRIBUTION) 
	      {
		  i = 2;
		  break;
	      }
	    }
	  }
	}
	if (intfc_new != NULL)
	    delete_interface(intfc_new);
	if (i == 2)
	{
	  (void) printf("WARNING in propagate_points_tangentially(), "
			"2nd scatter_front() failed\n"
			"MODIFY_TIME_STEP_NODE\n");
	  *dt_frac *= TIME_STEP_REDUCTION_FACTOR(front->interf);
	  DEBUG_LEAVE(propagate_points_tangentially)
	  return MODIFY_TIME_STEP; 
	}
	stop_clock("scatter_front");
	DEBUG_LEAVE(propagate_points_tangentially)
	return GOOD_STEP;
}		/*end propagate_points_tangentially*/

/*
*			unfold_interface_section():
*
*	This function called when tan_point_propagate() fails.  The presumed
*	reason for this failure is the the interface near the point being
*	propagated is folded,  i.e. the projection of the interface onto
*	the tangent plane at the point is tangled.  This function assumes
*	that the interfaces corresponding to s_tmp and s_new are copies of
*	one another so that the triangle and points lists are in exact
*	one-to-one correspondence with IDENTICAL positions and that the
*	surface tri lists are synchronized.  The procedure is to perform a
*	local redistribution of the interface sections about the points.
*	The procedure is as follows:
*
*	1. Search for very small tris about the points,  merge this tris with
*	   their neighbors.
*	2. If no small tris are found,  then delete the points tmpp and newp
*	   from their respective interfaces.
*
*/

LOCAL	void unfold_interface_section(
	POINT   *tmpp,
	POINT   *newp,
	TRI     *tri_tmp,
	TRI     *tri_new,
	SURFACE *s_tmp,
	SURFACE *s_new)
{
	INTERFACE    *intfc_tmp = s_tmp->interface;
	INTERFACE    *intfc_new = s_new->interface;
	TRI          **tris;
	boolean         repeat, modified;
	double        lbar, lmin;
	int          i, imin, v;
	int          nt, nt_tmp, nt_new;
	static TRI   **tris_tmp = NULL, **tris_new = NULL;
	static double *l;
	static int   max_num_tris = 0;
	const double  ltol = 0.01;

	debug_print("unfold","Entered unfold_interface_section\n");
	if (debugging("unfold") && debugging("consistency"))
	{
	    if (!consistent_interface(intfc_tmp))
	    {
	        screen("ERROR in unfold_interface_section(), "
		       "intfc_tmp is inconsistent\n");
	        clean_up(ERROR);
	    }
	    if (!consistent_interface(intfc_new))
	    {
	        screen("ERROR in unfold_interface_section(), "
		       "intfc_new is inconsistent\n");
	        clean_up(ERROR);
	    }
	}
	modified = NO;
	do
	{
	  repeat = NO;
	  nt_tmp = set_tri_list_around_point(tmpp,tri_tmp,&tris,intfc_tmp);
	  if (nt_tmp > max_num_tris)
	  {
	    if (tris_tmp)
	      free_these(3,tris_tmp,tris_new,l);
	    max_num_tris = 2*nt_tmp+1;
	    uni_array(&tris_tmp,max_num_tris,sizeof(TRI*));
	    uni_array(&tris_new,max_num_tris,sizeof(TRI*));
	    uni_array(&l,max_num_tris,FLOAT);
	  }
	  for (i = 0; i < nt_tmp; ++i)
	    tris_tmp[i] = tris[i];
	  tris_tmp[nt_tmp] = NULL;
	  nt_new = set_tri_list_around_point(newp,tri_new,&tris,intfc_new);
	  if (nt_new != nt_tmp)
	  {
	    screen("ERROR in unfold_interface_section(), "
		   "inconsisent tri lists on tmp and new interfaces\n");
	    clean_up(ERROR);
	  }
	  nt = nt_tmp;
	  for (i = 0; i < nt; ++i)
	    tris_new[i] = tris[i];
	  tris_new[nt] = NULL;

	  for (i = 0, lbar = 0.0; i < nt; ++i)
	  {
	    v = Vertex_of_point(tris_tmp[i],tmpp);
	    l[i] = length_of_tri_side(tris_tmp[i],Next_m3(v));
	    lbar += l[i];
	  }
	  lbar /= nt;
	  if (debugging("unfold"))
	  {
	      (void) printf("lbar = %"FFMT"\n",lbar);
	      for (i = 0; i < nt; ++i)
	      {
	          (void) printf("l[%d] = %"FFMT", l[%d]/lbar = %"FFMT"\n",
				i,l[i],i,l[i]/lbar);
	      }
	  }

	  imin = -1;
	  lmin = HUGE_VAL;
	  for (i = 0; i < nt; ++i)
	  {
	    if (l[i] < ltol*lbar)
	    {
	      if (l[i] < lmin)
	      {
		lmin = l[i];
		imin = i;
	      }
	    }
	  }
	  if (imin != -1)
	  {
	    if (debugging("unfold"))
	    {
	      (void) printf("deleting tri_tmp(%llu) side %d\n",
			    (long long unsigned int)tri_number(tri_tmp,intfc_tmp),imin);
	      print_tri(tri_tmp,intfc_tmp);
	      (void) printf("deleting tri_new(%llu) side %d\n",
			    (long long unsigned int)tri_number(tri_new,intfc_new),imin);
	      print_tri(tri_new,intfc_new);
	    }
	    if (!delete_side_of_tri(tris_new[imin],s_new,imin) ||
	        !delete_side_of_tri(tris_tmp[imin],s_tmp,imin))
	    {
		screen("ERROR in unfold_interface_section(), "
			"delete_side_of_tri() failed\n");
		clean_up(ERROR);
	    }
	    repeat = YES;
	    modified = YES;
	  }
	} while(repeat);
	
	if (!modified)
	{
	  /* Are the tri lists closed? */
	  if ((nt_tmp >= 3) &&
	      (Next_tri_at_vertex(tris_tmp[0],tmpp) == tris_tmp[nt_tmp-1]) &&
	      (Prev_tri_at_vertex(tris_tmp[nt_tmp-1],tmpp) == tris_tmp[0]) &&
	      (Next_tri_at_vertex(tris_new[0],newp) == tris_new[nt_new-1]) &&
	      (Prev_tri_at_vertex(tris_new[nt_new-1],newp) == tris_new[0]))
	  {
	    if (debugging("unfold"))
	    {
		(void) printf("deleting vertex tmpp(%llu) of tri_tmp(%llu)\n",
			      (long long unsigned int)point_number(tmpp),
			      (long long unsigned int)tri_number(tri_tmp,intfc_tmp));
		(void) printf("deleting vertex newp(%llu) of tri_new(%llu)\n",
			      (long long unsigned int)point_number(newp),
			      (long long unsigned int)tri_number(tri_new,intfc_tmp));
	    }
	    if (!delete_vertex_of_tri(tmpp,tri_tmp,s_tmp) ||
	        !delete_vertex_of_tri(newp,tri_new,s_new))
	    {
	      screen("ERROR in unfold_interface_section(), "
		     "delete_vertex_of_tri() failed\n");
	      clean_up(ERROR);
	    }
	  }
	  else
	  {
	    double BBL[3], BBU[3];

	    set_tri_list_bounding_box(tris_tmp,nt_tmp,BBL,BBU,NO,YES);
	    set_tri_list_bounding_box(tris_new,nt_new,BBL,BBU,YES,YES);

	    screen("ERROR in unfold_interface_section(), can't unfold\n");
	    (void) printf("Boundary_point(tmpp) = %d\n",Boundary_point(tmpp));
	    (void) printf("nt_tmp = %d\n",nt_tmp);
	    for (i = 0; i < nt_tmp; ++i)
	    {
	      (void) printf("tris_tmp[%d]\n",i);
	      print_tri(tris_tmp[i],intfc_tmp);
	    }
	    gview_plot_triangle_list("","tris_tmp",tris_tmp,nt_tmp,
		                     0.1,0.0,0.0,0.9,0.0,0.0,0.5,BBL,BBU);
	    gview_plot_vertices("","tmpp",&tmpp,1,BBL,BBU);

	    (void) printf("Boundary_point(newp) = %d\n",Boundary_point(newp));
	    (void) printf("nt_new = %d\n",nt_new);
	    for (i = 0; i < nt_new; ++i)
	    {
	      (void) printf("tris_new[%d]\n",i);
	      print_tri(tris_new[i],intfc_new);
	    }
	    gview_plot_triangle_list("","tris_new",tris_new,nt_new,
		                     0.1,0.0,0.0,0.9,0.0,0.0,0.5,BBL,BBU);
	    gview_plot_vertices("","newp",&newp,1,BBL,BBU);
	    clean_up(ERROR);
	  }
	}
	debug_print("unfold","Left unfold_interface_section\n");
}		/*end unfold_interface_section*/

#define   MAX_NUM_CURVE   500

LOCAL	void detach_and_propagate_curves(
	Front	*front,
	Front	*newfront,
	POINTER wave,
	double	dt)
{
	INTERFACE    *sav_intfc = current_interface();
	SURFACE	     **s, *olds, *news;
	CURVE	     **new_curves[MAX_NUM_CURVE];
	CURVE	     **c, **cc, *newc, *pc;
	P_LINK	     *p_table;	
	int	     i, num_c, p_size;
	boolean	     found = NO;

	DEBUG_ENTER(detach_and_propagate_curves)

	/*store all paired curves */
	i = 0;
	for (c =  front->interf->curves, 
	     cc = newfront->interf->curves;
	     c && *c &&
	     cc && *cc; 
	     c++, cc++, i++)
	{
	    if(i >= MAX_NUM_CURVE)
	    {
		printf("ERROR detach_and_propagate_curves, too many curves.\n");
		clean_up(ERROR);
	    }
	    
	    new_curves[i] = NULL;
	    if(!add_to_pointers(*cc, &new_curves[i]))
	    {
	        printf("ERROR detach_and_propagate_curves, "
		       "can not add into new_curves.\n");
		clean_up(ERROR);
	    }
	}

	/*detach newfront surface */
	set_current_interface(newfront->interf);
	
	for(s=newfront->interf->surfaces; s && *s; s++)
	{
	    /*ASSUME the moving surface type */
	    if (wave_type(*s) >= FIRST_PHYSICS_WAVE_TYPE&&
                wave_type(*s) != ELASTIC_BOUNDARY)
	    {
	        olds = *s;
		news = detach_one_surface(olds);
	        if(olds != news)
		    found = YES;
		break;
	    }
	}

	if(!found)
	{
	    for(i=0, c=front->interf->curves; c && *c; c++, i++)
	        curve_propagate(front, wave, *c, new_curves[i][0], dt);
	    
	    set_current_interface(sav_intfc);
	    DEBUG_LEAVE(detach_and_propagate_curves)
	    return;
	}
	
	/*printf("#detach surfce  %d to %d.\n", olds, news); */
	/*construct the olds and news pairs */
	num_c = 0;
	for(i = 0; i < 2; i++)
	{
	    c = i==0 ? news->pos_curves : news->neg_curves;
	    for(; c && *c; c++)
	        num_c++;
	}

	p_size = 4*num_c + 1;
	uni_array(&p_table,p_size,sizeof(P_LINK));
	reset_hash_table(p_table,p_size);
	
	for(i = 0; i < 2; i++)
	{
	    c =  i==0 ? olds->pos_curves : olds->neg_curves; 
	    cc = i==0 ? news->pos_curves : news->neg_curves; 

	    for( ; c  && *c && 
	           cc && *cc;  c++, cc++)
	    {
	        add_to_hash_table((POINTER)(*c),(POINTER)(*cc), p_table,p_size);
	    }
	}

	i = 0;
	for(c = front->interf->curves; c && *c; c++, i++)
	{
	    newc = new_curves[i][0];
	    pc = (CURVE*) find_from_hash_table((POINTER)newc,p_table,p_size);
	    if(pc != NULL)
	    {
	        if(curve_type(pc) != NEUMANN_CURVE_P)
		    delete_from_pointers(newc, &new_curves[i]);
		add_to_pointers(pc, &new_curves[i]);
	    }
	}
	
	free(p_table);
	/*delete the previous attached surface and associated curves */
	delete_scn(olds);

	if(debugging("detach_and_propagate_curves"))
	{
	    printf("#old curves\n");
	    for(c = front->interf->curves; c && *c; c++)
	        printf("#%p  %d\n", (void*)*c, curve_type(*c));
	
	    printf("#new curves\n");
	    for(c = newfront->interf->curves; c && *c; c++)
	        printf("#%p  %d\n", (void*)*c, curve_type(*c));

	    printf("#prop curves\n");
	    for(i=0, c = front->interf->curves; c && *c; c++, i++)
	    {
	        printf("#c  %p  %3d   ", (void*)*c, curve_type(*c));
	        for(cc = new_curves[i]; cc && *cc; cc++)
		    printf("|cc %p %3d  ", (void*)*cc, curve_type(*cc));
		printf("|\n");
	    }
	}

	set_current_interface(sav_intfc);

	for(i=0, c = front->interf->curves; c && *c; c++, i++)
	{
	    for(cc = new_curves[i]; cc && *cc; cc++)
	        curve_propagate(front, wave, *c, *cc, dt);
	}

	DEBUG_LEAVE(detach_and_propagate_curves)
}		/*end propagate_curve_points*/


LOCAL	void propagate_curve_points(
	Front	*front,
	Front	*newfront,
	POINTER wave,
	double	dt)
{
  	CURVE	**new_curves;
	
	DEBUG_ENTER(propagate_curve_points)
	start_clock("curve_propagate");

	detach_and_propagate_curves(front, newfront, wave, dt);

	order_interface(newfront->interf);
        for(new_curves = newfront->interf->curves; 
	    new_curves && *new_curves; new_curves++)
            reorder_curve_link_list(*new_curves);

	stop_clock("curve_propagate");
	
	DEBUG_LEAVE(propagate_curve_points)
}		/*end propagate_curve_points*/

LOCAL 	int	propagate_node_points(
	Front	*front,
	Front	*newfront,
	POINTER	wave,
	double	dt,
	double	*dt_frac)
{
	NODE		**on, *oldn, **nn, *newn;
	NODE_FLAG	flag;
	int		step_status, node_status;
	RPROBLEM	*rp;
	DEBUG_ENTER(propagate_node_points)

	start_clock("node_propagate");
	set_to_next_node_only(flag);
	step_status = GOOD_STEP;
	for (on = front->interf->nodes, nn = newfront->interf->nodes;
	     on && *on && nn && *nn && step_status==GOOD_STEP;
	     ++on , ++nn)
	{
	    oldn = *on;
	    newn = *nn;
	    node_status = (*front->node_propagate)(front,wave,oldn,newn,
						   &rp,dt,dt_frac,
						   flag,NULL);
	    switch(node_status)	
	    {
	    case GOOD_NODE:
	        break;
	    case ERROR_NODE:
	    default:
	        print_node_status("WARNING in propagate_node_points(), "
			          "node propagate failed with node_status ",
		                  node_status,"\n");
		step_status = ERROR_IN_STEP;
		break;
	    }
	}
	stop_clock("node_propagate");
	if (front->pp_grid)
	    step_status = syncronize_time_step_status(step_status,front->pp_grid);
	DEBUG_LEAVE(propagate_node_points)
	return step_status; 
}		/*end propagate_node_points*/
	   

LOCAL void debug_propagate_3d_front(
	Front 	*newfront)
{
	if (debugging("propagation_3d"))
	    gview_plot_interface("final_intfc_bs",newfront->interf);

	if (debugging("check_intersect"))
	{
	    CROSS  *cross;
	    boolean intersection_status;

	    intersection_status = intersections(newfront->interf,&cross,YES);
	    if (pp_min_status(intersection_status) == NO)
	    {
	    	screen("ERROR in debug_propagate_3d_front(), "
		       "intersections() failed.\n");
	    	gview_plot_interface("final_intfc_bs",newfront->interf);
	    	clean_up(ERROR);
	    }
	    intersection_status = (cross != NULL) ? NO : YES;
	    if (pp_min_status(intersection_status) == NO)
	    {
	    	screen("ERROR in debug_propagate_3d_front(), "
	    	       "before fscatter\n");
	    	gview_plot_interface("final_intfc_bs",newfront->interf);
	    	clean_up(ERROR);
	    }
	}

	if (debugging("check_tri_pt"))
	{
	    (void) printf("after propagating and before scatter_front()\n");
	    if (!consistent_interface(newfront->interf))
	    {
		screen("ERROR in debug_propagate_3d_front() "
		       "inconsistent interface pointers\n");
		clean_up(ERROR);
	    }
	}
}		/*end debug_propagate_3d_front*/

LOCAL void propagate_surface_points(
        Front           *front,
        Front           *newfront,
        POINTER         wave,
        double           dt,
        double           *V)
{
        INTERFACE               *intfc_old = front->interf;
        INTERFACE               *intfc_new = newfront->interf;
        HYPER_SURF              *oldhs, *newhs;
        HYPER_SURF_ELEMENT      *oldhse, *newhse;
        POINT                   *oldp, *newp;
        DEBUG_ENTER(propagate_surface_points)

        start_clock("propagate_surface_points");

	(void) next_point(intfc_old,NULL,NULL,NULL);
        (void) next_point(intfc_new,NULL,NULL,NULL);
        while (next_point(intfc_old,&oldp,&oldhse,&oldhs) &&
             next_point(intfc_new,&newp,&newhse,&newhs))
        {
	    /*if(Boundary_point(newp))*/
	    if (Boundary_point(newp) && wave_type(oldhs) != 
			FIRST_PHYSICS_WAVE_TYPE)
	    {
		ft_assign(left_state(newp),left_state(oldp),front->sizest);
		ft_assign(right_state(newp),right_state(oldp),front->sizest);
	        continue;
	    }
	    point_propagate(front,wave,oldp,newp,oldhse,oldhs,dt,V);
	}

        stop_clock("propagate_surface_points");
        DEBUG_LEAVE(propagate_surface_points)
}               /*end propagate_surface_points*/

LOCAL int struct_advance_front3d(
	double		dt,
	double		*dt_frac,
	Front		*front,
	Front		**newfront,
	POINTER		wave)
{
	static const char *fname = "struct_advance_front3d";
	boolean	   has_tracked_surfaces;
	int	   status;
	double	   V[MAXD];

	if (debugging("trace"))
	    (void) printf("Entering struct_advance_front3d()\n");
	*newfront = copy_front(front);
	has_tracked_surfaces = (front->interf->surfaces != NULL) ? YES : NO;
	if (pp_max_status(has_tracked_surfaces) == NO)
	{
	    set_size_of_intfc_state(size_of_state(front->interf));
	    set_copy_intfc_states(YES);
	    (*newfront)->interf = pp_copy_interface(front->interf);
	    status = ((*newfront)->interf != NULL) ? GOOD_STEP : ERROR_IN_STEP;
	    return return_advance_front(front,newfront,status,fname);
	}

		/* Initialize Newfront */

	start_clock("copy_interface");
	set_size_of_intfc_state(size_of_state(front->interf));
	set_copy_intfc_states(NO);
	(*newfront)->interf = pp_copy_interface(front->interf);
	if ((*newfront)->interf == NULL)
	{
	    (void) printf("WARNING in advance_3d_front(), "
	                  "unable to copy interface\n");
	    return return_advance_front(front,newfront,ERROR_IN_STEP,fname);
	}
	stop_clock("copy_interface");

		/* Propagate points on surfaces */

	start_clock("propagate");

	start_clock("coupled_propagate");
	set_copy_intfc_states(YES);
	init_intfc_curvature3d(front, front->interf);
	set_propagation_limits(front,*newfront);

	if (front->_point_propagate != NULL)
	    propagate_surface_points(front,*newfront,wave,dt,V);
	if (front->curve_propagate != NULL)
	    propagate_curve_points(front,*newfront,wave,dt);
	if (front->node_propagate != NULL)
	    propagate_node_points(front,*newfront,wave,dt,dt_frac);
	stop_clock("coupled_propagate");

	start_clock("scatter_front");
	if (!scatter_front(*newfront))
	{
	    stop_clock("scatter_front");
            (void) printf("ERROR in struct_advance_front3d(), "
                          "scatter_front() failed\n");
	    stop_clock("scatter_front");
	    clean_up(ERROR);
	}
	else
	{
	    stop_clock("scatter_front");
	    status = GOOD_STEP;
	}

	start_clock("interior_propagate");
	init_intfc_curvature3d(*newfront,(*newfront)->interf);
	if (front->interior_propagate != NULL)
	    (*front->interior_propagate)(*newfront,dt);
	stop_clock("interior_propagate");

	start_clock("scatter_front");
	if (!scatter_front(*newfront))
	{
	    stop_clock("scatter_front");
            (void) printf("ERROR in struct_advance_front3d(), "
                          "scatter_front() failed\n");
	    stop_clock("scatter_front");
	    clean_up(ERROR);
	}
	else
	{
	    stop_clock("scatter_front");
	    status = GOOD_STEP;
	}

	stop_clock("propagate");

	if (status == GOOD_STEP)
	    init_intfc_curvature3d(*newfront,(*newfront)->interf);
	
	return return_advance_front(front,newfront,status,fname);
}	/* end struct_advance_front3d */

LOCAL int simple_advance_front3d(
	double		dt,
	double		*dt_frac,
	Front		*front,
	Front		**newfront,
	POINTER		wave)
{
	static const char *fname = "simple_advance_front3d";
	boolean	   has_tracked_surfaces;
	int	   status;
	double	   V[MAXD];
	INTERFACE *intfc_old,*intfc_new;

	if (debugging("trace"))
	    (void) printf("Entering simple_advance_front3d()\n");
	*newfront = copy_front(front);
	has_tracked_surfaces = (front->interf->surfaces != NULL) ? YES : NO;
	if (pp_max_status(has_tracked_surfaces) == NO)
	{
	    set_size_of_intfc_state(size_of_state(front->interf));
	    set_copy_intfc_states(YES);
	    (*newfront)->interf = pp_copy_interface(front->interf);
	    status = ((*newfront)->interf != NULL) ? GOOD_STEP : ERROR_IN_STEP;
	    return return_advance_front(front,newfront,status,fname);
	}

		/* Initialize Newfront */

	start_clock("copy_interface");
	set_size_of_intfc_state(size_of_state(front->interf));
	set_copy_intfc_states(NO);
	(*newfront)->interf = pp_copy_interface(front->interf);
	if ((*newfront)->interf == NULL)
	{
	    (void) printf("WARNING in advance_3d_front(), "
	                  "unable to copy interface\n");
	    return return_advance_front(front,newfront,ERROR_IN_STEP,fname);
	}
	stop_clock("copy_interface");

		/* Propagate points on surfaces */

	intfc_old = front->interf;
	intfc_new = (*newfront)->interf;
	if (front->_surface_propagate != NULL)
	{
	    SURFACE **olds,**news;
	    for (olds = intfc_old->surfaces, news = intfc_new->surfaces;
		 olds && *olds; ++olds, ++news)
		(*front->_surface_propagate)(front,wave,*olds,*news,dt);
	}
	if (front->_curve_propagate != NULL)
	{
	    CURVE **oldc,**newc;
	    for (oldc = intfc_old->curves, newc = intfc_new->curves;
		 oldc && *oldc; ++oldc, ++newc)
		(*front->_curve_propagate)(front,wave,*oldc,*newc,dt);
	}
	if (front->_node_propagate != NULL)
	{
	    NODE **oldn,**newn;
	    for (oldn = intfc_old->nodes, newn = intfc_new->nodes;
		 oldn && *oldn; ++oldn, ++newn)
		(*front->_node_propagate)(front,wave,*oldn,*newn,dt);
	}

	start_clock("propagate");

	start_clock("interior_propagate");
	if (front->interior_propagate != NULL)
	    (*front->interior_propagate)(*newfront,dt);
	stop_clock("interior_propagate");

	stop_clock("propagate");
	status = GOOD_STEP;

	if (debugging("trace"))
	    (void) printf("Entering simple_advance_front3d()\n");
	return return_advance_front(front,newfront,status,fname);
}	/* end simple_advance_front3d */

LOCAL int spherical_advance_front3d(
	double		dt,
	double		*dt_frac,
	Front		*front,
	Front		**newfront,
	POINTER		wave)
{
	INTERFACE       *intfc_old;
	static const char *fname = "spherical_advance_front3d";
	double     V[MAXD];
	int status;

	if (debugging("trace"))
	    (void) printf("Entering spherical_advance_front3d()\n");

	intfc_old = front->interf;
	*newfront = copy_front(front);
	set_size_of_intfc_state(size_of_state(front->interf));
	set_copy_intfc_states(NO);

	start_clock("init_intfc_curvature3d");
	init_intfc_curvature3d(front,front->interf);
	stop_clock("init_intfc_curvature3d");

	start_clock("copy_interface");
	(*newfront)->interf = pp_copy_interface(front->interf);
        if ((*newfront)->interf == NULL)
        {
            (void) printf("WARNING in spherical_advance_front3d(), "
                          "unable to copy interface\n");
            return return_advance_front(front,newfront,ERROR_IN_STEP,fname);
        }
	stop_clock("copy_interface");

	if (front->_point_propagate != NULL)
        {
            boolean tri_list_status;
            start_clock("normal_propagate");

            tri_list_status =
                make_tri_comp_lists(intfc_old) == FUNCTION_FAILED ? NO : YES;

            if (pp_min_status(tri_list_status) == NO)
            {
                stop_clock("normal_propagate");
                screen("ERROR in propagate_3d_front(), "
                       "make_tri_comp_lists() failed\n");
                clean_up(ERROR);
            }
	    /* Set Default Propagation Limits */

            set_propagation_limits(front,*newfront);

            start_clock("init_propagate");
            init_propagate(front);
            stop_clock("init_propagate");

	    if (front->_point_propagate != NULL)
                propagate_surface_points(front,*newfront,wave,dt,V);
            else if (front->surface_propagate != NULL)
                surface_propagate(front,*newfront,wave,dt,V);
	    if (front->curve_propagate != NULL)
                propagate_curve_points(front,*newfront,wave,dt);
	    if (front->node_propagate != NULL)
            {
                status = propagate_node_points(front,*newfront,
                                                    wave,dt,dt_frac);
                if (status != GOOD_STEP)
                {
                    return return_advance_front(front,newfront,status,fname);;
                }
            }
	    start_clock("init_intfc_curvature3d");
	    init_intfc_curvature3d(*newfront,(*newfront)->interf);
	    stop_clock("init_intfc_curvature3d");
            stop_clock("normal_propagate");
        }
	
	resolve_front_collision(*newfront);

	if (debugging("trace"))
	    (void) printf("Leaving spherical_advance_front3d()\n");
}	/* end spherical_advance_front3d */

LOCAL int resolve_front_collision(
	Front *front)
{
	int i,j,num_spheres;
	SURFACE **s,**surfs;
	double **centers;
	INTERFACE *intfc = front->interf;
	int dim = Dimension(intfc);
	double min_dist = HUGE;

	printf("Entering resolve_front_collision()\n");
	num_spheres = 0;
	intfc_surface_loop(intfc,s)
	{
	    if (wave_type(*s) != ICE_PARTICLE_BOUNDARY)
		continue;
	    num_spheres++;
	}
	FT_MatrixMemoryAlloc((POINTER*)&centers,num_spheres,MAXD,
			sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&surfs,num_spheres,sizeof(SURFACE*));
	i = 0;
	intfc_surface_loop(intfc,s)
	{
	    if (wave_type(*s) != ICE_PARTICLE_BOUNDARY)
		continue;
	    surfs[i] = *s;
	    for (j = 0; j < dim; ++j)
	    	centers[i][j] = center_of_mass(*s)[j];
	    ++i;
	}
	/*
	printf("num_spheres = %d\n",num_spheres);
	for (i = 0; i < num_spheres; ++i)
	    printf("centers[%d] = %f %f %f\n",i,centers[i][0],centers[i][1],
					centers[i][2]);
	for (i = 0; i < num_spheres; ++i)
	    printf("radius[%d] = %f\n",i,spherical_radius(surfs[i]));
	*/
	for (i = 0; i < num_spheres; ++i)
	for (j = i+1; j < num_spheres; ++j)
	{
	    double d = distance_between_positions(centers[i],centers[j],dim);
	    d -= spherical_radius(surfs[i]);
	    d -= spherical_radius(surfs[j]);
	    /*
	    printf("touch_disntance[%d,%d] = %f\n",i,j,d);
	    */
	    if (min_dist > d) min_dist = d;
	}
	/*
	printf("min_dist = %f\n",min_dist);
	*/
	if (min_dist < 0.0)
	{
	    printf("Collision found, code needed!\n");
	    clean_up(0);
	}
	printf("Leaving resolve_front_collision()\n");
}	/* end resolve_front_collision */

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
*				example16.c:
*
*		User initialization example for Front Package:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*	
*	This example shows a circle in a multiple vortex field. It 
*	demonstrates the reversibility of the front tracking method.
*
*/

#include <FronTier.h>

	/*  Function Declarations */
static void test_propagate(Front*);
static int test_vortex_vel(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,
	                       HYPER_SURF*,double*);

char *in_name,*restart_state_name,*restart_name,*out_name;
boolean RestartRun;
int RestartStep;
boolean binary = YES;


int main(int argc, char **argv)
{
	static Front front;
	static RECT_GRID comp_grid;
	static F_BASIC_DATA f_basic;
	static LEVEL_FUNC_PACK level_func_pack;
	static VELO_FUNC_PACK velo_func_pack;
	VORTEX_PARAMS vortex_params; /* velocity function parameters */
	Locstate  sl;

	FT_Init(argc,argv,&f_basic);
	f_basic.dim = 2;

	/* Initialize basic computational data */

	f_basic.L[0] = 0.0;	f_basic.L[1] = 0.0;
	f_basic.U[0] = 1.0;	f_basic.U[1] = 1.0;
	f_basic.gmax[0] = 200;	f_basic.gmax[1] = 200;		//myex grid size
	f_basic.boundary[0][0] = f_basic.boundary[0][1] = PERIODIC_BOUNDARY;
	f_basic.boundary[1][0] = f_basic.boundary[1][1] = PERIODIC_BOUNDARY;
	f_basic.size_of_intfc_state = 0;

        in_name                 = f_basic.in_name;
        restart_state_name      = f_basic.restart_state_name;
        out_name                = f_basic.out_name;
        restart_name            = f_basic.restart_name;
        RestartRun              = f_basic.RestartRun;
        RestartStep             = f_basic.RestartStep;

        sprintf(restart_name,"%s/intfc-ts%s",restart_name,
                        right_flush(RestartStep,7));
        if (pp_numnodes() > 1)
            sprintf(restart_name,"%s-nd%s",restart_name, 
                                right_flush(pp_mynode(),4));

	FT_StartUp(&front,&f_basic);
	add_to_debug("free_end_node");

	if (!RestartRun)
	{
	    /* Initialize interface through level function */

	    level_func_pack.neg_component = 2;
	    level_func_pack.pos_component = 2;
	    level_func_pack.wave_type = FIRST_PHYSICS_WAVE_TYPE;

	    level_func_pack.func_params = NULL;
	    level_func_pack.func = NULL;

	    level_func_pack.num_points = 251;		//myex num points
	    level_func_pack.is_closed_curve = NO;

	    FT_MatrixMemoryAlloc((POINTER*)&level_func_pack.point_array,
	    			level_func_pack.num_points,
				2,sizeof(double));
	    int i;
	    for (i = 0; i < level_func_pack.num_points; ++i)
	    {
	    	double phi = i*PI/(double)(level_func_pack.num_points-1);
	    	level_func_pack.point_array[i][0] = 0.5 + 0.20*cos(phi);
	    	level_func_pack.point_array[i][1] = 0.35 + 0.20*sin(phi);
	    }

	    FT_InitIntfc(&front,&level_func_pack);
	    if (f_basic.dim < 3)
                FT_ClipIntfcToSubdomain(&front);
	}
	CURVE **c;
	for (c = front.interf->curves; c && *c; ++c)
	{
	    if (negative_component(*c) == positive_component(*c))
	    {
		node_type((*c)->start) = MONO_COMP_NODE;
		node_type((*c)->end) = MONO_COMP_NODE;
	    }
	}
	xgraph_2d_intfc("init_intfc",front.interf);

	/* Initialize velocity field function */

	vortex_params.dim = 2;
	vortex_params.type[0] = 'M';
	vortex_params.cos_time = 0;
	vortex_params.cen[0] = 0.5;
	vortex_params.cen[1] = 0.25;
	vortex_params.rad = 0.15;

	velo_func_pack.func_params = (POINTER)&vortex_params;
	velo_func_pack.func = test_vortex_vel;
	velo_func_pack.point_propagate = fourth_order_point_propagate;

	FT_InitVeloFunc(&front,&velo_func_pack);

	/* Propagate the front */

	test_propagate(&front);

	clean_up(0);
	return 0;
}

static  void test_propagate(
        Front *front)
{
        double CFL;
	VORTEX_PARAMS *vparams = (VORTEX_PARAMS*)front->vparams;
	char gname[100] = "immersed_intfc.gif";

	front->max_time = 1.0;
	front->max_step = 1000000;
	front->print_time_interval = 1.0;
	front->movie_frame_interval = 0.02;
	vparams->time = 0.5*front->max_time;

        CFL = Time_step_factor(front);
	Frequency_of_redistribution(front,GENERAL_WAVE) = 2;

	printf("CFL = %f\n",CFL);
	printf("Frequency_of_redistribution(front,GENERAL_WAVE) = %d\n",
		Frequency_of_redistribution(front,GENERAL_WAVE));

	if (!RestartRun)
	{
            FT_RedistMesh(front);
	    FT_ResetTime(front);

	    // Always output the initial interface.
	    FT_Save(front,out_name);
            FT_AddMovieFrame(front,out_name,binary);

	    // This is a virtual propagation to get maximum front 
	    // speed to determine the first time step.

	    FT_Propagate(front);
            FT_SetTimeStep(front);
	    FT_SetOutputCounter(front);
	}
	else
	{
            FT_SetOutputCounter(front);
	}

	FT_TimeControlFilter(front);

        for (;;)
        {
	    /* Propagating interface for time step dt */

	    FT_Propagate(front);
	    FT_AddTimeStepToCounter(front);

	    //Next time step determined by maximum speed of previous
	    //step, assuming the propagation is hyperbolic and
	    //is not dependent on second order derivatives of
	    //the interface such as curvature, and etc.

            FT_SetTimeStep(front);

	    /* Output section */

            printf("\ntime = %f   step = %5d   next dt = %f\n",
                        front->time,front->step,front->dt);
            fflush(stdout);

            if (FT_IsSaveTime(front))
		FT_Save(front,out_name);
            if (FT_IsMovieFrameTime(front))
                FT_AddMovieFrame(front,out_name,binary);

            if (FT_TimeLimitReached(front))
                    break;

	    FT_TimeControlFilter(front);
        }
        (void) delete_interface(front->interf);
}       /* end test_propagate */

/********************************************************************
 *	Sample (vortex) velocity function for the front    *
 ********************************************************************/

static int test_vortex_vel(
        POINTER params,
        Front *front,
        POINT *p,
        HYPER_SURF_ELEMENT *hse,
        HYPER_SURF *hs,
        double *vel)
{
        VORTEX_PARAMS *vortex_params;
        int i, dim;
        double coeff,coeff2,xtemp,ytemp;
        char *type;
        double *coords = Coords(p);
        vortex_params = (VORTEX_PARAMS*)params;
	double x,y,z;

	dim = vortex_params->dim;
        type = vortex_params->type;

        coeff2 = 1.0;
        if(vortex_params->time > 0 && front->time >= vortex_params->time)
            coeff2 = -1.0;
	x = coords[0];
	y = coords[1];
	if (dim == 3)
	    z = coords[2];

	if (dim == 2)
	{
	    switch(type[0])
	    {
	        case 'm':
		case 'M':
		    //four vortex work
		    xtemp = 4*PI*(x+0.5);
		    ytemp = 4*PI*(y+0.5);
		    vel[0] = coeff2*sin(xtemp)*sin(ytemp); 
		    vel[1] = coeff2*cos(xtemp)*cos(ytemp); 
		    //end four vortex work
		    break;
		case 's':
		case 'S':
		    //2D single vortex motion
		    vel[0] = -coeff2*sin(PI*x)*sin(PI*x)*sin(2*PI*y);
		    vel[1] = coeff2*sin(2*PI*x)*sin(PI*y)*sin(PI*y);
		    break;
		default:
		    screen("Undefined vortex type!\n");
		    clean_up(ERROR);
	    }
	}
	else if (dim == 3)
	{
	    switch(type[0])
	    {
	    case 'm':
	    case 'M':
		//double vortex work
	        vel[0] = coeff2*2*sin(PI*x)*sin(PI*x)*sin(2*PI*y)*sin(2*PI*z);
	        vel[1] = -coeff2*sin(2*PI*x)*sin(PI*y)*sin(PI*y)*sin(2*PI*z);
	        vel[2] = -coeff2*sin(2*PI*x)*sin(2*PI*y)*sin(PI*z)*sin(PI*z);
		break;
	    case 's':
	    case 'S':
		//shearing flow motion
		vel[0] = coeff2*sin(PI*x)*sin(PI*x)*sin(2*PI*y);
		vel[1] = -coeff2*sin(2*PI*x)*sin(PI*y)*sin(PI*y);
		vel[2] = coeff2*sqr(1-2*sqrt((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)));
	    default:
	        screen("Undefined time dependency type!\n");
	        clean_up(ERROR);
	    }
	}    
}       /* end vortex_vel */


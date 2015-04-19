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
*				example.c:
*
*		User initialization example for Front Package:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*/

#include <FronTier.h>

	/*  Function Declarations */
static void propagation_driver(Front*);
static double sphere_func(POINTER,double*);
static int test_curvature_vel(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,
			HYPER_SURF*,double*);

char *in_name,*restart_state_name,*restart_name,*out_name;
boolean RestartRun;
int RestartStep;

/********************************************************************
 *	Velocity function parameters for the front	 	    *
 ********************************************************************/

typedef struct {
	int dim;
	double coeff;
	double epsilon;
} TEST_CURV_PARAMS;

int main(int argc, char **argv)
{
	static Front front;
	static RECT_GRID comp_grid;
	static F_BASIC_DATA f_basic;
	static LEVEL_FUNC_PACK level_func_pack;
	static VELO_FUNC_PACK velo_func_pack;
	RECT_BOX_PARAMS rect_box_params;
	TEST_CURV_PARAMS curv_params; /* velocity function parameters */
	Locstate  sl;

	f_basic.dim = 3;	
	FT_Init(argc,argv,&f_basic);

	/* Initialize basic computational data */

	f_basic.L[0] = 0.0;	f_basic.L[1] = 0.0; 	f_basic.L[2] = 0.0;
	f_basic.U[0] = 1.0;	f_basic.U[1] = 1.0; 	f_basic.U[2] = 1.0;
	f_basic.gmax[0] = 64;	f_basic.gmax[1] = 64; f_basic.gmax[2] = 64;
	f_basic.boundary[0][0] = f_basic.boundary[0][1] = DIRICHLET_BOUNDARY;
	f_basic.boundary[1][0] = f_basic.boundary[1][1] = DIRICHLET_BOUNDARY;
	f_basic.boundary[2][0] = f_basic.boundary[2][1] = DIRICHLET_BOUNDARY;
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

	if (!RestartRun)
	{
	    /* Initialize interface through level function */
	    rect_box_params.dim = 3;
	    rect_box_params.center[0] = 0.5;
	    rect_box_params.center[1] = 0.5;
	    rect_box_params.center[2] = 0.5;
	    rect_box_params.length[0] = 0.4;
	    rect_box_params.length[1] = 0.3;
	    rect_box_params.length[2] = 0.2;

	    level_func_pack.neg_component = 1;
	    level_func_pack.pos_component = 2;
	    level_func_pack.func_params = (POINTER)&rect_box_params;
	    level_func_pack.func = rect_box_func;
	    level_func_pack.wave_type = FIRST_PHYSICS_WAVE_TYPE;

	    FT_InitIntfc(&front,&level_func_pack);
	}

	/* Initialize velocity field function */

	curv_params.dim = 3;
	curv_params.coeff = 0.1;
	curv_params.epsilon = 0.0001;

	velo_func_pack.func_params = (POINTER)&curv_params;
	velo_func_pack.func = test_curvature_vel;

	FT_InitFrontVeloFunc(&front,&velo_func_pack);

	/* For geometry-dependent velocity, use first 
	* order point propagation function, higher order
	* propagation requires surface propagate, currently
	* in writing, not yet in use. The following override
	* the assigned fourth_order_point_propagate.
	*/

	PointPropagationFunction(&front) = first_order_point_propagate;

	/* Propagate the front */

	propagation_driver(&front);

	clean_up(0);
	return 0;
}

static  void propagation_driver(
        Front *front)
{
        double CFL;

	front->max_time = 3; 
	front->max_step = 2;
	front->print_time_interval = 0.6;
	front->movie_frame_interval = 0.1;

        CFL = Time_step_factor(front) = 0.1;

	printf("CFL = %f\n",CFL);
	printf("Frequency_of_redistribution(front,GENERAL_WAVE) = %d\n",
		Frequency_of_redistribution(front,GENERAL_WAVE));

	if (!RestartRun)
	{
            FT_RedistMesh(front);
	    FT_ResetTime(front);

	    // Always output the initial interface.
	    FT_Save(front);
            FT_Draw(front);

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
		FT_Save(front);
            if (FT_IsDrawTime(front))
                FT_Draw(front);

            if (FT_TimeLimitReached(front))
                    break;

	    FT_TimeControlFilter(front);
	}
        (void) delete_interface(front->interf);
}       /* end propagation_driver */

static int test_curvature_vel(
        POINTER params,
        Front *front,
        POINT *p,
        HYPER_SURF_ELEMENT *hse,
        HYPER_SURF *hs,
        double *vel)
{
	TEST_CURV_PARAMS *curv_params = (TEST_CURV_PARAMS*) params;
        int i;
        double coeff,epsilon,eps;
        double kappa;
        double nor[MAXD];

        coeff = curv_params->coeff;
        epsilon = curv_params->epsilon;

        normal(p,hse,hs,nor,front);
        kappa = fabs(mean_curvature_at_point(p,hse,hs,front));

        for (i = 0; i < curv_params->dim; ++i)
        {
            vel[i] = nor[i]*(coeff - epsilon*kappa);
        }
}


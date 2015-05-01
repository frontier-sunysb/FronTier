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

char *restart_name;
boolean RestartRun;
int RestartStep;

int main(int argc, char **argv)
{
	static Front front;
	static F_BASIC_DATA f_basic;
	static LEVEL_FUNC_PACK level_func_pack;
	NORV_PARAMS curv_params; /* velocity function parameters */
	SURFACE *surf;

	f_basic.dim = 3;	
	FT_Init(argc,argv,&f_basic);

	/* Initialize basic computational data */

	f_basic.L[0] = 0.0;	f_basic.L[1] = 0.0; 	f_basic.L[2] = 0.0;
	f_basic.U[0] = 1.0;	f_basic.U[1] = 0.5; 	f_basic.U[2] = 0.5;
	f_basic.gmax[0] = 100;	f_basic.gmax[1] = 50; f_basic.gmax[2] = 50;
	f_basic.boundary[0][0] = f_basic.boundary[0][1] = DIRICHLET_BOUNDARY;
	f_basic.boundary[1][0] = f_basic.boundary[1][1] = DIRICHLET_BOUNDARY;
	f_basic.boundary[2][0] = f_basic.boundary[2][1] = DIRICHLET_BOUNDARY;
	f_basic.size_of_intfc_state = 0;

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
	    double x0,x1; // x-coordinates of the dumbbell centers
	    double y,z;   // y and z coordinate of the dumbbell centers
	    double R;	  // radius of the connecting cylunder
	    double r;	  // radius of the two end spheres

	    /* Initialize domain */
	    level_func_pack.pos_component = 2;
	    FT_InitIntfc(&front,&level_func_pack);

	    /* Initialize interior surface */
	    x0 = 0.25;		x1 = 0.75;
	    y = 0.25;
	    z = 0.25;
	    R = 0.15;
	    r = 0.075;
	    FT_MakeDumbBellSurf(&front,
			x0,x1,y,z,R,r,
			1,2,		// negative and positive components
			FIRST_PHYSICS_WAVE_TYPE, // wave type
			&surf);
	}

	/* Initialize velocity field function */

	curv_params.dim = 3;
	curv_params.coeff = -0.1;
	curv_params.epsilon = -0.0001;

	front.vfunc = NULL;
        FT_InitSurfVeloFunc(surf,
			"curvature_dependent_velocity",
                        (POINTER)&curv_params,
                        curvature_vel);

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
}

static  void propagation_driver(
        Front *front)
{
        double CFL;

	front->max_time = 1.0; 
	front->max_step = 1000;
	front->print_time_interval = 0.1;
	front->movie_frame_interval = 0.02;

        CFL = Time_step_factor(front) = 0.2;

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

	    FT_PrintTimeStamp(front);
            if (FT_IsSaveTime(front))
		FT_Save(front);
            if (FT_IsDrawTime(front))
                FT_Draw(front);

            if (FT_TimeLimitReached(front))
                    break;

	    FT_TimeControlFilter(front);
	}
}       /* end propagation_driver */

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
	TRANS_PARAMS trans_params; /* velocity function parameters */
	double center[MAXD],radius[MAXD];
	SURFACE *surf;

	f_basic.dim = 3;	
	FT_Init(argc,argv,&f_basic);

	/* Initialize basic computational data */

	f_basic.L[0] = 0.0;	f_basic.L[1] = 0.0; 	f_basic.L[2] = 0.0;
	f_basic.U[0] = 1.0;	f_basic.U[1] = 1.0; 	f_basic.U[2] = 1.0;
	f_basic.gmax[0] = 64;	f_basic.gmax[1] = 64; f_basic.gmax[2] = 64;
	f_basic.boundary[0][0] = f_basic.boundary[0][1] = PERIODIC_BOUNDARY;
	f_basic.boundary[1][0] = f_basic.boundary[1][1] = PERIODIC_BOUNDARY;
	f_basic.boundary[2][0] = f_basic.boundary[2][1] = PERIODIC_BOUNDARY;
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

	/* Initialize domain */
	level_func_pack.pos_component = 2;
	FT_InitIntfc(&front,&level_func_pack);

	/* Initialize interface through level function */
	center[0] = 0.5;
	center[1] = 0.5;
	center[2] = 0.5;
	radius[0] = 0.3;
	radius[1] = 0.2;
	radius[2] = 0.1;
	FT_MakeEllipticSurf(&front,center,radius,
                        1,2,    // negative and positive components
                        FIRST_PHYSICS_WAVE_TYPE,
                        1,      // refinement level
                        &surf);

	/* Initialize velocity field function */

	trans_params.dim = 3;
	trans_params.vel[0] = 0.1;
	trans_params.vel[1] = 0.22;
	trans_params.vel[2] = 0.153;

	front.vfunc = NULL;
        FT_InitSurfVeloFunc(surf,
			"translation_velocity",
                        (POINTER)&trans_params,
                        translation_vel);

	/* Propagate the front */

	PointPropagationFunction(&front) = fourth_order_point_propagate;
	propagation_driver(&front);

	clean_up(0);
}	/* end main */

static  void propagation_driver(
        Front *front)
{
	double CFL;
	front->max_time = 5; 
	front->max_step = 1000;
	front->print_time_interval = 0.5;
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

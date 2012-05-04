/************************************************************************************
FronTier is a set of libraries that implements differnt types of Front Traking algorithms.
Front Tracking is a numerical method for the solution of partial differential equations 
whose solutions have discontinuities.  


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

******************************************************************************/


/*
*				test_map_3d.c:
*
*		User initialization example for Front Package:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*	
*	This example shows how to use an array of points to make a curve.
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

/********************************************************************
 *      Velocity function parameters for the front                  *
 ********************************************************************/

typedef struct {
        double i1,i2;
        double cen1[2],cen2[2];
} DOUBLE_VORTEX_PARAMS;

int main(int argc, char **argv)
{
	static Front front;
	static RECT_GRID comp_grid;
	static F_BASIC_DATA f_basic;
	static LEVEL_FUNC_PACK level_func_pack;
	static VELO_FUNC_PACK velo_func_pack;
        DOUBLE_VORTEX_PARAMS dv_params; /* velocity function parameters */

	f_basic.dim = 2;
	FT_Init(argc,argv,&f_basic);

	/* Initialize basic computational data */

	f_basic.L[0] = 0.0;	f_basic.L[1] = 0.0;
	f_basic.U[0] = 1.0;	f_basic.U[1] = 1.0;
	f_basic.gmax[0] = 200;	f_basic.gmax[1] = 200;
	f_basic.boundary[0][0] = f_basic.boundary[0][1] = DIRICHLET_BOUNDARY;
	f_basic.boundary[1][0] = f_basic.boundary[1][1] = DIRICHLET_BOUNDARY;
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

	    level_func_pack.neg_component = 1;
	    level_func_pack.pos_component = 2;
	    level_func_pack.func_params = NULL;
	    level_func_pack.func = NULL;
	    level_func_pack.num_points = 500;
	    level_func_pack.is_closed_curve = YES;
	    level_func_pack.wave_type = FIRST_PHYSICS_WAVE_TYPE;

	    FT_MatrixMemoryAlloc((POINTER*)&level_func_pack.point_array,500,
	    				2,sizeof(double));
	    int i;
	    for (i = 0; i < 500; ++i)
	    {
	        double phi = i*2.0*PI/500.0;
	        level_func_pack.point_array[i][0] = 0.5 + 0.3*cos(phi);
	        level_func_pack.point_array[i][1] = 0.5 + 0.3*sin(phi);
	    }
	    FT_InitIntfc(&front,&level_func_pack);
	    if (f_basic.dim < 3)
                FT_ClipIntfcToSubdomain(&front);
	}

	/* Initialize velocity field function */

        dv_params.cen1[0] = 0.25;
        dv_params.cen1[1] = 0.50;
        dv_params.cen2[0] = 0.75;
        dv_params.cen2[1] = 0.50;
        dv_params.i1 = -0.5;
        dv_params.i2 =  0.5;

        velo_func_pack.func_params = (POINTER)&dv_params;
        velo_func_pack.func = test_vortex_vel;
        velo_func_pack.point_propagate = fourth_order_point_propagate;

        FT_InitVeloFunc(&front,&velo_func_pack);

        /* Propagate the front */

        test_propagate(&front);

	clean_up(0);
}

static  void test_propagate(
        Front *front)
{
        double CFL;

	front->max_time = 2.0;
	front->max_step = 400;
	front->print_time_interval = 1.0;
	front->movie_frame_interval = 0.02;

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
 *	Sample (circle) velocity function for the front    *
 ********************************************************************/

static int test_vortex_vel(
	POINTER params,
	Front *front,
	POINT *p,
	HYPER_SURF_ELEMENT *hse,
	HYPER_SURF *hs,
	double *vel)
{
	DOUBLE_VORTEX_PARAMS *dv_params = (DOUBLE_VORTEX_PARAMS*)params;
	double *coords = Coords(p);
	double d1,d2;
	double s1,s2;
	double *cen1 = dv_params->cen1;
	double *cen2 = dv_params->cen2;
	double dx1,dy1;
	double dx2,dy2;

	dx1 = coords[0] - cen1[0]; 
	dy1 = coords[1] - cen1[1];
	dx2 = coords[0] - cen2[0]; 
	dy2 = coords[1] - cen2[1];

	d1 = sqrt(sqr(dx1) + sqr(dy1));
	d2 = sqrt(sqr(dx2) + sqr(dy2));

	s1 = dv_params->i1/2.0/PI/d1;
	s2 = dv_params->i2/2.0/PI/d2;

	vel[0] =  s1*dy1/d1 + s2*dy2/d2;
	vel[1] = -s1*dx1/d1 - s2*dx2/d2;
}	/* end test_vortex_vel */

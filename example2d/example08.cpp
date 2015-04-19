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

#include <vector>
#include <FronTier.h>

	/*  Function Declarations */
static void test_propagate(Front*);
static double cap_h_func(POINTER,double*);
static int test_curvature_vel(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*,double*);


char *in_name,*restart_state_name,*restart_name,*out_name;
boolean RestartRun;
int RestartStep;

/********************************************************************
 *	Level function parameters for the initial interface 	    *
 ********************************************************************/

typedef struct {
        double x0;
        double y0;
        double h;
        double w;
} H_PARAMS;


typedef struct {
	int dim;
	double coeff;
	double epsilon;
} TEST_CURV_PARAMS;

/********************************************************************
 *	Velocity function parameters for the front	 	    *
 ********************************************************************/


int main(int argc, char **argv)
{
	static Front front;
	static RECT_GRID comp_grid;
	static F_BASIC_DATA f_basic;
	H_PARAMS h_params;
	static LEVEL_FUNC_PACK level_func_pack;
	static VELO_FUNC_PACK velo_func_pack;
	TEST_CURV_PARAMS curv_params; /* velocity function parameters */
	Locstate  sl;

	FT_Init(argc,argv,&f_basic);
	f_basic.dim = 2;

	/* Initialize basic computational data */

	f_basic.L[0] = 0.0;	f_basic.L[1] = 0.0;
	f_basic.U[0] = 1.0;	f_basic.U[1] = 1.0;
	f_basic.gmax[0] = 100;	f_basic.gmax[1] = 100;
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

	if (!RestartRun)
	{
	    /* Initialize interface through level function */
	    h_params.x0 = 0.5;
            h_params.y0 = 0.5;
            h_params.h = 0.4;
            h_params.w = 0.1;

	    level_func_pack.neg_component = 1;
	    level_func_pack.pos_component = 2;
	    level_func_pack.func_params = (POINTER)&h_params;
	    level_func_pack.func = cap_h_func;
	    level_func_pack.wave_type = FIRST_PHYSICS_WAVE_TYPE;

	    FT_InitIntfc(&front,&level_func_pack);
	    if (f_basic.dim < 3)
                FT_ClipIntfcToSubdomain(&front);
	}

	/* Initialize velocity field function */

	curv_params.dim = 2;
	curv_params.coeff = 0;
	curv_params.epsilon = 0.01;

	velo_func_pack.func_params = (POINTER)&curv_params;
	velo_func_pack.func = test_curvature_vel;
	velo_func_pack.point_propagate = first_order_point_propagate;

	FT_InitFrontVeloFunc(&front,&velo_func_pack);
	
        /* For geometry-dependent velocity, use first
        * order point propagation function, higher order
        * propagation requires surface propagate, currently
        * in writing, not yet in use. The following override
        * the assigned fourth_order_point_propagate.
        */

        front._point_propagate = first_order_point_propagate;

	/* Propagate the front */

	test_propagate(&front);

	clean_up(0);
	return 0;
}

static  void test_propagate(
        Front *front)
{
        double CFL;

	front->max_time = 1.0; 
	front->max_step = 100000;
	front->print_time_interval = 1.0;
	front->movie_frame_interval = 0.01;

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
}       /* end test_propagate */

/********************************************************************
 *	Sample (H shape) level function for the initial interface    *
 ********************************************************************/

static double cap_h_func(
        POINTER func_params,
        double *coords)
{
        H_PARAMS *h_params = (H_PARAMS*)func_params;
        double x0,y0,h,w;
        double dist[3];
        int i,imin;
        double r1,r2,dmin;

        x0 = h_params->x0;
        y0 = h_params->y0;
        h  = h_params->h;
        w  = h_params->w;

        dmin = HUGE;
        r1 = fabs(coords[0]-x0)/(h/2.0);
        r2 = fabs(coords[1]-y0)/(w/2.0);
        dist[0] = std::max(r1,r2) - 1;

        r1 = fabs(coords[0]-x0+(h-w)/2.0)/(w/2.0);
        r2 = fabs(coords[1]-y0)/(h/2.0);
        dist[1] = std::max(r1,r2) - 1;

        r1 = fabs(coords[0]-x0-(h-w)/2.0)/(w/2.0);
        r2 = fabs(coords[1]-y0)/(h/2.0);
        dist[2] = std::max(r1,r2) - 1;

        for (i = 0; i < 3; ++i)
        {
            if(dist[i]<dmin)
            {
                dmin = dist[i];
            }
        }

        return dmin;
}       /* end cap_h_func */

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

        GetFrontNormal(p,hse,hs,nor,front);
	GetFrontCurvature(p,hse,hs,&kappa,front);

        for (i = 0; i < curv_params->dim; ++i)
        {
            vel[i] = nor[i]*(coeff - epsilon*kappa);
        }
}


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
*				example0.c:
*
*		User initialization example for Front Package:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*	
*	This example shows a circle in a double vortex field. It demonstrates
*	the resolution of the front tracking method.
*
*/

#include <vector>
#include <FronTier.h>

	/*  Function Declarations */
static void test_propagate(Front*);
static int rotation_vel(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,
	                       HYPER_SURF*,double*);
static int trans_vel(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,
	                       HYPER_SURF*,double*);
static void print_vol_frac(Front*);

char *in_name,*restart_state_name,*restart_name,*out_name;
boolean RestartRun;
int RestartStep;

/********************************************************************
 *	Velocity function parameters for the front	 	    *
 ********************************************************************/

typedef struct {
        double omega_0;          /* angular velocity */
        double domega_dr;
        double cen[2];
} ROTATION_VEL_PARAMS;

typedef struct {
	double vel[2];
} TRANSLATION_VEL_PARAMS;

int main(int argc, char **argv)
{
	static Front front;
	static RECT_GRID comp_grid;
	static F_BASIC_DATA f_basic;
	static LEVEL_FUNC_PACK level_func_pack;
	static VELO_FUNC_PACK velo_func_pack;
	RECT_BOX_PARAMS box_params;	/* level function parameters */
	ROTATION_VEL_PARAMS rv_params; /* velocity function parameters */
	TRANSLATION_VEL_PARAMS tv_params; /* velocity function parameters */
	Locstate  sl;

	FT_Init(argc,argv,&f_basic);

	/* Initialize basic computational data */

	f_basic.dim = 2;
	f_basic.L[0] = 0.0;     f_basic.L[1] = 0.0;
        f_basic.U[0] = 1.0;     f_basic.U[1] = 1.0;
        f_basic.gmax[0] = 128;  f_basic.gmax[1] = 128;
        f_basic.boundary[0][0] = f_basic.boundary[0][1] = PERIODIC_BOUNDARY;
        f_basic.boundary[1][0] = f_basic.boundary[1][1] = PERIODIC_BOUNDARY;

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

	    box_params.dim = 2;
	    box_params.center[0] = 0.5;
	    box_params.center[1] = 0.5;
	    box_params.length[0] = 0.7;
	    box_params.length[1] = 0.3;

	    level_func_pack.neg_component = 1;
	    level_func_pack.pos_component = 2;
	    level_func_pack.func_params = (POINTER)&box_params;
	    level_func_pack.func = rect_box_func;
	    level_func_pack.wave_type = FIRST_PHYSICS_WAVE_TYPE;
	    FT_InitIntfc(&front,&level_func_pack);
	    FT_ClipIntfcToSubdomain(&front);
	}

	/* Initialize velocity field function */

	/* Use rotational velocity 
        rv_params.cen[0] = 0.5;
        rv_params.cen[1] = 0.5;
        rv_params.omega_0 = -2.0*PI;
        rv_params.domega_dr = 0.0;

	velo_func_pack.func_params = (POINTER)&rv_params;
	velo_func_pack.func = rotation_vel;
	*/
	tv_params.vel[0] = 0.5;
	tv_params.vel[1] = 0.7;
	velo_func_pack.func_params = (POINTER)&tv_params;
	velo_func_pack.func = trans_vel;
	

	velo_func_pack.point_propagate = fourth_order_point_propagate;
	FT_InitFrontVeloFunc(&front,&velo_func_pack);

	/* Propagate the front */

	test_propagate(&front);

	clean_up(0);
	return 0;
}

static  void test_propagate(
        Front *front)
{
        double CFL;

	front->max_time = 3;
	front->max_step = 10;
	front->print_time_interval = 2.0;
	front->movie_frame_interval = 0.05;

        CFL = Time_step_factor(front);
	Frequency_of_redistribution(front,GENERAL_WAVE) = 1000;

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
	    print_vol_frac(front);

	    //Next time step determined by maximum speed of previous
	    //step, assuming the propagation is hyperbolic and
	    //is not dependent on second order derivatives of
	    //the interface such as curvature, and etc.

            FT_SetTimeStep(front);

            printf("\ntime = %f   step = %5d   next dt = %f\n",
                        front->time,front->step,front->dt);
            fflush(stdout);

            if (FT_IsSaveTime(front))
                FT_Save(front);
            if (FT_IsDrawTime(front))
                FT_Draw(front);

            if (FT_TimeLimitReached(front))
                    break;

	    /* Output section, next dt may be modified */

	    FT_TimeControlFilter(front);
        }
        (void) delete_interface(front->interf);
}       /* end test_propagate */

/********************************************************************
 *	Sample (rotation) velocity function for the front    *
 ********************************************************************/

static int rotation_vel(
	POINTER params,
	Front *front,
	POINT *p,
	HYPER_SURF_ELEMENT *hse,
	HYPER_SURF *hs,
	double *vel)
{
	ROTATION_VEL_PARAMS *rv_params = (ROTATION_VEL_PARAMS*)params;
	double *coords = Coords(p);
	double V,xcomp,ycomp;
	double rad;
	double *cen = rv_params->cen;
	double omega_0 = rv_params->omega_0;
	double domega_dr = rv_params->domega_dr;
	double dx,dy;

	dx = coords[0] - cen[0]; 
	dy = coords[1] - cen[1];

	rad = sqrt(sqr(dx) + sqr(dy));
	if (rad == 0.0)
        {
            vel[0] = vel[1] = 0.0;
            return 1;
        }
	xcomp = fabs(coords[1]-cen[0])/rad;
        ycomp = fabs(coords[0]-cen[1])/rad;
        V = rad*(omega_0 + domega_dr*rad);
        if (coords[0]-cen[0] >= 0.0 && coords[1]-cen[1] >= 0.0) 
        {
            vel[0] = -V*xcomp;
            vel[1] =  V*ycomp;
        }
        else if (coords[0]-cen[0] <= 0.0 && coords[1]-cen[1] >= 0.0) 
        {
            vel[0] = -V*xcomp;
            vel[1] = -V*ycomp;
        }
        else if (coords[0]-cen[0] <= 0.0 && coords[1]-cen[1] <= 0.0) 
        {
            vel[0] =  V*xcomp;
            vel[1] = -V*ycomp;
        }
        else if (coords[0]-cen[0] >= 0.0 && coords[1]-cen[1] <= 0.0) 
        {
            vel[0] =  V*xcomp;
            vel[1] =  V*ycomp;
        }
}	/* end rotation_vel */

static int trans_vel(
	POINTER params,
	Front *front,
	POINT *p,
	HYPER_SURF_ELEMENT *hse,
	HYPER_SURF *hs,
	double *vel)
{
	TRANSLATION_VEL_PARAMS *tv_params = (TRANSLATION_VEL_PARAMS*)params;
	vel[0] = tv_params->vel[0];
	vel[1] = tv_params->vel[1];
}	/* trans_vel */

static void print_vol_frac(Front *front)
{
	static FILE *vfile;
	double vol_frac;

	if (vfile == NULL)
	{
	    vfile = fopen("vol_frac","w");
	    fprintf(vfile,"\"Volume fraction vs. time\"\n");
	}

	FT_MakeGridIntfc(front);
	vol_frac = FT_ComputeTotalVolumeFraction(front,1);
	fprintf(vfile,"%f %f\n",front->time,vol_frac);
	fflush(vfile);
	FT_FreeGridIntfc(front);
}	/* end print_vol_frac */

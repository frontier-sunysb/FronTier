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
*				example1.c:
*
*		User initialization example for Front Package:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	This is an example of three circles moving in a double vortex
*	velocity field.
*
*/

#include <FronTier.h>

	/*  Function Declarations */
static void test_propagate(Front*);
static double crystal_curve(POINTER,double*);
static int crystal_vel(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,
	                       HYPER_SURF*,double*);

char *in_name,*restart_state_name,*restart_name,*out_name;
boolean RestartRun;
int RestartStep;

/********************************************************************
 *	Level function parameters for the initial interface 	    *
 ********************************************************************/

typedef struct {
        int num_cir;
        double **cen;
        double *rad;
} TMC_PARAMS;


/********************************************************************
 *	Velocity function parameters for the front	 	    *
 ********************************************************************/

typedef struct {
	double i1,i2;
        double cen1[2],cen2[2];
} CRYSTAL_PARAMS;

int main(int argc, char **argv)
{
	static Front front;
	static RECT_GRID comp_grid;
	static F_BASIC_DATA f_basic;
	static LEVEL_FUNC_PACK level_func_pack;
	static VELO_FUNC_PACK velo_func_pack;
	CRYSTAL_PARAMS dv_params; /* velocity function parameters */
	TMC_PARAMS mc_params;
	Locstate  sl;

	FT_Init(argc,argv,&f_basic);
	f_basic.dim = 2;

	/* Initialize basic computational data */

	f_basic.L[0] = 0.0;	f_basic.L[1] = 0.0;
	f_basic.U[0] = 1.0;	f_basic.U[1] = 1.0;
	f_basic.gmax[0] = 100;	f_basic.gmax[1] = 100;
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

	    mc_params.num_cir = 3;
            FT_VectorMemoryAlloc((POINTER*)&mc_params.rad,mc_params.num_cir,
	    					FLOAT);
            FT_MatrixMemoryAlloc((POINTER*)&mc_params.cen,mc_params.num_cir,
	    					2,FLOAT);
	    mc_params.cen[0][0] = 0.3;
	    mc_params.cen[0][1] = 0.3;
	    mc_params.cen[1][0] = 0.7;
	    mc_params.cen[1][1] = 0.3;
	    mc_params.cen[2][0] = 0.5;
	    mc_params.cen[2][1] = 0.7;
	    mc_params.rad[0] = 0.1;
	    mc_params.rad[1] = 0.1;
	    mc_params.rad[2] = 0.1;

	    level_func_pack.neg_component = 1;
	    level_func_pack.pos_component = 2;
	    level_func_pack.func_params = (POINTER)&mc_params;
	    level_func_pack.func = crystal_curve;
	    level_func_pack.wave_type = FIRST_PHYSICS_WAVE_TYPE;

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
	velo_func_pack.func = crystal_vel;
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

	front->max_time = 0.4;
	front->max_step = 10000;
	front->print_time_interval = 1.0;
	front->movie_frame_interval = 0.02;

        CFL = Time_step_factor(front);

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
 *	Sample (circle) velocity function for the front    *
 ********************************************************************/

static int crystal_vel(
	POINTER params,
	Front *front,
	POINT *p,
	HYPER_SURF_ELEMENT *hse,
	HYPER_SURF *hs,
	double *vel)
{
	CRYSTAL_PARAMS *dv_params = (CRYSTAL_PARAMS*)params;
	double *coords = Coords(p);
	
	double dx = 1.0/45.0;
	double circ[2];
	double nor[2];
	int i,count;
	count = 0;
	for(i=0;i<20;i++)
	{
	     circ[0] = coords[0]+dx*sin( i*(1.0/20.0)*2.0*3.1415);
	     circ[1] = coords[1]+dx*cos( i*(1.0/20.0)*2.0*3.1415);
	     if(component(circ,front->interf) == 2)
	     count++;

	 }
	GetFrontNormal(p,hse,hs,nor,front);
	for (i = 0; i < 2; ++i)
	{
	    vel[i] = nor[i]*(((1.0*count)/20.0));
	    
	    // this slows the growth more
	    if(count <8)
	    vel[i]=0;
	    //vel[i] = nor[i]*exp(20.0/(1.0*count));
	}


}	/* end double_vortex_vel */

LOCAL double crystal_curve(
        POINTER func_params,
        double *coords)
{

	double dist, theta;

	dist =   sqrt(sqr(coords[0]-0.5) + sqr(coords[1]-0.5));
	theta = asin((coords[1]-0.5)/dist);
	return dist - .15 - 0.04*sin(5.0*theta) - .02*sin(10.0*theta);
}       /* end crystal_curve */

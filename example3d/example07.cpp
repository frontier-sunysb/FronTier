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
#include<stdlib.h>

	/*  Function Declarations */
static void propagation_driver(Front*);
static double msphere_func(POINTER,double*);
static int norm_vel_func(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,
	                       HYPER_SURF*,double*);

char *in_name,*restart_state_name,*restart_name,*out_name;
boolean RestartRun;
int RestartStep;
#if defined(__MPI__)
int subdomains[MAXD];
#endif /* defined(__MPI__) */

/********************************************************************
 *	Level function parameters for the initial interface 	    *
 ********************************************************************/
//multiple shape interface
typedef struct {
	double **cent;
	double *r;
        double cir_num;
} MSPHERE_PARAMS;

/********************************************************************
 *	Velocity function parameters for the front	 	    *
 ********************************************************************/

struct _TNORV_PARAMS
{
        int dim;
        double coeff;
};
typedef struct _TNORV_PARAMS TNORV_PARAMS;

int main(int argc, char **argv)
{
	static Front front;
	static RECT_GRID comp_grid;
	static F_BASIC_DATA f_basic;
	MSPHERE_PARAMS d_params;
	static LEVEL_FUNC_PACK level_func_pack;
	static VELO_FUNC_PACK velo_func_pack;
	TNORV_PARAMS norm_params; /* velocity function parameters */
	Locstate  sl;
	int i;

	f_basic.dim = 3;	
	FT_Init(argc,argv,&f_basic);

	/* Initialize basic computational data */

	f_basic.L[0] = 0;	f_basic.L[1] = 0; 	f_basic.L[2] =0; 
	f_basic.U[0] = 4.0;	f_basic.U[1] = 4.0; 	f_basic.U[2] = 4.0;
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
	    d_params.cir_num = 6;

	    FT_MatrixMemoryAlloc((POINTER*)&d_params.cent,(int)d_params.cir_num,(int)f_basic.dim,FLOAT);
	    FT_VectorMemoryAlloc((POINTER*)&d_params.r,(int)d_params.cir_num,FLOAT);
	
	    d_params.cent[0][0] = 1.4; 
	    d_params.cent[0][1] = 2.0; 
	    d_params.cent[0][2] = 2.6; 
	    d_params.r[0] = 0.2;
	    d_params.cent[1][0] = 2.6; 
	    d_params.cent[1][1] = 2.0; 
	    d_params.cent[1][2] = 2.6; 
	    d_params.r[1] = 0.2;
	    d_params.cent[2][0] = 1.4; 
	    d_params.cent[2][1] = 2.0; 
	    d_params.cent[2][2] = 1.4; 
	    d_params.r[2] = 0.2;
	    d_params.cent[3][0] = 2.6; 
	    d_params.cent[3][1] = 2.0; 
	    d_params.cent[3][2] = 1.4; 
	    d_params.r[3] = 0.2;
	    d_params.cent[4][0] = 2.0; 
	    d_params.cent[4][1] = 1.4; 
	    d_params.cent[4][2] = 2.0; 
	    d_params.r[4] = 0.2;
	    d_params.cent[5][0] = 2.0; 
	    d_params.cent[5][1] = 2.6; 
	    d_params.cent[5][2] = 2.0; 
	    d_params.r[5] = 0.2;
	
	    level_func_pack.neg_component = 1;
	    level_func_pack.pos_component = 2;
	    level_func_pack.func_params = (POINTER)&d_params;
	    level_func_pack.func = msphere_func;
	    level_func_pack.wave_type = FIRST_PHYSICS_WAVE_TYPE;

	    FT_InitIntfc(&front,&level_func_pack);
	}
	
	// check initial interface	
	/*
		print_front_output(&front,out_name);
		clean_up(0);
	*/
	
	/* Initialize velocity field function */

	norm_params.dim = 3; 
	norm_params.coeff = 0.5;

	velo_func_pack.func_params = (POINTER)&norm_params;
	velo_func_pack.func = norm_vel_func;

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

	front->max_time = 5.0; 
	front->max_step = 6;
	front->print_time_interval = 0.05;
	front->movie_frame_interval = 0.01;

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

/********************************************************************
 *	Sample (dummbell 3D) level function for the initial interface    *
 ********************************************************************/

static double msphere_func(
        POINTER func_params,
        double *coords)
{
        MSPHERE_PARAMS *d_params = (MSPHERE_PARAMS*)func_params;
	double **cent = d_params->cent; 
	double *r = d_params->r;
	int cir_num = d_params->cir_num;
	int i;
	double dist,dmin = HUGE;

	for(i = 0; i < cir_num; i++){
	   dist = sqrt( sqr(coords[0] - cent[i][0]) + sqr(coords[1] - cent[i][1]) + sqr(coords[2] - cent[i][2]) ) - r[i];
	   
	   if ( dist < dmin){
	 	dmin = dist;
	   }
	}

	return dmin;
}       /* end plate_func */

static int norm_vel_func(
        POINTER params,
        Front *front,
        POINT *p,
        HYPER_SURF_ELEMENT *hse,
        HYPER_SURF *hs,
        double *vel)
{
        TNORV_PARAMS *norv_params;
        int i;
        double coeff;
        double nor[MAXD];
                                                                                
        norv_params = (TNORV_PARAMS*)params;
        coeff = norv_params->coeff;
                                                                                
        GetFrontNormal(p,hse,hs,nor,front);
        for (i = 0; i < norv_params->dim; ++i)
        {
            vel[i] = nor[i]*coeff;
        }
}       /* end normal_vel_func */

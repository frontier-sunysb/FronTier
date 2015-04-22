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
static double plate_func(POINTER,double*);
static int test_rotation_vel(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*,double*);

char *in_name,*restart_state_name,*restart_name,*out_name;
boolean RestartRun;
int RestartStep;

/********************************************************************
 *	Level function parameters for the initial interface 	    *
 ********************************************************************/
//Plate shape interface
typedef struct {
        double x;
        double y;
	double z;
	double w1;
        double w2;
        double h;
} PLATE_PARAMS;


typedef struct
{       
       int dim;
       double x_0,y_0;
       double omega_0;
       double grad;
} TEST_ROTATION_PARAMS;

/********************************************************************
 *	Velocity function parameters for the front	 	    *
 ********************************************************************/


int main(int argc, char **argv)
{
	static Front front;
	static RECT_GRID comp_grid;
	static F_BASIC_DATA f_basic;
	PLATE_PARAMS d_params;
	static LEVEL_FUNC_PACK level_func_pack;
	static VELO_FUNC_PACK velo_func_pack;
	TEST_ROTATION_PARAMS cr_params; /* velocity function parameters */
	Locstate  sl;

	f_basic.dim = 3;	
	FT_Init(argc,argv,&f_basic);

	/* Initialize basic computational data */

	f_basic.L[0] = 0;	f_basic.L[1] = 0; 	f_basic.L[2] =2.0; 
	f_basic.U[0] = 5.0;	f_basic.U[1] = 1.0; 	f_basic.U[2] = 3.0;
	f_basic.gmax[0] = 100;	f_basic.gmax[1] = 20; f_basic.gmax[2] = 20;
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
	    d_params.x = 2.5;
	    d_params.y = 0.5;
            d_params.z = 2.5;
            d_params.w1 = 4.0;
            d_params.w2 = 0.8;
            d_params.h = 0.1;

	    level_func_pack.neg_component = 1;
	    level_func_pack.pos_component = 2;
	    level_func_pack.func_params = (POINTER)&d_params;
	    level_func_pack.func = plate_func;
	    level_func_pack.wave_type = FIRST_PHYSICS_WAVE_TYPE;

	    FT_InitIntfc(&front,&level_func_pack);
	}

	/* Initialize velocity field function */

	cr_params.dim = 3;
	cr_params.x_0 = 0.5;
	cr_params.y_0 = 2.5;
	cr_params.omega_0 = 2.0;
	cr_params.grad = 0.0;

	velo_func_pack.func_params = (POINTER)&cr_params;
	velo_func_pack.func = test_rotation_vel;

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

	front->max_time = 2.0; 
	front->max_step = 534;
	front->print_time_interval = 0.02;
	front->movie_frame_interval = 0.004;

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
 *	Sample (plate 3D) level function for the initial interface    *
 ********************************************************************/

static double plate_func(
        POINTER func_params,
        double *coords)
{
        PLATE_PARAMS *d_params = (PLATE_PARAMS*)func_params;
	double x,y,z,w1,w2,h,r0,r1,rx0,rx1,ry0,ry1;
	double arg;

        x = d_params->x;
        y = d_params->y;
	z = d_params->z;
        w1  = d_params->w1;
	w2 = d_params->w2;
        h  = d_params->h;
	
	if ((coords[2] - z) >=0 ){
	    r0 = (coords[2] - z)*w1 / h;
	    r1 = (coords[2] - z)*w2 / h;
	} else {
	    r0 = -(coords[2] - z)*w1 / h;
	    r1 = -(coords[2] - z)*w2 / h;
	}
	if ((coords[0] - x) >=0 ){
	    rx0 = (coords[0] - x)*w2 / w1;
	    rx1 = (coords[0] - x)*h / w1;
	} else {
	    rx0 = -(coords[0] - x)*2 / w1;
	    rx1 = -(coords[0] - x)*h / w1;
	}
	if ((coords[1] - y) >=0 ){
	    ry0 = (coords[1] - y)*w1 / w2;
	    ry1 = (coords[1] - y)*h / w2;
	} else {
	    ry0 = -(coords[1] - y)*w1 / w2;
	    ry1 = -(coords[1] - y)*h / w2;
	}
	if ((x - r0) <= coords[0] && coords[0] <= (x + r0) )
	{
	    if ((y - r1) <= coords[1] && coords[1] <= (y + r1) ){
	    	if (coords[2] >= z) 
		    arg = sqrt(h*(coords[2] - z)/2.0) - h/2.0;
		else                
		    arg = sqrt(-h*(coords[2] - z)/2.0) - h/2.0;
	    }
        }
	if ((y - rx0) <= coords[1] && coords[1] <= (y + rx0) )
	{
	    if ((z - rx1) <= coords[2] && coords[2] <= (z + rx1) ){
	    	if (coords[0] >= x) 
		    arg = sqrt(w1*(coords[0] - x)/2.0) - w1/2.0;
		else                
		    arg = sqrt(-w1*(coords[0] - x)/2.0) - w1/2.0;
	    }
        }
	if ((x - ry0) <= coords[0] && coords[0] <= (x + ry0) )
	{
	    if ((z - ry1) <= coords[2] && coords[2] <= (z + ry1) ){
	    	if (coords[1] >= y) 
		    arg = sqrt(w2*(coords[1] - y)/2.0) - w2/2.0;
		else                
		    arg = sqrt(-w2*(coords[1] - y)/2.0) - w2/2.0;
	    }
        }
	return arg;
}       /* end plate_func */

static int test_rotation_vel(
	POINTER params,
	Front *front,
	POINT *p,
	HYPER_SURF_ELEMENT *hse,
	HYPER_SURF *hs,
	double *vel)
{
	TEST_ROTATION_PARAMS *c_params;
	double x0,y0;
	double omega_0,domega_dr;
	double rad,V;
	double xcomp,ycomp;
	double *coords = Coords(p);

	c_params = (TEST_ROTATION_PARAMS*)params;

	x0 = c_params->x_0;
	y0 = c_params->y_0;
	omega_0 = c_params->omega_0;
	domega_dr = c_params->grad;

	rad = sqrt(sqr(coords[1] - x0) + sqr(coords[2] - y0));
	if (rad == 0.0)
	{
	    vel[0] = vel[1] = vel[2] = 0.0;
	    return 1;
	}
	xcomp = fabs(coords[2]-y0)/rad;
	ycomp = fabs(coords[1]-x0)/rad;
	V = rad*(omega_0 + domega_dr*rad)*(1+coords[0]*20);

	vel[0] = 0.0;
	if (coords[1]-x0 >= 0.0 && coords[2]-y0 >= 0.0) /*1st quadrant*/ 
	{
	    vel[1] = -V*xcomp;
	    vel[2] =  V*ycomp;
	}
	else if (coords[1]-x0 <= 0.0 && coords[2]-y0 >= 0.0) /*2nd quadrant*/ 
	{
	    vel[1] = -V*xcomp;
	    vel[2] = -V*ycomp;
	}
	else if (coords[1]-x0 <= 0.0 && coords[2]-y0 <= 0.0) /*3rd quadrant*/ 
	{
	    vel[1] =  V*xcomp;
	    vel[2] = -V*ycomp;
	}
	else if (coords[1]-x0 >= 0.0 && coords[2]-y0 <= 0.0) /*4th quadrant*/ 
	{
	    vel[1] =  V*xcomp;
	    vel[2] =  V*ycomp;
	}
	vel[0] = 0;
}	/* end shear_motion_vel */

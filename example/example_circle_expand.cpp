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
static int  test_vel(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,
			HYPER_SURF*,double*);
static void computeError(Front*);

char *in_name,*restart_state_name,*restart_name,*out_name;
boolean RestartRun;
int RestartStep;
boolean binary = YES;
/********************************************************************
 *	Level function parameters for the initial interface 	    *
 ********************************************************************/

typedef struct {
	        /* equation for line is x^2/a^2 + y^2/b^2 = 1 */
        double C; 
} RADIAL_VEL_PARAMS;


int main(int argc, char **argv)
{
	static Front front;
	static RECT_GRID comp_grid;
	static F_BASIC_DATA f_basic;
	static LEVEL_FUNC_PACK level_func_pack;
	static VELO_FUNC_PACK velo_func_pack;
	RADIAL_VEL_PARAMS rv_params; /* velocity function parameters */
	Locstate  sl;
	static int NUM = 200;

	FT_Init(argc,argv,&f_basic);
	f_basic.dim = 2;

	/* Initialize basic computational data */

	f_basic.L[0] = 0.0;	f_basic.L[1] = 0.0;
	f_basic.U[0] = 1.0;	f_basic.U[1] = 1.0;
	f_basic.gmax[0] = NUM;	f_basic.gmax[1] = NUM;		//myex grid size
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
	add_to_debug("high_order_reist");

	if (!RestartRun)
	{
	    /* Initialize interface through level function */

	    level_func_pack.neg_component = 1;
	    level_func_pack.pos_component = 2;
	    level_func_pack.wave_type = FIRST_PHYSICS_WAVE_TYPE;

	    level_func_pack.func_params = NULL;
	    level_func_pack.func = NULL;

	    level_func_pack.num_points = f_basic.gmax[0];	       
					//myex num points
	    level_func_pack.is_closed_curve = YES;

	    FT_MatrixMemoryAlloc((POINTER*)&level_func_pack.point_array,
	    			level_func_pack.num_points,
				2,sizeof(double));
	    int i;
	    for (i = 0; i < level_func_pack.num_points; ++i)
	    {
	    	double phi = i*2.0*PI/(double)level_func_pack.num_points;
	    	level_func_pack.point_array[i][0] = 0.5 + 0.15*cos(phi);
	    	level_func_pack.point_array[i][1] = 0.5 + 0.15*sin(phi);
	    }

	    FT_InitIntfc(&front,&level_func_pack);
	    if (f_basic.dim < 3)
                FT_ClipIntfcToSubdomain(&front);
	}
	Frequency_of_redistribution(&front,GENERAL_WAVE) = 
			level_func_pack.num_points/20;

	/* Initialize velocity field function */

	rv_params.C = 0.2;
	velo_func_pack.func_params = (POINTER)&rv_params;
	velo_func_pack.func = test_vel;
	velo_func_pack.point_propagate = first_order_point_propagate;

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

	front->max_time = 1.0;
	front->max_step = 1000000;
	front->print_time_interval = 1.0;
	front->movie_frame_interval = 0.02;
	vparams->time = 0.5*front->max_time;

        CFL = Time_step_factor(front);

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


static int test_vel(
        POINTER params,
        Front *front,
        POINT *p,
        HYPER_SURF_ELEMENT *hse,
        HYPER_SURF *hs,
        double *vel)
{
        RADIAL_VEL_PARAMS *rv_params = (RADIAL_VEL_PARAMS*)params;
	double C = rv_params->C;
        int i, dim = front->rect_grid->dim;
	double r = 0.0;
	double speed;
	double *coords = Coords(p);

	speed = C*front->rect_grid->h[0];
	for (i = 0; i < dim; ++i)
	    r += sqr(coords[i] - 0.5);
	r = sqrt(r);
	for (i = 0; i < dim; ++i)
	    vel[i] = C*(coords[i] - 0.5)/r;
}       /* end vortex_vel */

static void computeError(Front *front)
{
	INTERFACE *intfc = front->interf;
	CURVE **c;
	RADIAL_VEL_PARAMS *vparams = (RADIAL_VEL_PARAMS*)front->vparams;
	double C = vparams->C;
	double *points,*p;
	int i,num_points,num_curves;
	int dim = intfc->dim;
	double time = front->time;
	double L1_error,L2_error,Li_error;
	double r,R;
	static FILE *efile;

	if (efile == NULL)
	    efile = fopen("errorFile","w");
	R = 0.15 + C*time;
	num_curves = NumOfCurves(intfc);
	FT_VectorMemoryAlloc((POINTER*)&c,num_curves,sizeof(CURVE*));
	ArrayOfCurves(intfc,c);
	for (i = 0; i < num_curves; ++i)
	{
	    if (wave_type(c[i]) < FIRST_PHYSICS_WAVE_TYPE)
		continue;
	    num_points = NumOfCurvePoints(c[i]);
	    printf("circle radius is %24.18g\t", R);
	    printf("number of points: %d\t",num_points);
	    printf("radius/points is: %24.18g\n",R/num_points);
	    FT_VectorMemoryAlloc((POINTER*)&points,num_points*dim,sizeof(double));
	    ArrayOfCurvePoints(c[i],points);
	    break;
	}
	L1_error = 0.0; L2_error = 0.0; Li_error = 0.0;
	for (i = 0; i < num_points; ++i)
	{
	    p = points + i*dim;
	    r = sqr(p[0] - 0.5) + sqr(p[1] - 0.5);
	    r = sqrt(r);
	    L1_error += fabs(r-R);
	    L2_error += sqr(fabs(r-R));
	    if (Li_error<fabs(r-R))
	      Li_error = fabs(r-R);
	}
	L1_error /= (double)num_points;
	L2_error = sqrt(L2_error/(double)num_points);
	/*printf("number of curves: %d \n",num_curves);*/
	printf("time:%24.18g  L1:%24.18g  L2: %24.18g  Li:%24.18g\n",time,L1_error, L2_error,Li_error);
	fprintf(efile,"%24.18g   %24.18g  %24.18g  %24.18g\n",time,L1_error, L2_error,Li_error);
	fflush(efile);
	FT_FreeThese(2,c,points);
}	/* end computeError */

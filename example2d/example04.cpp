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
*				example2.c:
*
*		User initialization example for Front Package:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	This is example of three circles all moving the a normal velocity.
*	Bifurcation occurs when they meet each other. FronTier solves
*	the bifurcation automatically.
*
*/

#include <FronTier.h>

	/*  Function Declarations */
static void test_propagate(Front*);
static int norm_vel_func(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,
	                       HYPER_SURF*,double*);
static void map_output_interface(Front*,char*);

char *in_name,*restart_state_name,*restart_name,*out_name;
boolean RestartRun;
int RestartStep;

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
	static LEVEL_FUNC_PACK level_func_pack;
	static VELO_FUNC_PACK velo_func_pack;
	TNORV_PARAMS norm_params; /* velocity function parameters */
	MC_PARAMS mc_params;
	Locstate  sl;

	FT_Init(argc,argv,&f_basic);
	f_basic.dim = 2;

	/* Initialize basic computational data */

	f_basic.L[0] = 0.0;	f_basic.L[1] = 0.0;
	f_basic.U[0] = 1.0;	f_basic.U[1] = 1.0;
	f_basic.gmax[0] = 40;	f_basic.gmax[1] = 40;
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

	    mc_params.dim = 2;
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
	    level_func_pack.func = multi_circle_func;
	    level_func_pack.wave_type = FIRST_PHYSICS_WAVE_TYPE;

	    FT_InitIntfc(&front,&level_func_pack);
	    if (f_basic.dim < 3)
                FT_ClipIntfcToSubdomain(&front);
	}

	/* Initialize velocity field function */

	norm_params.dim = 2;
	norm_params.coeff = 0.1;

	velo_func_pack.func_params = (POINTER)&norm_params;
	velo_func_pack.func = norm_vel_func;
	velo_func_pack.point_propagate = first_order_point_propagate;

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

	front->max_time = 3.0;
	front->max_step = 10000;
	front->print_time_interval = 2.0;
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
        double curvature;
        double nor[MAXD];
                                                                                
        norv_params = (TNORV_PARAMS*)params;
        coeff = norv_params->coeff;
                                                                                
        GetFrontNormal(p,hse,hs,nor,front);
        for (i = 0; i < norv_params->dim; ++i)
        {
            vel[i] = nor[i]*coeff;
        }
}       /* end normal_vel_func */

static void map_output_interface(
	Front *front,
	char *out_name)
{
	INTERFACE *intfc = front->interf;
        FILE *out_file;
        char filename[100];
	int i,j,k,step;
	step = front->step;
	sprintf(filename,"%s.ts%s",out_name,right_flush(step,5));
	out_file = fopen(filename,"w");
	
	int num_curves = I_NumOfIntfcCurves(intfc);
	int dim = Dimension(intfc);
	int num_nodes = I_NumOfIntfcNodes(intfc);
	int num_bonds = I_NumOfIntfcBonds(intfc);
	int num_points = I_NumOfIntfcPoints(intfc);

	CURVE **curves;
	curves = (CURVE**)malloc(num_curves*sizeof(CURVE*));

	I_ArrayOfIntfcCurves(intfc,curves);

	fprintf(out_file,"Interface at step %d\n\n",step);
	fprintf(out_file,"Dimension = %d\n",dim);
	fprintf(out_file,"Number of Curves = %d\n",num_curves);
	fprintf(out_file,"Number of Nodes = %d\n",num_nodes);
	fprintf(out_file,"Number of Interface Bonds = %d\n",num_bonds);
	fprintf(out_file,"Number of Interface Points = %d\n",num_points);
	for (i = 0; i < num_curves; ++i)
	{
	    double *coords;

	    num_bonds = I_NumOfCurveBonds(curves[i]);
	    num_points = I_NumOfCurvePoints(curves[i]);
	    coords = (double*)malloc(num_points*dim*sizeof(double));

	    ArrayOfCurvePoints(curves[i],coords);

	    fprintf(out_file,"Number of Bonds on Curve %d = %d\n",
	    			i+1,num_bonds);
	    fprintf(out_file,"Number of Points on Curve %d = %d\n\n",
	    			i+1,num_points);
	    for (j = 0; j < num_points; ++j)
	    {
	    	for (k = 0; k < dim; ++k)
		    fprintf(out_file,"%f  ",coords[dim*j+k]);
	    	fprintf(out_file,"\n");
	    }
	    fprintf(out_file,"\n\n");
	}
	fclose(out_file);
}	/* end map_output */

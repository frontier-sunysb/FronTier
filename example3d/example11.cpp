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
*				example09.c:
*
*		User initialization example for FronTier++ Package:
*		This example checks component (index) of surfaces.
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*/

#include <FronTier.h>

/*  Function Declarations */
typedef struct{
	long num_bubbles; 	     /*number of bubbles in domain*/
	double **center;	     /*center of the bubbles*/
	double *radius;		     /*radius of the bubbles*/
} BUBBLE_PARAMS;

static void initBubbles(Front*,const BUBBLE_PARAMS*);
static void readBubbleParams(Front*, BUBBLE_PARAMS*);
static void setRestartBubbleVelocity(Front*);
static void propagation_driver(Front*);

char *restart_name,*in_name,*out_name;
boolean RestartRun;
int RestartStep;

int main(int argc, char **argv)
{
	static Front front;
	static F_BASIC_DATA f_basic;
	static LEVEL_FUNC_PACK level_func_pack;
	BUBBLE_PARAMS params;

	FT_Init(argc,argv,&f_basic);

	in_name 		= f_basic.in_name;
	out_name 		= f_basic.out_name;
	restart_name            = f_basic.restart_name;
        RestartRun              = f_basic.RestartRun;
        RestartStep             = f_basic.RestartStep;

        sprintf(restart_name,"%s/intfc-ts%s",restart_name,
                        right_flush(RestartStep,7));
        if (pp_numnodes() > 1)
        {
            sprintf(restart_name,"%s-nd%s",restart_name,
                        right_flush(pp_mynode(),4));
        }

	FT_ReadSpaceDomain(in_name,&f_basic);
	FT_StartUp(&front,&f_basic);
	FT_InitDebug(in_name);

	level_func_pack.pos_component = 2; /* default component*/

	if (!RestartRun)
	{
	    FT_InitIntfc(&front,&level_func_pack);
	    readBubbleParams(&front,&params);
	    initBubbles(&front,&params);
	}
	else
	    setRestartBubbleVelocity(&front);
	PointPropagationFunction(&front) = first_order_point_propagate;
	FT_Draw(&front);


	/* Propagate the front */
        propagation_driver(&front);
	clean_up(0);
}

static  void propagation_driver(
        Front *front)
{
        double CFL;

	FT_ReadTimeControl(in_name,front);
        CFL = Time_step_factor(front);

        printf("CFL = %f\n",CFL);
        printf("Frequency_of_redistribution(front,GENERAL_WAVE) = %d\n",
                Frequency_of_redistribution(front,GENERAL_WAVE));

        if (!RestartRun)
        {
            FT_RedistMesh(front);
            FT_ResetTime(front);

            /* Always output the initial interface.*/
            FT_Save(front);
            FT_Draw(front);

            /* This is a virtual propagation to get maximum front 
               speed to determine the first time step.*/

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

            /*Next time step determined by maximum speed of previous
              step, assuming the propagation is hyperbolic and
	      is not dependent on second order derivatives of
              the interface such as curvature, and etc.*/

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

}

static void initBubbles(
	Front *front,
	const BUBBLE_PARAMS* params)
{
	int i,j,dim;
	int w_type;
	int neg_comp,pos_comp;
	static NORV_PARAMS curv_params;
	SURFACE *surf;
	dim = FT_Dimension();

	curv_params.dim = 3;
        curv_params.coeff = 0.04;
        curv_params.epsilon = 0.001;
        front->vfunc = NULL;

	neg_comp = 1;
	pos_comp = 2;
	w_type = FIRST_PHYSICS_WAVE_TYPE;

	int count = 0, num_iter = 0;
	for (i = 0; i < params->num_bubbles; ++i)
	{
	    FT_MakeSphericalSurf(front,params->center[i],params->radius[i],
			neg_comp,pos_comp,w_type,2,&surf);
	    if (FT_CheckSurfCompConsistency(front,surf) == NO)
	    {
		(void) printf("Inconsistency, bubbles are intersecting!\n");
		clean_up(ERROR);
	    }
	    FT_InitSurfVeloFunc(surf,"curvature_dependent_velocity",
				(POINTER)&curv_params,curvature_vel);
	}
	FT_SetGlobalIndex(front);
}	/* end initBubbles */

static void setRestartBubbleVelocity(
	Front *front)
{
	static NORV_PARAMS curv_params;
	SURFACE **s;

	curv_params.dim = 3;
        curv_params.coeff = 0.04;
        curv_params.epsilon = 0.001;
        front->vfunc = NULL;
	intfc_surface_loop(front->interf,s)
	{
	    if (strcmp((*s)->vfunc_name,"curvature_dependent_velocity") == 0)
	    {
	    	FT_InitSurfVeloFunc((*s),"curvature_dependent_velocity",
				(POINTER)&curv_params,curvature_vel);
	    }
	}
}	/* end setRestartBubbleVelocity */

static void readBubbleParams(
	Front* front, 
	BUBBLE_PARAMS* params)
{
	char msg[200]; 
        char *inname = InName(front);
        FILE *infile = fopen(inname,"r");
	int  i,num_bubbles;
	double **center;
	double *radius;

	sprintf(msg,"Enter number of bubbles: ");
	CursorAfterString(infile,msg);
        fscanf(infile,"%d",&num_bubbles);
        (void) printf("%d\n",num_bubbles);
	params->num_bubbles = num_bubbles;
	FT_VectorMemoryAlloc((POINTER*)&radius,num_bubbles,sizeof(double));
	FT_MatrixMemoryAlloc((POINTER*)&center,num_bubbles,3,sizeof(double));
	params->center = center;
	params->radius = radius;
	for (i = 0; i < num_bubbles; ++i)
	{
	    sprintf(msg,"Enter center coordinate of bubble %d: ",i+1);
            CursorAfterString(infile,msg);
            fscanf(infile,"%lf %lf %lf",center[i],center[i]+1,center[i]+2);
            (void) printf("%f %f %f\n",center[i][0],center[i][1],center[i][2]);
	    sprintf(msg,"Enter radius of bubble %d: ",i+1);
            CursorAfterString(infile,msg);
            fscanf(infile,"%lf",radius+i);
            (void) printf("%f\n",radius[i]);
	}
	fclose(infile);
}	/* end readBubbleParams */

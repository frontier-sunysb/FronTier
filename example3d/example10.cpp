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
	unsigned short int seeds[3]; /*seed for rand number generation*/
	long num_bubbles; 	     /*number of bubbles in domain*/
	double mean_rad;	     /*mean radius of the bubbles*/
	double dev_rad;		     /*deviation of the radius*/
} RAND_BUBBLE_PARAMS;

static void initRandomBubbles(Front*,const RAND_BUBBLE_PARAMS*);
static void readRandBubbleParams(Front*, RAND_BUBBLE_PARAMS*);
static void propagation_driver(Front*);

char *restart_name,*in_name,*out_name;
boolean RestartRun;
int RestartStep;

int main(int argc, char **argv)
{
	static Front front;
	static F_BASIC_DATA f_basic;
	static LEVEL_FUNC_PACK level_func_pack;
	RAND_BUBBLE_PARAMS params;

	FT_Init(argc,argv,&f_basic);

	in_name 		= f_basic.in_name;
	out_name 		= f_basic.out_name;
	restart_name            = f_basic.restart_name;
        RestartRun              = f_basic.RestartRun;
        RestartStep             = f_basic.RestartStep;

	FT_ReadSpaceDomain(in_name,&f_basic);
	FT_StartUp(&front,&f_basic);
	FT_InitDebug(in_name);

	level_func_pack.pos_component = 2; /* default component*/
	FT_InitIntfc(&front,&level_func_pack);

	readRandBubbleParams(&front,&params);
	initRandomBubbles(&front,&params);
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

static void initRandomBubbles(
	Front *front,
	const RAND_BUBBLE_PARAMS* params)
{
	double R,x,center[MAXD],radius[MAXD];
	int neg_comp,pos_comp,w_type,dim,j;
	double *L = front->rect_grid->L;
	double *U = front->rect_grid->U;
	double *h = front->rect_grid->h;
	double tol = 5*h[0];
	static NORV_PARAMS curv_params;
	GAUSS_PARAMS gauss_params;
        UNIFORM_PARAMS uniform_params;
	SURFACE *surf;
	dim = FT_Dimension();

	curv_params.dim = 3;
        curv_params.coeff = 0.01;
        curv_params.epsilon = -0.001;
        front->vfunc = NULL;

	gauss_params.mu = params->mean_rad;
        gauss_params.sigma = params->dev_rad;
        uniform_params.a = 0.0;
        uniform_params.b = 1.0;

	neg_comp = 1;
	pos_comp = 2;
	w_type = FIRST_PHYSICS_WAVE_TYPE;

	int count = 0, num_iter = 0;
	while (count < params->num_bubbles)
	{
		num_iter ++;
		if (num_iter > 2*params->num_bubbles)
		{
			(void) printf("Too many times, "
				      "try to reduce numbers or radius\n");
			clean_up(ERROR);
		}
		R = gauss_center_limit((POINTER)&gauss_params,
			   const_cast<RAND_BUBBLE_PARAMS*>(params)->seeds);	
		for (j = 0; j < dim; ++j)
            	{
                    x = dist_uniform((POINTER)&uniform_params,
			   const_cast<RAND_BUBBLE_PARAMS*>(params)->seeds);
                    center[j] = L[j] + R + tol + x*(U[j] - L[j] - 2*(R+tol) );
		    radius[j] = R;
                }
		
		FT_MakeEllipticSurf(front,center,radius,neg_comp,pos_comp,
			w_type,2,&surf);
		FT_InitSurfVeloFunc(surf,(POINTER)&curv_params,curvature_vel);
		if (FT_CheckSurfCompConsistency(front,surf) == NO)
		{
			/*delete surface*/
			delete_surface(surf);
			continue;
		}
		count ++;
	}
	FT_SetGlobalIndex(front);
	PointPropagationFunction(front) = first_order_point_propagate;
}	/* end initInteriorSurfaces */

static void readRandBubbleParams(Front* front, RAND_BUBBLE_PARAMS* params)
{
	char msg[200]; 
        char *inname = InName(front);
        FILE *infile = fopen(inname,"r");
	int  num_bubbles;
	unsigned short int seeds[3];
	double r_bar, sigma;

	sprintf(msg,"Enter number of bubbles: ");
	CursorAfterString(infile,msg);
        fscanf(infile,"%d",&num_bubbles);
        (void) printf("%d\n",num_bubbles);
	params->num_bubbles = num_bubbles;

	sprintf(msg,"Enter mean radius of bubbles: ");
        CursorAfterString(infile,msg);
        fscanf(infile,"%lf",&r_bar);
        (void) printf("%f\n",r_bar);
	params->mean_rad = r_bar;

        sprintf(msg,"Enter standard deviation of radius: ");
        CursorAfterString(infile,msg);
        fscanf(infile,"%lf",&sigma);
        (void) printf("%f\n",sigma);	
	params->dev_rad = sigma;

	sprintf(msg,"Enter seeds for random number generation: ");
	CursorAfterString(infile,msg);
	for (int j = 0; j < 3; ++j)
        {
            fscanf(infile,"%hu ",&seeds[j]);
            (void) printf("%hu ",seeds[j]);
	    params->seeds[j] = seeds[j];
        }
	(void) printf("\n");
	fclose(infile);
}

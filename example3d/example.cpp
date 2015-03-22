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
*				example.c:
*
*		User initialization example for Front Package:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*/

#include <FronTier.h>

	/*  Function Declarations */
static void initInteriorSurfaces(Front*);
static void initVelocityField(Front*);
static void test_propagate(Front*);

char *restart_name,*in_name,*out_name;
boolean RestartRun;

int main(int argc, char **argv)
{
	static Front front;
	static F_BASIC_DATA f_basic;
	static LEVEL_FUNC_PACK level_func_pack;

	FT_Init(argc,argv,&f_basic);

	/* Initialize basic computational data */

	f_basic.size_of_intfc_state = 0;
        in_name                 = f_basic.in_name;
        out_name                = f_basic.out_name;
        RestartRun              = f_basic.RestartRun;

	FT_ReadSpaceDomain(in_name,&f_basic);
	FT_StartUp(&front,&f_basic);

        if (!RestartRun)
	{
	    level_func_pack.pos_component = 2;
	    FT_InitIntfc(&front,&level_func_pack);
	}
	initInteriorSurfaces(&front);
	FT_Draw(&front,out_name,YES);
	clean_up(0);

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
}

static  void test_propagate(
        Front *front)
{
        double CFL;
        CFL = Time_step_factor(front);

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
	FT_PrintTimeStamp(front);

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

            if (FT_IsSaveTime(front))
		FT_Save(front);
            if (FT_IsDrawTime(front))
                FT_Draw(front);

            if (FT_TimeLimitReached(front))
	    {
		FT_PrintTimeStamp(front);
                break;
	    }

	    FT_TimeControlFilter(front);
	    FT_PrintTimeStamp(front);
	}
	FT_FreeMainIntfc(front);
}       /* end test_propagate */

static void initInteriorSurfaces(
	Front *front)
{
	char *inname = InName(front);
	char string[100];
	FILE *infile = fopen(inname,"r");
	double center[MAXD],radius[MAXD];
	int neg_comp,pos_comp,w_type;
	SURFACE *surf;

	neg_comp = 1;
	pos_comp = 2;
	w_type = FIRST_PHYSICS_WAVE_TYPE;
	CursorAfterString(infile,"Enter interior surface type: ");
	fscanf(infile,"%s",string);
	printf("%s\n",string);
	switch (string[0])
	{
	case 's':
	case 'S':
	    center[0] = center[1] = center[2] = 0.5;
	    radius[0] = radius[1] = radius[2] = 0.3;
	    FT_MakeEllipticSurf(front,center,radius,neg_comp,pos_comp,
			w_type,2,&surf);
	    break;
	case 'e':
	case 'E':
	    center[0] = center[1] = center[2] = 0.5;
	    radius[0] = 0.2;
	    radius[1] = 0.4;
	    radius[2] = 0.3;
	    FT_MakeEllipticSurf(front,center,radius,neg_comp,pos_comp,
			w_type,2,&surf);
	    break;
	default:
	    printf("Unknown surface type!\n");
	    clean_up(ERROR);
	}
}	/* end initInteriorSurfaces */

static void initVelocityField(
	Front *front)
{
}	/* end initVelocityField */


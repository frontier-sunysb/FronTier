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
static void initInteriorSurfaces(Front*);
static void initVelocityField(Front*);
static void propagation_driver(Front*);

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
	FT_Draw(&front);
	clean_up(0);

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
}

static  void propagation_driver(
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
}       /* end propagation_driver */

/*
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
}*/	/* end initInteriorSurfaces */
static void initInteriorSurfaces(
	Front *front)
{
	char *inname = InName(front);
	char string[100];
	FILE *infile = fopen(inname,"r");
	double center[MAXD],radius[MAXD],edge[MAXD];
	double nor[MAXD],pt[MAXD],center1[MAXD];
	double x0,x1,y,z,R,RR,r,slope,edgel,height,height1;
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
	    if(strcmp(string,"sphere") == 0 || 
	    strcmp(string,"SPHERE") == 0)
	    {
		center[0] = center[1] = center[2] = 0.5;
	    	radius[0] = radius[1] = radius[2] = 0.3;
	    	FT_MakeEllipticSurf(front,center,radius,neg_comp,
		pos_comp,w_type,2,&surf);
	    }
	    if(strcmp(string,"stellatedoctahedron") == 0 
	       || strcmp(string,"STELLATEDOCTAHEDRON") == 0)
            {
                center[0] = center[1] = center[2] = 1;
                edgel = 2.828;
                FT_MakeStellatedOctahedronSurf(front,center,
		edgel,neg_comp,pos_comp,w_type,&surf);
            }
	    break;
	case 'e':
	case 'E':
	    center[0] = center[1] = center[2] = 0.5;
	    radius[0] = 0.5;
	    radius[1] = 0.3;
 	    radius[2] = 0.3;
	    FT_MakeEllipticSurf(front,center,radius,neg_comp,
	    pos_comp,w_type,2,&surf);
	    break; 
	case 'c':
	case 'C':
  	    if(strcmp(string,"cylinder") == 0 
	       || strcmp(string,"CYLINDER") == 0)
	    {
	    	center[0] = center[1] = center[2] = 0; 
	    	radius[0] = 0.3; 
	    	height = 1;  
	    	FT_MakeCylinderSurf(front,center,radius[0],
		height,neg_comp,pos_comp,
                        w_type,&surf);
	    }
	    if(strcmp(string,"cuboid") == 0 
	       || strcmp(string,"CUBOID") == 0)
            {
		center[0] = center[1] = center[2] = 0.5;
            	edge[0] = edge[1] = edge[2] = 0.3;
            	FT_MakeCuboidSurf(front,center,edge,neg_comp,
		pos_comp,w_type,&surf);
	    }
	    if(strcmp(string,"cone") == 0 
	       || strcmp(string,"CONE") == 0)
            {
            	center[0] = center[1] = 0.5;
            	center[2] = 0.0;
            	slope = 0.5;
            	height = 0.5;
            	FT_MakeConeSurf(front,center,slope,height,
		neg_comp,pos_comp,w_type,&surf);
	    }
	    if(strcmp(string,"cylindercrosscylinder") == 0 
	       || strcmp(string,"CYLINDERCROSSCYLINDER") == 0)
	    {
             	center[0] = center[1] = center[2] = 0.5;
             	center1[0] = center1[1] = center1[2] = 0.5;
             	r = 0.2;
             	R = 0.1;
             	height = 0.4;
             	height1 = 0.4;
             	FT_MakeCrossCylinderSurf(front,center,
		center1,r,R,height,height1,neg_comp,pos_comp,
		w_type,&surf);
	    }
	break;	
	case 'd':
	case 'D':
	    x0 = 0.25;
	    x1 = 0.75;
	    y = z = 0.25;
	    R = 0.15;
            r = 0.075;
	    FT_MakeDumbBellSurf(front,x0,x1,y,z,R,r,neg_comp,
	    pos_comp,w_type,&surf);
	break;
	case 'p':
	case 'P':
	    if(strcmp(string,"projectile") == 0 
	       || strcmp(string,"PROJECTILE") == 0)
	    {
		center[0] = center[1] = center[2] = 0.5; 
		R = 0.075; 
	    	r = 0.15;  
	    	height = 0.2;  
	    	FT_MakeProjectileSurf(front,center,R,r,height,
		neg_comp,pos_comp,w_type,&surf);
	    }
	    if(strcmp(string,"platform") == 0 
	       || strcmp(string,"PLATFORM") == 0)
	    {
		center[0] = center[1] = 0.5;
		center[2] = 0.2; 
		R = 0.49;  
		height = 0.4;
		slope = 0.5;
		FT_MakePlatformSurf(front,center,R,height,slope,
		neg_comp,pos_comp,w_type,&surf);
	    } 
 	    if(strcmp(string,"plane") == 0 
	       || strcmp(string,"PLANE") == 0)
	    {
             	nor[0] = -1.0;
            	nor[1] = nor[2] = 1.0;
             	pt[0] = pt[1] = pt[2] = 0.5;
             	FT_MakePlaneSurf(front,nor,pt,NO,neg_comp,
		pos_comp,w_type,&surf);
	    }
	break;
        case 't':
 	case 'T':
             center[0] = center[1] = center[2] = 1;
	     edgel = 2;
	     FT_MakeTetrahedronSurf(front,center,edgel,neg_comp,
	     pos_comp,w_type,&surf);
	break;
	case 'b':
	case 'B':
	     center[0] = center[1] = 0.5;
	     center[2] = 0.8;
	     R = 0.4;
	     RR = 0.38;
	     r = 0.1;
	     height = 0.3;
	     height1 = 0.05;
	     FT_MakeBowlSurf(front,center,R,r,RR,height,height1,
	     neg_comp,pos_comp,w_type,&surf);
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


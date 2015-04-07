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
*	User initialization example for Front Package:
*	Copyright 1999 by The University at Stony Brook, 
*	All rights reserved.
*
*/

#include <FronTier.h>

int main(int argc, char **argv)
{
	static Front front;
	static F_BASIC_DATA f_basic;
	static LEVEL_FUNC_PACK level_func_pack;
	SURFACE *surf;
	// for sphere
	double center[MAXD];
	double radii[MAXD];
	// for dumbbell
	double x0,x1,y0,z0,R,r;
	// for projectile
	double h;
	// for test of rotation
	double phi,theta;
	// for cuboid
	double edge[MAXD];
	// for cylinder
	double radius,height;
	// for cone
	double slope;
	// for tetrahedron
	double tedge;

	FT_Init(argc,argv,&f_basic);

	/* Initialize basic computational data */

	f_basic.size_of_intfc_state = 0;

	FT_ReadSpaceDomain(f_basic.in_name,&f_basic);
	FT_StartUp(&front,&f_basic);
	FT_InitDebug(f_basic.in_name);

	level_func_pack.pos_component = 2;

	FT_InitIntfc(&front,&level_func_pack);

	// Test making elliptic surface
	center[0] = center[1] = center[2] = 0.5;
	radii[0] = 0.4;
	radii[1] = 0.3;
	radii[2] = 0.3;
	FT_MakeEllipticSurf(&front,center,radii,
		exterior_component(front.interf),2,FIRST_PHYSICS_WAVE_TYPE,
		1,&surf);
	FT_Save(&front);
	delete_surface(surf);
	FT_AddTimeStepToCounter(&front);

	// Test making dumbbell surface
	x0 = 0.25;
	x1 = 0.75;
	y0 = 0.5;
	z0 = 0.5;
	R = 0.15;
	r = 0.075;
	FT_MakeDumbBellSurf(&front,x0,x1,y0,z0,R,r,
		exterior_component(front.interf),2,FIRST_PHYSICS_WAVE_TYPE,
		&surf);
	FT_Save(&front);
	FT_AddTimeStepToCounter(&front);
	// Test rotate the dumbbell surface
	phi = 30.0/180.0*PI;
	theta = 40.0/180.0*PI;
	FT_RotateSurface(surf,center,phi,theta);
	FT_Save(&front);
	FT_AddTimeStepToCounter(&front);
	delete_surface(surf);

	// Test making projectile surface
	center[0] = center[1] = center[2] = 0.5;
	R = 0.2;
	r = 0.3;
	h = 0.2;
	FT_MakeProjectileSurf(&front,center,R,r,h,
		exterior_component(front.interf),2,FIRST_PHYSICS_WAVE_TYPE,
                &surf);
	FT_Save(&front);
	FT_AddTimeStepToCounter(&front);
	// Test rotate the projectile surface
	phi = 30.0/180.0*PI;
	theta = 40.0/180.0*PI;
	FT_RotateSurface(surf,center,phi,theta);
	FT_Save(&front);
	FT_AddTimeStepToCounter(&front);
	delete_surface(surf);

	// Test making cuboid surface
	center[0] = center[1] = center[2] = 0.5;
	edge[0] = 0.3;
	edge[1] = 0.1;
	edge[2] = 0.1;
	FT_MakeCuboidSurf(&front,center,edge,
		exterior_component(front.interf),2,FIRST_PHYSICS_WAVE_TYPE,
		&surf);
	FT_Save(&front);
	delete_surface(surf);
	FT_AddTimeStepToCounter(&front);

	// Test making cylinder surface
	center[0] = center[1] = center[2] = 0.5;
	radius = 0.2;
	height = 0.3;
	FT_MakeCylinderSurf(&front,center,radius,height,
		exterior_component(front.interf),2,FIRST_PHYSICS_WAVE_TYPE,
		&surf);
	FT_Save(&front);
	delete_surface(surf);
	FT_AddTimeStepToCounter(&front);

	// Test making cone surface
	center[1] = center[2] = 0.5; 
	center[0] = 0.25;
	slope = 0.6;
	height = 0.5;
	FT_MakeConeSurf(&front,center,slope,height,
		exterior_component(front.interf),2,FIRST_PHYSICS_WAVE_TYPE,
                &surf);
	FT_Save(&front);
	delete_surface(surf);
	FT_AddTimeStepToCounter(&front);

	//Test making tetrahedron surface
	center[0] = center[1] = center[2] = 0.5;
        tedge = 0.4;
        FT_MakeTetrahedronSurf(&front,center,tedge,
                exterior_component(front.interf),2,FIRST_PHYSICS_WAVE_TYPE,
                &surf);
        FT_Save(&front);
        delete_surface(surf);

	clean_up(0);
}

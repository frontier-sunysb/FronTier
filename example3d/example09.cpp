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
static void initInteriorSurfaces(Front*);

int main(int argc, char **argv)
{
	static Front front;
	static F_BASIC_DATA f_basic;
	static LEVEL_FUNC_PACK level_func_pack;

	FT_Init(argc,argv,&f_basic);

	/* Hard coded domain settint */

	f_basic.dim = 3;
        f_basic.L[0] = 0.0;     f_basic.L[1] = 0.0;     f_basic.L[2] = 0.0;
        f_basic.U[0] = 1.0;     f_basic.U[1] = 1.0;     f_basic.U[2] = 1.0;
        f_basic.gmax[0] = 100;  f_basic.gmax[1] = 100; f_basic.gmax[2] = 100;
        f_basic.boundary[0][0] = f_basic.boundary[0][1] = DIRICHLET_BOUNDARY;
        f_basic.boundary[1][0] = f_basic.boundary[1][1] = DIRICHLET_BOUNDARY;
        f_basic.boundary[2][0] = f_basic.boundary[2][1] = DIRICHLET_BOUNDARY;
	f_basic.size_of_intfc_state = 0;


	FT_StartUp(&front,&f_basic);

	level_func_pack.pos_component = 2; // default component
	FT_InitIntfc(&front,&level_func_pack);
	initInteriorSurfaces(&front);
	FT_Draw(&front);
	clean_up(0);
}

static void initInteriorSurfaces(
	Front *front)
{
	double center[MAXD],radius[MAXD];
	int neg_comp,pos_comp,w_type;
	SURFACE *surf;

	neg_comp = 1;
	pos_comp = 2;
	w_type = FIRST_PHYSICS_WAVE_TYPE;

	center[0] = center[1] = center[2] = 0.4;
	radius[0] = radius[1] = radius[2] = 0.1;
	FT_MakeEllipticSurf(front,center,radius,neg_comp,pos_comp,
			w_type,2,&surf);
	if (FT_CheckSurfCompConsistency(front,surf) == YES)
	{
	    printf("Surface component is consistent!\n");	
	}
	center[0] = center[1] = center[2] = 0.6;
	radius[0] = radius[1] = radius[2] = 0.1;
	FT_MakeEllipticSurf(front,center,radius,neg_comp,pos_comp,
			w_type,2,&surf);
	if (FT_CheckSurfCompConsistency(front,surf) == YES)
	{
	    printf("Surface component is consistent!\n");	
	}
}	/* end initInteriorSurfaces */

/************************************************************************************
FronTier is a set of libraries that implements differnt types of Front Traking algorithms.
Front Tracking is a numerical method for the solution of partial differential equations 
whose solutions have discontinuities.  


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
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

******************************************************************************/


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

static void map_output_interface(Front*);

#if defined(__MPI__)
int subdomains[MAXD];
#endif /* defined(__MPI__) */


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
	MC_PARAMS mc_params;
	static int count = 0;
	char redist_name[100];

	/* Initialize basic computational data */

	f_basic.dim = 2;
	FT_Init(argc,argv,&f_basic);

	f_basic.L[0] = 0.0;	f_basic.L[1] = 0.0;
	f_basic.U[0] = 1.0;	f_basic.U[1] = 1.0;
	f_basic.gmax[0] = 40;	f_basic.gmax[1] = 40;
	f_basic.boundary[0][0] = f_basic.boundary[0][1] = DIRICHLET_BOUNDARY;
	f_basic.boundary[1][0] = f_basic.boundary[1][1] = DIRICHLET_BOUNDARY;
	f_basic.size_of_intfc_state = 0;

	FT_StartUp(&front,&f_basic);

	/* Initialize interface through level function */

	mc_params.dim = 2;
	mc_params.num_cir = 1;
        FT_VectorMemoryAlloc((POINTER*)&mc_params.rad,mc_params.num_cir,
					FLOAT);
        FT_MatrixMemoryAlloc((POINTER*)&mc_params.cen,mc_params.num_cir,
					2,FLOAT);
	mc_params.cen[0][0] = 0.5;
	mc_params.cen[0][1] = 0.5;
	mc_params.rad[0] = 0.35;

	level_func_pack.neg_component = 1;
	level_func_pack.pos_component = 2;
	level_func_pack.func_params = (POINTER)&mc_params;
	level_func_pack.func = multi_circle_func;
	level_func_pack.wave_type = FIRST_PHYSICS_WAVE_TYPE;

	FT_InitIntfc(&front,&level_func_pack);
	    if (f_basic.dim < 3)
                FT_ClipIntfcToSubdomain(&front);

	/* Original mesh from marching cube method */
	/* Can be viewed by the xgraph:
	 * xgraph -P -bg white file.xg */

        sprintf(restart_name,"%s/intfc-ts%s",restart_name,
                        right_flush(RestartStep,7));
        if (pp_numnodes() > 1)
            sprintf(restart_name,"%s-nd%s",restart_name, 
                                right_flush(pp_mynode(),4));
	xgraph_interface_curves(".",redist_name,front.interf,XY_PLANE);

	/* Default redistribute */
	FT_RedistMesh(&front);
	sprintf(redist_name,"redist-%d",right_flush(count++,3));
#if defined(__MPI__)
	sprintf(redist_name,"%s-nd%s",redist_name,right_flush(pp_mynode(),4));
#endif /* defined(__MPI__) */
	xgraph_interface_curves(".",redist_name,front.interf,XY_PLANE);

	/* Finer: spacing = 0.1 of rectangular mesh size */
	FrontSetSpacing(&front,0.1);
	FT_RedistMesh(&front);
	sprintf(redist_name,"redist-%d",right_flush(count++,3));
#if defined(__MPI__)
	sprintf(redist_name,"%s-nd%s",redist_name,right_flush(pp_mynode(),4));
#endif /* defined(__MPI__) */
	xgraph_interface_curves(".",redist_name,front.interf,XY_PLANE);

	/* Coarser: spacing = 3.0 of rectangular mesh size */
	FrontSetSpacing(&front,3.0);
	FT_RedistMesh(&front);
	sprintf(redist_name,"redist-%d",right_flush(count++,3));
#if defined(__MPI__)
	sprintf(redist_name,"%s-nd%s",redist_name,right_flush(pp_mynode(),4));
#endif /* defined(__MPI__) */
	xgraph_interface_curves(".",redist_name,front.interf,XY_PLANE);

	clean_up(0);
}

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

static char *out_name;

int main(int argc, char **argv)
{
	static Front front;
	static RECT_GRID comp_grid;
	static F_BASIC_DATA f_basic;
	static LEVEL_FUNC_PACK level_func_pack;
	MC_PARAMS mc_params;

	FT_Init(argc,argv,&f_basic);
	f_basic.dim = 2;

	/* Initialize basic computational data */

	f_basic.L[0] = 0.0;	f_basic.L[1] = 0.0;
	f_basic.U[0] = 1.0;	f_basic.U[1] = 1.0;
	f_basic.gmax[0] = 40;	f_basic.gmax[1] = 40;
	f_basic.boundary[0][0] = f_basic.boundary[0][1] = DIRICHLET_BOUNDARY;
	f_basic.boundary[1][0] = f_basic.boundary[1][1] = DIRICHLET_BOUNDARY;
	f_basic.size_of_intfc_state = 0;

        out_name                = f_basic.out_name;

	FT_StartUp(&front,&f_basic);

	/* Initialize interface through level function */

	mc_params.dim = 2;
	mc_params.num_cir = 3;
        FT_VectorMemoryAlloc((POINTER*)&mc_params.rad,mc_params.num_cir,FLOAT);
        FT_MatrixMemoryAlloc((POINTER*)&mc_params.cen,mc_params.num_cir,2,FLOAT);
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

	map_output_interface(&front);

	clean_up(0);
}

static void map_output_interface(
	Front *front)
{
	INTERFACE *intfc = front->interf;
        FILE *out_file;
	int i,j,k,step;
	step = front->step;
	char filename[100];

	out_file = fopen(filename,"w");
	
	int num_curves = FT_NumOfIntfcCurves(intfc);
	int dim = Dimension(intfc);
	int num_nodes = FT_NumOfIntfcNodes(intfc);
	int num_bonds = FT_NumOfIntfcBonds(intfc);
	int num_points = FT_NumOfIntfcPoints(intfc);

	CURVE **curves;
	curves = (CURVE**)malloc(num_curves*sizeof(CURVE*));

	FT_ArrayOfIntfcCurves(intfc,curves);

	fprintf(out_file,"Interface at step %d\n\n",step);
	fprintf(out_file,"Dimension = %d\n",dim);
	fprintf(out_file,"Number of Curves = %d\n",num_curves);
	fprintf(out_file,"Number of Nodes = %d\n",num_nodes);
	fprintf(out_file,"Number of Interface Bonds = %d\n",num_bonds);
	fprintf(out_file,"Number of Interface Points = %d\n",num_points);
	for (i = 0; i < num_curves; ++i)
	{
	    double *coords;

	    num_bonds = FT_NumOfCurveBonds(curves[i]);
	    num_points = FT_NumOfCurvePoints(curves[i]);
	    coords = (double*)malloc(num_points*dim*sizeof(double));

	    ArrayOfCurvePoints(curves[i],coords);

	    fprintf(out_file,"Number of Bonds on Curve %d = %d\n",
	    			i+1,num_bonds);
	    fprintf(out_file,"Number of Points on Curve %d = %d\n\n",
	    			i+1,num_points);
	    for (j = 0; j < num_points; ++j)
	    {
	    	for (k = 0; k < dim; ++k)
		    fprintf(out_file,"%lf  ",coords[dim*j+k]);
	    	fprintf(out_file,"\n");
	    }
	    fprintf(out_file,"\n\n");
	}
	fclose(out_file);
}	/* end map_output */

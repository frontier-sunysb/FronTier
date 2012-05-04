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
*				test_map_3d.c:
*
*		User initialization example for Front Package:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	This is example of three circles all moving the a normal velocity.
*	Bifurcation occurs when they meet each other. FronTier solves
*	the bifurcation automatically.
*
*	This example shows how to convert the interface into arrays
*	of floats. Functions introduced in this example include:
*
*	Dimension(intfc);
*
* int     NumOfPoints(SURFACE *s);
* int     NumOfPoints(INTERFACE *intfc);
* int     NumOfNodes(INTERFACE *intfc);
* int     NumOfBonds(INTERFACE *intfc);
* int     NumOfCurves(INTERFACE *intfc);
* int     NumOfSurfaces(INTERFACE *intfc);
* int     NumOfTris(SURFACE *s);
* int     NumOfTris(INTERFACE *intfc);

* void    ArrayOfPoints(INTERFACE *intfc, double *coords);
* void    ArrayOfTri(SURFACE *surface, TRI **tris);
* void    ArrayOfTri(SURFACE *surface, double *coords, int *vertices_index);
* void    ArrayOfTri(INTERFACE *intfc, TRI **tris);
* void    ArrayOfTri(INTERFACE *intfc, double *coords, int *vertex_indices);
*
*	See their use in function map_output_interface().
*
*/

#include <FronTier.h>

	/*  Function Declarations */

static void  map_output_interface(Front*);

char *out_name;

/********************************************************************
 *	Level function parameters for the initial interface 	    *
 ********************************************************************/

typedef struct {
        int num_cir;
        double **cen;
        double *rad;
} TMC_PARAMS;

int main(int argc, char **argv)
{
	static Front front;
	static RECT_GRID comp_grid;
	static F_BASIC_DATA f_basic;
	static LEVEL_FUNC_PACK level_func_pack;
	TMC_PARAMS mc_params;

	f_basic.dim = 3;	
	FT_Init(argc,argv,&f_basic);

	/* Initialize basic computational data */

	f_basic.L[0] = 0.0;	f_basic.L[1] = 0.0;	f_basic.L[2] = 0.0;
	f_basic.U[0] = 1.0;	f_basic.U[1] = 1.0;	f_basic.U[2] = 1.0;
	f_basic.gmax[0] = 40;	f_basic.gmax[1] = 40;	f_basic.gmax[2] = 40;
	f_basic.boundary[0][0] = f_basic.boundary[0][1] = DIRICHLET_BOUNDARY;
	f_basic.boundary[1][0] = f_basic.boundary[1][1] = DIRICHLET_BOUNDARY;
	f_basic.boundary[2][0] = f_basic.boundary[2][1] = DIRICHLET_BOUNDARY;
	f_basic.size_of_intfc_state = 0;

        out_name                = f_basic.out_name;

	FT_StartUp(&front,&f_basic);

	/* Start testing */

	front.step = 0;


	/* Initialize interface through level function */

	mc_params.num_cir = 1;
        FT_VectorMemoryAlloc((POINTER*)&mc_params.rad,mc_params.num_cir,FLOAT);
        FT_MatrixMemoryAlloc((POINTER*)&mc_params.cen,mc_params.num_cir,3,FLOAT);
	mc_params.cen[0][0] = 0.4;
	mc_params.cen[0][1] = 0.4;
	mc_params.cen[0][2] = 0.4;
	mc_params.rad[0] = 0.3;

	level_func_pack.neg_component = 1;
	level_func_pack.pos_component = 2;
	level_func_pack.func_params = (POINTER)&mc_params;
	level_func_pack.func = multi_circle_func;
	level_func_pack.wave_type = FIRST_PHYSICS_WAVE_TYPE;
	
	FT_InitIntfc(&front,&level_func_pack);

	map_output_interface(&front);

	clean_up(0);
}

static void SaveAsTecplot(
	char *filename, 
	int nP, 
	double *coords, 
	int nTri, 
	int *vIndices)
{
	FILE *hfile = fopen(filename, "w");
	fprintf(hfile, "TITLE = %s \n", filename);
	fprintf(hfile, "VARIABLES = \"X\" \"Y\" \"Z\" \n");
	fprintf(hfile, "ZONE N=%d, E=%d, F=FEPOINT, ET=TETRAHEDRON \n", 
					nP, nTri);
	for(int i=0; i<nP; i++)
	    fprintf(hfile,"%f %f %f \n",coords[i*3+0],coords[i*3+1],
	    				coords[i*3+2]);
	for(int i=0; i<nTri; i++)
	    fprintf(hfile, "%d %d %d %d", vIndices[i*3+0], vIndices[i*3+1],
					vIndices[i*3+2], vIndices[i*3+0]);
	fclose(hfile);
		
}

static void map_output_interface(
	Front *front)
{
	INTERFACE *intfc = front->interf;
        char filename[100];
	int i,j,step;
	step = front->step;
			
	int num_points = NumOfIntfcPoints(intfc);
	int num_tris = NumOfIntfcTris(intfc);
	int dim = 3;
	double *coords = (double*)malloc(num_points*dim*sizeof(double));
	int *vertex_indices = (int*)malloc(num_tris*3*sizeof(int));
	TRI **tris = (TRI**)malloc(num_tris*sizeof(TRI*));
	const char *out_name = "output";
	
	// the interface
	printf("Number of interface points = %d\n",num_points);
	printf("Number of interface triangles = %d\n\n",num_tris);

	ArrayOfIntfcTris_FT(intfc,tris);	
	printf("Print first 10 triangles:\n\n");
	for (i = 0; i < 10; i++)
	    print_tri(tris[i],intfc);
	printf("\n\n");

	sprintf(filename,"%s-%d.intfc.plt",out_name,step);
	ArrayOfIntfcTris(intfc,coords,vertex_indices);	
	SaveAsTecplot(filename,num_points,coords,num_tris,vertex_indices);
	
	int num_surfaces = NumOfSurfaces(intfc);
	printf("Number of surfaces = %d\n",num_surfaces);
	SURFACE **surfaces;
	surfaces = (SURFACE**)malloc(num_surfaces*sizeof(SURFACE*));
	
	ArrayOfSurfaces(intfc, surfaces);
	// print the tris of each surface
	for(i = 0; i < num_surfaces; i++)
	{
	    sprintf(filename,"%s-%d.surface-%d.plt",out_name,step,i);
	    num_points = NumOfSurfPoints(surfaces[i]);
	    num_tris   = NumOfSurfTris(surfaces[i]);
	    printf("Number of points on surface %d = %d\n",
	    			i+1,num_points);
	    printf("Number of tris on surface %d = %d\n",
	    			i+1,num_tris);

	    ArrayOfSurfTris_FT(surfaces[i],tris);	
	    printf("Print first 10 triangles:\n\n");
	    for (j = 0; j < 10; j++)
	        print_tri(tris[i],intfc);
	    printf("\n\n");

	    ArrayOfSurfTris(surfaces[i],coords,vertex_indices);
	    SaveAsTecplot(filename,num_points,coords,num_tris,vertex_indices);
	}
}	/* end map_output */

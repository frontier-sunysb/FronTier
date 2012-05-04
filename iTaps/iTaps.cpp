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
*				iTaps.c:
*
*		User initialization example for Front Package:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*	
*	This example shows a circle in a double vortex field. It demonstrates
*	the resolution of the front tracking method.
*
*/

#include <vector>
#include <FronTier.h>

	/*  Function Declarations */
static void iMesh_test(Front*);

char *in_name,*restart_state_name,*restart_name,*out_name;
boolean RestartRun;
int RestartStep;
boolean binary = YES;

/********************************************************************
 *	Level function parameters for the initial interface 	    *
 ********************************************************************/

typedef struct {
        /* equation for line is x^2/a^2 + y^2/b^2 = 1 */
        double x0;
        double y0;         
	double r;
        double w;         
	double h;
} DISK_PARAMS;


int main(int argc, char **argv)
{
	static Front front;
	static RECT_GRID comp_grid;
	static F_BASIC_DATA f_basic;
	static LEVEL_FUNC_PACK level_func_pack;
	DISK_PARAMS disk_params;	/* level function parameters */

	f_basic.dim = 2;	
	FT_Init(argc,argv,&f_basic);

	/* Initialize basic computational data */

	f_basic.L[0] = 0.0;	f_basic.L[1] = 0.0;
	f_basic.U[0] = 1.0;	f_basic.U[1] = 1.0;
	f_basic.gmax[0] = 128;	f_basic.gmax[1] = 128;
	f_basic.boundary[0][0] = f_basic.boundary[0][1] = PERIODIC_BOUNDARY;
	f_basic.boundary[1][0] = f_basic.boundary[1][1] = PERIODIC_BOUNDARY;
	f_basic.size_of_intfc_state = 0;

        in_name                 = f_basic.in_name;
        restart_state_name      = f_basic.restart_state_name;
        out_name                = f_basic.out_name;
        restart_name            = f_basic.restart_name;
        RestartRun              = f_basic.RestartRun;
        RestartStep             = f_basic.RestartStep;

        sprintf(restart_name,"%s.ts%s",restart_name,right_flush(RestartStep,7));
#if defined(__MPI__)
        sprintf(restart_name,"%s-nd%s",restart_name,right_flush(pp_mynode(),4));
#endif /* defined(__MPI__) */

	FT_StartUp(&front,&f_basic);

	if (!RestartRun)
	{
	    /* Initialize interface through level function */

	    disk_params.x0 = 0.5;
	    disk_params.y0 = 0.5;
	    disk_params.r = 0.3;
	    disk_params.w = 0.01;
	    disk_params.h = 0.4;

	    level_func_pack.neg_component = 1;
	    level_func_pack.pos_component = 2;
	    level_func_pack.func_params = (POINTER)&disk_params;
	    level_func_pack.func = slotted_disk_func;
	    level_func_pack.wave_type = FIRST_PHYSICS_WAVE_TYPE;
	    FT_InitIntfc(&front,&level_func_pack);
	}
        redistribute(&front,YES,NO);
        FT_AddMovieFrame(&front,out_name,binary);
	FT_Save(&front,out_name);
	
#ifdef IMESH
	iMesh_test(&front);
#endif /*ifdef IMESH */

	clean_up(0);
}	/* end main */

#ifdef IMESH
#define		MESH(m)		(iMesh_Instance*)(m)

static void iMesh_test(Front *front)
{
	iMesh_Instance instance;
	char options[10],name[100];
	int dim,err,name_len,options_len;
	iBase_EntitySetHandle entity_set_handle, entity_set_handle2;
	iBase_TagHandle tag_handle,tag_handle2;
	iBase_EntityHandle entity_handle;
	INTERFACE *intfc = front->interf;
	CURVE **c,*curve,*curve2,*curve3;
	POINT **ps, *p;
	BOND *b;
	FILE *file;
	int is_contained,tag_value_allocated,tag_value_size;

	sprintf(options,"2");
	options_len = strlen(options);
	iMesh_newMesh(options,&instance,&err,options_len);
	print_error_enum(err);

	iMesh_getGeometricDimension(instance,&dim,&err);
	printf("Dimension = %d\n",dim);
	print_error_enum(err);

	sprintf(name,"%s/intfc-ts0000000",front->out_name);
	name_len = strlen(name);
	iMesh_load(instance,entity_set_handle,
			name,options,&err,name_len,0);
	print_error_enum(err);

	sprintf(name,"intfc_save");
	name_len = strlen(name);
	iMesh_save(instance,entity_set_handle,
			name,options,&err,name_len,0);
	print_error_enum(err);

	iMesh_dtor(instance,&err);
	print_error_enum(err);

	iMesh_createEntSet(instance,NO,&entity_set_handle,&err);
	print_error_enum(err);
	iMesh_destroyEntSet(instance,entity_set_handle,&err);
	print_error_enum(err);
	for (c = intfc->curves; c && *c; ++c)
	    if (wave_type(*c) == FIRST_PHYSICS_WAVE_TYPE)
		curve = *c;
	//print_curve(curve);
	file = fopen("the_curve","w");
	xgraph_curve(file,curve,XY_PLANE);
	b = curve->first;
	p = b->start;

	/* Test entity set as interface */
	printf("\nTesting entity set as intfc\n");
	iMesh_createEntSet(instance,NO,&entity_set_handle,&err);
	print_error_enum(err);
	iMesh_createTag(instance,"intfc",sizeof(POINTER),iBase_BYTES,
		&tag_handle,&err,5);
	print_error_enum(err);
	iMesh_setEntSetData(instance,entity_set_handle,tag_handle,
			(const char*)intfc,sizeof(POINTER),&err);
	print_error_enum(err);

	entity_handle = entityOfBond(b);
	iMesh_isEntContained(instance,entity_set_handle,entity_handle,
			&is_contained,&err);
	print_error_enum(err);
	printf("Bond b is_contained = %s\n",(is_contained ? "YES" : "NO"));
	iMesh_deleteEnt(instance,entity_handle,&err);
	print_error_enum(err);

	entity_handle = entityOfPoint(p);
	iMesh_isEntContained(instance,entity_set_handle,entity_handle,
			&is_contained,&err);
	print_error_enum(err);
	printf("Point p is_contained = %s\n",(is_contained ? "YES" : "NO"));
	iMesh_deleteEnt(instance,entity_handle,&err);
	print_error_enum(err);
	iMesh_destroyEntSet(instance,entity_set_handle,&err);
	print_error_enum(err);
	

	/* Test entity set as curve */
	printf("\nTesting entity set as curve\n");
	iMesh_createEntSet(instance,NO,&entity_set_handle,&err);
	print_error_enum(err);
	iMesh_createTag(instance,"curve",sizeof(POINTER),iBase_BYTES,
		&tag_handle,&err,5);
	print_error_enum(err);
	iMesh_setEntSetData(instance,entity_set_handle,tag_handle,
			(const char*)curve,sizeof(POINTER),&err);
	print_error_enum(err);

	entity_handle = entityOfBond(b);
	iMesh_isEntContained(instance,entity_set_handle,entity_handle,
			&is_contained,&err);
	print_error_enum(err);
	printf("Bond b is_contained = %s\n",(is_contained ? "YES" : "NO"));
	iMesh_deleteEnt(instance,entity_handle,&err);
	print_error_enum(err);

	entity_handle = entityOfPoint(p);
	iMesh_isEntContained(instance,entity_set_handle,entity_handle,
			&is_contained,&err);
	print_error_enum(err);
	printf("Point p is_contained = %s\n",(is_contained ? "YES" : "NO"));
	iMesh_deleteEnt(instance,entity_handle,&err);
	print_error_enum(err);
	printf("Test step 3: curve = %d  intfc = %d\n",curve,intfc);

	/*Test iMesh_getEntSetData*/
	tag_value_allocated = YES;
	tag_value_size = sizeof(POINTER);
	iMesh_getEntSetData(instance,entity_set_handle,tag_handle,
			(char**)&curve2,&tag_value_allocated,&tag_value_size,
			&err);
	print_error_enum(err);
	printf("Test step 4: curve = %d  intfc = %d\n",curve,intfc);
	iMesh_destroyEntSet(instance,entity_set_handle,&err);
	print_error_enum(err);

	/*Test EntitySet Contained*/
	printf("\nTesting entity set as point_set\n");
	iMesh_createEntSet(instance,NO,&entity_set_handle,&err);
	print_error_enum(err);
	iMesh_createTag(instance,"intfc",sizeof(POINTER),iBase_BYTES,
		&tag_handle,&err,5);
	print_error_enum(err);
	printf("Test step 5: curve = %d  intfc = %d\n",curve,intfc);
	iMesh_setEntSetData(instance,entity_set_handle,tag_handle,
			(const char*)intfc,sizeof(POINTER),&err);
	print_error_enum(err);

//	iMesh_destroyEntSet(instance,entity_set_handle,&err);
//	ps = intfc->points;
//	if (ps == NULL)
//	    printf("the ps is NULL\n");
//	else
//	    printf("After getting the point set\n");
//	iMesh_createEntSet(instance,NO,&entity_set_handle2,&err);
//	print_error_enum(err);
//	iMesh_createTag(instance,"points",sizeof(POINTER),iBase_BYTES,
//		&tag_handle2,
//		&err,6);
//	print_error_enum(err);
//	iMesh_setEntSetData(instance,entity_set_handle2,tag_handle2,
//		(const char*)*ps,sizeof(POINTER),&err);
//	print_error_enum(err);
//	iMesh_isEntSetContained(instance,entity_set_handle,
//		entity_set_handle2,&is_contained,&err);
//	print_error_enum(err);

	printf("\nTesting entity set as curve_set\n");
//	iMesh_createEntSet(instance,NO,&entity_set_handle2,&err);
//	print_error_enum(err);
//	iMesh_createTag(instance,"curve",sizeof(POINTER),iBase_BYTES,
//		&tag_handle2,&err,5);
//	print_error_enum(err);
//	iMesh_setEntSetData(instance,entity_set_handle2,tag_handle2,
//		(const char*)*c,sizeof(POINTER),&err);
//	print_error_enum(err);
//	printf("the type for hand is %d\n", 
//		((FT_ESET_HANDLE *)entity_set_handle)->type);
//	printf("the type for hand_contained is %d\n", 
//		((FT_ESET_HANDLE *)entity_set_handle2)->type);
//	iMesh_isEntSetContained(instance,entity_set_handle,
//		entity_set_handle2,&is_contained,&err);
//	print_error_enum(err);
	
	iMesh_createEntSet(instance,NO,&entity_set_handle2,&err);
	print_error_enum(err);
	printf("Test step 6: curve = %d  intfc = %d\n",curve,intfc);
	iMesh_createTag(instance,"curve",sizeof(POINTER),iBase_BYTES,
		&tag_handle2,&err,5);
	print_error_enum(err);
	printf("Test 10: curve = %d\n",curve);
	iMesh_setEntSetData(instance,entity_set_handle2,tag_handle2,
		(const char*)curve,sizeof(POINTER),&err);
	print_error_enum(err);

	/*
	for (c = intfc->curves; c && *c; ++c)
	    if (wave_type(*c) != FIRST_PHYSICS_WAVE_TYPE)
		 curve2 = *c;
        iMesh_createEntSet(instance,NO,&entity_set_handle2,&err);
	print_error_enum(err);
	iMesh_createTag(instance,"curve",sizeof(POINTER),iBase_BYTES,
		&tag_handle2,&err,5);
	print_error_enum(err);
	iMesh_setEntSetData(instance,entity_set_handle2,tag_handle2,
		(const char*)curve2,sizeof(POINTER),&err);
	print_error_enum(err);

	printf("the type for entity_set_to add is %d\n", 
		((FT_ESET_HANDLE*)entity_set_handle2)->type);
	printf("the type for entity_set_handle is %d\n",
		((FT_ESET_HANDLE*)entity_set_handle)->type);
	iMesh_addEntSet(instance,entity_set_handle2,entity_set_handle,
		&err);
	print_error_enum(err);

	*/
	printf("Test step 7: curve = %d  intfc = %d\n",curve,intfc);
	iMesh_isEntSetContained(instance,entity_set_handle,
		entity_set_handle2,&is_contained,&err);
	print_error_enum(err);
	printf("curve is_contained = %s\n",
		(is_contained ? "YES" : "NO")); 	
}	/* end iMesh_test */
#endif /*ifdef IMESH */

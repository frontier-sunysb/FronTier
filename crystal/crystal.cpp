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
*				crystal.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*/

#include "crystal.h"
#include "crystal_basic.h"

#define		MAX_NUM_VERTEX_IN_CELL		20

	/*  Local Application Function Declarations */

static void 	solute_main_driver(Front*,C_CARTESIAN&);

char *in_name,*restart_state_name,*restart_name,*out_name;
boolean RestartRun;
boolean ReadFromInput;
int RestartStep;

int main(int argc, char **argv)
{
	static Front front;
	static F_BASIC_DATA f_basic;
	static LEVEL_FUNC_PACK level_func_pack;
	static VELO_FUNC_PACK velo_func_pack;
	static CRT_PARAMS cRparams;

	C_CARTESIAN       c_cartesian(front);

	FT_Init(argc,argv,&f_basic);

	/* Initialize basic computational data */

	f_basic.size_of_intfc_state = sizeof(STATE);
	
	//Initialize Petsc
	PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
	//if (debugging("trace")) printf("Passed PetscInitialize()\n");

	in_name      		= f_basic.in_name;
	restart_state_name      = f_basic.restart_state_name;
        out_name     		= f_basic.out_name;
        restart_name 		= f_basic.restart_name;
        RestartRun   		= f_basic.RestartRun;
        ReadFromInput   	= f_basic.ReadFromInput;
	RestartStep 		= f_basic.RestartStep;

	sprintf(restart_state_name,"%s/state.ts%s",restart_name,
			right_flush(RestartStep,7));
	sprintf(restart_name,"%s/intfc-ts%s",restart_name,
			right_flush(RestartStep,7));
#if defined(__MPI__)
	if (pp_numnodes() > 1)
	{
            sprintf(restart_name,"%s-nd%s",restart_name,
				right_flush(pp_mynode(),4));
            sprintf(restart_state_name,"%s-nd%s",restart_state_name,
				right_flush(pp_mynode(),4));
	}
#endif /* defined(__MPI__) */
	if (!ReadFromInput)
	{
	    (void) printf("ERROR: Input file needed!\n");
	    clean_up(ERROR);
	}
	FT_ReadSpaceDomain(in_name,&f_basic);

	FT_StartUp(&front,&f_basic);
	FT_InitDebug(in_name);

	cRparams.dim = f_basic.dim;
	front.extra2 = (POINTER)&cRparams;

	if (!RestartRun)
	{
	    /* Initialize interface through level function */
	    setInitialIntfc(&front,&level_func_pack,in_name);
	    FT_InitIntfc(&front,&level_func_pack);
	}
	else
	    read_restart_params(f_basic.dim,in_name,&front);

	read_crt_dirichlet_bdry_data(in_name,&front,f_basic);

	if (f_basic.dim == 2)
	    FT_ClipIntfcToSubdomain(&front);
	if (debugging("trace"))
	{
	    if (f_basic.dim == 2)
	    {
		char xg_name[100];
		sprintf(xg_name,"init_intfc-%d",pp_mynode());
		xgraph_2d_intfc(xg_name,front.interf);
	    }
	    else if (f_basic.dim == 3)
	    {
		char dname[100];
		sprintf(dname,"init_intfc-%s",right_flush(pp_mynode(),4));
		gview_plot_color_interface(dname,front.interf,YES);
	    }
	}

	read_crystal_params(in_name,&cRparams);
        if (debugging("trace")) printf("Passed read_crystal_params()\n");
	FT_ReadTimeControl(in_name,&front);

	/* Initialize velocity field function */

	
	velo_func_pack.func_params = (POINTER)&cRparams;
	velo_func_pack.func = NULL;

	FT_InitVeloFunc(&front,&velo_func_pack);

        c_cartesian.initMesh();
        c_cartesian.initMovieVariables();
	if (debugging("sample_solute"))
            c_cartesian.initSampleSolute(in_name);

        front._point_propagate = crystal_point_propagate;

	if (RestartRun)
	{
	    c_cartesian.readFrontInteriorStates(restart_state_name);
	}
	else
	{
	    initFrontStates(&front);
	    c_cartesian.setInitialCondition();
	}
	cRparams.field->vel = NULL;

	/* Propagate the front */

	solute_main_driver(&front,c_cartesian);

	PetscFinalize();
	clean_up(0);
}

static  void solute_main_driver(
        Front *front,
	C_CARTESIAN &c_cartesian)
{
        double CFL;
        int  i,dim = front->rect_grid->dim;
	CRT_PARAMS *cRparams = (CRT_PARAMS*)front->extra2;
	double h_min;
	double frac_dim,radius;
	FILE *Radius_file,*FracDim_file,*FracDim_vs_Damkl;
	boolean bdry_reached = NO;
	double **growth_data = NULL;
	int count = 0;
	double D = cRparams->D;
	double k = cRparams->k;
	Radius_file = FracDim_file = FracDim_vs_Damkl = NULL;
	SEED_PARAMS *s_params = (SEED_PARAMS*)cRparams->func_params;

	
	if (debugging("trace"))
	    printf("Entering solute_main_driver()\n");
	Curve_redistribution_function(front) = expansion_redistribute;
        CFL = Time_step_factor(front);

	h_min = front->rect_grid->h[0];
	for (i = 1; i < dim; ++i)
	    h_min = std::min(h_min,front->rect_grid->h[i]);

        FT_Draw(front);
	if (dim == 2 && front->vtk_movie_var)
            c_cartesian.vtk_plot_concentration2d(out_name);
	if (dim == 1)
            c_cartesian.oneDimPlot(out_name);
	if (debugging("trace"))
	    printf("Passed FT_Draw()\n");

	if (dim == 1)
	{
	    bdry_reached = NO;
	    record_1d_growth(front,&growth_data,&count);
	}
	else if (s_params != NULL)      //Interface_type != SEED
	{
	    bdry_reached = fractal_dimension(front,*s_params,&frac_dim,
						&radius);
	}
	if (pp_mynode() == 0 && dim != 1)
	{
	    char fname[200];
	    sprintf(fname,"%s/radius",out_name);
	    Radius_file = fopen(fname,"w");
	    (void) printf("Radius_file = %p\n",(POINTER)Radius_file);
	    sprintf(fname,"%s/FracDim",out_name);
	    FracDim_file = fopen(fname,"w");
	    sprintf(fname,"%s/DimDam",out_name);
	    FracDim_vs_Damkl = fopen(fname,"w");

	    fprintf(Radius_file,"\"Crystal radius\"\n\n");
	    fprintf(Radius_file,"%f  %f\n",front->time,radius);
	    fprintf(FracDim_file,"\"Fractal dimension\"\n\n");
	    fprintf(FracDim_file,"%f  %f\n",front->time,frac_dim);
	    fprintf(FracDim_vs_Damkl,"\"Fractal dimension vs Damkhl\"\n\n");
	    fprintf(FracDim_vs_Damkl,"%f  %f\n",k*radius/D,frac_dim);
	    fflush(Radius_file);
	    fflush(FracDim_file);
	    fflush(FracDim_vs_Damkl);
	}

        if (!RestartRun)
        {
	    if (debugging("sample"))
	    {
		cRparams->max_solute = -HUGE;	
		cRparams->min_solute =  HUGE;	
	    }
	    FT_ResetTime(front);
	    FT_Propagate(front);
	    c_cartesian.solve(front->dt);
	    if (debugging("sample"))
	    {
		printf("Front: max_solute = %f  min_solute = %f\n",
			cRparams->max_solute,cRparams->min_solute);
	    }
	    FT_SetTimeStep(front);
	    c_cartesian.setAdvectionDt();
	    front->dt = std::min(front->dt,CFL*c_cartesian.max_dt);
            FT_Save(front);

	    if (dim == 1)
	    	c_cartesian.oneDimPlot(out_name);
	    FT_SetOutputCounter(front);

        }
        else
        {
	    FT_SetOutputCounter(front);
        }
	if (debugging("trace"))
	    printf("Passed second restart check()\n");

	FT_TimeControlFilter(front);
	FT_PrintTimeStamp(front);

        for (;;)
        {
	    if (debugging("sample"))
	    {
		cRparams->max_solute = -HUGE;	
		cRparams->min_solute =  HUGE;	
	    }

	    FT_Propagate(front);

	    if (debugging("sample"))
	    {
		printf("Front: max_solute = %f  min_solute = %f\n",
			cRparams->max_solute,cRparams->min_solute);
	    }
	    c_cartesian.solve(front->dt);

	    if (debugging("trace"))
	    {
		print_storage("At end of time step","trace");
	    }

	    FT_AddTimeStepToCounter(front);

	    FT_SetTimeStep(front);
	    if (debugging("step_size"))
                printf("Time step from FrontHypTimeStep(): %f\n",front->dt);
	    front->dt = std::min(front->dt,CFL*c_cartesian.max_dt);
	    if (debugging("step_size"))
                printf("Time step from c_cartesian.max_dt(): %f\n",front->dt);

            if (FT_IsSaveTime(front))
	    {
		// Front standard output
                FT_Save(front);
		// Problem specific output
		c_cartesian.printFrontInteriorStates(out_name);
		if (debugging("trace")) printf("Passed Crystal "
                                "printFrontInteriorStates()\n");
	    }
            if (FT_IsDrawTime(front))
	    {
		// Front standard output
                FT_Draw(front);
		if (dim == 2 && front->vtk_movie_var)
                    c_cartesian.vtk_plot_concentration2d(out_name);
		if (dim == 1)
	    	    c_cartesian.oneDimPlot(out_name);
	    }
	    if (dim == 1)
	    {
	    	bdry_reached = NO;
		record_1d_growth(front,&growth_data,&count);
	    }
	    else if (s_params != NULL)
	    {
	    	bdry_reached = fractal_dimension(front,*s_params,&frac_dim,
					&radius);
	    }
	    if (debugging("trace"))
		printf("After auxilary output()\n");

	    if (bdry_reached) front->time_limit_reached = YES;
	    if (s_params != NULL)
	    {
		if (pp_mynode() == 0 && dim != 1 && !s_params->grow_from_floor)
	        {
		    fprintf(Radius_file,"%f  %f\n",front->time,radius);
		    fprintf(FracDim_file,"%f  %f\n",front->time,frac_dim);
		    fprintf(FracDim_vs_Damkl,"%f  %f\n",k*radius/D,frac_dim);
		    fflush(Radius_file);
		    fflush(FracDim_file);
		    fflush(FracDim_vs_Damkl);
	        }
	    }

            if (FT_TimeLimitReached(front))
	    {
	    	FT_PrintTimeStamp(front);
		if (dim == 1)
		    plot_growth_data(out_name,growth_data,count);
                break;
	    }
	    if (debugging("trace"))
		printf("After time output()\n");
	    /* Output section, next dt may be modified */

	    FT_TimeControlFilter(front);
	    FT_PrintTimeStamp(front);
        }
	if (s_params != NULL)
	{
	    if (pp_mynode() == 0 && dim != 1 && !s_params->grow_from_floor)
	    {
	        fclose(Radius_file);
	        fclose(FracDim_file);
	        fclose(FracDim_vs_Damkl);
	    }
        } 
}       /* end solute_main_driver */

void solute_read_front_states(
	FILE *infile,
	Front *front)
{
        INTERFACE *intfc = front->interf;
        STATE *sl,*sr;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
        double x;


        /* Initialize states at the interface */
        next_output_line_containing_string(infile,"Interface solute states:");
        next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
            FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
            fscanf(infile,"%lf",&x);
            sl->solute = x;
            fscanf(infile,"%lf",&x);
            sr->solute = x;
        }
}	/* end solute_read_front_states */

void solute_print_front_states(
	FILE *outfile,
	Front *front)
{
	INTERFACE *intfc = front->interf;
	STATE *sl,*sr;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;

        /* Initialize states at the interface */
        fprintf(outfile,"Interface solute states:\n");
        int count = 0;
        next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
            count++;
            FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
            fprintf(outfile,"%24.18g %24.18g\n",sl->solute,sr->solute);
        }

}	/* end solute_print_front_states */

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
*				crystal.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*/

#include <iFluid.h>
#include <crystal.h>
#include "subsurf.h"

	/*  Local Application Function Declarations */

static void 	subsurf_main_driver(Front*,C_CARTESIAN&,
			Incompress_Solver_Smooth_Basis*);
static void	solute_point_propagate(Front*,POINTER,POINT*,POINT*,	
			HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
static void     print_area_density(Front*,char*);

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
	static IF_PARAMS iFparams;
	static CRT_PARAMS cRparams;

	FT_Init(argc,argv,&f_basic);
	f_basic.size_of_intfc_state = sizeof(STATE);

	/* Construct the c_cartesian and l_cartesian class for crystal and
	 * incompressible fluid computation.
        */
	
	//Initialize Petsc
	PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
	//if (debugging("trace")) printf("Passed PetscInitialize()\n");


	C_CARTESIAN  c_cartesian(front);
	Incompress_Solver_Smooth_Basis *l_cartesian = NULL;
	if(f_basic.dim == 2)
	    l_cartesian = new Incompress_Solver_Smooth_2D_Cartesian(front);
	else if(f_basic.dim == 3)
	    l_cartesian = new Incompress_Solver_Smooth_3D_Cartesian(front);

	/* Initialize basic computational data */

	in_name      		= f_basic.in_name;
	restart_state_name      = f_basic.restart_state_name;
        out_name     		= f_basic.out_name;
        restart_name 		= f_basic.restart_name;
        RestartRun   		= f_basic.RestartRun;
        ReadFromInput   	= f_basic.ReadFromInput;
	RestartStep 		= f_basic.RestartStep;

	sprintf(restart_state_name,"%s/state.ts%s",restart_name,
			right_flush(RestartStep,7));
	sprintf(restart_name,"%s/intfc-ts%s",
			restart_name,right_flush(RestartStep,7));
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
	iFparams.dim = f_basic.dim;
	front.extra1 = (POINTER)&iFparams;
	front.extra2 = (POINTER)&cRparams;

	if (!RestartRun)
	{
	    if (f_basic.dim == 1)
	    {
		(void) printf("Subsurface problem has no 1D version!\n");
		clean_up(ERROR);
	    }

	    /* Initialize interface through level function */
	    setInitialIntfc(&front,&level_func_pack,in_name);

	    FT_InitIntfc(&front,&level_func_pack);
	    read_ss_dirichlet_bdry_data(in_name,&front,f_basic);
	    if (f_basic.dim != 3)
	    	FT_ClipIntfcToSubdomain(&front);
	}
	else
	{
	    read_ss_dirichlet_bdry_data(in_name,&front,f_basic);
	    read_restart_params(f_basic.dim,in_name,&front);
	}

	read_crystal_params(in_name,&cRparams);
	read_fluid_params(in_name,&iFparams);
	FT_ReadTimeControl(in_name,&front);

	/* Initialize velocity field function */

	velo_func_pack.func_params = (POINTER)&iFparams;
	velo_func_pack.func = NULL;

	FT_InitFrontVeloFunc(&front,&velo_func_pack);
        c_cartesian.initMesh();
        c_cartesian.initMovieVariables();
        l_cartesian->initMesh();
        l_cartesian->initMovieVariables();
	l_cartesian->findStateAtCrossing = ifluid_find_state_at_crossing;
	if (debugging("sample_velocity"))
            l_cartesian->initSampleVelocity(in_name);

        /* For geometry-dependent velocity, use first
        * order point propagation function, higher order
        * propagation requires surface propagate, currently
        * in writing, not yet in use. The following override
        * the assigned fourth_order_point_propagate.
        */

	if (debugging("sample_solute"))
            c_cartesian.initSampleSolute(in_name);
        front._point_propagate = solute_point_propagate;

	if (RestartRun)
	{
	    c_cartesian.readFrontInteriorStates(restart_state_name);
	    FT_FreeGridIntfc(&front);
	    l_cartesian->readFrontInteriorStates(restart_state_name);
	}
	else
	{
	    initFrontStates(&front);
	    c_cartesian.setInitialCondition();
	    FT_FreeGridIntfc(&front);
	    init_fluid_state_func(&front,l_cartesian);
	    l_cartesian->setInitialCondition();
	}
	cRparams.field->vel = iFparams.field->vel;

	/* Propagate the front */

	subsurf_main_driver(&front,c_cartesian,l_cartesian);

	PetscFinalize();
	clean_up(0);
}

static  void subsurf_main_driver(
        Front *front,
	C_CARTESIAN &c_cartesian,
	Incompress_Solver_Smooth_Basis *l_cartesian)
{
        int  dim = front->rect_grid->dim;
	boolean bdry_reached = NO;
	double CFL;

	if (debugging("trace"))
	    printf("Entering subsurf_main_driver()\n");

	Curve_redistribution_function(front) = expansion_redistribute;
        CFL = Time_step_factor(front);

	if (!RestartRun)
	{
	    FT_ResetTime(front);
	    FT_Propagate(front);
	    c_cartesian.solve(front->dt);
	    l_cartesian->solve(front->dt);
	    c_cartesian.timeStepAnalysis(NO);

	    FT_SetTimeStep(front);
	    c_cartesian.setAdvectionDt();
	    front->dt = std::min(front->dt,CFL*c_cartesian.max_dt);
	    l_cartesian->setAdvectionDt();
	    front->dt = std::min(front->dt,CFL*l_cartesian->max_dt);

	    FT_Save(front);
            c_cartesian.printFrontInteriorStates();
            l_cartesian->printFrontInteriorStates(out_name);

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
	    FT_Propagate(front);
	    c_cartesian.solve(front->dt);
	    l_cartesian->solve(front->dt);
	    FT_AddTimeStepToCounter(front);

	    FT_SetTimeStep(front);
	    front->dt = std::min(front->dt,CFL*c_cartesian.max_dt);
	    front->dt = std::min(front->dt,CFL*l_cartesian->max_dt);

            if (FT_IsSaveTime(front))
	    {
                FT_Save(front);
		c_cartesian.printFrontInteriorStates();
		l_cartesian->printFrontInteriorStates(out_name);
	    }
            if (FT_IsDrawTime(front))
	    {
                FT_Draw(front);
		c_cartesian.crystalDraw();
	    }

            if (FT_TimeLimitReached(front) || bdryReached(front))
	    {
		c_cartesian.timeStepAnalysis(YES);
		FT_PrintTimeStamp(front);
		FT_Draw(front);
		FT_Save(front);
		c_cartesian.printFrontInteriorStates();
		l_cartesian->printFrontInteriorStates(out_name);
                break;
	    }

	    if (FT_IsDrawTime(front))
                c_cartesian.timeStepAnalysis(YES);
	    else
	    	c_cartesian.timeStepAnalysis(NO);
	    FT_TimeControlFilter(front);
	    FT_PrintTimeStamp(front);
        }
}       /* end subsurf_main_driver */

static  void solute_point_propagate(
        Front *front,
        POINTER wave,
        POINT *oldp,
        POINT *newp,
        HYPER_SURF_ELEMENT *oldhse,
        HYPER_SURF         *oldhs,
        double              dt,
        double              *V)
{
	ifluid_point_propagate(front,wave,oldp,newp,oldhse,oldhs,dt,V);
	crystal_point_propagate(front,wave,oldp,newp,oldhse,oldhs,dt,V);
}	/* end solute_point_propagate */

static void print_area_density(
        Front *front,
        char *out_name)
{
        INTERFACE *intfc = front->interf;
        CURVE **c, *curve;
        double total_length, vol_frac, area_density;
        char fname[100];
        static FILE *lfile;
        static FILE *vfile;
        static FILE *afile;

        total_length = 0.0;
        vol_frac = 0.0;
        area_density = 0.0;

        if (lfile == NULL)
	{
            sprintf(fname,"%s/total_length",out_name);
            lfile = fopen(fname,"w");
            fprintf(lfile,"\"Total length vs. time\"\n");
        }

        if (vfile == NULL)
        {
            sprintf(fname,"%s/vol_frac",out_name);
            vfile = fopen(fname,"w");
            fprintf(vfile,"\"Volume fraction vs. time\"\n");
        }

        if (afile == NULL)
        {
            sprintf(fname,"%s/area_density",out_name);
            afile = fopen(fname,"w");
            fprintf(afile,"\"Area density vs. time\"\n");
        }

        curve = NULL;
        for (c = intfc->curves; c && *c; ++c)
        {
            if (wave_type(*c) == GROWING_BODY_BOUNDARY)
            {
                curve = *c;
                total_length += curve_length(curve);
            }
        }

        vol_frac = FT_ComputeTotalVolumeFraction(front,CRYSTAL_COMP);

        area_density = total_length/vol_frac;

        fprintf(lfile,"%f %f\n",front->time,total_length);
        fprintf(vfile,"%f %f\n",front->time,vol_frac);
        fprintf(afile,"%f %f\n",front->time,area_density);
        fflush(lfile);
        fflush(vfile);
        fflush(afile);
}       /* end print_area_density */

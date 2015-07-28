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

#include "crystal.h"

#define		MAX_NUM_VERTEX_IN_CELL		20

	/*  Local Application Function Declarations */

static void 	crystal_driver(Front*,C_CARTESIAN&);

char *in_name,*restart_state_name,*restart_name;
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
	
	/*Initialize Petsc */
	PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);

	in_name      		= f_basic.in_name;
	restart_state_name      = f_basic.restart_state_name;
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
	if (debugging("init"))
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
	FT_ReadTimeControl(in_name,&front);

	/* Initialize velocity field function */

	
	velo_func_pack.func_params = (POINTER)&cRparams;
	velo_func_pack.func = NULL;

	FT_InitFrontVeloFunc(&front,&velo_func_pack);

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
	cRparams.field->vel = NULL; 	// No convection

	/* Propagate the front */
	crystal_driver(&front,c_cartesian);

	PetscFinalize();
	clean_up(0);
}

static  void crystal_driver(
        Front *front,
	C_CARTESIAN &c_cartesian)
{
        double CFL;

	Curve_redistribution_function(front) = expansion_redistribute;
        CFL = Time_step_factor(front);

        c_cartesian.crystalDraw();

        if (!RestartRun)
        {
	    FT_ResetTime(front);
	    FT_Propagate(front);
	    c_cartesian.solve(front->dt);
	    c_cartesian.timeStepAnalysis(NO);
	    FT_SetTimeStep(front);
	    c_cartesian.setAdvectionDt();
	    front->dt = std::min(front->dt,CFL*c_cartesian.max_dt);
            FT_Save(front);
	    c_cartesian.printFrontInteriorStates();
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
	    FT_AddTimeStepToCounter(front);

	    FT_SetTimeStep(front);
	    front->dt = std::min(front->dt,CFL*c_cartesian.max_dt);

            if (FT_IsSaveTime(front))
	    {
                FT_Save(front);
		c_cartesian.printFrontInteriorStates();
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
                break;
	    }
            if (FT_IsDrawTime(front))
	    	c_cartesian.timeStepAnalysis(YES);
	    else
	    	c_cartesian.timeStepAnalysis(NO);
	    FT_TimeControlFilter(front);
	    FT_PrintTimeStamp(front);
        }
}       /* end crystal_driver */

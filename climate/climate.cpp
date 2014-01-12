/*********************************************************************
FronTier is a set of libraries that implements differnt types of Front 
Traking algorithms. Front Tracking is a numerical method for the solution 
of partial differential equations whose solutions have discontinuities.  


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

**********************************************************************/


/*
*				climate.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*/

#include "../iFluid/iFluid.h"
#include "../iFluid/ifluid_basic.h"
#include "climate.h"

#define		MAX_NUM_VERTEX_IN_CELL		20
	/*  Local Application Function Declarations */

static void	melting_flow_driver(Front*,CARTESIAN*,
			VCARTESIAN*, Incompress_Solver_Smooth_Basis *);
static int      rgbody_vel(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,
                        HYPER_SURF*,double*);
static double 	temperature_func(double*,COMPONENT,double);
static void 	read_movie_options(char*,PARAMS*);
static void	melt_flow_point_propagate(Front*,POINTER,POINT*,POINT*,
			HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);

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
	static PARAMS eqn_params;
	static IF_PARAMS iFparams;
	int dim;

	FT_Init(argc,argv,&f_basic);

        PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);

	CARTESIAN *t_cartesian = new CARTESIAN(front);
	VCARTESIAN *v_cartesian = new VCARTESIAN(front);
	Incompress_Solver_Smooth_Basis *l_cartesian = NULL;
	if (f_basic.dim == 2)
	    l_cartesian = new Incompress_Solver_Smooth_2D_Cartesian(front);
	else if (f_basic.dim == 3)
	    l_cartesian = new Incompress_Solver_Smooth_3D_Cartesian(front);

	/* Initialize basic computational data */
	f_basic.size_of_intfc_state = sizeof(STATE);

	in_name      		= f_basic.in_name;
	restart_state_name      = f_basic.restart_state_name;
        out_name     		= f_basic.out_name;
        restart_name 		= f_basic.restart_name;
        RestartRun   		= f_basic.RestartRun;
        ReadFromInput   	= f_basic.ReadFromInput;
	RestartStep 		= f_basic.RestartStep;
	dim	 		= f_basic.dim;

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
        FT_ReadSpaceDomain(in_name,&f_basic);

	FT_StartUp(&front,&f_basic);
	FT_InitDebug(in_name);

	eqn_params.dim = f_basic.dim;
	iFparams.dim = f_basic.dim;
	read_movie_options(in_name,&eqn_params);
	front.extra1 = (POINTER)&iFparams;
	front.extra2 = (POINTER)&eqn_params;
	read_CL_prob_type(&front);
	readPhaseParams(&front);
        read_iFparams(in_name,&iFparams);
	if (!RestartRun)
	{
	    if(eqn_params.no_droplets == NO)
	    {
		printf("Initializing droplets\n");
		level_func_pack.pos_component = LIQUID_COMP2;
	    	FT_InitIntfc(&front,&level_func_pack);
                initWaterDrops(&front);
	    	if (debugging("trace")) printf("Passed init water droplets()\n");
	    }
	    else
	    {
	        printf("No droplets contained\n");
	        level_func_pack.func_params = NULL;
                level_func_pack.func = NULL;
	        level_func_pack.pos_component = LIQUID_COMP2;
	        level_func_pack.wave_type = -1; 
	        FT_InitIntfc(&front,&level_func_pack);
	    }
	    read_melt_dirichlet_bdry_data(in_name,&front,f_basic);
	    if (f_basic.dim != 3)
                FT_ClipIntfcToSubdomain(&front);
	    if (debugging("trace"))
                printf("Passed FT_ClipIntfcToSubdomain()\n");
	}
	else
	{
	    readWaterDropsParams(&front,restart_state_name);
	    read_melt_dirichlet_bdry_data(in_name,&front,f_basic);
	}

	front._compute_force_and_torque = compute_ice_particle_force;
	FT_ReadTimeControl(in_name,&front);

	if(eqn_params.prob_type == PARTICLE_TRACKING)
	    gv_plot_scatter(&front);
	/* Initialize velocity field function */

	velo_func_pack.func_params = (POINTER)&iFparams;
	velo_func_pack.func = rgbody_vel;
	velo_func_pack.point_propagate = melt_flow_point_propagate;

	FT_InitVeloFunc(&front,&velo_func_pack);


        t_cartesian->initMesh();
	v_cartesian->field = t_cartesian->field;
        v_cartesian->initMesh();

	l_cartesian->initMesh();
	l_cartesian->findStateAtCrossing = ifluid_find_state_at_crossing;
	if (RestartRun)
	{
	    v_cartesian->readFrontInteriorState(restart_state_name);
	    t_cartesian->readFrontInteriorState(restart_state_name);
	    FT_FreeGridIntfc(&front);
	    l_cartesian->readFrontInteriorStates(restart_state_name);
	}
	else
	{
	    t_cartesian->setInitialCondition();
	    v_cartesian->setInitialCondition();
	    if (debugging("trace")) 
                printf("Passed cartesian setInitialCondition()\n");
	    FT_FreeGridIntfc(&front);

            init_fluid_state_func(&front,l_cartesian);
            l_cartesian->setInitialCondition();
            if (debugging("trace"))
                printf("Passed iFluid setInitialCondition()\n");
	}

	eqn_params.field->vel = iFparams.field->vel;
	eqn_params.field->pres = iFparams.field->pres;
	iFparams.field->temperature = eqn_params.field->temperature;
	printf("After reading data:\n");
	v_cartesian->checkField(); /*output variables to run-output*/

	FT_InitVeloFunc(&front,&velo_func_pack);

	if (debugging("trace")) printf("Passed FT_InitVeloFunc()\n");

	FT_SetGlobalIndex(&front);
	
	/* Propagate the front */
	melting_flow_driver(&front,t_cartesian, v_cartesian, l_cartesian);

	PetscFinalize();
	clean_up(0);
}

static  void melting_flow_driver(
        Front *front,
	CARTESIAN *t_cartesian,
	VCARTESIAN *v_cartesian,
	Incompress_Solver_Smooth_Basis *l_cartesian)
{
        double CFL;
        int  dim = front->rect_grid->dim;
	IF_PARAMS *iFparams;
	PARAMS *eqn_params;
        static LEVEL_FUNC_PACK level_func_pack;

	if (debugging("trace"))
	    printf("Entering melting_flow_driver()\n");
	Curve_redistribution_function(front) = expansion_redistribute;
        CFL = Time_step_factor(front);

	iFparams = (IF_PARAMS*)front->extra1;
	eqn_params = (PARAMS*)front->extra2;

	front->hdf_movie_var = NULL;

        if (!RestartRun)
        {
	    FT_ResetTime(front);
            FT_SetOutputCounter(front);
            /* Front standard output*/
            FT_Save(front,out_name);
	    if (dim == 1)
	    	t_cartesian->oneDimPlot(out_name);
	    else if (dim == 2)
		t_cartesian->vtk_plot_temperature2d(out_name);
	    else
	    {
		t_cartesian->vtk_plot3d(out_name,"temperature");
		v_cartesian->vtk_plot3d(out_name,"vapor");
		v_cartesian->vtk_plot3d(out_name,"supersat");
		l_cartesian->vtk_plot_scalar(out_name,"pres");
                if (eqn_params->prob_type == PARTICLE_TRACKING)
                    vtk_plot_scatter(front);
	    }


	    FT_ResetTime(front);
	    FT_Save(front,out_name);
            t_cartesian->printFrontInteriorState(out_name);
            v_cartesian->printFrontInteriorState(out_name);
            l_cartesian->printFrontInteriorStates(out_name);
	    printDropletsStates(front,out_name);

	    if(eqn_params->droplets_fixed == NO)
	    {
	        CondensationPreAdvance(front);
	        FrontPreAdvance(front);
	    }

	    FT_Propagate(front);
            if(eqn_params->prob_type == PARTICLE_TRACKING)
		ParticlePropagate(front);

	    v_cartesian->solve(front->dt);
	    printf("passing solving vapor\n\n");
	    t_cartesian->solve(front->dt);
	    printf("passing solving temperature\n\n");

	    l_cartesian->solve(front->dt);

	    FT_SetTimeStep(front);
	    t_cartesian->setAdvectionDt();
	    v_cartesian->setAdvectionDt();
	    l_cartesian->setAdvectionDt();
	    front->dt = std::min(front->dt,CFL*t_cartesian->m_dt);
	    front->dt = std::min(front->dt,CFL*v_cartesian->m_dt);
	    front->dt = std::min(front->dt,CFL*l_cartesian->max_dt);

            l_cartesian->initMovieVariables();
            t_cartesian->initMovieVariables();
            v_cartesian->initMovieVariables();
	    if (eqn_params->prob_type == PARTICLE_TRACKING)
		vtk_plot_scatter(front);
	    FT_AddMovieFrame(front,out_name,YES);
        }
        else
	{
	    FT_SetOutputCounter(front);
            t_cartesian->initMovieVariables();
            v_cartesian->initMovieVariables();
            if (eqn_params->prob_type == PARTICLE_TRACKING)
                vtk_plot_scatter(front);
            FT_AddMovieFrame(front,out_name,YES);
	}

	FT_TimeControlFilter(front);

        for (;;)
        {
	    if(eqn_params->droplets_fixed == NO)
	    {
	        CondensationPreAdvance(front);
	        FrontPreAdvance(front);
	    }
	    FT_Propagate(front);
            if(eqn_params->prob_type == PARTICLE_TRACKING)
		ParticlePropagate(front);

	    t_cartesian->solve(front->dt);
	    v_cartesian->solve(front->dt);
	    printf("Passed solving vapor and temperature equations\n");

	    l_cartesian->solve(front->dt);
	    printf("Passed solving NS equations\n");
	    
	    FT_AddTimeStepToCounter(front);
	    FT_SetTimeStep(front);
	    front->dt = FT_Min(front->dt,CFL*t_cartesian->m_dt);
	    front->dt = FT_Min(front->dt,CFL*v_cartesian->m_dt);
	    front->dt = FT_Min(front->dt,CFL*l_cartesian->max_dt);

            printf("\ntime = %10.9f   step = %7d   dt = %10.9f\n",
                        front->time,front->step,front->dt);
            fflush(stdout);

	    /*For checking the result*/
	    v_cartesian->checkField(); /*output variables to run-output*/
	    v_cartesian->recordField(out_name,"supersat");
	    v_cartesian->recordField(out_name,"vapor");
	    v_cartesian->recordField(out_name,"temperature");
	    v_cartesian->recordMomentum();
            if (FT_IsSaveTime(front))
	    {
		// Front standard output
                FT_Save(front,out_name);

		v_cartesian->recordRadius(out_name);
		// Problem specific output
		t_cartesian->printFrontInteriorState(out_name);
		v_cartesian->printFrontInteriorState(out_name);
		l_cartesian->printFrontInteriorStates(out_name);
	        printDropletsStates(front,out_name);
	    }
            if (FT_IsMovieFrameTime(front))
	    {
		// Front standard output
		if(eqn_params->prob_type == PARTICLE_TRACKING)
		    vtk_plot_scatter(front);
                FT_AddMovieFrame(front,out_name,YES);

		if (dim == 1)
	    	    t_cartesian->oneDimPlot(out_name);
		else if (dim == 2)
		    t_cartesian->vtk_plot_temperature2d(out_name);
		else
		{
		        t_cartesian->vtk_plot3d(out_name,"temperature");
		        v_cartesian->vtk_plot3d(out_name,"vapor");
		        v_cartesian->vtk_plot3d(out_name,"supersat");
		        l_cartesian->vtk_plot_scalar(out_name,"pres");
		}
	    }

            if (FT_TimeLimitReached(front))
	    {
		if(eqn_params->prob_type == PARTICLE_TRACKING)
                    vtk_plot_scatter(front);
	    	FT_AddMovieFrame(front,out_name,YES);
		v_cartesian->plotMomentum(out_name);
		v_cartesian->recordRadius(out_name);
                break;
	    }
	    /* Output section, next dt may be modified */

	    FT_TimeControlFilter(front);
        }
}       /* end melting_flow_driver */

static void melt_flow_point_propagate(
	Front *front,
        POINTER wave,
        POINT *oldp,
        POINT *newp,
        HYPER_SURF_ELEMENT *oldhse,
        HYPER_SURF         *oldhs,
        double              dt,
        double              *V)
{
	PARAMS *eqn_params = (PARAMS*)front->extra2;

        ifluid_point_propagate(front,wave,oldp,oldp,oldhse,oldhs,dt,V);
        melting_point_propagate(front,wave,oldp,newp,oldhse,oldhs,dt,V);
}       /* end melt_flow_point_propagate */

static int rgbody_vel(
        POINTER params,
        Front *front,
        POINT *p,
        HYPER_SURF_ELEMENT *hse,
        HYPER_SURF *hs,
        double *vel)
{
        int i,dim;
        dim  = front->rect_grid->dim;
        for (i = 0; i < dim; ++i)
        {
            vel[i] = center_of_mass_velo(hs)[i];
        }
        return YES;
}       /* end rgbody_vel */

static	double temperature_func(
	double *coords,
	COMPONENT comp,
	double T_l)
{
	if (comp != LIQUID_COMP2) return 0.0;
	return T_l;
}	/* end temperature_init_func */

static void read_movie_options(
        char *inname,
        PARAMS *params)
{
        static MOVIE_OPTION *movie_option;
        FILE *infile = fopen(inname,"r");
        char string[100];

        FT_ScalarMemoryAlloc((POINTER*)&movie_option,sizeof(MOVIE_OPTION));
        params->movie_option = movie_option;

        if (params->dim == 1) 
	{
	    movie_option->plot_temperature = YES;
	    return;
	}

        CursorAfterString(infile,"Type y to make movie of temperature:");
        fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
        if (string[0] == 'Y' || string[0] == 'y')
            movie_option->plot_temperature = YES;
        CursorAfterString(infile,"Type y to make movie of vapor mixing ratio:");
        fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
        if (string[0] == 'Y' || string[0] == 'y')
            movie_option->plot_vapor = YES;

	/* Default: not plot cross sectional variables */
        movie_option->plot_cross_section[0] = NO;
        movie_option->plot_cross_section[1] = NO;
        movie_option->plot_cross_section[2] = NO;
        fclose(infile);
}       /* end read_movie_options */

extern  double getStateTemperature(
        POINTER state)
{
        STATE *T_state = (STATE*)state;
        return T_state->temperature;
}       /* end getStateTemperature */

extern  double getStateVapor(
        POINTER state)
{
        STATE *T_state = (STATE*)state;
        return T_state->vapor;
}       /* end getStateVapor */

extern  double getStateSuper(
        POINTER state)
{
        STATE *T_state = (STATE*)state;
        return T_state->supersat;
}       /* end getStateSuper */

extern  void assignStateTemperature(
	double T,
        POINTER state)
{
        STATE *T_state = (STATE*)state;
        T_state->temperature = T;
}       /* end assignStateTemperature */

extern  void assignStateVapor(
	double T,
        POINTER state)
{
        STATE *T_state = (STATE*)state;
        T_state->vapor = T;
}       /* end assignStateVapor */

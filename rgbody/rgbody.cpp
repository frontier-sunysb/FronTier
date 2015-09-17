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

#include <iFluid.h>
#include <rgbody.h>

static void fluid_driver(Front*,Incompress_Solver_Smooth_Basis*);
static int rgbody_vel(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,
                        HYPER_SURF*,double*);

static void prompt_and_make_rigid_body(Front*,FILE*,POINTER,
		double(*func)(POINTER,double*),COMPONENT,COMPONENT,int);

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
	IF_PROB_TYPE prob_type;
	RG_PARAMS *rg_params;

	/* Initialize basic computational data */

	FT_Init(argc,argv,&f_basic);
	f_basic.size_of_intfc_state = sizeof(STATE);

	//Initialize Petsc
        PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);

	/*Construct Incompress Solver cartesian*/

	Incompress_Solver_Smooth_Basis *cartesian = NULL;
	if(f_basic.dim == 2)
	    cartesian = new Incompress_Solver_Smooth_2D_Cartesian(front);
	else if(f_basic.dim == 3)
	    cartesian = new Incompress_Solver_Smooth_3D_Cartesian(front);

        in_name                 = f_basic.in_name;
        restart_state_name      = f_basic.restart_state_name;
        out_name                = f_basic.out_name;
        restart_name            = f_basic.restart_name;
        RestartRun              = f_basic.RestartRun;
        RestartStep             = f_basic.RestartStep;

        sprintf(restart_state_name,"%s/state.ts%s",restart_name,
                        right_flush(RestartStep,7));
        sprintf(restart_name,"%s/intfc-ts%s",restart_name,
			right_flush(RestartStep,7));
	if (pp_numnodes() > 1)
	{
            sprintf(restart_name,"%s-nd%s",restart_name,
				right_flush(pp_mynode(),4));
            sprintf(restart_state_name,"%s-nd%s",restart_state_name,
                        	right_flush(pp_mynode(),4));
	}

	FT_ReadSpaceDomain(in_name,&f_basic);
	FT_StartUp(&front,&f_basic);
	FT_InitDebug(in_name);

	if (debugging("trace")) printf("Passed FT_StartUp()\n");
	iFparams.dim = f_basic.dim;
	front.extra1 = (POINTER)&iFparams;

	read_rg_prob_type(in_name,&prob_type);
	read_iFparams(in_name,&iFparams);

	if (!RestartRun)
	{
	    level_func_pack.neg_component = SOLID_COMP;
	    level_func_pack.pos_component = LIQUID_COMP2;
	    level_func_pack.func = NULL;
	    level_func_pack.point_array = NULL;
	    level_func_pack.func_params = iFparams.level_func_params;
	    if (f_basic.dim == 3)
	    	level_func_pack.set_3d_bdry = YES;

	    init_moving_bodies(&front,&level_func_pack,in_name,prob_type);
	    FT_InitIntfc(&front,&level_func_pack);
	    if (debugging("trace"))
	    {
		char test_name[100];
	    	printf("Passed FT_InitIntfc()\n");
		switch (f_basic.dim)
                {
                case 2:
                    sprintf(test_name,"init_intfc-%d.xg",pp_mynode());
                    xgraph_2d_intfc(test_name,front.interf);
                    break;
                case 3:
                    sprintf(test_name,"init_intfc-%d.xg",pp_mynode());
                    gview_plot_interface("gv-init",front.interf);
                    break;
                }
	    }
	    FT_PromptSetMixedTypeBoundary2d(in_name,&front);
	    read_iF_dirichlet_bdry_data(in_name,&front,f_basic);
	    if (f_basic.dim < 3)
	    	FT_ClipIntfcToSubdomain(&front);
	}
	else
	    restart_set_dirichlet_bdry_function(&front);

	front._compute_force_and_torque = ifluid_compute_force_and_torque;

	/* Initialize velocity field function */

	velo_func_pack.func_params = (POINTER)cartesian;
	velo_func_pack.func = rgbody_vel;
	velo_func_pack.point_propagate = ifluid_point_propagate;
	FT_InitFrontVeloFunc(&front,&velo_func_pack);
	if (debugging("trace"))
	    printf("Passed FT_InitFrontVeloFunc()\n");

	cartesian->initMesh();
	cartesian->initMovieVariables();
	if (debugging("sample_velocity"))
            cartesian->initSampleVelocity(in_name);
	init_fluid_state_func(cartesian,prob_type);
	if (debugging("trace"))
	    printf("Passed cartesian.initMesh()\n");
	cartesian->findStateAtCrossing = ifluid_find_state_at_crossing;

	if (RestartRun)
	{
	    cartesian->readFrontInteriorStates(restart_state_name);
	    FT_ScalarMemoryAlloc((POINTER*)&rg_params, sizeof(RG_PARAMS));
	    prompt_for_rigid_body_params(f_basic.dim,in_name,rg_params);
	    front.extra3 = (POINTER)rg_params;
	}
	else
	    cartesian->setInitialCondition();

	/* Propagate the front */

	fluid_driver(&front,cartesian);

	PetscFinalize();
	clean_up(0);
}

static  void fluid_driver(
        Front *front,
	Incompress_Solver_Smooth_Basis *cartesian)
{
        double CFL;
	RG_PARAMS *rg_params = (RG_PARAMS*)front->extra3;

	Curve_redistribution_function(front) = full_redistribute;

	FT_ReadTimeControl(in_name,front);
	CFL = Time_step_factor(front);

	if (!RestartRun)
	{
	    FT_RedistMesh(front);
	    FT_ResetTime(front);

	    if (debugging("trace"))
		printf("Calling initial FT_Propagate()\n");
	    FrontPreAdvance(front);
            FT_Propagate(front);
	    if (!rg_params->no_fluid)
	    {
	    	if (debugging("trace")) printf("Begin calling solve()\n");
	    	cartesian->solve(front->dt);
	    	if (debugging("trace")) printf("Passed solve()\n");
	    }
	    record_moving_body_data(out_name,front);
            FT_Draw(front);
	    FT_SetOutputCounter(front);
	    FT_SetTimeStep(front);
            front->dt = std::min(front->dt,CFL*cartesian->max_dt);
        }
        else
	    FT_SetOutputCounter(front);

	FT_TimeControlFilter(front);
	FT_PrintTimeStamp(front);

	if (debugging("trace"))
	{
	    printf("CFL = %f\n",CFL);
	    printf("Frequency_of_redistribution(front,GENERAL_WAVE) = %d\n",
			Frequency_of_redistribution(front,GENERAL_WAVE));
	}

	if (debugging("trace")) printf("Before time loop\n");
        for (;;)
        {
            /* Propagating interface for time step dt */

	    if (debugging("trace")) printf("Begin a time step\n");
	    FrontPreAdvance(front);
	    if (debugging("trace")) printf("Passed FrontPreAdvance()\n");
            FT_Propagate(front);

	    if (!rg_params->no_fluid)
	    {
	    	if (debugging("trace")) printf("Begin calling solve()\n");
	    	cartesian->solve(front->dt);
	    	if (debugging("trace")) printf("Passed solve()\n");
	    }

	    FT_AddTimeStepToCounter(front);
				
            //Next time step determined by maximum speed of previous
            //step, assuming the propagation is hyperbolic and
            //is not dependent on second order derivatives of
            //the interface such as curvature, and etc.

	    FT_SetTimeStep(front);
	    if (debugging("step_size"))
		printf("Time step from FrontHypTimeStep(): %f\n",front->dt);
            front->dt = std::min(front->dt,CFL*cartesian->max_dt);
	    if (debugging("step_size"))
                printf("Time step from l_cartesian->max_dt(): %f\n",front->dt);
	
            /* Output section */

	    record_moving_body_data(out_name,front);

            if (FT_IsSaveTime(front))
	    {
            	FT_Save(front);
		cartesian->printFrontInteriorStates(out_name);
	    }
            if (FT_IsDrawTime(front))
	    {
            	FT_Draw(front);
	    }

            if (FT_TimeLimitReached(front))
	    {
		FT_PrintTimeStamp(front);
                break;
	    }

	    FT_TimeControlFilter(front);
	    FT_PrintTimeStamp(front);
        }
	if (debugging("trace")) printf("After time loop\n");
}       /* end fluid_driver */


static int rgbody_vel(
	POINTER params,
	Front *front,
	POINT *p,
	HYPER_SURF_ELEMENT *hse,
	HYPER_SURF *hs,
	double *vel)
{
	int i,dim; 
	double omega,crds_at_com[MAXD];

	dim  = front->rect_grid->dim;
	for (i = 0; i < dim; ++i)
	{
	    crds_at_com[i] = Coords(p)[i] - rotation_center(hs)[i];
            vel[i] = center_of_mass_velo(hs)[i];
	}
	switch(dim)
	{
	case 2:
	    omega = angular_velo(hs);
	    vel[0] += omega*crds_at_com[1];
	    vel[1] += -omega*crds_at_com[0];
	    break;
	case 3:
//	    omega = angular_velo(hs);
//	    vel[0] += omega*crds_at_com[2];
//	    vel[2] += -omega*crds_at_com[0];
	    if (motion_type(hs) == ROTATION ||
		motion_type(hs) == PRESET_ROTATION ||
		motion_type(hs) == FREE_MOTION)
	    {
		vel[0] += -p_angular_velo(hs)[2] * crds_at_com[1]
                          +p_angular_velo(hs)[1] * crds_at_com[2];
                vel[1] +=  p_angular_velo(hs)[2] * crds_at_com[0]
                          -p_angular_velo(hs)[0] * crds_at_com[2];
                vel[2] += -p_angular_velo(hs)[1] * crds_at_com[0]
                          +p_angular_velo(hs)[0] * crds_at_com[1];
	    }
	}
	return YES;
}	/* end rgbody_vel */

extern void record_moving_body_data(
	char *out_name,
	Front *front)
{
	int i,j,num_moving_body;
        INTERFACE *intfc = front->interf;
        CURVE **c;
        SURFACE **s;
	static boolean first = YES;
	static FILE **torque_files,**omega_files;
	static FILE **force_files,**com_files;
	static FILE **pa_velo_files,**euler_files;
	static double **torque,*omega;
	static double **force,**com_velo;
	static double **pa_velo,**euler;
	int dim = intfc->dim;
	char fname[256];

	num_moving_body = -1;
	switch (dim)
	{
	case 2:
            for (c = intfc->curves; c && *c; ++c)
            {
            	if (wave_type(*c) == MOVABLE_BODY_BOUNDARY)
                    if (body_index(*c) > num_moving_body)
                    	num_moving_body = body_index(*c);
            }
	    break;
	case 3:
            for (s = intfc->surfaces; s && *s; ++s)
            {
            	if (wave_type(*s) == MOVABLE_BODY_BOUNDARY)
                    if (body_index(*s) > num_moving_body)
                    	num_moving_body = body_index(*s);
            }
	    break;
	}
        num_moving_body++;
        pp_global_imax(&num_moving_body,1);
        if (num_moving_body == 0) return;

        if (first)
        {
            FT_VectorMemoryAlloc((POINTER*)&omega,num_moving_body,FLOAT);
            FT_MatrixMemoryAlloc((POINTER*)&torque,num_moving_body,MAXD,FLOAT);
            FT_MatrixMemoryAlloc((POINTER*)&force,num_moving_body,MAXD,FLOAT);
            FT_MatrixMemoryAlloc((POINTER*)&com_velo,num_moving_body,
					MAXD,FLOAT);
	    FT_MatrixMemoryAlloc((POINTER*)&pa_velo,num_moving_body,
                                        MAXD,FLOAT);
            FT_MatrixMemoryAlloc((POINTER*)&euler,num_moving_body,
                                        MAXD+1,FLOAT);
	    if (pp_mynode() == 0)
	    {
            	FT_VectorMemoryAlloc((POINTER*)&torque_files,num_moving_body,
					sizeof(FILE*));
            	FT_VectorMemoryAlloc((POINTER*)&omega_files,num_moving_body,
					sizeof(FILE*));
            	FT_VectorMemoryAlloc((POINTER*)&force_files,num_moving_body,
					sizeof(FILE*));
            	FT_VectorMemoryAlloc((POINTER*)&com_files,num_moving_body,
					sizeof(FILE*));
		FT_VectorMemoryAlloc((POINTER*)&pa_velo_files,
                                        num_moving_body,sizeof(FILE*));
                FT_VectorMemoryAlloc((POINTER*)&euler_files,
                                        num_moving_body,sizeof(FILE*));
            	for (i = 0; i < num_moving_body; ++i)
		{
		    sprintf(fname,"%s/cen-of_mass-%d",out_name,i);
		    com_files[i] = fopen(fname,"w");
		    sprintf(fname,"%s/omega-%d",out_name,i);
		    omega_files[i] = fopen(fname,"w");
		    sprintf(fname,"%s/force-%d",out_name,i);
		    force_files[i] = fopen(fname,"w");
		    sprintf(fname,"%s/torque-%d",out_name,i);
		    torque_files[i] = fopen(fname,"w");
		    sprintf(fname,"%s/pa_velo-%d",out_name,i);
                    pa_velo_files[i] = fopen(fname,"w");
                    sprintf(fname,"%s/euler-%d",out_name,i);
                    euler_files[i] = fopen(fname,"w");
		}
	    }
        }
	switch (dim)
	{
	case 2:
	    for (c = intfc->curves; c && *c; ++c)
            {
            	if (wave_type(*c) == MOVABLE_BODY_BOUNDARY)
            	{
                    i = body_index(*c);
                    FrontForceAndTorqueOnHs(front,Hyper_surf(*c),front->dt,
					force[i],torque[i]);
                    omega[i] = angular_velo(*c);
                    for (j = 0; j < dim; ++j)
                    	com_velo[i][j] = center_of_mass_velo(*c)[j];
            	}
            }
	    break;
	case 3:
	    for (s = intfc->surfaces; s && *s; ++s)
            {
            	if (wave_type(*s) == MOVABLE_BODY_BOUNDARY)
            	{
                    i = body_index(*s);
                    FrontForceAndTorqueOnHs(front,Hyper_surf(*s),front->dt,
					force[i],torque[i]);
                    omega[i] = angular_velo(*s);
                    for (j = 0; j < dim; ++j)
                    	com_velo[i][j] = center_of_mass_velo(*s)[j];
		    for (j = 0; j < dim; ++j)
                        pa_velo[i][j] = p_angular_velo(*s)[j];
                    for (j = 0; j < 4; ++j)
                        euler[i][j] = euler_params(*s)[j];
            	}
            }
	    break;
	}
        for (i = 0; i < num_moving_body; ++i)
        {
            pp_global_sum(force[i],dim);
            pp_global_sum(torque[i],dim);
        }
        if (pp_mynode() != 0) return;

	if (first)
        {
            first = NO;
            for (i = 0; i < num_moving_body; ++i)
            {
                fprintf(torque_files[i],"\"Torque of body %d\"\n",i+1);
                fprintf(force_files[i],"\"Total force on body %d\"\n",i+1);
                fprintf(omega_files[i],"\"Angular velocity of body %d\"\n",i+1);                fprintf(com_files[i],"\"COM velocity of body %d\"\n",i+1);
            }
        }
        for (i = 0; i < num_moving_body; ++i)
        {
            fprintf(omega_files[i],"%f  %f\n",front->time,omega[i]);
            fprintf(torque_files[i],"%f  ",front->time);
            fprintf(force_files[i],"%f  ",front->time);
            fprintf(com_files[i],"%f  ",front->time);
	    fprintf(pa_velo_files[i],"%f \t",front->time);
            fprintf(euler_files[i],"%f \t",front->time);
            for (j = 0; j < dim; ++j)
            {
            	fprintf(torque_files[i],"%f  ",torque[i][j]);
                fprintf(force_files[i],"%f  ",force[i][j]);
                fprintf(com_files[i],"%f  ",com_velo[i][j]);
		fprintf(pa_velo_files[i],"%f \t",pa_velo[i][j]);
            }
            for (j = 0; j < 4; ++j)
                fprintf(euler_files[i],"%f \t",euler[i][j]);
            fprintf(torque_files[i],"\n");
            fprintf(force_files[i],"\n");
            fprintf(com_files[i],"\n");
	    fprintf(pa_velo_files[i],"\n");
            fprintf(euler_files[i],"\n");

            fflush(torque_files[i]);
            fflush(omega_files[i]);
            fflush(force_files[i]);
            fflush(com_files[i]);
	    fflush(pa_velo_files[i]);
            fflush(euler_files[i]);
        }
}	/* end record_moving_body_data */

extern void read_rg_prob_type(
	char *inname,
	IF_PROB_TYPE *prob_type)
{
	char string[100];
	FILE *infile = fopen(inname,"r");

	*prob_type = ERROR_TYPE;
	(void) printf("Available problem types are:\n");
        (void) printf("\tFLUID_SOLID_CIRCLE\n");
        (void) printf("\tFLUID_RIGID_BODY\n");
        (void) printf("\tOPEN_ROTOR\n");
        (void) printf("\tPRESSURE_PUMP\n");
        (void) printf("\tROTOR_ONE_FLUID\n");
        (void) printf("\tROTOR_TWO_FLUID\n");
        (void) printf("\tWINDMILL_2D\n");
        (void) printf("\tWINDMILL_3D\n");
        (void) printf("\tBEE_3D\n");
        (void) printf("\tHELICOPTER_3D\n");
        (void) printf("\tFLUID_SOLID_CONE\n");
	CursorAfterString(infile,"Enter problem type:");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	if (string[0] == 'F' || string[0] == 'f')
	{
            if (string[6] == 'S' || string[6] == 's')
	    {
		if (string[13] == 'i' || string[13] == 'I')
                    *prob_type = FLUID_SOLID_CIRCLE;
                else if (string[13] == 'o' || string[13] == 'O')
                    *prob_type = FLUID_SOLID_CONE;
            }
            else if (string[6] == 'R' || string[6] == 'r')
                *prob_type = FLUID_RIGID_BODY;
	}
        else if (string[0] == 'R' || string[0] == 'r')
	{
            if (string[6] == 'O' || string[6] == 'o')
                *prob_type = ROTOR_ONE_FLUID;
            else if (string[6] == 'T' || string[6] == 't')
                *prob_type = ROTOR_TWO_FLUID;
	}
        else if (string[0] == 'W' || string[0] == 'w')
	{
            if (string[9] == '2')
                *prob_type = WINDMILL_2D;
            else if (string[9] == '3')
                *prob_type = WINDMILL_3D;
	}
        else if (string[0] == 'B' || string[0] == 'b')
	{
            *prob_type = BEE_3D;
	}
        else if (string[0] == 'H' || string[0] == 'h')
	{
	    if (string[1] == 'E' || string[1] == 'e')
                *prob_type = HELICOPTER_3D;
            else if (string[1] == 'U' || string[1] == 'u')
                *prob_type = HUMAN_BODY_3D;
	}
	else if (string[0] == 'O' || string[0] == 'o')
            *prob_type = OPEN_ROTOR;
        else if (string[0] == 'P' || string[0] == 'p')
            *prob_type = PRESSURE_PUMP;

	assert(*prob_type != ERROR_TYPE);
	fclose(infile);
}	/* end read_rg_prob_type */

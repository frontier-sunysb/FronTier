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
*				melting.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*/

#include "../iFluid/iFluid.h"
#include "../iFluid/ifluid_basic.h"
#include "climate.h"

#define		MAX_NUM_VERTEX_IN_CELL		20

typedef struct {
        double center[3];
        double radius;
} TEST_SPHERE_PARAMS;

	/*  Local Application Function Declarations */

static void	melting_flow_driver(Front*,TEST_SPHERE_PARAMS,CARTESIAN*,
			VCARTESIAN*, Incompress_Solver_Smooth_Basis *);
static int      rgbody_vel(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,
                        HYPER_SURF*,double*);
static void	record_1d_growth(Front*,double***,int*);
static void 	plot_growth_data(char*,double**,int);
static double 	temperature_func(double*,COMPONENT,double);
static double 	crystal_curve(POINTER,double*);
static double 	sphere_surf(POINTER,double*);
static boolean 	fractal_dimension(Front*,TEST_SPHERE_PARAMS,double*,double*);
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
	TEST_SPHERE_PARAMS s_params;
	PARAMS eqn_params;
	IF_PARAMS iFparams;
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
	readPhaseParams(&front);
        read_fluid_params(&front);

	if (!RestartRun)
	{
	    if(!eqn_params.no_droplets)
	    {
		printf("Initializing droplets\n");
		level_func_pack.pos_component = LIQUID_COMP;
	    	FT_InitIntfc(&front,&level_func_pack);
                initWaterDrops(&front);
	    	if (debugging("trace")) printf("Passed init water droplets()\n");
	    }
	    else
	    {
	        printf("No droplets contained\n");
	        level_func_pack.func_params = NULL;
                level_func_pack.func = NULL;
	        level_func_pack.pos_component = LIQUID_COMP;
	        level_func_pack.neg_component = SOLID_COMP;
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
	    read_melt_dirichlet_bdry_data(in_name,&front,f_basic);

	front._compute_force_and_torque = compute_ice_particle_force;
        Time_step_factor(&front);
	FT_ReadTimeControl(in_name,&front);

	/* Initialize velocity field function */

	velo_func_pack.func_params = (POINTER)&iFparams;
	velo_func_pack.func = rgbody_vel;
	velo_func_pack.point_propagate = melt_flow_point_propagate;

	FT_InitVeloFunc(&front,&velo_func_pack);


        t_cartesian->initMesh();
	v_cartesian->field = t_cartesian->field;
        v_cartesian->initMesh();

	if (!eqn_params.no_fluid)
	{
	    l_cartesian->initMesh();
	    l_cartesian->findStateAtCrossing = ifluid_find_state_at_crossing;
	}
	if (debugging("sample_velocity"))
	    l_cartesian->initSampleVelocity(in_name);
	if (debugging("sample_temperature"))
            t_cartesian->initSampleTemperature(in_name);

	if (RestartRun)
	{
	    FT_ParallelExchIntfcBuffer(&front);
	    t_cartesian->readFrontInteriorState(restart_state_name);
	    v_cartesian->readFrontInteriorState(restart_state_name);
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
	    if (!eqn_params.no_fluid)
	    {
            	init_fluid_state_func(&front,l_cartesian);
            	l_cartesian->setInitialCondition();
	    }
            if (debugging("trace"))
                printf("Passed iFluid setInitialCondition()\n");
	}

	if (!eqn_params.no_fluid)
	{
	    eqn_params.field->vel = iFparams.field->vel;
	    eqn_params.field->pres = iFparams.field->pres;
	}
	FT_InitVeloFunc(&front,&velo_func_pack);

	if (debugging("trace")) printf("Passed FT_InitVeloFunc()\n");

        /* For geometry-dependent velocity, use first
        * order point propagation function, higher order
        * propagation requires surface propagate, currently
        * in writing, not yet in use. The following override
        * the assigned fourth_order_point_propagate.
        */
	FT_SetGlobalIndex(&front);
	front._point_propagate = melt_flow_point_propagate; 
	
	/* Propagate the front */
	melting_flow_driver(&front,s_params,t_cartesian, v_cartesian, l_cartesian);

	PetscFinalize();
	clean_up(0);
}

static  void melting_flow_driver(
        Front *front,
	TEST_SPHERE_PARAMS s_params,
	CARTESIAN *t_cartesian,
	VCARTESIAN *v_cartesian,
	Incompress_Solver_Smooth_Basis *l_cartesian)
{
        double CFL;
        int  dim = front->rect_grid->dim;
	IF_PARAMS *iFparams;
	PARAMS *eqn_params;
	double frac_dim,radius;
        static LEVEL_FUNC_PACK level_func_pack;
	boolean bdry_reached = NO;
	double **growth_data = NULL;
	int count = 0;

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
            // Front standard output
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
	    }

	    if (dim == 1)
	    {
	    	bdry_reached = NO;
		record_1d_growth(front,&growth_data,&count);
	    }
	    else
	    	bdry_reached = fractal_dimension(front,s_params,&frac_dim,
						&radius);
	    // Problem specific output

	    FT_ResetTime(front);
	    FT_Save(front,out_name);
            t_cartesian->printFrontInteriorState(out_name);
            v_cartesian->printFrontInteriorState(out_name);
	    if (eqn_params->no_fluid == NO)
            	l_cartesian->printFrontInteriorStates(out_name);

	    CondensationPreAdvance(front);
	    FrontPreAdvance(front);
	    FT_Propagate(front);

	    v_cartesian->solve(front->dt);
	    printf("passing solving vapor\n\n");
	    t_cartesian->solve(front->dt);
	    printf("passing solving temperature\n\n");
	    if (eqn_params->no_fluid == NO)
	    	l_cartesian->solve(front->dt);

	    FT_SetTimeStep(front);
	    t_cartesian->setAdvectionDt();
	    v_cartesian->setAdvectionDt();
	    front->dt = std::min(front->dt,CFL*t_cartesian->m_dt);
	    front->dt = std::min(front->dt,CFL*v_cartesian->m_dt);
	    if (eqn_params->no_fluid == NO)
	    {
	    	l_cartesian->setAdvectionDt();
	    	front->dt = std::min(front->dt,CFL*l_cartesian->max_dt);
	    }

	    if (eqn_params->no_fluid == NO)
            	l_cartesian->initMovieVariables();
            t_cartesian->initMovieVariables();
            v_cartesian->initMovieVariables();
	    FT_AddMovieFrame(front,out_name,YES);
        }
        else
	{
	    FT_SetOutputCounter(front);
            t_cartesian->initMovieVariables();
            v_cartesian->initMovieVariables();
            FT_AddMovieFrame(front,out_name,YES);
	}

	FT_TimeControlFilter(front);

        for (;;)
        {
	    CondensationPreAdvance(front);
	    FrontPreAdvance(front);
	    FT_Propagate(front);

	    t_cartesian->solve(front->dt);
	    v_cartesian->solve(front->dt);
	    printf("Passed solving AD equations\n");
	    if (eqn_params->no_fluid == NO)
	    	l_cartesian->solve(front->dt);
	    printf("Passed solving NS equations\n");
	    
	    FT_AddTimeStepToCounter(front);
	    FT_SetTimeStep(front);
	    front->dt = FT_Min(front->dt,CFL*t_cartesian->m_dt);
	    front->dt = FT_Min(front->dt,CFL*v_cartesian->m_dt);
	    if (eqn_params->no_fluid == NO)
	    	front->dt = FT_Min(front->dt,CFL*l_cartesian->max_dt);

            printf("\ntime = %10.9f   step = %7d   dt = %10.9f\n",
                        front->time,front->step,front->dt);
            fflush(stdout);

	    /*For checking the result*/
	    v_cartesian->checkField(); /*output some variables to run-output*/
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
	    	if (eqn_params->no_fluid == NO)
		    l_cartesian->printFrontInteriorStates(out_name);
	    }
            if (FT_IsMovieFrameTime(front))
	    {
		// Front standard output
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
	    if (dim == 1)
	    {
	    	bdry_reached = NO;
		record_1d_growth(front,&growth_data,&count);
	    }
	    else
	    	bdry_reached = fractal_dimension(front,s_params,&frac_dim,
					&radius);

	    if (bdry_reached) front->time_limit_reached = YES;
            if (FT_TimeLimitReached(front))
	    {
		if (dim == 1)
		    plot_growth_data(out_name,growth_data,count);
	    	FT_AddMovieFrame(front,out_name,YES);
		v_cartesian->plotMomentum(out_name);
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
	if (eqn_params->no_fluid == NO)
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


LOCAL double crystal_curve(
        POINTER func_params,
        double *coords)
{

        TEST_SPHERE_PARAMS *s_params = (TEST_SPHERE_PARAMS*)func_params;
        double dist, theta;
	double *cen = s_params->center;
	double radius = s_params->radius;

        dist =   sqrt(sqr(coords[0]-cen[0]) + sqr(coords[1]-cen[1]));
        theta = asin(fabs(coords[1]-0.5)/dist);
	if (coords[0]-0.5 < 0 && coords[1]-0.5 > 0)
	    theta = PI - theta;
	else if (coords[0]-0.5 < 0 && coords[1]-0.5 < 0)
	    theta = PI + theta;
	else if (coords[0]-0.5 > 0 && coords[1]-0.5 < 0)
	    theta = 2*PI - theta;
        return dist - radius + .003*sin(6.0*theta);
}       /* end crystal_curve */

static	double temperature_func(
	double *coords,
	COMPONENT comp,
	double T_l)
{
	if (comp != LIQUID_COMP) return 0.0;
	return T_l;
}	/* end temperature_init_func */

static double sphere_surf(
        POINTER func_params,
        double *coords)
{
        TEST_SPHERE_PARAMS *s_params = (TEST_SPHERE_PARAMS*)func_params;
	double x0,y0,z0,R;
	double distance;

        x0 = s_params->center[0];
        y0 = s_params->center[1];
        z0 = s_params->center[2];
	R = s_params->radius;

	distance = sqrt(sqr(coords[0] - x0) + sqr(coords[1] - y0) +
			sqr(coords[2] - z0)) - R;

        return distance;

}       /* end sphere_surf */

static boolean fractal_dimension(
	Front *front,
	TEST_SPHERE_PARAMS s_params,
	double *frac_dim,
	double *radius)
{
	double coords[MAXD],*center = s_params.center;
	double dist,r_sqr,r_max,r_min = s_params.radius;
	INTERFACE *grid_intfc,*intfc = front->interf;
	int i,j,k,*gmax,dim = intfc->dim;
	POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
	struct Table *T;
	RECT_GRID *grid;
	double *L,*U,*h;
	COMPONENT *gr_comp,comp;
	int N,Nc;
	boolean grid_intfc_made = NO;
	double ratio;
	boolean bdry_reached;
	double margin[MAXD];

	/* Find maximum radius of crystal growth */
	r_max = 0.0;
	grid = computational_grid(intfc);
	L = grid->GL;
	U = grid->GU;
	for (i = 0; i < dim; ++i)
	    margin[i] = 0.01*(U[i] - L[i]);
	bdry_reached = NO;
	next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
	    if (wave_type(hs) < FIRST_PHYSICS_WAVE_TYPE)
	    	continue;
	    r_sqr = 0.0;
	    for (i = 0; i < dim; ++i)
	    {
		if (Coords(p)[i] >= U[i] - margin[i] ||
		    Coords(p)[i] <= L[i] + margin[i])
		    bdry_reached = YES;
	    	r_sqr += sqr(Coords(p)[i] - center[i]);
	    }
	    if (r_max < r_sqr) r_max = r_sqr;
	}
	r_max = sqrt(r_max);
	*radius = r_max;
#if defined (__MPI__)
	pp_global_max(radius,1);
#endif /* defined (__MPI__) */

	/* Preparation for counting */

        if (front->grid_intfc == NULL)
        {
            FT_MakeGridIntfc(front);
            grid_intfc = front->grid_intfc;
            grid_intfc_made = YES;
        }
        else
            grid_intfc = front->grid_intfc;
        grid = &topological_grid(grid_intfc);
        gmax = grid->gmax;
        L = grid->L;
        h = grid->h;
	T = table_of_interface(grid_intfc);
        gr_comp = T->components;

	/* Start counting */
	N = Nc = 0;
	switch (dim)
	{
	case 2:
	    for (i = 0; i <= gmax[0]; ++i)
            for (j = 0; j <= gmax[1]; ++j)
            {
	    	comp = gr_comp[d_index2d(i,j,gmax)];
		coords[0] = L[0] + i*h[0];
                coords[1] = L[1] + j*h[1];
		dist = sqrt(sqr(coords[0] - center[0]) +
			    sqr(coords[1] - center[1]));
	    	if (dist > r_min && dist < r_max)
		{
		    ++N;
		    if (comp == SOLID_COMP)
		    	++Nc;
		}
	    }
	    break;
	case 3:
	    for (i = 0; i <= gmax[0]; ++i)
            for (j = 0; j <= gmax[1]; ++j)
	    for (k = 0; k <= gmax[2]; ++k)
            {
	    	comp = gr_comp[d_index3d(i,j,k,gmax)];
		coords[0] = L[0] + i*h[0];
                coords[1] = L[1] + j*h[1];
		coords[2] = L[2] + k*h[2];
		dist = sqrt(sqr(coords[0] - center[0]) +
			    sqr(coords[1] - center[1]) +
			    sqr(coords[2] - center[2]));
	    	if (dist > r_min && dist < r_max)
		{
		    ++N;
		    if (comp == SOLID_COMP)
		    	++Nc;
		}
	    }
	}
#if defined (__MPI__)
	pp_global_isum(&N,1);
	pp_global_isum(&Nc,1);
#endif /* defined (__MPI__) */
	if (grid_intfc_made)
	    FT_FreeGridIntfc(front);
	ratio = ((double)N)/((double)Nc);
	*frac_dim = (double)dim + log(ratio)/log(h[0]);
	return bdry_reached;
}	/* end fractal_dimension */
	
static void	record_1d_growth(
	Front *front,
	double ***growth_data,
	int *count)
{
	INTERFACE *intfc = front->interf;
	int i,j;
	POINT *p;
	HYPER_SURF *hs;
	HYPER_SURF_ELEMENT *hse;
	static int total_length = 0;
	STATE *sl,*sr;
	
	if (*count >= total_length)
	{
	    static double **tmp_data;
	    total_length += 1000;
	    FT_MatrixMemoryAlloc((POINTER*)&tmp_data,total_length,3,sizeof(double));
	    for (i = 0; i < *count; ++i)
	    for (j = 0; j < 3; ++j)
	    	tmp_data[i][j] = (*growth_data)[i][j];
	    FT_FreeThese(1,*growth_data);
	    *growth_data = tmp_data;
	}

	(*growth_data)[*count][0] = front->time;

	/* Initialize states at the interface */
	next_point(intfc,NULL,NULL,NULL);
	while (next_point(intfc,&p,&hse,&hs))
	{
	    if (is_bdry(p)) continue;
	    FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
	    if (negative_component(hs) == LIQUID_COMP)
	    {
		(*growth_data)[*count][1] = Coords(p)[0];
		(*growth_data)[*count][2] = sl->temperature;
	    }
	    else if (positive_component(hs) == LIQUID_COMP)
	    {
		(*growth_data)[*count][1] = Coords(p)[0];
		(*growth_data)[*count][2] = sr->temperature;
	    }
	}
	*count += 1;
}	/* end record_1d_growth */

static void plot_growth_data(
	char out_name[100],
	double **growth_data,
	int count)
{
	char fname[100];
	FILE *ofile;
	int i;

	sprintf(fname,"%s/posn.xg",out_name);
	ofile = fopen(fname,"w");
	fprintf(ofile,"\"Interface position vs. time\"\n");
	for (i = 0; i < count; ++i)
	    fprintf(ofile,"%f  %f\n",growth_data[i][0],growth_data[i][1]);
	fclose(ofile);

	sprintf(fname,"%s/solt.xg",out_name);
	ofile = fopen(fname,"w");
	fprintf(ofile,"\"Interface temperature vs. time\"\n");
	for (i = 0; i < count; ++i)
	    fprintf(ofile,"%f  %f\n",growth_data[i][0],growth_data[i][2]);
	fclose(ofile);
}	/* end plot_growth_data */

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

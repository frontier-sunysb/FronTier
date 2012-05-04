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
*				crystal.c:
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

#include "melting.h"
#include "melting_basic.h"

/********************************************************************
 *	Level function parameters for the initial interface 	    *
 ********************************************************************/


/********************************************************************
 *	Velocity function parameters for the front	 	    *
 ********************************************************************/
#define		MAX_NUM_VERTEX_IN_CELL		20

typedef struct {
        double center[3];
        double radius;
} TEST_SPHERE_PARAMS;

	/*  Local Application Function Declarations */

static void 	melting_driver(Front*,TEST_SPHERE_PARAMS,CARTESIAN&);
static void	record_1d_growth(Front*,double***,int*);
static void 	plot_growth_data(char*,double**,int);
static double 	temperature_func(double*,COMPONENT,double);
static double 	crystal_curve(POINTER,double*);
static double 	sphere_surf(POINTER,double*);
static boolean 	fractal_dimension(Front*,TEST_SPHERE_PARAMS,double*,double*);
static void 	read_movie_options(char*,PARAMS*);

char *in_name,*restart_state_name,*restart_name,*out_name;
boolean RestartRun;
boolean ReadFromInput;
int RestartStep;
boolean binary = YES;

int main(int argc, char **argv)
{
	static Front front;
	static F_BASIC_DATA f_basic;
	static LEVEL_FUNC_PACK level_func_pack;
	static VELO_FUNC_PACK velo_func_pack;
	TEST_SPHERE_PARAMS s_params;
	PARAMS eqn_params;
	SEED_PARAMS seed_params;
	int dim;

	CARTESIAN       cartesian(front);

	FT_Init(argc,argv,&f_basic);

        PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);

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
        sprintf(restart_name,"%s-nd%s",restart_name,right_flush(pp_mynode(),4));
        sprintf(restart_state_name,"%s-nd%s",restart_state_name,
			right_flush(pp_mynode(),4));
#endif /* defined(__MPI__) */

        FT_ReadSpaceDomain(in_name,&f_basic);

	FT_StartUp(&front,&f_basic);
	FT_InitDebug(in_name);

	readPhaseParams(in_name,&eqn_params);
	read_movie_options(in_name,&eqn_params);

	if (!RestartRun)
	{
	    /* Initialize interface through level function */

	    initPhaseIntfc(in_name,dim,&level_func_pack,&eqn_params);
	    if (debugging("trace")) printf("Passed initPhaseIntfc()\n");

	    FT_InitIntfc(&front,&level_func_pack);
	}
	if (debugging("trace")) printf("Passed FT_InitIntfc()\n");

        Time_step_factor(&front);
	FT_ReadTimeControl(in_name,&front);

	/* Initialize velocity field function */

	front.extra2 = (POINTER)&eqn_params;

	eqn_params.dim = f_basic.dim;

	velo_func_pack.func_params = (POINTER)&eqn_params;
	velo_func_pack.func = NULL;

        cartesian.initMesh();
	if (RestartRun)
	{
	    FT_ParallelExchIntfcBuffer(&front);
	    cartesian.readFrontInteriorState(restart_state_name);
	}
	else
	{
	    cartesian.setInitialCondition();
	}
	if (debugging("trace"))
	    printf("Passed state initialization\n");


	FT_InitVeloFunc(&front,&velo_func_pack);

	if (debugging("trace")) printf("Passed FT_InitVeloFunc()\n");

        /* For geometry-dependent velocity, use first
        * order point propagation function, higher order
        * propagation requires surface propagate, currently
        * in writing, not yet in use. The following override
        * the assigned fourth_order_point_propagate.
        */

        front._point_propagate = melting_point_propagate;

	/* Propagate the front */

	melting_driver(&front,s_params,cartesian);

	PetscFinalize();
	clean_up(0);
}

static  void melting_driver(
        Front *front,
	TEST_SPHERE_PARAMS s_params,
	CARTESIAN &cartesian)
{
        double CFL;
        int  dim = front->rect_grid->dim;
	PARAMS *eqn_params;
	double frac_dim,radius;
	FILE *Radius_file,*FracDim_file;
	boolean bdry_reached = NO;
	double **growth_data = NULL;
	int count = 0;
	Radius_file = FracDim_file = NULL;

	if (debugging("trace"))
	    printf("Entering melting_driver()\n");
	Curve_redistribution_function(front) = expansion_redistribute;
        CFL = Time_step_factor(front);

	eqn_params = (PARAMS*)front->vparams;

	front->hdf_movie_var = NULL;
        cartesian.initMovieVariables();
        if (debugging("trace"))
            printf("Passed initMovieVariables()\n");

        if (!RestartRun)
        {
	    FT_ResetTime(front);
            FT_SetOutputCounter(front);
            // Front standard output
            FT_Save(front,out_name);
            FT_AddMovieFrame(front,out_name,binary);
	    if (dim == 1)
	    	cartesian.oneDimPlot(out_name);

	    if (dim == 1)
	    {
	    	bdry_reached = NO;
		record_1d_growth(front,&growth_data,&count);
	    }
	    else
	    	bdry_reached = fractal_dimension(front,s_params,&frac_dim,
						&radius);
	    if (pp_mynode() == 0 && dim != 1)
	    {
		char fname[200];
		sprintf(fname,"%s/radius",out_name);
	    	Radius_file = fopen(fname,"w");
		sprintf(fname,"%s/FracDim",out_name);
	    	FracDim_file = fopen(fname,"w");

		fprintf(Radius_file,"\"Crystal radius\"\n\n");
		fprintf(Radius_file,"%f  %f\n",front->time,radius);
		fprintf(FracDim_file,"\"Fractal dimension\"\n\n");
		fprintf(FracDim_file,"%f  %f\n",front->time,frac_dim);
		fflush(Radius_file);
		fflush(FracDim_file);
	    }

	    // Problem specific output

	    FT_Propagate(front);
	    FT_SetTimeStep(front);
	    cartesian.setAdvectionDt();
	    front->dt = std::min(front->dt,CFL*cartesian.m_dt);
        }
        else
	    FT_SetOutputCounter(front);
	if (debugging("trace"))
	    printf("Passed second restart check()\n");

	FT_TimeControlFilter(front);

        for (;;)
        {
	    if (debugging("trace"))
		printf("Before FrontAdvance()\n");
	    FT_Propagate(front);
	    if (debugging("trace"))
		printf("After FrontAdvance(), before solve()\n");
	    cartesian.solve(front->dt);
	    if (debugging("trace"))
		printf("After solve()\n");

	    FT_AddTimeStepToCounter(front);

	    FT_SetTimeStep(front);
	    front->dt = std::min(front->dt,CFL*cartesian.m_dt);

            printf("\ntime = %f   step = %7d   dt = %f\n",
                        front->time,front->step,front->dt);
            fflush(stdout);

            if (FT_IsSaveTime(front))
	    {
		// Front standard output
                FT_Save(front,out_name);
		// Problem specific output
		cartesian.printFrontInteriorState(out_name);
	    }
            if (FT_IsMovieFrameTime(front))
	    {
		// Front standard output
        	cartesian.initMovieVariables();
                FT_AddMovieFrame(front,out_name,binary);
		if (dim == 1)
	    	    cartesian.oneDimPlot(out_name);
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
	    if (pp_mynode() == 0 && dim != 1)
	    {
		fprintf(Radius_file,"%f  %f\n",front->time,radius);
		fprintf(FracDim_file,"%f  %f\n",front->time,frac_dim);
		fflush(Radius_file);
		fflush(FracDim_file);
	    }

            if (FT_TimeLimitReached(front))
	    {
		if (dim == 1)
		    plot_growth_data(out_name,growth_data,count);
                break;
	    }
	    /* Output section, next dt may be modified */

	    FT_TimeControlFilter(front);
        }
	if (pp_mynode() == 0 && dim != 1)
	{
	    fclose(Radius_file);
	    fclose(FracDim_file);
	}
}       /* end melting_driver */

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

extern  double getStateTemperature(
        POINTER state)
{
        STATE *T_state = (STATE*)state;
        return T_state->temperature;
}       /* end getStateTemperature */


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

        CursorAfterString(infile,"Type y to make movie of solute:");
        fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
        if (string[0] == 'Y' || string[0] == 'y')
            movie_option->plot_temperature = YES;

	/* Default: not plot cross sectional variables */
        movie_option->plot_cross_section[0] = NO;
        movie_option->plot_cross_section[1] = NO;
        movie_option->plot_cross_section[2] = NO;
        if (params->dim == 3)
        {
            if (CursorAfterStringOpt(infile,
		"Type y to make yz cross section movie:"))
	    {
            	fscanf(infile,"%s",string);
	    	(void) printf("%s\n",string);
            	if (string[0] == 'Y' || string[0] == 'y')
                    movie_option->plot_cross_section[0] = YES;
	    }
            if (CursorAfterStringOpt(infile,
		"Type y to make xz cross section movie:"))
	    {
            	fscanf(infile,"%s",string);
	    	(void) printf("%s\n",string);
            	if (string[0] == 'Y' || string[0] == 'y')
                    movie_option->plot_cross_section[1] = YES;
	    }
            if (CursorAfterStringOpt(infile,
		"Type y to make xy cross section movie:"))
	    {
            	fscanf(infile,"%s",string);
	    	(void) printf("%s\n",string);
            	if (string[0] == 'Y' || string[0] == 'y')
                    movie_option->plot_cross_section[2] = YES;
	    }
        }
        fclose(infile);
}       /* end read_movie_options */

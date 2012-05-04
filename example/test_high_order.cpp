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
*				example16.c:
*
*		User initialization example for Front Package:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*	
*	This example shows a circle in a multiple vortex field. It 
*	demonstrates the reversibility of the front tracking method.
*
*/

#include <FronTier.h>

	/*  Function Declarations */
static void test_propagate(Front*);
static int  test_vel(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,
			HYPER_SURF*,double*);
static void computeError(Front*,int);
static void adjustIntfcPoints(INTERFACE *intfc,MC_PARAMS);
static void setPhysicsHsOrder(INTERFACE*,int);

char *in_name,*restart_state_name,*restart_name,*out_name;
boolean RestartRun;
int RestartStep;
boolean binary = YES;
/********************************************************************
 *	Level function parameters for the initial interface 	    *
 ********************************************************************/

typedef struct {
	        /* equation for line is x^2/a^2 + y^2/b^2 = 1 */
        double C; 
} RADIAL_VEL_PARAMS;


int main(int argc, char **argv)
{
	static Front front;
	static RECT_GRID comp_grid;
	static F_BASIC_DATA f_basic;
	static LEVEL_FUNC_PACK level_func_pack;
	static VELO_FUNC_PACK velo_func_pack;
	RADIAL_VEL_PARAMS rv_params; /* velocity function parameters */
	MC_PARAMS mc_params;
	int i;

	f_basic.dim = 2;	// default
	FT_Init(argc,argv,&f_basic);

	/* Initialize basic computational data */

	for (i = 0; i < f_basic.dim; ++i)
	{
	    f_basic.L[i] = 0.0;	
	    f_basic.U[i] = 1.0;	
	    f_basic.gmax[i] = 100;	
	    f_basic.boundary[i][0] = f_basic.boundary[i][1] = PERIODIC_BOUNDARY;
	}
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
	add_to_debug("high_order_redist");

	/* Initialize interface through level function */

	level_func_pack.neg_component = 1;
	level_func_pack.pos_component = 2;
	level_func_pack.wave_type = FIRST_PHYSICS_WAVE_TYPE;
	mc_params.num_cir = 1;
	FT_VectorMemoryAlloc((POINTER*)&mc_params.rad,mc_params.num_cir,FLOAT);
        FT_MatrixMemoryAlloc((POINTER*)&mc_params.cen,mc_params.num_cir,2,FLOAT);
	for (i = 0; i < f_basic.dim; ++i)
	{
            mc_params.cen[0][i] = 0.5;
	}
        mc_params.rad[0] = 0.15;
	mc_params.dim = f_basic.dim;
	level_func_pack.func_params = (POINTER)&mc_params;
        level_func_pack.func = multi_circle_func;

	FT_InitIntfc(&front,&level_func_pack);
	adjustIntfcPoints(front.interf,mc_params);
	if (f_basic.dim == 3)
	{
	    gview_plot_interface("gview",front.interf);
	}

	Frequency_of_redistribution(&front,GENERAL_WAVE) = 
			f_basic.gmax[0]/20;

	/* Initialize velocity field function */

	rv_params.C = 0.2;
	velo_func_pack.func_params = (POINTER)&rv_params;
	velo_func_pack.func = test_vel;
	velo_func_pack.point_propagate = first_order_point_propagate;

	FT_InitVeloFunc(&front,&velo_func_pack);

	/* Propagate the front */

	test_propagate(&front);

	clean_up(0);
	return 0;
}

static  void test_propagate(
        Front *front)
{
        int ip,im,status,count;
        Front *newfront;
        double dt,dt_frac,CFL;
        boolean is_print_time, is_movie_time, time_limit_reached;
        char s[10];
        double fcrds[MAXD];
        int  dim = front->rect_grid->dim;
	VORTEX_PARAMS *vparams = (VORTEX_PARAMS*)front->vparams;
	int redist_order;

	redist_order = 2;

	front->max_time = 1.0;
	front->max_step = 11;
	front->print_time_interval = 1.0;
	front->movie_frame_interval = 0.02;
	vparams->time = 0.5*front->max_time;

        CFL = Time_step_factor(front);

        printf("dim = %d\n", dim);
	printf("CFL = %f\n",CFL);
	printf("Frequency_of_redistribution(front,GENERAL_WAVE) = %d\n",
		Frequency_of_redistribution(front,GENERAL_WAVE));

	setPhysicsHsOrder(front->interf,redist_order);
	computeError(front,redist_order);
	if (!RestartRun)
	{
            redistribute(front,YES,NO);

            front->time = 0.0;
            front->dt = 0.0;
	    front->step = 0;

	    // Always output the initial interface.
	    FT_Save(front,out_name);
            FT_AddMovieFrame(front,out_name,binary);
            ip = im = 1;

	    // This is a virtual propagation to get maximum front 
	    // speed to determine the first time step.

            status = FrontAdvance(front->dt,&dt_frac,front,&newfront,
                                (POINTER)NULL);
            front->dt = CFL*FrontHypTimeStep(front); 
	}
	else
	{
	    ip = (int)(front->time/front->print_time_interval + 1.0);
            im = (int)(front->time/front->movie_frame_interval + 1.0);
	}

	front->dt = FrontOutputTimeControl(front,
			&is_movie_time,&is_print_time,
			&time_limit_reached,&im,&ip);

        for (;;)
        {
	    /* Propagating interface for time step dt */

            status = FrontAdvance(front->dt,&dt_frac,front,&newfront,
                                (POINTER)NULL);
            assign_interface_and_free_front(front,newfront);

            ++front->step;
            front->time += front->dt;
	    computeError(front,redist_order);

	    //Next time step determined by maximum speed of previous
	    //step, assuming the propagation is hyperbolic and
	    //is not dependent on second order derivatives of
	    //the interface such as curvature, and etc.

            front->dt = CFL*FrontHypTimeStep(front); 

	    /* Output section */

            printf("\ntime = %f   step = %5d   next dt = %f\n",
                        front->time,front->step,front->dt);
            fflush(stdout);

            if (is_print_time || time_limit_reached)
                print_front_output(front,out_name);
            if (is_movie_time || time_limit_reached)
                show_front_output(front,out_name,binary);

            if (time_limit_reached)
                    break;

	    /* Time and step control section */

	    front->dt = FrontOutputTimeControl(front,
			&is_movie_time,
			&is_print_time,
			&time_limit_reached,
			&im,&ip);
        }
        (void) delete_interface(front->interf);
}       /* end test_propagate */

/********************************************************************
 *	Sample (vortex) velocity function for the front    *
 ********************************************************************/


static int test_vel(
        POINTER params,
        Front *front,
        POINT *p,
        HYPER_SURF_ELEMENT *hse,
        HYPER_SURF *hs,
        double *vel)
{
        RADIAL_VEL_PARAMS *rv_params = (RADIAL_VEL_PARAMS*)params;
	double C = rv_params->C;
        int i, dim = front->rect_grid->dim;
	double r = 0.0;
	double speed;
	double *coords = Coords(p);

	speed = C*front->rect_grid->h[0];
	for (i = 0; i < dim; ++i)
	    r += sqr(coords[i] - 0.5);
	r = sqrt(r);
	for (i = 0; i < dim; ++i)
	    vel[i] = C*(coords[i] - 0.5)/r;
}       /* end vortex_vel */

static void computeError(
	Front *front,
	int redist_order)
{
	INTERFACE *intfc = front->interf;
	RADIAL_VEL_PARAMS *vparams = (RADIAL_VEL_PARAMS*)front->vparams;
	double C = vparams->C;
	int i,num_points,dim = intfc->dim;
	double time = front->time;
	double L1_error,L2_error,Li_error;
	double r,R;
	static FILE *L1File,*L2File,*LiFile;
	char fname[100];
	POINT   *p;
        HYPER_SURF_ELEMENT *hse;
        HYPER_SURF         *hs;


	if (L1File == NULL)
	{
	    sprintf(fname,"L1-error-order-%d",redist_order);
	    L1File = fopen(fname,"w");
	    sprintf(fname,"L2-error-order-%d",redist_order);
	    L2File = fopen(fname,"w");
	    sprintf(fname,"Li-error-order-%d",redist_order);
	    LiFile = fopen(fname,"w");
	    fprintf(L1File,"\"L1 Error\"\n");
	    fprintf(L2File,"\"L2 Error\"\n");
	    fprintf(LiFile,"\"Li Error\"\n");
	}

	R = 0.15 + C*time;

	L1_error = L2_error = Li_error = 0.0;
	num_points = 0;
	next_point(intfc,NULL, NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
	    if (wave_type(hs) < FIRST_PHYSICS_WAVE_TYPE)
		continue;
	    r = 0.0;
	    for (i = 0; i < dim; ++i)
		r += sqr(Coords(p)[i] - 0.5);
	    r = sqrt(r);
	    L1_error += fabs(r-R);
	    L2_error += sqr(r-R);
	    Li_error = max(fabs(r-R),Li_error);
	    num_points++;
	}
	L1_error /= (double)num_points;
	L2_error /= (double)num_points;
	L2_error = sqrt(L2_error);

	fprintf(L1File,"%24.18g  %24.18g\n",time,L1_error);
	fflush(L1File);
	fprintf(L2File,"%24.18g  %24.18g\n",time,L2_error);
	fflush(L2File);
	fprintf(LiFile,"%24.18g  %24.18g\n",time,Li_error);
	fflush(LiFile);
}	/* end computeError */

static void adjustIntfcPoints(
	INTERFACE *intfc,
	MC_PARAMS mc_params)
{
	int i,dim = intfc->dim;
	double *center = mc_params.cen[0];
	double r,exact_radius = mc_params.rad[0];
	double crds[i];
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
	next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
	    if (wave_type(hs) < FIRST_PHYSICS_WAVE_TYPE) continue;
	    r = 0.0;
	    for (i = 0; i < dim; ++i)
	    {
		crds[i] = Coords(p)[i] - center[i];
		r += sqr(crds[i]);
	    }
	    r = sqrt(r);
	    for (i = 0; i < dim; ++i)
		Coords(p)[i] = crds[i]*exact_radius/r + center[i];
	}
}	/* end adjustIntfcPoints */

static void setPhysicsHsOrder(
	INTERFACE *intfc,
	int order)
{
	CURVE **c;
	SURFACE **s;

	switch (intfc->dim)
	{
	case 2:
	    for (c = intfc->curves; c && *c; ++c)
	    	if (wave_type(*c) >= FIRST_PHYSICS_WAVE_TYPE)
		    (*c)->redist_order = order;
	case 3:
	    for (s = intfc->surfaces; s && *s; ++s)
	    	if (wave_type(*s) >= FIRST_PHYSICS_WAVE_TYPE)
		    (*s)->redist_order = order;
	}
}	/* end setPhysicsHsOrder */

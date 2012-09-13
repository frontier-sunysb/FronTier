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

#include <iFluid.h>
#include <airfoil.h>

	/*  Function Declarations */
static void airfoil_driver(Front*,Incompress_Solver_Smooth_Basis*);
static void print_hyper_surf_quality(Front*);
static void optimizeElasticMesh(Front*);
static void zero_state(COMPONENT,double*,L_STATE&,int,IF_PARAMS*);
static void reset_front_velocity(Front*);

char *in_name,*restart_state_name,*restart_name,*out_name;
boolean RestartRun;
boolean ReSetTime;
int RestartStep;
boolean binary = YES;
int constrained_propagate;

/********************************************************************
 *	Level function parameters for the initial interface 	    *
 ********************************************************************/

static void xgraph_front(Front*,char*);

int main(int argc, char **argv)
{
	static Front front;
	static F_BASIC_DATA f_basic;
	static LEVEL_FUNC_PACK level_func_pack;
	IF_PARAMS 	iFparams;
	AF_PARAMS 	af_params;

	FT_Init(argc,argv,&f_basic);
	f_basic.size_of_intfc_state = sizeof(STATE);

	//Initialize Petsc before FrontStartUP
        PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
        //if (debugging("trace")) printf("Passed PetscInitialize()\n");

	/*Construct Incompress Solver l_cartesian*/

	Incompress_Solver_Smooth_Basis *l_cartesian = NULL;
	if(f_basic.dim == 2)
	    l_cartesian = new Incompress_Solver_Smooth_2D_Cartesian(front);
	else if(f_basic.dim == 3)
	    l_cartesian = new Incompress_Solver_Smooth_3D_Cartesian(front);

	/* Initialize basic computational data */

        in_name                 = f_basic.in_name;
        restart_state_name      = f_basic.restart_state_name;
        out_name                = f_basic.out_name;
        restart_name            = f_basic.restart_name;
        RestartRun              = f_basic.RestartRun;
        RestartStep             = f_basic.RestartStep;
        ReSetTime             	= f_basic.ReSetTime;

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

	if (debugging("trace")) 
	    (void) printf("Passed FT_StartUp()\n");

        iFparams.dim = f_basic.dim;
        front.extra1 = (POINTER)&iFparams;
        front.extra2 = (POINTER)&af_params;
        read_iFparams(in_name,&iFparams);
        if (debugging("trace")) 
	    (void) printf("Passed read_iFparams()\n");
        read_iF_movie_options(in_name,&iFparams);
        if (debugging("trace")) 
	    (void) printf("Passed read_iF_movie_options()\n");


	setInitialIntfc(&front,&level_func_pack,in_name);
	if (!RestartRun)
	{
	    if (f_basic.dim == 3)
                level_func_pack.set_3d_bdry = YES;
	    FT_InitIntfc(&front,&level_func_pack);
	    if (f_basic.dim == 3 && debugging("trace"))
	    {
		gview_plot_interface("ginit",front.interf);
	    }
	    read_iF_dirichlet_bdry_data(in_name,&front,f_basic);
	    if (f_basic.dim < 3)
            	FT_ClipIntfcToSubdomain(&front);
	}
	else
	{
	    setRestartAirfoilIntfc(&front,&level_func_pack);
	    read_iF_dirichlet_bdry_data(in_name,&front,f_basic);
	}

	/* Time control */
	FT_ReadTimeControl(in_name,&front);

	/* Initialize velocity field function */

	initVelocityFunc(in_name,&front);

	if (!RestartRun)
	    optimizeElasticMesh(&front);
	    
	l_cartesian->findStateAtCrossing = af_find_state_at_crossing;
	l_cartesian->getInitialState = zero_state;
	l_cartesian->initMesh();
	if (debugging("sample_velocity"))
            l_cartesian->initSampleVelocity(in_name);
        if (debugging("trace"))
            (void) printf("Passed l_cartesian->initMesh()\n");
        if (RestartRun)
            l_cartesian->readFrontInteriorStates(restart_state_name);
        else
	{
            l_cartesian->setInitialCondition();
	    if (!RestartRun || ReSetTime)
	    	reset_front_velocity(&front);
	}
        if (debugging("trace"))
            (void) printf("Passed state initialization()\n");

	/* Propagate the front */

	if (!RestartRun)
	    set_equilibrium_mesh(&front);
	airfoil_driver(&front,l_cartesian);

	clean_up(0);
}

static  void airfoil_driver(
        Front *front,
	Incompress_Solver_Smooth_Basis *l_cartesian)
{
        double CFL;
        int  dim = front->rect_grid->dim;
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;

        CFL = Time_step_factor(front);
	Frequency_of_redistribution(front,GENERAL_WAVE) = 100;

	(void) printf("Frequency_of_redistribution(front,GENERAL_WAVE) = %d\n",
		Frequency_of_redistribution(front,GENERAL_WAVE));

	if (!RestartRun || ReSetTime)
	{
	    FT_ResetTime(front);

	    // Always output the initial interface.
	    if (dim == 2)
	    {
	    	xgraph_front(front,out_name);
	    }

	    if (debugging("trace"))
                (void) printf("Calling FT_Save()\n");
		FT_Save(front,out_name);
            if (debugging("trace"))
                (void) printf("Calling printFrontInteriorStates()\n");
                l_cartesian->printFrontInteriorStates(out_name);
	    if (debugging("trace"))
                (void) printf("Calling initMovieVariable()\n");
            l_cartesian->initMovieVariables();
            if (debugging("trace"))
                (void) printf("Calling FT_AddMovieFrame()\n");
            FT_AddMovieFrame(front,out_name,binary);
	    FT_Propagate(front);
	    if (!af_params->no_fluid)
	    {
            	if (debugging("trace")) printf("Calling ifluid solve()\n");
            	l_cartesian->solve(front->dt);
            	if (debugging("trace")) printf("Passed ifluid solve()\n");
	    }
	    print_airfoil_stat(front,out_name);

            FT_SetOutputCounter(front);
	    FT_SetTimeStep(front);
	    front->dt = std::min(front->dt,CFL*l_cartesian->max_dt);
	}
	else
	{
	    FT_SetOutputCounter(front);
	}

	FT_TimeControlFilter(front);

	if (debugging("step_size"))
                (void) printf("Time step from start: %f\n",front->dt);
        for (;;)
        {
	    /* Propagating interface for time step dt */

	    if (debugging("trace"))
                (void) printf("Before FT_Propagate()\n");

	    if (!af_params->no_fluid)
	    {
	    	coating_mono_hyper_surf(front);
	    	l_cartesian->applicationSetComponent();
	    }
            FT_Propagate(front);
	    if (!af_params->no_fluid)
	    {
	    	coating_mono_hyper_surf(front);
	    	l_cartesian->applicationSetStates();
	    }
            if (debugging("trace")) printf("Passed FT_Propagate()\n");

	    if (!af_params->no_fluid)
	    {
            	if (debugging("trace")) printf("Calling ifluid solve()\n");
            	l_cartesian->solve(front->dt);
            	if (debugging("trace")) printf("Passed ifluid solve()\n");
	    }
	    else
		l_cartesian->max_dt = HUGE;
	    if (debugging("trace"))
            {
                (void) printf("After solve()\n");
                (void) print_storage("at end of time step","trace");
            }

	    FT_AddTimeStepToCounter(front);

	    //Next time step determined by maximum speed of previous
	    //step, assuming the propagation is hyperbolic and
	    //is not dependent on second order derivatives of
	    //the interface such as curvature, and etc.

	    FT_SetTimeStep(front);
            if (debugging("step_size"))
                (void) printf("Time step from FrontHypTimeStep(): %f\n",
					front->dt);
            front->dt = std::min(front->dt,CFL*l_cartesian->max_dt);
            if (debugging("step_size"))
                (void) printf("Time step from l_cartesian->max_dt(): %f\n",
					front->dt);

	    /* Output section */

            (void) printf("\ntime = %f   step = %5d   next dt = %f\n",
                        front->time,front->step,front->dt);
            fflush(stdout);
	    print_airfoil_stat(front,out_name);

            if (FT_IsSaveTime(front))
	    {
		FT_Save(front,out_name);
                l_cartesian->printFrontInteriorStates(out_name);
	    }
	    if (debugging("trace"))
                (void) printf("After print output()\n");
            if (FT_IsMovieFrameTime(front))
	    {
		if (debugging("trace"))
                    (void) printf("Calling initMovieVariable()\n");
                l_cartesian->initMovieVariables();
                if (debugging("trace"))
                    (void) printf("Calling FT_AddMovieFrame()\n");
                FT_AddMovieFrame(front,out_name,binary);
	    }

            if (FT_TimeLimitReached(front))
                    break;

	    /* Time and step control section */

	    FT_TimeControlFilter(front);
        }
        (void) delete_interface(front->interf);
}       /* end airfoil_driver */

extern void test_tan_curve_propagate(
	Front           *fr,
        Front           *newfr,
        INTERFACE       *intfc,
        CURVE           *oldc,
        CURVE           *newc,
        double           dt)
{
	char xname[100];
	FILE *xfile;
	BOND *b,*oldb;
	double *dist_l,*dist_r;
	double **tan_l,**tan_r;
	POINT *p,*p_nb,**pts;
	int i,n,np,k;
	int dim = Dimension(intfc);
	double lambda;
	double **dp;
	double total_length,fixed_length;
	
	if (wave_type(newc) != ELASTIC_BOUNDARY)
	    return;

	if (debugging("airfoil"))
	    printf("Entering test_tan_curve_propagate()\n");
	lambda = 0.0;

	np = FT_NumOfCurvePoints(newc);
	FT_VectorMemoryAlloc((POINTER*)&dist_l,np,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&dist_r,np,sizeof(double));
	FT_MatrixMemoryAlloc((POINTER*)&tan_l,np,MAXD,sizeof(double));
	FT_MatrixMemoryAlloc((POINTER*)&tan_r,np,MAXD,sizeof(double));
	FT_MatrixMemoryAlloc((POINTER*)&dp,np,MAXD,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&pts,np,sizeof(POINT*));

	//calculate the total length of the newc 
	total_length = 0.0;
	for (b = newc->first; b != NULL; b = b->next)
	    total_length += bond_length(b);
	printf("Entering: %24.18g\t", total_length);

	b = newc->first;
	p = pts[0] = b->start;
	p_nb = b->end;
	dist_l[0] = 0.0;
	dist_r[0] = bond_length(b) - b->length0;
	for (i = 0; i < dim; ++i)
	{
	    tan_r[0][i] = Coords(p_nb)[i] - Coords(p)[i];
	    dp[0][i] = lambda*dist_r[0]*tan_r[0][i];
	}

	n = 1;
	for (b = newc->first; b != newc->last; b = b->next)
	{
	    p = pts[n] = b->end; 	
	    p_nb = b->start;
	    dist_l[n] = bond_length(b) - b->length0;
	    for (i = 0; i < dim; ++i)
		tan_l[n][i] = Coords(p_nb)[i] - Coords(p)[i];
	    p_nb = b->next->end;
	    dist_r[n] = bond_length(b->next) - b->next->length0;
	    for (i = 0; i < dim; ++i)
	    {
		tan_r[n][i] = Coords(p_nb)[i] - Coords(p)[i];
	    	dp[n][i] = lambda*(dist_l[n]*tan_l[n][i] +
				   dist_r[n]*tan_r[n][i]);
	    }
	    ++n;
	}
	b = newc->last;
	p = pts[n] = b->end;
	p_nb = b->start;
	dist_l[n] = bond_length(b) - b->length0;
	dist_r[n] = 0.0;
	for (i = 0; i < dim; ++i)
	{
	    tan_l[n][i] = Coords(p_nb)[i] - Coords(p)[i];
	    dp[n][i] = lambda*dist_l[n]*tan_l[n][i];
	}

	k = np / 2;
	b = newc->first;
	for (n = 1; n < k; ++n)
	    b = b->next;
	for (n = 1; n <= k; ++n)
	{
	    for (i = 0; i < dim; ++i)
		lambda += tan_l[k-n+1][i] * tan_l[k-n+1][i];
	    lambda = sqrt(lambda);
	    for (i = 0; i < dim; ++i)
		Coords(pts[k-n])[i] = Coords(pts[k-n+1])[i] + 
				b->length0*tan_l[k-n+1][i]/lambda;
	    lambda = 0.0;
	    b = b->prev;
	}

	b = newc->first;
	for (n = 0; n < k; ++n)
	    b = b->next;
	for (n = 1; n <= k; ++n)
	{
	    for (i = 0; i < dim; ++i)
		lambda += tan_r[k+n-1][i] * tan_r[k+n-1][i];
	    lambda = sqrt(lambda);
	    for (i = 0; i < dim; ++i)
	    	Coords(pts[k+n])[i] = Coords(pts[k+n-1])[i] + 
			b->length0*tan_r[k+n-1][i]/lambda;
	    lambda = 0.0;
	    b = b->next;
	}

	fixed_length = total_length = 0.0;
	for (b = newc->first, oldb = oldc->first; b != NULL; 
		b = b->next, oldb = oldb->next)
	{
	    set_bond_length(b,dim);
	    total_length += bond_length(b);
	    fixed_length += b->length0;
	}

	printf("Exiting: %24.18g\t", total_length);

	if (fr->step%50 == 0)
	{
	    sprintf(xname,"curves-%d",fr->step);
	    xfile = fopen(xname,"w");
	    xgraph_curve(xfile,newc,XY_PLANE);
	    xgraph_curve(xfile,newc,XY_PLANE);
	    fclose(xfile);
	}
	if (debugging("airfoil"))
	{
	    printf("total_length = %f  fixed_length = %f\n",
			total_length,fixed_length);
	    printf("Leaving test_tan_curve_propagate()\n");
	}
	FT_FreeThese(6,dist_l,dist_r,tan_l,tan_r,dp,pts);
}	/* end test_tan_curve_propagate */

static void xgraph_front(
	Front *front,
	char *outname)
{
	char fname[100];
	sprintf(fname,"%s/intfc-%s",outname,right_flush(front->step,4));
	xgraph_2d_intfc(fname,front->interf);
}	/* end xgraph_front */

static void print_hyper_surf_quality(
	Front *front)
{
	INTERFACE *intfc = front->interf;
	int dim = Dimension(intfc);
	CURVE **c,*curve;
	SURFACE **s,*surf;
	BOND *bond;
	TRI *tri;
	double max_area,min_area,max_length,min_length;
	double scaled_area,len[MAXD];
	int i;
	RECT_GRID *gr = &topological_grid(intfc);
	double *h = gr->h;

	switch (dim)
	{
	case 2:
	    max_length = 0.0;
	    min_length = HUGE;
	    for (c = intfc->curves; c && *c; ++c)
	    {
		if (wave_type(*c) != ELASTIC_BOUNDARY)
		    continue;
		curve = *c;
		for (bond = curve->first; bond != NULL; bond = bond->next)
		{
		    if (scaled_bond_length(bond,h,dim) > max_length)
			max_length = scaled_bond_length(bond,h,dim);
		    if (scaled_bond_length(bond,h,dim) < min_length)
			min_length = scaled_bond_length(bond,h,dim);
		}
	    }
	    (void) printf("\n\nElastic surface quality:\n");
	    (void) printf("min_length = %f\n",min_length);
	    (void) printf("max_length = %f\n",max_length);
	    (void) printf("\n\n");
	    break;
	case 3:
	    max_length = 0.0;
	    min_length = HUGE;
	    for (c = intfc->curves; c && *c; ++c)
	    {
		if (hsbdry_type(*c) != MONO_COMP_HSBDRY)
		    continue;
		curve = *c;
		for (bond = curve->first; bond != NULL; bond = bond->next)
		{
		    if (scaled_bond_length(bond,h,dim) > max_length)
			max_length = scaled_bond_length(bond,h,dim);
		    if (scaled_bond_length(bond,h,dim) < min_length)
			min_length = scaled_bond_length(bond,h,dim);
		}
	    }
	    (void) printf("\n\nElastic curve quality:\n");
	    (void) printf("min_scaled_length = %14.10f\n",min_length);
	    (void) printf("max_scaled_length = %14.10f\n",max_length);

	    max_length = 0.0;
	    min_length = HUGE;
	    for (c = intfc->curves; c && *c; ++c)
	    {
		if (hsbdry_type(*c) != STRING_HSBDRY)
		    continue;
		curve = *c;
		for (bond = curve->first; bond != NULL; bond = bond->next)
		{
		    if (scaled_bond_length(bond,h,dim) > max_length)
			max_length = scaled_bond_length(bond,h,dim);
		    if (scaled_bond_length(bond,h,dim) < min_length)
			min_length = scaled_bond_length(bond,h,dim);
		}
	    }
	    (void) printf("\nElastic string quality:\n");
	    (void) printf("min_scaled_length = %14.10f\n",min_length);
	    (void) printf("max_scaled_length = %14.10f\n",max_length);

	    max_area = max_length = 0.0;
	    min_area = min_length = HUGE;
	    for (s = intfc->surfaces; s && *s; ++s)
	    {
		if (wave_type(*s) != ELASTIC_BOUNDARY)
		    continue;
		surf = *s;
		for (tri = first_tri(surf); !at_end_of_tri_list(tri,surf); 
				tri = tri->next)
		{
		    scaled_tri_params(tri,h,&scaled_area,len);
		    if (scaled_area > max_area) max_area = scaled_area;
		    if (scaled_area < min_area) min_area = scaled_area;
		    for (i = 0; i < 3; ++i)
		    {
			if (len[i] > max_length) 
			    max_length = len[i];
			if (len[i] < min_length) 
			    min_length = len[i];
		    }
		}
	    }
	    (void) printf("\nElastic surface quality:\n");
	    (void) printf("min_scaled_area = %14.10f\n",min_area);  
	    (void) printf("max_scaled_area = %14.10f\n",max_area); 
	    (void) printf("min_scaled_tri_side = %14.10f\n",sqrt(min_length));
	    (void) printf("max_scaled_tri_side = %14.10f\n",sqrt(max_length));
	    (void) printf("\n\n");
	    break;
	}
}	/* end print_hyper_surf_quality */

static void zero_state(
        COMPONENT comp,
        double *coords,
        L_STATE& state,
        int dim,
        IF_PARAMS *iFparams)
{
        int i;
        for (i = 0; i < dim; ++i)
            state.m_U[i] = 0.0;
}       /* end zero_state */

static void reset_front_velocity(Front *front)
{
	INTERFACE *intfc = front->interf;
	POINT *p;
	HYPER_SURF_ELEMENT *hse;
	HYPER_SURF *hs;
	STATE *sl,*sr;
	int i,dim = front->rect_grid->dim;

	next_point(intfc,NULL,NULL,NULL);
	while (next_point(intfc,&p,&hse,&hs))
	{
	    FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
	    for (i = 0; i < dim; ++i)
	    {
		p->vel[i] = 0.0;
		sl->vel[i] = sr->vel[i] = 0.0;
		sl->Impct[i] = sr->Impct[i] = 0.0;
	    }
	}
}	/* end reset_velocity */

static void optimizeElasticMesh(
	Front *front)
{
	INTERFACE *intfc = front->interf;
	RECT_GRID *gr = computational_grid(intfc);
	boolean nothing_done;
	int i,status;
	CURVE **c,*curve;
	SURFACE **s,*surf;
	SCALED_REDIST_PARAMS scaled_redist_params;
	int old_string_pts,new_string_pts,old_canopy_pts,new_canopy_pts;
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	int num_opt_round;

	if (debugging("optimize_intfc"))
	{
	    (void) printf("Quality of mesh before optimization:\n");
	    print_hyper_surf_quality(front);
	    (void) printf("Checking consistency of interface\n");
	    consistent_interface(front->interf);
	    (void) printf("Checking completed\n");
	    gview_plot_interface("gview-before-optimize",intfc);
	    if (debugging("no_optimize"))
		clean_up(0);
	}
	if (debugging("no_optimize")) return;
	if (gr->dim == 2) return;

	num_opt_round = af_params->num_opt_round;
	scaled_redist_params.min_scaled_bond_length = 0.45;
	scaled_redist_params.max_scaled_bond_length = 1.05;

	scaled_redist_params.min_scaled_tri_area = 0.1083;
	scaled_redist_params.max_scaled_tri_area = 0.4330;
	scaled_redist_params.min_scaled_side_length = 0.45;
	scaled_redist_params.max_scaled_side_length = 1.05;
	scaled_redist_params.aspect_tol = 3.0;

	old_string_pts = old_canopy_pts = 0;
	for (s = intfc->surfaces; s && *s; ++s)
	    if (wave_type(*s) == ELASTIC_BOUNDARY)
		old_canopy_pts += FT_NumOfSurfPoints(*s);
	for (c = intfc->curves; c && *c; ++c)
	    if (hsbdry_type(*c) == STRING_HSBDRY)
		old_string_pts += FT_NumOfCurvePoints(*c) - 2;

	for (i = 0; i < num_opt_round; ++i)
	{
	    status = YES;
	    if (debugging("optimize_intfc"))
		(void) printf("Optimization round %d\n",i);
	    for (c = intfc->curves; c && *c; ++c)
	    {
	    	if (hsbdry_type(*c) != MONO_COMP_HSBDRY &&
		    hsbdry_type(*c) != STRING_HSBDRY)
		    continue;
	    	curve = *c;
	    	nothing_done = FT_OptimizeCurveMesh(front,curve,
				scaled_redist_params);
	    	status *= (int)nothing_done;
	    }
	    for (s = intfc->surfaces; s && *s; ++s)
	    {
	    	if (wave_type(*s) != ELASTIC_BOUNDARY)
		    continue;
	    	surf = *s;
	    	nothing_done = FT_OptimizeSurfMesh(front,surf,
				scaled_redist_params);
	    	status *= (int)nothing_done;
	    }
	    FT_ParallelExchIntfcBuffer(front);
	    if (debugging("optimize_intfc"))
	    {
		(void) printf("Quality of mesh after %d-th round:\n",i);
	    	print_hyper_surf_quality(front);
		(void) printf("Checking consistency of interface\n");
		consistent_interface(front->interf);
		(void) printf("After checking\n");
	    }
	    if (status) break;
	}

	new_string_pts = new_canopy_pts = 0;
	for (s = intfc->surfaces; s && *s; ++s)
	    if (wave_type(*s) == ELASTIC_BOUNDARY)
		new_canopy_pts += FT_NumOfSurfPoints(*s);
	for (c = intfc->curves; c && *c; ++c)
	    if (hsbdry_type(*c) == STRING_HSBDRY)
		new_string_pts += FT_NumOfCurvePoints(*c) - 2;
	af_params->m_s *= (double) old_canopy_pts/new_canopy_pts;
	af_params->m_c *= (double) old_string_pts/new_string_pts;
	if (debugging("optimize_intfc"))
	{
	    gview_plot_interface("gview-after-optimize",intfc);
	    clean_up(0);
	}
}	/* end optimizeElasticMesh */

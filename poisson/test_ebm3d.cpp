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
*		test_ebm3d.c
* modified from example3d/example01.c by robert shuqiang wang
*
*		User initialization example for Front Package:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*/

#include <FronTier.h>
#include <solver_petsc.h>
#include <ebm3d.h>
#include <mpi.h>

	/*  Function Declarations */
static void test_propagate(Front*);
static double sphere_func(POINTER,double*);
static int test_curvature_vel(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,
			HYPER_SURF*,double*);

char *in_name,*restart_state_name,*restart_name,*out_name;
boolean RestartRun;
int RestartStep;

/********************************************************************
 *	Level function parameters for the initial interface 	    *
 ********************************************************************/
//Dummbell shape interface, two spheres of radius R centered at (x0,y,z)
//and (x1,y,z) connected by a cylinder of radius rr along x0->x1. Assume
//x0<x1;
typedef struct {
        double center[3];
	double radius;
} TEST_SPHERE_PARAMS;


typedef struct {
	int dim;
	double coeff;
	double epsilon;
} TEST_CURV_PARAMS;

/********************************************************************
 *	Velocity function parameters for the front	 	    *
 ********************************************************************/


int main(int argc, char **argv)
{
	static Front front;
	static RECT_GRID comp_grid;
	F_BASIC_DATA f_basic;
	TEST_SPHERE_PARAMS s_params;
	static LEVEL_FUNC_PACK level_func_pack;
	VELO_FUNC_PACK velo_func_pack;
	TEST_CURV_PARAMS curv_params; /* velocity function parameters */
	Locstate  sl;

        //MPI_Init(&argc, &argv);
        pp_init(&argc, &argv);
	PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
             
        
	f_basic.dim = 3;	
        
        // added by robert wang
        // it happens that F_BASIC_DATA.subdomains[] are not set which cause the
        // the program to crash later.
        for(int i=0; i<f_basic.dim; i++)
            f_basic.subdomains[i] = 1;
        
	

	/* Initialize basic computational data */

	int size = 10;
	
	f_basic.L[0] =-1;	f_basic.L[1] =-1; 	f_basic.L[2] =-1;
	f_basic.U[0] = 1.0;	f_basic.U[1] = 1.0; 	f_basic.U[2] = 1.0;
	f_basic.gmax[0] = size;	f_basic.gmax[1] = size; f_basic.gmax[2] = size;
        
        FT_Init(argc,argv,&f_basic);        

	f_basic.boundary[0][0] = f_basic.boundary[0][1] = DIRICHLET_BOUNDARY;
	f_basic.boundary[1][0] = f_basic.boundary[1][1] = DIRICHLET_BOUNDARY;
	f_basic.boundary[2][0] = f_basic.boundary[2][1] = DIRICHLET_BOUNDARY;
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
	    s_params.center[0] = 0.5;
	    s_params.center[1] = 0.5;
	    s_params.center[2] = 0.5;
	    s_params.radius = 0.2;

	    level_func_pack.neg_component = 1;
	    level_func_pack.pos_component = 2;
	    level_func_pack.func_params = (POINTER)&s_params;
	    level_func_pack.func = sphere_func;

	    FT_InitIntfc(&front,&level_func_pack);
	}

	/* Initialize velocity field function */

/*
        curv_params.dim = 3;
        curv_params.coeff = 0.1;
        curv_params.epsilon = 0.0001;

        velo_func_pack.func_params = (POINTER)&curv_params;
        velo_func_pack.func = test_curvature_vel;
*/

	//FT_InitVeloFunc(&front,&velo_func_pack);

	/* For geometry-dependent velocity, use first 
	* order point propagation function, higher order
	* propagation requires surface propagate, currently
	* in writing, not yet in use. The following override
	* the assigned fourth_order_point_propagate.
	*/

	//front._point_propagate = first_order_point_propagate;

	/* Propagate the front */

	//test_propagate(&front);
        
        printf("hello, world!\n");
        
        EBM3D_CELL ebm_cell;
/*
        ebm_cell.test();
        clean_up(0);
        return 0;
*/
        
        EBM3D_LAPLACE laplace;
        PETSc solver;
        solver.SetTol(10e-9);
        EBM3D_CARTESIAN ebm_cartesian(&front, &solver, &laplace);        
        //ebm_cartesian.debug_get4PointsSolutionInterpCoeffs();
	//ebm_cartesian.debug_get10PointsSolutionInterpCoeffs();
        //ebm_cartesian.debug_get6PointsSolutionInterpCoeffs_Plane();

	//printf("solve the matrix\n");
        ebm_cartesian.solve();
	
        //ebm_cartesian.saveInterface_Tecplot("intfc.plt");
        //ebm_cartesian.saveReconstructedInterface_Tecplot("reconstructed.plt");
	
	//ebm_cartesian.saveReconstructedInterface_Tecplot();      // for debugging with the intfc normal.
        //ebm_cartesian.saveStates_Tecplot(-1, -1, f_basic.gmax[2]/2);
	//ebm_cartesian.saveStates_Tecplot();
	ebm_cartesian.saveStates_VTK();
        //ebm_cartesian.saveComponent_Tecplot("components.plt");

        
        printf("------------------ end of output ------------------\n\n");
        
	clean_up(0);
	return 0;
}

/*
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

	front->max_time = 3; 
	front->max_step = 1000;
	front->print_time_interval = 0.6;
	front->movie_frame_interval = 0.1;

        CFL = Time_step_factor(front) = 0.1;

        printf("dim = %d\n", dim);
	printf("CFL = %f\n",CFL);
	printf("Frequency_of_redistribution(front,GENERAL_WAVE) = %d\n",
		Frequency_of_redistribution(front,GENERAL_WAVE));

	if (!RestartRun)
	{
            redistribute(front,YES,NO);

            front->time = 0.0;
            front->dt = 0.0;
	    front->step = 0;

	    // Always output the initial interface.
	    FT_Save(front,out_name);
            FT_AddMovieFrame(front,out_name);
            ip = im = 1;

	    // This is a virtual propagation to get maximum front 
	    // speed to determine the first time step.

            status = FrontAdvance(front->dt,&dt_frac,front,&newfront,
                                (POINTER)NULL);
            front->dt = CFL*FrontHypTimeStep(front); 
	    front->dt = FrontOutputTimeControl(front,
			&is_movie_time,
			&is_print_time,
			&time_limit_reached,
			&im,&ip);
	}
	else
	{
	    ip = (int)(front->time/front->print_time_interval + 1.0);
            im = (int)(front->time/front->movie_frame_interval + 1.0);
	    time_limit_reached = NO;
	}

        for (;;)
        {
	    // Propagating interface for time step dt 

            status = FrontAdvance(front->dt,&dt_frac,front,&newfront,
                                (POINTER)NULL);
            assign_interface_and_free_front(front,newfront);

            ++front->step;
            front->time += front->dt;

	    //Next time step determined by maximum speed of previous
	    //step, assuming the propagation is hyperbolic and
	    //is not dependent on second order derivatives of
	    //the interface such as curvature, and etc.

            front->dt = CFL*FrontHypTimeStep(front); 

            printf("\ntime = %f   step = %5d   next dt = %f\n",
                        front->time,front->step,front->dt);
            fflush(stdout);

            if (is_print_time || time_limit_reached)
                print_front_output(front,out_name);
            if (is_movie_time || time_limit_reached)
                show_front_output(front,out_name);

            if (time_limit_reached)
                    break;

	    // Output section, next dt may be modified 

	    front->dt = FrontOutputTimeControl(front,
			&is_movie_time,
			&is_print_time,
			&time_limit_reached,
			&im,&ip);
        }
        (void) delete_interface(front->interf);
}       // end test_propagate 
*/

/********************************************************************
 *	Sample (dummbell 3D) level function for the initial interface    *
 ********************************************************************/

static double sphere_func(
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

}       /* end sphere_func */

static int test_curvature_vel(
        POINTER params,
        Front *front,
        POINT *p,
        HYPER_SURF_ELEMENT *hse,
        HYPER_SURF *hs,
        double *vel)
{
	TEST_CURV_PARAMS *curv_params = (TEST_CURV_PARAMS*) params;
        int i;
        double coeff,epsilon,eps;
        double kappa;
        double nor[MAXD];

        coeff = curv_params->coeff;
        epsilon = curv_params->epsilon;

        normal(p,hse,hs,nor,front);
        kappa = fabs(mean_curvature_at_point(p,hse,hs,front));

        for (i = 0; i < curv_params->dim; ++i)
        {
            vel[i] = nor[i]*(coeff - epsilon*kappa);
        }
}


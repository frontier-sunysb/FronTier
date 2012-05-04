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
 *          test_ebm2d.c
 * This file is mofidied from poison.c to test the embedded boundary method 
 * implemented in EBM2D.h/EBM2D.c.
 *
 * see also poison.c
 */

#include <ebm2d.h>
#include <frontier_ppgrid.h>
#include <geometry.h>

#include <FronTier.h>
// the following is used for "vector"
#undef _VMALLOC_H
#include <FronTier/util/vmalloc.h>              

#include <poisson.h>
#include <solver_petsc.h>
#include <mpi.h>

/********************************************************************
 *	Level function parameters for the initial interface 	    *
 ********************************************************************/

	/*  Local Application Function Declarations */

struct _EXP_PARAMS {
	double K[2];
};
typedef struct _EXP_PARAMS EXP_PARAMS;

static double	exp_func(POINTER,double*);
static double	sin_func(POINTER,double*);
static void	exp_func_flux(POINTER,double*,double*);
static void	sin_func_flux(POINTER,double*,double*);
static double	exp_solution(POINTER,double*);
static double	sin_solution(POINTER,double*);

static double 	circle_curve(POINTER,double*);
static void 	show_intfc_solute_states(INTERFACE*);
static void 	print_vtk_solution(const char*,POINTER,Front*);
static void 	print_error_analysis(const char*,POINTER,Front*);

char *in_name,*restart_state_name,*restart_name,*out_name;
boolean RestartRun;
int RestartStep;

int main(int argc, char **argv)
{
	static Front front;
	static F_BASIC_DATA f_basic;
	static LEVEL_FUNC_PACK level_func_pack;
	EXP_PARAMS sparams; /* function parameters */
	LAPLACE laplace;
	int i,j,num_nodes,*gmax;
	RECT_GRID *grid;
	char solute_name[100];
	//const char *dirname = "/nfs/user02/linli/FronTier/poison/";
        
	//pp_init(&argc, &argv);
	//PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
        
 	/* Initialize basic computational data */

	f_basic.dim = 2;	// Dimension
	for(int i=0; i<f_basic.dim; i++)
            f_basic.subdomains[i] = 1;
        
	/* Global Computational Domain */

	f_basic.L[0] = -1.0;	f_basic.L[1] = -1.0;
	f_basic.U[0] = 1.0;	f_basic.U[1] = 1.0;

	/* Global Rectangular Mesh */

	f_basic.gmax[0] = 10;	f_basic.gmax[1] = 10;
	//f_basic.gmax[0] = 40;	f_basic.gmax[1] = 40;
               
        FT_Init(argc,argv,&f_basic);

	/* Boundary conditions */

	f_basic.boundary[0][0] = f_basic.boundary[0][1] = NEUMANN_BOUNDARY;
	f_basic.boundary[1][0] = f_basic.boundary[1][1] = NEUMANN_BOUNDARY;

	/* Set front state size */

	f_basic.size_of_intfc_state = sizeof(STATE);

        in_name                 = f_basic.in_name;
        restart_state_name      = f_basic.restart_state_name;
        out_name                = f_basic.out_name;
        restart_name            = f_basic.restart_name;
        RestartRun              = f_basic.RestartRun;
        RestartStep             = f_basic.RestartStep;

        sprintf(restart_state_name,"%s-state.ts%s",restart_name,
                        right_flush(RestartStep,7));
        sprintf(restart_name,"%s.ts%s",restart_name,right_flush(RestartStep,7));
//#if defined(__MPI__)
        sprintf(restart_name,"%s-nd%s",restart_name,right_flush(pp_mynode(),4));
        sprintf(restart_state_name,"%s-nd%s",restart_state_name,
                        right_flush(pp_mynode(),4));
//#endif /* defined(__MPI__) */

	FT_StartUp(&front,&f_basic);
	PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL); 

	/* Initialize interface through level function */

	level_func_pack.neg_component = 1;
	level_func_pack.pos_component = 2;
	level_func_pack.func = circle_curve;

	printf("Calling FT_InitIntfc()\n");
	FT_InitIntfc(&front,&level_func_pack);

	printf("Calling redistribute()\n");
        redistribute(&front,YES,NO);

	/* Initialize velocity field function */

	front.interf->elliptic_comp = 2;
	//laplace.sparams = (POINTER)NULL;
	//laplace.sfunc     = sin_func;
	//laplace.flux_func = sin_func_flux;
	//laplace.solution  = sin_solution;
	//laplace.sparams   = (POINTER)&sparams;
	//laplace.sfunc     = exp_func;
	//laplace.flux_func = exp_func_flux;
	//laplace.solution  = exp_solution;

	/* Make grid interface */

	printf("Calling make_grid_intfc()\n");
	front.grid_intfc = make_grid_intfc(front.interf,
			EXPANDED_DUAL_GRID,NULL);

	/* This should be merged into above */

	
	printf("Calling make_emb_grid_intfc()\n");
	front.emb_grid_intfc = make_emb_grid_intfc(front.interf);

             
	//grid = &topological_grid(front.grid_intfc);
	//gmax = grid->gmax;
	//num_nodes = 1;
	//for (i = 0; i < 2; ++i)
        //    num_nodes *= (gmax[i] + 1);

	//laplace.solute = (double*)malloc(num_nodes*sizeof(double));
        //laplace.solute = new double[num_nodes];
	//laplace.D = 1.0;

	//embed_bdry_poison((POINTER)&laplace,&front);

	//scatter_top_grid_float_array(laplace.solute,&front);

	//print_vtk_solution(dirname,(POINTER)&laplace,&front);
	//print_error_analysis(dirname,(POINTER)&laplace,&front);
        
        //free(laplace.solute);
        
    
	//FronTierPPGrid pp;
	//pp.test(&front);
	//clean_up(0);

	//GEOMETRY geometry;
	//geometry.test();
	//clean_up(0);
	
	printf("hello, world! \n");
        PETSc solver;
        solver.SetTol(10e-9);
        EBM2D_LAPLACE mylaplace;
        EBM2D_CARTESIAN ebm_cartesian(&front, &solver, &mylaplace);
        
        ebm_cartesian.solve();
        ebm_cartesian.saveInterface_Tecplot("intfc.plt");
        ebm_cartesian.saveStates_Tecplot();
	ebm_cartesian.saveStates_VTK();
	
	//ebm_cartesian.saveDerivatives_Tecplot();
	//ebm_cartesian.saveReconstructedIntfc_Tecplot("reconstructed.plt",-1,-1);
        //ebm_cartesian.saveInterfaceEdges_Tecplot("edges.plt");
	//ebm_cartesian.saveComponent_Tecplot();
	
	
	
        //delete [] laplace.solute;
        
        
        
        
        clean_up(0);
        MPI_Finalize();
        printf("bye bye! \n");
        return 0;
}

LOCAL double circle_curve(
        POINTER func_params,
        double *coords)
{

        double dist, theta;

	//return 1.0;
        //dist =   sqrt(sqr(coords[0]-0.5) + sqr(coords[1]-0.5));
        dist =   sqrt(sqr(coords[0]-0.0) + sqr(coords[1]-0.0));
        
        //return dist - sqrt(0.1);
        return dist - 0.5;
}       /* end circle_curve */

static  void print_vtk_solution(
	const char *dirname,
        POINTER wave,
        Front *front)
{
	LAPLACE   *laplace = (LAPLACE*)wave;
        INTERFACE *grid_intfc = front->grid_intfc;
        RECT_GRID *grid;
        int i,j,ic,*gmax,icoords[MAXD];
        double *h,*L,coords[MAXD];
        double  s;
        Table *T;
        COMPONENT comp,*gr_comp;
        FILE *vtkfile;
	double *solute = laplace->solute;

        char filenm[100], buff[100];
        grid = &topological_grid(grid_intfc);
        T = table_of_interface(grid_intfc);
        gmax = grid->gmax;
        h = grid->h;
        L = grid->L;
	if(pp_numnodes() > 1)
	{
	    sprintf(filenm,"%s/solute-%s_nd-",dirname,
	    		right_flush(front->step,7),right_flush(pp_mynode(),4));
	    sprintf(buff,"%s.vtk",right_flush(pp_mynode(),4));
	    strcat(filenm, buff);
	}    
	else    
            sprintf(filenm,"%s/solute-%s.vtk",dirname,
	    		right_flush(front->step,7)); 
        vtkfile=fopen(filenm,"w");
        (void) fprintf(vtkfile,"# vtk DataFile Version 2.0\n");
        (void) fprintf(vtkfile,"comment line\n");
        (void) fprintf(vtkfile,"ASCII\n");
        (void) fprintf(vtkfile,"DATASET STRUCTURED_POINTS\n");
        (void) fprintf(vtkfile,"DIMENSIONS %d %d 1\n",gmax[0]-1,gmax[1]-1);
        (void) fprintf(vtkfile,"SPACING %f %f 0\n",h[0],h[1]);
        (void) fprintf(vtkfile,"ORIGIN %f %f 0\n",L[0],L[1]);
        (void) fprintf(vtkfile,"POINT_DATA %d\n",(gmax[0]-1)*(gmax[1]-1));
        (void) fprintf(vtkfile,"SCALARS solute float 1 \n");
        (void) fprintf(vtkfile,"LOOKUP_TABLE default \n");

        gr_comp = T->components; 

	// printing the vtk data
        for (j = 1; j < gmax[1]; ++j)
        {
	    for (i = 1; i < gmax[0]; ++i)
            {
                icoords[0] = i;
                icoords[1] = j;
                ic = d_index2d(i,j,gmax);
                comp = gr_comp[ic];
                s = solute[ic];
                if (comp == 1)
                {
                    s = 0.0;
                }   
                fprintf(vtkfile,"%f ",s);
            }
            fprintf(vtkfile," \n");
        }
 
	if (pp_mynode() == 1 && pp_numnodes() > 1) 
	//Only add file names to visit files via the first node.	
	{
	    for(i = 0; i < pp_numnodes(); i++)
	    {
	        char name[100];
		sprintf(name, "data%s/solute-%s_nd-",
				right_flush(front->step,7), 
				right_flush(front->step,7));
		sprintf(buff,"%s.vtk\n",right_flush(i,4)); 
	        strcat(name,buff);
	    }
        }
	fclose(vtkfile);
} 

static  void print_error_analysis(
	const char *dirname,
        POINTER wave,
        Front *front)
{
/*	LAPLACE   *laplace = (LAPLACE*)wave;
        INTERFACE *grid_intfc = front->grid_intfc;
        RECT_GRID *grid,*comp_grid;
        int i,j,ic,*gmax;
        double  s,exact_solu;
        Table *T;
        COMPONENT comp,*gr_comp;
        FILE *err_file;
        char filenm[100], buff[100];
	double L_1,L_2,L_I;
	double *h,*L,coords[MAXD];
	double *l_1,*l_2,*l_i,NumOfPts;
	int *num_mesh_pts,mn;
	int num_nodes;
	int xmax,ymax;
	int ix,iy;		// indices of expanded grid 
	double *solute = laplace->solute;
	double (*solution)(POINTER,double*) = laplace->solution;
	POINTER sparams = laplace->sparams;

	num_nodes = pp_numnodes();
	vector(&l_1,num_nodes,sizeof(double));
	vector(&l_2,num_nodes,sizeof(double));
	vector(&l_i,num_nodes,sizeof(double));
	vector(&num_mesh_pts,num_nodes,sizeof(int));
	for (i = 0; i < num_nodes; ++i)
	{
	    l_1[i] = l_2[i] = l_i[i] = 0.0;
	    num_mesh_pts[i] = 0;
	}

	comp_grid = front->rect_grid;
	xmax = comp_grid->gmax[0];
	ymax = comp_grid->gmax[1];
	L = comp_grid->L;
	h = comp_grid->h;

        grid = &topological_grid(grid_intfc);
        T = table_of_interface(grid_intfc);
        gmax = grid->gmax;
	if (num_nodes > 1)
	{
	    sprintf(filenm,"%s/solute-%s_nd-",dirname,
	    		right_flush(front->step,7),right_flush(pp_mynode(),4));
	    sprintf(buff,"%s.err",right_flush(pp_mynode(),4));
	    strcat(filenm, buff);
	}    
	else    
            sprintf(filenm,"%s/solute-%s.err",dirname,
	    		right_flush(front->step,7)); 
        err_file=fopen(filenm,"w");

        gr_comp = T->components; 

	fprintf(err_file,"Computational grid: %d X %d\n",
			front->rect_grid->gmax[0],
			front->rect_grid->gmax[1]);

	mn = pp_mynode();
	for (i = 0; i < xmax; ++i)
        {
            for (j = 0; j < ymax; ++j)
            {
		coords[0] = L[0] + (i+0.5)*h[0];
                coords[1] = L[1] + (j+0.5)*h[1];
		if (comp_grid->lbuf[0] == 0)
		    ix = i + 1;
		else
		    ix = i + comp_grid->lbuf[0];
		if (comp_grid->lbuf[1] == 0)
		    iy = j + 1;
		else
		    iy = j + comp_grid->lbuf[1];
                ic = d_index2d(ix,iy,gmax);
                s = solute[ic];
		exact_solu = solution(sparams,coords);
                fprintf(err_file,"(%12.8f %12.8f):  %f   %f\n",
				coords[0],coords[1],s,exact_solu);
	    	
	    	l_1[mn] += fabs(s - exact_solu);
	    	l_2[mn] += sqr(s - exact_solu);
		if (l_i[mn] < fabs(s - exact_solu))
		    l_i[mn] = fabs(s - exact_solu);
	    	num_mesh_pts[mn]++;
            }
            fprintf(err_file," \n");
        }
	if (num_nodes > 1)
	{
	    pp_global_max(l_1,num_nodes);
	    pp_global_max(l_2,num_nodes);
	    pp_global_max(l_i,num_nodes);
	    pp_global_imax(num_mesh_pts,num_nodes);
	}
	if (mn == 0)
	{
	    L_1 = L_2 = L_I = NumOfPts = 0.0;
	    for (i = 0; i < num_nodes; ++i)
	    {
		L_1 += l_1[i];
		L_2 += l_2[i];
		NumOfPts += (double)num_mesh_pts[i];
		if (L_I < l_i[i]) L_I = l_i[i];
	    }
	    L_1 /= NumOfPts;
	    L_2 = sqrt(L_2/NumOfPts);
            fprintf(err_file,"The L-1 norm of error: %g\n",L_1);
            fprintf(err_file,"The L-2 norm of error: %g\n",L_2);
            fprintf(err_file,"The L-I norm of error: %g\n",L_I);
	}
 
	fclose(err_file); */
} 

/*****************************************************************
 *  	The following is a set of test functions. Each set       *
 *  	consists of three functions: solution phi, source rho    *
 *  	and flux F and we have assumed:                          *
 *  	               \Delta phi = rho                          *
 *  	and                                                      *
 *  	               rho = div (F)                             *
 *****************************************************************/

LOCAL double exp_func(
        POINTER func_params,
        double *coords)
{
	EXP_PARAMS *exp_params = (EXP_PARAMS*)func_params;
	double *K = exp_params->K;
	double x,y,zz;

	x = coords[0];
	y = coords[1];
	zz = exp(K[0]*x*x+K[1]*y*y)*(4.0*K[0]*K[0]*x*x+
			4.0*K[1]*K[1]*y*y+2.0*(K[0]+K[1]));

	return zz;
}	/* end exp_func */

LOCAL void exp_func_flux(
        POINTER func_params,
        double *coords,
	double *flux)
{
	EXP_PARAMS *exp_params = (EXP_PARAMS*)func_params;
	double *K = exp_params->K;
	double x,y,zz;
	int i;

	printf("Entering exp_func_flux()\n");
	x = coords[0];
	y = coords[1];
	for (i = 0; i < 2; ++i)
	    flux[i] = 2.0*K[i]*coords[i]*exp(K[0]*coords[0]*coords[0]
	    		+ K[1]*coords[1]*coords[1]);
}	/* end exp_func_flux */

LOCAL double exp_solution(
	POINTER func_params,
        double *coords)
{
	EXP_PARAMS *exp_params = (EXP_PARAMS*)func_params;
	double *K = exp_params->K;
	double x,y,zz;

	x = coords[0];
	y = coords[1];
	zz = exp(K[0]*x*x+K[1]*y*y);
	return zz;
}	/* end exp_solution */

LOCAL double sin_func(
        POINTER func_params,
        double *coords)
{
	double x,y,zz;

	x = coords[0];
	y = coords[1];
	zz = cos(2.0*PI*x) + cos(2.0*PI*y);

	return zz;
}	/* end sin_func */

LOCAL void sin_func_flux(
        POINTER func_params,
        double *coords,
	double *flux)
{
	flux[0] = sin(2.0*PI*coords[0])/2.0/PI;
	flux[1] = sin(2.0*PI*coords[1])/2.0/PI;
}	/* end sin_func_flux */

LOCAL double sin_solution(
        POINTER func_params,
        double *coords)
{
	double x,y,zz;

	x = coords[0];
	y = coords[1];
	zz = -(cos(2.0*PI*x) + cos(2.0*PI*y))/4.0/PI/PI;

	return zz;
}	/* end sin_solution */

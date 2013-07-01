/****************************************************************
FronTier is a set of libraries that implements differnt types of 
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
*****************************************************************/


/*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*/

#include <iFluid.h>
#include <airfoil.h>

/*  Local Function Declarations */
static void spring_driver(Front*);
static void initSpringPropagation(Front*);
static void initCurvePropagation(Front*);
static void initNodePropagation(Front*);
static void spring_surface_propagate(Front*,POINTER,SURFACE*,SURFACE*,double);
static void spring_curve_propagate(Front*,POINTER,CURVE*,CURVE*,double);
static void spring_node_propagate(Front*,POINTER,NODE*,NODE*,double);
static void gviewSurfaceStrain(Front*);
static void gviewSurfaceStress(Front*);
static void naturalStressOfTri(TRI*,double);
static void vtkPlotSurfaceStress(Front*);

char *in_name,*restart_state_name,*restart_name,*out_name;
boolean RestartRun;
boolean ReSetTime;
int RestartStep;
boolean binary = YES;
int constrained_propagate;

int main(int argc, char **argv)
{
	static Front front;
	static F_BASIC_DATA f_basic;
	static LEVEL_FUNC_PACK level_func_pack;
	AF_PARAMS 	af_params;

	FT_Init(argc,argv,&f_basic);
	f_basic.size_of_intfc_state = sizeof(STATE);

	//Initialize Petsc before FrontStartUP
        PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);

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

        front.extra2 = (POINTER)&af_params;

	level_func_pack.pos_component = LIQUID_COMP2;
	if (!RestartRun)
	{
	    FT_InitIntfc(&front,&level_func_pack);
	    if (f_basic.dim < 3)
            	FT_ClipIntfcToSubdomain(&front);
	    else
		initSpringModel(&front);
	}
	setMotionParams(&front);

	/* Time control */
	FT_ReadTimeControl(in_name,&front);

	if (!RestartRun)
	{
	    optimizeElasticMesh(&front);
	    set_equilibrium_mesh(&front);
	}
	FT_SetGlobalIndex(&front);
	initSpringPropagation(&front);
	    
        if (RestartRun)
	{
	    readAfExtraDada(&front,restart_state_name);
	}

	if (!RestartRun || ReSetTime)
	    resetFrontVelocity(&front);

	/* Propagate the front */

	spring_driver(&front);

	clean_up(0);
}

static  void spring_driver(
        Front *front)
{
        double CFL;
        int  dim = front->rect_grid->dim;
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;

        CFL = Time_step_factor(front);
	Tracking_algorithm(front) = SIMPLE_TRACKING;

	if (!RestartRun || ReSetTime)
	{
	    FT_ResetTime(front);

	    // Always output the initial interface.
	    FT_Save(front,out_name);

	    printAfExtraDada(front,out_name);

            FT_AddMovieFrame(front,out_name,binary);

	    FT_Propagate(front);
	    print_airfoil_stat(front,out_name);

            FT_SetOutputCounter(front);
	    FT_SetTimeStep(front);
	}
	else
	{
	    FT_SetOutputCounter(front);
	}
	FT_TimeControlFilter(front);
        (void) printf("\ntime = %20.14f   step = %5d   next dt = %20.14f\n",
                        front->time,front->step,front->dt);
	
        for (;;)
        {
	    /* Propagating interface for time step dt */

            FT_Propagate(front);

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

	    /* Output section */

	    print_airfoil_stat(front,out_name);

            if (FT_IsSaveTime(front))
	    {
		FT_Save(front,out_name);
	    	printAfExtraDada(front,out_name);
		gviewSurfaceStrain(front);
		gviewSurfaceStress(front);
	    }
            if (FT_IsMovieFrameTime(front))
	    {
                FT_AddMovieFrame(front,out_name,binary);
	    }

            if (FT_TimeLimitReached(front))
	    {
            	(void) printf("\ntime = %20.14f   step = %5d   ",
				front->time,front->step);
		(void) printf("next dt = %20.14f\n",front->dt);
                break;
	    }

	    /* Time and step control section */

	    FT_TimeControlFilter(front);
	    print_storage("trace","after time loop");

            (void) printf("\ntime = %20.14f   step = %5d   next dt = %20.14f\n",
                        front->time,front->step,front->dt);
            fflush(stdout);
        }
        (void) delete_interface(front->interf);
}       /* end spring_driver */

static void initSpringPropagation(
	Front *front)
{
	Tracking_algorithm(front) = SIMPLE_TRACKING;
	front->_surface_propagate = spring_surface_propagate;
	front->_curve_propagate = spring_curve_propagate;
	initCurvePropagation(front);
	initNodePropagation(front);

	front->interior_propagate = fourth_order_elastic_surf_propagate;;
}	/* end initSpringPropagation */

struct _NODE_VEL_PARAMS {
	double dir[MAXD];
        double v0;
	double *time;
        double stop_time;
};
typedef struct _NODE_VEL_PARAMS NODE_VEL_PARAMS;

struct _CURVE_VEL_PARAMS {
	double dir[MAXD];
        double v0;
	double *time;
        double stop_time;
};
typedef struct _CURVE_VEL_PARAMS CURVE_VEL_PARAMS;

static void node_vel_func(
	POINTER vparams,
	double *vel)
{
	NODE_VEL_PARAMS *nvparams = (NODE_VEL_PARAMS*)vparams;	
	int i;
	for (i = 0; i < 3; ++i)	
	    vel[i] = 0.0;
	if (*(nvparams->time) >= nvparams->stop_time)
	    return;
	for (i = 0; i < 3; ++i)	
	{
	    vel[i] = nvparams->v0*nvparams->dir[i];
	}
}	/* end node_vel_func */

static void initCurvePropagation(
	Front *front)
{
	static CURVE_VEL_PARAMS *vparams;
	CURVE *curves[20];
	INTERFACE *intfc = front->interf;
	SURFACE **s,*surf;
	CURVE **c;
	char string[200];
	int i,num_curves = 0;
	FILE *infile;
	static AF_NODE_EXTRA node_extra;

	infile = fopen(InName(front),"r");
	if (!CursorAfterStringOpt(infile,"Enter yes to set curve motion:"))
	{
	    fclose(infile);
	    return;
	}
	else
	{
	    fscanf(infile,"%s",string);
	    (void) printf("%s\n",string);
	    if (string[0] != 'y' && string[0] != 'Y')
	    {
	    	fclose(infile);
		return;
	    }
	}

	front->_curve_propagate = spring_curve_propagate;
	node_extra.af_node_type = PRESET_NODE;
	intfc_surface_loop(intfc,s)
	{
	    if (wave_type(*s) == ELASTIC_BOUNDARY)
	    {
		surf = *s;
		surf_pos_curve_loop(surf,c)
		{
		    if (!pointer_in_list((POINTER)(*c),num_curves,
					(POINTER*)curves))
			curves[num_curves++] = *c;
		}
		surf_neg_curve_loop(surf,c)
		{
		    if (!pointer_in_list((POINTER)(*c),num_curves,
					(POINTER*)curves))
			curves[num_curves++] = *c;
		}
	    }
	}
	FT_VectorMemoryAlloc((POINTER*)&vparams,num_curves,
			sizeof(CURVE_VEL_PARAMS));
	for (i = 0; i < num_curves; ++i)
	{
	    (void) printf("For curve from (%f %f %f) to (%f %f %f)\n",
			Coords(curves[i]->first->start)[0],
			Coords(curves[i]->first->start)[1],
			Coords(curves[i]->first->start)[2],
			Coords(curves[i]->last->end)[0],
			Coords(curves[i]->last->end)[1],
			Coords(curves[i]->last->end)[2]);
	    CursorAfterString(infile,"Type yes to set the curve motion:");
	    fscanf(infile,"%s",string);
	    (void) printf("%s\n",string);
	    if (string[0] != 'y' && string[0] != 'Y')
		continue;
	    CursorAfterString(infile,"Enter velocity direction:");
	    fscanf(infile,"%lf %lf %lf",&vparams[i].dir[0],&vparams[i].dir[1],
				&vparams[i].dir[2]);
	    (void) printf("%f %f %f\n",vparams[i].dir[0],vparams[i].dir[1],
				vparams[i].dir[2]);
	    CursorAfterString(infile,"Enter speed:");
	    fscanf(infile,"%lf",&vparams[i].v0);
	    (void) printf("%f\n",vparams[i].v0);
	    CursorAfterString(infile,"Enter stop time:");
	    fscanf(infile,"%lf",&vparams[i].stop_time);
	    (void) printf("%f\n",vparams[i].stop_time);
	    curves[i]->vparams = (POINTER)&vparams[i];
	    curves[i]->vfunc = node_vel_func;
	    vparams[i].time = &front->time;
	    hsbdry_type(curves[i]) = PRESET_CURVE;
	    curves[i]->start->extra = &node_extra;
	    curves[i]->end->extra = &node_extra;
	}
}	/* end initCurvePropagation */

static void initNodePropagation(
	Front *front)
{
	static NODE_VEL_PARAMS *vparams;
	static AF_NODE_EXTRA node_extra;
	NODE **n,*nodes[20];
	INTERFACE *intfc = front->interf;
	SURFACE **s,*surf;
	CURVE **c;
	int i,num_nodes = 0;
	FILE *infile = fopen(InName(front),"r");

	intfc_node_loop(intfc,n)
	    sorted((*n)->posn) = NO;
	intfc_curve_loop(intfc,c)
	{
	    if (hsbdry_type(*c) == PRESET_CURVE) 
		sorted((*c)->start->posn) = sorted((*c)->end->posn) = YES;
	}
	intfc_surface_loop(intfc,s)
	{
	    if (wave_type(*s) == ELASTIC_BOUNDARY)
	    {
		surf = *s;
		surf_pos_curve_loop(surf,c)
		{
		    if (!pointer_in_list((POINTER)(*c)->start,num_nodes,
				(POINTER*)nodes) && !sorted((*c)->start->posn))
			nodes[num_nodes++] = (*c)->start;
		    if (!pointer_in_list((POINTER)(*c)->end,num_nodes,
				(POINTER*)nodes) && !sorted((*c)->end->posn))
			nodes[num_nodes++] = (*c)->end;
		}
		surf_neg_curve_loop(surf,c)
		{
		    if (!pointer_in_list((POINTER)(*c)->start,num_nodes,
				(POINTER*)nodes) && !sorted((*c)->start->posn))
			nodes[num_nodes++] = (*c)->start;
		    if (!pointer_in_list((POINTER)(*c)->end,num_nodes,
				(POINTER*)nodes) && !sorted((*c)->end->posn))
			nodes[num_nodes++] = (*c)->end;
		}
	    }
	}
	if (num_nodes == 0)
	{
	    fclose(infile);
	    return;
	}
	front->_node_propagate = spring_node_propagate;
	FT_VectorMemoryAlloc((POINTER*)&vparams,num_nodes,
			sizeof(NODE_VEL_PARAMS));
	node_extra.af_node_type = PRESET_NODE;
	for (i = 0; i < num_nodes; ++i)
	{
	    (void) printf("For node at (%f %f %f)\n",Coords(nodes[i]->posn)[0],
			Coords(nodes[i]->posn)[1],Coords(nodes[i]->posn)[2]);
	    CursorAfterString(infile,"Enter velocity direction:");
	    fscanf(infile,"%lf %lf %lf",&vparams[i].dir[0],&vparams[i].dir[1],
				&vparams[i].dir[2]);
	    (void) printf("%f %f %f\n",vparams[i].dir[0],vparams[i].dir[1],
				vparams[i].dir[2]);
	    CursorAfterString(infile,"Enter speed:");
	    fscanf(infile,"%lf",&vparams[i].v0);
	    (void) printf("%f\n",vparams[i].v0);
	    CursorAfterString(infile,"Enter stop time:");
	    fscanf(infile,"%lf",&vparams[i].stop_time);
	    (void) printf("%f\n",vparams[i].stop_time);
	    nodes[i]->vparams = (POINTER)&vparams[i];
	    nodes[i]->vfunc = node_vel_func;
	    nodes[i]->extra = &node_extra;
	    vparams[i].time = &front->time;
	}
}	/* end initNodePropagation */

static void spring_surface_propagate(
	Front *front,
	POINTER wave,
	SURFACE *olds,
	SURFACE *news,
	double dt)
{
}	/* end spring_surface_propagate */

static void spring_curve_propagate(
	Front *front,
	POINTER wave,
	CURVE *oldc,
	CURVE *newc,
	double dt)
{
	double vel[MAXD],s;
	int i,dim;
	POINT *oldp,*newp;
	BOND *oldb,*newb;

	if (debugging("trace"))
	    (void) printf("Entering spring_curve_propagate()\n");
	if (hsbdry_type(oldc) != PRESET_CURVE || oldc->vfunc == NULL) 
	{
	    if (debugging("trace"))
	    	(void) printf("Leaving spring_curve_propagate()\n");
	    return;
	}

	dim = FT_Dimension();
	(*oldc->vfunc)(oldc->vparams,vel);

	oldb = oldc->first;
	newb = newc->first;
	oldp = oldb->start;
	newp = newb->start;
	for (i = 0; i < dim; ++i)
	{
	    Coords(newp)[i] = Coords(oldp)[i] + dt*vel[i];
	    set_max_front_speed(i,fabs(vel[i]),NULL,Coords(newp),front);
	}
	s = mag_vector(vel,dim);
        set_max_front_speed(dim,s,NULL,Coords(newp),front);
	for (oldb = oldc->first, newb = newc->first; oldb != NULL;
		oldb = oldb->next, newb = newb->next)
	{
	    oldp = oldb->end;	newp = newb->end;
	    for (i = 0; i < dim; ++i)
	    {
	    	Coords(newp)[i] = Coords(oldp)[i] + dt*vel[i];
	    	set_max_front_speed(i,fabs(vel[i]),NULL,Coords(newp),front);
	    }
	    s = mag_vector(vel,dim);
            set_max_front_speed(dim,s,NULL,Coords(newp),front);
	}
	if (debugging("trace"))
	    (void) printf("Leaving spring_curve_propagate()\n");
}	/* end spring_curve_propagate */

static void spring_node_propagate(
	Front *front,
	POINTER wave,
	NODE *oldn,
	NODE *newn,
	double dt)
{
	double vel[MAXD],s;
	int i,dim;

	if (oldn->vfunc == NULL) return;
	dim = FT_Dimension();
	(*oldn->vfunc)(oldn->vparams,vel);
	for (i = 0; i < dim; ++i)
	{
	    Coords(newn->posn)[i] = Coords(oldn->posn)[i] + dt*vel[i];
	    set_max_front_speed(i,fabs(vel[i]),NULL,Coords(newn->posn),front);
	}
	s = mag_vector(vel,dim);
        set_max_front_speed(dim,s,NULL,Coords(newn->posn),front);
}	/* end spring_node_propagate */

static void gviewSurfaceStrain(
	Front *front)
{
}	/* end gviewSurfaceStrain */

static void gviewSurfaceStress(
	Front *front)
{
	char *outname = OutName(front);
	INTERFACE *intfc = front->interf;
	SURFACE **s;
	TRI *tri;
	POINT *p;
	double f[MAXD];
	char dirname[200];
	int i,j;
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	double ks = af_params->ks;
	double max_color = -HUGE;
	double min_color = HUGE;
	int n,N;
	double *color;
	
	n = 0;
	intfc_surface_loop(intfc,s)
	{
	    if (Boundary(*s)) continue;
	    surf_tri_loop(*s,tri)
	    {
		n++;
		naturalStressOfTri(tri,ks);
		if (max_color < tri->color)
		    max_color = tri->color;
		if (min_color > tri->color)
		    min_color = tri->color;
	    }
	}
	FT_VectorMemoryAlloc((POINTER*)&color,n,sizeof(double));
	printf("min_color = %f  max_color = %f\n",min_color,max_color);
	n = 0;
	intfc_surface_loop(intfc,s)
	{
	    if (Boundary(*s)) continue;
	    surf_tri_loop(*s,tri)
		tri->color = log(tri->color-min_color+1);
	}
	/* Smoothing loop */
	intfc_surface_loop(intfc,s)
	{
	    if (Boundary(*s)) continue;
	    I_SmoothSurfColor(*s,3);
	}
	sprintf(dirname,"%s/%s-ts%s",outname,"gv.stress",
			right_flush(front->step,7));
	gview_plot_color_scaled_interface(dirname,intfc);
	vtkPlotSurfaceStress(front);
}	/* end gviewSurfaceStress */

static void naturalStressOfTri(
	TRI *tri,
	double ks)
{
	double tau[3];
	double sigma[3];
	int i,j;
	double len,len0;
	double vec[3];
	double s[3],c[3];
	double b1,b2,arg,sigma1,sigma2;

	for (i = 0; i < 3; ++i)
	{
	    len0 = tri->side_length0[i];
	    for (j = 0; j < 3; ++j)
	    {
		vec[j] = Coords(Point_of_tri(tri)[(i+1)%3])[j] -
			Coords(Point_of_tri(tri)[i])[j];
	    }
	    len = Mag3d(vec);
	    tau[i] = ks*(len - len0);
	    c[i] = vec[0]/len;
	    s[i] = vec[1]/len;
	}
	// Convert to Cartesian tensor
	for (i = 0; i < 3; ++i)
	{
	    sigma[i] = sqr(c[i])*tau[0] + sqr(s[i])*tau[1] + s[i]*c[i]*tau[2];
	}
	// Diagonalize the stress tensor for principal directions
	b1 = -(sigma[0] + sigma[1]);
	b2 = sigma[0]*sigma[1] - 0.25*sqr(sigma[2]);
	arg = sqr(b1) - 4.0*b2;
	sigma1 = 0.5*(-b1 + sqrt(arg));
	sigma2 = 0.5*(-b1 - sqrt(arg));
	// Use von Mises stress as a measure
	tri->color = sqrt(sqr(sigma1) + sqr(sigma2) - sigma1*sigma2);
}	/* end naturalStressOfTri */

static void vtkPlotSurfaceStress(
	Front *front)
{
	char *outname = OutName(front);
	INTERFACE *intfc = front->interf;
	SURFACE **s;
	TRI *tri;
	POINT *p;
	char dirname[200],fname[200];
	int i,j;
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	int n,N;
	double *color;
	FILE *vfile;
	int num_tri;
	
	n = 0;
	sprintf(dirname,"%s/%s%s",outname,"vtk.ts",
			right_flush(front->step,7));
	if (!create_directory(dirname,NO))
        {
            printf("Cannot create directory %s\n",dirname);
            clean_up(ERROR);
        }
	sprintf(fname,"%s/%s",dirname,"stress.vtk");
	vfile = fopen(fname,"w");
	fprintf(vfile,"# vtk DataFile Version 3.0\n");
        fprintf(vfile,"Surface stress\n");
        fprintf(vfile,"ASCII\n");
        fprintf(vfile,"DATASET UNSTRUCTURED_GRID\n");

	num_tri = 0;

	intfc_surface_loop(intfc,s)
	{
	    if (Boundary(*s)) continue;
	    num_tri += (*s)->num_tri;
	}
	fprintf(vfile,"POINTS %d double\n", 3*num_tri);
	intfc_surface_loop(intfc,s)
	{
	    if (Boundary(*s)) continue;
	    surf_tri_loop(*s,tri)
	    {
		for (i = 0; i < 3; ++i)
		{
		    p = Point_of_tri(tri)[i];
		    fprintf(vfile,"%f %f %f\n",Coords(p)[0],Coords(p)[1],
						Coords(p)[2]);
		}
	    }
	}
	fprintf(vfile,"CELLS %i %i\n",num_tri,4*num_tri);
	n = 0;
	intfc_surface_loop(intfc,s)
	{
	    if (Boundary(*s)) continue;
	    surf_tri_loop(*s,tri)
	    {
		fprintf(vfile,"3 %i %i %i\n",3*n,3*n+1,3*n+2);
		n++;
	    }
	}
	fprintf(vfile, "CELL_TYPES %i\n",num_tri);
        intfc_surface_loop(intfc,s)
	{
            if (Boundary(*s)) continue;
            surf_tri_loop(*s,tri)
            {
                fprintf(vfile,"5\n");
            }
        }

	fprintf(vfile, "CELL_DATA %i\n", num_tri);
        fprintf(vfile, "SCALARS von_Mises_stress double\n");
        fprintf(vfile, "LOOKUP_TABLE default\n");
	intfc_surface_loop(intfc,s)
	{
	    if (Boundary(*s)) continue;
	    surf_tri_loop(*s,tri)
	    {
		fprintf(vfile,"%f\n",tri->color);
	    }
	}
}	/* end vtkSurfaceStress */


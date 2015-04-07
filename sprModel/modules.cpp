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
#include <airfoil.h>

static void initSpringModel2d(Front*);
static void initSpringModel3d(Front*);
static void cutToRectangle(FILE*,SURFACE*);
static void cutToTriangle(FILE*,SURFACE*);
static void cutToEllipse(FILE*,SURFACE*);
static void cutToCross(FILE*,SURFACE*);
static void cutToWing(FILE*,SURFACE*);
static void findClosestBdryPoint(double*,POINT**,CURVE**,SURFACE*);

extern void initSpringModel(Front *front)
{
	int dim = front->rect_grid->dim;

	if (debugging("trace"))
	    (void) printf("Entering initSpringModel()\n");
	switch (dim)
	{
	case 2:
	    initSpringModel2d(front);
	    break;
	case 3:
	    initSpringModel3d(front);
	    break;
	default:
	    (void) printf("Unknown dimension %d\n",dim);
	    clean_up(ERROR);
	}

	if (debugging("set_module"))
	    gview_plot_interface("module-step-1",front->interf);

	if (debugging("trace"))
	    (void) printf("Leaving initSpringModel()\n");
}	/* end initSpringModel */

static void initSpringModel2d(Front *front)
{
	int i,num_canopy;
	FILE *infile = fopen(InName(front),"r");

	if (debugging("trace"))
	    (void) printf("Entering initSpringModel2d()\n");

	if (debugging("set_module"))
	    gview_plot_interface("module-step-1",front->interf);

	if (debugging("trace"))
	    (void) printf("Leaving initSpringModel2d()\n");
	clean_up(0);
}	/* end initSpringModel */

static void initSpringModel3d(Front *front)
{
	int i;
	FILE *infile = fopen(InName(front),"r");
	double plane_nor[MAXD],plane_pt[MAXD];
        double height;
        double *L = front->rect_grid->L;
        double *U = front->rect_grid->U;
        COMPONENT amb_comp = front->interf->default_comp;
	SURFACE *surf;
	CURVE **c;
	char string[200];
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;

	if (debugging("trace"))
	    (void) printf("Entering initSpringModel3d()\n");

	CursorAfterString(infile,"Enter the height of the plane:");
        fscanf(infile,"%lf",&height);
        (void) printf("%f\n",height);

        for (i = 0; i < 2; ++i)
        {
            plane_nor[i] = 0.0;
            plane_pt[i] = 0.5*(L[i] + U[i]);
        }
        plane_nor[2] = 1.0;
        plane_pt[2] = height;
        FT_MakePlaneSurf(front,plane_nor,plane_pt,NO,amb_comp+1,amb_comp,
                        ELASTIC_BOUNDARY,&surf);
	negative_component(surf) = amb_comp;
	
	(void) printf("Available types of canopy boundaries are:\n");
	(void) printf("\tRectanglar (R)\n");
	(void) printf("\tTriangle (T)\n");
	(void) printf("\tElliptic (E)\n");
	(void) printf("\tCross (X)\n");
	(void) printf("\tWing (W)\n");
	CursorAfterString(infile,"Enter type of canopy boundary:");
        fscanf(infile,"%s",string);
        (void) printf("%s\n",string);
	switch (string[0])
	{
	case 'R':
	case 'r':
	    cutToRectangle(infile,surf);
	    break;
	case 'T':
	case 't':
	    cutToTriangle(infile,surf);
	    break;
	case 'E':
	case 'e':
	    cutToEllipse(infile,surf);
	    break;
	case 'X':
	case 'x':
	    cutToCross(infile,surf);
	    gview_plot_interface("test",surf->interface);
	    // Pending for debugging
	    clean_up(0);
	    break;
	case 'W':
	case 'w':
	    cutToWing(infile,surf);
	    break;
	}
	surf_pos_curve_loop(surf,c)
	{
	    if (is_closed_curve(*c))
	    {
		CursorAfterString(infile,
			"For closed curve, type yes if it is fixed:");
		fscanf(infile,"%s",string);
		(void) printf("%s\n",string);
		if (string[0] == 'Y' || string[0] == 'y')
		    hsbdry_type(*c) = FIXED_HSBDRY;
		else
		    hsbdry_type(*c) = MONO_COMP_HSBDRY;
	    }
	    else
	    {
		(void) printf("Curve from (%f %f) to (%f %f)\n",
				Coords((*c)->start->posn)[0],
				Coords((*c)->start->posn)[1],
				Coords((*c)->end->posn)[0],
				Coords((*c)->end->posn)[1]);
		CursorAfterString(infile,
			"For this curve, type yes if it is fixed:");
		fscanf(infile,"%s",string);
		(void) printf("%s\n",string);
		if (string[0] == 'Y' || string[0] == 'y')
		    hsbdry_type(*c) = FIXED_HSBDRY;
		else
		    hsbdry_type(*c) = MONO_COMP_HSBDRY;
	    }
	}
	surf_neg_curve_loop(surf,c)
	{
	    if (is_closed_curve(*c))
	    {
		CursorAfterString(infile,
			"For closed curve, type yes if it is fixed:");
		fscanf(infile,"%s",string);
		(void) printf("%s\n",string);
		if (string[0] == 'Y' || string[0] == 'y')
		    hsbdry_type(*c) = FIXED_HSBDRY;
		else
		    hsbdry_type(*c) = MONO_COMP_HSBDRY;
	    }
	    else
	    {
		(void) printf("Curve from (%f %f) to (%f %f)\n",
				Coords((*c)->start->posn)[0],
				Coords((*c)->start->posn)[1],
				Coords((*c)->end->posn)[0],
				Coords((*c)->end->posn)[1]);
		CursorAfterString(infile,
			"For this curve, type yes if it is fixed:");
		fscanf(infile,"%s",string);
		(void) printf("%s\n",string);
		if (string[0] == 'Y' || string[0] == 'y')
		    hsbdry_type(*c) = FIXED_HSBDRY;
		else
		    hsbdry_type(*c) = MONO_COMP_HSBDRY;
	    }
	}
	af_params->spring_model = MODEL1;
	if (CursorAfterStringOpt(infile,
            "Entering type of spring model: "))
        {
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            switch (string[0])
            {
            case '1':
                af_params->spring_model = MODEL1;
                break;
            case '2':
                af_params->spring_model = MODEL2;
                break;
            case '3':
                af_params->spring_model = MODEL3;
                break;
            default:
                break;
            }
        }
	fclose(infile);

	if (debugging("set_module"))
	    gview_plot_interface("module-step-1",front->interf);

	if (debugging("trace"))
	    (void) printf("Leaving initSpringModel3d()\n");
}	/* end initSpringModel */

static void cutToRectangle(
	FILE *infile,
	SURFACE *surf)
{
	static PLANE_PARAMS plane_params;
	double **insert_coords;
	double L[MAXD],U[MAXD];
	int i;
	double *N = plane_params.N;
	double *P = plane_params.P;
	POINT *pt;
	CURVE *curve;

	CursorAfterString(infile,"Enter lower bounds of the rectangle:");
	fscanf(infile,"%lf %lf",&L[0],&L[1]);
	(void) printf("%f %f\n",L[0],L[1]);
	CursorAfterString(infile,"Enter upper bounds of the rectangle:");
	fscanf(infile,"%lf %lf",&U[0],&U[1]);
	(void) printf("%f %f\n",U[0],U[1]);

	for (i = 0; i < 3; ++i)
	    P[i] = N[i] = 0.0;
	P[0] = L[0]; 	N[0] = 1.0;
	FT_CutSurfBdry(surf,plane_constr_func,(POINTER)&plane_params,NULL,0,2);

	for (i = 0; i < 3; ++i)
	    P[i] = N[i] = 0.0;
	P[1] = L[1]; 	N[1] = 1.0;
	FT_CutSurfBdry(surf,plane_constr_func,(POINTER)&plane_params,NULL,0,2);

	for (i = 0; i < 3; ++i)
	    P[i] = N[i] = 0.0;
	P[0] = U[0]; 	N[0] = -1.0;
	FT_CutSurfBdry(surf,plane_constr_func,(POINTER)&plane_params,NULL,0,2);

	for (i = 0; i < 3; ++i)
	    P[i] = N[i] = 0.0;
	P[1] = U[1]; 	N[1] = -1.0;
	FT_CutSurfBdry(surf,plane_constr_func,(POINTER)&plane_params,NULL,0,2);

	FT_InstallSurfEdge(surf,MONO_COMP_HSBDRY);

	FT_MatrixMemoryAlloc((POINTER*)&insert_coords,4,MAXD,sizeof(double));
	insert_coords[0][0] = L[0];
	insert_coords[0][1] = L[1];
	insert_coords[1][0] = L[0];
	insert_coords[1][1] = U[1];
	insert_coords[2][0] = U[0];
	insert_coords[2][1] = U[1];
	insert_coords[3][0] = U[0];
	insert_coords[3][1] = L[1];
	findClosestBdryPoint(insert_coords[0],&pt,&curve,surf);
	I_MoveNodeToPoint(pt,curve);
	for (i = 1; i < 4; ++i)
	{
	    findClosestBdryPoint(insert_coords[i],&pt,&curve,surf);
	    I_SplitCurve(pt,curve);
	}

	FT_FreeThese(1,insert_coords);
}	/* end cutToRectangle */

static void findClosestBdryPoint(
	double *coords,
	POINT **pt,
	CURVE **curve,
	SURFACE *surf)
{
	CURVE **c;
	BOND *b;
	POINT *p;
	double min_dist = HUGE;
	double dist;
	int dim = FT_Dimension();

	*pt = NULL;
	surf_pos_curve_loop(surf,c)
	{
	    curve_bond_loop(*c,b)
	    {
		p = b->start;
		dist = distance_between_positions(coords,Coords(p),dim);
		if (min_dist > dist)
		{
		    min_dist = dist;
		    *pt = p;
		    *curve = *c;
		}
		if (b->next == NULL)
		{
		    p = b->end;
		    dist = distance_between_positions(coords,Coords(p),dim);
		    if (min_dist > dist)
		    {
		    	min_dist = dist;
		    	*pt = p;
		    	*curve = *c;
		    }
		}
	    }
	}
	surf_neg_curve_loop(surf,c)
	{
	    curve_bond_loop(*c,b)
	    {
		p = b->start;
		dist = distance_between_positions(coords,Coords(p),dim);
		if (min_dist > dist)
		{
		    min_dist = dist;
		    *pt = p;
		    *curve = *c;
		}
		if (b->next == NULL)
		{
		    p = b->end;
		    dist = distance_between_positions(coords,Coords(p),dim);
		    if (min_dist > dist)
		    {
		    	min_dist = dist;
		    	*pt = p;
		    	*curve = *c;
		    }
		}
	    }
	}
}	/* end findClosestBdryPoint */

static void cutToEllipse(
	FILE *infile,
	SURFACE *surf)
{
}	/* end cutToEllipse */

static void cutToCross(
	FILE *infile,
	SURFACE *surf)
{
	static CROSS_CONSTR_PARAMS constr_params;
	static PLANE_PARAMS plane_params;
	double **insert_coords;
	double L[MAXD],U[MAXD];
	int i;
	double *N = plane_params.N;
	double *P = plane_params.P;

	(void) printf("Input the two crossing rectangles\n");
	CursorAfterString(infile,"Enter lower bounds of first rectangle:");
	fscanf(infile,"%lf %lf",&constr_params.L1[0],&constr_params.L1[1]);
	(void) printf("%f %f\n",constr_params.L1[0],constr_params.L1[1]);
	CursorAfterString(infile,"Enter upper bounds of first rectangle:");
	fscanf(infile,"%lf %lf",&constr_params.U1[0],&constr_params.U1[1]);
	(void) printf("%f %f\n",constr_params.U1[0],constr_params.U1[1]);
	CursorAfterString(infile,"Enter lower bounds of second rectangle:");
	fscanf(infile,"%lf %lf",&constr_params.L2[0],&constr_params.L2[1]);
	(void) printf("%f %f\n",constr_params.L2[0],constr_params.L2[1]);
	CursorAfterString(infile,"Enter upper bounds of second rectangle:");
	fscanf(infile,"%lf %lf",&constr_params.U2[0],&constr_params.U2[1]);
	(void) printf("%f %f\n",constr_params.U2[0],constr_params.U2[1]);
	FT_MatrixMemoryAlloc((POINTER*)&insert_coords,4,MAXD,sizeof(double));

	/* Cut to rectangle first to avoid corner problem */
	for (i = 0; i < 2; ++i)
	{
	    L[i] = (constr_params.L1[i] < constr_params.L2[i]) ? 
			constr_params.L1[i] : constr_params.L2[i];
	    U[i] = (constr_params.U1[i] > constr_params.U2[i]) ? 
			constr_params.U1[i] : constr_params.U2[i];
	}
	for (i = 0; i < 3; ++i)
	    P[i] = N[i] = 0.0;
	P[0] = L[0]; 	N[0] = 1.0;
	FT_CutSurfBdry(surf,plane_constr_func,(POINTER)&plane_params,NULL,0,2);

	for (i = 0; i < 3; ++i)
	    P[i] = N[i] = 0.0;
	P[1] = L[1]; 	N[1] = 1.0;
	FT_CutSurfBdry(surf,plane_constr_func,(POINTER)&plane_params,NULL,0,2);

	for (i = 0; i < 3; ++i)
	    P[i] = N[i] = 0.0;
	P[0] = U[0]; 	N[0] = -1.0;
	FT_CutSurfBdry(surf,plane_constr_func,(POINTER)&plane_params,NULL,0,2);

	for (i = 0; i < 3; ++i)
	    P[i] = N[i] = 0.0;
	P[1] = U[1]; 	N[1] = -1.0;
	FT_CutSurfBdry(surf,plane_constr_func,(POINTER)&plane_params,NULL,0,2);

	insert_coords[0][0] = FT_Max(constr_params.L1[0],constr_params.L2[0]);
	insert_coords[0][1] = FT_Max(constr_params.L1[1],constr_params.L2[1]);
	insert_coords[1][0] = FT_Min(constr_params.U1[0],constr_params.U2[0]);
	insert_coords[1][1] = FT_Min(constr_params.U1[1],constr_params.U2[1]);
	insert_coords[2][0] = FT_Max(constr_params.L1[0],constr_params.L2[0]);
	insert_coords[2][1] = FT_Min(constr_params.U1[1],constr_params.U2[1]);
	insert_coords[3][0] = FT_Min(constr_params.U1[0],constr_params.U2[0]);
	insert_coords[3][1] = FT_Max(constr_params.L1[1],constr_params.L2[1]);

	FT_CutSurfBdry(surf,xoss_constr_func,(POINTER)&constr_params,
			insert_coords,4,2);
	FT_InstallSurfEdge(surf,MONO_COMP_HSBDRY);
	FT_FreeThese(1,insert_coords);
}	/* end cutToCross */

static void cutToWing(
	FILE *infile,
	SURFACE *surf)
{
}	/* end cutToWing */

static void cutToTriangle(
	FILE *infile,
	SURFACE *surf)
{
	static PLANE_PARAMS plane_params;
	double **insert_coords;
	double L[MAXD],U[MAXD];
	int i;
	double *N = plane_params.N;
	double *P = plane_params.P;
	POINT *pt;
	CURVE *curve;

	FT_MatrixMemoryAlloc((POINTER*)&insert_coords,3,MAXD,sizeof(double));
	(void) printf("Triangle vertices on the x-y plane in counter-clock ");
	(void) printf("orientation\n");
	CursorAfterString(infile,"Enter x-y coordinates of vertex 1:");
	fscanf(infile,"%lf %lf",&insert_coords[0][0],&insert_coords[0][1]);
	(void) printf("%f %f\n",insert_coords[0][0],insert_coords[0][1]);
	CursorAfterString(infile,"Enter x-y coordinates of vertex 2:");
	fscanf(infile,"%lf %lf",&insert_coords[1][0],&insert_coords[1][1]);
	(void) printf("%f %f\n",insert_coords[1][0],insert_coords[1][1]);
	CursorAfterString(infile,"Enter x-y coordinates of vertex 3:");
	fscanf(infile,"%lf %lf",&insert_coords[2][0],&insert_coords[2][1]);
	(void) printf("%f %f\n",insert_coords[2][0],insert_coords[2][1]);

	for (i = 0; i < 3; ++i)
	    P[i] = N[i] = 0.0;
	P[0] = insert_coords[0][0];
	P[1] = insert_coords[0][1];
	N[0] = P[1] - insert_coords[1][1];
	N[1] = - P[0] + insert_coords[1][0];
	FT_CutSurfBdry(surf,plane_constr_func,(POINTER)&plane_params,NULL,0,2);

	for (i = 0; i < 3; ++i)
	    P[i] = N[i] = 0.0;
	P[0] = insert_coords[1][0];
	P[1] = insert_coords[1][1];
	N[0] = P[1] - insert_coords[2][1];
	N[1] = - P[0] + insert_coords[2][0];
	FT_CutSurfBdry(surf,plane_constr_func,(POINTER)&plane_params,NULL,0,2);

	for (i = 0; i < 3; ++i)
	    P[i] = N[i] = 0.0;
	P[0] = insert_coords[2][0];
	P[1] = insert_coords[2][1];
	N[0] = P[1] - insert_coords[0][1];
	N[1] = - P[0] + insert_coords[0][0];
	FT_CutSurfBdry(surf,plane_constr_func,(POINTER)&plane_params,NULL,0,2);

	FT_InstallSurfEdge(surf,MONO_COMP_HSBDRY);

	findClosestBdryPoint(insert_coords[0],&pt,&curve,surf);
	I_MoveNodeToPoint(pt,curve);
	for (i = 1; i < 3; ++i)
	{
	    findClosestBdryPoint(insert_coords[i],&pt,&curve,surf);
	    I_SplitCurve(pt,curve);
	}

	FT_FreeThese(1,insert_coords);
}	/* end cutToTriangle */


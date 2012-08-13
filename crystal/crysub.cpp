
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
*				crysub.c:
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

#include "crystal.h"
#include "crystal_basic.h"

/********************************************************************
 *	Level function parameters for the initial interface 	    *
 ********************************************************************/


/********************************************************************
 *	Velocity function parameters for the front	 	    *
 ********************************************************************/
#define		MAX_NUM_VERTEX_IN_CELL		20

	/*  Local Application Function Declarations */

static void	neumann_point_propagate(Front*,POINTER,POINT*,POINT*,	
			HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
static void	dirichlet_point_propagate(Front*,POINTER,POINT*,POINT*,	
			HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
static void	reaction_point_propagate(Front*,POINTER,POINT*,POINT*,	
			HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
static void 	euler_forward_scheme(CRT_PARAMS*,double,double,double,double,
				double,double*);
static void 	euler_backward_scheme(CRT_PARAMS*,double,double,double,double,
				double,double*);
static void 	middle_point_scheme(CRT_PARAMS*,double,double,double,double,
				double,double*);
static void 	constant_state(CRT_PARAMS*,double,double,double,double,
				double,double*);
static void 	setInitialIntfc1d(Front*,LEVEL_FUNC_PACK*,char*);
static void 	setInitialIntfc2d(Front*,LEVEL_FUNC_PACK*,char*);
static void 	setInitialIntfc3d(Front*,LEVEL_FUNC_PACK*,char*);

extern  double getStateSolute(
        POINTER state)
{
        STATE *solute_state = (STATE*)state;
        return solute_state->solute;
}       /* end getStateSolute */

extern  void crystal_point_propagate(
        Front *front,
        POINTER wave,
        POINT *oldp,
        POINT *newp,
        HYPER_SURF_ELEMENT *oldhse,
        HYPER_SURF         *oldhs,
        double              dt,
        double              *V)
{
	switch(wave_type(oldhs))
	{
	case SUBDOMAIN_BOUNDARY:
            return;
	case NEUMANN_BOUNDARY:
	    neumann_point_propagate(front,wave,oldp,newp,oldhse,
                                        oldhs,dt,V);
	    return;
	case DIRICHLET_BOUNDARY:
	    dirichlet_point_propagate(front,wave,oldp,newp,oldhse,
                                        oldhs,dt,V);
	    return;
	case GROWING_BODY_BOUNDARY:
	    reaction_point_propagate(front,wave,oldp,newp,oldhse,
                                        oldhs,dt,V);
	    return;
	}
}	/* end crystal_point_propagate */

static  void neumann_point_propagate(
        Front *front,
        POINTER wave,
        POINT *oldp,
        POINT *newp,
        HYPER_SURF_ELEMENT *oldhse,
        HYPER_SURF         *oldhs,
        double              dt,
        double              *V)
{
        int i, dim = front->rect_grid->dim;
	double nor[MAXD];
        double p1[MAXD],p2[MAXD];
        double *p0 = Coords(oldp);
        double dn,*h = front->rect_grid->h;
        double s0,s1,s2;
        STATE *sl,*sr,*state;
	CRT_PARAMS *cRparams = (CRT_PARAMS*)front->extra2;
        double *solute = cRparams->field->solute;

        for (i = 0; i < dim; ++i)
            Coords(newp)[i] = Coords(oldp)[i];
       	FT_GetStatesAtPoint(oldp,oldhse,oldhs,(POINTER*)&sl,(POINTER*)&sr);
       	state = (negative_component(oldhs) == SOLUTE_COMP) ? sl : sr;
	s0 = state->solute;
       	FT_NormalAtPoint(oldp,front,nor,SOLUTE_COMP);
       	dn = grid_size_in_direction(nor,h,dim);
       	for (i = 0; i < dim; ++i)
       	    p1[i] = p0[i] + nor[i]*dn;
	FT_IntrpStateVarAtCoords(front,SOLUTE_COMP,p1,solute,
			getStateSolute,&s1,&s0);
       	for (i = 0; i < dim; ++i)
       	    p2[i] = p1[i] + nor[i]*dn;
	FT_IntrpStateVarAtCoords(front,SOLUTE_COMP,p2,solute,
			getStateSolute,&s2,&s1);
	state = (negative_component(oldhs) == SOLUTE_COMP) ? 
		(STATE*)left_state(newp) : (STATE*)right_state(newp);
	state->solute = (4.0*s1 - s2)/3.0;
        return;
}       /* neumann_point_propagate */

static  void dirichlet_point_propagate(
        Front *front,
        POINTER wave,
        POINT *oldp,
        POINT *newp,
        HYPER_SURF_ELEMENT *oldhse,
        HYPER_SURF         *oldhs,
        double              dt,
        double              *V)
{
        int i, dim = front->rect_grid->dim;
        STATE *sl,*sr,*state,*bstate;
	CRT_PARAMS *cRparams = (CRT_PARAMS*)front->extra2;

        for (i = 0; i < dim; ++i)
            Coords(newp)[i] = Coords(oldp)[i];
	if (boundary_state(oldhs) != NULL)
	{
	    bstate = (STATE*)boundary_state(oldhs);
       	    FT_GetStatesAtPoint(oldp,oldhse,oldhs,(POINTER*)&sl,(POINTER*)&sr);
	    state =  (STATE*)left_state(newp);
	    state->solute = bstate->solute;
	    if (cRparams->max_solute < state->solute)
		cRparams->max_solute = state->solute;
	    if (cRparams->min_solute > state->solute)
		cRparams->min_solute = state->solute;
	    state =  (STATE*)right_state(newp);
	    state->solute = bstate->solute;
	    if (cRparams->max_solute < state->solute)
		cRparams->max_solute = state->solute;
	    if (cRparams->min_solute > state->solute)
		cRparams->min_solute = state->solute;
	}
	else
	{
       	    FT_GetStatesAtPoint(oldp,oldhse,oldhs,(POINTER*)&sl,(POINTER*)&sr);
	    state =  (STATE*)left_state(newp);
	    state->solute = sl->solute;
	    if (cRparams->max_solute < state->solute)
		cRparams->max_solute = state->solute;
	    if (cRparams->min_solute > state->solute)
		cRparams->min_solute = state->solute;
	    state =  (STATE*)right_state(newp);
	    state->solute = sr->solute;
	    if (cRparams->max_solute < state->solute)
		cRparams->max_solute = state->solute;
	    if (cRparams->min_solute > state->solute)
		cRparams->min_solute = state->solute;
	}
        return;
}       /* dirichlet_point_propagate */


static  void reaction_point_propagate(
        Front *front,
        POINTER wave,
        POINT *oldp,
        POINT *newp,
        HYPER_SURF_ELEMENT *oldhse,
        HYPER_SURF         *oldhs,
        double              dt,
        double              *V)
{
        double nor_speed,vel[MAXD];
        int i, dim = front->rect_grid->dim;
	double nor[MAXD];
        double p1[MAXD];
        double *p0 = Coords(oldp);
        double dn,*h = front->rect_grid->h;
        double s0,s1,grad_s,ans;
        STATE *sl,*sr,*state;
	CRT_PARAMS *cRparams = (CRT_PARAMS*)front->extra2;
        double *solute = cRparams->field->solute;
        double D = cRparams->D;
        double k = cRparams->k;
        double C_eq = cRparams->C_eq;
        double rho_s = cRparams->rho_s;
	double kappa;
	POINT_PROP_SCHEME point_prop_scheme = cRparams->point_prop_scheme;
	static void (*reaction_scheme)(CRT_PARAMS*,double,double,double,
				double,double,double*);
	static boolean first = YES;
	static double max_nor_speed = 0.0;
	REACTION_TYPE reaction_type = cRparams->reaction_type;

	if (first)
	{
	    first = NO;
	    switch (point_prop_scheme)
	    {
	    case EXPLICIT_EULER:
		reaction_scheme = euler_forward_scheme;
	    	break;
	    case IMPLICIT_EULER:
		reaction_scheme = euler_backward_scheme;
	    	break;
	    case MIDDLE_POINT:
		reaction_scheme = middle_point_scheme;
	    	break;
	    case CONSTANT_STATE:
		reaction_scheme = constant_state;
	    	break;
	    default:
		(void) printf("ERROR: Unknow reaction scheme!\n");
		clean_up(ERROR);
	    } 
	}
        FT_NormalAtPoint(oldp,front,nor,SOLUTE_COMP);
	/*Not yet working properly*/
	if (cRparams->add_curvature)
            FT_CurvatureAtPoint(oldp,front,&kappa);
	else
	    kappa = 0.0;

        dn = grid_size_in_direction(nor,h,dim);
        for (i = 0; i < dim; ++i)
            p1[i] = p0[i] + nor[i]*dn;

        FT_GetStatesAtPoint(oldp,oldhse,oldhs,(POINTER*)&sl,(POINTER*)&sr);
        state = (negative_component(oldhs) == SOLUTE_COMP) ? sl : sr;
        s0 = state->solute;

	FT_IntrpStateVarAtCoords(front,SOLUTE_COMP,p1,solute,
				getStateSolute,&s1,&state->solute);

        grad_s = (s1 - s0)/dn;
	switch (reaction_type)
	{
	case DEPOSITION_ONLY:
	    if (point_prop_scheme == CONSTANT_STATE)
	    	nor_speed = std::max(0.0,D*grad_s/rho_s);
	    else
	    	nor_speed = std::max(0.0,k*(s0 - C_eq)/rho_s);
	    break;
	case DISSOLUTION_ONLY:
	    if (point_prop_scheme == CONSTANT_STATE)
	    	nor_speed = std::min(0.0,D*grad_s/rho_s);
	    else
	    	nor_speed = std::min(0.0,k*(s0 - C_eq)/rho_s);
	    break;
	case REVERSIBLE_REACTION:
	    if (point_prop_scheme == CONSTANT_STATE)
	    	nor_speed = D*grad_s/rho_s;
	    else
	    	nor_speed = k*(s0 - C_eq)/rho_s;
	    break;
	default:
	    (void) printf("Unknow reaction type!\n");
	    clean_up(ERROR);
	}
        for (i = 0; i < dim; ++i)
        {
            vel[i] = nor[i]*nor_speed;
            Coords(newp)[i] = Coords(oldp)[i] + dt*vel[i];
        }
	if (max_nor_speed < nor_speed)
	{
	    max_nor_speed = nor_speed;
	    if (debugging("step_size"))
	    	(void) printf("Scaled max_prop_spacing = %f\n",
				max_nor_speed/dn*dt);
	}
	reaction_scheme(cRparams,dt,dn,kappa,s0,s1,&ans);

	/* Update the state of the new interface point */

	state = (negative_component(oldhs) == CRYSTAL_COMP) ? 
			(STATE*)left_state(newp) : (STATE*)right_state(newp);
        state->solute = rho_s;

	state = (negative_component(oldhs) == SOLUTE_COMP) ? 
			(STATE*)left_state(newp) : (STATE*)right_state(newp);
        s0 = state->solute = ans;

	if (cRparams->max_solute < state->solute)
		cRparams->max_solute = state->solute;
	if (cRparams->min_solute > state->solute)
		cRparams->min_solute = state->solute;
	nor_speed = 0.0;
        for (i = 0; i < dim; ++i)
        {
            vel[i] = nor[i]*k*(s0 - C_eq)/rho_s;
	    nor_speed += sqr(vel[i]);
            FT_RecordMaxFrontSpeed(i,fabs(vel[i]),NULL,Coords(newp),front);
        }
	FT_RecordMaxFrontSpeed(dim,sqrt(nor_speed),NULL,Coords(newp),front);
	if (the_point(oldp))
	{
	    printf("oldp = %f %f %f\n",Coords(oldp)[0],Coords(oldp)[1],
					Coords(oldp)[2]);
	    printf("newp = %f %f %f\n",Coords(newp)[0],Coords(newp)[1],
					Coords(newp)[2]);
	}
}       /* reaction_point_propagate */

static 	void euler_forward_scheme(
	CRT_PARAMS *cRparams,
	double dt,
	double dn,
	double kappa,
	double s0,
	double s1,
	double *ans)
{
        double D = cRparams->D;
        double k = cRparams->k;
        double C_eq = cRparams->C_eq;

	*ans = s0 + 2.0*dt/dn*(D*(s1 - s0)/dn - k*(s0 - C_eq));
}	/* end euler_forward_scheme */

static 	void euler_backward_scheme(
	CRT_PARAMS *cRparams,
	double dt,
	double dn,
	double kappa,
	double s0,
	double s1,
	double *ans)
{
        double D = cRparams->D;
        double k = cRparams->k;
        double C_eq = cRparams->C_eq;
	double c1,c2;

	c1 = 0.5*D*dt*(1.0/dn + kappa)/dn;
	c2 = 0.5*k*dt/dn;
	*ans = (s0 + c1*(s1 - s0) - c2*(s0 - C_eq) + c1*s1 + c2*C_eq)
			/(1.0 + c1 + c2);
}	/* end euler_backward_scheme */

static 	void middle_point_scheme(
	CRT_PARAMS *cRparams,
	double dt,
	double dn,
	double kappa,
	double s0,
	double s1,
	double *ans)
{
        double D = cRparams->D;
        double k = cRparams->k;
        double C_eq = cRparams->C_eq;
	double c1,c2;

	c1 = D*dt*(2.0/dn + kappa)/dn;
	c2 = 2.0*k*dt/dn;
	*ans = (s0 + c1*s1 + c2*C_eq)/(1.0 + c1 + c2);
}	/* end middle_point_scheme */

static 	void constant_state(
	CRT_PARAMS *cRparams,
	double dt,
	double dn,
	double kappa,
	double s0,
	double s1,
	double *ans)
{
        double C_eq = cRparams->C_eq;

	*ans = C_eq;
}	/* end middle_point_scheme */

extern	void	read_crystal_params(
	char *in_name,
	CRT_PARAMS *cRparams)
{
	FILE *infile;
	char string[200];

	infile = fopen(in_name,"r");
	CursorAfterString(infile,"Choose reaction type");
	(void) printf("\nAvailable reaction types are:\n");
	(void) printf("\tDEPOSITION_ONLY\n");
	(void) printf("\tDISSOLUTION_ONLY\n");
	(void) printf("\tREVERSIBLE_REACTION\n");
        CursorAfterString(infile,"Enter reaction type:");
        fscanf(infile,"%s",string);
	switch (string[0])
	{
	case 'd':
	case 'D':
	    if (string[1] == 'e' || string[1] == 'E')
		cRparams->reaction_type = DEPOSITION_ONLY;
	    else if (string[1] == 'i' || string[1] == 'I')
		cRparams->reaction_type = DISSOLUTION_ONLY;
	    else
	    {
		(void) printf("Unknown reaction type!\n");
		clean_up(ERROR);
	    }
	    break;
	case 'r':
	case 'R':
		cRparams->reaction_type = REVERSIBLE_REACTION;
	    break;
	default:
	    (void) printf("Unknown reaction type!\n");
	    clean_up(ERROR);
	}

	CursorAfterString(infile,"Diffusion coefficient:");
	fscanf(infile,"%lf",&cRparams->D);
	(void) printf("%f\n",cRparams->D);
	CursorAfterString(infile,"Growth rate:");
	fscanf(infile,"%lf",&cRparams->k);
	(void) printf("%f\n",cRparams->k);
	CursorAfterString(infile,"Equilibrium concentration:");
	fscanf(infile,"%lf",&cRparams->C_eq);
	(void) printf("%f\n",cRparams->C_eq);
	CursorAfterString(infile,"Ambient concentration:");
	fscanf(infile,"%lf",&cRparams->C_0);
	(void) printf("%f\n",cRparams->C_0);
	CursorAfterString(infile,"Crystal density:");
	fscanf(infile,"%lf",&cRparams->rho_s);
	(void) printf("%f\n",cRparams->rho_s);
	cRparams->gap = 0.1;
	CursorAfterStringOpt(infile,"Initial gap:");
	fscanf(infile,"%lf",&cRparams->gap);
	(void) printf("%f\n",cRparams->gap);
	cRparams->num_scheme = UNSPLIT_IMPLICIT;	/* default */
	cRparams->pde_order = 1;	/* default */
	if (CursorAfterStringOpt(infile,"Choose PDE scheme"))
	{
	    CursorAfterString(infile,"Enter scheme:");
	    fscanf(infile,"%s",string);
	    (void) printf("%s\n",string);
	    if ((string[0] == 'E' || string[0] == 'e') &&
	    	(string[1] == 'X' || string[1] == 'x')) 
	    	cRparams->num_scheme = UNSPLIT_EXPLICIT;
	    else if ((string[0] == 'I' || string[0] == 'i') &&
	    	(string[1] == 'M' || string[1] == 'm')) 
	    	cRparams->num_scheme = UNSPLIT_IMPLICIT;
	    else if ((string[0] == 'C' || string[0] == 'c') &&
	    	(string[1] == 'N' || string[1] == 'n')) 
	    	cRparams->num_scheme = CRANK_NICOLSON;
	}
	if (cRparams->num_scheme == UNSPLIT_IMPLICIT)
	{
	    if (CursorAfterStringOpt(infile,"Choose order of PDE scheme:"))
	    {
	    	fscanf(infile,"%d",&cRparams->pde_order);
		(void) printf("%d\n",cRparams->pde_order);
	    }
	}
	CursorAfterString(infile,"Choose point propagation scheme");
	(void) printf("\n");
        CursorAfterString(infile,"Enter scheme:");
        fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
        if ((string[0] == 'E' || string[0] == 'e') &&
            (string[1] == 'X' || string[1] == 'x'))
            cRparams->point_prop_scheme = EXPLICIT_EULER;
        else if ((string[0] == 'I' || string[0] == 'i') &&
            (string[1] == 'M' || string[1] == 'm'))
            cRparams->point_prop_scheme = IMPLICIT_EULER;
        else if ((string[0] == 'M' || string[0] == 'm') &&
            (string[1] == 'P' || string[1] == 'p'))
            cRparams->point_prop_scheme = MIDDLE_POINT;
        else if ((string[0] == 'C' || string[0] == 'c') &&
            (string[1] == 'S' || string[1] == 's'))
            cRparams->point_prop_scheme = CONSTANT_STATE;
        else
        {
            printf("Unknown point propagation scheme!\n");
            clean_up(ERROR);
        }
	CursorAfterString(infile,"Enter yes to add curvature effect:");
        fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	if (string[0] == 'Y' || string[0] == 'y')
            cRparams->add_curvature = YES;
	else
            cRparams->add_curvature = NO;
	fclose(infile);
}

extern boolean fractal_dimension(
	Front *front,
	SEED_PARAMS s_params,
	double *frac_dim,
	double *radius)
{
	double coords[MAXD],*center;
	double dist,r_sqr,r_max,r_min = s_params.seed_radius;
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
	boolean crystal_exist = NO;
	double ratio;
	boolean bdry_reached = NO;
	double margin[MAXD];

	if (s_params.grow_from_floor || s_params.grow_from_ceiling)
	    return bdry_reached;
	center = s_params.space_center[0];
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
	    crystal_exist = YES;
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
	crystal_exist = pp_max_status(crystal_exist);
#endif /* defined (__MPI__) */
	if (!crystal_exist)
	    return NO;

	/* Preparation for counting */

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
		    if (comp == CRYSTAL_COMP)
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
		    if (comp == CRYSTAL_COMP)
		    	++Nc;
		}
	    }
	}
#if defined (__MPI__)
	pp_global_isum(&N,1);
	pp_global_isum(&Nc,1);
#endif /* defined (__MPI__) */
	ratio = ((double)N)/((double)Nc);
	*frac_dim = (double)dim + log(ratio)/log(h[0]);
	return pp_max_status(bdry_reached);
}	/* end fractal_dimension */
	
extern void	record_1d_growth(
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
	    if (negative_component(hs) == SOLUTE_COMP)
	    {
		(*growth_data)[*count][1] = Coords(p)[0];
		(*growth_data)[*count][2] = sl->solute;
	    }
	    else if (positive_component(hs) == SOLUTE_COMP)
	    {
		(*growth_data)[*count][1] = Coords(p)[0];
		(*growth_data)[*count][2] = sr->solute;
	    }
	}
	*count += 1;
}	/* end record_1d_growth */

extern void plot_growth_data(
	char out_name[100],
	double **growth_data,
	int count)
{
	char fname[100];
	FILE *ofile;
	int i;

	sprintf(fname,"%s-posn.xg",out_name);
	ofile = fopen(fname,"w");
	fprintf(ofile,"\"Interface position vs. time\"\n");
	for (i = 0; i < count; ++i)
	    fprintf(ofile,"%f  %f\n",growth_data[i][0],growth_data[i][1]);
	fclose(ofile);

	sprintf(fname,"%s-solt.xg",out_name);
	ofile = fopen(fname,"w");
	fprintf(ofile,"\"Interface solute vs. time\"\n");
	for (i = 0; i < count; ++i)
	    fprintf(ofile,"%f  %f\n",growth_data[i][0],growth_data[i][2]);
	fclose(ofile);
}	/* end plot_growth_data */

extern	void read_seed_params(
	int dim,
	FILE *infile,
	SEED_PARAMS *s_params)
{
	char s[100];
	int i,j;

	s_params->dim = dim;

	s_params->num_floor_seeds = 0;
	s_params->num_ceiling_seeds = 0;
	s_params->num_space_seeds = 0;
	CursorAfterString(infile,"Are there seeds grow from floor:");
	fscanf(infile,"%s",s);
	(void) printf("%s\n",s);
	if (s[0] == 'y' || s[0] == 'Y')
	    s_params->grow_from_floor = YES;
	else
	    s_params->grow_from_floor = NO;
	if (s_params->grow_from_floor == YES)
	{
	    CursorAfterString(infile,"Enter floor level:");
	    fscanf(infile,"%lf",&s_params->floor_level);
	    (void) printf("%f\n",s_params->floor_level);
	    CursorAfterString(infile,"Enter number of floor seeds:");
	    fscanf(infile,"%d",&s_params->num_floor_seeds);
	    (void) printf("%d\n",s_params->num_floor_seeds);
	    FT_MatrixMemoryAlloc((POINTER*)&s_params->floor_center,
			s_params->num_floor_seeds,MAXD,sizeof(double));
	    for (i = 0; i < s_params->num_floor_seeds; ++i)
	    {
	    	sprintf(s,"Enter center coordinates of floor seed %d:",i+1);
	    	CursorAfterString(infile,s);
		for (j = 0; j < dim-1; ++j)
		{
		    fscanf(infile,"%lf ",&s_params->floor_center[i][j]);
	    	    (void) printf("%f ",s_params->floor_center[i][j]);
		}
		(void) printf("\n");
		s_params->floor_center[i][dim-1] = s_params->floor_level;
	    }
	}

	CursorAfterString(infile,"Are there seeds grow from ceiling:");
	fscanf(infile,"%s",s);
	(void) printf("%s\n",s);
	if (s[0] == 'y' || s[0] == 'Y')
	    s_params->grow_from_ceiling = YES;
	else
	    s_params->grow_from_ceiling = NO;
	if (s_params->grow_from_ceiling == YES)
	{
	    CursorAfterString(infile,"Enter ceiling level:");
	    fscanf(infile,"%lf",&s_params->ceiling_level);
	    (void) printf("%f\n",s_params->ceiling_level);
	    CursorAfterString(infile,"Enter number of ceiling seeds:");
	    fscanf(infile,"%d",&s_params->num_ceiling_seeds);
	    (void) printf("%d\n",s_params->num_ceiling_seeds);
	    FT_MatrixMemoryAlloc((POINTER*)&s_params->ceiling_center,
			s_params->num_ceiling_seeds,MAXD,sizeof(double));
	    for (i = 0; i < s_params->num_ceiling_seeds; ++i)
	    {
	    	sprintf(s,"Enter center coordinates of ceiling seed %d:",i+1);
	    	CursorAfterString(infile,s);
		for (j = 0; j < dim-1; ++j)
		{
		    fscanf(infile,"%lf ",&s_params->ceiling_center[i][j]);
	    	    (void) printf("%f ",s_params->ceiling_center[i][j]);
		}
		(void) printf("\n");
		s_params->ceiling_center[i][dim-1] = s_params->ceiling_level;
	    }
	}

	CursorAfterString(infile,"Are there seeds grow from space:");
	fscanf(infile,"%s",s);
	(void) printf("%s\n",s);
	if (s[0] == 'y' || s[0] == 'Y')
	    s_params->grow_from_space = YES;
	else
	    s_params->grow_from_space = NO;
	if (s_params->grow_from_space == YES)
	{
	    CursorAfterString(infile,"Enter number of space seeds:");
	    fscanf(infile,"%d",&s_params->num_space_seeds);
	    (void) printf("%d\n",s_params->num_space_seeds);
	    FT_MatrixMemoryAlloc((POINTER*)&s_params->space_center,
			s_params->num_space_seeds,MAXD,sizeof(double));
	    for (i = 0; i < s_params->num_space_seeds; ++i)
	    {
	    	sprintf(s,"Enter center coordinates of space seed %d:",i+1);
	    	CursorAfterString(infile,s);
		for (j = 0; j < dim; ++j)
		{
		    fscanf(infile,"%lf ",&s_params->space_center[i][j]);
	    	    (void) printf("%f ",s_params->space_center[i][j]);
		}
		(void) printf("\n");
	    }
	    s_params->add_space_seed_pert = NO;		// default
	    CursorAfterString(infile,"Enter yes to add perturbation:");
	    fscanf(infile,"%s",s);
	    (void) printf("%s\n",s);
	    if (s[0] == 'y' || s[0] == 'Y')
	    {
		s_params->add_space_seed_pert = YES;
	    	CursorAfterString(infile,"Enter number of period:");
		fscanf(infile,"%lf ",&s_params->nu);
	    	(void) printf("%f\n",s_params->nu);
	    	CursorAfterString(infile,"Enter amplitude:");
		fscanf(infile,"%lf ",&s_params->amp);
	    	(void) printf("%f\n",s_params->amp);
	    	CursorAfterString(infile,"Enter phase shift:");
		fscanf(infile,"%lf ",&s_params->phase);
	    	(void) printf("%f\n",s_params->phase);
	    }
	}

	CursorAfterString(infile,"Enter radius of seeds:");
	fscanf(infile,"%lf",&s_params->seed_radius);
	(void) printf("%f\n",s_params->seed_radius);
}	/* end read_seed_params */

extern void read_crt_dirichlet_bdry_data(
	char *inname,
	Front *front,
	F_BASIC_DATA f_basic)
{
	char msg[100],s[100];
	int i,dim = front->rect_grid->dim;
	FILE *infile = fopen(inname,"r");
	STATE state;
	HYPER_SURF *hs;

	for (i = 0; i < dim; ++i)
	{
	    if (f_basic.boundary[i][0] == DIRICHLET_BOUNDARY)
	    {
		hs = NULL;
		if (rect_boundary_type(front->interf,i,0) == DIRICHLET_BOUNDARY)
		    hs = FT_RectBoundaryHypSurf(front->interf,DIRICHLET_BOUNDARY,
						i,0);
		sprintf(msg,"For lower boundary in %d-th dimension",i);
		CursorAfterString(infile,msg);
		(void) printf("\n");
		CursorAfterString(infile,"Enter type of Dirichlet boundary:");
		fscanf(infile,"%s",s);
		(void) printf("%s\n",s);
		switch (s[0])
		{
		case 'c':			// Constant state
		case 'C':
		    CursorAfterString(infile,"Enter solute concentration:");
		    fscanf(infile,"%lf",&state.solute);
		    (void) printf("%f\n",state.solute);
		    FT_SetDirichletBoundary(front,NULL,NULL,NULL,
					(POINTER)&state,hs);
		    break;
		default: 
		    printf("ERROR: Dirichlet type %s not implemented\n",s);
		    clean_up(ERROR);
		}
	    }
	    if (f_basic.boundary[i][1] == DIRICHLET_BOUNDARY)
	    {
		hs = NULL;
		if (rect_boundary_type(front->interf,i,1) == DIRICHLET_BOUNDARY)
		    hs = FT_RectBoundaryHypSurf(front->interf,DIRICHLET_BOUNDARY,
						i,1);
		sprintf(msg,"For upper boundary in %d-th dimension",i);
		CursorAfterString(infile,msg);
		(void) printf("\n");
		CursorAfterString(infile,"Enter type of Dirichlet boundary:");
		fscanf(infile,"%s",s);
		(void) printf("%s\n",s);
		switch (s[0])
		{
		case 'c':			// Constant state
		case 'C':
		    CursorAfterString(infile,"Enter solute concentration:");
		    fscanf(infile,"%lf",&state.solute);
		    (void) printf("%f\n",state.solute);
		    FT_SetDirichletBoundary(front,NULL,NULL,NULL,
					(POINTER)&state,hs);
		    break;
		default: 
		    printf("ERROR: Dirichlet type %s not implemented\n",s);
		    clean_up(ERROR);
		}
	    }
	}
	fclose(infile);
}	/* end read_crt_dirichlet_bdry_data */


static void initFrontStates(
	Front *front)
{
	INTERFACE *intfc = front->interf;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
        STATE *sl,*sr;
	CRT_PARAMS *cRparams = (CRT_PARAMS*)front->extra2;
        double rho_s = cRparams->rho_s;

	next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
            FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);

            if (positive_component(hs) == SOLUTE_COMP)
                sr->solute = cRparams->C_eq;
            else if (positive_component(hs) == CRYSTAL_COMP)
                sr->solute = rho_s;
            else
                sr->solute = 0.0;
            if (negative_component(hs) == SOLUTE_COMP)
                sl->solute = cRparams->C_eq;
            else if (negative_component(hs) == CRYSTAL_COMP)
                sl->solute = rho_s;
            else
                sl->solute = 0.0;
        }
}	/* end initFrontStates */

extern void setInitialIntfc(
	Front *front,
	LEVEL_FUNC_PACK *level_func_pack,
	char *inname)
{
	int dim = front->rect_grid->dim;
	level_func_pack->neg_component = CRYSTAL_COMP;
        level_func_pack->pos_component = SOLUTE_COMP;
        level_func_pack->wave_type = GROWING_BODY_BOUNDARY;

	switch(dim)
	{
	case 1:
	    setInitialIntfc1d(front,level_func_pack,inname);
	    break;
	case 2:
	    setInitialIntfc2d(front,level_func_pack,inname);
	    break;
	case 3:
	    setInitialIntfc3d(front,level_func_pack,inname);
	    break;
	}

}	/* end setInitialIntfc */

static void setInitialIntfc1d(
	Front *front,
	LEVEL_FUNC_PACK *level_func_pack,
	char *inname)
{
	FILE *infile = fopen(inname,"r");
        char string[256];
	int dim = front->rect_grid->dim;
	CRT_PARAMS *cRparams = (CRT_PARAMS*)front->extra2;
	static double **point;

	FT_MatrixMemoryAlloc((POINTER*)&point,1,1,FLOAT);
        level_func_pack->num_points = 1;
        level_func_pack->point_array = point;

	CursorAfterString(infile,"Enter phase transition point coordinate:");
	fscanf(infile,"%lf",&point[0][0]);
	(void) printf("%f\n",point[0][0]);

	fclose(infile);
}	/* end setInitialIntfc1d */

static void setInitialIntfc2d(
	Front *front,
	LEVEL_FUNC_PACK *level_func_pack,
	char *inname)
{
	FILE *infile = fopen(inname,"r");
        char string[256];
	RECT_GRID *gr = front->rect_grid;
	int dim = gr->dim;
	static SEED_PARAMS s_params;
	static TRIANGLE_PARAMS t_params;
	static RECTANGLE_PARAMS r_params;
	CRT_PARAMS *cRparams = (CRT_PARAMS*)front->extra2;

	(void) printf("Available initial interface types are:\n");
	(void) printf("Seed (S)\n");
	(void) printf("Triangle (T)\n");
	(void) printf("Rectangle (R)\n");
	CursorAfterString(infile,"Enter initial interface type: ");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	switch (string[0])
	{
	case 's':
	case 'S':
	    read_seed_params(2,infile,&s_params);
	    level_func_pack->func_params = (POINTER)&s_params;
            level_func_pack->func = seed_func;
	    cRparams->func_params = (POINTER)&s_params;
	    cRparams->func = seed_func;
	    break;
	case 't':
	case 'T':
	    CursorAfterString(infile,"Enter three points for triangle:");
            CursorAfterString(infile,"Enter the first point x0, and y0: ");
            fscanf(infile,"%lf %lf\n",&t_params.x[0],&t_params.y[0]);
            (void) printf("%f %f\n",t_params.x[0],t_params.y[0]);
            CursorAfterString(infile,"Enter the second point x1, and y1: ");
            fscanf(infile,"%lf %lf\n",&t_params.x[1],&t_params.y[1]);
            (void) printf("%f %f\n",t_params.x[1],t_params.y[1]);
            CursorAfterString(infile,"Enter the third point x2, and y2: ");
            fscanf(infile,"%lf %lf\n",&t_params.x[2],&t_params.y[2]);
            (void) printf("%f %f\n",t_params.x[2],t_params.y[2]);
	    
	    level_func_pack->func_params = (POINTER)&t_params;
	    level_func_pack->func = triangle_func;
	    cRparams->func_params = NULL;
	    cRparams->func = NULL;
	    break;
	case 'r':
	case 'R':
	    CursorAfterString(infile,"Enter the left down point first\n");
	    CursorAfterString(infile,"Enter the first point x0, and y0: ");
	    fscanf(infile,"%lf %lf\n",&r_params.x0,&r_params.y0);
	    (void) printf("%f %f\n",r_params.x0,r_params.y0);
	    CursorAfterString(infile,"Enter the length of horizontal: ");
	    fscanf(infile,"%lf \n",&r_params.a);
	    (void) printf("%f \n",r_params.a);
	    CursorAfterString(infile,"Enter the length of vertical: ");
	    fscanf(infile,"%lf \n",&r_params.b);
	    (void) printf("%f \n",r_params.b);

	    level_func_pack->func_params = (POINTER)&r_params;
	    level_func_pack->func = rectangle_func;
	    cRparams->func_params = NULL;
	    cRparams->func = NULL;
	    break;
	default: 
	    (void) printf("Unknow type of initial interface!\n");
	    clean_up(ERROR);
	}

	fclose(infile);
}	/* end setInitialIntfc2d */

static void setInitialIntfc3d(
	Front *front,
	LEVEL_FUNC_PACK *level_func_pack,
	char *inname)
{
	FILE *infile = fopen(inname,"r");
        char string[256];
	int dim = front->rect_grid->dim;
	static SEED_PARAMS s_params;
	static RECT_BOX_PARAMS r_params;
	CRT_PARAMS *cRparams = (CRT_PARAMS*)front->extra2;

	CursorAfterString(infile,"Enter initial interface type: ");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	switch (string[0])
	{
	case 's':
	case 'S':
	    read_seed_params(3,infile,&s_params);
	    level_func_pack->func_params = (POINTER)&s_params;
            level_func_pack->func = seed_func;
	    level_func_pack->set_3d_bdry = YES;
	    cRparams->func_params = (POINTER)&s_params;
	    cRparams->func = seed_func;
	    break;
	case 'r':
	case 'R':
	    r_params.dim = 3;
	    CursorAfterString(infile,"Enter center of rectangular box:");
	    fscanf(infile,"%lf %lf %lf\n",&r_params.center[0],
				&r_params.center[1],&r_params.center[2]);
	    (void) printf("%f %f %f\n",r_params.center[0],r_params.center[1],
				r_params.center[2]);
	    CursorAfterString(infile,"Enter length of rectangular box:");
            fscanf(infile,"%lf %lf %lf\n",&r_params.length[0],
                                &r_params.length[1],&r_params.length[2]);
            (void) printf("%f %f %f\n",r_params.length[0],
                                r_params.length[1],r_params.length[2]);
	    level_func_pack->func_params = (POINTER)&r_params;
	    level_func_pack->func = rect_box_func;
	    level_func_pack->set_3d_bdry = YES;
	    cRparams->func_params = NULL;
	    cRparams->func = NULL;
	    break;
	default: 
	    (void) printf("Unknow type of initial interface!\n");
	    clean_up(ERROR);
	}

	fclose(infile);
}	/* end setInitialIntfc3d */

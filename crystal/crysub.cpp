
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

enum _FUEL_SAMPLE {
    	Void = 1,
	Metal,
	MetOx,
	Oxide
};
typedef enum _FUEL_SAMPLE FUEL_SAMPLE;

struct _SAMPLE_BDRY_PARAMS {
    	int bdry_type;
};
typedef struct _SAMPLE_BDRY_PARAMS SAMPLE_BDRY_PARAMS;

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
static double   level_sample_func(POINTER,double*);

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
        double rho_s;
	double kappa;
	POINT_PROP_SCHEME point_prop_scheme = cRparams->point_prop_scheme;
	static void (*reaction_scheme)(CRT_PARAMS*,double,double,double,
				double,double,double*);
	static boolean first = YES;
	static double max_nor_speed = 0.0;
	REACTION_TYPE reaction_type = cRparams->reaction_type;
	double CFL;

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
	if (cRparams->crystal_dens_func != NULL)
	    rho_s = (*cRparams->crystal_dens_func)(
			cRparams->crystal_dens_params,Coords(oldp));
	else
	    rho_s = cRparams->rho_s;
	
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

	CFL = Time_step_factor(front);
	if (fabs(nor_speed*dt) > CFL*dn)
	{
	    (void) printf("WARNING: speed higher than allowed by CFL\n");
	    (void) printf("CFL*dn/dt = %f  nor_speed = %f\n",
				CFL*dn/dt,nor_speed);
	    nor_speed = nor_speed*CFL*dn/dt/fabs(nor_speed);
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
        if (cRparams->crystal_dens_func != NULL)
	    rho_s = (*cRparams->crystal_dens_func)(
		    	cRparams->crystal_dens_params,Coords(newp));
	else
	    rho_s = cRparams->rho_s;

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

extern	void read_crystal_params(
	char *in_name,
	CRT_PARAMS *cRparams)
{
	FILE *infile;
	char string[200],s[100];
	static HALF_MOON_PARAMS half_moon_params;
	static FUEL_SAMPLE_PARAMS fuel_sample_params;

	infile = fopen(in_name,"r");
	CursorAfterString(infile,"Choose reaction type");
	(void) printf("\nAvailable reaction types are:\n");
	(void) printf("\tDEPOSITION_ONLY\n");
	(void) printf("\tDISSOLUTION_ONLY\n");
	(void) printf("\tREVERSIBLE_REACTION\n");
        CursorAfterString(infile,"Enter reaction type:");
        fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
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
	cRparams->crystal_dens_func = NULL;
	cRparams->crystal_dens_params = NULL;
	if (cRparams->reaction_type == DISSOLUTION_ONLY)
	{
	    if (CursorAfterStringOpt(infile,
		"Enter yes for varying crystal density:"))
	    {
	      fscanf(infile,"%s",s);
	      (void) printf("%s\n",s);
	      if (s[0] == 'y' || s[0] == 'Y')
	      {
		CursorAfterString(infile,
			"Enter crystal density function type:");
		fscanf(infile,"%s",string);
		(void) printf("%s\n",string);
		switch (string[0])
		{
		case 'h':
		case 'H':
		    CursorAfterString(infile,
			"Enter direction of the cut-line:");
		    fscanf(infile,"%d",&half_moon_params.dir);
		    (void) printf("%d\n",half_moon_params.dir);
		    CursorAfterString(infile,
			"Enter coordinate of the cut-line:");
		    fscanf(infile,"%lf",&half_moon_params.cutline);
		    (void) printf("%f\n",half_moon_params.cutline);
		    CursorAfterString(infile,
			"Enter lower and upper density:");
		    fscanf(infile,"%lf %lf",&half_moon_params.lower_dens,
				&half_moon_params.upper_dens);
		    (void) printf("%f %f\n",half_moon_params.lower_dens,
				half_moon_params.upper_dens);
		    cRparams->crystal_dens_func = half_moon_density;
		    cRparams->crystal_dens_params = (POINTER)&half_moon_params;
		    break;
		case 'p':
		case 'P':
		    cRparams->crystal_dens_func = perturbed_density;
		    cRparams->crystal_dens_params = NULL; 
		    break;
		case 'f':
		case 'F':
		    CursorAfterString(infile,
			"Enter crystal density of three samples:");
		    fscanf(infile,"%lf %lf %lf",&fuel_sample_params.rho_1,
			    			&fuel_sample_params.rho_2,
						&fuel_sample_params.rho_3);
		    (void) printf("%f %f %f\n", fuel_sample_params.rho_1,
			    			fuel_sample_params.rho_2,
						fuel_sample_params.rho_3);
		    cRparams->crystal_dens_func = fuel_sample_density;
		    cRparams->crystal_dens_params = 
					(POINTER)&fuel_sample_params;
		    break;
		default:
		    (void) printf("Unknow crystal density function!\n");
		    clean_up(ERROR);
		}
	      }
	      else
	      {
		CursorAfterString(infile,"Crystal density:");
		fscanf(infile,"%lf",&cRparams->rho_s);
		(void) printf("%f\n",cRparams->rho_s);
	      }
	    }
	}
	else
	{
	    CursorAfterString(infile,"Crystal density:");
	    fscanf(infile,"%lf",&cRparams->rho_s);
	    (void) printf("%f\n",cRparams->rho_s);
	}
	cRparams->gap = 0.0;     //No gap for dissolution
	CursorAfterStringOpt(infile,"Initial gap:");
	fscanf(infile,"%lf",&cRparams->gap);
	(void) printf("%f\n",cRparams->gap);
	cRparams->num_scheme = UNSPLIT_IMPLICIT;	/* default */
	cRparams->pde_order = 2;	/* default */
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
	int i_surf;

	for (i = 0; i < dim; ++i)
	{
	    if (f_basic.boundary[i][0] == DIRICHLET_BOUNDARY)
	    {
		hs = NULL;
		i_surf = 2*i;
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
		    FT_InsertDirichletBoundary(front,NULL,NULL,NULL,
					(POINTER)&state,hs,i_surf);
		    break;
		default: 
		    printf("ERROR: Dirichlet type %s not implemented\n",s);
		    clean_up(ERROR);
		}
	    }
	    if (f_basic.boundary[i][1] == DIRICHLET_BOUNDARY)
	    {
		hs = NULL;
		i_surf = 2*i + 1;
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
		    FT_InsertDirichletBoundary(front,NULL,NULL,NULL,
					(POINTER)&state,hs,i_surf);
		    break;
		default: 
		    printf("ERROR: Dirichlet type %s not implemented\n",s);
		    clean_up(ERROR);
		}
	    }
	}
	fclose(infile);
}	/* end read_crt_dirichlet_bdry_data */


extern void initFrontStates(
	Front *front)
{
	INTERFACE *intfc = front->interf;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
        STATE *sl,*sr;
	CRT_PARAMS *cRparams = (CRT_PARAMS*)front->extra2;
        double rho_s;

	next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
	    if (cRparams->crystal_dens_func != NULL)
	    	rho_s = (*cRparams->crystal_dens_func)(
			cRparams->crystal_dens_params,Coords(p));
	    else
		rho_s = cRparams->rho_s;
            
	    FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);

            if (positive_component(hs) == SOLUTE_COMP)
	    {
		if (wave_type(hs) == GROWING_BODY_BOUNDARY)
		    sr->solute = cRparams->C_0; 
	    	    //no gap for dissolution, C_eq for precipitation
		else
		    sr->solute = cRparams->C_0;
	    }
            else if (positive_component(hs) == CRYSTAL_COMP)
                sr->solute = rho_s;
            else
                sr->solute = 0.0;
            
	    if (negative_component(hs) == SOLUTE_COMP)
	    {
		if (wave_type(hs) == GROWING_BODY_BOUNDARY)
		    sl->solute = cRparams->C_0;	
	            //no gap for dissolution, C_eq for precipitation
		else
		    sl->solute = cRparams->C_0;
	    }
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
        char string[256],s[200];
	RECT_GRID *gr = front->rect_grid;
	int dim = gr->dim;
	static SEED_PARAMS s_params;
	static TRIANGLE_PARAMS t_params;
	static RECTANGLE_PARAMS r_params;
	static SLOTTED_CIRCLE_PARAMS sc_params;
	static SAMPLE_BDRY_PARAMS sb_params;
	CRT_PARAMS *cRparams = (CRT_PARAMS*)front->extra2;

	CursorAfterString(infile,"Enter initial interface type: ");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	switch (string[0])
	{
	case 's':
	case 'S':
	    if (string[1] == 'e' || string[1] == 'E')
	    {
	        read_seed_params(2,infile,&s_params);
	    	level_func_pack->func_params = (POINTER)&s_params;
            	level_func_pack->func = seed_func;
	    	cRparams->func_params = (POINTER)&s_params;
	    	cRparams->func = seed_func;
	    }
	    else if (string[1] == 'l' || string[1] == 'L')
	    {
		CursorAfterString(infile,"Enter center coordinates of the circle:");
		fscanf(infile,"%lf %lf\n",&sc_params.x0,&sc_params.y0);
		(void) printf("%f %f\n",sc_params.x0,sc_params.y0);
		CursorAfterString(infile,"Enter the radius of the circle:");
		fscanf(infile,"%lf\n",&sc_params.r);
		(void) printf("%f\n",sc_params.r);
		CursorAfterString(infile,"Enter the width and depth of the notch:");
		fscanf(infile,"%lf %lf\n",&sc_params.w,&sc_params.h);
		(void) printf("%f %f\n",sc_params.w,sc_params.h);
		sc_params.add_pert = NO;     		//default
		CursorAfterString(infile,"Enter yes to add perturbation:");
		fscanf(infile,"%s",s);
		(void) printf("%s\n",s);
		if (s[0] == 'y' || s[0] == 'Y')
		{
                    sc_params.add_pert = YES;
                    CursorAfterString(infile,"Enter number of period:");
                    fscanf(infile,"%lf ",&sc_params.nu);
                    (void) printf("%f\n",sc_params.nu);
                    CursorAfterString(infile,"Enter amplitude:");
                    fscanf(infile,"%lf ",&sc_params.amp);
                    (void) printf("%f\n",sc_params.amp);
                    CursorAfterString(infile,"Enter phase shift:");
                    fscanf(infile,"%lf ",&sc_params.phase);
                    (void) printf("%f\n",sc_params.phase);
		}
		
		level_func_pack->func_params = (POINTER)&sc_params;
		level_func_pack->func = slotted_circle_func;
		cRparams->func_params = NULL;
		cRparams->func = NULL;
	    }
	    else
	    {
		(void) printf("Unknow type of initial interface!\n");
		clean_up(ERROR);
	    }
	    break;
	case 'o':               // ORNL sample problem
	case 'O':
	    CursorAfterString(infile,"Enter sample boundary type:");
	    fscanf(infile,"%s",s);
	    (void) printf("%s\n",s);
	    switch (s[0])
	    {
	    case 'r':
	    case 'R':
		sb_params.bdry_type = REFLECTION_BOUNDARY;
		break;
	    case 'd':
	    case 'D':
		sb_params.bdry_type = DIRICHLET_BOUNDARY;
		break;
	    default:
		(void) printf("Unknown sample boundary %s\n",string);
		clean_up(ERROR);
	    }
	    level_func_pack->func_params = (POINTER)&sb_params;
	    level_func_pack->func = level_sample_func;
	    cRparams->func_params = NULL;
	    cRparams->func = NULL;
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
	case 'f':
	case 'F':
	    CursorAfterString(infile,"Enter center coordinates of the circle:");
	    fscanf(infile,"%lf %lf\n",&sc_params.x0,&sc_params.y0);
	    (void) printf("%f %f\n",sc_params.x0,sc_params.y0);
	    CursorAfterString(infile,"Enter the radius of the circle:");
	    fscanf(infile,"%lf\n",&sc_params.r);
	    (void) printf("%f\n",sc_params.r);
	    CursorAfterString(infile,"Enter the width and depth of four slots:");
	    fscanf(infile,"%lf %lf\n",&sc_params.w,&sc_params.h);
	    (void) printf("%f %f\n",sc_params.w,sc_params.h);
	    sc_params.add_pert = NO;

	    level_func_pack->func_params = (POINTER)&sc_params;
	    level_func_pack->func = four_slotted_circle_func;
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

extern double half_moon_density(
	POINTER params,
	double *coords)
{
	HALF_MOON_PARAMS *hm_params = (HALF_MOON_PARAMS*)params;
	int dir = hm_params->dir;
	double cutline = hm_params->cutline;

	return (coords[dir] < cutline) ? hm_params->lower_dens :
					 hm_params->upper_dens;
}	/* end half_moon_density */

extern double perturbed_density(
	POINTER params,
	double *coords)
{
	double x = coords[0] - 2.0;
	double y = coords[1] - 3.0;
	
	double average,perturbed;

	average = 5*sin(310*x*2.5+(1.0/6.0)*PI)*sin(207*y*2.5+(1.0/5.0)*PI)+8;
	perturbed = 10*sin(100000*PI*x*2.5+(1.0/4.0)*PI)*sin(1000011*y*2.5+(3.0/10.0)*PI);
	if (perturbed < 0.0 && average + perturbed > 0.0)
	{
	    if (average + perturbed > 0.6 && average + perturbed < 9.5)
	    {
		return 13.0;
	    }
	    else if (average + perturbed >= 9.5 && average + perturbed < 11.8)
	    {
		return 0.5;
	    }
	    else
     	    {
		return average + perturbed;
	    }
	}
	else
	{
	    if (average > 0.6 && average < 9.8)
	    {
		return 13.0;
	    }
	    else if (average >= 9.8 && average < 11.8)
	    {	
		return 0.5;
	    }
	    else
	    {
		return average;
	    }
	}
}

struct _QUAD {
        double *verts[4];
	FUEL_SAMPLE quad_sample;
};
typedef struct _QUAD QUAD;

static boolean is_in_quad(double*,QUAD*);
static void vert_icoords(double*,double*,double*,int*);
extern boolean sample_func(POINTER,double*,QUAD*);

extern double fuel_sample_density(
	POINTER params,
	double *coords)
{
        QUAD qad;
	FUEL_SAMPLE_PARAMS *fparams = (FUEL_SAMPLE_PARAMS*)params;
	if (!sample_func(NULL,coords,&qad))
	    return fparams->rho_2;
	switch (qad.quad_sample)
	{
	case Metal:
	    return fparams->rho_2;				
	case Void:
	case Oxide:		
	    return fparams->rho_2;	        
	case MetOx:		
	    return fparams->rho_3;				
	}
}       /* end fuel_sample_density */

static double level_sample_func(
	POINTER func_params,
	double *coords)
{
        QUAD qad;
	if (!sample_func(func_params,coords,&qad)) return 1.0;
	if (qad.quad_sample == Metal) return 1.0;
	else return -1.0;
}       /* end level_sample_func */

extern boolean sample_func(
	POINTER params,
	double *coords,
	QUAD *qad)
{
        int i,nq;
	static boolean first = YES;
	static double **verts;
	static int **quads,**quad_store;
	static double L[MAXD],U[MAXD],h[MAXD];
	static int **num_quads,****mesh_quads;
	static int gmax[MAXD];
	static double epsilon;
	int **pq,icoords[MAXD];
	double crds[MAXD];
	SAMPLE_BDRY_PARAMS *sb_params = (SAMPLE_BDRY_PARAMS*)params;

	if (first)
        {
            FILE *infile;
            char string[200];
            double x[MAXD];
            boolean recorded;
            int j,k,nv,count;
            int quad[4];
            int **icoords_quad;						
	    
	    if (debugging("trace"))
	        (void) printf("Entering sample_func() first time\n");
            first = NO;
            for (i = 0; i < 2; ++i)
            {
                L[i] =  HUGE;
                U[i] = -HUGE;
            }
            infile = fopen("sample/fuel-test1.mesh","r");
            fgetstring(infile,"vertices = [");
            nv = 0;
            while (YES)
            {
                nv++;
                fgetstring(infile,"[");
                fscanf(infile,"%lf",&x[0]);
                string[0] = getc(infile);
                fscanf(infile,"%lf",&x[1]);
                fscanf(infile,"%s",string);
                if (L[0] > x[0]) L[0] = x[0];
                if (L[1] > x[1]) L[1] = x[1];
                if (U[0] < x[0]) U[0] = x[0];
                if (U[1] < x[1]) U[1] = x[1];
                if (string[0] == '}')
                    break;
                else
                    fscanf(infile,"%*s %*s %*s");
            }

	    fgetstring(infile,"elements = [");
	    nq = 0;
	    while (YES)
	    {
	    	nq++;
		fgetstring(infile,"[");
		for (i = 0; i < 4; ++i)
		{
		    fscanf(infile,"%d",&quad[i]);
		    fgetstring(infile,",");
		}
		fgetstring(infile,"\"");
		fscanf(infile,"%s",string);
		if (strncmp(string,"MetOx",5) == 0)
		{
		    if (string[7] != ',') break;
		}
		else if (strncmp(string,"Oxide",5) == 0)
		{
		    if (string[7] != ',') break;
		}
		else if (strncmp(string,"Metal",5) == 0)
		{
		    if (string[7] != ',') break;
		}
		else if (strncmp(string,"Void",4) == 0)
		{
		    if (string[6] != ',') break;
		}
		fscanf(infile,"%*s %*s %*s");
	    }
	    if (debugging("fuel_sample"))
	    {
	        (void) printf("L = %f %f  U = %f %f\n",L[0],L[1],U[0],U[1]);
	        (void) printf("Number of vertices = %d\n",nv);
	        (void) printf("Number of quads = %d\n",nq);
	        clean_up(0);
	    }
	    FT_MatrixMemoryAlloc((POINTER*)&verts,nv,MAXD,sizeof(double));
	    FT_MatrixMemoryAlloc((POINTER*)&quads,nq,5,sizeof(int));

	    rewind(infile);
	    nv = 0;
	    fgetstring(infile,"vertices = [");
	    while (YES)
	    {
	        fgetstring(infile,"[");
		fscanf(infile,"%lf",&verts[nv][0]);
		string[0] = getc(infile);
		fscanf(infile,"%lf",&verts[nv][1]);
		fscanf(infile,"%s",string);
		if (string[0] == '}')
		    break;
		else
		    fscanf(infile,"%*s %*s %*s");
		nv++;
	    }
	    nv++;
	
	    for (i = 0; i < 2; ++i)
	    {
	    	gmax[i] = 200;
		h[i] = (U[i] - L[i])/200.0;
	    }
	    epsilon = 0.00000001*h[0];
	    FT_MatrixMemoryAlloc((POINTER*)&num_quads,gmax[0],gmax[1],
		    			sizeof(int));
	    FT_MatrixMemoryAlloc((POINTER*)&icoords_quad,4,MAXD,sizeof(int));

	    nq = count = 0;
	    fgetstring(infile,"elements = [");
	    while (YES)
	    {
		fgetstring(infile,"[");
		for (i = 0; i < 4; ++i)
		{
		    fscanf(infile,"%d",&quads[nq][i]);
		    fgetstring(infile,",");
		    vert_icoords(L,h,verts[quads[nq][i]],icoords_quad[i]);
		    for (j = 0; j < 2; ++j)
		    {
		        if (icoords_quad[i][j] < 0)
			    icoords_quad[i][j] = 0;
			if (icoords_quad[i][j] >= gmax[j])
			    icoords_quad[i][j] = gmax[j] - 1;
		    }
		}
		for (i = 0; i < 4; ++i)
		{
		    recorded = NO;
		    for (j = 0; j < i; ++j)
		    {
		        if (icoords_quad[i][0] == icoords_quad[j][0] &&
			    icoords_quad[i][1] == icoords_quad[j][1])
			    recorded = YES;
		    }
		    if (!recorded)
		    {
		        num_quads[icoords_quad[i][0]][icoords_quad[i][1]]++;
			count++;
		    }
		}
		fgetstring(infile,"\"");
		fscanf(infile,"%s",string);
		if (strncmp(string,"MetOx",5) == 0)
		{
		    quads[nq][4] = (int)MetOx;
		    if (string[7] != ',') break;
		}
		else if (strncmp(string,"Oxide",5) == 0)
		{
		    quads[nq][4] = (int)Oxide;
		    if (string[7] != ',') break;
		}
		else if (strncmp(string,"Metal",5) == 0)
		{
		    quads[nq][4] = (int)Metal;
		    if (string[7] != ',') break;
		}
		else if (strncmp(string,"Void",4) == 0)
		{
		    quads[nq][4] = (int)Void;
		    if (string[6] != ',') break;
		}
		fscanf(infile,"%*s %*s %*s");
		nq++;
	    }
	    nq++;
	    FT_VectorMemoryAlloc((POINTER*)&quad_store,count,sizeof(int*));
	    FT_MatrixMemoryAlloc((POINTER*)&mesh_quads,gmax[0],gmax[1],
		    			sizeof(int**));
	    pq = quad_store;
	    count = 0;
	    for (i = 0; i < gmax[0]; ++i)
	    for (j = 0; j < gmax[1]; ++j)
	    {
	        mesh_quads[i][j] = pq;
		pq += num_quads[i][j];
		num_quads[i][j] = 0;
	    }
	    for (i = 0; i < nq; ++i)
	    {
	        for (j = 0; j < 4; ++j)
		{
		    vert_icoords(L,h,verts[quads[i][j]],icoords_quad[j]);
		    for (k = 0; k < 2; ++k)
		    {
		        if (icoords_quad[j][k] < 0)
			    icoords_quad[j][k] = 0;
			if (icoords_quad[j][k] >= gmax[k])
			    icoords_quad[j][k] = gmax[k] - 1;
		    }
		    recorded = NO;
		    for (k = 0; k < j; ++k)
		    {
		        if (icoords_quad[j][0] == icoords_quad[k][0] &&
			    icoords_quad[j][1] == icoords_quad[k][1])
			    recorded = YES;
		    }
		    if (!recorded)
		    {
		        int ii,jj;
			ii = icoords_quad[j][0];
			jj = icoords_quad[j][1];
			mesh_quads[ii][jj][num_quads[ii][jj]++] = quads[i];
			count++;
		    }
		}
	    }
	    if (debugging("gview_sample"))
	    {
	        FILE *gfile = fopen("sample.list","w");
		fprintf(gfile,"{ LIST \n");
		fprintf(gfile,"{ OFF\n");
		fprintf(gfile,"%d %d %d\n",nv,nq,0);
		for (i = 0; i < nv; ++i)
		    fprintf(gfile,"%f %f %f\n",verts[i][0],verts[i][1],0.0);
		for (i = 0; i < nq; ++i)
		{
		    fprintf(gfile,"4 %d %d %d %d ",quads[i][0],quads[i][1],
			    		quads[i][2],quads[i][3]);
		    switch (quads[i][4])
		    {
		    case Void:
			fprintf(gfile,"%f %f %f %f\n",1.0,1.0,1.0,0.5);
			break;
		    case Metal:
			fprintf(gfile,"%f %f %f %f\n",1.0,0.0,0.0,0.5);
			break;
		    case MetOx:
			fprintf(gfile,"%f %f %f %f\n",0.0,1.0,0.0,0.5);
			break;
		    case Oxide:
			fprintf(gfile,"%f %f %f %f\n",0.0,0.0,1.0,0.5);
			break;
		    }
		}
		fprintf(gfile,"}\n}\n");
		fclose(gfile);
	    }
	}

	/* For possible reflection */
	for (i = 0; i < 2; ++i)
	{
	    crds[i] = coords[i];
	    if (coords[i] < L[i])
	    {
	        if (sb_params == NULL)
		    return NO;
		else if (sb_params->bdry_type == DIRICHLET_BOUNDARY)
		    return NO;
		else if (sb_params->bdry_type == REFLECTION_BOUNDARY)
		    crds[i] = L[i] + (L[i] - coords[i]);
	    }
	    if (coords[i] > U[i])
	    {
	        if (sb_params == NULL)
		    return NO;
		else if (sb_params->bdry_type == DIRICHLET_BOUNDARY)
		    return NO;
		else if (sb_params->bdry_type == REFLECTION_BOUNDARY)
		    crds[i] = U[i] - (coords[i] - U[i]);
	    }
	    if (coords[i] == L[i]) crds[i] = L[i] + epsilon;
	    if (coords[i] == U[i]) crds[i] = U[i] - epsilon;
	}
	vert_icoords(L,h,crds,icoords);
	for (i = 0; i < 2; ++i)
	{
	    if (icoords[i] < 0)
		icoords[i] = 0;
	    if (icoords[i] >= gmax[i])
		icoords[i] = gmax[i] - 1;
	}
	nq = num_quads[icoords[0]][icoords[1]];
	pq = mesh_quads[icoords[0]][icoords[1]];
	for (i = 0; i < nq; ++i)
	{
	    qad->verts[0] = verts[pq[i][0]];
	    qad->verts[1] = verts[pq[i][1]];
	    qad->verts[2] = verts[pq[i][2]];
	    qad->verts[3] = verts[pq[i][3]];
	    if (is_in_quad(crds,qad))
	    {
	        qad->quad_sample = (FUEL_SAMPLE)pq[i][4];
		if (debugging("sample_func"))
		{
		    (void) printf("coords = %f %f\n",coords[0],coords[1]);
		    (void) printf("icoords = %d %d\n",icoords[0],icoords[1]);
		    (void) printf("Number of quad in cell: %d\n",nq);
		    (void) printf("In quads[%d]:\n",i);
		    (void) printf("Vertices:\n");
		    (void) printf("qad.verts[0] = %f %f\n",
			    		qad->verts[0][0],qad->verts[0][1]);
		    (void) printf("qad.verts[1] = %f %f\n",
			    		qad->verts[1][0],qad->verts[1][1]);
		    (void) printf("qad.verts[2] = %f %f\n",
			    		qad->verts[2][0],qad->verts[2][1]);
		    (void) printf("qad.verts[3] = %f %f\n",
			    		qad->verts[3][0],qad->verts[3][1]);
		    (void) printf("Quad sample: %d\n",qad->quad_sample);
		}
		return YES;
	    }
	}
	if (debugging("sample_func"))
	{
	    (void) printf("No quad found!\n");
	}
	return NO;
}       /* end sample_func */

static boolean is_in_quad(
	double *coords,
	QUAD *quad)
{
	double **verts = quad->verts;
	double v1[MAXD],v2[MAXD],c0,cr;
	int i,j;

	for (j = 0; j < 2; ++j)
	{
	    v1[j] = coords[j] - verts[0][j];
	    v2[j] = verts[1][j] - verts[0][j];
	}
	Cross2d(v1,v2,c0);
	for (i = 1; i < 4; ++i)
	{
	    for (j = 0; j < 2; ++j)
	    {
	        v1[j] = coords[j] - verts[i][j];
		v2[j] = verts[(i+1)%4][j] - verts[i][j];
	    }
	    Cross2d(v1,v2,cr);
	    if (!same_sign(c0,cr)) return NO;
	}
	return YES;
}       /* end is_in_quad */

static void vert_icoords(
	double *L,
	double *h,
	double *vert,
	int *icoords)
{
	int i;
	for (i = 0; i < 2; ++i)
	{
	    icoords[i] = irint(floor((vert[i] - L[i])/h[i]));
	}
}       /* end vert_icoords */

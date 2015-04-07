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

#include "inverse.h"
#define eps_plus  80
#define eps_minus  2

static int tangential_direction(int i,double *T,double *N,int dim);

extern double sourceFunc(
	POINTER params,
	int D, 
	double *P)
{
	CIM_PARAMS *cim_params = (CIM_PARAMS*)params;

	if (cim_params->jump_type == EXACT_JUMP)
	    return exact_source(params,D,P);
	else
	{
	    double r = sqrt(sqr(P[0])+sqr(P[1]));
	    double v = exact_eps(D)*(9.0*r-1.0-16.0*r*r);
	    return v;
	}
}

extern double jumpEpsGradDotNorm(
	POINTER params,
	int D, 
	double *N, 
	double *P)
{
	CIM_PARAMS *cim_params = (CIM_PARAMS*)params;

	switch (cim_params->jump_type)
	{
	case EXACT_JUMP:
    	    return exact_jump_eps_gradu_dot_n(params,D,N,P);
	case CONSTANT_JUMP:
	    return cim_params->jump_eps_grad_u_dot_n*D;
	default:
	    (void) printf("Case not implemented!\n");
	    clean_up(ERROR);
	}
}	/* end jumpEpsGradDotNorm */

extern double jumpU(
	POINTER params,
	int D, 
	double *P)
{
	CIM_PARAMS *cim_params = (CIM_PARAMS*)params;

	switch (cim_params->jump_type)
	{
	case EXACT_JUMP:
    	    return exact_jump_u(params,D,P);
	case CONSTANT_JUMP:
	    return cim_params->jump_u*D;
	default:
	    (void) printf("Case not implemented!\n");
	    clean_up(ERROR);
	}
}	/* end jumpU */

extern double jumpGradDotTan(
	POINTER params,
	int D, 
	int i, 
	double *N, 
	double *P)
{
	CIM_PARAMS *cim_params = (CIM_PARAMS*)params;
	switch (cim_params->jump_type)
	{
	case EXACT_JUMP:
	    return exact_jump_gradu_dot_t(params,D,i,N,P);
	case CONSTANT_JUMP:
	    return cim_params->jump_grad_u_dot_t;
	default:
	    (void) printf("Case not implemented!\n");
	    clean_up(ERROR);
	}
}	/* end jumpGradDotTan */

extern double getStateU(POINTER state)
{
	STATE *st = (STATE*)state;
	return st->u;
}

extern void assignStateU(
	double u,
	POINTER state)
{
	STATE *st = (STATE*)state;
	st->u = u;
}

extern int cim_find_state_at_crossing(
	Front *front,
	int *icoords,
        GRID_DIRECTION dir,
        int comp,
        POINTER *state,
        HYPER_SURF **hs,
        double *crx_coords)
{
	boolean status;
	STATE *st;
	int D = (comp == 2) ? 1 : -1;
	CIM_PARAMS *cim_params = (CIM_PARAMS*)front->extra1;
	INTERFACE *grid_intfc;

	status = FT_StateStructAtGridCrossing(front,grid_intfc,icoords,dir,
				comp,state,hs,crx_coords);
        if (status == NO) return NO_PDE_BOUNDARY;
	st = (STATE*)*state;

        if (wave_type(*hs) == FIRST_PHYSICS_WAVE_TYPE) 
	    return NO_PDE_BOUNDARY;
	else if (wave_type(*hs) == DIRICHLET_BOUNDARY)
	{
	    if (cim_params->jump_type == EXACT_JUMP)
		st->u = exact_solution((POINTER)cim_params,D,crx_coords);
            return DIRICHLET_PDE_BOUNDARY;
	}
	
}       /* cim_find_state_at_crossing */

extern void initCimIntfcParams(
	char *inname,
	Front *front,
	LEVEL_FUNC_PACK *level_func_pack)
{
	CIM_PARAMS *cim_params = (CIM_PARAMS*)front->extra1;	
	static CIRCLE_PARAMS circle_params;
	int i,dim = front->rect_grid->dim;
	FILE *infile = fopen(inname,"r");

        level_func_pack->wave_type = FIRST_PHYSICS_WAVE_TYPE;
        level_func_pack->neg_component = 1;
        level_func_pack->pos_component = 2;
	switch (cim_params->intfc_num)
	{
	case 1:
	    circle_params.dim = dim;
	    CursorAfterString(infile,"Enter center of the circle:");
	    for (i = 0; i < dim; ++i)
	    {
		fscanf(infile,"%lf ",&circle_params.cen[i]);
		(void) printf(" %f",circle_params.cen[i]);
            	circle_params.cen[i] = 0.0;
	    }
	    (void) printf("\n");
	    CursorAfterString(infile,"Enter radius of the circle:");
	    fscanf(infile,"%lf",&circle_params.R);
	    (void) printf("%f\n",circle_params.R);
            circle_params.add_plan_surf = NO;
            circle_params.add_perturbation = NO;
            level_func_pack->func = level_circle_func;
            level_func_pack->func_params = (POINTER)&circle_params;
	    break;
	case 2:
	    circle_params.dim = dim;
            level_func_pack->func = intfc_func_case2;
            level_func_pack->func_params = NULL;
	    break;
	case 3:
	    circle_params.dim = dim;
            level_func_pack->func = intfc_func_case3;
            level_func_pack->func_params = NULL;
	    break;
	case 4:
	    circle_params.dim = dim;
            level_func_pack->func = intfc_func_case4;
            level_func_pack->func_params = NULL;
	    break;
	case 5:
	    circle_params.dim = dim;
            level_func_pack->func = intfc_func_case5;
            level_func_pack->func_params = NULL;
	    break;
	case 6:
	    circle_params.dim = dim;
            level_func_pack->func = intfc_func_case6;
            level_func_pack->func_params = NULL;
	    break;
	case 7:
	    circle_params.dim = dim;
            level_func_pack->func = intfc_func_case7;
            level_func_pack->func_params = NULL;
	    break;
	case 8:
	    circle_params.dim = dim;
            level_func_pack->func = intfc_func_case8;
            level_func_pack->func_params = NULL;
	    break;
	case 9:
	    circle_params.dim = dim;
            level_func_pack->func = intfc_func_case9;
            level_func_pack->func_params = NULL;
	    break;
	case 10:
	    circle_params.dim = dim;
            level_func_pack->func = intfc_func_case10;
            level_func_pack->func_params = NULL;
	    break;
	default:
	    (void) printf("Unknown interface case number!\n");
	    clean_up(ERROR);
	}
}	/* end initCimIntfcParams */

extern void print_front_states(
        FILE *outfile,
        Front *front)
{
        INTERFACE *intfc = front->interf;
        STATE *sl,*sr;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;

        /* Initialize states at the interface */
        fprintf(outfile,"Interface u states:\n");
        int count = 0;
        next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
            count++;
            FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
            fprintf(outfile,"%24.18g %24.18g\n",sl->u,sr->u);
        }

}       /* end print_front_states */

void read_front_states(
	FILE *infile,
	Front *front)
{
        INTERFACE *intfc = front->interf;
        STATE *sl,*sr;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
        double x;

        /* Initialize states at the interface */
        next_output_line_containing_string(infile,"Interface u states:");
        next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
            FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
            fscanf(infile,"%lf",&x);
            sl->u = x;
            fscanf(infile,"%lf",&x);
            sr->u = x;
        }
}	/* end solute_read_front_states */

/***************************************************************
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
****************************************************************/

#include "iFluid.h"
#include "ifluid_basic.h"

	/*  Function Declarations */
static void neumann_point_propagate(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
static void dirichlet_point_propagate(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
static void contact_point_propagate(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
static void rgbody_point_propagate(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
static void iF_flowThroughBoundaryState2d(double*,HYPER_SURF*,Front*,POINTER,
        		POINTER);
static void iF_flowThroughBoundaryState3d(double*,HYPER_SURF*,Front*,POINTER,
        		POINTER);
static void get_time_dependent_params(int,FILE*,POINTER*);
static void addToEnergyFlux(RECT_GRID*,HYPER_SURF*,double*,double*,int,int,
			boolean);

static double (*getStateVel[3])(POINTER) = {getStateXvel,getStateYvel,
                                        getStateZvel};

static void ifluid_compute_force_and_torque2d(Front*,HYPER_SURF*,double,
                        double*,double*);
static void ifluid_compute_force_and_torque3d(Front*,HYPER_SURF*,double,
                        double*,double*);
static boolean force_on_hse(HYPER_SURF_ELEMENT*,HYPER_SURF*,RECT_GRID*,double*,
                                        double*,double*,boolean);
static boolean force_on_hse2d(HYPER_SURF_ELEMENT*,HYPER_SURF*,RECT_GRID*,
                                        double*,double*,double*,boolean);
static boolean force_on_hse3d(HYPER_SURF_ELEMENT*,HYPER_SURF*,RECT_GRID*,
                                        double*,double*,double*,boolean);
static double intrp_between(double,double,double,double,double);

extern double getStatePres(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->pres;
}	/* end getStatePres */

extern double getStatePhi(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->phi;
}	/* end getStatePres */

extern double getStateVort(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->vort;
}	/* end getStateVort */

extern double getStateXvel(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->vel[0];
}	/* end getStateXvel */

extern double getStateYvel(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->vel[1];
}	/* end getStateYvel */

extern double getStateZvel(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->vel[2];
}	/* end getStateZvel */

extern double getStateXimp(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->impulse[0];
}	/* end getStateXimp */

extern double getStateYimp(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->impulse[1];
}	/* end getStateYimp */

extern double getStateZimp(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->impulse[2];
}	/* end getStateZimp */

extern void read_iF_dirichlet_bdry_data(
	char *inname,
	Front *front,
	F_BASIC_DATA f_basic)
{
	char msg[100],s[100];
	int i,k,dim = front->rect_grid->dim;
	FILE *infile = fopen(inname,"r");
	static STATE *state;
	POINTER func_params;
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
		    FT_ScalarMemoryAlloc((POINTER*)&state,sizeof(STATE));
		    CursorAfterString(infile,"Enter velocity:");
		    for (k = 0; k < dim; ++k)
		    {
			fscanf(infile,"%lf",&state->vel[k]);
			(void) printf("%f ",state->vel[k]);
		    }
		    (void) printf("\n");
		    CursorAfterString(infile,"Enter pressure:");
		    fscanf(infile,"%lf",&state->pres);
		    (void) printf("%f\n",state->pres);
		    state->phi = getPhiFromPres(front,state->pres);
		    FT_InsertDirichletBoundary(front,NULL,NULL,
					NULL,(POINTER)state,hs,i_surf);
		    break;
		case 'f':			// Flow through state
		case 'F':
		    FT_InsertDirichletBoundary(front,iF_flowThroughBoundaryState,
					"flowThroughBoundaryState",NULL,
					NULL,hs,i_surf);
		    break;
		case 't':			// Flow through state
		case 'T':
		    get_time_dependent_params(dim,infile,&func_params);
		    FT_InsertDirichletBoundary(front,iF_timeDependBoundaryState,
					"iF_timeDependBoundaryState",
					func_params,NULL,hs,i_surf);
		    break;
		default:
		    (void) printf("Unknown Dirichlet boundary!\n");
		    clean_up(ERROR);
		}
	    }
            if (f_basic.boundary[i][1] == DIRICHLET_BOUNDARY)
	    {
		hs = NULL;
		i_surf = 2*i + 1;
                if (rect_boundary_type(front->interf,i,1) == DIRICHLET_BOUNDARY)                    hs = FT_RectBoundaryHypSurf(front->interf,DIRICHLET_BOUNDARY,
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
		    FT_ScalarMemoryAlloc((POINTER*)&state,sizeof(STATE));
		    CursorAfterString(infile,"Enter velocity:");
		    for (k = 0; k < dim; ++k)
		    {
			fscanf(infile,"%lf ",&state->vel[k]);
			(void) printf("%f ",state->vel[k]);
		    }
		    (void) printf("\n");
		    CursorAfterString(infile,"Enter pressure:");
		    fscanf(infile,"%lf",&state->pres);
		    (void) printf("%f\n",state->pres);
		    state->phi = getPhiFromPres(front,state->pres);
		    FT_InsertDirichletBoundary(front,NULL,NULL,
					NULL,(POINTER)state,hs,i_surf);
		    break;
		case 'f':			// Flow through state
		case 'F':
		    FT_InsertDirichletBoundary(front,iF_flowThroughBoundaryState,
					"flowThroughBoundaryState",NULL,
					NULL,hs,i_surf);
		    break;
		case 't':			// Flow through state
		case 'T':
		    get_time_dependent_params(dim,infile,&func_params);
		    FT_InsertDirichletBoundary(front,iF_timeDependBoundaryState,
					"iF_timeDependBoundaryState",
					func_params,NULL,hs,i_surf);
		    break;
		default:
		    (void) printf("Unknown Dirichlet boundary!\n");
		    clean_up(ERROR);
		}
	    }
	}
	fclose(infile);
}	/* end read_iF_dirichlet_bdry_data */

extern void restart_set_dirichlet_bdry_function(Front *front)
{
	INTERFACE *intfc = front->interf;
	int i;
	BOUNDARY_STATE  *bstate;
	const char *s;
	for (i = 0; i < num_bstates(intfc); ++i)
	{
	    bstate = bstate_list(intfc)[i];
	    if (bstate == NULL) continue;
	    s = bstate->_boundary_state_function_name;
	    if (s == NULL) continue;
	    if (strcmp(s,"flowThroughBoundaryState") == 0)
            	bstate->_boundary_state_function = iF_flowThroughBoundaryState;
	}
}	/* end restart_set_dirichlet_bdry_function */

extern void iF_timeDependBoundaryState(
        double          *p0,
        HYPER_SURF      *hs,
        Front           *front,
        POINTER         params,
        POINTER         state)
{
	TIME_DEPENDENT_PARAMS *td_params = (TIME_DEPENDENT_PARAMS*)params;
	STATE *iF_state = (STATE*)state;
	int i,dim = front->rect_grid->dim;
	double time = front->time;
	double *T = td_params->T;
	double omega = td_params->omega;
	double phase = td_params->phase;
	static int step = 0;

	switch (td_params->td_type)
	{
	case CONSTANT:
	    for (i = 0; i < dim; ++i)
	    	iF_state->vel[i] = td_params->v_base[i];
	    iF_state->pres = td_params->p_base;
	    break;
	case PULSE_FUNC:
	    if (time <= T[0])
	    {
	    	for (i = 0; i < dim; ++i)
	    	    iF_state->vel[i] = td_params->v_base[i];
	    	iF_state->pres = td_params->p_base;
	    }
	    else if (time > T[0] && time <= T[1])
	    {
	    	for (i = 0; i < dim; ++i)
	    	    iF_state->vel[i] = td_params->v_base[i] + 
				(time - T[0])/(T[1] - T[0])*
				(td_params->v_peak[i] - td_params->v_base[i]);
	    	iF_state->pres = td_params->p_base + 
				(time - T[0])/(T[1] - T[0])*
				(td_params->p_peak - td_params->p_base);
	    }
	    else if (time > T[1] && time <= T[2])
	    {
	    	for (i = 0; i < dim; ++i)
	    	    iF_state->vel[i] = td_params->v_peak[i];
	    	iF_state->pres = td_params->p_peak;
	    }
	    else if (time > T[2] && time <= T[3])
	    {
	    	for (i = 0; i < dim; ++i)
	    	    iF_state->vel[i] = td_params->v_peak[i] + 
				(time - T[2])/(T[3] - T[2])*
				(td_params->v_tail[i] - td_params->v_peak[i]);
	    	iF_state->pres = td_params->p_peak + 
				(time - T[2])/(T[3] - T[2])*
				(td_params->p_tail - td_params->p_peak);
	    }
	    if (time > T[3])
	    {
	    	for (i = 0; i < dim; ++i)
	    	    iF_state->vel[i] = td_params->v_tail[i];
	    	iF_state->pres = td_params->p_tail;
	    }
	    break;
	case SINE_FUNC:
	    for (i = 0; i < dim; ++i)
	    	iF_state->vel[i] = td_params->v_base[i] +
			td_params->v_amp[i]*sin(omega*time + phase);
	    iF_state->pres = td_params->p_base + td_params->p_amp*
				sin(omega*time + phase);
	    break;
	default:
	    (void) printf("In iF_timeDependBoundaryState(), unknown type!\n");
	    clean_up(ERROR);
	}
	if (debugging("time_depend_bdry"))
	{
	    if (step != front->step)
	    {
	    	printf("time = %f  vel = %f %f  p = %f\n",time,iF_state->vel[0],
				iF_state->vel[1],iF_state->pres);
	    	step = front->step;
	    }
	}
}	/* end iF_timeDependBoundaryState */

extern void iF_flowThroughBoundaryState(
        double          *p0,
        HYPER_SURF      *hs,
        Front           *front,
        POINTER         params,
        POINTER         state)
{
	switch (front->rect_grid->dim)
	{
	case 2:
	    return iF_flowThroughBoundaryState2d(p0,hs,front,params,state);
	case 3:
	    return iF_flowThroughBoundaryState3d(p0,hs,front,params,state);
	}
}	/* end iF_flowThroughBoundaryState */

static void iF_flowThroughBoundaryState3d(
        double          *p0,
        HYPER_SURF      *hs,
        Front           *front,
        POINTER         params,
        POINTER         state)
{
	Tan_stencil **tsten;
	Nor_stencil *nsten;
	FLOW_THROUGH_PARAMS *ft_params = (FLOW_THROUGH_PARAMS*)params;
	POINT *oldp = ft_params->oldp;
	COMPONENT comp = ft_params->comp;
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	IF_FIELD *field = iFparams->field;
	double dir[MAXD];
	double u[3];		/* velocity in the sweeping direction */
	double v[3][MAXD];	/* velocity in the orthogonal direction */
	double pres[3];		/* pressure stencil */
	double f_u;		/* u flux in the sweeping direction */
	double f_v[MAXD];	/* v flux in the orthogonal direction */
	double f_pres;		/* pressure flux */
	double dn,dt = front->dt;
	STATE *oldst,*newst = (STATE*)state;
	STATE  **sts;
	POINTER sl,sr;
	int i,j,k,dim = front->rect_grid->dim;
	int nrad = 2;

	if (debugging("flow_through"))
	    printf("Entering iF_flowThroughBoundaryState3d()\n");

	FT_GetStatesAtPoint(oldp,oldp->hse,oldp->hs,&sl,&sr);
	oldst = NULL;
	if (comp == negative_component(hs))  
	    oldst = (STATE*)sl;
	else 
	    oldst = (STATE*)sr;

	nsten = FT_CreateNormalStencil(front,oldp,comp,nrad);
	for (i = 0; i < dim; ++i)
	    dir[i] = nsten->nor[i];
	dn = FT_GridSizeInDir(dir,front);

	u[1] = 0.0;
	for (j = 0; j < 2; ++j)
	{
	    pres[j] = oldst->pres;
	}
	for (i = 0; i < dim; ++i)
	{
	    double vtmp;
	    FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],
			field->vel[i],getStateVel[i],&vtmp,&oldst->vel[i]);
	    u[1] += vtmp*dir[i];
	    newst->vel[i] = vtmp;
	}
	FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],field->pres,
                            getStatePres,&pres[2],&oldst->pres);

	f_pres = linear_flux(u[1],pres[0],pres[1],pres[2]);
	newst->pres = oldst->pres - dt/dn*f_pres;
	if (debugging("flow_through"))
	{
	    (void) print_Nor_stencil(front,nsten);
	    (void) printf("new velocity after normal prop: %f %f %f\n",
			newst->vel[0],newst->vel[1],newst->vel[2]);
	}
	
	tsten = FrontGetTanStencils(front,oldp,nrad);
	if (tsten == NULL) return;

	for (k = 0; k < dim-1; ++k)
	{
	    for (i = 0; i < dim; ++i)
	    	dir[i] = tsten[k]->dir[i];
	    dn = FT_GridSizeInDir(dir,front);

	    if (comp == negative_component(hs))  
	    	sts = (STATE**)tsten[k]->leftst;
	    else 
	    	sts = (STATE**)tsten[k]->rightst;

	    if (debugging("flow_through"))
	    {
	    	(void) printf("Ambient component: %d\n",comp);
	    	(void) printf("Tangential grid size = %f\n",dn);
		(void) printf("For direction %d\n",k);
	    	(void) print_Tan_stencil(front,tsten[k]);
	    }

	    for (j = 0; j < 3; ++j)
	    	u[j] = 0.0;
	    for (j = 0; j < 3; ++j)
	    {
	    	pres[j] = sts[j-1]->pres;
	    	for (i = 0; i < dim; ++i)
		    u[j] += sts[j-1]->vel[i]*dir[i];
	    	for (i = 0; i < dim; ++i)
	    	{
		    v[j][i] = sts[j-1]->vel[i] - u[j]*dir[i];
	    	}
	    }

	    f_u = burger_flux(u[0],u[1],u[2]);
	    for (i = 0; i < dim; ++i)
	    {
	    	f_v[i] = linear_flux(u[1],v[0][i],v[1][i],v[2][i]);
	    }
	    f_pres = linear_flux(u[1],pres[0],pres[1],pres[2]);

	    for (i = 0; i < dim; ++i)
	    {
	    	newst->vel[i] += - dt/dn*(f_u*dir[i] + f_v[i]) ;
	    }
	    newst->pres += - dt/dn*f_pres;
	}
	if (debugging("flow_through"))
	{
	    (void) printf("State after tangential sweep:\n");
	    (void) print_general_vector("Velocity: ",newst->vel,dim,"\n");
	    (void) printf("Pressure: %f\n",newst->pres);
	}
}	/* end iF_flowThroughBoundaryState3d */

static void iF_flowThroughBoundaryState2d(
        double          *p0,
        HYPER_SURF      *hs,
        Front           *front,
        POINTER         params,
        POINTER         state)
{
	Tan_stencil **tsten;
	Nor_stencil *nsten;
	FLOW_THROUGH_PARAMS *ft_params = (FLOW_THROUGH_PARAMS*)params;
	POINT *oldp = ft_params->oldp;
	COMPONENT comp = ft_params->comp;
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	IF_FIELD *field = iFparams->field;
	double dir[MAXD];
	double u[3];		/* velocity in the sweeping direction */
	double v[3][MAXD];	/* velocity in the orthogonal direction */
	double vort[3];		/* vorticity stencil */
	double pres[3];		/* pressure stencil */
	double f_u;		/* u flux in the sweeping direction */
	double f_v[MAXD];	/* v flux in the orthogonal direction */
	double f_vort;		/* vort flux */
	double f_pres;		/* pressure flux */
	double dn,dt = front->dt;
	STATE *oldst,*newst = (STATE*)state;
	STATE  **sts;
	POINTER sl,sr;
	int i,j,dim = front->rect_grid->dim;
	int nrad = 2;
	
	if (debugging("flow_through"))
	    printf("Entering iF_flowThroughBoundaryState2d()\n");

	FT_GetStatesAtPoint(oldp,oldp->hse,oldp->hs,&sl,&sr);
	nsten = FT_CreateNormalStencil(front,oldp,comp,nrad);
	dn = FT_GridSizeInDir(nsten->nor,front);
	if (debugging("flow_through"))
	{
	    (void) printf("Normal grid size = %f\n",dn);
	    (void) print_Nor_stencil(front,nsten);
	}

	if (comp == negative_component(hs))  
	    oldst = (STATE*)sl;
	else 
	    oldst = (STATE*)sr;

	u[1] = 0.0;
	for (j = 0; j < 2; ++j)
	{
	    vort[j] = oldst->vort;
	    pres[j] = oldst->pres;
	}
	for (i = 0; i < dim; ++i)
	{
	    double vtmp;
	    FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],
			field->vel[i],getStateVel[i],&vtmp,&oldst->vel[i]);
	    u[1] += vtmp*dir[i];
	    newst->vel[i] = vtmp;
	}

	FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],field->vort,
                            getStateVort,&vort[2],&oldst->vort);
	FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],field->pres,
                            getStatePres,&pres[2],&oldst->pres);

	f_vort = linear_flux(u[1],vort[0],vort[1],vort[2]);
	f_pres = linear_flux(u[1],pres[0],pres[1],pres[2]);

	newst->vort = oldst->vort - dt/dn*f_vort;
	newst->pres = oldst->pres - dt/dn*f_pres;
	
	tsten = FrontGetTanStencils(front,oldp,nrad);

	if (debugging("flow_through"))
	{
	    (void) printf("Ambient component: %d\n",comp);
	    (void) printf("Tangential grid size = %f\n",dn);
	    (void) print_Tan_stencil(front,tsten[0]);
	}

	if (comp == negative_component(hs))  
	    sts = (STATE**)tsten[0]->leftst;
	else 
	    sts = (STATE**)tsten[0]->rightst;

	for (i = 0; i < dim; ++i)
	    dir[i] = tsten[0]->dir[i];
	dn = FT_GridSizeInDir(dir,front);

	for (j = 0; j < 3; ++j)
	    u[j] = 0.0;
	for (j = 0; j < 3; ++j)
	{
	    vort[j] = sts[j-1]->vort;
	    pres[j] = sts[j-1]->pres;
	    for (i = 0; i < dim; ++i)
	    {
		u[j] += sts[j-1]->vel[i]*dir[i];
	    }
	    for (i = 0; i < dim; ++i)
	    {
		v[j][i] = sts[j-1]->vel[i] - dir[i]*u[j];
	    }
	}

	f_u = burger_flux(u[0],u[1],u[2]);
	for (i = 0; i < dim; ++i)
	    f_v[i] = linear_flux(u[1],v[0][i],v[1][i],v[2][i]);
	f_vort = linear_flux(u[1],vort[0],vort[1],vort[2]);
	f_pres = linear_flux(u[1],pres[0],pres[1],pres[2]);

	for (i = 0; i < dim; ++i)
	    newst->vel[i] += - dt/dn*(f_u*dir[i] + f_v[i]) ;
	newst->vort += - dt/dn*f_vort;
	newst->pres += - dt/dn*f_pres;
	
	if (debugging("flow_through"))
	{
	    (void) printf("State after tangential sweep:\n");
	    (void) print_general_vector("Velocity: ",newst->vel,dim,"\n");
	    (void) printf("Vorticity: %f\n",newst->vort);
	    (void) printf("Pressure: %f\n",newst->pres);
	}
}       /* end iF_flowThroughBoundaryState2d */

extern void ifluid_point_propagate(
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
	case MOVABLE_BODY_BOUNDARY:
	case ICE_PARTICLE_BOUNDARY:
	    return rgbody_point_propagate(front,wave,oldp,newp,oldhse,
					oldhs,dt,V);
	case NEUMANN_BOUNDARY:
	case GROWING_BODY_BOUNDARY:
	    return neumann_point_propagate(front,wave,oldp,newp,oldhse,
					oldhs,dt,V);
	case DIRICHLET_BOUNDARY:
	    return dirichlet_point_propagate(front,wave,oldp,newp,oldhse,
					oldhs,dt,V);
	default:
	    return contact_point_propagate(front,wave,oldp,newp,oldhse,
					oldhs,dt,V);
	}
}       /* ifluid_point_propagate */

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
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	IF_FIELD *field = iFparams->field;
        int i, dim = front->rect_grid->dim;
	double *m_pre = field->pres;
	double *m_phi = field->phi;
	double *m_vor = field->vort;
	STATE *oldst,*newst;
	POINTER sl,sr;
	COMPONENT comp;
	double nor[MAXD],p1[MAXD];
	double *p0 = Coords(oldp);
	double dn,*h = front->rect_grid->h;

	FT_GetStatesAtPoint(oldp,oldhse,oldhs,&sl,&sr);
	if (ifluid_comp(negative_component(oldhs)))
	{
	    comp = negative_component(oldhs);
	    oldst = (STATE*)sl;
	    newst = (STATE*)left_state(newp);
	}
	else if (ifluid_comp(positive_component(oldhs)))
	{
	    comp = positive_component(oldhs);
	    oldst = (STATE*)sr;
	    newst = (STATE*)right_state(newp);
	}
	FT_NormalAtPoint(oldp,front,nor,comp);

	dn = grid_size_in_direction(nor,h,dim);
	for (i = 0; i < dim; ++i)
	    p1[i] = p0[i] + nor[i]*dn;

	for (i = 0; i < dim; ++i)
	{
            Coords(newp)[i] = Coords(oldp)[i];
	    newst->vel[i] = 0.0;
            FT_RecordMaxFrontSpeed(i,0.0,NULL,Coords(newp),front);
	}
	FT_IntrpStateVarAtCoords(front,comp,p1,m_pre,
			getStatePres,&newst->pres,&oldst->pres);
	/*
	FT_IntrpStateVarAtCoords(front,comp,p1,m_phi,
			getStatePhi,&newst->phi,&oldst->phi);
	*/
	if (dim == 2)
	{
	    FT_IntrpStateVarAtCoords(front,comp,p1,m_vor,
			getStateVort,&newst->vort,&oldst->vort);
	}
	FT_RecordMaxFrontSpeed(dim,0.0,NULL,Coords(newp),front);
        return;
}	/* end neumann_point_propagate */

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
        double speed;
        int i, dim = front->rect_grid->dim;
	STATE *newst = NULL;
	STATE *bstate;
	COMPONENT comp;

	if (debugging("dirichlet_bdry"))
	{
	    (void) printf("Entering dirichlet_point_propagate()\n");
	    (void) print_general_vector("oldp:  ",Coords(oldp),dim,"\n");
	}

	if (ifluid_comp(negative_component(oldhs)))
	{
	    newst = (STATE*)left_state(newp);
	    comp = negative_component(oldhs);
	}
	else if (ifluid_comp(positive_component(oldhs)))
	{
	    newst = (STATE*)right_state(newp);
	    comp = positive_component(oldhs);
	}
	if (newst == NULL) return;	// node point

	if (boundary_state(oldhs) != NULL)
	{
	    bstate = (STATE*)boundary_state(oldhs);
            for (i = 0; i < dim; ++i)
	    {
	    	newst->vel[i] = bstate->vel[i];
		FT_RecordMaxFrontSpeed(i,fabs(newst->vel[i]),NULL,Coords(newp),
					front);
	    }
	    speed = mag_vector(newst->vel,dim);
	    FT_RecordMaxFrontSpeed(dim,speed,NULL,Coords(newp),front);
            newst->pres = bstate->pres;
	    newst->phi = bstate->phi;
            newst->vort = 0.0;

	    if (debugging("dirichlet_bdry"))
	    {
		(void) printf("Preset boundary state:\n");
		(void) print_general_vector("Velocity: ",newst->vel,dim,"\n");
		(void) printf("Pressure: %f\n",newst->pres);
		(void) printf("Vorticity: %f\n",newst->vort);
	    }
	}
	else if (boundary_state_function(oldhs))
	{
	    if (strcmp(boundary_state_function_name(oldhs),
		       "flowThroughBoundaryState") == 0)
	    {
		FLOW_THROUGH_PARAMS ft_params;
		oldp->hse = oldhse;
		oldp->hs = oldhs;
	    	ft_params.oldp = oldp;
	    	ft_params.comp = comp;
	    	(*boundary_state_function(oldhs))(Coords(oldp),oldhs,front,
			(POINTER)&ft_params,(POINTER)newst);	
	    }
	    else if (strcmp(boundary_state_function_name(oldhs),
		       "iF_timeDependBoundaryState") == 0)
	    {
		TIME_DEPENDENT_PARAMS *td_params = (TIME_DEPENDENT_PARAMS*)
				boundary_state_function_params(oldhs);
	    	(*boundary_state_function(oldhs))(Coords(oldp),oldhs,front,
			(POINTER)td_params,(POINTER)newst);	
	    }
            for (i = 0; i < dim; ++i)
		FT_RecordMaxFrontSpeed(i,fabs(newst->vel[i]),NULL,Coords(newp),
					front);
	    speed = mag_vector(newst->vel,dim);
	    FT_RecordMaxFrontSpeed(dim,speed,NULL,Coords(newp),front);
	}
	if (debugging("dirichlet_bdry"))
	    (void) printf("Leaving dirichlet_point_propagate()\n");
        return;
}	/* end dirichlet_point_propagate */

static  void contact_point_propagate(
        Front *front,
        POINTER wave,
        POINT *oldp,
        POINT *newp,
        HYPER_SURF_ELEMENT *oldhse,
        HYPER_SURF         *oldhs,
        double              dt,
        double              *V)
{
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	IF_FIELD *field = iFparams->field;
        double vel[MAXD],s;
        int i, dim = front->rect_grid->dim;
	double *m_pre = field->pres;
	double *m_vor = field->vort;
	double *p0;
	STATE *oldst,*newst;
	POINTER sl,sr;
	double pres,vort;

        (*front->vfunc)(front->vparams,front,oldp,oldhse,oldhs,vel);
        for (i = 0; i < dim; ++i)
        {
            Coords(newp)[i] = Coords(oldp)[i] + dt*vel[i];
            FT_RecordMaxFrontSpeed(i,fabs(vel[i]),NULL,Coords(newp),front);
        }

	FT_GetStatesAtPoint(oldp,oldhse,oldhs,&sl,&sr);
	oldst = (STATE*)sl;
	p0 = Coords(newp);
	FT_IntrpStateVarAtCoords(front,-1,p0,m_pre,getStatePres,&pres,
				&oldst->pres);
	if (dim == 2)
	{
	    FT_IntrpStateVarAtCoords(front,-1,p0,m_vor,getStateVort,&vort,
				&oldst->vort);
	}

	newst = (STATE*)left_state(newp);
	newst->vort = vort;
	newst->pres = pres;
	for (i = 0; i < dim; ++i)
	{
	    newst->vel[i] = vel[i];
	}
	newst = (STATE*)right_state(newp);
	newst->vort = vort;
	newst->pres = pres;
	for (i = 0; i < dim; ++i)
	{
	    newst->vel[i] = vel[i];
	}

	s = mag_vector(vel,dim);
	FT_RecordMaxFrontSpeed(dim,s,NULL,Coords(newp),front);
}	/* end contact_point_propagate */

static  void rgbody_point_propagate(
        Front *front,
        POINTER wave,
        POINT *oldp,
        POINT *newp,
        HYPER_SURF_ELEMENT *oldhse,
        HYPER_SURF         *oldhs,
        double              dt,
        double              *V)
{
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	IF_FIELD *field = iFparams->field;
        double vel[MAXD];
        int i, dim = front->rect_grid->dim;
	double dn,*h = front->rect_grid->h;
	double *m_pre = field->pres;
	double *m_vor = field->vort;
	double nor[MAXD],p1[MAXD];
	double *p0 = Coords(oldp);
	STATE *oldst,*newst;
	POINTER sl,sr;
	COMPONENT comp;

	FT_GetStatesAtPoint(oldp,oldhse,oldhs,&sl,&sr);
	if (ifluid_comp(negative_component(oldhs)))
	{
	    comp = negative_component(oldhs);
	    oldst = (STATE*)sl;
	    newst = (STATE*)left_state(newp);
	}
	else if (ifluid_comp(positive_component(oldhs)))
	{
	    comp = positive_component(oldhs);
	    oldst = (STATE*)sr;
	    newst = (STATE*)right_state(newp);
	}
	FT_NormalAtPoint(oldp,front,nor,comp);

	dn = grid_size_in_direction(nor,h,dim);
	for (i = 0; i < dim; ++i)
	    p1[i] = p0[i] + nor[i]*dn;

        if (wave_type(oldhs) == MOVABLE_BODY_BOUNDARY)
        {
            double omega_dt,crds_com[MAXD];
            omega_dt = angular_velo(oldhs)*dt;
            for (i = 0; i < dim; ++i)
            {
                vel[i] = center_of_mass_velo(oldhs)[i];
                crds_com[i] = Coords(oldp)[i] +dt*vel[i] -
                        center_of_mass(oldhs)[i];
            }
            vel[0] += -angular_velo(oldhs)*crds_com[1]*cos(omega_dt) -
                     angular_velo(oldhs)*crds_com[0]*sin(omega_dt);
            vel[1] +=  angular_velo(oldhs)*crds_com[0]*cos(omega_dt) -
                     angular_velo(oldhs)*crds_com[1]*sin(omega_dt);
            for (i = 0; i < dim; ++i)
            {
                Coords(newp)[i] = Coords(oldp)[i] + dt*vel[i];
                newst->vel[i] = vel[i];
                FT_RecordMaxFrontSpeed(i,fabs(vel[i]),NULL,Coords(newp),front);
            }
        }
        else
        {
            fourth_order_point_propagate(front,NULL,oldp,newp,oldhse,
                                    oldhs,dt,vel);
        }
	for (i = 0; i < dim; ++i) newst->vel[i] = vel[i];
	FT_IntrpStateVarAtCoords(front,comp,p1,m_pre,
			getStatePres,&newst->pres,&oldst->pres);
	if (dim == 2)
	{
	    FT_IntrpStateVarAtCoords(front,comp,p1,m_vor,
			getStateVort,&newst->vort,&oldst->vort);
	}
        return;
}	/* end rgbody_point_propagate */

extern void fluid_print_front_states(
	FILE *outfile,
	Front *front)
{
	INTERFACE *intfc = front->interf;
        STATE *sl,*sr;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
	int dim = intfc->dim;

	fprintf(outfile,"Interface ifluid states:\n");
        next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
            FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
            fprintf(outfile,"%24.18g %24.18g\n",getStatePres(sl),
                                getStatePres(sr));
            fprintf(outfile,"%24.18g %24.18g\n",getStatePhi(sl),
                                getStatePhi(sr));
            if (dim == 2)
            {
                fprintf(outfile,"%24.18g %24.18g\n",getStateXvel(sl),
                                getStateXvel(sr));
                fprintf(outfile,"%24.18g %24.18g\n",getStateYvel(sl),
                                getStateYvel(sr));
                fprintf(outfile,"%24.18g %24.18g\n",getStateVort(sl),
                                getStateVort(sr));
                fprintf(outfile,"%24.18g %24.18g\n",getStateXimp(sl),
                                getStateXimp(sr));
                fprintf(outfile,"%24.18g %24.18g\n",getStateYimp(sl),
                                getStateYimp(sr));
            }
            if (dim == 3)
            {
                fprintf(outfile,"%24.18g %24.18g\n",getStateXvel(sl),
                                getStateXvel(sr));
                fprintf(outfile,"%24.18g %24.18g\n",getStateYvel(sl),
                                getStateYvel(sr));
                fprintf(outfile,"%24.18g %24.18g\n",getStateZvel(sl),
                                getStateZvel(sr));
                fprintf(outfile,"%24.18g %24.18g\n",getStateXimp(sl),
                                getStateXimp(sr));
                fprintf(outfile,"%24.18g %24.18g\n",getStateYimp(sl),
                                getStateYimp(sr));
                fprintf(outfile,"%24.18g %24.18g\n",getStateZimp(sl),
                                getStateZimp(sr));
            }
        }
}	/* end fluid_print_front_states */

extern void fluid_read_front_states(
	FILE *infile,
	Front *front)
{
	INTERFACE *intfc = front->interf;
        STATE *sl,*sr;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
        STATE *lstate,*rstate;
	int dim = intfc->dim;

	next_output_line_containing_string(infile,"Interface ifluid states:");
        next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
            FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
            lstate = (STATE*)sl;        rstate = (STATE*)sr;
            fscanf(infile,"%lf %lf",&lstate->pres,&rstate->pres);
            fscanf(infile,"%lf %lf",&lstate->phi,&rstate->phi);
            fscanf(infile,"%lf %lf",&lstate->vel[0],&rstate->vel[0]);
            fscanf(infile,"%lf %lf",&lstate->vel[1],&rstate->vel[1]);
            if (dim == 2)
                fscanf(infile,"%lf %lf",&lstate->vort,&rstate->vort);
            if (dim == 3)
                fscanf(infile,"%lf %lf",&lstate->vel[2],&rstate->vel[2]);
            fscanf(infile,"%lf %lf",&lstate->impulse[0],&rstate->impulse[0]);
            fscanf(infile,"%lf %lf",&lstate->impulse[1],&rstate->impulse[1]);
            if (dim == 3)
            	fscanf(infile,"%lf %lf",&lstate->impulse[2],&rstate->impulse[2]);
        }
}	/* end fluid_read_front_states */

extern void read_iFparams(
	char *inname,
	IF_PARAMS *iFparams)
{
	char string[100];
	FILE *infile = fopen(inname,"r");
	int i,dim = iFparams->dim;

	/* defaults numerical schemes */
	iFparams->num_scheme.projc_method = SIMPLE;
	iFparams->num_scheme.advec_method = WENO;
	iFparams->num_scheme.ellip_method = SIMPLE_ELLIP;

	CursorAfterString(infile,"Enter projection type:");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	switch (string[0])
	{
	case 'S':
	case 's':
	    iFparams->num_scheme.projc_method = SIMPLE;
	    break;
	case 'B':
	case 'b':
	    iFparams->num_scheme.projc_method = BELL_COLELLA;
	    break;
	case 'K':
	case 'k':
	    iFparams->num_scheme.projc_method = KIM_MOIN;
	    break;
	case 'P':
	case 'p':
	    iFparams->num_scheme.projc_method = PEROT_BOTELLA;
	}
	assert(iFparams->num_scheme.projc_method != ERROR_PROJC_SCHEME);
	(void) printf("The default advection order is WENO-Runge-Kutta 4\n");
	iFparams->adv_order = 4;
	if (CursorAfterStringOpt(infile,"Enter advection order:"))
	{
	    fscanf(infile,"%d",&iFparams->adv_order);
	    (void) printf("%d\n",iFparams->adv_order);
	}

	(void) printf("Available elliptic methods are:\n");
	(void) printf("\tSimple elliptic (S)\n");
	(void) printf("\tCIM elliptic (C)\n");
	(void) printf("\tDouble elliptic (DB)\n");
	(void) printf("\tDual elliptic (DU)\n");
	if (CursorAfterStringOpt(infile,"Enter elliptic method:"))
	{
	    fscanf(infile,"%s",string);
	    (void) printf("%s\n",string);
	    switch (string[0])
	    {
	    case 'S':
	    case 's':
	    	iFparams->num_scheme.ellip_method = SIMPLE_ELLIP;
	    	break;
	    case 'c':
	    case 'C':
	    	iFparams->num_scheme.ellip_method = CIM_ELLIP;
	    	break;
	    case 'd':
	    case 'D':
		if (string[1] == 'b' || string[1] == 'B')
	    	    iFparams->num_scheme.ellip_method = DOUBLE_ELLIP;
		else if (string[1] == 'u' || string[1] == 'U')
	    	    iFparams->num_scheme.ellip_method = DUAL_ELLIP;
	    	break;
	    }
	}

	for (i = 0; i < dim; ++i) iFparams->U_ambient[i] = 0.0;
        if (CursorAfterStringOpt(infile,"Enter fluid ambient velocity:"))
        {
            for (i = 0; i < dim; ++i)
            {
                fscanf(infile,"%lf ",&iFparams->U_ambient[i]);
                (void) printf("%f ",iFparams->U_ambient[i]);
            }
            (void) printf("\n");
        }
	iFparams->ub_speed = HUGE;
        if (CursorAfterStringOpt(infile,"Enter upper bound for speed:"))
	{
            fscanf(infile,"%lf ",&iFparams->ub_speed);
            (void) printf("%f\n",iFparams->ub_speed);
	}
	iFparams->total_div_cancellation = NO;
        if (CursorAfterStringOpt(infile,	
		"Enter yes to use total divergence cancellation:"))
	{
	    fscanf(infile,"%s",string);
	    (void) printf("%s\n",string);
	    if (string[0] == 'y' || string[0] == 'Y')
	    	iFparams->total_div_cancellation = YES;
	}
        if (CursorAfterStringOpt(infile,
		"Enter density and viscosity of the fluid:"))
        {
            fscanf(infile,"%lf %lf",&iFparams->rho2,&iFparams->mu2);
            (void) printf("%f %f\n",iFparams->rho2,iFparams->mu2);
	}
	iFparams->use_eddy_visc = NO;
        if (CursorAfterStringOpt(infile,
		"Enter yes to use eddy viscosity:"))
        {
	    fscanf(infile,"%s",string);
	    (void) printf("%s\n",string);
	    if (string[0] == 'y' || string[0] == 'Y')
	    {
	    	iFparams->use_eddy_visc = YES;
        	CursorAfterString(infile,
			"Enter maximum distance for eddy viscosity:");
            	fscanf(infile,"%lf",&iFparams->ymax);
            	(void) printf("%f\n",iFparams->ymax);
	    }
	}
        if (CursorAfterStringOpt(infile,"Enter gravity:"))
        {
            for (i = 0; i < dim; ++i)
            {
                fscanf(infile,"%lf ",&iFparams->gravity[i]);
                (void) printf("%f ",iFparams->gravity[i]);
            }
            (void) printf("\n");
        }

	fclose(infile);
}	/* end read_iFparams */

extern boolean isDirichletPresetBdry(
	Front *front,
	int *icoords,
	GRID_DIRECTION dir,
	COMPONENT comp)
{
	HYPER_SURF *hs;
	POINTER intfc_state;
	double crx_coords[MAXD];
	INTERFACE *grid_intfc = front->grid_intfc;

	if (!FT_StateStructAtGridCrossing(front,grid_intfc,icoords,dir,
                                comp,&intfc_state,&hs,crx_coords))
	    return NO;
	if (wave_type(hs) != DIRICHLET_BOUNDARY)
	    return NO;
	if (boundary_state(hs) != NULL)
	    return NO;
	return YES;
}	/* end isDirichletPresetBdry */

static void get_time_dependent_params(
	int dim,
	FILE *infile,
	POINTER *params)
{
	static TIME_DEPENDENT_PARAMS *td_params;
	char string[100];
	int i;

	FT_ScalarMemoryAlloc((POINTER*)&td_params,
			sizeof(TIME_DEPENDENT_PARAMS));
	CursorAfterString(infile,"Enter type of time-dependent function:");
	fscanf(infile,"%s",string);
	(void) printf(" %s\n",string);
	switch (string[0])
	{
	case 'C':
	case 'c':
	    td_params->td_type = CONSTANT;
	    CursorAfterString(infile,"Enter base velocity:");
	    for (i = 0; i < dim; ++i)
	    {
		fscanf(infile,"%lf ",&td_params->v_base[i]);
		(void) printf("%f ",td_params->v_base[i]);
	    }
	    (void) printf("\n");
	    CursorAfterString(infile,"Enter base pressure:");
	    fscanf(infile,"%lf ",&td_params->p_base);
	    (void) printf("%f\n",td_params->p_base);
	    break;
	case 'P':
	case 'p':
	    td_params->td_type = PULSE_FUNC;
	    CursorAfterString(infile,"Enter base velocity:");
	    for (i = 0; i < dim; ++i)
	    {
		fscanf(infile,"%lf ",&td_params->v_base[i]);
		(void) printf("%f ",td_params->v_base[i]);
	    }
	    (void) printf("\n");
	    CursorAfterString(infile,"Enter base pressure:");
	    fscanf(infile,"%lf ",&td_params->p_base);
	    (void) printf("%f\n",td_params->p_base);
	    CursorAfterString(infile,"Enter peak velocity:");
	    for (i = 0; i < dim; ++i)
	    {
		fscanf(infile,"%lf ",&td_params->v_peak[i]);
		(void) printf("%f ",td_params->v_peak[i]);
	    }
	    (void) printf("\n");
	    CursorAfterString(infile,"Enter peak pressure:");
	    fscanf(infile,"%lf ",&td_params->p_peak);
	    (void) printf("%f\n",td_params->p_peak);
	    CursorAfterString(infile,"Enter tail velocity:");
	    for (i = 0; i < dim; ++i)
	    {
		fscanf(infile,"%lf ",&td_params->v_tail[i]);
		(void) printf("%f ",td_params->v_tail[i]);
	    }
	    (void) printf("\n");
	    CursorAfterString(infile,"Enter tail pressure:");
	    fscanf(infile,"%lf ",&td_params->p_tail);
	    (void) printf("%f\n",td_params->p_tail);
	    CursorAfterString(infile,"Enter time to rise:");
	    fscanf(infile,"%lf ",&td_params->T[0]);
	    (void) printf("%f\n",td_params->T[0]);
	    CursorAfterString(infile,"Enter time to reach peak:");
	    fscanf(infile,"%lf ",&td_params->T[1]);
	    (void) printf("%f\n",td_params->T[1]);
	    CursorAfterString(infile,"Enter time to fall:");
	    fscanf(infile,"%lf ",&td_params->T[2]);
	    (void) printf("%f\n",td_params->T[2]);
	    CursorAfterString(infile,"Enter time to reach tail:");
	    fscanf(infile,"%lf ",&td_params->T[3]);
	    (void) printf("%f\n",td_params->T[3]);
	    break;
	case 'S':
	case 's':
	    td_params->td_type = SINE_FUNC;
	    CursorAfterString(infile,"Enter base velocity:");
	    for (i = 0; i < dim; ++i)
	    {
		fscanf(infile,"%lf ",&td_params->v_base[i]);
		(void) printf("%f ",td_params->v_base[i]);
	    }
	    (void) printf("\n");
	    CursorAfterString(infile,"Enter base pressure:");
	    fscanf(infile,"%lf ",&td_params->p_base);
	    (void) printf("%f\n",td_params->p_base);
	    CursorAfterString(infile,"Enter velocity amplitude:");
	    for (i = 0; i < dim; ++i)
	    {
		fscanf(infile,"%lf ",&td_params->v_amp[i]);
		(void) printf("%f ",td_params->v_amp[i]);
	    }
	    (void) printf("\n");
	    CursorAfterString(infile,"Enter pressure amplitude:");
	    fscanf(infile,"%lf ",&td_params->p_amp);
	    (void) printf("%f\n",td_params->p_amp);
	    CursorAfterString(infile,"Enter oscillation period:");
	    fscanf(infile,"%lf ",&td_params->omega);
	    (void) printf("%f\n",td_params->omega);
	    td_params->omega = 2.0*PI/td_params->omega;
	    CursorAfterString(infile,"Enter initial phase:");
	    fscanf(infile,"%lf ",&td_params->phase);
	    (void) printf("%f\n",td_params->phase);
	    td_params->phase *= PI/180.0;
	    break;
	default:
	    (void) printf("Unknown type of time-dependent function!\n");
	    clean_up(ERROR);
	}

	*params = (POINTER)td_params;	
}	/* end get_time_dependent_params */

extern void recordBdryEnergyFlux(
	Front *front,
	char *out_name)
{
	RECT_GRID *rgr = front->rect_grid;
	INTERFACE *intfc = front->interf;
	HYPER_SURF *hs;
	double Energy_in,Energy_out;
	int dir,side;
	int dim = rgr->dim;
	static FILE *ein_file,*eout_file;
	char file_name[100];

	if (ein_file == NULL && pp_mynode() == 0)
	{
	    sprintf(file_name,"%s/in_energy.xg",out_name);
	    ein_file = fopen(file_name,"w");
	    fprintf(ein_file,"\"Energy influx\" vs. time\n");
	    sprintf(file_name,"%s/out_energy.xg",out_name);
	    eout_file = fopen(file_name,"w");
	    fprintf(eout_file,"\"Energy outflux\" vs. time\n");
	}
	Energy_in = Energy_out = 0.0;
	for (dir = 0; dir < dim; ++dir)
	for (side = 0; side < 2; ++side)
	{
	    hs = FT_RectBoundaryHypSurf(intfc,DIRICHLET_BOUNDARY,dir,side);
	    if (hs == NULL) continue;
	    if (boundary_state(Hyper_surf(hs)))
	    {
		addToEnergyFlux(rgr,hs,&Energy_in,&Energy_out,dir,side,YES);
	    }
	    else if (boundary_state_function(Hyper_surf(hs)))
	    {
		if (strcmp(boundary_state_function_name(hs),
                       "iF_timeDependBoundaryState") == 0)
		{
		    addToEnergyFlux(rgr,hs,&Energy_in,&Energy_out,dir,side,YES);
		}
		else if (strcmp(boundary_state_function_name(hs),
                       "flowThroughBoundaryState") == 0)
		{
		    addToEnergyFlux(rgr,hs,&Energy_in,&Energy_out,dir,side,NO);
		}
	    }
	}
	pp_global_sum(&Energy_in,1);
	pp_global_sum(&Energy_out,1);

	if (pp_mynode() == 0)
	{
	    fprintf(ein_file,"%f %f\n",front->time,Energy_in);
	    fprintf(eout_file,"%f %f\n",front->time,Energy_out);
	}
}	/* end recordBdryEnergyFlux */


static void addToEnergyFlux(
	RECT_GRID *rgr,
	HYPER_SURF *hs,
	double *Energy_in,
	double *Energy_out,
	int dir,
	int side,
	boolean is_influx)
{
	int i,dim = rgr->dim;
	double *L = rgr->L;	
	double *U = rgr->U;	
	CURVE *c;
	SURFACE *s;
	BOND *b = NULL;
	TRI *t;
	double ave_coord,engy_flux,vel;
	boolean is_outside_hse,is_left_side;
	STATE *sl,*sr,*state;
	POINT *pts[MAXD];

	is_left_side = (ifluid_comp(negative_component(hs))) ? YES : NO;
	switch (dim)
	{
	case 2:
	    c = Curve_of_hs(hs);
	    for (b = c->first; b != NULL; b = b->next)
	    {
		is_outside_hse = NO;
		pts[0] = b->start;
		pts[1] = b->end;
		for (i = 1; i < dim; ++i)
		{
		    ave_coord = (Coords(pts[0])[(dir+i)%dim] +
				 Coords(pts[1])[(dir+i)%dim])/2.0;
		    if (ave_coord < L[(dir+i)%dim] ||
			ave_coord > U[(dir+i)%dim])
			is_outside_hse = YES;
		}
		if (is_outside_hse) continue;
		engy_flux = 0.0;
		for (i = 0; i < 2; ++i)
		{
		    FT_GetStatesAtPoint(pts[i],Hyper_surf_element(b),hs,
				(POINTER*)&sl,(POINTER*)&sr);
		    state = (is_left_side) ? sl : sr;
		    if (is_influx)
		    	vel = (side == 0) ? state->vel[dir] : -state->vel[dir];
		    else
		    	vel = (side == 1) ? state->vel[dir] : -state->vel[dir];
		    //engy_flux += 0.5*state->dens*(sqr(state->vel[0]) +
		    engy_flux += 0.5*(sqr(state->vel[0]) +
					sqr(state->vel[1]))*vel;
		}
		engy_flux *= bond_length(b)/2.0;
		if (is_influx)
		    *Energy_in += engy_flux;
		else
		    *Energy_out += engy_flux;
	    }
	    break;
	case 3:
	    s = Surface_of_hs(hs);
	    for (t = first_tri(s); !at_end_of_tri_list(t,s); t = t->next)
	    {
		is_outside_hse = NO;
		for (i = 0; i < 3; ++i)
		    pts[i] = Point_of_tri(t)[i];
		for (i = 1; i < dim; ++i)
		{
		    ave_coord = (Coords(pts[0])[(dir+i)%dim] +
				 Coords(pts[1])[(dir+i)%dim] +
				 Coords(pts[2])[(dir+i)%dim])/3.0;
		    if (ave_coord < L[(dir+i)%dim] ||
			ave_coord > U[(dir+i)%dim])
			is_outside_hse = YES;
		}
		if (is_outside_hse) continue;
		engy_flux = 0.0;
		for (i = 0; i < 3; ++i)
		{
		    FT_GetStatesAtPoint(pts[i],Hyper_surf_element(t),hs,
				(POINTER*)&sl,(POINTER*)&sr);
		    state = (is_left_side) ? sl : sr;
		    if (is_influx)
		    	vel = (side == 0) ? state->vel[dir] : -state->vel[dir];
		    else
		    	vel = (side == 1) ? state->vel[dir] : -state->vel[dir];
		    //engy_flux += 0.5*state->dens*(sqr(state->vel[0]) +
		    engy_flux += 0.5*(sqr(state->vel[0]) +
					sqr(state->vel[1]) +
					sqr(state->vel[2]))*vel;
		}
		engy_flux *= bond_length(b)/3.0;
		if (is_influx)
		    *Energy_in += engy_flux;
		else
		    *Energy_out += engy_flux;
	    }
	    break;
	}
}	/* end addToEnergyFlux */

extern double p_jump(
	POINTER params,
	int D,
	double *coords)
{
	return 0.0;
}	/* end p_jump */

extern double grad_p_jump_n(
	POINTER params,
	int D,
	double *N,
	double *coords)
{
	return 0.0;
}	/* end grad_p_jump_n */

extern double grad_p_jump_t(
	POINTER params,
	int D,
	int i,
	double *N,
	double *coords)
{
	return 0.0;
}	/* end grad_p_jump_t */

extern double getPhiFromPres(
        Front *front,
        double pres)
{
        IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
        switch (iFparams->num_scheme.projc_method)
        {
        case BELL_COLELLA:
            return 0.0;
        case KIM_MOIN:
            return 0.0;
        case SIMPLE:
        case PEROT_BOTELLA:
            return pres;
        default:
            (void) printf("Unknown projection type\n");
            clean_up(0);
        }
}       /* end getPhiFromPres */

extern double getPressure(
        Front *front,
        double *coords,
        double *base_coords)
{
        INTERFACE *intfc = front->interf;
        int i,dim = Dimension(intfc);
        POINT *p0;
        double pres,pres0;
        IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
        double *g = iFparams->gravity;
        double rho = iFparams->rho2;
        boolean hyper_surf_found = NO;

        return 0.0;
        pres0 = 1.0;
        if (dim == 2)
        {
            CURVE **c;
            for (c = intfc->curves; c && *c; ++c)
            {
                if (wave_type(*c) == DIRICHLET_BOUNDARY &&
                    boundary_state(*c) != NULL)
                {
                    p0 = (*c)->first->start;
                    pres0 = getStatePres(boundary_state(*c));
                    hyper_surf_found = YES;
                }
            }
        }
        else if (dim == 3)
        {
            SURFACE **s;
            for (s = intfc->surfaces; s && *s; ++s)
            {
                if (wave_type(*s) == DIRICHLET_BOUNDARY &&
                    boundary_state(*s) != NULL)
                {
                    p0 = Point_of_tri(first_tri(*s))[0];
                    pres0 = getStatePres(boundary_state(*s));
                    hyper_surf_found = YES;
                }
            }
        }
        pres = pres0;
        if (hyper_surf_found)
        {
            for (i = 0; i < dim; ++i)
                pres -= rho*(coords[i] - Coords(p0)[i])*g[i];
        }
        else if (base_coords != NULL)
        {
            for (i = 0; i < dim; ++i)
                pres -= rho*(coords[i] - Coords(p0)[i])*g[i];
        }
        return pres;
}       /* end getPressure */

extern int ifluid_find_state_at_dual_crossing(
	Front *front,
	int *icoords,
	GRID_DIRECTION dir,
	int comp,
	POINTER *state,
	HYPER_SURF **hs,
	double *crx_coords)
{
	boolean status;
	INTERFACE *grid_intfc = front->comp_grid_intfc;
	status = FT_StateStructAtGridCrossing(front,grid_intfc,icoords,dir,
				comp,state,hs,crx_coords);
	if (status == NO) 
	    return NO_PDE_BOUNDARY;
	if (wave_type(*hs) == FIRST_PHYSICS_WAVE_TYPE) 
	    return NO_PDE_BOUNDARY;
	if (wave_type(*hs) == NEUMANN_BOUNDARY) 
	    return CONST_V_PDE_BOUNDARY;
	if (wave_type(*hs) == MOVABLE_BODY_BOUNDARY) 
	    return CONST_V_PDE_BOUNDARY;
	if (wave_type(*hs) == ICE_PARTICLE_BOUNDARY) 
	    return CONST_V_PDE_BOUNDARY;
	if (wave_type(*hs) == DIRICHLET_BOUNDARY) 
	{
	    if (boundary_state(*hs))
	    	return CONST_V_PDE_BOUNDARY;
	    else
	    	return CONST_P_PDE_BOUNDARY;
	}
}	/* ifluid_find_state_at_crossing */

extern  void ifluid_compute_force_and_torque(
        Front *fr,
        HYPER_SURF *hs,
        double dt,
        double *force,
        double *torque)
{
        switch (fr->rect_grid->dim)
        {
        case 2:
            return ifluid_compute_force_and_torque2d(fr,hs,dt,force,torque);
        case 3:
            return ifluid_compute_force_and_torque3d(fr,hs,dt,force,torque);
        }
}       /* end ifluid_compute_force_and_torque */

static  void ifluid_compute_force_and_torque2d(
        Front *front,
        HYPER_SURF *hs,
        double dt,
        double *force,
        double *torque)
{
        RECT_GRID *gr = computational_grid(front->interf);
        double f[MAXD],rr[MAXD];
        double t,pres;
        double area[MAXD],posn[MAXD];
        BOND *b;
        boolean pos_side;
        int i,dim = gr->dim;
        IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
        double *gravity = iFparams->gravity;
        CURVE *curve = Curve_of_hs(hs);

        if (debugging("rigid_body"))
            (void) printf("Entering ifluid_compute_force_and_torque2d()\n");

        if (ifluid_comp(negative_component(curve)))
            pos_side = NO;
        else
            pos_side = YES;

        for (i = 0; i < dim; ++i)
        {
            force[i] = 0.0;
        }
        *torque = 0.0;
        for (b = curve->first; b != NULL; b = b->next)
        {
            if (force_on_hse(Hyper_surf_element(b),Hyper_surf(curve),gr,
                        &pres,area,posn,pos_side))
            {
                for (i = 0; i < dim; ++i)
                {
                    f[i] = pres*area[i];
                    rr[i] = posn[i] - rotation_center(curve)[i];
                    //rr[i] = 0.5*(Coords(b->start)[i] + Coords(b->end)[i])
                                //- rotation_center(curve)[i];
                    force[i] += f[i];
                }
                Cross2d(rr,f,t);
                *torque += t;
            }
        }
         /* Add gravity to the total force */
        if (motion_type(curve) != ROTATION)
        {
            for (i = 0; i < dim; ++i)
                force[i] += gravity[i]*total_mass(curve);
        }
        if (debugging("rigid_body"))
        {
            (void) printf("Leaving ifluid_compute_force_and_torque2d()\n");
            (void) printf("total_force = %f %f\n",force[0],force[1]);
            (void) printf("torque = %f\n",*torque);
        }
}       /* end ifluid_compute_force_and_torque2d */

#define         MAX_TRI_FOR_INTEGRAL            100
static  void ifluid_compute_force_and_torque3d(
        Front *front,
        HYPER_SURF *hs,
        double dt,
        double *force,
        double *torque)
{
        RECT_GRID *gr = computational_grid(front->interf);
        double f[MAXD],rr[MAXD];
        double t[MAXD],tdir,pres;
        double area[MAXD],posn[MAXD];
        TRI *tri;
        boolean pos_side;
        int i,dim = gr->dim;
        IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
        double *gravity = iFparams->gravity;
        SURFACE *surface = Surface_of_hs(hs);

        if (ifluid_comp(negative_component(surface)))
            pos_side = NO;
        else
            pos_side = YES;

        for (i = 0; i < dim; ++i)
        {
            force[i] = 0.0;
            torque[i] = 0.0;
        }
        for (tri = first_tri(surface); !at_end_of_tri_list(tri,surface);
                        tri = tri->next)
        {
            if (force_on_hse(Hyper_surf_element(tri),Hyper_surf(surface),gr,
                        &pres,area,posn,pos_side))
            {
                for (i = 0; i < dim; ++i)
                {
                    f[i] = pres*area[i];
                    force[i] += f[i];
                    rr[i] = posn[i] - rotation_center(surface)[i];
                }
                Cross3d(rr,f,t);
                tdir = Dot3d(t,(rotation_direction(hs)));
                for (i = 0; i < dim; ++i)
                {
                    t[i] = tdir*rotation_direction(hs)[i];
                    torque[i] += t[i];
                }
            }
        }
         /* Add gravity to the total force */
        if (motion_type(surface) != ROTATION)
        {
            for (i = 0; i < dim; ++i)
                force[i] += gravity[i]*total_mass(surface);
        }
        if (debugging("rigid_body"))
        {
            printf("In ifluid_compute_force_and_torque3d()\n");
            printf("total_force = %f %f %f\n",force[0],force[1],force[2]);
            printf("torque = %f %f %f\n",torque[0],torque[1],torque[2]);
        }
}       /* end ifluid_compute_force_and_torque3d */

static boolean force_on_hse(
        HYPER_SURF_ELEMENT *hse,        /* Bond (2D) or tri (3D) */
        HYPER_SURF *hs,                 /* Curve (2D) or surface (3D) */
        RECT_GRID *gr,                  /* Rectangular grid */
        double *pres,           /* Average pressure */
        double *area,           /* Area as a vector, pointing onto body */
        double *posn,           /* Position of the pressure */
        boolean pos_side)       /* Is the body on the positive side of hs? */
{
        int dim = gr->dim;
        switch (dim)
        {
        case 2:
            return force_on_hse2d(hse,hs,gr,pres,area,posn,pos_side);
        case 3:
            return force_on_hse3d(hse,hs,gr,pres,area,posn,pos_side);
        default:
            return NO;
        }

}       /* end force_on_hse */

static boolean force_on_hse2d(
        HYPER_SURF_ELEMENT *hse,
        HYPER_SURF *hs,
        RECT_GRID *gr,
        double *pres,
        double *area,
        double *posn,
        boolean pos_side)
{
        double crds1[MAXD],crds2[MAXD];
        double p1,p2;
        Locstate s1,s2;
        BOND *b = Bond_of_hse(hse);
        CURVE *c = Curve_of_hs(hs);
        double *L = gr->L;
        double *U = gr->U;
        int i;

        /* Get pressure at two end points of the bond */
        if (b->start == c->start->posn)
            s1 = pos_side ? right_start_state(c) : left_start_state(c);
        else
            s1 = pos_side ? right_state(b->start) : left_state(b->start);
        if (b->end == c->end->posn)
            s2 = pos_side ? right_end_state(c) : left_end_state(c);
        else
            s2 = pos_side ? right_state(b->end) : left_state(b->end);

        p1 = getStatePres(s1);  p2 = getStatePres(s2);
        for (i = 0; i < 2; ++i)
        {
            crds1[i] = Coords(b->start)[i];
            crds2[i] = Coords(b->end)[i];
        }

        /* Cut and interpolate if one end is outside the domain */
        for (i = 0; i < 2; ++i)
        {
            if (crds1[i] <= L[i])
            {
                if (crds2[i] <= L[i]) return NO; // both ends out
                else
                {
                    crds1[(i+1)%2] = intrp_between(crds1[i],crds2[i],L[i],
                                crds1[(i+1)%2],crds2[(i+1)%2]);
                    p1 = intrp_between(crds1[i],crds2[i],L[i],p1,p2);
                    crds1[i] = L[i];
                }
            }
            if (crds1[i] >= U[i])
            {
                if (crds2[i] >= U[i]) return NO; // both ends out
                else
                {
                    crds1[(i+1)%2] = intrp_between(crds1[i],crds2[i],U[i],
                                crds1[(i+1)%2],crds2[(i+1)%2]);
                    p1 = intrp_between(crds1[i],crds2[i],U[i],p1,p2);
                    crds1[i] = U[i];
                }
            }
        }
        for (i = 0; i < 2; ++i)
        {
            if (crds2[i] <= L[i])
            {
                if (crds1[i] <= L[i]) return NO; // both ends out
                else
                {
                    crds2[(i+1)%2] = intrp_between(crds1[i],crds2[i],L[i],
                                crds1[(i+1)%2],crds2[(i+1)%2]);
                    p2 = intrp_between(crds1[i],crds2[i],L[i],p1,p2);
                    crds2[i] = L[i];
                }
            }
            if (crds2[i] >= U[i])
            {
                if (crds1[i] >= U[i]) return NO; // both ends out
                else
                {
                    crds2[(i+1)%2] = intrp_between(crds1[i],crds2[i],U[i],
                                crds1[(i+1)%2],crds2[(i+1)%2]);
                    p2 = intrp_between(crds1[i],crds2[i],U[i],p1,p2);
                    crds2[i] = U[i];
                }
            }
        }
        area[0] = pos_side ? crds1[1] - crds2[1] : crds2[1] - crds1[1];
        area[1] = pos_side ? crds2[0] - crds1[0] : crds1[0] - crds2[0];
        *pres = 0.5*(p1 + p2);
        posn[0] = 0.5*(crds1[0] + crds2[0]);
        posn[1] = 0.5*(crds1[1] + crds2[1]);
        return YES;
}       /* end force_on_hse2d */

static boolean force_on_hse3d(
        HYPER_SURF_ELEMENT *hse,
        HYPER_SURF *hs,
        RECT_GRID *gr,
        double *pres,
        double *area,
        double *posn,
        boolean pos_side)
{
        TRI *t = Tri_of_hse(hse);
        POINT *point;
        Locstate sl,sr;
        int i,j,dim = gr->dim;

        *pres = 0.0;
        for (i = 0; i < 3; ++i)
            posn[i] = 0.0;
        for (i = 0; i < 3; ++i)
        {
            point = Point_of_tri(t)[i];
            for (j = 0; j < dim; ++j)
                posn[j] += Coords(point)[j];
            FT_GetStatesAtPoint(point,hse,hs,&sl,&sr);
            if (pos_side)
                *pres += getStatePres(sr);
            else
                *pres += getStatePres(sl);
        }
        *pres /= 3.0;
        for (i = 0; i < dim; ++i)
        {
            area[i] = pos_side ? -Tri_normal(t)[i] : Tri_normal(t)[i];
            posn[i] /= 3.0;
        }
        /* Need to treat subdomain boundary */
        return YES;
}       /* end force_on_hse3d */

static double intrp_between(
        double x1,
        double x2,
        double x,
        double y1,
        double y2)
{
        double y;
        if (x1 == x2) return y1;
        y = y1 + (y2 - y1)/(x2 - x1)*(x - x1);
        return y;
}

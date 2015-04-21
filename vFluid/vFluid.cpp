/****************************************************************
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
*****************************************************************/

#include "vcartsn.h"

	/*  Function Declarations */
static void vfluid_driver(Front*,CARTESIAN*);
static void promptForDirichletBdryState(FILE*,Front*,HYPER_SURF**,int);
static void vF_flowThroughBoundaryState(double*,HYPER_SURF*,Front*,POINTER,
				POINTER);

char *in_name,*restart_state_name,*restart_name;
boolean RestartRun;
boolean ReSetTime;
int RestartStep;

int main(int argc, char **argv)
{
	static Front front;
	static F_BASIC_DATA f_basic;
	static LEVEL_FUNC_PACK level_func_pack;
	static PARAMS params;
	static VELO_FUNC_PACK velo_func_pack;

	FT_Init(argc,argv,&f_basic);
	f_basic.dim = 2;
	f_basic.size_of_intfc_state = sizeof(STATE);

	/* Initialize basic computational data */

        in_name                 = f_basic.in_name;
        restart_state_name      = f_basic.restart_state_name;
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

	/* Initialize solver for scalar field*/
        CARTESIAN *cartesian = new CARTESIAN(front);

	level_func_pack.pos_component = LIQUID_COMP;
	if (!RestartRun)
	{
	    FT_InitIntfc(&front,&level_func_pack);
	    read_vF_dirichlet_bdry_data(&front);
	}
	else
	{
	}

	/* Time control */
	FT_ReadTimeControl(in_name,&front);


	cartesian->initMesh();

        if (!RestartRun)
	{
	    cartesian->initFlowState();
	}
        else
	{
	}


        cartesian->initDrawVariables();

	/* Propagate the front */

	vfluid_driver(&front,cartesian);

	clean_up(0);
}

static  void vfluid_driver(
        Front *front,
	CARTESIAN *cartesian)
{
        double CFL;
        int  dim = front->rect_grid->dim;

        CFL = Time_step_factor(front);
	Tracking_algorithm(front) = STRUCTURE_TRACKING;

	if (!RestartRun || ReSetTime)
	{
	    FT_ResetTime(front);

	    FT_Save(front);
            FT_Draw(front);

	    FT_Propagate(front);

            cartesian->solve(front->dt);

            FT_SetOutputCounter(front);
	    FT_SetTimeStep(front);
	    front->dt = std::min(front->dt,CFL*cartesian->max_dt);
	}
	else
	{
	    FT_SetOutputCounter(front);
	}
	FT_TimeControlFilter(front);
	FT_PrintTimeStamp(front);
	
        for (;;)
        {
	    /* Propagating interface for time step dt */

            cartesian->solve(front->dt);

	    FT_AddTimeStepToCounter(front);

	    //Next time step determined by maximum speed of previous
	    //step, assuming the propagation is hyperbolic and
	    //is not dependent on second order derivatives of
	    //the interface such as curvature, and etc.

	    FT_SetTimeStep(front);
            front->dt = std::min(front->dt,CFL*cartesian->max_dt);

	    /* Output section */

            if (FT_IsSaveTime(front))
	    {
		FT_Save(front);
	    }
            if (FT_IsDrawTime(front))
	    {
                FT_Draw(front);
	    }

            if (FT_TimeLimitReached(front))
	    {
		FT_PrintTimeStamp(front);
                break;
	    }

	    /* Time and step control section */

	    FT_TimeControlFilter(front);
	    FT_PrintTimeStamp(front);
        }
}       /* end vfluid_driver */

extern double getStateXvel(POINTER state)
{
        STATE *fstate = (STATE*)state;
        return fstate->vel[0];
}       /* end getStateXvel */

extern double getStateYvel(POINTER state)
{
        STATE *fstate = (STATE*)state;
        return fstate->vel[1];
}       /* end getStateYvel */

extern double getStateEta(POINTER state)
{
        STATE *fstate = (STATE*)state;
        return fstate->eta;
}       /* end getStateEta */


extern void read_vF_dirichlet_bdry_data(
	Front *front)
{
	char msg[100];
	int i,j,dim = front->rect_grid->dim;
	FILE *infile = fopen(InName(front),"r");
	INTERFACE *intfc = front->interf;
	HYPER_SURF *hs;
	int i_hs = 0;

	(void) printf("Available type of Dirichlet boundary include: \n");
	(void) printf("\tConstant state (C)\n");
	(void) printf("\tFlow through (F)\n");
	for (i = 0; i < dim; ++i)
	for (j = 0; j < 2; ++j)
	{
	    if (rect_boundary_type(intfc,i,j) == DIRICHLET_BOUNDARY)
	    {
		hs = NULL;
	        if (rect_boundary_type(front->interf,i,j) == DIRICHLET_BOUNDARY)
		    hs = FT_RectBoundaryHypSurf(front->interf,
					DIRICHLET_BOUNDARY,i,j);
		if (j == 0)
		    sprintf(msg,"For lower boundary in %d-th dimension",i);
		else
		    sprintf(msg,"For upper boundary in %d-th dimension",i);
		CursorAfterString(infile,msg);
		(void) printf("\n");
		promptForDirichletBdryState(infile,front,&hs,i_hs);
		i_hs++;
	    }
	    else if (rect_boundary_type(intfc,i,j) == MIXED_TYPE_BOUNDARY)
            {
		HYPER_SURF **hss;
		int k,nhs;
                hss = FT_MixedBoundaryHypSurfs(intfc,i,j,DIRICHLET_BOUNDARY,
                                        &nhs);
                printf("Number of Dirichlet boundaries on dir %d side %d: %d\n",
                                        i,j,nhs);
                if (dim == 2)
                {
                    for (k = 0; k < nhs; ++k)
                    {
                        CURVE *c = Curve_of_hs(hss[k]);
                        (void) printf("Curve %d start and end at: ",k+1);
                        (void) printf("(%f %f)->(%f %f)\n",
                                  Coords(c->start->posn)[0],
                                  Coords(c->start->posn)[1],
                                  Coords(c->end->posn)[0],
                                  Coords(c->end->posn)[1]);
                        promptForDirichletBdryState(infile,front,hss+k,i_hs);
                        i_hs++;
                    }
                }
            }
	}
	fclose(infile);
}	/* end read_vF_dirichlet_bdry_data */

static void promptForDirichletBdryState(
	FILE *infile,
	Front *front,
	HYPER_SURF **hs,
        int i_hs)
{
	static STATE *state;
	char s[100];
	POINTER func_params;
	int dim = FT_Dimension();
	int k;

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
	    CursorAfterString(infile,"Enter eta:");
	    fscanf(infile,"%lf",&state->eta);
	    (void) printf("%f\n",state->eta);
	    FT_InsertDirichletBoundary(front,NULL,NULL,
			NULL,(POINTER)state,*hs,i_hs);
	    break;
	case 'f':			// Flow through state
	case 'F':
	    FT_InsertDirichletBoundary(front,
			vF_flowThroughBoundaryState,"flowThroughBoundaryState",
			NULL,NULL,*hs,i_hs);
	    break;
	default:
	    (void) printf("Unknown Dirichlet boundary!\n");
	    clean_up(ERROR);
	}
}	/* end promptForDirichletBdryState */

static void vF_flowThroughBoundaryState(
        double          *p0,
        HYPER_SURF      *hs,
        Front           *front,
        POINTER         ft_params,
        POINTER         state)
{
	Tan_stencil **tsten;
	Nor_stencil *nsten;
	PARAMS *params = (PARAMS*)front->extra1;
	VF_FIELD *field = params->field;
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
	    printf("Entering iF_flowThroughBoundaryState()\n");

	/*
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
	*/
}       /* end vF_flowThroughBoundaryState */

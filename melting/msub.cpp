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
#include "melting.h"

        /*  Function Declarations */
static void neumann_point_propagate(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
static void dirichlet_point_propagate(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
static void contact_point_propagate(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
static void zero_state(COMPONENT,double*,IF_FIELD*,int,int,IF_PARAMS*);

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
        double vel[MAXD],s;
        int i, dim = front->rect_grid->dim;
        double dn,*h = front->rect_grid->h;
        double *m_pre = iFparams->field->pres;
        double *m_vor = iFparams->field->vort;
        double nor[MAXD],tan[MAXD],p1[MAXD];
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
        tan[0] = -nor[1];       tan[1] = nor[0];

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
        }
        else
        {
            for (i = 0; i < dim; ++i)
                vel[i] = 0.0;
        }
        for (i = 0; i < dim; ++i)
        {
            Coords(newp)[i] = Coords(newp)[i] + dt*vel[i];
            newst->vel[i] = vel[i];
            FT_RecordMaxFrontSpeed(i,fabs(vel[i]),NULL,Coords(newp),front);
        }
        FT_IntrpStateVarAtCoords(front,comp,p1,m_pre,
                        getStatePres,&newst->pres,&oldst->pres);
        if (dim == 2)
        {
            FT_IntrpStateVarAtCoords(front,comp,p1,m_vor,
                        getStateVort,&newst->vort,&oldst->vort);
        }
        s = mag_vector(vel,dim);
        FT_RecordMaxFrontSpeed(dim,s,NULL,Coords(newp),front);
        return;
}       /* end neumann_point_propagate */

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
        STATE *newst;
        STATE *bstate;
        FLOW_THROUGH_PARAMS ft_params;
        COMPONENT comp;

        if (debugging("dirichlet_bdry"))
        {
            printf("Entering dirichlet_point_propagate()\n");
            print_general_vector("oldp:  ",Coords(oldp),dim,"\n");
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
        if (newst == NULL) return;      // node point
	
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
            newst->vort = 0.0;

            if (debugging("dirichlet_bdry"))
            {
                printf("Preset boundary state:\n");
                print_general_vector("Velocity: ",newst->vel,dim,"\n");
                printf("Pressure: %f\n",newst->pres);
                printf("Vorticity: %f\n",newst->vort);
            }
        }
        else if (boundary_state_function(oldhs))
        {
            ft_params.oldp = oldp;
            ft_params.comp = comp;
            (*boundary_state_function(oldhs))(Coords(oldp),oldhs,front,
                        (POINTER)&ft_params,(POINTER)newst);
	}
        if (debugging("dirichlet_bdry"))
            printf("Leaving dirichlet_point_propagate()\n");
        return;
}       /* end dirichlet_point_propagate */

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
        double vel[MAXD],s;
        int i, dim = front->rect_grid->dim;
        double *m_pre = iFparams->field->pres;
        double *m_vor = iFparams->field->vort;
        double *p0;
        STATE *oldst,*newst;
        POINTER sl,sr;
        double pres,vort;

        (*front->vfunc)(front->vparams,front,oldp,oldhse,oldhs,vel);
        for (i = 0; i < dim; ++i)
        {
            Coords(newp)[i] = Coords(newp)[i] + dt*vel[i];
            FT_RecordMaxFrontSpeed(i,fabs(vel[i]),NULL,Coords(newp),front);
        }

        FT_GetStatesAtPoint(oldp,oldhse,oldhs,&sl,&sr);
        oldst = (STATE*)sl;
        p0 = Coords(newp);
        FT_IntrpStateVarAtCoords(front,-1,p0,m_pre,getStatePres,&pres,
                                &oldst->pres);
        FT_IntrpStateVarAtCoords(front,-1,p0,m_vor,getStateVort,&vort,
                                &oldst->vort);

        newst = (STATE*)left_state(newp);
        newst->vort = vort;
        newst->pres = pres;
        for (i = 0; i < dim; ++i)
            newst->vel[i] = vel[i];
        newst = (STATE*)right_state(newp);
        newst->vort = vort;
        newst->pres = pres;
        for (i = 0; i < dim; ++i)
            newst->vel[i] = vel[i];

        s = mag_vector(vel,dim);
        FT_RecordMaxFrontSpeed(dim,s,NULL,Coords(newp),front);
}       /* end contact_point_propagate */

void read_melt_dirichlet_bdry_data(
        char *inname,
        Front *front,
        F_BASIC_DATA f_basic)
{
        char msg[100],s[100];
        int i,k,dim = front->rect_grid->dim;
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
                case 'c':                       // Constant state
                case 'C':
                    CursorAfterString(infile,"Enter velocity:");
                    for (k = 0; k < dim; ++k)
                    {
                        fscanf(infile,"%lf",&state.vel[k]);
                        (void) printf("%f ",state.vel[k]);
                    }
		    (void) printf("\n");
                    CursorAfterString(infile,"Enter pressure:");
                    fscanf(infile,"%lf",&state.pres);
                    (void) printf("%f\n",state.pres);
                    CursorAfterString(infile,"Enter fluid temperature:");
                    fscanf(infile,"%lf",&state.temperature);
                    (void) printf("%f\n",state.temperature);
                    FT_InsertDirichletBoundary(front,NULL,NULL,NULL,
                                (POINTER)&state,hs,i_surf);
                    break;
		case 'f':                       // Flow through state
                case 'F':
                    FT_InsertDirichletBoundary(front,melt_flowThroughBoundaryState,
                                "flowThroughBoundaryState",NULL,NULL,hs,i_surf);
                    break;
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
                case 'c':                       // Constant state
                case 'C':
                    CursorAfterString(infile,"Enter velocity:");
                    for (k = 0; k < dim; ++k)
                    {
                        fscanf(infile,"%lf ",&state.vel[k]);
                        (void) printf("%f ",state.vel[k]);
                    }
                    (void) printf("\n");
                    CursorAfterString(infile,"Enter pressure:");
                    fscanf(infile,"%lf",&state.pres);
                    (void) printf("%f\n",state.pres);
                    CursorAfterString(infile,"Enter fluid temperature:");
                    fscanf(infile,"%lf",&state.temperature);
                    (void) printf("%f\n",state.temperature);
                    FT_InsertDirichletBoundary(front,NULL,NULL,NULL,
                                (POINTER)&state,hs,i_surf);
                    break;
                case 'f':                       // Flow through state
                case 'F':
                    FT_InsertDirichletBoundary(front,melt_flowThroughBoundaryState,
                                "flowThroughBoundaryState",NULL,NULL,hs,i_surf);
                    break;
                }
            }
	}
        fclose(infile);
}       /* end read_melt_dirichlet_bdry_data */

static double (*getStateVel[MAXD])(POINTER) = {getStateXvel,getStateYvel,
                                        getStateZvel};

extern void melt_flowThroughBoundaryState(
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
        PARAMS *eqn_params = (PARAMS*)front->extra2;
        IF_FIELD *iF_field = iFparams->field;
        PHASE_FIELD *pH_field = eqn_params->field;
        double dir[MAXD];
        double u[3];            /* velocity in the sweeping direction */
        double v[3][MAXD];      /* velocity in the orthogonal direction */
        double vort[3];         /* vorticity stencil */
        double pres[3];         /* pressure stencil */
        double temp[3];         /* temperature stencil */
        double f_u;             /* u flux in the sweeping direction */
        double f_v[MAXD];       /* v flux in the orthogonal direction */
        double f_vort;          /* vort flux */
        double f_pres;          /* pressure flux */
        double f_temp;          /* temperature flux */
        double dn,dt = front->dt;
        STATE *newst = (STATE*)state;
        STATE  **sts;
        int i,j,dim = front->rect_grid->dim;
        int nrad = 2;

        if (debugging("flow_through"))
            printf("Entering melt_flowThroughBoundaryState()\n");
	
	tsten = FrontGetTanStencils(front,oldp,nrad);
        for (i = 0; i < dim; ++i)
            dir[i] = tsten[0]->dir[i];
        dn = FT_GridSizeInDir(dir,front);

        if (comp == negative_component(hs))
            sts = (STATE**)tsten[0]->leftst;
        else
            sts = (STATE**)tsten[0]->rightst;

        if (debugging("flow_through"))
        {
            printf("Ambient component: %d\n",comp);
            printf("hs = %p  oldp->hs = %p\n",(POINTER)hs,(POINTER)oldp->hs);
            printf("Time step = %f  Tangential grid size = %f\n",dt,dn);
            printf("Tangential direction: ");
            for (j = 0; j < dim; ++j)
                printf("%f ",tsten[0]->dir[j]);
            printf("\n");
            printf("Tan_stencil at point p(%f %f)\n",Coords(oldp)[0],
                                Coords(oldp)[1]);
            printf("Left points:\n");
            for (i = 0; i < nrad; ++i)
            {
                for (j = 0; j < dim; ++j)
                    printf("%f ",Coords(tsten[0]->p[-i])[j]);
                printf("\n");
            }
            printf("Right points:\n");
            for (i = 0; i < nrad; ++i)
            {
                for (j = 0; j < dim; ++j)
                    printf("%f ",Coords(tsten[0]->p[i])[j]);
                printf("\n");
            }
        }

	for (j = 0; j < 3; ++j)
            u[j] = 0.0;
        for (j = 0; j < 3; ++j)
        {
            vort[j] = sts[j-1]->vort;
            pres[j] = sts[j-1]->pres;
            temp[j] = sts[j-1]->temperature;
            for (i = 0; i < dim; ++i)
            {
                u[j] += sts[j-1]->vel[i]*dir[i];
                v[j][i] = sts[j-1]->vel[i]*(1.0 - dir[i]);
            }
        }

        f_u = burger_flux(u[0],u[1],u[2]);
        for (i = 0; i < dim; ++i)
            f_v[i] = linear_flux(u[1],v[0][i],v[1][i],v[2][i]);
        f_vort = linear_flux(u[1],vort[0],vort[1],vort[2]);
        f_pres = linear_flux(u[1],pres[0],pres[1],pres[2]);
        f_temp = linear_flux(u[1],temp[0],temp[1],temp[2]);

        for (i = 0; i < dim; ++i)
            newst->vel[i] = sts[0]->vel[i] - dt/dn*(
                f_u*dir[i] + f_v[i]) ;
        newst->vort = sts[0]->vort - dt/dn*f_vort;
        newst->pres = sts[0]->pres - dt/dn*f_pres;
        newst->temperature = sts[0]->temperature - dt/dn*f_temp;

        nsten = FT_CreateNormalStencil(front,oldp,comp,nrad);
        for (i = 0; i < dim; ++i)
            dir[i] = nsten->nor[i];
        dn = FT_GridSizeInDir(dir,front);

        if (debugging("flow_through"))
        {
            printf("Time step = %f  Normal grid size = %f\n",dt,dn);
            printf("Normal direction: ");
            for (j = 0; j < dim; ++j)
                printf("%f ",nsten->nor[j]);
            printf("\n");
            printf("Nor_stencil at point p(%f %f)\n",Coords(oldp)[0],
                                Coords(oldp)[1]);
            printf("Nor_stencil:\n");
            for (i = 0; i < nrad; ++i)
            {
                for (j = 0; j < dim; ++j)
                    printf("%f ",nsten->pts[i][j]);
                printf("\n");
            }
        }
	
	for (j = 0; j < 3; ++j)
            u[j] = 0.0;
        for (j = 0; j < 2; ++j)
        {
            vort[j] = sts[0]->vort;
            pres[j] = sts[0]->pres;
            temp[j] = sts[0]->temperature;
            for (i = 0; i < dim; ++i)
            {
                u[j] += sts[0]->vel[i]*dir[i];
                v[j][i] = sts[0]->vel[i]*(1.0 - dir[i]);
            }
        }
        for (i = 0; i < dim; ++i)
        {
            double vtmp;
            FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],
                        iF_field->vel[i],getStateVel[i],&vtmp,&sts[0]->vel[i]);
            u[2] += vtmp*dir[i];
            v[2][i] = vtmp*(1.0 - dir[i]);
        }
        if (dim == 2)
        {
            FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],
                        iF_field->vort,getStateVort,&vort[2],&sts[1]->vort);
        }
        FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],iF_field->pres,
                            getStatePres,&pres[2],&sts[1]->pres);
        FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],pH_field->temperature,
                            getStateTemperature,&temp[2],&sts[1]->temperature);

        f_u = burger_flux(u[0],u[1],u[2]);
        for (i = 0; i < dim; ++i)
            f_v[i] = linear_flux(u[1],v[0][i],v[1][i],v[2][i]);
        f_vort = linear_flux(u[1],vort[0],vort[1],vort[2]);
        f_pres = linear_flux(u[1],pres[0],pres[1],pres[2]);
        f_temp = linear_flux(u[1],temp[0],temp[1],temp[2]);

        for (i = 0; i < dim; ++i)
            newst->vel[i] += - dt/dn*(f_u*dir[i] + f_v[i]) ;
        newst->vort += - dt/dn*f_vort;
        newst->pres += - dt/dn*f_pres;
        newst->temperature += - dt/dn*f_temp;
	
	if (debugging("flow_through"))
        {
            printf("flow through boundary state:\n");
            print_general_vector("Velocity: ",newst->vel,dim,"\n");
            printf("Pressure: %f\n",newst->pres);
            printf("Vorticity: %f\n",newst->vort);
            printf("Temperature: %f\n",newst->temperature);
        }
}       /*end melt_flowThroughBoundaryState */

void init_fluid_state_func(
	Front *front,
        Incompress_Solver_Smooth_Basis *l_cartesian)
{
	l_cartesian->getInitialState = zero_state;
}	/* end init_fluid_state_func */

static void zero_state(
        COMPONENT comp,
        double *coords,
        IF_FIELD *field,
        int index,
        int dim,
        IF_PARAMS *iFparams)
{
        int i;
        for (i = 0; i < dim; ++i)
            field->vel[i][index] = 0.0;
        field->pres[index] = 0.0;
}       /* end zero_state */

extern void read_fluid_params(
	Front *front)
{
	char string[100];
	PARAMS *eqn_params = (PARAMS*)front->extra2;
	char *inname = InName(front);
	FILE *infile = fopen(inname,"r");
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;

	if (eqn_params->no_fluid == YES) return;
	/* defaults numerical schemes */
        iFparams->num_scheme.projc_method = SIMPLE;
        iFparams->num_scheme.advec_method = WENO;
        iFparams->num_scheme.ellip_method = SIMPLE_ELLIP;

	CursorAfterString(infile,"Enter density and viscosity of the fluid:");
	fscanf(infile,"%lf %lf",&iFparams->rho2,&iFparams->mu2);
	(void) printf("%f %f\n",iFparams->rho2,iFparams->mu2);
	iFparams->m_comp1 = SOLID_COMP;
	iFparams->m_comp2 = LIQUID_COMP2;
        CursorAfterString(infile,"Enter surface tension:");
        fscanf(infile,"%lf",&iFparams->surf_tension);
	(void) printf("%f\n",iFparams->surf_tension);
        CursorAfterString(infile,"Enter factor of smoothing radius:");
        fscanf(infile,"%lf",&iFparams->smoothing_radius);
	(void) printf("%f\n",iFparams->smoothing_radius);
	iFparams->num_scheme.projc_method = ERROR_PROJC_SCHEME;
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
	(void) printf("The default advection order is WENO-Runge-Kutta 4");
        iFparams->adv_order = 4;
        if (CursorAfterStringOpt(infile,"Enter advection order:"))
        {
            fscanf(infile,"%d",&iFparams->adv_order);
            (void) printf("%d\n",iFparams->adv_order);
        }
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
            }
        }

	fclose(infile);
}	/* end read_fluid_params */

extern double getStatePres(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->pres;
}	/* end getStatePres */

extern double getStatePhi(POINTER state)
{
        STATE *fstate = (STATE*)state;
        return fstate->phi;
}       /* end getStatePres */

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

extern double getStateMu(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->mu;
}	/* end getStateMu */

extern double getStateTemp(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->temperature;
}	/* end getStateTemp */

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
            if (dim == 2)
            {
                fprintf(outfile,"%24.18g %24.18g\n",getStateXvel(sl),
                                getStateXvel(sr));
                fprintf(outfile,"%24.18g %24.18g\n",getStateYvel(sl),
                                getStateYvel(sr));
                fprintf(outfile,"%24.18g %24.18g\n",getStateVort(sl),
                                getStateVort(sr));
            }
            if (dim == 3)
            {
                fprintf(outfile,"%24.18g %24.18g\n",getStateXvel(sl),
                                getStateXvel(sr));
                fprintf(outfile,"%24.18g %24.18g\n",getStateYvel(sl),
                                getStateYvel(sr));
                fprintf(outfile,"%24.18g %24.18g\n",getStateZvel(sl),
                                getStateZvel(sr));
            }
        }
}       /* end fluid_print_front_states */

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
            fscanf(infile,"%lf %lf",&lstate->vel[0],&rstate->vel[0]);
            fscanf(infile,"%lf %lf",&lstate->vel[1],&rstate->vel[1]);
            if (dim == 2)
                fscanf(infile,"%lf %lf",&lstate->vort,&rstate->vort);
            if (dim == 3)
                fscanf(infile,"%lf %lf",&lstate->vel[2],&rstate->vel[2]);
        }
}       /* end fluid_read_front_states */

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
	case NEUMANN_BOUNDARY:
	case MOVABLE_BODY_BOUNDARY:
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
}       /* end isDirichletPresetBdry */

extern double p_jump(
        POINTER params,
        int D,
        double *coords)
{
        return 0.0;
}       /* end p_jump */

extern double grad_p_jump_n(
        POINTER params,
        int D,
        double *N,
        double *coords)
{
        return 0.0;
}       /* end grad_p_jump_n */

extern double grad_p_jump_t(
        POINTER params,
        int D,
        int i,
        double *N,
        double *coords)
{
        return 0.0;
}       /* end grad_p_jump_t */

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
            return NEUMANN_PDE_BOUNDARY;
	if (wave_type(*hs) == GROWING_BODY_BOUNDARY)
            return NEUMANN_PDE_BOUNDARY;
        if (wave_type(*hs) == DIRICHLET_BOUNDARY)
            return DIRICHLET_PDE_BOUNDARY;
}     /*ifluid_find_state_at_crossing */


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


#include <iFluid.h>
#include "climate.h"
#include <time.h>
        /*  Function Declarations */
static void neumann_point_propagate(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
static void dirichlet_point_propagate(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
static void contact_point_propagate(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
static void zero_state(COMPONENT,double*,IF_FIELD*,int,int,IF_PARAMS*);
static void rand_state(COMPONENT,double*,IF_FIELD*,int,int,IF_PARAMS*);
static void Taylor_Green_state(COMPONENT,double*,IF_FIELD*,int,int,IF_PARAMS*);
static void Fourier_state(COMPONENT,double*,IF_FIELD*,int,int,IF_PARAMS*);
static  void compute_ice_particle_force2d_velo(Front*,HYPER_SURF*,double,double*,double*);
static  void compute_ice_particle_force3d_velo(Front*,HYPER_SURF*,double,double*,double*);
static  void compute_ice_particle_force2d_pres(Front*,HYPER_SURF*,double,double*,double*);
static  void compute_ice_particle_force3d_pres(Front*,HYPER_SURF*,double,double*,double*);
static boolean force_on_hse(HYPER_SURF_ELEMENT*,HYPER_SURF*,RECT_GRID*,double*,
                                        double*,double*,boolean);
static boolean force_on_hse2d(HYPER_SURF_ELEMENT*,HYPER_SURF*,RECT_GRID*,
                                        double*,double*,double*,boolean);
static boolean force_on_hse3d(HYPER_SURF_ELEMENT*,HYPER_SURF*,RECT_GRID*,
                                        double*,double*,double*,boolean);
static double intrp_between(double,double,double,double,double);
static double compute_supersat_3d(Front*,SURFACE*);
static double compute_supersat_2d(Front*,CURVE*);
static boolean supersat_on_hse(HYPER_SURF_ELEMENT*,HYPER_SURF*,
			RECT_GRID*,double*,boolean);
static boolean supersat_on_hse2d(HYPER_SURF_ELEMENT*,HYPER_SURF*,
			RECT_GRID*,double*,boolean);
static boolean supersat_on_hse3d(HYPER_SURF_ELEMENT*,HYPER_SURF*,
			RECT_GRID*,double*,boolean);

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
        PARAMS *eqn_params = (PARAMS*)front->extra2;
        switch(eqn_params->init_state)
        {
            case RAND_STATE:
                l_cartesian->getInitialState = rand_state;
                break;
            case TAYLOR_STATE:
                l_cartesian->getInitialState = Taylor_Green_state;
                break;
	    case ZERO_STATE:
		l_cartesian->getInitialState = zero_state;
	        break;
            default:
                l_cartesian->getInitialState = zero_state;
        }
}	/* end init_fluid_state_func */

static void Taylor_Green_state(
        COMPONENT comp,
        double *coords,
        IF_FIELD *field,
        int index,
        int dim,
        IF_PARAMS *iFparams)
{
	int i;
	double tcoords[MAXD];
	double **vel = field->vel;
	for (i = 0; i < dim; i++)
	    tcoords[i] = 2*PI*coords[i] - PI;
	switch (dim)
	{
	    case 2:
	        vel[0][index] = sin(tcoords[0]) * cos(tcoords[1]);
	        vel[1][index] = -cos(tcoords[0]) * sin(tcoords[1]);
		break;
	    case 3:
		vel[0][index] = sin(tcoords[0])*cos(tcoords[1])*cos(tcoords[2]);
		vel[1][index] = -cos(tcoords[0])*sin(tcoords[1])*cos(tcoords[2]);
		vel[2][index] = 0;
		break;
	    default:
		printf("Unknown dim = %d\n",dim);
		clean_up(ERROR);
	}
}       /* end Taylor_Green_state */

static void rand_state(
        COMPONENT comp,
        double *coords,
        IF_FIELD *field,
        int index,
        int dim,
        IF_PARAMS *iFparams)
{
	short unsigned int seed[3] = {time(NULL)-index,
				      time(NULL),
				      time(NULL)+index};
	//short unsigned int seed[3] = {index,index-10,index+10};
	GAUSS_PARAMS gauss_params;
	double r_bar = 0;
	double sigma = 0.5;
        int i;
	gauss_params.mu = r_bar;
	gauss_params.sigma = sigma;
        for (i = 0; i < dim; ++i)
            field->vel[i][index] = gauss_center_limit((POINTER)&gauss_params,seed);
}       /* end rand_state */

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

extern double jumpEpsGradDotNorm(
        POINTER params,
        int D,
        double *N,
        double *P)
{
        printf("Entering jumpEpsGradDotNorm(), to be written\n");
        clean_up(0);
}       /* end jumpEpsGradDotNorm */

extern double jumpT(
        POINTER params,
        int D,
        double *P)
{
        return 0.0;
}       /* end jumpT */

extern double jumpGradDotTan(
        POINTER params,
        int D,
        int i,
        double *N,
        double *P)
{
        return 0.0;
}       /* end jumpGradDotTan */

extern  void compute_ice_particle_force(
        Front *fr,
        HYPER_SURF *hs,
        double dt,
        double *force,
        double *torque)
{
        switch (fr->rect_grid->dim)
        {
        case 2:
            return compute_ice_particle_force2d_velo(fr,hs,dt,force,torque);
        case 3:
            return compute_ice_particle_force3d_velo(fr,hs,dt,force,torque);
        }
}       /* end compute_ice_particle_force */

static  void compute_ice_particle_force2d_pres(
        Front *fr,
        HYPER_SURF *hs,
        double dt,
        double *force,
        double *torque)
{
        RECT_GRID *gr = computational_grid(fr->interf);
        double f[MAXD],rr[MAXD];
        double t,pres;
        double area[MAXD],posn[MAXD];
        BOND *b;
        boolean pos_side;
        int i,dim = gr->dim;
        IF_PARAMS *iFparams = (IF_PARAMS*)fr->extra1;
        double *gravity = iFparams->gravity;
        CURVE *curve = Curve_of_hs(hs);

        if (debugging("rigid_body"))
            (void) printf("Entering compute_ice_particle_force2d()\n");

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
                    rr[i] = 0.5*(Coords(b->start)[i] + Coords(b->end)[i])
                                - rotation_center(curve)[i];
                    force[i] += f[i];
                }
		printf("posn:(%f,%f), area=(%f,%f), pres=%f\n",posn[0],posn[1],area[0],area[1],pres);
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
            (void) printf("Leaving compute_ice_particle_force2d()\n");
            (void) printf("total_force = %f %f\n",force[0],force[1]);
            (void) printf("torque = %f\n",*torque);
        }
}

#define         MAX_TRI_FOR_INTEGRAL            100

static  void compute_ice_particle_force3d_pres(
        Front *fr,
        HYPER_SURF *hs,
        double dt,
        double *force,
        double *torque)
{
        RECT_GRID *gr = computational_grid(fr->interf);
        double *L = gr->L;
        double *U = gr->U;
        double f[MAXD];
        double pres;
        double area[MAXD],posn[MAXD];
        double tri_center[MAXD];
        TRI *tri;
        POINT *p;
        boolean pos_side;
        int i,j,dim = gr->dim;
        IF_PARAMS *iFparams = (IF_PARAMS*)fr->extra1;
        double *gravity = iFparams->gravity;
        SURFACE *surface = Surface_of_hs(hs);
        boolean out_domain_tri;

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
            out_domain_tri = NO;
            for (i = 0; i < dim; ++i)
                tri_center[i] = 0;
            for (j = 0; j < 3; ++j)
            {
                p = Point_of_tri(tri)[j];
                for (i = 0; i < dim; ++i)
                    tri_center[i] += Coords(p)[i];
            }
            for (i = 0; i < dim; ++i)
            {
                tri_center[i] /= 3.0;
                if (tri_center[i] <= L[i] || tri_center[i] > U[i])
                    out_domain_tri = YES;
            }
            if (out_domain_tri == YES)
            {
                continue;
            }

            if (force_on_hse(Hyper_surf_element(tri),Hyper_surf(surface),gr,
                        &pres,area,posn,pos_side))
            {
                for (i = 0; i < dim; ++i)
                {
                    f[i] = pres*area[i];
                    force[i] += f[i];
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
            printf("In compute_ice_particle_force3d()\n");
            printf("body_index = %d\n",body_index(hs));
            printf("total_force = %f %f %f\n\n",force[0],force[1],force[2]);
        }
}

static  void compute_ice_particle_force2d_velo(
        Front *fr,
        HYPER_SURF *hs,
        double dt,
        double *force,
        double *torque)
{
        RECT_GRID *gr = FT_GridIntfcTopGrid(fr);
        int *gmax = FT_GridIntfcTopGmax(fr);
        double cvel,*coords;
        int i,index,dim = gr->dim;
	int ic[MAXD];
        IF_PARAMS *iFparams = (IF_PARAMS*)fr->extra1;
        double **vel = iFparams->field->vel;
        double *gravity = iFparams->gravity;
        PARAMS *eqn_params = (PARAMS*)fr->extra2;
        CURVE *curve = Curve_of_hs(hs);
	BOND  *b;
	POINT *p;

        /*For computing finite respone time*/
        double R        = spherical_radius(curve);/*droplet radius*/
        double rho_l    = eqn_params->rho_l;/*water droplet density*/
        double rho_0    = iFparams->rho2;/*fluid density*/
        double mu       = iFparams->mu2;/*viscosity*/
        double tau_p    = 2 * rho_l*R*R/(9*rho_0*mu);/*response time*/

        if (debugging("rigid_body"))
            (void) printf("Entering compute_ice_particle_force2d()\n");

        for (i = 0; i < dim; ++i)
        {
            force[i] = 0.0;
        }
        *torque = 0.0;

	int max_index;
	double max_vel = 0;
	double vel_mag;
	for (b = curve->first; b != NULL; b=b->next)
	{
	    p = b->end;
	    rect_in_which(Coords(p),ic,gr);
	    index = d_index(ic,gmax,dim);
	    vel_mag = sqr(vel[0][index])+sqr(vel[1][index]);
	    if(vel_mag > max_vel && vel_mag < HUGE)
	    {
		max_vel = vel_mag;
		max_index = index;
	    }
	}
	
        for(i = 0; i < dim; ++i)
        {
            cvel   = center_of_mass_velo(curve)[i];
            force[i]   = (vel[i][max_index]-cvel)/tau_p
			  + gravity[i];
            force[i]  *= total_mass(curve);
        }
        if (debugging("rigid_body"))
        {
            (void) printf("Leaving compute_ice_particle_force2d()\n");
            (void) printf("body_index = %d\n",body_index(hs));
            (void) printf("Flow velo:[%f,%f]\n",vel[0][max_index],vel[1][max_index]);
            (void) printf("cent velo:[%f,%f]\n",center_of_mass_velo(curve)[0],
                                                center_of_mass_velo(curve)[1]);
        }
}       /* end compute_ice_particle_force2d */

static  void compute_ice_particle_force3d_velo(
        Front *fr,
        HYPER_SURF *hs,
        double dt,
        double *force,
        double *torque)
{
        RECT_GRID *gr = FT_GridIntfcTopGrid(fr);
        int *gmax = FT_GridIntfcTopGmax(fr);
        double cvel,center[MAXD];
	int ic[MAXD];
        int i,index,dim = gr->dim;
        IF_PARAMS *iFparams = (IF_PARAMS*)fr->extra1;
        double **vel = iFparams->field->vel;
        PARAMS *eqn_params = (PARAMS*)fr->extra2;
        double *gravity = iFparams->gravity;
        SURFACE *surface = Surface_of_hs(hs);
	TRI *tri;
	POINT *p;

        /*For computing finite respone time*/
        double R = spherical_radius(surface);/*droplet radius*/
        double rho_l    = eqn_params->rho_l;/*water droplet density*/
        double rho_0      = iFparams->rho2;/*fluid density*/
        double mu       = iFparams->mu2;/*viscosity*/
        double tau_p    = 2 * rho_l*R*R/(9*rho_0*mu);/*response time*/

        if (debugging("rigid_body"))
            (void) printf("Entering compute_ice_particle_force2d()\n");

        for (i = 0; i < dim; ++i)
        {
            force[i] = 0.0;
        }
        *torque = 0.0;

        int max_index;
        double max_vel = 0;
	double vel_mag;
	for(tri = first_tri(surface); !at_end_of_tri_list(tri,surface);
                        tri = tri->next)
        {
	    for(i = 0; i < 3; ++i)
	    {
	        p = Point_of_tri(tri)[i];
                rect_in_which(Coords(p),ic,gr);
                index = d_index(ic,gmax,dim);
		vel_mag = sqr(vel[0][index])+sqr(vel[1][index]) 
                                     +sqr(vel[2][index]);
                if(vel_mag > max_vel && vel_mag < HUGE)
                {
                    max_vel = vel_mag;
                    max_index = index;
                }
	    }
        }

        /*get information for center of mass*/
        for(i = 0; i < dim; ++i)
        {
            force[i]   = (vel[i][max_index]-cvel)/tau_p
                         + gravity[i];
            force[i]  *= total_mass(surface);
        }
        if (debugging("rigid_body"))
        {
            printf("In compute_ice_particle_force3d()\n");
            (void) printf("body_index = %d\n",body_index(hs));
            (void) printf("Flow velo:[%f,%f,%f] at point(%f,%f,%f)\n",
			vel[0][max_index],vel[1][max_index],vel[2][max_index],
			Coords(p)[0],Coords(p)[1],Coords(p)[2]);
            (void) printf("cent velo:[%f,%f,%f]\n",center_of_mass_velo(surface)[0],
                                                   center_of_mass_velo(surface)[1],
                                                   center_of_mass_velo(surface)[2]);
        }
}       /* end compute_ice_particle_force3d */

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
	double v1[MAXD],V2[MAXD];
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
                if (crds2[i] <= L[i]) return NO; /*both ends out*/
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
                if (crds2[i] >= U[i]) return NO; /* both ends out*/
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


        area[0] = fabs(crds1[1] - crds2[1]);
        area[1] = fabs(crds1[0] - crds2[0]);
	if(crds2[1] > crds1[1])
	    area[0] *= -1;
	if(crds2[0] < crds1[0])
	    area[1] *= -1;
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

extern void ParticlePropagate(Front *fr)
{
        RECT_GRID *gr = FT_GridIntfcTopGrid(fr);
        RECT_GRID *rect_grid = fr->rect_grid;
        IF_PARAMS *iFparams = (IF_PARAMS*)fr->extra1;
	PARAMS *eqn_params = (PARAMS*)fr->extra2;
	PARTICLE* particle_array = eqn_params->particle_array;
	double **vel = iFparams->field->vel;
	double *supersat = eqn_params->field->supersat;
	double *gravity = iFparams->gravity;
        int *gmax = FT_GridIntfcTopGmax(fr);
	int i, j, index, dim = gr->dim;
	double T;
	int ic[MAXD];
	double *center;
	double s; /*restore local supersaturation*/
	double *cvel; /*center velocity for droplets*/
	double a;  /*acceleration*/
	double dt = fr->dt;

        /*computing finite respone time*/
        double rho_0    = iFparams->rho2;/*fluid density*/
        double mu       = iFparams->mu2;/*viscosity*/
	double R, rho, tau_p, delta_R;
	double R_max = 0;
	double R_min = HUGE;

	for (i = 0; i < eqn_params->num_drops; i++)
	{
            /*computing finite respone time*/
            R        = particle_array[i].radius;/*droplet radius*/
            rho      = particle_array[i].rho;/*water droplet density*/
            tau_p    = 2 * rho*R*R/(9*rho_0*mu);/*response time*/

	    if (R == 0)
	    {
		R_min = 0;
	        continue;
	    }
	    /*find index at coords*/
	    center = particle_array[i].center;
	    rect_in_which(center,ic,gr);
	    index = d_index(ic,gmax,dim);
	    cvel = particle_array[i].vel;
	    /*compute radius for particle[i]*/
	    s = supersat[index];
	    delta_R = R*R+2*eqn_params->K*s*dt;
	    if(delta_R < 0)
		R = 0;
	    else
	        R = sqrt(delta_R);

	    particle_array[i].radius = R;
	    /*save max and min radius*/
	    if(R > R_max)
		R_max = R;
	    if(R < R_min)
		R_min = R;

	    /*compute velocity for particle[i] with implicit method*/
	    for(j = 0; j < dim; ++j)
            {
		/*compute velocity*/
		cvel[j] += vel[j][index]/tau_p + gravity[j];
		cvel[j] /= (1+1/tau_p); 
	  	/*compute center of particle[i]*/
		center[j] += cvel[j]*dt;

		/*handle periodic drops*/
		T = rect_grid->U[j]-rect_grid->L[j];	
		if (center[j] > rect_grid->U[j])
		    center[j] = rect_grid->L[j]+fmod(center[j],T);
		if (center[j] < rect_grid->L[j])
		    center[j] = rect_grid->U[j]+fmod(center[j],T);
		if(isnan(center[j]))
		{
		    printf("center[%d]=nan, T = %f, domain=[%f,%f]\n",
				j,T,rect_grid->L[j],rect_grid->U[j]);
		    clean_up(ERROR);
		}
	    }

	    if (debugging("particles"))
	    {
	        printf("\nDrop[%d]:\n",i);
	        printf("Supersat: %f\n",s);
	        printf("Condensation rate: %20.19f\n",eqn_params->K);
	        printf("delta_R = %f\n",delta_R);
	        printf("dt = %f\n",dt);
	        printf("Radius:%15.14f\n",R);
	        printf("center:[%f,%f,%f]\n", center[0],center[1],center[2]);
	        printf("Flow_vel[%f,%f,%f]\nc_vel[%f,%f,%f]\n",
		    vel[0][index],vel[1][index],vel[2][index],
		    cvel[0],cvel[1],cvel[2]);
	        printf("Response time = %f\n",tau_p);
	    }
	}
	printf("max radius = %20.14f, min radius = %20.14f\n",R_max,R_min);
}

extern void CondensationPreAdvance(Front *fr)
{
        INTERFACE *intfc = fr->interf;
	PARAMS *eqn_params = (PARAMS*)fr->extra2;
        CURVE **c;
        SURFACE **s;
	double cond_speed;
	double avg_supersat;
        switch (fr->rect_grid->dim)
        {
        case 2:
            for (c = intfc->curves; c && *c; ++c)
            {
                if(wave_type(*c) == ICE_PARTICLE_BOUNDARY)
                {
		    compute_supersat_2d(fr,*c);
		    /*update radius*/
                    avg_supersat = f_hyper_surf(*c)->_surface_tension;
                    cond_speed = eqn_params->K*avg_supersat
                                        /spherical_radius(*c);
                    spherical_radius(*c) += cond_speed * (fr->dt);
		    /*update mass*/
		    total_mass(*c) = eqn_params->rho_l*PI*4.0/3.0*
				spherical_radius(*c)*spherical_radius(*c)*
				spherical_radius(*c);
		}
            }
	    break;
        case 3:
            for (s = intfc->surfaces; s && *s; ++s)
            {
                if(wave_type(*s) == ICE_PARTICLE_BOUNDARY)
                {
		    compute_supersat_3d(fr,*s);
                    /*update radius*/
		    avg_supersat = f_hyper_surf(*s)->_surface_tension;
		    cond_speed = eqn_params->K*avg_supersat
					/spherical_radius(*s);
		    spherical_radius(*s) += cond_speed * (fr->dt);
		    /*update total mass*/
                    total_mass(*s) = eqn_params->rho_l*PI*4.0/3.0*
                                spherical_radius(*s)*spherical_radius(*s)*
                                spherical_radius(*s);
		    printf("Radius=%19.18f, average supersat=%f, condensation speed=%f\n",
			    spherical_radius(*s),avg_supersat,cond_speed);
		}
            }
	    break;
	default:
	    printf("Unkown dim = %i\n",fr->rect_grid->dim);
	    clean_up(ERROR);
        }
}

static double compute_supersat_3d(Front *fr,SURFACE *surf)
{
        RECT_GRID *gr = computational_grid(fr->interf);
        double avg_supersat,supersat;
        TRI *tri;
        boolean pos_side;
        int count,i,dim = gr->dim;
        PARAMS *eqn_params = (PARAMS*)fr->extra2;

        if(negative_component(surf) == LIQUID_COMP2)
            pos_side = NO;
        else
            pos_side = YES;
        avg_supersat = 0;
        count =0;
        for (tri = first_tri(surf); !at_end_of_tri_list(tri,surf);
			tri = tri->next)
        {
            if(supersat_on_hse(Hyper_surf_element(tri),Hyper_surf(surf),gr,
                        &supersat,pos_side))
            {
                avg_supersat += supersat;
                count ++;
            }
        }
        if(count != 0)
            avg_supersat = avg_supersat/count;
        else
            printf("WARNING: no supersat found\n");
        /*use _surface_tension as a temporary space*/
        f_hyper_surf(surf)->_surface_tension = avg_supersat;
}

static double compute_supersat_2d(Front *fr,CURVE *curve)
{
        RECT_GRID *gr = computational_grid(fr->interf);
        double avg_supersat,supersat;
        BOND *b;
        boolean pos_side;
        int count,i,dim = gr->dim;
        PARAMS *eqn_params = (PARAMS*)fr->extra2;

        if(negative_component(curve) == LIQUID_COMP2)
            pos_side = NO;
        else
            pos_side = YES;
        avg_supersat = 0;
        count =0;
        for (b = curve->first; b!=NULL; b=b->next)
        {
            if(supersat_on_hse(Hyper_surf_element(b),Hyper_surf(curve),gr,
                        &supersat,pos_side))
            {
                avg_supersat += supersat;
                count ++;
            }
        }
        if(count != 0)
            avg_supersat = avg_supersat/count;
        else
            printf("WARNING: no supersat found\n");
	/*use _surface_tension as a temporary space*/
        f_hyper_surf(curve)->_surface_tension = avg_supersat;
}

static boolean supersat_on_hse(
        HYPER_SURF_ELEMENT *hse,
        HYPER_SURF *hs,
        RECT_GRID *gr,
        double *supersat,
        boolean pos_side)
{
        int dim = gr->dim;
        switch (dim)
        {
        case 2:
            return supersat_on_hse2d(hse,hs,gr,supersat,pos_side);
        case 3:
            return supersat_on_hse3d(hse,hs,gr,supersat,pos_side);
        default:
            return NO;
        }
}

static boolean supersat_on_hse3d(
        HYPER_SURF_ELEMENT *hse,
        HYPER_SURF *hs,
        RECT_GRID *gr,
        double *supersat,
        boolean pos_side)
{
        double crds1[MAXD], crds2[MAXD];
        double sup1,sup2;
        Locstate s1,s2;
	POINTER sl,sr;
        SURFACE *s = Surface_of_hs(hs);
        TRI  *t = Tri_of_hse(hse);
	POINT *point;
	int i,dim = gr->dim;

	*supersat = 0;
	for (i = 0; i < 3; ++i)
	{
	    point = Point_of_tri(t)[i];
	    FT_GetStatesAtPoint(point,hse,hs,&sl,&sr);
	    if(pos_side)
		*supersat += getStateSuper(sr);
	    else 
		*supersat += getStateSuper(sl);
	}
	*supersat /= 3.0;
        return YES;
}

static boolean supersat_on_hse2d(
        HYPER_SURF_ELEMENT *hse,
        HYPER_SURF *hs,
        RECT_GRID *gr,
        double *supersat,
        boolean pos_side)
{
        double crds1[MAXD], crds2[MAXD];
        double sup1,sup2;
        int i;
        Locstate s1,s2;
        CURVE *c = Curve_of_hs(hs);
        BOND  *b = Bond_of_hse(hse);

        if (b->start == c->start->posn)
            s1 = pos_side ? right_start_state(c) : left_start_state(c);
        else
            s1 = pos_side ? right_state(b->start) : left_state(b->start);
        if (b->end == c->end->posn)
            s2 = pos_side ? right_end_state(c) : left_end_state(c);
        else
            s2 = pos_side ? right_state(b->end) : left_state(b->end);
        sup1 = getStateSuper(s1); sup2 = getStateSuper(s2);
        for (i = 0; i < 2; ++i)
        {
            crds1[i] = Coords(b->start)[i];
            crds2[i] = Coords(b->end)[i];
        }
        *supersat = 0.5*(sup1+sup2);
        return YES;
}       /*end supersat_on_hse2d*/

extern void printDropletsStates(Front* front, char* outname)
{
        char filename[100];
        FILE *outfile;
	PARAMS* eqn_params = (PARAMS*)front->extra2;
	int dim = Dimension(front->interf);
        PARTICLE* particle_array = eqn_params->particle_array;
	int i,j,num_drops = eqn_params->num_drops;
            
        sprintf(filename,"%s/state.ts%s",outname,
                        right_flush(front->step,7));
        sprintf(filename,"%s-drops",filename);
        outfile = fopen(filename,"w");
	for (i = 0; i < num_drops; i++)
	{
	    fprintf(outfile,"%24.18g\n",particle_array[i].radius);
	    for (j = 0; j < dim; j++)
		fprintf(outfile,"%24.18g ",particle_array[i].center[j]);
	    fprintf(outfile,"\n");
	    for (j = 0; j < dim; j++)
		fprintf(outfile,"%24.18g ",particle_array[i].vel[j]);
	    fprintf(outfile,"\n");
	}
	fclose(outfile);
}

extern void gv_plot_scatter(Front* front)
{
	PARAMS* eqn_params = (PARAMS*)front->extra2;
	PARTICLE* particle_array = eqn_params->particle_array;
	int dim = front->rect_grid->dim;
	int i, j, n = eqn_params->num_drops;
	double *coords;
	char *outname = OutName(front);
	char fname[256];
	FILE* file;
	
	sprintf(fname,"%s/particle.off",outname);	
	file = fopen(fname,"w");
	fprintf(file,"appearance {linewidth 3}\n");
	fprintf(file,"VECT\n");
	fprintf(file,"%d %d %d\n\n",n,n,n);
	/*num of polyline, num of point, num of colors*/

        for (i = 0; i < n; i++)
        {
            fprintf(file,"1 ");/*num of points for each polyline*/
        }
	fprintf(file,"\n\n");

        for (i = 0; i < n; i++)
        {
            fprintf(file,"1 ");/*num of colors for each polyline*/
        }
        fprintf(file,"\n\n");
	
	for (i = 0; i < n; i++)
	{
	    coords = particle_array[i].center;
	    if (dim == 2)
	    {
	            fprintf(file,"%f %f 0\n",coords[0],coords[1]);/*coord for each points*/
	    }
	    else if (dim == 3)
	    {
		fprintf(file,"%f %f %f\n",coords[0],coords[1],coords[2]);
	    }
	}
	fprintf(file,"\n");
	for (i = 0; i < n; i++)
	    fprintf(file,"0 0 1 1\n");/*color for each point*/
	fclose(file);
}

extern void vtk_plot_scatter(Front* front)
{
        PARAMS* eqn_params = (PARAMS*)front->extra2;
        PARTICLE* particle_array = eqn_params->particle_array;
        int dim = front->rect_grid->dim;
        int i, j, count, n;
        double *coords;
        char *outname = OutName(front);
        char fname[256];
        FILE* file;

	count =0;
	for (i = 0; i < eqn_params->num_drops; i++)
	{
	    if (particle_array[i].radius != 0)
		count ++;
	}
	n = count; /*number of droplets with positive radius*/

        sprintf(fname,"%s/vtk/vtk.ts%s/",outname,
				right_flush(front->step,7));
	if (!create_directory(fname,NO))
        {
            printf("Cannot create directory %s\n",fname);
            clean_up(ERROR);
        }
	sprintf(fname,"%s/particle.vtk",fname);
        file = fopen(fname,"w");
        fprintf(file,"# vtk DataFile Version 3.0\n");
        fprintf(file,"%s\n","particles");
        fprintf(file,"ASCII\n");

	fprintf(file,"DATASET UNSTRUCTURED_GRID\n");
	fprintf(file,"POINTS %d FLOAT\n",n);

        for (i = 0; i < eqn_params->num_drops; i++)
        {
	    if (particle_array[i].radius == 0)
		continue;
            coords = particle_array[i].center;
            if (dim == 2)
                fprintf(file,"%f %f 0\n",coords[0],coords[1]);
            else if (dim == 3)
                fprintf(file,"%f %f %f\n",coords[0],coords[1],coords[2]);
        }
	fprintf(file,"POINT_DATA %d\n",n);
	fprintf(file,"SCALARS radius FLOAT\n");
	fprintf(file,"LOOKUP_TABLE default\n");
	for (i = 0; i < eqn_params->num_drops; i++)
        {
            if (particle_array[i].radius == 0)
                continue;
            fprintf(file,"%20.14f\n",particle_array[i].radius);
        }
	printf("%d number of droplets contained\n",n);
	fclose(file);
}

extern void vtk_plot_sample_traj(Front* front)
{
        PARAMS* eqn_params = (PARAMS*)front->extra2;
        int dim = front->rect_grid->dim;
        int i, j;
        double *coords;
        char *outname = OutName(front);
        char fname[256];
        FILE* file;
	/*static array for preserving trajectory*/
	static double traj[MAXD][MAX_STEP];
	static int step = 0;

	for(i = 0; i < dim; i++)
	{
	    traj[i][step] = eqn_params->particle_array[0].center[i];
	}
	if(step == 0)
	{  
	    step++;
	    return;
	}
        sprintf(fname,"%s/vtk/vtk.ts%s/",outname,
				right_flush(front->step,7));
	if (!create_directory(fname,NO))
        {
            printf("Cannot create directory %s\n",fname);
            clean_up(ERROR);
        }
	sprintf(fname,"%s/sample_particle.vtk",fname);
        file = fopen(fname,"w");
        fprintf(file,"# vtk DataFile Version 3.0\n");
        fprintf(file,"%s\n","sample_particles");
        fprintf(file,"ASCII\n");

	fprintf(file,"DATASET UNSTRUCTURED_GRID\n");
	fprintf(file,"POINTS %d FLOAT\n",step+1);

	for (i = 0; i < step+1; i++)
	{
            if (dim == 2)
                fprintf(file,"%f %f 0\n",traj[0][i],traj[1][i]);
            else if (dim == 3)
                fprintf(file,"%f %f %f\n",traj[0][i],traj[1][i],traj[2][i]);
	}
	fprintf(file,"CELLS %d %d\n",step,3*step);
	for (i = 0; i < step; i++)
	{
	    fprintf(file,"2 %d %d\n",i,i+1);
	}
	fprintf(file,"CELL_TYPES %d\n",step);
	for (i = 0; i < step; i++)
        {
            fprintf(file,"3\n");
        }
	step ++;
	fclose(file);
}

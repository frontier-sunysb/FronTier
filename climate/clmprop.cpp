/***************************************************************
FronTier is a set of libraries that implements differnt types of 
Front Traking algorithms. Front Tracking is a numerical method 
for the solution of partial differential equations whose solutions 
have discontinuities.  


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

****************************************************************/

#include <iFluid.h>
#include "climate.h"
static void condensation_point_propagate(Front*,POINTER,POINT*,POINT*,	
		HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
static void neumann_point_propagate(Front*,POINTER,POINT*,POINT*,	
		HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
static void dirichlet_point_propagate(Front*,POINTER,POINT*,POINT*,	
		HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);

extern  void melting_point_propagate(
        Front *front,
        POINTER wave,
        POINT *oldp,
        POINT *newp,
        HYPER_SURF_ELEMENT *oldhse,
        HYPER_SURF         *oldhs,
        double              dt,
        double              *V)
{
	PARAMS *eqn_params = (PARAMS*)front->extra2;
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
	case ICE_PARTICLE_BOUNDARY:
	    condensation_point_propagate(front,wave,oldp,newp,oldhse,
                                        oldhs,dt,V);
	    return;
        }
}	/* end melting_point_propagate */

static  void condensation_point_propagate(
        Front *front,
        POINTER wave,
        POINT *oldp,
        POINT *newp,
        HYPER_SURF_ELEMENT *oldhse,
        HYPER_SURF         *oldhs,
        double              dt,
        double              *V)
{
        if (debugging("point_propagate"))
            (void) printf("Entering condensation_point_propagate()\n");
	RECT_GRID *gr = computational_grid(front->interf);
	int i, index;
	int dim = gr->dim;
	int icoords[MAXD];
        double nor[MAXD], vel[MAXD], p1[MAXD];
        double *p0 = Coords(oldp);
	double R, supersat, center[MAXD];
	double dn, *h = front->rect_grid->h;
	double speed;
	double relative_dist;
	/*For updating the state*/
	PARAMS *eqn_params = (PARAMS*)front->extra2;
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	double *m_pre = iFparams->field->pres;
	double *m_vor = iFparams->field->vort;
	double *m_temp = eqn_params->field->temperature;
	double *m_vap = eqn_params->field->vapor;
	double *m_sup = eqn_params->field->supersat;

	double K = eqn_params->K;
	CURVE *curve = Curve_of_hs(oldhs);
	SURFACE *surface = Surface_of_hs(oldhs);
	STATE *state, *newst, *oldst;
	POINTER sl,sr;
	COMPONENT comp;

	FT_GetStatesAtPoint(oldp,oldhse,oldhs,&sl,&sr);
	if (negative_component(oldhs) == LIQUID_COMP)
        {
            comp = negative_component(oldhs);
            oldst = (STATE*)sl;
            newst = (STATE*)left_state(newp);
        }
        else if (positive_component(oldhs) == LIQUID_COMP)
        {
            comp = positive_component(oldhs);
            oldst = (STATE*)sr;
            newst = (STATE*)right_state(newp);
        }
        GetFrontNormal(oldp,oldhse,oldhs,nor,front);
	
        dn = grid_size_in_direction(nor,h,dim);
        for (i = 0; i < dim; ++i)
            p1[i] = p0[i] + nor[i]*dn;
	/*use _surface_tension as a temperary storage*/
	supersat = f_hyper_surf(oldhs)->_surface_tension;
	if (dim == 2)
	    R = spherical_radius(curve);
	else if (dim == 3)
	    R = spherical_radius(surface);
	else
	    printf("Unknown dim = %d\n",dim);
	speed =  K * supersat/R;
	for (i = 0; i < dim; ++i)
	{	
	    vel[i] = speed*nor[i];
	    Coords(newp)[i] = Coords(oldp)[i] + dt*vel[i];
	    newst->vel[i] = vel[i];
	    FT_RecordMaxFrontSpeed(i,fabs(vel[i]),NULL,Coords(newp),front);
	}
	FT_IntrpStateVarAtCoords(front,comp,p1,m_pre,
                        getStatePres,&newst->pres,&oldst->pres);
        FT_IntrpStateVarAtCoords(front,comp,p1,m_vap,
                        getStateVapor,&newst->vapor,&oldst->vapor);
        FT_IntrpStateVarAtCoords(front,comp,p1,m_sup,
                        getStateSuper,&newst->supersat,&oldst->supersat);
        FT_IntrpStateVarAtCoords(front,comp,p1,m_temp,
                        getStateTemperature,&newst->temperature,&oldst->temperature);
        relative_dist = fabs(speed*dt/dn);

        if (*front->max_scaled_propagation < relative_dist)
        {
            *front->max_scaled_propagation = relative_dist;
            if (relative_dist > 0.5)
            {
                (void) printf("WARNING: propagation too large!\n");
                (void) printf("relative distance = %f\n",relative_dist);
            }
            for (i = 0; i < dim; ++i)
                front->max_prop_point[i] = Coords(oldp)[i];
        }
	
        if (debugging("point_propagate"))
        {
            (void) printf("oldp = %f %f %f\n",Coords(oldp)[0],Coords(oldp)[1],
                        Coords(oldp)[2]);
            (void) printf("newp = %f %f %f\n",Coords(newp)[0],Coords(newp)[1],
                        Coords(newp)[2]);
            (void) printf("relative propagate = %f %f %f\n",
                        (Coords(newp)[0]-Coords(oldp)[0])/h[0],
                        (Coords(newp)[1]-Coords(oldp)[1])/h[1],
                        (Coords(newp)[2]-Coords(oldp)[2])/h[2]);
            (void) printf("relative dist = %f\n",relative_dist);
            (void) printf("Leaving condensation_point_propagate()\n");
        }
}

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
        PARAMS *params = (PARAMS*)front->extra2;
        double *Temp = params->field->temperature;
        COMPONENT ext_comp = exterior_component(front->interf);
        COMPONENT phase_comp;

	phase_comp = (negative_component(oldhs) == ext_comp) ?
                positive_component(oldhs) : negative_component(oldhs);
        for (i = 0; i < dim; ++i)
            Coords(newp)[i] = Coords(oldp)[i];
        FT_GetStatesAtPoint(oldp,oldhse,oldhs,(POINTER*)&sl,(POINTER*)&sr);
        state = (negative_component(oldhs) == phase_comp) ? sl : sr;
        s0 = state->temperature;
        FT_NormalAtPoint(oldp,front,nor,phase_comp);
        dn = grid_size_in_direction(nor,h,dim);
        for (i = 0; i < dim; ++i)
            p1[i] = p0[i] + nor[i]*dn;
        FT_IntrpStateVarAtCoords(front,phase_comp,p1,Temp,
                        getStateTemperature,&s1,&s0);
        for (i = 0; i < dim; ++i)
            p2[i] = p1[i] + nor[i]*dn;
        FT_IntrpStateVarAtCoords(front,phase_comp,p2,Temp,
                        getStateTemperature,&s2,&s1);
        state = (negative_component(oldhs) == phase_comp) ?
                (STATE*)left_state(newp) : (STATE*)right_state(newp);
        state->temperature = (4.0*s1 - s2)/3.0;
        return;
}	/* end neumann_point_propagate */

static void dirichlet_point_propagate(
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
        PARAMS *params = (PARAMS*)front->extra2;

	for (i = 0; i < dim; ++i)
            Coords(newp)[i] = Coords(oldp)[i];
        if (boundary_state(oldhs) != NULL)
        {
            bstate = (STATE*)boundary_state(oldhs);
            FT_GetStatesAtPoint(oldp,oldhse,oldhs,(POINTER*)&sl,(POINTER*)&sr);
            state =  (STATE*)left_state(newp);
            state->temperature = bstate->temperature;
            if (params->max_temperature < state->temperature)
                params->max_temperature = state->temperature;
            if (params->min_temperature > state->temperature)
                params->min_temperature = state->temperature;
            state =  (STATE*)right_state(newp);
            state->temperature = bstate->temperature;
            if (params->max_temperature < state->temperature)
                params->max_temperature = state->temperature;
            if (params->min_temperature > state->temperature)
                params->min_temperature = state->temperature;
        }
	else
	{
            FT_GetStatesAtPoint(oldp,oldhse,oldhs,(POINTER*)&sl,(POINTER*)&sr);
            state =  (STATE*)left_state(newp);
            if (params->max_temperature < state->temperature)
                params->max_temperature = state->temperature;
            if (params->min_temperature > state->temperature)
                params->min_temperature = state->temperature;
            state =  (STATE*)right_state(newp);
            if (params->max_temperature < state->temperature)
                params->max_temperature = state->temperature;
            if (params->min_temperature > state->temperature)
                params->min_temperature = state->temperature;
        }
        return;
}       /* end dirichlet_point_propagate */


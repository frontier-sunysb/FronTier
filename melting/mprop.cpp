/***************************************************************
FronTier is a set of libraries that implements different types of 
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
#include "melting.h"
static void ice_point_propagate_secondorder(Front*,POINTER,POINT*,POINT*,
                HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
static void ice_point_propagate(Front*,POINTER,POINT*,POINT*,	
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
        case GROWING_BODY_BOUNDARY:
	    if (eqn_params->num_scheme == UNSPLIT_EXPLICIT_CIM)
            	ice_point_propagate_secondorder(front,wave,oldp,newp,oldhse,
                                        oldhs,dt,V);
	    else
            	ice_point_propagate(front,wave,oldp,newp,oldhse,oldhs,dt,V);
            return;
        }
}	/* end melting_point_propagate */

static  void ice_point_propagate_secondorder(
        Front *front,
        POINTER wave,
        POINT *oldp,
        POINT *newp,
        HYPER_SURF_ELEMENT *oldhse,
        HYPER_SURF         *oldhs,
        double              dt,
        double              *V)
{
        double vel[MAXD];
        int i, dim = front->rect_grid->dim;
	double nor_l[MAXD],nor_s[MAXD];
        double p1[MAXD];
        double *p0 = Coords(oldp);
        double dn,*h = front->rect_grid->h;
	double grad_s,grad_l;
        double s0,s1;
        STATE *sl,*sr,*state;
	PARAMS *eqn_params = (PARAMS*)front->extra2;
	double *temperature = eqn_params->field->temperature;
        double rho_s = eqn_params->rho[0];
        double rho_l = eqn_params->rho[1];
        double k_s = eqn_params->k[0];
        double k_l = eqn_params->k[1];
        double Cp_s = eqn_params->Cp[0];
        double Cp_l = eqn_params->Cp[1];
        double L = eqn_params->L[0];
	double H;
	double relative_dist;
	double speed;

	int    m, icoord[MAXD], ic0[MAXD], ic1[MAXD], ic2[MAXD], id0, id1, id2;
	int    icL[MAXD], icS[MAXD], gmin[MAXD], s[MAXD];
	double P[MAXD], gd0[MAXD], gd1[MAXD], gd2[MAXD], T0, T1, T2;
	int       *top_gmax = FT_GridIntfcTopGmax(front);
	RECT_GRID *top_grid = FT_GridIntfcTopGrid(front);
        COMPONENT *top_comp = FT_GridIntfcTopComp(front);
	double coords[MAXD], crx_coords[MAXD], temp_nb[2];
	double a, b[MAXD], r, rL_min, rS_min, alpha[MAXD];
	double TxxL, TxxR, TxL, TxR, TxS;
	double Txx[MAXD], Txy[MAXD][MAXD], Tx[MAXD];
	boolean fr_crx_grid_seg;
	const GRID_DIRECTION dir[3][2] =
                {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};

	if (debugging("point_propagate"))
	    (void) printf("Entering ice_point_propagate_secondorder()\n");

	for(i = 0; i < dim; ++i) gmin[i] = 0;
	
	rect_in_which(p0, ic0, top_grid);
	id0 = d_index(ic0,top_gmax,dim);
	coords_of_grid_point(ic0,gd0,top_grid);
	// find minimum distance 
	rL_min = 1e5; rS_min = 1e5;
	r = 0;
	for(i = 0; i < dim; ++i) r += (gd0[i]-p0[i])*(gd0[i]-p0[i]);
	if(top_comp[id0] == LIQUID_COMP && r < rL_min)
	{
	    for(i = 0; i < dim; ++i) icL[i] = ic0[i];
	    rL_min = r;
	}
	if(top_comp[id0] == SOLID_COMP && r < rS_min)
	{
	    for(i = 0; i < dim; ++i) icS[i] = ic0[i];
	    rS_min = r;
	}
	for(i = 0; i < dim; ++i) ic1[i] = ic0[i] + 1;
	id1 = d_index(ic1,top_gmax,dim);
	coords_of_grid_point(ic1,gd1,top_grid);
	r = 0;
	for(i = 0; i < dim; ++i) r += (gd1[i]-p0[i])*(gd1[i]-p0[i]);
	if(top_comp[id1] == LIQUID_COMP && r < rL_min)
	{
	    for(i = 0; i < dim; ++i) icL[i] = ic1[i];
	    rL_min = r;
	}
	if(top_comp[id1] == SOLID_COMP && r < rS_min)
	{
	    for(i = 0; i < dim; ++i) icS[i] = ic1[i];
	    rS_min = r;
	}
	for(i = 0; i < dim; ++i) 
	{
	    for(m = 0; m < dim; ++m)
	        ic2[m] = (m==i) ? ic0[m] + 1 : ic0[m];
	    id2 = d_index(ic2,top_gmax,dim);
	    coords_of_grid_point(ic2,gd2,top_grid);
	    r = 0;
	    for(m = 0; m < dim; ++m) r += (gd2[m]-p0[m])*(gd2[m]-p0[m]);
	    if(top_comp[id2] == LIQUID_COMP && r < rL_min)
	    {
	        for(i = 0; i < dim; ++i) icL[i] = ic2[i];
	        rL_min = r;
    	    }
	    if(top_comp[id2] == SOLID_COMP && r < rS_min)
	    {
	        for(i = 0; i < dim; ++i) icS[i] = ic2[i];
	        rS_min = r;
	    }
	    for(m = 0; m < dim; ++m)
	        ic2[m] = (m==i) ? ic1[m] - 1 : ic0[m];
	    id2 = d_index(ic2,top_gmax,dim);
	    coords_of_grid_point(ic2,gd2,top_grid);
	    r = 0;
	    for(m = 0; m < dim; ++m) r += (gd2[m]-p0[m])*(gd2[m]-p0[m]);
	    if(top_comp[id2] == LIQUID_COMP && r < rL_min)
	    {
	        for(i = 0; i < dim; ++i) icL[i] = ic2[i];
	        rL_min = r;
    	    }
	    if(top_comp[id2] == SOLID_COMP && r < rS_min)
	    {
	        for(i = 0; i < dim; ++i) icS[i] = ic2[i];
	        rS_min = r;
	    }
	}
	// SOLID part
	for(i = 0; i < dim; ++i) ic0[i] = icS[i];
	id0 = d_index(ic0,top_gmax,dim);
	coords_of_grid_point(ic0,gd0,top_grid);
	for(i = 0; i < dim; ++i)
	{
	    T0 = temperature[id0];
	    Txx[i] = 0;
	    Tx[i]  = 0;
	    s[i]   = 0;
	    for(m = 0; m < 2; ++m) 
	    {
		fr_crx_grid_seg = FT_StateVarAtGridCrossing(front,
                                    ic0,dir[i][m],SOLID_COMP,
				    getStateTemperature,&temp_nb[m],crx_coords);
		if(fr_crx_grid_seg) 
		{
		    s[i] = 2*m-1;
		    a = alpha[i] = fabs(crx_coords[i] - gd0[i])/top_grid->h[i];
		    next_ip_in_dir(ic0,dir[i][1-m],ic1,gmin,top_gmax);
		    id1 = d_index(ic1, top_gmax, dim);
		    T1  = temperature[id1]-temp_nb[m];
		    next_ip_in_dir(ic1,dir[i][1-m],ic2,gmin,top_gmax);
		    id2 = d_index(ic2, top_gmax, dim);
		    T2  = temperature[id2]-temp_nb[m];
		    Txx[i] += (1-a)*T2 + 2*(a*a+a-1)/(1+a)*T1 - (1+a)*
				(T0-temp_nb[m]) - (T1+temp_nb[m] - T0);
		}
		else
		{
		    next_ip_in_dir(ic0,dir[i][m],ic1,gmin,top_gmax);
		    id1 = d_index(ic1, top_gmax, dim);
		    Txx[i] += (temperature[id1]-temperature[id0]);
		    if(m == 0)
		    	Tx[i]  += -(temperature[id1]-temperature[id0]);
		    else
		    	Tx[i]  += (temperature[id1]-temperature[id0]);
		}
	    }
	    if(s[i] == 0) 
	    {
		Tx[i] *= 0.5;
		b[i] = 0;
	    } 
	    else 
		b[i] = 0.5*s[i]*top_grid->h[i];
	    for(m = i+1; m < dim; ++m) {
		Txy[i][m] = 0.0; // not implemented!
	    }
	}
	for(i = 0; i < dim; ++i) 
	{
	    for(m = 0; m < dim; ++m) 
	    {
		if(m == i) 
		    Tx[i] += (p0[i] - b[i])*Txx[i]/top_grid->h[i];
		else
		    Tx[i] += (p0[i])*Txy[i][m]/top_grid->h[i];
	    }
	    Tx[i] /= top_grid->h[i];
	}
        GetFrontNormal(oldp,oldhse,oldhs,nor_s,front);
        dn = grid_size_in_direction(nor_s,h,dim);
        if (negative_component(oldhs) == SOLID_COMP)
            for (i = 0; i < dim; ++i)
                nor_s[i] *= -1.0;
	grad_s = 0.0;
	for(i = 0; i < dim; ++i)  
	    grad_s += Tx[i]*nor_s[i];
	// LIQUID part
	for(i = 0; i < dim; ++i) ic0[i] = icL[i];
	id0 = d_index(ic0,top_gmax,dim);
	coords_of_grid_point(ic0,gd0,top_grid);
	for(i = 0; i < dim; ++i)
	{
	    T0 = temperature[id0];
	    Txx[i] = 0;
	    Tx[i]  = 0;
	    s[i]   = 0;
	    for(m = 0; m < 2; ++m) 
	    {
		fr_crx_grid_seg = FT_StateVarAtGridCrossing(front,
                                    ic0,dir[i][m],LIQUID_COMP,
				    getStateTemperature,&temp_nb[m],crx_coords);
		if(fr_crx_grid_seg) 
		{
		    s[i] = 2*m-1;
		    a = alpha[i] = fabs(crx_coords[i] - gd0[i])/top_grid->h[i];
		    next_ip_in_dir(ic0,dir[i][1-m],ic1,gmin,top_gmax);
		    id1 = d_index(ic1, top_gmax, dim);
		    T1  = temperature[id1]-temp_nb[m];
		    next_ip_in_dir(ic1,dir[i][1-m],ic2,gmin,top_gmax);
		    id2 = d_index(ic2, top_gmax, dim);
		    T2  = temperature[id2]-temp_nb[m];
		    Txx[i] += (1-a)*T2 + 2*(a*a+a-1)/(1+a)*T1 - 
				(1+a)*(T0-temp_nb[m]) - (T1+temp_nb[m] - T0);
		}
		else
		{
		    next_ip_in_dir(ic0,dir[i][m],ic1,gmin,top_gmax);
		    id1 = d_index(ic1, top_gmax, dim);
		    Txx[i] += (temperature[id1]-temperature[id0]);
		    if(m == 0)
		    	Tx[i]  += -(temperature[id1]-temperature[id0]);
		    else
		    	Tx[i]  += (temperature[id1]-temperature[id0]);
		}
	    }
	    if(s[i] == 0) 
	    {
		Tx[i] *= 0.5;
		b[i] = 0;
	    } 
	    else 
		b[i] = 0.5*s[i]*top_grid->h[i];
	    for(m = i+1; m < dim; ++m) {
		Txy[i][m] = Txy[m][i] = 0.0; // not implemented!
	    }
	}
	for(i = 0; i < dim; ++i) 
	{
	    for(m = 0; m < dim; ++m) 
	    {
		if(m == i) 
		    Tx[i] += (p0[i] - b[i])*Txx[i]/top_grid->h[i];
		else
		    Tx[i] += (p0[i])*Txy[i][m]/top_grid->h[i];
	    }
	    Tx[i] /= top_grid->h[i];
	}
        GetFrontNormal(oldp,oldhse,oldhs,nor_l,front);
        if (negative_component(oldhs) == LIQUID_COMP)
            for (i = 0; i < dim; ++i)
                nor_l[i] *= -1.0;
	grad_l = 0.0;
	for(i = 0; i < dim; ++i)  
	    grad_l += Tx[i]*nor_l[i];

	H = Cp_l*k_l*grad_l + Cp_s*k_s*grad_s;
	speed = H/L/rho_s;
        for (i = 0; i < dim; ++i)
        {
	    if (H > 0)
		vel[i] = speed*nor_s[i];
	    else
		vel[i] = speed*nor_s[i];
	}

        for (i = 0; i < dim; ++i)
        {
            Coords(newp)[i] = Coords(oldp)[i] + dt*vel[i];
            FT_RecordMaxFrontSpeed(i,fabs(vel[i]),NULL,Coords(newp),front);
        }
	
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

	/* Update the state of the new interface point */
	state = (STATE*)left_state(newp);
	if (negative_component(oldhs) == LIQUID_COMP)
	    state->temperature = eqn_params->Ti[0];			
	else if (negative_component(oldhs) == SOLID_COMP)
            state->temperature = eqn_params->Ti[0];

	state = (STATE*)right_state(newp);
	if (positive_component(oldhs) == LIQUID_COMP)
	    state->temperature = eqn_params->Ti[0];
	else if (positive_component(oldhs) == SOLID_COMP)
            state->temperature = eqn_params->Ti[0];

	if (debugging("point_propagate"))
        {
	    (void) printf("oldp = %f %f %f\n",Coords(oldp)[0],Coords(oldp)[1],
    			Coords(oldp)[2]);
            (void) printf("newp = %f %f %f\n",Coords(newp)[0],Coords(newp)[1],
			Coords(newp)[2]);
	    (void) printf("H = %f  vel = %f %f %f\n",H,vel[0],vel[1],vel[2]);
	    (void) printf("relative propagate = %f %f %f\n",
			(Coords(newp)[0]-Coords(oldp)[0])/h[0],
			(Coords(newp)[1]-Coords(oldp)[1])/h[1],
			(Coords(newp)[2]-Coords(oldp)[2])/h[2]);
	    (void) printf("relative dist = %f\n",relative_dist); 
	    (void) printf("Leaving ice_point_propagate_secondorder()\n");
        }
	if (dim > 1)
	    clean_up(0);
}       /* ice_point_propagate_secondorder */

static  void ice_point_propagate(
        Front *front,
        POINTER wave,
        POINT *oldp,
        POINT *newp,
        HYPER_SURF_ELEMENT *oldhse,
        HYPER_SURF         *oldhs,
        double              dt,
        double              *V)
{
        double vel[MAXD];
        int i, dim = front->rect_grid->dim;
	double nor_l[MAXD],nor_s[MAXD];
        double p1[MAXD];
        double *p0 = Coords(oldp);
        double dn,*h = front->rect_grid->h;
	double grad_s,grad_l;
        double s0,s1;
        STATE *sl,*sr,*state;
	PARAMS *eqn_params = (PARAMS*)front->extra2;
	double *temperature = eqn_params->field->temperature;
        double rho_s = eqn_params->rho[0];
        double rho_l = eqn_params->rho[1];
        double k_s = eqn_params->k[0];
        double k_l = eqn_params->k[1];
        double Cp_s = eqn_params->Cp[0];
        double Cp_l = eqn_params->Cp[1];
        double L = eqn_params->L[0];
	double H;
	double relative_dist;
	double speed;

	if (debugging("point_propagate"))
	    (void) printf("Entering ice_point_propagate()\n");

        GetFrontNormal(oldp,oldhse,oldhs,nor_l,front);
        if (negative_component(oldhs) == LIQUID_COMP)
            for (i = 0; i < dim; ++i)
                nor_l[i] *= -1.0;
        dn = grid_size_in_direction(nor_l,h,dim);
        for (i = 0; i < dim; ++i)
            p1[i] = p0[i] + nor_l[i]*dn;
        FT_GetStatesAtPoint(oldp,oldhse,oldhs,(POINTER*)&sl,(POINTER*)&sr);
        state = (negative_component(oldhs) == LIQUID_COMP) ? sl : sr;
        s0 = state->temperature;
	FT_IntrpStateVarAtCoords(front,LIQUID_COMP,p1,temperature,
				getStateTemperature,&s1,(double*)state);
        grad_l = (s1 - s0)/dn;

        GetFrontNormal(oldp,oldhse,oldhs,nor_s,front);
        if (negative_component(oldhs) == SOLID_COMP)
            for (i = 0; i < dim; ++i)
                nor_s[i] *= -1.0;
        dn = grid_size_in_direction(nor_s,h,dim);
        for (i = 0; i < dim; ++i)
            p1[i] = p0[i] + nor_s[i]*dn;
        FT_GetStatesAtPoint(oldp,oldhse,oldhs,(POINTER*)&sl,(POINTER*)&sr);
        state = (negative_component(oldhs) == SOLID_COMP) ? sl : sr;
        s0 = state->temperature;
	FT_IntrpStateVarAtCoords(front,SOLID_COMP,p1,temperature,
				getStateTemperature,&s1,(double*)state);

        grad_s = (s1 - s0)/dn;
	H = Cp_l*k_l*grad_l + Cp_s*k_s*grad_s;
	speed = H/L/rho_s;
        for (i = 0; i < dim; ++i)
        {
	    if (H > 0)
		vel[i] = speed*nor_s[i];
	    else
		vel[i] = speed*nor_s[i];
	}

        for (i = 0; i < dim; ++i)
        {
            Coords(newp)[i] = Coords(oldp)[i] + dt*vel[i];
            FT_RecordMaxFrontSpeed(i,fabs(vel[i]),NULL,Coords(newp),front);
        }
	
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

	/* Update the state of the new interface point */
	state = (STATE*)left_state(newp);
	if (negative_component(oldhs) == LIQUID_COMP)
	    state->temperature = eqn_params->Ti[0];			
	else if (negative_component(oldhs) == SOLID_COMP)
            state->temperature = eqn_params->Ti[0];

	state = (STATE*)right_state(newp);
	if (positive_component(oldhs) == LIQUID_COMP)
	    state->temperature = eqn_params->Ti[0];
	else if (positive_component(oldhs) == SOLID_COMP)
            state->temperature = eqn_params->Ti[0];

	if (debugging("point_propagate"))
        {
	    (void) printf("oldp = %f %f %f\n",Coords(oldp)[0],Coords(oldp)[1],
    			Coords(oldp)[2]);
            (void) printf("newp = %f %f %f\n",Coords(newp)[0],Coords(newp)[1],
			Coords(newp)[2]);
	    (void) printf("H = %f  vel = %f %f %f\n",H,vel[0],vel[1],vel[2]);
	    (void) printf("relative propagate = %f %f %f\n",
			(Coords(newp)[0]-Coords(oldp)[0])/h[0],
			(Coords(newp)[1]-Coords(oldp)[1])/h[1],
			(Coords(newp)[2]-Coords(oldp)[2])/h[2]);
	    (void) printf("relative dist = %f\n",relative_dist); 
	    (void) printf("Leaving ice_point_propagate()\n");
        }
}       /* ice_point_propagate */

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


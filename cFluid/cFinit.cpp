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

#include <cFluid.h>

static void getRTState(STATE*,EQN_PARAMS*,double*,COMPONENT);
static void getRMState(STATE*,EQN_PARAMS*,double*,COMPONENT);
static void getBubbleState(STATE*,EQN_PARAMS*,double*,COMPONENT);
static void getAmbientState(STATE*,EQN_PARAMS*,double*,COMPONENT);
static void getBlastState(STATE*,EQN_PARAMS*,double*,COMPONENT);
static void getShockSineWaveState(STATE*,EQN_PARAMS*,double*,COMPONENT);
static void getAccuracySineWaveState(STATE*,EQN_PARAMS*,double*,COMPONENT);
static void behind_state(int,double,double*,int,STATE*,STATE*);
static double intfcPertHeight(FOURIER_POLY*,double*);
static double getStationaryVelocity(EQN_PARAMS*);

void G_CARTESIAN::initSinePertIntfc(
	LEVEL_FUNC_PACK *level_func_pack,
	char *inname)
{
	static FOURIER_POLY *level_func_params;
	FILE *infile = fopen(inname,"r");
	int i,j,num_modes;
	char mesg[100];
	EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
	static double	L[3], U[3];

	FT_ScalarMemoryAlloc((POINTER*)&level_func_params,sizeof(FOURIER_POLY));
	dim = level_func_params->dim = front->rect_grid->dim;
	
	ft_assign(L, front->rect_grid->L, 3*DOUBLE);
	ft_assign(U, front->rect_grid->U, 3*DOUBLE);
	
	level_func_params->L = L;
	level_func_params->U = U;

	level_func_pack->neg_component = GAS_COMP1;
	level_func_pack->pos_component = GAS_COMP2;
	CursorAfterString(infile,"Enter mean position of fluid interface:");
	fscanf(infile,"%lf",&level_func_params->z0);
	(void) printf("%f\n",level_func_params->z0);
	CursorAfterString(infile,"Enter number of sine modes:");
	fscanf(infile,"%d",&num_modes);
	(void) printf("%d\n",num_modes);
	level_func_params->num_modes = num_modes;
	FT_MatrixMemoryAlloc((POINTER*)&level_func_params->nu,num_modes,
				dim-1,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&level_func_params->phase,num_modes,
				sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&level_func_params->A,num_modes,
				sizeof(double));
	for (i = 0; i < num_modes; ++i)
	{
	    sprintf(mesg,"Enter frequency of mode %d:",i+1);
	    CursorAfterString(infile,mesg);
	    for (j = 0; j < dim-1; ++j)
	    {
	    	fscanf(infile,"%lf",&level_func_params->nu[i][j]);
		(void) printf("%f ",level_func_params->nu[i][j]);
	    }
	    (void) printf("\n");
	    sprintf(mesg,"Enter amplitude of mode %d:",i+1);
	    CursorAfterString(infile,mesg);
	    fscanf(infile,"%lf",&level_func_params->A[i]);
	    (void) printf("%f\n",level_func_params->A[i]);
	    sprintf(mesg,"Enter phase of mode %d:",i+1);
	    CursorAfterString(infile,mesg);
	    fscanf(infile,"%lf",&level_func_params->phase[i]);
	    (void) printf("%f\n",level_func_params->phase[i]);
	}

	eqn_params->level_func_params = (POINTER)level_func_params;
	level_func_pack->func_params = (POINTER)level_func_params;
	level_func_pack->func = level_wave_func;
	level_func_pack->wave_type = FIRST_PHYSICS_WAVE_TYPE;
	fclose(infile);
}	/* end initRayleiTaylor */

static double intfcPertHeight(
	FOURIER_POLY *wave_params,
	double *coords)
{
	double arg,z,k,phase;
	int i,j,num_modes,dim;
	double *L = wave_params->L;
	double *U = wave_params->U;

	dim = wave_params->dim;
	num_modes = wave_params->num_modes;
	z = wave_params->z0;

	for (i = 0; i < num_modes; ++i)
	{
	    arg = 0.0;
	    for (j = 0; j < dim-1; ++j)
	    {
		k = wave_params->nu[i][j]*2.0*PI/(U[j]-L[j]);
		arg += k*coords[j];
	    }
	    phase = wave_params->phase[i]*PI/180.0;
	    arg -= phase;
	    z += wave_params->A[i]*sin(arg);
	}
	return z;
}	/* end intfcPertHeight */

void G_CARTESIAN::setRayleiTaylorParams(char *inname)
{
	int i;
	FILE *infile = fopen(inname,"r");
	EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
	double		pinf,einf,gamma;
	char s[100], str[256];

	sprintf(str, "Enter gamma, pinf, einf of the fluid with comp %d:", 
	             GAS_COMP1);
	CursorAfterString(infile,str);
	fscanf(infile,"%lf %lf %lf",&gamma,&pinf,&einf);
	(void) printf("%f %f %f\n",gamma,pinf,einf);
	(eqn_params->eos[GAS_COMP1]).gamma = gamma;
	(eqn_params->eos[GAS_COMP1]).pinf = pinf;
	(eqn_params->eos[GAS_COMP1]).einf = einf;
	
	sprintf(str, "Enter gamma, pinf, einf of the fluid with comp %d:", 
		     GAS_COMP2);
	CursorAfterString(infile,str);
	fscanf(infile,"%lf %lf %lf",&gamma,&pinf,&einf);
	(void) printf("%f %f %f\n",gamma,pinf,einf);
	(eqn_params->eos[GAS_COMP2]).gamma = gamma;
	(eqn_params->eos[GAS_COMP2]).pinf = pinf;
	(eqn_params->eos[GAS_COMP2]).einf = einf;

	CursorAfterString(infile,"Enter density of top fluid:");
	fscanf(infile,"%lf",&eqn_params->rho2);
	(void) printf("%f\n",eqn_params->rho2);
	CursorAfterString(infile,"Enter density of bottom fluid:");
	fscanf(infile,"%lf",&eqn_params->rho1);
	(void) printf("%f\n",eqn_params->rho1);
	CursorAfterString(infile,"Enter gravity:");
	for (i = 0; i < dim; ++i)
	{
	    fscanf(infile,"%lf",&eqn_params->gravity[i]);
	    (void) printf("%f ",eqn_params->gravity[i]);
	}
	(void) printf("\n");
	CursorAfterString(infile,"Enter pressure at interface:");
	fscanf(infile,"%lf",&eqn_params->p0);
	(void) printf("%f\n",eqn_params->p0);
	CursorAfterString(infile,"Type yes to track the interface:");
	fscanf(infile,"%s",s);
	(void) printf("%s\n",s);
	if (s[0] == 'y' || s[0] == 'Y')
	    eqn_params->tracked = YES;
	else
	    eqn_params->tracked = NO;
	fclose(infile);
}	/* end initRayleiTaylorParams */

void G_CARTESIAN::initRayleiTaylorStates()
{
	int i,j,k,l,index;
	EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
	double coords[MAXD];
	COMPONENT comp;
	STATE *sl,*sr,state;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
	INTERFACE *intfc = front->interf;
	double *dens = field.dens;
	double *engy = field.engy;
	double *pres = field.pres;
	double **momn = field.momn;

        next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
	    FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
	    getRTState(sl,eqn_params,Coords(p),negative_component(hs));
	    getRTState(sr,eqn_params,Coords(p),positive_component(hs));
	}
	FT_MakeGridIntfc(front);
	setDomain();

	switch (dim)
	{
	case 2:
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index2d(i,j,top_gmax);
		comp = top_comp[index];
		getRectangleCenter(index,coords);
	    	getRTState(&state,eqn_params,coords,comp);
		dens[index] = state.dens;
		pres[index] = state.pres;
		engy[index] = state.engy;
		for (l = 0; l < dim; ++l)
		    momn[l][index] = state.momn[l];
	    }
	    break;
	case 3:
	    for (k = imin[2]; k <= imax[2]; ++k)
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index3d(i,j,k,top_gmax);
		comp = top_comp[index];
		getRectangleCenter(index,coords);
	    	getRTState(&state,eqn_params,coords,comp);
		dens[index] = state.dens;
		pres[index] = state.pres;
		engy[index] = state.engy;
		for (l = 0; l < dim; ++l)
		    momn[l][index] = state.momn[l];
	    }
	    break;
	}
	scatMeshStates();
}	/* end initRayleiTaylorStates */

void G_CARTESIAN::setRichtmyerMeshkovParams(char *inname)
{
	FILE *infile = fopen(inname,"r");
	EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
	char s[100], str[256];
	double	pinf, einf, gamma;

	sprintf(str,"Enter gamma, pinf, einf of the fluid with comp %d:", 
	             GAS_COMP1);
	CursorAfterString(infile,str);
	fscanf(infile,"%lf %lf %lf",&gamma,&pinf,&einf);
	(void) printf("%f %f %f\n",gamma,pinf,einf);
	(eqn_params->eos[GAS_COMP1]).gamma = gamma;
	(eqn_params->eos[GAS_COMP1]).pinf = pinf;
	(eqn_params->eos[GAS_COMP1]).einf = einf;
	
	sprintf(str,"Enter gamma, pinf, einf of the fluid with comp %d:", 
		     GAS_COMP2);
	CursorAfterString(infile,str);
	fscanf(infile,"%lf %lf %lf",&gamma,&pinf,&einf);
	(void) printf("%f %f %f\n",gamma,pinf,einf);
	(eqn_params->eos[GAS_COMP2]).gamma = gamma;
	(eqn_params->eos[GAS_COMP2]).pinf = pinf;
	(eqn_params->eos[GAS_COMP2]).einf = einf;

	CursorAfterString(infile,"Enter density of top fluid:");
	fscanf(infile,"%lf",&eqn_params->rho2);
	(void) printf("%f\n",eqn_params->rho2);
	CursorAfterString(infile,"Enter density of bottom fluid:");
	fscanf(infile,"%lf",&eqn_params->rho1);
	(void) printf("%f\n",eqn_params->rho1);
	CursorAfterString(infile,"Enter pressure at interface:");
	fscanf(infile,"%lf",&eqn_params->p0);
	(void) printf("%f\n",eqn_params->p0);
	CursorAfterString(infile,"Enter Mach number of shock:");
	fscanf(infile,"%lf",&eqn_params->Mach_number);
	(void) printf("%f\n",eqn_params->Mach_number);
	CursorAfterString(infile,"Enter position of shock:");
	fscanf(infile,"%lf",&eqn_params->shock_position);
	(void) printf("%f\n",eqn_params->shock_position);
	CursorAfterString(infile,"Enter direction of shock:");
	fscanf(infile,"%d",&eqn_params->shock_dir);
	(void) printf("%d\n",eqn_params->shock_dir);
	CursorAfterString(infile,"Type yes to track the interface:");
	fscanf(infile,"%s",s);
	(void) printf("%s\n",s);
	if (s[0] == 'y' || s[0] == 'Y')
	    eqn_params->tracked = YES;
	else
	    eqn_params->tracked = NO;
	if (CursorAfterStringOpt(infile,"Type yes for stationary contact: "))
        {
	    fscanf(infile,"%s",s);
	    (void) printf("%s\n",s);
            if (s[0] == 'y' || s[0] == 'Y')
	    {
                eqn_params->contact_stationary = YES;
		eqn_params->contact_vel = getStationaryVelocity(eqn_params);
	    }
            else
                eqn_params->contact_stationary = NO;
        }
	fclose(infile);
}	/* end setRayleiTaylorParams */

void G_CARTESIAN::initRichtmyerMeshkovStates()
{
	int i,j,k,l,index;
	EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
	double coords[MAXD];
	COMPONENT comp;
	STATE *sl,*sr,state;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
	INTERFACE *intfc = front->interf;
	double *dens = field.dens;
	double *engy = field.engy;
	double *pres = field.pres;
	double **momn = field.momn;

        next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
	    FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
	    getRMState(sl,eqn_params,Coords(p),negative_component(hs));
	    getRMState(sr,eqn_params,Coords(p),positive_component(hs));
	}
	FT_MakeGridIntfc(front);
	setDomain();

	switch (dim)
	{
	case 2:
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index2d(i,j,top_gmax);
		comp = top_comp[index];
		getRectangleCenter(index,coords);
	    	getRMState(&state,eqn_params,coords,comp);
		dens[index] = state.dens;
		pres[index] = state.pres;
		engy[index] = state.engy;
		for (l = 0; l < dim; ++l)
		    momn[l][index] = state.momn[l];
	    }
	    break;
	case 3:
	    for (k = imin[2]; k <= imax[2]; ++k)
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index3d(i,j,k,top_gmax);
		comp = top_comp[index];
		getRectangleCenter(index,coords);
	    	getRMState(&state,eqn_params,coords,comp);
		dens[index] = state.dens;
		pres[index] = state.pres;
		engy[index] = state.engy;
		for (l = 0; l < dim; ++l)
		    momn[l][index] = state.momn[l];
	    }
	    break;
	}
	scatMeshStates();
}	/* end initRichtmyerMeshkovStates */

//EOS dep
static void getRTState(
	STATE *state,
	EQN_PARAMS *eqn_params,
	double *coords,
	COMPONENT comp)
{
	FOURIER_POLY	*wave_params;
	EOS_PARAMS	*eos;
	double		z_intfc;
	double 		rho1 = eqn_params->rho1;
	double 		rho2 = eqn_params->rho2;
	double 		p0 = eqn_params->p0;
	double 		*g = eqn_params->gravity;
	double 		dz, gz, c2, gamma;
	int    		i,dim;
	double		tmp;

	eos = &(eqn_params->eos[comp]);
	state->eos = eos;
	gamma = eos->gamma;

	wave_params = (FOURIER_POLY*)eqn_params->level_func_params;
	dim = wave_params->dim;
	z_intfc = intfcPertHeight(wave_params,coords);
	dz = coords[dim-1] - z_intfc;
	gz = g[dim-1];

	/* Constant density */
	for (i = 0; i < dim; ++i)
	    state->momn[i] = 0.0;
	switch (comp)
	{
	case GAS_COMP1:
	    c2 = gamma*(p0+eos->pinf)/rho1;
	    tmp = exp(gamma*gz/c2*dz);
	    state->dens = rho1*tmp;
	    state->pres = state->dens*c2/gamma - eos->pinf;
	    state->engy = EosInternalEnergy(state);
	    break;
	case GAS_COMP2:
	    c2 = gamma*(p0+eos->pinf)/rho2;
	    tmp = exp(gamma*gz/c2*dz);
	    state->dens = rho2*tmp;
	    state->pres = state->dens*c2/gamma - eos->pinf;
	    state->engy = EosInternalEnergy(state);
	    break;
	case EXT_COMP:
	    state->dens = 0.0;
	    state->pres = 0.0;
	    state->engy = 0.0;
	    break;
	default:
	    printf("ERROR: Unknown component %d in getRTState()!\n",comp);
	    clean_up(ERROR);
	}
}	/* end getRTState */

static void getRMState(
	STATE *state,
	EQN_PARAMS *eqn_params,
	double *coords,
	COMPONENT comp)
{
	FOURIER_POLY *wave_params;
	EOS_PARAMS	*eos;
	double rho1 = eqn_params->rho1;
	double rho2 = eqn_params->rho2;
	double p0 = eqn_params->p0;
	double shock_position = eqn_params->shock_position;
	double Mach_number = eqn_params->Mach_number;
	double shock_speed;
	double csp = eqn_params->contact_vel;
	int shock_dir = eqn_params->shock_dir;
	int i,dim;
 
	if (debugging("rm_state"))
	    printf("Entering getRMState(), coords = %f %f\n",
				coords[0],coords[1]);
	wave_params = (FOURIER_POLY*)eqn_params->level_func_params;
	dim = wave_params->dim;

	/* Constant density */
	for (i = 0; i < dim; ++i)
	    state->vel[i] = state->momn[i] = 0.0;
	state->dim = dim;
	eos = &(eqn_params->eos[comp]);
	state->eos = eos;
	
	switch (comp)
	{
	case GAS_COMP1:
	    state->dens = rho1;
	    state->pres = p0;
	    state->engy = EosInternalEnergy(state);
	    break;
	case GAS_COMP2:
	    state->dens = rho2;
	    state->pres = p0;
	    state->engy = EosInternalEnergy(state);
	    break;
	case EXT_COMP:
	    state->dens = 0.0;
	    state->pres = 0.0;
	    state->engy = 0.0;
	    break;
	default:
	    printf("ERROR: Unknown component %d in getRTState()!\n",comp);
	    clean_up(ERROR);
	}
	if (debugging("rm_state"))
	{
	    printf("Before calling behind_state()\n");
	    printf("state = %f %f %f\n",state->dens,state->pres,
					state->vel[0]);
	}
	if ((shock_dir ==  1 && coords[dim-1] < shock_position) ||
	    (shock_dir == -1 && coords[dim-1] > shock_position))
	{
	    behind_state(SHOCK_MACH_NUMBER,Mach_number,
			&shock_speed,shock_dir,state,state);	
	    state->engy = EosEnergy(state);
	    if (debugging("rm_state"))
	    {
	    	printf("After calling behind_state()\n");
	    	printf("state = %f %f %f\n",state->dens,state->pres,
			state->vel[0]);
	    }
	}
	state->vel[dim-1] -= csp;
	state->momn[dim-1] = state->vel[dim-1]*state->dens;
	state->engy = EosEnergy(state);

}	/* end getRMState */

static void behind_state(
	int		which_parameter,
	double		parameter,
	double		*shock_speed,
	int		shock_dir,
	STATE		*ahead_state,
	STATE		*behind_state)
{
	double		r0, p0, u0;		/* ahead state */
	double		r1, p1, u1;		/* behind state */
	double		U;			/* shock speed */
	double		M0n;			/* shock mack number,
						   relative to ahead flow */
	double		M0nsq;			/* steady normal ahead Mach
						   number squared */
	int		dim;

	dim = ahead_state->dim;
	r0  = ahead_state->dens;
	p0  = ahead_state->pres;
	u0  = ahead_state->vel[dim-1]*shock_dir;

	switch(which_parameter)
	{
	case SHOCK_MACH_NUMBER:
	    M0n = parameter;
	    *shock_speed = U = u0 + M0n*EosSoundSpeed(ahead_state);
	    M0nsq = sqr(M0n);
	    p1 = EosMaxBehindShockPres(M0nsq,ahead_state);
	    u1 =  u0 + (p0 - p1) / (r0*(u0 - U)); 
	    r1 = r0*((u0 - U)/(u1 - U));
	    if (debugging("rm_state"))
	    {
		printf("M0n = %f  shock_speed = %f\n",M0n,*shock_speed);
		printf("p1 = %f  u1 = %f  r1 = %f\n",p1,u1,r1);
	    }
	    break;
	default:
	    screen("ERROR in behind_state(), "
	           "unknown parameter %d\n",which_parameter);
	    clean_up(ERROR);
	}	
	behind_state->dens = r1;
	behind_state->pres = p1;
	behind_state->vel[dim-1] = u1*shock_dir;
	behind_state->momn[dim-1] = r1*u1*shock_dir;
}		/*end behind_state */


void G_CARTESIAN::initCirclePlaneIntfc(
	LEVEL_FUNC_PACK *level_func_pack,
	char *inname)
{
	EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
	FILE *infile = fopen(inname,"r");
	static CIRCLE_PARAMS *circle_params;
	int i,dim;
	PROB_TYPE prob_type = eqn_params->prob_type;

	FT_ScalarMemoryAlloc((POINTER*)&circle_params,sizeof(CIRCLE_PARAMS));
        circle_params->dim = dim = front->rect_grid->dim;
        circle_params->add_plan_surf = NO;
        CursorAfterString(infile,"Enter the center of the circle:");
        for (i = 0; i < dim; ++i)
	{
            fscanf(infile,"%lf",&circle_params->cen[i]);
            (void) printf("%f ",circle_params->cen[i]);
	}
	(void) printf("\n");
        CursorAfterString(infile,"Enter radius of the circle:");
        fscanf(infile,"%lf",&circle_params->R);
        (void) printf("%f\n",circle_params->R);

        if (prob_type == BUBBLE_SURFACE)
        {
            CursorAfterString(infile,"Enter height of the surface:");
            fscanf(infile,"%lf",&circle_params->H);
            (void) printf("%f\n",circle_params->H);
            circle_params->add_plan_surf = YES;
        }
	level_func_pack->func_params = (POINTER)circle_params;
	eqn_params->level_func_params = (POINTER)circle_params;

	switch (prob_type)
	{
	case TWO_FLUID_BUBBLE:
        case BUBBLE_SURFACE:
            level_func_pack->neg_component = GAS_COMP1;
            level_func_pack->pos_component = GAS_COMP2;
            level_func_pack->func = level_circle_func;
            level_func_pack->wave_type = FIRST_PHYSICS_WAVE_TYPE;
            break;
        case FLUID_SOLID_CIRCLE:
            level_func_pack->neg_component = SOLID_COMP;
            level_func_pack->pos_component = GAS_COMP1;
            level_func_pack->func = level_circle_func;
	    level_func_pack->wave_type = MOVABLE_BODY_BOUNDARY;
            break;
	default:
	    (void) printf("ERROR Wrong type in initCirclePlaneIntfc()\n");
	    clean_up(ERROR);
	}
	fclose(infile);	
}	/* end initCirclePlaneIntfc */


void G_CARTESIAN::initBubbleStates()
{
	int i,j,k,l,index;
	EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
	double coords[MAXD];
	COMPONENT comp;
	STATE *sl,*sr,state;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
	INTERFACE *intfc = front->interf;
	double *dens = field.dens;
	double *engy = field.engy;
	double *pres = field.pres;
	double **momn = field.momn;

        next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
	    FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
	    getBubbleState(sl,eqn_params,Coords(p),negative_component(hs));
	    getBubbleState(sr,eqn_params,Coords(p),positive_component(hs));
	}
	FT_MakeGridIntfc(front);
	setDomain();

	switch (dim)
	{
	case 2:
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index2d(i,j,top_gmax);
		comp = top_comp[index];
		getRectangleCenter(index,coords);
	    	getBubbleState(&state,eqn_params,coords,comp);
		dens[index] = state.dens;
		pres[index] = state.pres;
		engy[index] = state.engy;
		for (l = 0; l < dim; ++l)
		    momn[l][index] = state.momn[l];
	    }
	    break;
	case 3:
	    for (k = imin[2]; k <= imax[2]; ++k)
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index3d(i,j,k,top_gmax);
		comp = top_comp[index];
		getRectangleCenter(index,coords);
	    	getBubbleState(&state,eqn_params,coords,comp);
		dens[index] = state.dens;
		pres[index] = state.pres;
		engy[index] = state.engy;
		for (l = 0; l < dim; ++l)
		    momn[l][index] = state.momn[l];
	    }
	    break;
	}
	scatMeshStates();
}	/* end initRayleiTaylorStates */

static void getBubbleState(
	STATE *state,
	EQN_PARAMS *eqn_params,
	double *coords,
	COMPONENT comp)
{
	static CIRCLE_PARAMS	*circle_params;
	EOS_PARAMS	*eos;
	double		z0;
	double 		p1 = eqn_params->p1;
	double 		p2 = eqn_params->p2;
	double 		rho1 = eqn_params->rho1;
	double 		rho2 = eqn_params->rho2;
	double 		*g = eqn_params->gravity;
	double 		dz, gz, c2, gamma;
	int    		i,dim;
	double		tmp;

	eos = &(eqn_params->eos[comp]);
	state->eos = eos;
	gamma = eos->gamma;

	circle_params = (CIRCLE_PARAMS*)eqn_params->level_func_params;
	dim = circle_params->dim;
	z0 = circle_params->cen[dim-1];
	dz = coords[dim-1] - z0;
	gz = g[dim-1];

	/* Constant density */
	for (i = 0; i < dim; ++i)
	    state->momn[i] = 0.0;
	switch (comp)
	{
	case GAS_COMP1:
	    c2 = gamma*(p1+eos->pinf)/rho1;
	    tmp = exp(gamma*gz/c2*dz);
	    state->dens = rho1*tmp;
	    state->pres = state->dens*c2/gamma - eos->pinf;
	    state->engy = EosInternalEnergy(state);
	    break;
	case GAS_COMP2:
	    c2 = gamma*(p2+eos->pinf)/rho2;
	    tmp = exp(gamma*gz/c2*dz);
	    state->dens = rho2*tmp;
	    state->pres = state->dens*c2/gamma - eos->pinf;
	    state->engy = EosInternalEnergy(state);
	    break;
	case EXT_COMP:
	    state->dens = 0.0;
	    state->pres = 0.0;
	    state->engy = 0.0;
	    break;
	default:
	    printf("ERROR: Unknown component %d in getRTState()!\n",comp);
	    clean_up(ERROR);
	}
}	/* end getBubbleState */

void G_CARTESIAN::setBubbleParams(char *inname)
{
	int i;
	FILE *infile = fopen(inname,"r");
	EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
	double		pinf,einf,gamma;
	char s[100], str[256];

	sprintf(str, "Enter gamma, pinf, einf of the fluid with comp %d:", 
	             GAS_COMP1);
	CursorAfterString(infile,str);
	fscanf(infile,"%lf %lf %lf",&gamma,&pinf,&einf);
	(void) printf("%f %f %f\n",gamma,pinf,einf);
	(eqn_params->eos[GAS_COMP1]).gamma = gamma;
	(eqn_params->eos[GAS_COMP1]).pinf = pinf;
	(eqn_params->eos[GAS_COMP1]).einf = einf;
	
	sprintf(str, "Enter gamma, pinf, einf of the fluid with comp %d:", 
		     GAS_COMP2);
	CursorAfterString(infile,str);
	fscanf(infile,"%lf %lf %lf",&gamma,&pinf,&einf);
	(void) printf("%f %f %f\n",gamma,pinf,einf);
	(eqn_params->eos[GAS_COMP2]).gamma = gamma;
	(eqn_params->eos[GAS_COMP2]).pinf = pinf;
	(eqn_params->eos[GAS_COMP2]).einf = einf;

	CursorAfterString(infile,"Enter density and pressure inside bubble:");
	fscanf(infile,"%lf %lf",&eqn_params->rho2,&eqn_params->p2);
	(void) printf("%f %f\n",eqn_params->rho2,eqn_params->p2);
	CursorAfterString(infile,"Enter density and pressure outside bubble:");
	fscanf(infile,"%lf %lf",&eqn_params->rho1,&eqn_params->p1);
	(void) printf("%f %f\n",eqn_params->rho1,eqn_params->p1);
	CursorAfterString(infile,"Enter gravity:");
	for (i = 0; i < dim; ++i)
	{
	    fscanf(infile,"%lf",&eqn_params->gravity[i]);
	    (void) printf("%f ",eqn_params->gravity[i]);
	}
	(void) printf("\n");
	CursorAfterString(infile,"Type yes to track the interface:");
	fscanf(infile,"%s",s);
	(void) printf("%s\n",s);
	if (s[0] == 'y' || s[0] == 'Y')
	    eqn_params->tracked = YES;
	else
	    eqn_params->tracked = NO;
	fclose(infile);
}	/* end setBubbleParams */

void G_CARTESIAN::setImplosionParams(char *inname)
{
	int i;
	FILE *infile = fopen(inname,"r");
	EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
	double		pinf,einf,gamma;
	char s[100], str[256];

	sprintf(str, "Enter gamma, pinf, einf of the fluid with comp %d:", 
	             GAS_COMP1);
	CursorAfterString(infile,str);
	fscanf(infile,"%lf %lf %lf",&gamma,&pinf,&einf);
	(void) printf("%f %f %f\n",gamma,pinf,einf);
	(eqn_params->eos[GAS_COMP1]).gamma = gamma;
	(eqn_params->eos[GAS_COMP1]).pinf = pinf;
	(eqn_params->eos[GAS_COMP1]).einf = einf;
	
	sprintf(str, "Enter gamma, pinf, einf of the fluid with comp %d:", 
		     GAS_COMP2);
	CursorAfterString(infile,str);
	fscanf(infile,"%lf %lf %lf",&gamma,&pinf,&einf);
	(void) printf("%f %f %f\n",gamma,pinf,einf);
	(eqn_params->eos[GAS_COMP2]).gamma = gamma;
	(eqn_params->eos[GAS_COMP2]).pinf = pinf;
	(eqn_params->eos[GAS_COMP2]).einf = einf;

	CursorAfterString(infile,"Enter density and pressure inside contact:");
	fscanf(infile,"%lf %lf",&eqn_params->rho2,&eqn_params->p2);
	(void) printf("%f %f\n",eqn_params->rho2,eqn_params->p2);
	CursorAfterString(infile,"Enter density and pressure outside contact:");
	fscanf(infile,"%lf %lf",&eqn_params->rho1,&eqn_params->p1);
	(void) printf("%f %f\n",eqn_params->rho1,eqn_params->p1);
	CursorAfterString(infile,"Enter density and pressure at the boundary:");
	fscanf(infile,"%lf %lf",&eqn_params->rho0,&eqn_params->p0);
	(void) printf("%f %f\n",eqn_params->rho0,eqn_params->p0);
	CursorAfterString(infile,"Enter gravity:");
	for (i = 0; i < dim; ++i)
	{
	    fscanf(infile,"%lf",&eqn_params->gravity[i]);
	    (void) printf("%f ",eqn_params->gravity[i]);
	}
	(void) printf("\n");
	CursorAfterString(infile,"Type yes to track the interface:");
	fscanf(infile,"%s",s);
	(void) printf("%s\n",s);
	if (s[0] == 'y' || s[0] == 'Y')
	    eqn_params->tracked = YES;
	else
	    eqn_params->tracked = NO;
	fclose(infile);
}	/* end setBubbleParams */

void G_CARTESIAN::setMTFusionParams(char *inname)
{
	int i;
	FILE *infile = fopen(inname,"r");
	EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
	double		pinf,einf,gamma;
	char s[100], str[256];

	sprintf(str, "Enter gamma, pinf, einf of the fluid with comp %d:", 
	             GAS_COMP1);
	CursorAfterString(infile,str);
	fscanf(infile,"%lf %lf %lf",&gamma,&pinf,&einf);
	(void) printf("%f %f %f\n",gamma,pinf,einf);
	(eqn_params->eos[GAS_COMP1]).gamma = gamma;
	(eqn_params->eos[GAS_COMP1]).pinf = pinf;
	(eqn_params->eos[GAS_COMP1]).einf = einf;
	
	sprintf(str, "Enter gamma, pinf, einf of the fluid with comp %d:", 
		     GAS_COMP2);
	CursorAfterString(infile,str);
	fscanf(infile,"%lf %lf %lf",&gamma,&pinf,&einf);
	(void) printf("%f %f %f\n",gamma,pinf,einf);
	(eqn_params->eos[GAS_COMP2]).gamma = gamma;
	(eqn_params->eos[GAS_COMP2]).pinf = pinf;
	(eqn_params->eos[GAS_COMP2]).einf = einf;

	CursorAfterString(infile,"Enter density and pressure inside contact:");
	fscanf(infile,"%lf %lf",&eqn_params->rho2,&eqn_params->p2);
	(void) printf("%f %f\n",eqn_params->rho2,eqn_params->p2);
	CursorAfterString(infile,"Enter density and pressure outside contact:");
	fscanf(infile,"%lf %lf",&eqn_params->rho1,&eqn_params->p1);
	(void) printf("%f %f\n",eqn_params->rho1,eqn_params->p1);
	CursorAfterString(infile,"Enter gravity:");
	for (i = 0; i < dim; ++i)
	{
	    fscanf(infile,"%lf",&eqn_params->gravity[i]);
	    (void) printf("%f ",eqn_params->gravity[i]);
	}
	(void) printf("\n");
	CursorAfterString(infile,"Type yes to track the interface:");
	fscanf(infile,"%s",s);
	(void) printf("%s\n",s);
	if (s[0] == 'y' || s[0] == 'Y')
	    eqn_params->tracked = YES;
	else
	    eqn_params->tracked = NO;
	fclose(infile);
}	/* end setMTFusionParams */

void G_CARTESIAN::initImplosionIntfc(
	LEVEL_FUNC_PACK *level_func_pack,
	char *inname)
{
	FILE *infile = fopen(inname,"r");
	static CIRCLE_PARAMS *circle_params;
	COMPONENT neg_comp,pos_comp;
	int i,dim,num_segs;
	double  (*func)(POINTER,double*);
        POINTER func_params;
        CURVE **wall,**contact;
	char s[100];

	FT_ScalarMemoryAlloc((POINTER*)&circle_params,sizeof(CIRCLE_PARAMS));
        circle_params->dim = dim = front->rect_grid->dim;
        circle_params->add_plan_surf = NO;
        CursorAfterString(infile,"Enter the center of implosion:");
        for (i = 0; i < dim; ++i)
	{
            fscanf(infile,"%lf",&circle_params->cen[i]);
            (void) printf("%f ",circle_params->cen[i]);
	}
	(void) printf("\n");
        CursorAfterString(infile,"Enter radius of the wall:");
        fscanf(infile,"%lf",&circle_params->R);
        (void) printf("%f\n",circle_params->R);

	func_params = (POINTER)circle_params;
	func = level_circle_func;
	neg_comp = GAS_COMP1;
	pos_comp = SOLID_COMP;
	wall = (CURVE**)FT_CreateLevelHyperSurfs(front->rect_grid,front->interf,
                    	neg_comp,pos_comp,func,func_params,DIRICHLET_BOUNDARY,
			&num_segs);

	CursorAfterString(infile,
                        "Enter radius of the two fluid interface:");
        fscanf(infile,"%lf",&circle_params->R);
        (void) printf("%f\n",circle_params->R);
	CursorAfterString(infile,"Enter yes to add perturbation:");
        fscanf(infile,"%s",s);
	(void) printf("%s\n",s);
	if (s[0] == 'y' || s[0] == 'Y')
	{
	    FOURIER_POLY *fpoly;
            circle_params->add_perturbation = YES;
	    FT_ScalarMemoryAlloc((POINTER*)&circle_params->fpoly,
				sizeof(FOURIER_POLY));
	    fpoly = circle_params->fpoly;
	    CursorAfterString(infile,"Enter number of Fourier modes:");
            fscanf(infile,"%d",&fpoly->num_modes);
            (void) printf("%d\n",fpoly->num_modes);
	    FT_VectorMemoryAlloc((POINTER*)&fpoly->A,fpoly->num_modes,
				sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&fpoly->phase,fpoly->num_modes,
				sizeof(double));
	    FT_MatrixMemoryAlloc((POINTER*)&fpoly->nu,1,fpoly->num_modes,
				sizeof(double));
	    for (i = 0; i < fpoly->num_modes; ++i)
	    {
		sprintf(s,"For %d-th mode",i);
	    	CursorAfterString(infile,s);
		sprintf(s,"enter frequency, amplitude and phase:");
	    	CursorAfterString(infile,s);
		fscanf(infile,"%lf %lf %lf",&fpoly->nu[0][i],&fpoly->A[i],
					&fpoly->phase[i]);
		(void) printf("%f %f %f\n",fpoly->nu[0][i],fpoly->A[i],
					fpoly->phase[i]);
	    }
	}
        neg_comp = GAS_COMP2;
        pos_comp = GAS_COMP1;
	contact = (CURVE**)FT_CreateLevelHyperSurfs(front->rect_grid,
			front->interf,neg_comp,pos_comp,func,func_params,
			FIRST_PHYSICS_WAVE_TYPE,&num_segs);

	level_func_pack->func = NULL;
	level_func_pack->point_array = NULL;

	fclose(infile);	
}	/* end initImplosionIntfc */

void G_CARTESIAN::initMTFusionIntfc(
	LEVEL_FUNC_PACK *level_func_pack,
	char *inname)
{
	FILE *infile = fopen(inname,"r");
	static CIRCLE_PARAMS *circle_params;
	COMPONENT neg_comp,pos_comp;
	int i,dim,num_segs;
	double  (*func)(POINTER,double*);
        POINTER func_params;
        CURVE **wall,**contact;
	char s[100];

	FT_ScalarMemoryAlloc((POINTER*)&circle_params,sizeof(CIRCLE_PARAMS));
        circle_params->dim = dim = front->rect_grid->dim;
        circle_params->add_plan_surf = NO;
	circle_params->add_perturbation = NO;
        CursorAfterString(infile,"Enter the center of the sphere:");
        for (i = 0; i < dim; ++i)
	{
            fscanf(infile,"%lf",&circle_params->cen[i]);
            (void) printf("%f ",circle_params->cen[i]);
	}
	(void) printf("\n");
        CursorAfterString(infile,"Enter radius of the sphere:");
        fscanf(infile,"%lf",&circle_params->R);
        (void) printf("%f\n",circle_params->R);

	func_params = (POINTER)circle_params;
	func = level_circle_func;
	neg_comp = GAS_COMP1;
	pos_comp = SOLID_COMP;
	wall = (CURVE**)FT_CreateLevelHyperSurfs(front->rect_grid,front->interf,
                    	neg_comp,pos_comp,func,func_params,DIRICHLET_BOUNDARY,
			&num_segs);

	CursorAfterString(infile,
                        "Enter radius of the two fluid interface:");
        fscanf(infile,"%lf",&circle_params->R);
        (void) printf("%f\n",circle_params->R);
	CursorAfterString(infile,"Enter yes to add perturbation:");
        fscanf(infile,"%s",s);
	(void) printf("%s\n",s);
	if (s[0] == 'y' || s[0] == 'Y')
	{
	    FOURIER_POLY *fpoly;
            circle_params->add_perturbation = YES;
	    FT_ScalarMemoryAlloc((POINTER*)&circle_params->fpoly,
				sizeof(FOURIER_POLY));
	    fpoly = circle_params->fpoly;
	    CursorAfterString(infile,"Enter number of Fourier modes:");
            fscanf(infile,"%d",&fpoly->num_modes);
            (void) printf("%d\n",fpoly->num_modes);
	    FT_VectorMemoryAlloc((POINTER*)&fpoly->A,fpoly->num_modes,
				sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&fpoly->phase,fpoly->num_modes,
				sizeof(double));
	    FT_MatrixMemoryAlloc((POINTER*)&fpoly->nu,1,fpoly->num_modes,
				sizeof(double));
	    for (i = 0; i < fpoly->num_modes; ++i)
	    {
		sprintf(s,"For %d-th mode",i);
	    	CursorAfterString(infile,s);
		sprintf(s,"enter frequency, amplitude and phase:");
	    	CursorAfterString(infile,s);
		fscanf(infile,"%lf %lf %lf",&fpoly->nu[0][i],&fpoly->A[i],
					&fpoly->phase[i]);
		(void) printf("%f %f %f\n",fpoly->nu[0][i],fpoly->A[i],
					fpoly->phase[i]);
	    }
	}
        neg_comp = GAS_COMP2;
        pos_comp = GAS_COMP1;
	contact = (CURVE**)FT_CreateLevelHyperSurfs(front->rect_grid,
			front->interf,neg_comp,pos_comp,func,func_params,
			FIRST_PHYSICS_WAVE_TYPE,&num_segs);

	level_func_pack->func = NULL;
	level_func_pack->point_array = NULL;
	if (debugging("trace"))
	{
	    char dirname[100];
	    sprintf(dirname,"init_intfc-%d",pp_mynode());
	    if (dim == 2)
		xgraph_2d_intfc(dirname,front->interf);
	    else if (dim == 3)
		gview_plot_interface(dirname,front->interf);
	}

	fclose(infile);	
}	/* end initMTFusionIntfc */

void G_CARTESIAN::initProjectileIntfc(
	LEVEL_FUNC_PACK *level_func_pack,
	char *inname)
{
	FILE *infile = fopen(inname,"r");
	static PROJECTILE_PARAMS *proj_params;
	static RECTANGLE_PARAMS *rparams;
	COMPONENT neg_comp,pos_comp;
	int i,dim,num_segs;
	double  (*func)(POINTER,double*);
        POINTER func_params;
        CURVE **wall,**projectile;
	double gun_length,gun_thickness,gap;
	RECT_GRID *rgr = front->rect_grid;
	RG_PARAMS rgb_params;
	char string[100];

	FT_ScalarMemoryAlloc((POINTER*)&proj_params,sizeof(PROJECTILE_PARAMS));
	FT_ScalarMemoryAlloc((POINTER*)&rparams,sizeof(RECTANGLE_PARAMS));
        proj_params->dim = dim = rgr->dim;
	gap = 0.001*rgr->h[1];
        CursorAfterString(infile,"Enter the center of projectile:");
        for (i = 0; i < dim; ++i)
	{
            fscanf(infile,"%lf",&proj_params->cen[i]);
            (void) printf("%f ",proj_params->cen[i]);
	}
	(void) printf("\n");
        CursorAfterString(infile,"Enter the radius of the projectile:");
        fscanf(infile,"%lf",&proj_params->R);
        (void) printf("%f\n",proj_params->R);
        CursorAfterString(infile,"Enter the head height of the projectile:");
        fscanf(infile,"%lf",&proj_params->r);
        (void) printf("%f\n",proj_params->r);
        CursorAfterString(infile,"Enter the butt height of the projectile:");
        fscanf(infile,"%lf",&proj_params->h);
        (void) printf("%f\n",proj_params->h);
	proj_params->R -= 0.5*gap;

	func_params = (POINTER)proj_params;
	func = projectile_func;
	neg_comp = SOLID_COMP;
	pos_comp = GAS_COMP1;
	projectile = (CURVE**)FT_CreateLevelHyperSurfs(front->rect_grid,
			front->interf,neg_comp,pos_comp,func,func_params,
			MOVABLE_BODY_BOUNDARY,&num_segs);

	if (CursorAfterStringOpt(infile,"Enter yes if there is gun:"))
	{
	    fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
	    if (string[0] == 'Y' || string[0] == 'y')
	    {
            	CursorAfterString(infile,"Enter the gun open position:");
            	fscanf(infile,"%lf",&gun_length);
            	(void) printf("%f\n",gun_length);
            	CursorAfterString(infile,"Enter the gun wall thickness:");
            	fscanf(infile,"%lf",&gun_thickness);
            	(void) printf("%f\n",gun_thickness);
	    	gun_thickness -= 0.5*gap;

	    	rparams->x0 = rgr->L[0] - rgr->h[0];
	    	rparams->y0 = proj_params->cen[1] + 
			(proj_params->R + gap);
	    	rparams->a = gun_length*2.0;
	    	rparams->b = gun_thickness;
		adjust_rectangle_params(rparams,front->rect_grid);

	    	func = rectangle_func;
	    	func_params = (POINTER)rparams;

	    	neg_comp = SOLID_COMP;
	    	pos_comp = GAS_COMP1;
	    	wall = (CURVE**)FT_CreateLevelHyperSurfs(front->rect_grid,
			front->interf,
                    	neg_comp,pos_comp,func,func_params,NEUMANN_BOUNDARY,
			&num_segs);

	    	rparams->y0 = proj_params->cen[1] - 
			(proj_params->R + gap) - gun_thickness;
		adjust_rectangle_params(rparams,front->rect_grid);
	    	wall = (CURVE**)FT_CreateLevelHyperSurfs(front->rect_grid,
			front->interf,
                    	neg_comp,pos_comp,func,func_params,NEUMANN_BOUNDARY,
			&num_segs);
	    }
	}

	level_func_pack->func = NULL;
	level_func_pack->point_array = NULL;

	fclose(infile);	
}	/* end initProjectileIntfc */

void G_CARTESIAN::initImplosionStates()
{
	int i,j,k,l,index;
	EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
	double coords[MAXD];
	COMPONENT comp;
	STATE *sl,*sr,state;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
	INTERFACE *intfc = front->interf;
	double *dens = field.dens;
	double *engy = field.engy;
	double *pres = field.pres;
	double **momn = field.momn;

        next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
	    FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
	    getAmbientState(sl,eqn_params,Coords(p),negative_component(hs));
	    getAmbientState(sr,eqn_params,Coords(p),positive_component(hs));
	}
	FT_MakeGridIntfc(front);
	setDomain();

	switch (dim)
	{
	case 2:
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index2d(i,j,top_gmax);
		comp = top_comp[index];
		getRectangleCenter(index,coords);
	    	getAmbientState(&state,eqn_params,coords,comp);
		dens[index] = state.dens;
		pres[index] = state.pres;
		engy[index] = state.engy;
		for (l = 0; l < dim; ++l)
		    momn[l][index] = state.momn[l];
	    }
	    break;
	case 3:
	    for (k = imin[2]; k <= imax[2]; ++k)
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index3d(i,j,k,top_gmax);
		comp = top_comp[index];
		getRectangleCenter(index,coords);
	    	getAmbientState(&state,eqn_params,coords,comp);
		dens[index] = state.dens;
		pres[index] = state.pres;
		engy[index] = state.engy;
		for (l = 0; l < dim; ++l)
		    momn[l][index] = state.momn[l];
	    }
	    break;
	}
	scatMeshStates();
}	/* end initImplosionStates */

void G_CARTESIAN::initMTFusionStates()
{
	int i,j,k,l,index;
	EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
	double coords[MAXD];
	COMPONENT comp;
	STATE *sl,*sr,state;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
	INTERFACE *intfc = front->interf;
	double *dens = field.dens;
	double *engy = field.engy;
	double *pres = field.pres;
	double **momn = field.momn;

        next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
	    FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
	    getAmbientState(sl,eqn_params,Coords(p),negative_component(hs));
	    getAmbientState(sr,eqn_params,Coords(p),positive_component(hs));
	}
	FT_MakeGridIntfc(front);
	setDomain();

	switch (dim)
	{
	case 2:
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index2d(i,j,top_gmax);
		comp = top_comp[index];
		getRectangleCenter(index,coords);
	    	getAmbientState(&state,eqn_params,coords,comp);
		dens[index] = state.dens;
		pres[index] = state.pres;
		engy[index] = state.engy;
		for (l = 0; l < dim; ++l)
		    momn[l][index] = state.momn[l];
	    }
	    break;
	case 3:
	    for (k = imin[2]; k <= imax[2]; ++k)
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index3d(i,j,k,top_gmax);
		comp = top_comp[index];
		getRectangleCenter(index,coords);
	    	getAmbientState(&state,eqn_params,coords,comp);
		dens[index] = state.dens;
		pres[index] = state.pres;
		engy[index] = state.engy;
		for (l = 0; l < dim; ++l)
		    momn[l][index] = state.momn[l];
	    }
	    break;
	}
	scatMeshStates();
}	/* end initMTFusionStates */

static void getAmbientState(
	STATE *state,
	EQN_PARAMS *eqn_params,
	double *coords,
	COMPONENT comp)
{
	EOS_PARAMS	*eos;
	double rho1 = eqn_params->rho1;
	double rho2 = eqn_params->rho2;
	double p1 = eqn_params->p1;
	double p2 = eqn_params->p2;
	double *v1 = eqn_params->v1;
	double *v2 = eqn_params->v2;
	int i,dim;
 
	if (debugging("ambient"))
	    printf("Entering getAmbientState(), coords = %f %f\n",
				coords[0],coords[1]);
	dim = eqn_params->dim;

	/* Constant density */
	for (i = 0; i < dim; ++i)
	    state->vel[i] = state->momn[i] = 0.0;
	state->dim = dim;
	eos = &(eqn_params->eos[comp]);
	state->eos = eos;
	
	switch (comp)
	{
	case GAS_COMP1:
	    state->dens = rho1;
	    state->pres = p1;
	    for (i = 0; i < dim; ++i)
	    {
	    	state->vel[i] = v1[i];
	    	state->momn[i] = rho1*v1[i];
	    }
	    state->engy = EosEnergy(state);
	    break;
	case GAS_COMP2:
	    state->dens = rho2;
	    state->pres = p2;
	    for (i = 0; i < dim; ++i)
	    {
	    	state->vel[i] = v2[i];
	    	state->momn[i] = rho2*v2[i];
	    }
	    state->engy = EosEnergy(state);
	    break;
	case SOLID_COMP:
	    state->dens = 0.0;
	    state->pres = 0.0;
	    state->engy = 0.0;
	    break;
	default:
	    printf("ERROR: Unknown component %d in getAmbientState()!\n",
				comp);
	    clean_up(ERROR);
	}
	if (debugging("ambient_state"))
	    (void) printf("Leaving getAmbientState(): state = %d %f %f %f\n",
			comp,state->dens,state->pres,state->engy);
}	/* end getAmbientState */

void G_CARTESIAN::setProjectileParams(char *inname)
{
	int i;
	FILE *infile = fopen(inname,"r");
	EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
	double		pinf,einf,gamma;
	char str[256];

	sprintf(str, "Enter gamma, pinf, einf of ambient air:");
	CursorAfterString(infile,str);
	fscanf(infile,"%lf %lf %lf",&gamma,&pinf,&einf);
	(void) printf("%f %f %f\n",gamma,pinf,einf);
	(eqn_params->eos[GAS_COMP1]).gamma = gamma;
	(eqn_params->eos[GAS_COMP1]).pinf = pinf;
	(eqn_params->eos[GAS_COMP1]).einf = einf;
	
	CursorAfterString(infile,"Enter density and pressure of ambient air:");
	fscanf(infile,"%lf %lf",&eqn_params->rho1,&eqn_params->p1);
	(void) printf("%f %f\n",eqn_params->rho1,eqn_params->p1);
	CursorAfterString(infile,"Enter gravity:");
	for (i = 0; i < dim; ++i)
	{
	    fscanf(infile,"%lf",&eqn_params->gravity[i]);
	    (void) printf("%f ",eqn_params->gravity[i]);
	}
	(void) printf("\n");
	fclose(infile);
}	/* end setProjectileParams */

void G_CARTESIAN::initProjectileStates()
{
	int i,j,k,l,index;
	EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
	double coords[MAXD];
	COMPONENT comp;
	STATE *sl,*sr,state;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
	INTERFACE *intfc = front->interf;
	double *dens = field.dens;
	double *engy = field.engy;
	double *pres = field.pres;
	double **momn = field.momn;

        next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
	    FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
	    getAmbientState(sl,eqn_params,Coords(p),negative_component(hs));
	    getAmbientState(sr,eqn_params,Coords(p),positive_component(hs));
	}
	FT_MakeGridIntfc(front);
	setDomain();

	switch (dim)
	{
	case 2:
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index2d(i,j,top_gmax);
		comp = top_comp[index];
		getRectangleCenter(index,coords);
	    	getAmbientState(&state,eqn_params,coords,comp);
		dens[index] = state.dens;
		pres[index] = state.pres;
		engy[index] = state.engy;
		for (l = 0; l < dim; ++l)
		    momn[l][index] = state.momn[l];
	    }
	    break;
	case 3:
	    for (k = imin[2]; k <= imax[2]; ++k)
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index3d(i,j,k,top_gmax);
		comp = top_comp[index];
		getRectangleCenter(index,coords);
	    	getAmbientState(&state,eqn_params,coords,comp);
		dens[index] = state.dens;
		pres[index] = state.pres;
		engy[index] = state.engy;
		for (l = 0; l < dim; ++l)
		    momn[l][index] = state.momn[l];
	    }
	    break;
	}
	scatMeshStates();
}	/* end initProjectileStates */

void G_CARTESIAN::initRiemannProb(
	LEVEL_FUNC_PACK *level_func_pack,
	char *inname)
{
	eqn_params = (EQN_PARAMS*)front->extra1;
	FILE *infile = fopen(inname,"r");
	double x;
	char string[100];

	level_func_pack->neg_component = GAS_COMP1;
	level_func_pack->pos_component = GAS_COMP2;
	level_func_pack->func = NULL;
	level_func_pack->func_params = NULL;
	level_func_pack->point_array = NULL;

	switch (eqn_params->prob_type)
	{
	case RIEMANN_PROB:
	    level_func_pack->num_points = 1;
	    FT_MatrixMemoryAlloc((POINTER*)&level_func_pack->point_array,
			1,MAXD,sizeof(double));
	    CursorAfterString(infile,"Enter position of interface:");
	    fscanf(infile,"%lf",&x);
	    (void) printf("%f\n",x);
	    level_func_pack->point_array[0][0] = x;
	    level_func_pack->wave_type = FIRST_PHYSICS_WAVE_TYPE;
	    break;
	case ONED_BLAST:
	case ONED_SSINE:
	case ONED_ASINE:
	    break;
	default:
	    (void) printf("Unknow problem type\n");
	    clean_up(ERROR);
	}
	CursorAfterString(infile,"Type yes to track the interface:");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	if (string[0] == 'y' || string[0] == 'Y')
            eqn_params->tracked = YES;
        else
            eqn_params->tracked = NO;
	fclose(infile);
}	/* end initRiemannProb */

void G_CARTESIAN::setRiemProbParams(char *inname)
{
	FILE *infile = fopen(inname,"r");
	EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
	double		pinf,einf,gamma;
	char str[256];

	sprintf(str, "Enter gamma, pinf, einf of ambient air:");
	CursorAfterString(infile,str);
	fscanf(infile,"%lf %lf %lf",&gamma,&pinf,&einf);
	(void) printf("%f %f %f\n",gamma,pinf,einf);
	(eqn_params->eos[GAS_COMP1]).gamma = gamma;
	(eqn_params->eos[GAS_COMP1]).pinf = pinf;
	(eqn_params->eos[GAS_COMP1]).einf = einf;
	(eqn_params->eos[GAS_COMP2]).gamma = gamma;
	(eqn_params->eos[GAS_COMP2]).pinf = pinf;
	(eqn_params->eos[GAS_COMP2]).einf = einf;
	
	CursorAfterString(infile,"Enter left and right density:");
	fscanf(infile,"%lf %lf",&eqn_params->rho1,&eqn_params->rho2);
	(void) printf("%f %f\n",eqn_params->rho1,eqn_params->rho2);
	CursorAfterString(infile,"Enter left and right pressure:");
	fscanf(infile,"%lf %lf",&eqn_params->p1,&eqn_params->p2);
	(void) printf("%f %f\n",eqn_params->p1,eqn_params->p2);
	CursorAfterString(infile,"Enter left and right velocity:");
	fscanf(infile,"%lf %lf",&eqn_params->v1[0],&eqn_params->v2[0]);
	(void) printf("%f %f\n",eqn_params->v1[0],eqn_params->v2[0]);
	fclose(infile);
}	/* end setRiemProbParams */

void G_CARTESIAN::setOnedParams(char *inname)
{
	FILE *infile = fopen(inname,"r");
	EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
	double		pinf,einf,gamma;
	char str[256];

	sprintf(str, "Enter gamma, pinf, einf of ambient air:");
        CursorAfterString(infile,str);
        fscanf(infile,"%lf %lf %lf",&gamma,&pinf,&einf);
        (void) printf("%f %f %f\n",gamma,pinf,einf);
        (eqn_params->eos[GAS_COMP1]).gamma = gamma;
        (eqn_params->eos[GAS_COMP1]).pinf = pinf;
        (eqn_params->eos[GAS_COMP1]).einf = einf;
        (eqn_params->eos[GAS_COMP2]).gamma = gamma;
        (eqn_params->eos[GAS_COMP2]).pinf = pinf;
        (eqn_params->eos[GAS_COMP2]).einf = einf;
}	/* end setOnedProbParams */

void G_CARTESIAN::initRiemProbStates()
{
	int i,l,index;
	EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
	double coords[MAXD];
	COMPONENT comp;
	STATE *sl,*sr,state;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
	INTERFACE *intfc = front->interf;
	double *dens = field.dens;
	double *engy = field.engy;
	double *pres = field.pres;
	double **momn = field.momn;

        next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
	    FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
	    getAmbientState(sl,eqn_params,Coords(p),negative_component(hs));
	    getAmbientState(sr,eqn_params,Coords(p),positive_component(hs));
	}
	FT_MakeGridIntfc(front);
	setDomain();

	switch (dim)
	{
	case 1:
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index1d(i,top_gmax);
		comp = top_comp[index];
		getRectangleCenter(index,coords);
	    	getAmbientState(&state,eqn_params,coords,comp);
		dens[index] = state.dens;
		pres[index] = state.pres;
		engy[index] = state.engy;
		for (l = 0; l < dim; ++l)
		    momn[l][index] = state.momn[l];
	    }
	    break;
	case 2:
	case 3:
	default:
	    (void) printf("initRiemProbStates() not for case dim = %d\n",dim);
	    clean_up(ERROR);
	}
	scatMeshStates();
}	/* end initRiemProbStates */

void G_CARTESIAN::initBlastWaveStates()
{
	int i,l,index;
	EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
	double coords[MAXD];
	COMPONENT comp;
	STATE *sl,*sr,state;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
	INTERFACE *intfc = front->interf;
	double *dens = field.dens;
	double *engy = field.engy;
	double *pres = field.pres;
	double **momn = field.momn;

        next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
	    FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
	    getBlastState(sl,eqn_params,Coords(p),negative_component(hs));
	    getBlastState(sr,eqn_params,Coords(p),positive_component(hs));
	}
	FT_MakeGridIntfc(front);
	setDomain();

	switch (dim)
	{
	case 1:
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index1d(i,top_gmax);
		comp = top_comp[index];
		getRectangleCenter(index,coords);
	    	getBlastState(&state,eqn_params,coords,comp);
		dens[index] = state.dens;
		pres[index] = state.pres;
		engy[index] = state.engy;
		for (l = 0; l < dim; ++l)
		    momn[l][index] = state.momn[l];
	    }
	    break;
	case 2:
	case 3:
	default:
	    (void) printf("initBlastWaveStates() not for case dim = %d\n",dim);
	    clean_up(ERROR);
	}
	scatMeshStates();
}	/* end initBlastWaveStates */

void G_CARTESIAN::initShockSineWaveStates()
{
	int i,l,index;
	EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
	double coords[MAXD];
	COMPONENT comp;
	STATE *sl,*sr,state;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
	INTERFACE *intfc = front->interf;
	double *dens = field.dens;
	double *engy = field.engy;
	double *pres = field.pres;
	double **momn = field.momn;

        next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
	    FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
	    getShockSineWaveState(sl,eqn_params,Coords(p),
				negative_component(hs));
	    getShockSineWaveState(sr,eqn_params,Coords(p),
				positive_component(hs));
	}
	FT_MakeGridIntfc(front);
	setDomain();

	switch (dim)
	{
	case 1:
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index1d(i,top_gmax);
		comp = top_comp[index];
		getRectangleCenter(index,coords);
	    	getShockSineWaveState(&state,eqn_params,coords,comp);
		dens[index] = state.dens;
		pres[index] = state.pres;
		engy[index] = state.engy;
		for (l = 0; l < dim; ++l)
		    momn[l][index] = state.momn[l];
	    }
	    break;
	case 2:
	case 3:
	default:
	    (void) printf("initShockSineWaveStates() not for case dim = %d\n",
				dim);
	    clean_up(ERROR);
	}
	scatMeshStates();
}	/* end initShockSineWaveStates */

void G_CARTESIAN::initAccuracySineWaveStates()
{
	int i,l,index;
	EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
	double coords[MAXD];
	COMPONENT comp;
	STATE *sl,*sr,state;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
	INTERFACE *intfc = front->interf;
	double *dens = field.dens;
	double *engy = field.engy;
	double *pres = field.pres;
	double **momn = field.momn;

        next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
	    FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
	    getAccuracySineWaveState(sl,eqn_params,Coords(p),
				negative_component(hs));
	    getAccuracySineWaveState(sr,eqn_params,Coords(p),
				positive_component(hs));
	}
	FT_MakeGridIntfc(front);
	setDomain();

	switch (dim)
	{
	case 1:
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index1d(i,top_gmax);
		comp = top_comp[index];
		getRectangleCenter(index,coords);
	    	getAccuracySineWaveState(&state,eqn_params,coords,comp);
		dens[index] = state.dens;
		pres[index] = state.pres;
		engy[index] = state.engy;
		for (l = 0; l < dim; ++l)
		    momn[l][index] = state.momn[l];
	    }
	    break;
	case 2:
	case 3:
	default:
	    (void) printf("initAcuracySineWaveStates() not for case dim = %d\n",
				dim);
	    clean_up(ERROR);
	}
	scatMeshStates();
}	/* end initAccuracySineWaveStates */

static void getBlastState(
	STATE *state,
	EQN_PARAMS *eqn_params,
	double *coords,
	COMPONENT comp)
{
	EOS_PARAMS	*eos;
	eos = &(eqn_params->eos[comp]);
	state->eos = eos;
	state->dim = 1;
	if (coords[0] < 0.1)
	{
	    state->dens = 1.0;
	    state->vel[0] = 0.0;
	    state->pres = 1000.0;
	}
	else if (coords[0] > 0.9)
	{
	    state->dens = 1.0;
	    state->vel[0] = 0.0;
	    state->pres = 100.0;
	}
	else
	{
	    state->dens = 1.0;
	    state->vel[0] = 0.0;
	    state->pres = 0.01;
	}
	state->momn[0] = state->vel[0] * state->dens;
	state->engy = EosInternalEnergy(state);
}	/* end getBlastState */

static void getShockSineWaveState(
	STATE *state,
	EQN_PARAMS *eqn_params,
	double *coords,
	COMPONENT comp)
{
	EOS_PARAMS	*eos;
	eos = &(eqn_params->eos[comp]);
	state->eos = eos;
	state->dim = 1;
	if (coords[0] < -4.0)
	{
	    state->dens = 3.857143;
	    state->momn[0] = 2.629329;
	    state->pres = 10.333333;
	}
	else
	{
	    state->dens = 1.0 + 0.2*sin(5.0*coords[0]);
	    state->momn[0] = 0.0;
	    state->pres = 1.0;
	}
	state->vel[0] = state->momn[0]/state->dens;
	state->engy = EosEnergy(state);
}	/* end getShockSineWaveState */

static void getAccuracySineWaveState(
	STATE *state,
	EQN_PARAMS *eqn_params,
	double *coords,
	COMPONENT comp)
{
	EOS_PARAMS	*eos;
	eos = &(eqn_params->eos[comp]);
	state->eos = eos;
	state->dim = 1;
	state->dens = 1.0 + 0.2 * sin(PI*coords[0]);
	state->vel[0] = 0.7;
	state->pres = 1.0;
	state->momn[0]= state->dens * state->vel[0];
	state->engy = EosEnergy(state);
}	/* end getAccuracySineWaveState */

void G_CARTESIAN::initRectPlaneIntfc(
        LEVEL_FUNC_PACK *level_func_pack,
        char *inname)
{
        EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
        FILE *infile = fopen(inname,"r");
        static RECT_BOX_PARAMS *rect_params;
        int i,dim;
        PROB_TYPE prob_type = eqn_params->prob_type;

        FT_ScalarMemoryAlloc((POINTER*)&rect_params,sizeof(RECT_BOX_PARAMS));
        rect_params->dim = dim = front->rect_grid->dim;
        CursorAfterString(infile,"Enter the center of the rectangle:");
        for (i = 0; i < dim; ++i)
        {
            fscanf(infile,"%lf",&rect_params->center[i]);
            (void) printf("%f ",rect_params->center[i]);
        }
        (void) printf("\n");
        CursorAfterString(infile,"Enter lengths of the rectangle:");
        for (i = 0; i < dim; ++i)
        {
            fscanf(infile,"%lf",&rect_params->length[i]);
            (void) printf("%f\n",rect_params->length[i]);
        }
        (void) printf("\n");

        level_func_pack->func_params = (POINTER)rect_params;

        switch (prob_type)
        {
        case FLUID_SOLID_RECT:
            level_func_pack->neg_component = SOLID_COMP;
            level_func_pack->pos_component = GAS_COMP1;
            level_func_pack->func = rect_box_func;
            level_func_pack->wave_type = MOVABLE_BODY_BOUNDARY;
            break;
        default:
            (void) printf("ERROR: entering wrong initialization function\n");
            clean_up(ERROR);
        }
        fclose(infile);
}       /* end initRectPlaneIntfc */

void G_CARTESIAN::initTrianglePlaneIntfc(
        LEVEL_FUNC_PACK *level_func_pack,
        char *inname)
{
        EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
        FILE *infile = fopen(inname,"r");
        static TRIANGLE_PARAMS *tri_params;
        int i,dim;
        char msg[100];
        PROB_TYPE prob_type = eqn_params->prob_type;

        FT_ScalarMemoryAlloc((POINTER*)&tri_params,sizeof(TRIANGLE_PARAMS));

        CursorAfterString(infile,"Triangle is specified by three vertices");
        (void) printf("\n");
        for (i = 0; i < 3; ++i)
        {
            sprintf(msg,"Enter coordinates of point %d:",i+1);
            CursorAfterString(infile,msg);
            fscanf(infile,"%lf %lf",&tri_params->x[i],&tri_params->y[i]);
            (void) printf("%f %f\n",tri_params->x[i],tri_params->y[i]);
        }

        level_func_pack->func_params = (POINTER)tri_params;

        switch (prob_type)
        {
        case FLUID_SOLID_TRIANGLE:
            level_func_pack->neg_component = SOLID_COMP;
            level_func_pack->pos_component = GAS_COMP1;
            level_func_pack->func = triangle_func;
            //level_func_pack->wave_type = NEUMANN_BOUNDARY;
            level_func_pack->wave_type = MOVABLE_BODY_BOUNDARY;
            break;
        default:
            (void) printf("ERROR: entering wrong initialization function\n");
            clean_up(ERROR);
        }
        fclose(infile);
}       /* end initTrianglePlaneIntfc */

void G_CARTESIAN::initCylinderPlaneIntfc(
        LEVEL_FUNC_PACK *level_func_pack,
        char *inname)
{
        EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
        FILE *infile = fopen(inname,"r");
        static CYLINDER_PARAMS *cylinder_params;
        int i;
        PROB_TYPE prob_type = eqn_params->prob_type;

        FT_ScalarMemoryAlloc((POINTER*)&cylinder_params,sizeof(CYLINDER_PARAMS));
        CursorAfterString(infile,"Enter the center of the cylinder:");
        for (i = 0; i < 3; ++i)
        {
            fscanf(infile,"%lf",&cylinder_params->center[i]);
            (void) printf("%f ",cylinder_params->center[i]);
        }
        (void) printf("\n");
        CursorAfterString(infile,"Enter radius of the cylinder:");
        fscanf(infile,"%lf",&cylinder_params->radius);
        (void) printf("%f\n",cylinder_params->radius);
        CursorAfterString(infile,"Enter height of the cylinder:");
        fscanf(infile,"%lf",&cylinder_params->height);
        (void) printf("%f\n",cylinder_params->height);

        level_func_pack->func_params = (POINTER)cylinder_params;

        switch (prob_type)
        {
        case FLUID_SOLID_CYLINDER:
            level_func_pack->neg_component = SOLID_COMP;
            level_func_pack->pos_component = GAS_COMP1;
            level_func_pack->func = cylinder_func;
            //level_func_pack->wave_type = NEUMANN_BOUNDARY;
            level_func_pack->wave_type = MOVABLE_BODY_BOUNDARY;
            break;
        default:
            (void) printf("ERROR: entering wrong initialization function\n");
            clean_up(ERROR);
        }
        fclose(infile);
}       /* end initCylinderPlaneIntfc */

extern  void prompt_for_rigid_body_params(
        int dim,
        char *inname,
        RG_PARAMS *rgb_params)
{
        int i;
        char msg[100],s[100],ss[100];
        FILE *infile = fopen(inname,"r");
        boolean is_preset_motion = NO;
        double mag_dir;

        if (debugging("rgbody"))
            (void) printf("Enter prompt_for_rigid_body_params()\n");

        rgb_params->dim = dim;
        CursorAfterString(infile,"Type yes if motion is preset: ");
        fscanf(infile,"%s",s);
        (void) printf("%s\n",s);
        if (s[0] == 'y' || s[0] == 'Y')
        {
            (void) printf("Available preset motion types are:\n");
            (void) printf("\tPRESET_TRANSLATION\n");
            (void) printf("\tPRESET_ROTATION\n");
            (void) printf("\tPRESET_MOTION (general)\n");
            CursorAfterString(infile,"Enter type of preset motion: ");
            fscanf(infile,"%s",s);
            (void) printf("%s\n",s);
            switch(s[7])
            {
            case 'M':
                rgb_params->motion_type = PRESET_MOTION;
                break;
            case 'C':
                rgb_params->motion_type = PRESET_COM_MOTION;
                break;
            case 'T':
                rgb_params->motion_type = PRESET_TRANSLATION;
                break;
            case 'R':
                rgb_params->motion_type = PRESET_ROTATION;
                break;
            default:
                (void) printf("Unknow type of preset motion!\n");
                clean_up(ERROR);
            }
        }
        else
        {
            (void) printf("Available dynamic motion types are:\n");
            (void) printf("\tFREE_MOTION:\n");
            (void) printf("\tCOM_MOTION (center of mass):\n");
            (void) printf("\tTRANSLATION:\n");
            (void) printf("\tROTATION:\n");
            CursorAfterString(infile,"Enter type of dynamic motion: ");
            fscanf(infile,"%s",s);
            (void) printf("%s\n",s);
            switch(s[0])
            {
            case 'F':
                rgb_params->motion_type = FREE_MOTION;
                break;
            case 'C':
                rgb_params->motion_type = COM_MOTION;
                break;
            case 'T':
                rgb_params->motion_type = TRANSLATION;
                break;
            case 'R':
                rgb_params->motion_type = ROTATION;
                break;
            default:
                (void) printf("Unknow type of motion!\n");
                clean_up(ERROR);
            }
        }
        if (rgb_params->motion_type == TRANSLATION ||
            rgb_params->motion_type == PRESET_TRANSLATION)
        {
            mag_dir = 0.0;
            CursorAfterString(infile,"Enter the direction of motion:");
            for (i = 0; i < dim; ++i)
            {
                fscanf(infile,"%lf",&rgb_params->translation_dir[i]);
                (void) printf("%f ",rgb_params->translation_dir[i]);
                mag_dir += sqr(rgb_params->translation_dir[i]);
            }
            (void) printf("\n");
            mag_dir = sqrt(mag_dir);
            for (i = 0; i < dim; ++i)
                rgb_params->translation_dir[i] /= mag_dir;
        }
        if (rgb_params->motion_type == FREE_MOTION ||
            rgb_params->motion_type == COM_MOTION ||
            rgb_params->motion_type == TRANSLATION)
        {
            sprintf(msg,"Enter the total mass for rigid body:");
            CursorAfterString(infile,msg);
            fscanf(infile,"%lf",&rgb_params->total_mass);
            (void) printf("%f\n",rgb_params->total_mass);
        }
        if (rgb_params->motion_type == FREE_MOTION ||
            rgb_params->motion_type == COM_MOTION ||
            rgb_params->motion_type == TRANSLATION ||
            rgb_params->motion_type == PRESET_MOTION ||
            rgb_params->motion_type == PRESET_TRANSLATION)
        {
            sprintf(msg,"Enter the initial center of mass for rigid body:");
            CursorAfterString(infile,msg);
            for (i = 0; i < dim; ++i)
            {
                fscanf(infile,"%lf",&rgb_params->center_of_mass[i]);
                (void) printf("%f ",rgb_params->center_of_mass[i]);
            }
            (void) printf("\n");
            sprintf(msg,"Enter the initial center of mass velocity:");
            CursorAfterString(infile,msg);
            for (i = 0; i < dim; ++i)
            {
                fscanf(infile,"%lf",&rgb_params->cen_of_mass_velo[i]);
                (void) printf("%f ",rgb_params->cen_of_mass_velo[i]);
            }
            (void) printf("\n");
        }
        if (rgb_params->motion_type == PRESET_ROTATION)
        {
            /* 2D preset rotation is always about the z-axis */
            /* 3D preset rotation axis along rotation_dir */
            if (dim == 3)
            {
                mag_dir = 0.0;
                CursorAfterString(infile,"Enter the direction of rotation:");
                for (i = 0; i < dim; ++i)
                {
                    fscanf(infile,"%lf",&rgb_params->rotation_dir[i]);
                    (void) printf("%f ",rgb_params->rotation_dir[i]);
                    mag_dir += sqr(rgb_params->rotation_dir[i]);
                }
                (void) printf("\n");
                mag_dir = sqrt(mag_dir);
                for (i = 0; i < dim; ++i)
                    rgb_params->rotation_dir[i] /= mag_dir;
                /* initialize the euler parameters */
                rgb_params->euler_params[0] = 1.0;
                for (i = 1; i < 4; ++i)
                    rgb_params->euler_params[i] = 0.0;
            }
            /* Center of axis is the coordinate of a point on the axis */
            CursorAfterString(infile,"Enter rotation center:");
            for (i = 0; i < dim; ++i)
            {
                fscanf(infile,"%lf",&rgb_params->rotation_cen[i]);
                (void) printf("%f ",rgb_params->rotation_cen[i]);
            }
            (void) printf("\n");
            CursorAfterString(infile,"Enter preset angular velocity:");
            fscanf(infile,"%lf",&rgb_params->angular_velo);
            (void) printf("%f\n",rgb_params->angular_velo);
            if (dim == 3)
            {
                /* used to update the maximum speed in 3D cases */
                for (i = 0; i < dim; ++i)
                    rgb_params->p_angular_velo[i] = rgb_params->angular_velo
                                        * rgb_params->rotation_dir[i];
            }
        }
        if (rgb_params->motion_type == ROTATION)
        {
            if (CursorAfterStringOpt(infile,
                "Type yes if rigid body will rotate about an point:"))
            {
                fscanf(infile,"%s",s);
                (void) printf("%s\n",s);
                if (s[0] == 'y' || s[0] == 'Y')
                {
                    sprintf(msg,"Enter rotation center:");
                    CursorAfterString(infile,msg);
                    for (i = 0; i < dim; ++i)
                    {
                        fscanf(infile,"%lf",&rgb_params->rotation_cen[i]);
                        (void) printf("%f ",rgb_params->rotation_cen[i]);
                    }
                    (void) printf("\n");
                }
            }
            if (CursorAfterStringOpt(infile,
                "Type yes if rigid body will rotate about an axis:"))
            {
                fscanf(infile,"%s",s);
                (void) printf("%s\n",s);
                if (s[0] == 'y' || s[0] == 'Y')
                {
                    /* For 2D, it is always about the z-axis */
                    if (dim == 3)
                    {
                        sprintf(msg,"Enter direction of the axis:");
                        CursorAfterString(infile,msg);
                        for (i = 0; i < dim; ++i)
                        {
                            fscanf(infile,"%lf",&rgb_params->rotation_dir[i]);
                            (void) printf("%f ",rgb_params->rotation_dir[i]);
                            mag_dir += sqr(rgb_params->rotation_dir[i]);
                        }
                        mag_dir = sqrt(mag_dir);
                        for (i = 0; i < dim; ++i)
                            rgb_params->rotation_dir[i] /= mag_dir;
                        (void) printf("\n");
                    }
                }
            }
        }
        if (rgb_params->motion_type == FREE_MOTION ||
            rgb_params->motion_type == ROTATION)
        {
            CursorAfterString(infile,"Enter the moment of inertial: ");
            if (dim == 2)
            {
                fscanf(infile,"%lf",&rgb_params->moment_of_inertial);
                (void) printf("%f\n",rgb_params->moment_of_inertial);
            }
            else if (dim == 3)
            {
                for (i = 0; i < dim; ++i)
                {
                    fscanf(infile,"%lf",&rgb_params->p_moment_of_inertial[i]);
                    (void) printf("%f ",rgb_params->p_moment_of_inertial[i]);
                }
                (void) printf("\n");
            }
            CursorAfterString(infile,"Enter initial angular velocity: ");
            if (dim == 2)
            {
                fscanf(infile,"%lf",&rgb_params->angular_velo);
                (void) printf("%f\n",rgb_params->angular_velo);
            }
            else if (dim == 3)
            {
                for (i = 0; i < dim; ++i)
                {
                    fscanf(infile,"%lf",&rgb_params->p_angular_velo[i]);
                    (void) printf("%f ",rgb_params->p_angular_velo[i]);
                }
                (void) printf("\n");
                /* initialize the euler parameters */
                rgb_params->euler_params[0] = 1.0;
                for (i = 1; i < 4; ++i)
                    rgb_params->euler_params[i] = 0.0;
            }
        }

        if (debugging("rgbody"))
            (void) printf("Leaving prompt_for_rigid_body_params()\n");
}       /* end prompt_for_rigid_body_params */

extern void set_rgbody_params(
        RG_PARAMS rg_params,
        HYPER_SURF *hs)
{
        int i,dim = rg_params.dim;
        total_mass(hs) = rg_params.total_mass;
        mom_inertial(hs) = rg_params.moment_of_inertial;
        angular_velo(hs) = rg_params.angular_velo;
        motion_type(hs) = rg_params.motion_type;
        surface_tension(hs) = 0.0;
        for (i = 0; i < dim; ++i)
        {
            center_of_mass(hs)[i] = rg_params.center_of_mass[i];
            center_of_mass_velo(hs)[i] =
                                rg_params.cen_of_mass_velo[i];
            rotation_center(hs)[i] =
                                rg_params.rotation_cen[i];
            translation_dir(hs)[i] = rg_params.translation_dir[i];
            if (dim == 3)
            {
                rotation_direction(hs)[i] = rg_params.rotation_dir[i];
                p_mom_inertial(hs)[i] = rg_params.p_moment_of_inertial[i];
                p_angular_velo(hs)[i] = rg_params.p_angular_velo[i];
            }
        }
        if (dim == 3)
        {
            for (i = 0; i < 4; i++)
                euler_params(hs)[i] = rg_params.euler_params[i];
        }
}       /* end set_rgbody_params */

static double getStationaryVelocity(
	EQN_PARAMS *eqn_params)
{
	RIEMANN_INPUT input;
	RIEMANN_SOLN riem_soln;
	int shock_dir = eqn_params->shock_dir;
	STATE state;
	double shock_speed,contact_speed;
	double Mach_number = eqn_params->Mach_number;
	int dim = eqn_params->dim;

	input.left_state.d = eqn_params->rho1;
	input.right_state.d = eqn_params->rho2;
	input.left_state.gamma = eqn_params->eos[GAS_COMP1].gamma;
	input.right_state.gamma = eqn_params->eos[GAS_COMP2].gamma;

	input.left_state.p = input.right_state.p = eqn_params->p0;
	input.left_state.u = input.right_state.u = 0.0;
	state.dim = dim;
	if (shock_dir == 1)
	{
	    state.dens = eqn_params->rho1;
	    state.pres = eqn_params->p0;
	    state.eos = &(eqn_params->eos[GAS_COMP1]);
	    state.engy = EosInternalEnergy(&state);
	}
	else if (shock_dir == -1)
	{
	    state.dens = eqn_params->rho2;
	    state.pres = eqn_params->p0;
	    state.eos = &(eqn_params->eos[GAS_COMP2]);
	    state.engy = EosInternalEnergy(&state);
	}
	behind_state(SHOCK_MACH_NUMBER,Mach_number,&shock_speed,					shock_dir,&state,&state);

	if (shock_dir == 1)
        {
	    input.left_state.d = state.dens;
	    input.left_state.p = state.pres;
	    input.left_state.u = state.vel[dim-1];
	}
        else if (shock_dir == -1)
        {
	    input.right_state.d = state.dens;
	    input.right_state.p = state.pres;
	    input.right_state.u = state.vel[dim-1];
	}
	RiemannSolution(input,&riem_soln);
	contact_speed = riem_soln.contact.speed_contact;
	return contact_speed;
}	/* end getStationaryVelocity */

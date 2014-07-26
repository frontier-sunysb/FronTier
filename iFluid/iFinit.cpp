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

#include <iFluid.h>
#include <ifluid_basic.h>

static void TG_state(COMPONENT,double*,IF_FIELD*,int,int,IF_PARAMS*);
static void kh_state(COMPONENT,double*,IF_FIELD*,int,int,IF_PARAMS*);
static void zero_state(COMPONENT,double*,IF_FIELD*,int,int,IF_PARAMS*);
static void ambient_state(COMPONENT,double*,IF_FIELD*,int,int,IF_PARAMS*);
static void random_state(COMPONENT,double*,IF_FIELD*,int,int,IF_PARAMS *);
static void initRayleiTaylorIntfc(Front*,LEVEL_FUNC_PACK*,char*);
static void initKHIntfc(Front*,LEVEL_FUNC_PACK*,char*);
static void initCirclePlaneIntfc(Front*,LEVEL_FUNC_PACK*,char*,IF_PROB_TYPE);
static void initRectPlaneIntfc(Front*,LEVEL_FUNC_PACK*,char*,IF_PROB_TYPE);
static void initTrianglePlaneIntfc(Front*,LEVEL_FUNC_PACK*,char*,IF_PROB_TYPE);
static void initChannelFlow(Front*,LEVEL_FUNC_PACK*,char*);
static void initCylinderPlaneIntfc(Front*,LEVEL_FUNC_PACK*,char*,IF_PROB_TYPE);

extern void setInitialIntfc(
	Front *front,
        LEVEL_FUNC_PACK *level_func_pack,
        char *inname,
	IF_PROB_TYPE prob_type)
{
        IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	iFparams->m_comp1 = LIQUID_COMP1;
	iFparams->m_comp2 = LIQUID_COMP2;
        switch (prob_type)
        {
	case TWO_FLUID_BUBBLE:
        case BUBBLE_SURFACE:
            initCirclePlaneIntfc(front,level_func_pack,inname,prob_type);
            break;
        case FLUID_SOLID_CIRCLE:
	    iFparams->m_comp1 = SOLID_COMP;
            initCirclePlaneIntfc(front,level_func_pack,inname,prob_type);
            break;
        case FLUID_SOLID_RECT:
	    iFparams->m_comp1 = SOLID_COMP;
            initRectPlaneIntfc(front,level_func_pack,inname,prob_type);
            break;
        case FLUID_SOLID_TRIANGLE:
	    iFparams->m_comp1 = SOLID_COMP;
            initTrianglePlaneIntfc(front,level_func_pack,inname,prob_type);
            break;
        case TWO_FLUID_RT:
            initRayleiTaylorIntfc(front,level_func_pack,inname);
            break;
        case TWO_FLUID_KH:
            initKHIntfc(front,level_func_pack,inname);
            break;
	case TAYLOR_GREEN_VORTEX:
        case CHANNEL_FLOW:
	    iFparams->m_comp1 = SOLID_COMP;
            initChannelFlow(front,level_func_pack,inname);
            break;
        case FLUID_SOLID_CYLINDER:
            iFparams->m_comp1 = SOLID_COMP;
            initCylinderPlaneIntfc(front,level_func_pack,inname,prob_type);
            break;
	default:
	    (void) printf("In setInitialIntfc unknown type: %d\n",prob_type);
        }
}       /* end setInitialIntfc */

static void initRayleiTaylorIntfc(
	Front *front,
	LEVEL_FUNC_PACK *level_func_pack,
	char *inname)
{
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	static FOURIER_POLY level_func_params;
	static double L[MAXD],U[MAXD];
	FILE *infile = fopen(inname,"r");
	int i,j,dim,num_modes;
	char mesg[100];

	dim = level_func_params.dim = front->rect_grid->dim;
	level_func_params.L = L;
	level_func_params.U = U;
	for (i = 0; i < dim; ++i)
	{
	    level_func_params.L[i] = front->rect_grid->L[i];
	    level_func_params.U[i] = front->rect_grid->U[i];
	}

	level_func_pack->neg_component = LIQUID_COMP1;
        level_func_pack->pos_component = LIQUID_COMP2;
	level_func_pack->wave_type = FIRST_PHYSICS_WAVE_TYPE;
	CursorAfterString(infile,"Enter mean position of fluid interface:");
	fscanf(infile,"%lf",&level_func_params.z0);
	(void) printf("%f\n",level_func_params.z0);
	CursorAfterString(infile,"Enter number of sine modes:");
	fscanf(infile,"%d",&num_modes);
	(void) printf("%d\n",num_modes);
	level_func_params.num_modes = num_modes;
	FT_MatrixMemoryAlloc((POINTER*)&level_func_params.nu,num_modes,
			dim-1,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&level_func_params.phase,num_modes,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&level_func_params.A,num_modes,sizeof(double));
	for (i = 0; i < num_modes; ++i)
	{
	    sprintf(mesg,"Enter frequency of mode %d:",i+1);
	    CursorAfterString(infile,mesg);
	    for (j = 0; j < dim; ++j)
	    {
	    	fscanf(infile,"%lf",&level_func_params.nu[i][j]);
		(void) printf("%f ",level_func_params.nu[i][j]);
	    }
	    (void) printf("\n");
	    sprintf(mesg,"Enter amplitude of mode %d:",i+1);
	    CursorAfterString(infile,mesg);
	    fscanf(infile,"%lf",&level_func_params.A[i]);
	    (void) printf("%f\n",level_func_params.A[i]);
	    sprintf(mesg,"Enter phase of mode %d:",i+1);
	    CursorAfterString(infile,mesg);
	    fscanf(infile,"%lf",&level_func_params.phase[i]);
	    (void) printf("%f\n",level_func_params.phase[i]);
	}
	CursorAfterString(infile,"Enter density and viscosity of fluid 1:");
        fscanf(infile,"%lf %lf",&iFparams->rho1,&iFparams->mu1);
	(void) printf("%f %f\n",iFparams->rho1,iFparams->mu1);
        CursorAfterString(infile,"Enter density and viscosity of fluid 2:");
        fscanf(infile,"%lf %lf",&iFparams->rho2,&iFparams->mu2);
	(void) printf("%f %f\n",iFparams->rho2,iFparams->mu2);
	CursorAfterString(infile,"Enter gravity:");
        for (i = 0; i < dim; ++i)
	{
            fscanf(infile,"%lf",&iFparams->gravity[i]);
	    (void) printf("%f ",iFparams->gravity[i]);
	}
	(void) printf("\n");
	CursorAfterString(infile,"Enter surface tension:");
        fscanf(infile,"%lf",&iFparams->surf_tension);
	(void) printf("%f\n",iFparams->surf_tension);
	CursorAfterString(infile,"Enter factor of smoothing radius:");
        fscanf(infile,"%lf",&iFparams->smoothing_radius);
	(void) printf("%f\n",iFparams->smoothing_radius);

	level_func_pack->func_params = (POINTER)&level_func_params;
	level_func_pack->func = level_wave_func;
	fclose(infile);
}	/* end initRayleiTaylor */

static void initKHIntfc(
	Front *front,
	LEVEL_FUNC_PACK *level_func_pack,
	char *inname)
{
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	static FOURIER_POLY level_func_params;
	static double L[MAXD],U[MAXD];
	FILE *infile = fopen(inname,"r");
	int i,j,dim,num_modes;
	char mesg[100];

	dim = level_func_params.dim = front->rect_grid->dim;
	level_func_params.L = L;
	level_func_params.U = U;
	for (i = 0; i < dim; ++i)
	{
	    level_func_params.L[i] = front->rect_grid->L[i];
	    level_func_params.U[i] = front->rect_grid->U[i];
	}

	level_func_pack->neg_component = LIQUID_COMP1;
        level_func_pack->pos_component = LIQUID_COMP2;
	level_func_pack->wave_type = FIRST_PHYSICS_WAVE_TYPE;
	CursorAfterString(infile,"Enter mean position of fluid interface:");
	fscanf(infile,"%lf",&level_func_params.z0);
	CursorAfterString(infile,"Enter number of sine modes:");
	fscanf(infile,"%d",&num_modes);
	level_func_params.num_modes = num_modes;
	FT_MatrixMemoryAlloc((POINTER*)&level_func_params.nu,num_modes,
			dim-1,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&level_func_params.phase,num_modes,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&level_func_params.A,num_modes,sizeof(double));
	for (i = 0; i < num_modes; ++i)
	{
	    sprintf(mesg,"Enter frequency of mode %d:",i+1);
	    CursorAfterString(infile,mesg);
	    for (j = 0; j < dim; ++j)
	    	fscanf(infile,"%lf",&level_func_params.nu[i][j]);
	    sprintf(mesg,"Enter amplitude of mode %d:",i+1);
	    CursorAfterString(infile,mesg);
	    fscanf(infile,"%lf",&level_func_params.A[i]);
	    sprintf(mesg,"Enter phase of mode %d:",i+1);
	    CursorAfterString(infile,mesg);
	    fscanf(infile,"%lf",&level_func_params.phase[i]);
	}
	CursorAfterString(infile,"Enter density and viscosity of fluid 1:");
        fscanf(infile,"%lf %lf",&iFparams->rho1,&iFparams->mu1);
        CursorAfterString(infile,"Enter density and viscosity of fluid 2:");
        fscanf(infile,"%lf %lf",&iFparams->rho2,&iFparams->mu2);
        CursorAfterString(infile,"Enter bottom horizontal velocity:");
	for (j = 0; j < dim-1; ++j)
            fscanf(infile,"%lf",&iFparams->U1[j]);
        CursorAfterString(infile,"Enter top horizontal velocity:");
	for (j = 0; j < dim-1; ++j)
            fscanf(infile,"%lf",&iFparams->U2[j]);

	CursorAfterString(infile,"Enter gravity:");
        for (i = 0; i < dim; ++i)
            fscanf(infile,"%lf",&iFparams->gravity[i]);
	CursorAfterString(infile,"Enter surface tension:");
        fscanf(infile,"%lf",&iFparams->surf_tension);
	CursorAfterString(infile,"Enter factor of smoothing radius:");
        fscanf(infile,"%lf",&iFparams->smoothing_radius);

	level_func_pack->func_params = (POINTER)&level_func_params;
	level_func_pack->func = level_wave_func;
	fclose(infile);
}	/* end initRayleiTaylor */

static void initCirclePlaneIntfc(
	Front *front,
	LEVEL_FUNC_PACK *level_func_pack,
	char *inname,
	IF_PROB_TYPE prob_type)
{
        IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	FILE *infile = fopen(inname,"r");
	static CIRCLE_PARAMS *circle_params;
	int i,dim;

	iFparams = (IF_PARAMS*)front->extra1;
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

	switch (prob_type)
	{
	case TWO_FLUID_BUBBLE:
        case BUBBLE_SURFACE:
            level_func_pack->neg_component = LIQUID_COMP1;
            level_func_pack->pos_component = LIQUID_COMP2;
            level_func_pack->func = level_circle_func;
            level_func_pack->wave_type = FIRST_PHYSICS_WAVE_TYPE;
	    CursorAfterString(infile,"Enter density and viscosity of fluid 1:");
            fscanf(infile,"%lf %lf",&iFparams->rho1,&iFparams->mu1);
            (void) printf("%f %f\n",iFparams->rho1,iFparams->mu1);
            CursorAfterString(infile,"Enter density and viscosity of fluid 2:");
            fscanf(infile,"%lf %lf",&iFparams->rho2,&iFparams->mu2);
            (void) printf("%f %f\n",iFparams->rho2,iFparams->mu2);
            break;
        case FLUID_SOLID_CIRCLE:
	    iFparams->m_comp1 = SOLID_COMP;
            level_func_pack->neg_component = SOLID_COMP;
            level_func_pack->pos_component = LIQUID_COMP2;
            level_func_pack->func = level_circle_func;
            level_func_pack->wave_type = MOVABLE_BODY_BOUNDARY;
            CursorAfterString(infile,
			"Enter density and viscosity of the fluid:");
            fscanf(infile,"%lf %lf",&iFparams->rho2,&iFparams->mu2);
            (void) printf("%f %f\n",iFparams->rho2,iFparams->mu2);
            break;
	default:
	    (void) printf("ERROR: entering wrong initialization function\n");
	    clean_up(ERROR);
	}
	CursorAfterString(infile,"Enter gravity:");
        for (i = 0; i < dim; ++i)
            fscanf(infile,"%lf\n",&iFparams->gravity[i]);
	CursorAfterString(infile,"Enter surface tension:");
        fscanf(infile,"%lf",&iFparams->surf_tension);
        (void) printf("%f\n",iFparams->surf_tension);
	CursorAfterString(infile,"Enter factor of smoothing radius:");
        fscanf(infile,"%lf",&iFparams->smoothing_radius);
        (void) printf("%f\n",iFparams->smoothing_radius);
	if (prob_type == TWO_FLUID_BUBBLE &&
	    iFparams->surf_tension != 0.0)	
	{
	    double mag_g = mag_vector(iFparams->gravity,dim);
	    double rho_b = iFparams->rho1;
	    double rho_o = iFparams->rho2;
	    double mu_b = iFparams->mu1;
	    double mu_o = iFparams->mu2;
	    double sigma = iFparams->surf_tension;
	    double de = 2.0*circle_params->R;
	    double V_morton = mag_g*sqr(mu_o)*sqr(mu_o)/rho_o/sqr(sigma)/sigma;
	    double V_eotvos = rho_o*mag_g*sqr(de)/sigma;
	    printf("Bubble simulation\n");
	    printf("Morton number   = %g\n",V_morton);
	    printf("Eotvos number   = %g\n",V_eotvos);
	    printf("Density ratio   = %g\n",rho_o/rho_b);
	    printf("Viscosity ratio = %g\n",mu_o/mu_b);
	}
	
	fclose(infile);	
}	/* end initCirclePlaneIntfc */

extern void init_fluid_state_func(
	Incompress_Solver_Smooth_Basis *l_cartesian,
	IF_PROB_TYPE prob_type)
{
	switch (prob_type)
	{
	case TWO_FLUID_KH:
	    l_cartesian->getInitialState = kh_state;
	    break;
	case TAYLOR_GREEN_VORTEX:
            l_cartesian->getInitialState = TG_state;
            break;
	case RANDOM_FLOW:
	    l_cartesian->getInitialState = random_state;
	    break;
	default:
	    l_cartesian->getInitialState = ambient_state;
	}
}	/* end init_fluid_state_func */

static void TG_state(
        COMPONENT comp,
        double *coords,
        IF_FIELD *field,
        int index,
        int dim,
        IF_PARAMS *iFparams)
{
        int i;
        double tcoords[MAXD], a[MAXD] = {2*PI, 2*PI, 2*PI}, A[MAXD] = {1., -1., 0.};
        double **vel = field->vel;
        switch (dim)
        {
            case 2:
                for (i = 0; i < dim; i++)
                    tcoords[i] = coords[i];
                vel[0][index] = A[0]*sin(tcoords[0]) * cos(tcoords[1]);
                vel[1][index] = A[1]*cos(tcoords[0]) * sin(tcoords[1]);
                break;
            case 3:
                for (i = 0; i < dim; i++)
                    tcoords[i] = a[i] * coords[i];
                vel[0][index] = A[0]*cos(tcoords[0])*sin(tcoords[1])*sin(tcoords[2]);
                vel[1][index] = A[1]*sin(tcoords[0])*cos(tcoords[1])*sin(tcoords[2]);
                vel[2][index] = A[2]*sin(tcoords[0])*sin(tcoords[1])*cos(tcoords[2]);
                break;
            default:
                printf("Unknown dim = %d\n",dim);
                clean_up(ERROR);
        }

}

static void kh_state(
	COMPONENT comp,
	double *coords,
	IF_FIELD *field,
	int index,
	int dim,
	IF_PARAMS *iFparams)
{
	int i;
	double **vel = field->vel;

	for (i = 0; i < dim-1; ++i)
	{
	    vel[i][index] = (comp == LIQUID_COMP1) ? iFparams->U1[i] :
				iFparams->U2[i];
	}
	vel[dim-1][index] = 0.0;
}	/* end kh_state */

static void zero_state(
	COMPONENT comp,
	double *coords,
	IF_FIELD *field,
	int index,
	int dim,
	IF_PARAMS *iFparams)
{
	int i;
	double **vel = field->vel;
	for (i = 0; i < dim; ++i)
	    vel[i][index] = 0.0;
}	/* end zero_state */

static void random_state(
        COMPONENT comp,
        double *coords,
        IF_FIELD *field,
        int index,
        int dim,
        IF_PARAMS *iFparams)
{
        short unsigned int seed[3] = {2,72,7172};
        int i;
        double **vel = field->vel;
        for (i = 0; i < dim; ++i)
        {
            vel[i][index] = (erand48(seed) - 0.5) * 10;
        }
}       /* end random_state */

static void ambient_state(
	COMPONENT comp,
	double *coords,
	IF_FIELD *field,
	int index,
	int dim,
	IF_PARAMS *iFparams)
{
	int i;
	double *U_ambient = iFparams->U_ambient;
	double **vel = field->vel;
	for (i = 0; i < dim; ++i)
	{
	    if (ifluid_comp(comp))
	    	vel[i][index] = U_ambient[i];
	    else
	    	vel[i][index] = 0.0;
	}
}	/* end ambient_state */

extern void read_iF_prob_type(
	char *inname,
	IF_PROB_TYPE *prob_type)
{
	char string[100];
	FILE *infile = fopen(inname,"r");

	*prob_type = ERROR_TYPE;
	CursorAfterString(infile,"Enter problem type:");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	if (string[0] == 'T' || string[0] == 't')
	{
	    if (string[10] == 'B' || string[10] == 'b')
	    	*prob_type = TWO_FLUID_BUBBLE;
	    else if (string[10] == 'R' || string[10] == 'r')
	    	*prob_type = TWO_FLUID_RT;
	    else if (string[10] == 'K' || string[10] == 'k')
	    	*prob_type = TWO_FLUID_KH;
	    else if (string[7] == 'G' || string[7] == 'g')
                *prob_type = TAYLOR_GREEN_VORTEX;
	} 
	else if (string[0] == 'F' || string[0] == 'f')
	{
            if (string[6] == 'S' || string[6] == 's')
	    {
		if (string[12] == 'C' || string[12] == 'c')
                {
                    if (string[13] == 'I' || string[13] == 'i')
                        *prob_type = FLUID_SOLID_CIRCLE;
                    else if (string[13] == 'Y' || string[13] == 'y')
                        *prob_type = FLUID_SOLID_CYLINDER;
                }
		else if (string[12] == 'R' || string[12] == 'r')
                    *prob_type = FLUID_SOLID_RECT;
		else if (string[12] == 'T' || string[12] == 't')
                    *prob_type = FLUID_SOLID_TRIANGLE;
	    }
            else if (string[6] == 'R' || string[6] == 'r')
                *prob_type = FLUID_RIGID_BODY;
	}
	else if (string[0] == 'B' || string[0] == 'b')
	{
	    *prob_type = BUBBLE_SURFACE;
	}
        else if (string[0] == 'R' || string[0] == 'r')
	{
            if (string[6] == 'O' || string[6] == 'o')
                *prob_type = ROTOR_ONE_FLUID;
            else if (string[6] == 'T' || string[6] == 't')
                *prob_type = ROTOR_TWO_FLUID;
	}
	else if (string[0] == 'C' || string[0] == 'c')
	{
            *prob_type = CHANNEL_FLOW;
	}

	assert(*prob_type != ERROR_TYPE);
	fclose(infile);
}	/* end read_iFparams */

static void initChannelFlow(
	Front *front,
	LEVEL_FUNC_PACK *level_func_pack,
	char *inname)
{
        IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	FILE *infile = fopen(inname,"r");
	int i,dim = front->rect_grid->dim;

	iFparams = (IF_PARAMS*)front->extra1;

	// No interface 
        level_func_pack->neg_component = LIQUID_COMP2;
        level_func_pack->pos_component = LIQUID_COMP2;
        level_func_pack->func = NULL;
        CursorAfterString(infile,"Enter density and viscosity of the fluid:");
        fscanf(infile,"%lf %lf",&iFparams->rho2,&iFparams->mu2);
        (void) printf("%f %f\n",iFparams->rho2,iFparams->mu2);
	CursorAfterString(infile,"Enter gravity:");
        for (i = 0; i < dim; ++i)
            fscanf(infile,"%lf",&iFparams->gravity[i]);
	
	fclose(infile);	
}	/* end initChannelFlow */

static void initRectPlaneIntfc(
	Front *front,
	LEVEL_FUNC_PACK *level_func_pack,
	char *inname,
	IF_PROB_TYPE prob_type)
{
        IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	FILE *infile = fopen(inname,"r");
	static RECT_BOX_PARAMS *rect_params;
	int i,dim;

	iFparams = (IF_PARAMS*)front->extra1;
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
	    iFparams->m_comp1 = SOLID_COMP;
            level_func_pack->neg_component = SOLID_COMP;
            level_func_pack->pos_component = LIQUID_COMP2;
            level_func_pack->func = rect_box_func;
            level_func_pack->wave_type = MOVABLE_BODY_BOUNDARY;
            //level_func_pack->wave_type = NEUMANN_BOUNDARY;
            CursorAfterString(infile,
			"Enter density and viscosity of the fluid:");
            fscanf(infile,"%lf %lf",&iFparams->rho2,&iFparams->mu2);
            (void) printf("%f %f\n",iFparams->rho2,iFparams->mu2);
            break;
	default:
	    (void) printf("ERROR: entering wrong initialization function\n");
	    clean_up(ERROR);
	}
	CursorAfterString(infile,"Enter gravity:");
        for (i = 0; i < dim; ++i)
            fscanf(infile,"%lf\n",&iFparams->gravity[i]);
	CursorAfterString(infile,"Enter surface tension:");
        fscanf(infile,"%lf",&iFparams->surf_tension);
        (void) printf("%f\n",iFparams->surf_tension);
	CursorAfterString(infile,"Enter factor of smoothing radius:");
        fscanf(infile,"%lf",&iFparams->smoothing_radius);
        (void) printf("%f\n",iFparams->smoothing_radius);
	fclose(infile);	
}	/* end initRectPlaneIntfc */

static void initTrianglePlaneIntfc(
	Front *front,
	LEVEL_FUNC_PACK *level_func_pack,
	char *inname,
	IF_PROB_TYPE prob_type)
{
        IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	FILE *infile = fopen(inname,"r");
	static TRIANGLE_PARAMS *tri_params;
	int i,dim;
	char msg[100];

	iFparams = (IF_PARAMS*)front->extra1;
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
	    iFparams->m_comp1 = SOLID_COMP;
            level_func_pack->neg_component = SOLID_COMP;
            level_func_pack->pos_component = LIQUID_COMP2;
            level_func_pack->func = triangle_func;
            level_func_pack->wave_type = MOVABLE_BODY_BOUNDARY;
            CursorAfterString(infile,
			"Enter density and viscosity of the fluid:");
            fscanf(infile,"%lf %lf",&iFparams->rho2,&iFparams->mu2);
            (void) printf("%f %f\n",iFparams->rho2,iFparams->mu2);
            break;
	default:
	    (void) printf("ERROR: entering wrong initialization function\n");
	    clean_up(ERROR);
	}
	CursorAfterString(infile,"Enter gravity:");
        for (i = 0; i < dim; ++i)
            fscanf(infile,"%lf\n",&iFparams->gravity[i]);
	CursorAfterString(infile,"Enter surface tension:");
        fscanf(infile,"%lf",&iFparams->surf_tension);
        (void) printf("%f\n",iFparams->surf_tension);
	CursorAfterString(infile,"Enter factor of smoothing radius:");
        fscanf(infile,"%lf",&iFparams->smoothing_radius);
        (void) printf("%f\n",iFparams->smoothing_radius);
	fclose(infile);	
}	/* end initTrianglePlaneIntfc */

static void initCylinderPlaneIntfc(
        Front *front,
        LEVEL_FUNC_PACK *level_func_pack,
        char *inname,
        IF_PROB_TYPE prob_type)
{
        IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
        FILE *infile = fopen(inname,"r");
        static CYLINDER_PARAMS *cylinder_params;
        int i;

        iFparams = (IF_PARAMS*)front->extra1;
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
            iFparams->m_comp1 = SOLID_COMP;
            level_func_pack->neg_component = SOLID_COMP;
            level_func_pack->pos_component = LIQUID_COMP2;
            level_func_pack->func = cylinder_func;
            level_func_pack->wave_type = MOVABLE_BODY_BOUNDARY;
            //level_func_pack->wave_type = NEUMANN_BOUNDARY;
            CursorAfterString(infile,
                        "Enter density and viscosity of the fluid:");
            fscanf(infile,"%lf %lf",&iFparams->rho2,&iFparams->mu2);
            (void) printf("%f %f\n",iFparams->rho2,iFparams->mu2);
            break;
        default:
            (void) printf("ERROR: entering wrong initialization function\n");
            clean_up(ERROR);
        }
        CursorAfterString(infile,"Enter gravity:");
        for (i = 0; i < 3; ++i)
            fscanf(infile,"%lf\n",&iFparams->gravity[i]);
        CursorAfterString(infile,"Enter surface tension:");
        fscanf(infile,"%lf",&iFparams->surf_tension);
        (void) printf("%f\n",iFparams->surf_tension);
        CursorAfterString(infile,"Enter factor of smoothing radius:");
        fscanf(infile,"%lf",&iFparams->smoothing_radius);
        (void) printf("%f\n",iFparams->smoothing_radius);

        fclose(infile);
}       /* end initCylinderPlaneIntfc */

extern  void prompt_for_rigid_body_params(
        int dim,
        char *inname,
        RG_PARAMS *rgb_params)
{
        int i;
        char msg[100],s[100];
        FILE *infile = fopen(inname,"r");

        if (debugging("rgbody"))
            (void) printf("Enter prompt_for_rigid_body_params()\n");

        rgb_params->dim = dim;
        sprintf(msg,"Enter the total mass for rigid body:");
        CursorAfterString(infile,msg);
        fscanf(infile,"%lf",&rgb_params->total_mass);
        (void) printf("%f\n",rgb_params->total_mass);
        sprintf(msg,"Enter the center of mass for rigid body:");
        CursorAfterString(infile,msg);
        for (i = 0; i < dim; ++i)
        {
            fscanf(infile,"%lf",&rgb_params->center_of_mass[i]);
            (void) printf("%f ",rgb_params->center_of_mass[i]);
        }
        (void) printf("\n");
        CursorAfterString(infile,
                "Type yes if rigid body will only rotate about an axis:");
        fscanf(infile,"%s",s);
        (void) printf("%s\n",s);
        if (s[0] == 'y' || s[0] == 'Y')
        {
            if (dim == 3)
            {
                double mag_dir = 0.0;
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

            sprintf(msg,"Enter center of the axis:");
            CursorAfterString(infile,msg);
            for (i = 0; i < dim; ++i)
            {
                fscanf(infile,"%lf",&rgb_params->rotation_cen[i]);
                (void) printf("%f ",rgb_params->rotation_cen[i]);
            }
            (void) printf("\n");

            sprintf(msg,"Enter the moment of inertial about the axis:");
            CursorAfterString(infile,msg);
            fscanf(infile,"%lf",&rgb_params->moment_of_inertial);
            (void) printf("%f\n",rgb_params->moment_of_inertial);

            CursorAfterString(infile,
                        "Type yes if angular velocity is preset: ");
            fscanf(infile,"%s",s);
            (void) printf("%s\n",s);
            if (s[0] == 'y' || s[0] == 'Y')
            {
                rgb_params->motion_type = PRESET_MOTION;
                CursorAfterString(infile,"Enter preset angular velocity: ");
            }
            else
            {
                rgb_params->motion_type = ROTATION;
                CursorAfterString(infile,"Enter initial angular velocity: ");
            }
            fscanf(infile,"%lf",&rgb_params->angular_velo);
            (void) printf("%f\n",rgb_params->angular_velo);
        }
        else
        {
            sprintf(msg,"Enter the moment of inertial about center of mass:");
            CursorAfterString(infile,msg);
            fscanf(infile,"%lf",&rgb_params->moment_of_inertial);
            (void) printf("%f\n",rgb_params->moment_of_inertial);

            rgb_params->motion_type = FREE_MOTION;
            CursorAfterString(infile,
                        "Type yes if you want vertical motion only?: ");
            fscanf(infile,"%s",s);
            (void) printf("%s\n",s);
            if (s[0] == 'y' || s[0] == 'Y')
                rgb_params->motion_type = VERTICAL_MOTION;
            CursorAfterString(infile,
                        "Type yes if you want horizontal motion only?: ");
            fscanf(infile,"%s",s);
            (void) printf("%s\n",s);
            if (s[0] == 'y' || s[0] == 'Y')
                rgb_params->motion_type = HORIZONTAL_MOTION;

            sprintf(msg,"Enter the initial center of mass velocity:");
            CursorAfterString(infile,msg);
            for (i = 0; i < dim; ++i)
            {
                fscanf(infile,"%lf",&rgb_params->cen_of_mass_velo[i]);
                (void) printf("%f ",rgb_params->cen_of_mass_velo[i]);
            }
            (void) printf("\n");

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
            if (dim == 3)
                rotation_direction(hs)[i] =
                                rg_params.rotation_dir[i];
        }
}       /* end set_rgbody_params */

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
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
****************************************************************/

#include <iFluid.h>
#include "melting.h"

#define		MAX_NUM_VERTEX_IN_CELL		20

typedef struct {
        double center[3];
        double radius;
} TEST_SPHERE_PARAMS;

	/*  Local Application Function Declarations */

static void	initPhaseIntfc1d(char*,LEVEL_FUNC_PACK*,PARAMS*);
static void	initPhaseIntfc2d(char*,LEVEL_FUNC_PACK*,PARAMS*);
static void	initPhaseIntfc3d(char*,LEVEL_FUNC_PACK*,PARAMS*);
static void	read_seed_params(int,char*,SEED_PARAMS*);

static double crystal_curve(
        POINTER func_params,
        double *coords)
{

        TEST_SPHERE_PARAMS *s_params = (TEST_SPHERE_PARAMS*)func_params;
        double dist, theta;
	double *cen = s_params->center;
	double radius = s_params->radius;

        dist =   sqrt(sqr(coords[0]-cen[0]) + sqr(coords[1]-cen[1]));
        theta = asin(fabs(coords[1]-0.5)/dist);
	if (coords[0]-0.5 < 0 && coords[1]-0.5 > 0)
	    theta = PI - theta;
	else if (coords[0]-0.5 < 0 && coords[1]-0.5 < 0)
	    theta = PI + theta;
	else if (coords[0]-0.5 > 0 && coords[1]-0.5 < 0)
	    theta = 2*PI - theta;
        return dist - radius + .003*sin(6.0*theta);
}       /* end crystal_curve */

extern void readPhaseParams(
	Front *front)
{
	FILE *infile;
	char scheme[200];
	int i,num_phases;
	char string[200];
	char *in_name = InName(front);
	PARAMS *eqn_params = (PARAMS*)front->extra2;

	infile = fopen(in_name,"r");

	CursorAfterString(infile,"Enter number of phases:");
	fscanf(infile,"%d",&num_phases);
	eqn_params->num_phases = num_phases;
	(void) printf(" %d\n",num_phases);
	FT_VectorMemoryAlloc((POINTER*)&eqn_params->T0,num_phases,
					sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&eqn_params->rho,num_phases,
					sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&eqn_params->Cp,num_phases,
					sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&eqn_params->k,num_phases,
					sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&eqn_params->Ti,num_phases-1,
					sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&eqn_params->L,num_phases-1,
					sizeof(double));
	for (i = 0; i < num_phases; ++i)
	{
	    sprintf(string,"Enter ambient temperature of phase %d:",i+1);
	    CursorAfterString(infile,string);
	    fscanf(infile,"%lf",&eqn_params->T0[i]);
	    (void) printf("%f\n",eqn_params->T0[i]);
	    sprintf(string,"Enter density of phase %d:",i+1);
	    CursorAfterString(infile,string);
	    fscanf(infile,"%lf",&eqn_params->rho[i]);
	    (void) printf("%f\n",eqn_params->rho[i]);
	    sprintf(string,"Enter specific heat of phase %d:",i+1);
	    CursorAfterString(infile,string);
	    fscanf(infile,"%lf",&eqn_params->Cp[i]);
	    (void) printf("%f\n",eqn_params->Cp[i]);
	    sprintf(string,"Enter thermal conductivity of phase %d:",i+1);
	    CursorAfterString(infile,string);
	    fscanf(infile,"%lf",&eqn_params->k[i]);
	    (void) printf("%f\n",eqn_params->k[i]);
	    if (i != num_phases-1)
	    {
	        sprintf(string,"Enter melting temperature of interface %d:",
						i+1);
	        CursorAfterString(infile,string);
	    	fscanf(infile,"%lf",&eqn_params->Ti[i]);
	        (void) printf("%f\n",eqn_params->Ti[i]);
	        sprintf(string,"Enter latent heat at interface %d:",i+1);
	        CursorAfterString(infile,string);
	    	fscanf(infile,"%lf",&eqn_params->L[i]);
	        (void) printf("%f\n",eqn_params->L[i]);
	    }
	}

	eqn_params->num_scheme = UNSPLIT_IMPLICIT;  // default
	eqn_params->pde_order = 2; //default

	if (CursorAfterStringOpt(infile,"Choose numerical scheme"))
	{
	    (void) printf("\n");
	    CursorAfterString(infile,"Enter scheme:");
	    fscanf(infile,"%s",scheme);
	    (void) printf("%s\n",scheme);
	    if ((scheme[0] == 'E' || scheme[0] == 'e') &&
	    	(scheme[1] == 'X' || scheme[1] == 'x')) 
	    	eqn_params->num_scheme = UNSPLIT_EXPLICIT;
	    else if ((scheme[0] == 'E' || scheme[0] == 'e') &&
	    	(scheme[1] == 'C' || scheme[1] == 'c')) 
	    	eqn_params->num_scheme = UNSPLIT_EXPLICIT_CIM;
	    else if ((scheme[0] == 'I' || scheme[0] == 'i') &&
	    	(scheme[1] == 'M' || scheme[1] == 'm')) 
	    	eqn_params->num_scheme = UNSPLIT_IMPLICIT;
	    else if ((scheme[0] == 'C' || scheme[0] == 'c') &&
	    	(scheme[1] == 'N' || scheme[1] == 'n')) 
	    	eqn_params->num_scheme = CRANK_NICOLSON;
	    else if ((scheme[0] == 'I' || scheme[0] == 'i') &&
	    	(scheme[1] == 'C' || scheme[1] == 'c')) 
	    	eqn_params->num_scheme = UNSPLIT_IMPLICIT_CIM;
	}
	
	if (eqn_params->num_scheme == UNSPLIT_IMPLICIT)
	{
	    if (CursorAfterStringOpt(infile,"Choose order of PDE scheme:"))
	    {
	    	fscanf(infile,"%d",&eqn_params->pde_order);
		(void) printf("%d\n",eqn_params->pde_order);
	    }
	}
	eqn_params->no_fluid = YES;
	if (CursorAfterStringOpt(infile,"Enter yes to turn on fluid solver:"))
	{
	    fscanf(infile,"%s",string);
	    (void) printf("%s\n",string);
	    if (string[0] == 'y' || string[0] == 'Y')
		eqn_params->no_fluid = NO;
	}
	fclose(infile);
}

extern void initPhaseIntfc(
	char *in_name,
	int dim,
	LEVEL_FUNC_PACK *level_func_pack,
	PARAMS *eqn_params)
{
	switch (dim)
	{
	case 1:
	    return initPhaseIntfc1d(in_name,level_func_pack,eqn_params);
	case 2:
            return initPhaseIntfc2d(in_name,level_func_pack,eqn_params);
        case 3:
            return initPhaseIntfc3d(in_name,level_func_pack,eqn_params);
	}
}	/* end initPhaseIntfc */

static void initPhaseIntfc1d(
	char *in_name,
	LEVEL_FUNC_PACK *level_func_pack,
	PARAMS *eqn_params)
{
	FILE *infile;
	double **points;
	int i,num_phases = eqn_params->num_phases;
	COMPONENT neg_comp,pos_comp;
	POINT *p;
	char string[200];
	
	infile = fopen(in_name,"r");

	FT_MatrixMemoryAlloc((POINTER*)&points,num_phases-1,MAXD,
				sizeof(double));
	for (i = 0; i < num_phases-1; ++i)
	{
	    sprintf(string,"Enter position of interface %d:",i+1);
	    CursorAfterString(infile,string);
            fscanf(infile,"%lf",&points[i][0]);
            (void) printf("%f\n",points[i][0]);
	    neg_comp = (eqn_params->T0[i] < eqn_params->Ti[i]) ?
                                SOLID_COMP : LIQUID_COMP;
            pos_comp = (eqn_params->T0[i] > eqn_params->Ti[i]) ?
                                SOLID_COMP : LIQUID_COMP;
            p = make_point(points[i],neg_comp,pos_comp);
            wave_type(Hyper_surf(p)) = GROWING_BODY_BOUNDARY;
	}
}	/* end initPhaseIntfc1d */

static void initPhaseIntfc2d(
	char *in_name,
	LEVEL_FUNC_PACK *level_func_pack,
	PARAMS *eqn_params)
{
	static SEED_PARAMS seed_params;
	read_seed_params(2,in_name,&seed_params);
	if (seed_params.num_floor_seeds == 0 &&
            seed_params.num_ceiling_seeds == 0 &&
            seed_params.num_space_seeds == 0)
        {
            level_func_pack->func_params = NULL;
            level_func_pack->func = NULL;
        }
        else
        {
            level_func_pack->func_params = (POINTER)&seed_params;
            level_func_pack->func = seed_func;
        }
	level_func_pack->neg_component = SOLID_COMP;
        level_func_pack->pos_component = LIQUID_COMP;
        level_func_pack->wave_type = GROWING_BODY_BOUNDARY;

}	/* end initPhaseIntfc2d */

static void initPhaseIntfc3d(
	char *in_name,
	LEVEL_FUNC_PACK *level_func_pack,
	PARAMS *eqn_params)
{
	FILE *infile = fopen(in_name,"r");	
	char string[256];
	static SEED_PARAMS seed_params;
	static CUBOID_PARAMS cuboid_params;
	static TETRAHEDRON_PARAMS tetrahedron_params;
	
	CursorAfterString(infile,"Enter initial interface type: ");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	switch (string[0])
	{
	case 's':
	case 'S':
	    read_seed_params(3,in_name,&seed_params);
	    if (seed_params.num_floor_seeds == 0 &&
            	seed_params.num_ceiling_seeds == 0 &&
            	seed_params.num_space_seeds == 0)
      	    {
        	level_func_pack->func_params = NULL;
            	level_func_pack->func = NULL;
       	    }
            else
            {
                level_func_pack->func_params = (POINTER)&seed_params;
                level_func_pack->func = seed_func;
            }
	    level_func_pack->neg_component = SOLID_COMP;
            level_func_pack->pos_component = LIQUID_COMP;
            level_func_pack->wave_type = GROWING_BODY_BOUNDARY;
	    level_func_pack->set_3d_bdry = YES;
	    break;
	case 'c':
	case 'C':
	    CursorAfterString(infile,"Enter center of Cuboid: ");
	    fscanf(infile,"%lf %lf %lf\n",&cuboid_params.center[0],
				&cuboid_params.center[1],&cuboid_params.center[2]);
	    (void) printf("%f %f %f\n",cuboid_params.center[0],
				cuboid_params.center[1],cuboid_params.center[2]);
	    CursorAfterString(infile,"Enter edge of Cuboid: ");
	    fscanf(infile,"%lf %lf %lf\n",&cuboid_params.edge[0],
				&cuboid_params.edge[1],&cuboid_params.edge[2]);
	    (void) printf("%f %f %f\n",cuboid_params.edge[0],
				cuboid_params.edge[1],cuboid_params.edge[2]);
	    level_func_pack->func_params = (POINTER)&cuboid_params;
	    level_func_pack->func = cuboid_func;
	    level_func_pack->neg_component = SOLID_COMP;
            level_func_pack->pos_component = LIQUID_COMP;
            level_func_pack->wave_type = GROWING_BODY_BOUNDARY;
            level_func_pack->set_3d_bdry = YES;
	    break;
	case 't':
	case 'T':
	    CursorAfterString(infile,"Enter center of Tetrahedron: ");
	    fscanf(infile,"%lf %lf %lf\n",&tetrahedron_params.center[0],
				&tetrahedron_params.center[1],
				&tetrahedron_params.center[2]);
	    (void) printf("%f %f %f\n",tetrahedron_params.center[0],
				tetrahedron_params.center[1],
				tetrahedron_params.center[2]);
	    CursorAfterString(infile,"Enter edge of Tetrahedron: ");
	    fscanf(infile,"%lf\n",&tetrahedron_params.edge);
	    (void) printf("%f\n",tetrahedron_params.edge);
	    level_func_pack->func_params = (POINTER)&tetrahedron_params;
            level_func_pack->func = tetrahedron_func;
            level_func_pack->neg_component = SOLID_COMP;
            level_func_pack->pos_component = LIQUID_COMP;
            level_func_pack->wave_type = GROWING_BODY_BOUNDARY;
            level_func_pack->set_3d_bdry = YES; 
	    break;
	default:
	    (void) printf("Unknow type of initial interface!\n");
	    clean_up(ERROR);	    
	}

	fclose(infile);
}	/* end initPhaseIntfc3d */

static	void read_seed_params(
	int dim,
	char *in_name,
	SEED_PARAMS *s_params)
{
	FILE *infile;
	char s[100];
	int i,j;

	infile = fopen(in_name,"r");	
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
		    CursorAfterString(infile,"Enter temperature:");
		    fscanf(infile,"%lf",&state.temperature);
		    (void) printf("%f\n",state.temperature);
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
		    CursorAfterString(infile,"Enter temperature:");
		    fscanf(infile,"%lf",&state.temperature);
		    (void) printf("%f\n",state.temperature);
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

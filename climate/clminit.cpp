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
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
****************************************************************/

#include <iFluid.h>
#include "climate.h"

#define		MAX_NUM_VERTEX_IN_CELL		20

typedef struct {
        double center[3];
        double radius;
} TEST_SPHERE_PARAMS;

	/*  Local Application Function Declarations */

static void initRandomDrops(Front*,double**,double*,int*,int*,
                                double,double);

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
        fscanf(infile,"%d",&eqn_params->num_phases);
	num_phases = eqn_params->num_phases;
	(void) printf("%d phases are included\n",num_phases);
	FT_VectorMemoryAlloc((POINTER*)&eqn_params->T0,num_phases,
					sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&eqn_params->Ti,num_phases-1,
					sizeof(double));
	for (i = 0; i < num_phases; ++i)
	{
	    sprintf(string,"Enter ambient temperature of phase %d:",i+1);
	    CursorAfterString(infile,string);
	    fscanf(infile,"%lf",&eqn_params->T0[i]);
	    (void) printf("%f\n",eqn_params->T0[i]);
	    sprintf(string,"Enter diffusivity of phase %d:",i+1);
	    CursorAfterString(infile,string);
	    fscanf(infile,"%lf",&eqn_params->D);
	    (void) printf("%f\n",eqn_params->D);
            if (i != num_phases-1)
            {
                sprintf(string,"Enter melting temperature of interface %d:",
                                                i+1);
                CursorAfterString(infile,string);
                fscanf(infile,"%lf",&eqn_params->Ti[i]);
                (void) printf("%f\n",eqn_params->Ti[i]);
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

        eqn_params->no_droplets = YES;
        if (CursorAfterStringOpt(infile,"Enter yes to turn on droplets:"))
        {
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'y' || string[0] == 'Y')
                eqn_params->no_droplets = NO;
        }

	fclose(infile);
}

extern void initWaterDrops(
	Front *front)
{
	char string[100],msg[200];
        char *inname = InName(front);
        FILE *infile = fopen(inname,"r");
        int i,j,l,dir,num_drops;
        int *gindex;
        double **center,*radius;
        double r_bar,sigma;
        CURVE *curve;
        SURFACE *surf,*psurf;
        double *L = front->rect_grid->L;
        double *U = front->rect_grid->U;
        double T[MAXD];
        int dim = front->rect_grid->dim;
        int w_type;
        IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	PARAMS *eqn_params  = (PARAMS*)front->extra2;
        double drop_dens;

        (void) printf("Water phase state can be\n");
        (void) printf("\tIce Particle (I)\n");
        (void) printf("\tLiquid Water Drop (W)\n");
        CursorAfterString(infile,"Enter phase state of water drop:");
        fscanf(infile,"%s",string);
        (void) printf("%s\n",string);
        switch (string[0])
        {
        case 'I':
        case 'i':
            w_type = ICE_PARTICLE_BOUNDARY;
            iFparams->m_comp1 = SOLID_COMP;
            iFparams->m_comp2 = LIQUID_COMP;
            break;
        case 'W':
        case 'w':
            (void) printf("Liquid water state not yet implemented\n");
            clean_up(ERROR);
            break;
        default:
            (void) printf("Unknow phase state of water\n");
            clean_up(ERROR);
        }
        CursorAfterString(infile,"Enter number of water drops:");
        fscanf(infile,"%d",&num_drops);
        (void) printf("%d\n",num_drops);
        CursorAfterString(infile,"Enter density of water drops:");
        fscanf(infile,"%lf",&drop_dens);
        (void) printf("%f\n",drop_dens);
	eqn_params->rho_l = drop_dens;
        CursorAfterString(infile,"Enter coefficient for condensation:");
        fscanf(infile,"%lf",&eqn_params->K);
        (void) printf("%20.19f\n",eqn_params->K);

        FT_VectorMemoryAlloc((POINTER*)&gindex,8*num_drops,sizeof(int));
        FT_VectorMemoryAlloc((POINTER*)&radius,8*num_drops,sizeof(double));
        FT_MatrixMemoryAlloc((POINTER*)&center,8*num_drops,MAXD,sizeof(double));

        (void) printf("Two methods for initialization:\n");
        (void) printf("\tPrompt initialization (P)\n");
        (void) printf("\tRandom initialization (R)\n");
        CursorAfterString(infile,"Enter method of initialization:");
        fscanf(infile,"%s",string);
        (void) printf("%s\n",string);
        switch (string[0])
        {
        case 'p':
        case 'P':
            for (i = 0; i < num_drops; ++i)
            {
                sprintf(msg,"Enter center of water drop %d:",i+1);
                CursorAfterString(infile,msg);
                for (j = 0; j < dim; ++j)
                {
                    fscanf(infile,"%lf ",&center[i][j]);
                    (void) printf("%f ",center[i][j]);
                }
                (void) printf("\n");
                sprintf(msg,"Enter radius of water drop %d:",i+1);
                CursorAfterString(infile,msg);
                fscanf(infile,"%lf",&radius[i]);
                (void) printf("%f\n",radius[i]);
            }
            break;
        case 'r':
        case 'R':
            sprintf(msg,"Enter mean radius of water drop:");
            CursorAfterString(infile,msg);
            fscanf(infile,"%lf",&r_bar);
            (void) printf("%f\n",r_bar);
            sprintf(msg,"Enter standard deviation of radius:");
            CursorAfterString(infile,msg);
            fscanf(infile,"%lf",&sigma);
            (void) printf("%f\n",sigma);
            initRandomDrops(front,center,radius,gindex,&num_drops,r_bar,sigma);
            break;
        default:
            (void) printf("Unknown option for initialization!\n");
            clean_up(ERROR);
        }
        for (i = 0; i < num_drops; ++i)
        {
            double radii[MAXD];
            for (j = 0; j < dim; ++j)
                radii[j] = radius[i];
            if (dim == 2)
            {
                FT_MakeEllipticCurve(front,center[i],radii,SOLID_COMP,
                        LIQUID_COMP,w_type,1,&curve);
                Gindex(curve) = gindex[i] + 10;
		node_type(curve->start) = CLOSED_NODE;
		body_index(curve) = gindex[i];
		total_mass(curve) = drop_dens*PI*4.0/3.0*
				radius[i]*radius[i]*radius[i];
		spherical_radius(curve) = radius[i];
		for(l = 0; l < dim; ++l)
		{
		    center_of_mass(curve)[l] = center[i][l];
		    center_of_mass_velo(curve)[l] = 0;
		}
		motion_type(curve) = COM_MOTION;
            }
            else if (dim == 3)
            {
                SURFACE **s;
                FT_MakeEllipticSurf(front,center[i],radii,SOLID_COMP,
                        LIQUID_COMP,w_type,1,&surf);
                Gindex(surf) = gindex[i] + 10;
                body_index(surf) = gindex[i];
                total_mass(Hyper_surf(surf)) = drop_dens*PI*4.0/3.0*
                                radius[i]*radius[i]*radius[i];
                for (l = 0; l < dim; ++l)
                {
                    center_of_mass(surf)[l] = center[i][l];
                    center_of_mass_velo(surf)[l] = 0.0;
                }
                motion_type(Hyper_surf(surf)) = COM_MOTION;
                spherical_radius(Hyper_surf(surf)) = radius[i];
                intfc_surface_loop(front->interf,s)
                {
                    if (*s == surf) continue;
                    if (Gindex(surf) == Gindex(*s))
                        I_AddTwoSurfaces(*s,surf);
                }
            }
        }
        FT_FreeThese(3,radius,gindex,center);
        if (debugging("init_intfc"))
        {
            if (dim == 2)
                xgraph_2d_intfc("test.xg",front->interf);
            else if (dim == 3)
                gview_plot_interface("init_intfc",front->interf);
        }
}       /* end initWaterDrops */

static void initRandomDrops(
        Front *front,
        double **center,
        double *radius,
        int *gindex,
        int *num_drops,
        double r_bar,
        double sigma)
{
        int i,j,ii,jj,n,dir,side;
        GAUSS_PARAMS gauss_params;
        UNIFORM_PARAMS uniform_params;
        unsigned short int xsubi[3];
        double x,dist,R;
        int dim = FT_Dimension();
        double *L = front->rect_grid->L;
        double *U = front->rect_grid->U;
        double *h = front->rect_grid->h;
        double T[MAXD];
        boolean periodic_pair_passed;
        int np;                 // number of periodic image
        int n0,num_d0;          // number of true drops (without periodics)
        double **pcenter;       // centers of periodic image
        double min_h;

        xsubi[0] = 10;
        xsubi[1] = 100;
        xsubi[2] = 1000;
        FT_MatrixMemoryAlloc((POINTER*)&pcenter,8,MAXD,sizeof(double));
        num_d0 = *num_drops;

        gauss_params.mu = r_bar;
        gauss_params.sigma = sigma;
        uniform_params.a = 0.0;
        uniform_params.b = 1.0;

        min_h = h[0];
        for (i = 0; i < dim; ++i)
        {
            T[i] = U[i] - L[i];
            if (min_h > h[i]) min_h = h[i];
        }

        n = n0 = 0;
        for (i = 0; i < 10*num_d0; ++i)
        {
            for (j = 0; j < dim; ++j)
            {
                x = dist_uniform((POINTER)&uniform_params,xsubi);
                center[n][j] = L[j] + x*(U[j] - L[j]);
            }
            R = radius[n] = gauss_center_limit((POINTER)&gauss_params,xsubi);
            if (R < 2*min_h) continue;
            for (j = 0; j < n; ++j)
            {
                dist = distance_between_positions(center[j],center[n],dim);
                if (dist < (radius[j] + radius[n]))
                    break;
            }
            if (j < n) continue;

            for (jj = 0; jj < dim; ++jj)
                pcenter[0][jj] = center[n][jj];
            np = 1;
            periodic_pair_passed = YES;
            for (dir = 0; dir < dim; ++dir)
            {
                if (FT_BoundaryType(dir,0) == PERIODIC_BOUNDARY)
                {
                    for (ii = 0; ii < np; ++ii)
                    {
                        for (jj = 0; jj < dim; ++jj)
                            pcenter[np+ii][jj] = pcenter[ii][jj];
                        if (pcenter[ii][dir] > 0.5*(L[dir] + U[dir]))
                            pcenter[np+ii][dir] -= T[dir];
                        else
                            pcenter[np+ii][dir] += T[dir];
                    }
                    for (ii = 0; ii < np; ++ii)
                    for (jj = 0; jj < n; ++jj)
                    {
                        dist = distance_between_positions(pcenter[np+ii],
                                        center[jj],dim);
                        if (dist < (radius[jj] + radius[n]))
                                periodic_pair_passed = NO;
                    }
                    if (periodic_pair_passed == NO) break;
                    np *= 2;
                }
            }
            if (periodic_pair_passed == NO) continue;
            else
            {
                for (ii = 0; ii < np; ++ii)
                {
                    for (jj = 0; jj < dim; ++jj)
                    {
                        center[n][jj] = pcenter[ii][jj];
                        gindex[n] = n0;
                        radius[n] = R;
                    }
                    n++;
                }
            }
            n0++;
            if (n0 == num_d0) break;
        }
        *num_drops = n;
        FT_FreeThese(1,pcenter);
}       /* end initRandomDrops */

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

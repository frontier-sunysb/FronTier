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

#include <cFluid.h>

extern void read_cFluid_params(
	char *inname,
	EQN_PARAMS *eqn_params)
{
	char string[100];
	FILE *infile = fopen(inname,"r");
	int i,dim = eqn_params->dim;

	CursorAfterString(infile,"Enter numerical scheme for interior solver:");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	switch (string[0])
	{
	case 'T':
	case 't':
	    switch (string[4])
	    {
	    case '1':
		eqn_params->num_scheme = TVD_FIRST_ORDER;
		break;
	    case '2':
		eqn_params->num_scheme = TVD_SECOND_ORDER;
		break;
	    case '4':
		eqn_params->num_scheme = TVD_FOURTH_ORDER;
		break;
	    default:
		printf("Numerical scheme %s not implemented!\n",string);
		clean_up(ERROR);
	    }
	    break;
	case 'W':
	case 'w':
	    switch (string[5])
	    {
	    case '1':
		eqn_params->num_scheme = WENO_FIRST_ORDER;
		break;
	    case '2':
		eqn_params->num_scheme = WENO_SECOND_ORDER;
		break;
	    case '4':
		eqn_params->num_scheme = WENO_FOURTH_ORDER;
		break;
	    default:
		printf("Numerical scheme %s not implemented!\n",string);
		clean_up(ERROR);
	    }
	    break;
	default:
	    printf("Numerical scheme %s not implemented!\n",string);
	    clean_up(ERROR);
	}
	CursorAfterString(infile,"Enter order of point propagator:");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	switch (string[0])
	{
	case '1':
	    eqn_params->point_prop_scheme = FIRST_ORDER;
	    break;
	case '2':
	    eqn_params->point_prop_scheme = SECOND_ORDER;
	    break;
	case '4':
	    eqn_params->point_prop_scheme = FOURTH_ORDER;
	    break;
	default:
	    printf("Point propagator order %s not implemented!\n",string);
	    clean_up(ERROR);
	}

	eqn_params->use_base_soln = NO;
	if (CursorAfterStringOpt(infile,
		"Enter yes for comparison with base data:"))
	{
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'y' || string[0] == 'Y')
            {
            	eqn_params->use_base_soln = YES;
	    	CursorAfterString(infile,"Enter base directory name:");
            	fscanf(infile,"%s",eqn_params->base_dir_name);
            	(void) printf("%s\n",eqn_params->base_dir_name);
            	CursorAfterString(infile,"Enter number of comparing steps:");
            	fscanf(infile,"%d",&eqn_params->num_step);
            	(void) printf("%d\n",eqn_params->num_step);
            	FT_VectorMemoryAlloc((POINTER*)&eqn_params->steps,
                                eqn_params->num_step,sizeof(int));
            	for (i = 0; i < eqn_params->num_step; ++i)
            	{
                    sprintf(string,"Enter index of step %d:",i+1);
                    CursorAfterString(infile,string);
                    fscanf(infile,"%d",&eqn_params->steps[i]);
                    (void) printf("%d\n",eqn_params->steps[i]);
            	}
            	FT_ScalarMemoryAlloc((POINTER*)&eqn_params->f_basic,
                                sizeof(F_BASIC_DATA));
            	eqn_params->f_basic->dim = dim;
	    }
	}

	assert(eqn_params->prob_type != ERROR_TYPE);
	fclose(infile);

	if (eqn_params->use_base_soln == YES)
	    FT_ReadComparisonDomain(inname,eqn_params->f_basic);
}	/* end read_cFluid_params */

extern void read_movie_options(
	Front *front)
{
	static MOVIE_OPTION *movie_option;
	FILE *infile;
	char string[100];
	EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
	char *inname;

	inname = InName(front);
	infile = fopen(inname,"r");
	FT_ScalarMemoryAlloc((POINTER*)&movie_option,sizeof(MOVIE_OPTION));
	eqn_params->movie_option = movie_option;
	movie_option->set_bounds = NO;	// default
	if (CursorAfterStringOpt(infile,"Type y to set movie bounds:"))
	{
	    fscanf(infile,"%s",string);
	    (void) printf("%s\n",string);
	    if (string[0] == 'Y' || string[0] == 'y')
		movie_option->set_bounds = YES;
	}
	front->hdf_cut_frame = NO;
	if (CursorAfterStringOpt(infile,"Type y to add cut frame:"))
	{
	    fscanf(infile,"%s",string);
	    (void) printf("%s\n",string);
	    if (string[0] == 'Y' || string[0] == 'y')
	    {
		front->hdf_cut_frame = YES;
		CursorAfterString(infile,"Enter x-bounds for cut frame:");
	    	fscanf(infile,"%lf %lf",&front->cut_L[0],&front->cut_U[0]);
		(void) printf("%f %f\n",front->cut_L[0],front->cut_U[0]);
		CursorAfterString(infile,"Enter y-bounds for cut frame:");
	    	fscanf(infile,"%lf %lf",&front->cut_L[1],&front->cut_U[1]);
		(void) printf("%f %f\n",front->cut_L[1],front->cut_U[1]);
		if (CursorAfterStringOpt(infile,"Enter resolution level:"))
		{
		    fscanf(infile,"%d",&front->resolution_level);
		    (void) printf("%d\n",front->resolution_level);
		}
	    }
	}

	CursorAfterString(infile,"Type y to make movie of density:");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	if (string[0] == 'Y' || string[0] == 'y')
	{
	    movie_option->plot_dens = YES;
	    if (movie_option->set_bounds)
	    {
		CursorAfterString(infile,"Enter min and max density:");
		fscanf(infile,"%lf %lf",&movie_option->min_dens,
				&movie_option->max_dens);
		(void) printf("%f %f\n",movie_option->min_dens,
				movie_option->max_dens);
	    }
	}
	CursorAfterString(infile,"Type y to make movie of pressure:");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	if (string[0] == 'Y' || string[0] == 'y')
	{
	    movie_option->plot_pres = YES;
	    if (movie_option->set_bounds)
	    {
		CursorAfterString(infile,"Enter min and max pressure:");
		fscanf(infile,"%lf %lf",&movie_option->min_pres,
				&movie_option->max_pres);
		(void) printf("%f %f\n",movie_option->min_pres,
				movie_option->max_pres);
	    }
	}
	if (eqn_params->dim != 1)
	{
	    CursorAfterString(infile,"Type y to make movie of vorticity:");
	    fscanf(infile,"%s",string);
	    (void) printf("%s\n",string);
	    if (string[0] == 'Y' || string[0] == 'y')
	    {
	    	movie_option->plot_vort = YES;
	    }
	}
	CursorAfterString(infile,"Type y to make movie of velocity:");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	if (string[0] == 'Y' || string[0] == 'y')
	{
	    movie_option->plot_velo = YES;
	    if (movie_option->set_bounds)
	    {
		CursorAfterString(infile,"Enter min and max x-velocity:");
		fscanf(infile,"%lf %lf",&movie_option->min_xvel,
				&movie_option->max_xvel);
		(void) printf("%f %f\n",movie_option->min_xvel,
				movie_option->max_xvel);
		CursorAfterString(infile,"Enter min and max y-velocity:");
		fscanf(infile,"%lf %lf",&movie_option->min_yvel,
				&movie_option->max_yvel);
		(void) printf("%f %f\n",movie_option->min_yvel,
				movie_option->max_yvel);
	    }
	}
	CursorAfterString(infile,"Type y to make movie of mach number:");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	if (string[0] == 'Y' || string[0] == 'y')
	{
	    movie_option->plot_mach = YES;
	    if (movie_option->set_bounds)
	    {
		CursorAfterString(infile,"Enter min and max mach number:");
		fscanf(infile,"%lf %lf",&movie_option->min_mach,
				&movie_option->max_mach);
		(void) printf("%f %f\n",movie_option->min_mach,
				movie_option->max_mach);
	    }
	}

	if (eqn_params->dim == 3)
	{
	    CursorAfterString(infile,"Type y to make yz cross section movie:");
	    fscanf(infile,"%s",string);
	    (void) printf("%s\n",string);
	    if (string[0] == 'Y' || string[0] == 'y')
		movie_option->plot_cross_section[0] = YES;
	    CursorAfterString(infile,"Type y to make xz cross section movie:");
	    fscanf(infile,"%s",string);
	    (void) printf("%s\n",string);
	    if (string[0] == 'Y' || string[0] == 'y')
		movie_option->plot_cross_section[1] = YES;
	    CursorAfterString(infile,"Type y to make xy cross section movie:");
	    fscanf(infile,"%s",string);
	    (void) printf("%s\n",string);
	    if (string[0] == 'Y' || string[0] == 'y')
		movie_option->plot_cross_section[2] = YES;
	}
	fclose(infile);
}	/* end read_movie_options */

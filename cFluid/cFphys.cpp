#include <cFluid.h>

extern void read_cFluid_params(
	char *inname,
	EQN_PARAMS *eqn_params)
{
	char string[100];
	FILE *infile = fopen(inname,"r");
	int i,dim = eqn_params->dim;

	eqn_params->prob_type = ERROR_TYPE;
	eqn_params->tracked = YES;		// Default
	CursorAfterString(infile,"Enter problem type:");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	if (string[0] == 'T' || string[0] == 't')
	{
	    if (string[10] == 'B' || string[10] == 'b')
	    	eqn_params->prob_type = TWO_FLUID_BUBBLE;
	    else if (string[10] == 'R' || string[10] == 'r')
	    {
		if (string[11] == 'T' || string[11] == 't')
	    	    eqn_params->prob_type = TWO_FLUID_RT;
		else if (string[11] == 'M' || string[11] == 'm')
	    	    eqn_params->prob_type = TWO_FLUID_RM;
	    }
	} 
	else if (string[0] == 'F' || string[0] == 'f')
	    eqn_params->prob_type = FLUID_SOLID_CIRCLE;
	else if (string[0] == 'B' || string[0] == 'b')
	    eqn_params->prob_type = BUBBLE_SURFACE;
	else if (string[0] == 'I' || string[0] == 'i')
	    eqn_params->prob_type = IMPLOSION;
	else if (string[0] == 'P' || string[0] == 'p')
	    eqn_params->prob_type = PROJECTILE;
	else if (string[0] == 'R' || string[0] == 'r')
	    eqn_params->prob_type = RIEMANN_PROB;
	else if (string[0] == 'M' || string[0] == 'm')
	    eqn_params->prob_type = MT_FUSION;
	else if (string[0] == 'O' || string[0] == 'o')
	{
	    if (string[5] == 'S' || string[5] == 's')
	    	eqn_params->prob_type = ONED_SSINE;
	    else if (string[5] == 'B' || string[5] == 'b')
	    	eqn_params->prob_type = ONED_BLAST;
	}
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
            	eqn_params->use_base_soln = YES;
	}
	if (eqn_params->use_base_soln == YES)
        {
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

	assert(eqn_params->prob_type != ERROR_TYPE);
	fclose(infile);

	if (eqn_params->use_base_soln == YES)
	    FT_ReadComparisonDomain(inname,eqn_params->f_basic);
}	/* end read_cFluid_params */

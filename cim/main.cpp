#include "cim.h"
#include "solver.h"

static boolean cim_driver(Front*,C_CARTESIAN&);
static void initCimTestParams(char*,Front*);

char *in_name,*restart_state_name,*restart_name,*out_name;
boolean RestartRun;
boolean ReSetTime;
int RestartStep;
boolean binary = YES;

static void CIM_flowThroughBoundaryState(double*,HYPER_SURF*,Front*,
			POINTER,POINTER);
static void CIM_timeDependBoundaryState(double*,HYPER_SURF*,Front*,
                        POINTER,POINTER);
static void read_dirichlet_bdry_data(char*,Front*,F_BASIC_DATA);

int main(int argc, char **argv)
{
	static Front front;
        static F_BASIC_DATA f_basic;
        static LEVEL_FUNC_PACK level_func_pack;
	CIRCLE_PARAMS circle_params;
	C_CARTESIAN       c_cartesian(front);
	static CIM_PARAMS cim_params;

	FT_Init(argc,argv,&f_basic);
	f_basic.size_of_intfc_state = sizeof(STATE);

        PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);

	/* Initialize basic computational data */

        in_name                 = f_basic.in_name;
        restart_state_name      = f_basic.restart_state_name;
        out_name                = f_basic.out_name;
        restart_name            = f_basic.restart_name;
        RestartRun              = f_basic.RestartRun;
        RestartStep             = f_basic.RestartStep;
        ReSetTime               = f_basic.ReSetTime;

        sprintf(restart_state_name,"%s/state.ts%s",restart_name,
                        right_flush(RestartStep,7));
        sprintf(restart_name,"%s/intfc-ts%s",restart_name,
                        right_flush(RestartStep,7));
        if (pp_numnodes() > 1)
        {
            sprintf(restart_name,"%s-nd%s",restart_name,
                                right_flush(pp_mynode(),4));
            sprintf(restart_state_name,"%s-nd%s",restart_state_name,
                                right_flush(pp_mynode(),4));
        }

        FT_ReadSpaceDomain(in_name,&f_basic);
        FT_StartUp(&front,&f_basic);
        FT_InitDebug(in_name);
	front.extra1 = (POINTER)&cim_params;
	initCimTestParams(in_name,&front);
	if (debugging("trace"))
            (void) printf("Passed FT_StartUp()\n");

	initCimIntfcParams(in_name,&front,&level_func_pack);

	c_cartesian.w_type = level_func_pack.wave_type;
	c_cartesian.neg_comp = level_func_pack.neg_component;
	c_cartesian.pos_comp = level_func_pack.pos_component;

	FT_InitIntfc(&front,&level_func_pack);
	read_dirichlet_bdry_data(in_name,&front,f_basic);
        if (f_basic.dim == 2)
	{
	    char xg_name[100];
            FT_ClipIntfcToSubdomain(&front);
	    sprintf(xg_name,"%s/init_intfc-%d",out_name,pp_mynode());
            xgraph_2d_intfc(xg_name,front.interf);
	}
	else if (f_basic.dim == 3)
	{
	    char gv_name[100];
	    sprintf(gv_name,"%s/init_intfc-%d",out_name,pp_mynode());
	    gview_plot_interface(gv_name,front.interf);
	}

	cim_driver(&front,c_cartesian);
}

static boolean cim_driver(
	Front *front,
	C_CARTESIAN &c_cartesian)
{
	c_cartesian.solve();
	c_cartesian.initMovieVariables();
	FT_AddMovieFrame(front,out_name,NO);
	viewTopVariable(front,c_cartesian.numer_soln,NO,0.0,0.0,
			(char*)"test_dir",(char*)"solution");
	return YES;
}	/* end cim_driver */

static void read_dirichlet_bdry_data(
	char *inname,
	Front *front,
	F_BASIC_DATA f_basic)
{
	char msg[100],s[100];
	int i,dim = front->rect_grid->dim;
	FILE *infile = fopen(inname,"r");
	static STATE *state;
	POINTER func_params;
	HYPER_SURF *hs;

	for (i = 0; i < dim; ++i)
	{
	    if (f_basic.boundary[i][0] == DIRICHLET_BOUNDARY)
	    {
		hs = NULL;
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
		    FT_ScalarMemoryAlloc((POINTER*)&state,sizeof(STATE));
		    CursorAfterString(infile,"Enter boundary u:");
		    fscanf(infile,"%lf",&state->u);
		    (void) printf("%f\n",state->u);
		    FT_SetDirichletBoundary(front,NULL,NULL,
					NULL,(POINTER)state,hs);
		    break;
		case 'f':			// Flow through state
		case 'F':
		    FT_SetDirichletBoundary(front,CIM_flowThroughBoundaryState,
					"flowThroughBoundaryState",NULL,
					NULL,hs);
		    break;
		case 't':			// Flow through state
		case 'T':
		    //get_time_dependent_params(dim,infile,&func_params);
		    FT_SetDirichletBoundary(front,CIM_timeDependBoundaryState,
					"CIM_timeDependBoundaryState",
					func_params,NULL,hs);
		    break;
		default:
		    (void) printf("Unknown Dirichlet boundary!\n");
		    clean_up(ERROR);
		}
	    }
            if (f_basic.boundary[i][1] == DIRICHLET_BOUNDARY)
	    {
		hs = NULL;
                if (rect_boundary_type(front->interf,i,1) == DIRICHLET_BOUNDARY)                    hs = FT_RectBoundaryHypSurf(front->interf,DIRICHLET_BOUNDARY,
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
		    FT_ScalarMemoryAlloc((POINTER*)&state,sizeof(STATE));
		    CursorAfterString(infile,"Enter boundary u:");
		    fscanf(infile,"%lf",&state->u);
		    (void) printf("%f\n",state->u);
		    FT_SetDirichletBoundary(front,NULL,NULL,
					NULL,(POINTER)state,hs);
		    break;
		case 'f':			// Flow through state
		case 'F':
		    FT_SetDirichletBoundary(front,CIM_flowThroughBoundaryState,
					"flowThroughBoundaryState",NULL,
					NULL,hs);
		    break;
		case 't':			// Flow through state
		case 'T':
		    //get_time_dependent_params(dim,infile,&func_params);
		    FT_SetDirichletBoundary(front,CIM_timeDependBoundaryState,
					"CIM_timeDependBoundaryState",
					func_params,NULL,hs);
		    break;
		default:
		    (void) printf("Unknown Dirichlet boundary!\n");
		    clean_up(ERROR);
		}
	    }
	}
	fclose(infile);
}	/* end read_dirichlet_bdry_data */


static void CIM_flowThroughBoundaryState(
        double          *p0,
        HYPER_SURF      *hs,
        Front           *front,
        POINTER         params,
        POINTER         state)
{
	printf("Entering CIM_flowThroughBoundaryState()\n");
}	/* CIM_flowThroughBoundaryState */

static void CIM_timeDependBoundaryState(
        double          *p0,
        HYPER_SURF      *hs,
        Front           *front,
        POINTER         params,
        POINTER         state)
{
}	/* end CIM_timeDependBoundaryState */

static void initCimTestParams(
	char *inname,
	Front *front)
{
	CIM_PARAMS *cim_params = (CIM_PARAMS*)front->extra1;
	FILE *infile = fopen(inname,"r");
	int i,dim = front->rect_grid->dim;

	cim_params->dim = dim;
	for (i = 0; i < dim; ++i)
	    cim_params->h[i] = front->rect_grid->h[i];
	CursorAfterString(infile,"Enter run case number:");
	fscanf(infile,"%d",&cim_params->Run_case);
	(void) printf("%d\n",cim_params->Run_case);
	if (dim == 2)
	{
	    if (cim_params->Run_case == 8)
	    {
		(void) printf("Run case %d not implemented in 2D!\n",
				cim_params->Run_case);
		clean_up(ERROR);
	    }
	}
	else if (dim == 3)
	{
	    if (cim_params->Run_case == 3 ||
		cim_params->Run_case == 4 ||
		cim_params->Run_case == 5 ||
		cim_params->Run_case == 6 ||
		cim_params->Run_case == 7) 
	    {
		(void) printf("Run case %d not implemented in 3D!\n",
				cim_params->Run_case);
		clean_up(ERROR);
	    }
	}

	CursorAfterString(infile,"Enter interface number:");
	fscanf(infile,"%d",&cim_params->intfc_num);
	(void) printf("%d\n",cim_params->intfc_num);
	fclose(infile);
}	/* end initCimTestParams */

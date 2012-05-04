#include "cim.h"
#include "solver.h"

static boolean cim_driver(Front*,C_CARTESIAN&);

char *in_name,*restart_state_name,*restart_name,*out_name;
boolean RestartRun;
boolean ReSetTime;
int RestartStep;
boolean binary = YES;

int main(int argc, char **argv)
{
	static Front front;
        static F_BASIC_DATA f_basic;
        static LEVEL_FUNC_PACK level_func_pack;
	CIRCLE_PARAMS circle_params;
	C_CARTESIAN       c_cartesian(front);

	FT_Init(argc,argv,&f_basic);
	f_basic.size_of_intfc_state = 0;

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
	if (debugging("trace"))
            (void) printf("Passed FT_StartUp()\n");

	circle_params.dim = 2;
	circle_params.R = 0.5;
	circle_params.cen[0] = 0.0;
	circle_params.cen[1] = 0.0;
	circle_params.add_plan_surf = NO;
	circle_params.add_perturbation = NO;
	level_func_pack.func = level_circle_func;
	level_func_pack.func_params = (POINTER)&circle_params;
	level_func_pack.wave_type = FIRST_PHYSICS_WAVE_TYPE;
	level_func_pack.neg_component = 1;
        level_func_pack.pos_component = 2;
	FT_InitIntfc(&front,&level_func_pack);
        if (f_basic.dim < 3)
	{
	    char xg_name[100];
            FT_ClipIntfcToSubdomain(&front);
	    sprintf(xg_name,"init_intfc-%d",pp_mynode());
            xgraph_2d_intfc(xg_name,front.interf);
	}

	c_cartesian.initMesh();
	cim_driver(&front,c_cartesian);
}

static boolean cim_driver(
	Front *front,
	C_CARTESIAN &c_cartesian)
{
	c_cartesian.cim_solve();
	return YES;
}	/* end cim_driver */


#include <iFluid.h>
#include <airfoil.h>

extern void initParachuteModules(Front *front)
{
	int i,num_canopy;
	FILE *infile = fopen(InName(front),"r");
	SURFACE *surf;

	if (debugging("trace"))
	    (void) printf("Entering initParachuteModules()\n");

	if (debugging("set_module"))
	    gview_plot_interface("module-step-1",front->interf);

	CursorAfterString(infile,"Enter number of canopy surfaces:");
        fscanf(infile,"%d",&num_canopy);
        (void) printf("%d\n",num_canopy);

	for (i = 0; i < num_canopy; ++i)
	{
	    CgalCanopySurface(infile,front,&surf);
	}

	if (debugging("trace"))
	    (void) printf("Leaving initParachuteModules()\n");
}	/* end initParachuteModules */

extern void initParachuteDefault(
	Front *front)
{
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	af_params->is_parachute_system = YES;
	af_params->num_opt_round = 20;
        af_params->spring_model = MODEL1;
	af_params->gore_len_fac = 1.0;
}	/* end initParachuteDefault */

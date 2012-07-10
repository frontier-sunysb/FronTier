#include <iFluid.h>
#include <airfoil.h>

typedef struct {
        double N[MAXD];         /* normal of the plane */
        double P[MAXD];         /* a point on the plane */
        double cen[MAXD];       /* center of the vent */
        double radius;          /* radius of the vent */
} CONSTR_PARAMS;

typedef struct {
	int dim;
        double L[MAXD];         /* Lower bounds of box */
        double U[MAXD];         /* Upper bounds of box */
} RECT_CONSTR_PARAMS;

static void zero_state(COMPONENT,double*,L_STATE&,int,IF_PARAMS*);
static void setInitialIntfc2d(Front*,LEVEL_FUNC_PACK*,char*);
static void setInitialIntfc3d(Front*,LEVEL_FUNC_PACK*,char*);
static boolean parachute_constr_func(POINTER,double*);
static boolean rect_constr_func(POINTER,double*);
static boolean cross_constr_func(POINTER,double*);
static boolean circle_constr_func(POINTER,double*);
static boolean install_string_and_rotate(INTERFACE*,SURFACE*,POINTER,int);
static boolean install_strings(INTERFACE*,SURFACE*,POINTER,int);
static boolean change_mono_boundary(INTERFACE*,SURFACE*,POINTER,int);
static void rotate_point(POINT*,double*,double,double,boolean);
static void set_side_curves(double*,double*,SURFACE*,CURVE**,
			CURVE**,CURVE**,CURVE**);

extern void setInitialIntfc(
        Front *front,
        LEVEL_FUNC_PACK *level_func_pack,
        char *inname)
{
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	FILE *infile = fopen("inname","w");
	char string[100];

	level_func_pack->wave_type = ELASTIC_BOUNDARY;
	iFparams->m_comp1 = SOLID_COMP;
        iFparams->m_comp2 = LIQUID_COMP2;
	if (CursorAfterStringOpt(infile,
            "Entering yes to set wave type to FIRST_PHYSICS_WAVE_TYPE: "))
	{
	    fscanf(infile,"%s",string);
	    (void) printf("%s\n",string);
	    if (string[0] == 'y' || string[0] == 'Y')
		level_func_pack->wave_type = FIRST_PHYSICS_WAVE_TYPE;
	}
	fclose(infile);

	switch (front->rect_grid->dim)
	{
	case 2:
	    return setInitialIntfc2d(front,level_func_pack,inname);
	case 3:
	    return setInitialIntfc3d(front,level_func_pack,inname);
	}
}	/* end setInitialIntfc */

static void setInitialIntfc3d(
        Front *front,
        LEVEL_FUNC_PACK *level_func_pack,
        char *inname)
{
	char string[100];
	FILE *infile = fopen(inname,"r");
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	int i;
	static ELLIP_PARAMS ellip_params;
	static PLANE_PARAMS plane_params;
	static STRING_PARAMS *string_params;
	double *cen,*rad;
	static CONSTR_PARAMS constr_params;
	static BDRY_PARAMS bdry_params;
	static RECT_CONSTR_PARAMS rect_constr_params;
	static CIRCLE_PARAMS circle_constr_params;
	int num_canopy;

	level_func_pack->set_3d_bdry = YES;
	level_func_pack->neg_component = LIQUID_COMP2;
        level_func_pack->pos_component = LIQUID_COMP2;	
	level_func_pack->func_params = NULL;
        level_func_pack->func = NULL;
	af_params->is_parachute_system = NO;
	af_params->spring_model = MODEL1;	// default
	CursorAfterString(infile,"Enter number of canopy surfaces:");
	fscanf(infile,"%d",&num_canopy);
	(void) printf("%d\n",num_canopy);
	level_func_pack->num_mono_hs = num_canopy;
	FT_VectorMemoryAlloc((POINTER*)&string_params,num_canopy,
			sizeof(STRING_PARAMS));

	(void) printf("Choices of initial surface are:\n");
	(void) printf("\tEllipsoid (E)\n");
	(void) printf("\tPlane (P)\n");
	(void) printf("\tNone (N)\n");
	CursorAfterString(infile,"Enter initial canopy surface type:");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	switch (string[0])
	{
	case 'E':
	case 'e':
	    cen = ellip_params.cen;
	    rad = ellip_params.rad;
	    CursorAfterString(infile,
			"Enter center coordinate of the ellipsoid:");
	    fscanf(infile,"%lf %lf %lf",&cen[0],&cen[1],&cen[2]);
	    (void) printf("%f %f %f\n",cen[0],cen[1],cen[2]);
	    CursorAfterString(infile,
			"Enter three radii of the ellipsoid:");
	    fscanf(infile,"%lf %lf %lf",&rad[0],&rad[1],&rad[2]);
	    (void) printf("%f %f %f\n",rad[0],rad[1],rad[2]);
	    level_func_pack->is_mono_hs = YES;
	    level_func_pack->func = ellipsoid_func;
	    level_func_pack->func_params = (POINTER)&ellip_params;
	    level_func_pack->constr_func = parachute_constr_func;
	    level_func_pack->constr_params = (POINTER)&constr_params;
	    constr_params.N[0] = constr_params.N[1] = 0.0;
	    constr_params.N[2] = 1.0;
	    constr_params.P[0] = cen[0];
	    constr_params.P[1] = cen[1];
	    CursorAfterString(infile,
			"Enter the height of canopy boundary:");
	    fscanf(infile,"%lf",&constr_params.P[2]);
	    (void) printf("%f\n",constr_params.P[2]);
	    constr_params.radius = 0.0;
	    CursorAfterString(infile,"Enter yes to cut a vent on canopy:");
	    fscanf(infile,"%s",string);
	    (void) printf("%s\n",string);
	    if (string[0] == 'y' || string[0] == 'Y')
	    {
	    	constr_params.cen[0] = cen[0];
	    	constr_params.cen[1] = cen[1];
	    	CursorAfterString(infile,"Enter radius of the vent:");
	    	fscanf(infile,"%lf",&constr_params.radius);
	    	(void) printf("%f\n",constr_params.radius);
	    }
	    CursorAfterString(infile,"Enter yes to attach strings to canopy:");
	    fscanf(infile,"%s",string);
	    (void) printf("%s\n",string);
	    if (string[0] == 'y' || string[0] == 'Y')
	    {
	    	level_func_pack->attach_string = YES;
	    	level_func_pack->string_func = install_string_and_rotate;
	    	level_func_pack->string_params = (POINTER)string_params;
		af_params->is_parachute_system = YES;
		for (i = 0; i < num_canopy; ++i)
		{
		    string_params[i].cen[0] = cen[0];
		    string_params[i].cen[1] = cen[1];
	    	    CursorAfterString(infile,"Enter number of chords:");
	    	    fscanf(infile,"%d",&string_params[i].num_strings);
	    	    (void) printf("%d\n",string_params[i].num_strings);
	    	    CursorAfterString(infile,"Enter start angle of chord:");
	    	    fscanf(infile,"%lf",&string_params[i].start_angle);
	    	    (void) printf("%f\n",string_params[i].start_angle);
	    	    CursorAfterString(infile,"Enter initial position of load:");
	    	    fscanf(infile,"%lf %lf %lf",
					&string_params[i].coords_load[0],
					&string_params[i].coords_load[1],
					&string_params[i].coords_load[2]);
	    	    (void) printf("%f %f %f\n",
					string_params[i].coords_load[0],
					string_params[i].coords_load[1],
					string_params[i].coords_load[2]);
	    	    CursorAfterString(infile,"Enter rotation angles:");
	    	    fscanf(infile,"%lf %lf",&string_params[i].theta,
					&string_params[i].phi);
	    	    (void) printf("%f %f\n",string_params[i].theta,
					string_params[i].phi);
		    string_params[i].theta *= PI/180.0;
		    string_params[i].phi *= PI/180.0;
		}
	    }
	    else
		level_func_pack->attach_string = NO;
	    break;
	case 'T':
	case 't':
	    cen = ellip_params.cen;
	    rad = ellip_params.rad;
	    CursorAfterString(infile,
			"Enter vertex coordinate of the paraboloid:");
	    fscanf(infile,"%lf %lf %lf",&cen[0],&cen[1],&cen[2]);
	    (void) printf("%f %f %f\n",cen[0],cen[1],cen[2]);
	    CursorAfterString(infile,
			"Enter coefficients of the paraboloid:");
	    fscanf(infile,"%lf %lf",&rad[0],&rad[1]);
	    (void) printf("%f %f\n",rad[0],rad[1]);
	    rad[0] = -rad[0]; rad[1] = -rad[1];
	    level_func_pack->is_mono_hs = YES;
	    level_func_pack->func = paraboloid_func;
	    level_func_pack->func_params = (POINTER)&ellip_params;
	    level_func_pack->constr_func = parachute_constr_func;
	    level_func_pack->constr_params = (POINTER)&constr_params;
	    constr_params.N[0] = constr_params.N[1] = 0.0;
	    constr_params.N[2] = 1.0;
	    constr_params.P[0] = cen[0];
	    constr_params.P[1] = cen[1];
	    CursorAfterString(infile,
			"Enter the height of canopy boundary:");
	    fscanf(infile,"%lf",&constr_params.P[2]);
	    (void) printf("%f\n",constr_params.P[2]);
	    constr_params.radius = 0.0;
	    CursorAfterString(infile,"Enter yes to cut a vent on canopy:");
	    fscanf(infile,"%s",string);
	    (void) printf("%s\n",string);
	    if (string[0] == 'y' || string[0] == 'Y')
	    {
	    	constr_params.cen[0] = cen[0];
	    	constr_params.cen[1] = cen[1];
	    	CursorAfterString(infile,"Enter radius of the vent:");
	    	fscanf(infile,"%lf",&constr_params.radius);
	    	(void) printf("%f\n",constr_params.radius);
	    }
	    CursorAfterString(infile,"Enter yes to attach strings to canopy:");
	    fscanf(infile,"%s",string);
	    (void) printf("%s\n",string);
	    if (string[0] == 'y' || string[0] == 'Y')
	    {
	    	level_func_pack->attach_string = YES;
	    	level_func_pack->string_func = install_string_and_rotate;
	    	level_func_pack->string_params = (POINTER)string_params;
		af_params->is_parachute_system = YES;
		for (i = 0; i < num_canopy; ++i)
		{
		    string_params[i].cen[0] = cen[0];
		    string_params[i].cen[1] = cen[1];
	    	    CursorAfterString(infile,"Enter number of chords:");
	    	    fscanf(infile,"%d",&string_params[i].num_strings);
	    	    (void) printf("%d\n",string_params[i].num_strings);
	    	    CursorAfterString(infile,"Enter start angle of chord:");
	    	    fscanf(infile,"%lf",&string_params[i].start_angle);
	    	    (void) printf("%f\n",string_params[i].start_angle);
	    	    CursorAfterString(infile,"Enter initial position of load:");
	    	    fscanf(infile,"%lf %lf %lf",
					&string_params[i].coords_load[0],
					&string_params[i].coords_load[1],
					&string_params[i].coords_load[2]);
	    	    (void) printf("%f %f %f\n",
					string_params[i].coords_load[0],
					string_params[i].coords_load[1],
					string_params[i].coords_load[2]);
	    	    CursorAfterString(infile,"Enter rotation angles:");
	    	    fscanf(infile,"%lf %lf",&string_params[i].theta,
					&string_params[i].phi);
	    	    (void) printf("%f %f\n",string_params[i].theta,
					string_params[i].phi);
		    string_params[i].theta *= PI/180.0;
		    string_params[i].phi *= PI/180.0;
		}
	    }
	    else
		level_func_pack->attach_string = NO;
	    break;

	case 'P':
	case 'p':
	    level_func_pack->func_params = (POINTER)&plane_params;
	    level_func_pack->func = plane_func;
	    level_func_pack->is_mono_hs = YES;
	    rect_constr_params.dim = 3;
	    plane_params.P[0] = plane_params.P[1] = 0.0;
	    CursorAfterString(infile,
			"Enter the height of the plane:");
	    fscanf(infile,"%lf",&plane_params.P[2]);
	    (void) printf("%f\n",plane_params.P[2]);
	    plane_params.N[0] = plane_params.N[1] = 0.0;
	    plane_params.N[2] = 1.0;
            CursorAfterString(infile,"Enter type of plane edge:");
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'c' || string[0] == 'c')
            {
	    	level_func_pack->constr_params = (POINTER)&circle_constr_params;
	    	level_func_pack->constr_func = circle_constr_func;
	    	cen = circle_constr_params.cen;
		CursorAfterString(infile,"Enter circle center:");
	        fscanf(infile,"%lf %lf",&circle_constr_params.cen[0],
				&circle_constr_params.cen[1]);
	        (void) printf("%f %f\n",circle_constr_params.cen[0],
				circle_constr_params.cen[1]);
		CursorAfterString(infile,"Enter circle radius:");
	        fscanf(infile,"%lf",&circle_constr_params.R);
	        (void) printf("%f\n",circle_constr_params.R);

		CursorAfterString(infile,
				"Enter yes to attach strings to canopy:");
		fscanf(infile,"%s",string);
		(void) printf("%s\n",string);
                if (string[0] == 'y' || string[0] == 'Y')
                {
                    af_params->is_parachute_system = YES;
                    level_func_pack->attach_string = YES;
                    level_func_pack->string_func = install_string_and_rotate;
                    level_func_pack->string_params = (POINTER)string_params;
                    for (i = 0; i < num_canopy; ++i)
                    {
                        string_params[i].cen[0] = cen[0];
                        string_params[i].cen[1] = cen[1];
                        CursorAfterString(infile,"Enter number of chords:");
                        fscanf(infile,"%d",&string_params[i].num_strings);
                        (void) printf("%d\n",string_params[i].num_strings);
                        CursorAfterString(infile,"Enter start angle of chord:");
                        fscanf(infile,"%lf",&string_params[i].start_angle);
                        (void) printf("%f\n",string_params[i].start_angle);
                        CursorAfterString(infile,
					"Enter initial position of load:");
                        fscanf(infile,"%lf %lf %lf",
                                        &string_params[i].coords_load[0],
                                        &string_params[i].coords_load[1],
                                        &string_params[i].coords_load[2]);
                        (void) printf("%f %f %f\n",
                                        string_params[i].coords_load[0],
                                        string_params[i].coords_load[1],
                                        string_params[i].coords_load[2]);
                        CursorAfterString(infile,"Enter rotation angles:");
                        fscanf(infile,"%lf %lf",&string_params[i].theta,
                                        &string_params[i].phi);
                        (void) printf("%f %f\n",string_params[i].theta,
                                        string_params[i].phi);
                        string_params[i].theta *= PI/180.0;
                        string_params[i].phi *= PI/180.0;
                    }		
		}
		else if (CursorAfterStringOpt(infile,
				"Enter yes to change canopy boundary:"))
		{
		    fscanf(infile,"%s",string);
		    (void) printf("%s\n",string);
                    if (string[0] == 'y' || string[0] == 'Y')
		    {
                    	level_func_pack->attach_string = YES;
                    	level_func_pack->string_func = change_mono_boundary;
	    	    	level_func_pack->string_params = NULL;
		    }
		}
            }
            else if (string[0] == 'r' || string[0] == 'R')
            {
	    	level_func_pack->constr_params = (POINTER)&rect_constr_params;
	    	level_func_pack->constr_func = rect_constr_func;
	        CursorAfterString(infile,
			"Enter the box lower boundary:");
	        fscanf(infile,"%lf %lf %lf",&rect_constr_params.L[0],
			&rect_constr_params.L[1],&rect_constr_params.L[2]);
	        (void) printf("%f %f %f\n",rect_constr_params.L[0],
			rect_constr_params.L[1],rect_constr_params.L[2]);
	        CursorAfterString(infile,
			"Enter the box upper boundary:");
	        fscanf(infile,"%lf %lf %lf",&rect_constr_params.U[0],
			&rect_constr_params.U[1],&rect_constr_params.U[2]);
	        (void) printf("%f %f %f\n",rect_constr_params.U[0],
			rect_constr_params.U[1],rect_constr_params.U[2]);
		if (CursorAfterStringOpt(infile,
			"Enter yes to change canopy boundary:"))
		{
		    fscanf(infile,"%s",string);
		    (void) printf("%s\n",string);
                    if (string[0] == 'y' || string[0] == 'Y')
		    {
                    	level_func_pack->attach_string = YES;
                    	level_func_pack->string_func = change_mono_boundary;
			bdry_params.L[0] = rect_constr_params.L[0];
			bdry_params.L[1] = rect_constr_params.L[1];
			bdry_params.U[0] = rect_constr_params.U[0];
			bdry_params.U[1] = rect_constr_params.U[1];
			for (i = 0; i < 2; ++i)
			{
			    bdry_params.lower_bdry[i] = 
				bdry_params.upper_bdry[i] = NO;	

			    sprintf(string,"For direction %d",i);
			    CursorAfterString(infile,string);
		    	    (void) printf("%s\n",string);

			    CursorAfterString(infile,
				"Enter yes to fix lower boundary:");
		    	    fscanf(infile,"%s",string);
		    	    (void) printf("%s\n",string);
                    	    if (string[0] == 'y' || string[0] == 'Y')
				bdry_params.lower_bdry[i] = YES;	
			    else
			    {
				bdry_params.lower_side[i] = NO_LOAD;
				if (CursorAfterStringOpt(infile,
				   "Enter load type:"))
				{
		    		    fscanf(infile,"%s",string);
		    		    (void) printf("%s\n",string);
				    switch(string[0])
				    {
				    case 'n':
				    case 'N':
					break;
				    case 'f':
				    case 'F':
					bdry_params.lower_side[i] = FREE_LOAD;
					break;
				    case 'r':
				    case 'R':
					bdry_params.lower_side[i] = RIGID_LOAD;
					break;
				    }
			    	    CursorAfterString(infile,
					"Enter total mass:");
		    		    fscanf(infile,"%lf",
					&bdry_params.lower_mass[i]);
		    		    (void) printf("%f\n",
					bdry_params.lower_mass[i]);
				}
			    }

			    CursorAfterString(infile,
				"Enter yes to fix upper boundary:");
		    	    fscanf(infile,"%s",string);
		    	    (void) printf("%s\n",string);
                    	    if (string[0] == 'y' || string[0] == 'Y')
				bdry_params.upper_bdry[i] = YES;	
			    else
			    {
				bdry_params.upper_side[i] = NO_LOAD;
				if (CursorAfterStringOpt(infile,
				   "Enter load type:"))
				{
		    		    fscanf(infile,"%s",string);
		    		    (void) printf("%s\n",string);
				    switch(string[0])
				    {
				    case 'n':
				    case 'N':
					break;
				    case 'f':
				    case 'F':
					bdry_params.upper_side[i] = FREE_LOAD;
					break;
				    case 'r':
				    case 'R':
					bdry_params.upper_side[i] = RIGID_LOAD;
					break;
				    }
			    	    CursorAfterString(infile,
					"Enter total mass:");
		    		    fscanf(infile,"%lf",
					&bdry_params.upper_mass[i]);
		    		    (void) printf("%f\n",
					bdry_params.upper_mass[i]);
				}
			    }
			}
	    	    	level_func_pack->string_params = (POINTER)&bdry_params;
		    }
		}
	    }
	    else if (string[0] == 'x' || string[0] =='X')
	    {
		level_func_pack->constr_params = (POINTER)&rect_constr_params;
	    	level_func_pack->constr_func = cross_constr_func;
	        CursorAfterString(infile,
			"Enter the box lower boundary:");
	        fscanf(infile,"%lf %lf %lf",&rect_constr_params.L[0],
			&rect_constr_params.L[1],&rect_constr_params.L[2]);
	        (void) printf("%f %f %f\n",rect_constr_params.L[0],
			rect_constr_params.L[1],rect_constr_params.L[2]);
	        CursorAfterString(infile,
			"Enter the box upper boundary:");
	        fscanf(infile,"%lf %lf %lf",&rect_constr_params.U[0],
			&rect_constr_params.U[1],&rect_constr_params.U[2]);
	        (void) printf("%f %f %f\n",rect_constr_params.U[0],
			rect_constr_params.U[1],rect_constr_params.U[2]);

		CursorAfterString(infile,
			"Enter yes to attach strings to canopy:");
		fscanf(infile,"%s",string);
		(void) printf("%s\n",string);
		if (string[0] == 'y' || string[0] == 'Y')
		{
		    af_params->is_parachute_system = YES;
	    	    level_func_pack->attach_string = YES;
	    	    level_func_pack->string_func = install_strings;
	    	    level_func_pack->string_params = (POINTER)string_params;
		    for (i = 0; i < num_canopy; ++i)
		    {
	    		CursorAfterString(infile,"Enter number of chords:");
	    		fscanf(infile,"%d",&string_params[i].num_strings);
	    		(void) printf("%d\n",string_params[i].num_strings);
	    		CursorAfterString(infile,
					"Enter initial position of load:");
	    		fscanf(infile,"%lf %lf %lf",
					&string_params[i].coords_load[0],
					&string_params[i].coords_load[1],
					&string_params[i].coords_load[2]);
	    		(void) printf("%f %f %f\n",
					string_params[i].coords_load[0],
					string_params[i].coords_load[1],
					string_params[i].coords_load[2]);
			string_params[i].cen[0] = 
					string_params[i].coords_load[0];
			string_params[i].cen[1] = 
					string_params[i].coords_load[1];
			string_params[i].L[0] = rect_constr_params.L[0];
			string_params[i].L[1] = rect_constr_params.L[1];
			string_params[i].U[0] = rect_constr_params.U[0];
			string_params[i].U[1] = rect_constr_params.U[1];
		    }
		}
	    }
	    break;
	case 'N':
	case 'n':
	    break;
	}
	af_params->pert_params.pert_type = NO_PERT;
	if (CursorAfterStringOpt(infile,
	    "Entering perturbation type: "))
	{
	    fscanf(infile,"%s",string);
	    (void) printf("%s\n",string);
	    switch (string[0])
	    {
	    case 'n':
	    case 'N':
		break;
	    case 'p':
	    case 'P':
		af_params->pert_params.pert_type = PARALLEL_RAND_PERT;
		break;
	    case 'o':
	    case 'O':
		af_params->pert_params.pert_type = ORTHOGONAL_RAND_PERT;
		break;
	    case 'r':
	    case 'R':
		af_params->pert_params.pert_type = RADIAL_PERT;
	    	CursorAfterString(infile,"Enter perturbation center:");
	        fscanf(infile,"%lf %lf",&af_params->pert_params.cen[0],
				&af_params->pert_params.cen[1]);
		(void) printf("%f %f\n",af_params->pert_params.cen[0],
				af_params->pert_params.cen[1]);
	    	CursorAfterString(infile,"Enter perturbation radius:");
	        fscanf(infile,"%lf",&af_params->pert_params.pert_radius);
		(void) printf("%f\n",af_params->pert_params.pert_radius);
	    	CursorAfterString(infile,"Enter perturbation amplitude:");
	        fscanf(infile,"%lf",&af_params->pert_params.pert_amp);
		(void) printf("%f\n",af_params->pert_params.pert_amp);
		break;
	    case 'l':
	    case 'L':
		af_params->pert_params.pert_type = LINEAR_PERT;
	    	CursorAfterString(infile,"Enter perturbation direction:");
	        fscanf(infile,"%d",&af_params->pert_params.dir);
		(void) printf("%d\n",af_params->pert_params.dir);
	    	CursorAfterString(infile,"Enter perturbation center:");
	        fscanf(infile,"%lf",&af_params->pert_params.x0);
		(void) printf("%f\n",af_params->pert_params.x0);
	    	CursorAfterString(infile,"Enter perturbation lower end:");
	        fscanf(infile,"%lf",&af_params->pert_params.xl);
		(void) printf("%f\n",af_params->pert_params.xl);
	    	CursorAfterString(infile,"Enter perturbation upper end:");
	        fscanf(infile,"%lf",&af_params->pert_params.xu);
		(void) printf("%f\n",af_params->pert_params.xu);
	    	CursorAfterString(infile,"Enter perturbation amplitude:");
	        fscanf(infile,"%lf",&af_params->pert_params.pert_amp);
		(void) printf("%f\n",af_params->pert_params.pert_amp);
		break;
	    case 's':
	    case 'S':
		af_params->pert_params.pert_type = SINE_PERT;
		break;
	    }
	}
}	/* end setInitialIntfc3d */

static void setInitialIntfc2d(
        Front *front,
        LEVEL_FUNC_PACK *level_func_pack,
        char *inname)
{
	char string[100];
	FILE *infile = fopen(inname,"r");
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	int *gmax = front->rect_grid->gmax;
        int i,np;
	double x,y,phi,dx,dy,xL,yL,xU,yU,height;
	double cen[2],rad[2],Amp,mu,phi0;

        level_func_pack->is_closed_curve = NO;
	level_func_pack->neg_component = LIQUID_COMP2;
        level_func_pack->pos_component = LIQUID_COMP2;	
        level_func_pack->wave_type = ELASTIC_BOUNDARY;

        level_func_pack->num_points = np = (int)1.25*gmax[0];  /* default */  

	level_func_pack->func_params = NULL;
        level_func_pack->func = NULL;
	af_params->is_parachute_system = NO;
	af_params->spring_model = MODEL1;	// default

	if (CursorAfterStringOpt(infile,"Enter number of point: "))
	{
	    fscanf(infile,"%d",&np);
	    (void) printf("%d\n",np);
	    level_func_pack->num_points = np;
	}
	FT_MatrixMemoryAlloc((POINTER*)&level_func_pack->point_array,np,
				2,sizeof(double));

	(void) printf("Choices of initial surface are:\n");
	(void) printf("\tLine (L)\n");
	(void) printf("\tEllipsoid (E)\n");
	(void) printf("\tSine curve (S)\n");
	CursorAfterString(infile,"Enter initial curve type:");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	switch (string[0])
	{
	case 'L':
	case 'l':
	    CursorAfterString(infile,
                        "Enter coordinates of line start point:");
            fscanf(infile,"%lf %lf",&xL,&yL);
            (void) printf("%f %f\n",xL,yL);
            CursorAfterString(infile,
                        "Enter coordinates of line end point:");
            fscanf(infile,"%lf %lf",&xU,&yU);
            (void) printf("%f %f\n",xU,yU);
            dx = (xU - xL)/(np-1);
            dy = (yU - yL)/(np-1);
            for (i = 0; i < np; ++i)
            {
                level_func_pack->point_array[i][0] = xL + i*dx;
                level_func_pack->point_array[i][1] = yL + i*dy;
            }
	    break;
	case 'E':
	case 'e':
	    CursorAfterString(infile,"Enter center of the ellipse:");
	    fscanf(infile,"%lf %lf",&cen[0],&cen[1]);
	    (void) printf("%f %f\n",cen[0],cen[1]);
	    CursorAfterString(infile,"Enter radii of the ellipse:");
	    fscanf(infile,"%lf %lf",&rad[0],&rad[1]);
	    (void) printf("%f %f\n",rad[0],rad[1]);
	    for (i = 0; i < np; ++i)
            {
                phi = i*PI/(double)(np-1);
                level_func_pack->point_array[i][0] = cen[0] + rad[0]*cos(phi);
                level_func_pack->point_array[i][1] = cen[1] + rad[1]*sin(phi);
            }
	    break;
	case 'S':
	case 's':
	    CursorAfterString(infile,"Enter vertical coordinate of sine wave:");
	    fscanf(infile,"%lf",&height);
	    (void) printf("%f\n",height);
	    CursorAfterString(infile,
			"Enter horizontal end coordinates of plane:");
	    fscanf(infile,"%lf %lf",&xL,&xU);
	    (void) printf("%f %f\n",xL,xU);
	    CursorAfterString(infile,"Enter amplitude of sine wave:");
	    fscanf(infile,"%lf",&Amp);
	    (void) printf("%f\n",Amp);
	    CursorAfterString(infile,"Enter number of period of sine wave:");
	    fscanf(infile,"%lf",&mu);
	    (void) printf("%f\n",mu);
	    CursorAfterString(infile,"Enter start phase of sine wave:");
	    fscanf(infile,"%lf",&phi0);
	    (void) printf("%f\n",phi0);
	    dx = (xU-xL)/(double)(np-1);
	    for (i = 0; i < np; ++i)
            {
                x = xL + i*dx;
                phi = 2.0*mu*PI*x - phi0;
                y = height + Amp*sin(phi);
                level_func_pack->point_array[i][0] = x;
                level_func_pack->point_array[i][1] = y;
            }
	    break;
	}
	CursorAfterString(infile,"Enter string start node type:");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	switch (string[0])
	{
	case 'f':
	case 'F':
	    if (string[1] == 'i' || string[1] == 'I')
	    	af_params->start_type = FIXED_END;
	    else if (string[1] == 'r' || string[1] == 'R')
	    	af_params->start_type = FREE_END;
	    else
	    {
	    	(void) printf("Unknow start node type\n");
	    	clean_up(ERROR);
	    }
	    break;
	case 'l':
	case 'L':
	    af_params->start_type = LOADED_END;
	    break;
	default:
	    (void) printf("Unknow start node type\n");
	    clean_up(ERROR);
	}
	CursorAfterString(infile,"Enter string end node type:");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	switch (string[0])
	{
	case 'f':
	case 'F':
	    if (string[1] == 'i' || string[1] == 'I')
	    	af_params->end_type = FIXED_END;
	    else if (string[1] == 'r' || string[1] == 'R')
	    	af_params->end_type = FREE_END;
	    else
	    {
	    	(void) printf("Unknow end node type\n");
	    	clean_up(ERROR);
	    }
	    break;
	case 'l':
	case 'L':
	    af_params->end_type = LOADED_END;
	    break;
	default:
	    (void) printf("Unknow end node type\n");
	    clean_up(ERROR);
	}

	af_params->pert_params.pert_type = NO_PERT;
	if (CursorAfterStringOpt(infile,
	    "Entering perturbation type: "))
	{
	    fscanf(infile,"%s",string);
	    (void) printf("%s\n",string);
	    switch (string[0])
	    {
	    case 'n':
	    case 'N':
		break;
	    case 'p':
	    case 'P':
		af_params->pert_params.pert_type = PARALLEL_RAND_PERT;
		break;
	    case 'o':
	    case 'O':
		af_params->pert_params.pert_type = ORTHOGONAL_RAND_PERT;
		break;
	    case 's':
	    case 'S':
		af_params->pert_params.pert_type = SINE_PERT;
		break;
	    }
	}
	af_params->pert_params.pert_amp = 0.2;
	if (CursorAfterStringOpt(infile,
	    "Entering perturbation amplitude: "))
	{
	    fscanf(infile,"%lf",&af_params->pert_params.pert_amp);
	    (void) printf("%f\n",af_params->pert_params.pert_amp);
	}
	if (CursorAfterStringOpt(infile,
	    "Entering type of spring model: "))
	{
	    fscanf(infile,"%s",string);
	    (void) printf("%s\n",string);
	    switch (string[0])
	    {
	    case '1':
		af_params->spring_model = MODEL1;
		break;
	    case '2':
		af_params->spring_model = MODEL2;
		break;
	    case '3':
		af_params->spring_model = MODEL3;
		break;
	    default:
		break;
	    }
	}
	fclose(infile);

}	/* end setInitialIntfc2d */


static boolean circle_constr_func(
        POINTER params,
        double *coords)
{
        CIRCLE_PARAMS *circle_constr_params = (CIRCLE_PARAMS*)params;
	double *cen = circle_constr_params->cen;
	double R = circle_constr_params->R;
	double r;
	
	r = sqrt(sqr(coords[0] - cen[0]) + sqr(coords[1] - cen[1]));
	if ( r <= R) return YES;
	else return NO;
}	/* end circle_constr_func */

static boolean rect_constr_func(
        POINTER params,
        double *coords)
{
        RECT_CONSTR_PARAMS *rect_constr_params = (RECT_CONSTR_PARAMS*)params;
	int i;
	int dim = rect_constr_params->dim;
	double *L = rect_constr_params->L;
	double *U = rect_constr_params->U;
	for (i = 0; i < dim-1; ++i)
	{
	    if (coords[i] < L[i] || coords[i] > U[i])
		return NO;
	}
	return YES;
}	/* end rect_constr_func */

static boolean cross_constr_func(
        POINTER params,
        double *coords)
{
        RECT_CONSTR_PARAMS *rect_constr_params = (RECT_CONSTR_PARAMS*)params;
	int i;
	int dim = rect_constr_params->dim;
	double *L = rect_constr_params->L;
	double *U = rect_constr_params->U;
	double LL[3],UU[3];
	int counter1,counter2;

	LL[0] = L[1]; LL[1] = L[0]; LL[2] = L[2];
	UU[0] = U[1]; UU[1] = U[0]; UU[2] = U[2];

	counter1 = counter2 = 0;
	for (i = 0; i < dim-1; ++i)
        {
            if (coords[i] > L[i] && coords[i] < U[i])
		counter1++;
        }
	for (i = 0; i < dim-1; ++i)
	{
	    if (coords[i] > LL[i] && coords[i] < UU[i])
		counter2++;
	}
	if (counter1 == 2 || counter2 == 2)
	    return YES;
        return NO;
}	/* end cross_constr_func */


static boolean parachute_constr_func(
        POINTER params,
        double *coords)
{
        CONSTR_PARAMS *constr_params = (CONSTR_PARAMS*)params;
        double *N,v[MAXD],dist;
        int i;
	double *cen = constr_params->cen;

        N = constr_params->N;
        for (i = 0; i < 3; ++i)
            v[i] = coords[i] - constr_params->P[i];
        dist = Dot3d(N,v);
        if (dist < 0.0) return NO;
	dist = sqrt(sqr(coords[0] - cen[0]) + sqr(coords[1] - cen[1]));
	return (dist > constr_params->radius) ? YES : NO;
}       /* end parachute_constr_func */

static boolean install_string_and_rotate(
	INTERFACE *intfc,
	SURFACE *surf,
	POINTER params,
	int ip)
{
	CURVE **c,*canopy_bdry,**string_curves;
	POINT **pts,**string_pts;
	BOND **bonds,**string_bonds;
	STRING_PARAMS *tmp_params = (STRING_PARAMS*)params;
	STRING_PARAMS *string_params = tmp_params+ip;
	double *cen = string_params->cen,*cload,coords[MAXD];
	double ave_radius_sqr,max_radius_sqr;
	int i,j,k,num_curves,num_points,num_bonds,num_strings;
	int nb;
	double *string_angle,start_angle,d_angle;
	double theta1,theta2,d1,d2,rot_theta,rot_phi;
	NODE *nload,**string_nodes;
	boolean node_moved;
	AF_NODE_EXTRA *extra;
	double spacing,dir[MAXD],*h = computational_grid(intfc)->h;
	BOND *b;

	num_strings = string_params->num_strings;
	start_angle = string_params->start_angle;
	rot_theta = string_params->theta;
	rot_phi   = string_params->phi;
	cload = string_params->coords_load;
	d_angle = 2*PI/num_strings;

	canopy_bdry = NULL;
	nload = make_node(Point(cload));
	FT_ScalarMemoryAlloc((POINTER*)&extra,sizeof(AF_NODE_EXTRA));
	extra->af_node_type = LOAD_NODE;
	nload->extra = (POINTER)extra;

	FT_VectorMemoryAlloc((POINTER*)&string_angle,num_strings,
				sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&string_pts,num_strings,
				sizeof(POINT*));
	FT_VectorMemoryAlloc((POINTER*)&string_bonds,num_strings,
				sizeof(BOND*));
	FT_VectorMemoryAlloc((POINTER*)&string_nodes,num_strings,
				sizeof(NODE*));
	FT_VectorMemoryAlloc((POINTER*)&string_curves,num_strings,
				sizeof(CURVE*));
	for (i = 0; i < num_strings; ++i)
	{
	    string_angle[i] = start_angle + i*d_angle;
	    if (string_angle[i] > 2.0*PI) 
		string_angle[i] -= 2.0*PI;
	}

	num_curves = FT_NumOfSurfCurves(surf);
	FT_VectorMemoryAlloc((POINTER*)&c,num_curves,sizeof(CURVE*));
	FT_ArrayOfSurfCurves(surf,c);

	max_radius_sqr = 0.0;
	for (i = 0; i < num_curves; ++i)
	{
	    if (hsbdry_type(c[i]) != MONO_COMP_HSBDRY)
		continue;
	    num_points = FT_NumOfCurvePoints(c[i]);
	    FT_VectorMemoryAlloc((POINTER*)&pts,num_points,sizeof(POINT*));
	    FT_ArrayOfCurvePoints(c[i],pts);
	    ave_radius_sqr = 0.0;
	    for (j = 0; j < num_points; ++j)
	    {
		ave_radius_sqr += sqr(Coords(pts[j])[0] - cen[0]) +
			      sqr(Coords(pts[j])[1] - cen[1]);
	    }
	    ave_radius_sqr /= (double)num_points;
	    if (ave_radius_sqr > max_radius_sqr)
	    {
		max_radius_sqr = ave_radius_sqr;
		canopy_bdry = c[i];
	    }
	    FT_FreeThese(1,pts);
	}
	FT_FreeThese(1,c);

	num_bonds = FT_NumOfCurveBonds(canopy_bdry);
	FT_VectorMemoryAlloc((POINTER*)&bonds,num_bonds,sizeof(BOND*));
	FT_ArrayOfCurveBonds(canopy_bdry,bonds);
	for (i = 0; i < num_bonds; ++i)
	{
	    theta1 = plane_angle(cen,Coords(bonds[i]->start));
	    theta2 = plane_angle(cen,Coords(bonds[i]->end));
	    if (fabs(theta1 - theta2) > PI) 
	    {
		if (theta2 > theta1) theta2 -= 2*PI;
		else theta1 -= 2*PI;
	    }
	    for (j = 0; j < num_strings; ++j)
	    {
		if (within_interval(theta1,theta2,string_angle[j]) ||
		    within_interval(theta1,theta2,string_angle[j]-2*PI))
		{
	    	    d1 = distance_between_positions(cen,								Coords(bonds[i]->start),2);
	    	    d2 = distance_between_positions(cen,								Coords(bonds[i]->end),2);
		    d1 = 0.5*(d1 + d2);
		    coords[0] = cen[0] + d1*cos(string_angle[j]);
		    coords[1] = cen[1] + d1*sin(string_angle[j]);
		    coords[2] = 0.5*(Coords(bonds[i]->start)[2] + 
					Coords(bonds[i+1]->end)[2]);
		    string_pts[j] = Point(coords);
		    string_bonds[j] = bonds[i];
		} 
	    }
	}
	FT_FreeThese(1,bonds);

	for (i = 0; i < num_strings; ++i)
	{
	    insert_point_in_bond(string_pts[i],string_bonds[i],canopy_bdry);
	    string_bonds[i] = string_bonds[i]->next;
	}
	node_moved = NO;
	for (i = 0; i < num_strings; ++i)
	{
	    canopy_bdry = FT_CurveOfPoint(intfc,string_pts[i]);
	    if (is_closed_curve(canopy_bdry) && !node_moved)
	    {
		move_closed_loop_node(canopy_bdry,string_bonds[i]);
	        string_nodes[i] = FT_NodeOfPoint(intfc,string_pts[i]);
	    	FT_ScalarMemoryAlloc((POINTER*)&extra,sizeof(AF_NODE_EXTRA));
	    	extra->af_node_type = STRING_NODE;
	    	string_nodes[i]->extra = (POINTER)extra;
		node_moved = YES;
		continue;
	    }
	    c = split_curve(string_pts[i],string_bonds[i],canopy_bdry,0,0,0,0);
	    string_nodes[i] = FT_NodeOfPoint(intfc,string_pts[i]);
	    FT_ScalarMemoryAlloc((POINTER*)&extra,sizeof(AF_NODE_EXTRA));
	    extra->af_node_type = STRING_NODE;
	    string_nodes[i]->extra = (POINTER)extra;
	}
	for (i = 0; i < num_strings; ++i)
	{
	    string_curves[i] = make_curve(0,0,string_nodes[i],nload);
	    hsbdry_type(string_curves[i]) = STRING_HSBDRY;
	    spacing = separation(string_nodes[i]->posn,nload->posn,3);
	    for (j = 0; j < 3; ++j)
		dir[j] = (Coords(nload->posn)[j] - 
			Coords(string_nodes[i]->posn)[j])/spacing;
	    nb = (int)spacing/(1.1*h[0]);
	    spacing /= (double)nb;
	    b = string_curves[i]->first;
	    for (j = 1; j < nb; ++j)
	    {
	    	for (k = 0; k < 3; ++k)
		    coords[k] = Coords(string_nodes[i]->posn)[k] + 
					j*dir[k]*spacing;
		insert_point_in_bond(Point(coords),b,string_curves[i]);
		b = b->next;
	    }
	}
	FT_FreeThese(5,string_angle,string_pts,string_bonds,string_nodes,
				string_curves);

	num_points = FT_NumOfSurfPoints(surf);
	FT_VectorMemoryAlloc((POINTER*)&pts,num_points,sizeof(POINT*));
	FT_ArrayOfSurfPoints(surf,pts);
	rotate_point(nload->posn,cload,rot_phi,rot_theta,YES);
	for (i = 0; i < num_points; ++i)
	    rotate_point(pts[i],cload,rot_phi,rot_theta,NO);

	FT_FreeThese(1,pts);
	return YES;
}	/* end install_string_and_rotate */

static boolean install_strings(
	INTERFACE *intfc,
	SURFACE *surf,
	POINTER params,
	int ip)
{
	CURVE **c,*canopy_bdry,**string_curves;
	POINT **pts,**string_pts;
	BOND **bonds,**string_bonds;
	STRING_PARAMS *tmp_params = (STRING_PARAMS*)params;
	STRING_PARAMS *string_params = tmp_params+ip;
	double *cen = string_params->cen,*cload,coords[MAXD];
	double ave_radius_sqr,max_radius_sqr;
	int i,j,k,num_curves,num_points,num_bonds,num_strings;
	int nb;
	double *string_angle,start_angle,d_angle;
	double theta1,theta2,d1,d2,rot_theta,rot_phi;
	NODE *nload,**string_nodes;
	boolean node_moved;
	AF_NODE_EXTRA *extra;
	double spacing,dir[MAXD],*h = computational_grid(intfc)->h;
	BOND *b;
	double *U,*L, width;

	num_strings = string_params->num_strings;
	start_angle = string_params->start_angle;
	rot_theta = string_params->theta;
	rot_phi   = string_params->phi;
	cload = string_params->coords_load;
	d_angle = 2*PI/num_strings;
	U = string_params->U;
	L = string_params->L;
	width = U[1]-L[1];
	printf("width=%f\n",width);

	canopy_bdry = NULL;
	nload = make_node(Point(cload));
	FT_ScalarMemoryAlloc((POINTER*)&extra,sizeof(AF_NODE_EXTRA));
	extra->af_node_type = LOAD_NODE;
	nload->extra = (POINTER)extra;

	FT_VectorMemoryAlloc((POINTER*)&string_angle,num_strings,
				sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&string_pts,num_strings,
				sizeof(POINT*));
	FT_VectorMemoryAlloc((POINTER*)&string_bonds,num_strings,
				sizeof(BOND*));
	FT_VectorMemoryAlloc((POINTER*)&string_nodes,num_strings,
				sizeof(NODE*));
	FT_VectorMemoryAlloc((POINTER*)&string_curves,num_strings,
				sizeof(CURVE*));

	for (i = 0; i < 5; ++i)
	{
	    coords[0] = L[1] + i*(width/4.0); coords[1] = L[0];
	    string_angle[i] = plane_angle(cen,coords);
	    coords[0] = U[0]; coords[1] = L[1] + i*(width/4.0);
	    string_angle[i+5] = plane_angle(cen,coords);
	    coords[0] = L[1] + i*(width/4.0); coords[1] = U[0];
	    string_angle[i+10] = plane_angle(cen,coords);
	    coords[0] = L[0]; coords[1] = L[1] + i*(width/4.0);
	    string_angle[i+15] = plane_angle(cen,coords);
	}

	for (i = 0; i < num_strings; ++i)
        {
            if (string_angle[i] > 2.0*PI) 
        	string_angle[i] -= 2.0*PI;
        }

	num_curves = FT_NumOfSurfCurves(surf);
	FT_VectorMemoryAlloc((POINTER*)&c,num_curves,sizeof(CURVE*));
	FT_ArrayOfSurfCurves(surf,c);

	max_radius_sqr = 0.0;
	for (i = 0; i < num_curves; ++i)
	{
	    if (hsbdry_type(c[i]) != MONO_COMP_HSBDRY)
		continue;
	    num_points = FT_NumOfCurvePoints(c[i]);
	    FT_VectorMemoryAlloc((POINTER*)&pts,num_points,sizeof(POINT*));
	    FT_ArrayOfCurvePoints(c[i],pts);
	    ave_radius_sqr = 0.0;
	    for (j = 0; j < num_points; ++j)
	    {
		ave_radius_sqr += sqr(Coords(pts[j])[0] - cen[0]) +
			      sqr(Coords(pts[j])[1] - cen[1]);
	    }
	    ave_radius_sqr /= (double)num_points;
	    if (ave_radius_sqr > max_radius_sqr)
	    {
		max_radius_sqr = ave_radius_sqr;
		canopy_bdry = c[i];
	    }
	    FT_FreeThese(1,pts);
	}
	FT_FreeThese(1,c);

	num_bonds = FT_NumOfCurveBonds(canopy_bdry);
	FT_VectorMemoryAlloc((POINTER*)&bonds,num_bonds,sizeof(BOND*));
	FT_ArrayOfCurveBonds(canopy_bdry,bonds);
	for (i = 0; i < num_bonds; ++i)
	{
	    theta1 = plane_angle(cen,Coords(bonds[i]->start));
	    theta2 = plane_angle(cen,Coords(bonds[i]->end));
	    if (fabs(theta1 - theta2) > PI) 
	    {
		if (theta2 > theta1) theta2 -= 2*PI;
		else theta1 -= 2*PI;
	    }
	    for (j = 0; j < num_strings; ++j)
	    {
		if (within_interval(theta1,theta2,string_angle[j]) ||
		    within_interval(theta1,theta2,string_angle[j]-2*PI))
		{
	    	    d1 = distance_between_positions(cen,								Coords(bonds[i]->start),2);
	    	    d2 = distance_between_positions(cen,								Coords(bonds[i]->end),2);
		    d1 = 0.5*(d1 + d2);
		    coords[0] = cen[0] + d1*cos(string_angle[j]);
		    coords[1] = cen[1] + d1*sin(string_angle[j]);
		    coords[2] = 0.5*(Coords(bonds[i]->start)[2] + 
					Coords(bonds[i+1]->end)[2]);
		    string_pts[j] = Point(coords);
		    string_bonds[j] = bonds[i];
		} 
	    }
	}
	FT_FreeThese(1,bonds);

	for (i = 0; i < num_strings; ++i)
	{
	    insert_point_in_bond(string_pts[i],string_bonds[i],canopy_bdry);
	    string_bonds[i] = string_bonds[i]->next;
	}
	node_moved = NO;
	for (i = 0; i < num_strings; ++i)
	{
	    canopy_bdry = FT_CurveOfPoint(intfc,string_pts[i]);
	    if (is_closed_curve(canopy_bdry) && !node_moved)
	    {
		move_closed_loop_node(canopy_bdry,string_bonds[i]);
	        string_nodes[i] = FT_NodeOfPoint(intfc,string_pts[i]);
	    	FT_ScalarMemoryAlloc((POINTER*)&extra,sizeof(AF_NODE_EXTRA));
	    	extra->af_node_type = STRING_NODE;
	    	string_nodes[i]->extra = (POINTER)extra;
		node_moved = YES;
		continue;
	    }
	    c = split_curve(string_pts[i],string_bonds[i],canopy_bdry,0,0,0,0);
	    string_nodes[i] = FT_NodeOfPoint(intfc,string_pts[i]);
	    FT_ScalarMemoryAlloc((POINTER*)&extra,sizeof(AF_NODE_EXTRA));
	    extra->af_node_type = STRING_NODE;
	    string_nodes[i]->extra = (POINTER)extra;
	}
	for (i = 0; i < num_strings; ++i)
	{
	    string_curves[i] = make_curve(0,0,string_nodes[i],nload);
	    hsbdry_type(string_curves[i]) = STRING_HSBDRY;
	    spacing = separation(string_nodes[i]->posn,nload->posn,3);
	    for (j = 0; j < 3; ++j)
		dir[j] = (Coords(nload->posn)[j] - 
			Coords(string_nodes[i]->posn)[j])/spacing;
	    nb = (int)spacing/(1.1*h[0]);
	    spacing /= (double)nb;
	    b = string_curves[i]->first;
	    for (j = 1; j < nb; ++j)
	    {
	    	for (k = 0; k < 3; ++k)
		    coords[k] = Coords(string_nodes[i]->posn)[k] + 
					j*dir[k]*spacing;
		insert_point_in_bond(Point(coords),b,string_curves[i]);
		b = b->next;
	    }
	}
	FT_FreeThese(5,string_angle,string_pts,string_bonds,string_nodes,
				string_curves);
	return YES;
}	/* end install_strings */

static boolean change_mono_boundary(
	INTERFACE *intfc,
	SURFACE *surf,
	POINTER params,
	int ip)
{

	if (params != NULL)
	{
	    BDRY_PARAMS *bdry_params = (BDRY_PARAMS*)params;
	    double *L,*U;
	    boolean *lower_bdry = bdry_params->lower_bdry;
	    boolean *upper_bdry = bdry_params->upper_bdry;
	    CURVE *cside00,*cside01,*cside10,*cside11;
	    L = bdry_params->L;
	    U = bdry_params->U;
	    set_side_curves(L,U,surf,&cside00,&cside01,&cside10,&cside11);
	    if (lower_bdry[0] == YES)
		hsbdry_type(cside00) = FIXED_HSBDRY;
	    if (upper_bdry[0] == YES)
		hsbdry_type(cside01) = FIXED_HSBDRY;
	    else
	    {
		static C_PARAMS c_params;
		int npts = FT_NumOfCurvePoints(cside01);
		c_params.load_type = bdry_params->upper_side[0];
		c_params.load_mass = bdry_params->upper_mass[0];
		c_params.point_mass = bdry_params->upper_mass[0]/npts;
		c_params.dir = 0;
		cside01->extra = (POINTER)&c_params;
	    }
	    if (lower_bdry[1] == YES)
		hsbdry_type(cside10) = FIXED_HSBDRY;
	    if (upper_bdry[1] == YES)
		hsbdry_type(cside11) = FIXED_HSBDRY;
	}
	else
	{
	    /* Circular surface, only one boundary */
	    CURVE **c;
	    for (c = surf->pos_curves; c && *c; ++c)
	    {
		if (hsbdry_type(*c) == MONO_COMP_HSBDRY)
		    hsbdry_type(*c) = FIXED_HSBDRY;
	    }
	    for (c = surf->neg_curves; c && *c; ++c)
	    {
		if (hsbdry_type(*c) == MONO_COMP_HSBDRY)
		    hsbdry_type(*c) = FIXED_HSBDRY;
	    }
	}
	return YES;
}	/* end change_mono_boundary */

static void rotate_point(
	POINT *p,
	double *center,			// Rotation center
	double phi,			// Polar angle
	double theta,			// Azimuthal angle
	boolean first)
{
	static double roty_matrix[3][3];
	static double rotz_matrix[3][3];
	double pt[3];
	int i,j;

	if (first == YES)
	{
	    /* Set up rotation matrix */
	    roty_matrix[0][0] = cos(theta);
	    roty_matrix[0][1] = 0.0;
	    roty_matrix[0][2] = sin(theta);
	    roty_matrix[1][0] = 0.0;
	    roty_matrix[1][1] = 1.0;
	    roty_matrix[1][2] = 0.0;
	    roty_matrix[2][0] = -sin(theta);
	    roty_matrix[2][1] = 0.0;
	    roty_matrix[2][2] = cos(theta);
	    
	    rotz_matrix[0][0] = cos(phi);
	    rotz_matrix[0][1] = -sin(phi);
	    rotz_matrix[0][2] = 0.0;
	    rotz_matrix[1][0] = sin(phi);
	    rotz_matrix[1][1] = cos(phi);
	    rotz_matrix[1][2] = 0.0;
	    rotz_matrix[2][0] = 0.0;
	    rotz_matrix[2][1] = 0.0;
	    rotz_matrix[2][2] = 1.0;
	}
	for (i = 0; i < 3; i++)
	    Coords(p)[i] -= center[i];
	for (i = 0; i < 3; i++)
	{
	    pt[i] = 0.0; 
	    for (j = 0; j < 3; j++)
	    {
		pt[i] += roty_matrix[i][j]*Coords(p)[j];
	    }
	}
	for (i = 0; i < 3; i++)
	{
	    Coords(p)[i] = 0.0;
	    for (j = 0; j < 3; j++)
	    {
		Coords(p)[i] += rotz_matrix[i][j]*pt[j];
	    }
	}
	for (i = 0; i < 3; i++)
	    Coords(p)[i] += center[i];
}	/* end rotate_point */

static void zero_state(
        COMPONENT comp,
        double *coords,
        L_STATE& state,
        int dim,
        IF_PARAMS *af_params)
{
        int i;
        for (i = 0; i < dim; ++i)
            state.m_U[i] = 0.0;
        state.m_U[1] = 0.0;
}       /* end zero_state */

extern void setRestartAirfoilIntfc(
	Front *front,
        LEVEL_FUNC_PACK *level_func_pack)
{
	INTERFACE *intfc = front->interf;
	AF_NODE_EXTRA *extra;	
	NODE **n,*payload_node,*string_node;
	CURVE **c;
	boolean is_payload_node;
	STRING_PARAMS *string_params;

	string_params = (STRING_PARAMS*)level_func_pack->string_params;

	payload_node = NULL;
	for (n = intfc->nodes; n && *n; ++n)
	{
	    is_payload_node = YES;
	    for (c = (*n)->in_curves; c && *c; ++c)
	    {
		if (hsbdry_type(*c) != STRING_HSBDRY)
		{
		    is_payload_node = NO;
		    break;
		}
	    } 
	    if (is_payload_node == NO) continue;
	    for (c = (*n)->out_curves; c && *c; ++c)
	    {
		if (hsbdry_type(*c) != STRING_HSBDRY)
		{
		    is_payload_node = NO;
		    break;
		}
	    } 
	    if (is_payload_node == YES)
	    {
		payload_node = *n;
		break;
	    }
	}
	if (payload_node != NULL)
	{
	    FT_ScalarMemoryAlloc((POINTER*)&extra,sizeof(AF_NODE_EXTRA));
	    extra->af_node_type = LOAD_NODE;
	    payload_node->extra = (POINTER)extra;
	}
	for (c = payload_node->in_curves; c && *c; ++c)
	{
	    string_node = (*c)->start;
	    FT_ScalarMemoryAlloc((POINTER*)&extra,sizeof(AF_NODE_EXTRA));
	    extra->af_node_type = STRING_NODE;
	    string_node->extra = (POINTER)extra;
	}
	for (c = payload_node->out_curves; c && *c; ++c)
	{
	    string_node = (*c)->end;
	    FT_ScalarMemoryAlloc((POINTER*)&extra,sizeof(AF_NODE_EXTRA));
	    extra->af_node_type = STRING_NODE;
	    string_node->extra = (POINTER)extra;
	}
}	/* end setRestartAirfoilIntfc */

static void closest_point_on_curve(
	POINT **p_closest,
	BOND **b_closest,
	double *crds,
	CURVE *c)
{
	BOND *b;
	double d,dmin;
	POINT *p,*pmin;
	BOND *bmin;
	pmin = p = c->first->start;
	dmin = d = distance_between_positions(crds,Coords(p),2);
	bmin = c->first;
	for (b = c->first; b != NULL; b = b->next)
	{
	    p = b->end;
	    d = distance_between_positions(crds,Coords(p),2);
	    if (d < dmin)
	    {
		dmin = d;
		pmin = p;
		bmin = b;
	    }
	}
	*p_closest = pmin;
	*b_closest = bmin;
}	/* end closest_point_on_curve */

static boolean point_on_curve(
	POINT *p,
	BOND **b,
	CURVE *c)
{
	BOND *bond;
	if (p == c->first->start)
	{
	    *b = c->first;
	    return YES;
	}
	for (bond = c->first; bond != NULL; bond = bond->next)
	{
	    if (p == bond->end)
	    {
		*b = bond;
		return YES;
	    }
	}
	return NO;
}	/* end point_on_curve */

static void set_side_curves(
	double *L,
	double *U,
	SURFACE *s,
	CURVE **cside00,
	CURVE **cside01,
	CURVE **cside10,
	CURVE **cside11)
{
	BOND *b00,*b01,*b10,*b11;
	POINT *p00,*p01,*p10,*p11;
	double crds[MAXD];
	CURVE **c,*cbdry,*curves[4];
	int i,nc;
	for (c = s->pos_curves; c && *c; ++c)
	{
	    if (hsbdry_type(*c) == MONO_COMP_HSBDRY)
		cbdry = *c;
	}
	for (c = s->neg_curves; c && *c; ++c)
	{
	    if (hsbdry_type(*c) == MONO_COMP_HSBDRY)
		cbdry = *c;
	}

	crds[0] = L[0];	crds[1] = L[1];
	closest_point_on_curve(&p00,&b00,crds,cbdry);
	crds[0] = U[0];	crds[1] = L[1];
	closest_point_on_curve(&p01,&b01,crds,cbdry);
	crds[0] = L[0];	crds[1] = U[1];
	closest_point_on_curve(&p10,&b10,crds,cbdry);
	crds[0] = U[0];	crds[1] = U[1];
	closest_point_on_curve(&p11,&b11,crds,cbdry);

	change_node_of_closed_curve(p00,cbdry);

	c = split_curve(p01,b01,cbdry,0,0,0,0);
	nc = 0;
	curves[nc++] = c[0];
	curves[nc++] = c[1];

	for (i = 0; i < nc; ++i)
	{
	    if (point_on_curve(p10,&b10,curves[i]))
	    {
		c = split_curve(p10,b10,curves[i],0,0,0,0);
		curves[i] = c[0];
	    }
	}
	curves[nc++] = c[1];
	for (i = 0; i < nc; ++i)
	{
	    if (point_on_curve(p11,&b11,curves[i]))
	    {
		c = split_curve(p11,b11,curves[i],0,0,0,0);
		curves[i] = c[0];
	    }
	}
	curves[nc++] = c[1];
	*cside00 = *cside01 = *cside10 = *cside11 = NULL;
	for (i = 0; i < nc; ++i)
	{
	    if (point_on_curve(p00,&b00,curves[i]) &&
		point_on_curve(p10,&b10,curves[i]))
		*cside00 = curves[i];
	    if (point_on_curve(p01,&b01,curves[i]) &&
		point_on_curve(p11,&b11,curves[i]))
		*cside01 = curves[i];
	    if (point_on_curve(p00,&b00,curves[i]) &&
		point_on_curve(p01,&b01,curves[i]))
		*cside10 = curves[i];
	    if (point_on_curve(p10,&b10,curves[i]) &&
		point_on_curve(p11,&b11,curves[i]))
		*cside11 = curves[i];
	}
}	/* end set_side_curves */

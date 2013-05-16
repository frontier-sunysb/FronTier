#include <iFluid.h>
#include <climate.h>

static void prompt_for_rigid_body_params(int,char*,RG_PARAMS*);
static void set_rgbody_params(RG_PARAMS,HYPER_SURF*);
static void zero_state(COMPONENT,double*,IF_FIELD*,int,int,IF_PARAMS*);
static void initWaterDrops2d(Front*);
static void initWaterDrops3d(Front*);

extern void init_fluid_state_func(
	Incompress_Solver_Smooth_Basis *cartesian)
{
	cartesian->getInitialState = zero_state;
}	/* end init_fluid_state_func */

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
}       /* end zero_state */

static	void prompt_for_rigid_body_params(
	int dim,
	char *inname,
	RG_PARAMS *rgb_params)
{
	int i;
	char msg[100],s[100];
	FILE *infile = fopen(inname,"r");

	if (debugging("rgbody")) 
	    (void) printf("Enter prompt_for_rigid_body_params()\n");

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

            CursorAfterString(infile,
			"Type yes if you want vertical motion only?: ");
	    fscanf(infile,"%s",s);
	    (void) printf("%s\n",s);
            if (s[0] == 'y' || s[0] == 'Y')
	    	rgb_params->motion_type = VERTICAL_MOTION;
            else
	    	rgb_params->motion_type = FREE_MOTION;

	    sprintf(msg,"Enter the initial center of mass velocity:");
	    CursorAfterString(infile,msg);
	    for (i = 0; i < dim; ++i)
	    {
	    	fscanf(infile,"%lf",&rgb_params->cen_of_mass_velo[i]);
	    	(void) printf("%f ",rgb_params->cen_of_mass_velo[i]);
	    }

        }
	    
	if (debugging("rgbody")) 
	    (void) printf("Leaving prompt_for_rigid_body_params()\n");
}	/* end prompt_for_rigid_body_params */

static void set_rgbody_params(
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
}	/* end set_rgbody_params */

extern void initWaterDrops(
	Front *front)
{
	int dim = front->rect_grid->dim;
	printf("Entering initWaterDrops()\n");
	switch (dim)
	{
	case 2:
	    initWaterDrops2d(front);
	    clean_up(0);
	    return;
	case 3:
	    initWaterDrops3d(front);
	    return;
	}
}	/* end initWaterDrops */

static void initWaterDrops2d(
	Front *front)
{
	char string[100],msg[200];
	char *inname = InName(front);
	FILE *infile = fopen(inname,"r");
	int i,num_drops;
	double **center,*radius;
	CURVE *curve;

	(void) printf("Two methods for initialization:\n");
	(void) printf("\tPrompt initialization (P)\n");
	(void) printf("\tRandom initialization (R)\n");
	CursorAfterString(infile,"Enter method of initialization:");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	CursorAfterString(infile,"Enter number of water drops:");
	fscanf(infile,"%d",&num_drops);
	(void) printf("%d\n",num_drops);
	FT_VectorMemoryAlloc((POINTER*)&radius,num_drops,sizeof(double));
	FT_MatrixMemoryAlloc((POINTER*)&center,num_drops,MAXD,sizeof(double));
	switch (string[0])
	{
	case 'p':
	case 'P':
	    for (i = 0; i < num_drops; ++i)
	    {
		sprintf(msg,"Enter center of water drop %d:",i+1);
		CursorAfterString(infile,msg);
		fscanf(infile,"%lf %lf",&center[i][0],&center[i][1]);
		(void) printf("%f %f\n",center[i][0],center[i][1]);
		sprintf(msg,"Enter radius of water drop %d:",i+1);
		CursorAfterString(infile,msg);
		fscanf(infile,"%lf",&radius[i]);
		(void) printf("%f\n",radius[i]);
	    }
	    break;
	case 'd':
	case 'D':
	    break;
	default:
	    (void) printf("Unknown option for initialization!\n");
	    clean_up(ERROR);
	}
	for (i = 0; i < num_drops; ++i)
	{
	    double radii[2];
	    radii[0] = radii[1] = radius[i];
	    FT_MakeEllipticCurve(front,center[i],radii,SOLID_COMP,
			LIQUID_COMP2,PARTICLE_BOUNDARY,&curve);
	}
	xgraph_2d_intfc("test.xg",front->interf);
}	/* end initWaterDrops2d */

static void initWaterDrops3d(
	Front *front)
{
}	/* end initWaterDrops3d */

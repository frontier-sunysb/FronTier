#include <iFluid.h>
#include <climate.h>

static void prompt_for_rigid_body_params(int,char*,RG_PARAMS*);
static void set_rgbody_params(RG_PARAMS,HYPER_SURF*);
static void zero_state(COMPONENT,double*,IF_FIELD*,int,int,IF_PARAMS*);
static void initRandomDrops(Front*,double**,double*,int,
				double,double);

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
	char string[100],msg[200];
	char *inname = InName(front);
	FILE *infile = fopen(inname,"r");
	int i,j,dir,num_drops;
	double **center,*radius;
	double r_bar,sigma;
	CURVE *curve;
	SURFACE *surf,*psurf;
	double *L = front->rect_grid->L;
	double *U = front->rect_grid->U;
	double T[MAXD];
	int dim = front->rect_grid->dim;
	int w_type;

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
	FT_VectorMemoryAlloc((POINTER*)&radius,num_drops,sizeof(double));
	FT_MatrixMemoryAlloc((POINTER*)&center,num_drops,MAXD,sizeof(double));

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
	    initRandomDrops(front,center,radius,num_drops,r_bar,sigma);
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
	    	FT_MakeEllipticCurve(front,center[i],radii,SOLID_COMP,
			LIQUID_COMP2,w_type,1,&curve);
	    else if (dim == 3)
	    	FT_MakeEllipticSurf(front,center[i],radii,SOLID_COMP,
			LIQUID_COMP2,w_type,1,&surf);
	    for (dir = 0; dir < dim; ++dir)
	    {
		continue;
		if (FT_BoundaryType(dir,0) == PERIODIC_BOUNDARY)
		{
		    printf("Before copying:\n");
		    printf("number of surfs = %d\n",
				FT_NumOfIntfcSurfaces(surf->interface));
		    psurf = I_CopySurface(surf);
		    printf("After copying:\n");
		    printf("number of surfs = %d\n",
				FT_NumOfIntfcSurfaces(surf->interface));
	    	    for (j = 0; j < dim; ++j)
		    	T[j] = 0.0;
		    if (center[i][dir] > 0.5*(L[dir] + U[dir]))
			T[dir] = L[dir] - U[dir];
		    else
			T[dir] = U[dir] - L[dir];
		    I_ShiftSurface(psurf,T);
		}
	    }
	}
	if (debugging("init_intfc"))
	{
	    if (dim == 2)
	    	xgraph_2d_intfc("test.xg",front->interf);
	    else if (dim == 3)
		gview_plot_interface("init_intfc",front->interf);
	}
}	/* end initWaterDrops */

static void initRandomDrops(
	Front *front,
	double **center,
	double *radius,
	int num_drops,
	double r_bar,
	double sigma)
{
	int i,j,n,dir,side;
	GAUSS_PARAMS gauss_params;
	UNIFORM_PARAMS uniform_params;
	unsigned short int xsubi[3];
	double x,dist;
	int dim = FT_Dimension();
	double *L = front->rect_grid->L;
	double *U = front->rect_grid->U;
	double T[MAXD];
	boolean periodic_pair_passed;

	xsubi[0] = 10;
	xsubi[1] = 100;
	xsubi[2] = 1000;

	n = 0;
	gauss_params.mu = r_bar;
	gauss_params.sigma = sigma;
	uniform_params.a = 0.0;
	uniform_params.b = 1.0;

	for (i = 0; i < dim; ++i)
	    T[i] = U[i] - L[i];
	for (i = 0; i < num_drops*3; ++i)
	{
	    for (j = 0; j < dim; ++j)
	    {
	    	x = dist_uniform((POINTER)&uniform_params,xsubi);
	    	center[n][j] = L[j] + x*(U[j] - L[j]);
	    }
	    radius[n] = gauss_center_limit((POINTER)&gauss_params,xsubi);
	    for (j = 0; j < n; ++j)
	    {
		dist = distance_between_positions(center[j],center[n],dim);
		if (dist < (radius[j] + radius[n]))
		    break;
	    }
	    if (j < n) continue;
	    periodic_pair_passed = YES;
	    for (dir = 0; dir < dim; ++dir)
	    {
		if (FT_BoundaryType(dir,0) == PERIODIC_BOUNDARY)
		{
		    // Test periodic image
		    double pcen[MAXD];
		    for (j = 0; j < dim; ++j)
			pcen[j] = center[n][j];
		    if (center[n][dir] > 0.5*(L[dir] + U[dir]))
			pcen[dir] -= T[dir];
		    else
			pcen[dir] += T[dir];
	    	    for (j = 0; j < n; ++j)
	    	    {
			dist = distance_between_positions(center[j],pcen,dim);
			if (dist < (radius[j] + radius[n]))
				periodic_pair_passed = NO;
		    }
		}
	    }
	    if (periodic_pair_passed == NO) continue;
	    n++;
	    if (n == num_drops) break;
	}
}	/* end initRandomDrops */

#include <iFluid.h>
#include <climate.h>

static void prompt_for_rigid_body_params(int,char*,RG_PARAMS*);
static void set_rgbody_params(RG_PARAMS,HYPER_SURF*);
static void zero_state(COMPONENT,double*,IF_FIELD*,int,int,IF_PARAMS*);
static void ambient_state(COMPONENT,double*,IF_FIELD*,int,int,IF_PARAMS*);
static void initRandomDrops(Front*,double**,double*,int*,int*,
				double,double);

extern void init_fluid_state_func(
	Incompress_Solver_Smooth_Basis *cartesian)
{
	cartesian->getInitialState = ambient_state;
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

static void ambient_state(
        COMPONENT comp,
        double *coords,
        IF_FIELD *field,
        int index,
        int dim,
        IF_PARAMS *iFparams)
{
        int i;
        double **vel = field->vel;
	double *U_ambient = iFparams->U_ambient;
        for (i = 0; i < dim; ++i)
            vel[i][index] = U_ambient[i];
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
	    iFparams->m_comp2 = LIQUID_COMP2;
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
	(void) printf("%d\n",drop_dens);
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
			LIQUID_COMP2,w_type,1,&curve);
		Gindex(curve) = gindex[i] + 10;
	    }
	    else if (dim == 3)
	    {
		SURFACE **s;
	    	FT_MakeEllipticSurf(front,center[i],radii,SOLID_COMP,
			LIQUID_COMP2,w_type,1,&surf);
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
	CursorAfterString(infile,
		"Enter density and viscosity of ambient fluid:");
        fscanf(infile,"%lf %lf",&iFparams->rho2,&iFparams->mu2);
        (void) printf("%f %f\n",iFparams->rho2,iFparams->mu2);
	CursorAfterString(infile,"Enter gravity:");
        for (i = 0; i < dim; ++i)
        {
            fscanf(infile,"%lf",&iFparams->gravity[i]);
            (void) printf("%f ",iFparams->gravity[i]);
        }
        (void) printf("\n");
	FT_FreeThese(3,radius,gindex,center);
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
	int np;			// number of periodic image
	int n0,num_d0;		// number of true drops (without periodics)
	double **pcenter;	// centers of periodic image
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
}	/* end initRandomDrops */

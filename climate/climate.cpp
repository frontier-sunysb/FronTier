#include <iFluid.h>
#include <climate.h>

static void fluid_driver(Front*,Incompress_Solver_Smooth_Basis*);
static int rgbody_vel(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,
                        HYPER_SURF*,double*);

static void prompt_and_make_rigid_body(Front*,FILE*,POINTER,
		double(*func)(POINTER,double*),COMPONENT,COMPONENT,int);
static boolean force_on_hse(HYPER_SURF_ELEMENT*,HYPER_SURF*,RECT_GRID*,double*,
                                        double*,double*,boolean);
static boolean force_on_hse2d(HYPER_SURF_ELEMENT*,HYPER_SURF*,RECT_GRID*,
					double*,double*,double*,boolean);
static boolean force_on_hse3d(HYPER_SURF_ELEMENT*,HYPER_SURF*,RECT_GRID*,
					double*,double*,double*,boolean);
static double intrp_between(double,double,double,double,double);
static void compute_ice_particle_force(Front*,HYPER_SURF*,double,
			double*,double*);
static void compute_ice_particle_force2d(Front*,HYPER_SURF*,double,
			double*,double*);
static void compute_ice_particle_force3d(Front*,HYPER_SURF*,double,
			double*,double*);
static void gview_particle_trajectory(Front*,boolean);

char *in_name,*restart_state_name,*restart_name,*out_name;
boolean RestartRun;
boolean ReadFromInput;
int RestartStep;
boolean binary = NO;

int main(int argc, char **argv)
{
	static Front front;
	static F_BASIC_DATA f_basic;
	static LEVEL_FUNC_PACK level_func_pack;
	static VELO_FUNC_PACK velo_func_pack;
	IF_PARAMS iFparams;

	/* Initialize basic computational data */

	FT_Init(argc,argv,&f_basic);
	f_basic.size_of_intfc_state = sizeof(STATE);

	//Initialize Petsc
        PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);

	/*Construct Incompress Solver cartesian*/

	Incompress_Solver_Smooth_Basis *cartesian = NULL;
	if(f_basic.dim == 2)
	    cartesian = new Incompress_Solver_Smooth_2D_Cartesian(front);
	else if(f_basic.dim == 3)
	    cartesian = new Incompress_Solver_Smooth_3D_Cartesian(front);

        in_name                 = f_basic.in_name;
        restart_state_name      = f_basic.restart_state_name;
        out_name                = f_basic.out_name;
        restart_name            = f_basic.restart_name;
        RestartRun              = f_basic.RestartRun;
        RestartStep             = f_basic.RestartStep;

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

	if (debugging("trace")) printf("Passed FT_StartUp()\n");
	iFparams.dim = f_basic.dim;
	front.extra1 = (POINTER)&iFparams;

	read_iFparams(in_name,&iFparams);
	read_iF_movie_options(in_name,&iFparams);

	if (!RestartRun)
	{
	    level_func_pack.pos_component = LIQUID_COMP2;
	    FT_InitIntfc(&front,&level_func_pack);
	    read_iF_dirichlet_bdry_data(in_name,&front,f_basic);
	    if (f_basic.dim < 3)
	    	FT_ClipIntfcToSubdomain(&front);
	    initWaterDrops(&front);
	}
	else
	    restart_set_dirichlet_bdry_function(&front);

	front._compute_force_and_torque = compute_ice_particle_force;

	/* Initialize velocity field function */

	velo_func_pack.func_params = (POINTER)cartesian;
	velo_func_pack.func = rgbody_vel;
	velo_func_pack.point_propagate = ifluid_point_propagate;
	FT_InitVeloFunc(&front,&velo_func_pack);
	if (debugging("trace"))
	    printf("Passed FT_InitVeloFunc()\n");

	cartesian->initMesh();
	if (debugging("sample_velocity"))
            cartesian->initSampleVelocity(in_name);
	init_fluid_state_func(cartesian);
	if (debugging("trace"))
	    printf("Passed cartesian.initMesh()\n");
	cartesian->findStateAtCrossing = ifluid_find_state_at_crossing;

	if (RestartRun)
	    cartesian->readFrontInteriorStates(restart_state_name);
	else
	    cartesian->setInitialCondition();

	/* Propagate the front */

	fluid_driver(&front,cartesian);

	PetscFinalize();
	clean_up(0);
}

static  void fluid_driver(
        Front *front,
	Incompress_Solver_Smooth_Basis *cartesian)
{
        double CFL;

	Curve_redistribution_function(front) = full_redistribute;

	FT_ReadTimeControl(in_name,front);
	CFL = Time_step_factor(front);

	if (!RestartRun)
	{
	    FT_RedistMesh(front);
	    FT_ResetTime(front);
	    gview_particle_trajectory(front,NO);

	    if (debugging("trace"))
		printf("Calling initial FT_Propagate()\n");
	    FrontPreAdvance(front);
            FT_Propagate(front);
	    cartesian->solve(front->dt);
	    FT_SetOutputCounter(front);
	    FT_SetTimeStep(front);
        }
        else
	    FT_SetOutputCounter(front);

	FT_TimeControlFilter(front);
        printf("\ntime = %f   step = %5d   next dt = %f\n",
                        front->time,front->step,front->dt);
        fflush(stdout);

	if (debugging("trace"))
	{
	    printf("CFL = %f\n",CFL);
	    printf("Frequency_of_redistribution(front,GENERAL_WAVE) = %d\n",
			Frequency_of_redistribution(front,GENERAL_WAVE));
	}

	if (debugging("trace")) printf("Before time loop\n");
        for (;;)
        {
            /* Propagating interface for time step dt */

	    if (debugging("trace")) printf("Begin a time step\n");
	    FrontPreAdvance(front);
	    if (debugging("trace")) printf("Passed FrontPreAdvance()\n");
            FT_Propagate(front);

	    if (debugging("trace")) printf("Begin calling solve()\n");
	    cartesian->solve(front->dt);
	    if (debugging("trace")) printf("Passed solve()\n");

	    FT_AddTimeStepToCounter(front);
				
            //Next time step determined by maximum speed of previous
            //step, assuming the propagation is hyperbolic and
            //is not dependent on second order derivatives of
            //the interface such as curvature, and etc.

	    FT_SetTimeStep(front);
	    if (debugging("step_size"))
		printf("Time step from FrontHypTimeStep(): %f\n",front->dt);
            front->dt = std::min(front->dt,CFL*cartesian->max_dt);
	    if (debugging("step_size"))
                printf("Time step from l_cartesian->max_dt(): %f\n",front->dt);
	
            /* Output section */

	    gview_particle_trajectory(front,NO);
            if (FT_IsSaveTime(front))
	    {
            	FT_Save(front,out_name);
		cartesian->printFrontInteriorStates(out_name);
	    }
            if (FT_IsMovieFrameTime(front))
	    {
	    	if (debugging("trace")) 
		    printf("Calling compMovieVariables()\n");
	        cartesian->initMovieVariables();
	    	if (debugging("trace")) 
		    printf("Calling FT_AddMovieFrame()\n");
            	FT_AddMovieFrame(front,out_name,binary);
	    }

            if (FT_TimeLimitReached(front))
	    {
            	printf("\ntime = %f   step = %5d   next dt = %f\n",
                        front->time,front->step,front->dt);
            	fflush(stdout);
                break;
	    }

	    FT_TimeControlFilter(front);

            printf("\ntime = %f   step = %5d   next dt = %f\n",
                        front->time,front->step,front->dt);
            fflush(stdout);
        }
	gview_particle_trajectory(front,YES);
	if (debugging("trace")) printf("After time loop\n");
}       /* end fluid_driver */


static int rgbody_vel(
	POINTER params,
	Front *front,
	POINT *p,
	HYPER_SURF_ELEMENT *hse,
	HYPER_SURF *hs,
	double *vel)
{
	int i,dim; 

	dim  = front->rect_grid->dim;
	for (i = 0; i < dim; ++i)
	{
            vel[i] = center_of_mass_velo(hs)[i];
	}
	return YES;
}	/* end rgbody_vel */

/* 	Function computes velocity of center of mass and
 *  	angular velocity of a regid body, must be a closed curve. 
*/

extern	void compute_ice_particle_force(
	Front *fr,
	HYPER_SURF *hs,
	double dt,
	double *force,
	double *torque)
{
	switch (fr->rect_grid->dim)
	{
	case 2:
	    return compute_ice_particle_force2d(fr,hs,dt,force,torque);
	case 3:
	    return compute_ice_particle_force3d(fr,hs,dt,force,torque);
	}
}	/* end compute_ice_particle_force */

static	void compute_ice_particle_force2d(
	Front *fr,
	HYPER_SURF *hs,
	double dt,
	double *force,
	double *torque)
{
	RECT_GRID *gr = computational_grid(fr->interf);
	double f[MAXD],rr[MAXD];
	double t,pres;
	double area[MAXD],posn[MAXD];
	BOND *b;
	boolean pos_side;
	int i,dim = gr->dim;
	IF_PARAMS *iFparams = (IF_PARAMS*)fr->extra1;
	double *gravity = iFparams->gravity;
	CURVE *curve = Curve_of_hs(hs);

	if (debugging("rigid_body"))
	    (void) printf("Entering compute_ice_particle_force2d()\n");

	if (ifluid_comp(negative_component(curve)))
	    pos_side = NO;
	else 
	    pos_side = YES;

	for (i = 0; i < dim; ++i)
	{
	    force[i] = 0.0;
	}
	*torque = 0.0;
	for (b = curve->first; b != NULL; b = b->next)
	{
	    if (force_on_hse(Hyper_surf_element(b),Hyper_surf(curve),gr,
			&pres,area,posn,pos_side))
	    {
	    	for (i = 0; i < dim; ++i)
	    	{
		    f[i] = pres*area[i];
	    	    rr[i] = 0.5*(Coords(b->start)[i] + Coords(b->end)[i])
				- rotation_center(curve)[i];
	    	    force[i] += f[i];
	    	}
	    	Cross2d(rr,f,t);
	    	*torque += t;
	    }
	}
	 /* Add gravity to the total force */
	if (motion_type(curve) != ROTATION)
	{
	    for (i = 0; i < dim; ++i)
	    	force[i] += gravity[i]*total_mass(curve);
	}
	if (debugging("rigid_body"))
	{
	    (void) printf("Leaving compute_ice_particle_force2d()\n");
	    (void) printf("total_force = %f %f\n",force[0],force[1]);
	    (void) printf("torque = %f\n",*torque);
	}
}	/* end compute_ice_particle_force2d */

#define         MAX_TRI_FOR_INTEGRAL            100
static	void compute_ice_particle_force3d(
	Front *fr,
	HYPER_SURF *hs,
	double dt,
	double *force,
	double *torque)
{
	RECT_GRID *gr = computational_grid(fr->interf);
	double *L = gr->L;
	double *U = gr->U;
	double f[MAXD];
	double pres;
	double area[MAXD],posn[MAXD];
	double tri_center[MAXD];
	TRI *tri;
	POINT *p;
	boolean pos_side;
	int i,j,dim = gr->dim;
	IF_PARAMS *iFparams = (IF_PARAMS*)fr->extra1;
	double *gravity = iFparams->gravity;
	SURFACE *surface = Surface_of_hs(hs);
	boolean out_domain_tri;

	if (ifluid_comp(negative_component(surface)))
	    pos_side = NO;
	else 
	    pos_side = YES;

	for (i = 0; i < dim; ++i)
	{
	    force[i] = 0.0;
	    torque[i] = 0.0;
	}
	for (tri = first_tri(surface); !at_end_of_tri_list(tri,surface); 
			tri = tri->next)
	{
	    // Eliminate periodic duplicates
	    out_domain_tri = NO;
	    for (i = 0; i < dim; ++i)
		tri_center[i] = 0;
	    for (j = 0; j < 3; ++j)
	    {
		p = Point_of_tri(tri)[j];
	    	for (i = 0; i < dim; ++i)
		    tri_center[i] += Coords(p)[i];
	    }
	    for (i = 0; i < dim; ++i)
	    {
		tri_center[i] /= 3.0;
		if (tri_center[i] <= L[i] || tri_center[i] > U[i])
		    out_domain_tri = YES;
	    }
	    if (out_domain_tri == YES) 
	    {
		continue;
	    }

	    if (force_on_hse(Hyper_surf_element(tri),Hyper_surf(surface),gr,
			&pres,area,posn,pos_side))
	    {
	    	for (i = 0; i < dim; ++i)
	    	{
		    f[i] = pres*area[i];
	    	    force[i] += f[i];
		}
	    }
	}
	 /* Add gravity to the total force */
	if (motion_type(surface) != ROTATION)
	{
	    for (i = 0; i < dim; ++i)
	    	force[i] += gravity[i]*total_mass(surface);
	}
	if (debugging("rigid_body"))
	{
	    printf("In compute_ice_particle_force3d()\n");
	    printf("body_index = %d\n",body_index(hs));
	    printf("total_force = %f %f %f\n\n",force[0],force[1],force[2]);
	}
}	/* end compute_ice_particle_force3d */


static boolean force_on_hse(
	HYPER_SURF_ELEMENT *hse,	/* Bond (2D) or tri (3D) */
	HYPER_SURF *hs,			/* Curve (2D) or surface (3D) */
	RECT_GRID *gr,			/* Rectangular grid */
	double *pres,		/* Average pressure */
	double *area,		/* Area as a vector, pointing onto body */
	double *posn,		/* Position of the pressure */
	boolean pos_side)	/* Is the body on the positive side of hs? */
{
	int dim = gr->dim;
	switch (dim)
	{
	case 2: 
	    return force_on_hse2d(hse,hs,gr,pres,area,posn,pos_side);
	case 3: 
	    return force_on_hse3d(hse,hs,gr,pres,area,posn,pos_side);
	default:
	    return NO;
	}
	
}	/* end force_on_hse */

static boolean force_on_hse2d(
	HYPER_SURF_ELEMENT *hse,
	HYPER_SURF *hs,
	RECT_GRID *gr,
	double *pres,
	double *area,
	double *posn,
	boolean pos_side)
{
	double crds1[MAXD],crds2[MAXD];
	double p1,p2;
	Locstate s1,s2;
	BOND *b = Bond_of_hse(hse);
	CURVE *c = Curve_of_hs(hs);
	double *L = gr->L;
	double *U = gr->U;
	int i;
	
	/* Get pressure at two end points of the bond */
	if (b->start == c->start->posn)
	    s1 = pos_side ? right_start_state(c) : left_start_state(c);
	else
	    s1 = pos_side ? right_state(b->start) : left_state(b->start);
	if (b->end == c->end->posn)
	    s2 = pos_side ? right_end_state(c) : left_end_state(c);
	else
	    s2 = pos_side ? right_state(b->end) : left_state(b->end);

	p1 = getStatePres(s1);	p2 = getStatePres(s2);
	for (i = 0; i < 2; ++i)
	{
	    crds1[i] = Coords(b->start)[i];
	    crds2[i] = Coords(b->end)[i];
	}

	/* Cut and interpolate if one end is outside the domain */
	for (i = 0; i < 2; ++i)
	{
	    if (crds1[i] <= L[i])
	    {
		if (crds2[i] <= L[i]) return NO; // both ends out
		else
		{
		    crds1[(i+1)%2] = intrp_between(crds1[i],crds2[i],L[i],
				crds1[(i+1)%2],crds2[(i+1)%2]);
		    p1 = intrp_between(crds1[i],crds2[i],L[i],p1,p2);
		    crds1[i] = L[i];
		}
	    }
	    if (crds1[i] >= U[i])
	    {
		if (crds2[i] >= U[i]) return NO; // both ends out
		else
		{
		    crds1[(i+1)%2] = intrp_between(crds1[i],crds2[i],U[i],
				crds1[(i+1)%2],crds2[(i+1)%2]);
		    p1 = intrp_between(crds1[i],crds2[i],U[i],p1,p2);
		    crds1[i] = U[i];
		}
	    }
	}
	for (i = 0; i < 2; ++i)
	{
	    if (crds2[i] <= L[i])
	    {
		if (crds1[i] <= L[i]) return NO; // both ends out
		else
		{
		    crds2[(i+1)%2] = intrp_between(crds1[i],crds2[i],L[i],
				crds1[(i+1)%2],crds2[(i+1)%2]);
		    p2 = intrp_between(crds1[i],crds2[i],L[i],p1,p2);
		    crds2[i] = L[i];
		}
	    }
	    if (crds2[i] >= U[i])
	    {
		if (crds1[i] >= U[i]) return NO; // both ends out
		else
		{
		    crds2[(i+1)%2] = intrp_between(crds1[i],crds2[i],U[i],
				crds1[(i+1)%2],crds2[(i+1)%2]);
		    p2 = intrp_between(crds1[i],crds2[i],U[i],p1,p2);
		    crds2[i] = U[i];
		}
	    }
	}
	area[0] = pos_side ? crds1[1] - crds2[1] : crds2[1] - crds1[1];
	area[1] = pos_side ? crds2[0] - crds1[0] : crds1[0] - crds2[0];
	*pres = 0.5*(p1 + p2);
	posn[0] = 0.5*(crds1[0] + crds2[0]);
	posn[1] = 0.5*(crds1[1] + crds2[1]);
	return YES;
}	/* end force_on_hse2d */

static boolean force_on_hse3d(
	HYPER_SURF_ELEMENT *hse,
	HYPER_SURF *hs,
	RECT_GRID *gr,
	double *pres,
	double *area,
	double *posn,
	boolean pos_side)
{
        TRI *t = Tri_of_hse(hse);
	POINT *point;
	Locstate sl,sr;
        int i,j,dim = gr->dim;

	*pres = 0.0;
	for (i = 0; i < 3; ++i)
	    posn[i] = 0.0;
	for (i = 0; i < 3; ++i)
	{
	    point = Point_of_tri(t)[i];
	    for (j = 0; j < dim; ++j)
		posn[j] += Coords(point)[j];
	    FT_GetStatesAtPoint(point,hse,hs,&sl,&sr);
	    if (pos_side)
		*pres += getStatePres(sr);
	    else
		*pres += getStatePres(sl);
	}
	*pres /= 3.0;
	for (i = 0; i < dim; ++i)
	{
	    area[i] = pos_side ? -Tri_normal(t)[i] : Tri_normal(t)[i];
	    posn[i] /= 3.0;
	}
	/* Need to treat subdomain boundary */
	return YES;
}	/* end force_on_hse3d */

static double intrp_between(
	double x1,
	double x2,
	double x,
	double y1,
	double y2)
{
	double y;
	if (x1 == x2) return y1;
	y = y1 + (y2 - y1)/(x2 - x1)*(x - x1);
	return y;
}

#define		MAX_TIME_STEPS		5000

static void gview_particle_trajectory(
	Front *front,
	boolean end_of_run)
{
	INTERFACE *intfc = front->interf;
	SURFACE **s;
	static int num_particles = 0;
	static int *pindex,max_index,min_index;
	static boolean first = YES;
	static double ***com;
	int step = front->step;
	static int num_steps;
	int i,ip,dim = Dimension(intfc);

	if (first == YES)
	{
	    max_index = 0;
	    min_index = 10000;
	    intfc_surface_loop(intfc,s)
	    {
		if (wave_type(*s) == ICE_PARTICLE_BOUNDARY)
		{
		    if (max_index < body_index(*s))
			max_index = body_index(*s);
		    if (min_index > body_index(*s))
			min_index = body_index(*s);
		}
	    }
	    pp_global_imax(&max_index,1);
	    pp_global_imin(&min_index,1);
	    num_particles = max_index - min_index + 1;
	    tri_array(&com,MAX_TIME_STEPS,num_particles,MAXD,sizeof(double));
	    uni_array(&pindex,num_particles,sizeof(int));
	    num_steps = MAX_TIME_STEPS;
	}	
	if (step >= num_steps)
	{
	    FT_FreeThese(1,com);
	    num_steps += MAX_TIME_STEPS;
	    tri_array(&com,num_steps,num_particles,MAXD,sizeof(double));
	}
	intfc_surface_loop(intfc,s)
	{
	    if (wave_type(*s) == ICE_PARTICLE_BOUNDARY)
	    {
		ip = body_index(*s);
		for (i = 0; i < dim; ++i)
		    com[step][ip][i] = center_of_mass(*s)[i];
	    }
	}
}	/* end gview_particle_trajectory */

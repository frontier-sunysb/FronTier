#include <iFluid.h>
#include <rgbody.h>

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
static void ifluid_compute_force_and_torque2d(Front*,HYPER_SURF*,double,
			double*,double*);
static void ifluid_compute_force_and_torque3d(Front*,HYPER_SURF*,double,
			double*,double*);

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
	RG_PROB_TYPE prob_type;

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

	read_rg_prob_type(in_name,&prob_type);
	read_iFparams(in_name,&iFparams);

	if (!RestartRun)
	{
	    level_func_pack.neg_component = SOLID_COMP;
	    level_func_pack.pos_component = LIQUID_COMP2;
	    level_func_pack.func = NULL;
	    level_func_pack.point_array = NULL;
	    level_func_pack.func_params = iFparams.level_func_params;
	    if (f_basic.dim == 3)
	    	level_func_pack.set_3d_bdry = YES;

	    init_moving_bodies(&front,&level_func_pack,in_name,prob_type);
	    FT_InitIntfc(&front,&level_func_pack);
	    if (debugging("trace"))
	    	printf("Passed FT_InitIntfc()\n");
	    read_iF_dirichlet_bdry_data(in_name,&front,f_basic);
	    if (f_basic.dim < 3)
	    	FT_ClipIntfcToSubdomain(&front);
	}
	else
	    restart_set_dirichlet_bdry_function(&front);

	front._compute_force_and_torque = ifluid_compute_force_and_torque;

	/* Initialize velocity field function */

	velo_func_pack.func_params = (POINTER)cartesian;
	velo_func_pack.func = rgbody_vel;
	velo_func_pack.point_propagate = ifluid_point_propagate;
	FT_InitVeloFunc(&front,&velo_func_pack);
	if (debugging("trace"))
	    printf("Passed FT_InitVeloFunc()\n");

	cartesian->initMesh();
	cartesian->initMovieVariables();
	if (debugging("sample_velocity"))
            cartesian->initSampleVelocity(in_name);
	init_fluid_state_func(cartesian,prob_type);
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

	    if (debugging("trace"))
		printf("Calling initial FT_Propagate()\n");
	    FrontPreAdvance(front);
            FT_Propagate(front);
	    //cartesian->solve(front->dt);
	    record_moving_body_data(out_name,front);
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

	    record_moving_body_data(out_name,front);

            if (FT_IsSaveTime(front))
	    {
            	FT_Save(front,out_name);
		cartesian->printFrontInteriorStates(out_name);
	    }
            if (FT_IsMovieFrameTime(front))
	    {
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
	double omega,crds_at_com[MAXD];

	dim  = front->rect_grid->dim;
	for (i = 0; i < dim; ++i)
	{
	    crds_at_com[i] = Coords(p)[i] - rotation_center(hs)[i];
            vel[i] = center_of_mass_velo(hs)[i];
	}
	switch(dim)
	{
	case 2:
	    omega = angular_velo(hs);
	    vel[0] += omega*crds_at_com[1];
	    vel[1] += -omega*crds_at_com[0];
	    break;
	case 3:
	    omega = angular_velo(hs);
            vel[0] += omega*crds_at_com[2];
            vel[2] += -omega*crds_at_com[0];
	}
	return YES;
}	/* end rgbody_vel */

/* 	Function computes velocity of center of mass and
 *  	angular velocity of a regid body, must be a closed curve. 
*/

extern	void ifluid_compute_force_and_torque(
	Front *fr,
	HYPER_SURF *hs,
	double dt,
	double *force,
	double *torque)
{
	switch (fr->rect_grid->dim)
	{
	case 2:
	    return ifluid_compute_force_and_torque2d(fr,hs,dt,force,torque);
	case 3:
	    return ifluid_compute_force_and_torque3d(fr,hs,dt,force,torque);
	}
}	/* end ifluid_compute_force_and_torque */

static	void ifluid_compute_force_and_torque2d(
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
	    (void) printf("Entering ifluid_compute_force_and_torque2d()\n");

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
	    (void) printf("Leaving ifluid_compute_force_and_torque2d()\n");
	    (void) printf("total_force = %f %f\n",force[0],force[1]);
	    (void) printf("torque = %f\n",*torque);
	}
}	/* end ifluid_compute_force_and_torque2d */

#define         MAX_TRI_FOR_INTEGRAL            100
static	void ifluid_compute_force_and_torque3d(
	Front *fr,
	HYPER_SURF *hs,
	double dt,
	double *force,
	double *torque)
{
	RECT_GRID *gr = computational_grid(fr->interf);
	double f[MAXD],rr[MAXD];
	double t[MAXD],tdir,pres;
	double area[MAXD],posn[MAXD];
	TRI *tri;
	boolean pos_side;
	int i,dim = gr->dim;
	IF_PARAMS *iFparams = (IF_PARAMS*)fr->extra1;
	double *gravity = iFparams->gravity;
	SURFACE *surface = Surface_of_hs(hs);

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
	    if (force_on_hse(Hyper_surf_element(tri),Hyper_surf(surface),gr,
			&pres,area,posn,pos_side))
	    {
	    	for (i = 0; i < dim; ++i)
	    	{
		    f[i] = pres*area[i];
	    	    force[i] += f[i];
		    rr[i] = posn[i] - rotation_center(surface)[i];
		}
		Cross3d(rr,f,t);
		tdir = Dot3d(t,(rotation_direction(hs)));
	    	for (i = 0; i < dim; ++i)
		{
		    t[i] = tdir*rotation_direction(hs)[i];
		    torque[i] += t[i];
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
	    printf("In ifluid_compute_force_and_torque3d()\n");
	    printf("total_force = %f %f %f\n",force[0],force[1],force[2]);
	    printf("torque = %f %f %f\n",torque[0],torque[1],torque[2]);
	}
}	/* end ifluid_compute_force_and_torque3d */


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

extern void record_moving_body_data(
	char *out_name,
	Front *front)
{
	int i,j,num_moving_body;
        INTERFACE *intfc = front->interf;
        CURVE **c;
        SURFACE **s;
	static boolean first = YES;
	static FILE **torque_files,**omega_files;
	static FILE **force_files,**com_files;
	static double *torque,*omega;
	static double **force,**com_velo;
	int dim = intfc->dim;
	char fname[256];

	num_moving_body = -1;
	switch (dim)
	{
	case 2:
            for (c = intfc->curves; c && *c; ++c)
            {
            	if (wave_type(*c) == MOVABLE_BODY_BOUNDARY)
                    if (body_index(*c) > num_moving_body)
                    	num_moving_body = body_index(*c);
            }
	    break;
	case 3:
            for (s = intfc->surfaces; s && *s; ++s)
            {
            	if (wave_type(*s) == MOVABLE_BODY_BOUNDARY)
                    if (body_index(*s) > num_moving_body)
                    	num_moving_body = body_index(*s);
            }
	    break;
	}
        num_moving_body++;
        pp_global_imax(&num_moving_body,1);
        if (num_moving_body == 0) return;

        if (first)
        {
            FT_VectorMemoryAlloc((POINTER*)&torque,num_moving_body,FLOAT);
            FT_VectorMemoryAlloc((POINTER*)&omega,num_moving_body,FLOAT);
            FT_MatrixMemoryAlloc((POINTER*)&force,num_moving_body,MAXD,FLOAT);
            FT_MatrixMemoryAlloc((POINTER*)&com_velo,num_moving_body,
					MAXD,FLOAT);
	    if (pp_mynode() == 0)
	    {
            	FT_VectorMemoryAlloc((POINTER*)&torque_files,num_moving_body,
					sizeof(FILE*));
            	FT_VectorMemoryAlloc((POINTER*)&omega_files,num_moving_body,
					sizeof(FILE*));
            	FT_VectorMemoryAlloc((POINTER*)&force_files,num_moving_body,
					sizeof(FILE*));
            	FT_VectorMemoryAlloc((POINTER*)&com_files,num_moving_body,
					sizeof(FILE*));
            	for (i = 0; i < num_moving_body; ++i)
		{
		    sprintf(fname,"%s/cen-of_mass-%d",out_name,i);
		    com_files[i] = fopen(fname,"w");
		    sprintf(fname,"%s/omega-%d",out_name,i);
		    omega_files[i] = fopen(fname,"w");
		    sprintf(fname,"%s/force-%d",out_name,i);
		    force_files[i] = fopen(fname,"w");
		    sprintf(fname,"%s/torque-%d",out_name,i);
		    torque_files[i] = fopen(fname,"w");
		}
	    }
        }
	switch (dim)
	{
	case 2:
	    for (c = intfc->curves; c && *c; ++c)
            {
            	if (wave_type(*c) == MOVABLE_BODY_BOUNDARY)
            	{
                    i = body_index(*c);
                    FrontForceAndTorqueOnHs(front,Hyper_surf(*c),front->dt,
					force[i],&torque[i]);
                    omega[i] = angular_velo(*c);
                    for (j = 0; j < dim; ++j)
                    	com_velo[i][j] = center_of_mass_velo(*c)[j];
            	}
            }
	    break;
	case 3:
	    for (s = intfc->surfaces; s && *s; ++s)
            {
            	if (wave_type(*s) == MOVABLE_BODY_BOUNDARY)
            	{
                    i = body_index(*s);
                    FrontForceAndTorqueOnHs(front,Hyper_surf(*s),front->dt,
					force[i],&torque[i]);
                    omega[i] = angular_velo(*s);
                    for (j = 0; j < dim; ++j)
                    	com_velo[i][j] = center_of_mass_velo(*s)[j];
            	}
            }
	    break;
	}
        for (i = 0; i < num_moving_body; ++i)
        {
            pp_global_sum(force[i],dim);
            pp_global_sum(&torque[i],1);
        }
        if (pp_mynode() != 0) return;

	if (first)
        {
            first = NO;
            for (i = 0; i < num_moving_body; ++i)
            {
                fprintf(torque_files[i],"\"Torque of body %d\"\n",i+1);
                fprintf(force_files[i],"\"Total force on body %d\"\n",i+1);
                fprintf(omega_files[i],"\"Angular velocity of body %d\"\n",i+1);                fprintf(com_files[i],"\"COM velocity of body %d\"\n",i+1);
            }
        }
        for (i = 0; i < num_moving_body; ++i)
        {
            fprintf(torque_files[i],"%f  %f\n",front->time,torque[i]);
            fprintf(omega_files[i],"%f  %f\n",front->time,omega[i]);
            fprintf(force_files[i],"%f  ",front->time);
            fprintf(com_files[i],"%f  ",front->time);
            for (j = 0; j < dim; ++j)
            {
                fprintf(force_files[i],"%f  ",force[i][j]);
                fprintf(com_files[i],"%f  ",com_velo[i][j]);
            }
            fprintf(force_files[i],"\n");
            fprintf(com_files[i],"\n");

            fflush(torque_files[i]);
            fflush(omega_files[i]);
            fflush(force_files[i]);
            fflush(com_files[i]);
        }
}	/* end record_moving_body_data */

extern void read_rg_prob_type(
	char *inname,
	RG_PROB_TYPE *prob_type)
{
	char string[100];
	FILE *infile = fopen(inname,"r");

	*prob_type = ERROR_TYPE;
	CursorAfterString(infile,"Enter problem type:");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	if (string[0] == 'F' || string[0] == 'f')
	{
            if (string[6] == 'S' || string[6] == 's')
                *prob_type = FLUID_SOLID_CIRCLE;
            else if (string[6] == 'R' || string[6] == 'r')
                *prob_type = FLUID_RIGID_BODY;
	}
        else if (string[0] == 'R' || string[0] == 'r')
	{
            if (string[6] == 'O' || string[6] == 'o')
                *prob_type = ROTOR_ONE_FLUID;
            else if (string[6] == 'T' || string[6] == 't')
                *prob_type = ROTOR_TWO_FLUID;
	}
        else if (string[0] == 'W' || string[0] == 'w')
	{
            if (string[9] == '2')
                *prob_type = WINDMILL_2D;
            else if (string[9] == '3')
                *prob_type = WINDMILL_3D;
	}
        else if (string[0] == 'B' || string[0] == 'b')
	{
            *prob_type = BEE_3D;
	}
        else if (string[0] == 'H' || string[0] == 'h')
	{
            *prob_type = HELICOPTER_3D;
	}

	assert(*prob_type != ERROR_TYPE);
	fclose(infile);
}	/* end read_rg_prob_type */

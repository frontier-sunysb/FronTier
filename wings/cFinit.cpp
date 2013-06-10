#include <cFluid.h>

static void initWing(Front*);

extern void initModules(
	Front *front)
{
	char *inname = InName(front);
	FILE *infile = fopen(inname,"r");
	int dim = Dimension(front->interf);
	char string[100];

	CursorAfterString(infile,"Enter module type:");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	fclose(infile);
	switch (dim)
	{
	case 2:
	    switch (string[0])
	    {
	    case 'w':
	    case 'W':
		initWing(front);
		break;
	    default:
		(void) printf("Unknown 2D module type: %s\n",string);
	    }
	    break;
	case 3:
	    switch (string[0])
	    {
	    default:
		(void) printf("Unknown 3D module type: %s\n",string);
	    }
	default:
	    (void) printf("Unknown dimension: %d\n",dim);
	    clean_up(ERROR);
	}
	if (debugging("init_intfc"))
	    gview_plot_interface("init_intfc",front->interf);
}	/* end initModules */

static void initWing(
	Front *front)
{
	
	char *inname = InName(front);
	FILE *infile = fopen(inname,"r");
	char string[100];
	FILE *wfile;
	int i,num_points;
	double **point_array;
	int w_type;
	CURVE *curve;
	double com[MAXD];	// Center of mass
	double alpha;		// Pitch angle
	POINT *p;
	BOND *b;

	if (debugging("trace"))
	    (void) printf("Entering initWing()\n");
	CursorAfterString(infile,"Enter input file name for wing:");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	wfile = fopen(string,"r");
	fscanf(wfile,"%d",&num_points);
	FT_MatrixMemoryAlloc((POINTER*)&point_array,num_points,MAXD,
				sizeof(double));
	for (i = 0; i < num_points; ++i)
	{
	    fscanf(wfile,"%lf %lf %lf",&point_array[i][0],&point_array[i][1],
					&point_array[i][2]);
	}
	if (point_array[0][0] == point_array[num_points-1][0] &&
	    point_array[0][1] == point_array[num_points-1][1] &&
	    point_array[0][2] == point_array[num_points-1][2]) 
	    num_points--;
	w_type = NEUMANN_BOUNDARY;
	if (CursorAfterStringOpt(infile,"Enter yes to set wing in motion:"))
	{
	    fscanf(infile,"%s",string);
	    (void) printf("%s\n",string);
	    if (string[0] == 'y' || string[0] == 'Y')
		w_type = MOVABLE_BODY_BOUNDARY;
	}
	curve = FT_MakePointArrayCurve(front,num_points,point_array,
		exterior_component(front->interf),GAS_COMP1,YES,w_type);
	if (w_type == MOVABLE_BODY_BOUNDARY)
	{
	    motion_type(curve) = PRESET_MOTION;
	}
	CursorAfterString(infile,"Enter center of mass:");
	fscanf(infile,"%lf %lf",&com[0],&com[1]);
	(void) printf("%f %f\n",com[0],com[1]);
	for (i = 0; i < 2; ++i)
	    center_of_mass(Hyper_surf(curve))[i] = com[i];
	CursorAfterString(infile,"Enter pitch angle in degree:");
	fscanf(infile,"%lff",&alpha);
	(void) printf("%f\n",alpha);
	alpha *= -2.0*PI/360.0;
	p = curve->first->start;
	rotate_point_with_angle(p,com,alpha,YES);
	for (b = curve->first; b != NULL; b = b->next)
	{
	    if (b == curve->last && is_closed_curve(curve))
		continue;
	    p = b->end;
	    rotate_point_with_angle(p,com,alpha,NO);
	}

	if (debugging("trace"))
	{
	    (void) printf("Number of point on the wing: %d\n",num_points);
	    (void) printf("Interface after initWing()\n");
	    print_interface(front->interf);
	    (void) printf("Leaving initWing()\n");
	}
	fclose(infile);
	fclose(wfile);
}	/* end initWing */

extern void initStateParams(
	Front *front)
{
	char *inname = InName(front);
	EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
	char str[200];
	double pinf,einf,gamma;

	FILE *infile = fopen(inname,"r");
	sprintf(str, "Enter gamma, pinf, einf of the fluid with comp %d:",
                     GAS_COMP1);
        CursorAfterString(infile,str);
        fscanf(infile,"%lf %lf %lf",&gamma,&pinf,&einf);
        (void) printf("%f %f %f\n",gamma,pinf,einf);
        (eqn_params->eos[GAS_COMP1]).gamma = gamma;
        (eqn_params->eos[GAS_COMP1]).pinf = pinf;
        (eqn_params->eos[GAS_COMP1]).einf = einf;
	sprintf(str, "Enter density of fluid with comp %d:",GAS_COMP1);
	CursorAfterString(infile,str);
        fscanf(infile,"%lf",&eqn_params->rho1);
        (void) printf("%f\n",eqn_params->rho1);
	sprintf(str, "Enter pressure of fluid with comp %d:",GAS_COMP1);
	CursorAfterString(infile,str);
        fscanf(infile,"%lf",&eqn_params->p0);
        (void) printf("%f\n",eqn_params->p0);
	sprintf(str, "Enter velocity of fluid with comp %d:",GAS_COMP1);
	CursorAfterString(infile,str);
        fscanf(infile,"%lf %lf",&eqn_params->v1[0],&eqn_params->v1[1]);
        (void) printf("%f %f\n",eqn_params->v1[0],eqn_params->v1[1]);

}	/* end initStateParams */

void G_CARTESIAN::initFrontInteriorStates()
{
	int i,j,k,l,index;
        EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
        double coords[MAXD];
        COMPONENT comp;
        STATE *sl,*sr,state;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
        INTERFACE *intfc = front->interf;
        double *dens = field.dens;
        double *engy = field.engy;
        double *pres = field.pres;
        double **momn = field.momn;
	EOS_PARAMS      *eos;

	eos = &(eqn_params->eos[GAS_COMP1]);
	next_point(intfc,NULL,NULL,NULL);
	while (next_point(intfc,&p,&hse,&hs))
        {
            FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
	    if (negative_component(hs) == GAS_COMP1)
	    {
		sl->eos = eos;
		sl->dens = eqn_params->rho1;
		sl->pres = eqn_params->p0;
		sl->momn[0] = eqn_params->rho1*eqn_params->v1[0];
		sl->momn[1] = eqn_params->rho1*eqn_params->v1[1];
		sl->engy = EosInternalEnergy(sl);
	    }
	    else if (positive_component(hs) == GAS_COMP1)
	    {
		sr->eos = eos;
		sr->dens = eqn_params->rho1;
		sr->pres = eqn_params->p0;
		sr->momn[0] = eqn_params->rho1*eqn_params->v1[0];
		sr->momn[1] = eqn_params->rho1*eqn_params->v1[1];
		sr->engy = EosInternalEnergy(sr);
	    }
        }
        FT_MakeGridIntfc(front);
	setDomain();
	state.eos = eos;
	switch (dim)
	{
	case 2:
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index2d(i,j,top_gmax);
		comp = top_comp[index];
		if (comp != GAS_COMP1) continue;
		state.momn[0] = momn[0][index] = 
				eqn_params->rho1*eqn_params->v1[0];
		state.momn[1] = momn[1][index] = 
				eqn_params->rho1*eqn_params->v1[1];
		state.dens = dens[index] = eqn_params->rho1;
		state.pres = pres[index] = eqn_params->p0;
		engy[index] = EosInternalEnergy(&state);
	    }
	    break;
	case 3:
	    for (k = imin[2]; k <= imax[2]; ++k)
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index3d(i,j,k,top_gmax);
		comp = top_comp[index];
		if (comp != GAS_COMP1) continue;
		for (l = 0; l < dim; ++l)
		    momn[l][index] = 0.0;
		state.dens = dens[index] = eqn_params->rho1;
		state.pres = pres[index] = eqn_params->p0;
		engy[index] = EosInternalEnergy(&state);
	    }
	    break;
	}
	scatMeshStates();
}	/* end initInteriorStates */
	
extern void initMotionParams(
	Front *front)
{
	CURVE *moving_curve = NULL;
	char *inname = InName(front);
	char str[200];
	FILE *infile;
	static MOTION_PARAMS motion_params;
	CURVE **c;
	int i;

	for (c = front->interf->curves; c && *c; ++c)
	{
	    if (wave_type(*c) == MOVABLE_BODY_BOUNDARY)
	    {
		moving_curve = *c;
	    }
	}
	if (moving_curve == NULL) return;

	infile = fopen(inname,"r");
	(void) printf("%f\n",motion_params.amplitude);
	(void) printf("Motion of wing application are preset\n");
	(void) printf("\tMotion types are\n");
	(void) printf("\tOscillatory (O)\n");
	(void) printf("\tRotation (R)\n");
	CursorAfterString(infile,"Enter angular motion type:");
	fscanf(infile,"%s",str);
	(void) printf("%s\n",str);
	switch (str[0])
	{
	case 'O':
	case 'o':
	    motion_params.wing_motion_type = OSCILLATORY_MOTION;
	    CursorAfterString(infile,"Enter amplitude of angular velocity:");
	    fscanf(infile,"%lf",&motion_params.amplitude);
	    (void) printf("%f\n",motion_params.amplitude);
	    CursorAfterString(infile,"Enter period of angular velocity:");
	    fscanf(infile,"%lf",&motion_params.period);
	    (void) printf("%f\n",motion_params.period);
	    break;
	case 'R':
	case 'r':
	    motion_params.wing_motion_type = ROTATIONAL_MOTION;
	    CursorAfterString(infile,"Enter angular velocity:");
	    fscanf(infile,"%lf",&motion_params.angular_velo);
	    (void) printf("%f\n",motion_params.angular_velo);
	    break;
	default:
	    (void) printf("Unknown angular motion!\n");
	    clean_up(ERROR);
	}
	for (c = front->interf->curves; c && *c; ++c)
	{
	    if (wave_type(*c) == MOVABLE_BODY_BOUNDARY)
	    {
		(*c)->extra = (POINTER)&motion_params;
	    }
	}
	fclose(infile);
}	/* end initMotionParams */
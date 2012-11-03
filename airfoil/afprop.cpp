#include <iFluid.h>
#include <airfoil.h>

static double (*getStateVel[3])(POINTER) = {getStateXvel,getStateYvel,
                                        getStateZvel};
static double (*getStateVort3d[3])(POINTER) = {getStateXvort,getStateYvort,
                                        getStateZvort};
static SURFACE *canopy_of_string_node(NODE*);
static void convert_to_point_mass(Front*,AF_PARAMS*);
static void airfoil_curve_propagate(Front*,POINTER,CURVE*,CURVE*,double);
static void curve_prop_without_interaction(Front*,POINTER,CURVE*,CURVE*,double);
static void curve_prop_with_interaction(Front*,POINTER,CURVE*,CURVE*,double);
static void curve_point_prop_with_interaction(Front*,POINT*,POINT*,BOND*,
					double);
static	int arrayOfMonoHsbdry(INTERFACE*,CURVE**);
static	int arrayOfGoreHsbdry(INTERFACE*,CURVE**);
static 	int getGoreNodes(INTERFACE*,NODE**);


typedef struct {
	double cen[MAXD];
	double v0;
	double stop_time;
} VERTICAL_PARAMS;


extern void initVelocityFunc(
	char *inname,
	Front *front)
{
	static VELO_FUNC_PACK velo_func_pack;
	static VORTEX_PARAMS *vortex_params; /* velocity function parameters */
        static BIPOLAR_PARAMS *dv_params;
	static VERTICAL_PARAMS *vert_params;
	static TOROIDAL_PARAMS *toro_params;
	static PARABOLIC_PARAMS *para_params;
	static SINGULAR_PARAMS *sing_params;
	FILE *infile = fopen(inname,"r");
	int i,dim = front->rect_grid->dim;
	char string[100];
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;

	af_params->no_fluid = NO;
        if (CursorAfterStringOpt(infile,
            "Entering yes to turn off fluid solver: "))
        {
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'y' || string[0] == 'Y')
                af_params->no_fluid = YES;
        }
	if (af_params->no_fluid == YES)
	{
	    front->curve_propagate = airfoil_curve_propagate;
	    velo_func_pack.point_propagate = airfoil_point_propagate;
            CursorAfterString(infile,"Enter velocity function: ");
	    fscanf(infile,"%s",string);
	    (void) printf("%s\n",string);
	    switch (string[0])
            {
            case 'r':
            case 'R':
	    	FT_ScalarMemoryAlloc((POINTER*)&vortex_params,
				sizeof(VORTEX_PARAMS));
            	front->max_time = 0.4;
            	front->movie_frame_interval = 0.02;
            	vortex_params->dim = 2;
            	vortex_params->type[0] = 'M';
            	vortex_params->cos_time = 0;
            	vortex_params->cen[0] = 0.5;
            	vortex_params->cen[1] = 0.25;
            	vortex_params->rad = 0.15;
            	vortex_params->time = 0.5*front->max_time;
            	velo_func_pack.func_params = (POINTER)vortex_params;
            	velo_func_pack.func = vortex_vel;
            	break;
            case 'd':
            case 'D':
	    	FT_ScalarMemoryAlloc((POINTER*)&dv_params,
				sizeof(BIPOLAR_PARAMS));
            	dv_params->cen1[0] = 0.25;
            	dv_params->cen1[1] = 0.25;
            	dv_params->cen2[0] = 0.75;
            	dv_params->cen2[1] = 0.25;
            	dv_params->i1 = -0.5;
            	dv_params->i2 =  0.5;
            	velo_func_pack.func_params = (POINTER)dv_params;
            	velo_func_pack.func = double_vortex_vel;
            	break;
            case 'v':
            case 'V':
	    	FT_ScalarMemoryAlloc((POINTER*)&vert_params,
				sizeof(VERTICAL_PARAMS));
		CursorAfterString(infile,"Enter center velocity:");
        	fscanf(infile,"%lf",&vert_params->v0);
        	(void) printf("%f\n",vert_params->v0);
		CursorAfterString(infile,"Enter stop motion time:");
        	fscanf(infile,"%lf",&vert_params->stop_time);
        	(void) printf("%f\n",vert_params->stop_time);
		CursorAfterString(infile,"Enter center of vertical motion:");
        	fscanf(infile,"%lf %lf",&vert_params->cen[0],
				&vert_params->cen[1]);
            	velo_func_pack.func_params = (POINTER)vert_params;
            	velo_func_pack.func = vertical_velo;
            	break;
            case 't':
            case 'T':
	    	FT_ScalarMemoryAlloc((POINTER*)&toro_params,
				sizeof(TOROIDAL_PARAMS));
		CursorAfterString(infile,"Enter center of toroidal motion:");
        	fscanf(infile,"%lf %lf %lf",&toro_params->tcen[0],
				&toro_params->tcen[1],&toro_params->tcen[2]);
        	(void) printf("%f %f %f\n",toro_params->tcen[0],
				toro_params->tcen[1],toro_params->tcen[2]);
		CursorAfterString(infile,"Enter distance to poloidal center:");
        	fscanf(infile,"%lf",&toro_params->R0);
        	(void) printf("%f\n",toro_params->R0);
		CursorAfterString(infile,"Enter velocity magnitude:");
        	fscanf(infile,"%lf",&toro_params->v0);
        	(void) printf("%f\n",toro_params->v0);
		CursorAfterString(infile,"Enter stop motion time:");
        	fscanf(infile,"%lf",&toro_params->stop_time);
        	(void) printf("%f\n",toro_params->stop_time);
            	velo_func_pack.func_params = (POINTER)toro_params;
            	velo_func_pack.func = toroidal_velo;
            	break;
            case 'p':
            case 'P':
	    	FT_ScalarMemoryAlloc((POINTER*)&para_params,
				sizeof(PARABOLIC_PARAMS));
		CursorAfterString(infile,"Enter center of parabolic velocity:");
        	fscanf(infile,"%lf %lf",&para_params->cen[0],
				&para_params->cen[1]);
        	(void) printf("%f %f\n",para_params->cen[0],
				para_params->cen[1]);
		CursorAfterString(infile,"Enter center velocity:");
        	fscanf(infile,"%lf",&para_params->v0);
        	(void) printf("%f\n",para_params->v0);
		CursorAfterString(infile,"Enter downward concavity:");
        	fscanf(infile,"%lf",&para_params->a);
        	(void) printf("%f\n",para_params->a);
		CursorAfterString(infile,"Enter stop motion time:");
        	fscanf(infile,"%lf",&para_params->stop_time);
        	(void) printf("%f\n",para_params->stop_time);
            	velo_func_pack.func_params = (POINTER)para_params;
            	velo_func_pack.func = parabolic_velo;
            	break;
            case 's':
            case 'S':
	    	FT_ScalarMemoryAlloc((POINTER*)&sing_params,
				sizeof(SINGULAR_PARAMS));
		CursorAfterString(infile,"Enter center of velocity:");
        	fscanf(infile,"%lf %lf",&sing_params->cen[0],
				&sing_params->cen[1]);
        	(void) printf("%f %f\n",sing_params->cen[0],
				sing_params->cen[1]);
		CursorAfterString(infile,"Enter center velocity:");
        	fscanf(infile,"%lf",&sing_params->v0);
        	(void) printf("%f\n",sing_params->v0);
		CursorAfterString(infile,"Enter radius of center:");
        	fscanf(infile,"%lf",&sing_params->R);
        	(void) printf("%f\n",sing_params->R);
		CursorAfterString(infile,"Enter stop motion time:");
        	fscanf(infile,"%lf",&sing_params->stop_time);
        	(void) printf("%f\n",sing_params->stop_time);
            	velo_func_pack.func_params = (POINTER)sing_params;
            	velo_func_pack.func = singular_velo;
            	break;
            case 'z':
            case 'Z':
            	velo_func_pack.func_params = NULL;
            	velo_func_pack.func = zero_velo;
            	break;
	    default:
		(void) printf("Unknown velocity function, use zero_velo()\n");
            	velo_func_pack.func_params = NULL;
            	velo_func_pack.func = zero_velo;
		break;
            }	
	}
	else
	{
	    if (dim == 3)
	    	front->curve_propagate = airfoil_curve_propagate;
	    velo_func_pack.point_propagate = airfoil_point_propagate;
            velo_func_pack.func_params = NULL;
            velo_func_pack.func = airfoil_velo;
	}
	FT_InitVeloFunc(front,&velo_func_pack);
	if (af_params->no_fluid == YES || af_params->is_parachute_system == NO)
	{
            CursorAfterString(infile,"Enter tangential propagator:");
	    fscanf(infile,"%s",string);
	    (void) printf("%s\n",string);
	    if (dim == 2)
	    {
	    	switch (string[0])
	    	{
	    	case 'n':
	    	case 'N':
	    	    front->tan_curve_propagate = NULL;
	    	    break;
	    	case 'f':
	    	case 'F':
	    	    front->tan_curve_propagate 
				= fixed_length_tan_curve_propagate;
	    	    break;
	    	case 'e':
	    	case 'E':
	    	    if (string[1] == '2')
	    	    	front->tan_curve_propagate 
				= second_order_elastic_curve_propagate;
	    	    else
	    	    	front->tan_curve_propagate 
				= fourth_order_elastic_curve_propagate;
	    	    break;
	    	default:
		    (void) printf("Unknown tangential propagator!\n");
		    clean_up(ERROR);
	    	}
	    }
	    else if (dim == 3)
	    {
	    	switch (string[0])
	    	{
	    	case 'n':
	    	case 'N':
	    	    front->tan_curve_propagate = NULL;
	    	    break;
	    	case 'e':
	    	case 'E':
	    	    if (string[1] == '2')
	    	    	front->interior_propagate 
				= second_order_elastic_surf_propagate;
	    	    else
	    	    	front->interior_propagate 
				= fourth_order_elastic_surf_propagate;
	    	    break;
	    	case 'p':
	    	case 'P':
	    	    front->interior_propagate 
				= fourth_order_elastic_set_propagate;
	    	    break;
	    	default:
		    (void) printf("Unknown tangential propagator!\n");
		    clean_up(ERROR);
	    	}
	    }
	}
	else
	    front->interior_propagate = fourth_order_elastic_set_propagate;

	CursorAfterString(infile,"Enter gravity:");
        for (i = 0; i < dim; ++i)
        {
            fscanf(infile,"%lf",&af_params->gravity[i]);
            (void) printf("%f ",af_params->gravity[i]);
	    iFparams->gravity[i] = af_params->gravity[i];
        }
        (void) printf("\n");
	CursorAfterString(infile,"Enter payload:");
        fscanf(infile,"%lf",&af_params->payload);
        (void) printf("%f\n",af_params->payload);
	if (af_params->no_fluid == NO)
	{
	    CursorAfterString(infile,
                        "Enter density and viscosity of the fluid:");
            fscanf(infile,"%lf %lf",&iFparams->rho2,&iFparams->mu2);
            (void) printf("%f %f\n",iFparams->rho2,iFparams->mu2);
            CursorAfterString(infile,"Enter surface tension:");
            fscanf(infile,"%lf",&iFparams->surf_tension);
            (void) printf("%f\n",iFparams->surf_tension);
            CursorAfterString(infile,"Enter factor of smoothing radius:");
            fscanf(infile,"%lf",&iFparams->smoothing_radius);
            (void) printf("%f\n",iFparams->smoothing_radius);
            CursorAfterString(infile,"Enter porosity of canopy:");
            fscanf(infile,"%lf",&af_params->gamma);
            (void) printf("%f\n",af_params->gamma);
            CursorAfterString(infile,"Enter area density of canopy:");
            fscanf(infile,"%lf",&af_params->area_dens);
            (void) printf("%f\n",af_params->area_dens);
	}
	af_params->n_tan = 1;
	CursorAfterString(infile,"Enter tangential sub step number:");
	fscanf(infile,"%d",&af_params->n_tan);
	(void) printf("%d\n",af_params->n_tan);
	if (dim == 3)
	{
	    CursorAfterString(infile,"Enter fabric spring constant:");
            fscanf(infile,"%lf",&af_params->ks);
            (void) printf("%f\n",af_params->ks);
            CursorAfterString(infile,"Enter fabric friction constant:");
            fscanf(infile,"%lf",&af_params->lambda_s);
            (void) printf("%f\n",af_params->lambda_s);
            CursorAfterString(infile,"Enter fabric point mass:");
            fscanf(infile,"%lf",&af_params->m_s);
            (void) printf("%f\n",af_params->m_s);
	    if (af_params->attach_gores == YES)
	    {
		CursorAfterString(infile,"Enter gore spring constant:");
        	fscanf(infile,"%lf",&af_params->kg);
        	(void) printf("%f\n",af_params->kg);
        	CursorAfterString(infile,"Enter gore friction constant:");
        	fscanf(infile,"%lf",&af_params->lambda_g);
        	(void) printf("%f\n",af_params->lambda_g);
        	CursorAfterString(infile,"Enter gore point mass:");
        	fscanf(infile,"%lf",&af_params->m_g);
        	(void) printf("%f\n",af_params->m_g);
	    }
	}
	CursorAfterString(infile,"Enter string spring constant:");
        fscanf(infile,"%lf",&af_params->kl);
        (void) printf("%f\n",af_params->kl);
        CursorAfterString(infile,"Enter string friction constant:");
        fscanf(infile,"%lf",&af_params->lambda_l);
        (void) printf("%f\n",af_params->lambda_l);
        CursorAfterString(infile,"Enter string point mass:");
        fscanf(infile,"%lf",&af_params->m_l);
        (void) printf("%f\n",af_params->m_l);
	//convert_to_point_mass(front,af_params);
	fclose(infile);
}	/* end initVelocityFunc */

extern int zero_velo(
	POINTER params,
        Front *front,
        POINT *p,
        HYPER_SURF_ELEMENT *hse,
        HYPER_SURF *hs,
        double *vel)
{
	vel[0] = vel[1] = vel[2] = 0.0;
	return YES;
}	/* end zero_velo */

extern int vertical_velo(
        POINTER params,
        Front *front,
        POINT *p,
        HYPER_SURF_ELEMENT *hse,
        HYPER_SURF *hs,
        double *vel)
{
	VERTICAL_PARAMS *vert_params = (VERTICAL_PARAMS*)params;
        double *coords = Coords(p);
	double dist,v_vert;
	double v0 = vert_params->v0;
	double stop_time = vert_params->stop_time;

	dist = sqrt(sqr(coords[0] - 0.5) + sqr(coords[1] - 0.5));
	v_vert = (0.15 - dist)/0.15;
	vel[0] = vel[1] = 0.0;
	if (front->time < stop_time)
	    vel[2] = v_vert*v0;
	else
	    vel[2] = 0.0;
	return YES;
}       /* end vertical_velo */


extern void elastic_point_propagate(
        Front *front,
        POINTER wave,
        POINT *oldp,
        POINT *newp,
        HYPER_SURF_ELEMENT *oldhse,
        HYPER_SURF         *oldhs,
        double              dt,
        double              *V)
{
	STATE *newsl,*newsr;
	STATE *sl,*sr;
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
        IF_FIELD *field = iFparams->field;
	int i, dim = front->rect_grid->dim;
	double *vort = field->vort;
	double **vort3d = field->vort3d;
	double **vel = field->vel;
	double *pres = field->pres;
	COMPONENT base_comp = positive_component(oldhs);
	double pp[MAXD],pm[MAXD],nor[MAXD],h;
	double area_dens = af_params->area_dens;
	double left_nor_speed,right_nor_speed;
	double dv[MAXD];

	if (af_params->no_fluid)
	{
	    fourth_order_point_propagate(front,wave,oldp,newp,oldhse,
				oldhs,dt,V);
	    ft_assign(left_state(newp),left_state(oldp),front->sizest);
	    ft_assign(right_state(newp),right_state(oldp),front->sizest);
	    return;
	}

	FT_GetStatesAtPoint(oldp,oldhse,oldhs,(POINTER*)&sl,(POINTER*)&sr);
	newsl = (STATE*)left_state(newp);
	newsr = (STATE*)right_state(newp);

	FT_NormalAtPoint(oldp,front,nor,NO_COMP);
	h = FT_GridSizeInDir(nor,front);
	for (i = 0; i < dim; ++i)
	{
	    pm[i] = Coords(oldp)[i] - h*nor[i];
	    pp[i] = Coords(oldp)[i] + h*nor[i];
	}

        FT_IntrpStateVarAtCoords(front,base_comp-1,pm,pres,
			getStatePres,&newsl->pres,&sl->pres);
        FT_IntrpStateVarAtCoords(front,base_comp+1,pp,pres,
			getStatePres,&newsr->pres,&sr->pres);

	/* Impulse is incremented by the fluid pressure force */
	for (i = 0; i < dim; ++i)
	{
	    dv[i] = 0.0;
	    if (debugging("rigid_canopy"))
		dv[i] = 0.0;
	    else if (front->step > 5)
		dv[i] = (sl->pres - sr->pres)*nor[i]*dt/area_dens;
	    newsr->Impct[i] = newsl->Impct[i] = sl->Impct[i] + dv[i];
	    newsr->vel[i] = newsl->vel[i] = sl->vel[i];
	}

	/* Interpolating vorticity for the hyper surface point */
	if (dim == 2)
        {
            FT_IntrpStateVarAtCoords(front,base_comp-1,pm,vort,
				getStateVort,&newsl->vort,&sl->vort);
            FT_IntrpStateVarAtCoords(front,base_comp+1,pp,vort,
				getStateVort,&newsr->vort,&sr->vort);
        }
	else if (dim == 3)
	{
	    for (i = 0; i < dim; ++i)
	    {
            	FT_IntrpStateVarAtCoords(front,base_comp-1,pm,vort3d[i],
			getStateVort3d[i],&newsl->vort3d[i],&sl->vort3d[i]);
            	FT_IntrpStateVarAtCoords(front,base_comp+1,pp,vort3d[i],
			getStateVort3d[i],&newsr->vort3d[i],&sr->vort3d[i]);
	    }
	}
}       /* elastic_point_propagate */

/*	Given string node, the function finds the corresponding
*	canopy surface.
*/

static SURFACE *canopy_of_string_node(NODE *n)
{
	SURFACE *canopy,**s;
	CURVE *c,**curves;
	int i,nc;
	boolean canopy_found = NO;

	canopy = NULL;
	nc = FT_NumOfNodeCurves(n);
	FT_VectorMemoryAlloc((POINTER*)&curves,nc,sizeof(CURVE*));
	FT_ArrayOfNodeCurves(n,curves);

	for (i = 0; i < nc; ++i)
	{
	    c = curves[i];
	    for (s = c->pos_surfaces; s && *s; ++s)
	    {
		if (wave_type(*s) == ELASTIC_BOUNDARY)
		{
		    canopy_found = YES;
		    canopy = *s;
		    break;
		}
	    }
	   if (canopy_found) break;
	    for (s = c->neg_surfaces; s && *s; ++s)
	    {
		if (wave_type(*s) == ELASTIC_BOUNDARY)
		{
		    canopy_found = YES;
		    canopy = *s;
		    break;
		}
	    }
	}
	return (canopy_found == YES) ? canopy : NULL;
}	/* end canopy_of_string_node */

extern int toroidal_velo(
        POINTER params,
        Front *front,
        POINT *p,
        HYPER_SURF_ELEMENT *hse,
        HYPER_SURF *hs,
        double *vel)
{
	TOROIDAL_PARAMS *toro_params = (TOROIDAL_PARAMS*)params;
        double *coords = Coords(p);
	int i,dim = front->rect_grid->dim;
	double v0 = toro_params->v0;
	double stop_time = toro_params->stop_time;
	double tcoords[2]; 		/* toroidal coords */
	double pvel[2];		/* projected velocity */
	double *tcen = toro_params->tcen; /* toroidal center of vel */
	double R0 = toro_params->R0; /* radial dist of poloidal center */
	double d1,d2;
	double s1,s2;
	double dx1,dy1;
	double dx2,dy2;
	double dtol = 0.000001*front->rect_grid->h[0];

	if (front->time >= stop_time)
	{
	    for (i = 0; i < dim; ++i)
		vel[i] = 0.0;
	    return YES;
	}
	/* Project 3D to 2D */
	tcoords[0] = sqrt(sqr(coords[0] - tcen[0]) + sqr(coords[1] - tcen[1]));
	tcoords[1] = coords[2] - tcen[2];

	dx1 = tcoords[0] - R0;
	dx2 = tcoords[0] + R0;
	dy1 = dy2 = tcoords[1];

	d1 = sqr(dx1) + sqr(dy1);
	d2 = sqr(dx2) + sqr(dy2);
	s1 = v0;
	s2 = -v0;
	if (d1 < dtol || d2 < dtol)
	    pvel[0] = pvel[1] = 0.0;
	else if (front->time < stop_time)
	{
	    pvel[0] =  s1*dy1/d1 + s2*dy2/d2;
	    pvel[1] = -s1*dx1/d1 - s2*dx2/d2;
	} 
	else
	    pvel[0] = pvel[1] = 0.0;
	    
	/* give it back to 3D */
	vel[0] = pvel[0]*(coords[0]-tcen[0])/tcoords[0];
	vel[1] = pvel[0]*(coords[1]-tcen[1])/tcoords[0];
	vel[2] = pvel[1];
	return YES;
}       /* end toroidal_velo */


extern int parabolic_velo(
        POINTER params,
        Front *front,
        POINT *p,
        HYPER_SURF_ELEMENT *hse,
        HYPER_SURF *hs,
        double *vel)
{
	PARABOLIC_PARAMS *para_params = (PARABOLIC_PARAMS*)params;
        double *coords = Coords(p);
	int i,dim = front->rect_grid->dim;
	double v0 = para_params->v0;
	double a = para_params->a;
	double *cen = para_params->cen;
	double R_sqr = 0.0;

	for (i = 0; i < dim-1; ++i)
	{
	    R_sqr += sqr(coords[i] - cen[i]);
	    vel[i] = 0.0;
	}
	vel[dim-1] = -0.5*a*R_sqr + v0;
	return YES;
}	/* end parabolic_velo */

extern int singular_velo(
        POINTER params,
        Front *front,
        POINT *p,
        HYPER_SURF_ELEMENT *hse,
        HYPER_SURF *hs,
        double *vel)
{
	SINGULAR_PARAMS *para_params = (SINGULAR_PARAMS*)params;
        double *coords = Coords(p);
	int i,dim = front->rect_grid->dim;
	double v0 = para_params->v0;
	double R = para_params->R;
	double *cen = para_params->cen;
	double r = 0.0;
	for (i = 0; i < dim-1; ++i)
	{
	    r += sqr(coords[i] - cen[i]);
	    vel[i] = 0.0;
	}
	r = sqrt(r);
	if (r < R)
	    vel[dim-1] = v0;
	else
	    vel[dim-1] = 0.0;
	return YES;
}	/* end sigular_velo */

extern void airfoil_point_propagate(
        Front *front,
        POINTER wave,
        POINT *oldp,
        POINT *newp,
        HYPER_SURF_ELEMENT *oldhse,
        HYPER_SURF         *oldhs,
        double              dt,
        double              *V)
{
        if (wave_type(oldhs) == ELASTIC_BOUNDARY)
            return elastic_point_propagate(front,wave,oldp,newp,oldhse,oldhs,
                                        dt,V);
        else
            return ifluid_point_propagate(front,wave,oldp,newp,oldhse,oldhs,
                                        dt,V);
}       /* airfoil_point_propagate */

extern void fourth_order_elastic_set_propagate(
	Front           *newfr,
        double           fr_dt)
{
	INTERFACE *old_intfc = newfr->interf;
	INTERFACE *new_intfc = newfr->interf;
	NODE *oldn,*newn;	/* old and new payload nodes */
	AF_NODE_EXTRA *extra;
	CURVE **old_str_curves,**old_mono_hsbdry,**old_gore_hsbdry;
	CURVE **new_str_curves,**new_mono_hsbdry,**new_gore_hsbdry;
	NODE **old_str_nodes,**old_gore_nodes;
	NODE **new_str_nodes,**new_gore_nodes;
	SURFACE *news;
	int i,num_str_curves,num_mono_hsbdry,num_gore_hsbdry;
	static PARACHUTE_SET old_geom_set;
	static PARACHUTE_SET new_geom_set;
	NODE **n;
	int num_gore_nodes = numOfGoreNodes(old_intfc);

	if (debugging("trace"))
	    (void) printf("Entering fourth_order_elastic_set_propagate()\n");
	for (n = old_intfc->nodes; n && *n; ++n)
	{
	    extra = (AF_NODE_EXTRA*)(*n)->extra;
	    if (extra == NULL || extra->af_node_type != LOAD_NODE) 
		continue;
	    oldn = *n;
	    break;
	}
	for (n = new_intfc->nodes; n && *n; ++n)
	{
	    extra = (AF_NODE_EXTRA*)(*n)->extra;
	    if (extra == NULL || extra->af_node_type != LOAD_NODE) 
		continue;
	    newn = *n;
	    break;
	}

	old_geom_set.load_node = oldn;
	new_geom_set.load_node = newn;

	num_str_curves = FT_NumOfNodeCurves(oldn);

	if (num_str_curves != FT_NumOfNodeCurves(newn))
	{
	    (void) printf("ERROR in af_node_propagate(): "
			  "old and new number of string curves not equal!\n");
	    clean_up(ERROR);
	}
	FT_VectorMemoryAlloc((POINTER*)&old_str_curves,num_str_curves,
			sizeof(CURVE*));
	FT_VectorMemoryAlloc((POINTER*)&new_str_curves,num_str_curves,
			sizeof(CURVE*));
	FT_VectorMemoryAlloc((POINTER*)&old_str_nodes,num_str_curves,
			sizeof(NODE*));
	FT_VectorMemoryAlloc((POINTER*)&new_str_nodes,num_str_curves,
			sizeof(NODE*));
	FT_VectorMemoryAlloc((POINTER*)&old_gore_nodes,num_gore_nodes,
			sizeof(NODE*));
	FT_VectorMemoryAlloc((POINTER*)&new_gore_nodes,num_gore_nodes,
			sizeof(NODE*));
	FT_ArrayOfNodeCurves(oldn,old_str_curves);
	FT_ArrayOfNodeCurves(newn,new_str_curves);

	old_geom_set.num_strings = num_str_curves;
	new_geom_set.num_strings = num_str_curves;

	old_geom_set.string_node = old_str_nodes;
	new_geom_set.string_node = new_str_nodes;

	old_geom_set.string_curves = old_str_curves;
	new_geom_set.string_curves = new_str_curves;

	for (i = 0; i < num_str_curves; ++i)
	{
	    old_str_nodes[i] = (old_str_curves[i]->start == oldn) ? 
			old_str_curves[i]->end : old_str_curves[i]->start;
	    new_str_nodes[i] = (new_str_curves[i]->start == newn) ? 
			new_str_curves[i]->end : new_str_curves[i]->start;
	}
	news = canopy_of_string_node(new_str_nodes[0]);

	num_gore_nodes  = numOfGoreNodes(old_intfc);
	num_gore_hsbdry = numOfGoreHsbdry(old_intfc);
	num_mono_hsbdry = numOfMonoHsbdry(old_intfc);

	FT_VectorMemoryAlloc((POINTER*)&old_mono_hsbdry,num_mono_hsbdry,
			sizeof(CURVE*));
	FT_VectorMemoryAlloc((POINTER*)&new_mono_hsbdry,num_mono_hsbdry,
			sizeof(CURVE*));
	num_mono_hsbdry = arrayOfMonoHsbdry(old_intfc,old_mono_hsbdry);
	num_mono_hsbdry = arrayOfMonoHsbdry(new_intfc,new_mono_hsbdry);
	FT_VectorMemoryAlloc((POINTER*)&old_gore_hsbdry,num_gore_hsbdry,
			sizeof(CURVE*));
	FT_VectorMemoryAlloc((POINTER*)&new_gore_hsbdry,num_gore_hsbdry,
			sizeof(CURVE*));
	num_gore_hsbdry = arrayOfGoreHsbdry(old_intfc,old_gore_hsbdry);
	num_gore_hsbdry = arrayOfGoreHsbdry(new_intfc,new_gore_hsbdry);

	old_geom_set.num_mono_hsbdry = num_mono_hsbdry;
	new_geom_set.num_mono_hsbdry = num_mono_hsbdry;
	old_geom_set.num_gore_hsbdry = num_gore_hsbdry;
	new_geom_set.num_gore_hsbdry = num_gore_hsbdry;
	old_geom_set.num_gore_nodes = num_gore_nodes;
	new_geom_set.num_gore_nodes = num_gore_nodes;

	old_geom_set.mono_hsbdry = old_mono_hsbdry;
	new_geom_set.mono_hsbdry = new_mono_hsbdry;
	old_geom_set.gore_hsbdry = old_gore_hsbdry;
	new_geom_set.gore_hsbdry = new_gore_hsbdry;
	getGoreNodes(old_intfc,old_gore_nodes);
	getGoreNodes(new_intfc,new_gore_nodes);
	old_geom_set.gore_nodes = old_gore_nodes;
	new_geom_set.gore_nodes = new_gore_nodes;
	old_geom_set.canopy = news;
	new_geom_set.canopy = news;
	old_geom_set.front = newfr;
	new_geom_set.front = newfr;

	if (debugging("para_set"))
	{
	    (void) printf("The parachute contains:\n");
	    (void) printf("%d string chord curves\n",num_str_curves);
	    (void) printf("%d canopy boundary curves\n",num_mono_hsbdry);
	    (void) printf("%d canopy gore curves\n",num_gore_hsbdry);
	}
	fourth_order_parachute_propagate(newfr,&old_geom_set,&new_geom_set);
	FT_FreeThese(10,old_str_curves,new_str_curves,
		       old_str_nodes,new_str_nodes,
		       old_mono_hsbdry,new_mono_hsbdry,
		       old_gore_hsbdry,new_gore_hsbdry,
		       old_gore_nodes,new_gore_nodes);
	
	if (debugging("trace"))
	    (void) printf("Leaving fourth_order_elastic_set_propagate()\n");
}	/* end fourth_order_elastic_set_propagate() */

static void convert_to_point_mass(
	Front *front,
	AF_PARAMS *af_params)
{
	INTERFACE *intfc;
	int num_str_pts,num_fabric_pts;
	SURFACE **s;
	int dim = Dimension(intfc);
	
	switch (dim)
	{
	case 2:
	    printf("In convert_to_point_mass(): 2D code needed!\n");
	    clean_up(ERROR);
	case 3:
	    for (s = intfc->surfaces; s && *s; ++s)
	    {
		if (wave_type(*s) == ELASTIC_BOUNDARY)
		{
		}
	    }
	}
}	/* end convert_to_point_mass */

extern void airfoil_curve_propagate(
        Front *front,
        POINTER wave,
	CURVE *oldc,
	CURVE *newc,
        double dt)
{
	int dim = front->rect_grid->dim;

	if (dim != 3) return;
	switch (hsbdry_type(oldc))
	{
	case STRING_HSBDRY:
	case MONO_COMP_HSBDRY:
	    return curve_prop_without_interaction(front,wave,oldc,newc,dt);
	case GORE_HSBDRY:
	    return curve_prop_with_interaction(front,wave,oldc,newc,dt);
	default:
	    return;
	}
}	/* end airfoil_curve_propagate */

static void curve_prop_without_interaction(
        Front *front,
        POINTER wave,
	CURVE *oldc,
	CURVE *newc,
        double dt)
{
	BOND *oldb,*newb;
	POINT *oldp,*newp;

	oldp = oldc->start->posn;
	newp = newc->start->posn;
	ft_assign(left_state(newp),left_state(oldp),front->sizest);
	ft_assign(right_state(newp),right_state(oldp),front->sizest);

	oldp = oldc->end->posn;
	newp = newc->end->posn;
	ft_assign(left_state(newp),left_state(oldp),front->sizest);
	ft_assign(right_state(newp),right_state(oldp),front->sizest);

	for (oldb = oldc->first, newb = newc->first; oldb != oldc->last;
		oldb = oldb->next, newb = newb->next)
	{
	    oldp = oldb->end;
	    newp = newb->end;
	    ft_assign(left_state(newp),left_state(oldp),front->sizest);
	    ft_assign(right_state(newp),right_state(oldp),front->sizest);
	}
}	/* end curve_prop_without_interaction */

static void curve_prop_with_interaction(
        Front *front,
        POINTER wave,
	CURVE *oldc,
	CURVE *newc,
        double dt)
{
	BOND *oldb,*newb;
	POINT *oldp,*newp;

	if (debugging("interact_curve"))
	{
	    (void) printf("Entering curve_prop_with_interaction()\n");
	}
	oldp = oldc->start->posn;
	newp = newc->start->posn;
	ft_assign(left_state(newp),left_state(oldp),front->sizest);
	ft_assign(right_state(newp),right_state(oldp),front->sizest);

	oldp = oldc->end->posn;
	newp = newc->end->posn;
	ft_assign(left_state(newp),left_state(oldp),front->sizest);
	ft_assign(right_state(newp),right_state(oldp),front->sizest);

	for (oldb = oldc->first, newb = newc->first; oldb != oldc->last;
		oldb = oldb->next, newb = newb->next)
	{
	    oldp = oldb->end;
	    newp = newb->end;
	    curve_point_prop_with_interaction(front,oldp,newp,oldb,dt);
	}
	if (debugging("interact_curve"))
	{
	    (void) printf("Leaving curve_prop_with_interaction()\n");
	}
}	/* end curve_prop_with_interaction */

static void curve_point_prop_with_interaction(
	Front *front,
	POINT *oldp,
	POINT *newp,
	BOND *oldb,
	double dt)
{
	BOND_TRI **btris;
	HYPER_SURF_ELEMENT *oldhse;
	HYPER_SURF         *oldhs;
	STATE *sl,*sr,*newsl,*newsr;
	double mag_nor,branch_nor[MAXD],nor[MAXD];
	double pm[MAXD],pp[MAXD],h;
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	IF_FIELD *field = iFparams->field;
	double **vel = field->vel;
        double *pres = field->pres;
	double area_dens = af_params->area_dens;
	double left_nor_speed,right_nor_speed,dv;
	COMPONENT base_comp;
	int i;

	sl = (STATE*)left_state(oldp);		
	sr = (STATE*)right_state(oldp);
	newsl = (STATE*)left_state(newp);	
	newsr = (STATE*)right_state(newp);

	for (i = 0; i < 3; ++i) nor[i] = 0.0;
	for (btris = Btris(oldb); btris && *btris; ++btris)
	{
	    oldp->hse = oldhse = Hyper_surf_element((*btris)->tri);
	    oldp->hs = oldhs = Hyper_surf((*btris)->surface);
	    FT_NormalAtPoint(oldp,front,branch_nor,NO_COMP);
	    base_comp = positive_component(oldhs);
	    for (i = 0; i < 3; ++i) nor[i] += branch_nor[i];
	}
	mag_nor = Mag3d(nor);
	for (i = 0; i < 3; ++i) nor[i] /= mag_nor;
	h = FT_GridSizeInDir(nor,front);
	for (i = 0; i < 3; ++i)
	{
	    pm[i] = Coords(oldp)[i] - h*nor[i];
            pp[i] = Coords(oldp)[i] + h*nor[i];
	}
	FT_IntrpStateVarAtCoords(front,base_comp-1,pm,pres,
                        getStatePres,&newsl->pres,&sl->pres);
        FT_IntrpStateVarAtCoords(front,base_comp+1,pp,pres,
                        getStatePres,&newsr->pres,&sr->pres);
	for (i = 0; i < 3; ++i)
        {
            FT_IntrpStateVarAtCoords(front,base_comp-1,pm,vel[i],
                        getStateVel[i],&newsl->vel[i],&sl->vel[i]);
            FT_IntrpStateVarAtCoords(front,base_comp+1,pp,vel[i],
                        getStateVel[i],&newsr->vel[i],&sr->vel[i]);
        }
        left_nor_speed = Dot3d(newsl->vel,nor);
        right_nor_speed = Dot3d(newsr->vel,nor);
	for (i = 0; i < 3; ++i)
        {
            newsl->vel[i] -= left_nor_speed*nor[i];
            newsr->vel[i] -= right_nor_speed*nor[i];
        }
	/* Impulse is incremented by the fluid pressure force */
        for (i = 0; i < 3; ++i)
        {
	    dv = 0.0;
	    if (front->step > 5)
	    	dv = (sl->pres - sr->pres)*nor[i]*dt/area_dens;
	    if (debugging("rigid_canopy"))
	    	dv = 0.0;
	    newsr->Impct[i] = newsl->Impct[i] = sl->Impct[i] + dv;
	}
}	/* end curve_point_prop_with_interaction */

static boolean is_gore_node(NODE*);

extern int numOfMonoHsbdry(
	INTERFACE *intfc)
{
	CURVE **c;
	int nc = 0;
	for (c = intfc->curves; c && *c; ++c)
	{
	    if (hsbdry_type(*c) == MONO_COMP_HSBDRY) nc++;
	} 
	return nc;
}	/* end numOfMonoBdry */

extern int numOfGoreHsbdry(
	INTERFACE *intfc)
{
	CURVE **c;
	int nc = 0;
	for (c = intfc->curves; c && *c; ++c)
	{
	    if (hsbdry_type(*c) == GORE_HSBDRY) nc++;
	} 
	return nc;
}	/* end numOfMonoBdry */

static int arrayOfMonoHsbdry(
	INTERFACE *intfc,
	CURVE **mono_curves)
{
	CURVE **c;
	int nc = 0;
	for (c = intfc->curves; c && *c; ++c)
	{
	    if (hsbdry_type(*c) == MONO_COMP_HSBDRY) 
	    {
		mono_curves[nc] = *c;
		nc++;
	    }
	} 
	return nc;
}	/* end arrayOfMonoBdry */

static int arrayOfGoreHsbdry(
	INTERFACE *intfc,
	CURVE **gore_curves)
{
	CURVE **c;
	int nc = 0;
	for (c = intfc->curves; c && *c; ++c)
	{
	    if (hsbdry_type(*c) == GORE_HSBDRY) 
	    {
		gore_curves[nc] = *c;
		nc++;
	    }
	} 
	return nc;
}	/* end arrayOfGoreBdry */

extern int numOfGoreNodes(
	INTERFACE *intfc)
{
	NODE **n;
	int num_gore_nodes = 0;
	AF_NODE_EXTRA *extra;

	for (n = intfc->nodes; n && *n; ++n)
	{
	    if ((*n)->extra == NULL)
		continue;
	    extra = (AF_NODE_EXTRA*)(*n)->extra;
	    if (extra->af_node_type == GORE_NODE)
		num_gore_nodes++;
	}
	return num_gore_nodes;
}	/* numOfGoreNodes */

static boolean is_gore_node(
	NODE *node)
{
	CURVE **c;
	AF_NODE_EXTRA *extra;

	if ((node)->extra == NULL)
	    return NO;
	extra = (AF_NODE_EXTRA*)(node)->extra;
	if (extra->af_node_type == GORE_NODE)
	    return YES;
	else 
	    return NO;
}	/* end is_gore_node */

static int getGoreNodes(
	INTERFACE *intfc,
	NODE **gore_nodes)
{
	NODE **n;
	int num_nodes = 0;

	for (n = intfc->nodes; n && *n; ++n)
	{
	    if (is_gore_node(*n))
		gore_nodes[num_nodes++] = *n;
	}
	return num_nodes;
}	/* getGoreNodes */


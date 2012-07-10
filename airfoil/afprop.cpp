#include <iFluid.h>
#include <airfoil.h>

static double (*getStateVel[3])(POINTER) = {getStateXvel,getStateYvel,
                                        getStateZvel};
static double (*getStateVort3d[3])(POINTER) = {getStateXvort,getStateYvort,
                                        getStateZvort};
static SURFACE *canopy_of_string_node(NODE*);
static void convert_to_point_mass(Front*,AF_PARAMS*);
static void airfoil_curve_propagate(Front*,POINTER,CURVE*,CURVE*,double);

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
	    	    	front->tan_surface_propagate 
				= second_order_elastic_surf_propagate;
	    	    else
	    	    	front->tan_surface_propagate 
				= fourth_order_elastic_surf_propagate;
	    	    break;
	    	case 'p':
	    	case 'P':
	    	    front->tan_surface_propagate 
				= fourth_order_elastic_set_propagate;
	    	    break;
	    	default:
		    (void) printf("Unknown tangential propagator!\n");
		    clean_up(ERROR);
	    	}
	    }
	}
	else
	    front->tan_surface_propagate = fourth_order_elastic_set_propagate;

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
	}
	CursorAfterString(infile,"Enter string spring constant:");
        fscanf(infile,"%lf",&af_params->kl);
        (void) printf("%f\n",af_params->kl);
        CursorAfterString(infile,"Enter string friction constant:");
        fscanf(infile,"%lf",&af_params->lambda_c);
        (void) printf("%f\n",af_params->lambda_c);
        CursorAfterString(infile,"Enter string point mass:");
        fscanf(infile,"%lf",&af_params->m_c);
        (void) printf("%f\n",af_params->m_c);
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
	double dv;
	double r;
	int j;
	static double Z[20],R[21];
	static int np[20],total_np;
	static int step = 0;
	static boolean first = YES;
	static double min_lift,max_lift,ave_lift;
	static double Ave_lift;

	/*TMP*/
	if (debugging("ave_lift"))
	{
	    if (first)
	    {
	    	first = NO;
	    	for (i = 0; i < 21; ++i)
	    	{
		    R[i] = i*5.0/20.0;
		    printf("R[%d] = %f\n",i,R[i]);
	    	}
	    	for (i = 0; i < 20; ++i)
	    	{
		    Z[i] = 0.0;
		    np[i] = 0;
	    	}
	    	max_lift = -HUGE;
	    	min_lift = HUGE;
	    	ave_lift = 0.0;
	    	total_np = 0;
	    	if (front->step == 2)
		    Ave_lift = ave_lift;
	    	else
		    Ave_lift = 0.0;
	    }
	    if (step != front->step)
	    {
	    	step = front->step;
	    	printf("Impact vs. radius:\n");
	    	for (i = 0; i < 20; ++i)
	    	{
		    Z[i] /= (double)np[i];
		    printf("%f  %f\n",R[i]+0.125,Z[i]);
		    Z[i] = 0.0;
		    np[i] = 0;
	    	}
	    	ave_lift /= (double)total_np;
	    	printf("Max lift = %f  Min lift = %f  Ave lift = %f\n",
			max_lift,min_lift,ave_lift);
	    	max_lift = -HUGE;
	    	min_lift = HUGE;
	    	ave_lift = 0.0;
	    	total_np = 0;
	    }
	}

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
	for (i = 0; i < dim; ++i)
	{
            FT_IntrpStateVarAtCoords(front,base_comp-1,pm,vel[i],
			getStateVel[i],&newsl->vel[i],&sl->vel[i]);
            FT_IntrpStateVarAtCoords(front,base_comp+1,pp,vel[i],
			getStateVel[i],&newsr->vel[i],&sr->vel[i]);
	}
	left_nor_speed = scalar_product(newsl->vel,nor,dim);
	right_nor_speed = scalar_product(newsr->vel,nor,dim);
	for (i = 0; i < dim; ++i)
	{
	    newsl->vel[i] -= left_nor_speed*nor[i];
	    newsr->vel[i] -= right_nor_speed*nor[i];
	}

	r = 0.0;
	for (i = 0; i < dim-1; ++i)
	    r += sqr(Coords(oldp)[i] - 7.0);
	r = sqrt(r);

	/* Impulse is incremented by the fluid pressure force */
	for (i = 0; i < dim; ++i)
	{
	    dv = 0.0;
	    if (debugging("rigid_canopy"))
		dv = 0.0;
	    else if (debugging("ave_lift") && i == 2)
		dv = Ave_lift*dt;
	    else if (front->step > 15)
		dv = (sl->pres - sr->pres)*nor[i]*dt/area_dens;
	    newsr->Impct[i] = newsl->Impct[i] = sl->Impct[i] + dv;
	    if (i == 2)
	    {
	    	for (j = 0; j < 20; ++j)
		{
		    if (R[j] < r && r < R[j+1])
		    {
			Z[j] += (sl->pres - sr->pres)*nor[i]/area_dens;
			np[j]++;
		    }
		}
		if (max_lift < (sl->pres - sr->pres)*nor[i]/area_dens)
		    max_lift = (sl->pres - sr->pres)*nor[i]/area_dens;
		if (min_lift > (sl->pres - sr->pres)*nor[i]/area_dens)
		    min_lift = (sl->pres - sr->pres)*nor[i]/area_dens;
		ave_lift += (sl->pres - sr->pres)*nor[i]/area_dens;
		total_np++;
	    }
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

extern int af_node_propagate(
	Front *front,
	POINTER wave,
	NODE *oldn,
	NODE *newn,
	RPROBLEM **rp,
	double dt,
	double *dt_frac,
	NODE_FLAG flag,
	POINTER user)
{
	AF_NODE_EXTRA *extra = (AF_NODE_EXTRA*)oldn->extra;
	CURVE **old_str_curves,**old_elas_bdry_curves;
	CURVE **new_str_curves,**new_elas_bdry_curves;
	NODE **old_str_nodes,**new_str_nodes;
	SURFACE *olds,*news;
	int i,num_str_curves,num_elas_bdry;
	PARACHUTE_SET old_geom_set;
	PARACHUTE_SET new_geom_set;

	if (extra == NULL) return GOOD_NODE;
	if (extra->af_node_type != LOAD_NODE) return GOOD_NODE;

	if (debugging("trace"))
	{
	    (void) printf("Entering af_node_propagate()\n");
	    (void) printf("posn = %f %f %f\n",Coords(oldn->posn)[0],
				Coords(oldn->posn)[1],Coords(oldn->posn)[2]);
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
	olds = canopy_of_string_node(old_str_nodes[0]);
	news = canopy_of_string_node(new_str_nodes[0]);

	printf("olds = %p  news = %p\n",(POINTER)olds,(POINTER)news);
	num_elas_bdry = FT_NumOfSurfCurves(olds);
	if (num_str_curves != FT_NumOfSurfCurves(news))
	{
	    (void) printf("ERROR in af_node_propagate(): "
			  "old and new number of canopy bady not equal!\n");
	    clean_up(ERROR);
	}
	FT_VectorMemoryAlloc((POINTER*)&old_elas_bdry_curves,num_elas_bdry,
			sizeof(CURVE*));
	FT_VectorMemoryAlloc((POINTER*)&new_elas_bdry_curves,num_elas_bdry,
			sizeof(CURVE*));
	FT_ArrayOfSurfCurves(olds,old_elas_bdry_curves);
	FT_ArrayOfSurfCurves(news,new_elas_bdry_curves);

	old_geom_set.num_elas_bdry_curves = num_elas_bdry;
	new_geom_set.num_elas_bdry_curves = num_elas_bdry;
	old_geom_set.elas_bdry_curves = old_elas_bdry_curves;
	new_geom_set.elas_bdry_curves = new_elas_bdry_curves;
	old_geom_set.canopy = olds;
	new_geom_set.canopy = news;
	old_geom_set.front = front;
	new_geom_set.front = front;

	fourth_order_parachute_propagate(front,&old_geom_set,&new_geom_set);
	FT_FreeThese(6,old_str_curves,new_str_curves,old_str_nodes,
		new_str_nodes,old_elas_bdry_curves,new_elas_bdry_curves);
	
	if (debugging("trace"))
	    (void) printf("Leaving af_node_propagate()\n");
	return GOOD_NODE;
}	/* end af_node_propagate */

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
	Front *front,
	Front           *newfr,
        INTERFACE       *intfc,
        SURFACE         *olds,
        SURFACE         *news,
        double           fr_dt)
{
	INTERFACE *old_intfc = news->interface;
	INTERFACE *new_intfc = news->interface;
	NODE *oldn,*newn;	/* old and new payload nodes */
	AF_NODE_EXTRA *extra;
	CURVE **old_str_curves,**old_elas_bdry_curves;
	CURVE **new_str_curves,**new_elas_bdry_curves;
	NODE **old_str_nodes,**new_str_nodes;
	int i,dim,num_str_curves,num_elas_bdry;
	PARACHUTE_SET old_geom_set;
	PARACHUTE_SET new_geom_set;
	NODE **n;

	if (wave_type(news) != ELASTIC_BOUNDARY) return;
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
	dim = front->rect_grid->dim;
	/*
	for (i = 0; i < dim; ++i)
	    old_geom_set.V_refs[i] = new_geom_set.V_refs[i] =
			oldn->posn->vel[i];
	*/
	for (i = 0; i < dim; ++i)
	    old_geom_set.V_refs[i] = new_geom_set.V_refs[i] = 0.0;

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

	num_elas_bdry = FT_NumOfSurfCurves(news);
	FT_VectorMemoryAlloc((POINTER*)&old_elas_bdry_curves,num_elas_bdry,
			sizeof(CURVE*));
	FT_VectorMemoryAlloc((POINTER*)&new_elas_bdry_curves,num_elas_bdry,
			sizeof(CURVE*));
	FT_ArrayOfSurfCurves(news,old_elas_bdry_curves);
	FT_ArrayOfSurfCurves(news,new_elas_bdry_curves);

	old_geom_set.num_elas_bdry_curves = num_elas_bdry;
	new_geom_set.num_elas_bdry_curves = num_elas_bdry;
	old_geom_set.elas_bdry_curves = old_elas_bdry_curves;
	new_geom_set.elas_bdry_curves = new_elas_bdry_curves;
	old_geom_set.canopy = news;
	new_geom_set.canopy = news;
	old_geom_set.front = newfr;
	new_geom_set.front = newfr;

	fourth_order_parachute_propagate(newfr,&old_geom_set,&new_geom_set);
	FT_FreeThese(6,old_str_curves,new_str_curves,old_str_nodes,
		new_str_nodes,old_elas_bdry_curves,new_elas_bdry_curves);
	
	if (debugging("trace"))
	    (void) printf("Leaving fourth_order_elastic_set_propagate()\n");
}	/* end fourth_order_elastic_set_propagate() */

static void convert_to_point_mass(
	Front *front,
	AF_PARAMS *af_params)
{
	INTERFACE *intfc;
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
        double              dt)
{
	BOND *oldb,*newb;
	POINT *oldp,*newp;
	int dim = front->rect_grid->dim;

	if (debugging("string"))
	    (void) printf("Entering airfoil_curve_propagate()\n");
	if (dim != 3) return;
	if (hsbdry_type(oldc) != STRING_HSBDRY &&
	    hsbdry_type(oldc) != MONO_COMP_HSBDRY) 
	    return;

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
	if (debugging("string"))
	    (void) printf("Leaving airfoil_curve_propagate()\n");
}	/* end airfoil_curve_propagate */


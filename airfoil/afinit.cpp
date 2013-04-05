#include <iFluid.h>
#include <airfoil.h>

static void zero_state(COMPONENT,double*,IF_FIELD*,int,int,IF_PARAMS*);
static void setInitialIntfcAF2d(Front*,LEVEL_FUNC_PACK*,char*);
static void setInitialIntfcAF3d(Front*,LEVEL_FUNC_PACK*,char*);

extern void setInitialIntfcAF(
        Front *front,
        LEVEL_FUNC_PACK *level_func_pack,
        char *inname)
{
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	FILE *infile = fopen(inname,"r");
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
	    return setInitialIntfcAF2d(front,level_func_pack,inname);
	case 3:
	    return setInitialIntfcAF3d(front,level_func_pack,inname);
	}
}	/* end setInitialIntfcAF */

static void setInitialIntfcAF3d(
        Front *front,
        LEVEL_FUNC_PACK *level_func_pack,
        char *inname)
{
	char string[100];
	FILE *infile = fopen(inname,"r");
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	int num_canopy;

	level_func_pack->set_3d_bdry = YES;
	level_func_pack->neg_component = LIQUID_COMP2;
        level_func_pack->pos_component = LIQUID_COMP2;	
	level_func_pack->func_params = NULL;
        level_func_pack->func = NULL;
	af_params->is_parachute_system = NO;
	af_params->cut_vent = NO;
	af_params->num_opt_round = 20;
        af_params->spring_model = MODEL1;	// default
	af_params->attach_gores = NO;		// default
	af_params->gore_len_fac = 1.0;		// default
	CursorAfterString(infile,"Enter number of canopy surfaces:");
	fscanf(infile,"%d",&num_canopy);
	(void) printf("%d\n",num_canopy);
	level_func_pack->num_mono_hs = num_canopy;

	(void) printf("Choices of initial surface are:\n");
	(void) printf("\tEllipsoid (E)\n");
	(void) printf("\tPlane (P)\n");
	(void) printf("\tT-10 (T)\n");
	(void) printf("\tNone (N)\n");
	CursorAfterString(infile,"Enter initial canopy surface type:");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	switch (string[0])
	{
	case 'E':
	case 'e':
	    initEllipticSurf(infile,front,level_func_pack);
	    break;
	case 'T':
	case 't':
	    initParabolicSurf(infile,front,level_func_pack);
	    break;
	case 'P':
	case 'p':
	    initPlaneSurf(infile,front,level_func_pack);
	    break;
	case 'N':
	case 'n':
	    break;
	}
	af_params->pert_params.pert_type = NO_PERT;
	(void) printf("Available perturbation types are:\n");
	(void) printf("\tNo perturbation (N)\n");
	(void) printf("\tParallel random perturbation (P)\n");
	(void) printf("\tOrthogonal random perturbation (O)\n");
	(void) printf("\tRadial perturbation (R)\n");
	(void) printf("\tLinear perturbation (L)\n");
	(void) printf("\tSine perturbation (S)\n");
	(void) printf("\tDefault is no perturbation\n");
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
	if (CursorAfterStringOpt(infile,
            "Entering number of canopy optimization rounds: "))
        {
            fscanf(infile,"%d",&af_params->num_opt_round);
            (void) printf("%d\n",af_params->num_opt_round);
        }
}	/* end setInitialIntfcAF3d */

static void setInitialIntfcAF2d(
        Front *front,
        LEVEL_FUNC_PACK *level_func_pack,
        char *inname)
{
	char string[100];
	FILE *infile = fopen(inname,"r");
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	int *gmax = front->rect_grid->gmax;
	double *L = front->rect_grid->L;
	double *U = front->rect_grid->U;
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

}	/* end setInitialIntfcAF2d */


static void zero_state(
        COMPONENT comp,
        double *coords,
	IF_FIELD *field,
	int index,
        int dim,
        IF_PARAMS *af_params)
{
        int i;
        for (i = 0; i < dim; ++i)
            field->vel[i][index] = 0.0;
        field->vel[1][index] = 0.0;
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
	boolean is_gore_node;
	boolean gore_curve_attached;
	STRING_PARAMS *string_params;

	string_params = (STRING_PARAMS*)level_func_pack->string_params;

	payload_node = NULL;
	for (n = intfc->nodes; n && *n; ++n)
	{
	    gore_curve_attached = NO;
	    is_payload_node = YES;
	    is_gore_node = YES;
	    for (c = (*n)->in_curves; c && *c; ++c)
	    {
		if (hsbdry_type(*c) != STRING_HSBDRY)
		    is_payload_node = NO;
		if (hsbdry_type(*c) != MONO_COMP_HSBDRY &&
			 hsbdry_type(*c) != GORE_HSBDRY)
		    is_gore_node = NO;
		if (hsbdry_type(*c) == GORE_HSBDRY)
		    gore_curve_attached = YES;
	    } 
	    for (c = (*n)->out_curves; c && *c; ++c)
	    {
		if (hsbdry_type(*c) != STRING_HSBDRY)
		    is_payload_node = NO;
		if (hsbdry_type(*c) != MONO_COMP_HSBDRY &&
			 hsbdry_type(*c) != GORE_HSBDRY)
		    is_gore_node = NO;
		if (hsbdry_type(*c) == GORE_HSBDRY)
		    gore_curve_attached = YES;
	    } 
	    if (!gore_curve_attached) is_gore_node = NO;
	    if (is_payload_node == YES)
	    {
		payload_node = *n;
	    }
	    if (is_gore_node == YES)
	    {
	    	FT_ScalarMemoryAlloc((POINTER*)&extra,sizeof(AF_NODE_EXTRA));
	    	extra->af_node_type = GORE_NODE;
	    	(*n)->extra = (POINTER)extra;
	    }
	}
	if (payload_node != NULL)
	{
	    FT_ScalarMemoryAlloc((POINTER*)&extra,sizeof(AF_NODE_EXTRA));
	    extra->af_node_type = LOAD_NODE;
	    payload_node->extra = (POINTER)extra;
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
	}
}	/* end setRestartAirfoilIntfc */

extern void printAfExtraDada(
	Front *front,
	char *out_name)
{
	INTERFACE *intfc = front->interf;
        STATE *sl,*sr;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
	int i,dim = intfc->dim;
	FILE *outfile;
	char filename[200];
	CURVE **c;
	NODE **n;
	BOND *b;

	sprintf(filename,"%s/state.ts%s",out_name,
                        right_flush(front->step,7));
#if defined(__MPI__)
        if (pp_numnodes() > 1)
            sprintf(filename,"%s-nd%s",filename,right_flush(pp_mynode(),4));
#endif /* defined(__MPI__) */
        sprintf(filename,"%s-afdata",filename);
        outfile = fopen(filename,"w");

	fprintf(outfile,"\nAirfoil extra front state data:\n");

	next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
            if (wave_type(hs) != ELASTIC_BOUNDARY) continue;
            FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
            for (i = 0; i < dim; ++i)
                fprintf(outfile,"%24.18g %24.18g\n",sl->impulse[i],sr->impulse[i]);
            for (i = 0; i < dim; ++i)
                fprintf(outfile,"%24.18g ",p->vel[i]);
	    fprintf(outfile,"\n");
        }
	for (c = intfc->curves; c && *c; ++c)
	{
	    b = (*c)->first;	p = b->start;
	    sl = (STATE*)left_state(p);
	    sr = (STATE*)right_state(p);
            for (i = 0; i < dim; ++i)
                fprintf(outfile,"%24.18g ",p->vel[i]);
	    fprintf(outfile,"\n");
            fprintf(outfile,"%24.18g %24.18g\n",sl->pres,sr->pres);
            for (i = 0; i < dim; ++i)
                fprintf(outfile,"%24.18g ",sl->impulse[i]);
	    fprintf(outfile,"\n");
            for (i = 0; i < dim; ++i)
                fprintf(outfile,"%24.18g ",sr->impulse[i]);
	    fprintf(outfile,"\n");
	    for (b = (*c)->first; b != NULL; b = b->next)
	    {
		p = b->end;
	    	sl = (STATE*)left_state(p);
	    	sr = (STATE*)right_state(p);
            	for (i = 0; i < dim; ++i)
                    fprintf(outfile,"%24.18g ",p->vel[i]);
            	fprintf(outfile,"%24.18g %24.18g\n",sl->pres,sr->pres);
	    	fprintf(outfile,"\n");
            	for (i = 0; i < dim; ++i)
                    fprintf(outfile,"%24.18g ",sl->impulse[i]);
	    	fprintf(outfile,"\n");
            	for (i = 0; i < dim; ++i)
                    fprintf(outfile,"%24.18g ",sr->impulse[i]);
	    	fprintf(outfile,"\n");
	    }
	}
	for (n = intfc->nodes; n && *n; ++n)
	{
	    p = (*n)->posn;
	    sl = (STATE*)left_state(p);
	    sr = (STATE*)right_state(p);
            for (i = 0; i < dim; ++i)
                fprintf(outfile,"%24.18g ",p->vel[i]);
            fprintf(outfile,"%24.18g %24.18g\n",sl->pres,sr->pres);
	    fprintf(outfile,"\n");
            for (i = 0; i < dim; ++i)
                fprintf(outfile,"%24.18g ",sl->impulse[i]);
	    fprintf(outfile,"\n");
            for (i = 0; i < dim; ++i)
                fprintf(outfile,"%24.18g ",sr->impulse[i]);
	    fprintf(outfile,"\n");
	}
	fprintf(outfile,"\nCurve extra data:\n");
	for (c = intfc->curves; c && *c; ++c)
	{
	    C_PARAMS *c_params = (C_PARAMS*)(*c)->extra;
	    if (c_params == NULL)
                fprintf(outfile,"curve extra: no\n");
	    else
	    {
                fprintf(outfile,"curve extra: yes\n");
                fprintf(outfile,"point_mass = %24.18g\n",c_params->point_mass);
                fprintf(outfile,"load_mass = %24.18g\n",c_params->load_mass);
                fprintf(outfile,"load_type = %d\n",c_params->load_type);
                fprintf(outfile,"dir = %d\n",c_params->dir);
	    }
	}
	fprintf(outfile,"\nNode extra data:\n");
	for (n = intfc->nodes; n && *n; ++n)
	{
	    AF_NODE_EXTRA *n_params = (AF_NODE_EXTRA*)(*n)->extra;
	    if (n_params == NULL)
                fprintf(outfile,"node extra: no\n");
	    else
	    {
                fprintf(outfile,"node extra: yes\n");
                fprintf(outfile,"af_node_type = %d\n",n_params->af_node_type);
	    }
	}
}	/* end printAfExtraDada */

extern void readAfExtraDada(
	Front *front,
	char *restart_name)
{
	INTERFACE *intfc = front->interf;
        STATE *sl,*sr;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
	int i,dim = intfc->dim;
	FILE *infile;
	char filename[200];
	CURVE **c;
	NODE **n;
	BOND *b;
	char string[100];

        sprintf(filename,"%s-afdata",restart_name);
        infile = fopen(filename,"r");

	printf("filename = %s\n",filename);
	next_output_line_containing_string(infile,
		"Airfoil extra front state data:");

	next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
            if (wave_type(hs) != ELASTIC_BOUNDARY) continue;
            FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
            for (i = 0; i < dim; ++i)
                fscanf(infile,"%lf %lf\n",&sl->impulse[i],&sr->impulse[i]);
            for (i = 0; i < dim; ++i)
                fscanf(infile,"%lf ",&p->vel[i]);
	    fscanf(infile,"\n");
        }
	for (c = intfc->curves; c && *c; ++c)
	{
	    b = (*c)->first;	p = b->start;
	    sl = (STATE*)left_state(p);
	    sr = (STATE*)right_state(p);
            for (i = 0; i < dim; ++i)
                fscanf(infile,"%lf ",&p->vel[i]);
            fscanf(infile,"%lf %lf",&sl->pres,&sr->pres);
	    fscanf(infile,"\n");
            for (i = 0; i < dim; ++i)
                fscanf(infile,"%lf ",&sl->impulse[i]);
            for (i = 0; i < dim; ++i)
                fscanf(infile,"%lf ",&sr->impulse[i]);
	    fscanf(infile,"\n");
	    for (b = (*c)->first; b != NULL; b = b->next)
	    {
		p = b->end;
	    	sl = (STATE*)left_state(p);
	    	sr = (STATE*)right_state(p);
            	for (i = 0; i < dim; ++i)
                    fscanf(infile,"%lf ",&p->vel[i]);
            	fscanf(infile,"%lf %lf",&sl->pres,&sr->pres);
	    	fscanf(infile,"\n");
            	for (i = 0; i < dim; ++i)
               	    fscanf(infile,"%lf ",&sl->impulse[i]);
            	for (i = 0; i < dim; ++i)
                    fscanf(infile,"%lf ",&sr->impulse[i]);
	    }
	}
	for (n = intfc->nodes; n && *n; ++n)
	{
	    p = (*n)->posn;
	    sl = (STATE*)left_state(p);
	    sr = (STATE*)right_state(p);
            for (i = 0; i < dim; ++i)
                fscanf(infile,"%lf ",&p->vel[i]);
            fscanf(infile,"%lf %lf",&sl->pres,&sr->pres);
	    fscanf(infile,"\n");
            for (i = 0; i < dim; ++i)
                fscanf(infile,"%lf ",&sl->impulse[i]);
            for (i = 0; i < dim; ++i)
                fscanf(infile,"%lf ",&sr->impulse[i]);
	}
	next_output_line_containing_string(infile,"Curve extra data:");
	for (c = intfc->curves; c && *c; ++c)
	{
	    C_PARAMS *c_params;
	    next_output_line_containing_string(infile,"curve extra:");
            fscanf(infile,"%s",string);
	    if (string[0] == 'n') continue;
	    FT_ScalarMemoryAlloc((POINTER*)&c_params,sizeof(C_PARAMS));
	    next_output_line_containing_string(infile,"point_mass = ");
            fscanf(infile,"%lf",&c_params->point_mass);
	    next_output_line_containing_string(infile,"load_mass = ");
            fscanf(infile,"%lf",&c_params->load_mass);
	    next_output_line_containing_string(infile,"load_type = ");
            fscanf(infile,"%d",&c_params->load_type);
	    next_output_line_containing_string(infile,"dir = ");
            fscanf(infile,"%d",&c_params->dir);
	    (*c)->extra = (POINTER)c_params;
	}
	next_output_line_containing_string(infile,"Node extra data:");
	for (n = intfc->nodes; n && *n; ++n)
	{
	    AF_NODE_EXTRA *n_params;
	    next_output_line_containing_string(infile,"node extra:");
            fscanf(infile,"%s",string);
	    if (string[0] == 'n') continue;
	    FT_ScalarMemoryAlloc((POINTER*)&n_params,sizeof(AF_NODE_EXTRA));
	    next_output_line_containing_string(infile,"af_node_type =");
            fscanf(infile,"%d",&n_params->af_node_type);
	}
}	/* end readAfExtraDada */


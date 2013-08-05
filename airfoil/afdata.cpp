#include <iFluid.h>
#include <airfoil.h>

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
                fprintf(outfile,"%24.18g %24.18g\n",sl->impulse[i],
					sr->impulse[i]);
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
	fclose(outfile);
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
	    fgetstring(infile,"curve extra:");
            fscanf(infile,"%s",string);
	    if (string[0] == 'n') continue;
	    FT_ScalarMemoryAlloc((POINTER*)&c_params,sizeof(C_PARAMS));
	    fgetstring(infile,"point_mass = ");
            fscanf(infile,"%lf",&c_params->point_mass);
	    fgetstring(infile,"load_mass = ");
            fscanf(infile,"%lf",&c_params->load_mass);
	    fgetstring(infile,"load_type = ");
            fscanf(infile,"%d",&c_params->load_type);
	    fgetstring(infile,"dir = ");
            fscanf(infile,"%d",&c_params->dir);
	    (*c)->extra = (POINTER)c_params;
	}
	next_output_line_containing_string(infile,"Node extra data:");
	for (n = intfc->nodes; n && *n; ++n)
	{
	    AF_NODE_EXTRA *n_params;
	    fgetstring(infile,"node extra:");
            fscanf(infile,"%s",string);
	    if (string[0] == 'n') continue;
	    FT_ScalarMemoryAlloc((POINTER*)&n_params,sizeof(AF_NODE_EXTRA));
	    fgetstring(infile,"af_node_type =");
            fscanf(infile,"%d",&n_params->af_node_type);
	    (*n)->extra = (POINTER)n_params;
	}
}	/* end readAfExtraDada */

extern void printHyperSurfQuality(
	Front *front)
{
	INTERFACE *intfc = front->interf;
	int dim = Dimension(intfc);
	CURVE **c,*curve;
	SURFACE **s,*surf;
	BOND *bond;
	TRI *tri;
	double max_area,min_area,max_length,min_length;
	double scaled_area,len[MAXD];
	int i;
	RECT_GRID *gr = &topological_grid(intfc);
	double *h = gr->h;

	switch (dim)
	{
	case 2:
	    max_length = 0.0;
	    min_length = HUGE;
	    for (c = intfc->curves; c && *c; ++c)
	    {
		if (wave_type(*c) != ELASTIC_BOUNDARY)
		    continue;
		curve = *c;
		for (bond = curve->first; bond != NULL; bond = bond->next)
		{
		    if (scaled_bond_length(bond,h,dim) > max_length)
			max_length = scaled_bond_length(bond,h,dim);
		    if (scaled_bond_length(bond,h,dim) < min_length)
			min_length = scaled_bond_length(bond,h,dim);
		}
	    }
	    (void) printf("\n\nElastic surface quality:\n");
	    (void) printf("min_length = %f\n",min_length);
	    (void) printf("max_length = %f\n",max_length);
	    (void) printf("\n\n");
	    break;
	case 3:
	    max_length = 0.0;
	    min_length = HUGE;
	    for (c = intfc->curves; c && *c; ++c)
	    {
		if (hsbdry_type(*c) != MONO_COMP_HSBDRY &&
		    hsbdry_type(*c) != GORE_HSBDRY)
		    continue;
		curve = *c;
		for (bond = curve->first; bond != NULL; bond = bond->next)
		{
		    if (scaled_bond_length(bond,h,dim) > max_length)
			max_length = scaled_bond_length(bond,h,dim);
		    if (scaled_bond_length(bond,h,dim) < min_length)
			min_length = scaled_bond_length(bond,h,dim);
		}
	    }
	    (void) printf("\n\nElastic curve quality:\n");
	    (void) printf("min_scaled_length = %14.10f\n",min_length);
	    (void) printf("max_scaled_length = %14.10f\n",max_length);

	    max_length = 0.0;
	    min_length = HUGE;
	    for (c = intfc->curves; c && *c; ++c)
	    {
		if (hsbdry_type(*c) != STRING_HSBDRY)
		    continue;
		curve = *c;
		for (bond = curve->first; bond != NULL; bond = bond->next)
		{
		    if (scaled_bond_length(bond,h,dim) > max_length)
			max_length = scaled_bond_length(bond,h,dim);
		    if (scaled_bond_length(bond,h,dim) < min_length)
			min_length = scaled_bond_length(bond,h,dim);
		}
	    }
	    (void) printf("\nElastic string quality:\n");
	    (void) printf("min_scaled_length = %14.10f\n",min_length);
	    (void) printf("max_scaled_length = %14.10f\n",max_length);

	    max_area = max_length = 0.0;
	    min_area = min_length = HUGE;
	    for (s = intfc->surfaces; s && *s; ++s)
	    {
		if (wave_type(*s) != ELASTIC_BOUNDARY)
		    continue;
		surf = *s;
		for (tri = first_tri(surf); !at_end_of_tri_list(tri,surf); 
				tri = tri->next)
		{
		    scaled_tri_params(tri,h,&scaled_area,len);
		    if (scaled_area > max_area) max_area = scaled_area;
		    if (scaled_area < min_area) min_area = scaled_area;
		    for (i = 0; i < 3; ++i)
		    {
			if (len[i] > max_length) 
			    max_length = len[i];
			if (len[i] < min_length) 
			    min_length = len[i];
		    }
		}
	    }
	    (void) printf("\nElastic surface quality:\n");
	    (void) printf("min_scaled_area = %14.10f\n",min_area);  
	    (void) printf("max_scaled_area = %14.10f\n",max_area); 
	    (void) printf("min_scaled_tri_side = %14.10f\n",sqrt(min_length));
	    (void) printf("max_scaled_tri_side = %14.10f\n",sqrt(max_length));
	    (void) printf("\n\n");
	    break;
	}
}	/* end printHyperSurfQuality */

extern void optimizeElasticMesh(
	Front *front)
{
	INTERFACE *intfc = front->interf;
	RECT_GRID *gr = computational_grid(intfc);
	boolean nothing_done;
	int i,status;
	CURVE **c,*curve;
	SURFACE **s,*surf;
	SCALED_REDIST_PARAMS scaled_redist_params;
	int old_string_pts,new_string_pts,old_canopy_pts,new_canopy_pts;
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	int num_opt_round;

	if (debugging("trace"))
	    (void) printf("Entering optimizeElasticMesh()\n");
	if (debugging("optimize_intfc"))
	{
	    (void) printf("Quality of mesh before optimization:\n");
	    printHyperSurfQuality(front);
	    (void) printf("Checking consistency of interface\n");
	    consistent_interface(front->interf);
	    (void) printf("Checking completed\n");
	    gview_plot_interface("gview-before-optimize",intfc);
	    if (debugging("no_optimize"))
		clean_up(0);
	}
	if (debugging("no_optimize")) return;
	if (gr->dim == 2) return;

	num_opt_round = af_params->num_opt_round;
	scaled_redist_params.min_scaled_bond_length = 0.45/2.0;
	scaled_redist_params.max_scaled_bond_length = 1.05/2.0;

	scaled_redist_params.min_scaled_tri_area = 0.1083;
	scaled_redist_params.max_scaled_tri_area = 0.4330;
	scaled_redist_params.min_scaled_tri_area = 0.1083/2.0;
	scaled_redist_params.max_scaled_tri_area = 0.4330/2.0;
	scaled_redist_params.min_scaled_side_length = 0.45/2.0;
	scaled_redist_params.max_scaled_side_length = 1.05/2.0;
	scaled_redist_params.aspect_tol = 3.0;

	old_string_pts = old_canopy_pts = 0;
	for (s = intfc->surfaces; s && *s; ++s)
	    if (wave_type(*s) == ELASTIC_BOUNDARY)
		old_canopy_pts += FT_NumOfSurfPoints(*s);
	for (c = intfc->curves; c && *c; ++c)
	    if (hsbdry_type(*c) == STRING_HSBDRY)
		old_string_pts += FT_NumOfCurvePoints(*c) - 2;

	printf("num_opt_round = %d\n",num_opt_round);
	num_opt_round = 20;
	for (i = 0; i < num_opt_round; ++i)
	{
	    status = YES;
	    if (debugging("optimize_intfc"))
		(void) printf("Optimization round %d\n",i);
	    for (c = intfc->curves; c && *c; ++c)
	    {
	    	if (hsbdry_type(*c) != MONO_COMP_HSBDRY &&
		    hsbdry_type(*c) != STRING_HSBDRY &&
		    hsbdry_type(*c) != GORE_HSBDRY)
		    continue;
	    	curve = *c;
	    	nothing_done = FT_OptimizeCurveMesh(front,curve,
				scaled_redist_params);
	    	status *= (int)nothing_done;
	    }
	    for (s = intfc->surfaces; s && *s; ++s)
	    {
	    	if (wave_type(*s) != ELASTIC_BOUNDARY)
		    continue;
	    	surf = *s;
	    	nothing_done = FT_OptimizeSurfMesh(front,surf,
				scaled_redist_params);
	    	status *= (int)nothing_done;
	    }
	    FT_ParallelExchIntfcBuffer(front);
	    if (debugging("optimize_intfc"))
	    {
		(void) printf("Quality of mesh after %d-th round:\n",i);
	    	printHyperSurfQuality(front);
		(void) printf("Checking consistency of interface\n");
		consistent_interface(front->interf);
		(void) printf("After checking\n");
	    }
	    if (status) break;
	}

	new_string_pts = new_canopy_pts = 0;
	for (s = intfc->surfaces; s && *s; ++s)
	    if (wave_type(*s) == ELASTIC_BOUNDARY)
		new_canopy_pts += FT_NumOfSurfPoints(*s);
	for (c = intfc->curves; c && *c; ++c)
	    if (hsbdry_type(*c) == STRING_HSBDRY)
		new_string_pts += FT_NumOfCurvePoints(*c) - 2;
	if (debugging("optimize_intfc"))
	{
	    gview_plot_interface("gview-after-optimize",intfc);
	    //clean_up(0);
	}
	if (debugging("trace"))
	    (void) printf("Leaving optimizeElasticMesh()\n");
}	/* end optimizeElasticMesh */

extern void modifyInitialization(
	Front *front)
{
	char *inname = InName(front);
	FILE *infile = fopen(inname,"r");
	char string[200];
	double disp[MAXD],center[MAXD];
	double phi,theta;
	INTERFACE *intfc = front->interf;
	
	if (CursorAfterStringOpt(infile,
            "Entering yes to modify initialization:"))
        {
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] != 'y' && string[0] != 'Y')
		return;
        }
	CursorAfterString(infile,
		"Enter yes for translation of interior interface:");
        fscanf(infile,"%s",string);
        (void) printf("%s\n",string);
	if (string[0] != 'y' || string[0] != 'Y')
	{
	    CursorAfterString(infile,"Enter displacement of translation:");
            fscanf(infile,"%lf %lf %lf",disp,disp+1,disp+2);
            (void) printf("%f %f %f\n",disp[0],disp[1],disp[2]);
	    I_TransInteriorIntfcPoints(intfc,disp);
	}
	CursorAfterString(infile,
		"Enter yes for rotation of interior interface:");
        fscanf(infile,"%s",string);
        (void) printf("%s\n",string);
	if (string[0] != 'y' || string[0] != 'Y')
	{
	    CursorAfterString(infile,"Enter center of rotation:");
            fscanf(infile,"%lf %lf %lf",center,center+1,center+2);
            (void) printf("%f %f %f\n",center[0],center[1],center[2]);
	    CursorAfterString(infile,"Enter azimuthal and polar angles:");
            fscanf(infile,"%lf %lf",&phi,&theta);
            (void) printf("%f %f\n",phi,theta);
	    theta *= PI/180.0;
	    phi *= PI/180.0;
	    I_SphericalRotateInteriorIntfcPoints(intfc,center,phi,theta);
	}
}	/* end modifyInitialization */

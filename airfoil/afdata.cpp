/***************************************************************
FronTier is a set of libraries that implements differnt types of 
Front Traking algorithms. Front Tracking is a numerical method for 
the solution of partial differential equations whose solutions have 
discontinuities.  

Copyright (C) 1999 by The University at Stony Brook. 

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
****************************************************************/

#include <iFluid.h>
#include <airfoil.h>

static void naturalStressOfTri(TRI*,double);
static void singleCanopyModification(Front*);
static void bifurcateCanopyModification(Front*);
static void copyParachuteSet(PARACHUTE_SET,PARACHUTE_SET*);
static void rotateParachuteSet(PARACHUTE_SET*,double*,double,double);

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
            fscanf(infile,"%d",(int*)&c_params->load_type);
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
            fscanf(infile,"%d",(int*)&n_params->af_node_type);
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
	}
	if (debugging("trace"))
	    (void) printf("Leaving optimizeElasticMesh()\n");
}	/* end optimizeElasticMesh */

extern void modifyInitialization(
	Front *front)
{
	char *inname = InName(front);
	FILE *infile = fopen(inname,"r");
	char string[200],input_string[200];
	boolean bifurcate_initialization;
	
	if (!CursorAfterStringOpt(infile,
            "Entering yes to modify initialization:"))
	{
	    fclose(infile);
	    return;
	}
	else
        {
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] != 'y' && string[0] != 'Y')
	    {
	    	fclose(infile);
		return;
	    }
        }
	bifurcate_initialization = NO;
	if (CursorAfterStringOpt(infile,
            "Enter yes to bifurcate initialization:"))
	{
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'y' || string[0] == 'Y')
		bifurcate_initialization = YES;
        }
    	fclose(infile);
	if (bifurcate_initialization)
	    bifurcateCanopyModification(front);
	else
	    singleCanopyModification(front);
}	/* end modifyInitialization */

extern void setStressColor(
	Front *front)
{
	INTERFACE *intfc = front->interf;
	SURFACE **s;
	TRI *tri;
	POINT *p;
	double f[MAXD];
	int i,j;
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	double ks = af_params->ks;
	double max_color = -HUGE;
	double min_color = HUGE;
	int n,N;
	double *color;
	
	n = 0;
	intfc_surface_loop(intfc,s)
	{
	    if (Boundary(*s)) continue;
	    surf_tri_loop(*s,tri)
	    {
		n++;
		naturalStressOfTri(tri,ks);
		if (max_color < tri->color)
		    max_color = tri->color;
		if (min_color > tri->color)
		    min_color = tri->color;
	    }
	}
	FT_VectorMemoryAlloc((POINTER*)&color,n,sizeof(double));
	n = 0;
	intfc_surface_loop(intfc,s)
	{
	    if (Boundary(*s)) continue;
	    surf_tri_loop(*s,tri)
		tri->color = log(tri->color-min_color+1);
	}
	/* Smoothing loop */
	intfc_surface_loop(intfc,s)
	{
	    if (Boundary(*s)) continue;
	    I_SmoothSurfColor(*s,3);
	}
}	/* end setStressColor */

extern void initMovieStress(
	char *inname,
	Front *front)
{
	char string[100];
	FILE *infile = fopen(inname,"r");
	front->print_gview_color = NO;
	if (CursorAfterStringOpt(infile,
            "Type y to plot surface stress: "))
        {
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'y' || string[0] == 'Y')
	    {
		front->print_gview_color = YES;
		FT_AddVtkIntfcMovieVariable(front,"VMSTRESS");
		return;
	    }
        }
        fclose(infile);
}	/* end initMovieStress */

static void naturalStressOfTri(
	TRI *tri,
	double ks)
{
	double tau[3];
	double sigma[3];
	int i,j;
	double len,len0;
	double vec[3];
	double s[3],c[3];
	double b1,b2,arg,sigma1,sigma2;

	for (i = 0; i < 3; ++i)
	{
	    len0 = tri->side_length0[i];
	    for (j = 0; j < 3; ++j)
	    {
		vec[j] = Coords(Point_of_tri(tri)[(i+1)%3])[j] -
			Coords(Point_of_tri(tri)[i])[j];
	    }
	    len = Mag3d(vec);
	    tau[i] = ks*(len - len0);
	    c[i] = vec[0]/len;
	    s[i] = vec[1]/len;
	}
	// Convert to Cartesian tensor
	for (i = 0; i < 3; ++i)
	{
	    sigma[i] = sqr(c[i])*tau[0] + sqr(s[i])*tau[1] + s[i]*c[i]*tau[2];
	}
	// Diagonalize the stress tensor for principal directions
	b1 = -(sigma[0] + sigma[1]);
	b2 = sigma[0]*sigma[1] - 0.25*sqr(sigma[2]);
	arg = sqr(b1) - 4.0*b2;
	sigma1 = 0.5*(-b1 + sqrt(arg));
	sigma2 = 0.5*(-b1 - sqrt(arg));
	// Use von Mises stress as a measure
	tri->color = sqrt(sqr(sigma1) + sqr(sigma2) - sigma1*sigma2);
}	/* end naturalStressOfTri */

extern void poisson_ratio(
	Front *front)
{
        INTERFACE *intfc = front->interf;
        SURFACE **s;
	CURVE **c;
	BOND *b;
	double x_min, y_min, x_max, y_max;

        intfc_surface_loop(intfc,s)
        {
            if (Boundary(*s)) continue;
	    surf_pos_curve_loop(*s,c)
	    {
		b = (*c)->first;
		x_min = x_max = b->start->_coords[0];
		y_min = y_max = b->start->_coords[1];
		for (b = (*c)->first; b != NULL; b = b->next)
		{
		    if (x_min > b->end->_coords[0])
			x_min = b->end->_coords[0];

		    if (x_max < b->end->_coords[0])
			x_max = b->end->_coords[0];

		    if (y_min > b->end->_coords[1])
			y_min = b->end->_coords[1];

		    if (y_max < b->end->_coords[1])
			y_max = b->end->_coords[1];
		}
		printf("curve x-d boundary: x_min: %f\t x_max: %f\n",
					x_min,x_max);
		printf("curve y-d boundary: y_min: %f\t y_max: %f\n\n",
					y_min,y_max);
	    }
	}
}	/* poisson_ratio */

static void singleCanopyModification(
	Front *front)
{
	char *inname = InName(front);
	FILE *infile = fopen(inname,"r");
	char string[200],input_string[200];
	double disp[MAXD],center[MAXD];
	double phi,theta;
	INTERFACE *intfc = front->interf;
	double L[MAXD],U[MAXD];
	int i,dim,gmax[MAXD];

	dim = FT_Dimension();
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
	if (CursorAfterStringOpt(infile,
            "Entering yes to modify computational grid:"))
        {
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'y' || string[0] == 'Y')
	    {
		for (i = 0; i < dim; ++i)
        	{
	            sprintf(input_string,
				"New domain limit in %d-th dimension:",i);
	            CursorAfterString(infile,input_string);
	            fscanf(infile,"%lf %lf",&L[i],&U[i]);
	            (void) printf("%f %f\n",L[i],U[i]);
        	}
		CursorAfterString(infile,"New computational grid:");
        	for (i = 0; i < dim; ++i)
        	{
	            fscanf(infile,"%d",&gmax[i]);
		    (void) printf("%d ",gmax[i]);
        	}
        	(void) printf("\n");
		FT_ResetDomainAndGrid(front,L,U,gmax);
	    }
        }
	fclose(infile);
}	/* end singleCanopyModification */

static void bifurcateCanopyModification(
	Front *front)
{
	char *inname = InName(front);
	FILE *infile = fopen(inname,"r");
	char string[200],input_string[200];
	PARACHUTE_SET *parachute_set;
	int num_canopy;
	double disp[MAXD],center[MAXD];
	double phi,theta;
	INTERFACE *intfc = front->interf;
	double L[MAXD],U[MAXD];
	int i,dim,gmax[MAXD];

	dim = FT_Dimension();

	CursorAfterString(infile,
		"Enter total number of canopy sets:");
        fscanf(infile,"%d",&num_canopy);
        (void) printf("%d\n",num_canopy);
	FT_VectorMemoryAlloc((POINTER*)&parachute_set,num_canopy,
			sizeof(PARACHUTE_SET));
	/* Get the original set */
	assembleParachuteSet(front,&parachute_set[0],2);
	for (i = 1; i < num_canopy; ++i)
	{
	    copyParachuteSet(parachute_set[0],&parachute_set[i]);
	}
	for (i = 0; i < num_canopy; ++i)
	{
	    sprintf(string,"For canopy set %d",i+1);
	    CursorAfterString(infile,string);
	    (void) printf("\n");

	    	CursorAfterString(infile,
		"Enter yes for translation of canopy set:");
            fscanf(infile,"%s",input_string);
            (void) printf("%s\n",input_string);
	    if (input_string[0] == 'y' || input_string[0] == 'Y')
	    {
	    	CursorAfterString(infile,"Enter displacement of translation:");
            	fscanf(infile,"%lf %lf %lf",disp,disp+1,disp+2);
            	(void) printf("%f %f %f\n",disp[0],disp[1],disp[2]);
	    }
	    CursorAfterString(infile,
		"Enter yes for rotation of canopy set:");
            fscanf(infile,"%s",input_string);
            (void) printf("%s\n",input_string);
	    if (input_string[0] == 'y' || input_string[0] == 'Y')
	    {
	    	CursorAfterString(infile,"Enter center of rotation:");
            	fscanf(infile,"%lf %lf %lf",center,center+1,center+2);
            	(void) printf("%f %f %f\n",center[0],center[1],center[2]);
	    	CursorAfterString(infile,"Enter azimuthal and polar angles:");
            	fscanf(infile,"%lf %lf",&phi,&theta);
            	(void) printf("%f %f\n",phi,theta);
	    	theta *= PI/180.0;
	    	phi *= PI/180.0;
	    	rotateParachuteSet(&parachute_set[i],center,phi,theta);
	    }
	}
	InstallNewLoadNode(front,num_canopy);

	if (CursorAfterStringOpt(infile,
            "Entering yes to modify computational grid:"))
        {
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'y' || string[0] == 'Y')
	    {
		for (i = 0; i < dim; ++i)
        	{
	            sprintf(input_string,
				"New domain limit in %d-th dimension:",i);
	            CursorAfterString(infile,input_string);
	            fscanf(infile,"%lf %lf",&L[i],&U[i]);
	            (void) printf("%f %f\n",L[i],U[i]);
        	}
		CursorAfterString(infile,"New computational grid:");
        	for (i = 0; i < dim; ++i)
        	{
	            fscanf(infile,"%d",&gmax[i]);
		    (void) printf("%d ",gmax[i]);
        	}
        	(void) printf("\n");
		FT_ResetDomainAndGrid(front,L,U,gmax);
	    }
        }
	fclose(infile);
}	/* end bifurcateCanopyModification */

static void copyParachuteSet(
	PARACHUTE_SET orig_set,
	PARACHUTE_SET *copy_set)
{
	int i,j,ns,nc,nn;
	INTERFACE *cur_intfc = current_interface();
	Front *front = orig_set.front;
	NODE *start,*end;
	CURVE **c,**pos_curves,**neg_curves;
	AF_NODE_EXTRA *extra;

	set_current_interface(front->interf);

	ns = copy_set->num_surfs = orig_set.num_surfs;
	nc = copy_set->num_curves = orig_set.num_curves;
	nn = copy_set->num_nodes = orig_set.num_nodes;

	for (i = 0; i < nn; ++i)
	{
	    copy_set->nodes[i] = copy_node(orig_set.nodes[i]);
	    if (is_load_node(orig_set.nodes[i]))
	    {
		copy_set->load_node = copy_set->nodes[i];
        	FT_ScalarMemoryAlloc((POINTER*)&extra,sizeof(AF_NODE_EXTRA));
		copy_set->load_node->extra = extra;
		extra->af_node_type = LOAD_NODE;
	    }
	}
	for (i = 0; i < nc; ++i)
	{
	    start = end = NULL;
	    for (j = 0; j < nn; ++j)
	    {
		if (orig_set.curves[i]->start == orig_set.nodes[j])
		    start = copy_set->nodes[j];
		if (orig_set.curves[i]->end == orig_set.nodes[j])
		    end = copy_set->nodes[j];
	    }
	    if (start == NULL || end == NULL)
	    {
		printf("Cannot find start or end node\n");
		clean_up(ERROR);
	    }
	    copy_set->curves[i] = copy_curve(orig_set.curves[i],start,end);
	}
	for (i = 0; i < ns; ++i)
	{
	    pos_curves = neg_curves = NULL;
	    for (j = 0; j < nc; ++j)
	    {
		surf_pos_curve_loop(orig_set.surfs[i],c)
		{
		    if (*c == orig_set.curves[j])
			unique_add_to_pointers(copy_set->curves[j],&pos_curves);
		}
		surf_neg_curve_loop(orig_set.surfs[i],c)
		{
		    if (*c == orig_set.curves[j])
			unique_add_to_pointers(copy_set->curves[j],&neg_curves);
		}
	    }
	    copy_set->surfs[i] = copy_surface(orig_set.surfs[i],pos_curves,
				neg_curves,YES);
	}
}	/* end copyParachuteSet */

static void rotateParachuteSet(
	PARACHUTE_SET *parachute_set,
	double *center,
	double phi,
	double theta)
{
	int i;
	for (i = 0; i < parachute_set->num_surfs; ++i)
	{
	    I_SphericalRotateInteriorSurfPoints(parachute_set->surfs[i],
					center,phi,theta);
	}
	for (i = 0; i < parachute_set->num_curves; ++i)
	{
	    I_SphericalRotateInteriorCurvePoints(parachute_set->curves[i],
					center,phi,theta);
	}
	for (i = 0; i < parachute_set->num_nodes; ++i)
	{
	    I_SphericalRotatePoint(parachute_set->nodes[i]->posn,
					center,phi,theta,NO);
	}
}	/* end rotateParachuteSet */

extern void InstallNewLoadNode(
	Front *front,
	int num_canopy)
{
	INTERFACE *intfc = front->interf;
	FILE *infile = fopen(InName(front),"r");
	NODE **n, *sec_nload, *nload;
	CURVE **string_curves;
	AF_NODE_EXTRA *extra;
	BOND *bond;
	double center[MAXD],newload[MAXD],dir[MAXD],coords[MAXD];
	double spacing,*h = front->rect_grid->h;
	int i,j,k,nb;
        INTERFACE *cur_intfc;

        if (CursorAfterStringOpt(infile,"Enter new load position:"))
	{
            fscanf(infile,"%lf %lf %lf",newload,newload+1,newload+2);
            (void) printf("%f %f %f\n",newload[0],newload[1],newload[2]);
	}
	else
	{
	    fclose(infile);
	    return;
	}

        CursorAfterString(infile,"Enter center of rotation:");
        fscanf(infile,"%lf %lf %lf",center,center+1,center+2);
        (void) printf("%f %f %f\n",center[0],center[1],center[2]);

        cur_intfc = current_interface();
	set_current_interface(intfc);
	FT_VectorMemoryAlloc((POINTER*)&string_curves,num_canopy+1,
                                sizeof(CURVE*));
	sec_nload = make_node(Point(center));
        FT_ScalarMemoryAlloc((POINTER*)&extra,sizeof(AF_NODE_EXTRA));
        extra->af_node_type = SEC_LOAD_NODE;
        sec_nload->extra = (POINTER)extra;

	i = 0;
	intfc_node_loop(intfc,n)
	{
	    extra = (AF_NODE_EXTRA*)((*n)->extra);
	    if (extra == NULL)
		continue;
	    if (extra->af_node_type == LOAD_NODE)
	    {
		extra->af_node_type = THR_LOAD_NODE;
		string_curves[i] = make_curve(0,0,(*n),sec_nload);
		hsbdry_type(string_curves[i]) = STRING_HSBDRY;
		spacing = separation((*n)->posn,sec_nload->posn,3);
		for (j = 0; j < 3; ++j)
                    dir[j] = (Coords(sec_nload->posn)[j] -
                        Coords((*n)->posn)[j])/spacing;
		nb = (int)spacing/(0.25*h[0]);
		spacing /= (double)nb;
		bond = string_curves[i]->first;
		bond->length0 = spacing;
		for (j = 1; j < nb; ++j)
		{
		    for (k = 0; k < 3; ++k)
                    	coords[k] = Coords((*n)->posn)[k] +
                                        j*dir[k]*spacing;
                    insert_point_in_bond(Point(coords),bond,string_curves[i]);
                    bond = bond->next;
		    bond->length0 = spacing;
		}
		i++;
	    }
	}	

	nload = make_node(Point(newload));
        FT_ScalarMemoryAlloc((POINTER*)&extra,sizeof(AF_NODE_EXTRA));
        extra->af_node_type = LOAD_NODE;
        nload->extra = (POINTER)extra;

	string_curves[num_canopy] = make_curve(0,0,sec_nload,nload);
	hsbdry_type(string_curves[num_canopy]) = STRING_HSBDRY;
	spacing = separation(sec_nload->posn,nload->posn,3);
	for (j = 0; j < 3; ++j)
            dir[j] = (Coords(nload->posn)[j] -
                        Coords(sec_nload->posn)[j])/spacing;
	nb = (int)spacing/(0.25*h[0]);
	spacing /= (double)nb;
	bond = string_curves[num_canopy]->first;
	bond->length0 = spacing;
	for (j = 1; j < nb; ++j)
	{
	    for (k = 0; k < 3; ++k)
               	coords[k] = Coords(sec_nload->posn)[k] +
        		j*dir[k]*spacing;
            insert_point_in_bond(Point(coords),bond,string_curves[i]);
            bond = bond->next;
	    bond->length0 = spacing;
	}
	set_current_interface(cur_intfc);
	fclose(infile);
	FT_FreeThese(1,string_curves);
}	/* end InstallNewLoadNode */

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


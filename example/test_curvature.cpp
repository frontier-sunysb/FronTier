#include <FronTier.h>

char *in_name,*restart_state_name,*restart_name,*out_name;
boolean RestartRun;
int RestartStep;
boolean binary = YES;

/********************************************************************
 *	Level function parameters for the initial interface 	    *
 ********************************************************************/

typedef struct {
	int dim;
	double center[MAXD];
        double R; 
} CIRCLE_PARAMS;

	/*  Function Declarations */
static void set_geom_functions(Front*);
static void test_normal(POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*,double*,Front*);
static void reset_hyper_sphere_points(Front*,CIRCLE_PARAMS);
static void print_error_norms(Front*,char*,CIRCLE_PARAMS);
static double test_curvature(POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*,Front*);
static double level_circle_func(POINTER,double*);

int main(int argc, char **argv)
{
	static Front front;
	static RECT_GRID comp_grid;
	static F_BASIC_DATA f_basic;
	static LEVEL_FUNC_PACK level_func_pack;
	CIRCLE_PARAMS circle_params;	/* level function parameters */
	int i,dim;

	FT_Init(argc,argv,&f_basic);
	dim = f_basic.dim;

        f_basic.size_of_intfc_state = 0;

        in_name                 = f_basic.in_name;
        restart_state_name      = f_basic.restart_state_name;
        out_name                = f_basic.out_name;
        restart_name            = f_basic.restart_name;
        RestartRun              = f_basic.RestartRun;
        RestartStep             = f_basic.RestartStep;

        sprintf(restart_name,"%s.ts%s",restart_name,right_flush(RestartStep,7));
#if defined(__MPI__)
        sprintf(restart_name,"%s-nd%s",restart_name,right_flush(pp_mynode(),4));
#endif /* defined(__MPI__) */

	FT_ReadSpaceDomain(in_name,&f_basic);
	FT_StartUp(&front,&f_basic);

	    /* Initialize interface through level function */

	for (i = 0; i < dim; ++i)
	    circle_params.center[i] = 0.5;
	circle_params.R = 0.3;
	circle_params.dim = dim;

	level_func_pack.neg_component = 1;
	level_func_pack.pos_component = 2;
	level_func_pack.func_params = (POINTER)&circle_params;
	level_func_pack.func = level_circle_func;
	level_func_pack.wave_type = FIRST_PHYSICS_WAVE_TYPE;

	FT_InitIntfc(&front,&level_func_pack);
	reset_hyper_sphere_points(&front,circle_params);
	//set_geom_functions(&front);
	print_error_norms(&front,out_name,circle_params);

	clean_up(0);
}

/********************************************************************
 *	Hyper sphere level function for the initial interface    *
 ********************************************************************/

static double level_circle_func(
        POINTER func_params,
        double *coords)
{
	CIRCLE_PARAMS *circle_params = (CIRCLE_PARAMS*)func_params;
	double *center,R,dist;
	int i,dim;

	center = circle_params->center;
	R  = circle_params->R;
	dim  = circle_params->dim;

	dist = 0.0;
	for (i = 0; i < dim; ++i)
	    dist += sqr(coords[i] - center[i]);	
	dist = sqrt(dist) - R;
	return dist;
}	/* end level_circle_func */

static void reset_hyper_sphere_points(
	Front *front,
	CIRCLE_PARAMS circle_params)
{
	POINT              *p;
        HYPER_SURF_ELEMENT *hse;
        HYPER_SURF         *hs;
	INTERFACE	   *intfc = front->interf;
	double	 	   *center = circle_params.center;
	double	 	   r, R = circle_params.R;
	int		   i,dim = circle_params.dim;

	next_point(intfc,NULL,NULL,NULL);
	while (next_point(intfc,&p,&hse,&hs))
	{
	    if (wave_type(hs) < FIRST_PHYSICS_WAVE_TYPE) continue;
	    r = 0.0;
	    for (i = 0; i < dim; ++i)
		r += sqr(Coords(p)[i] - center[i]);
	    r = sqrt(r);
	    for (i = 0; i < dim; ++i)
		Coords(p)[i] = R/r*(Coords(p)[i] - center[i]) + center[i];
	}
}	/* end adjust_hyper_sphere_points */

static void set_geom_functions(Front *front)
{
	INTERFACE *intfc = front->interf;
	interface_normal(intfc) = test_normal;
	interface_curvature(intfc) = test_curvature;
}

static void print_error_norms(
	Front *front,
	char *out_name,
	CIRCLE_PARAMS circle_params)
{
	POINT              *p;
        HYPER_SURF_ELEMENT *hse;
        HYPER_SURF         *hs;
	INTERFACE	   *intfc = front->interf;
	double	 	   *center = circle_params.center;
	double	 	   r, R = circle_params.R;
	int		   i,dim = circle_params.dim;
	double		   nor[MAXD],exact_nor[MAXD];
	double		   kappa, exact_kappa = 1.0/R;
	double		   N1_curvature,N2_curvature,Ni_curvature;
	double		   N1_norm[MAXD],N2_norm[MAXD],Ni_norm[MAXD];
	double		   N;
	FILE		   *outfile = fopen(out_name,"w");

	N1_curvature = N2_curvature = Ni_curvature = 0.0;
	for (i = 0; i < dim; ++i)
	    N1_norm[i] = N2_norm[i] = Ni_norm[i] = 0.0;
	N = 0.0;

	next_point(intfc,NULL,NULL,NULL);
	fprintf(outfile,"Computed Curvature      Exact Curvature\n");
	while (next_point(intfc,&p,&hse,&hs))
	{
	    if (wave_type(hs) < FIRST_PHYSICS_WAVE_TYPE) continue;
	    r = 0.0;
	    for (i = 0; i < dim; ++i)
		r += sqr(Coords(p)[i] - center[i]);
	    r = sqrt(r);
	    for (i = 0; i < dim; ++i)
		exact_nor[i] = (Coords(p)[i] - center[i])/r;

	    GetFrontNormal(p,hse,hs,nor,front);
	    GetFrontCurvature(p,hse,hs,&kappa,front);

	    fprintf(outfile,"%18.14f   %18.14f\n",kappa,exact_kappa);
	    N1_curvature += fabs(kappa - exact_kappa);
	    N2_curvature += sqr(kappa - exact_kappa);
	    if (Ni_curvature < fabs(kappa - exact_kappa)) 
		Ni_curvature = fabs(kappa - exact_kappa);
	    for (i = 0; i < dim; ++i)
	    {
		N1_norm[i] = fabs(nor[i] - exact_nor[i]);
		N2_norm[i] = sqr(nor[i] - exact_nor[i]);
	    	if (Ni_norm[i] < fabs(nor[i] - exact_nor[i])) 
		    Ni_norm[i] = fabs(nor[i] - exact_nor[i]);
	    }
	    N += 1.0;
	}
	N1_curvature = N1_curvature/N;
	N2_curvature = sqrt(N2_curvature/N);
	for (i = 0; i < dim; ++i)
	{
	    N1_norm[i] = N1_norm[i]/N;
	    N2_norm[i] = sqrt(N2_norm[i]/N);
	}
	
	fprintf(outfile,"\n\n");
	fprintf(outfile,"Error norms for normal:\n");
	fprintf(outfile,"L-1 norm:  ");
	for (i = 0; i < dim; ++i)
	    fprintf(outfile,"%24.18g  ",N1_norm[i]);
	fprintf(outfile,"\n");
	fprintf(outfile,"L-2 norm:  ");
	for (i = 0; i < dim; ++i)
	    fprintf(outfile,"%24.18g  ",N2_norm[i]);
	fprintf(outfile,"\n");
	fprintf(outfile,"L-I norm:  ");
	for (i = 0; i < dim; ++i)
	    fprintf(outfile,"%24.18g  ",Ni_norm[i]);
	fprintf(outfile,"\n\n");
	fprintf(outfile,"Error norms for curvature:\n");
	fprintf(outfile,"L-1 norm:   %24.18g\n",N1_curvature);
	fprintf(outfile,"L-2 norm:   %24.18g\n",N2_curvature);
	fprintf(outfile,"L-I norm:   %24.18g\n",Ni_curvature);
}	/* end print_error_norms */

static void test_normal(
        POINT              *p,
        HYPER_SURF_ELEMENT *hse,
        HYPER_SURF         *hs,
        double              *nor,
        Front              *front)
{
	/*
	switch (front->rect_grid->dim)
	{
	case 2:
	    return test_nromal2d();
	case 3:
	    return test_nromal3d();
	}
	*/
}	/* end test_normal */

static double test_curvature(
        POINT              *p,
        HYPER_SURF_ELEMENT *hse,
        HYPER_SURF         *hs,
        Front              *fr)
{
	/*
	switch (front->rect_grid->dim)
	{
	case 2:
	    return test_curvature2d();
	case 3:
	    return test_curvature3d();
	}
	*/
}	/* end test_curvature */

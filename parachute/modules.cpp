#include <iFluid.h>
#include <airfoil.h>

static void initParachuteDefault(Front*);

static void initSingleModule(FILE*,Front*);

static void initCanopySurface(FILE*,Front*,SURFACE**);
static void initEllipticSurface(FILE*,Front*,SURFACE**);
static void initParabolicSurface(FILE*,Front*,SURFACE**);
static void initFlatSurface(FILE*,Front*,SURFACE**);

static void cutToCanopyType(FILE*,SURFACE*);
static void cutToRectangle(FILE*,SURFACE*);
static void cutToEllipse(FILE*,SURFACE*);
static void cutToCross(FILE*,SURFACE*);
static void cutToWing(FILE*,SURFACE*);

static void installStringChord(FILE*,Front*,SURFACE*);
static void installParallelString(FILE*,Front*,SURFACE*);
static void installAngularString(FILE*,Front*,SURFACE*);

extern void initParachuteModules(Front *front)
{
	int i,num_canopy;
	FILE *infile = fopen(InName(front),"r");

	if (debugging("trace"))
	    (void) printf("Entering initParachuteModules()\n");

	if (debugging("set_module"))
	    gview_plot_interface("module-step-1",front->interf);

	initParachuteDefault(front);

	CursorAfterString(infile,"Enter number of canopy surfaces:");
        fscanf(infile,"%d",&num_canopy);
        (void) printf("%d\n",num_canopy);

	for (i = 0; i < num_canopy; ++i)
	{
	    initSingleModule(infile,front);
	}

	if (debugging("trace"))
	    (void) printf("Leaving initParachuteModules()\n");
}	/* end initParachuteModules */

static void initSingleModule(
	FILE *infile,
	Front *front)
{
	SURFACE *surf;
	initCanopySurface(infile,front,&surf);
	if (debugging("set_module"))
	    gview_plot_interface("module-step-2",front->interf);

	cutToCanopyType(infile,surf);
	if (debugging("set_module"))
	    gview_plot_interface("module-step-3",front->interf);

	installStringChord(infile,front,surf);
	if (debugging("set_module"))
	    gview_plot_interface("module-step-4",front->interf);
	
}	/* end initSingleModule */

static void initCanopySurface(
	FILE *infile,
	Front *front,
	SURFACE **surf)
{
	char string[200];
        (void) printf("Available canopy surface types are:\n");
        (void) printf("\tFLAT (F)\n");
        (void) printf("\tPARABOLIC (P)\n");
        (void) printf("\tELLIPTIC (E)\n");
	CursorAfterString(infile,"Enter canopy surface type:");
        fscanf(infile,"%s",string);
        (void) printf("%s\n",string);
	switch (string[0])
        {
	case 'E':
        case 'e':
            initEllipticSurface(infile,front,surf);
            break;
        case 'P':
        case 'p':
            initParabolicSurface(infile,front,surf);
            break;
        case 'F':
        case 'f':
            initFlatSurface(infile,front,surf);
            break;
	default:
	    (void) printf("Unknown canopy surface type\n");
	    clean_up(ERROR);
	}
}	/* end initCanopySurface */

static void initEllipticSurface(
	FILE *infile,
	Front *front,
	SURFACE **surf)
{
}	/* end initEllipticSurface */

static void initParabolicSurface(
	FILE *infile,
	Front *front,
	SURFACE **surf)
{
}	/* end initParabolicSurface */

static void initFlatSurface(
	FILE *infile,
	Front *front,
	SURFACE **surf)
{
	double plane_nor[MAXD],plane_pt[MAXD];
	double height;
	double *L = front->rect_grid->L;
	double *U = front->rect_grid->U;
	int i;
	COMPONENT amb_comp = front->interf->default_comp;

	CursorAfterString(infile,"Enter the height of the plane:");
	fscanf(infile,"%lf",&height);
        (void) printf("%f\n",height);	

	for (i = 0; i < 2; ++i)
	{
	    plane_nor[i] = 0.0;
	    plane_pt[i] = 0.5*(L[i] + U[i]);
	}
	plane_nor[2] = 1.0;
	plane_pt[2] = height;
	FT_MakePlaneSurf(front,plane_nor,plane_pt,NO,amb_comp+1,amb_comp,
			ELASTIC_BOUNDARY,surf);
}	/* end initFlatSurface */

static void cutToCanopyType(
	FILE *infile,
	SURFACE *surf)
{
	char string[200];

	(void) printf("Available types of canopy boundaries are:\n");
	(void) printf("\tCircular (C)\n");
	(void) printf("\tRectanglar (R)\n");
	(void) printf("\tElliptic (E)\n");
	(void) printf("\tCross (X)\n");
	(void) printf("\tWing (W)\n");
	CursorAfterString(infile,"Enter type of canopy boundary:");
        fscanf(infile,"%s",string);
        (void) printf("%s\n",string);
	switch (string[0])
	{
	case 'C':
	case 'c':
	    cutToEllipse(infile,surf);
	    break;
	case 'R':
	case 'r':
	    cutToRectangle(infile,surf);
	    break;
	case 'E':
	case 'e':
	    cutToEllipse(infile,surf);
	    break;
	case 'X':
	case 'x':
	    cutToCross(infile,surf);
	    break;
	case 'W':
	case 'w':
	    cutToWing(infile,surf);
	    break;
	}
}	/* end cutToCanopyType */

static void cutToRectangle(
	FILE *infile,
	SURFACE *surf)
{
}	/* end cutToRectangle */

static void cutToEllipse(
	FILE *infile,
	SURFACE *surf)
{
}	/* end cutToEllipse */

static void cutToCross(
	FILE *infile,
	SURFACE *surf)
{
	static CROSS_CONSTR_PARAMS constr_params;
	static PLANE_PARAMS plane_params;
	double **insert_coords;
	double L[MAXD],U[MAXD];
	int i;
	double *N = plane_params.N;
	double *P = plane_params.P;

	(void) printf("Input the two crossing rectangles\n");
	CursorAfterString(infile,"Enter lower bounds of first rectangle:");
	fscanf(infile,"%lf %lf",&constr_params.L1[0],&constr_params.L1[1]);
	(void) printf("%f %f\n",constr_params.L1[0],constr_params.L1[1]);
	CursorAfterString(infile,"Enter upper bounds of first rectangle:");
	fscanf(infile,"%lf %lf",&constr_params.U1[0],&constr_params.U1[1]);
	(void) printf("%f %f\n",constr_params.U1[0],constr_params.U1[1]);
	CursorAfterString(infile,"Enter lower bounds of second rectangle:");
	fscanf(infile,"%lf %lf",&constr_params.L2[0],&constr_params.L2[1]);
	(void) printf("%f %f\n",constr_params.L2[0],constr_params.L2[1]);
	CursorAfterString(infile,"Enter upper bounds of second rectangle:");
	fscanf(infile,"%lf %lf",&constr_params.U2[0],&constr_params.U2[1]);
	(void) printf("%f %f\n",constr_params.U2[0],constr_params.U2[1]);
	FT_MatrixMemoryAlloc((POINTER*)&insert_coords,4,MAXD,sizeof(double));

	/* Cut to rectangle first to avoid corner problem */
	for (i = 0; i < 2; ++i)
	{
	    L[i] = (constr_params.L1[i] < constr_params.L2[i]) ? 
			constr_params.L1[i] : constr_params.L2[i];
	    U[i] = (constr_params.U1[i] > constr_params.U2[i]) ? 
			constr_params.U1[i] : constr_params.U2[i];
	}
	for (i = 0; i < 3; ++i)
	    P[i] = N[i] = 0.0;
	P[0] = L[0]; 	N[0] = 1.0;
	FT_CutSurfBdry(surf,plane_constr_func,(POINTER)&plane_params,NULL,0,2);

	for (i = 0; i < 3; ++i)
	    P[i] = N[i] = 0.0;
	P[1] = L[1]; 	N[1] = 1.0;
	FT_CutSurfBdry(surf,plane_constr_func,(POINTER)&plane_params,NULL,0,2);

	for (i = 0; i < 3; ++i)
	    P[i] = N[i] = 0.0;
	P[0] = U[0]; 	N[0] = -1.0;
	FT_CutSurfBdry(surf,plane_constr_func,(POINTER)&plane_params,NULL,0,2);

	for (i = 0; i < 3; ++i)
	    P[i] = N[i] = 0.0;
	P[1] = U[1]; 	N[1] = -1.0;
	FT_CutSurfBdry(surf,plane_constr_func,(POINTER)&plane_params,NULL,0,2);

	insert_coords[0][0] = FT_Max(constr_params.L1[0],constr_params.L2[0]);
	insert_coords[0][1] = FT_Max(constr_params.L1[1],constr_params.L2[1]);
	insert_coords[1][0] = FT_Min(constr_params.U1[0],constr_params.U2[0]);
	insert_coords[1][1] = FT_Min(constr_params.U1[1],constr_params.U2[1]);
	insert_coords[2][0] = FT_Max(constr_params.L1[0],constr_params.L2[0]);
	insert_coords[2][1] = FT_Min(constr_params.U1[1],constr_params.U2[1]);
	insert_coords[3][0] = FT_Min(constr_params.U1[0],constr_params.U2[0]);
	insert_coords[3][1] = FT_Max(constr_params.L1[1],constr_params.L2[1]);

	FT_CutSurfBdry(surf,xoss_constr_func,(POINTER)&constr_params,
			insert_coords,4,2);
	FT_InstallSurfEdge(surf,MONO_COMP_HSBDRY);
	FT_FreeThese(1,insert_coords);
}	/* end cutToCross */

static void cutToWing(
	FILE *infile,
	SURFACE *surf)
{
}	/* end cutToWing */

static void installStringChord(
	FILE *infile,
	Front *front,
	SURFACE *surf)
{
	char string[100];

	CursorAfterString(infile,"Enter yes to attach strings to canopy:");
        fscanf(infile,"%s",string);
        (void) printf("%s\n",string);
	if (string[0] != 'Y' && string[0] != 'y')
	    return;
	(void) printf("Available string chord configurations are\n");
	(void) printf("\tParallel (P)\n");
	(void) printf("\tAngular (A)\n");
	CursorAfterString(infile,"Enter string configuration type:");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	switch (string[0])
	{
	case 'P':
	case 'p':
	    installParallelString(infile,front,surf);
	    break;
	case 'A':
	case 'a':
	    installAngularString(infile,front,surf);
	    break;
	default:
	    (void) printf("Unknow string configuration type\n");
	    clean_up(ERROR);
	}
}	/* end installStringChord */

static void installParallelString(
	FILE *infile,
	Front *front,
	SURFACE *surf)
{
	INTERFACE *intfc = front->interf;
	CURVE **c,*canopy_bdry,**string_curves;
	POINT **pts,**string_pts;
	BOND **bonds;
	double cen[MAXD],cload[MAXD],coords[MAXD];
	double ave_radius_sqr,max_radius_sqr;
	int i,j,k,num_curves,num_points,num_bonds,num_strings,num_pairs[2];
	int nb;
	double *string_angle,start_angle,d_angle;
	double xs,xe,ys,ye;	// start and end coords for string in each dir
	double xl,xu,yl,yu;	// lower and upper side coords for the string
	double theta1,theta2,d1,d2,rot_theta,rot_phi;
	NODE *nload,**string_nodes;
	boolean node_moved;
	AF_NODE_EXTRA *extra;
	double spacing,dir[MAXD],*h = computational_grid(intfc)->h;
	BOND *b;
	double *U,*L, width;

	printf("Entering installParallelString()\n");
	if (CursorAfterStringOpt(infile,"Enter yes to attach gores on canopy:"))
        {
	}
	CursorAfterString(infile,"Enter number of paris in x-direction:");
        fscanf(infile,"%d",&num_pairs[0]);
        (void) printf("%d\n",num_pairs[0]);
	CursorAfterString(infile,"Enter start and end x-coordinates:");
        fscanf(infile,"%lf %lf",&xs,&xe);
        (void) printf("%f %f\n",xs,xe);
	CursorAfterString(infile,"Enter y position of lower and upper sides:");
        fscanf(infile,"%lf %lf",&yl,&yu);
        (void) printf("%f %f\n",yl,yu);
	
	CursorAfterString(infile,"Enter number of paris in y-direction:");
        fscanf(infile,"%d",&num_pairs[1]);
        (void) printf("%d\n",num_pairs[1]);
	CursorAfterString(infile,"Enter start and end y-coordinates:");
        fscanf(infile,"%lf %lf",&ys,&ye);
        (void) printf("%f %f\n",ys,ye);
	CursorAfterString(infile,"Enter x position of lower and upper sides:");
        fscanf(infile,"%lf %lf",&xl,&xu);
        (void) printf("%f %f\n",xl,xu);

	CursorAfterString(infile,"Enter initial position of load:");
        fscanf(infile,"%lf %lf %lf",&cload[0],&cload[1],&cload[2]);
        (void) printf("%f %f %f\n",cload[0],cload[1],cload[2]);
	cen[0] = cload[0];	cen[1] = cload[1];
	
	num_strings = 2*(num_pairs[0] + num_pairs[1]);

	canopy_bdry = NULL;
	nload = make_node(Point(cload));
	FT_ScalarMemoryAlloc((POINTER*)&extra,sizeof(AF_NODE_EXTRA));
	extra->af_node_type = LOAD_NODE;
	nload->extra = (POINTER)extra;

	FT_VectorMemoryAlloc((POINTER*)&string_angle,num_strings,
				sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&string_pts,num_strings,
				sizeof(POINT*));
	FT_VectorMemoryAlloc((POINTER*)&string_nodes,num_strings,
				sizeof(NODE*));
	FT_VectorMemoryAlloc((POINTER*)&string_curves,num_strings,
				sizeof(CURVE*));

	width = (xe - xs)/(num_pairs[0] - 1.0);
	for (i = 0; i < num_pairs[0]; ++i)
	{
	    coords[0] = xs + i*width;
	    coords[1] = yl;
	    string_angle[2*i] = plane_angle(cen,coords);
	    coords[1] = yu;
	    string_angle[2*i + 1] = plane_angle(cen,coords);
	}
	width = (ye - ys)/(num_pairs[1] - 1.0);
	for (i = num_pairs[0]; i < num_pairs[0]+num_pairs[1]; ++i)
	{
	    coords[0] = xl;
	    coords[1] = ys + (i - num_pairs[0])*width;
	    string_angle[2*i] = plane_angle(cen,coords);
	    coords[0] = xu;
	    string_angle[2*i + 1] = plane_angle(cen,coords);
	}

	for (i = 0; i < num_strings; ++i)
        {
            if (string_angle[i] > 2.0*PI) 
        	string_angle[i] -= 2.0*PI;
        }

	num_curves = FT_NumOfSurfCurves(surf);
	FT_VectorMemoryAlloc((POINTER*)&c,num_curves,sizeof(CURVE*));
	FT_ArrayOfSurfCurves(surf,c);

	max_radius_sqr = 0.0;
	for (i = 0; i < num_curves; ++i)
	{
	    if (hsbdry_type(c[i]) != MONO_COMP_HSBDRY)
		continue;
	    num_points = FT_NumOfCurvePoints(c[i]);
	    FT_VectorMemoryAlloc((POINTER*)&pts,num_points,sizeof(POINT*));
	    FT_ArrayOfCurvePoints(c[i],pts);
	    ave_radius_sqr = 0.0;
	    for (j = 0; j < num_points; ++j)
	    {
		ave_radius_sqr += sqr(Coords(pts[j])[0] - cen[0]) +
			      sqr(Coords(pts[j])[1] - cen[1]);
	    }
	    ave_radius_sqr /= (double)num_points;
	    if (ave_radius_sqr > max_radius_sqr)
	    {
		max_radius_sqr = ave_radius_sqr;
		canopy_bdry = c[i];
	    }
	    FT_FreeThese(1,pts);
	}
	FT_FreeThese(1,c);

	num_bonds = FT_NumOfCurveBonds(canopy_bdry);
	FT_VectorMemoryAlloc((POINTER*)&bonds,num_bonds,sizeof(BOND*));
	FT_ArrayOfCurveBonds(canopy_bdry,bonds);
	for (i = 0; i < num_bonds; ++i)
	{
	    theta1 = plane_angle(cen,Coords(bonds[i]->start));
	    theta2 = plane_angle(cen,Coords(bonds[i]->end));
	    if (fabs(theta1 - theta2) > PI) 
	    {
		if (theta2 > theta1) theta2 -= 2*PI;
		else theta1 -= 2*PI;
	    }
	    for (j = 0; j < num_strings; ++j)
	    {
		if (within_interval(theta1,theta2,string_angle[j]) ||
		    within_interval(theta1,theta2,string_angle[j]-2*PI))
		{
	    	    d1 = distance_between_positions(cen,
					Coords(bonds[i]->start),2);
	    	    d2 = distance_between_positions(cen,
					Coords(bonds[i]->end),2);
		    d1 = 0.5*(d1 + d2);
		    coords[0] = cen[0] + d1*cos(string_angle[j]);
		    coords[1] = cen[1] + d1*sin(string_angle[j]);
		    coords[2] = 0.5*(Coords(bonds[i]->start)[2] + 
					Coords(bonds[i+1]->end)[2]);
		    string_pts[j] = Point(coords);
	    	    insert_point_in_bond(string_pts[j],bonds[i],canopy_bdry);
		} 
	    }
	}
	FT_FreeThese(1,bonds);

	node_moved = NO;
	for (i = 0; i < num_strings; ++i)
	{
	    canopy_bdry = FT_CurveOfPoint(intfc,string_pts[i],&b);
	    if (is_closed_curve(canopy_bdry) && !node_moved)
	    {
		move_closed_loop_node(canopy_bdry,b);
	        string_nodes[i] = FT_NodeOfPoint(intfc,string_pts[i]);
	    	FT_ScalarMemoryAlloc((POINTER*)&extra,sizeof(AF_NODE_EXTRA));
	    	extra->af_node_type = STRING_NODE;
	    	string_nodes[i]->extra = (POINTER)extra;
		node_moved = YES;
		continue;
	    }
	    c = split_curve(string_pts[i],b,canopy_bdry,0,0,0,0);
	    string_nodes[i] = FT_NodeOfPoint(intfc,string_pts[i]);
	    FT_ScalarMemoryAlloc((POINTER*)&extra,sizeof(AF_NODE_EXTRA));
	    extra->af_node_type = STRING_NODE;
	    string_nodes[i]->extra = (POINTER)extra;
	}
	for (i = 0; i < num_strings; ++i)
	{
	    string_curves[i] = make_curve(0,0,string_nodes[i],nload);
	    hsbdry_type(string_curves[i]) = STRING_HSBDRY;
	    spacing = separation(string_nodes[i]->posn,nload->posn,3);
	    for (j = 0; j < 3; ++j)
		dir[j] = (Coords(nload->posn)[j] - 
			Coords(string_nodes[i]->posn)[j])/spacing;
	    nb = (int)spacing/(1.1*h[0]);
	    spacing /= (double)nb;
	    b = string_curves[i]->first;
	    for (j = 1; j < nb; ++j)
	    {
	    	for (k = 0; k < 3; ++k)
		    coords[k] = Coords(string_nodes[i]->posn)[k] + 
					j*dir[k]*spacing;
		insert_point_in_bond(Point(coords),b,string_curves[i]);
		b = b->next;
	    }
	}
	FT_FreeThese(4,string_angle,string_pts,string_nodes,string_curves);
}	/* end installParallelString */

static void installAngularString(
	FILE *infile,
	Front *front,
	SURFACE *surf)
{
}	/* end installAngularString */

static void initParachuteDefault(
	Front *front)
{
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	af_params->is_parachute_system = YES;
	af_params->num_opt_round = 20;
        af_params->spring_model = MODEL1;
	af_params->gore_len_fac = 1.0;
}	/* end initParachuteDefault */

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

typedef struct {
        int dim;
        double cen[MAXD];
	double ang[MAXD];
        double rad;
} ARC_PARAMS;



static void initSingleModule(Front*);
static void initMultiModule(Front*,int);
static void initRigidBody(Front*);
static void MergeTwoIntfc(INTERFACE*,INTERFACE*);
static void CopyNodeInfo(INTERFACE*,INTERFACE*);
static void modifyCanopySet(FILE*,Front*,SURFACE*);
static boolean curve_of_boundary_hs(CURVE*);
static boolean line_seg_func(POINTER,double,double*);
static boolean arc_func(POINTER,double,double*);
static void rgb_init(Front*,RG_PARAMS);


extern void init2DModules(Front *front)
{
	int i,num_canopy;
	FILE *infile = fopen(InName(front),"r");
	SURFACE *surf;
	RG_PARAMS rgb_params;

	if (debugging("trace"))
	    (void) printf("Entering init2DModules()\n");

	if (debugging("set_module"))
	    gview_plot_interface("module-step-1",front->interf);

	num_canopy = 0;
	if (CursorAfterStringOpt(infile,"Enter number of canopy curve:"))
	{
            fscanf(infile,"%d",&num_canopy);
	    fclose(infile);
            (void) printf("%d\n",num_canopy);
	}

	if (num_canopy == 1)
	    initSingleModule(front);
	else if (num_canopy > 1)
	    initMultiModule(front,num_canopy);

	initRigidBody(front);
	rgb_init(front,rgb_params);

	if (debugging("trace"))
	    (void) printf("Leaving init2DModules()\n");
}	/* end init2DModules */

extern void initParachuteDefault(
	Front *front)
{
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	af_params->is_parachute_system = YES;
	af_params->num_opt_round = 20;
        af_params->spring_model = MODEL1;
	af_params->gore_len_fac = 1.0;
}	/* end initParachuteDefault */

static void initSingleModule(
        Front *front)
{
	FILE *infile = fopen(InName(front),"r");
	SURFACE *surf;

	/*
	CgalCanopySurface(infile,front,&surf);
	InstallNewLoadNode(front,1);
	*/
	fclose(infile);
}	/* end initSingleModule */

static void initMultiModule(
	Front *front,
	int num_canopy)
{
	INTERFACE *intfc = front->interf;
	SURFACE **surfs;
	FILE *infile = fopen(InName(front),"r");
	double center[MAXD];
	double phi,theta;
	char string[10];
	int i;
        INTERFACE *cur_intfc;
        cur_intfc = current_interface();

	/*
	FT_VectorMemoryAlloc((POINTER*)&surfs,num_canopy,sizeof(SURFACE*));

	for (i = 0; i < num_canopy; ++i)
	{
	    CgalCanopySurface(infile,front,&surfs[i]);
	    modifyCanopySet(infile,front,surfs[i]);
	}

	InstallNewLoadNode(front,num_canopy);
	set_current_interface(cur_intfc);
	FT_FreeThese(1,surfs);
	*/
}	/* end initMultiModule */

static void CopyNodeInfo(
	INTERFACE *intfc,
	INTERFACE *newintfc)
{
	AF_NODE_EXTRA *extra, *tmp;
        INTERFACE *cur_intfc;
	NODE **n,**newn;

        cur_intfc = current_interface();
	set_current_interface(newintfc);
	intfc_node_loop(intfc,n)
	{
	    intfc_node_loop(newintfc,newn)
	    {
		if (fabs(Coords((*n)->posn)[0]-Coords((*newn)->posn)[0]) < 1e-6
		 && fabs(Coords((*n)->posn)[1]-Coords((*newn)->posn)[1]) < 1e-6
		 && fabs(Coords((*n)->posn)[2]-Coords((*newn)->posn)[2]) < 1e-6)
		{
		    tmp = (AF_NODE_EXTRA*)((*n)->extra);
		    if (tmp != NULL)
		    {
            	        FT_ScalarMemoryAlloc((POINTER*)&extra,
				sizeof(AF_NODE_EXTRA));
            	    	extra->af_node_type = tmp->af_node_type;
            	        (*newn)->extra = (POINTER)extra;
		    }
		    else
			(*newn)->extra = NULL;
		    break;
		}
	    }
	}
	set_current_interface(cur_intfc);
}	/* end CopyNodeInfo */


static void MergeTwoIntfc(
	INTERFACE *intfc,
	INTERFACE *tmp_intfc)
{
	P_LINK    *p_table;
	int p_size;
	SURFACE   **tmp_s,*news;
	CURVE **tmp_c,*newc;
	NODE **tmp_n,*newn;
	AF_NODE_EXTRA *extra, *tmp_ex;
        INTERFACE *cur_intfc;

        cur_intfc = current_interface();
	set_current_interface(intfc);

	p_size = 4*(tmp_intfc->num_points) + 1;
        uni_array(&p_table,p_size,sizeof(P_LINK));
        reset_hash_table(p_table,p_size);

	for (tmp_s = tmp_intfc->surfaces; tmp_s && *tmp_s; ++tmp_s)
	{
	    if (is_bdry_hs(Hyper_surf(*tmp_s)))
		continue;
	    news = copy_buffer_surface(*tmp_s,p_table,p_size);
	    Hyper_surf_index(news) = Hyper_surf_index((*tmp_s));
	    
	}
        intfc_curve_loop(tmp_intfc,tmp_c)
        {
            if (curve_of_boundary_hs(*tmp_c))
                continue;
	    newc = matching_curve(*tmp_c,p_table,p_size);
	    hsbdry_type(newc) = hsbdry_type(*tmp_c);
	}
	intfc_node_loop(tmp_intfc,tmp_n)
	{
	    tmp_ex = (AF_NODE_EXTRA*)((*tmp_n)->extra);
	    if (NULL == tmp_ex)
		continue;
	    newn = matching_node((*tmp_n),p_table,p_size);
	    FT_ScalarMemoryAlloc((POINTER*)&extra,sizeof(AF_NODE_EXTRA));
	    extra->af_node_type = tmp_ex->af_node_type;
	    newn->extra = (POINTER)extra;	    
	}
	free(p_table);
	set_current_interface(cur_intfc);
}	/* end MergeTwoIntfc */

static boolean curve_of_boundary_hs(
        CURVE *c)
{
        SURFACE **s;
        curve_pos_surf_loop(c,s)
        {
            if (Boundary_hs(Hyper_surf(*s)))
                return YES;
        }
        curve_neg_surf_loop(c,s)
        {
            if (Boundary_hs(Hyper_surf(*s)))
                return YES;
        }
        return NO;
}       /* end curve_of_boundary_hs */

static void modifyCanopySet(
	FILE *infile,
	Front *front,
	SURFACE *canopy)
{
	char string[200];
	double displacement[MAXD];	// Translation displacement
	double center[MAXD];		// Center of rotation
	double phi,theta;		// Spherical angles of rotation
	int i;
	int nc,nn;			// nc, nn: number of curves and nodes;
	CURVE **c,*curves[500];
	NODE **n,*nodes[500];
	TRI *tri;
	BOND *b;
	POINT *p;

        if (CursorAfterStringOpt(infile,
            "Entering yes to modify initialization:"))
        {
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] != 'y' && string[0] != 'Y')
                return;
        }
	else
	    return;

        CursorAfterString(infile,
                "Enter yes for rotation of canopy:");
        fscanf(infile,"%s",string);
        (void) printf("%s\n",string);
        if (string[0] == 'y' || string[0] == 'Y')
        {
            CursorAfterString(infile,"Enter center of rotation:");
            fscanf(infile,"%lf %lf",center,center+1);
            (void) printf("%f %f\n",center[0],center[1]);
            CursorAfterString(infile,"Enter azimuthal and polar angles:");
            fscanf(infile,"%lf %lf",&phi,&theta);
            (void) printf("%f %f\n",phi,theta);
            theta *= PI/180.0;
            phi *= PI/180.0;
        }

	/* Assemble curves and nodes */
	nc = nn = 0;
	surf_pos_curve_loop(canopy,c)
	{
	    if (!pointer_in_list(*c,nc,(POINTER*)curves))
	    {
		curves[nc++] = *c;
		if (!pointer_in_list((*c)->start,nn,(POINTER*)nodes))
		    nodes[nn++] = (*c)->start;
		if (!pointer_in_list((*c)->end,nn,(POINTER*)nodes))
		    nodes[nn++] = (*c)->end;
	    }
	}
	surf_neg_curve_loop(canopy,c)
	{
	    if (!pointer_in_list(*c,nc,(POINTER*)curves))
	    {
		curves[nc++] = *c;
		if (!pointer_in_list((*c)->start,nn,(POINTER*)nodes))
		    nodes[nn++] = (*c)->start;
		if (!pointer_in_list((*c)->end,nn,(POINTER*)nodes))
		    nodes[nn++] = (*c)->end;
	    }
	}
	for (i = 0; i < nn; ++i)
	{
	    node_in_curve_loop(nodes[i],c)
	    {
		if (!pointer_in_list(*c,nc,(POINTER*)curves))
		{
		    curves[nc++] = *c;
		    if (!pointer_in_list((*c)->start,nn,(POINTER*)nodes))
		    	nodes[nn++] = (*c)->start;
		    if (!pointer_in_list((*c)->end,nn,(POINTER*)nodes))
		    	nodes[nn++] = (*c)->end;
		}
	    }
	    node_out_curve_loop(nodes[i],c)
	    {
		if (!pointer_in_list(*c,nc,(POINTER*)curves))
		{
		    curves[nc++] = *c;
		    if (!pointer_in_list((*c)->start,nn,(POINTER*)nodes))
		    	nodes[nn++] = (*c)->start;
		    if (!pointer_in_list((*c)->end,nn,(POINTER*)nodes))
		    	nodes[nn++] = (*c)->end;
		}
	    }
	}
	/* Systematically rotate all points */
	I_SphericalRotateInteriorSurfPoints(canopy,center,phi,theta);
	for (i = 0; i < nc; ++i)
	    I_SphericalRotateInteriorCurvePoints(curves[i],center,phi,theta);
	for (i = 0; i < nn; ++i)
	    I_SphericalRotatePoint(nodes[i]->posn,center,phi,theta,NO);
}	/* end modifyCanopySet */

static void initRigidBody(
	Front *front)
{
	FILE *infile = fopen(InName(front),"r");
	char string[100];
	double cen[MAXD];
	double radius,radii[MAXD];
	double len,dist;
	int w_type, x_type;
	int i, num, k, dim = FT_Dimension();
	int neg_comp,pos_comp;
	CURVE *curve;
	LINE_SEG_PARAMS line_params;
	ARC_PARAMS arc_params;

	if (CursorAfterStringOpt(infile,"Enter yes to add rigid body:"))
	{
	    fscanf(infile,"%s",string);
	    (void) printf("%s\n",string);
	    if (string[0] != 'y' && string[0] != 'Y')
		return;
	}
	else
	    return;

	w_type = MOVABLE_BODY_BOUNDARY;
	(void) printf("Rigid body can be fixed (F) or Movable (M)\n");
	(void) printf("The default is Movable (M)\n");
	if (CursorAfterStringOpt(infile,"Type yes if the rigid body is fixed:"))
	{
	    fscanf(infile,"%s",string);
	    (void) printf("%s\n",string);
	    if (string[0] == 'y' || string[0] == 'Y')
		w_type = NEUMANN_BOUNDARY;
	}

	(void) printf("Available type of rigid body include:\n");
	(void) printf("\tCircle        (c)\n");
	(void) printf("\tLine segment  (l)\n");
	(void) printf("\tLine segments (n)\n");
	(void) printf("\tCircular arc  (a)\n");
	(void) printf("\tTwo arcs      (t)\n");
	(void) printf("\tTwo segments  (s)\n");

	CursorAfterString(infile,"Enter type of rigid body:");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	switch (string[0])
	{
	case 'c':
	case 'C':
	    CursorAfterString(infile,"Enter center of the circle:");
	    fscanf(infile,"%lf %lf",cen,cen+1);
	    (void) printf("%f %f\n",cen[0],cen[1]);
	    CursorAfterString(infile,"Enter radius of the circle:");
	    fscanf(infile,"%lf",&radius);
	    (void) printf("%f\n",radius);
	    for (i = 0; i < dim; ++i) radii[i] = radius;
	    neg_comp = SOLID_COMP;
	    pos_comp = LIQUID_COMP2;
	    FT_MakeEllipticCurve(front,cen,radii,neg_comp,pos_comp,w_type,
				2.0,&curve);
	    node_type(curve->start) = CLOSED_NODE;
	    break;
	case 'L':
	case 'l':
	    CursorAfterString(infile,"Enter start coordinate:");
	    fscanf(infile,"%lf %lf",line_params.coords_start,
				line_params.coords_start+1);
	    (void) printf("%f %f\n",line_params.coords_start[0],
				line_params.coords_start[1]);
	    CursorAfterString(infile,"Enter end coordinate:");
	    fscanf(infile,"%lf %lf",line_params.coords_end,
				line_params.coords_end+1);
	    (void) printf("%f %f\n",line_params.coords_end[0],
				line_params.coords_end[1]);
	    neg_comp = LIQUID_COMP2;
	    pos_comp = LIQUID_COMP2;
	    line_params.dim = 2;
	    curve = FT_MakeParametricCurve(front,neg_comp,pos_comp,w_type,
				line_seg_func,(POINTER)&line_params,2,NO);
	    node_type(curve->start) = node_type(curve->end) = FIXED_NODE;
	    break;
	case 'N':
        case 'n':
            CursorAfterString(infile,"Enter number of segments:");
            fscanf(infile,"%d",&num);
            (void) printf("%d\n",num);
            CursorAfterString(infile,"Enter length of segments:");
            fscanf(infile,"%lf",&len);
            (void) printf("%f\n",len);
            CursorAfterString(infile,"Enter distance between segments:");
            fscanf(infile,"%lf",&dist);
	    (void) printf("%f\n",dist);
            CursorAfterString(infile,
				"Enter start coordinates for first segment:");
            fscanf(infile,"%lf %lf",line_params.coords_start,
                                line_params.coords_start+1);
            (void) printf("%f %f\n",line_params.coords_start[0],
                                line_params.coords_start[1]);
	    for (k = 0; k < num; k++)
	    {
	    	line_params.coords_end[0] =  line_params.coords_start[0] + len;
	    	line_params.coords_end[1] = line_params.coords_start[1];
            	neg_comp = LIQUID_COMP2;
            	pos_comp = LIQUID_COMP2;
            	line_params.dim = 2;
            	curve = FT_MakeParametricCurve(front,neg_comp,pos_comp,w_type,
                                line_seg_func,(POINTER)&line_params,2,NO);
            	node_type(curve->start) = node_type(curve->end) = FIXED_NODE;
	    	line_params.coords_start[0] = line_params.coords_end[0] + dist;            
            }
	    break;
	case 'A':
        case 'a':
            CursorAfterString(infile,"Enter center coordinates:");
            fscanf(infile,"%lf %lf",arc_params.cen,
                                arc_params.cen+1);
            (void) printf("%f %f\n",arc_params.cen[0],
                                arc_params.cen[1]);
            CursorAfterString(infile,"Enter radius:");
            fscanf(infile,"%lf",&(arc_params.rad));
            (void) printf("%f\n",arc_params.rad);
	    CursorAfterString(infile,"Enter start end angles in degrees:");
            fscanf(infile,"%lf %lf",arc_params.ang,
                                arc_params.ang+1);
            (void) printf("%f %f\n",arc_params.ang[0],
                                arc_params.ang[1]);

            neg_comp = LIQUID_COMP2;
            pos_comp = LIQUID_COMP2;
            arc_params.dim = 3;
            curve = FT_MakeParametricCurve(front,neg_comp,pos_comp,w_type,
                                arc_func,(POINTER)&arc_params,3,NO);
            node_type(curve->start) = node_type(curve->end) = FIXED_NODE;
            break;
	case 'S':
        case 's':
            CursorAfterString(infile,
				"Enter start coordinates of first segment:");
            fscanf(infile,"%lf %lf",line_params.coords_start,
				line_params.coords_start+1);
            (void) printf("%f %f\n",line_params.coords_start[0],
                                line_params.coords_start[1]);
            CursorAfterString(infile,"Enter end coordinates of first segment:");
            fscanf(infile,"%lf %lf",line_params.coords_end,
                                line_params.coords_end+1);
            (void) printf("%f %f\n",line_params.coords_end[0],
                                line_params.coords_end[1]);
            neg_comp = LIQUID_COMP2;
            pos_comp = LIQUID_COMP2;
            line_params.dim = 2;
            curve = FT_MakeParametricCurve(front,neg_comp,pos_comp,w_type,
                                line_seg_func,(POINTER)&line_params,2,NO);
            node_type(curve->start) = node_type(curve->end) = FIXED_NODE;

            CursorAfterString(infile,
				"Enter start coordinates of second segment:");
            fscanf(infile,"%lf %lf",line_params.coords_start,
                                line_params.coords_start+1);
            (void) printf("%f %f\n",line_params.coords_start[0],
                                line_params.coords_start[1]);
            CursorAfterString(infile,
				"Enter end coordinates of second segment:");
            fscanf(infile,"%lf %lf",line_params.coords_end,
                                line_params.coords_end+1);
            (void) printf("%f %f\n",line_params.coords_end[0],
                                line_params.coords_end[1]);
            neg_comp = LIQUID_COMP2;
            pos_comp = LIQUID_COMP2;
            line_params.dim = 2;
            curve = FT_MakeParametricCurve(front,neg_comp,pos_comp,w_type,
                                line_seg_func,(POINTER)&line_params,2,NO);
            node_type(curve->start) = node_type(curve->end) = FIXED_NODE;
            break;
	case 'T':
        case 't':
            CursorAfterString(infile,"Enter center coordinates for first arc:");
            fscanf(infile,"%lf %lf",arc_params.cen,
                                arc_params.cen+1);
            (void) printf("%f %f\n",arc_params.cen[0],
                                arc_params.cen[1]);
            CursorAfterString(infile,"Enter radius for first arc:");
            fscanf(infile,"%lf",&(arc_params.rad));
            (void) printf("%f\n",arc_params.rad);
            CursorAfterString(infile,"Enter start end angles for first arc:");
            fscanf(infile,"%lf %lf",arc_params.ang,
                                arc_params.ang+1);
            (void) printf("%f %f\n",arc_params.ang[0],
                                arc_params.ang[1]);
	    
	    neg_comp = LIQUID_COMP2;
            pos_comp = LIQUID_COMP2;
            arc_params.dim = 3;
            curve = FT_MakeParametricCurve(front,neg_comp,pos_comp,w_type,
                                arc_func,(POINTER)&arc_params,3,NO);
            node_type(curve->start) = node_type(curve->end) = FIXED_NODE;
		
	    CursorAfterString(infile,
				"Enter center coordinates for second arc:");
            fscanf(infile,"%lf %lf",arc_params.cen,
                                arc_params.cen+1);
            (void) printf("%f %f\n",arc_params.cen[0],
                                arc_params.cen[1]);
	    CursorAfterString(infile,"Enter radius for second arc:");
            fscanf(infile,"%lf",&(arc_params.rad));
            (void) printf("%f\n",arc_params.rad);
            CursorAfterString(infile,"Enter start end angles for second arc:");
            fscanf(infile,"%lf %lf",arc_params.ang,
                                arc_params.ang+1);
            (void) printf("%f %f\n",arc_params.ang[0],
                                arc_params.ang[1]);

            neg_comp = LIQUID_COMP2;
            pos_comp = LIQUID_COMP2;
            arc_params.dim = 3;
            curve = FT_MakeParametricCurve(front,neg_comp,pos_comp,w_type,
                                arc_func,(POINTER)&arc_params,3,NO);
            node_type(curve->start) = node_type(curve->end) = FIXED_NODE;
            break;

	default:
	    (void) printf("Unknow type of rigid body!\n");
	    clean_up(ERROR);
	}

	fclose(infile);
}	/* end initRigidBody */

static boolean line_seg_func(
	POINTER params,
	double t,
	double *coords)
{
	LINE_SEG_PARAMS *l_params = (LINE_SEG_PARAMS*)params;
	double *coords_start = l_params->coords_start;
	double *coords_end = l_params->coords_end;
	int i,dim = l_params->dim;
	for (i = 0; i < dim; ++i)
	{
	    if (coords_end[i] == coords_start[i])
		coords[i] = coords_start[i];
	    else
	    	coords[i] = coords_start[i] + t*(coords_end[i] - 
				coords_start[i]);
	}
	return YES;
}	/* end line_seg_func */


static boolean arc_func(
        POINTER params,
        double t,
        double *coords)
{
        ARC_PARAMS *l_params = (ARC_PARAMS*)params;
        double *cen = l_params->cen;
        double ang[2];
        double rad = l_params->rad;
        int i,dim = l_params->dim;
	
	ang[0] = PI*l_params->ang[0]/180;
	ang[1] = PI*l_params->ang[1]/180;
        coords[0] = cen[0] + rad*cos(ang[1] - ang[1]*t*(ang[1]-ang[0])/ang[1]);
        coords[1] = cen[1] + rad*sin(ang[1] - ang[1]*t*(ang[1]-ang[0])/ang[1]);
        return YES;
}       /* end arc_func */

extern void initInnerBoundary(
	Front *front,
	LEVEL_FUNC_PACK *level_func_pack)
{
	char string[100];
	FILE *infile = fopen(InName(front),"r");
	double cen[MAXD],R;
	static CIRCLE_PARAMS circle_params;

	if (!CursorAfterStringOpt(infile,"Enter yes to add inner boundary: "))
	    return;
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	if (string[0] != 'y' && string[0] != 'Y')
	    return;

	level_func_pack->pos_component = SOLID_COMP;
	level_func_pack->neg_component = LIQUID_COMP2;
	(void) printf("Available inner boundary types are\n");
	(void) printf("\tCircle (c)\n");
	CursorAfterString(infile,"Entering inner boundary type: ");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	switch(string[0])
	{
	case 'c':
	case 'C':
	    circle_params.dim = 2;
	    CursorAfterString(infile,"Entering center of the circle: ");
	    fscanf(infile,"%lf %lf",circle_params.cen,circle_params.cen+1);
	    (void) printf("%f %f\n",circle_params.cen[0],circle_params.cen[1]);
	    CursorAfterString(infile,"Entering radius of the circle: ");
	    fscanf(infile,"%lf",&circle_params.R);
	    (void) printf("%f\n",circle_params.R);
	    level_func_pack->func_params = (POINTER)&circle_params;
	    level_func_pack->func = level_circle_func;
	    level_func_pack->wave_type = NEUMANN_BOUNDARY;
	    break;
	default:
	    (void) printf("Unknown inner boundary type\n");
	    clean_up(ERROR);
	}
}	/* end initInnerBoundary */

extern  void prompt_for_rigid_body_params(
        int dim,
        char *inname,
        RG_PARAMS *rgb_params)
{
        int i;
        char msg[100],s[100],ss[100];
        FILE *infile = fopen(inname,"r");
        boolean is_preset_motion = NO;
        double mag_dir;

        if (debugging("rgbody"))
            (void) printf("Enter prompt_for_rigid_body_params()\n");

        rgb_params->dim = dim;
        CursorAfterString(infile,"Type yes if motion is preset: ");
        fscanf(infile,"%s",s);
        (void) printf("%s\n",s);
        if (s[0] == 'y' || s[0] == 'Y')
        {
            (void) printf("Available preset motion types are:\n");
            (void) printf("\tPRESET_TRANSLATION\n");
            (void) printf("\tPRESET_ROTATION\n");
            (void) printf("\tPRESET_MOTION (general)\n");
            CursorAfterString(infile,"Enter type of preset motion: ");
            fscanf(infile,"%s",s);
            (void) printf("%s\n",s);
            switch(s[7])
            {
            case 'M':
                rgb_params->motion_type = PRESET_MOTION;
                break;
            case 'C':
                rgb_params->motion_type = PRESET_COM_MOTION;
                break;
            case 'T':
                rgb_params->motion_type = PRESET_TRANSLATION;
                break;
            case 'R':
                rgb_params->motion_type = PRESET_ROTATION;
                break;
            default:
                (void) printf("Unknow type of preset motion!\n");
                clean_up(ERROR);
            }
        }
        else
        {
            (void) printf("Available dynamic motion types are:\n");
            (void) printf("\tFREE_MOTION:\n");
            (void) printf("\tCOM_MOTION (center of mass):\n");
            (void) printf("\tTRANSLATION:\n");
            (void) printf("\tROTATION:\n");
            CursorAfterString(infile,"Enter type of dynamic motion: ");
            fscanf(infile,"%s",s);
            (void) printf("%s\n",s);
            switch(s[0])
            {
            case 'F':
                rgb_params->motion_type = FREE_MOTION;
                break;
            case 'C':
                rgb_params->motion_type = COM_MOTION;
                break;
            case 'T':
                rgb_params->motion_type = TRANSLATION;
                break;
            case 'R':
                rgb_params->motion_type = ROTATION;
                break;
            default:
                (void) printf("Unknow type of motion!\n");
                clean_up(ERROR);
            }
        }
        if (rgb_params->motion_type == TRANSLATION ||
            rgb_params->motion_type == PRESET_TRANSLATION)
        {
            mag_dir = 0.0;
            CursorAfterString(infile,"Enter the direction of motion:");
            for (i = 0; i < dim; ++i)
            {
                fscanf(infile,"%lf",&rgb_params->translation_dir[i]);
                (void) printf("%f ",rgb_params->translation_dir[i]);
                mag_dir += sqr(rgb_params->translation_dir[i]);
            }
            (void) printf("\n");
            mag_dir = sqrt(mag_dir);
            for (i = 0; i < dim; ++i)
                rgb_params->translation_dir[i] /= mag_dir;
        }
        if (rgb_params->motion_type == FREE_MOTION ||
            rgb_params->motion_type == COM_MOTION ||
            rgb_params->motion_type == TRANSLATION)
        {
            sprintf(msg,"Enter the total mass for rigid body:");
            CursorAfterString(infile,msg);
            fscanf(infile,"%lf",&rgb_params->total_mass);
            (void) printf("%f\n",rgb_params->total_mass);
        }
        if (rgb_params->motion_type == FREE_MOTION ||
            rgb_params->motion_type == COM_MOTION ||
            rgb_params->motion_type == TRANSLATION ||
            rgb_params->motion_type == PRESET_MOTION ||
            rgb_params->motion_type == PRESET_TRANSLATION)
        {
            sprintf(msg,"Enter the initial center of mass for rigid body:");
            CursorAfterString(infile,msg);
            for (i = 0; i < dim; ++i)
            {
                fscanf(infile,"%lf",&rgb_params->center_of_mass[i]);
                (void) printf("%f ",rgb_params->center_of_mass[i]);
            }
            (void) printf("\n");
            sprintf(msg,"Enter the initial center of mass velocity:");
            CursorAfterString(infile,msg);
            for (i = 0; i < dim; ++i)
            {
                fscanf(infile,"%lf",&rgb_params->cen_of_mass_velo[i]);
                (void) printf("%f ",rgb_params->cen_of_mass_velo[i]);
            }
            (void) printf("\n");
        }
        if (rgb_params->motion_type == PRESET_ROTATION)
        {
            /* 2D preset rotation is always about the z-axis */
            /* Center of axis is the coordinate of a point on the axis */
            CursorAfterString(infile,"Enter rotation center:");
            for (i = 0; i < dim; ++i)
            {
                fscanf(infile,"%lf",&rgb_params->rotation_cen[i]);
                (void) printf("%f ",rgb_params->rotation_cen[i]);
            }
            (void) printf("\n");
            CursorAfterString(infile,"Enter preset angular velocity:");
            fscanf(infile,"%lf",&rgb_params->angular_velo);
            (void) printf("%f\n",rgb_params->angular_velo);
        }
        if (rgb_params->motion_type == ROTATION)
        {
            CursorAfterString(infile,"Enter rotation center:");
            for (i = 0; i < dim; ++i)
            {
                fscanf(infile,"%lf",&rgb_params->rotation_cen[i]);
                (void) printf("%f ",rgb_params->rotation_cen[i]);
            }
            (void) printf("\n");
        }
        if (rgb_params->motion_type == FREE_MOTION ||
            rgb_params->motion_type == ROTATION)
        {
            CursorAfterString(infile,"Enter the moment of inertial: ");
            fscanf(infile,"%lf",&rgb_params->moment_of_inertial);
            (void) printf("%f\n",rgb_params->moment_of_inertial);
            CursorAfterString(infile,"Enter initial angular velocity: ");
            fscanf(infile,"%lf",&rgb_params->angular_velo);
            (void) printf("%f\n",rgb_params->angular_velo);
        }

        if (debugging("rgbody"))
            (void) printf("Leaving prompt_for_rigid_body_params()\n");
}       /* end prompt_for_rigid_body_params */

extern void set_rgbody_params(
        RG_PARAMS rg_params,
        HYPER_SURF *hs)
{
        int i,dim = rg_params.dim;
        total_mass(hs) = rg_params.total_mass;
        mom_inertial(hs) = rg_params.moment_of_inertial;
        angular_velo(hs) = rg_params.angular_velo;
        motion_type(hs) = rg_params.motion_type;
        surface_tension(hs) = 0.0;
        for (i = 0; i < dim; ++i)
        {
            center_of_mass(hs)[i] = rg_params.center_of_mass[i];
            center_of_mass_velo(hs)[i] =
                                rg_params.cen_of_mass_velo[i];
            rotation_center(hs)[i] =
                                rg_params.rotation_cen[i];
            translation_dir(hs)[i] = rg_params.translation_dir[i];
        }
}       /* end set_rgbody_params */

static void rgb_init(Front *front,
        RG_PARAMS rgb_params)
{
        CURVE **c;
        SURFACE **s;

        if (debugging("trace"))
            printf("Entering rgb_init() \n");
        if (FT_Dimension() == 1) return;
        else if (FT_Dimension() == 2)
        {
            for (c = front->interf->curves; c && *c; ++c)
            {
                if (wave_type(*c) == MOVABLE_BODY_BOUNDARY)
                {
                    prompt_for_rigid_body_params(front->f_basic->dim,
                                front->f_basic->in_name,&rgb_params);
                    body_index(*c) = 0;
                    set_rgbody_params(rgb_params,Hyper_surf(*c));
                }
            }
        }
        else
        {
            for (s = front->interf->surfaces; s && *s; ++s)
            {
                if (wave_type(*s) == MOVABLE_BODY_BOUNDARY)
                {
                    prompt_for_rigid_body_params(front->f_basic->dim,
                                front->f_basic->in_name,&rgb_params);
                    body_index(*s) = 0;
                    set_rgbody_params(rgb_params,Hyper_surf(*s));
                }
            }
        }
        if (debugging("trace"))
            printf("Leaving rgb_init() \n");
}       /* end rgb_init */

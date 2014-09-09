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

static void initSingleModule(Front*);
static void initMultiModule(Front*,int);
static void initRigidBody(Front*);
static void MergeTwoIntfc(INTERFACE*,INTERFACE*);
static void CopyNodeInfo(INTERFACE*,INTERFACE*);
static void modifyCanopySet(FILE*,Front*,SURFACE*);
static boolean curve_of_boundary_hs(CURVE*);
static boolean line_seg_func(POINTER,double,double*);

extern void init2DModules(Front *front)
{
	int i,num_canopy;
	FILE *infile = fopen(InName(front),"r");
	SURFACE *surf;

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
	int w_type;
	int i,dim = FT_Dimension();
	int neg_comp,pos_comp;
	CURVE *curve;
	LINE_SEG_PARAMS line_params;

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
	(void) printf("\tCircle (c)\n");
	(void) printf("\tLine segment (l)\n");
	(void) printf("\tCirculer arc (a)\n");
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
	case 'A':
	case 'a':
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
	    	coords[i] = coords_start[i] + 
			t*(coords_end[i] - coords_start[i]);
	}
	return YES;
}	/* end line_seg_func */

/***************************************************************
FronTier is a set of libraries that implements different types of 
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
#include <rgbody.h>

static void make_rigid_body(Front*,char*,POINTER,
		double(*func)(POINTER,double*),COMPONENT,COMPONENT,int);
static void make_rigid_body2d(Front*,char*,POINTER,
		double(*func)(POINTER,double*),COMPONENT,COMPONENT,int);
static void make_rigid_body3d(Front*,char*,POINTER,
		double(*func)(POINTER,double*),COMPONENT,COMPONENT,int);
static POINTER get_propeller_params(FILE *);
static void initRotorIntfc(Front*,LEVEL_FUNC_PACK*,char*,IF_PROB_TYPE);
static void initWindMillIntfc3d(Front*,LEVEL_FUNC_PACK*,char*,IF_PROB_TYPE);
static void initBeeIntfc3d(Front*,LEVEL_FUNC_PACK*,char*,IF_PROB_TYPE);
static void initApacheIntfc3d(Front*,LEVEL_FUNC_PACK*,char*,IF_PROB_TYPE);
static void initConeIntfc(Front*,LEVEL_FUNC_PACK*,char*,IF_PROB_TYPE);
static void initPressurePump(Front*,LEVEL_FUNC_PACK*,char*,IF_PROB_TYPE);
static void initHumanBody3d(Front*,LEVEL_FUNC_PACK*,char*,IF_PROB_TYPE);
/*TMP*/
static	void insert_surface_tris_into_another(SURFACE*,SURFACE*);

static void initConeIntfc(
        Front *front,
        LEVEL_FUNC_PACK *level_func_pack,
        char *inname,
        IF_PROB_TYPE prob_type)
{
        IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
        FILE *infile = fopen(inname,"r");
        COMPONENT neg_comp;
        COMPONENT pos_comp;
        char msg[100],s[100];
        double  (*func)(POINTER,double*);
        POINTER func_params;
        static CONE_PARAMS cone_params;
        int dim = 3;
	static RG_PARAMS rg_params;

	rg_params.no_fluid = YES;
	front->extra3 = (POINTER)&rg_params;
        CursorAfterString(infile,"Enter the center of the cone:");
        fscanf(infile,"%lf %lf %lf",&cone_params.center[0],
                        &cone_params.center[1],&cone_params.center[2]);
        (void) printf("%f %f %f\n",cone_params.center[0],
                        cone_params.center[1],cone_params.center[2]);
        CursorAfterString(infile,"Enter the slope of the cone:");
        fscanf(infile,"%lf",&cone_params.slope);
        (void) printf("%f \n",cone_params.slope);
        CursorAfterString(infile,"Enter the height of the cone:");
        fscanf(infile,"%lf",&cone_params.height);
        (void) printf("%f \n",cone_params.height);

        func_params = (POINTER)&cone_params;
        func = cone_func;
        neg_comp = LIQUID_COMP2;
        pos_comp = SOLID_COMP;
        make_rigid_body(front,inname,func_params,func,
                                pos_comp,neg_comp,0);
        fclose(infile);
        FT_ParallelExchIntfcBuffer(front);
}

static void initRotorIntfc(
	Front *front,
	LEVEL_FUNC_PACK *level_func_pack,
	char *inname,
	IF_PROB_TYPE prob_type)
{
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	FILE *infile = fopen(inname,"r");
        COMPONENT neg_comp;
        COMPONENT pos_comp;
        char msg[100],s[100];
        double  (*func)(POINTER,double*);
        POINTER func_params;
	static CIRCLE_PARAMS circle_params;
	CURVE **wall,**contact;
	int num_segs;

	circle_params.dim = 2;
	circle_params.add_plan_surf = NO;
	if (prob_type != OPEN_ROTOR)
        {
	    CursorAfterString(infile,"Enter center of the cylinder:");
	    fscanf(infile,"%lf  %lf",&circle_params.cen[0],
				&circle_params.cen[1]);
	    (void) printf("%f %f\n",circle_params.cen[0],circle_params.cen[1]);
	    CursorAfterString(infile,"Enter radius of the cylinder:");
	    fscanf(infile,"%lf",&circle_params.R);
	    (void) printf("%f\n",circle_params.R);
	    func_params = (POINTER)&circle_params;
            func = level_circle_func;
	    neg_comp = (prob_type == ROTOR_ONE_FLUID) ?
				LIQUID_COMP2 : LIQUID_COMP1;
	    pos_comp = SOLID_COMP;
	    wall = (CURVE**)FT_CreateLevelHyperSurfs(front->rect_grid,
			front->interf,neg_comp,pos_comp,func,func_params,
			NEUMANN_BOUNDARY,&num_segs);
	}
	else
        {
            pos_comp = SOLID_COMP;
            neg_comp = LIQUID_COMP2;
        }
	if (prob_type == ROTOR_TWO_FLUID)
	{
	    CursorAfterString(infile,
			"Enter radius of the two fluid interface:");
	    fscanf(infile,"%lf",&circle_params.R);
	    (void) printf("%f\n",circle_params.R);
	    neg_comp = LIQUID_COMP2;
	    pos_comp = LIQUID_COMP1;
	    contact = (CURVE**)FT_CreateLevelHyperSurfs(front->rect_grid,
			front->interf,neg_comp,pos_comp,func,func_params,
			FIRST_PHYSICS_WAVE_TYPE,&num_segs);
	    sprintf(msg,"Enter density and viscosity of fluid 1:");
	    CursorAfterString(infile,msg);
	    fscanf(infile,"%lf %lf",&iFparams->rho1,&iFparams->mu1);
	    (void) printf("%f %f\n",iFparams->rho1,iFparams->mu1);
	    sprintf(msg,"Enter density and viscosity of fluid 2:");
	    CursorAfterString(infile,msg);
	    fscanf(infile,"%lf %lf",&iFparams->rho2,&iFparams->mu2);
	    (void) printf("%f %f\n",iFparams->rho2,iFparams->mu2);
	    neg_comp = LIQUID_COMP2;
	    pos_comp = SOLID_COMP;
	}
	else
	{
	    sprintf(msg,"Enter density and viscosity of the fluid:");
	    CursorAfterString(infile,msg);
	    fscanf(infile,"%lf %lf",&iFparams->rho2,&iFparams->mu2);
	    (void) printf("%f %f\n",iFparams->rho2,iFparams->mu2);
	}

	sprintf(msg,"Enter the shape of the rotor:");
	CursorAfterString(infile,msg);
	fscanf(infile,"%s",s);
	(void) printf("%s\n",s);
	switch(s[0])
	{
	case 'p':
        case 'P':
	    func_params = get_propeller_params(infile);
            func = propeller_func;
	    make_rigid_body(front,inname,func_params,func,
				pos_comp,neg_comp,0);
	    break;
	}
	fclose(infile);
	FT_ParallelExchIntfcBuffer(front);
}	/* end initRotorIntfc */

static void initWindMillIntfc3d(
	Front *front,
	LEVEL_FUNC_PACK *level_func_pack,
	char *inname,
	IF_PROB_TYPE prob_type)
{
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	FILE *infile = fopen(inname,"r");
	INTERFACE *intfc = front->interf;
        RECT_GRID *gr = front->rect_grid;
        COMPONENT neg_comp;
        COMPONENT pos_comp;
        int dim = gr->dim;
        char msg[100];
	SURFACE *wing,*rotor,*support;
	char vtk_name[100],string[100];
	static RG_PARAMS rgb_params;

	sprintf(string,"Enter the vtk file name for wing:");
	CursorAfterString(infile,string);
	fscanf(infile,"%s",vtk_name);
	(void) printf("%s\n",vtk_name);
	neg_comp = SOLID_COMP;
	pos_comp = LIQUID_COMP2;
	read_vtk_surface(intfc,neg_comp,pos_comp,vtk_name,&wing);
	wave_type(wing) = MOVABLE_BODY_BOUNDARY;
	rgb_params.no_fluid = NO;
	front->extra3 = (POINTER)&rgb_params;
	if (debugging("windmill"))
	{
	    (void) printf("Wing surface read from vtk file %s\n",vtk_name);
	    (void) gview_plot_interface("g-wm-wing",intfc);
	}

	sprintf(string,"Enter the vtk file name for rotor:");
        CursorAfterString(infile,string);
        fscanf(infile,"%s",vtk_name);
        (void) printf("%s\n",vtk_name);
        neg_comp = SOLID_COMP;
        pos_comp = LIQUID_COMP2;
        read_vtk_surface(intfc,neg_comp,pos_comp,vtk_name,&rotor);
        wave_type(rotor) = MOVABLE_BODY_BOUNDARY;
        if (debugging("windmill"))
        {
            (void) printf("Rotor surface read from vtk file %s\n",vtk_name);
            (void) gview_plot_interface("g-wm-rotor",intfc);
        }
	insert_surface_tris_into_another(wing,rotor);
	delete_surface(rotor);

	prompt_for_rigid_body_params(dim,inname,&rgb_params);
	rgb_params.dim = 3;
	set_rgbody_params(rgb_params,Hyper_surf(wing));

	if (debugging("windmill"))
            (void) printf("Passed prompt_for_rigid_body()\n");

	sprintf(string,"Enter the vtk file name for support structure:");
	CursorAfterString(infile,string);
	fscanf(infile,"%s",vtk_name);
	(void) printf("%s\n",vtk_name);
	neg_comp = SOLID_COMP;
	pos_comp = LIQUID_COMP2;
	read_vtk_surface(intfc,neg_comp,pos_comp,vtk_name,&support);
	wave_type(support) = NEUMANN_BOUNDARY;
	if (debugging("windmill"))
	{
	    (void) printf("Support surface read from vtk file %s\n",vtk_name);
	    (void) gview_plot_interface("g-wm-full",intfc);
	    if (consistent_interface(intfc))
	    {
		(void) printf("Interface read from vtk files is consistent\n");
		(void) print_interface(intfc);
	    }
	    else
	    {
		(void) printf("Interface read from vtk files not consistent\n");
		clean_up(ERROR);
	    }
	}
	sprintf(msg,"Enter density and viscosity of the fluid:");
	CursorAfterString(infile,msg);
	fscanf(infile,"%lf %lf",&iFparams->rho2,&iFparams->mu2);
	(void) printf("%f %f\n",iFparams->rho2,iFparams->mu2);

	fclose(infile);
	FT_ParallelExchIntfcBuffer(front);
}	/* end initWindMillIntfc3d */

static void initFluidRgb(
	Front *front,
	LEVEL_FUNC_PACK *level_func_pack,
	char *inname,
	IF_PROB_TYPE prob_type)
{
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	FILE *infile = fopen(inname,"r");
        COMPONENT neg_comp;
        COMPONENT pos_comp;
        char msg[100],s[100];
        double  (*func)(POINTER,double*);
        POINTER func_params;

	sprintf(msg,"Enter density and viscosity of the fluid:");
	CursorAfterString(infile,msg);
	fscanf(infile,"%lf %lf",&iFparams->rho2,&iFparams->mu2);
	(void) printf("%f %f\n",iFparams->rho2,iFparams->mu2);

	neg_comp = LIQUID_COMP2;
	pos_comp = SOLID_COMP;

	sprintf(msg,"Enter the shape of the rigid body:");
        CursorAfterString(infile,msg);
        fscanf(infile,"%s",s);
        (void) printf("%s\n",s);
        switch(s[0])
        {
        case 'p':
        case 'P':
	    func_params = get_propeller_params(infile);
            func = propeller_func;
	    make_rigid_body(front,inname,func_params,func,
				pos_comp,neg_comp,0);
	    break;
	}
	fclose(infile);
	FT_ParallelExchIntfcBuffer(front);
}	/* end initFluidRgb by JDKim */

extern void init_moving_bodies(
        Front *front,
        LEVEL_FUNC_PACK *level_func_pack,
        char *inname,
        IF_PROB_TYPE prob_type)
{
        IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
        iFparams->m_comp2 = LIQUID_COMP2;
	switch (prob_type)
        {
	case ROTOR_ONE_FLUID:
	case OPEN_ROTOR:
            iFparams->m_comp1 = SOLID_COMP;
	case ROTOR_TWO_FLUID:
	    initRotorIntfc(front,level_func_pack,inname,prob_type);
	    break;
	case FLUID_RIGID_BODY:
            iFparams->m_comp1 = SOLID_COMP;
	    initFluidRgb(front,level_func_pack,inname,prob_type);
	    break;
	case WINDMILL_3D:
            iFparams->m_comp1 = SOLID_COMP;
	    initWindMillIntfc3d(front,level_func_pack,inname,prob_type);
	    break;
	case BEE_3D:
            iFparams->m_comp1 = SOLID_COMP;
	    initBeeIntfc3d(front,level_func_pack,inname,prob_type);
	    break;
	case HELICOPTER_3D:
            iFparams->m_comp1 = SOLID_COMP;
	    initApacheIntfc3d(front,level_func_pack,inname,prob_type);
	    break;
	case FLUID_SOLID_CONE:
            iFparams->m_comp1 = SOLID_COMP;
            initConeIntfc(front,level_func_pack,inname,prob_type);
            break;
	case PRESSURE_PUMP:
            iFparams->m_comp1 = SOLID_COMP;
            initPressurePump(front,level_func_pack,inname,prob_type);
            break;
	case HUMAN_BODY_3D:
            iFparams->m_comp1 = SOLID_COMP;
            initHumanBody3d(front,level_func_pack,inname,prob_type);
            break;
	default:
	    (void) printf("ERROR: wrong type in init_moving_bodies!\n");
	    clean_up(ERROR);
	}
}	/* init_moving_bodies */

static void make_rigid_body(
	Front *front,
	char *inname,
	POINTER func_params,
	double  (*func)(POINTER,double*),
	COMPONENT neg_comp,
	COMPONENT pos_comp,
	int id)
{
	int dim = front->rect_grid->dim;
	switch (dim)
	{
	case 2:
	    make_rigid_body2d(front,inname,func_params,func,neg_comp,
				pos_comp,id);
	    return;
	case 3:
	    make_rigid_body3d(front,inname,func_params,func,neg_comp,
				pos_comp,id);
	    return;
	}
}	/* end make_rigid_body */

static void make_rigid_body2d(
	Front *front,
	char *inname,
	POINTER func_params,
	double  (*func)(POINTER,double*),
	COMPONENT neg_comp,
	COMPONENT pos_comp,
	int id)
{
	INTERFACE *intfc = front->interf;
	RECT_GRID *gr = front->rect_grid;
	int j,dim = intfc->dim;
	int num_segs;
	CURVE **curves;
	static RG_PARAMS rg_params;

	rg_params.no_fluid = NO;
	front->extra3 = (POINTER)&rg_params;
	curves = (CURVE**)FT_CreateLevelHyperSurfs(gr,intfc,neg_comp,pos_comp,
			func,func_params,MOVABLE_BODY_BOUNDARY,&num_segs);
	if (curves == NULL || num_segs == 0) return;
	prompt_for_rigid_body_params(dim,inname,&rg_params);

	rg_params.dim = 2;
	for (j = 0; j < num_segs; ++j)
	{
	    body_index(curves[j]) = id;
	    set_rgbody_params(rg_params,Hyper_surf(curves[j]));
	}
}	/* end make_rigid_body2d */

static void make_rigid_body3d(
	Front *front,
	char *inname,
	POINTER func_params,
	double  (*func)(POINTER,double*),
	COMPONENT neg_comp,
	COMPONENT pos_comp,
	int id)
{
	INTERFACE *intfc = front->interf;
	RECT_GRID *gr = front->rect_grid;
	int i,dim = intfc->dim;
	int num_segs;
	SURFACE **ps;
	static RG_PARAMS rg_params;

	rg_params.no_fluid = NO;
	front->extra3 = (POINTER)&rg_params;
	ps = (SURFACE**)FT_CreateLevelHyperSurfs(gr,intfc,neg_comp,pos_comp,
			func,func_params,MOVABLE_BODY_BOUNDARY,&num_segs);
	if (ps == NULL || num_segs == 0) return;
	prompt_for_rigid_body_params(dim,inname,&rg_params);
	rg_params.dim = 3;
	for (i = 0; i < num_segs; ++i)
	{
	    body_index(ps[i]) = id;
	    set_rgbody_params(rg_params,Hyper_surf(ps[i]));
	}
}	/* end make_rigid_body3d */

static POINTER get_propeller_params(FILE *infile)
{
        static PROPELLER_PARAMS pparams;

        CursorAfterString(infile,"Enter the center:");
        fscanf(infile,"%lf %lf",&pparams.x[0],&pparams.y[0]);
	(void) printf("%f %f\n",pparams.x[0],pparams.y[0]);
        CursorAfterString(infile,"Enter the number of wings:");
        fscanf(infile,"%d",&pparams.NofW);
	(void) printf("%d\n",pparams.NofW);
        CursorAfterString(infile,"Enter the inner radius:");
        fscanf(infile,"%lf",&pparams.r[0]);
	(void) printf("%f\n",pparams.r[0]);
        CursorAfterString(infile,"Enter the outer radius:");
        fscanf(infile,"%lf",&pparams.r[1]);
	(void) printf("%f\n",pparams.r[1]);
	if (debugging("init_rigid_body"))
	{
	    printf("Propeller shaped rigid body\n");
	    printf("Center: %f %f\n",pparams.x[0],pparams.y[0]);
	    printf("Number of wings: %d\n",pparams.NofW);
	    printf("Inner radius: %f\n",pparams.r[0]);
	    printf("Outer radius: %f\n",pparams.r[1]);
	}
        return (POINTER)&pparams;
}	/* end get_propeller_params */

static	void insert_surface_tris_into_another(
	SURFACE *s1,
	SURFACE *s2)
{
	int i,n,num_tri = s2->num_tri;
	TRI *tri,**tris;
	FT_VectorMemoryAlloc((POINTER*)&tris,num_tri,sizeof(TRI));
	n = 0;
	for (tri = first_tri(s2); !at_end_of_tri_list(tri,s2); tri = tri->next)
	    tris[n++] = tri;
	for (i = 0; i < n; ++i)
	    insert_tri_at_tail_of_list(tris[i],s1);
	FT_FreeThese(1,tris);
	
}	/* end insert_surface_tris_into_another */

static void initBeeIntfc3d(
	Front *front,
	LEVEL_FUNC_PACK *level_func_pack,
	char *inname,
	IF_PROB_TYPE prob_type)
{
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	FILE *infile = fopen(inname,"r");
	INTERFACE *intfc = front->interf;
        RECT_GRID *gr = front->rect_grid;
        COMPONENT neg_comp;
        COMPONENT pos_comp;
        int i,dim = gr->dim;
        char msg[100];
	SURFACE *body_parts[22];
	char string[100];
	char **vtk_name;
	static RG_PARAMS rgb_params;
	char fname[200];

	printf("Entering initBeeIntfc3D()\n");
	FT_MatrixMemoryAlloc((POINTER*)&vtk_name,22,100,sizeof(char));
	sprintf(vtk_name[0],"vtk-bee/ANTENNA01.vtk");
	sprintf(vtk_name[1],"vtk-bee/ANTENNA02.vtk");
	sprintf(vtk_name[2],"vtk-bee/BACK.vtk");
	sprintf(vtk_name[3],"vtk-bee/BODY.vtk");
	sprintf(vtk_name[4],"vtk-bee/HEAD.vtk");
	sprintf(vtk_name[5],"vtk-bee/LEFT-EYE.vtk");
	sprintf(vtk_name[6],"vtk-bee/LOWERLEG01.vtk");
	sprintf(vtk_name[7],"vtk-bee/LOWERLEG02.vtk");
	sprintf(vtk_name[8],"vtk-bee/LOWERLEG03.vtk");
	sprintf(vtk_name[9],"vtk-bee/LOWERLEG04.vtk");
	sprintf(vtk_name[10],"vtk-bee/LOWERLEG05.vtk");
	sprintf(vtk_name[11],"vtk-bee/LOWERLEG06.vtk");
	sprintf(vtk_name[12],"vtk-bee/RIGHT-EYE.vtk");
	sprintf(vtk_name[13],"vtk-bee/UPPERLEG01.vtk");
	sprintf(vtk_name[14],"vtk-bee/UPPERLEG02.vtk");
	sprintf(vtk_name[15],"vtk-bee/UPPERLEG03.vtk");
	sprintf(vtk_name[16],"vtk-bee/UPPERLEG04.vtk");
	sprintf(vtk_name[17],"vtk-bee/UPPERLEG05.vtk");
	sprintf(vtk_name[18],"vtk-bee/UPPERLEG06.vtk");
	sprintf(vtk_name[19],"vtk-bee/WING1.vtk");
	sprintf(vtk_name[20],"vtk-bee/WING2.vtk");

	rgb_params.no_fluid = NO;
	front->extra3 = (POINTER)&rgb_params;
	for (i = 0; i < 21; ++i)
	{
	    neg_comp = SOLID_COMP;
	    pos_comp = LIQUID_COMP2;
	    read_vtk_surface(intfc,neg_comp,pos_comp,vtk_name[i],
				&body_parts[i]);
	    wave_type(body_parts[i]) = MOVABLE_BODY_BOUNDARY;
	}
	if (debugging("bee"))
	{
	    sprintf(fname,"g-bee");
	    (void) gview_plot_interface(fname,intfc);
	}

	//insert_surface_tris_into_another(wing,rotor);
	//delete_surface(rotor);

	//prompt_for_rigid_body_params(dim,inname,&rgb_params);
	//rgb_params.dim = 3;
	//set_rgbody_params(rgb_params,Hyper_surf(wing));

	printf("More code needs to be added!\n");
	printf("Exiting initBeeIntfc3d()\n");
	clean_up(0);
	sprintf(msg,"Enter density and viscosity of the fluid:");
	CursorAfterString(infile,msg);
	fscanf(infile,"%lf %lf",&iFparams->rho2,&iFparams->mu2);
	(void) printf("%f %f\n",iFparams->rho2,iFparams->mu2);

	fclose(infile);
	FT_ParallelExchIntfcBuffer(front);
}	/* end initBeeIntfc3d */

static void initApacheIntfc3d(
	Front *front,
	LEVEL_FUNC_PACK *level_func_pack,
	char *inname,
	IF_PROB_TYPE prob_type)
{
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	FILE *infile = fopen(inname,"r");
	INTERFACE *intfc = front->interf;
        RECT_GRID *gr = front->rect_grid;
        COMPONENT neg_comp;
        COMPONENT pos_comp;
        int i,dim = gr->dim;
        char msg[100];
	SURFACE *body_parts[22];
	char string[100];
	char **vtk_name;
	static RG_PARAMS rgb_params;
	char fname[200];

	printf("Entering initApacheIntfc3d()\n");
	FT_MatrixMemoryAlloc((POINTER*)&vtk_name,3,100,sizeof(char));
	sprintf(vtk_name[0],"vtk-apache/APACHE.vtk");
	sprintf(vtk_name[1],"vtk-apache/APACROT2.vtk");

	rgb_params.no_fluid = NO;
	front->extra3 = (POINTER)&rgb_params;
	for (i = 0; i < 2; ++i)
	{
	    neg_comp = SOLID_COMP;
	    pos_comp = LIQUID_COMP2;
	    read_vtk_surface(intfc,neg_comp,pos_comp,vtk_name[i],
				&body_parts[i]);
	    wave_type(body_parts[i]) = MOVABLE_BODY_BOUNDARY;
	}
	if (debugging("apache"))
	{
	    sprintf(fname,"g-apache");
	    (void) gview_plot_interface(fname,intfc);
	}

	//insert_surface_tris_into_another(wing,rotor);
	//delete_surface(rotor);

	//prompt_for_rigid_body_params(dim,inname,&rgb_params);
	//rgb_params.dim = 3;
	//set_rgbody_params(rgb_params,Hyper_surf(wing));

	printf("More code needs to be added!\n");
	printf("Exiting initApacheIntfc3d()\n");
	clean_up(0);
	sprintf(msg,"Enter density and viscosity of the fluid:");
	CursorAfterString(infile,msg);
	fscanf(infile,"%lf %lf",&iFparams->rho2,&iFparams->mu2);
	(void) printf("%f %f\n",iFparams->rho2,iFparams->mu2);

	fclose(infile);
	FT_ParallelExchIntfcBuffer(front);
}	/* end initApacheIntfc3d */

static void initPressurePump(
	Front *front,
	LEVEL_FUNC_PACK *level_func_pack,
	char *inname,
	IF_PROB_TYPE prob_type)
{
	FILE *infile = fopen(inname,"r");
	int i,num_nodes;
	double **node_coords;
	int neg_comp,pos_comp;
	CURVE *pump;
	static RG_PARAMS rg_params;
	int dim = FT_Dimension();

	printf("Entering initPressurePump()\n");
	rg_params.no_fluid = NO;
	front->extra3 = (POINTER)&rg_params;
	CursorAfterString(infile,"Enter number of node point of pump: ");
	fscanf(infile,"%d",&num_nodes);
	(void) printf("%d\n",num_nodes);
	FT_MatrixMemoryAlloc((POINTER*)&node_coords,num_nodes,MAXD,
				sizeof(double));
	CursorAfterString(infile,"Enter coordinates of node points: ");
	(void) printf("\n");
	for (i = 0; i < num_nodes; ++i)
	{
	    fscanf(infile,"%lf %lf",&node_coords[i][0],&node_coords[i][1]);
	    (void) printf("%f %f\n",node_coords[i][0],node_coords[i][1]);
	}
	neg_comp = SOLID_COMP;
	pos_comp = LIQUID_COMP2;
	pump = FT_MakeNodeArrayCurve(front,num_nodes,node_coords,neg_comp,
			pos_comp,NO,0.75,MOVABLE_BODY_BOUNDARY);

	fclose(infile);

	prompt_for_rigid_body_params(dim,inname,&rg_params);
	set_rgbody_params(rg_params,Hyper_surf(pump));
}	/* end initPressurePump */

static void initHumanBody3d(
        Front *front,
        LEVEL_FUNC_PACK *level_func_pack,
        char *inname,
        IF_PROB_TYPE prob_type)
{
        IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
        FILE *infile = fopen(inname,"r");
        INTERFACE *intfc = front->interf;
        RECT_GRID *gr = front->rect_grid;
        COMPONENT neg_comp;
        COMPONENT pos_comp;
        int dim = gr->dim;
        char msg[100];
        SURFACE *human_body;
        char vtk_name[100],string[100];
        static RG_PARAMS rgb_params;

        sprintf(string,"Enter the vtk file name for human body:");
        CursorAfterString(infile,string);
        fscanf(infile,"%s",vtk_name);
        (void) printf("%s\n",vtk_name);
        neg_comp = SOLID_COMP;
        pos_comp = LIQUID_COMP2;
        read_vtk_surface(intfc,neg_comp,pos_comp,vtk_name,&human_body);
        wave_type(human_body) = MOVABLE_BODY_BOUNDARY;
        rgb_params.no_fluid = YES;
        front->extra3 = (POINTER)&rgb_params;
        if (debugging("human_body"))
        {
            (void) printf("Hunam body surface read from vtk file %s\n",
                                        vtk_name);
            (void) gview_plot_interface("g-human-body",intfc);
        }

        prompt_for_rigid_body_params(dim,inname,&rgb_params);
        rgb_params.dim = 3;
        set_rgbody_params(rgb_params,Hyper_surf(human_body));

        if (debugging("human_body"))
            (void) printf("Passed prompt_for_rigid_body()\n");

        sprintf(msg,"Enter density and viscosity of the fluid:");
        CursorAfterString(infile,msg);
        fscanf(infile,"%lf %lf",&iFparams->rho2,&iFparams->mu2);
        (void) printf("%f %f\n",iFparams->rho2,iFparams->mu2);

        fclose(infile);
        FT_ParallelExchIntfcBuffer(front);
}       /* end initWindMillIntfc3d */

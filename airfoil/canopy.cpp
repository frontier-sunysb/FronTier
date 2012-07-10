/************************************************************************************
FronTier is a set of libraries that implements differnt types of Front Traking algorithms.
Front Tracking is a numerical method for the solution of partial differential equations 
whose solutions have discontinuities.  


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
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

******************************************************************************/

#include <iFluid.h>
#include <airfoil.h>
#include "solver.h"

static void spring_force_at_point1(double*,POINT*,TRI*,SURFACE*,double);
static void spring_force_at_point2(double*,POINT*,TRI*,SURFACE*,double);
static boolean is_pore(Front*,HYPER_SURF_ELEMENT*,double*);
static void compute_canopy_accel(PARACHUTE_SET*,double**,double**,double**);
static void compute_string_accel(PARACHUTE_SET*,double**,double**,double**);
static void assign_canopy_field(PARACHUTE_SET*,double**,double**);
static void assign_string_field(PARACHUTE_SET*,double**,double**);
static void propagate_canopy(PARACHUTE_SET*,double**);
static void propagate_string(PARACHUTE_SET*,double**);
static void coating_mono_hyper_surf2d(Front*);
static void coating_mono_hyper_surf3d(Front*);
static void compute_total_canopy_force2d(Front*,double*,double*);
static void compute_total_canopy_force3d(Front*,double*,double*);
static void compute_center_of_mass_velo(PARACHUTE_SET*);
static void set_canopy_velocity(PARACHUTE_SET*,double**);

#define 	MAX_NUM_RING1		30

static void spring_force_at_point1(
	double *f,
	POINT *p,
	TRI *tri,
	SURFACE *surf,
	double ks)
{
	TRI *tris[MAX_NUM_RING1];
	int i,j,k,nt;
	POINT *p_nb;
	double length0,length,dir[3];
	
	PointAndFirstRingTris(p,Hyper_surf_element(tri),Hyper_surf(surf),
				&nt,tris);
	for (k = 0; k < 3; ++k) f[k] = 0.0;
	for (i = 0; i < nt; ++i)
	{
	    for (j = 0; j < 3; ++j)
	    {
		if (Point_of_tri(tris[i])[j] == p)
		{
		    length0 = tris[i]->side_length0[j];
		    p_nb = Point_of_tri(tris[i])[(j+1)%3];
		    length = separation(p,p_nb,3);
	    	    for (k = 0; k < 3; ++k)
		    {
			dir[k] = (Coords(p_nb)[k] - Coords(p)[k])/length;
			f[k] += ks*(length - length0)*dir[k];
		    }
		    if (is_side_bdry(tris[i],(j+2)%3))
		    {
			(void) printf("Detect boundary "
				"in spring_force_at_point1()\n");
			clean_up(ERROR);
		    }
		}
	    }
	}
}	/* end spring_force_at_point1 */

extern boolean is_string_node(NODE *n)
{
	AF_NODE_EXTRA *af_node_extra;
	if (n->extra == NULL) return NO;
	af_node_extra = (AF_NODE_EXTRA*)n->extra;
	if (af_node_extra->af_node_type == STRING_NODE) return YES;
	return NO;
}	/* end is_string_node */

extern boolean is_load_node(NODE *n)
{
	AF_NODE_EXTRA *af_node_extra;
	if (n->extra == NULL) return NO;
	af_node_extra = (AF_NODE_EXTRA*)n->extra;
	if (af_node_extra->af_node_type == LOAD_NODE) return YES;
	return NO;
}	/* end is_load_node */

extern void coating_mono_hyper_surf(
	Front *front)
{
	int dim = front->rect_grid->dim;
	switch (dim)
	{
	case 2:
	    coating_mono_hyper_surf2d(front);
	    return;
	case 3:
	    coating_mono_hyper_surf3d(front);
	    return;
	}
}	/* end coating_mono_hyper_surf */

static void coating_mono_hyper_surf2d(
	Front *front)
{
	INTERFACE *grid_intfc = front->grid_intfc;
	RECT_GRID *top_grid = &topological_grid(grid_intfc);
	struct Table *T = table_of_interface(grid_intfc);
	COMPONENT *top_comp = T->components;
        COMPONENT          comp;
        INTERFACE          *intfc = front->interf;
	double 		   *L = top_grid->L;
	double 		   *h = top_grid->h;
	double             coords[MAXD];
        double             t[MAXD],p[MAXD],vec[MAXD];
	const double 	   *nor;
	CURVE **c,*immersed_curve;
        HYPER_SURF_ELEMENT *hse;
        HYPER_SURF         *hs;
	BOND *b;
	COMPONENT base_comp;
	int i,index,nb,index_nb,*top_gmax = top_grid->gmax;
	int dim = top_grid->dim;
	int icoords[MAXD],icn[MAXD],smin[MAXD],smax[MAXD];
	GRID_DIRECTION dir[4] = {WEST,EAST,SOUTH,NORTH};

	if (debugging("trace"))
	    (void) printf("Entering coating_mono_hyper_surf2d()\n");

	immersed_curve = NULL;
	for (c = grid_intfc->curves; c && *c; ++c)
	{
	    if (wave_type(*c) == ELASTIC_BOUNDARY)
	    {
		immersed_curve = *c;
		comp = base_comp = negative_component(*c);
		hs = Hyper_surf(immersed_curve);
		break;
	    }
	}
	if (immersed_curve == NULL)
	{
	    (void) printf("ERROR: In coating_mono_hyper_surf3d()"
			 " no immersed_curve found!\n");
	    clean_up(ERROR);
	}

	for (icoords[0] = 1; icoords[0] < top_gmax[0]; ++icoords[0])
	for (icoords[1] = 1; icoords[1] < top_gmax[1]; ++icoords[1])
	{
	    index = d_index(icoords,top_gmax,dim);
	    for (i = 0; i < dim; ++i)
		coords[i] = L[i] + icoords[i]*h[i];
	    if (nearest_interface_point_within_range(coords,comp,grid_intfc,
			NO_BOUNDARIES,hs,p,t,&hse,&hs,3))
	    {
		if (wave_type(hs) != ELASTIC_BOUNDARY) continue;
		b = Bond_of_hse(hse);
	    	t[1] = Coords(b->start)[0] - Coords(b->end)[0]; // t is normal
	    	t[0] = Coords(b->end)[1] - Coords(b->start)[1];
	    	for (i = 0; i < dim; ++i)
		    vec[i] = coords[i] - p[i];
	    	if (scalar_product(vec,t,dim) > 0.0)
		    top_comp[index] = base_comp + 1;
	    	else
		    top_comp[index] = base_comp - 1;
	    }
	}
	negative_component(immersed_curve) = base_comp - 1;
	positive_component(immersed_curve) = base_comp + 1;
	if (debugging("coat_comp"))
	{
	    for (icoords[1] = 1; icoords[1] < top_gmax[1]; ++icoords[1])
	    {
	    	for (icoords[0] = 1; icoords[0] < top_gmax[0]; ++icoords[0])
	    	{
		    index = d_index(icoords,top_gmax,dim);
		    (void) printf("%d",top_comp[index]);
	    	}
	    	(void) printf("\n");
	    }
	}
	if (debugging("trace"))
	    (void) printf("Leaving coating_mono_hyper_surf2d()\n");
}	/* end coating_mono_hyper_surf2d */

static void coating_mono_hyper_surf3d(
	Front *front)
{
	INTERFACE *grid_intfc = front->grid_intfc;
	RECT_GRID *top_grid = &topological_grid(grid_intfc);
	struct Table *T = table_of_interface(grid_intfc);
	COMPONENT *top_comp = T->components;
        COMPONENT          comp;
        INTERFACE          *intfc = front->interf;
	double 		   *L = top_grid->L;
	double 		   *h = top_grid->h;
	double             coords[MAXD];
        double             t[MAXD],p[MAXD],vec[MAXD];
	const double 	   *nor;
	SURFACE **s,*immersed_surf;
        HYPER_SURF_ELEMENT *hse;
        HYPER_SURF         *hs;
	COMPONENT base_comp;
	int i,index,nb,index_nb,*top_gmax = top_grid->gmax;
	int dim = top_grid->dim;
	int icoords[MAXD],icn[MAXD],smin[MAXD],smax[MAXD];
	GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};

	if (debugging("trace"))
	    (void) printf("Entering coating_mono_hyper_surf3d()\n");
	immersed_surf = NULL;
	for (s = grid_intfc->surfaces; s && *s; ++s)
	{
	    if (wave_type(*s) == ELASTIC_BOUNDARY)
	    {
		immersed_surf = *s;
		hs = Hyper_surf(*s);
		comp = base_comp = negative_component(*s);
		break;
	    }
	}
	if (immersed_surf == NULL)
	{
	    (void) printf("ERROR: In coating_mono_hyper_surf3d()"
			 " No immersed_surf found!\n");
	    clean_up(ERROR);
	}

	for (icoords[0] = 1; icoords[0] < top_gmax[0]; ++icoords[0])
	for (icoords[1] = 1; icoords[1] < top_gmax[1]; ++icoords[1])
	for (icoords[2] = 1; icoords[2] < top_gmax[2]; ++icoords[2])
	{
	    index = d_index(icoords,top_gmax,dim);
	    for (i = 0; i < dim; ++i)
		coords[i] = L[i] + icoords[i]*h[i];
	    if (nearest_interface_point_within_range(coords,comp,grid_intfc,
			NO_BOUNDARIES,hs,p,t,&hse,&hs,3))
	    {
	    	nor = Tri_normal(Tri_of_hse(hse));
	    	for (i = 0; i < dim; ++i)
		    vec[i] = coords[i] - p[i];
	    	if (scalar_product(vec,nor,dim) > 0.0)
		    top_comp[index] = base_comp + 1;
	    	else
		    top_comp[index] = base_comp - 1;
	    }
	}
	negative_component(immersed_surf) = base_comp - 1;
	positive_component(immersed_surf) = base_comp + 1;
	if (debugging("coat_comp"))
	{
	    icoords[0] = top_gmax[0]/2;
	    for (icoords[2] = 0; icoords[2] <= top_gmax[2]; ++icoords[2])
	    {
	    	for (icoords[1] = 0; icoords[1] <= top_gmax[1]; ++icoords[1])
	    	{
		    index = d_index(icoords,top_gmax,dim);
		    printf("%d",top_comp[index]);
	    	}
	    	printf("\n");
	    }
	}
	if (debugging("immersed_surf") && front->step%1 == 0)
	{
	    int icrd_nb[MAXD],index_nb,n;
	    POINTER l_state,u_state;
	    double crx_coords[MAXD],crx_nb[MAXD];
	    static double *pl,*pu,*vz,*x;
	    FILE *pfile;
	    char pname[200];

	    n = 0;
	    if (pu == NULL)
	    {
	    	FT_VectorMemoryAlloc((POINTER*)&pu,top_gmax[1],sizeof(double));
	    	FT_VectorMemoryAlloc((POINTER*)&pl,top_gmax[1],sizeof(double));
	    	FT_VectorMemoryAlloc((POINTER*)&vz,top_gmax[1],sizeof(double));
	    	FT_VectorMemoryAlloc((POINTER*)&x,top_gmax[1],sizeof(double));
	    }
	    icrd_nb[0] = icoords[0] = top_gmax[0]/2;
	    for (icoords[1] = 0; icoords[1] <= top_gmax[1]; ++icoords[1])
	    {
	    	icrd_nb[1] = icoords[1];
	        for (icoords[2] = 2; icoords[2] < top_gmax[2]-1; ++icoords[2])
		{
		    index = d_index(icoords,top_gmax,dim);
	    	    icrd_nb[2] = icoords[2] + 1;
		    index_nb = d_index(icrd_nb,top_gmax,dim);
		    if (top_comp[index] != top_comp[index_nb] &&
		  	FT_StateStructAtGridCrossing(front,icoords,UPPER,
                                top_comp[index],&l_state,&hs,crx_coords) &&
			FT_StateStructAtGridCrossing(front,icrd_nb,LOWER,
                                top_comp[index_nb],&u_state,&hs,crx_coords))
		    {
			pl[n] = getStatePres(l_state);
			pu[n] = getStatePres(u_state);
			vz[n] = getStateZvel(l_state);
			x[n] = crx_coords[1];
			n++;
		    }
		}
	    }
	    sprintf(pname,"cpres-%d.xg",front->step);
	    pfile = fopen(pname,"w");
	    fprintf(pfile,"\"Lower pressure\"\n");
	    for (i = 0; i < n; ++i)
		fprintf(pfile,"%f %f\n",x[i],pl[i]);
	    fprintf(pfile,"\n\n\"Upper pressure\"\n");
	    for (i = 0; i < n; ++i)
		fprintf(pfile,"%f %f\n",x[i],pu[i]);
	    fprintf(pfile,"\n\n\"Pressure difference\"\n");
	    for (i = 0; i < n; ++i)
		fprintf(pfile,"%f %f\n",x[i],pl[i]-pu[i]);
	    fclose(pfile);
	    sprintf(pname,"cvelz-%d.xg",front->step);
	    pfile = fopen(pname,"w");
	    fprintf(pfile,"\"Z-velocity\"\n");
	    for (i = 0; i < n; ++i)
		fprintf(pfile,"%f %f\n",x[i],vz[i]);
	    fclose(pfile);
	}
	if (debugging("trace"))
	    (void) printf("Leaving coating_mono_hyper_surf3d()\n");
}	/* end coating_mono_hyper_surf3d */

extern void compute_total_canopy_force(
	Front *front,
	double *pos_force,
	double *neg_force)
{
	int dim = front->rect_grid->dim;
	switch (dim)
	{
	case 2:
	    compute_total_canopy_force2d(front,pos_force,neg_force);
	    return;
	case 3:
	    compute_total_canopy_force3d(front,pos_force,neg_force);
	    return;
	}
}	/* end compute_total_canopy_force */

static void compute_total_canopy_force2d(
	Front *front,
        double *pos_force,
        double *neg_force)
{
		
}	/* end compute_total_canopy_force2d */

static void compute_total_canopy_force3d(
	Front *front,
        double *pos_force,
        double *neg_force)
{
	TRI *tri;
	SURFACE *surf;
	INTERFACE *intfc = front->interf;
	POINT *p;
	STATE *sl,*sr;
	double pres_p,pres_m;
	double area[MAXD];
	int i;
	static FILE *pfile;

	if (debugging("trace"))
	    (void) printf("Entering compute_total_canopy_force3d()\n");
	if (pfile == NULL)
	{
	    pfile = fopen("payload","w");
	    fprintf(pfile,"\"Net lift vs time\"\n");
	}
	for (i = 0; i < 3; ++i)
	    pos_force[i] = neg_force[i] = 0.0;
	next_tri(intfc,NULL,NULL);
	while (next_tri(intfc,&tri,&surf))
	{
	    if (wave_type(surf) != ELASTIC_BOUNDARY)
		continue; 
	    pres_p = pres_m = 0.0;
	    for (i = 0; i < 3; ++i)
	    {
		p = Point_of_tri(tri)[i];
		FT_GetStatesAtPoint(p,Hyper_surf_element(tri),Hyper_surf(surf),
				(POINTER*)&sl,(POINTER*)&sr);
		pres_m += sl->pres;
		pres_p += sr->pres;
		area[i] = Tri_normal(tri)[i];
	    }
	    for (i = 0; i < 3; ++i)
	    {
		pos_force[i] -= pres_p*area[i]/3.0;
		neg_force[i] += pres_m*area[i]/3.0;
	    }
	}
	if (debugging("trace"))
	    (void) printf("Leaving compute_total_canopy_force3d()\n");
}	/* end compute_total_canopy_force3d */

extern int airfoil_velo(
        POINTER params,
        Front *front,
        POINT *p,
        HYPER_SURF_ELEMENT *hse,
        HYPER_SURF *hs,
        double *vel)
{
	/*TMP to be written*/
        vel[0] = vel[1] = vel[2] = 0.0;
	return YES;
}       /* end airfoil_velo */

extern int af_find_state_at_crossing(
        Front *front,
        int *icoords,
        GRID_DIRECTION dir,
        int comp,
        POINTER *state,
        HYPER_SURF **hs,
        double *crx_coords)
{
        boolean status;
	HYPER_SURF_ELEMENT *hse;

        status = FT_StateStructAtGridCrossing2(front,icoords,dir,comp,state,hs,
                                        &hse,crx_coords);
        if (status == NO) 
	    return NO_PDE_BOUNDARY;
        if (wave_type(*hs) == FIRST_PHYSICS_WAVE_TYPE) 
	    return NO_PDE_BOUNDARY;
	if (wave_type(*hs) == ELASTIC_BOUNDARY && 
		is_pore(front,hse,crx_coords))
	    return NO_PDE_BOUNDARY;
	if (wave_type(*hs) == DIRICHLET_BOUNDARY)
	    return DIRICHLET_PDE_BOUNDARY;
        return NEUMANN_PDE_BOUNDARY;
}       /* af_find_state_at_crossing */

static boolean is_pore(
	Front *front,
	HYPER_SURF_ELEMENT *hse,
	double *crx_coords)
{
	INTERFACE *intfc = front->grid_intfc;
	SURFACE **s;
	TRI *tri;
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	double gamma = af_params->gamma;
	static int current_step = -1;

	if (front->rect_grid->dim != 3) return NO;
	if (gamma == 0.0) return NO;
	if (front->step != current_step)
	{
	    double R;
	    current_step = front->step;
	    for (s = intfc->surfaces; s && *s; ++s)
	    {
		if (wave_type(*s) != ELASTIC_BOUNDARY) continue;
		for (tri = first_tri(*s); !at_end_of_tri_list(tri,*s);
				tri = tri->next)
		{
		    R = drand48();
		    if (R < gamma) Tri_index(tri) = 0;
		    else Tri_index(tri) = 1;
		}
	    }
	}
	tri = Tri_of_hse(hse);
	return (Tri_index(tri) == 0) ? YES : NO;
}	/* end is_pore */


extern void fourth_order_parachute_propagate(
	Front *fr,
        PARACHUTE_SET *old_geom_set,
        PARACHUTE_SET *new_geom_set)
{
	static int size = 0;
	static double **x_old,**x_new,**v_old,**v_new,**f_old,**f_new;
        static double **x_mid,**v_mid,**f_mid,**v_end;
	double lambda_s,m_s,lambda_c,m_c;
	AF_PARAMS *af_params = (AF_PARAMS*)fr->extra2;
	int i,j,num_pts;
	int n,n_tan = af_params->n_tan;
	double fr_dt = fr->dt;
	double dt = fr_dt/(double)n_tan;
	int num_strings = old_geom_set->num_strings;
	int n_cps,n_sps;
	double dt_tol;
	double xcom[MAXD],vcom[MAXD];
	double coeff1,coeff2;

	n_cps = FT_NumOfSurfPoints(old_geom_set->canopy);  /* canopy pts */
	n_sps = 1;				/* load node */
	for (i = 0; i < num_strings; ++i)	/* string interior pts */
	    n_sps += FT_NumOfCurvePoints(old_geom_set->string_curves[i]) - 2;

	if (debugging("trace"))
	    (void) printf("Entering fourth_order_parachute_propagate()\n");
	if (debugging("step_size"))
	{
	    double *spfr = Spfr(fr);
	    printf("Before fourth_order_parachute_propagate()\n");
	    for (i = 0; i <= 3; ++i)
		printf("Max front speed(%d) = %f\n",i,spfr[i]);
	}

	num_pts = n_cps + n_sps;

	new_geom_set->ks = old_geom_set->ks = af_params->ks;
	new_geom_set->lambda_s = old_geom_set->lambda_s = af_params->lambda_s;
	new_geom_set->m_s = old_geom_set->m_s = af_params->m_s;
	new_geom_set->kl = old_geom_set->kl = af_params->kl;
	new_geom_set->lambda_c = old_geom_set->lambda_c = af_params->lambda_c;
	new_geom_set->m_c = old_geom_set->m_c = af_params->m_c;
	new_geom_set->n_cps = old_geom_set->n_cps = n_cps;
	new_geom_set->n_sps = old_geom_set->n_sps = n_sps;
	dt_tol = sqrt((af_params->m_s)/(af_params->ks))/10.0;
	if (dt_tol > sqrt((af_params->m_c)/(af_params->kl))/10.0)
	    dt_tol = sqrt((af_params->m_c)/(af_params->kl))/10.0;
	if (debugging("step_size"))
	{
	    (void) printf("ks = %f  m_s = %f  lambda_s = %f\n",
			new_geom_set->ks,
			new_geom_set->m_s,
			new_geom_set->lambda_s);
	    (void) printf("kl = %f  m_c = %f  lambda_c = %f\n",
			new_geom_set->kl,
			new_geom_set->m_c,
			new_geom_set->lambda_c);
	}

	if (dt > dt_tol)
	{
	    n_tan = (int)(fr_dt/dt_tol);
	    dt = fr_dt/(double)n_tan;
	}
	if (debugging("step_size"))
	{
	    (void) printf("fr_dt = %f  dt_tol = %20.14f  dt = %20.14f\n",
				fr_dt,dt_tol,dt);
	    (void) printf("Number of tangential sub-steps = %d\n",n_tan);
	}
	new_geom_set->dt = old_geom_set->dt = dt;

	if (size < num_pts)
	{
	    size = num_pts;
	    if (v_old != NULL)
	    {
		FT_FreeThese(9,v_old,v_new,x_old,x_new,f_old,f_new,
				v_mid,x_mid,f_mid);
	    }
	    FT_MatrixMemoryAlloc((POINTER*)&x_old,size,3,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&v_old,size,3,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&f_old,size,3,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&x_mid,size,3,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&v_mid,size,3,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&f_mid,size,3,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&x_new,size,3,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&v_new,size,3,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&f_new,size,3,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&v_end,size,3,sizeof(double));
	}

	compute_canopy_accel(old_geom_set,f_old,x_old,v_old);
	compute_string_accel(old_geom_set,f_old,x_old,v_old);

	for (n = 0; n < n_tan; ++n)
	{
	    for (i = 0; i < size; ++i)
            for (j = 0; j < 3; ++j)
            {
                x_new[i][j] = x_old[i][j] + dt*v_old[i][j]/6.0;
                v_new[i][j] = v_old[i][j] + dt*f_old[i][j]/6.0;
                x_mid[i][j] = x_old[i][j] + 0.5*v_old[i][j]*dt;
                v_mid[i][j] = v_old[i][j] + 0.5*f_old[i][j]*dt;
            }
            assign_canopy_field(new_geom_set,x_mid,v_mid);
            assign_string_field(new_geom_set,x_mid,v_mid);
	    compute_canopy_accel(new_geom_set,f_mid,x_mid,v_mid);
	    compute_string_accel(new_geom_set,f_mid,x_mid,v_mid);

	    for (i = 0; i < size; ++i)
            for (j = 0; j < 3; ++j)
            {
                x_new[i][j] += dt*v_mid[i][j]/3.0;
                v_new[i][j] += dt*f_mid[i][j]/3.0;
                x_mid[i][j] = x_old[i][j] + 0.5*v_mid[i][j]*dt;
                v_mid[i][j] = v_old[i][j] + 0.5*f_mid[i][j]*dt;
            }
            assign_canopy_field(new_geom_set,x_mid,v_mid);
            assign_string_field(new_geom_set,x_mid,v_mid);
	    compute_canopy_accel(new_geom_set,f_mid,x_mid,v_mid);
	    compute_string_accel(new_geom_set,f_mid,x_mid,v_mid);

	    for (i = 0; i < size; ++i)
            for (j = 0; j < 3; ++j)
            {
                x_new[i][j] += dt*v_mid[i][j]/3.0;
                v_new[i][j] += dt*f_mid[i][j]/3.0;
                x_mid[i][j] = x_old[i][j] + v_mid[i][j]*dt;
                v_mid[i][j] = v_old[i][j] + f_mid[i][j]*dt;
            }
            assign_canopy_field(new_geom_set,x_mid,v_mid);
            assign_string_field(new_geom_set,x_mid,v_mid);
	    compute_canopy_accel(new_geom_set,f_mid,x_mid,v_mid);
	    compute_string_accel(new_geom_set,f_mid,x_mid,v_mid);

	    for (i = 0; i < size; ++i)
            for (j = 0; j < 3; ++j)
            {
                x_new[i][j] += dt*v_mid[i][j]/6.0;
                v_new[i][j] += dt*f_mid[i][j]/6.0;
            }

	    propagate_canopy(new_geom_set,x_new);
	    propagate_string(new_geom_set,x_new);
            assign_canopy_field(new_geom_set,x_new,v_new);
            assign_string_field(new_geom_set,x_new,v_new);
	    if (n != n_tan-1)
	    {
	    	compute_canopy_accel(new_geom_set,f_old,x_old,v_old);
	    	compute_string_accel(new_geom_set,f_old,x_old,v_old);
	    }
	    else
	    {
	    	if (dt == 0.0)
		{
	    	    for (i = 0; i < size; ++i)
            	    for (j = 0; j < 3; ++j)
		    	v_end[i][j] = 0.0;
		}
		else
		{
	    	    for (i = 0; i < size; ++i)
            	    for (j = 0; j < 3; ++j)
		    	v_end[i][j] = (x_new[i][j] - x_old[i][j])/dt;
		}
	    }
	}
	compute_center_of_mass_velo(new_geom_set);
	set_canopy_velocity(new_geom_set,v_end);

	if (debugging("trace"))
	    (void) printf("Leaving fourth_order_parachute_propagate()\n");
}	/* end fourth_order_parachute_propagate */

static void compute_canopy_accel(
	PARACHUTE_SET *geom_set,
	double **f,
	double **x,
	double **v)
{
	int i,n = 0;
	int ns,nb;
	int n_start,n_end;
	Front *fr = geom_set->front;
	AF_PARAMS *af_params = (AF_PARAMS*)fr->extra2;
	void (*compute_node_accel)(PARACHUTE_SET*,NODE*,double**,
                                double**,double **,int*);
        void (*compute_curve_accel)(PARACHUTE_SET*,CURVE*,double**,
                                double**,double **,int*);
        void (*compute_surf_accel)(PARACHUTE_SET*,SURFACE*,double**,
                                double**,double **,int*);

        switch (af_params->spring_model)
        {
        case MODEL1:
            compute_surf_accel = compute_surf_accel1;
            compute_curve_accel = compute_curve_accel1;
            compute_node_accel = compute_node_accel1;
            break;
        case MODEL2:
            compute_surf_accel = compute_surf_accel2;
            compute_curve_accel = compute_curve_accel2;
            compute_node_accel = compute_node_accel2;
            break;
        case MODEL3:
        default:
            (void) printf("Model function not implemented yet!\n");
            clean_up(ERROR);
        }

	ns = geom_set->num_strings;
	nb = geom_set->num_elas_bdry_curves;

	if (debugging("string_chord") || debugging("rigid_canopy") || debugging("ave_lift"))
	    n_start = n;
	compute_surf_accel(geom_set,geom_set->canopy,f,x,v,&n);
	for (i = 0; i < ns; ++i)
	    compute_node_accel(geom_set,geom_set->string_node[i],f,x,v,&n);
	for (i = 0; i < nb; ++i)
	{
	    compute_curve_accel(geom_set,geom_set->elas_bdry_curves[i],f,x,v,&n);
	    if (is_closed_curve(geom_set->elas_bdry_curves[i]))
		compute_node_accel(geom_set,geom_set->elas_bdry_curves[i]->start,
				f,x,v,&n);	
	}
	if (debugging("string_chord") || debugging("rigid_canopy") || debugging("ave_lift"))
	{
	    n_end = n;
	    for (i = n_start; i < n_end; ++i)
	    {
		f[i][0] = f[i][1] = f[i][2] = 0.0;
		if (debugging("string_chord"))
		    v[i][0] = v[i][1] = v[i][2] = 0.0;
	    }
	}
}	/* end compute_canopy_accel */

static void compute_string_accel(
	PARACHUTE_SET *geom_set,
	double **f,
	double **x,
	double **v)
{
	int i,n;
	Front *fr = geom_set->front;
	int dim = fr->rect_grid->dim;
	IF_PARAMS *iFparams = (IF_PARAMS*)fr->extra1;
	double g = iFparams->gravity[dim-1];
	int ns = geom_set->num_strings;
	AF_PARAMS *af_params = (AF_PARAMS*)fr->extra2;
	void (*compute_node_accel)(PARACHUTE_SET*,NODE*,double**,
                                double**,double **,int*);
        void (*compute_curve_accel)(PARACHUTE_SET*,CURVE*,double**,
                                double**,double **,int*);

        switch (af_params->spring_model)
        {
        case MODEL1:
            compute_curve_accel = compute_curve_accel1;
            compute_node_accel = compute_node_accel1;
            break;
        case MODEL2:
            compute_curve_accel = compute_curve_accel2;
            compute_node_accel = compute_node_accel2;
            break;
        case MODEL3:
            compute_curve_accel = compute_curve_accel3;
            compute_node_accel = compute_node_accel3;
            break;
        default:
            (void) printf("Model function not implemented yet!\n");
            clean_up(ERROR);
        }

	if (debugging("trace"))
	    (void) printf("Entering compute_string_accel()\n");

	n = geom_set->n_cps;
	compute_node_accel(geom_set,geom_set->load_node,f,x,v,&n);
	for (i = 0; i < ns; ++i)
	{
	    compute_curve_accel(geom_set,geom_set->string_curves[i],f,x,v,&n);
	}

	if (debugging("trace"))
	    (void) printf("Leaving compute_string_accel()\n");
}	/* end  compute_string_accel */

static void assign_canopy_field(
	 PARACHUTE_SET *geom_set,
        double **x,
        double **v)
{
	int n = 0;
	int i,ns,nb;
	ns = geom_set->num_strings;
	nb = geom_set->num_elas_bdry_curves;

	assign_surf_field(geom_set->canopy,x,v,&n);
	for (i = 0; i < ns; ++i)
	    assign_node_field(geom_set->string_node[i],x,v,&n);
	for (i = 0; i < nb; ++i)
	{
	    assign_curve_field(geom_set->elas_bdry_curves[i],x,v,&n);
	    if (is_closed_curve(geom_set->elas_bdry_curves[i]))
		 assign_node_field(geom_set->elas_bdry_curves[i]->start,
			x,v,&n);
	}
}	/* end assign_canopy_field */

static void assign_string_field(
	 PARACHUTE_SET *geom_set,
        double **x,
        double **v)
{
	int n;
	int i,ns;

	ns = geom_set->num_strings;
	n = geom_set->n_cps;

	if (debugging("trace"))
	    (void) printf("Entering assign_parachute_field()\n");

	assign_node_field(geom_set->load_node,x,v,&n);
	for (i = 0; i < ns; ++i)
	    assign_curve_field(geom_set->string_curves[i],x,v,&n);

	if (debugging("trace"))
	    (void) printf("Leaving assign_parachute_field()\n");
}       /* end assign_string_field */

extern void assign_node_field(
	NODE *node,
	double **x,
	double **v,
	int *n)
{
	int i,dim = Dimension(node->interface);
	CURVE **c;

	for (i = 0; i < dim; ++i)
	{
	    Coords(node->posn)[i] = x[*n][i];
	    node->posn->vel[i] = v[*n][i];
	}
	for (c = node->out_curves; c && *c; ++c)
	    set_bond_length((*c)->first,dim);
	for (c = node->in_curves; c && *c; ++c)
	    set_bond_length((*c)->last,dim);
	(*n)++;
}	/* end assign_node_field */
	
extern void assign_curve_field(
	CURVE *curve,
	double **x,
	double **v,
	int *n)
{
	int i,j,dim = Dimension(curve->interface);
	BOND *b;

	i = *n;
	for (b = curve->first; b != curve->last; b = b->next)
	{
	    for (j = 0; j < dim; ++j)
	    {
	    	Coords(b->end)[j] = x[i][j];
	    	b->end->vel[j] = v[i][j];
	    }
	    set_bond_length(b,dim);
	    i++;
	}
	set_bond_length(curve->first,dim);
	set_bond_length(curve->last,dim);
	*n = i;
}	/* end assign_curve_field */
	
extern void assign_surf_field(
	SURFACE *surf,
	double **x,
	double **v,
	int *n)
{
	int i,j,k;
	TRI *tri;
	POINT *p;
	STATE *sl,*sr;

	unsort_surf_point(surf);
	i = *n;
	for (tri = first_tri(surf); !at_end_of_tri_list(tri,surf); 
			tri = tri->next)
	{
	    for (j = 0; j < 3; ++j)
	    {
		p = Point_of_tri(tri)[j];
		if (sorted(p) || Boundary_point(p)) continue;
		sl = (STATE*)left_state(p);
		sr = (STATE*)right_state(p);
		for (k = 0; k < 3; ++k)
		{
		    Coords(p)[k] = x[i][k];
		    p->vel[k] = v[i][k];
		}
		sorted(p) = YES;
	    	++i;
	    }
	}
	*n = i;
}	/* end assign_surf_field */
	
extern void compute_surf_accel1(
	PARACHUTE_SET *geom_set,
	SURFACE *surf,
	double **f,
	double **x,
	double **v,
	int *n)
{
	int j,k;
	TRI *tri;
	POINT *p;
	int dim = 3;
	double ks = geom_set->ks;
	double m_s = geom_set->m_s;
	double lambda_s = geom_set->lambda_s;

	unsort_surf_point(surf);
	for (tri = first_tri(surf); !at_end_of_tri_list(tri,surf); 
			tri = tri->next)
	{
	    for (j = 0; j < 3; ++j)
	    {
		p = Point_of_tri(tri)[j];
		if (sorted(p) || Boundary_point(p)) continue;
		for (k = 0; k < dim; ++k)
		{
		    x[*n][k] = Coords(p)[k];
		    v[*n][k] = p->vel[k];
		}
	    	spring_force_at_point1(f[*n],p,tri,surf,ks);
		for (k = 0; k < dim; ++k)
		{
		    f[*n][k] -= lambda_s*(v[*n][k]);
		    f[*n][k] /= m_s;
		}
		sorted(p) = YES;
	    	++(*n);
	    }
	}
}	/* end compute_surf_accel1 */

extern void compute_curve_accel1(
	PARACHUTE_SET *geom_set,
	CURVE *curve,
	double **f,
	double **x,
	double **v,
	int *n)
{
	int i,j;
	double x_diff;
	BOND *b;
	double dir[MAXD],len0,vect[MAXD];
	int dim = Dimension(curve->interface);
	double kc,m_c,lambda_c;
	double *V_refs = geom_set->V_refs;

	if (dim == 3)
	{
	    if (hsbdry_type(curve) == STRING_HSBDRY)
	    {
	    	kc = geom_set->kl;
	    	m_c = geom_set->m_c;
	    	lambda_c = geom_set->lambda_c;
	    }
	    else
	    {
	    	kc = geom_set->ks;
	    	m_c = geom_set->m_s;
	    	lambda_c = geom_set->lambda_s;
	    }
	}
	else
	{
	    kc = geom_set->kl;
	    m_c = geom_set->m_c;
	    lambda_c = geom_set->lambda_c;
	}
	i = *n;
	for (b = curve->first; b != curve->last; b = b->next)
	{
	    for (j = 0; j < dim; ++j)
	    {
	    	x[i][j] = Coords(b->end)[j];
	    	v[i][j] = b->end->vel[j];
		f[i][j] = -lambda_c*(v[i][j] - V_refs[j])/m_c;
	    }
	    i++;
	}

	i = *n;
	for (b = curve->first; b != NULL; b = b->next)
	{
	    len0 = bond_length0(b);
	    x_diff = bond_length(b) - len0; 
	    for (j = 0; j < dim; ++j)
	    {
		dir[j] = (Coords(b->end)[j] - Coords(b->start)[j])
				/bond_length(b);
		vect[j] = x_diff*dir[j];
		if (b != curve->first)
		{
	    	    f[i-1][j]   += kc*vect[j]/m_c;
		}
		if (b != curve->last)
		{
	    	    f[i][j] -= kc*vect[j]/m_c;
		}
	    }
	    if (b != curve->last) i++;
	}

	if (dim == 3)
	{
	    POINT *p,*p_nb;
	    BOND_TRI **btris;
	    TRI **tris;
	    int j,k,side,nt;
	    double length0,length;
	    i = *n;
	    for (b = curve->first; b != curve->last; b = b->next)
	    {
		p = b->end;
		for (btris = Btris(b); btris && *btris; ++btris)
		{
		    nt = FT_FirstRingTrisAroundPoint(p,(*btris)->tri,&tris);
		    for (j = 0; j < nt; ++j)
		    {
			for (side = 0; side < 3; ++side)
			{
			    if (p == Point_of_tri(tris[j])[side])
			    {
				if (is_side_bdry(tris[j],side))
				    continue;
				p_nb = Point_of_tri(tris[j])[(side+1)%3];
				length0 = tris[j]->side_length0[side];
				length = separation(p,p_nb,3);
				for (k = 0; k < 3; ++k)
                        	{
                            	    dir[k] = (Coords(p_nb)[k] - 
						Coords(p)[k])/length;
                            	    f[i][k] += kc*(length - length0)*
					dir[k]/m_c;
                        	}
			    }
			}
		    }
		}
		i++;
	    }
	}
	*n = i;
}	/* end compute_curve_accel1 */

extern void compute_node_accel1(
	PARACHUTE_SET *geom_set,
	NODE *node,
	double **f,
	double **x,
	double **v,
	int *n)
{
	CURVE **c;
	BOND *b;
	double x_diff,len0,len,dir[MAXD],vect[MAXD];
	POINT *p,*p_nb;
	int i,j,dim = Dimension(node->interface);
	double ks = geom_set->ks;
	double kl = geom_set->kl;
	double m_s = geom_set->m_s;
	double m_c = geom_set->m_c;
	double *V_refs = geom_set->V_refs;
	double lambda_s = geom_set->lambda_s;
	double lambda_c = geom_set->lambda_c;
	double m_l;

	if (dim == 3)
	{
	    AF_NODE_EXTRA *extra = (AF_NODE_EXTRA*)node->extra;
	    if (extra != NULL && extra->af_node_type == LOAD_NODE)
	    {
	    	Front *front = geom_set->front;
	    	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;

	    	m_l = af_params->payload;
	    }
	}

	for (i = 0; i < dim; ++i)
	{
	    x[*n][i] = Coords(node->posn)[i];
	    v[*n][i] = node->posn->vel[i];
	    f[*n][i] = 0.0;
	}
	for (c = node->out_curves; c && *c; ++c)
	{
	    b = (*c)->first;
	    len0 = bond_length0(b);
	    x_diff = bond_length(b) - len0; 
	    for (j = 0; j < dim; ++j)
	    {
		dir[j] = (Coords(b->end)[j] - Coords(b->start)[j])
				/bond_length(b);
		vect[j] = x_diff*dir[j];
		if (dim == 3)
		{
		    if (is_load_node(node) == YES)
		    	f[*n][j]   += kl*vect[j]/m_l;
		    else if (hsbdry_type(*c) == STRING_HSBDRY)
	    	    	f[*n][j]   += kl*vect[j]/m_s;
		    else if (hsbdry_type(*c) == MONO_COMP_HSBDRY)
	    	    	f[*n][j]   += ks*vect[j]/m_s;
		}
		else
		    f[*n][j]   += kl*vect[j]/m_c;
	    }
	}
	for (c = node->in_curves; c && *c; ++c)
	{
	    b = (*c)->last;
	    len0 = bond_length0(b);
	    x_diff = bond_length(b) - len0; 
	    for (j = 0; j < dim; ++j)
	    {
		dir[j] = (Coords(b->end)[j] - Coords(b->start)[j])
				/bond_length(b);
		vect[j] = x_diff*dir[j];
		if (dim == 3)
		{
		    if (is_load_node(node) == YES)
		    	f[*n][j]   -= kl*vect[j]/m_l;
		    else if (hsbdry_type(*c) == STRING_HSBDRY)
	    	    	f[*n][j]   -= kl*vect[j]/m_s;
		    else if (hsbdry_type(*c) == MONO_COMP_HSBDRY)
	    	    	f[*n][j]   -= ks*vect[j]/m_s;
		}
		else
		    f[*n][j]   -= kl*vect[j]/m_c;
	    }
	}
	if (dim == 3)
	{
	    BOND_TRI **btris;
	    TRI **tris;
	    int k,side,nt,ns;
	    SURFACE *out_surfs[10];

	    ns = 0;
	    for (c = node->out_curves; c && *c; ++c)
	    {
		b = (*c)->first;
		p = b->start;
		for (btris = Btris(b); btris && *btris; ++btris)
		{
		    out_surfs[ns++] = (*btris)->surface;
		    nt = FT_FirstRingTrisAroundPoint(p,(*btris)->tri,&tris);
		    for (j = 0; j < nt; ++j)
		    {
			for (side = 0; side < 3; ++side)
			{
			    if (p == Point_of_tri(tris[j])[side])
			    {
				if (is_side_bdry(tris[j],side))
				    continue;
				p_nb = Point_of_tri(tris[j])[(side+1)%3];
				len0 = tris[j]->side_length0[side];
				len = separation(p,p_nb,3);
	    			x_diff = len - len0; 
				for (k = 0; k < 3; ++k)
                        	{
                            	    dir[k] = (Coords(p_nb)[k] - 
						Coords(p)[k])/len;
                            	    f[*n][k] += ks*x_diff*dir[k]/m_s;
                        	}
			    }
			}
		    }
		}
	    }
	    for (c = node->in_curves; c && *c; ++c)
	    {
		if (is_closed_curve(*c)) continue;
		b = (*c)->last;
		p = b->end;
		for (btris = Btris(b); btris && *btris; ++btris)
		{
		    boolean duplicate_surf = NO;
		    for (j = 0; j < ns; ++j)
			if ((*btris)->surface == out_surfs[j])
			    duplicate_surf = YES;
		    if (duplicate_surf == YES) continue;
		    nt = FT_FirstRingTrisAroundPoint(p,(*btris)->tri,&tris);
		    for (j = 0; j < nt; ++j)
		    {
			for (side = 0; side < 3; ++side)
			{
			    if (p == Point_of_tri(tris[j])[side])
			    {
				if (is_side_bdry(tris[j],side))
				    continue;
				p_nb = Point_of_tri(tris[j])[(side+1)%3];
				len0 = tris[j]->side_length0[side];
				len = separation(p,p_nb,3);
				x_diff = len - len0;
				for (k = 0; k < 3; ++k)
                        	{
                            	    dir[k] = (Coords(p_nb)[k] - 
						Coords(p)[k])/len;
                            	    f[*n][k] += ks*x_diff*dir[k]/m_s;
                        	}
			    }
			}
		    }
		}
	    }
	    if (!is_load_node(node))
	    {
	    	for (i = 0; i < 3; ++i)
	    	    f[*n][i] -= lambda_s*(v[*n][i] - V_refs[i])/m_s;
	    }
	}
	else
	{
	    for (i = 0; i < 3; ++i)
	    	f[*n][i] -= lambda_c*(v[*n][i] - V_refs[i])/m_c;
	}
	(*n)++;
}	/* end compute_node_accel1 */

static void propagate_canopy(
	PARACHUTE_SET *geom_set,
	double **x)
{
	int i,j;
	TRI *tri;
	POINT *p;
	STATE *sl,*sr;
	HYPER_SURF_ELEMENT *hse;
        HYPER_SURF         *hs;
	double *v;
	Front *front = geom_set->front;
	SURFACE *canopy = geom_set->canopy;
	CURVE *curve;
	NODE *node;
	BOND *b;
	double dt = geom_set->dt;
	int dim = front->rect_grid->dim;
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	double *g = iFparams->gravity;
	int ns,n = 0;

	if (debugging("string_chord") || debugging("folding"))
	    return;

	if (debugging("trace"))
	    (void) printf("Entering propagate_canopy()\n");

	unsort_surf_point(canopy);
	hs = Hyper_surf(canopy);
	for (tri = first_tri(canopy); !at_end_of_tri_list(tri,canopy); 
			tri = tri->next)
	{
	    hse = Hyper_surf_element(tri);
	    for (i = 0; i < 3; ++i)
	    {
		p = Point_of_tri(tri)[i];
		if (sorted(p) || Boundary_point(p)) continue;
		FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
		for (j = 0; j < 3; ++j)
		{
	    	    x[n][j] += (sl->Impct[j] + 0.5*g[j]*dt)*dt;
	    	    sr->Impct[j] = sl->Impct[j] = sl->Impct[j] + g[j]*dt; 
		}
		sorted(p) = YES;
	    	++n;
	    }
	}
	ns = geom_set->num_strings;
	for (i = 0; i < ns; ++i)
	{
	    node = geom_set->string_node[i];
	    sl = (STATE*)left_state(node->posn);
	    sr = (STATE*)right_state(node->posn);
	    for (j = 0; j < 3; ++j)
            {
	    	x[n][j] += (sl->Impct[j] + 0.5*g[j]*dt)*dt;
	    	sr->Impct[j] = sl->Impct[j] = sl->Impct[j] + g[j]*dt; 
            }
	    ++n;
	}
	for (i = 0; i < ns; ++i)
        {
	    curve = geom_set->elas_bdry_curves[i];
	    for (b = curve->first; b != curve->last; b = b->next)
	    {
		p = b->end;
		sl = (STATE*)left_state(p);
		sr = (STATE*)right_state(p);
		for (j = 0; j < dim; ++j)
		{
	    	    x[n][j] += (sl->Impct[j] + 0.5*g[j]*dt)*dt;
	    	    sr->Impct[j] = sl->Impct[j] = sl->Impct[j] + g[j]*dt; 
		}
		++n;
	    }
	    if (is_closed_curve(curve))
	    {
	    	node = curve->start;
	    	sl = (STATE*)left_state(node->posn);
	    	sr = (STATE*)right_state(node->posn);
	    	for (j = 0; j < 3; ++j)
            	{
	    	    x[n][j] += (sl->Impct[j] + 0.5*g[j]*dt)*dt;
	    	    sr->Impct[j] = sl->Impct[j] = sl->Impct[j] + g[j]*dt; 
            	}
	    	++n;
	    }
	}
	if (debugging("trace"))
	    (void) printf("Leaving propagate_canopy()\n");
}	/* end propagate_canopy */

static void propagate_string(
	PARACHUTE_SET *geom_set,
	double **x)
{
	int i,j;
	POINT *p;
	CURVE *c;
	BOND *b;
	STATE *sl,*sr;
	Front *front = geom_set->front;
	NODE *load_node = geom_set->load_node;
	double dt = geom_set->dt;
	int dim = front->rect_grid->dim;
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	double *g = iFparams->gravity;
	double area_dens = af_params->area_dens;
	double ext_force[MAXD];
	int n,ns;

	if (debugging("trace"))
	    (void) printf("Entering propagate_string()\n");
	n = geom_set->n_cps;
	ns = geom_set->num_strings;

	sl = (STATE*)left_state(load_node->posn);
	sr = (STATE*)right_state(load_node->posn);
	for (j = 0; j < dim; ++j)
	{
	    x[n][j] += (sl->Impct[j] + 0.5*g[j]*dt)*dt;
	    sr->Impct[j] = sl->Impct[j] = sl->Impct[j] + g[j]*dt; 
	}
	++n;
	if (debugging("folding"))
	    return;

	for (i = 0; i < ns; ++i)
	{
	    c = geom_set->string_curves[i];
	    for (b = c->first; b != c->last; b = b->next)
	    {
		p = b->end;
		sl = (STATE*)left_state(p);
		sr = (STATE*)right_state(p);
		for (j = 0; j < dim; ++j)
		{
	    	    x[n][j] += (sl->Impct[j] + 0.5*g[j]*dt)*dt;
	    	    sr->Impct[j] = sl->Impct[j] = sl->Impct[j] + g[j]*dt; 
		}
		++n;
	    }
	}
	if (debugging("trace"))
	    (void) printf("Leaving propagate_string()\n");
}	/* end propagate_string */

static void compute_center_of_mass_velo(
	PARACHUTE_SET *geom_set)
{
	int i,j;
	TRI *tri;
	POINT *p;
	STATE *state;
	Front *front = geom_set->front;
	SURFACE *canopy = geom_set->canopy;
	NODE *node;
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	double area_dens = af_params->area_dens;
	double xt[MAXD],vt[MAXD],xcan[MAXD],vcan[MAXD],xload[MAXD],vload[MAXD];
	double area,mass_canopy,payload;
	double *xcom,*vcom;

	if (debugging("trace"))
	    (void) printf("Entering compute_center_of_mass_velo()\n");

	for (j = 0; j < 3; ++j)
	    vcan[j] = 0.0;
	area = mass_canopy = 0.0;
	for (tri = first_tri(canopy); !at_end_of_tri_list(tri,canopy); 
			tri = tri->next)
	{
	    for (j = 0; j < 3; ++j)
	    {
		vt[j] = 0.0;
		xt[j] = 0.0;
	    }
	    for (i = 0; i < 3; ++i)
	    {
		p = Point_of_tri(tri)[i];
		state = (STATE*)left_state(p);
		for (j = 0; j < 3; ++j)
		{
		    vt[j] += state->vel[j]/3.0;
		    xt[j] += Coords(p)[j]/3.0;
		}
	    }
	    for (j = 0; j < 3; ++j)
	    {
		vcan[j] += vt[j]*tri_area(tri);
		xcan[j] += xt[j]*tri_area(tri);
	    }
	    area += tri_area(tri);
	}
	mass_canopy += area_dens*area;
	for (j = 0; j < 3; ++j)
	{
	    vcan[j] /= area;
	    xcan[j] /= area;
	}

	node = geom_set->load_node;
	state = (STATE*)left_state(node->posn);
	for (j = 0; j < 3; ++j)
	{
	    vload[j] = state->vel[j];
	    xload[j] = Coords(node->posn)[j];
	}
	payload = af_params->payload;

	xcom = center_of_mass(Hyper_surf(canopy));
	vcom = center_of_mass_velo(Hyper_surf(canopy));
	for (j = 0; j < 3; ++j)
	{
	    vcom[j] = (vcan[j]*mass_canopy + vload[j]*payload)/
				(mass_canopy + payload);
	    xcom[j] = (xcan[j]*mass_canopy + xload[j]*payload)/
				(mass_canopy + payload);
	}
	if (debugging("trace"))
	    (void) printf("Leaving compute_center_of_mass_velo()\n");
}	/* end compute_center_of_mass_velo */

static void set_canopy_velocity(
	PARACHUTE_SET *geom_set,
	double **v)
{
	int i,j;
	TRI *tri;
	BOND *b;
	POINT *p;
	BOND_TRI **btris;
	STATE *sl,*sr;
	HYPER_SURF_ELEMENT *hse;
        HYPER_SURF         *hs;
	Front *front = geom_set->front;
	SURFACE *canopy = geom_set->canopy;
	CURVE **c,*curve;
	NODE *node,*load_node = geom_set->load_node;
	int dim = front->rect_grid->dim;
	double nor[MAXD],nor_speed,max_speed;
	double *vel;
	int n,ns,nb;

	if (debugging("trace"))
	    (void) printf("Entering set_canopy_velocity()\n");

	if (debugging("string_chord") || debugging("folding"))
	    return;
	if (debugging("step_size"))
	{
	    double *spfr = Spfr(front);
	    (void) printf("Before set_canopy_velocity()\n");
	    for (i = 0; i < dim; ++i)
            {
                (void) printf("front: spfr[%d] %g\n",i,spfr[i]);
            }
            (void) printf("front: spfr[%d] %g\n",i,spfr[i]);
	}

	unsort_surf_point(canopy);
	hs = Hyper_surf(canopy);
	max_speed = 0.0;
	n = 0;
	for (tri = first_tri(canopy); !at_end_of_tri_list(tri,canopy); 
			tri = tri->next)
	{
	    hse = Hyper_surf_element(tri);
	    for (i = 0; i < 3; ++i)
	    {
		p = Point_of_tri(tri)[i];
		if (sorted(p) || Boundary_point(p)) continue;
		FT_NormalAtPoint(p,front,nor,NO_COMP);
		FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
		vel = v[n];
		nor_speed = scalar_product(vel,nor,dim);
		for (j = 0; j < 3; ++j)
		{
		    sl->vel[j] = nor_speed*nor[j];
		    sr->vel[j] = nor_speed*nor[j];
	    	    FT_RecordMaxFrontSpeed(j,sl->vel[j],NULL,Coords(p),front);
		}
	    	FT_RecordMaxFrontSpeed(3,Mag3d(sl->vel),NULL,Coords(p),front);
		if (max_speed < Mag3d(sl->vel)) max_speed = Mag3d(sl->vel);
		sorted(p) = YES;
		n++;
	    }
	}

	ns = geom_set->num_strings;
	for (i = 0; i < ns; ++i)
	{
	    node = geom_set->string_node[i];
	    for (c = node->out_curves; c && *c; ++c)
            {
		if (hsbdry_type(*c) != MONO_COMP_HSBDRY) continue;
                b = (*c)->first;
                p = b->start;
                btris = Btris(b);
                p->hse = hse = Hyper_surf_element((*btris)->tri);
                p->hs = hs = Hyper_surf((*btris)->surface);
                FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
		FT_NormalAtPoint(p,front,nor,NO_COMP);
		vel = v[n];
		nor_speed = scalar_product(vel,nor,dim);
		if (max_speed < fabs(nor_speed)) max_speed = fabs(nor_speed);
                for (j = 0; j < dim; ++j)
		    sl->vel[j] = sr->vel[j] = nor_speed*nor[j];
            }
            for (c = node->in_curves; c && *c; ++c)
            {
		if (hsbdry_type(*c) != MONO_COMP_HSBDRY) continue;
                b = (*c)->last;
                p = b->end;
                btris = Btris(b);
                p->hse = hse = Hyper_surf_element((*btris)->tri);
                p->hs = hs = Hyper_surf((*btris)->surface);
                FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
		FT_NormalAtPoint(p,front,nor,NO_COMP);
		vel = v[n];
		nor_speed = scalar_product(vel,nor,dim);
		if (max_speed < fabs(nor_speed)) max_speed = fabs(nor_speed);
                for (j = 0; j < dim; ++j)
		    sl->vel[j] = sr->vel[j] = nor_speed*nor[j];
            }
	    n++;
	}
	nb = geom_set->num_elas_bdry_curves;
	for (i = 0; i < nb; ++i)
	{
	    curve = geom_set->elas_bdry_curves[i];
	    for (b = curve->first; b != curve->last; b = b->next)
            {
            	p = b->end;
            	btris = Btris(b);
            	if (btris && *btris)
            	{
                    p->hse = hse = Hyper_surf_element((*btris)->tri);
                    p->hs = hs = Hyper_surf((*btris)->surface);
                    FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
		    FT_NormalAtPoint(p,front,nor,NO_COMP);
		    vel = v[n];
		    nor_speed = scalar_product(vel,nor,dim);
                    for (j = 0; j < dim; ++j)
		    	sl->vel[j] = sr->vel[j] = nor_speed*nor[j];
            	}
            	n++;
            }
	}


	if (debugging("step_size"))
	{
	    double *spfr = Spfr(front);
	    (void) printf("In set_canopy_velocity(): max_speed = %f\n",
				max_speed);
	    (void) printf("After set_canopy_velocity()\n");
	    for (i = 0; i < dim; ++i)
            {
                (void) printf("front: spfr[%d] %g\n",i,spfr[i]);
            }
            (void) printf("front: spfr[%d] %g\n",i,spfr[i]);
	}
	n = geom_set->n_cps;
	sl = (STATE*)left_state(load_node->posn);
	sr = (STATE*)right_state(load_node->posn);
	vel = v[n];
	for (j = 0; j < 3; ++j)
	{
	    sl->vel[j] = vel[j];
	    sr->vel[j] = vel[j];
	}
	if (debugging("trace"))
	    (void) printf("Leaving set_canopy_velocity()\n");
}	/* end set_canopy_velocity */

void Incompress_Solver_Smooth_Basis::applicationSetComponent(void)
{
	int i,icrd[MAXD],ic;
        int size = (int)cell_center.size();

        // cell center components
        for (i = 0; i < size; i++)
        {
            cell_center[i].comp =
                        getComponent(cell_center[i].icoords);
        }
	if (debugging("set_shifted_states"))
	{
	    printf("Sample component in applicationSetComponent()\n");
	    if (dim == 3)
	    {
            	icrd[0] = top_gmax[0]/2;
            	for (icrd[2] = 0; icrd[2] <= top_gmax[2]; ++icrd[2])
            	{
                    for (icrd[1] = 0; icrd[1] <= top_gmax[1]; ++icrd[1])
                    {
                        ic = d_index(icrd,top_gmax,dim);
                        printf("%d",top_comp[ic]);
                    }
                    printf("\n");
            	}
	    }
        }
}	/* end applicationSetComponent */

static void print_state(POINTER state);
void Incompress_Solver_Smooth_Basis::applicationSetStates(void)
{
	double coords[MAXD];
	int *icoords;
	int i,j,size = (int)cell_center.size();
	int id;
	STATE state;
	int ave_comp;
	double p_intfc[MAXD],t[MAXD];
	HYPER_SURF_ELEMENT *hse;
	HYPER_SURF *hs;
	double dist;
	
	setDomain();
	for (i = 0; i < size; i++)
        {
            icoords = cell_center[i].icoords;
            if (cell_center[i].comp != -1 &&
                cell_center[i].comp != top_comp[i])
            {
		for (j = 0; j < dim; ++j)
		    coords[j] = top_L[j] + icoords[j]*top_h[j];
		id = d_index(icoords,top_gmax,dim);
		if (fabs(cell_center[i].comp - top_comp[i]) != 2)
		    continue;

		if (debugging("set_crossed_state"))
		{
		    double r;
		    printf("\n");
		    printf("Shifted component:\n");
		    printf("icoords = %d %d %d\n",icoords[0],icoords[1],
					icoords[2]);
		    printf("old comp = %d  new comp = %d\n",
					cell_center[i].comp,top_comp[i]);
		    r = sqrt(sqr(coords[0] - 7.0) + sqr(coords[1] - 7.0));
		    printf("Radius = %f\n",r);
		}

		ave_comp = (cell_center[i].comp + top_comp[i])/2;
		if (!FT_FindNearestIntfcPointInRange(front,ave_comp,
			coords,p_intfc,t,&hse,&hs,2))
		    continue;

		dist = 0.0;
		for (j = 0; j < dim; ++j)
		    dist += sqr(coords[j] - p_intfc[j]);
		dist = sqrt(dist);
		if (debugging("set_crossed_state"))
		{
		    printf("coords  = %f %f %f\n",coords[0],coords[1],
					coords[2]);
		    printf("p_intfc = %f %f %f\n",p_intfc[0],p_intfc[1],
					p_intfc[2]);
		}
		if (dist > top_h[0]*Time_step_factor(front))
		{
		    if (debugging("set_crossed_state"))
			printf("external point: dist = %f\n",dist);
		    continue;
		}

		FrontNearestIntfcState(front,coords,ave_comp,(POINTER)&state);

		if (debugging("set_crossed_state"))
		{
		    printf("Old velocity  : %f %f %f\n",
				cell_center[id].m_state.m_U[0],
				cell_center[id].m_state.m_U[1],
				cell_center[id].m_state.m_U[2]);
		    printf("Intfc velocity: %f %f %f\n",state.vel[0],
			state.vel[1],state.vel[2]);
		    printf("Old pressure   = %f  Old phi   = %f\n",
				cell_center[id].m_state.m_P,
				cell_center[id].m_state.m_phi);
		    printf("Intfc pressure = %f  Intfc phi = %f\n",
				state.pres,state.phi);
		}
		for (j = 0; j < dim; ++j)
		    cell_center[id].m_state.m_U[j] = state.vel[j];
	    }
        }
	FT_FreeGridIntfc(front);
	FT_MakeGridIntfc(front);
}	/* end applicationSetStates */

extern void compute_node_accel2(
	PARACHUTE_SET *geom_set,
	NODE *node,
	double **f,
	double **x,
	double **v,
	int *n)
{
	CURVE **c;
	BOND *b;
	double x_diff,len0,len,dir[MAXD],vect[MAXD];
	POINT *p,*p_nb;
	int i,j,dim = Dimension(node->interface);
	double ks = geom_set->ks;
	double kl = geom_set->kl;
	double m_s = geom_set->m_s;
	double m_c = geom_set->m_c;
	double *V_refs = geom_set->V_refs;
	double lambda_s = geom_set->lambda_s;
	double lambda_c = geom_set->lambda_c;
	double m_l;

	if (dim == 3)
	{
	    AF_NODE_EXTRA *extra = (AF_NODE_EXTRA*)node->extra;
	    if (extra != NULL && extra->af_node_type == LOAD_NODE)
	    {
	    	Front *front = geom_set->front;
	    	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;

	    	m_l = af_params->payload;
	    }
	}

	for (i = 0; i < dim; ++i)
	{
	    x[*n][i] = Coords(node->posn)[i];
	    v[*n][i] = node->posn->vel[i];
	    f[*n][i] = 0.0;
	}
	for (c = node->out_curves; c && *c; ++c)
	{
	    b = (*c)->first;
	    len0 = bond_length0(b);
	    for (j = 0; j < dim; ++j)
	    {
		vect[j] = Coords(b->end)[j] - Coords(b->start)[j] -
				len0*b->dir0[j];
		if (dim == 3)
		{
		    if (is_load_node(node) == YES)
		    	f[*n][j]   += kl*vect[j]/m_l;
		    else if (hsbdry_type(*c) == STRING_HSBDRY)
	    	    	f[*n][j]   += kl*vect[j]/m_s;
		    else if (hsbdry_type(*c) == MONO_COMP_HSBDRY)
	    	    	f[*n][j]   += ks*vect[j]/m_s;
		}
		else
		    f[*n][j]   += kl*vect[j]/m_c;
	    }
	}
	for (c = node->in_curves; c && *c; ++c)
	{
	    b = (*c)->last;
	    len0 = bond_length0(b);
	    for (j = 0; j < dim; ++j)
	    {
		vect[j] = Coords(b->end)[j] - Coords(b->start)[j] -
				len0*b->dir0[j];
		if (dim == 3)
		{
		    if (is_load_node(node) == YES)
		    	f[*n][j]   -= kl*vect[j]/m_l;
		    else if (hsbdry_type(*c) == STRING_HSBDRY)
	    	    	f[*n][j]   -= kl*vect[j]/m_s;
		    else if (hsbdry_type(*c) == MONO_COMP_HSBDRY)
	    	    	f[*n][j]   -= ks*vect[j]/m_s;
		}
		else
		    f[*n][j]   -= kl*vect[j]/m_c;
	    }
	}
	if (dim == 3)
	{
	    BOND_TRI **btris;
	    TRI **tris;
	    int k,side,nt,ns;
	    SURFACE *out_surfs[10];

	    ns = 0;
	    for (c = node->out_curves; c && *c; ++c)
	    {
		b = (*c)->first;
		p = b->start;
		for (btris = Btris(b); btris && *btris; ++btris)
		{
		    out_surfs[ns++] = (*btris)->surface;
		    nt = FT_FirstRingTrisAroundPoint(p,(*btris)->tri,&tris);
		    for (j = 0; j < nt; ++j)
		    {
			for (side = 0; side < 3; ++side)
			{
			    if (p == Point_of_tri(tris[j])[side])
			    {
				if (is_side_bdry(tris[j],side))
				    continue;
				p_nb = Point_of_tri(tris[j])[(side+1)%3];
				len0 = tris[j]->side_length0[side];
				len = separation(p,p_nb,3);
	    			x_diff = len - len0; 
				for (k = 0; k < 3; ++k)
                        	{
                            	    dir[k] = (Coords(p_nb)[k] - 
						Coords(p)[k])/len;
                            	    vect[k] = (Coords(p_nb)[k] - Coords(p)[k])
					- len0*tris[j]->side_dir0[side][k];
                            	    //f[*n][k] += ks*x_diff*dir[k]/m_s;
                            	    f[*n][k] += ks*vect[k]/m_s;
                        	}
			    }
			}
		    }
		}
	    }
	    for (c = node->in_curves; c && *c; ++c)
	    {
		if (is_closed_curve(*c)) continue;
		b = (*c)->last;
		p = b->end;
		for (btris = Btris(b); btris && *btris; ++btris)
		{
		    boolean duplicate_surf = NO;
		    for (j = 0; j < ns; ++j)
			if ((*btris)->surface == out_surfs[j])
			    duplicate_surf = YES;
		    if (duplicate_surf == YES) continue;
		    nt = FT_FirstRingTrisAroundPoint(p,(*btris)->tri,&tris);
		    for (j = 0; j < nt; ++j)
		    {
			for (side = 0; side < 3; ++side)
			{
			    if (p == Point_of_tri(tris[j])[side])
			    {
				if (is_side_bdry(tris[j],side))
				    continue;
				p_nb = Point_of_tri(tris[j])[(side+1)%3];
				len0 = tris[j]->side_length0[side];
				len = separation(p,p_nb,3);
				x_diff = len - len0;
				for (k = 0; k < 3; ++k)
                        	{
                            	    dir[k] = (Coords(p_nb)[k] - 
						Coords(p)[k])/len;
                            	    vect[k] = (Coords(p_nb)[k] - Coords(p)[k])
					- len0*tris[j]->side_dir0[side][k];
                            	    f[*n][k] += ks*vect[k]/m_s;
                        	}
			    }
			}
		    }
		}
	    }
	    if (!is_load_node(node))
	    {
	    	for (i = 0; i < 3; ++i)
	    	    f[*n][i] -= lambda_s*(v[*n][i] - V_refs[i])/m_s;
	    }
	}
	else
	{
	    for (i = 0; i < 3; ++i)
	    	f[*n][i] -= lambda_c*(v[*n][i] - V_refs[i])/m_c;
	}
	(*n)++;
}	/* end compute_node_accel2 */

extern void compute_curve_accel2(
	PARACHUTE_SET *geom_set,
	CURVE *curve,
	double **f,
	double **x,
	double **v,
	int *n)
{
	int i,j;
	double x_diff;
	BOND *b;
	double dir[MAXD],len0,vect[MAXD];
	int dim = Dimension(curve->interface);
	double kc,m_c,lambda_c;
	double *V_refs = geom_set->V_refs;

	if (dim == 3)
	{
	    if (hsbdry_type(curve) == STRING_HSBDRY)
	    {
	    	kc = geom_set->kl;
	    	m_c = geom_set->m_c;
	    	lambda_c = geom_set->lambda_c;
	    }
	    else
	    {
	    	kc = geom_set->ks;
	    	m_c = geom_set->m_s;
	    	lambda_c = geom_set->lambda_s;
	    }
	}
	else
	{
	    kc = geom_set->kl;
	    m_c = geom_set->m_c;
	    lambda_c = geom_set->lambda_c;
	}
	i = *n;
	for (b = curve->first; b != curve->last; b = b->next)
	{
	    for (j = 0; j < dim; ++j)
	    {
	    	x[i][j] = Coords(b->end)[j];
	    	v[i][j] = b->end->vel[j];
		f[i][j] = -lambda_c*(v[i][j] - V_refs[j])/m_c;
	    }
	    i++;
	}

	i = *n;
	for (b = curve->first; b != NULL; b = b->next)
	{
	    len0 = bond_length0(b);
	    for (j = 0; j < dim; ++j)
	    {
		vect[j] = Coords(b->end)[j] - Coords(b->start)[j] -
				len0*b->dir0[j];
		if (b != curve->first)
		{
	    	    f[i-1][j]   += kc*vect[j]/m_c;
		}
		if (b != curve->last)
		{
	    	    f[i][j] -= kc*vect[j]/m_c;
		}
	    }
	    if (b != curve->last) i++;
	}

	if (dim == 3)
	{
	    POINT *p,*p_nb;
	    BOND_TRI **btris;
	    TRI **tris;
	    int j,k,side,nt;
	    double length0,length;
	    i = *n;
	    for (b = curve->first; b != curve->last; b = b->next)
	    {
		p = b->end;
		for (btris = Btris(b); btris && *btris; ++btris)
		{
		    nt = FT_FirstRingTrisAroundPoint(p,(*btris)->tri,&tris);
		    for (j = 0; j < nt; ++j)
		    {
			for (side = 0; side < 3; ++side)
			{
			    if (p == Point_of_tri(tris[j])[side])
			    {
				if (is_side_bdry(tris[j],side))
				    continue;
				p_nb = Point_of_tri(tris[j])[(side+1)%3];
				length0 = tris[j]->side_length0[side];
				length = separation(p,p_nb,3);
				for (k = 0; k < 3; ++k)
                        	{
                            	    dir[k] = (Coords(p_nb)[k] - 
						Coords(p)[k])/length;
				    vect[k] = Coords(p_nb)[k] - Coords(p)[k]
					- length0*tris[j]->side_dir0[side][k];
                            	    //f[i][k] += kc*(length - length0)*
					//dir[k]/m_c;
                            	    f[i][k] += kc*vect[k]/m_c;
                        	}
			    }
			}
		    }
		}
		i++;
	    }
	}
	*n = i;
}	/* end compute_curve_accel2 */

extern void compute_surf_accel2(
	PARACHUTE_SET *geom_set,
	SURFACE *surf,
	double **f,
	double **x,
	double **v,
	int *n)
{
	int j,k;
	TRI *tri;
	POINT *p;
	int dim = 3;
	double ks = geom_set->ks;
	double m_s = geom_set->m_s;
	double lambda_s = geom_set->lambda_s;

	unsort_surf_point(surf);
	for (tri = first_tri(surf); !at_end_of_tri_list(tri,surf); 
			tri = tri->next)
	{
	    for (j = 0; j < 3; ++j)
	    {
		p = Point_of_tri(tri)[j];
		if (sorted(p) || Boundary_point(p)) continue;
		for (k = 0; k < dim; ++k)
		{
		    x[*n][k] = Coords(p)[k];
		    v[*n][k] = p->vel[k];
		}
	    	spring_force_at_point2(f[*n],p,tri,surf,ks);
		for (k = 0; k < dim; ++k)
		{
		    f[*n][k] -= lambda_s*(v[*n][k]);
		    f[*n][k] /= m_s;
		}
		sorted(p) = YES;
	    	++(*n);
	    }
	}
}	/* end compute_surf_accel2 */

static void spring_force_at_point2(
	double *f,
	POINT *p,
	TRI *tri,
	SURFACE *surf,
	double ks)
{
	TRI *tris[MAX_NUM_RING1];
	int i,j,k,nt;
	POINT *p_nb;
	double length0,length,dir[3],vect[3];
	
	PointAndFirstRingTris(p,Hyper_surf_element(tri),Hyper_surf(surf),
				&nt,tris);
	for (k = 0; k < 3; ++k) f[k] = 0.0;
	for (i = 0; i < nt; ++i)
	{
	    for (j = 0; j < 3; ++j)
	    {
		if (Point_of_tri(tris[i])[j] == p)
		{
		    length0 = tris[i]->side_length0[j];
		    p_nb = Point_of_tri(tris[i])[(j+1)%3];
		    length = separation(p,p_nb,3);
	    	    for (k = 0; k < 3; ++k)
		    {
			dir[k] = (Coords(p_nb)[k] - Coords(p)[k])/length;
			vect[k] = (Coords(p_nb)[k] - Coords(p)[k]) -
				length0*tris[i]->side_dir0[j][k];
			f[k] += ks*vect[k];
		    }
		    if (is_side_bdry(tris[i],(j+2)%3))
		    {
			(void) printf("Detect boundary "
				"in spring_force_at_point1()\n");
			clean_up(ERROR);
		    }
		}
	    }
	}
}	/* end spring_force_at_point2 */

extern void compute_node_accel3(
	PARACHUTE_SET *geom_set,
	NODE *node,
	double **f,
	double **x,
	double **v,
	int *n)
{
	CURVE **c;
	BOND *b;
	double x_diff,len0,len,dir[MAXD],vect[MAXD];
	POINT *p,*p_nb;
	int i,j,dim = Dimension(node->interface);
	double ks = geom_set->ks;
	double kl = geom_set->kl;
	double m_s = geom_set->m_s;
	double m_c = geom_set->m_c;
	double *V_refs = geom_set->V_refs;
	double lambda_s = geom_set->lambda_s;
	double lambda_c = geom_set->lambda_c;
	double m_l;

	if (dim == 3)
	{
	    AF_NODE_EXTRA *extra = (AF_NODE_EXTRA*)node->extra;
	    if (extra != NULL && extra->af_node_type == LOAD_NODE)
	    {
	    	Front *front = geom_set->front;
	    	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;

	    	m_l = af_params->payload;
	    }
	}

	for (i = 0; i < dim; ++i)
	{
	    x[*n][i] = Coords(node->posn)[i];
	    v[*n][i] = node->posn->vel[i];
	    f[*n][i] = 0.0;
	}
	for (c = node->out_curves; c && *c; ++c)
	{
	    b = (*c)->first;
	    len0 = bond_length0(b);
	    for (j = 0; j < dim; ++j)
	    {
		vect[j] = len0*((Coords(b->end)[j] - Coords(b->start)[j])
				/bond_length(b) - b->dir0[j]);
		if (dim == 3)
		{
		    if (is_load_node(node) == YES)
		    	f[*n][j]   += kl*vect[j]/m_l;
		    else if (hsbdry_type(*c) == STRING_HSBDRY)
	    	    	f[*n][j]   += kl*vect[j]/m_s;
		    else if (hsbdry_type(*c) == MONO_COMP_HSBDRY)
	    	    	f[*n][j]   += ks*vect[j]/m_s;
		}
		else
		    f[*n][j]   += kl*vect[j]/m_c;
	    }
	}
	for (c = node->in_curves; c && *c; ++c)
	{
	    b = (*c)->last;
	    len0 = bond_length0(b);
	    for (j = 0; j < dim; ++j)
	    {
		vect[j] = len0*((Coords(b->end)[j] - Coords(b->start)[j])
				/bond_length(b) - b->dir0[j]);
		if (dim == 3)
		{
		    if (is_load_node(node) == YES)
		    	f[*n][j]   -= kl*vect[j]/m_l;
		    else if (hsbdry_type(*c) == STRING_HSBDRY)
	    	    	f[*n][j]   -= kl*vect[j]/m_s;
		    else if (hsbdry_type(*c) == MONO_COMP_HSBDRY)
	    	    	f[*n][j]   -= ks*vect[j]/m_s;
		}
		else
		    f[*n][j]   -= kl*vect[j]/m_c;
	    }
	}
	if (dim == 3)
	{
	    BOND_TRI **btris;
	    TRI **tris;
	    int k,side,nt,ns;
	    SURFACE *out_surfs[10];

	    ns = 0;
	    for (c = node->out_curves; c && *c; ++c)
	    {
		b = (*c)->first;
		p = b->start;
		for (btris = Btris(b); btris && *btris; ++btris)
		{
		    out_surfs[ns++] = (*btris)->surface;
		    nt = FT_FirstRingTrisAroundPoint(p,(*btris)->tri,&tris);
		    for (j = 0; j < nt; ++j)
		    {
			for (side = 0; side < 3; ++side)
			{
			    if (p == Point_of_tri(tris[j])[side])
			    {
				if (is_side_bdry(tris[j],side))
				    continue;
				p_nb = Point_of_tri(tris[j])[(side+1)%3];
				len0 = tris[j]->side_length0[side];
				len = separation(p,p_nb,3);
	    			x_diff = len - len0; 
				for (k = 0; k < 3; ++k)
                        	{
                            	    dir[k] = (Coords(p_nb)[k] - 
						Coords(p)[k])/len;
                            	    vect[k] = len0*(dir[k]
					- tris[j]->side_dir0[side][k]);
                            	    f[*n][k] += ks*vect[k]/m_s;
                        	}
			    }
			}
		    }
		}
	    }
	    for (c = node->in_curves; c && *c; ++c)
	    {
		if (is_closed_curve(*c)) continue;
		b = (*c)->last;
		p = b->end;
		for (btris = Btris(b); btris && *btris; ++btris)
		{
		    boolean duplicate_surf = NO;
		    for (j = 0; j < ns; ++j)
			if ((*btris)->surface == out_surfs[j])
			    duplicate_surf = YES;
		    if (duplicate_surf == YES) continue;
		    nt = FT_FirstRingTrisAroundPoint(p,(*btris)->tri,&tris);
		    for (j = 0; j < nt; ++j)
		    {
			for (side = 0; side < 3; ++side)
			{
			    if (p == Point_of_tri(tris[j])[side])
			    {
				if (is_side_bdry(tris[j],side))
				    continue;
				p_nb = Point_of_tri(tris[j])[(side+1)%3];
				len0 = tris[j]->side_length0[side];
				len = separation(p,p_nb,3);
				x_diff = len - len0;
				for (k = 0; k < 3; ++k)
                        	{
                            	    dir[k] = (Coords(p_nb)[k] - 
						Coords(p)[k])/len;
                            	    vect[k] = len0*(dir[k]
					- tris[j]->side_dir0[side][k]);
                            	    f[*n][k] += ks*vect[k]/m_s;
                        	}
			    }
			}
		    }
		}
	    }
	    if (!is_load_node(node))
	    {
	    	for (i = 0; i < 3; ++i)
	    	    f[*n][i] -= lambda_s*(v[*n][i] - V_refs[i])/m_s;
	    }
	}
	else
	{
	    for (i = 0; i < 3; ++i)
	    	f[*n][i] -= lambda_c*(v[*n][i] - V_refs[i])/m_c;
	}
	(*n)++;
}	/* end compute_node_accel3 */

extern void compute_curve_accel3(
	PARACHUTE_SET *geom_set,
	CURVE *curve,
	double **f,
	double **x,
	double **v,
	int *n)
{
	int i,j;
	double x_diff;
	BOND *b;
	double dir[MAXD],len0,vect[MAXD];
	int dim = Dimension(curve->interface);
	double kc,m_c,lambda_c;
	double *V_refs = geom_set->V_refs;

	if (dim == 3)
	{
	    if (hsbdry_type(curve) == STRING_HSBDRY)
	    {
	    	kc = geom_set->kl;
	    	m_c = geom_set->m_c;
	    	lambda_c = geom_set->lambda_c;
	    }
	    else
	    {
	    	kc = geom_set->ks;
	    	m_c = geom_set->m_s;
	    	lambda_c = geom_set->lambda_s;
	    }
	}
	else
	{
	    kc = geom_set->kl;
	    m_c = geom_set->m_c;
	    lambda_c = geom_set->lambda_c;
	}
	i = *n;
	for (b = curve->first; b != curve->last; b = b->next)
	{
	    for (j = 0; j < dim; ++j)
	    {
	    	x[i][j] = Coords(b->end)[j];
	    	v[i][j] = b->end->vel[j];
		f[i][j] = -lambda_c*(v[i][j] - V_refs[j])/m_c;
	    }
	    i++;
	}

	i = *n;
	for (b = curve->first; b != NULL; b = b->next)
	{
	    len0 = bond_length0(b);
	    for (j = 0; j < dim; ++j)
	    {
		vect[j] = len0*((Coords(b->end)[j] - Coords(b->start)[j])
				/bond_length(b) - b->dir0[j]);
		if (b != curve->first)
		{
	    	    f[i-1][j]   += kc*vect[j]/m_c;
		}
		if (b != curve->last)
		{
	    	    f[i][j] -= kc*vect[j]/m_c;
		}
	    }
	    if (b != curve->last) i++;
	}
	printf("f[20] = %f %f\n",f[20][0],f[20][1]);

	if (dim == 3)
	{
	    POINT *p,*p_nb;
	    BOND_TRI **btris;
	    TRI **tris;
	    int j,k,side,nt;
	    double length0,length;
	    i = *n;
	    for (b = curve->first; b != curve->last; b = b->next)
	    {
		p = b->end;
		for (btris = Btris(b); btris && *btris; ++btris)
		{
		    nt = FT_FirstRingTrisAroundPoint(p,(*btris)->tri,&tris);
		    for (j = 0; j < nt; ++j)
		    {
			for (side = 0; side < 3; ++side)
			{
			    if (p == Point_of_tri(tris[j])[side])
			    {
				if (is_side_bdry(tris[j],side))
				    continue;
				p_nb = Point_of_tri(tris[j])[(side+1)%3];
				length0 = tris[j]->side_length0[side];
				length = separation(p,p_nb,3);
				for (k = 0; k < 3; ++k)
                        	{
                            	    dir[k] = (Coords(p_nb)[k] - 
						Coords(p)[k])/length;
				    vect[k] = length0*(dir[k]
					- tris[j]->side_dir0[side][k]);
                            	    f[i][k] += kc*vect[k]/m_c;
                        	}
			    }
			}
		    }
		}
		i++;
	    }
	}
	*n = i;
}	/* end compute_curve_accel3 */


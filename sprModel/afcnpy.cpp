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
static void compute_total_canopy_force2d(Front*,double*,double*);
static void compute_total_canopy_force3d(Front*,double*,double*);
static void compute_center_of_mass_velo(PARACHUTE_SET*);
static void set_canopy_velocity(PARACHUTE_SET*,double**);
static boolean curve_in_pointer_list(CURVE*,CURVE**);

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
	
	if (is_registered_point(surf,p))
	{
	    for (i = 0; i < 3; ++i)
		f[i] = 0.0;
	    return;
	}
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
	if (wave_type(*hs) == ELASTIC_BOUNDARY) //TMP
	    return CONST_V_PDE_BOUNDARY;
	if (wave_type(*hs) == DIRICHLET_BOUNDARY)
	{
	    if (boundary_state(*hs))
	    	return CONST_V_PDE_BOUNDARY;
	    else
	    	return CONST_P_PDE_BOUNDARY;
	}
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
        PARACHUTE_SET *new_geom_set)
{
	static int size = 0;
	static double **x_old,**x_new,**v_old,**v_new;
	static double **x_pos,**v_pos;
	static double **accel;
	double lambda_s,m_s,lambda_l,m_l;
	AF_PARAMS *af_params = (AF_PARAMS*)fr->extra2;
	int i,j,num_pts;
	int n,n_tan = af_params->n_tan;
	double fr_dt = fr->dt;
	double dt = fr_dt/(double)n_tan;
	int num_strings = new_geom_set->num_strings;
	int n_cps,n_sps;
	double dt_tol;
	double xcom[MAXD],vcom[MAXD];
	double coeff1,coeff2;
	static SPRING_VERTEX *sv;
	static boolean first = YES;

	start_clock("set_data");

	n_cps = FT_NumOfSurfPoints(new_geom_set->canopy);  /* canopy pts */
	n_sps = 1;				/* load node */
	for (i = 0; i < num_strings; ++i)	/* string interior pts */
	{
	    n_sps += FT_NumOfCurvePoints(new_geom_set->string_curves[i]) - 2;
	}

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

	new_geom_set->ks = af_params->ks;
	new_geom_set->lambda_s = af_params->lambda_s;
	new_geom_set->m_s = af_params->m_s;

	new_geom_set->kl = af_params->kl;
	new_geom_set->lambda_l = af_params->lambda_l;
	new_geom_set->m_l = af_params->m_l;

	new_geom_set->kg = af_params->kg;
	new_geom_set->lambda_g = af_params->lambda_g;
	new_geom_set->m_g = af_params->m_g;

	new_geom_set->n_cps = n_cps;
	new_geom_set->n_sps = n_sps;
	dt_tol = sqrt((af_params->m_s)/(af_params->ks))/10.0;
	if (af_params->m_l != 0.0 &&
	    dt_tol > sqrt((af_params->m_l)/(af_params->kl))/10.0)
	    dt_tol = sqrt((af_params->m_l)/(af_params->kl))/10.0;
	if (af_params->m_g != 0.0 &&
	    dt_tol > sqrt((af_params->m_g)/(af_params->kg))/10.0)
	    dt_tol = sqrt((af_params->m_g)/(af_params->kg))/10.0;
	if (debugging("step_size"))
	{
	    (void) printf("Input surface parameters:\n");
	    (void) printf("ks = %f  m_s = %f  lambda_s = %f\n",
			new_geom_set->ks,
			new_geom_set->m_s,
			new_geom_set->lambda_s);
	    (void) printf("Input string parameters:\n");
	    (void) printf("kl = %f  m_l = %f  lambda_l = %f\n",
			new_geom_set->kl,
			new_geom_set->m_l,
			new_geom_set->lambda_l);
	    (void) printf("Input gore parameters:\n");
	    (void) printf("kg = %f  m_g = %f  lambda_g = %f\n",
			new_geom_set->kg,
			new_geom_set->m_g,
			new_geom_set->lambda_g);
	    (void) printf("n_cps = %d  n_sps = %d  num_pts = %d\n",
			n_cps,n_sps,num_pts);
	    (void) printf("number of strings = %d\n",
			new_geom_set->num_strings);
	    (void) printf("number of mono bdry = %d\n",
			new_geom_set->num_mono_hsbdry);
	    (void) printf("number of gore bdry = %d\n",
			new_geom_set->num_gore_hsbdry);
	    (void) printf("number of gore node = %d\n",
			new_geom_set->num_gore_nodes);
	}

	if (dt > dt_tol)
	{
	    n_tan = (int)(fr_dt/dt_tol);
	    dt = fr_dt/(double)n_tan;
	}
	(void) printf("\nfr_dt = %f  dt_tol = %20.14f  dt = %20.14f\n",
				fr_dt,dt_tol,dt);
	(void) printf("Number of interior sub-steps = %d\n\n",n_tan);

	new_geom_set->dt = dt;

	if (size < num_pts)
	{
	    size = num_pts;
	    if (v_old != NULL)
	    {
		FT_FreeThese(8,v_old,v_new,x_old,x_new,sv,x_pos,v_pos,accel);
	    }
	    FT_VectorMemoryAlloc((POINTER*)&x_pos,size,sizeof(double*));
	    FT_VectorMemoryAlloc((POINTER*)&v_pos,size,sizeof(double*));
	    FT_VectorMemoryAlloc((POINTER*)&sv,size,sizeof(SPRING_VERTEX));

	    FT_MatrixMemoryAlloc((POINTER*)&x_old,size,3,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&v_old,size,3,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&x_new,size,3,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&v_new,size,3,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&accel,size,3,sizeof(double));
	}

	count_canopy_spring_neighbors(new_geom_set,sv);
	count_string_spring_neighbors(new_geom_set,sv);
	if (first)
	{
	    set_spring_vertex_memory(sv,size);
	    first = NO;
	}
	set_canopy_spring_vertex(new_geom_set,x_pos,v_pos,sv);
	set_string_spring_vertex(new_geom_set,x_pos,v_pos,sv);
	stop_clock("set_data");

	start_clock("compute_data");
	for (i = 0; i < size; ++i)
	{
	    compute_spring_accel1(sv[i],accel[i],3);
	}
	for (i = 0; i < size; ++i)
        for (j = 0; j < 3; ++j)
	{
	    x_old[i][j] = x_pos[i][j];
	    v_old[i][j] = v_pos[i][j];
	}

	for (n = 0; n < n_tan; ++n)
	{
	    for (i = 0; i < size; ++i)
            for (j = 0; j < 3; ++j)
            {
                x_new[i][j] = x_old[i][j] + dt*v_old[i][j]/6.0;
                v_new[i][j] = v_old[i][j] + dt*accel[i][j]/6.0;
                x_pos[i][j] = x_old[i][j] + 0.5*v_old[i][j]*dt;
                v_pos[i][j] = v_old[i][j] + 0.5*accel[i][j]*dt;
            }
	    for (i = 0; i < size; ++i)
	    {
	    	compute_spring_accel1(sv[i],accel[i],3);
	    }

	    for (i = 0; i < size; ++i)
            for (j = 0; j < 3; ++j)
            {
                x_new[i][j] += dt*v_pos[i][j]/3.0;
                v_new[i][j] += dt*accel[i][j]/3.0;
                x_pos[i][j] = x_old[i][j] + 0.5*v_pos[i][j]*dt;
                v_pos[i][j] = v_old[i][j] + 0.5*accel[i][j]*dt;
            }
	    for (i = 0; i < size; ++i)
	    {
	    	compute_spring_accel1(sv[i],accel[i],3);
	    }

	    for (i = 0; i < size; ++i)
            for (j = 0; j < 3; ++j)
            {
                x_new[i][j] += dt*v_pos[i][j]/3.0;
                v_new[i][j] += dt*accel[i][j]/3.0;
                x_pos[i][j] = x_old[i][j] + v_pos[i][j]*dt;
                v_pos[i][j] = v_old[i][j] + accel[i][j]*dt;
            }
	    for (i = 0; i < size; ++i)
	    {
	    	compute_spring_accel1(sv[i],accel[i],3);
	    }

	    for (i = 0; i < size; ++i)
            for (j = 0; j < 3; ++j)
            {
                x_new[i][j] += dt*v_pos[i][j]/6.0;
                v_new[i][j] += dt*accel[i][j]/6.0;
            }
	    propagate_canopy(new_geom_set,x_new);
	    propagate_string(new_geom_set,x_new);
	    for (i = 0; i < size; ++i)
            for (j = 0; j < 3; ++j)
            {
		x_pos[i][j] = x_new[i][j];
		v_pos[i][j] = v_new[i][j];
	    }
	    if (n != n_tan-1)
	    {
		for (i = 0; i < size; ++i)
        	for (j = 0; j < 3; ++j)
		{
	    	    x_old[i][j] = x_pos[i][j];
	    	    v_old[i][j] = v_pos[i][j];
		}
	    	for (i = 0; i < size; ++i)
		{
	    	    compute_spring_accel1(sv[i],accel[i],3);
		}
	    }
	}
	stop_clock("compute_data");

	compute_center_of_mass_velo(new_geom_set);
	set_canopy_velocity(new_geom_set,v_new);

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
	int ns,nbc,ngc,ng;
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
	if (debugging("canopy"))
	    (void) printf("Entering compute_canopy_accel()\n");

	ns = geom_set->num_strings;
	nbc = geom_set->num_mono_hsbdry;
	ngc = geom_set->num_gore_hsbdry;
	ng = geom_set->num_gore_nodes;

	if (debugging("string_chord") || debugging("rigid_canopy") || 
	    debugging("ave_lift"))
	    n_start = n;
	compute_surf_accel(geom_set,geom_set->canopy,f,x,v,&n);
	for (i = 0; i < ng; ++i)
            compute_node_accel(geom_set,geom_set->gore_nodes[i],f,x,v,&n);
	for (i = 0; i < ns; ++i)
	    compute_node_accel(geom_set,geom_set->string_node[i],f,x,v,&n);
	for (i = 0; i < ngc; ++i)
	    compute_curve_accel(geom_set,geom_set->gore_hsbdry[i],f,x,v,&n);
	for (i = 0; i < nbc; ++i)
	{
	    compute_curve_accel(geom_set,geom_set->mono_hsbdry[i],f,x,
				v,&n);
	    if (is_closed_curve(geom_set->mono_hsbdry[i]))
		compute_node_accel(geom_set,
				geom_set->mono_hsbdry[i]->start,f,x,v,&n);	
	}
	if (debugging("string_chord") || debugging("rigid_canopy") 
			|| debugging("ave_lift"))
	{
	    n_end = n;
	    for (i = n_start; i < n_end; ++i)
	    {
		f[i][0] = f[i][1] = f[i][2] = 0.0;
		if (debugging("string_chord"))
		    v[i][0] = v[i][1] = v[i][2] = 0.0;
	    }
	}
	if (debugging("canopy"))
	    (void) printf("Leaving compute_canopy_accel()\n");
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

	if (debugging("canopy"))
	    (void) printf("Entering compute_string_accel()\n");

	n = geom_set->n_cps;
	compute_node_accel(geom_set,geom_set->load_node,f,x,v,&n);
	for (i = 0; i < ns; ++i)
	{
	    compute_curve_accel(geom_set,geom_set->string_curves[i],f,x,v,&n);
	}

	if (debugging("canopy"))
	    (void) printf("Leaving compute_string_accel()\n");
}	/* end  compute_string_accel */

static void assign_canopy_field(
	PARACHUTE_SET *geom_set,
        double **x,
        double **v)
{
	int n = 0;
	int i,ns,nbc,ngc,ng;

	if (debugging("canopy"))
	    (void) printf("Entering assign_canopy_field()\n");
	ng = geom_set->num_gore_nodes;
	ns = geom_set->num_strings;
	nbc = geom_set->num_mono_hsbdry;
	ngc = geom_set->num_gore_hsbdry;

	assign_surf_field(geom_set->canopy,x,v,&n);
	for (i = 0; i < ng; ++i)
            assign_node_field(geom_set->gore_nodes[i],x,v,&n);
	for (i = 0; i < ns; ++i)
	    assign_node_field(geom_set->string_node[i],x,v,&n);
	for (i = 0; i < ngc; ++i)
	    assign_curve_field(geom_set->gore_hsbdry[i],x,v,&n);
	for (i = 0; i < nbc; ++i)
	{
	    assign_curve_field(geom_set->mono_hsbdry[i],x,v,&n);
	    if (is_closed_curve(geom_set->mono_hsbdry[i]))
		 assign_node_field(geom_set->mono_hsbdry[i]->start,
			x,v,&n);
	}
	if (debugging("canopy"))
	    (void) printf("Leaving assign_canopy_field()\n");
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

	if (debugging("canopy"))
	    (void) printf("Entering assign_string_field()\n");

	assign_node_field(geom_set->load_node,x,v,&n);
	for (i = 0; i < ns; ++i)
	    assign_curve_field(geom_set->string_curves[i],x,v,&n);

	if (debugging("canopy"))
	    (void) printf("Leaving assign_string_field()\n");
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
	double dir[MAXD],len,len0,vect[MAXD];
	int dim = Dimension(curve->interface);
	double kl,m_l,lambda_l;

	if (dim == 3)
	{
	    if (hsbdry_type(curve) == STRING_HSBDRY)
	    {
	    	kl = geom_set->kl;
	    	m_l = geom_set->m_l;
	    	lambda_l = geom_set->lambda_l;
	    }
	    else if (hsbdry_type(curve) == GORE_HSBDRY)
	    {
	    	kl = geom_set->kg;
	    	m_l = geom_set->m_g;
	    	lambda_l = geom_set->lambda_g;
	    }
	    else
	    {
	    	kl = geom_set->ks;
	    	m_l = geom_set->m_s;
	    	lambda_l = geom_set->lambda_s;
	    }
	}
	else
	{
	    kl = geom_set->kl;
	    m_l = geom_set->m_l;
	    lambda_l = geom_set->lambda_l;
	}
	i = *n;
	for (b = curve->first; b != curve->last; b = b->next)
	{
	    x[i] = Coords(b->end);
	    v[i] = b->end->vel;
	    for (j = 0; j < dim; ++j)
	    {
		f[i][j] = -lambda_l*v[i][j]/m_l;
	    }
	    i++;
	}

	i = *n;
	for (b = curve->first; b != NULL; b = b->next)
	{
	    len = separation(b->start,b->end,dim);
	    len0 = bond_length0(b);
	    x_diff = len - len0; 
	    for (j = 0; j < dim; ++j)
	    {
		dir[j] = (Coords(b->end)[j] - Coords(b->start)[j])
				/len;
		vect[j] = x_diff*dir[j];
		if (b != curve->first)
		{
	    	    f[i-1][j]   += kl*vect[j]/m_l;
		}
		if (b != curve->last)
		{
	    	    f[i][j] -= kl*vect[j]/m_l;
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
	    double ks = geom_set->ks;
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
                            	    f[i][k] += ks*(length - length0)*
					dir[k]/m_l;
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
	INTERFACE *intfc = geom_set->front->interf;
	int i,j,dim = Dimension(intfc);
	double ks = geom_set->ks;
	double kl = geom_set->kl;
	double kg = geom_set->kg;
	double mass;
	double lambda_s = geom_set->lambda_s;
	double lambda_l = geom_set->lambda_l;
	double lambda_g = geom_set->lambda_g;

	if (dim == 3)
	{
	    AF_NODE_EXTRA *extra = (AF_NODE_EXTRA*)node->extra;
	    if (extra != NULL)
	    {
		if (extra->af_node_type == LOAD_NODE)
		{
	    	    Front *front = geom_set->front;
	    	    AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
		    mass = af_params->payload;
		}
		else if (extra->af_node_type == GORE_NODE)
                    mass = geom_set->m_g;
		else if (extra->af_node_type == STRING_NODE)
                    mass = geom_set->m_s;
	    }
	    else
                mass = geom_set->m_s;
	}
	else
            mass = geom_set->m_l;

	x[*n] = Coords(node->posn);
	v[*n] = node->posn->vel;
	for (i = 0; i < dim; ++i)
	{
	    f[*n][i] = 0.0;
	}
	for (c = node->out_curves; c && *c; ++c)
	{
	    b = (*c)->first;
	    len = separation(b->start,b->end,dim);
	    len0 = bond_length0(b);
	    x_diff = len - len0; 
	    for (j = 0; j < dim; ++j)
	    {
		dir[j] = (Coords(b->end)[j] - Coords(b->start)[j])
				/len;
		vect[j] = x_diff*dir[j];
		if (dim == 3)
		{
		    if (is_load_node(node) == YES)
		    	f[*n][j]   += kl*vect[j]/mass;
		    else if (hsbdry_type(*c) == STRING_HSBDRY)
	    	    	f[*n][j]   += kl*vect[j]/mass;
		    else if (hsbdry_type(*c) == MONO_COMP_HSBDRY)
	    	    	f[*n][j]   += ks*vect[j]/mass;
		    else if (hsbdry_type(*c) == GORE_HSBDRY)
	    	    	f[*n][j]   += kg*vect[j]/mass;
		}
		else
		    f[*n][j]   += kl*vect[j]/mass;
	    }
	}
	for (c = node->in_curves; c && *c; ++c)
	{
	    if (curve_in_pointer_list(*c,node->out_curves) && 
		!is_closed_curve(*c)) 
		continue;
	    b = (*c)->last;
	    len = separation(b->start,b->end,dim);
	    len0 = bond_length0(b);
	    x_diff = len - len0; 
	    for (j = 0; j < dim; ++j)
	    {
		dir[j] = (Coords(b->end)[j] - Coords(b->start)[j])
				/len;
		vect[j] = x_diff*dir[j];
		if (dim == 3)
		{
		    if (is_load_node(node) == YES)
		    	f[*n][j]   -= kl*vect[j]/mass;
		    else if (hsbdry_type(*c) == STRING_HSBDRY)
	    	    	f[*n][j]   -= kl*vect[j]/mass;
		    else if (hsbdry_type(*c) == MONO_COMP_HSBDRY)
	    	    	f[*n][j]   -= ks*vect[j]/mass;
		    else if (hsbdry_type(*c) == GORE_HSBDRY)
	    	    	f[*n][j]   -= kg*vect[j]/mass;
		}
		else
		    f[*n][j]   -= kl*vect[j]/mass;
	    }
	}
	if (dim == 3)
	{
	    BOND_TRI **btris;
	    TRI **tris,*tri_list[500];
	    int k,side,nt,num_tris;
	    TRI *tri;

	    num_tris = 0;
	    p = node->posn;
	    for (c = node->out_curves; c && *c; ++c)
	    {
		b = (*c)->first;
		for (btris = Btris(b); btris && *btris; ++btris)
		{
		    nt = FT_FirstRingTrisAroundPoint(p,(*btris)->tri,&tris);
		    for (j = 0; j < nt; ++j)
		    {
			if (!pointer_in_list((POINTER)tris[j],num_tris,
					(POINTER*)tri_list))
			    tri_list[num_tris++] = tris[j];
		    }
		}
	    }
	    for (c = node->in_curves; c && *c; ++c)
	    {
		b = (*c)->last;
		for (btris = Btris(b); btris && *btris; ++btris)
		{
		    nt = FT_FirstRingTrisAroundPoint(p,(*btris)->tri,&tris);
		    for (j = 0; j < nt; ++j)
		    {
			if (!pointer_in_list((POINTER)tris[j],num_tris,
					(POINTER*)tri_list))
			    tri_list[num_tris++] = tris[j];
		    }
		}
	    }
	    for (i = 0; i < num_tris; ++i)
	    {
		tri = tri_list[i];
		for (side = 0; side < 3; ++side)
		{
		    if (p == Point_of_tri(tri)[side])
		    {
			if (is_side_bdry(tri,side))
			    continue;
			p_nb = Point_of_tri(tri)[(side+1)%3];
			len0 = tri->side_length0[side];
			len = separation(p,p_nb,3);
    			x_diff = len - len0; 
			for (k = 0; k < 3; ++k)
                       	{
                       	    dir[k] = (Coords(p_nb)[k] - 
					Coords(p)[k])/len;
                       	    f[*n][k] += ks*x_diff*dir[k]/mass;
                       	}
		    }
		}
	    }
	    if (!is_load_node(node))
	    {
	    	for (i = 0; i < 3; ++i)
	    	    f[*n][i] -= lambda_s*v[*n][i]/mass;
	    }
	}
	else
	{
	    for (i = 0; i < 3; ++i)
	    	f[*n][i] -= lambda_l*v[*n][i]/mass;
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
	int ng,ngc,ns,n = 0;

	if (debugging("string_chord") || debugging("folding"))
	    return;

	if (debugging("canopy"))
	    (void) printf("Entering propagate_canopy()\n");

	hs = Hyper_surf(canopy);
	propagate_surface(geom_set,canopy,x,&n);
	ng = geom_set->num_gore_nodes;
	for (i = 0; i < ng; ++i)
	{
	    node = geom_set->gore_nodes[i];
	    propagate_node(geom_set,node,x,&n);
	}
	ns = geom_set->num_strings;
	for (i = 0; i < ns; ++i)
	{
	    node = geom_set->string_node[i];
	    propagate_node(geom_set,node,x,&n);
	}
	ngc = geom_set->num_gore_hsbdry;
	for (i = 0; i < ngc; ++i)
        {
	    curve = geom_set->gore_hsbdry[i];
	    propagate_curve(geom_set,curve,x,&n);
	}
	for (i = 0; i < ns; ++i)
        {
	    curve = geom_set->mono_hsbdry[i];
	    propagate_curve(geom_set,curve,x,&n);
	    if (is_closed_curve(curve))
	    {
	    	node = curve->start;
	    	propagate_node(geom_set,node,x,&n);
	    }
	}
	if (debugging("canopy"))
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

	if (debugging("canopy"))
	    (void) printf("Entering propagate_string()\n");
	n = geom_set->n_cps;
	ns = geom_set->num_strings;

	propagate_node(geom_set,load_node,x,&n);
	if (debugging("folding"))
	    return;

	for (i = 0; i < ns; ++i)
	{
	    c = geom_set->string_curves[i];
	    propagate_curve(geom_set,c,x,&n);
	}
	if (debugging("canopy"))
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

	if (debugging("canopy"))
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
	if (debugging("canopy"))
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
	int n,ng,ngc,ns,nbc;
	double crds_max[MAXD];
	int gindex_max;

	if (debugging("canopy"))
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
		    sl->vel[j] = sl->impulse[j] + nor_speed*nor[j];
		    sr->vel[j] = sr->impulse[j] + nor_speed*nor[j];
	    	    FT_RecordMaxFrontSpeed(j,sl->vel[j],NULL,Coords(p),front);
		}
	    	FT_RecordMaxFrontSpeed(3,Mag3d(sl->vel),NULL,Coords(p),front);
		if (max_speed < Mag3d(sl->vel)) 
		{
		    max_speed = Mag3d(sl->vel);
		    gindex_max = Gindex(p);
		    for (j = 0; j < 3; ++j)
			crds_max[j] = Coords(p)[j];
		}
		sorted(p) = YES;
		n++;
	    }
	}

	ng = geom_set->num_gore_nodes;
	for (i = 0; i < ng; ++i)
	{
	    node = geom_set->gore_nodes[i];
		if (the_point(node->posn))
		{
		    printf("n = %d\n",n);
		    printf("position 1\n");
		    clean_up(0);
		}
	    for (c = node->out_curves; c && *c; ++c)
            {
		if (hsbdry_type(*c) != GORE_HSBDRY) continue;
                b = (*c)->first;
                p = b->start;
		for (btris = Btris(b); btris && *btris; ++btris)
		{
                    p->hse = hse = Hyper_surf_element((*btris)->tri);
                    p->hs = hs = Hyper_surf((*btris)->surface);
                    FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
		    FT_NormalAtPoint(p,front,nor,NO_COMP);
		    vel = v[n];
		    nor_speed = scalar_product(vel,nor,dim);
		    if (max_speed < fabs(nor_speed)) 
		    {
			max_speed = fabs(nor_speed);
		    	gindex_max = Gindex(p);
		    	for (j = 0; j < 3; ++j)
			    crds_max[j] = Coords(p)[j];
		    }
                    for (j = 0; j < dim; ++j)
		    	sl->vel[j] = sr->vel[j] = sl->impulse[j] + nor_speed*nor[j];
		}
            }
            for (c = node->in_curves; c && *c; ++c)
            {
		if (hsbdry_type(*c) != GORE_HSBDRY) continue;
                b = (*c)->last;
                p = b->end;
                btris = Btris(b);
		for (btris = Btris(b); btris && *btris; ++btris)
		{
                    p->hse = hse = Hyper_surf_element((*btris)->tri);
                    p->hs = hs = Hyper_surf((*btris)->surface);
                    FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
		    FT_NormalAtPoint(p,front,nor,NO_COMP);
		    vel = v[n];
		    nor_speed = scalar_product(vel,nor,dim);
		    if (max_speed < fabs(nor_speed)) 
		    {
			max_speed = fabs(nor_speed);
		    	gindex_max = Gindex(p);
		    	for (j = 0; j < 3; ++j)
			    crds_max[j] = Coords(p)[j];
		    }
                    for (j = 0; j < dim; ++j)
		    	sl->vel[j] = sr->vel[j] = 
				sl->impulse[j] + nor_speed*nor[j];
		}
            }
	    n++;
	}

	ns = geom_set->num_strings;
	for (i = 0; i < ns; ++i)
	{
	    node = geom_set->string_node[i];
		if (the_point(node->posn))
		{
		    printf("n = %d\n",n);
		    printf("position 2\n");
		    clean_up(0);
		}
	    for (c = node->out_curves; c && *c; ++c)
            {
		if (hsbdry_type(*c) != MONO_COMP_HSBDRY &&
		    hsbdry_type(*c) != GORE_HSBDRY) 
		    continue;
                b = (*c)->first;
                p = b->start;
		for (btris = Btris(b); btris && *btris; ++btris)
		{
                    p->hse = hse = Hyper_surf_element((*btris)->tri);
                    p->hs = hs = Hyper_surf((*btris)->surface);
                    FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
		    FT_NormalAtPoint(p,front,nor,NO_COMP);
		    vel = v[n];
		    nor_speed = scalar_product(vel,nor,dim);
		    if (max_speed < fabs(nor_speed)) 
		    {
		    	max_speed = fabs(nor_speed);
		    	gindex_max = Gindex(p);
		    	for (j = 0; j < 3; ++j)
			    crds_max[j] = Coords(p)[j];
		    }
                    for (j = 0; j < dim; ++j)
		    	sl->vel[j] = sr->vel[j] = 
				sl->impulse[j] + nor_speed*nor[j];
		}
            }
            for (c = node->in_curves; c && *c; ++c)
            {
		if (hsbdry_type(*c) != MONO_COMP_HSBDRY &&
		    hsbdry_type(*c) != GORE_HSBDRY) 
		    continue;
                b = (*c)->last;
                p = b->end;
		for (btris = Btris(b); btris && *btris; ++btris)
		{
                    p->hse = hse = Hyper_surf_element((*btris)->tri);
                    p->hs = hs = Hyper_surf((*btris)->surface);
                    FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
		    FT_NormalAtPoint(p,front,nor,NO_COMP);
		    vel = v[n];
		    nor_speed = scalar_product(vel,nor,dim);
		    if (max_speed < fabs(nor_speed)) 
		    {
		    	max_speed = fabs(nor_speed);
		    	gindex_max = Gindex(p);
		    	for (j = 0; j < 3; ++j)
			    crds_max[j] = Coords(p)[j];
		    }
                    for (j = 0; j < dim; ++j)
		    	sl->vel[j] = sr->vel[j] = 
				sl->impulse[j] + nor_speed*nor[j];
		    }
            }
	    n++;
	}

	ngc = geom_set->num_gore_hsbdry;
	for (i = 0; i < ngc; ++i)
	{
	    curve = geom_set->gore_hsbdry[i];
	    for (b = curve->first; b != curve->last; b = b->next)
            {
            	p = b->end;
		if (the_point(p))
		{
		    printf("n = %d\n",n);
		    printf("position 3\n");
		    clean_up(0);
		}
		for (btris = Btris(b); btris && *btris; ++btris)
            	{
                    p->hse = hse = Hyper_surf_element((*btris)->tri);
                    p->hs = hs = Hyper_surf((*btris)->surface);
                    FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
		    FT_NormalAtPoint(p,front,nor,NO_COMP);
		    vel = v[n];
		    nor_speed = scalar_product(vel,nor,dim);
                    for (j = 0; j < dim; ++j)
		    	sl->vel[j] = sr->vel[j] = 
				sl->impulse[j] + nor_speed*nor[j];
            	}
            	n++;
            }
	}

	nbc = geom_set->num_mono_hsbdry;
	for (i = 0; i < nbc; ++i)
	{
	    curve = geom_set->mono_hsbdry[i];
	    for (b = curve->first; b != curve->last; b = b->next)
            {
            	p = b->end;
		if (the_point(p))
		{
		    printf("n = %d\n",n);
		    printf("position 4\n");
		    clean_up(0);
		}
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
		    	sl->vel[j] = sr->vel[j] = sl->impulse[j] + nor_speed*nor[j];
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
	    (void) printf("Coordinate of max speed: %f %f %f\n",
				crds_max[0],crds_max[1],crds_max[2]);
	    (void) printf("Global index of point-max: %d\n",gindex_max);
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
	if (debugging("canopy"))
	    (void) printf("Leaving set_canopy_velocity()\n");
}	/* end set_canopy_velocity */

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
	double kg = geom_set->kg;
	double mass;
	double lambda_s = geom_set->lambda_s;
	double lambda_l = geom_set->lambda_l;
	double lambda_g = geom_set->lambda_g;

	if (dim == 3)
	{
	    AF_NODE_EXTRA *extra = (AF_NODE_EXTRA*)node->extra;
	    if (extra != NULL)
	    {
		if (extra->af_node_type == LOAD_NODE)
		{
	    	    Front *front = geom_set->front;
	    	    AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	    	    mass = af_params->payload;
		}
		else if (extra->af_node_type == GORE_NODE)
		    mass = geom_set->m_g;
		else if (extra->af_node_type == STRING_NODE)
		    mass = geom_set->m_s;
	    }
	    else
		mass = geom_set->m_s;
	}
	else
	    mass = geom_set->m_l;

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
		    	f[*n][j]   += kl*vect[j]/mass;
		    else if (hsbdry_type(*c) == STRING_HSBDRY)
	    	    	f[*n][j]   += kl*vect[j]/mass;
		    else if (hsbdry_type(*c) == MONO_COMP_HSBDRY)
	    	    	f[*n][j]   += ks*vect[j]/mass;
		    else if (hsbdry_type(*c) == GORE_HSBDRY)
	    	    	f[*n][j]   += kg*vect[j]/mass;
		}
		else
		    f[*n][j]   += kl*vect[j]/mass;
	    }
	}
	for (c = node->in_curves; c && *c; ++c)
	{
	    if (curve_in_pointer_list(*c,node->out_curves) && 
		!is_closed_curve(*c)) 
		continue;
	    b = (*c)->last;
	    len0 = bond_length0(b);
	    for (j = 0; j < dim; ++j)
	    {
		vect[j] = Coords(b->end)[j] - Coords(b->start)[j] -
				len0*b->dir0[j];
		if (dim == 3)
		{
		    if (is_load_node(node) == YES)
		    	f[*n][j]   -= kl*vect[j]/mass;
		    else if (hsbdry_type(*c) == STRING_HSBDRY)
	    	    	f[*n][j]   -= kl*vect[j]/mass;
		    else if (hsbdry_type(*c) == MONO_COMP_HSBDRY)
	    	    	f[*n][j]   -= ks*vect[j]/mass;
		    else if (hsbdry_type(*c) == GORE_HSBDRY)
	    	    	f[*n][j]   -= kg*vect[j]/mass;
		}
		else
		    f[*n][j]   -= kl*vect[j]/mass;
	    }
	}
	if (dim == 3)
	{
	    BOND_TRI **btris;
	    TRI **tris,*tri_list[500];
	    int k,side,nt,num_tris;
	    TRI *tri;

	    num_tris = 0;
	    p = node->posn;
	    for (c = node->out_curves; c && *c; ++c)
	    {
		b = (*c)->first;
		for (btris = Btris(b); btris && *btris; ++btris)
		{
		    nt = FT_FirstRingTrisAroundPoint(p,(*btris)->tri,&tris);
		    for (j = 0; j < nt; ++j)
		    {
			if (!pointer_in_list((POINTER)tris[j],num_tris,
					(POINTER*)tri_list))
			    tri_list[num_tris++] = tris[j];
		    }
		}
	    }
	    for (c = node->in_curves; c && *c; ++c)
	    {
		b = (*c)->last;
		for (btris = Btris(b); btris && *btris; ++btris)
		{
		    nt = FT_FirstRingTrisAroundPoint(p,(*btris)->tri,&tris);
		    for (j = 0; j < nt; ++j)
		    {
			if (!pointer_in_list((POINTER)tris[j],num_tris,
					(POINTER*)tri_list))
			    tri_list[num_tris++] = tris[j];
		    }
		}
	    }
	    for (i = 0; i < num_tris; ++i)
	    {
		tri = tri_list[i];
		for (side = 0; side < 3; ++side)
		{
		    if (p == Point_of_tri(tri)[side])
		    {
			if (is_side_bdry(tri,side))
			    continue;
			p_nb = Point_of_tri(tri)[(side+1)%3];
			len0 = tri->side_length0[side];
			len = separation(p,p_nb,3);
			x_diff = len - len0;
			for (k = 0; k < 3; ++k)
                       	{
                       	    dir[k] = tri->side_dir0[side][k]; 
                       	    vect[k] = (Coords(p_nb)[k] - Coords(p)[k])
					- len0*dir[k];
                       	    f[*n][k] += ks*vect[k]/mass;
                       	}
		    }
		}
	    }
	    if (!is_load_node(node))
	    {
	    	for (i = 0; i < 3; ++i)
	    	    f[*n][i] -= lambda_s*v[*n][i]/mass;
	    }
	}
	else
	{
	    for (i = 0; i < 3; ++i)
	    	f[*n][i] -= lambda_l*v[*n][i]/mass;
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
	double kl,m_l,lambda_l,ks,lambda_s;

	if (dim == 3)
	{
	    if (hsbdry_type(curve) == STRING_HSBDRY)
	    {
	    	kl = geom_set->kl;
	    	m_l = geom_set->m_l;
	    	lambda_l = geom_set->lambda_l;
	    }
	    else if (hsbdry_type(curve) == GORE_HSBDRY)
	    {
	    	kl = geom_set->kg;
	    	m_l = geom_set->m_g;
	    	lambda_l = geom_set->lambda_g;
	    }
	    else
	    {
	    	kl = geom_set->ks;
	    	m_l = geom_set->m_s;
	    	lambda_l = geom_set->lambda_s;
	    }
	    ks = geom_set->ks;
	    lambda_s = geom_set->lambda_s;
	}
	else
	{
	    kl = geom_set->kl;
	    m_l = geom_set->m_l;
	    lambda_l = geom_set->lambda_l;
	}
	i = *n;
	for (b = curve->first; b != curve->last; b = b->next)
	{
	    for (j = 0; j < dim; ++j)
	    {
	    	x[i][j] = Coords(b->end)[j];
	    	v[i][j] = b->end->vel[j];
		f[i][j] = -lambda_l*v[i][j]/m_l;
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
	    	    f[i-1][j]   += kl*vect[j]/m_l;
		}
		if (b != curve->last)
		{
	    	    f[i][j] -= kl*vect[j]/m_l;
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
				for (k = 0; k < 3; ++k)
                        	{
                            	    dir[k] = tris[j]->side_dir0[side][k]; 
				    vect[k] = Coords(p_nb)[k] - Coords(p)[k]
					- length0*dir[k];
                            	    f[i][k] += ks*vect[k]/m_l;
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
			dir[k] = tris[i]->side_dir0[j][k];
			vect[k] = (Coords(p_nb)[k] - Coords(p)[k]) -
				//length0*tris[i]->side_dir0[j][k];
				length0*dir[k];
			f[k] += ks*vect[k];
		    }
		    if (is_side_bdry(tris[i],(j+2)%3))
		    {
			(void) printf("Detect boundary "
				"in spring_force_at_point2()\n");
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
	double m_l = geom_set->m_l;
	double lambda_s = geom_set->lambda_s;
	double lambda_l = geom_set->lambda_l;
	double payload;

	if (dim == 3)
	{
	    AF_NODE_EXTRA *extra = (AF_NODE_EXTRA*)node->extra;
	    if (extra != NULL && extra->af_node_type == LOAD_NODE)
	    {
	    	Front *front = geom_set->front;
	    	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;

	    	payload = af_params->payload;
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
		    	f[*n][j]   += kl*vect[j]/payload;
		    else if (hsbdry_type(*c) == STRING_HSBDRY)
	    	    	f[*n][j]   += kl*vect[j]/m_s;
		    else if (hsbdry_type(*c) == MONO_COMP_HSBDRY)
	    	    	f[*n][j]   += ks*vect[j]/m_s;
		}
		else
		    f[*n][j]   += kl*vect[j]/m_l;
	    }
	}
	for (c = node->in_curves; c && *c; ++c)
	{
	    //if (curve_in_pointer_list(*c,node->out_curves)) continue;
	    b = (*c)->last;
	    len0 = bond_length0(b);
	    for (j = 0; j < dim; ++j)
	    {
		vect[j] = len0*((Coords(b->end)[j] - Coords(b->start)[j])
				/bond_length(b) - b->dir0[j]);
		if (dim == 3)
		{
		    if (is_load_node(node) == YES)
		    	f[*n][j]   -= kl*vect[j]/payload;
		    else if (hsbdry_type(*c) == STRING_HSBDRY)
	    	    	f[*n][j]   -= kl*vect[j]/m_s;
		    else if (hsbdry_type(*c) == MONO_COMP_HSBDRY)
	    	    	f[*n][j]   -= ks*vect[j]/m_s;
		}
		else
		    f[*n][j]   -= kl*vect[j]/m_l;
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
	   	//if (curve_in_pointer_list(*c,node->out_curves)) continue;
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
	    	    f[*n][i] -= lambda_s*v[*n][i]/m_s;
	    }
	}
	else
	{
	    for (i = 0; i < 3; ++i)
	    	f[*n][i] -= lambda_l*v[*n][i]/m_l;
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
	double kl,m_l,lambda_l;

	if (dim == 3)
	{
	    if (hsbdry_type(curve) == STRING_HSBDRY)
	    {
	    	kl = geom_set->kl;
	    	m_l = geom_set->m_l;
	    	lambda_l = geom_set->lambda_l;
	    }
	    else
	    {
	    	kl = geom_set->ks;
	    	m_l = geom_set->m_s;
	    	lambda_l = geom_set->lambda_s;
	    }
	}
	else
	{
	    kl = geom_set->kl;
	    m_l = geom_set->m_l;
	    lambda_l = geom_set->lambda_l;
	}
	i = *n;
	for (b = curve->first; b != curve->last; b = b->next)
	{
	    for (j = 0; j < dim; ++j)
	    {
	    	x[i][j] = Coords(b->end)[j];
	    	v[i][j] = b->end->vel[j];
		f[i][j] = -lambda_l*v[i][j]/m_l;
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
	    	    f[i-1][j]   += kl*vect[j]/m_l;
		}
		if (b != curve->last)
		{
	    	    f[i][j] -= kl*vect[j]/m_l;
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
                            	    f[i][k] += kl*vect[k]/m_l;
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

static boolean curve_in_pointer_list(
	CURVE *c,
	CURVE **c_list)
{
	CURVE **pc;
	if (c_list == NULL) return NO;
	for (pc = c_list; pc && *pc; ++pc)
	{
	    if (c == *pc) return YES;
	}
	return NO;
}	/* end curve_in_pointer_list */

extern boolean is_registered_point(
	SURFACE *surf,
	POINT *p)
{
	REGISTERED_PTS *rgp = (REGISTERED_PTS*)surf->extra;
	int i,num_pts;
	int *global_ids;
	
	if (rgp == NULL) return NO;

	num_pts = rgp->num_pts;
	global_ids = rgp->global_ids;
	for (i = 0; i < num_pts; ++i)
	{
	    if (Gindex(p) == global_ids[i])
		return YES;
	}
	return NO;
}	/* end is_registered_point */

extern void propagate_surface(
        PARACHUTE_SET *geom_set,
        SURFACE *surf,
        double **x,
        int *n)
{
        int i,j;
        TRI *tri;
        POINT *p;
        STATE *sl,*sr;
        HYPER_SURF_ELEMENT *hse;
        HYPER_SURF         *hs;
	double dt = geom_set->dt;
	Front *front = geom_set->front;
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	double *g = iFparams->gravity;

	hs = Hyper_surf(surf);
	unsort_surf_point(surf);
        for (tri = first_tri(surf); !at_end_of_tri_list(tri,surf);
                        tri = tri->next)
        {
            hse = Hyper_surf_element(tri);
            for (i = 0; i < 3; ++i)
            {
                p = Point_of_tri(tri)[i];
                if (sorted(p) || Boundary_point(p)) continue;
                FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
		if (is_registered_point(surf,p))
		{
                    for (j = 0; j < 3; ++j)
                    {
                        x[*n][j] += sl->impulse[j]*dt;
                        sr->impulse[j] = sl->impulse[j] = sl->impulse[j];
                    }
		}
		else
                {
                    for (j = 0; j < 3; ++j)
                    {
                        x[*n][j] += (sl->impulse[j] + 0.5*g[j]*dt)*dt;
                        sr->impulse[j] = sl->impulse[j] = 
					sl->impulse[j] + g[j]*dt;
                    }
                }
                sorted(p) = YES;
                ++(*n);
            }
        }
}       /* propagate_surface */

extern void propagate_node(
        PARACHUTE_SET *geom_set,
	NODE *node,
        double **x,
        int *n)
{
        int i,j;
        POINT *p;
        STATE *sl,*sr;
	double dt = geom_set->dt;
	Front *front = geom_set->front;
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	double *g = iFparams->gravity;
	int dim = front->rect_grid->dim;

        sl = (STATE*)left_state(node->posn);
        sr = (STATE*)right_state(node->posn);
        for (j = 0; j < dim; ++j)
        {
            x[*n][j] += (sl->impulse[j] + 0.5*g[j]*dt)*dt;
            sr->impulse[j] = sl->impulse[j] = sl->impulse[j] + g[j]*dt;
        }
        ++(*n);
}	/* end propagate_node */

extern void propagate_curve(
        PARACHUTE_SET *geom_set,
	CURVE *curve,
        double **x,
        int *n)
{
        int i,j;
        POINT *p;
	BOND *b;
        STATE *sl,*sr;
	double dt = geom_set->dt;
	Front *front = geom_set->front;
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	double *g = iFparams->gravity;
	int dim = front->rect_grid->dim;

	for (b = curve->first; b != curve->last; b = b->next)
        {
            p = b->end;
            sl = (STATE*)left_state(p);
            sr = (STATE*)right_state(p);
            for (j = 0; j < dim; ++j)
            {
                x[*n][j] += (sl->impulse[j] + 0.5*g[j]*dt)*dt;
                sr->impulse[j] = sl->impulse[j] = sl->impulse[j] + g[j]*dt;
            }
            ++(*n);
        }
}	/* end propagate_curve */

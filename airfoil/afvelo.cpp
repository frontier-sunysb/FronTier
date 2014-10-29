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
	double cen[MAXD];
	double v0;
	double stop_time;
} VERTICAL_PARAMS;

typedef struct {
	double v0[MAXD];
	double stop_time;
} RANDOMV_PARAMS;

typedef struct {
	int num_pts;
	int *global_ids;
	double vel[MAXD];
} FIXAREA_PARAMS;

static void initVelocityFunc(FILE*,Front*);
static int zero_velo(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*,
                                double*);
static int toroidal_velo(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*,
                                double*);
static int parabolic_velo(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*,
                                double*);
static int singular_velo(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*,
                                double*);
static int vertical_velo(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*,
                                double*);
static int random_velo(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*,
                                double*);
static int marker_velo(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*,
                                double*);
static void init_fixarea_params(Front*,FILE*,FIXAREA_PARAMS*);
static void init_fixpoint_params(Front*,FILE*,FIXAREA_PARAMS*);
static void convert_to_point_mass(Front*, AF_PARAMS*);
static void checkSetGoreNodes(INTERFACE*);
static void set_gore_node(NODE*);


extern void setMotionParams(
	Front *front)
{
	FILE *infile = fopen(InName(front),"r");
	int i,dim = front->rect_grid->dim;
	char string[100];
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	INTERFACE *intfc = front->interf;
	boolean status;

	af_params->no_fluid = NO;
	if (dim == 3 && numOfGoreHsbdry(intfc) != 0)
	{
	    af_params->attach_gores = YES;
	    checkSetGoreNodes(intfc);
	}
        if (CursorAfterStringOpt(infile,
            "Entering yes to turn off fluid solver: "))
        {
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'y' || string[0] == 'Y')
                af_params->no_fluid = YES;
        }
	if (af_params->no_fluid == YES)
	{
	    front->curve_propagate = airfoil_curve_propagate;
	    initVelocityFunc(infile,front);
	}
	else
	{
	    front->_point_propagate = airfoil_point_propagate;
	    if (dim == 3)
	    	front->curve_propagate = airfoil_curve_propagate;
	}

	if (af_params->no_fluid == YES || 
	    af_params->is_parachute_system == NO)
	{
            CursorAfterString(infile,"Enter interior propagator:");
	    fscanf(infile,"%s",string);
	    (void) printf("%s\n",string);
	    if (dim == 2)
	    {
	    	switch (string[0])
	    	{
	    	case 'n':
	    	case 'N':
	    	    front->tan_curve_propagate = NULL;
	    	    break;
	    	case 'f':
	    	case 'F':
	    	    front->tan_curve_propagate 
				= fixed_length_tan_curve_propagate;
	    	    break;
	    	case 'e':
	    	case 'E':
	    	    if (string[1] == '2')
	    	    	front->tan_curve_propagate 
				= second_order_elastic_curve_propagate;
	    	    else
		    {
	    	    	front->tan_curve_propagate 
				= fourth_order_elastic_curve_propagate;
#if defined(__GPU__)
            		if (CursorAfterStringOpt(infile,
				"Enter yes to use GPU solver:"))
			{
	    		    fscanf(infile,"%s",string);
	    		    (void) printf("%s\n",string);
			    if (string[0] == 'y' || string[0] == 'Y')
				af_params->use_gpu = YES;
			}
#endif
		    }
	    	    break;
	    	default:
		    (void) printf("Unknown interior propagator!\n");
		    clean_up(ERROR);
	    	}
	    }
	    else if (dim == 3)
	    {
	    	switch (string[0])
	    	{
	    	case 'n':
	    	case 'N':
	    	    front->interior_propagate = NULL;
	    	    break;
	    	case 'e':
	    	case 'E':
	    	    if (string[1] == '2')
	    	    	front->interior_propagate 
				= second_order_elastic_surf_propagate;
	    	    else
		    {
	    	    	front->interior_propagate 
				= fourth_order_elastic_surf_propagate;
				//= fourth_order_elastic_set_propagate;
#if defined(__GPU__)
            		if (CursorAfterStringOpt(infile,
				"Enter yes to use GPU solver:"))
			{
	    		    fscanf(infile,"%s",string);
	    		    (void) printf("%s\n",string);
			    if (string[0] == 'y' || string[0] == 'Y')
				af_params->use_gpu = YES;
			}
#endif
		    }
	    	    break;
	    	case 'p':
	    	case 'P':
	    	    front->interior_propagate 
				= fourth_order_elastic_set_propagate;
	    	    break;
	    	default:
		    (void) printf("Unknown interior propagator!\n");
		    clean_up(ERROR);
	    	}
	    }
	}
	else
	    front->interior_propagate = fourth_order_elastic_set_propagate;

	CursorAfterString(infile,"Enter gravity:");
        for (i = 0; i < dim; ++i)
        {
            fscanf(infile,"%lf",&af_params->gravity[i]);
            (void) printf("%f ",af_params->gravity[i]);
        }
        (void) printf("\n");
	if (CursorAfterStringOpt(infile,"Enter payload:"))
	{
            fscanf(infile,"%lf",&af_params->payload);
            (void) printf("%f\n",af_params->payload);
	}

	if (af_params->no_fluid == NO)
	{
	    CursorAfterString(infile,
                        "Enter density and viscosity of the fluid:");
            fscanf(infile,"%lf %lf",&iFparams->rho2,&iFparams->mu2);
            (void) printf("%f %f\n",iFparams->rho2,iFparams->mu2);
	    if (FT_FrontContainWaveType(front,CONTACT))
	    {
            	CursorAfterString(infile,"Enter surface tension:");
            	fscanf(infile,"%lf",&iFparams->surf_tension);
            	(void) printf("%f\n",iFparams->surf_tension);
	    }
	    if (FT_FrontContainWaveType(front,ELASTIC_BOUNDARY))
	    {
            	CursorAfterString(infile,"Enter porosity of canopy:");
            	fscanf(infile,"%lf",&af_params->gamma);
            	(void) printf("%f\n",af_params->gamma);
            	CursorAfterString(infile,"Enter area density of canopy:");
            	fscanf(infile,"%lf",&af_params->area_dens);
            	(void) printf("%f\n",af_params->area_dens);
	    }
            CursorAfterString(infile,"Enter factor of smoothing radius:");
            fscanf(infile,"%lf",&iFparams->smoothing_radius);
            (void) printf("%f\n",iFparams->smoothing_radius);
            for (i = 0; i < dim; ++i)
	    	iFparams->gravity[i] = af_params->gravity[i];
	}
	status = FT_FrontContainWaveType(front,ELASTIC_BOUNDARY);
	status = pp_max_status(status);
	if (!status) return;

	af_params->n_sub = 1;
	CursorAfterString(infile,"Enter interior sub step number:");
	fscanf(infile,"%d",&af_params->n_sub);
	(void) printf("%d\n",af_params->n_sub);
        af_params->use_total_mass = NO;
        if (CursorAfterStringOpt(infile,"Enter yes to use total mass:"))
        {
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'y' || string[0] == 'Y')
                af_params->use_total_mass = YES;
        }

	if (dim == 3)
	{
	    CursorAfterString(infile,"Enter fabric spring constant:");
            fscanf(infile,"%lf",&af_params->ks);
            (void) printf("%f\n",af_params->ks);
            CursorAfterString(infile,"Enter fabric friction constant:");
            fscanf(infile,"%lf",&af_params->lambda_s);
            (void) printf("%f\n",af_params->lambda_s);
            if (af_params->use_total_mass)
            {
                CursorAfterString(infile,"Enter fabric total mass:");
                fscanf(infile,"%lf",&af_params->total_canopy_mass);
                (void) printf("%f\n",af_params->total_canopy_mass);
            }
            else
	    {
	    	CursorAfterString(infile,"Enter fabric point mass:");
            	fscanf(infile,"%lf",&af_params->m_s);
            	(void) printf("%f\n",af_params->m_s);
	    }
	    af_params->m_g = 0.0;
            if (af_params->attach_gores == YES)
	    {
		CursorAfterString(infile,"Enter gore spring constant:");
        	fscanf(infile,"%lf",&af_params->kg);
        	(void) printf("%f\n",af_params->kg);
        	CursorAfterString(infile,"Enter gore friction constant:");
        	fscanf(infile,"%lf",&af_params->lambda_g);
        	(void) printf("%f\n",af_params->lambda_g);
                if (af_params->use_total_mass)
                {
                    CursorAfterString(infile,"Enter gore total mass:");
                    fscanf(infile,"%lf",&af_params->total_gore_mass);
                    (void) printf("%f\n",af_params->total_gore_mass);
                }
                else
                {
                    CursorAfterString(infile,"Enter gore point mass:");
                    fscanf(infile,"%lf",&af_params->m_g);
                    (void) printf("%f\n",af_params->m_g);
                }

	    }
	}
	CursorAfterString(infile,"Enter string spring constant:");
        fscanf(infile,"%lf",&af_params->kl);
        (void) printf("%f\n",af_params->kl);
        CursorAfterString(infile,"Enter string friction constant:");
        fscanf(infile,"%lf",&af_params->lambda_l);
        (void) printf("%f\n",af_params->lambda_l);
        if (af_params->use_total_mass)
        {
            CursorAfterString(infile,"Enter string total mass:");
            fscanf(infile,"%lf",&af_params->total_string_mass);
            (void) printf("%f\n",af_params->total_string_mass);
        }
        else
        {
            CursorAfterString(infile,"Enter string point mass:");
            fscanf(infile,"%lf",&af_params->m_l);
            (void) printf("%f\n",af_params->m_l);
        }
	af_params->num_smooth_layers = 1;
	if (CursorAfterStringOpt(infile,"Enter number of smooth layers:"))
	{
            fscanf(infile,"%d",&af_params->num_smooth_layers);
            (void) printf("%d\n",af_params->num_smooth_layers);
	}
        if (af_params->use_total_mass)
            convert_to_point_mass(front,af_params);
	fclose(infile);
}	/* end setMotionParams */

static void convert_to_point_mass(
        Front *front,
        AF_PARAMS *af_params)
{
        INTERFACE *intfc;
        int num_str_pts, num_fabric_pts, num_gore_pts;
        SURFACE **s;
        CURVE **c;
        intfc = front->interf;
        int dim = Dimension(intfc);

        switch (dim)
        {
        case 2:
            num_str_pts = 0;
            for (c = intfc->curves; c && *c; ++c)
            {
                if (wave_type(*c) == ELASTIC_BOUNDARY)
                    num_str_pts +=  FT_NumOfCurvePoints(*c);
            }
            num_str_pts -= 2; //ignore the boundary points
            af_params->m_l = af_params->total_string_mass/num_str_pts;
            printf("string total mass = %f\n",af_params->total_string_mass);
            printf("string point number = %d\n",num_str_pts);
            printf("string point mass = %f\n",af_params->m_l);
            break;
        case 3:
            num_str_pts = num_fabric_pts  = 0;
            for (s = intfc->surfaces; s && *s; ++s)
            {
                if (wave_type(*s) == ELASTIC_BOUNDARY)
                    num_fabric_pts += FT_NumOfSurfPoints(*s);
            }
            for (c = intfc->curves; c && *c; ++c)
            {
                if (hsbdry_type(*c) == STRING_HSBDRY)
                     num_str_pts += FT_NumOfCurvePoints(*c);
            }
            af_params->m_s = af_params->total_canopy_mass/num_fabric_pts;
            if (num_str_pts != 0)
                af_params->m_l = af_params->total_string_mass/num_str_pts;
            else
                af_params->m_l = 0.01;
            printf("fabric total mass = %f\n",af_params->total_canopy_mass);
            printf("fabric point number = %d\n",num_fabric_pts);
            printf("fabric point mass = %f\n",af_params->m_s);
        }

}       /* end convert_to_point_mass */


static void initVelocityFunc(
	FILE *infile,
	Front *front)
{
	static VELO_FUNC_PACK velo_func_pack;
	static VORTEX_PARAMS *vortex_params; /* velocity function parameters */
        static BIPOLAR_PARAMS *dv_params;
	static VERTICAL_PARAMS *vert_params;
	static RANDOMV_PARAMS *randv_params;
	static TOROIDAL_PARAMS *toro_params;
	static PARABOLIC_PARAMS *para_params;
	static SINGULAR_PARAMS *sing_params;
	static FIXAREA_PARAMS *fixarea_params;
	int i,dim = front->rect_grid->dim;
	char string[100];
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;

	if (af_params->no_fluid == YES)
	{
	    front->curve_propagate = airfoil_curve_propagate;
	    velo_func_pack.point_propagate = airfoil_point_propagate;
	    (void) printf("Available velocity functions are:\n");
	    (void) printf("\tVortex velocity (R)\n");
	    (void) printf("\tDouble vortex velocity (D)\n");
	    (void) printf("\tVertical velocity (V)\n");
	    (void) printf("\tToroidal velocity (T)\n");
	    (void) printf("\tParabolic velocity (P)\n");
	    (void) printf("\tSingular velocity (S)\n");
	    (void) printf("\tZero velocity (Z)\n");
	    (void) printf("\tFixed area velocity (FA)\n");
	    (void) printf("\tFixed point velocity (FP)\n");
	    (void) printf("\tFree fall velocity (FF)\n");
            CursorAfterString(infile,"Enter velocity function: ");
	    fscanf(infile,"%s",string);
	    (void) printf("%s\n",string);
	    switch (string[0])
            {
            case 'r':
            case 'R':
		if (string[1] == 'o' || string[1] == 'O')
		{
	    	    FT_ScalarMemoryAlloc((POINTER*)&vortex_params,
				sizeof(VORTEX_PARAMS));
            	    front->max_time = 0.4;
            	    front->movie_frame_interval = 0.02;
            	    vortex_params->dim = 2;
            	    vortex_params->type[0] = 'M';
            	    vortex_params->cos_time = 0;
            	    vortex_params->cen[0] = 0.5;
            	    vortex_params->cen[1] = 0.25;
            	    vortex_params->rad = 0.15;
            	    vortex_params->time = 0.5*front->max_time;
            	    velo_func_pack.func_params = (POINTER)vortex_params;
            	    velo_func_pack.func = vortex_vel;
		}
		else if (string[1] == 'a' || string[1] == 'A')
		{
	    	    FT_ScalarMemoryAlloc((POINTER*)&randv_params,
				sizeof(RANDOMV_PARAMS));
		    CursorAfterString(infile,"Enter random amplitude:");
		    for (i = 0; i < dim; ++i)
		    {
        	    	fscanf(infile,"%lf",&randv_params->v0[i]);
        	    	(void) printf("%fi ",randv_params->v0[i]);
		    }
        	    (void) printf("\n");
		    CursorAfterString(infile,"Enter stop motion time:");
        	    fscanf(infile,"%lf",&randv_params->stop_time);
        	    (void) printf("%f\n",randv_params->stop_time);
            	    velo_func_pack.func_params = (POINTER)randv_params;
            	    velo_func_pack.func = random_velo;
		}
		else
		{
		    (void) printf("ERROR: need either RO or RA\n");
		    clean_up(ERROR);
		}
            	break;
            case 'd':
            case 'D':
	    	FT_ScalarMemoryAlloc((POINTER*)&dv_params,
				sizeof(BIPOLAR_PARAMS));
            	dv_params->cen1[0] = 0.25;
            	dv_params->cen1[1] = 0.25;
            	dv_params->cen2[0] = 0.75;
            	dv_params->cen2[1] = 0.25;
            	dv_params->i1 = -0.5;
            	dv_params->i2 =  0.5;
            	velo_func_pack.func_params = (POINTER)dv_params;
            	velo_func_pack.func = double_vortex_vel;
            	break;
            case 'v':
            case 'V':
	    	FT_ScalarMemoryAlloc((POINTER*)&vert_params,
				sizeof(VERTICAL_PARAMS));
		CursorAfterString(infile,"Enter center velocity:");
        	fscanf(infile,"%lf",&vert_params->v0);
        	(void) printf("%f\n",vert_params->v0);
		CursorAfterString(infile,"Enter stop motion time:");
        	fscanf(infile,"%lf",&vert_params->stop_time);
        	(void) printf("%f\n",vert_params->stop_time);
		CursorAfterString(infile,"Enter center of vertical motion:");
        	fscanf(infile,"%lf %lf",&vert_params->cen[0],
				&vert_params->cen[1]);
            	velo_func_pack.func_params = (POINTER)vert_params;
            	velo_func_pack.func = vertical_velo;
            	break;
            case 't':
            case 'T':
	    	FT_ScalarMemoryAlloc((POINTER*)&toro_params,
				sizeof(TOROIDAL_PARAMS));
		CursorAfterString(infile,"Enter center of toroidal motion:");
        	fscanf(infile,"%lf %lf %lf",&toro_params->tcen[0],
				&toro_params->tcen[1],&toro_params->tcen[2]);
        	(void) printf("%f %f %f\n",toro_params->tcen[0],
				toro_params->tcen[1],toro_params->tcen[2]);
		CursorAfterString(infile,"Enter distance to poloidal center:");
        	fscanf(infile,"%lf",&toro_params->R0);
        	(void) printf("%f\n",toro_params->R0);
		CursorAfterString(infile,"Enter velocity magnitude:");
        	fscanf(infile,"%lf",&toro_params->v0);
        	(void) printf("%f\n",toro_params->v0);
		CursorAfterString(infile,"Enter stop motion time:");
        	fscanf(infile,"%lf",&toro_params->stop_time);
        	(void) printf("%f\n",toro_params->stop_time);
            	velo_func_pack.func_params = (POINTER)toro_params;
            	velo_func_pack.func = toroidal_velo;
            	break;
            case 'p':
            case 'P':
	    	FT_ScalarMemoryAlloc((POINTER*)&para_params,
				sizeof(PARABOLIC_PARAMS));
		CursorAfterString(infile,"Enter center of parabolic velocity:");
        	fscanf(infile,"%lf %lf",&para_params->cen[0],
				&para_params->cen[1]);
        	(void) printf("%f %f\n",para_params->cen[0],
				para_params->cen[1]);
		CursorAfterString(infile,"Enter center velocity:");
        	fscanf(infile,"%lf",&para_params->v0);
        	(void) printf("%f\n",para_params->v0);
		CursorAfterString(infile,"Enter downward concavity:");
        	fscanf(infile,"%lf",&para_params->a);
        	(void) printf("%f\n",para_params->a);
		CursorAfterString(infile,"Enter stop motion time:");
        	fscanf(infile,"%lf",&para_params->stop_time);
        	(void) printf("%f\n",para_params->stop_time);
            	velo_func_pack.func_params = (POINTER)para_params;
            	velo_func_pack.func = parabolic_velo;
            	break;
            case 's':
            case 'S':
	    	FT_ScalarMemoryAlloc((POINTER*)&sing_params,
				sizeof(SINGULAR_PARAMS));
		CursorAfterString(infile,"Enter center of velocity:");
        	fscanf(infile,"%lf %lf",&sing_params->cen[0],
				&sing_params->cen[1]);
        	(void) printf("%f %f\n",sing_params->cen[0],
				sing_params->cen[1]);
		CursorAfterString(infile,"Enter center velocity:");
        	fscanf(infile,"%lf",&sing_params->v0);
        	(void) printf("%f\n",sing_params->v0);
		CursorAfterString(infile,"Enter radius of center:");
        	fscanf(infile,"%lf",&sing_params->R);
        	(void) printf("%f\n",sing_params->R);
		CursorAfterString(infile,"Enter stop motion time:");
        	fscanf(infile,"%lf",&sing_params->stop_time);
        	(void) printf("%f\n",sing_params->stop_time);
            	velo_func_pack.func_params = (POINTER)sing_params;
            	velo_func_pack.func = singular_velo;
            	break;
            case 'z':
            case 'Z':
            	velo_func_pack.func_params = NULL;
            	velo_func_pack.func = zero_velo;
            	break;
            case 'f':
            case 'F':
	    	FT_ScalarMemoryAlloc((POINTER*)&fixarea_params,
				sizeof(FIXAREA_PARAMS));
		if (string[1] == 'a' || string[1] == 'A')
		    init_fixarea_params(front,infile,fixarea_params);
		else if (string[1] == 'p' || string[1] == 'P')
		    init_fixpoint_params(front,infile,fixarea_params);
		else if (string[1] == 'f' || string[1] == 'F')
		{
            	    velo_func_pack.func_params = NULL;
            	    velo_func_pack.func = NULL;
		}
            	velo_func_pack.func_params = (POINTER)fixarea_params;
            	velo_func_pack.func = marker_velo;
            	break;
	    default:
		(void) printf("Unknown velocity function, use zero_velo()\n");
            	velo_func_pack.func_params = NULL;
            	velo_func_pack.func = zero_velo;
		break;
            }	
	}
	FT_InitVeloFunc(front,&velo_func_pack);
}	/* end initVelocityFunc */

static int zero_velo(
	POINTER params,
        Front *front,
        POINT *p,
        HYPER_SURF_ELEMENT *hse,
        HYPER_SURF *hs,
        double *vel)
{
	vel[0] = vel[1] = vel[2] = 0.0;
	return YES;
}	/* end zero_velo */

static int random_velo(
	POINTER params,
        Front *front,
        POINT *p,
        HYPER_SURF_ELEMENT *hse,
        HYPER_SURF *hs,
        double *vel)
{
	RANDOMV_PARAMS *randv_params = (RANDOMV_PARAMS*)params;
	double *v0 = randv_params->v0;
	double stop_time = randv_params->stop_time;
	unsigned short int xsubi[3];

	if (front->time >= stop_time)
	{
	    vel[0] = vel[1] = vel[2] = 0.0;
	}
	else
	{
	    xsubi[0] = 7256;
	    xsubi[1] = 764;
	    xsubi[2] = 2163;
	    vel[0] = v0[0]*(2.0*erand48(xsubi) - 1.0);
	    vel[1] = v0[1]*(2.0*erand48(xsubi) - 1.0);
	    vel[2] = v0[2]*(2.0*erand48(xsubi) - 1.0);
	}
	return YES;
}	/* end random_velo */

static int vertical_velo(
        POINTER params,
        Front *front,
        POINT *p,
        HYPER_SURF_ELEMENT *hse,
        HYPER_SURF *hs,
        double *vel)
{
	VERTICAL_PARAMS *vert_params = (VERTICAL_PARAMS*)params;
        double *coords = Coords(p);
	double dist,v_vert;
	double v0 = vert_params->v0;
	double stop_time = vert_params->stop_time;

	dist = sqrt(sqr(coords[0] - 0.5) + sqr(coords[1] - 0.5));
	v_vert = (0.15 - dist)/0.15;
	vel[0] = vel[1] = 0.0;
	if (front->time < stop_time)
	    vel[2] = v_vert*v0;
	else
	    vel[2] = 0.0;
	return YES;
}       /* end vertical_velo */

static int toroidal_velo(
        POINTER params,
        Front *front,
        POINT *p,
        HYPER_SURF_ELEMENT *hse,
        HYPER_SURF *hs,
        double *vel)
{
	TOROIDAL_PARAMS *toro_params = (TOROIDAL_PARAMS*)params;
        double *coords = Coords(p);
	int i,dim = front->rect_grid->dim;
	double v0 = toro_params->v0;
	double stop_time = toro_params->stop_time;
	double tcoords[2]; 		/* toroidal coords */
	double pvel[2];		/* projected velocity */
	double *tcen = toro_params->tcen; /* toroidal center of vel */
	double R0 = toro_params->R0; /* radial dist of poloidal center */
	double d1,d2;
	double s1,s2;
	double dx1,dy1;
	double dx2,dy2;
	double dtol = 0.000001*front->rect_grid->h[0];

	if (front->time >= stop_time)
	{
	    for (i = 0; i < dim; ++i)
		vel[i] = 0.0;
	    return YES;
	}
	/* Project 3D to 2D */
	tcoords[0] = sqrt(sqr(coords[0] - tcen[0]) + sqr(coords[1] - tcen[1]));
	tcoords[1] = coords[2] - tcen[2];

	dx1 = tcoords[0] - R0;
	dx2 = tcoords[0] + R0;
	dy1 = dy2 = tcoords[1];

	d1 = sqr(dx1) + sqr(dy1);
	d2 = sqr(dx2) + sqr(dy2);
	s1 = v0;
	s2 = -v0;
	if (d1 < dtol || d2 < dtol)
	    pvel[0] = pvel[1] = 0.0;
	else if (front->time < stop_time)
	{
	    pvel[0] =  s1*dy1/d1 + s2*dy2/d2;
	    pvel[1] = -s1*dx1/d1 - s2*dx2/d2;
	} 
	else
	    pvel[0] = pvel[1] = 0.0;
	    
	/* give it back to 3D */
	vel[0] = pvel[0]*(coords[0]-tcen[0])/tcoords[0];
	vel[1] = pvel[0]*(coords[1]-tcen[1])/tcoords[0];
	vel[2] = pvel[1];
	return YES;
}       /* end toroidal_velo */


static int parabolic_velo(
        POINTER params,
        Front *front,
        POINT *p,
        HYPER_SURF_ELEMENT *hse,
        HYPER_SURF *hs,
        double *vel)
{
	PARABOLIC_PARAMS *para_params = (PARABOLIC_PARAMS*)params;
        double *coords = Coords(p);
	int i,dim = front->rect_grid->dim;
	double v0 = para_params->v0;
	double a = para_params->a;
	double *cen = para_params->cen;
	double R_sqr = 0.0;

	for (i = 0; i < dim-1; ++i)
	{
	    R_sqr += sqr(coords[i] - cen[i]);
	    vel[i] = 0.0;
	}
	vel[dim-1] = -0.5*a*R_sqr + v0;
	return YES;
}	/* end parabolic_velo */

static int singular_velo(
        POINTER params,
        Front *front,
        POINT *p,
        HYPER_SURF_ELEMENT *hse,
        HYPER_SURF *hs,
        double *vel)
{
	SINGULAR_PARAMS *para_params = (SINGULAR_PARAMS*)params;
        double *coords = Coords(p);
	int i,dim = front->rect_grid->dim;
	double v0 = para_params->v0;
	double R = para_params->R;
	double *cen = para_params->cen;
	double r = 0.0;
	for (i = 0; i < dim-1; ++i)
	{
	    r += sqr(coords[i] - cen[i]);
	    vel[i] = 0.0;
	}
	r = sqrt(r);
	if (r < R)
	    vel[dim-1] = v0;
	else
	    vel[dim-1] = 0.0;
	return YES;
}	/* end sigular_velo */

struct _SHAPE_PARAMS {
	int shape_id;
	double L[2];
	double U[2];
	double cen[2];
	double R[2];
};
typedef struct _SHAPE_PARAMS SHAPE_PARAMS;

static boolean within_shape(SHAPE_PARAMS,double*);

static void init_fixarea_params(
	Front *front,
	FILE *infile,
	FIXAREA_PARAMS *fixarea_params)
{
	char string[100];
	SHAPE_PARAMS sparams;
	int i,num_pts,dim = front->rect_grid->dim;
	static REGISTERED_PTS *registered_pts;
	POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
	INTERFACE *intfc = front->interf;
	SURFACE *surf;

	(void) printf("Available initial areas are:\n");
	(void) printf("\tRectangle (R)\n");
	(void) printf("\tEllipse (E)\n");
	CursorAfterString(infile,"Enter initial shape of fixed area:");
        fscanf(infile,"%s",string);
        (void) printf("%s\n",string);
	switch (string[0])
	{
	case 'r':
	case 'R':
	    sparams.shape_id = 0;
	    CursorAfterString(infile,"Enter rectangle lower bounds:");
            fscanf(infile,"%lf %lf",&sparams.L[0],&sparams.L[1]);
            (void) printf("%f %f\n",sparams.L[0],sparams.L[1]);
	    CursorAfterString(infile,"Enter rectangle upper bounds:");
            fscanf(infile,"%lf %lf",&sparams.U[0],&sparams.U[1]);
            (void) printf("%f %f\n",sparams.U[0],sparams.U[1]);
	    break;
	case 'e':
	case 'E':
	    sparams.shape_id = 1;
	    CursorAfterString(infile,"Enter center of ellipse:");
            fscanf(infile,"%lf %lf",&sparams.cen[0],&sparams.cen[1]);
            (void) printf("%f %f\n",sparams.cen[0],sparams.cen[1]);
	    CursorAfterString(infile,"Enter radii of ellipse:");
            fscanf(infile,"%lf %lf",&sparams.R[0],&sparams.R[1]);
            (void) printf("%f %f\n",sparams.R[0],sparams.R[1]);
	    break;
	}
	CursorAfterString(infile,"Enter area velocity:");
        for (i = 0; i < dim; ++i)
        {
            fscanf(infile,"%lf",&fixarea_params->vel[i]);
            (void) printf("%f ",fixarea_params->vel[i]);
        }
        (void) printf("\n");
	num_pts = 0;
	next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
            if (wave_type(hs) != ELASTIC_BOUNDARY) continue;
	    if (within_shape(sparams,Coords(p)))
	    {
		num_pts++;	
	    }
	}
	FT_VectorMemoryAlloc((POINTER*)&fixarea_params->global_ids,num_pts,
				sizeof(int));
	FT_ScalarMemoryAlloc((POINTER*)&registered_pts,sizeof(REGISTERED_PTS));
	FT_VectorMemoryAlloc((POINTER*)&registered_pts->global_ids,num_pts,
				sizeof(int));
	fixarea_params->num_pts = num_pts;
	registered_pts->num_pts = num_pts;

	num_pts = 0;
	next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
            if (wave_type(hs) != ELASTIC_BOUNDARY) continue;
	    if (within_shape(sparams,Coords(p)))
	    {
		fixarea_params->global_ids[num_pts] = Gindex(p);	
		registered_pts->global_ids[num_pts] = Gindex(p);	
		num_pts++;
	    	surf = Surface_of_hs(hs);
	    	surf->extra = (POINTER)registered_pts;
	    }
	}
}	/* end init_fixarea_params */

static boolean within_shape(
	SHAPE_PARAMS sparams,
	double *coords)
{
	double dist;
	switch (sparams.shape_id)
	{
	case 0:
	    if (coords[0] < sparams.L[0]) return NO;
	    if (coords[0] > sparams.U[0]) return NO;
	    if (coords[1] < sparams.L[1]) return NO;
	    if (coords[1] > sparams.U[1]) return NO;
	    return YES;
	case 1:
	    dist = sqr((coords[0] - sparams.cen[0])/sparams.R[0]) +
	           sqr((coords[1] - sparams.cen[1])/sparams.R[1]);
	    if (dist < 1.0) return YES;
	    else return NO;
	}
}	/* within_shape */

static void init_fixpoint_params(
	Front *front,
	FILE *infile,
	FIXAREA_PARAMS *fixarea_params)
{
	int i,dim = front->rect_grid->dim;
	static REGISTERED_PTS *registered_pts;
	POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
	INTERFACE *intfc = front->interf;
	SURFACE *surf;
	double coords[MAXD];
	double dist,min_dist;

	FT_VectorMemoryAlloc((POINTER*)&fixarea_params->global_ids,1,
				sizeof(int));
	FT_ScalarMemoryAlloc((POINTER*)&registered_pts,sizeof(REGISTERED_PTS));
	FT_VectorMemoryAlloc((POINTER*)&registered_pts->global_ids,1,
				sizeof(int));
	fixarea_params->num_pts = 1;
	registered_pts->num_pts = 1;

	CursorAfterString(infile,"Enter initial fixed point coordinates:");
        for (i = 0; i < dim; ++i)
        {
            fscanf(infile,"%lf",&coords[i]);
            (void) printf("%f ",coords[i]);
        }
        (void) printf("\n");
	CursorAfterString(infile,"Enter point velocity:");
        for (i = 0; i < dim; ++i)
        {
            fscanf(infile,"%lf",&fixarea_params->vel[i]);
            (void) printf("%f ",fixarea_params->vel[i]);
        }
        (void) printf("\n");

	min_dist = HUGE;
	next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
            if (wave_type(hs) != ELASTIC_BOUNDARY) continue;
	    dist = distance_between_positions(coords,Coords(p),dim);
	    if (dist < min_dist)
	    {
		min_dist = dist;
		registered_pts->global_ids[0] = Gindex(p);
		fixarea_params->global_ids[0] = Gindex(p);
	    	surf = Surface_of_hs(hs);
	    	surf->extra = (POINTER)registered_pts;
	    }
        }
}	/* end init_fixarea_params */

static int marker_velo(
        POINTER params,
        Front *front,
        POINT *p,
        HYPER_SURF_ELEMENT *hse,
        HYPER_SURF *hs,
        double *vel)
{
	FIXAREA_PARAMS *fixarea_params = (FIXAREA_PARAMS*)params;
	int i,j,dim = front->rect_grid->dim;
	int num_pts = fixarea_params->num_pts;
	int *global_ids = fixarea_params->global_ids;
	for (i = 0; i < num_pts; ++i)
	{
	    if (Gindex(p) == global_ids[i])
	    {
		for (j = 0; j < dim; ++j)
		    vel[j] = fixarea_params->vel[j];
		return YES;
	    }
	}
	for (j = 0; j < dim; ++j)
	    vel[j] = 0.0;
	return YES;
}	/* end marker_velo */

extern void resetFrontVelocity(Front *front)
{
	INTERFACE *intfc = front->interf;
	POINT *p;
	HYPER_SURF_ELEMENT *hse;
	HYPER_SURF *hs;
	STATE *sl,*sr;
	int i,dim = front->rect_grid->dim;
	CURVE **c;
	BOND *b;

	next_point(intfc,NULL,NULL,NULL);
	while (next_point(intfc,&p,&hse,&hs))
	{
	    FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
	    for (i = 0; i < dim; ++i)
	    {
		p->vel[i] = 0.0;
		sl->vel[i] = sr->vel[i] = 0.0;
		sl->impulse[i] = sr->impulse[i] = 0.0;
	    }
	}
	if (dim == 3)
	{
	    for (c = intfc->curves; c && *c; ++c)
	    {
		p = (*c)->start->posn;
		sl = (STATE*)left_state(p);
		sr = (STATE*)right_state(p);
	        for (i = 0; i < dim; ++i)
		{
		    p->vel[i] = 0.0;
		    sl->vel[i] = sr->vel[i] = 0.0;
		    sl->impulse[i] = sr->impulse[i] = 0.0;
		}
		for (b = (*c)->first; b != (*c)->last; b = b->next)
		{
		    p = b->end;
		    sl = (STATE*)left_state(p);
		    sr = (STATE*)right_state(p);
	            for (i = 0; i < dim; ++i)
		    {
		    	p->vel[i] = 0.0;
		    	sl->vel[i] = sr->vel[i] = 0.0;
		    	sl->impulse[i] = sr->impulse[i] = 0.0;
		    }
		}
		p = (*c)->end->posn;
		sl = (STATE*)left_state(p);
		sr = (STATE*)right_state(p);
	        for (i = 0; i < dim; ++i)
		{
		    p->vel[i] = 0.0;
		    sl->vel[i] = sr->vel[i] = 0.0;
		    sl->impulse[i] = sr->impulse[i] = 0.0;
		}
	    }
	}
}	/* end resetFrontVelocity */

static void checkSetGoreNodes(
	INTERFACE *intfc)
{
	CURVE **c;
	intfc_curve_loop(intfc,c)
	{
	    if (hsbdry_type(*c) == GORE_HSBDRY)
	    {
		set_gore_node((*c)->start);
		set_gore_node((*c)->end);
	    }
	}
}	/* end checkSetGoreNodes */


static void set_gore_node(
	NODE *n)
{
	static AF_NODE_EXTRA *extra;
	boolean is_gore_node = NO;
	CURVE **c;

	for (c = n->in_curves; c && *c; ++c)
	    if (hsbdry_type(*c) == STRING_HSBDRY)
		return;
	    else if (hsbdry_type(*c) == GORE_HSBDRY)
		is_gore_node = YES;
	for (c = n->out_curves; c && *c; ++c)
	    if (hsbdry_type(*c) == STRING_HSBDRY)
		return;
	    else if (hsbdry_type(*c) == GORE_HSBDRY)
		is_gore_node = YES;
	if (!is_gore_node) return;
	if (n->extra == NULL)
	{
	    if (extra ==  NULL)
	    {
	    	FT_ScalarMemoryAlloc((POINTER*)&extra,sizeof(AF_NODE_EXTRA));
		extra->af_node_type = GORE_NODE;
	    }
	    n->extra = (POINTER)extra;
	}
}	/* end set_gore_node */

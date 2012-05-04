#include <FronTier.h>

#include "curvature.h"

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
        double R[MAXD]; 
} TEST_CIRCLE_PARAMS;

	/*  Function Declarations */
static void set_geom_functions(char*,Front*);
static void test_normal(POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*,double*,Front*);

static void test_normal2d(POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*,Front*);
static void test_normal3d(POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*,Front*);

static double test_curvature2d(POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*,Front*);
static double test_curvature3d(POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*,Front*);

static double level_hyper_sphere(POINTER,double*);
static double level_torus_func(POINTER,double*);
static double level_ellipse_func(POINTER,double*);

static void reset_hyper_sphere_points(Front*,TEST_CIRCLE_PARAMS);
static void reset_torus_points(Front*,TEST_CIRCLE_PARAMS);
static void reset_ellipse_points(Front*,TEST_CIRCLE_PARAMS);

static void print_hyper_error(Front*,char*,TEST_CIRCLE_PARAMS);

static double test_curvature(POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*,Front*);

enum {
	TEST_HYPER_SPHERE		= 	1,
	TEST_ELLIPSE,
	TEST_TORUS
};
int method;

LOCAL   double   compute_curvature2d1(POINT*,POINT**,int);

LOCAL   void obtain_orthbases(double*,double*);

LOCAL void FrontFirstOrderNorm(Front*);

int main(int argc, char **argv)
{
	static Front front;
	static RECT_GRID comp_grid;
	static F_BASIC_DATA f_basic;
	static LEVEL_FUNC_PACK level_func_pack;
	TEST_CIRCLE_PARAMS circle_params;	/* level function parameters */
	int i,dim;
	FILE *infile;
	char s[100];

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
	circle_params.dim = dim;

	infile = fopen(in_name,"r");
	CursorAfterString(infile,"Enter test case:");
	fscanf(infile,"%s",s);
	switch (s[0])
	{
	case 'C':
	case 'S':
	    method = TEST_HYPER_SPHERE;
	    CursorAfterString(infile,"Enter center:");
	    for (i = 0; i < dim; ++i)
	    	fscanf(infile,"%lf",&circle_params.center[i]);
	    CursorAfterString(infile,"Enter radius:");
	    	fscanf(infile,"%lf",&circle_params.R[0]);
	    level_func_pack.func = level_hyper_sphere;
	    break;
	case 'E':
	    method = TEST_ELLIPSE;
            CursorAfterString(infile,"Enter center:");
            for (i = 0; i < dim; ++i)
                fscanf(infile,"%lf",&circle_params.center[i]);
            CursorAfterString(infile,"Enter radius:");
            for (i = 0; i < dim; ++i)
                fscanf(infile,"%lf",&circle_params.R[i]);
            
            printf("%f %f %f\n",circle_params.R[0],circle_params.R[1],circle_params.R[2]);    
            
            level_func_pack.func = level_ellipse_func;
	    break;
	case 'T':
	    method = TEST_TORUS;
            CursorAfterString(infile,"Enter center:");
            for (i = 0; i < dim; ++i)
                fscanf(infile,"%lf",&circle_params.center[i]);
            CursorAfterString(infile,"Enter radius:");
            for (i = 0; i < dim-1; ++i)
                fscanf(infile,"%lf",&circle_params.R[i]);
	    level_func_pack.func = level_torus_func;
	}

	    /* Initialize interface through level function */

	level_func_pack.neg_component = 1;
	level_func_pack.pos_component = 2;
	level_func_pack.func_params = (POINTER)&circle_params;

        level_func_pack.wave_type = FIRST_PHYSICS_WAVE_TYPE;
        
	FT_InitIntfc(&front,&level_func_pack);
        FT_Save(&front,out_name);
	switch (method)
	{
	case TEST_HYPER_SPHERE:
	    reset_hyper_sphere_points(&front,circle_params);
	    break;
	case TEST_ELLIPSE:
            reset_ellipse_points(&front,circle_params);
	    break;
	case TEST_TORUS:
	    reset_torus_points(&front,circle_params);
	    break;
	}
        set_geom_functions(in_name,&front);
	if(dim == 3) FrontFirstOrderNorm(&front);
	print_hyper_error(&front,out_name,circle_params);
	clean_up(0);
}

LOCAL void FrontFirstOrderNorm(Front *front)
{
        INTERFACE *intfc = front->interf;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
        POINT *p;
        next_point(intfc,NULL,NULL,NULL);

        while (next_point(intfc,&p,&hse,&hs))
        {
            normal(p,hse,hs,p->_nor0,front);
            //GetFrontNormal(p,hse,hs,p->_nor,front);
        }
}     /* end IntfcFirstOrderNorm */

/********************************************************************
*              Torus level function for the initial interface       *
 ********************************************************************/
static double level_torus_func(
        POINTER func_params,
        double *coords)
{
        TEST_CIRCLE_PARAMS *circle_params = (TEST_CIRCLE_PARAMS*)func_params;
        double *center,R,dist;
        int i,dim;
        double x0,y0,z0,x,y,z,r;
        double distance;
        
        x = coords[0];
        y = coords[1];
        z = coords[2];
        R = circle_params->R[0];
        r = circle_params->R[1];
        
        distance = sqrt(sqr(R-sqrt(sqr(x)+sqr(y))) +sqr(z)) - r;
        return distance;
          
}   /* end torus_func */


static void reset_torus_points(
        Front *front,
        TEST_CIRCLE_PARAMS circle_params)
{
        POINT              *p;
        HYPER_SURF_ELEMENT *hse;
        HYPER_SURF         *hs;
        INTERFACE          *intfc = front->interf;
        double             *center = circle_params.center;
        double             r, *R = circle_params.R;
        int                i,dim = circle_params.dim;
        
        double                   norm,x,y,z,x0,y0,cnt[3],nrm[3];

        r = R[1];
        next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
            if (wave_type(hs) < FIRST_PHYSICS_WAVE_TYPE) continue;
            
            x = Coords(p)[0];
            y = Coords(p)[1];
            z = Coords(p)[2];
            
            cnt[0] = x/sqrt(x*x + y*y)*R[0];
            cnt[1] = y/sqrt(x*x + y*y)*R[0];
            cnt[2] = 0;
            norm = sqrt((x - cnt[0])*(x - cnt[0]) + (y - cnt[1])*(y - cnt[1])
                                    +(z - cnt[2])*(z - cnt[2]));
            nrm[0] = (x - cnt[0])/norm;
            nrm[1] = (y - cnt[1])/norm;
            nrm[2] = (z - cnt[2])/norm;
            
            Coords(p)[0] = cnt[0] + r *nrm[0];
            Coords(p)[1] = cnt[1] + r *nrm[1];
            Coords(p)[2] = cnt[2] + r *nrm[2];
        }
}       /* end adjust_torus_points */

/********************************************************************
*              Ellipse & ellipsoid  level function for the initial interface       *
 ********************************************************************/
static double level_ellipse_func(
        POINTER func_params,
        double *coords)
{
        TEST_CIRCLE_PARAMS *circle_params = (TEST_CIRCLE_PARAMS*)func_params;
        double *center,*R,dist;
        int i,dim;

        center = circle_params->center;
        R  = circle_params->R;
        dim  = circle_params->dim;

        dist = 0.0;
        for (i = 0; i < dim; ++i)
            dist += sqr(coords[i] - center[i])/sqr(R[i]);

        dist = sqrt(dist) - 1;
        return dist;
}  /* end ellipse_func */

static void reset_ellipse_points(
        Front *front,
        TEST_CIRCLE_PARAMS circle_params)
{
        POINT              *p;
        HYPER_SURF_ELEMENT *hse;
        HYPER_SURF         *hs;
        INTERFACE          *intfc = front->interf;
        
        double             *center = circle_params.center;
        double             r, b, *R = circle_params.R; //a = 2*r; b = r;
        int                i,dim = circle_params.dim;
        double             norm,x,y,z,theta; /// z - - for ellipsoid
        
        b = R[1];
        next_point(intfc,NULL,NULL,NULL);
        if(dim == 2)
        {
            while (next_point(intfc,&p,&hse,&hs))
            {
                if (wave_type(hs) < FIRST_PHYSICS_WAVE_TYPE) continue;
                if (Boundary_point(p) || Boundary_hs(hs)) continue;

                x = Coords(p)[0];
                y = Coords(p)[1];
                theta = atan((y - center[1])/(x - center[0]));

                if((x - center[0]) >= 0.0 && (y - center[1]) >= 0.0)
                {
                   Coords(p)[0] = 2.0*b*cos(theta) + center[0];
                   Coords(p)[1] = b*sin(theta) + center[1];
                }
                if((x - center[0]) < 0.0 && (y - center[1]) < 0.0)
                {
                   Coords(p)[0] = -2.0*b*cos(theta) + center[0];
                   Coords(p)[1] = -b*sin(theta) + center[1];
                }
                if((x - center[0]) < 0.0 && (y - center[1]) > 0.0)
                {
                   Coords(p)[0]= -2.0*b*cos(theta) + center[0];
                   Coords(p)[1]= -b*sin(theta) + center[1];
                }
                if((x - center[0]) > 0.0 && (y - center[1]) < 0.0)
                {
                   Coords(p)[0] = 2.0*b*cos(theta) + center[0];
                   Coords(p)[1] = b*sin(theta) + center[1];
                }
             }
         }else
         {
             while (next_point(intfc,&p,&hse,&hs))
             {   
                 x = Coords(p)[0];
                 y = Coords(p)[1];
                 z = Coords(p)[2];
                 theta = sqr(x - center[0])/sqr(R[0]) 
                           + sqr(y - center[1])/sqr(R[1])     
                           + sqr(z - center[2])/sqr(R[2]) - 1.0;
                 Coords(p)[0] = sqrt(1 + theta)*x;
                 Coords(p)[1] = sqrt(1 + theta)*y;
                 Coords(p)[2] = sqrt(1 + theta)*z;
             } 
         } 
}       /* end adjust_ellipse_points */


/********************************************************************
 *	Hyper sphere level function for the initial interface    *
 ********************************************************************/

static double level_hyper_sphere(
        POINTER func_params,
        double *coords)
{
	TEST_CIRCLE_PARAMS *circle_params = (TEST_CIRCLE_PARAMS*)func_params;
	double *center,*R,dist;
	int i,dim;

	center = circle_params->center;
	R  = circle_params->R;
	dim  = circle_params->dim;

	dist = 0.0;
	for (i = 0; i < dim; ++i)
	    dist += sqr(coords[i] - center[i]);	
	dist = sqrt(dist) - R[0];
	return dist;
}	/* end level_hyper_sphere */

static void reset_hyper_sphere_points(
	Front *front,
	TEST_CIRCLE_PARAMS circle_params)
{
	POINT              *p;
        HYPER_SURF_ELEMENT *hse;
        HYPER_SURF         *hs;
	INTERFACE	   *intfc = front->interf;
	double	 	   *center = circle_params.center;
	double	 	   r, R = circle_params.R[0];
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
}	/* end reset_hyper_sphere_points */

static void set_geom_functions(
	char *in_name,
	Front *front)
{
	INTERFACE *intfc = front->interf;
	int dim = intfc->dim;
	FILE *infile = fopen(in_name,"r");
	TANGENT_METHOD t_method;
	NORMAL_METHOD n_method;
	CURVATURE_METHOD c_method;
	char s[100];

	CursorAfterString(infile,"Enter yes to test local function:");
	fscanf(infile,"%s",s);
	if (s[0] == 'Y' || s[0] == 'y')
	{
	    interface_normal(intfc) = test_normal;
	    interface_curvature(intfc) = test_curvature;
	    fclose(infile);
	    return;
	}
	if (dim == 2)
	{
	    CursorAfterString(infile,"Enter normal method:");
	    fscanf(infile,"%s",s);
	    if (s[0] == 'F' || s[0] == 'f')
	    {
		n_method = FIRST_ORDER_NORMAL;
	    	CursorAfterString(infile,"Enter tangent method:");
	    	fscanf(infile,"%s",s);
		if (s[1] == 'I' || s[1] == 'i')
		    t_method = LINEAR_SECANT;
		else if (s[1] == 'A' || s[1] == 'a')
		    t_method = LANGRANGIAN_INTERPOLANT;
		else
		{
		    screen("Unknown tangent method!\n");
		    clean_up(ERROR);
		}
		c_method = NORMAL_CURVATURE;
	    }
	    else if (s[0] == 'W' || s[0] == 'w')
	    {
		n_method = WLSP_NORMAL;
		t_method = WLSP_TANGENT;
		c_method = WLSP_CURVATURE;
	    }
	}
	else if (dim == 3)
	{
	    CursorAfterString(infile,"Enter normal method:");
	    if (s[0] == 'W' || s[0] == 'w')
	    {
		n_method = WLSP_NORMAL;
		c_method = WLSP_CURVATURE;
	    }
	    else if (s[0] == 'P' || s[0] == 'p')
	    {
		n_method = PLANE_FIT_NORMAL;
		c_method = NORMAL_CURVATURE;
	    }
	    else if (s[0] == 'S' || s[0] == 's')
	    {
		n_method = SINE_WEIGHTED_NORMAL;
		c_method = NORMAL_CURVATURE;
	    }
	    else if (s[0] == 'A' || s[0] == 'a')
	    {
		n_method = AREA_WEIGHTED_NORMAL;
		c_method = NORMAL_CURVATURE;
	    }
	}
	fclose(infile);
	FrontSetGeomVarMethods(front,t_method,n_method,c_method);
	return;
}

static void test_normal(
        POINT              *p,
        HYPER_SURF_ELEMENT *hse,
        HYPER_SURF         *hs,
        double              *nor,
        Front              *front)
{
	int i,dim = front->rect_grid->dim;
	switch (front->rect_grid->dim)
        {
        case 2:
            test_normal2d(p,hse,hs,front);
            break;
        case 3:
            test_normal3d(p,hse,hs,front);
            break;
        }
	for (i = 0; i < dim; ++i) 
	    nor[i] = p->_nor[i];
}	/* end test_normal */

static double test_curvature(
        POINT              *p,
        HYPER_SURF_ELEMENT *hse,
        HYPER_SURF         *hs,
        Front              *fr)
{
	switch (fr->rect_grid->dim)
	{
	case 2:
	    return test_curvature2d(p,hse,hs,fr);
	case 3:
	    return test_curvature3d(p,hse,hs,fr);
	}
	
}	/* end test_curvature */

////////////////Least square 2D ////////////////////

/*******************************************************************
 * Function for the curvature computation of seven points with LSQ *
 *******************************************************************/
/* Prototypes for local functions. */
/*************************************************************
 *
 *    FUNCTION: safeqr
 *
 *    Perform Householder QR factorization to compute
 *    [Q, A(1:n,1:n)] = qr(A(1:m,1:n),0); b = Q'*b; rnk=n;
 **************************************************************/
int safeqr(
        double *A,
        int *A_dim1,
        int *A_dim2,
        double *b,
        int *b_dim1,
        double tol)
{
        double v[13];
        int v_dim1;
        static double dtemp_0[91];
        int dtemp_0_dim1;
        double dtemp_1[13];
        int dtemp_1_dim1;
        int rnk_out, m, n, itemp_0, itemp_1;
        double vnrm2, dtemp_2, dtemp_3, dtemp_4, dtemp_5;
        int itemp_2, itemp_3, itemp_4, itemp_5, itemp_6;
        int itemp_7, itemp_8, itemp_9, itemp_10, itemp_11;
        int itemp_12, itemp_13, itemp_14, itemp_15, itemp_16;
        int itemp_17, itemp_18, itemp_19, itemp_20, itemp_21;
        int itemp_22, itemp_23, itemp_24, itemp_25, itemp_26;
        int itemp_27, i, j, i3;

        n = *A_dim2;
        m = *A_dim1;
        rnk_out = n;
        itemp_0 = n;
        for(i = 0; i <= itemp_0 - 1; i += 1)
        {
            ct_set_max(v_dim1, m - i, 0);
            for(j = 0; j <= m - 1 - i; j += 1)
            {
                v[j] = A[(i + j)*(*A_dim2) + i];
            }
            if(v[0] >= 0.0)
            {
                ct_set_max(itemp_11, m - i, 0);
                if(v_dim1 > 1)
                {
                    itemp_1 = itemp_11;
                }
                else
                {
                    itemp_1 = 1;
                }
                ct_set_max(itemp_12, m - i, 0);
                if(v_dim1 > 1)
                {
                    itemp_2 = itemp_12;
                }
                else
                {
                    itemp_2 = 1;
                }
                dtemp_3 = 0.0;

                for(j = 0; j <= itemp_1 - 1; j += 1)
                {
                    dtemp_3 = dtemp_3 + v[j]*v[j];
                }
                v[0] =  v[0] + sqrt(dtemp_3);
            }
            else
            {
                ct_set_max(itemp_13, m - i, 0);
                if (v_dim1 > 1)
                {
                    itemp_3 = itemp_13;
                }
                else
                {
                    itemp_3 = 1;
                }

                ct_set_max(itemp_14, m - i, 0);
                if (v_dim1 > 1)
                {
                    itemp_4 = itemp_14;
                }
                else
                {
                    itemp_4 = 1;
                }
                dtemp_4 = 0.0;
                for(j = 0; j <= itemp_3 - 1; j += 1)
                {
                    dtemp_4 = dtemp_4 + v[j]*v[j];
                }
                v[0] = v[0] - sqrt(dtemp_4);
            }
            ct_set_max(itemp_15, m - i, 0);
            if(v_dim1 > 1)
            {
                itemp_5 = itemp_15;
            }
            else
            {
                itemp_5 = 1;
            }
            ct_set_max(itemp_16, m - i, 0);
            if(v_dim1 > 1)
            {
                itemp_6 = itemp_16;
            }
            else
            {
                itemp_6 = 1;
            }
            dtemp_5 = 0.0;
            for(j = 0; j <= itemp_5 - 1; j += 1)
            {
                dtemp_5 = dtemp_5 + v[j]*v[j];
            }
            vnrm2 = sqrt(dtemp_5);
            if(vnrm2 > 0.0)
            {
                for(j = 0; j <= v_dim1 - 1; j += 1)
                {
                    v[j] = v[j]/vnrm2;
                }
            }
            for(j = 1 + i; j <= n; j += 1)
            {
                ct_set_max(itemp_17, m - i, 0);
                if(v_dim1 > 1)
                {
                    itemp_7 = itemp_17;
                }
                else
                {
                    itemp_7 = 1;
                }
                ct_set_max(itemp_18, m - i, 0);
                dtemp_2 = 0.0;
                for(i3 = 0; i3 <= itemp_7 - 1; i3 += 1)
                {
                    dtemp_2 = dtemp_2 + v[i3] * A[(i + i3)*(*A_dim2) - 1 + j];
                }
                ct_set_max(itemp_19, m - i, 0);
                ct_set_max(itemp_20, m - i, 0);
                if(v_dim1 > 1)
                {
                    itemp_8 = itemp_20;
                }
                else
                {
                    itemp_8 = 1;
                }
                if(itemp_19 == 1)
                {
                    dtemp_0_dim1 = itemp_8;
                    for(i3 = 0; i3 <= dtemp_0_dim1 - 1; i3 += 1)
                    {
                        dtemp_0[i3] = A[i*(*A_dim2) - 1 + j]
                                            - 2.0*v[i3]*dtemp_2;
                    }
                }
                else
                {
                    dtemp_0_dim1 = itemp_19;
                    if (itemp_8 == 1)
                    {
                        for(i3 = 0; i3 <= dtemp_0_dim1 - 1; i3 += 1)
                        {
                            dtemp_0[i3] = A[(i + i3)*(*A_dim2) - 1 + j]
                                                 - 2.0*v[0]*dtemp_2;
                        }
                    }
                    else
                    {
                        for(i3 = 0; i3 <= dtemp_0_dim1 - 1; i3 += 1)
                        {
                            dtemp_0[i3] = A[(i + i3)*(*A_dim2) - 1 + j]
                                                 - 2.0*v[i3] *dtemp_2;
                        }
                    }
                }
                ct_set_max(itemp_21, m - i, 0);
                itemp_27 = 0;
                if(dtemp_0_dim1 == 1)
                {
                    for(i3=0; i3 <= itemp_21 - 1; i3 += 1)
                    {
                        A[(i + i3)*(*A_dim2) - 1 + j] = dtemp_0[0];
                    }
                }
                else
                {
                    for(i3 = 0; i3 <= itemp_21 - 1; i3 += 1)
                    {
                        A[(i + i3)*(*A_dim2) - 1 + j] = dtemp_0[itemp_27];
                        itemp_27 = 1 + itemp_27;
                    }
                }
            }
            /* Estimate rank of matrix */
            if((fabs(A[i*(*A_dim2) + i]) < tol)&&(rnk_out == n))
            {
                rnk_out = i;
                if(rnk_out > 3)
                {
                    return rnk_out;
                }
            }
            ct_set_max(itemp_22, m - i, 0);
            if(v_dim1 > 1)
            {
                itemp_9 = itemp_22;
            }
            else
            {
                itemp_9 = 1;
            }
            ct_set_max(itemp_23, m - i, 0);
            dtemp_2 = 0.0;
            for(j = 0; j <= itemp_9 - 1; j += 1)
            {
                dtemp_2 = dtemp_2 + v[j]*b[i + j];
            }
            ct_set_max(itemp_24, m - i, 0);
            ct_set_max(itemp_25, m - i, 0);
            if(v_dim1 > 1)
            {
                itemp_10 = itemp_25;
            }
            else
            {
                itemp_10 = 1;
            }
            if (itemp_24 == 1)
            {
                dtemp_1_dim1 = itemp_10;
                for(j = 0; j <= dtemp_1_dim1 - 1; j += 1)
                {
                    dtemp_1[j] = b[i] - 2.0*v[j]*dtemp_2;
                }
            }
            else
            {
                dtemp_1_dim1 = itemp_24;
                if(itemp_10 == 1)
                {
                    for(j = 0; j <= dtemp_1_dim1 - 1; j += 1)
                    {
                        dtemp_1[j] =  b[i + j] - 2.0*v[0]*dtemp_2;
                    }
                }
                else
                {
                    for(j = 0; j <= dtemp_1_dim1 - 1; j += 1)
                    {
                        dtemp_1[j] =  b[i + j] - 2.0*v[j]*dtemp_2;
                    }
                }
            }
            ct_set_max(itemp_26, m - i, 0);
            itemp_27 = 0;
            if(dtemp_1_dim1 == 1)
            {
                for(j = 0; j <= itemp_26 - 1; j += 1)
                {
                    b[i + j] = dtemp_1[0];
                }
            }
            else
            {
                for(j = 0; j <= itemp_26 - 1; j += 1)
                {
                    b[i + j] = dtemp_1[itemp_27];
                    itemp_27 = 1 + itemp_27;
                }
            }
        }
        return rnk_out;
}

/*************************************************************
 *                FUNCTION: polyfit2d_lhf                    *
 *************************************************************/
void  polyfit2d_lhf(
        double normal_out[2],
        double *curvature_out,
        double *ngbpts,
        int ngbpts_dim1,
        int mid,
        int deg)
{
        static double newpts[26];
        static double AS[91];
        int AS_dim1, AS_dim2;
        double b[13];
        int b_dim1;
        double b_1[13];
        double ws[7];
        static double dtemp_0[26];
        static double dtemp_1[91];
        int npnts, ncols, rnk, n, i;
        double tng[2], nrm[2];
        double Q_rotation[2][2];
        double u_max, w, f_t, f_tt, ell;
        double dtemp_2, dtemp_3, dtemp_4, dtemp_5, dtemp_6;
        int j;
        double dtemp_7;
 	/*npnts is the number of points being fit */
        npnts = ngbpts_dim1;
        
        /*Calculate the tangent vector using edge-length weighting */

        tng[0] =  ngbpts[(mid - 1) << 1]
                 -  ngbpts[(mid - 2) << 1]
                 + ( ngbpts[mid << 1] -  ngbpts[(mid - 1) << 1]);
        tng[1] =  ngbpts[1 + ((mid - 1) << 1)]
                 -  ngbpts[1 + ((mid - 2) << 1)]
                 + ( ngbpts[1 + (mid << 1)] -  ngbpts[1 + ((mid - 1) << 1)]);
        dtemp_2 = sqrt(tng[0]*tng[0] +  tng[1]*tng[1]);
        tng[0] =  tng[0]/dtemp_2;
        tng[1] =  tng[1]/dtemp_2;
        nrm[0] = - tng[1];
        nrm[1] = tng[0];

        /*normal vector */
        Q_rotation[0][0] = tng[0];
        Q_rotation[1][0] = tng[1];
        Q_rotation[0][1] = nrm[0];
        Q_rotation[1][1] = nrm[1];

        /*the rotation matrix */
        /*Calculate the coordinates in the newly local uv-axis */

        for(i = 0; i <= npnts - 1; i += 1)
        {
            newpts[i << 1] = (ngbpts[i << 1] -  ngbpts[(mid << 1) - 2])
                            *Q_rotation[0][0] + (ngbpts[1 + (i << 1)]
                            -  ngbpts[(mid << 1) - 1])*Q_rotation[1][0];
            newpts[1 + (i << 1)] = (ngbpts[i << 1] -  ngbpts[(mid << 1) - 2])
                            *Q_rotation[0][1] + (ngbpts[1 + (i << 1)]
                            -  ngbpts[(mid << 1) - 1])*Q_rotation[1][1];
        }

        /* Calculate Vandermonde matrix */
        ncols = deg + 1;
        AS_dim1 = npnts;
        AS_dim2 = ncols;
        for(i = 0; i <= AS_dim1 - 1; i += 1)
        {
            AS[i*AS_dim2] = 1.0;
        }
        for(i = 0; i <= npnts - 1; i += 1)
        {
            dtemp_0[i] = newpts[i << 1];
        }
        for(i = 0; i <= AS_dim1 - 1; i += 1)
        {
            AS[1 + i*AS_dim2] = dtemp_0[i];
        }
        for(i = 3; i <= ncols; i += 1)
        {
            for(j = 0; j <= AS_dim1 - 1; j += 1)
            {
                dtemp_1[j] =  AS[ j*AS_dim2 - 2 + i]*newpts[j << 1];
            }
            for(j = 0; j <= AS_dim1 - 1; j += 1)
            {
                AS[j*AS_dim2 - 1 + i] = dtemp_1[j];
            }
        }

        /* Assign a weight to each columns of AS */
        u_max = -HUGE_VAL;
        b_dim1 = npnts;
        for(i = 0; i <= npnts - 1; i += 1)
        {
            u_max = ct_max(fabs(newpts[i << 1]), u_max);
            b[i] = newpts[1 + (i << 1)];
        }

        for(i = 0; i <= npnts - 1; i += 1)
        {
            dtemp_4 = fabs(newpts[i << 1] + 0.1*u_max);
            dtemp_5 = pow(dtemp_4, deg);
            w = 1.0;
            for(j = 0; j <= AS_dim2 - 1; j += 1)
            {
                AS[i*AS_dim2 + j] =  AS[i*AS_dim2 + j]*w;
            }
            b[i] = b[i]*w;
          }

        /* Rescale Vandermonde matrix */
        for(i = 0; i <= ncols - 1; i += 1)
        {
            dtemp_3 = 0.0;
            for(j = 0; j <= AS_dim1 - 1; j += 1)
            {
                dtemp_3 = dtemp_3 +  AS[j*AS_dim2 + i]*AS[j*AS_dim2 + i];
            }
            dtemp_7 = sqrt(dtemp_3);
            for(j = 0; j <= AS_dim1 - 1; j += 1)
            {
                AS[j*AS_dim2 + i] =  AS[j*AS_dim2 + i]/dtemp_7;
            }
            ws[i] = dtemp_7;
        }
        /* Solve using QR factorization.  */
        /* Afterwards, AS contains R and b contains Q'b */
        rnk = safeqr(AS, &AS_dim1, &AS_dim2, b, &b_dim1, 0.001);
        if ((rnk < ncols) && (rnk < 3))
            rnk = 3;

        for(i = 0; i <= b_dim1 - 1; i += 1)
        {
            b_1[i] = b[i];
        }
        /* Perform backward substitution to compute  */
        /* bs = triu(R(1:ncols,1:ncols))\bs; */
        n = rnk;
        for(i = n; i >= 1; i += -1)
        {
            for(j = 1 + i; j <= n; j += 1)
            {
                b_1[i - 1] =  b_1[i - 1] - AS[(i - 1)*AS_dim2 - 1 + j]
                                   *b_1[j - 1];
            }
            b_1[i - 1] =  b_1[i - 1]/AS[(i - 1)*AS_dim2 - 1 + i];
        }
        /*Calculate normal and curvature */
        f_t =  b_1[1]/ws[1];
        f_tt = 2.0*b_1[2]/ws[2];
        ell = sqrt(1.0 + f_t*f_t);
        dtemp_6 = -f_t;
        normal_out[0] =  Q_rotation[0][0]*(dtemp_6 / ell);
        normal_out[0] =  normal_out[0] +  Q_rotation[0][1]*(1.0 / ell);
        normal_out[1] =  Q_rotation[1][0]*(dtemp_6 / ell);
        normal_out[1] =  normal_out[1] +  Q_rotation[1][1]*(1.0 / ell);
        *curvature_out = f_tt/sqrt(ell*ell*ell);
        return;
}

/*******************************************************************
 *  Function for the curvature computation of six points with LSQ   *
 *******************************************************************/
LOCAL double compute_curvature2d1(
	POINT	*p,
        POINT   **pts,
	int	num_pts)

{
        double curvature;
        double normal[MAXD];
        static double **ngbpts;
        int i,j;

        if (ngbpts == NULL)
	    FT_MatrixMemoryAlloc((POINTER*)&ngbpts,num_pts,2,FLOAT);

        /* 2*nrad+1 points fitting */
	for (i = 0; i < num_pts; ++i)
        {
            ngbpts[i][0] = Coords(pts[i])[0];
            ngbpts[i][1] = Coords(pts[i])[1];
        }

        polyfit2d_lhf(normal,&curvature,*ngbpts,num_pts,num_pts/2+1,num_pts/2);
        
        /* seven points cubic*/

        p->_nor[0] = normal[0];
        p->_nor[1] = normal[1];
        p->curvature = curvature;
        return curvature;
}


static void test_normal2d(
        POINT           *p,
        HYPER_SURF_ELEMENT      *hse,
        HYPER_SURF              *hs,
        Front           *front)
{
        double           curvature = 0.0;
        CURVE           *c = Curve_of_hs(hs);
        BOND            *b = Bond_of_hse(hse);
	POINT		*pts[7];
        int             i, dim = c->interface->dim;

        if (b == NULL || c == NULL)
            return;
	if (!FrontGetPointChain(p,pts,5))
	{
	    screen("WARNING open curve too short for curvature!\n");
	    curvature = 0.0;
	}
        curvature = compute_curvature2d1(p,pts,5);
}              

static double test_curvature2d(
        POINT           *p,
        HYPER_SURF_ELEMENT      *hse,
        HYPER_SURF              *hs,
        Front           *front)
{
      return p->curvature;
}


/**********************************************************************
 *      This is a new function for Yanhong Zhao to experiment new     *
 *      algorithm for the calculation of curvature.    3D             *
 **********************************************************************/

LOCAL void obtain_orthbases(
        double *tang,
        double *nrm1st)
{
        double t[3],t2[3];
        double norm,k;
        int i;
        if(fabs(nrm1st[0]) > fabs(nrm1st[1])
                         && fabs(nrm1st[0]) > fabs(nrm1st[2]))
        {
            t[0] = 0.0;
            t[1] = 1.0;
            t[2] = 0.0;
        }
        else
        {
            t[0] = 1.0;
            t[1] = 0.0;
            t[2] = 0.0;
        }
        norm = 0.0;
        k = Dot3d(t,nrm1st);

        for(i = 0;i < 3; i++)
        {
            t[i] = t[i] - k*nrm1st[i];
        }
        norm = Mag3d(t);
        for(i = 0;i < 3; i++)
        {
            tang[i] = t[i]/norm;
        }

        Cross3d(tang,nrm1st,t2);
        for(i = 0; i < 3; i++)
            tang[i+3] = t2[i];
}

#define         MAX_RING_PTS            40

//This function is particularly for normal of first order
EXPORT  void  GetFrontFirstNormal(
        POINT              *p,
        HYPER_SURF_ELEMENT *hse,
        HYPER_SURF         *hs,
        double              *nor,
        Front              *fr)
{
        TRI *tri,**tris;
        POINT *pts_ring1[MAX_RING_PTS];
        POINT *pts_ring2[MAX_RING_PTS];
        POINT *pc,*pp; //current point
        int np1,np2,i,j,nt1,num;
        double nor_f[3],norm=0; ///normal of tri
        const double  *fnor;
        double t1[3],t2[3];
        INTERFACE *intfc = p->hs->interface;
        TRI **ptris,*t,*tris1[20];

        PointArrayRing2(p,hse,hs,&np1,&np2,pts_ring1,pts_ring2);

	tri = Tri_of_hse(p->hse);
        nt1 = set_tri_list_around_point(p,tri,&ptris,intfc);
        for(i = 0; i < 3; i++)
        {
            nor_f[i] = 0;
        }

        for(i = 0; i < nt1; i++)
        {
            t = ptris[i];
            fnor = Tri_normal(t);
            for(j = 0; j < 3; j++)
            {
                nor_f[j] = nor_f[j] + fnor[j];
            }
        }
        for(i = 0; i < 3; i++)
        {
            norm = norm + nor_f[i]*nor_f[i];
        }
        norm = sqrt(norm);
        for(i = 0; i < 3; i++)
        {
            p->_nor[i] = nor_f[i]/norm;
            nor[i] = p->_nor[i];
        }
}

//This function is particularly for normal of second order
EXPORT  void  GetFrontSecondNormal(
        POINT              *p,
        HYPER_SURF_ELEMENT *hse,
        HYPER_SURF         *hs,
        Front              *fr)
{
        TRI **tris;
        POINT *pts_ring1[MAX_RING_PTS];
        POINT *pts_ring2[MAX_RING_PTS];
        POINT *pc,*pt;
        int np1,np2;
        int i,j;

        double tang[6],nrm1st[3],nrm[3],nrm_t[3],nrm_v[50][3];
        double xs_v[50][3],xss_v[50][2],nrm_vs[3];
        double PT[3][3],xss_r[50][2];
        double ws_row[50],w;

        int degree = 2;
        double temp,ncoeffs[]={3,6,10,15,21,28};
        double bs[50][2];
        double x_t[3],coeffs[20];
        double cs[10][2]={0.0},cond;
        double curvature;
        int xs_dim, nrm_coor_dim;
        double xs[60][3],prcurvs[2],maxprdir[3];
        int deg;
        double nrm_coor[60][3];       
 
        pc = p; //point to current point
        PointArrayRing2(pc,hse,hs,&np1,&np2,pts_ring1,pts_ring2);

        for(i = 0; i < 3; i++)
        {
            nrm1st[i] = p->_nor0[i];
        }

        for(i = 0; i <3; i++)
            nrm1st[i] = -nrm1st[i];

        xs_dim = np1 + np2 + 1;
        nrm_coor_dim = np1 + np2 + 1;
        
        xs[0][0] = Coords(pc)[0];
        xs[0][1] = Coords(pc)[1];
        xs[0][2] = Coords(pc)[2];

        nrm_coor[0][0] = pc->_nor0[0];
        nrm_coor[0][1] = pc->_nor0[1];
        nrm_coor[0][2] = pc->_nor0[2];

        for(i = 0; i < np1; i++)
        {
            xs[i+1][0] = Coords(pts_ring1[i])[0];
            xs[i+1][1] = Coords(pts_ring1[i])[1];
            xs[i+1][2] = Coords(pts_ring1[i])[2];
            nrm_coor[i+1][0] = pts_ring1[i]->_nor0[0];
            nrm_coor[i+1][1] = pts_ring1[i]->_nor0[1];
            nrm_coor[i+1][2] = pts_ring1[i]->_nor0[2];
        
        }

        for(i = 0; i < np2; i++)
        {
            xs[i+np1+1][0] = Coords(pts_ring2[i])[0];
            xs[i+np1+1][1] = Coords(pts_ring2[i])[1];
            xs[i+np1+1][2] = Coords(pts_ring2[i])[2];
            nrm_coor[i+np1+1][0] = pts_ring2[i]->_nor0[0];
            nrm_coor[i+np1+1][1] = pts_ring2[i]->_nor0[1];
            nrm_coor[i+np1+1][2] = pts_ring2[i]->_nor0[2];            
        }

        polyfit_lhf3d(nrm,&deg,prcurvs, 
                           maxprdir,*xs,xs_dim,*nrm_coor,nrm_coor_dim,
                           2,1);


        for(i = 0; i < 3; i++)
        {
            pc->_nor[i] = nrm[i];
        }
        pc->curvature = fabs((prcurvs[0]+prcurvs[1])/2);
}

/*************************************************************
 *
 * FUNCTION: polyfit_lhf3d
 *
 *POLYFIT_LHF3D Compute normal, principal curvatures, and principal direction.
 * [NRM,DEG,PRCURVS,MAXPRDIR] = POLYFIT_LHF_SURF_POINT(XS, NRMS_COOR, DEGREE, STRIP)
 * Computes normal NRM, principal curvatures 
 * xs: xs(1,1:3) is the target point and  xs(2:end,1:3)  are neighbor points
 * nrm_coor: nrm_coor(1,1:3) is the normal at target point and  nrm_coor(2:end,1:3) are normals at neighbor points
 *************************************************************/

void  polyfit_lhf3d(
        double nrm_out[3], 
        int *deg_out, 
        double prcurvs_out[2], 
        double maxprdir_out[3], 
        double *xs, 
        int xs_dim1, 
        double *nrm_coor, 
        int nrm_coor_dim1, 
        int degree, 
        int strip)
{
        static double us[3075];
        int us_dim1;
        static double bs[1025];
        int bs_dim1;
        static double nrms_v[3075];
        int nrms_v_dim1;
        static double ws_row[1025];
        int ws_row_dim1;
        double cs[6];
        int cs_dim1;
        static double dtemp_0[3075];
        static double dtemp_1[3075];
        int dtemp_1_dim1;
        static double dtemp_2[3075];
        int strip_1, nverts, n, itemp_0, itemp_1;
        double absnrm[3], t1_1[3], t2[3], nrm_l[3], maxprdir_l[3];
        int t1[3];
        double grad[2];
        double P[3][3];
        double H[2][2];
        static int maxprdir_1[3] = 
                 {
                   0, 0, 0
                 };
        static int t1_2[3] = 
                 {
                   1, 0, 0
                 };
        static int t1_3[3] = 
                 {
                   0, 1, 0
                 };
        double dtemp_3[3];
        double dtemp_4, dtemp_5, dtemp_6, dtemp_7, dtemp_8;
        double dtemp_9, dtemp_10, dtemp_11, dtemp_12, dtemp_13;
        double dtemp_14, dtemp_15, dtemp_16, dtemp_17, dtemp_18;
        int itemp_2, i, j;
        double dtemp_19;

        strip_1 = strip;
        /* First, compute the rotation matrix */

        for (i = 0; i <= 2; i += 1) 
        {
            absnrm[i] = fabs(nrm_coor[i]);
            nrm_out[i] = nrm_coor[i];
        }
        if (( absnrm[0] >  absnrm[1])&&( absnrm[0] >  absnrm[2])) 
        {
            for (i = 0; i <= 2; i += 1) 
            {
                t1[i] = t1_3[i];
            }
        }
        else 
        {
            for (i = 0; i <= 2; i += 1) 
            {
                t1[i] = t1_2[i];
            }
        }
        dtemp_4 = 0.0;
        for (i = 0; i <= 2; i += 1) 
        {
            dtemp_4 = dtemp_4 + (((double) t1[i]))*nrm_out[i];
        }
        dtemp_5 = 0.0;
        for (i = 0; i <= 2; i += 1) 
        {
            dtemp_18 = (((double) t1[i])) - dtemp_4*nrm_out[i];
            dtemp_5 = dtemp_5 + dtemp_18*dtemp_18;
            t1_1[i] = dtemp_18;
        }
        dtemp_6 = sqrt(dtemp_5);
        for (i = 0; i <= 2; i += 1) 
        {
            t1_1[i] = t1_1[i]/dtemp_6;
        }
        /* CROSS: Eficient routine for computing cross product */

        t2[0] = nrm_out[1]*t1_1[2] - nrm_out[2]*t1_1[1];
        t2[1] = nrm_out[2]*t1_1[0] - nrm_out[0]*t1_1[2];
        t2[2] = nrm_out[0]*t1_1[1] - nrm_out[1]*t1_1[0];
        nverts = xs_dim1;
        /* Evaluate local coordinate system */

        us_dim1 = nverts - strip;
        if (nverts < strip) 
        {
            us_dim1 = 0;
        }
        memset(us,0,((nverts - strip) << 1)*8);
        ct_set_max(bs_dim1,nverts - strip,0);
        memset(bs,0,(nverts - strip)*8);
        us[0] = 0.0;
        us[1] = 0.0;
        /* compute weight based on angles of normals */
        nrms_v_dim1 = ct_max(nverts - 1, 0);
        for(i = 0; i <= nverts - 2; i += 1) 
        {
            dtemp_7 = 0.0;
            dtemp_8 = 0.0;
            dtemp_9 = 0.0;
            for(j = 0; j <= 2; j += 1) 
            {
                dtemp_19 = xs[(1 + i)*3 + j] - xs[j];
                dtemp_7 = dtemp_7 + dtemp_19*t1_1[j];
                dtemp_8 = dtemp_8 + dtemp_19*t2[j];
                dtemp_9 = dtemp_9 + dtemp_19*nrm_out[j];
            }
            us[2 - (strip << 1) + (i << 1)] = dtemp_7;
            us[3 - (strip << 1) + (i << 1)] = dtemp_8;
            bs[1 - strip + i] = dtemp_9;
            for(j = 0; j <= 3 - 1; j +=1 ) 
            {
                nrms_v[i*3 + j] = nrm_coor[(1 + i)*3 + j];
            }
        }
        if(!strip_1) 
        {
            strip_1 = 0;
            for (i = 0; i <= nrms_v_dim1 - 1; i += 1) 
            {
                dtemp_0[i] =  nrms_v[i*3]*nrm_out[0] 
                               + nrms_v[1 + i*3]*nrm_out[1];
            }
            itemp_0 = ct_max(nrms_v_dim1, 0);
            if (nrms_v_dim1 == 1) 
            {
                dtemp_1_dim1 = itemp_0;
                for (i=0; i<=itemp_0 - 1; i+=1) 
                {
                    dtemp_1[i] =  dtemp_0[0] + nrms_v[2 + i*3]*nrm_out[2];
                }
            }
            else 
            {
                dtemp_1_dim1 = nrms_v_dim1;
                if(itemp_0 == 1) 
                {
                    for (i = 0; i <= nrms_v_dim1 - 1; i += 1)
                    {
                        dtemp_1[i] = dtemp_0[i] + nrms_v[2]*nrm_out[2];
                    }
                }
                else 
                {
                    for (i = 0; i <= dtemp_1_dim1 - 1; i += 1) 
                    {
                        dtemp_1[i] = dtemp_0[i] + nrms_v[2 + i*3]*nrm_out[2];
                    }
                }
            }
            itemp_1 = 0;
            if(dtemp_1_dim1) 
            {
                itemp_1 = dtemp_1_dim1;
            }
            ws_row_dim1 = itemp_1 + 1;
            ws_row[0] = 1.0;
            for(i = 0; i <= dtemp_1_dim1 - 1; i += 1) 
            {
                dtemp_2[i] = ct_max(dtemp_1[i], 0.0);
                ws_row[1 + i] = dtemp_2[i];
            }
        }
        else 
        {
            ws_row_dim1 = nrms_v_dim1;
            for(i = 0; i <= nrms_v_dim1 - 1; i += 1) 
            {
                dtemp_17 = nrms_v[i*3]*nrm_out[0] + nrms_v[1 + i*3]*nrm_out[1] 
                             +  nrms_v[2 + i*3]*nrm_out[2];
                dtemp_10 = ct_max(dtemp_17, 0.0);
                ws_row[i] = dtemp_10;
            }
        }
        /* Compute the coefficients */

        *deg_out = eval_vander_bivar(us, us_dim1, bs, &bs_dim1,
                            degree, ws_row, ws_row_dim1, strip_1);
        /* Convert coefficients into normals and curvatures */

        if((*deg_out) == 1) 
        {
            n = 3 - strip_1;
        }
        else 
        {
            n = 6 - strip_1;
        }
        itemp_2 = ct_max(strip_1 - 1 + n, 0);
        if(bs_dim1 > 1) 
        {
            cs_dim1 = itemp_2;
        }
        else 
        {
            cs_dim1 = 1;
        }
        for(i = 0; i <= cs_dim1 - 1; i += 1) 
        {
            cs[i] = bs[1 - strip_1 + i];
        }
        grad[0] = cs[0];
        grad[1] = cs[1];
        dtemp_3[0] = - grad[0];
        dtemp_3[1] = - grad[1];
        dtemp_3[2] = 1.0;
        dtemp_5 = sqrt(1.0 + (grad[0]*grad[0] + grad[1]*grad[1]));
        for(i = 0; i <= 2; i += 1) 
        {
            nrm_l[i] = dtemp_3[i]/dtemp_5;
            P[i][0] = t1_1[i];
            P[i][1] = t2[i];
            P[i][2] = nrm_out[i];
        }
        /* nrm = P*nrm_l; */

        dtemp_11 = 0.0;
        dtemp_12 = 0.0;
        dtemp_13 = 0.0;
        for(i = 0; i <= 2; i += 1) 
        {
            dtemp_11 = dtemp_11 + P[0][i]*nrm_l[i];
            dtemp_12 = dtemp_12 + P[1][i]*nrm_l[i];
            dtemp_13 = dtemp_13 + P[2][i]*nrm_l[i];
        }
        nrm_out[0] = dtemp_11;
        nrm_out[1] = dtemp_12;
        nrm_out[2] = dtemp_13;
        if ((cs_dim1 >= 5)&&((*deg_out) > 1)) 
        {
            H[0][0] = cs[2];
            H[0][1] = cs[3];
            H[1][0] = cs[3];
            H[1][1] = cs[4];
            eval_curvature_lhf_surf(prcurvs_out, maxprdir_l, grad, H);
            /* maxprdir = P * maxprdir_l; */

            dtemp_14 = 0.0;
            dtemp_15 = 0.0;
            dtemp_16 = 0.0;
            for(i = 0; i <= 2; i += 1) 
            {
                dtemp_14 = dtemp_14 + P[0][i]*maxprdir_l[i];
                dtemp_15 = dtemp_15 + P[1][i]*maxprdir_l[i];
                dtemp_16 = dtemp_16 + P[2][i]*maxprdir_l[i];
            }
            maxprdir_out[0] = dtemp_14;
            maxprdir_out[1] = dtemp_15;
            maxprdir_out[2] = dtemp_16;
        }
        else 
        {
            prcurvs_out[0] = 0.0;
            prcurvs_out[1] = 0.0;
            for(i = 0; i <= 2; i += 1) 
            {
                maxprdir_out[i] = ((double) maxprdir_1[i]);
            }
        }
        return;
}

/*************************************************************
 *
 * FUNCTION: eval_curvature_lhf_surf
 *
 *EVAL_CURVATURE_LHF_SURF Compute principal curvature, principal direction 
 *and pseudo-inverse.
 * [CURVS,DIR,JINV] = EVAL_CURVATURE_LHF_SURF(GRAD,H) Computes principal 
 * curvature in 2x1 CURVS, principal direction of maximum curvature in 3x2 
 * DIR, and pseudo-inverse of J in 2x3 JINV.  Input arguments are the
 * gradient of the height function in 2x1 GRAD, and the Hessian of the
 * height function in 2x2 H with a local coordinate frame.
 * See also EVAL_CURVATURE_LHFINV_SURF, EVAL_CURVATURE_PARA_SURF
 *************************************************************/

void  eval_curvature_lhf_surf(
        double curvs_out[2], 
        double dir_out[3], 
        double grad[2], 
        double H[2][2])
{
        double grad_sqnorm, grad_norm, ell, ell2, ell3;
        double c, s, kH2, tmp, dtemp_0;
        double v[2], W1[2];
        double W[2][2];
        double U[3][2];
        double d1[2];
        int i;

        grad_sqnorm = grad[0]*grad[0] + grad[1]*grad[1];
        grad_norm = sqrt(grad_sqnorm);
        /* Compute key parameters */

        ell = sqrt(1.0 + grad_sqnorm);
        ell2 = 1.0 + grad_sqnorm;
        ell3 = ell*(1.0 + grad_sqnorm);
        if(grad_norm == 0.0)
        {
            c = 1.0;
            s = 0.0;
        }
        else 
        {
            c =  grad[0]/grad_norm;
            s =  grad[1]/grad_norm;
        }
        /* Compute mean curvature and Gaussian curvature */
        /* kH2 = (H(1,1)+H(2,2))/ell - grad*H*grad'/ell3; */
        /* kG =  (H(1,1)*H(2,2)-H(1,2)^2)/ell2^2; */
        /* Solve quadratic equation to compute principal curvatures */

        v[0] = c*H[0][0] + s*H[0][1];
        v[1] = c*H[0][1] + s*H[1][1];
        W1[0] = (v[0]*c + v[1]*s)/ell3;
        W1[1] = (v[0]*(-s) + v[1]*c)/ell2;
        W[0][0] = W1[0];
        W[0][1] = W1[1];
        W[1][0] = W1[1];
        W[1][1] = ((c*H[0][1] - s*H[0][0])*(-s) 
                       + (c*H[1][1] - s*H[0][1])*c)/ell;
        /* Lambda = eig(W); */

        kH2 =  W[0][0] +  W[1][1];
        tmp = sqrt((W[0][0] - W[1][1])*(W[0][0] - W[1][1]) 
                        + 4.0*W[0][1]*W[0][1]);
        if(kH2 > 0.0) 
        {
            curvs_out[0] = (kH2 + tmp)/2.0;
            curvs_out[1] = (kH2 - tmp)/2.0;
        }
        else 
        {
            curvs_out[0] = (kH2 - tmp)/2.0;
            curvs_out[1] = (kH2 + tmp)/2.0;
        }
        /* Compute principal directions, first with basis of left  */
        /* singular vectors of Jacobian */
        /* Compute principal directions in 3-D space */

        U[0][0] = c/ell;
        U[0][1] = -s;
        U[1][0] = s/ell;
        U[1][1] = c;
        U[2][0] = grad_norm/ell;
        U[2][1] = 0.0;
        if(curvs_out[0] ==  curvs_out[1]) 
        {
            for(i = 0; i <= 2; i += 1) 
            {
                dir_out[i] = U[i][0];
            }
        }
        else 
        {
            if(fabs(W[0][0] - curvs_out[1]) > fabs(W[0][0] - curvs_out[0])) 
            {
                d1[0] = W[0][0] - curvs_out[1];
                d1[1] = W[0][1];
            }
            else
            {
                d1[0] = - W[0][1];
                d1[1] = W[0][0] -  curvs_out[0];
            }
            dtemp_0 = sqrt(d1[0]*d1[0] + d1[1]*d1[1]);
            d1[0] = d1[0]/dtemp_0;
            d1[1] = d1[1]/dtemp_0;
            dir_out[0] = U[0][0]*d1[0] + U[0][1]*d1[1];
            dir_out[1] = U[1][0]*d1[0] + U[1][1]*d1[1];
            dir_out[2] = U[2][0]*d1[0] + U[2][1]*d1[1];
        }
        return;
}

/*************************************************************
 *
 * FUNCTION: eval_vander_bivar
 *
 *EVAL_VANDER_BIVAR Evaluate generalized Vandermonde matrix.
 * [BS,DEGREE] = EVAL_VANDER_BIVAR(US,BS,DEGREE,WS,STRIP) Evaluates
 * generalized Vandermonde matrix V, and solve V\BS. It supports up to
 * degree 6.
 * See also EVAL_VANDER_UNIVAR
 *************************************************************/

int eval_vander_bivar(
        double *us, 
        int us_dim1, 
        double *bs, 
        int *bs_dim1, 
        int degree, 
        double *ws, 
        int ws_dim1, 
        int strip)
{
        static double V[28672];
        int V_dim1, V_dim2;
        static double ws1[1024];
        static double ts[1024];
        int ncols;
        static double v[1024];
        int degree_1_out, npnts, jj, rnk, ncols_sub;
        double t2, t, vnrm2, s, w;
        int ncoeff, itemp_0, itemp_1, i, j;
        static int ncoeffs[6] = 
                   {
                      3, 6, 10, 15, 21, 28
                   };
        static double weights[10] = 
                      {
                         1., 1.0, 1.0, 2.0, 1.0, 2.0, 6.0, 2.0, 2.0, 6.0
                      };
        double dtemp_0, dtemp_1, dtemp_2, dtemp_3, dtemp_4;
        int i3;

        degree_1_out = degree;
        /* Determine degree of fitting */

        npnts = us_dim1;
        /* Determine degree of polynomial */

        ncols =  ncoeffs[degree - 1] - strip;
        while((ncoeffs[degree_1_out - 1] - strip > npnts)&&(degree_1_out > 1))
        {
            degree_1_out = degree_1_out - 1;
            ncols =  ncoeffs[degree_1_out - 1] - strip;
        }
        /* Construct matrix */

        V_dim1 = ct_max(npnts, 0);
        ct_set_max(V_dim2, ncols, 0);
        memset(V, 0, npnts * ncols * 8);
        /* Allocate the maximum size */

        for (i = 0; i <= npnts - 1; i += 1) 
        {
            V[i*V_dim2] = 1.0;
            jj = 2 - strip;
            V[i*V_dim2 - 1 + jj] = us[i << 1];
            jj = 1 + jj;
            V[i*V_dim2 - 1 + jj] = us[1 + (i << 1)];
            for(j = 2; j <= degree_1_out; j += 1) 
            {
                jj = 1 + jj;
                V[i*V_dim2 - 1 + jj] = V[i*V_dim2 - 1 + jj - j]*us[i << 1];
                for(i3 = 0; i3 <= j - 1; i3 += 1) 
                {
                    jj = 1 + jj;
                    V[i*V_dim2 - 1 + jj] = V[i*V_dim2 - 2 + jj - j]
                                                *us[1 + (i << 1)];
                }
            }
        }
        /* Scale rows to assign different weights to different points */

        if(degree_1_out > 1) 
        {
        /* Scale weights to be inversely proportional to distance */

            dtemp_1 = 0.0;
            for(i = 0; i <= us_dim1 - 1; i += 1) 
            {
                dtemp_4 = us[i << 1]*us[i << 1] 
                           + us[1 + (i << 1)]*us[1 + (i << 1)];
                dtemp_1 = dtemp_1 + dtemp_4;
                ws1[i] = dtemp_4;
            }
            for(i = 0; i <= us_dim1 - 1; i += 1) 
            {
                ws1[i] = ws1[i] + 0.01*(dtemp_1/(((double)npnts)));
            }
            if(degree_1_out < 4) 
            {
                for(i = 0; i <= npnts - 1; i += 1) 
                {
                    if(ws1[i] != 0.0) 
                    {
                        ws1[i] = ws[i]/sqrt(ws1[i]);
                    }
                    else 
                    {
                        ws1[i] = ws[i];
                    }
                }
            }
            else 
            {
                for(i = 0; i <= npnts - 1; i += 1) 
                {
                    if(ws1[i] != 0.0) 
                    {
                        ws1[i] = ws[i]/ws1[i];
                    }
                    else 
                    {
                        ws1[i] = ws[i];
                    }
                }
            }
        }
        else 
        {
            for(i = 0; i <= ws_dim1 - 1; i += 1) 
            {
                ws1[i] = ws[i];
            }
        }
        
        for(i = 0; i <= npnts - 1; i += 1) 
        {
            for(j = 0; j <= V_dim2 - 1; j += 1) 
            {
                V[i*V_dim2 + j] = V[i*V_dim2 + j]*ws1[i];
            }
            bs[i] = bs[i]*ws1[i];
        }
        /* Scale columns to reduce condition number */

        memset(ts, 0, ncols*8);
        for (i = 0; i <= ncols - 1; i += 1) 
        {
            /*NORM2_VEC Computes the 2-norm of a vector. */
            /* NORM2_VEC(V) Computes the 2-norm of column vector V, */
            /*       taking care not to cause unnecessary overflow. */
            /* See also SQNORM2_VEC */

            w = 0.0;
            for(j = 0; j <= V_dim1 - 1; j += 1) 
            {
                dtemp_2 = fabs(V[j*V_dim2 + i]);
                w = ct_max(dtemp_2, w);
                v[j] = V[j*V_dim2 + i];
            }
            s = 0.0;
            if(w == 0.0) 
            {
                /* W can be zero for max(0,nan,...) */
                /* adding all three entries together will make sure */
                /* NaN will not disappear. */

                for(j = 0; j <= V_dim1 - 1; j += 1) 
                {
                    s = s + v[j];
                }
            }
            else 
            {
                for(j = 0; j <= V_dim1 - 1; j += 1) 
                {
                    dtemp_0 = v[j]/w;
                    s = s + dtemp_0*dtemp_0;
                }
                s = w*sqrt(s);
            }
            dtemp_3 = s;
            if (dtemp_3 == 0.0) 
            {
                dtemp_3 = 1.0;
            }
            else 
            {
                for(j = 0; j <= V_dim1 - 1; j += 1) 
                {
                    V[j*V_dim2 + i] = V[j*V_dim2 + i]/dtemp_3;
                }
            }
            ts[i] = dtemp_3;
        }  
        /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
        /* Perform Householder QR factorization to compute */
        /* [Q, V(1:n,1:n)] = qr(V(1:m,1:n),0); bs = Q'*bs; rnk=n; */

        rnk = ncols;
        for(i = 0; i <= ncols - 1; i += 1) 
        {
            itemp_1 = npnts - i;
            if(npnts < i) 
            {
                itemp_1 = 0;
            }
            for(j = 0; j <= npnts - 1 - i; j += 1) 
            {
                v[j] = V[(i + j)*V_dim2 + i];
            }
            /* We don't need to worry about overflow
                                  , since V has been rescaled. */

            t2 = 0.0;
            for(j = 0; j <= itemp_1 - 1; j += 1) 
            {
                t2 = t2 + v[j]*v[j];
            }
            t = sqrt(t2);
            if(v[0] >= 0.0) 
            {
                vnrm2 = sqrt(2.0*(t2 + v[0]*t));
                v[0] = v[0] + t;
            }
            else 
            {
                vnrm2 = sqrt(2.0*(t2 - v[0]*t));
                v[0] = v[0] - t;
            }
            if(vnrm2 > 0.0) 
            {
                for(j = 0; j <= itemp_1 - 1; j += 1) 
                {
                    v[j] = v[j]/vnrm2;
                }
            }
            /* Optimized version for */
            /* V(k:npnts,k:ncols) = V(k:npnts,k:ncols) 
                                 - 2*v*(v'*V(k:npnts,k:ncols)); */

            for(j = 1 + i; j <= ncols; j += 1) 
            {
                t2 = 0.0;
                for(i3 = 0; i3 <= itemp_1 - 1; i3 += 1) 
                {
                    t2 = t2 + v[i3]*V[(i + i3)*V_dim2 - 1 + j];
                }
                t2 = t2 + t2;
                for(i3 = 0; i3 <= itemp_1 - 1; i3 += 1) 
                {
                    V[(i + i3)*V_dim2 - 1 + j] = V[(i + i3)*V_dim2 - 1 + j] 
                                                       - t2*v[i3];
                }
            }
            /* Estimate rank of matrix */

            if(( fabs(V[i*V_dim2 + i]) < 0.001)&&(rnk == ncols)) 
            {
                rnk = i;
                if(rnk > 3) 
                {
                    break;
                }
            }
            /* Optimized version for */
            /* bs(k:npnts,:) = bs(k:npnts,:) - 2*v*(v'*bs(k:npnts,:)); */

            t2 = 0.0;
            for(j = 0; j <= itemp_1 - 1; j += 1) 
            {
                t2 = t2 +  v[j]*bs[i + j];
            }
            t2 = t2 + t2;
            for(j = 0; j <= itemp_1 - 1; j += 1) 
            {
                bs[i + j] =  bs[i + j] - t2*v[j];
            }
        }
        /*%%%%%%%%%%%%%%%%%%%%% */

        ncols_sub = ncols;
        while((rnk < ncols_sub)&&(degree_1_out > 1)) 
        {
            degree_1_out = degree_1_out - 1;
            ncols_sub = ncoeffs[degree_1_out - 1] - strip;
        }
        /*%%%%%%%%%%%%%%%%%%%%% */
        /* Perform backward substitution to compute */
        /* bs = triu(R(1:ncols,1:ncols))\bs; */

        for(i = ncols_sub; i >= 1; i += -1) 
        {
            for(j = 1 + i; j <= ncols_sub; j += 1) 
            {
                bs[i - 1] = bs[i - 1] - V[(i - 1)*V_dim2 - 1 + j]*bs[j - 1];
            }
            if(V[(i - 1)*V_dim2 - 1 + i] == 0.0) 
            {
                /* Singular matrix */

                /*__CATALYTIC_warning("Matrix is singular.", 19);*/
	        printf("Matrix is singular.\n");
                /*#ok<WNTAG> */

                bs[i - 1] = 0.0;
            }
            else 
            {
                bs[i - 1] =  bs[i - 1] /  V[ (i - 1) * V_dim2 - 1 + i];
            }
        }
        /*%%%%%%%%%%%%%%%%%%%%% */
        /* Force weights to be double precision. */
        
        ncoeff = 6 - strip;
        itemp_0 = ct_min(ncols_sub, ncoeff);
        for(i = 0; i <= itemp_0 - 1; i += 1) 
        {
            bs[i] =  bs[i]*weights[strip + i]/ts[i];
        }
        return degree_1_out;
}

static void test_normal3d(
        POINT           *p,
        HYPER_SURF_ELEMENT      *hse,
        HYPER_SURF              *hs,
        Front           *front)
{

        GetFrontSecondNormal(p,hse,hs,front);
}

static double test_curvature3d(
        POINT           *p,
        HYPER_SURF_ELEMENT      *hse,
        HYPER_SURF              *hs,
        Front           *front)
{
        return p->curvature;
}

static void print_hyper_error(
	Front *front,
	char *out_name,
	TEST_CIRCLE_PARAMS circle_params)
{
	POINT              *p;
        HYPER_SURF_ELEMENT *hse;
        HYPER_SURF         *hs;
	INTERFACE	   *intfc = front->interf;
	double	 	   *center = circle_params.center;
	double	 	   r, *R =  circle_params.R;
	int		   i,dim = circle_params.dim;
	double		   nor[MAXD],exact_nor[MAXD];
	double		   kappa, exact_kappa;
	double		   N1_curvature,N2_curvature,Ni_curvature;
	double		   N1_norm[MAXD],N2_norm[MAXD],Ni_norm[MAXD];
        double             iN1_norm[MAXD],iN2_norm[MAXD],iNi_norm[MAXD];
        double             norm_n1, norm_n2, norm_ni;
	double		   N, ra, rb, a, b, kk;
        double             cosv,sinv,cosu,sinu,cos2v,cos2u;
	double             x1,x2,y1,y2,z1,z2;
        char		   fname[100];
	FILE		   *outfile;
        double             df,dd;
        
        double error,err[3],errnormal;
        double mean_k,x,y,z,norm,acosv,t;
        

	sprintf(fname,"%s-error",out_name);
	outfile = fopen(fname,"w");

	N1_curvature = N2_curvature = Ni_curvature = 0.0;
	for (i = 0; i < dim; ++i)
	{
            N1_norm[i] = N2_norm[i] = Ni_norm[i] = 0.0;
            iN1_norm[i] = iN2_norm[i] = iNi_norm[i] = 0.0;
        }    
	N = 0.0;
        norm_n1 = 0;
        norm_n2 = 0;
        norm_ni = 0;

	fprintf(outfile,"Computed Curvature      Exact Curvature\n");
	next_point(intfc,NULL,NULL,NULL);
	
        while (next_point(intfc,&p,&hse,&hs))
	{
	    if (wave_type(hs) < FIRST_PHYSICS_WAVE_TYPE) continue;
	    r = 0.0;
	    for (i = 0; i < dim; ++i)
		r += sqr(Coords(p)[i] - center[i]);
	    r = sqrt(r);
	    
            switch (method)
            {
               case TEST_HYPER_SPHERE:
                  for (i = 0; i < dim; ++i)
                     exact_nor[i] = (Coords(p)[i] - center[i])/r;
                     exact_kappa = 1.0/R[0];
                  break;
               case TEST_ELLIPSE:
                  if(dim == 2)
                  {
                      b = R[1];
                      kk = 4.0*(Coords(p)[1] - center[1])
                                    /(Coords(p)[0] - center[0]);
                      exact_nor[0] = 1.0/sqrt(1 + kk*kk);
                      if(Coords(p)[0] == center[0])
                      {
                         exact_nor[1] = 1.0;
                      }
                      else
                      {
                         exact_nor[1] = fabs(kk/sqrt(1 + kk*kk));
                      } 
                      ra = Coords(p)[0] - center[0];
                      rb = Coords(p)[1] - center[1];
                      exact_kappa = 1.0/(4*sqr(sqr(b))*(sqr(ra)
                                      /(16*sqr(sqr(b))) + sqr(rb)/sqr(sqr(b)))
                                         *sqrt(sqr(ra)/(16*sqr(sqr(b))) 
                                          + sqr(rb)/sqr(sqr(b))));
                  }
                  else
                  {   
                      x1 = 1; x2 = 0;
                      y1 = (x1 - Coords(p)[0])*(Coords(p)[1] - 0.5)/sqr(R[1])
                                 *sqr(R[0])/(Coords(p)[0] - 0.5) + Coords(p)[1];
                      y2 = (x2 - Coords(p)[0])*(Coords(p)[1] - 0.5)/sqr(R[1])
                                 *sqr(R[0])/(Coords(p)[0] - 0.5) + Coords(p)[1];
                      z1 = (x1 - Coords(p)[0])*(Coords(p)[2] - 0.5)/sqr(R[2])
                                 *sqr(R[0])/(Coords(p)[0] - 0.5) + Coords(p)[2];
                      z2 = (x2 - Coords(p)[0])*(Coords(p)[2] - 0.5)/sqr(R[2])
                                 *sqr(R[0])/(Coords(p)[0] - 0.5) + Coords(p)[2];
                      exact_nor[0] = (x1 - x2)/sqrt(sqr(x1 - x2) 
                                        + sqr(y1 - y2) + sqr(z1 - z2));
                      exact_nor[1] = (y1 - y2)/sqrt(sqr(x1 - x2) 
                                        + sqr(y1 - y2) + sqr(z1 - z2));
                      exact_nor[2] = (z1 - z2)/sqrt(sqr(x1 - x2) 
                                        + sqr(y1 - y2) + sqr(z1 - z2));
                 
                      cosv = (Coords(p)[2] - center[2])/R[2];
                      sinv = sqrt(1 - sqr(cosv));
                      cosu = (Coords(p)[0] - center[0])/(R[0]*sinv);
                      sinu = (Coords(p)[1] - center[1])/(R[1]*sinv);
                      cos2v = 2*sqr(cosv) - 1;
                      cos2u = 2*sqr(cosu) - 1;
          
                      exact_kappa = R[0]*R[1]*R[2]*(3*(sqr(R[0]) + sqr(R[1])) 
                                        + 2*sqr(R[2]) + (sqr(R[0]) + sqr(R[1])
                                        - 2*sqr(R[2]))*cos2v - 2*(sqr(R[0]) 
                                          - sqr(R[1]))*cos2u*sqr(sinv))
                                  /(8*(sqr(R[0])*sqr(R[1])*sqr(cosv) 
                                      + sqr(R[2])*(sqr(R[1])*sqr(cosu)
                                      + sqr(R[0])*sqr(sinu))*sqr(sinv))
                                    *sqrt(sqr(R[0])*sqr(R[1])*sqr(cosv) 
                                      + sqr(R[2])*(sqr(R[1])*sqr(cosu)                                                + sqr(R[0])*sqr(sinu))*sqr(sinv)));
                   }  
                   break;
               case TEST_TORUS:
                  ////////////calculate the exact normal for torus/////
                  x = Coords(p)[0];
                  y = Coords(p)[1];
                  z = Coords(p)[2];
                  acosv = sqrt(R[1]*R[1] - z*z)/R[0];
                  r = sqrt(x*x + y*y);
                  if(r < R[0]) acosv = -acosv;
                  t = 1.0 - R[0]/r;
                  err[0] = t*x;
                  err[1] = t*y;
                  err[2] = z;
                  norm = 0;
                  for(i = 0; i < 3; i++)
                  {
                     norm = norm + err[i]*err[i];
                  }
                  for (i = 0; i < dim; ++i)
                     exact_nor[i] = -err[i]/sqrt(norm);
                  exact_kappa = (1.0 + 2*acosv)/(2.0*R[1]*(1.0+acosv));
                  break;
            }
	    GetFrontNormal(p,hse,hs,nor,front);
	    GetFrontCurvature(p,hse,hs,&kappa,front);
             
            if(kappa < 0.0) 
                kappa = fabs(kappa);
                ///for test ellipse there is a point of negative kappa

            N1_curvature += fabs(kappa - exact_kappa);
            N2_curvature += sqr(kappa - exact_kappa);

            if ((Ni_curvature < fabs(kappa - exact_kappa))&&(kappa!=0))
                Ni_curvature = fabs(kappa - exact_kappa);

            fprintf(outfile,"%18.14f   %18.14f\n",kappa,exact_kappa);
              
            for (i = 0; i < dim; ++i)
	    {
                nor[i] = p->_nor[i];//new method
                iN1_norm[i] = fabs(fabs(nor[i]) - fabs(exact_nor[i]));
                iN2_norm[i] = sqr(fabs(nor[i]) - fabs(exact_nor[i]));
                
                N1_norm[i] += fabs(fabs(nor[i]) - fabs(exact_nor[i]));
		N2_norm[i] += sqr(fabs(nor[i]) - fabs(exact_nor[i]));
	    	if ((Ni_norm[i] < fabs(fabs(nor[i]) - fabs(exact_nor[i]))))
                 Ni_norm[i] = fabs(fabs(nor[i]) - fabs(exact_nor[i]));
            }

            if(dim == 3)
            {
                norm_n1 += iN1_norm[0] + iN1_norm[1] + iN1_norm[2];
                norm_n2 += iN2_norm[0] + iN2_norm[1] + iN2_norm[2]; 
                norm_ni = Ni_norm[0] + Ni_norm[1] + Ni_norm[2];
            }
            if(dim == 2)
            {
                norm_n1 += iN1_norm[0] + iN1_norm[1];
                norm_n2 += iN2_norm[0] + iN2_norm[1];
                norm_ni = Ni_norm[0] + Ni_norm[1];
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
        
        norm_n1 = norm_n1/N;
        norm_n2 = sqrt(norm_n2/N);
        
        fprintf(outfile,"N:%f\n",N);
	
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
        
        fprintf(outfile,"Error norms for normal:\n");
        fprintf(outfile,"L-1 norm:   %24.18g\n",norm_n1);
        fprintf(outfile,"L-2 norm:   %24.18g\n",norm_n2);
        fprintf(outfile,"L-I norm:   %24.18g\n",norm_ni);
        
        fprintf(outfile,"\n\n"); 
 
	fprintf(outfile,"Error norms for curvature:\n");
	fprintf(outfile,"L-1 norm:   %24.18g\n",N1_curvature);
	fprintf(outfile,"L-2 norm:   %24.18g\n",N2_curvature);
	fprintf(outfile,"L-I norm:   %24.18g\n",Ni_curvature);
}	/* end print_hyper_error */

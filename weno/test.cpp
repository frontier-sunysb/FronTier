/*
 * test.c
 *
 *  Created on: Jun 12, 2011
 *      Author: yli
 */

#include "weno.h"

static void readWenoParams(PARAMS*,char*);
static void setSolutionType(INIT_TYPE);

char *in_name,*out_name;

static void (*exact_soln)(double,double,double*,double*,int);
static void wave_func(double,double,double*,double*,int);
static void hamp_func(double,double,double*,double*,int);
static void square_func(double,double,double*,double*,int);

int main(int argc, char **argv) 
{
	F_BASIC_DATA f_basic; /* Only some basic features of f_basic is used */
        Front front;          /* Only some basic feature of front is used */
        double *x,*u_old,*u_new,*u_sol; /* mesh points and solutions */
        double *x_movie,*y_movie;
        double dx,dt;
        char movie_caption[100];
        char time_label[100];
        char xg_name[200];
        char gd_name[200];
        FILE *xg_file;
        double xmin,xmax,umin,umax,height;
	double L,U;
        double CFL;     /* CFL safety factor */
        int buffer_size = 3;
	int i,imin,imax,mesh_size,expanded_mesh_size;
	int boundary_type[2];
	PARAMS params;
	double wave_speed;

	FT_Init(argc,argv,&f_basic);
        in_name         = f_basic.in_name;
        out_name        = f_basic.out_name;

	FT_ReadSpaceDomain(in_name,&f_basic);
	L = f_basic.L[0];
        U = f_basic.U[0];
        mesh_size = f_basic.gmax[0];
	boundary_type[0] = f_basic.boundary[0][0];
        boundary_type[1] = f_basic.boundary[0][1];
        imin = buffer_size;
        imax = imin + mesh_size;

	FT_ReadTimeControl(in_name,&front);
        CFL = Time_step_factor(&front);
        FT_ResetTime(&front);
	readWenoParams(&params,in_name);
	setEquationType(params);
	setSolutionType(params.init_type);
	wave_speed = params.a;

	/* Allocating memory for arrays */
        expanded_mesh_size = mesh_size + 2*buffer_size;
        FT_VectorMemoryAlloc((POINTER*)&u_old,expanded_mesh_size,
                                        sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&u_new,expanded_mesh_size,
                                        sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&u_sol,expanded_mesh_size,
                                        sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&x,expanded_mesh_size,
                                        sizeof(double));
	dx = (U - L)/(mesh_size - 1);
        FT_TimeControlFilter(&front);

	/* Set the initial condition */
        for (i = 0; i < expanded_mesh_size; i++)
	{
            x[i] = L + (i - buffer_size)*dx;
	}
        exact_soln(wave_speed,0.0,x,u_old,expanded_mesh_size);
	
	/* Set frame margin for GD movie output */
        xmin = L;       xmax = U;
        umin = HUGE;    umax = -HUGE;
        for (i = 0; i < mesh_size; i++)
        {
            if (umin > u_old[i]) umin = u_old[i];
            if (umax < u_old[i]) umax = u_old[i];
        }
        height = umax - umin;
        umin -= 0.15*height;    umax += 0.15*height;

        /* Open and initialize GD movie output */
        sprintf(movie_caption,"u vs. x");
        sprintf(gd_name,"%s/soln.gif",out_name);
        gd_initplot(gd_name,movie_caption,xmin,xmax,umin,umax,2);

	/* Time loop */
	front.dt = dt = 0.4*dx/wave_speed;
        for (;;)
        {
            /* Advancing numerical solution */
	    printf("dx = %f  dt = %f\n",dx,dt);
	    Weno5(mesh_size,u_old+buffer_size,u_new+buffer_size,dx,front.dt);

            /* Swapping solution storage */
            for (i = imin; i < imax; i++)
            {
                u_old[i] = u_new[i];
            }

            /* Time and step control */
            FT_AddTimeStepToCounter(&front);
            (void) printf("\ntime = %20.14f   step = %5d   ",
                                front.time,front.step);

	    /* Movie frame */
            if (FT_IsMovieFrameTime(&front))
            {
                /* Numerical solution */
                x_movie = x+buffer_size;
                y_movie = u_old+buffer_size;
                gd_plotdata(mesh_size,x_movie,y_movie);

                /* Exact solution */
                exact_soln(wave_speed,front.time,x,u_sol,expanded_mesh_size);
                x_movie = x+buffer_size;
                y_movie = u_sol+buffer_size;
                gd_plotdata(mesh_size,x_movie,y_movie);

                /* Time label */
                sprintf(time_label,"Time = %6.3f",front.time);
                gd_plotframe(time_label);
            }
            /* Output date control */
            if (FT_IsSaveTime(&front))
            {
                /* Numerical solution */
		sprintf(xg_name,"%s/num_sol-%d.xg",out_name,front.ip);
                xg_file = fopen(xg_name,"w");
                fprintf(xg_file,"\"u vs. x\"\n");
                for (i = imin; i < imax; ++i)
                {
                    fprintf(xg_file,"%f  %f\n",x[i],u_old[i]);
                }
                fclose(xg_file);

                /* Exact solution */
                exact_soln(wave_speed,front.time,x,u_sol,expanded_mesh_size);
                sprintf(xg_name,"%s/exc-%d.xg",out_name,front.ip);
                xg_file = fopen(xg_name,"w");
                fprintf(xg_file,"\"u vs. x\"\n");
                for (i = imin; i < imax; ++i)
                {
                    fprintf(xg_file,"%f  %f\n",x[i],u_sol[i]);
                }
                fclose(xg_file);
            }

	    /* Termination control */
            if (FT_TimeLimitReached(&front))
            {
                front.dt = dt;
                FT_TimeControlFilter(&front); /* reduce time step for output */
                (void) printf("next dt = %20.14f\n",front.dt);
                break;
            }

            front.dt = dt;
            FT_TimeControlFilter(&front); /* reduce time step for output */
            (void) printf("next dt = %20.14f\n",front.dt);
        }
        gd_closeplot();
	free(x);
	free(u_old);
	free(u_new);
	free(u_sol);
	clean_up(0);
}

static void readWenoParams(
	PARAMS *params,
	char *inname)
{
	char string[100];
        FILE *infile = fopen(inname,"r");

	CursorAfterString(infile,"Enter type of equation:");
        fscanf(infile,"%s",string);
        (void) printf(" %s\n",string);
        switch (string[0])
        {
        case 'l':
        case 'L':
	    params->eqn_type = LINEAR_EQN;
	    break;
        case 'b':
        case 'B':
	    params->eqn_type = BURGER_EQN;
	    break;
	default:
	    (void) printf("Unknown equation type!\n");
	    clean_up(ERROR);
	}
	CursorAfterString(infile,"Enter type of initial condition:");
        fscanf(infile,"%s",string);
        (void) printf(" %s\n",string);
        switch (string[0])
        {
        case 's':
        case 'S':
	    params->init_type = SQUARE;
	    break;
        case 'h':
        case 'H':
	    params->init_type = HAMP;
	    break;
        case 'w':
        case 'W':
	    params->init_type = WAVE;
	    break;
	default:
	    (void) printf("Unknown initial condition type!\n");
	    clean_up(ERROR);
	}
	CursorAfterString(infile,"Enter coefficient a:");
        fscanf(infile,"%lf",&params->a);
        (void) printf(" %f\n",params->a);
	fclose(infile);
}	/* end readWenoParams */

static void wave_func(
        double a,
        double time,
        double *x,
        double *u,
        int mesh_size)
{
        int i;
        double arg;

        for (i = 0; i < mesh_size; ++i)
        {
            arg = x[i] - a*time;
            u[i] = sin(2.0*PI*(arg - 1.0));
        }
}

static void hamp_func(
        double a,
        double time,
        double *x,
        double *u,
        int mesh_size)
{
        int i;
        double arg;

        for (i = 0; i < mesh_size; ++i)
        {
            arg = x[i] - a*time;
            if (arg >= 1.0 && arg <= 2.0)
                u[i] = sin(PI*(arg - 1.0))*sin(PI*(arg - 1.0));
            else
                u[i] = 0.0;
        }
}

static void square_func(
        double a,
        double time,
        double *x,
        double *u,
        int mesh_size)
{
        int i;
        double arg;

        for (i = 0; i < mesh_size; ++i)
        {
            arg = x[i] - a*time;
	    
            if (arg >= 1.0 && arg <= 2.0)
                u[i] = 1.0;
            else
                u[i] = 0.0;
        }
}

static void setSolutionType(
        INIT_TYPE init_type)
{
        switch (init_type)
        {
        case WAVE:
	    exact_soln = wave_func;
            break;
        case HAMP:
	    exact_soln = hamp_func;
            break;
        case SQUARE:
	    exact_soln = square_func;
            break;
        default:
            (void) printf("Unknown equation type!\n");
            clean_up(ERROR);
        }
}       /* end setSolutionType */

/*
*			      wave.cpp:
*
*	Sample wave equation code for AMS-528, Numerical PDE.
*	Coded by Xiaolin Li
*
*/

#include <FronTier.h>

enum    _NUM_SCHEME {
        ERROR_SCHEME    =       -1,
	EXACT_SOLUTION_ONLY,
        LW_CENTRAL_EXPLICIT,
        LW_CENTRAL_IMPLICIT,
        CRANK_NICOLSON,
        MORE_SCHEMES
};
typedef enum _NUM_SCHEME NUM_SCHEME;

enum    _OPTION_TYPE {
	ERROR_TYPE = -1,
	EROUPEAN_CALL,
	EROUPEAN_PUT,
        MORE_OPTIONS
};
typedef enum _OPTION_TYPE OPTION_TYPE;

struct _BS_PARAMS {
        double r;               	/* Interest rate */
        double sigma;               	/* Volatility */
        double E;               	/* Strike price */
        double T;               	/* Expiry time */
        NUM_SCHEME num_scheme;  	/* Numerical scheme choice */
	int option_type;		/* Option type */
};
typedef struct _BS_PARAMS BS_PARAMS;

#define		Nfunc(x)		(0.5*(1.0 + erf((x))))

char *in_name,*out_name;
/* for exact solution */
static void exact_call_soln(BS_PARAMS,double,double*,double*,int);
static void exact_put_soln(BS_PARAMS,double,double*,double*,int);

/* functions for different schemes */
static void expexp_step(double,double,double*,double*,double*,int,BS_PARAMS);
static void expimp_step(double,double,double*,double*,double*,int,BS_PARAMS);
static void impimp_step(double,double,double*,double*,double*,int,BS_PARAMS);
static void tridiagonal_solver(int,double*,double*,double*,double*,double*);
static void init_bs_params(char*,BS_PARAMS*);

int main(int argc, char **argv)
{ 
	F_BASIC_DATA f_basic; /* Only some basic features of f_basic is used */
	Front front;	      /* Only some basic feature of front is used */
	double *x,*u_old,*u_new,*u_sol; /* mesh points and solutions */
	double *x_movie,*y_movie;
	double dx,dt;
	char movie_caption[100];
	char time_label[100];
	char xg_name[200];
	char gd_name[200];
	FILE *xg_file;
	double xmin,xmax,umin,umax,height;
	void (*exact_soln)(BS_PARAMS,double,double*,double*,int);
	double R,*u_stencil;
	double L,U;
	double CFL;	/* CFL safety factor */
	int buffer_size = 1;   /* if needed, can be changed to other numbers */
	int i,mesh_size,out_size;
	int boundary_type[2];
	double b;		/* diffusion coefficient */
	BS_PARAMS bs_params;

	/* I/O names */
	FT_Init(argc,argv,&f_basic);
	in_name		= f_basic.in_name;
	out_name	= f_basic.out_name;

	/* Get domain information */
	FT_ReadSpaceDomain(in_name,&f_basic);
	L = f_basic.L[0];
	U = f_basic.U[0];
	mesh_size = f_basic.gmax[0];
	boundary_type[0] = f_basic.boundary[0][0];
	boundary_type[1] = f_basic.boundary[0][1];

	/* Time control of the run and I/O */
	FT_ReadTimeControl(in_name,&front);
	init_bs_params(in_name,&bs_params);

	if (bs_params.option_type == EROUPEAN_PUT)
	    exact_soln = exact_put_soln;
	else if (bs_params.option_type == EROUPEAN_CALL)
	    exact_soln = exact_call_soln;
	out_size = (int)(2*bs_params.E/U*mesh_size);

	CFL = Time_step_factor(&front);
	FT_ResetTime(&front);

	/* Allocating memory for arrays */
	FT_VectorMemoryAlloc((POINTER*)&x,mesh_size+1,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&u_old,mesh_size+1,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&u_new,mesh_size+1,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&u_sol,mesh_size+1,sizeof(double));

	/* Initialization of states */

	dx = (U - L)/mesh_size;

	front.dt = dt = 0.25*dx*dx/sqr(bs_params.sigma*U);
	FT_TimeControlFilter(&front);

	/* Set the initial condition */
	for (i = 0; i <= mesh_size; i++)
	    x[i] = L + i*dx;
	exact_soln(bs_params,front.time,x,u_old,mesh_size);

	/* Assign numerical scheme */

	/* Set frame margin for GD movie output */
	xmin = L;	xmax = U*out_size/mesh_size;
	umin = HUGE;	umax = -HUGE;
	for (i = 0; i <= out_size; i++)
	{
	    if (umin > u_old[i]) umin = u_old[i];
	    if (umax < u_old[i]) umax = u_old[i];
	}
	height = umax - umin;
	umin -= 0.15*height;	umax += 0.15*height;

	/* Open and initialize GD movie output */
	sprintf(movie_caption,"u vs. x");
	sprintf(gd_name,"%s.gif",out_name);
	gd_initplot(gd_name,movie_caption,xmin,xmax,umin,umax,2);

	/* Time loop */
	for (;;)
	{
	    /* Advancing numerical solution */
	    u_old[0] = 0.0;
	    u_old[mesh_size] = u_old[mesh_size-1] + dx;
	    printf("u_old[%d] = %f\n",mesh_size,u_old[mesh_size]);
	    if (bs_params.num_scheme == LW_CENTRAL_EXPLICIT)
	    	expexp_step(front.dt,dx,x,u_old,u_new,mesh_size,bs_params);
	    else if (bs_params.num_scheme == LW_CENTRAL_IMPLICIT)
	    	expimp_step(front.dt,dx,x,u_old,u_new,mesh_size,bs_params);
	    else if (bs_params.num_scheme == CRANK_NICOLSON)
	    	impimp_step(front.dt,dx,x,u_old,u_new,mesh_size,bs_params);

	    /* Swapping solution storage */
	    for (i = 0; i <= mesh_size; i++)
	    {
		u_old[i] = u_new[i];
	    }

	    /* Time and step control */
	    FT_AddTimeStepToCounter(&front);
	    (void) printf("\ntime = %20.14f   step = %5d   ",
                        	front.time,front.step);

	    /* Update boundary condition */

	    /* Movie frame */
	    if (FT_IsMovieFrameTime(&front))
	    {
		/* Numerical solution */
		if (bs_params.num_scheme != EXACT_SOLUTION_ONLY)
		{
		    x_movie = x+buffer_size;
		    y_movie = u_old+buffer_size;
		    gd_plotdata(mesh_size,x_movie,y_movie);
		}

		/* Exact solution */
		exact_soln(bs_params,front.time,x,u_sol,mesh_size);
		x_movie = x;
		y_movie = u_sol;
		gd_plotdata(out_size+1,x_movie,y_movie);

		/* Time label */
		sprintf(time_label,"Time = %6.3f",front.time);
		gd_plotframe(time_label);
	    }
	    /* Output date control */
	    if (FT_IsSaveTime(&front))
	    {
		/* Numerical solution */
		if (bs_params.num_scheme != EXACT_SOLUTION_ONLY)
		{
		    sprintf(xg_name,"%s-num_sol-%d.xg",out_name,front.ip);
		    xg_file = fopen(xg_name,"w");
		    fprintf(xg_file,"\"u vs. x\"\n");
		    for (i = 0; i <= mesh_size; ++i)
		    {
		    	fprintf(xg_file,"%f  %f\n",x[i],u_old[i]);
		    }
		    fclose(xg_file);
		}

		/* Exact solution */
		exact_soln(bs_params,front.time,x,u_sol,mesh_size);
		sprintf(xg_name,"%s-exc-%d.xg",out_name,front.ip);
		xg_file = fopen(xg_name,"w");
		fprintf(xg_file,"\"u vs. x\"\n");
		for (i = 0; i <= mesh_size; ++i)
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
}

static void exact_call_soln(
	BS_PARAMS bs_params,
        double tau,
        double *s,
        double *c,
        int mesh_size)
{
	double r,sigma,E;
        double d1,d2,N1,N2;
	int i;

	r = bs_params.r;
	sigma = bs_params.sigma;
	E = bs_params.E;
        for (i = 0; i <= mesh_size; ++i)
        {
	    if (tau == 0.0)
	    {
		if (s[i] < E) 
		    N1 = N2 = 0.0;
		else 
		    N1 = N2 = 1.0;
	    }
	    else
	    {
	    	if (s[i] == 0.0) 
		    N1 = N2 = 0.0;
		else
		{
	    	    d1 = (log(s[i]/E) + (r + 0.5*sqr(sigma))*tau)/sigma
		    				/sqrt(tau);
	    	    d2 = (log(s[i]/E) + (r - 0.5*sqr(sigma))*tau)/sigma
		    				/sqrt(tau);
	    	    N1 = Nfunc(d1);	N2 = Nfunc(d2);
		}
	    }
	    c[i] = s[i]*N1 - E*exp(-r*tau)*N2;
	}
}

static void exact_put_soln(
	BS_PARAMS bs_params,
        double tau,
        double *s,
        double *p,
        int mesh_size)
{
	int i;
	double r = bs_params.r;
	double E = bs_params.E;
	exact_call_soln(bs_params,tau,s,p,mesh_size);
	for (i = 0; i <= mesh_size; ++i)
	{
	    p[i] -= s[i] - E*exp(-r*tau);
	}
}

static void tridiagonal_solver(
	int n,
	double *a,
	double *b,
	double *c,
	double *d,
	double *u)
{
	static double *A,*B;
	int i;
	double denom;

	if (A == NULL)
	{
	    FT_VectorMemoryAlloc((POINTER*)&A,n-1,sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&B,n-1,sizeof(double));
	}

	A[n-2] = -a[n-1]/b[n-1];
	B[n-2] =  d[n-1]/b[n-1];

	for (i = n-2; i > 0; --i)
	{
	    denom = A[i]*c[i] + b[i];
	    A[i-1] = -a[i]/denom;
	    B[i-1] = (d[i] - B[i]*c[i])/denom;
	}

	u[0] = (d[0] - B[0]*c[0])/(A[0]*c[0] + b[0]);
	for (i = 0; i < n-1; ++i)
	{
	    u[i+1] = A[i]*u[i] + B[i];
	}
}
	
static void init_bs_params(
	char *in_name,
	BS_PARAMS *bs_params)
{
	FILE *infile = fopen(in_name,"r");
	char string[100];

	CursorAfterString(infile,"Enter option type:");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	if (string[0] == 'e' || string[0] == 'E')
	{
	    if (string[9] == 'c' || string[9] == 'C')
	    	bs_params->option_type = EROUPEAN_CALL;
	    else if (string[9] == 'p' || string[9] == 'P')
	    	bs_params->option_type = EROUPEAN_PUT;
	}
	CursorAfterString(infile,"Enter interest rate:");
	fscanf(infile,"%lf",&bs_params->r);
	(void) printf(" %f\n",bs_params->r);
	CursorAfterString(infile,"Enter volatility:");
	fscanf(infile,"%lf",&bs_params->sigma);
	(void) printf(" %f\n",bs_params->sigma);
	CursorAfterString(infile,"Enter strike price:");
	fscanf(infile,"%lf",&bs_params->E);
	(void) printf(" %f\n",bs_params->E);
	CursorAfterString(infile,"Enter expiry time:");
	fscanf(infile,"%lf",&bs_params->T);
	(void) printf(" %f\n",bs_params->T);

	CursorAfterString(infile,"Enter numerical scheme:");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	if (string[0] == 'e' || string[0] == 'E')
	{
	    bs_params->num_scheme = EXACT_SOLUTION_ONLY;
	}
	else if (string[0] == 'c' || string[0] == 'C')
	{
	    bs_params->num_scheme = CRANK_NICOLSON;
	}
	else if (string[0] == 'l' || string[0] == 'L')
	{
	    if (string[11] == 'e' || string[11] == 'E')
	    	bs_params->num_scheme = LW_CENTRAL_EXPLICIT;
	    else if (string[11] == 'i' || string[11] == 'I')
	    {
	    	bs_params->num_scheme = LW_CENTRAL_IMPLICIT;
	    }
	}
	else
	{
	    (void) printf("Unknown numerical scheme!\n");
	    clean_up(ERROR);
	}
}	/* end init_bs_params */

static void expexp_step(
	double dt,
	double ds,
	double *s,
	double *v_old,
	double *v_new,
	int mesh_size,
	BS_PARAMS bs_params)
{
	int i;
	double r = bs_params.r;
	double sigma = bs_params.sigma;
	double lambda = dt/ds;

	for (i = 1; i < mesh_size; ++i)
	{
	    v_new[i] = v_old[i] + 0.5*lambda*r*s[i]*(v_old[i+1] - v_old[i-1])
	    		+ 0.5*sqr(lambda*r*s[i])*(v_old[i+1] - 2.0*v_old[i]
			+ v_old[i-1])
	  		+ 0.5*dt*sqr(sigma*s[i])*(v_old[i+1] - 2.0*v_old[i] + 
			v_old[i-1])/sqr(ds) - r*dt*v_old[i];
	}
}

static void impimp_step(
	double dt,
	double ds,
	double *s,
	double *u_old,
	double *u_new,
	int mesh_size,
	BS_PARAMS bs_params)
{
	int i,n;
	static double *a,*b,*c,*d;
	double *pu_old,*pu_new;

	if (a == NULL) /* Allocate memory for coefficients */
	{
	    FT_VectorMemoryAlloc((POINTER*)&a,n,sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&b,n,sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&c,n,sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&d,n,sizeof(double));
	}

	tridiagonal_solver(n,a,b,c,d,pu_new);
}

static void expimp_step(
	double dt,
	double ds,
	double *s,
	double *u_old,
	double *u_new,
	int mesh_size,
	BS_PARAMS bs_params)
{
	int i,n;
	static double *a,*b,*c,*d;
	double *pu_old,*pu_new;

	if (a == NULL) /* Allocate memory for coefficients */
	{
	    FT_VectorMemoryAlloc((POINTER*)&a,n,sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&b,n,sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&c,n,sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&d,n,sizeof(double));
	}

	tridiagonal_solver(n,a,b,c,d,pu_new);
}

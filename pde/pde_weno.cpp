/*
*				pde_weno.c:
*
*                        5-th order Weno method for 1D invisid Burgers Equation
*/

#include <FronTier.h>

char *out_name;


double   dx,dt;
double   valuemax;
double   L = -2.0;
double   U =  2.0;
double   CFLnumber = 0.5 ;
int      mesh_size = 401;
int      ghost_size = 3; // number of ghost points, related to the order of accuracy (fifth order)
int      shift_size = mesh_size + ghost_size;
int      extended_size = mesh_size + 2*ghost_size;

double   end_time = 10.0;	//static void wave_func(double,double,double*,double*,int);

/* functions for different schemes */


static void GetFlux_Weno (double *x, double *u_old, double *u_new, double *rhs);
/* Function for getting flux*/
static void Runge_Kutta (double *x, double *u_old, double *u_new, double *rhs);
/* Function for discretisize time (3-th order Runge_Kutta) */

int main(int argc, char **argv)
{
       	
	F_BASIC_DATA f_basic;
	double *x,*u_old,*u_new,*u_sol;
	double *rhs; //right hand side of the space-discretization equation
	double *x_draw;
	double tttime;
	int ieno = 2; //switch for choosing eno or weno
	int i,j,n,movie_interval;
	char movie_caption[100];
	char time_label[100];
	char xg_name[200];
	char gd_name[200];
	FILE *xg_file;
	double xmin, xmax, umin, umax, height;

	printf("The left boundary:");
	scanf("%lf", &L);

	printf("The right boundary:");
	scanf("%lf", &U);

	printf("The number of points:");
	scanf("%d", &mesh_size);


	shift_size = mesh_size + ghost_size;
	extended_size = mesh_size + 2*ghost_size;

	printf("CFL number:");
	scanf("%lf", &CFLnumber);

	printf("The end time:");
	scanf("%lf", &end_time);

	printf ("number of points: %d,  end time: %f,  CFL number: %f    \n", mesh_size, end_time, CFLnumber);

	/* Allocating memory for arrays */
	FT_VectorMemoryAlloc((POINTER*)&x,     extended_size,sizeof(double));  // extended interval
	FT_VectorMemoryAlloc((POINTER*)&u_old, extended_size,sizeof(double));  //old time step solution (extended)
	FT_VectorMemoryAlloc((POINTER*)&u_new, mesh_size,sizeof(double));  //new time step solution
	FT_VectorMemoryAlloc((POINTER*)&u_sol, mesh_size,sizeof(double));  //real solution array
	
	FT_VectorMemoryAlloc((POINTER*)&rhs,   mesh_size,sizeof(double)); //right hand side in R-K process
	FT_VectorMemoryAlloc((POINTER*)&x_draw, mesh_size, sizeof(double)); // for the movie x value





	/* Will take output file name from command line */
	FT_Init(argc,argv,&f_basic);
	out_name = f_basic.out_name;

	/* Initialization of states */

	dx = (U - L)/(mesh_size - 1);
	dt = CFLnumber*dx;
	movie_interval = (int)(0.15/dt);

	/* Set the initial condition */
	for (i = -ghost_size; i < shift_size; i++)
	    x[i + ghost_size] = L + i*dx; // desiganate the value of x[i] (including ghost points)
	for (i = 0; i < mesh_size; i++)
	    x_draw[i] = L + i*dx;

	for (i = 0; i < mesh_size; i++) // The intial value for function
	{
	  u_new[i] = 0.25 + 0.5*sin(PI*x_draw[i]);
	  u_old[i + ghost_size] = u_new[i];
	}
	for (i = 0; i <= ghost_size; i++)
	{

	  u_old[ghost_size - i] = u_old[ghost_size - i + mesh_size - 1];
	  u_old[ghost_size + mesh_size - 1 + i] = u_old[ghost_size + i];
	}

	//wave_func(1.0,0.0,x,u_old,mesh_size);

	/* Assign numerical scheme */

	//numerical_scheme = upwind;

	/* Set frame margin for GD output */

	xmin = L;	xmax = U;
	umin = HUGE;	umax = -HUGE;
	for (i = 0; i < mesh_size; i++)
	{
	    if (umin > u_new[i]) umin = u_new[i];
	    if (umax < u_new[i]) umax = u_new[i];
	}
	height = umax - umin;
	umin -= 0.5*height;	umax += 0.5*height;

	/* Open and initialize GD output */


/*just a test for GD package*/
	sprintf(movie_caption,"u vs. x");
	sprintf(gd_name,"%s.gif",out_name);
	gd_initplot(gd_name,movie_caption,xmin,xmax,umin,umax,2);
	gd_plotdata(mesh_size,x_draw,u_new);
	sprintf(time_label,"Time = 0.00000");
        gd_plotframe(time_label);


/*Main Time loop */

	n = 0; tttime = 0.0; 

	while (tttime <= end_time)
	{

		n++;
		Runge_Kutta(x, u_old, u_new, rhs);
		
		if (n%movie_interval == 0)
		{
			gd_plotdata(mesh_size,x_draw,u_new);
			sprintf(time_label, "Time=%6.3f",tttime);
			gd_plotframe(time_label);
		}
		tttime += dt;
	}
	gd_closeplot();

	return 0;
}


static void GetFlux_Weno (double *x, double *u_old, double *u_new, double *rhs)
{
	int i,j,k;
	//int all_size = mesh_size + 2*ghost_size;
	//int shift_size = mesh_size + ghost_size;
	double em;
	double fp;
	double epweno = 1.0E-06;

	double *f; // f(u[i]);
	double *flux; 
	double *dfp,*dfm;
	double cdx = (double)(1.0/dx);
	//printf("%f\n", cdx);
	double hh[5][3];
	double t1,t2,t3,tt1,tt2,tt3,s1,s2,s3,t0;

	
	FT_VectorMemoryAlloc((POINTER*)&f, extended_size, sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&flux, extended_size, sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&dfp, extended_size, sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&dfm, extended_size, sizeof(double));


	//for (i=0; i<=ghost_size; i++) //deal with boundary
	//{
	//	u_old[ghost_size-i]=u_old[ghost_size-i+mesh_size-1];
	//	u_old[ghost_size+mesh_size-1+i]=u_old[ghost_size+i];
	//}

	em = 0.0;

	for (i = 0; i < extended_size; i++)
	{
		f[i] = 0.5*u_old[i]*u_old[i];
	//	f[i] = 0.5 * u_old[i];
		fp = fabs(u_old[i]);
		if (em <= fp)
			em = fp;

	}
	valuemax = em;
	//printf("maxvalue=%f\n", valuemax);

	for (i = -ghost_size; i < shift_size - 1; i++)
	{
		j = i + ghost_size;
		dfp[j] = 0.5*( f[j+1] - f[j] + em*(u_old[j+1] - u_old[j]));
		dfm[j] = f[j+1] - f[j] - dfp[j];
	}

	for (i = -1; i < mesh_size; i++)
	{
		j = i + ghost_size;
		hh[1][1] = dfp[j-2];
		hh[2][1] = dfp[j-1];
		hh[3][1] = dfp[j];
		hh[4][1] = dfp[j+1];
		hh[1][2] = -dfm[j+2];
		hh[2][2] = -dfm[j+1];
		hh[3][2] = -dfm[j];
		hh[4][2] = -dfm[j-1];

		flux[j] = ( -f[j-1] + 7.0*(f[j] + f[j+1]) - f[j+2])/12.0;

		for (k = 1; k <= 2; k++)
		{
			t1 = hh[1][k] - hh[2][k];
			t2 = hh[2][k] - hh[3][k];
			t3 = hh[3][k] - hh[4][k];

			tt1 = 13.0*t1*t1 + 3.0*(hh[1][k] - 3.0*hh[2][k])*(hh[1][k] - 3.0*hh[2][k]);
			tt2 = 13.0*t2*t2 + 3.0*(hh[2][k] +  hh[3][k])*(hh[2][k] + hh[3][k]);
			tt3 = 13.0*t3*t3 + 3.0*(3.0*hh[3][k] - hh[4][k])*(3.0*hh[3][k] - hh[4][k]);

			tt1 = (epweno + tt1)*(epweno + tt1);
			tt2 = (epweno + tt2)*(epweno + tt2);
			tt3 = (epweno + tt3) * (epweno + tt3);

			s1 = tt2*tt3;
			s2 = 6.0*tt1*tt3;
			s3 = 3.0*tt1*tt2;
			
			t0 = 1.0/( s1 + s2 + s3);

			s1 = s1*t0;
			s3 = s3*t0;

			flux[j] = flux[j] + (s1*(t2 - t1) + (0.5*s3 - 0.25)*(t3 - t2))/3.0;

		}

		
	}


	for (i = 0; i < mesh_size; i++)
	{
		j = i + ghost_size;
		rhs[i] = ( flux[j-1] - flux[j] )*cdx;
	}

	free(f);free(flux);
	FT_FreeThese(2,dfp,dfm);


}



static void Runge_Kutta (double *x, double *u_old, double *u_new, double *rhs)
{
	double *u_mid1,*u_mid2;
	//double *u_mid3,*u_mid4;

	FT_VectorMemoryAlloc((POINTER*)&u_mid1, extended_size, sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&u_mid2, extended_size, sizeof(double));
	//FT_VectorMemoryAlloc((POINTER*)&u_mid3, mesh_size, sizeof(double));
	//FT_VectorMemoryAlloc((POINTER*)&u_mid4, mesh_size, sizeof(double));

	int i;

	GetFlux_Weno(x,u_old,u_new,rhs);

	dt = dx*CFLnumber/valuemax;

	for (i = 0; i < mesh_size; i++)
		u_mid1[i + ghost_size] = u_new[i] + dt*rhs[i];
	for (i = 0; i <= ghost_size; i++)//deal with boundary
	{
		u_mid1[ghost_size - i] = u_mid1[ghost_size - i + mesh_size - 1];
		u_mid1[ghost_size + mesh_size - 1 + i] = u_mid1[ghost_size + i];
	}
//////////////////////////////////////////////////////////////////////////////
	GetFlux_Weno(x, u_mid1, u_new, rhs);

	for (i = 0; i < mesh_size; i++)
		u_mid2[i + ghost_size] = 0.75*u_new[i] + 0.25*(u_mid1[i + ghost_size] + dt*rhs[i]);
	for (i = 0; i <= ghost_size; i++)//deal with boundary
	{
		u_mid2[ghost_size - i] = u_mid2[ghost_size - i + mesh_size - 1];
		u_mid2[ghost_size + mesh_size - 1 + i] = u_mid2[ghost_size + i];
	}
////////////////////////////////////////////////////////////////////////////////////////
	GetFlux_Weno(x,u_mid2,u_new,rhs);
	for (i = 0; i < mesh_size; i++)
		u_new[i] = (u_new[i] + 2*(u_mid2[i + ghost_size] + dt*rhs[i]))/3.0;
	for (i = 0; i < mesh_size; i++)
		u_old[i + ghost_size] = u_new[i];
//////////////////////////////////////////////////////////////////////////////////////////

	for (i = 0; i <= ghost_size; i++)//deal with boundary
	{
	  u_old[ghost_size - i] = u_old[ghost_size - i + mesh_size - 1];
	  u_old[ghost_size + mesh_size - 1 + i] = u_old[ghost_size + i];
	}

	FT_FreeThese(2,u_mid1,u_mid2);
}





















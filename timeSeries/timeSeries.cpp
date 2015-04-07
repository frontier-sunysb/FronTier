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

#include "timeSeries.h"

/* Not yet in use */
static void initXgraphPlots(char*,int*,double**);
static void makeXgraphPlot(char*,int,double*,double*,int,double);

static void makeSolnMovie(char*,double*,double*,int,double,double,double,
				double,double,char*);
static void readTimeParams(TIME_PARAMS*,char*);
static void getMovieBound(TIME_PARAMS*,double*,double*);

static void time_driver(Front*);
static double timeSeriesEvent(Front*);

char *in_name,*out_name;

int main(int argc, char **argv)
{
	static Front front;
	static F_BASIC_DATA f_basic;
	static TIME_PARAMS ts_params;

	FT_Init(argc,argv,&f_basic);
	f_basic.dim = 1;
	in_name                 = f_basic.in_name;
	out_name                = f_basic.out_name;

	FT_ReadTimeControl(in_name,&front);
	readTimeParams(&ts_params,in_name);
	front.extra1 = &ts_params;

	time_driver(&front);

	clean_up(0);
}

static void time_driver(
	Front *front)
{
	double *ts_time,*ts_events;
	int size;
	TIME_PARAMS *ts_params = (TIME_PARAMS*)front->extra1;
	double dt;
	double tmin,tmax,umin,umax;
	char caption[100];
	int n = 0;

	tmin = 0.0;
	tmax = front->max_time;
	getMovieBound(ts_params,&umin,&umax);
	sprintf(caption,"Time series");

	FT_ResetTime(front);
	front->dt = dt = 1.0/ts_params->T_density;
	size = (int)front->max_time/dt + 500;
	FT_VectorMemoryAlloc((POINTER*)&ts_time,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&ts_events,size,sizeof(double));

	FT_SetOutputCounter(front);
	FT_TimeControlFilter(front);
	FT_AddTimeStepToCounter(front);

	ts_time[n] = front->time;
	ts_events[n] = timeSeriesEvent(front);
	makeSolnMovie(out_name,ts_time,ts_events,n,tmin,tmax,umin,umax,
				front->time,caption);
	n++;

	for (;;)
	{

	    ts_time[n] = front->time;
	    ts_events[n] = timeSeriesEvent(front);

	    FT_AddTimeStepToCounter(front);
	    front->dt = dt;

	    if (FT_IsDrawTime(front))
            {
		makeSolnMovie(out_name,ts_time,ts_events,n,tmin,tmax,
				umin,umax,front->time,caption);
	    }
	    if (FT_IsSaveTime(front))
            {
		;
	    }
	    if (FT_TimeLimitReached(front))
            {
		break;
	    }
	    FT_TimeControlFilter(front);
	    n++;
	}
}	/* end time_driver */

static double timeSeriesEvent(
	Front *front)
{
	TIME_PARAMS *ts_params = (TIME_PARAMS*)front->extra1;
	POINTER pdf_params;
	unsigned short int *seeds;
	double (*random_func)(POINTER,unsigned short int*);

	random_func = ts_params->random_func;
	pdf_params = ts_params->pdf_params;
	seeds = ts_params->seeds;

	switch (ts_params->ts_type)
	{
	case WHITE_NOISE:
	    return (*random_func)(pdf_params,seeds);
	default:
	    clean_up(ERROR);
	}
}	/* end timeSeriesEvent */

static void makeXgraphPlot(
	char *out_name,
	int N,
	double *x,
	double *u,
	int l,
	double time)
{
	char xname[256];
	FILE *xfile;
	int i;

	sprintf(xname,"%s/soln-%d.xg",out_name,l);
	xfile = fopen(xname,"w");
	fprintf(xfile,"\"time = %6.3f\"\n",time);
	for (i = 0; i < N; ++i)
	    	fprintf(xfile,"%f %f\n",x[i],u[i]);
	fclose(xfile);
}	/* end makeXgraphPlot */

static void makeSolnMovie(
	char *out_name,
	double *x,
	double *u,
	int N,
	double xmin,
	double xmax,
	double umin,
	double umax,
	double time,
	char *caption)
{
	static boolean first = YES;
	char time_label[100];
	char gd_name[256];

	sprintf(gd_name,"%s/timeSeries.gif",out_name);

	sprintf(time_label,"t = %6.3f",time);
	if (first)
	{
	    first = NO;
	    gd_initplot(gd_name,caption,xmin,xmax,umin,umax,3);
	    gd_plotdata(N+1,x,u);
	    gd_plotframe(time_label);
	}
	else
	{
	    gd_appendplot(gd_name,caption,xmin,xmax,umin,umax,3);
	    gd_plotdata(N+1,x,u);
	    gd_plotframe(time_label);
	}
}	/* end makeSolnMovie */

static void initXgraphPlots(
	char *inname,
	int *N,
	double **time)
{
	char string[100];
	FILE *infile = fopen(inname,"r");
	static double *tstore;
	int i,n;
	double T0,T1,dT;

	CursorAfterString(infile,"Enter number of time frames:");
	fscanf(infile,"%d",&n);
        (void) printf(" %d\n",n);
	FT_VectorMemoryAlloc((POINTER*)&tstore,n,sizeof(double));
	CursorAfterString(infile,"Enter yes for non-uniform time sequence:");
	fscanf(infile,"%s",string);
        (void) printf(" %s\n",string);
	if (string[0] == 'y' || string[0] == 'Y')
	{
	    for (i = 0; i < n; ++i)
	    {
		sprintf(string,"Enter time for %d-th frame:",i);
		CursorAfterString(infile,string);
		fscanf(infile,"%lf",&tstore[i]);
		printf(" %f\n",tstore[i]);
	    }
	}
	else
	{
	    CursorAfterString(infile,"Enter start time:");
	    fscanf(infile,"%lf",&T0);
	    printf("%f\n",T0);
	    CursorAfterString(infile,"Enter end time:");
	    fscanf(infile,"%lf",&T1);
	    printf("%f\n",T1);
	    dT = (T1 - T0)/(n-1);
	    for (i = 0; i < n; ++i)
		tstore[i] = T0 + n*dT;
	}
	fclose(infile);

	*N = n;
	*time = tstore;
}	/* initXgraphPlots */

static void readTimeParams(
	TIME_PARAMS *params,
	char *inname)
{
	static GAUSS_PARAMS *gauss_params;
	static EXP_PARAMS *exp_params;
	static POWER_PARAMS *power_params;
	static UNIFORM_PARAMS *uniform_params;
	static STABLE_PARAMS *stable_params;
	char string[100];
        FILE *infile = fopen(inname,"r");

	CursorAfterString(infile,"Enter time series type:");
        fscanf(infile,"%s",string);
        (void) printf(" %s\n",string);
	switch (string[0])
	{
	case 'w':
	case 'W':
	    params->ts_type = WHITE_NOISE;
	    break;
	default:
	    (void) printf("Unknown time series type!\n");
	    clean_up(ERROR);
	}
	CursorAfterString(infile,"Enter name of distribution function:");
        fscanf(infile,"%s",string);
        (void) printf(" %s\n",string);
	switch (string[0])
	{
	case 'g':
	case 'G':
	    FT_ScalarMemoryAlloc((POINTER*)&gauss_params,sizeof(GAUSS_PARAMS));
	    switch (string[6])
	    {
	    case 'n':
	    case 'N':
	    	params->rand_type = GAUSS_NEWTON;
		params->random_func = gauss_newton;
		break;
	    case 'b':
	    case 'B':
	    	params->rand_type = GAUSS_BM;
		params->random_func = gauss_box_muller;
		break;
	    case 'c':
	    case 'C':
	    	params->rand_type = GAUSS_CL;
		params->random_func = gauss_center_limit;
		break;
	    }
	    CursorAfterString(infile,"Enter mathematical expectation:");
            fscanf(infile,"%lf",&gauss_params->mu);
            (void) printf(" %f\n",gauss_params->mu);
	    CursorAfterString(infile,"Enter standard deviation:");
            fscanf(infile,"%lf",&gauss_params->sigma);
            (void) printf(" %f\n",gauss_params->sigma);
	    params->pdf_params = (POINTER)gauss_params;
	    break;
	case 'e':
	case 'E':
	    FT_ScalarMemoryAlloc((POINTER*)&exp_params,sizeof(EXP_PARAMS));
	    params->rand_type = EXPONENTIAL;
	    CursorAfterString(infile,"Enter parameter lambda:");
            fscanf(infile,"%lf",&exp_params->lambda);
            (void) printf(" %f\n",exp_params->lambda);
	    params->pdf_params = (POINTER)exp_params;
	    params->random_func = dist_exponential;
	    break;
	case 'c':
	case 'C':
	    params->rand_type = CAUCHY;
	    FT_ScalarMemoryAlloc((POINTER*)&gauss_params,sizeof(GAUSS_PARAMS));
	    CursorAfterString(infile,"Enter mathematical expectation:");
            fscanf(infile,"%lf",&gauss_params->mu);
            (void) printf(" %f\n",gauss_params->mu);
	    CursorAfterString(infile,"Enter standard deviation:");
            fscanf(infile,"%lf",&gauss_params->sigma);
            (void) printf(" %f\n",gauss_params->sigma);
	    params->pdf_params = (POINTER)gauss_params;
	    params->random_func = dist_cauchy;
	    break;
	case 'p':
	case 'P':
	    FT_ScalarMemoryAlloc((POINTER*)&power_params,sizeof(POWER_PARAMS));
	    params->rand_type = POWER;
	    CursorAfterString(infile,"Enter the power:");
            fscanf(infile,"%d",&power_params->power);
            (void) printf(" %d\n",power_params->power);
	    params->pdf_params = (POINTER)power_params;
	    params->random_func = dist_power;
	    break;
	case 'u':
	case 'U':
	    FT_ScalarMemoryAlloc((POINTER*)&uniform_params,
				sizeof(UNIFORM_PARAMS));
	    params->rand_type = UNIFORM;
	    CursorAfterString(infile,"Enter the lower and upper bounds:");
            fscanf(infile,"%lf %lf",&uniform_params->a,&uniform_params->b);
            (void) printf(" %f %f\n",uniform_params->a,uniform_params->b);
	    params->pdf_params = (POINTER)uniform_params;
	    params->random_func = dist_uniform;
	    break;
	case 's':
	case 'S':
	    FT_ScalarMemoryAlloc((POINTER*)&stable_params,
				sizeof(STABLE_PARAMS));
	    params->rand_type = STABLE;
	    CursorAfterString(infile,"Enter alpha:");
            fscanf(infile,"%lf",&stable_params->alpha);
            (void) printf(" %f\n",stable_params->alpha);
	    CursorAfterString(infile,"Enter beta:");
            fscanf(infile,"%lf",&stable_params->beta);
            (void) printf(" %f\n",stable_params->beta);
	    CursorAfterString(infile,"Enter sigma:");
            fscanf(infile,"%lf",&stable_params->sigma);
            (void) printf(" %f\n",stable_params->sigma);
	    CursorAfterString(infile,"Enter mu:");
            fscanf(infile,"%lf",&stable_params->mu);
            (void) printf(" %f\n",stable_params->mu);
	    params->pdf_params = (POINTER)stable_params;
	    params->random_func = dist_stable;
	    break;
	case 'm':
	case 'M':
	    params->rand_type = MIDDLE;
	    params->random_func = dist_middle;
	    break;
	default:
	    (void) printf("Unknown random variable type!\n");
	    clean_up(ERROR);
	}

	CursorAfterString(infile,"Enter type of random seed:");
        fscanf(infile,"%s",string);
        (void) printf(" %s\n",string);
	switch (string[0])
	{
	case 'f':
	case 'F':
	    params->seed_type = FIXED_SEED;
	    break;
	case 'r':
	case 'R':
	    params->seed_type = RANDOM_SEED;
	    break;
	case 'i':
	case 'I':
	    params->seed_type = INPUT_SEED;
	    break;
	}
	if (params->seed_type == INPUT_SEED)
	{
	    CursorAfterString(infile,"Enter three random seeds:");
            fscanf(infile,"%hu %hu %hu",&params->seeds[0],&params->seeds[1],
					&params->seeds[2]);
            (void) printf(" %d %d %d\n",params->seeds[0],params->seeds[1],
					params->seeds[2]);
	}
	/* The following is for stock simulation */
	CursorAfterString(infile,"Enter incident density:");
        fscanf(infile,"%d",&params->T_density);
	(void) printf(" %d\n",params->T_density);

	fclose(infile);
}	/* end readTimeParams */

static void assignRandomSeeds(
	TIME_PARAMS *params)
{
	time_t seconds;

	switch (params->seed_type)
	{
	case INPUT_SEED:
	    return;
	case FIXED_SEED:
	    params->seeds[0] = 271;
	    params->seeds[1] = 6253;
	    params->seeds[2] = 176;
	    return;
	case RANDOM_SEED:
	    seconds = time(NULL);
	    params->seeds[0] = seconds%9829;
	    params->seeds[1] = seconds%9883;
	    params->seeds[2] = seconds%9743;
	    return;
	}
}	/* end assignRandomSeeds */

static void getMovieBound(
	TIME_PARAMS *ts_params,
	double *umin,
	double *umax)
{
	GAUSS_PARAMS *gauss_params;

	switch (ts_params->rand_type)
	{
	case GAUSS_NEWTON:
	case GAUSS_BM:
	case GAUSS_CL:
	    gauss_params = (GAUSS_PARAMS*)ts_params->pdf_params;
	    *umin = gauss_params->mu - 5.0*gauss_params->sigma;
	    *umax = gauss_params->mu + 5.0*gauss_params->sigma;
	    break;
	default:
	    (void) printf("In getMovieBound(), unknown random type!\n");
	    clean_up(ERROR);
	}
}	/* end getMovieBound */

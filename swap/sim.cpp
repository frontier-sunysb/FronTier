#include "swap.h"

static void SlopeBasedTrade(DATA_SET*);
static void Slope1BasedTrade(DATA_SET*);
static void IntegralSlopeTrade(DATA_SET*);

extern void InvestSimulation(
	DATA_SET *data)
{
	int i_method;
	DATA_SET *copy_data = CopyData(data);
	printf("Available trading methods are:\n");
	printf("\tSlope-based trading (0)\n");
	printf("\tLeast-Square-Slope-based trading (1)\n");
	printf("\tIntegral-Slope trading (2)\n");
	printf("Enter trading method: ");
	scanf("%d",&i_method);
	switch (i_method)
	{
	case 0:
	    SlopeBasedTrade(copy_data);
	    break;
	case 1:
	    Slope1BasedTrade(copy_data);
	    break;
	case 2:
            IntegralSlopeTrade(copy_data);
	    break;
	default:
	    printf("Unknown trading method\n");
	    return;
	}
}	/* end InvestSimulation */

static void SlopeBasedTrade(
	DATA_SET *data)
{
	int i,j,i1,i2,is,N;
	int M = data->num_assets-1;
	double LA,LB;
	double QA,QB,QC;
	double S0[5],S[5];
	double s0[5],s1[5],q1[5];
	double s0_max,s1_max;
	double *value;
	double time[MAX_TRACE];
	int num_swap = 0;
	int is1_max,is0_max;
	double total_value0,total_value;
	char string[100];
	double dt;
	int ns1;	/* number of point for secondary slope */
	boolean cash_swap = NO;
	double cash;

start_trade:

	N = data->num_backtrace - 5;
	if (N > MAX_TRACE - 5) N = MAX_TRACE - 5;
	printf("Enter time interval: ");
	scanf("%lf",&dt);
	printf("Enter number of point for secondary slope: ");
	scanf("%d",&ns1);
	printf("Enter yes to include cash swap: ");
	scanf("%s",string);
	cash_swap = NO;
	cash = 0.0;
	num_swap = 0;
	if (string[0] == 'y' || string[0] == 'Y')
	{
	    cash_swap = YES;
	}

	total_value0 = 0.0;
	is = data->num_days - N;
	for (i = 0; i < M; ++i)
	{
	    printf("For stock %s, enter starting number of shares: ",
				data->assets[i].asset_name);
	    scanf("%lf",S0+i);
	    total_value0 += S0[i]*data->assets[i].value[is];
	}
	printf("\nTotal starting value = %f\n",total_value0);
	total_value = 0.0;
	is = data->num_days - 1;
	for (i = 0; i < M; ++i)
	{
	    total_value += S0[i]*data->assets[i].value[is];
	}
	printf("Total number of steps = %d\n",N);
	printf("Total untraded end value = %f\n",total_value);
	printf("Change: %f percent\n",100*(total_value/total_value0 - 1.0));
	printf("\n");

	for (i = 0; i < N; ++i) 
	    time[i] = i*dt;

	for (j = 0; j < N; ++j)
	{
	    is = data->num_days - N + j - (ns1 - 1);
	    s0_max = s1_max = -HUGE;
	    for (i = 0; i < M; ++i)
	    {
	    	value = data->assets[i].norm_value + is;
		s0[i] = value[ns1-1] - value[ns1-2];
		S[i] = data->assets[i].value[is+(ns1-1)];
		LeastSquareLinear(time,value,ns1,&LA,&LB);
		LeastSquareQuadr(time,value,ns1,&QA,&QB,&QC);
		s1[i] = LA;
		q1[i] = QA;
		if (s0_max < s0[i]) 
		{
		    s0_max = s0[i];
		    is0_max = i;
		}
		if (s1_max < s1[i]) 
		{
		    s1_max = s1[i];
		    is1_max = i;
		}
	    }
	    for (i1 = 0; i1 < M; ++i1)
	    {
		if (cash_swap)
		{
		    if (S0[i1] != 0.0 && s0[is0_max] < 0.0 &&
			s1[i1] < 0.0)
		    {
		    	cash += S0[i1]*S[i1];
		    	S0[i1] = 0.0;
		    	num_swap++;
		    }
		    if (cash != 0.0 && s0[is0_max] > 0.0 && 
			i1 == is0_max && s1[i1] > 0.0)
		    {
			S0[i1] += cash/S[i1];
			cash = 0.0;
		    	num_swap++;
		    }
		}
		if (
		    S0[i1] != 0.0 
		    && i1 != is0_max
		    && s1[i1] < s1[is0_max]
		   )
		{
		    S0[is0_max] += S0[i1]*S[i1]/S[is0_max];
		    S0[i1] = 0.0;
		    num_swap++;
		}
	    }
	}
	total_value = 0.0;
	is = data->num_days - 1;
	for (i = 0; i < M; ++i)
	{
	    printf("For stock %s, end number of shares: %f\n",
				data->assets[i].asset_name,S0[i]);
	    total_value += S0[i]*data->assets[i].value[is];
	}
	total_value += cash;
	printf("Total traded end value = %f\n",total_value);
	printf("Change: %f percent\n",100*(total_value/total_value0 - 1.0));
	printf("Number of swaps = %d\n",num_swap);
	printf("Type yes to restart simulation: ");
	scanf("%s",string);
	if (string[0] == 'y' || string[0] == 'Y')
	    goto start_trade;
}	/* end SlopeBasedTrade */

static void Slope1BasedTrade(
	DATA_SET *data)
{
	int i,j,i1,i2,is,N;
	int M = data->num_assets-1;
	double LA,LB;
	double QA,QB,QC;
	double S0[5],S[5];
	double s0[5],s1[5],q1[5];
	double s0_max,s1_max;
	double *value;
	double time[MAX_TRACE];
	int num_swap = 0;
	int is1_max,is0_max;
	double total_value0,total_value;
	char string[100];
	double dt;
	int ns1;	/* number of point for secondary slope */
	boolean cash_swap = NO;
	double cash;

start_trade:

	N = data->num_backtrace - 5;
	if (N > MAX_TRACE - 5) N = MAX_TRACE - 5;
	printf("Enter time interval: ");
	scanf("%lf",&dt);
	printf("Enter number of point for secondary slope: ");
	scanf("%d",&ns1);
	printf("Enter yes to include cash swap: ");
	scanf("%s",string);
	cash_swap = NO;
	cash = 0.0;
	num_swap = 0;
	if (string[0] == 'y' || string[0] == 'Y')
	{
	    cash_swap = YES;
	}

	total_value0 = 0.0;
	is = data->num_days - N;
	for (i = 0; i < M; ++i)
	{
	    printf("For stock %s, enter starting number of shares: ",
				data->assets[i].asset_name);
	    scanf("%lf",S0+i);
	    total_value0 += S0[i]*data->assets[i].value[is];
	}
	printf("\nTotal starting value = %f\n",total_value0);
	total_value = 0.0;
	is = data->num_days - 1;
	for (i = 0; i < M; ++i)
	{
	    total_value += S0[i]*data->assets[i].value[is];
	}
	printf("Total untraded end value = %f\n",total_value);
	printf("Change: %f percent\n",100*(total_value/total_value0 - 1.0));
	printf("\n");

	for (i = 0; i < N; ++i) 
	    time[i] = i*dt;

	for (j = 0; j < N; ++j)
	{
	    is = data->num_days - N + j - (ns1 - 1);
	    s0_max = s1_max = -HUGE;
	    for (i = 0; i < M; ++i)
	    {
	    	value = data->assets[i].norm_value + is;
		s0[i] = value[ns1-1] - value[ns1-2];
		S[i] = data->assets[i].value[is+(ns1-1)];
		LeastSquareLinear(time,value,ns1,&LA,&LB);
		LeastSquareQuadr(time,value,ns1,&QA,&QB,&QC);
		s1[i] = LA;
		q1[i] = QA;
		if (s0_max < s0[i]) 
		{
		    s0_max = s0[i];
		    is0_max = i;
		}
		if (s1_max < s1[i]) 
		{
		    s1_max = s1[i];
		    is1_max = i;
		}
	    }
	    for (i1 = 0; i1 < M; ++i1)
	    {
		if (cash_swap)
		{
		    if (S0[i1] != 0.0 && s0[is0_max] < 0.0 &&
			s1[i1] < 0.0)
		    {
		    	cash += S0[i1]*S[i1];
		    	S0[i1] = 0.0;
		    	num_swap++;
		    }
		    if (cash != 0.0 && s0[is0_max] > 0.0 && 
			i1 == is0_max && s1[i1] > 0.0)
		    {
			S0[i1] += cash/S[i1];
			cash = 0.0;
		    	num_swap++;
		    }
		}
		if (
		    S0[i1] != 0.0 
		    && i1 != is1_max
		    //&& s1[i1] < s1[is0_max]
		   )
		{
		    S0[is1_max] += S0[i1]*S[i1]/S[is1_max];
		    S0[i1] = 0.0;
		    num_swap++;
		}
	    }
	}
	total_value = 0.0;
	is = data->num_days - 1;
	for (i = 0; i < M; ++i)
	{
	    printf("For stock %s, end number of shares: %f\n",
				data->assets[i].asset_name,S0[i]);
	    total_value += S0[i]*data->assets[i].value[is];
	}
	total_value += cash;
	printf("Total traded end value = %f\n",total_value);
	printf("Change: %f percent\n",100*(total_value/total_value0 - 1.0));
	printf("Number of swaps = %d\n",num_swap);
	printf("Type yes to restart simulation: ");
	scanf("%s",string);
	if (string[0] == 'y' || string[0] == 'Y')
	    goto start_trade;

} 	/* end Slope1BasedTrade */

#define		MAX_SIM_NUM		100
static void IntegralSlopeTrade(
	DATA_SET *data)
{
	int i,j,i1,i2,is,N;
	int M = data->num_assets;
	double LA,LB;
	double S0[MAX_SIM_NUM],S[MAX_SIM_NUM],nS[MAX_SIM_NUM];
	double s0[MAX_SIM_NUM],s1[MAX_SIM_NUM];
	double *value;
	double time[MAX_TRACE];
	int num_swap = 0;
	double total_value0,total_value;
	char string[100];
	double dt = 1.0;
	int ns1;	/* number of point for secondary slope */
	double dv;
	int itmp,isort[200];
	int count = 0;
	FILE *tfile;
	char tfile_name[100];

	printf("Type yes to adjust base: ");
	scanf("%s",string);
	if (string[0] == 'y' || string[0] == 'Y')
	    ContinueBaseAdjust(data);
start_trade:

	create_directory("trade",NO);
	sprintf(tfile_name,"trade/record-%d",count);
	tfile = fopen(tfile_name,"w");
	N = data->num_backtrace - 5;
	if (N > MAX_TRACE - 5) N = MAX_TRACE - 5;
	printf("Enter number of point for slope: ");
	scanf("%d",&ns1);
	printf("Enter normalized integral difference for trading: ");
	scanf("%lf",&dv);
	num_swap = 0;

	total_value0 = 0.0;
	is = data->num_days - N;
	for (i = 0; i < M; ++i)
	{
	    printf("For stock %s, enter starting number of shares: ",
				data->assets[i].asset_name);
	    scanf("%lf",S0+i);
	    total_value0 += S0[i]*data->assets[i].value[is];
	}
	printf("\nTotal starting value = %f\n",total_value0);
	fprintf(tfile,"\nTotal starting value = %f\n",total_value0);
	total_value = 0.0;
	is = data->num_days - 1;
	fprintf(tfile,"Starting assets:\n");
	for (i = 0; i < M; ++i)
	{
	    total_value += S0[i]*data->assets[i].value[is];
	    fprintf(tfile,"%5s   %f    %10.2f\n",
				data->assets[i].asset_name,
				S0[i],data->assets[i].value[is]);
	}
	printf("Total untraded end value = %f\n",total_value);
	printf("Change: %f percent\n",100*(total_value/total_value0 - 1.0));
	printf("\n");
	fprintf(tfile,"Total untraded end value = %f\n",total_value);
	fprintf(tfile,"Change: %f percent\n",
				100*(total_value/total_value0 - 1.0));
	fprintf(tfile,"\n");

	for (i = 0; i < N; ++i) 
	    time[i] = i*dt;

	for (j = 0; j < N; ++j)
	{
	    is = data->num_days - N + j - (ns1 - 1);
	    for (i = 0; i < M; ++i)
	    {
	    	value = data->assets[i].norm_value + is;
		s0[i] = value[ns1-1] - value[ns1-2];
		S[i] = data->assets[i].value[is+(ns1-1)];
		nS[i] = data->assets[i].norm_value[is+(ns1-1)];
		LeastSquareLinear(time,value,ns1,&LA,&LB);
		s1[i] = LA;
		isort[i] = i;
	    }
	    for (i1 = 0; i1 < M; ++i1)
	    for (i2 = 0; i2 < M; ++i2)
	    {
		if (nS[isort[i1]] < nS[isort[i2]])
	 	{
		    itmp = isort[i1];
		    isort[i1] = isort[i2];
		    isort[i2] = itmp;
		}
	    }
	    for (i1 = 0; i1 < M; ++i1)
	    for (i2 = 0; i2 < M; ++i2)
	    {
		if (i1 == i2) continue;
		if (
		    S0[i1] != 0.0 
		    && s1[i1] > s1[i2]
		    && 100.0*(nS[i1] - nS[i2]) > dv
		   )
		{
		    S0[i2] += S0[i1]*S[i1]/S[i2];
		    S0[i1] = 0.0;
		    num_swap++;
		}
	    }
	}
	total_value = 0.0;
	is = data->num_days - 1;
	fprintf(tfile,"Ending assets:\n");
	for (i = 0; i < M; ++i)
	{
	    printf("For stock %s, end number of shares: %f\n",
				data->assets[i].asset_name,S0[i]);
	    total_value += S0[i]*data->assets[i].value[is];
	    fprintf(tfile,"%5s   %f    %10.2f\n",
				data->assets[i].asset_name,
				S0[i],data->assets[i].value[is]);
	}
	fprintf(tfile,"Total traded end value = %f\n",total_value);
	fprintf(tfile,"Change: %f percent\n",
			100*(total_value/total_value0 - 1.0));
	fprintf(tfile,"Number of swaps = %d\n",num_swap);
	fclose(tfile);
	printf("Total traded end value = %f\n",total_value);
	printf("Change: %f percent\n",100*(total_value/total_value0 - 1.0));
	printf("Number of swaps = %d\n",num_swap);
	printf("Type yes to restart simulation: ");
	scanf("%s",string);
	if (string[0] == 'y' || string[0] == 'Y')
	{
	    count++;
	    goto start_trade;
	}
}	/* end IntegralSlopeTrade */

#include "swap.h"
#include <ctime>
#include <iostream>
#include <locale>

static void SlopeBasedTrade(MARKET_DATA*);
static void Slope1BasedTrade(MARKET_DATA*);
static void IntegralSlopeTrade(MARKET_DATA*);
static void PyramidTrade(MARKET_DATA*);
static void ComputePyramidParams(double*,double*,int*,int,double,double*,
				double*);
static void GetFiveStates(double*,double*,int*,int,double,double*,double*);
static void PrintDevLinearProfile(char**,double*,double*,double,double,
				int*,int);
static int GetAccountStates(PORTFOLIO*,STATE_INFO*);
static double CurrentBusinessTime();
static double AverageNormPrice(double*,int);
static double MedianNormPrice(double*,int);


extern void InvestSimulation(
	MARKET_DATA *data)
{
	int i_method;
	MARKET_DATA *copy_data = CopyMarketData(data);
	printf("Available trading methods are:\n");
	printf("\tSlope-based trading (0)\n");
	printf("\tLeast-Square-Slope-based trading (1)\n");
	printf("\tIntegral-Slope trading (2)\n");
	printf("\tPyramid trading (3)\n");
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
	case 3:
	    PyramidTrade(copy_data);
	    break;
	default:
	    printf("Unknown trading method\n");
	    return;
	}
	FreeMarketData(copy_data);
}	/* end InvestSimulation */

static void SlopeBasedTrade(
	MARKET_DATA *data)
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
	is = data->num_segs - N;
	for (i = 0; i < M; ++i)
	{
	    printf("For stock %s, enter starting number of shares: ",
				data->assets[i].asset_name);
	    scanf("%lf",S0+i);
	    total_value0 += S0[i]*data->assets[i].value[is];
	}
	printf("\nTotal starting value = %f\n",total_value0);
	total_value = 0.0;
	is = data->num_segs - 1;
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
	    is = data->num_segs - N + j - (ns1 - 1);
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
	is = data->num_segs - 1;
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
	MARKET_DATA *data)
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
	is = data->num_segs - N;
	for (i = 0; i < M; ++i)
	{
	    printf("For stock %s, enter starting number of shares: ",
				data->assets[i].asset_name);
	    scanf("%lf",S0+i);
	    total_value0 += S0[i]*data->assets[i].value[is];
	}
	printf("\nTotal starting value = %f\n",total_value0);
	total_value = 0.0;
	is = data->num_segs - 1;
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
	    is = data->num_segs - N + j - (ns1 - 1);
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
	is = data->num_segs - 1;
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
	MARKET_DATA *data)
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
	is = data->num_segs - N;
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
	is = data->num_segs - 1;
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
	    is = data->num_segs - N + j - (ns1 - 1);
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
	is = data->num_segs - 1;
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

static void PyramidTrade_old(
	MARKET_DATA *data)
{
	int i, j, k, i1, i2;
	int is, isf, N;
	int M = data->num_assets-1;
	
	double des_S0[100], ini_cash[100], old_value[100], new_invest[100];  
	double tmp_os[100], tmp_nS[100], percent[100];
	double linear_value[100];
	double *norm_price[100];
	
	char string[100];
	double total_value0, total_value;
	
	int itmp,isort[200];
	int num_swap = 0;
	double slope,b;		// Slope and intercept of linear profile

	printf("Type yes to adjust base: ");
	scanf("%s", string);
	if (string[0] == 'Y' || string[0] == 'y')
	   ContinueBaseAdjust(data);
	exit(0);

start_trade:
	printf("Enter slope of the profile: ");
	scanf("%lf",&slope);
	N = data->num_backtrace;
	if (N > MAX_TRACE) N = MAX_TRACE;
	is = data->num_segs - N;  /* starting date */
	isf = data->num_segs - 1; /* ending date */

	for (i = 0; i <= M; ++i)
	{
	    isort[i] = i;
	    tmp_os[i] = data->assets[i].value[is];
	    tmp_nS[i] =  data->assets[i].norm_value[is];
	    norm_price[i] = data->assets[i].norm_value;
	}
	for (i1 = 0; i1 < M; ++i1)
	{
	    for (i2 = i1 + 1; i2 <= M; ++i2)
	    {
	        if (tmp_nS[isort[i1]] < tmp_nS[isort[i2]]) 
	        {
		    itmp = isort[i1];
		    isort[i1] = isort[i2];
		    isort[i2] = itmp;
		} 
	    }
	}

	printf("Enter total investment value: ");
	scanf("%lf",&total_value0);
	/* Get each stock num according to investment value */
	for (i = 0; i <= M; ++i)
	{
	    des_S0[isort[i]] = total_value0/tmp_os[isort[i]]/(M+1);  
	}

	/* Calculate the intercept */
	b = total_value0;
	for (i = 0; i <= M; ++i)
	{
	    b -= slope*data->assets[isort[i]].norm_value[is];
	    percent[i] = ini_cash[i]/total_value0; 
	}
	b /= (M + 1);

	printf("\nTotal starting value = %f\n",total_value0);	
	printf("slope = %f  intercept = %f\n",slope,b);
	total_value = 0.0;
	for (i = 0; i <= M; ++i) 
	{
	    total_value += des_S0[isort[i]]*data->assets[isort[i]].value[isf];
	    linear_value[i] = slope*data->assets[i].norm_value[is] + b;
	}
	printf("Sorted starting prices:\n");
	for (i = 0; i <= M; ++i)
	{
	    printf("%4s: %f  %8.2f  %8f\n",data->assets[isort[i]].asset_name,
				data->assets[isort[i]].norm_value[is],
				linear_value[isort[i]],
				des_S0[isort[i]]*tmp_os[isort[i]]);
	}

	printf("Total untraded end value = %f\n",total_value);
	printf("Change: %f percent\n",100 * (total_value/total_value0 - 1.0));
	printf("\n");
	exit(0);


	for (i = is + 1; i <=  isf; i++) 
	{			
	    total_value = 0;    
	    for (j = 0; j <= M; j++)    
	    {
		tmp_nS[isort[j]] = data->assets[isort[j]].norm_value[i];
		tmp_os[isort[j]] = data->assets[isort[j]].value[i];
		old_value[j] = tmp_os[isort[j]] * des_S0[isort[j]];
		total_value += old_value[j];
	    }
	    /* order change due to Norm value */
	    for (j = 0; j < M; j++ )   
	    {
		for (k = j + 1; k <= M; k++)
		{
		    if (tmp_nS[isort[j]] <  tmp_nS[isort[k]])
		    {
		       itmp = isort[j];
		       isort[j] = isort[k];
		       isort[k] = itmp;
		    }
		}
	    }
	    /* each stock number des_S0[i] updates due to triangle holding */
	    for (j = 0; j <= M; j++) 
	    {
		new_invest[j] =  total_value * percent[j];
		if (new_invest[j] > old_value[j])
		{
		   des_S0[isort[j]] += (new_invest[j] - old_value[j])
				/tmp_os[isort[j]];
		   num_swap ++;
		} else if(new_invest[j] < old_value[j])
		{
		   des_S0[isort[j]] -= (old_value[j] - new_invest[j])
				/tmp_os[isort[j]]; 
		   num_swap ++;
		}
	    }
	}   /* end first for */
	total_value = 0;
	/* Get total value in the last day */
	for (i = 0; i <= M ; i++)
	{
	    total_value += des_S0[i] * data->assets[i].value[isf];  
	}
	printf("Total traded end value = %f\n",total_value);
	printf("Change: %f percent\n",100*(total_value/total_value0 - 1.0));
	printf("Number of swaps = %d\n",num_swap);
	printf("Type yes to restart simulation: ");
	scanf("%s",string);
	if (string[0] == 'y' || string[0] == 'Y')
	    goto start_trade;
} 	/* end	PyramidTrade */

extern void PrintCurrentLinearProfile(PORTFOLIO *account)
{
	double *price,*nprice,*bprice;
	int *num_shares;
	char **names;
	double polar_ratio,slope,b;
	int i,id,n,M,N;
	MARKET_DATA *data = account->data;

	N = data->num_assets;
	id = data->num_segs - 1;
	polar_ratio = account->polar_ratio;

	FT_VectorMemoryAlloc((POINTER*)&price,N,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&nprice,N,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&bprice,N,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&num_shares,N,sizeof(int));
	FT_VectorMemoryAlloc((POINTER*)&names,N,sizeof(char*));
	n = 0;
	for (i = 0; i < N; ++i)
	{
	    if (!account->data_map[i])
		continue;
	    bprice[n] = data->assets[i].base_value;
	    num_shares[n] = account->shares[i];
	    names[n] = data->assets[i].asset_name;
	    n++;
	}
	M = n;
	GetCurrentPrice(names,price,M);
	for (i = 0; i < M; ++i)
	    nprice[i] = price[i]/bprice[i];

	ComputePyramidParams(price,nprice,num_shares,M,polar_ratio,&slope,&b);
	PrintDevLinearProfile(names,price,nprice,slope,b,num_shares,M);
	FT_FreeThese(4,names,price,nprice,num_shares);
}	/* end PrintCurrentLinearProfile */

static void PrintDevLinearProfile(
	char **names,
	double *price,
	double *norm_price,
	double slope,		// Slope for A_i = slope*p_i + b
	double b,		// Intercept for A_i = slope*p_i + b
	int *num_shares,
	int M)
{
	int i,j,i1,i2;
	int itmp,isort[100];
	double linear_value[100];
	double current_value[100];
	double market_value = 0.0;
	double ave_nprice;
	int lin_num_shares[100];
	int N_c,N_l;	// Number of normalized shares 

	for (i = 0; i < M; ++i)
	{
	    isort[i] = i;
	}

	for (i1 = 0; i1 < M-1; ++i1)
	{
	    for (i2 = i1 + 1; i2 < M; ++i2)
	    {
	        if (norm_price[isort[i1]] < norm_price[isort[i2]]) 
	        {
		    itmp = isort[i1];
		    isort[i1] = isort[i2];
		    isort[i2] = itmp;
		} 
	    }
	}
	N_c = 0;
	for (i = 0; i < M; ++i)
	{
	    market_value += price[i]*num_shares[i];
	    N_c += irint(price[i]*num_shares[i]/norm_price[i]);
	}

	market_value = 0.0;
	N_l = 0;
	for (i = 0; i < M; ++i)
	{
	    linear_value[i] = slope*norm_price[i] + b;
	    lin_num_shares[i] = linear_value[i]/price[i];
	    current_value[i] = num_shares[i]*price[i];
	    N_l += irint(linear_value[i]/norm_price[i]);
	}
	ave_nprice = MedianNormPrice(norm_price,M);
	printf("ave_nprice = %f\n",ave_nprice);
	printf("Equity NPrice   L-Val    C-Val    D-Val  "
		"L-Share  C-Share  D-Share\n");
	for (i = 0; i < M; ++i)
	{
	    i1 = isort[i];
	    printf("  %4s  %5.3f  %6.0f   %6.0f   %6.0f   %6d   %6d   %6d\n",
				names[i1],norm_price[i1],
				linear_value[i1],current_value[i1],
				current_value[i1]-linear_value[i1],
				lin_num_shares[i1],num_shares[i1],
				num_shares[i1]-lin_num_shares[i1]);
			
	}
	printf("Average normalized price: %5.3f\n\n",ave_nprice);
}	/* end PrintDevLinearProfile */


extern void PrintDataStates(PORTFOLIO *account)
{
	int i,N;
	char string[100];
	double business_time;
	MARKET_DATA *data = account->data;
	int ie = account->eindex;
	int id = data->num_segs-1;
	int Ne;
	STATE_INFO sinfo;
	int *S = sinfo.S;

	N = data->num_assets;

	Ne = GetAccountStates(account,&sinfo);

	printf("  Current   State-0   State-1   State-2   State-3   State-4\n");
	printf("%9d %9d %9d %9d %9d %9d\n",S[5],S[0],S[1],S[2],S[3],S[4]);
	printf("  State-C %9d %9d %9d %9d %9d\n",S[0]-S[5],S[1]-S[5],S[2]-S[5],
					S[3]-S[5],S[4]-S[5]);
	printf("\n");
	printf("Equivalent to %d %s shares\n",Ne,data->assets[ie].asset_name);
	printf("\n");
	printf("Type yes to record states: ");
	scanf("%s",string);
	if (string[0] == 'y' || string[0] == 'Y')
	{
	    char fname[256];
	    FILE *sfile;
	    const char scolor[][100] = {"yellow","green","pink","blue",
					"red","aqua"};

	    business_time = CurrentBusinessTime();
	    create_directory("state",NO);
	    for (i = 0; i < 6; ++i)
	    {
		sprintf(fname,"state/%s-state-%d",account->account_name,i);
		sfile = fopen(fname,"r");
		if (sfile != NULL)
		{
		    fclose(sfile);
            	    sfile = fopen(fname,"a");
		}
		else
		{
            	    sfile = fopen(fname,"w");
		    fprintf(sfile,"Next\n");
		    fprintf(sfile,"color=%s\n",scolor[i]);
		    fprintf(sfile,"thickness = 1.5\n");
		}
		fprintf(sfile,"%8.2f  %9d\n",business_time,S[i]);
		fclose(sfile);
	    }
	    sprintf(fname,"Eshare/%s-%s",account->account_name,
				data->assets[ie].asset_name);
	    create_directory("Eshare",NO);
	    sfile = fopen(fname,"r");
	    if (sfile != NULL)
	    {
		fclose(sfile);
            	sfile = fopen(fname,"a");
	    }
	    else
	    {
            	sfile = fopen(fname,"w");
		fprintf(sfile,"Next\n");
		fprintf(sfile,"color=%s\n",data->assets[ie].color);
		fprintf(sfile,"thickness = 1.5\n");
	    }
	    fprintf(sfile,"%8.2f  %9d\n",business_time,Ne);
	}
}	/* end PrintDataStates */

static int GetAccountStates(
	PORTFOLIO *account,
	STATE_INFO *sinfo)
{
	double *price,*nprice,*bprice;
	double np_min,np_max,np_ave;
	double slope,b,pratio;
	char **names;
	char string[100];
	int i,j,id,n,M,N,idmin,idmax;
	int *num_shares;
	int V_total,C;
	int eindex,Ne,ne;
	int *isort,*ns_L,*ns_U;
	MARKET_DATA *data = account->data;

	N = data->num_assets;
	id = data->num_segs - 1;
	pratio = account->polar_ratio;
	eindex = account->eindex;

	FT_VectorMemoryAlloc((POINTER*)&price,N,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&nprice,N,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&bprice,N,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&num_shares,N,sizeof(int));
	FT_VectorMemoryAlloc((POINTER*)&names,N,sizeof(char*));
	FT_VectorMemoryAlloc((POINTER*)&ns_L,N,sizeof(int));
	FT_VectorMemoryAlloc((POINTER*)&ns_U,N,sizeof(int));
	FT_VectorMemoryAlloc((POINTER*)&isort,N,sizeof(int));
	n = 0;
	np_max = -HUGE;
	np_min = HUGE;
	C = V_total = 0.0;
	ne = -1;
	for (i = 0; i < N; ++i)
	{
	    if (!account->data_map[i])
		continue;
	    bprice[n] = data->assets[i].base_value;
	    num_shares[n] = account->shares[i];
	    names[n] = data->assets[i].asset_name;
	    if (eindex == i) ne = n;
	    n++;
	}
	M = n;

	GetCurrentPrice(names,price,M);
	for (i = 0; i < M; ++i)
	{
	    isort[i] = i;
	    nprice[i] = price[i]/bprice[i];
	    if (np_min > nprice[i]) 
		np_min = nprice[i];
	    if (np_max < nprice[i]) 
		np_max = nprice[i];
	    V_total += price[i]*num_shares[i];
	    C += irint(price[i]*num_shares[i]/nprice[i]);
	}

	np_ave = MedianNormPrice(nprice,M);
	for (i = 0; i < M-1; ++i)
        for (j = i+1; j < M; ++j)
        {
	    if (nprice[isort[i]] < nprice[isort[j]])
	    {
		int itmp = isort[j];
                isort[j] = isort[i];
                isort[i] = itmp;
	    }
	}
	strcpy(sinfo->stock_max,names[isort[0]]);
	strcpy(sinfo->stock_min,names[isort[M-1]]);
	sinfo->dnp_max = (nprice[isort[0]] - nprice[isort[1]])/np_ave;
	sinfo->dnp_min = (nprice[isort[M-2]] - nprice[isort[M-1]])/np_ave;

	sinfo->S[0] = irint(V_total/np_max);
	sinfo->S[2] = irint(V_total/np_ave);
	sinfo->S[4] = irint(V_total/np_min);
	ComputePyramidParams(price,nprice,num_shares,M,1/pratio,&slope,&b);
	sinfo->S[1] = sinfo->S[3] = 0.0;
	for (i = 0; i < M; ++i)
	{
	    ns_U[i] = irint((slope*nprice[i] + b)/nprice[i]);
	    sinfo->S[1] += ns_U[i];
	}
	ComputePyramidParams(price,nprice,num_shares,M,pratio,&slope,&b);
	sinfo->svalue = sinfo->dvalue = 0.0;
	for (i = 0; i < M; ++i)
	{
	    double Vl,Vc;
	    Vc = num_shares[i]*price[i];
	    Vl = slope*nprice[i] + b;
	    ns_L[i] = irint(Vl/nprice[i]);
	    sinfo->S[3] += ns_L[i];
	    if (Vc - Vl > sinfo->svalue)
	    {
		idmax = i;
		sinfo->svalue = Vc - Vl;
	    }
	    if (Vc - Vl < sinfo->dvalue)
	    {
		idmin = i;
		sinfo->dvalue = Vc - Vl;
	    }
	}
	strncpy(sinfo->sname,names[idmax],20);
	strncpy(sinfo->dname,names[idmin],20);
	sinfo->dnp = (nprice[idmax] - nprice[idmin])/np_ave;
	sinfo->C = C;
	if (0 <= ne && ne < N)
	    Ne = irint(V_total/price[ne]);
	else
	    Ne = 0;
	FT_FreeThese(4,names,price,nprice,num_shares);
	return Ne;
}	/* end GetAccountStates */

static void ComputePyramidParams(
	double *price,
	double *nprice,
	int *shares,
	int M,
	double m,
	double *slope,
	double *intercept)
{
	double a[2],b[2],c[2],x[2];
	int i,N;
	double pmin,pmax;
	int imin,imax;

	a[0] = 0.0;
	c[0] = c[1] = 0.0;
	pmin = HUGE;
	pmax = -HUGE;
	N = 0;
	for (i = 0; i < M; ++i)
	{
	    if (nprice[i] < pmin)
	    {
		pmin = nprice[i];
		imin = i;
	    }
	    if (nprice[i] > pmax)
	    {
		pmax = nprice[i];
		imax = i;
	    }
	    a[0] += nprice[i];
	    c[0] += price[i]*shares[i];
	    N++;
	}
	a[1] = nprice[imin] - m*nprice[imax];
	b[0] = (double)N;	
	b[1] = (1.0 - m);
	linear_22_equation(a,b,c,x);
	*slope = x[0];
	*intercept = x[1];
}	/* end ComputePyramidParams */

static double CurrentBusinessTime()
{
	std::time_t t = std::time(NULL);
	char string[100];
	double year,week,day,hour,minute;
	double business_time;

	strftime(string,100*sizeof(char),"%G",std::localtime(&t));
	year = (double)atoi(string);
	strftime(string,100*sizeof(char),"%U",std::localtime(&t));
	week = (double)atoi(string);
	strftime(string,100*sizeof(char),"%u",std::localtime(&t));
	day = (double)atoi(string);
	strftime(string,100*sizeof(char),"%H",std::localtime(&t));
	hour = (double)atoi(string);
	strftime(string,100*sizeof(char),"%M",std::localtime(&t));
	minute = (double)atoi(string);
	business_time = week*5.0 + day*8.0 + hour + minute/60.0 - 9.0;
	
	return business_time;
}	/* end CurrentBusinessTime */

static double AverageNormPrice(
	double *nprice,
	int M)
{
	int i,i1,i2;
	int itmp,isort[100];
	double ave_nprice = 0.0;

	for (i = 0; i < M; ++i)
	    isort[i] = i;

	for (i1 = 0; i1 < M-1; ++i1)
	{
	    for (i2 = i1 + 1; i2 < M; ++i2)
	    {
	        if (nprice[isort[i1]] < nprice[isort[i2]]) 
	        {
		    itmp = isort[i1];
		    isort[i1] = isort[i2];
		    isort[i2] = itmp;
		} 
	    }
	}
	for (i = 0; i < M-1; ++i)
	{
	    i1 = isort[i];
	    i2 = isort[i+1];
	    ave_nprice += 0.5*(nprice[i1] + nprice[i2])*
			fabs(nprice[i1] - nprice[i2]);
	}
	i1 = isort[0];
	i2 = isort[M-1];
	ave_nprice /= fabs(nprice[i1] - nprice[i2]);
	return ave_nprice;
}	/* end AverageNormPrice */

extern void ReportDataStates(
	PORTFOLIO **accounts,
	int num_accounts)
{
	int i,ie,N,Ne;
	int seconds;
	MARKET_DATA *data = accounts[0]->data;
	STATE_INFO sinfo;
	int *S = sinfo.S;
	std::time_t t;
	char date[100],hour[10],minute[10],string[100];

	N = data->num_assets;

	printf("Enter number of minutes for reporting: ");
	scanf("%d",&seconds);
	seconds *= 60;	
	while (YES)
	{
 	    t = std::time(NULL);
	    strftime(date,100*sizeof(char),"%D",std::localtime(&t));
	    printf("\n");
	    printf("========================%s==========",date);
	    printf("=================\n");
	    for (i = 0; i < num_accounts; ++i)
	    {
 	    	t = std::time(NULL);
	    	strftime(hour,100*sizeof(char),"%H",std::localtime(&t));
	    	strftime(minute,100*sizeof(char),"%M",std::localtime(&t));
	    	Ne = GetAccountStates(accounts[i],&sinfo);
		ie = accounts[i]->eindex;

		printf("Account %s at time %s:%s \n",accounts[i]->account_name,
				hour,minute);
	    	printf("  Current   State-0   State-1   State-2"
		   "   State-3   State-4\n");
	    	printf("%9d %9d %9d %9d %9d %9d\n",S[5],S[0],S[1],S[2],S[3],
					S[4]);
	    	printf("  State-C %9d %9d %9d %9d %9d\n",S[0]-S[5],S[1]-S[5],
					S[2]-S[5],S[3]-S[5],S[4]-S[5]);
	    	printf("Equivalent to %d %s shares\n",Ne,
					data->assets[ie].asset_name);
		printf("Top surplus (%4s): %8.2f | ",sinfo.sname,sinfo.svalue);
		printf("Top deficit (%4s): %8.2f\n",sinfo.dname,sinfo.dvalue);
		printf("Percetage diff between surplus and deficit: %5.2f\n",
				sinfo.dnp*100);
		printf("   Top stock %4s with gap to main body:%5.2f percent\n",
				sinfo.stock_max,sinfo.dnp_max*100);
		printf("Bottom stock %4s with gap to main body:%5.2f percent\n",
				sinfo.stock_min,sinfo.dnp_min*100);
	    	printf("\n");
	    }
	    sleep(seconds);
	}
}	/* end ReportDataStates */

static void GetFiveStates(
	double *price,
	double *nprice,
	int *shares,
	int M,
	double pratio,
	double *S,
	double *sL)
{
	double slope,b;
	int i;
	double V,V_total;
	double np_min,np_max,np_mid;
	double ns_L,ns_U;

	S[5] = V_total = 0.0;
	np_min = HUGE;
	np_max = -HUGE;
	for (i = 0; i < M; ++i)
	{
	    V_total += shares[i]*price[i];
	    S[5] += shares[i]*price[i]/nprice[i];
	    if (np_min > nprice[i]) np_min = nprice[i];
	    if (np_max < nprice[i]) np_max = nprice[i];
	}
	np_mid = MedianNormPrice(nprice,M);
	S[0] = irint(V_total/np_max);
	S[2] = irint(V_total/np_mid);
	S[4] = irint(V_total/np_min);
	ComputePyramidParams(price,nprice,shares,M,1/pratio,&slope,&b);
	S[1] = S[3] = 0.0;
	for (i = 0; i < M; ++i)
	{
	    V = slope*nprice[i] + b;
	    ns_U = irint(V/nprice[i]);
	    S[1] += ns_U;
	}
	ComputePyramidParams(price,nprice,shares,M,pratio,&slope,&b);
	for (i = 0; i < M; ++i)
	{
	    V = slope*nprice[i] + b;
	    ns_L = irint(V/nprice[i]);
	    sL[i] = V;
	    S[3] += ns_L;
	}
}	/* GetFiveStates */

static void PyramidTrade(
	MARKET_DATA *data)
{
	char string[100];
	char fname[256],dirname[256];
	int i,j,n,M = data->num_assets;
	int istart,iend,N;
	int *shares;
	double V0,Vf;		   // initial and final value
	double *price,*nprice,*vL; // price, normal price, linear value profile
	double S[6];		   // five states, S[5] is current state
	double pratio;		   // profile polarization ratio
	const char scolor[][100] = {"yellow","green","pink","blue",
                                        "red","aqua"};
	boolean *data_map;
	double **value,**norm_value;
	FILE *sfile[6];
	int imax,imin;
	double np_max,np_min;
	int num_trades;
	double lambda;
	
	FT_VectorMemoryAlloc((POINTER*)&shares,M,sizeof(int));
	FT_VectorMemoryAlloc((POINTER*)&price,M,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&nprice,M,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&vL,M,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&data_map,M,sizeof(boolean));
	FT_VectorMemoryAlloc((POINTER*)&value,M,sizeof(double*));
	FT_VectorMemoryAlloc((POINTER*)&norm_value,M,sizeof(double*));
	for (i = 0; i < data->num_assets; ++i)
	{
	    data_map[i] = NO;
	    printf("Type yes to include %4s: ",data->assets[i].asset_name);
	    scanf("%s",string);
	    if (string[0] == 'y') data_map[i] = YES;
	}

start_trade:
	ContinueBaseAdjust(data);
	printf("Enter polarization ratio of the profile: ");
	scanf("%lf",&pratio);
	printf("Enter exchange level parameter: ");
	scanf("%lf",&lambda);


	iend = data->num_segs - 1; /* ending segment */
	printf("Enter total trading segments (back from current): ");
	scanf("%d",&N);
	if (iend - N < data->num_backtrace) 
	    istart = data->num_backtrace + 1;
	else
	    istart = data->num_segs - N;  /* starting segment */

	printf("Enter total investment fund: ");
	scanf("%lf",&V0);
	np_max = -HUGE;
	np_min = HUGE;
	n = 0;
	for (i = 0; i < data->num_assets; ++i)
	{
	    if (!data_map[i]) continue;
	    value[n] = data->assets[i].value;
	    norm_value[n] = data->assets[i].norm_value;
	    price[n] = value[n][istart];
	    nprice[n] = norm_value[n][istart];
	    if (nprice[n] > np_max)
	    {
		np_max = nprice[n];
		imax = n;
	    }
	    if (nprice[n] < np_min)
	    {
		np_min = nprice[n];
		imin = n;
	    }
	    n++;
	}
	M = n;
	printf("Choose initial portfolio profile: ");
	printf("Avaliable profiles are:\n");
	printf("\tAll in lowest normal share price (0)\n");
	printf("\tLower linear profile (1)\n");
	printf("\tUniform profile (2)\n");
	printf("\tUpper linear profile (3)\n");
	printf("\tAll in highest normal share price (4)\n");
	printf("Enter choice: ");
	scanf("%s",string);
	for (i = 0; i < M; ++i)
	    shares[i] = 0;
	switch(string[0])
	{
	case '0':
	    shares[imin] = V0/price[imin];
	    break;
	case '1':
	    shares[imin] = V0/price[imin];
	    GetFiveStates(price,nprice,shares,M,pratio,S,vL);
	    S[5] = 0;
	    for (i = 0; i < M; ++i)
	    {
	    	shares[i] = vL[i]/price[i];
		S[5] += irint(vL[i]/nprice[i]);
	    }
	    break;
	case '2':
	    for (i = 0; i < M; ++i)
	    	shares[i] = V0/M/price[i];
	    break;
	case '3':
	    shares[imin] = V0/price[imin];
	    GetFiveStates(price,nprice,shares,M,1.0/pratio,S,vL);
	    S[5] = 0;
	    for (i = 0; i < M; ++i)
	    {
	    	shares[i] = vL[i]/price[i];
		S[5] += irint(vL[i]/nprice[i]);
	    }
	    break;
	case '4':
	    shares[imax] = V0/price[imax];
	    break;
	default:
	    printf("Unknow choice!\n");
	    goto start_trade;
	}

	Vf = 0.0;
	for (i = 0; i < M; ++i)
	{
	    Vf += shares[i]*value[i][iend];
	}

	printf("Total untraded end value = %f\n",Vf);
	printf("Change: %f percent\n",100 * (Vf/V0 - 1.0));
	printf("\n");

	create_directory("simulation/state",NO);
	for (i = 0; i < 6; i++)    
	{
	    sprintf(fname,"simulation/state/s%d",i);
	    sfile[i] = fopen(fname,"w");
	    fprintf(sfile[i],"Next\n");
            fprintf(sfile[i],"color=%s\n",scolor[i]);
            fprintf(sfile[i],"thickness = 1.5\n");
	}
	num_trades = 0;
	for (i = istart + 1; i <=  iend; i++) 
	{			
	    for (j = 0; j < M; j++)    
	    {
		price[j] = value[j][i];
		nprice[j] = norm_value[j][i];
	    }
	    GetFiveStates(price,nprice,shares,M,pratio,S,vL);
	    if (S[5] < S[3] - lambda*(S[3] - S[1]))
	    {
		num_trades++;
	    	for (j = 0; j < M; j++)    
		    shares[j] = vL[j]/price[j];
	    }
	    for (j = 0; j < 6; j++)    
		fprintf(sfile[j],"%4d  %f\n",i-istart,S[j]);
	}
	for (i = 0; i < 6; i++)    
	    fclose(sfile[i]);
	Vf = 0.0;
	for (i = 0; i < M ; i++)
	{
	    Vf += shares[i]*value[i][iend];  
	}
	printf("Total traded end value = %f\n",Vf);
	printf("Change: %f percent\n",100*(Vf/V0 - 1.0));
	printf("Number of trades = %d\n",num_trades);
	printf("Type yes to restart simulation: ");
	scanf("%s",string);
	if (string[0] == 'y' || string[0] == 'Y')
	    goto start_trade;
} 	/* end	PyramidTrade */

static double MedianNormPrice(
	double *nprice,
	int M)
{
	int i,i1,i2;
	int itmp,isort[100];
	double mid_nprice;

	for (i = 0; i < M; ++i)
	    isort[i] = i;

	for (i1 = 0; i1 < M-1; ++i1)
	{
	    for (i2 = i1 + 1; i2 < M; ++i2)
	    {
	        if (nprice[isort[i1]] < nprice[isort[i2]]) 
	        {
		    itmp = isort[i1];
		    isort[i1] = isort[i2];
		    isort[i2] = itmp;
		} 
	    }
	}
	if (M%2 == 1)
	{
	    i = (M - 1)/2;
	    mid_nprice = nprice[isort[i]];
	}
	else
	{
	    i1 = M/2 - 1;
	    i2 = M/2;
	    mid_nprice = 0.5*(nprice[isort[i1]] + nprice[isort[i2]]);
	}
	return mid_nprice;
}	/* end MedianNormPrice */


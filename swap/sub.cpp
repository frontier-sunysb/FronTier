#include "swap.h"

static void AddAssets(DATA_SET*);
static void DeleteAssets(DATA_SET*);
static void PrintRateOfChange(DATA_SET*);
static double AssetTimeIntegral(DATA_SET*,int,int,int);

extern double GetLeastSquare(
	DATA_SET *data,
	int m,
	int n,
	double *Q,
	double *L)
{
	int i,is,N;
	double *value;
	double *time;
	double QA,QB,QC;
	double LA,LB;
	FILE *xfile;
	char fname[100];
	int id = (n-1)/4;
	double base_value = data->assets[m].base_value;

	is = (data->num_days < n) ? 0 : data->num_days - n;
	N = data->num_days - is;
	value = data->assets[m].norm_value + is;
	FT_VectorMemoryAlloc((POINTER*)&time,N,sizeof(double));
	for (i = 0; i < N; ++i) time[i] = i*0.25;
	LeastSquareLinear(time,value,N,&LA,&LB);
	sprintf(fname,"xg/Fit-l-%d-%d.xg",m,id);
	xfile = fopen(fname,"w");
	fprintf(xfile,"color=green\n");
	fprintf(xfile,"thickness = 1.5\n");
	for (i = 0; i < N; ++i)
	{
	    fprintf(xfile,"%f %f\n",time[i],value[i]);
	}
	fprintf(xfile,"Next\n");
	fprintf(xfile,"color=red\n");
	fprintf(xfile,"thickness = 1.5\n");
	for (i = 0; i < N; ++i)
	{
	    fprintf(xfile,"%f %f\n",time[i],LA*time[i]+LB);
	}
	fclose(xfile);
	LeastSquareQuadr(time,value,N,&QA,&QB,&QC);
	sprintf(fname,"xg/Fit-q-%d-%d.xg",m,id);
	xfile = fopen(fname,"w");
	fprintf(xfile,"color=green\n");
	fprintf(xfile,"thickness = 1.5\n");
	for (i = 0; i < N; ++i)
	{
	    fprintf(xfile,"%f %f\n",time[i],value[i]);
	}
	fprintf(xfile,"Next\n");
	fprintf(xfile,"color=red\n");
	fprintf(xfile,"thickness = 1.5\n");
	for (i = 0; i < N; ++i)
	{
	    fprintf(xfile,"%f %f\n",time[i],QA*sqr(time[i])+QB*time[i]+QC);
	}
	fclose(xfile);
	FT_FreeThese(1,time);
	Q[0] = QA;	Q[1] = QB;	Q[2] = QC;
	L[0] = LA;	L[1] = LB;
}	/* end GetLeastSquare */

extern double LeastSquareQuadr(
	double *x,
	double *y,
	int N,
	double *A,
	double *B,
	double *C)
{
	double x0,x1,x2,x3,x4,yx0,yx1,yx2;
	double sx0,sx1,sx2,sx3,sx4,syx0,syx1,syx2;
	double a[3],b[3],c[3],d[3];
	double dinom;
	int i;

	x0 = x1 = x2 = x3 = x4 = yx0 = yx1 = yx2 = 0.0;
	sx0 = sx1 = sx2 = sx3 = sx4 = syx0 = syx1 = syx2 = 0.0;
	for (i = 0; i < N; ++i)
	{
	    x0 = 1.0;
	    x1 = x0*x[i];
	    x2 = x1*x[i];
	    x3 = x2*x[i];
	    x4 = x3*x[i];
	    yx0 = y[i];
	    yx1 = y[i]*x1;
	    yx2 = y[i]*x2;
	    sx0 += x0;
	    sx1 += x1;
	    sx2 += x2;
	    sx3 += x3;
	    sx4 += x4;
	    syx0 += yx0;
	    syx1 += yx1;
	    syx2 += yx2;
	}
	a[0] =               sx4;
	a[1] = b[0] =        sx3;
	a[2] = b[1] = c[0] = sx2;
	       b[2] = c[1] = sx1;
		      c[2] = sx0;
	d[0] = syx2;
	d[1] = syx1;
	d[2] = syx0;
	dinom = Det3d(a,b,c);
	*A = Det3d(d,b,c)/dinom;
	*B = Det3d(a,d,c)/dinom;
	*C = Det3d(a,b,d)/dinom;
}	/* end LeastSquareQuadr */

#define	Det2d(a,b) ((a)[0]*(b)[1] - (a)[1]*(b)[0]) 

extern double LeastSquareLinear(
	double *x,
	double *y,
	int N,
	double *A,
	double *B)
{
	double x0,x1,x2,yx0,yx1;
	double sx0,sx1,sx2,syx0,syx1;
	double a[2],b[2],c[2];
	double dinom;
	int i;

	x0 = x1 = x2 = yx0 = yx1 = 0.0;
	sx0 = sx1 = sx2 = syx0 = syx1 = 0.0;
	for (i = 0; i < N; ++i)
	{
	    x0 = 1.0;
	    x1 = x0*x[i];
	    x2 = x1*x[i];
	    yx0 = y[i];
	    yx1 = y[i]*x1;
	    sx0 += x0;
	    sx1 += x1;
	    sx2 += x2;
	    syx0 += yx0;
	    syx1 += yx1;
	}
	a[0] =        sx2;
	a[1] = b[0] = sx1;
	       b[1] = sx0;
	c[0] = syx1;
	c[1] = syx0;
	dinom = Det2d(a,b);
	*A = Det2d(c,b)/dinom;
	*B = Det2d(a,c)/dinom;
}	/* end LeastSquareLinear */

extern void ContinueBaseAdjust(
        DATA_SET *data)
{
	int time_p;    /* time period */
	int i,j,k,M = data->num_assets;
	double ave_value;
	int save_num_days = data->num_days;
	int backtrace;
	char fname[256],dirname[256];
	FILE **bfiles;
	double *base0;

	FT_VectorMemoryAlloc((POINTER*)&bfiles,M,sizeof(FILE*));
	FT_VectorMemoryAlloc((POINTER*)&base0,M,sizeof(double));

	printf("There are total %d days of data\n",data->num_days);
	printf("Enter backtrace for base update: ");
	scanf("%d", &backtrace);
	data->num_backtrace = backtrace;
	sprintf(dirname,"simulation/%s",data->account_name);
	create_directory(dirname,NO);
	for (i = 0; i < M; ++i)
	{
	    sprintf(fname,"%s/%s",dirname,data->assets[i].asset_name);
	    bfiles[i] = fopen(fname,"w");
	    fprintf(bfiles[i],"Next\n");
	    fprintf(bfiles[i],"color=%s\n",data->assets[i].color);
	    fprintf(bfiles[i],"thickness = 1.5\n");
	}

	for (j = 10; j <= save_num_days; j++)
	{
	    data->num_days = j;
	    AdjustBase(data);
	    for (i = 0; i < data->num_assets; i++)
	    {
		if (j == 10)
		    base0[i] = data->assets[i].base_value;
		k = data->num_days-1;
		data->assets[i].norm_value[k] =
			data->assets[i].value[k]/data->assets[i].base_value;
		fprintf(bfiles[i],"%f  %f\n",(double)j,
			data->assets[i].base_value/base0[i]);
	    }
	}
	for (i = 0; i < data->num_assets; i++)
	    fclose(bfiles[i]);
	FT_FreeThese(2,bfiles,base0);
	data->num_days = save_num_days;
}	/* end ContinueBaseAdjust */

extern void AdjustBase(
	DATA_SET *data)
{
	int i;
	int is,N;
	double I;

	is = (data->num_days < data->num_backtrace) ? 0 :
			data->num_days - data->num_backtrace;
	N = data->num_days - is;

	for (i = 0; i < data->num_assets; ++i)	
	{
	    I = AssetTimeIntegral(data,i,is,data->num_days);
	    data->assets[i].base_value = I/N;
	}
}	/* end AdjustBase */

extern void ModifyData(
	DATA_SET *data)
{
	int i,j,m,n;
	int istar;		/* center asset index */
	char string[100];
	double ave_value;
	double mfactor;

	printf("Type yes to change data set: ");
	scanf("%s",string);
	if (string[0] != 'y' && string[0] != 'Y')
	    return;

	printf("Type yes to add assets: ");
	scanf("%s",string);
	if (string[0] == 'y' || string[0] == 'Y')
	{
	    AddAssets(data);
	}

	printf("Type yes to delete assets: ");
	scanf("%s",string);
	if (string[0] == 'y' || string[0] == 'Y')
	{
	    DeleteAssets(data);
	}

	printf("Type yes to change number of back trace: ");
	scanf("%s",string);
	if (string[0] == 'y' || string[0] == 'Y')
	{
	    printf("Current back trace: %d\n",data->num_backtrace);
	    printf("Enter number of back trace: ");
	    scanf("%d",&data->num_backtrace);
	}

	printf("Type yes to change account profile slope: ");
	scanf("%s",string);
	if (string[0] == 'y' || string[0] == 'Y')
	{
	    printf("Current profile slope: %f\n",data->polar_ratio);
	    printf("Enter account profile slope: ");
	    scanf("%lf",&data->polar_ratio);
	}
   	PromptForDataMap(data);

	printf("Type yes to modify data base: ");
	scanf("%s",string);
	if (string[0] != 'y' && string[0] != 'Y')
	    return;
	AdjustBase(data);
}	/* end ModifyData */

static void AddAssets(
	DATA_SET *data)
{
	int i,j,num_add;
	ASSET *new_assets;
	int new_num_assets;

	printf("Current number of assets: %d\n",data->num_assets);
	printf("Enter number of assets to be added: ");
	scanf("%d",&num_add);

	new_num_assets = data->num_assets + num_add;
	FT_VectorMemoryAlloc((POINTER*)&new_assets,new_num_assets+1,
                                sizeof(ASSET));
	for (j = 0; j < new_num_assets+1; ++j)
        {
	    FT_VectorMemoryAlloc((POINTER*)&new_assets[j].value,
                                data->num_days+10,sizeof(double));
            FT_VectorMemoryAlloc((POINTER*)&new_assets[j].norm_value,
                                data->num_days+10,sizeof(double));
	    if (j < data->num_assets)
	    {
		new_assets[j].base_value = data->assets[j].base_value;
		strcpy(new_assets[j].asset_name,data->assets[j].asset_name);
		strcpy(new_assets[j].color,data->assets[j].color);
		for (i = 0; i < data->num_days; ++i)
		{
		    new_assets[j].value[i] = data->assets[j].value[i];
		    new_assets[j].norm_value[i] = data->assets[j].norm_value[i];
		}
	    }
	    else if (j < new_num_assets)
	    {
		printf("Enter name of the new asset: ");
		scanf("%s",new_assets[j].asset_name);
		printf("Enter color of the new asset: ");
		scanf("%s",new_assets[j].color);
		printf("Enter base value of the new asset: ");
		scanf("%lf",&new_assets[j].base_value);
		for (i = 0; i < data->num_days; ++i)
		{
		    new_assets[j].value[i] = new_assets[j].base_value;
		    new_assets[j].norm_value[i] = 1.0;
		}
	    }
	    else
	    {
		int im = data->num_assets;
		new_assets[j].base_value = data->assets[im].base_value;
		strcpy(new_assets[j].asset_name,data->assets[im].asset_name);
		strcpy(new_assets[j].color,data->assets[im].color);
		for (i = 0; i < data->num_days; ++i)
		{
		    new_assets[j].value[i] = data->assets[im].value[i];
		    new_assets[j].norm_value[i] = 
				data->assets[im].norm_value[i];
		}
	    }
	}

	data->assets = new_assets;
	data->num_assets = new_num_assets;
	FT_FreeThese(1,data->assets);
}	/* end AddAssets */

static void DeleteAssets(
        DATA_SET *data)
{
}       /* end DeleteAssets */

extern void PromptForDataMap(
	DATA_SET *data)
{
	int i;
	char string[100];

	printf("Type yes to change data map: ");
	scanf("%s",string);
	if (string[0] != 'y' && string[0] != 'Y') return;

	for (i = 0; i < data->num_assets; ++i)
	{
	    printf("Include asset %s ? ",data->assets[i].asset_name);
	    scanf("%s",string);
	    if (string[0] == 'y' || string[0] == 'Y')
	    	data->data_map[i] = YES;
	    else
	    	data->data_map[i] = NO;
	}
}	/* end PromptForDataMap */

extern void AddData(
	DATA_SET *data)
{
	int j,m,M;
	char string[100];
	double ave_value;
	char **stock_names;
	double *current_price;

	m = data->num_days;
	M = data->num_assets;

	FT_VectorMemoryAlloc((POINTER*)&stock_names,M,sizeof(char*));
	FT_VectorMemoryAlloc((POINTER*)&current_price,M,sizeof(double));
	for (j = 0; j < M; ++j)
	    stock_names[j] = data->assets[j].asset_name;
	GetCurrentPrice(stock_names,current_price,M);
	printf("Current prices:\n");
	for (j = 0; j < M; ++j)
	    printf("\t%5s: %f\n",stock_names[j],current_price[j]);
	printf("Type yes to add to database: ");
	scanf("%s",string);
	if (string[0] != 'y' && string[0] != 'Y')
	{
	    FT_FreeThese(2,stock_names,current_price);
	    return;
	}
	data->new_data = YES;

	ave_value = 0.0;
	for (j = 0; j < M; ++j)
	{
	    data->assets[j].value[m] = current_price[j];
	    if (data->assets[j].value[m] == -1)
		data->assets[j].value[m] = data->assets[j].value[m-1];
	    data->assets[j].norm_value[m] = data->assets[j].value[m]/
					data->assets[j].base_value;
	    ave_value += data->assets[j].norm_value[m];
	}
	data->assets[j].base_value = 1.0;
	data->assets[j].value[m] = ave_value/data->num_assets;
	data->assets[j].norm_value[m] = data->assets[j].value[m];
	data->num_days += 1;
	FT_FreeThese(2,stock_names,current_price);
}	/* end AddData */

extern void XgraphData(
	DATA_SET *data)
{
	char fname[100];
	FILE *xfile;
	int i,j,k,M,is;
	double time;
	double *value,ave_value;

	create_directory("xg",NO);
	M = data->num_assets;
	is = (data->num_days < data->num_backtrace) ? 0 :
		data->num_days - data->num_backtrace;
	for (j = 0; j < M; ++j)
	{
	    if (data->data_map[j] != YES)
	    	continue;
	    sprintf(fname,"xg/N-%d-.xg",j);
	    xfile = fopen(fname,"w");
	    fprintf(xfile,"color=%s\n",data->assets[j].color);
	    fprintf(xfile,"thickness = 1.5\n");
	    value = data->assets[j].norm_value;
	    for (i = is; i < data->num_days; ++i)
	    {
		if (value[i] == 0.0) continue;
		time = (double)i;
		fprintf(xfile,"%f %f\n",time,value[i]);
	    }
	    fclose(xfile);
	}
	for (j = 0; j < M-1; ++j)
	for (k = j+1; k < M; ++k)
	{
	    if (data->data_map[j] != YES || data->data_map[k] != YES)
	    	continue;
	    sprintf(fname,"xg/M-%d-%d.xg",j,k);
	    xfile = fopen(fname,"w");
	    fprintf(xfile,"color=green\n");
	    fprintf(xfile,"thickness = 1.5\n");
	    for (i = is; i < data->num_days; ++i)
	    {
		if (value[i] == 0.0) continue;
		time = (double)i;
		ave_value = 0.5*(data->assets[j].norm_value[i] +
				 data->assets[k].norm_value[i]);
		fprintf(xfile,"%f %f\n",time,ave_value);
	    }
	    fclose(xfile);
	}
	for (j = 0; j < M-1; ++j)
	for (k = j+1; k < M; ++k)
	{
	    if (data->data_map[j] != YES || data->data_map[k] != YES)
	    	continue;
	    sprintf(fname,"xg/D-%d-%d.xg",j,k);
	    xfile = fopen(fname,"w");

	    fprintf(xfile,"color=%s\n",data->assets[j].color);
	    fprintf(xfile,"thickness = 1.5\n");
	    for (i = is; i < data->num_days; ++i)
	    {
		if (value[i] == 0.0) continue;
		time = (double)i;
		ave_value = 0.5*(data->assets[j].norm_value[i] +
				 data->assets[k].norm_value[i]);
		fprintf(xfile,"%f %f\n",time,
			100*(data->assets[j].norm_value[i] - ave_value));
	    }
	    fprintf(xfile,"Next\n");
	    fprintf(xfile,"color=%s\n",data->assets[k].color);
	    fprintf(xfile,"thickness = 1.5\n");
	    for (i = is; i < data->num_days; ++i)
	    {
		if (value[i] == 0.0) continue;
		time = (double)i;
		ave_value = 0.5*(data->assets[j].norm_value[i] +
				 data->assets[k].norm_value[i]);
		fprintf(xfile,"%f %f\n",time,
			100*(data->assets[k].norm_value[i] - ave_value));
	    }
	    fclose(xfile);
	}
}	/* end XgraphData */

extern void TradeInfo(
	DATA_SET *data)
{
	int i,j;
	char string[100];
	boolean closed_trade_exist = NO;
	boolean partial_loop_exist = NO;

	SortTradeOrder(data);
	PrintOpenTrade(data);
	for (i = 0; i < data->num_trades; ++i)
	{
	    if (data->trades[i].status != OPEN)
	    {
	    	PrintClosedTradeLoop(NULL,data->trades[i],data);
		if (data->trades[i].status == CLOSED)
		    closed_trade_exist = YES;
		else if (data->trades[i].status == PARTIAL_CLOSED)
		    partial_loop_exist = YES;
	    }
	}
	if (closed_trade_exist)
	{
	    printf("Type yes to save and delete closed tradings: ");
	    scanf("%s",string);
	    if (string[0] == 'y' || string[0] == 'Y')
	    {
		SaveDeleteClosedTradeLoop(data);
		WriteAccountData(data);
	    }
	}
	if (partial_loop_exist)
	{
	    printf("Choose to save-close (s) or wrap (w) the loop: ");
	    scanf("%s",string);
	    if (string[0] == 's' || string[0] == 'S')
	    {
		SaveDeleteClosedTradeLoop(data);
		WriteAccountData(data);
	    }
	    else if (string[0] == 'w' || string[0] == 'W')
	    {
		WrapPartialTradeLoop(data);
		WriteAccountData(data);
	    }
	}
	printf("Type yes to print linear profile: ");
	scanf("%s",string);
	if (string[0] == 'y' || string[0] == 'Y')
	    PrintCurrentLinearProfile(data);
}	/* end TradeInfo */

extern void DataCompare(
	DATA_SET *data)
{
	int i1,i2,j;
	double nu1,nu2;
	double mu1,mu2;
	double norm_value1,norm_value2,dvalue;
	char string[100];
	int is,ie;
	char *stock_names[2];
	double prices[2];

	PrintAssetList(data);
	printf("Enter indices of two assets: ");
	scanf("%d %d",&i1,&i2);
	stock_names[0] = data->assets[i1].asset_name;
	stock_names[1] = data->assets[i2].asset_name;
	GetCurrentPrice(stock_names,prices,2);
	norm_value1 = 100.0*prices[0]/data->assets[i1].base_value;
	norm_value2 = 100.0*prices[1]/data->assets[i2].base_value;
	dvalue = norm_value1 - norm_value2;
	if (dvalue > 0.0)
	{
	    printf("\n%4s > %4s  %5.2f percent\n\n",stock_names[0],
				stock_names[1],dvalue);
	}
	else
	{
	    printf("\n%4s < %4s  %5.2f percent\n\n",stock_names[0],
				stock_names[1],-dvalue);
	}

	printf("Type yes to calculate amplification factor: ");
	scanf("%s",string);
	if (string[0] == 'y' || string[0] == 'Y')
	{
	    printf("Enter old values of two assets: ");
	    scanf("%lf %lf",&nu1,&nu2);
	    printf("Enter new values of two assets: ");
	    scanf("%lf %lf",&mu1,&mu2);
	    printf("The amplification factor = %f\n",nu1*mu2/nu2/mu1);
	}
}	/* end DataCompare */

extern void DataTrend(
	DATA_SET *data)
{
	int i,j,k,is,N;
	int M = data->num_assets;
	double QA,QB,QC,LA,LB,mLA,mLB;
	FILE *sfile0,*sfile1,*sfile2;
	FILE *qfile1,*qfile2;
	FILE *vfile0,*vfile1;
	FILE *mfile0,*mfile1;
	FILE *msfile0,*msfile1,*msfile2;
	double *value,time[MAX_TRACE],mean[MAX_TRACE];
	double value_ave,mean_ave[MAX_TRACE];
	double *pmean,*pmean_ave;
	double s0_max,s0[20],s1[20];
	int is0[20],is1[20];
	int is0_max;
	int iswap;
	char string[100];
	int num_data = 0;

	vfile0 = fopen("xg/v-0.xg","w");
	vfile1 = fopen("xg/v-1.xg","w");
	mfile0 = fopen("xg/m-0.xg","w");
	mfile1 = fopen("xg/m-1.xg","w");
	sfile0 = fopen("xg/s-0.xg","w");
	sfile1 = fopen("xg/s-1.xg","w");
	sfile2 = fopen("xg/s-2.xg","w");
	qfile1 = fopen("xg/c-1.xg","w");
	qfile2 = fopen("xg/c-2.xg","w");
	msfile0 = fopen("xg/ms-0.xg","w");
	msfile1 = fopen("xg/ms-1.xg","w");
	msfile2 = fopen("xg/ms-2.xg","w");

	N = data->num_backtrace - 12;
	if (N > MAX_TRACE - 12) N = MAX_TRACE - 12;

	for (i = 0; i < N; ++i) 
	    time[i] = i*0.25;
	is = data->num_days - N - 12;
	for (i = 0; i < N + 12; ++i) 
	{
	    mean[i] = 0.0;
	    mean_ave[i] = 0.0;
	    num_data = 0;
	    for (j = 0; j < M; ++j)
	    {
		if (!data->data_map[j]) continue;
		value = data->assets[j].norm_value + is + i;
		mean[i] += value[0];
		value_ave = (0.5*(value[0] + value[-3]) +
				  value[-1] + value[-2])/3.0;
		mean_ave[i] += value_ave;
		num_data++;
	    }
	    mean[i] /= num_data;
	    mean_ave[i] /= num_data;
	}

	for (i = 0; i < M; ++i)
	{
	    if (data->data_map[i] != YES)
	    	continue;
	    fprintf(sfile0,"color=%s\n",data->assets[i].color);
	    fprintf(sfile0,"thickness = 1.5\n");
	    fprintf(sfile1,"color=%s\n",data->assets[i].color);
	    fprintf(sfile1,"thickness = 1.5\n");
	    fprintf(sfile2,"color=%s\n",data->assets[i].color);
	    fprintf(sfile2,"thickness = 1.5\n");
	    fprintf(qfile1,"color=%s\n",data->assets[i].color);
	    fprintf(qfile1,"thickness = 1.5\n");
	    fprintf(qfile2,"color=%s\n",data->assets[i].color);
	    fprintf(qfile2,"thickness = 1.5\n");
	    fprintf(vfile0,"color=%s\n",data->assets[i].color);
	    fprintf(vfile0,"thickness = 1.5\n");
	    fprintf(vfile1,"color=%s\n",data->assets[i].color);
	    fprintf(vfile1,"thickness = 1.5\n");
	    fprintf(mfile0,"color=%s\n",data->assets[i].color);
	    fprintf(mfile0,"thickness = 1.5\n");
	    fprintf(mfile1,"color=%s\n",data->assets[i].color);
	    fprintf(mfile1,"thickness = 1.5\n");
	    fprintf(msfile0,"color=%s\n",data->assets[i].color);
	    fprintf(msfile0,"thickness = 1.5\n");
	    fprintf(msfile1,"color=%s\n",data->assets[i].color);
	    fprintf(msfile1,"thickness = 1.5\n");
	    fprintf(msfile2,"color=%s\n",data->assets[i].color);
	    fprintf(msfile2,"thickness = 1.5\n");

	    for (j = 0; j < N; ++j)
	    {
		is = data->num_days - N - 3;
		value = data->assets[i].norm_value + is + j;
		pmean = mean + j + 9;
		pmean_ave = mean_ave + j + 9;
		/* Value of the quater day */
		fprintf(vfile0,"%f  %f\n",time[j],100.0*value[3]);

		/* Average value of the day */
		value_ave = (0.5*(value[0] + value[3]) +
				value[1] + value[2])/3.0;
		fprintf(vfile1,"%f  %f\n",time[j],100.0*value_ave);
		fprintf(sfile0,"%f  %f\n",time[j],
				400.0*(value[3] - value[2]));
		if (j == N-1) s0[i] = 100.0*(value[3] - value[2])/0.25;

		fprintf(mfile0,"%f  %f\n",time[j],100.0*(value[3] - pmean[3]));
		fprintf(mfile1,"%f  %f\n",time[j],
				100.0*(value_ave - pmean_ave[3]));
		fprintf(msfile0,"%f  %f\n",time[j],
				400.0*((value[3] - value[2]) -
				       (pmean[3] - pmean[2])));

		/* Average slope of the day */
		LeastSquareLinear(time,value,4,&LA,&LB);
		LeastSquareLinear(time,pmean,4,&mLA,&mLB);
		fprintf(sfile1,"%f  %f\n",time[j],100.0*LA);
		fprintf(msfile1,"%f  %f\n",time[j],100.0*(LA - mLA));
		/* Average concavity of the day */
		LeastSquareQuadr(time,value,4,&QA,&QB,&QC);
		fprintf(qfile1,"%f  %f\n",time[j],100.0*QA);
		if (j == N-1) s1[i] = LA;

		/* Average slope of two days */
		is = data->num_days - N - 7;
		value = data->assets[i].norm_value + is + j;

		LeastSquareLinear(time,value,8,&LA,&LB);
		fprintf(sfile2,"%f  %f\n",time[j],100.0*LA);
		fprintf(msfile2,"%f  %f\n",time[j],100.0*(LA - mLA));
		/* Average concavity of two days */
		LeastSquareQuadr(time,value,8,&QA,&QB,&QC);
		fprintf(qfile2,"%f  %f\n",time[j],100.0*QA);

	    }
	    fprintf(vfile0,"Next\n");
	    fprintf(vfile1,"Next\n");
	    fprintf(sfile0,"Next\n");
	    fprintf(sfile1,"Next\n");
	    fprintf(sfile2,"Next\n");
	    fprintf(qfile1,"Next\n");
	    fprintf(qfile2,"Next\n");
	    fprintf(mfile0,"Next\n");
	    fprintf(mfile1,"Next\n");
	    fprintf(msfile0,"Next\n");
	    fprintf(msfile1,"Next\n");
	    fprintf(msfile2,"Next\n");
	}
	fprintf(vfile0,"Next\n");
	fprintf(vfile0,"color=aqua\n");
	fprintf(vfile0,"thickness = 2\n");
	for (j = 0; j < N; ++j)
	{
	    is = data->num_days - N - 3;
	    pmean = mean + j + 9;
	    fprintf(vfile0,"%f  %f\n",time[j],100.0*pmean[3]);
	}
	fclose(vfile0);
	fclose(vfile1);
	fclose(sfile0);
	fclose(sfile1);
	fclose(sfile2);
	fclose(qfile1);
	fclose(qfile2);
	fclose(mfile0);
	fclose(mfile1);
	fclose(msfile0);
	fclose(msfile1);
	fclose(msfile2);
}	/* end DataTrend */

extern DATA_SET *CopyData(
	DATA_SET *data)
{
	int i,j;
	DATA_SET *copy_data;
	double ave_value;
	
	FT_ScalarMemoryAlloc((POINTER*)&copy_data,sizeof(DATA_SET));
	copy_data->num_assets = data->num_assets;
	copy_data->num_days = data->num_days;
	copy_data->num_backtrace = data->num_backtrace;
	memcpy(copy_data->account_name,data->account_name,
				strlen(data->account_name)+1);
	FT_VectorMemoryAlloc((POINTER*)&copy_data->assets,data->num_assets+1,
				sizeof(ASSET));
	FT_VectorMemoryAlloc((POINTER*)&copy_data->data_map,data->num_assets+1,
				sizeof(boolean));
	FT_MatrixMemoryAlloc((POINTER*)&copy_data->date,data->num_days+10,100,
				sizeof(char));
	for (j = 0; j < data->num_assets+1; ++j)
	{
	    memcpy(copy_data->assets[j].asset_name,data->assets[j].asset_name,
				strlen(data->assets[j].asset_name)+1);
	    memcpy(copy_data->assets[j].color,data->assets[j].color,
				strlen(data->assets[j].color)+1);
	    copy_data->assets[j].base_value = data->assets[j].base_value;
	    FT_VectorMemoryAlloc((POINTER*)&copy_data->assets[j].value,
				data->num_days+10,sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&copy_data->assets[j].norm_value,
				data->num_days+10,sizeof(double));
	    copy_data->data_map[j] = data->data_map[j];
	}
	for (i = 0; i < data->num_days; ++i)
	{
	    	for (j = 0; j < data->num_assets+1; ++j)
		{
		    copy_data->assets[j].value[i] = data->assets[j].value[i];
		    copy_data->assets[j].norm_value[i] = 
			  data->assets[j].value[i]/data->assets[j].base_value;
		    ave_value += data->assets[j].norm_value[i];
		}
		ave_value = 0.0;
	    	for (j = 0; j < data->num_assets; ++j)
		    ave_value += data->assets[j].norm_value[i];
		copy_data->assets[j].value[i] = ave_value/data->num_assets;
		copy_data->assets[j].norm_value[i] = ave_value/data->num_assets;
	}
	return copy_data;
}	/* end CopyData */

extern void FreeData(
	DATA_SET *data)
{
	int i,j;
	for (j = 0; j < data->num_assets+1; ++j)
	    FT_FreeThese(2,data->assets[j].value,data->assets[j].norm_value);
	FT_FreeThese(3,data->assets,data->data_map,data->date);
	FT_FreeThese(1,data);
}	/* end FreeData */

static void PrintRateOfChange(
	DATA_SET *data)
{
	int i,n,M;
	double base_value,fit_value;
	double Q[3],L[2];

	n = data->num_days-1;
	M = data->num_assets;

	printf("\n");
	printf("          Immediate   "
	       		"     One-day     "
	       		"          Two-Day\n");
	printf("               L    "
	       		"   L      Q      M"
	       		"       L      Q      M\n");
	for (i = 0; i < M; ++i)
	{
	    if (data->data_map[i] != YES)
	        continue;
	    base_value = data->assets[i].base_value;
	    printf("%2d  %5s: ",i,data->assets[i].asset_name);
	    L[0] = (data->assets[i].norm_value[n] - 
			data->assets[i].norm_value[n-1])/0.25;
	    printf(" %6.2f ",100*L[0]);
	    GetLeastSquare(data,i,5,Q,L);
	    fit_value = base_value*(Q[0]+Q[1]+Q[2]);
	    printf("|%6.2f ",100*L[0]);
	    printf("%6.2f ",100*Q[0]);
	    printf("%6.2f ",fit_value);
	    GetLeastSquare(data,i,9,Q,L);
	    fit_value = base_value*(Q[0]*sqr(2.0)+Q[1]*2.0+Q[2]);
	    printf("|%6.2f ",100*L[0]);
	    printf("%6.2f ",100*Q[0]);
	    printf("%6.2f ",fit_value);
	    printf("\n");
	}
}	/* end PrintRateOfChange */

extern void ReadMarketData(
	DATA_SET *data)
{
	char *data_name = data->data_name;
	FILE *infile = fopen(data_name,"r");
	int i,j;
	
	if (infile == NULL)
	{
	    double ave_value = 0.0;
	    double ave_base_value = 0.0;
	    printf("Enter data set name: ");
	    scanf("%s",data->data_name);
	    sprintf(data_name,"%s",data->data_name);
	    printf("Enter number of assets: ");
	    scanf("%d",&data->num_assets);
	    data->num_days = 1;
	    FT_VectorMemoryAlloc((POINTER*)&data->assets,data->num_assets+1,
				sizeof(ASSET));
	    FT_VectorMemoryAlloc((POINTER*)&data->data_map,data->num_assets+1,
				sizeof(boolean));

	    for (j = 0; j < data->num_assets+1; ++j)
	    {
		FT_VectorMemoryAlloc((POINTER*)&data->assets[j].value,
				data->num_days+10,sizeof(double));
		FT_VectorMemoryAlloc((POINTER*)&data->assets[j].norm_value,
				data->num_days+10,sizeof(double));
		data->data_map[j] = YES;
		if (j < data->num_assets)
		{
		    printf("Enter name of asset %d: ",j+1);
		    scanf("%s",data->assets[j].asset_name);
		    printf("Enter color of asset %d: ",j+1);
		    scanf("%s",data->assets[j].color);
		    printf("Enter value of asset %d: ",j+1);
		    scanf("%lf",&data->assets[j].value[0]);
		    ave_value += data->assets[j].value[0];
		    printf("Enter base value of asset %d: ",j+1);
		    scanf("%lf",&data->assets[j].base_value);
		    data->assets[j].norm_value[0] = data->assets[j].value[0]/
				data->assets[j].base_value;
		    ave_base_value += data->assets[j].base_value;
		}
		else
		{
		    data->assets[j].value[0] = 1.0;
		    data->assets[j].base_value = 1.0;
		    data->assets[j].norm_value[0] = 1.0;
		    sprintf(data->assets[j].asset_name,"Mean");
		    sprintf(data->assets[j].color,"green");
		}
	    }
	}
	else
	{
	    double ave_value = 0.0;
	    fscanf(infile,"%*s");
	    fscanf(infile,"%d",&data->num_assets);
	    fscanf(infile,"%d",&data->num_days);
	    fscanf(infile,"%d",&data->num_backtrace);
	    fscanf(infile,"%d",&data->istar);
	    FT_VectorMemoryAlloc((POINTER*)&data->assets,data->num_assets+1,
				sizeof(ASSET));
	    FT_VectorMemoryAlloc((POINTER*)&data->data_map,data->num_assets+1,
				sizeof(boolean));
	    FT_MatrixMemoryAlloc((POINTER*)&data->date,data->num_days+10,100,
				sizeof(char));
	    for (j = 0; j < data->num_assets+1; ++j)
	    {
		fscanf(infile,"%s\n",data->assets[j].asset_name);
		fscanf(infile,"%s\n",data->assets[j].color);
		fscanf(infile,"%lf\n",&data->assets[j].base_value);
		FT_VectorMemoryAlloc((POINTER*)&data->assets[j].value,
				data->num_days+10,sizeof(double));
		FT_VectorMemoryAlloc((POINTER*)&data->assets[j].norm_value,
				data->num_days+10,sizeof(double));
	    }
	    for (i = 0; i < data->num_days; ++i)
	    {
	    	for (j = 0; j < data->num_assets+1; ++j)
		{
		    fscanf(infile,"%lf ",&data->assets[j].value[i]);
		    data->assets[j].norm_value[i] = data->assets[j].value[i]/
					data->assets[j].base_value;
		    ave_value += data->assets[j].norm_value[i];
		}
		ave_value = 0.0;
	    	for (j = 0; j < data->num_assets; ++j)
		    ave_value += data->assets[j].norm_value[i];
		data->assets[j].value[i] = ave_value/data->num_assets;
		data->assets[j].norm_value[i] = ave_value/data->num_assets;
	    }
	    fclose(infile);
	}
}	/* end ReadMarketData */

extern void ReadAccountData(
	DATA_SET *data)
{
	char *account_name = data->account_name;
	FILE *infile = fopen(account_name,"r");
	int i,j;
	
        FT_VectorMemoryAlloc((POINTER*)&data->shares,data->num_assets,
				sizeof(int));
        for (j = 0; j < data->num_assets; ++j)
	    data->shares[j] = 0;
	if (infile == NULL)
	{
            for (j = 0; j < data->num_assets; ++j)
		data->data_map[j] = YES;
	    return;
	}
	else
	{
	    if (fgetstring(infile,"Account Investment Polarization Ratio ="))
	    {
		fscanf(infile,"%lf",&data->polar_ratio);
	    }
	    if (fgetstring(infile,"Data Map"))
	    {
                for (j = 0; j < data->num_assets; ++j)
                    fscanf(infile,"%d",(int*)data->data_map+j);
	    }
	    if (fgetstring(infile,"Investment Portfolio"))
            {
                for (j = 0; j < data->num_assets; ++j)
                    fscanf(infile,"%d",data->shares+j);
            }
            FT_VectorMemoryAlloc((POINTER*)&data->trades,MAX_NUM_TRADE,
                                        sizeof(TRADE));
            if (fgetstring(infile,"Trading Record"))
            {
                fscanf(infile,"%d",&data->num_trades);
                for (j = 0; j < data->num_trades; ++j)
		{
                    fscanf(infile,"%d",&data->trades[j].num_stages);
                    fscanf(infile,"%d",(int*)&data->trades[j].status);
		    for (i = 0; i < data->trades[j].num_stages; ++i)
		    {
                    	fscanf(infile,"%d %d\n",
				&data->trades[j].index_sell[i],
                                &data->trades[j].index_buy[i]);
                    	fscanf(infile,"%d %d\n",
				&data->trades[j].num_shares_sell[i],
                                &data->trades[j].num_shares_buy[i]);
                    	fscanf(infile,"%lf %lf\n",
				&data->trades[j].price_sell[i],
                                &data->trades[j].price_buy[i]);
		    }
		}
            }
            if (fgetstring(infile,"Index of Share Equivalence"))
            {
                fscanf(infile,"%d\n",&data->eindex);
	    }
	    fclose(infile);
	}
}	/* end ReadAccountData */

extern void WriteMarketData(
	DATA_SET *data)
{
	FILE *outfile;	
	int i,j;
	char string[100];
	char *data_name = data->data_name;

	if (data->new_data == YES)
	{
	    printf("Save new data? ");
	    scanf("%s",string);
	    if (string[0] != 'y' && string[0] != 'Y')
	    	data->num_days--;
	    else
	    {
	    	printf("Adjust data base? ");
	    	scanf("%s",string);
	    	if (string[0] == 'y' || string[0] == 'Y')
		{
		    char dirname[256],fname[256];
		    FILE *bfile;
		    int m;
		    sprintf(dirname,"base/%s",data->data_name);
		    create_directory(dirname,NO);
		    AdjustBase(data);
	    	    for (m = 0; m < data->num_assets; ++m)	
	    	    {
			sprintf(fname,"%s/%s",dirname,
					data->assets[m].asset_name);
			if ((bfile = fopen(fname,"r")) == NULL)
			    bfile = fopen(fname,"w");
			else
			{
			    fclose(bfile);
			    bfile = fopen(fname,"a");
			}
			fprintf(bfile,"%f  %f\n",(double)data->num_days,
					data->assets[m].base_value);
			fclose(bfile);
	    	    }
		}
	    }
	}
	outfile = fopen(data_name,"w");
	fprintf(outfile,"%s\n",data->data_name);
	fprintf(outfile,"%d\n",data->num_assets);
	fprintf(outfile,"%d\n",data->num_days);
	fprintf(outfile,"%d\n",data->num_backtrace);
	fprintf(outfile,"%d\n",data->istar);
	for (j = 0; j < data->num_assets+1; ++j)
	{
	    fprintf(outfile,"%s\n",data->assets[j].asset_name);
	    fprintf(outfile,"%s\n",data->assets[j].color);
	    fprintf(outfile,"%f\n",data->assets[j].base_value);
	}
	for (i = 0; i < data->num_days; ++i)
	{
	    for (j = 0; j < data->num_assets+1; ++j)
	    {
		fprintf(outfile,"%f ",data->assets[j].value[i]);
	    }
	    fprintf(outfile,"\n");
	}
	fclose(outfile);
}	/* end WriteMarketData */

extern void WriteAccountData(
	DATA_SET *data)
{
	FILE *outfile;	
	int i,j;
	char string[100];
	char *account_name = data->account_name;

	outfile = fopen(account_name,"w");
	fprintf(outfile,"%s\n",data->account_name);
	fprintf(outfile,"Account Investment Polarization Ratio = %f\n",
				data->polar_ratio);
	fprintf(outfile,"Data Map\n");
	for (i = 0; i < data->num_assets; ++i)
	    fprintf(outfile,"%d\n",data->data_map[i]);
	if (data->shares)
	{
	    fprintf(outfile,"Investment Portfolio\n");
	    for (j = 0; j < data->num_assets; ++j)
		fprintf(outfile,"%d\n",data->shares[j]);
	}
	if (data->num_trades != 0)
	{
	    SortTradeOrder(data);
	    fprintf(outfile,"Trading Record\n");
	    fprintf(outfile,"%d\n",data->num_trades);
	    for (j = 0; j < data->num_trades; ++j)
	    {
	    	fprintf(outfile,"%d\n",data->trades[j].num_stages);
	    	fprintf(outfile,"%d\n",data->trades[j].status);
		for (i = 0; i < data->trades[j].num_stages; ++i)
		{
	    	    fprintf(outfile,"%d %d\n",data->trades[j].index_sell[i],
				data->trades[j].index_buy[i]);
	    	    fprintf(outfile,"%d %d\n",
				data->trades[j].num_shares_sell[i],
				data->trades[j].num_shares_buy[i]);
	    	    fprintf(outfile,"%f %f\n",data->trades[j].price_sell[i],
				data->trades[j].price_buy[i]);
		}
	    }
	}
        fprintf(outfile,"Index of Share Equivalence");
        fprintf(outfile,"%d\n",data->eindex);
	fclose(outfile);
}	/* end WriteAccountData */

extern void RankData(
	DATA_SET *data)
{
	int i,j;
	int n = data->num_days-1;
	int M = data->num_assets;
	double dvalue,d2value;
	char string[100];
	double Q[3],L[2];
	double fit_value,base_value;
	int *isort;
	double mean_norm_value = 0.0;
	int num_ave = 0;
	int is_first,is_last;
	DATA_SET *copy_data;

	copy_data = CopyData(data);
	printf("Number of backtrace: %d\n",copy_data->num_backtrace);

start_ranking:

	FT_VectorMemoryAlloc((POINTER*)&isort,M,sizeof(int));
	for (i = 0; i < M; ++i) 
	{
	    isort[i] = i;
	    num_ave++;
	    copy_data->assets[i].norm_value[n] =
			copy_data->assets[i].value[n]/
			copy_data->assets[i].base_value;
	    mean_norm_value += copy_data->assets[i].norm_value[n];
	}

	for (i = 0; i < M-1; ++i)
	for (j = i+1; j < M; ++j)
	{
	    if (copy_data->assets[isort[i]].norm_value[n] < 
		copy_data->assets[isort[j]].norm_value[n])
	    {
		int itmp = isort[j];
		isort[j] = isort[i];
		isort[i] = itmp;
	    }
	}
	is_first = is_last = -1;
	for (i = 0; i < M; ++i)
	{
	    if (copy_data->data_map[isort[i]] == YES && is_first == -1)
		is_first = i;
	    if (copy_data->data_map[isort[M-i-1]] == YES && is_last == -1)
		is_last = M-i-1;
	}
	printf("\t");
	for (i = 0; i < M; ++i)
	{
	    if (i == is_first) continue;
	    if (copy_data->data_map[isort[i]] != YES)
	    	continue;
	    printf("%6s  ",copy_data->assets[isort[i]].asset_name);
	}
	printf("\n");
	for (i = 0; i < M; ++i)
	{
	    if (copy_data->data_map[isort[i]] != YES)
	    	continue;
	    if (i == is_last) 
	    {
	    	printf("\n");
		break;
	    }
	    printf("%6s  ",copy_data->assets[isort[i]].asset_name);
	    for (j = 0; j <= i; ++j)
	    {
	    	if (j == is_first) continue;
	    	if (copy_data->data_map[isort[j]] != YES)
	    	    continue;
		printf("        ");
	    }
	    for (j = i+1; j < M; ++j)
	    {
	    	if (copy_data->data_map[isort[j]] != YES)
	    	    continue;
		double dvalue = (copy_data->assets[isort[i]].norm_value[n]
			- copy_data->assets[isort[j]].norm_value[n])*100.0;
		printf("%6.2f  ",dvalue);
	    }
	    printf("\n");
	}
	printf("Change back trace and re-rank: ");
	scanf("%s",string);
	if (string[0] == 'y')
	{
	    printf("Enter new number of backtrace: ");
	    scanf("%d",&copy_data->num_backtrace);
	    AdjustBase(copy_data);
	    goto start_ranking;
	}
	FT_FreeThese(1,isort);
	FreeData(copy_data);
}	/* end RankData */


extern void PrintAssetList(
	DATA_SET *data)
{
	int i;
	for (i = 0; i < data->num_assets; ++i)
	{
	    printf("Assets %d: %4s\n",i,data->assets[i].asset_name);
	}
}	/* end PrintAssetList */

extern void InitTrade(
	DATA_SET *data)
{
	int i1,i2,j;
	double nu1,nu2;
	double mu1,mu2;
	double norm_value1,norm_value2,dvalue;
	char string[100];
	int is,ie;
	char fname[200];
	int i;
	FILE *xfile;
	double time,ratio,v[2],v_ave;
	char *stock_names[2];

	PrintAssetList(data);
	printf("Enter indices of two assets: ");
	scanf("%d %d",&i1,&i2);
	is = (data->num_days < data->num_backtrace) ? 0 :
                	data->num_days - data->num_backtrace;
	ie = data->num_days;
	stock_names[0] = data->assets[i1].asset_name;
	stock_names[1] = data->assets[i2].asset_name;
	GetCurrentPrice(stock_names,v,2);
	create_directory("xg",NO);
	sprintf(fname,"xg/I-ratio.xg");
	xfile = fopen(fname,"w");
        fprintf(xfile,"color=red\n");
        fprintf(xfile,"thickness = 1.5\n");
	for (i = is; i < ie; ++i)
        {
            time = (double)i;
	    ratio = data->assets[i1].norm_value[i]/
				data->assets[i2].norm_value[i];
            fprintf(xfile,"%f %f\n",time,ratio);
        }
        fprintf(xfile,"%f %f\n",time+1.0,v[0]/data->assets[i1].base_value
				/v[1]*data->assets[i2].base_value);
        fprintf(xfile,"\nNext\n");
        fprintf(xfile,"color=green\n");
        fprintf(xfile,"thickness = 1.5\n");
	for (i = is; i < ie; ++i)
        {
            time = (double)(i-2);
	    ratio = 0.0;
	    for (j = i-3; j <= i; ++j)
	    	ratio += data->assets[i1].norm_value[j]/
				data->assets[i2].norm_value[j];
	    ratio /= 4.0;
            fprintf(xfile,"%f %f\n",time,ratio);
        }
	if (string[0] == 'y' || string[0] == 'Y')
	{
            time = (double)(i-2);
	    ratio = 0.0;
	    for (j = i-3; j < i; ++j)
	    	ratio += data->assets[i1].norm_value[j]/
				data->assets[i2].norm_value[j];
	    ratio += v[0]/data->assets[i1].base_value
                                /v[1]*data->assets[i2].base_value;
	    ratio /= 4.0;
            fprintf(xfile,"%f %f\n",time,ratio);
	}
        fprintf(xfile,"\nNext\n");
        fprintf(xfile,"color=blue\n");
        fprintf(xfile,"thickness = 1.5\n");
	for (i = is; i < ie; ++i)
        {
            time = (double)(i-6);
	    ratio = 0.0;
	    for (j = i-11; j <= i; ++j)
	    	ratio += data->assets[i1].norm_value[j]/
				data->assets[i2].norm_value[j];
	    ratio /= 12.0;
            fprintf(xfile,"%f %f\n",time,ratio);
        }
	if (string[0] == 'y' || string[0] == 'Y')
	{
            time = (double)(i-6);
	    ratio = 0.0;
	    for (j = i-11; j < i; ++j)
	    	ratio += data->assets[i1].norm_value[j]/
				data->assets[i2].norm_value[j];
	    ratio += v[0]/data->assets[i1].base_value
                                /v[1]*data->assets[i2].base_value;
	    ratio /= 12.0;
            fprintf(xfile,"%f %f\n",time,ratio);
	}
	fclose(xfile);

	is = data->num_days - data->num_backtrace;
	ie = data->num_days;
	sprintf(fname,"xg/Compare.xg");
	xfile = fopen(fname,"w");
        fprintf(xfile,"color=%s\n",data->assets[i1].color);
        fprintf(xfile,"thickness = 1.5\n");
	for (i = is; i < ie; ++i)
	{
	    time = (double)i;
	    v_ave = 0.5*(data->assets[i1].norm_value[i] +
			     	data->assets[i2].norm_value[i]);
	    fprintf(xfile,"%f  %f\n",time,
				data->assets[i1].norm_value[i]-v_ave);
	}
	if (string[0] == 'y' || string[0] == 'Y')
	{
	    v_ave = 0.5*(v[0]/data->assets[i1].base_value +
				v[1]/data->assets[i2].base_value);
	    fprintf(xfile,"%f  %f\n",time+1,
				v[0]/data->assets[i1].base_value-v_ave);
	}
        fprintf(xfile,"\nNext\n");
        fprintf(xfile,"color=%s\n",data->assets[i2].color);
        fprintf(xfile,"thickness = 1.5\n");
	for (i = is; i < ie; ++i)
	{
	    time = (double)i;
	    v_ave = 0.5*(data->assets[i1].norm_value[i] +
			     	data->assets[i2].norm_value[i]);
	    fprintf(xfile,"%f  %f\n",time,
				data->assets[i2].norm_value[i]-v_ave);
	}
	if (string[0] == 'y' || string[0] == 'Y')
	{
	    v_ave = 0.5*(v[0]/data->assets[i1].base_value +
				v[1]/data->assets[i2].base_value);
	    fprintf(xfile,"%f  %f\n",time+1,
				v[1]/data->assets[i2].base_value-v_ave);
	}
	fclose(xfile);
	printf("Comparing two assets:\n");
        printf("%d  %4s  color=%s\n",i1,data->assets[i1].asset_name,
				data->assets[i1].color);
        printf("%d  %4s  color=%s\n",i2,data->assets[i2].asset_name,
				data->assets[i2].color);
}	/* end InitTrade */

extern void PrintEquityIndex(DATA_SET *data)
{
	int i;
	int M = data->num_assets;
	for (i = 0; i < M/3; ++i)
	{
	    printf("%5s %2d   %5s %2d   %5s %2d\n",
			data->assets[3*i].asset_name,3*i,
			data->assets[3*i+1].asset_name,3*i+1,
			data->assets[3*i+2].asset_name,3*i+2);
	}
	for (i = 3*(M/3); i < M; ++i)
	    printf("%5s %2d   ",data->assets[i].asset_name,i);
	printf("\n");
}	/* end PrintEquityIndex */

extern void PrintAccountValue(
	DATA_SET *data)
{
	int i,M,n;
	double market_value = 0.0;
	char string[100];
	M = data->num_assets;
	n = data->num_days-1;
	for (i = 0; i < M; ++i)
	{
	    market_value += data->shares[i]*data->assets[i].value[n];
	}
	printf("Total account market value: %f\n",market_value);
	printf("Type yes to print shares: ");
	scanf("%s",string);
	if (string[0] == 'y')
	{
	    for (i = 0; i < M; ++i)
	    	printf("%4s:  %d\n",data->assets[i].asset_name,
				data->shares[i]);
	}
}	/* end PrintAccountValue */

static double AssetTimeIntegral(
	DATA_SET *data,
	int index,
	int istart,
	int iend)
{
	int i,j;
	double I = 0.0;

	for (i = istart+1; i < iend-1; ++i)
	    I += data->assets[index].value[i]; 
	I += 0.5*(data->assets[index].value[istart] +
		  data->assets[index].value[iend-1]);
	return I;
}	/* end AssetTimeIntegral */

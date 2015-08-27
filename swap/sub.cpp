#include "swap.h"

static double NormDiff(DATA_SET*,int,int);
static void AddAssets(DATA_SET*);
static void DeleteAssets(DATA_SET*);

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
	int istar;		  /* center asset index */
	int time_p;    /* time period */
	int i, j, k, m;
	double ave_value;
	DATA_SET * period_data;  /* time period DATA_SET */
	FT_ScalarMemoryAlloc((POINTER*)&period_data,sizeof(DATA_SET));

	printf("There are total %d days of data\n",data->num_days);
	printf("Enter number of days for base update: ");
	scanf("%d", &time_p);
	FT_VectorMemoryAlloc((POINTER*)&period_data->assets,data->num_assets+1,
					sizeof(ASSET));
	FT_MatrixMemoryAlloc((POINTER*)&period_data->date,data->num_days,100,
					sizeof(char));
	strcpy(period_data->data_name, "period_data");
	period_data->num_days = time_p;
	period_data->num_assets = data->num_assets;
	period_data->num_backtrace = period_data->num_days;
   	for (j = 0; j < data->num_assets+1; ++j)
    	{
	    FT_VectorMemoryAlloc((POINTER*)&period_data->assets[j].value,
				period_data->num_days,sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&period_data->assets[j].norm_value,
				period_data->num_days,sizeof(double));
	}
	for (i = 0; i < period_data->num_assets + 1; i++ )
	{
	    strcpy(period_data->assets[i].asset_name, 
				data->assets[i].asset_name);
	    period_data->assets[i].base_value = data->assets[i].base_value;
	    strcpy(period_data->assets[i].color, data->assets[i].color);
	}
	printf("Available assets are: \n");
   	for (j = 0; j < data->num_assets; ++j)
	    printf("\t%s (%d)\n",data->assets[j].asset_name,j);
	printf("Enter center asset index: ");
	scanf("%d",&istar);
	for (j = 0; j <= data->num_days - time_p; j++ )
	{
	    for (k = 0; k < time_p; k++)
	    {
		ave_value = 0;
		for (i = 0; i < period_data->num_assets; i++)
		{
		    period_data->assets[i].value[k] = 
				data->assets[i].value[j+k];
		    period_data->assets[i].norm_value[k] = 
				period_data->assets[i].value[k]/
				period_data->assets[i].base_value;
		    ave_value += period_data->assets[i].norm_value[k];
		}
		period_data->assets[i].value[k] = 
				ave_value/period_data->num_assets;
		period_data->assets[i].norm_value[k] = 
				period_data->assets[i].value[k];
	    }

	    for (m = 0; m < data->num_assets; ++m)	
	    {
		if (m == istar) continue;
		AdjustBase(period_data,istar,m);
	    }
	}
	for (j = 0; j < data->num_days; j++)
	{
	    ave_value = 0;
	    for (i = 0; i < period_data->num_assets ; i++)
	    {
		data->assets[i].base_value = period_data->assets[i].base_value;
		data->assets[i].norm_value[j] = data->assets[i].value[j]/
				data->assets[i].base_value;
		ave_value += data->assets[i].norm_value[j];
	    }
	}
	data->assets[i].value[j] = ave_value / data->num_assets;
	data->assets[i].norm_value[j] = data->assets[i].value[j];
}	/* end ContinueBaseAdjust */

extern void AdjustBase(
	DATA_SET *data,
        int istar,
        int iplnt)
{
	double xl,xu;
	double yl,yu;
	double save_base;
	int i,j,N = 10;
	int num_days = data->num_days;
	int is;
	double ave_value = 0.0;

	is = (data->num_days < data->num_backtrace) ? 0 :
			data->num_days - data->num_backtrace;
	for (i = is; i < num_days; ++i)
	    ave_value += data->assets[istar].value[i];

	data->assets[istar].base_value = ave_value/(num_days-is);
	//data->assets[istar].base_value = data->assets[istar].value[is];
	data->assets[iplnt].base_value = data->assets[iplnt].value[is];

	xl = xu = data->assets[iplnt].base_value;
	yl = yu = NormDiff(data,istar,iplnt);
	if (yl > 0.0)
	{
	    for (i = 0; i < N; ++i)
	    {
		data->assets[iplnt].base_value *= 2.0;
		if ((yl = NormDiff(data,istar,iplnt)) < 0.0)
		{
		    xl = data->assets[iplnt].base_value;
		    break;
		}
	    }
	}
	else
	{
	    for (i = 0; i < N; ++i)
	    {
		data->assets[iplnt].base_value /= 2.0;
		if ((yu = NormDiff(data,istar,iplnt)) > 0.0)
		{
		    xu = data->assets[iplnt].base_value;
		    break;
		}
	    }
	}
	if (i == 3)
	{
	    printf("Cannot find two starting base\n");
	    clean_up(ERROR);	
	}
	N += 22;
	for (i = 0; i < N; ++i)
	{
	    save_base = data->assets[iplnt].base_value;
	    data->assets[iplnt].base_value = 0.5*(xl + xu);
	    for (j = 0; j < num_days; ++j)
	    {
		if (data->assets[iplnt].value[j] == save_base)
		{
		    data->assets[iplnt].value[j] = 0.5*(xl + xu);
		    data->assets[iplnt].norm_value[j] = 1.0;
		}
	    }
	    if (NormDiff(data,istar,iplnt) > 0.0)
	    	xu = data->assets[iplnt].base_value;
	    else
	    	xl = data->assets[iplnt].base_value;
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

	printf("Type yes to modify data base: ");
	scanf("%s",string);
	if (string[0] != 'y' && string[0] != 'Y')
	    return;

	printf("Available modification types are:\n");
	printf("\tManual (M)\n");
	printf("\tAutomatic (A)\n");
	printf("Enter type of modification: ");
	scanf("%s",string);
	switch (string[0])
	{
	case 'm':
	case 'M':
	    for (j = 0; j < data->num_assets; ++j)
	    {
	    	printf("Enter modification factor of %s base value: ",
			data->assets[j].asset_name);
		scanf("%lf",&mfactor);
		data->assets[j].base_value *= mfactor;
		data->assets[j].value[0] *= mfactor;
	    }
	    for (i = 0; i < data->num_days; ++i)
	    {
		ave_value = 0.0;
	    	for (j = 0; j < data->num_assets; ++j)
		{
		    data->assets[j].norm_value[i] = data->assets[j].value[i]/
					data->assets[j].base_value;
		    ave_value += data->assets[j].norm_value[i];
		}
		data->assets[j].value[i] = ave_value/data->num_assets;
		data->assets[j].norm_value[i] = ave_value/data->num_assets;
	    }
	    break;
	case 'a':
	case 'A':
	    printf("Available assets are:\n");
	    for (j = 0; j < data->num_assets; ++j)
	    	printf("\t%s (%d)\n",data->assets[j].asset_name,j);
	    printf("Enter center asset index: ");
	    scanf("%d",&istar);
	    for (j = 0; j < data->num_assets; ++j)	
	    {
		if (j == istar) continue;
		AdjustBase(data,istar,j);
	    }
	    break;
	default:
	    printf("Unknow type of modification!\n");
	    clean_up(ERROR);
	}
	printf("Print normalized difference? ");
	scanf("%s",string);
	if (string[0] != 'y' && string[0] != 'Y')
	    return;
	else
	{
	    for (i = 0; i < data->num_assets-1; ++i)
	    for (j = i+1; j < data->num_assets; ++j)
		printf("Diff-%s-%s = %f\n",
			data->assets[i].asset_name,data->assets[j].asset_name,
			NormDiff(data,i,j));
	}
}	/* end ModifyData */

static double NormDiff(
	DATA_SET *data,
	int istar,
	int iplnt)
{
	int i,is;
	double area = 0.0;

	is = (data->num_days < data->num_backtrace) ? 0 :
			data->num_days - data->num_backtrace;
	data->assets[iplnt].norm_value[0] = 
		data->assets[iplnt].value[0]/data->assets[iplnt].base_value;
	for (i = is; i < data->num_days; ++i)
	{
	    data->assets[istar].norm_value[i] = 
		data->assets[istar].value[i]/data->assets[istar].base_value;
	    data->assets[iplnt].norm_value[i] = 
		data->assets[iplnt].value[i]/data->assets[iplnt].base_value;
	}
	for (i = is+1; i < data->num_days; ++i)
	{
	    area += 0.5*(data->assets[iplnt].norm_value[i-1] -
	    		 data->assets[istar].norm_value[i-1] +
		 	 data->assets[iplnt].norm_value[i] -
		 	 data->assets[istar].norm_value[i]);
	}
	return area;
}	/* end NormDiff */

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
	for (i = 0; i < data->num_assets; ++i)
	    data->data_map[i] = YES;
	printf("Print data information of all assets? ");
	scanf("%s",string);
	if (string[0] == 'y') return;
	for (i = 0; i < data->num_assets; ++i)
	{
	    printf("Include asset %s ? ",data->assets[i].asset_name);
	    scanf("%s",string);
	    if (string[0] != 'y')
	    	data->data_map[i] = NO;
	}
}	/* end PromptForDataMap */

extern void ReadData(
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
	    if (fgetstring(infile,"Investment Portfolio"))
            {
                FT_VectorMemoryAlloc((POINTER*)&data->shares,data->num_assets,
                                        sizeof(int));
                for (j = 0; j < data->num_assets; ++j)
                    fscanf(infile,"%d",data->shares+j);
            }
            if (fgetstring(infile,"Trading Record"))
            {
                FT_VectorMemoryAlloc((POINTER*)&data->trades,MAX_NUM_TRADE,
                                        sizeof(TRADE));
                fscanf(infile,"%d",&data->num_trades);
                for (j = 0; j < data->num_trades; ++j)
		{
                    fscanf(infile,"%d",(int*)&data->trades[j].closed);
                    fscanf(infile,"%d %d\n",data->trades[j].index_buy,
                                data->trades[j].index_buy+1);
                    fscanf(infile,"%d %d\n",data->trades[j].index_sell,
                                data->trades[j].index_sell+1);
                    fscanf(infile,"%d %d\n",data->trades[j].num_shares_buy,
                                data->trades[j].num_shares_buy+1);
                    fscanf(infile,"%d %d\n",data->trades[j].num_shares_sell,
                                data->trades[j].num_shares_sell+1);
                    fscanf(infile,"%lf %lf\n",data->trades[j].price_buy,
                                data->trades[j].price_buy+1);
                    fscanf(infile,"%lf %lf\n",data->trades[j].price_sell,
                                data->trades[j].price_sell+1);
		}
            }
	    fclose(infile);
	}
}	/* end ReadData */

extern void WriteData(
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
	}
	outfile = fopen(data_name,"w");
	fprintf(outfile,"%s\n",data->data_name);
	fprintf(outfile,"%d\n",data->num_assets);
	fprintf(outfile,"%d\n",data->num_days);
	fprintf(outfile,"%d\n",data->num_backtrace);
	printf("\n\n");
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
	if (data->shares)
	{
	    fprintf(outfile,"Investment Portfolio\n");
	    for (j = 0; j < data->num_assets; ++j)
		fprintf(outfile,"%d\n",data->shares[j]);
	}
	if (data->num_trades != 0)
	{
	    fprintf(outfile,"Trading Record\n");
	    fprintf(outfile,"%d\n",data->num_trades);
	    for (j = 0; j < data->num_trades; ++j)
	    {
	    	fprintf(outfile,"%d\n",data->trades[j].closed);
	    	fprintf(outfile,"%d %d\n",data->trades[j].index_buy[0],
				data->trades[j].index_buy[1]);
	    	fprintf(outfile,"%d %d\n",data->trades[j].index_sell[0],
				data->trades[j].index_sell[1]);
	    	fprintf(outfile,"%d %d\n",data->trades[j].num_shares_buy[0],
				data->trades[j].num_shares_buy[1]);
	    	fprintf(outfile,"%d %d\n",data->trades[j].num_shares_sell[0],
				data->trades[j].num_shares_sell[1]);
	    	fprintf(outfile,"%f %f\n",data->trades[j].price_buy[0],
				data->trades[j].price_buy[1]);
	    	fprintf(outfile,"%f %f\n",data->trades[j].price_sell[0],
				data->trades[j].price_sell[1]);
	    }
	}
	fclose(outfile);
}	/* end WriteData */

extern void AddData(
	DATA_SET *data)
{
	int i,j,m,n;
	char string[100];
	double ave_value;

	printf("Type yes to add new entry data : ");
	data->new_data = NO;
	scanf("%s",string);
	if (string[0] != 'y' && string[0] != 'Y')
	    return;
	data->new_data = YES;
	n = 1;
	m = data->num_days;
	for (i = 0; i < n; ++i)
	{
	    ave_value = 0.0;
	    for (j = 0; j < data->num_assets; ++j)
	    {
	    	printf("Enter value of %s for entry+%d: ",
			data->assets[j].asset_name,i+1);
		scanf("%lf",&data->assets[j].value[m+i]);
		if (data->assets[j].value[m+i] == -1)
		    data->assets[j].value[m+i] = data->assets[j].value[m+i-1];
		data->assets[j].norm_value[m+i] = data->assets[j].value[m+i]/
					data->assets[j].base_value;
		ave_value += data->assets[j].norm_value[m+i];
	    }
	    data->assets[j].base_value = 1.0;
	    data->assets[j].value[m+i] = ave_value/data->num_assets;
	    data->assets[j].norm_value[m+i] = data->assets[j].value[m+i];
	}
	data->num_days += n;
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

extern void DataInfo(
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

	printf("Type yes to get rate of change: ");
	scanf("%s",string);
	if (string[0] == 'y' || string[0] == 'Y')
	{

	    printf("\n");
	    printf("          Immediate   "
	       		"     One-day     "
	       		"          Two-Day\n");
	    printf("               L    "
	       		"   L      Q      M"
	       		"       L      Q      M\n");
	    n = data->num_days - 1;
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
	}

	printf("\n");
	FT_VectorMemoryAlloc((POINTER*)&isort,M,sizeof(int));
	for (i = 0; i < M; ++i) 
	{
	    isort[i] = i;
	    num_ave++;
	    mean_norm_value += data->assets[i].norm_value[n];
	}

	for (i = 0; i < M-1; ++i)
	for (j = i+1; j < M; ++j)
	{
	    if (data->assets[isort[i]].norm_value[n] < 
		data->assets[isort[j]].norm_value[n])
	    {
		int itmp = isort[j];
		isort[j] = isort[i];
		isort[i] = itmp;
	    }
	}
	is_first = is_last = -1;
	for (i = 0; i < M; ++i)
	{
	    if (data->data_map[isort[i]] == YES && is_first == -1)
		is_first = i;
	    if (data->data_map[isort[M-i-1]] == YES && is_last == -1)
		is_last = M-i-1;
	}
	printf("\t");
	for (i = 0; i < M; ++i)
	{
	    if (i == is_first) continue;
	    if (data->data_map[isort[i]] != YES)
	    	continue;
	    printf("%6s  ",data->assets[isort[i]].asset_name);
	}
	printf("\n");
	for (i = 0; i < M; ++i)
	{
	    if (data->data_map[isort[i]] != YES)
	    	continue;
	    if (i == is_last) 
	    {
	    	printf("\n");
		break;
	    }
	    printf("%6s  ",data->assets[isort[i]].asset_name);
	    for (j = 0; j <= i; ++j)
	    {
	    	if (j == is_first) continue;
	    	if (data->data_map[isort[j]] != YES)
	    	    continue;
		printf("        ");
	    }
	    for (j = i+1; j < M; ++j)
	    {
	    	if (data->data_map[isort[j]] != YES)
	    	    continue;
		double dvalue = (data->assets[isort[i]].norm_value[n]
			- data->assets[isort[j]].norm_value[n])*100.0;
		printf("%6.2f  ",dvalue);
	    }
	    printf("\n");
	}
	printf("Type to print profitable tradings: ");
	scanf("%s",string);
	if (string[0] == 'y' || string[0] == 'Y')
	    PrintProfitableTrade(data);
	printf("\nCurrent average normalized value: %7.4f\n",
			mean_norm_value/num_ave);
	FT_FreeThese(1,isort);
}	/* end DataInfo */

extern void DataCompare(
	DATA_SET *data)
{
	int i1,i2,j;
	double nu1,nu2;
	double mu1,mu2;
	double norm_value1,norm_value2,dvalue;
	char string[100];
	int is,ie;

	for (j = 0; j < data->num_assets; ++j)
	{
	    printf("Index of assets are: %d for %4s\n",j,
		data->assets[j].asset_name);
	}
	printf("Enter indices of two assets: ");
	scanf("%d %d",&i1,&i2);
	printf("Type yes to calculate normalized difference: ");
	scanf("%s",string);
	if (string[0] == 'y' || string[0] == 'Y')
	{
	    printf("Enter values of two assets: ");
	    scanf("%lf %lf",&nu1,&nu2);
	    norm_value1 = 100.0*nu1/data->assets[i1].base_value;
	    norm_value2 = 100.0*nu2/data->assets[i2].base_value;
	    dvalue = norm_value1 - norm_value2;
	    if (dvalue > 0.0)
	    {
	    	printf("\n%4s > %4s  %5.2f percent\n\n",
				data->assets[i1].asset_name,
				data->assets[i2].asset_name,dvalue);
	    }
	    else
	    {
	    	printf("\n%4s < %4s  %5.2f percent\n\n",
				data->assets[i1].asset_name,
				data->assets[i2].asset_name,-dvalue);
	    }
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
	printf("Type yes to plot initiation ratio: ");
	scanf("%s",string);
	if (string[0] == 'y' || string[0] == 'Y')
	{
	    char fname[200];
	    int i;
	    FILE *xfile;
	    double time,ratio;
	    is = (data->num_days < data->num_backtrace) ? 0 :
                	data->num_days - data->num_backtrace;
	    ie = data->num_days;
	    create_directory("xg",NO);
	    sprintf(fname,"xg/I-ratio.xg");
	    xfile = fopen(fname,"w");
            fprintf(xfile,"color=red\n");
            fprintf(xfile,"thickness = 1.5\n");
	    for (i = is; i < ie; ++i)
            {
                time = (double)i;
		ratio = data->assets[i1].value[i]/data->assets[i2].value[i];
                fprintf(xfile,"%f %f\n",time,ratio);
            }
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
	memcpy(copy_data->data_name,data->data_name,strlen(data->data_name)+1);
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

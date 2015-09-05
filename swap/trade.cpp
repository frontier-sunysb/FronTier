#include "swap.h"
#include <ctime>
#include <iostream>
#include <locale>

static void PrintOpenTrade(DATA_SET*);
static double AmpRatio(TRADE,DATA_SET*,int);
static double currentAmpRatio(TRADE,DATA_SET*,int);
static int StartIndex(TRADE);

#define 	MAX_NUM_TRADE		200

extern void InvestShares(
	DATA_SET *data)
{
	char string[100];
	int num_shares;
	int i;

	if (data->shares == NULL)
	{
	    FT_VectorMemoryAlloc((POINTER*)&data->shares,data->num_assets,
					sizeof(int));
	    for (i = 0; i < data->num_assets; ++i)
		data->shares[i] = 0;
	}
	printf("Type yes to add investment: ");
	scanf("%s",string);
	if (string[0] == 'y' || string[0] == 'Y')
	{
	    for (i = 0; i < data->num_assets; ++i)
	    {
		printf("%-4s current share number %d, add: ",
			data->assets[i].asset_name,data->shares[i]);
		scanf("%d",&num_shares);
		data->shares[i] += num_shares;
	    }
	    printf("Investment update:\n");
	    for (i = 0; i < data->num_assets; ++i)
		printf("%-4s: %d\n",data->assets[i].asset_name,
				data->shares[i]);
	}
}	/* end InvestShares */

extern void TradeShares(
	DATA_SET *data)
{
	boolean to_trade = YES;
	int i,n,is,ib,ns,nb;
	char string[100];

	printf("The current open trades are:\n");
	PrintOpenTrade(data);
	while (to_trade)
	{
	    printf("Type yes to trade: ");
	    scanf("%s",string);
	    if (string[0] != 'y' && string[0] != 'Y')
		break;
	    if (data->num_trades >= MAX_NUM_TRADE)
	    {
		printf("There are too many open trades!\n");
		break;
	    }
	    printf("Type yes if this is a new trade: ");
	    scanf("%s",string);
	    if (string[0] == 'y' || string[0] == 'Y')
	    {
		n = data->num_trades;
		for (i = 0; i < data->num_assets; ++i)
		{
	    	    printf("Index of assets are: %d for %4s\n",i,
			    data->assets[i].asset_name);
		}
	    	printf("Enter equity index, number of shares "
			"and price to sell: ");
		scanf("%d %d %lf",&data->trades[n].index_sell[0],
				&data->trades[n].num_shares_sell[0],
				&data->trades[n].price_sell[0]);
	    	printf("Enter equity index, number of shares "
			"and price to buy: ");
		scanf("%d %d %lf",&data->trades[n].index_buy[0],
				&data->trades[n].num_shares_buy[0],
				&data->trades[n].price_buy[0]);
		data->trades[n].status = OPEN;
		is = data->trades[n].index_sell[0];
		ns = data->trades[n].num_shares_sell[0];
		ib = data->trades[n].index_buy[0];
		nb = data->trades[n].num_shares_buy[0];
		data->shares[is] -= ns;
		data->shares[ib] += nb;
		(data->num_trades)++;
		data->trades[n].num_stages++;
	    }
	    else
	    {
		int nstage;
		PrintOpenTrade(data);
		printf("Select index of open trade: ");
		scanf("%d",&n);
		nstage = data->trades[n].num_stages;
	    	printf("Enter equity index, number of shares "
			"and price to sell: ");
		scanf("%d %d %lf",&data->trades[n].index_sell[nstage],
				&data->trades[n].num_shares_sell[nstage],
				&data->trades[n].price_sell[nstage]);
	    	printf("Enter equity index, number of shares "
			"and price to buy: ");
		scanf("%d %d %lf",&data->trades[n].index_buy[nstage],
				&data->trades[n].num_shares_buy[nstage],
				&data->trades[n].price_buy[nstage]);
		if (data->trades[n].index_buy[nstage-1] != 
		    data->trades[n].index_sell[nstage])
		{
		    printf("Invalid matching trade!\n");
		    continue;
		}
		is = data->trades[n].index_sell[nstage];
		ns = data->trades[n].num_shares_sell[nstage];
		ib = data->trades[n].index_buy[nstage];
		nb = data->trades[n].num_shares_buy[nstage];
		data->shares[is] -= ns;
		data->shares[ib] += nb;
		if (data->trades[n].index_buy[nstage] ==
		    data->trades[n].index_sell[0])
		    data->trades[n].status = CLOSED;
		else
		{
		    for (i = 1; i < data->trades[n].num_stages; ++i)
		    {
			if (data->trades[n].index_buy[nstage] ==
                    	    data->trades[n].index_sell[i])
			{
		    	    data->trades[n].status = PARTIAL_CLOSED;
			    break;
			}
		    }
		}
		data->trades[n].num_stages++;
	    }
	}
	printf("Investment update:\n");
	for (i = 0; i < data->num_assets; ++i)
	    printf("%-4s: %d\n",data->assets[i].asset_name,data->shares[i]);
}	/* end TradeShares */

static void PrintOpenTrade(
	DATA_SET *data)
{
	int i,j;
	int ib,is;
	int nb,ns;
	double pb,ps;
	double npb,nps;
	for (i = 0; i < data->num_trades; ++i)
	{
	    if (data->trades[i].status == CLOSED) continue;
	    ib = data->trades[i].index_buy[0];
	    is = data->trades[i].index_sell[0];
	    nb = data->trades[i].num_shares_buy[0];
	    ns = data->trades[i].num_shares_sell[0];
	    pb = data->trades[i].price_buy[0];
	    ps = data->trades[i].price_sell[0];
	    npb = data->trades[i].price_buy[0]/data->assets[ib].base_value;
	    nps = data->trades[i].price_sell[0]/data->assets[is].base_value;
	    printf("%2d: Sell  %-4s %d at %7.3f | ",i,
			data->assets[is].asset_name,ns,ps);
	    printf("Buy %-4s %d at %7.3f",
			data->assets[ib].asset_name,nb,pb);
	    printf("  A = %f\n",
			currentAmpRatio(data->trades[i],data,0));
	    for (j = 1; j < data->trades[i].num_stages; ++j)
	    {
	    	ib = data->trades[i].index_buy[j];
	    	is = data->trades[i].index_sell[j];
	    	nb = data->trades[i].num_shares_buy[j];
	    	ns = data->trades[i].num_shares_sell[j];
	    	pb = data->trades[i].price_buy[j];
	    	ps = data->trades[i].price_sell[j];
	    	printf("    Sell  %-4s %d at %7.3f | ",
			data->assets[is].asset_name,ns,ps);
	    	printf("Buy %-4s %d at %7.3f",
			data->assets[ib].asset_name,nb,pb);
	    	printf("  A = %6.4f\n",
			currentAmpRatio(data->trades[i],data,j));
	    }
	}
}	/* PrintOpenTrade */

extern void PrintProfitableTrade(
	DATA_SET *data)
{
	int i,j,ns,is,ib;
	double ratio;
	TRADE trade;

	for (i = 0; i < data->num_trades; ++i)
	{
	    if (data->trades[i].status == CLOSED) continue;
	    ft_assign(&trade,&data->trades[i],sizeof(TRADE));
	    if (data->trades[i].status == PARTIAL_CLOSED)
		trade.num_stages = ns = StartIndex(data->trades[i]);
	    else
	    	ns = data->trades[i].num_stages;
	    for (j = 0; j < ns; ++j)
	    {
	    	is = data->trades[i].index_buy[ns-1];
	    	ib = data->trades[i].index_sell[j];
	    	ratio = currentAmpRatio(trade,data,j);
		if (j > 0) printf("   ");
		else printf("%2d ",i);
	    	printf("Amp factor of open trade %4s(S)/%-4s(B) = %f\n",
				data->assets[is].asset_name,
				data->assets[ib].asset_name,ratio);
	    }
	}
}	/* end PrintProfitableTrade */

extern void PrintClosedTrade(
	FILE *outfile,
	TRADE trade,
	DATA_SET *data)
{
	int *is,*ib;
	int *ns,*nb;
	double *ps,*pb;
	double cash_gain;
	double net_gain;

	printf("Entering PrintClosedTrade()\n");
	is = trade.index_sell;
	ib = trade.index_buy;
	ns = trade.num_shares_sell;
	nb = trade.num_shares_buy;
	ps = trade.price_sell;
	pb = trade.price_buy;
	cash_gain = ns[0]*ps[0] + ns[1]*ps[1] - nb[0]*pb[0] - nb[1]*pb[1];
	net_gain = pb[1]*(nb[1]-ns[0]) + ps[1]*(nb[0]-ns[1]) + cash_gain;

	if (trade.status != CLOSED) return;

	if (outfile == NULL)
	{
	    if (is[0] != ib[1] || is[1] != ib[0])
	    {
	    	printf("Warning: not a matched trade!\n\n");
		clean_up(ERROR);
	    }
	    printf("\n");
	    printf(" Forward trade: Sell  %-4s %5d at %7.3f | ",
                        data->assets[is[0]].asset_name,ns[0],ps[0]);
            printf("Buy %-4s %5d at %7.3f\n",
                        data->assets[ib[0]].asset_name,nb[0],pb[0]);	    
	    printf("Backward trade: Sell  %-4s %5d at %7.3f | ",
                        data->assets[is[1]].asset_name,ns[1],ps[1]);
            printf("Buy %-4s %5d at %7.3f\n",
                        data->assets[ib[1]].asset_name,nb[1],pb[1]);	    
	    printf("Gain: %-4s %d; %-4s %d; Cash: %f Amp-ratio: %f\n",
			data->assets[is[0]].asset_name,
			nb[1]-ns[0],
			data->assets[is[1]].asset_name,
			nb[0]-ns[1],
			cash_gain,
			ps[0]*ps[1]/pb[0]/pb[1]);
	    printf("Net gain at present prices: $%f\n",net_gain);
	}
	else
	{
	    char date[200];
	    std::time_t t = std::time(NULL);
	    strftime(date,100*sizeof(char),"%D",std::localtime(&t));
	    if (is[0] != ib[1] || is[1] != ib[0])
	    {
	    	printf("Warning: not a matched trade!\n\n");
		return;
	    }
	    fprintf(outfile,"Date: %s\n",date);
	    fprintf(outfile," Forward trade: Sell  %-4s %5d at %7.3f | ",
                        data->assets[is[0]].asset_name,ns[0],ps[0]);
            fprintf(outfile,"Buy %-4s %5d at %7.3f\n",
                        data->assets[ib[0]].asset_name,nb[0],pb[0]);	    
	    fprintf(outfile,"Backward trade: Sell  %-4s %5d at %7.3f | ",
                        data->assets[is[1]].asset_name,ns[1],ps[1]);
            fprintf(outfile,"Buy %-4s %5d at %7.3f\n",
                        data->assets[ib[1]].asset_name,nb[1],pb[1]);	    
	    fprintf(outfile,"Gain: %-4s %d; %-4s %d; Cash: %f Amp-ratio: %f\n",
			data->assets[is[0]].asset_name,
			nb[1]-ns[0],
			data->assets[is[1]].asset_name,
			nb[0]-ns[1],
			cash_gain,
			ps[0]*ps[1]/pb[0]/pb[1]);
	    fprintf(outfile,"Net gain at present prices: $%f\n",net_gain);
	    fprintf(outfile,"\n");
	}
}	/* end PrintClosedTrade */

extern void SaveDeleteClosedTrade(
	DATA_SET *data)
{
	char fname[200];
	FILE *sfile;
	int i,j;

	create_directory("record",NO);
	sprintf(fname,"record/%s",data->data_name);
	sfile = fopen(fname,"r");
	if (sfile != NULL)
	{
	    fclose(sfile);
	    sfile = fopen(fname,"a");
	}
	else
	    sfile = fopen(fname,"w");
	for (i = 0; i < data->num_trades; ++i)
	{
	    if (!data->trades[i].status == CLOSED) continue;
	    PrintClosedTrade(sfile,data->trades[i],data);
	    for (j = i+1; j < data->num_trades; ++j)
	    	ft_assign(&data->trades[j-1],&data->trades[j],sizeof(TRADE));

	    data->num_trades--;
	    i--;
	}
	fclose(sfile);
}	/* end SaveDeleteClosedTrade */

static double AmpRatio(
	TRADE trade,
	DATA_SET *data,
	int istart)
{
	double ratio;
	int n = data->num_days-1;
	int ns = trade.num_stages-1;
	int ib = trade.index_sell[istart];
	int is = trade.index_buy[ns];
	int i;

	double st = data->assets[is].value[n];
	double bt = data->assets[ib].value[n];
	ratio = 1.0;
	for (i = istart; i < trade.num_stages; ++i)
	{
	    ratio *= trade.price_sell[i]/trade.price_buy[i];
	}
	ratio *= st/bt;
	return ratio;
}	/* end AmpRatio */

static double currentAmpRatio(
	TRADE trade,
	DATA_SET *data,
	int istart)
{
	TRADE cur_trade;
	int ns = trade.num_stages;
	int nd = data->num_days-1;
	int is,ib;

	ft_assign(&cur_trade,&trade,sizeof(TRADE));
	is = cur_trade.index_sell[ns] = trade.index_buy[ns-1];
	ib = cur_trade.index_buy[ns] = trade.index_sell[istart];
	cur_trade.price_sell[ns] = data->assets[is].value[nd];
	cur_trade.price_buy[ns] = data->assets[ib].value[nd];
	cur_trade.num_stages++;
	return AmpRatio(cur_trade,data,istart);
}	/* end currentAmpRatio */

static int StartIndex(		// Start index of partially closed trade
	TRADE trade)
{
	int i,N = trade.num_stages-1;
	switch (trade.status)
	{
	case OPEN:
	    return -1;
	case CLOSED:
	    return 0;
	case PARTIAL_CLOSED:
	    for (i = 0; i < N; ++i)
	    {
		if (trade.index_sell[i] == trade.index_buy[N])
		    break;
	    }
	    return (i == N) ? -1 : i;
	default:
	    printf("Unknown trade status!\n");
	    clean_up(ERROR);
	}
}	/* end StartIndex */
	
extern void PrintClosedTradeLoop(
	FILE *outfile,
	TRADE trade,
	DATA_SET *data)
{
	int *is,*ib;
	int *ns,*nb;
	int share_gain[MAX_NUM_STAGES];
	int i,istart,iend,inext,iday;
	double *ps,*pb;
	double cash_gain,net_gain,Amp_ratio;
	double current_ib_value;

	if (trade.status == OPEN) return;

	is = trade.index_sell;
	ib = trade.index_buy;
	ns = trade.num_shares_sell;
	nb = trade.num_shares_buy;
	ps = trade.price_sell;
	pb = trade.price_buy;
	istart = StartIndex(trade);
	iend = trade.num_stages - 1;
	iday = data->num_days - 1;

	cash_gain = 0.0;
	net_gain = 0.0;
	Amp_ratio = 1.0;
	for (i = istart; i <= iend; ++i)
	{
	    inext = (i == iend) ? istart : i+1;
	    current_ib_value = data->assets[ib[i]].value[iday];
	    cash_gain += ns[i]*ps[i] - nb[i]*pb[i];
	    share_gain[i] = nb[i] - ns[inext];
	    net_gain += current_ib_value*share_gain[i];
	    Amp_ratio *= ps[i]/pb[i];
	}

	if (outfile == NULL)
	{
	    printf("\nClosed sell loop:\n");
	    for (i = istart; i <= iend; ++i)
	    {
	    	printf("%d: Sell  %-4s %5d at %7.3f | ",i,
                        data->assets[is[i]].asset_name,ns[i],ps[i]);
            	printf("Buy %-4s %5d at %7.3f Gain %5d shares\n",
                        data->assets[ib[i]].asset_name,nb[i],pb[i],
			share_gain[i]);	    
	    }
	    printf("Cash gain: %f Net gain: %f Amp ratio: %f\n",
			cash_gain,net_gain,Amp_ratio);
	}
	else
	{
	    char date[200];
	    std::time_t t = std::time(NULL);
	    strftime(date,100*sizeof(char),"%D",std::localtime(&t));
	    fprintf(outfile,"Date: %s\n",date);
	    fprintf(outfile,"\nClosed sell loop:\n");
	    for (i = istart; i <= iend; ++i)
	    {
	    	fprintf(outfile,"%d: Sell  %-4s %5d at %7.3f | ",i,
                        data->assets[is[i]].asset_name,ns[i],ps[i]);
            	fprintf(outfile,"Buy %-4s %5d at %7.3f Gain %5d shares\n",
                        data->assets[ib[i]].asset_name,nb[i],pb[i],
			share_gain[i]);	    
	    }
	    fprintf(outfile,"Cash gain: %f Net gain: %f Amp ratio: %f\n",
			cash_gain,net_gain,Amp_ratio);
	}
}	/* end PrintPartialClosedTrade */

extern void SaveDeleteClosedTradeLoop(
	DATA_SET *data)
{
	char fname[200];
	FILE *sfile;
	int i,j;

	create_directory("record",NO);
	sprintf(fname,"record/%s",data->data_name);
	sfile = fopen(fname,"r");
	if (sfile != NULL)
	{
	    fclose(sfile);
	    sfile = fopen(fname,"a");
	}
	else
	    sfile = fopen(fname,"w");
	for (i = 0; i < data->num_trades; ++i)
	{
	    if (data->trades[i].status == OPEN) continue;
	    PrintClosedTradeLoop(sfile,data->trades[i],data);
	    if (data->trades[i].status == CLOSED)
	    {
	    	for (j = i+1; j < data->num_trades; ++j)
	    	    ft_assign(&data->trades[j-1],&data->trades[j],
				sizeof(TRADE));
	    	data->num_trades--;
	    	i--;
	    }
	    else if (data->trades[i].status == PARTIAL_CLOSED)
	    {
		data->trades[i].num_stages = StartIndex(data->trades[i]);
		data->trades[i].status = OPEN;
	    }
	}
	fclose(sfile);
}	/* end SaveDeleteClosedTradeLoop */


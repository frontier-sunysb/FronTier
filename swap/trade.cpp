#include "swap.h"

static void PrintOpenTrade(DATA_SET*);
static void PrintClosedTrade(FILE*,TRADE,DATA_SET*);

#define 	MAX_NUM_TRADE		200

extern void InvestShares(
	DATA_SET *data)
{
	char string[100];
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
		printf("Add number of shares to %s: ",
				data->assets[i].asset_name);
		scanf("%d",data->shares+i);
	    }
	}
}	/* end InvestShares */

extern void TradeShares(
	DATA_SET *data)
{
	boolean to_trade = YES;
	int i,n;
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
			"and price to buy: ");
		scanf("%d %d %lf",&data->trades[n].index_buy[0],
				&data->trades[n].num_shares_buy[0],
				&data->trades[n].price_buy[0]);
	    	printf("Enter equity index, number of shares "
			"and price to sell: ");
		scanf("%d %d %lf",&data->trades[n].index_sell[0],
				&data->trades[n].num_shares_sell[0],
				&data->trades[n].price_sell[0]);
		data->trades[n].closed = NO;
		(data->num_trades)++;
	    }
	    else
	    {
		PrintOpenTrade(data);
		printf("Select index of open trade: ");
		scanf("%d",&n);
	    	printf("Enter equity index, number of shares "
			"and price to buy: ");
		scanf("%d %d %lf",&data->trades[n].index_buy[1],
				&data->trades[n].num_shares_buy[1],
				&data->trades[n].price_buy[1]);
	    	printf("Enter equity index, number of shares "
			"and price to sell: ");
		scanf("%d %d %lf",&data->trades[n].index_sell[1],
				&data->trades[n].num_shares_sell[1],
				&data->trades[n].price_sell[1]);
		if (data->trades[n].index_buy[0] != 
		    data->trades[n].index_sell[1] ||
		    data->trades[n].index_buy[1] !=
                    data->trades[n].index_sell[0])
		{
		    printf("Invalid matching trade!\n");
		    continue;
		}
		data->trades[n].closed = YES;
	    }
	}
}	/* end TradeShares */

static void PrintOpenTrade(
	DATA_SET *data)
{
	int i;
	int ib,is;
	int nb,ns;
	double pb,ps;
	double npb,nps;
	for (i = 0; i < data->num_trades; ++i)
	{
	    if (data->trades[i].closed) continue;
	    ib = data->trades[i].index_buy[0];
	    is = data->trades[i].index_sell[0];
	    nb = data->trades[i].num_shares_buy[0];
	    ns = data->trades[i].num_shares_sell[0];
	    pb = data->trades[i].price_buy[0];
	    ps = data->trades[i].price_sell[0];
	    npb = data->trades[i].price_buy[0]/data->assets[ib].base_value;
	    nps = data->trades[i].price_sell[0]/data->assets[is].base_value;
	    printf("Trade %2d: Sell  %-4s %d at %7.3f | ",i,
			data->assets[is].asset_name,ns,ps);
	    printf("Buy %-4s %d at %7.3f\n",
			data->assets[ib].asset_name,nb,pb);
	    if (npb < nps)
		printf("\t  %-4s > %-4s %6.2f percent\n",
				data->assets[is].asset_name,
				data->assets[ib].asset_name,100*(nps-npb));
	    else
		printf("\t  %-4s < %-4s %6.2f percent\n",
				data->assets[is].asset_name,
				data->assets[ib].asset_name,100*(npb-nps));
	}
}	/* PrintOpenTrade */

extern void PrintProfitableTrade(
	DATA_SET *data)
{
	TRADE trade;
	int i,is,ib;
	int n = data->num_days - 1;
	double ps[2],pb[2];
	double ratio;

	for (i = 0; i < data->num_trades; ++i)
	{
	    if (data->trades[i].closed) continue;
	    ft_assign(&trade,&data->trades[i],sizeof(TRADE));
	    trade.closed = YES;
	    ib = data->trades[i].index_buy[0];
	    is = data->trades[i].index_sell[0];
	    trade.index_buy[1] = is;
	    trade.index_sell[1] = ib;
	    pb[0] = trade.price_buy[0];
	    ps[0] = trade.price_sell[0];
	    pb[1] = trade.price_buy[1] = data->assets[is].value[n];
	    ps[1] = trade.price_sell[1] = data->assets[ib].value[n];
	    trade.num_shares_sell[1] = trade.num_shares_buy[0];
	    trade.num_shares_buy[1] = (int)trade.num_shares_sell[1]*ps[1]/pb[1];
	    ratio = ps[0]*ps[1]/pb[0]/pb[1];
	    printf("Amplification factor of open trade %4s/%-4s= %f\n",
				data->assets[is].asset_name,
				data->assets[ib].asset_name,ratio);
	    if (ratio > 1.0)
	    	PrintClosedTrade(NULL,trade,data);
	}
}	/* end PrintProfitableTrade */

static void PrintClosedTrade(
	FILE *outfile,
	TRADE trade,
	DATA_SET *data)
{
	int *is,*ib;
	int *ns,*nb;
	double *ps,*pb;
	is = trade.index_sell;
	ib = trade.index_buy;
	ns = trade.num_shares_sell;
	nb = trade.num_shares_buy;
	ps = trade.price_sell;
	pb = trade.price_buy;

	if (outfile == NULL)
	{
	    printf(" Forward trade: Sell  %-4s %d at %7.3f | ",
                        data->assets[is[0]].asset_name,ns[0],ps[0]);
            printf("Buy %-4s %d at %7.3f\n",
                        data->assets[ib[0]].asset_name,nb[0],pb[0]);	    
	    printf("Backward trade: Sell  %-4s %d at %7.3f | ",
                        data->assets[is[1]].asset_name,ns[1],ps[1]);
            printf("Buy %-4s %d at %7.3f\n",
                        data->assets[ib[1]].asset_name,nb[1],pb[1]);	    
	    if (is[0] != ib[1] || is[1] != ib[0])
	    {
	    	printf("Warning: not a matched trade!\n\n");
		return;
	    }
	    printf("Gain: %-4s %d; %-4s %d; Cash: %f Amp-ratio: %f\n",
			data->assets[is[0]].asset_name,
			nb[1]-ns[0],
			data->assets[is[1]].asset_name,
			nb[0]-ns[1],
			ns[0]*ps[0]+ns[1]*ps[1]-nb[0]*pb[0]-nb[1]*pb[1],
			ps[0]*ps[1]/pb[0]/pb[1]);
	}
	else
	{
	    if (is[0] != ib[1] || is[1] != ib[0])
	    {
	    	printf("Warning: not a matched trade!\n\n");
		return;
	    }
	    fprintf(outfile," Forward trade: Sell  %-4s %d at %7.3f | ",
                        data->assets[is[0]].asset_name,ns[0],ps[0]);
            fprintf(outfile,"Buy %-4s %d at %7.3f\n",
                        data->assets[ib[0]].asset_name,nb[0],pb[0]);	    
	    fprintf(outfile,"Backward trade: Sell  %-4s %d at %7.3f | ",
                        data->assets[is[1]].asset_name,ns[1],ps[1]);
            fprintf(outfile,"Buy %-4s %d at %7.3f\n",
                        data->assets[ib[1]].asset_name,nb[1],pb[1]);	    
	    fprintf(outfile,"Gain: %-4s %d; %-4s %d; Cash: %f Amp-ratio: %f\n",
			data->assets[is[0]].asset_name,
			nb[1]-ns[0],
			data->assets[is[1]].asset_name,
			nb[0]-ns[1],
			ns[0]*ps[0]+ns[1]*ps[1]-nb[0]*pb[0]+nb[1]*pb[1],
			ps[0]*ps[1]/pb[0]/pb[1]);
	}
}	/* end PrintClosedTrade */

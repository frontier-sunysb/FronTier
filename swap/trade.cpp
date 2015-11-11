#include "swap.h"
#include <ctime>
#include <iostream>
#include <locale>

static double AmpRatio(TRADE,PORTFOLIO*,int);
static double currentAmpRatio(TRADE,PORTFOLIO*,int);
static int StartIndex(TRADE);
static void SplitTrade(PORTFOLIO*);
static void MergeTrades(PORTFOLIO*);
static void WrapMiddleTrades(PORTFOLIO*);
static void WrapTrade(TRADE*);

#define 	MAX_NUM_TRADE		200

extern void InvestShares(
	PORTFOLIO *account)
{
	MARKET_DATA *data = account->data;
	char string[100];
	int num_shares;
	int i;

	if (account->shares == NULL)
	{
	    FT_VectorMemoryAlloc((POINTER*)&account->shares,data->num_assets,
					sizeof(int));
	    for (i = 0; i < data->num_assets; ++i)
		account->shares[i] = 0;
	}
	for (i = 0; i < data->num_assets; ++i)
	    printf("%-4s current share number %d\n",
			data->assets[i].asset_name,account->shares[i]);
	printf("Type yes to add investment: ");
	scanf("%s",string);
	if (string[0] == 'y' || string[0] == 'Y')
	{
	    for (i = 0; i < data->num_assets; ++i)
	    {
		printf("%-4s add: ",data->assets[i].asset_name);
		scanf("%d",&num_shares);
		account->shares[i] += num_shares;
	    }
	    printf("Investment update:\n");
	    for (i = 0; i < data->num_assets; ++i)
		printf("%-4s: %d\n",data->assets[i].asset_name,
				account->shares[i]);
	}
}	/* end InvestShares */

extern boolean TradeShares(
	PORTFOLIO *account)
{
	boolean to_trade = YES;
	int i,n,is,ib,ns,nb,istart;
	char string[100];
	double ratio,ps,pb;
	boolean traded = NO;
	double prices[2];
	char *stock_names[2];
	MARKET_DATA *data = account->data;

	while (to_trade)
	{
	    if (account->num_trades >= MAX_NUM_TRADE)
	    {
		printf("There are too many open trades!\n");
		break;
	    }
	    printf("Type yes to print liear profile: ");
	    scanf("%s",string);
	    if (string[0] == 'y' || string[0] == 'Y')
		PrintCurrentLinearProfile(account);
	    printf("Type yes if this is a new trade: ");
	    scanf("%s",string);
	    if (string[0] == 'y' || string[0] == 'Y')
	    {
		double total_cash;
		n = account->num_trades;
		for (i = 0; i < data->num_assets; ++i)
		{
	    	    printf("Index of assets are: %d for %4s\n",i,
			    data->assets[i].asset_name);
		}
	    	printf("Enter equity indices for sell and buy: ");
		scanf("%d %d",&is,&ib);
		printf("Enter total cash value of the trade: ");
		scanf("%lf",&total_cash);
		stock_names[0] = data->assets[is].asset_name;
		stock_names[1] = data->assets[ib].asset_name;
		GetCurrentPrice(stock_names,prices,2);
		ps = prices[0];
		pb = prices[1];
		printf("Number of shares to sell and buy: %d %d\n",
				irint(total_cash/ps),irint(total_cash/pb));
	    	printf("Enter realized sell shares and price for %s: ",
				stock_names[0]);
		scanf("%d %lf",&account->trades[n].num_shares_sell[0],
				&account->trades[n].price_sell[0]);
	    	printf("Enter realized buy shares and price for %s: ",
				stock_names[1]);
		scanf("%d %lf",&account->trades[n].num_shares_buy[0],
				&account->trades[n].price_buy[0]);
		account->trades[n].status = OPEN;
		account->trades[n].index_sell[0] = is;
		account->trades[n].index_buy[0] = ib;
		ns = account->trades[n].num_shares_sell[0];
		nb = account->trades[n].num_shares_buy[0];
		account->shares[is] -= ns;
		account->shares[ib] += nb;
		(account->num_trades)++;
		account->trades[n].num_stages++;
		traded = YES;
	    }
	    else
	    {
		int nstage,sindex;
		PrintAssetList(data);
		printf("Enter sell index: ");
		scanf("%d",&sindex);
		PrintOpenTradeForIndex(account,sindex);
		printf("Select number of open trade: ");
		scanf("%d",&n);
		nstage = account->trades[n].num_stages;
		account->trades[n].index_sell[nstage] = sindex;
		printf("Enter buy index: ");
		scanf("%d",&account->trades[n].index_buy[nstage]);  
		if (account->trades[n].index_buy[nstage-1] != 
		    account->trades[n].index_sell[nstage])
		{
		    printf("Invalid matching trade!\n");
		    continue;
		}
		is = account->trades[n].index_sell[nstage];
		ib = account->trades[n].index_buy[nstage];
		stock_names[0] = data->assets[is].asset_name;
		stock_names[1] = data->assets[ib].asset_name;
		istart = -1;
		for (i = 0; i < nstage; ++i)
		{
		    if (account->trades[n].index_sell[i] == 
			account->trades[n].index_buy[nstage])
		    {
			istart = i;
			break;
		    }
		}
		ratio = 1.0;
		GetCurrentPrice(stock_names,prices,2);
		ps = prices[0];
		pb = prices[1];
		ratio *= ps/pb;
		if (istart != -1)
		{
		    const char *gain_loss;
		    for (i = istart; i < account->trades[n].num_stages; ++i)
        	    {
        	    	ratio *= account->trades[n].price_sell[i]/
				account->trades[n].price_buy[i];
        	    }
		    printf("Amp ratio: %f\n",ratio);
		    gain_loss = (ratio > 1.0) ? "Gain" : "Loss";
		    printf("Number of shares to sell and buy:\n");
		    printf("Option 1: %d %d  %s in %s\n",
			account->trades[n].num_shares_buy[nstage-1],
			irint(account->trades[n].num_shares_sell[istart]*ratio),
			gain_loss,data->assets[ib].asset_name);
		    printf("Option 2: %d %d  %s in %s\n",
			irint(account->trades[n].num_shares_buy[nstage-1]/ratio),
			account->trades[n].num_shares_sell[istart],
			gain_loss,data->assets[is].asset_name);
		    printf("Option 3: %d %d  %s in Cash\n",
			account->trades[n].num_shares_buy[nstage-1],
			account->trades[n].num_shares_sell[istart],gain_loss);
	    	    printf("Trade? ");
	    	    scanf("%s",string);
	    	    if (string[0] != 'y' && string[0] != 'Y')
			break;
		}
		else
		{
		    printf("Number of shares to sell and buy: %d %d\n",
			account->trades[n].num_shares_buy[nstage-1],
			irint(account->trades[n].num_shares_buy[nstage-1]*ratio));
		}
		is = account->trades[n].index_sell[nstage];
	    	printf("Enter realized shares "
			"and price sold for %s: ",stock_names[0]);
		scanf("%d %lf",&account->trades[n].num_shares_sell[nstage],
				&account->trades[n].price_sell[nstage]);
	    	printf("Enter realized shares "
			"and price bought for %s: ",stock_names[1]);
		scanf("%d %lf",&account->trades[n].num_shares_buy[nstage],
				&account->trades[n].price_buy[nstage]);
		ns = account->trades[n].num_shares_sell[nstage];
		ib = account->trades[n].index_buy[nstage];
		nb = account->trades[n].num_shares_buy[nstage];
		account->shares[is] -= ns;
		account->shares[ib] += nb;
		if (account->trades[n].index_buy[nstage] ==
		    account->trades[n].index_sell[0])
		{
		    printf("Loop is fully closed, type open to keep it open: ");
		    scanf("%s",string);
		    if (string[0] != 'o')
		    	account->trades[n].status = CLOSED;
		}
		else
		{
		    for (i = 1; i < account->trades[n].num_stages; ++i)
		    {
			if (account->trades[n].index_buy[nstage] ==
                    	    account->trades[n].index_sell[i])
			{
		    	    account->trades[n].status = PARTIAL_CLOSED;
			    break;
			}
		    }
		}
		account->trades[n].num_stages++;
		traded = YES;
	    }
	    printf("Type yes to continue trading: ");
	    scanf("%s",string);
	    if (string[0] != 'y' && string[0] != 'Y')
		break;
	}
	if (!traded) return NO;
	printf("Investment update:\n");
	for (i = 0; i < data->num_assets; ++i)
	    printf("%-4s: %d\n",data->assets[i].asset_name,account->shares[i]);
	return YES;
}	/* end TradeShares */

extern void PrintOpenTrade(
	PORTFOLIO *account)
{
	int i,j,is,ib,ns,nsell,nbuy;
	double ratio;
	TRADE trade;
	char **stock_names;
	double *prices;
	MARKET_DATA *data = account->data;
	int M = data->num_assets;
	int n = data->num_segs;

	FT_VectorMemoryAlloc((POINTER*)&stock_names,M,sizeof(char*));
	FT_VectorMemoryAlloc((POINTER*)&prices,M,sizeof(double));
	for (i = 0; i < M; ++i)
	    stock_names[i] = data->assets[i].asset_name;
	GetCurrentPrice(stock_names,prices,M);
	for (i = 0; i < M; ++i)
	    data->assets[i].price[n] = prices[i];
	data->num_segs++;

	for (i = 0; i < account->num_trades; ++i)
	{
	    if (account->trades[i].status == CLOSED) continue;
	    ft_assign(&trade,&account->trades[i],sizeof(TRADE));
	    if (account->trades[i].status == PARTIAL_CLOSED)
		trade.num_stages = ns = StartIndex(account->trades[i]);
	    else
	    	ns = account->trades[i].num_stages;
	    for (j = 0; j < ns; ++j)
	    {
	    	is = account->trades[i].index_buy[ns-1];
	    	nsell = account->trades[i].num_shares_buy[ns-1];
	    	ib = account->trades[i].index_sell[j];
	    	nbuy = account->trades[i].num_shares_sell[j];
	    	ratio = currentAmpRatio(trade,account,j);
		if (j > 0) printf("   ");
		else printf("%2d ",i);
	    	printf("Open trade (%d/%d) (Sell)%4s/%-4s (%4d/%4d) Amp = %f\n",
				is,ib,data->assets[is].asset_name,
				data->assets[ib].asset_name,nsell,nbuy,ratio);
	    }
	}
	FT_FreeThese(2,stock_names,prices);
	data->num_segs--;
}	/* end PrintOpenTrade */

static double AmpRatio(
	TRADE trade,
	PORTFOLIO *account,
	int istart)
{
	MARKET_DATA *data = account->data;
	double ratio;
	int n = data->num_segs-1;
	int ns = trade.num_stages-1;
	int ib = trade.index_sell[istart];
	int is = trade.index_buy[ns];
	int i;

	double st = data->assets[is].price[n];
	double bt = data->assets[ib].price[n];
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
	PORTFOLIO *account,
	int istart)
{
	MARKET_DATA *data = account->data;
	TRADE cur_trade;
	int ns = trade.num_stages;
	int nd = data->num_segs-1;
	int is,ib;
	char *stock_names[2];
	double prices[2];

	ft_assign(&cur_trade,&trade,sizeof(TRADE));
	is = cur_trade.index_sell[ns] = trade.index_buy[ns-1];
	ib = cur_trade.index_buy[ns] = trade.index_sell[istart];
	cur_trade.price_sell[ns] = data->assets[is].price[nd];
	cur_trade.price_buy[ns] = data->assets[ib].price[nd];
	cur_trade.num_stages++;
	return AmpRatio(cur_trade,account,istart);
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
	    for (i = 1; i < N; ++i)
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
	
extern double PrintClosedTradeLoop(
	FILE *outfile,
	TRADE trade,
	PORTFOLIO *account)
{
	MARKET_DATA *data = account->data;
	int *is,*ib;
	int *ns,*nb;
	int share_gain[MAX_NUM_STAGES];
	int i,istart,iend,inext,iday;
	double *ps,*pb;
	double cash_gain,net_gain,Amp_ratio;
	double current_ib_value;
	double fee = 0.0;			// trading cost

	if (trade.status == OPEN) return 0.0;

	is = trade.index_sell;
	ib = trade.index_buy;
	ns = trade.num_shares_sell;
	nb = trade.num_shares_buy;
	ps = trade.price_sell;
	pb = trade.price_buy;
	istart = StartIndex(trade);
	iend = trade.num_stages - 1;
	iday = data->num_segs - 1;

	cash_gain = 0.0;
	net_gain = 0.0;
	Amp_ratio = 1.0;
	for (i = istart; i <= iend; ++i)
	{
	    inext = (i == iend) ? istart : i+1;
	    current_ib_value = data->assets[ib[i]].price[iday];
	    cash_gain += ns[i]*ps[i] - nb[i]*pb[i];
	    share_gain[i] = nb[i] - ns[inext];
	    net_gain += current_ib_value*share_gain[i];
	    Amp_ratio *= ps[i]/pb[i];
	    fee += 20.0;
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
	    printf("Cash gain: %7.2f Share gain: %7.2f Total: %7.2f Amp: %f\n",
			cash_gain,net_gain,cash_gain+net_gain,Amp_ratio);
	}
	else
	{
	    char date[200];
	    std::time_t t = std::time(NULL);
	    strftime(date,100*sizeof(char),"%D",std::localtime(&t));
	    fprintf(outfile,"\nDate: %s\n",date);
	    fprintf(outfile,"Closed sell loop:\n");
	    for (i = istart; i <= iend; ++i)
	    {
	    	fprintf(outfile,"%d: Sell  %-4s %5d at %7.3f | ",i,
                        data->assets[is[i]].asset_name,ns[i],ps[i]);
            	fprintf(outfile,"Buy %-4s %5d at %7.3f Gain %5d shares\n",
                        data->assets[ib[i]].asset_name,nb[i],pb[i],
			share_gain[i]);	    
	    }
	    fprintf(outfile,"Cash gain: %7.2f Share gain: %7.2f ",
				cash_gain,net_gain);
	    fprintf(outfile,"Total: %7.2f Amp: %5.3f\n",
				cash_gain+net_gain-fee,Amp_ratio);
	}
	return cash_gain+net_gain-fee;
}	/* end PrintPartialClosedTrade */

extern void SaveDeleteClosedTradeLoop(
	PORTFOLIO *account)
{
	MARKET_DATA *data = account->data;
	char fname[200];
	FILE *sfile;
	int i,j;

	create_directory("record",NO);
	sprintf(fname,"record/%s",account->account_name);
	sfile = fopen(fname,"r");
	if (sfile != NULL)
	{
	    fclose(sfile);
	    sfile = fopen(fname,"a");
	}
	else
	    sfile = fopen(fname,"w");
	for (i = 0; i < account->num_trades; ++i)
	{
	    if (account->trades[i].status == OPEN) continue;
	    PrintClosedTradeLoop(sfile,account->trades[i],account);
	    if (account->trades[i].status == CLOSED)
	    {
	    	for (j = i+1; j < account->num_trades; ++j)
	    	    ft_assign(&account->trades[j-1],&account->trades[j],
				sizeof(TRADE));
	    	account->num_trades--;
	    	i--;
	    }
	    else if (account->trades[i].status == PARTIAL_CLOSED)
	    {
		account->trades[i].num_stages = StartIndex(account->trades[i]);
		account->trades[i].status = OPEN;
	    }
	}
	fclose(sfile);
}	/* end SaveDeleteClosedTradeLoop */

extern void SortTradeOrder(
	PORTFOLIO *account)
{
	int i,j;
	TRADE tmp_trade;

	if (account->num_trades == 0) return;
	for (i = 0; i < account->num_trades-1; ++i)
	{
            for (j = i+1; j < account->num_trades; ++j)
            {
                if (currentAmpRatio(account->trades[i],account,0) >
                    currentAmpRatio(account->trades[j],account,0))
                {
                    ft_assign(&tmp_trade,&account->trades[i],sizeof(TRADE));
                    ft_assign(&account->trades[i],&account->trades[j],
						sizeof(TRADE));
                    ft_assign(&account->trades[j],&tmp_trade,sizeof(TRADE));
                }
            }
	}
}	/* end SortTradeOrder */

extern void WrapPartialTradeLoop(
	PORTFOLIO *account)
{
	char fname[200];
	FILE *sfile;
	int i,j,is,ie;
	double last_cost,profit;

	for (i = 0; i < account->num_trades; ++i)
	{
	    if (account->trades[i].status == OPEN ||
		account->trades[i].status == CLOSED) 
		continue;
	    profit = PrintClosedTradeLoop(sfile,account->trades[i],account);
	    is = StartIndex(account->trades[i]);
	    ie = account->trades[i].num_stages - 1;
	    last_cost = account->trades[i].num_shares_buy[is-1]*
			account->trades[i].price_buy[is-1];
	    last_cost -= profit;
	    account->trades[i].num_shares_buy[is-1] = 
			account->trades[i].num_shares_buy[ie];
	    account->trades[i].price_buy[is-1] = last_cost/
			account->trades[i].num_shares_buy[is-1];
	    account->trades[i].num_stages = StartIndex(account->trades[i]);
	    account->trades[i].status = OPEN;
	}
}	/* end WrapPartialTradeLoop */

extern boolean ExperimentTrade(
	PORTFOLIO *account)
{
	MARKET_DATA *data = account->data;
	boolean to_trade = YES;
	int i,n,is,ib,ns,nb,istart;
	int iday = data->num_segs-1;
	char string[100];
	double ratio,ps,pb;
	boolean traded = NO;

	PrintCurrentLinearProfile(account);
	if (strcmp(account->account_name,"pswap-data") == 0 ||
	    strcmp(account->account_name,"etrade-data") == 0)
	{
	    printf("Cannot use real account to experiment!\n");
	    return NO;
	}
	printf("Trade on current prices? ");
    	scanf("%s",string);
	if (string[0] != 'y' && string[0] != 'Y')
	    return NO;

	while (to_trade)
	{
	    if (account->num_trades >= MAX_NUM_TRADE)
	    {
		printf("There are too many open trades!\n");
		break;
	    }
	    PrintOpenTrade(account);
	    printf("Type yes if this is a new trade: ");
	    scanf("%s",string);
	    if (string[0] == 'y' || string[0] == 'Y')
	    {
		double total_cash;
		n = account->num_trades;
		for (i = 0; i < data->num_assets; ++i)
		{
	    	    printf("Index of assets are: %d for %4s\n",i,
			    data->assets[i].asset_name);
		}
	    	printf("Enter equity indices for sell and buy: ");
		scanf("%d %d",&is,&ib);
		account->trades[n].index_sell[0] = is;
		account->trades[n].index_buy[0] = ib;
		printf("Enter total cash value of the trade: ");
		scanf("%lf",&total_cash);
		ps = data->assets[is].price[iday];
		pb = data->assets[ib].price[iday];
		printf("Number of shares to sell and buy: %d %d\n",
				irint(total_cash/ps),irint(total_cash/pb));
		account->trades[n].num_shares_sell[0] = irint(total_cash/ps);
		account->trades[n].num_shares_buy[0] = irint(total_cash/pb);
		account->trades[n].status = OPEN;
		account->trades[n].price_sell[0] = ps;
		account->trades[n].price_buy[0] = pb;
		account->shares[is] -= account->trades[n].num_shares_sell[0];
		account->shares[ib] += account->trades[n].num_shares_buy[0];
		(account->num_trades)++;
		account->trades[n].num_stages++;
		traded = YES;
	    }
	    else
	    {
		int nstage,sindex;
		PrintAssetList(data);
		printf("Enter sell index: ");
		scanf("%d",&sindex);
		PrintOpenTradeForIndex(account,sindex);
		printf("Select number of open trade: ");
		scanf("%d",&n);
		nstage = account->trades[n].num_stages;
		account->trades[n].index_sell[nstage] = sindex;
		printf("Enter buy index: ");
		scanf("%d",&account->trades[n].index_buy[nstage]);  
		if (account->trades[n].index_buy[nstage-1] != 
		    account->trades[n].index_sell[nstage])
		{
		    printf("Invalid matching trade!\n");
		    continue;
		}
		is = account->trades[n].index_sell[nstage];
		ib = account->trades[n].index_buy[nstage];
		istart = -1;
		for (i = 0; i < nstage; ++i)
		{
		    if (account->trades[n].index_sell[i] == 
			account->trades[n].index_buy[nstage])
		    {
			istart = i;
			break;
		    }
		}
		ratio = 1.0;
		ps = data->assets[is].price[iday];
		pb = data->assets[ib].price[iday];
		ratio *= ps/pb;
		printf("ratio = %f\n",ratio);
		if (istart != -1)
		{
		    const char *gain_loss;
		    for (i = istart; i < account->trades[n].num_stages; ++i)
        	    {
        	    	ratio *= account->trades[n].price_sell[i]/
				account->trades[n].price_buy[i];
        	    }
		    printf("Amp ratio: %f\n",ratio);
		    gain_loss = (ratio > 1.0) ? "Gain" : "Loss";
		    printf("Number of shares to sell and buy:\n");
		    printf("Option 1: %d %d  %s in %s\n",
			account->trades[n].num_shares_buy[nstage-1],
			irint(account->trades[n].num_shares_sell[istart]*ratio),
			gain_loss,data->assets[ib].asset_name);
		    printf("Option 2: %d %d  %s in %s\n",
			irint(account->trades[n].num_shares_buy[nstage-1]/ratio),
			account->trades[n].num_shares_sell[istart],
			gain_loss,data->assets[is].asset_name);
		    printf("Option 3: %d %d  %s in Cash\n",
			account->trades[n].num_shares_buy[nstage-1],
			account->trades[n].num_shares_sell[istart],gain_loss);
	    	    printf("Trade? ");
	    	    scanf("%s",string);
	    	    if (string[0] != 'y' && string[0] != 'Y')
			break;
		}
		else
		{
		    printf("Number of shares to sell and buy: %d %d\n",
			account->trades[n].num_shares_buy[nstage-1],
			irint(account->trades[n].num_shares_buy[nstage-1]*ratio));
		}
		is = account->trades[n].index_sell[nstage];
		ib = account->trades[n].index_buy[nstage];
	    	printf("Enter number of shares to sell and buy: ");
		scanf("%d %d",&ns,&nb);
		account->trades[n].price_sell[nstage] = ps;
		account->trades[n].price_buy[nstage] = pb;
		account->trades[n].num_shares_sell[nstage] = ns;
		account->trades[n].num_shares_buy[nstage] = nb;
		account->shares[is] -= ns;
		account->shares[ib] += nb;
		if (account->trades[n].index_buy[nstage] ==
		    account->trades[n].index_sell[0])
		{
		    printf("Loop is closed, type open to keep it open: ");
		    scanf("%s",string);
		    if (string[0] != 'o')
		    	account->trades[n].status = CLOSED;
		}
		else
		{
		    for (i = 1; i < account->trades[n].num_stages; ++i)
		    {
			if (account->trades[n].index_buy[nstage] ==
                    	    account->trades[n].index_sell[i])
			{
		    	    account->trades[n].status = PARTIAL_CLOSED;
			    break;
			}
		    }
		}
		account->trades[n].num_stages++;
		traded = YES;
	    }
	    printf("Type yes to continue trading: ");
	    scanf("%s",string);
	    if (string[0] != 'y' && string[0] != 'Y')
		break;
	    printf("New profile:\n");
	    PrintCurrentLinearProfile(account);
	}
	if (!traded) return NO;
	printf("Investment update:\n");
	for (i = 0; i < data->num_assets; ++i)
	    printf("%-4s: %d\n",data->assets[i].asset_name,account->shares[i]);
	return YES;
}	/* end TradeShares */

extern void PrintOpenTradeForIndex(
	PORTFOLIO *account,
	int index)
{
	MARKET_DATA *data = account->data;
	int i,j,is,ib,ns,nsell,nbuy;
	double ratio;
	TRADE trade;
	char **stock_names;
	double *prices;
	int M = data->num_assets;
	int n = data->num_segs;

	FT_VectorMemoryAlloc((POINTER*)&stock_names,M,sizeof(char*));
	FT_VectorMemoryAlloc((POINTER*)&prices,M,sizeof(double));
	for (i = 0; i < M; ++i)
	    stock_names[i] = data->assets[i].asset_name;
	GetCurrentPrice(stock_names,prices,M);
	for (i = 0; i < M; ++i)
	    data->assets[i].price[n] = prices[i];
	data->num_segs++;

	for (i = 0; i < account->num_trades; ++i)
	{
	    if (account->trades[i].status != OPEN) continue;
	    ft_assign(&trade,&account->trades[i],sizeof(TRADE));
	    ns = account->trades[i].num_stages;
	    is = account->trades[i].index_buy[ns-1];
	    if (is != index) continue;
	    for (j = 0; j < ns; ++j)
	    {
	    	nsell = account->trades[i].num_shares_buy[ns-1];
	    	ib = account->trades[i].index_sell[j];
	    	nbuy = account->trades[i].num_shares_sell[j];
	    	ratio = currentAmpRatio(trade,account,j);
		if (j > 0) printf("   ");
		else printf("%2d ",i);
	    	printf("Open trade (%d/%d) (Sell)%4s/%-4s (%4d/%4d) Amp = %f\n",
				is,ib,data->assets[is].asset_name,
				data->assets[ib].asset_name,nsell,nbuy,ratio);
	    }
	}
	FT_FreeThese(2,stock_names,prices);
	data->num_segs--;
}	/* end PrintOpenTrade */

extern void ReOrganizeTrade(
	PORTFOLIO *account)
{
	char string[100];
	printf("Available operations are: \n");
	printf("\tSplit trade (s)\n");
	printf("\tMerge trades (m)\n");
	printf("\tWrap up middle trades (w)\n");
	printf("Enter choice: ");
	scanf("%s",string);
	switch (string[0])
	{
	case 's':
	    SplitTrade(account);
	    return;
	case 'm':
	    MergeTrades(account);
	    return;
	case 'w':
	    WrapMiddleTrades(account);
	    return;
	default: 
	    printf("Unknown option!\n");
	    return;
	}
}	/* end FragmentTrade */

static void WrapMiddleTrades(
	PORTFOLIO *account)
{
	int i;
	PrintOpenTrade(account);
	printf("Selcte index of trade to wrap up: ");
	scanf("%d",&i);
	WrapTrade(&account->trades[i]);
}	/* end WrapMiddleTrades */	

static void WrapTrade(
	TRADE *trade)
{
	int i,nstage;
	double gain,last_buy;
	
	nstage = trade->num_stages;
	for (i = 0; i < nstage-1; ++i)
	{
	    gain +=  trade->num_shares_sell[i+1]*trade->price_sell[i+1] -
			trade->num_shares_buy[i]*trade->price_buy[i];
	    gain -= 20.0;	// trading cost
	}
	last_buy = trade->price_buy[nstage-1]*trade->num_shares_buy[nstage-1];
	last_buy -= gain;
	trade->index_buy[0] = trade->index_buy[nstage-1];
	trade->num_shares_buy[0] = trade->num_shares_buy[nstage-1];
	trade->price_buy[0] = last_buy/trade->num_shares_buy[0];
	trade->num_stages = 1;
}	/* end WrapTrade */

static void MergeTrades(
	PORTFOLIO *account)
{
	int i,j,i1,i2;
	int is1,ib1;
	int is2,ib2;
	int nstage1,nstage2;
	int num_shares;
	double buy_cost;

	PrintOpenTrade(account);
	printf("Selcte two indices of trades to merge: ");
	scanf("%d %d",&i1,&i2);
	nstage1 = account->trades[i1].num_stages;
	nstage2 = account->trades[i2].num_stages;

	is1 = account->trades[i1].index_sell[0];
	ib1 = account->trades[i1].index_buy[nstage1-1];
	is2 = account->trades[i2].index_sell[0];
	ib2 = account->trades[i2].index_buy[nstage2-1];
	if (is1 != is2 || ib1 != ib2)
	{
	    printf("Unmatched trades, cannot merge!\n");
	    return;
	}
	WrapTrade(&account->trades[i1]);
	WrapTrade(&account->trades[i2]);
	num_shares = account->trades[i1].num_shares_buy[0] +
			account->trades[i2].num_shares_buy[0];
	buy_cost = account->trades[i1].num_shares_buy[0]*
			account->trades[i1].price_buy[0] +
			account->trades[i2].num_shares_buy[0]*
			account->trades[i2].price_buy[0];
	account->trades[i1].num_shares_buy[0] = num_shares;
	account->trades[i1].price_buy[0] = buy_cost/num_shares;
	for (i = i2+1; i < account->num_trades; ++i)
	    ft_assign(&account->trades[i-1],&account->trades[i],sizeof(TRADE));
	account->num_trades--;
}	/* end MergeTrades */

static void SplitTrade(
	PORTFOLIO *account)
{
	int i,ns,index;
	double frac;
	TRADE *new_trade;

	PrintOpenTrade(account);
	printf("Selcte the index of trade to slit: ");
	scanf("%d",&index);
	printf("Enter fraction of new trade: ");
	scanf("%lf",&frac);
	if (frac > 0.7)
	{
	    printf("Fraction number too large, should be less than 0.7\n");
	    return;
	}
	new_trade = &account->trades[account->num_trades];
	new_trade->status = account->trades[index].status;
	ns = new_trade->num_stages = account->trades[index].num_stages;
	for (i = 0; i < ns; ++i)
	{
	    new_trade->price_buy[i] = account->trades[index].price_buy[i];
	    new_trade->price_sell[i] = account->trades[index].price_sell[i];
	    new_trade->index_buy[i] = account->trades[index].index_buy[i];
	    new_trade->index_sell[i] = account->trades[index].index_sell[i];
	    new_trade->num_shares_buy[i] = 
			irint(frac*account->trades[index].num_shares_buy[i]);
	    new_trade->num_shares_sell[i] = 
			irint(frac*account->trades[index].num_shares_sell[i]);
	    account->trades[index].num_shares_buy[i] -=
			new_trade->num_shares_buy[i];
	    account->trades[index].num_shares_sell[i] -=
			new_trade->num_shares_sell[i];
	}
	account->num_trades++;
}	/* end SplitTrade */

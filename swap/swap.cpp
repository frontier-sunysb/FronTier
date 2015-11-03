#include "swap.h"

int main(int argc, char **argv)
{
	MARKET_DATA *data;
	PORTFOLIO *account;
	boolean input_data = NO;
	boolean get_web_data = NO;
	char string[100];

	FT_Init(argc,argv,NULL);
	FT_ScalarMemoryAlloc((POINTER*)&data,sizeof(MARKET_DATA));
	FT_ScalarMemoryAlloc((POINTER*)&account,sizeof(PORTFOLIO));
	account->data = data;

	argc--;
	argv++;

	while (argc >= 1)
	{
	    if (argv[0][0] != '-')
	    {
	    	printf("Usage: swap [-i data-file] [-a account-file] [-w]\n");
		printf("For: existing input or getting data from web\n");
            	exit(1);
	    }
	    switch (argv[0][1]) {
	    case 'w':
                get_web_data = YES;
		argc -= 1;
                argv += 1;
                break;
	    case 'i':
		input_data = YES;
                zero_scalar(data->data_name,256);
                strcpy(data->data_name,argv[1]);
                argc -= 2;
                argv += 2;
                break;
	    case 'a':
		input_data = YES;
                zero_scalar(account->account_name,256);
                strcpy(account->account_name,argv[1]);
                argc -= 2;
                argv += 2;
                break;
	    }
	}

	if (input_data)
	{
	    ReadMarketData(data);
	    ReadAccountData(account);
	}
	else if (get_web_data)
        {
            GetStockFile(data);
            SaveStockFile(data);
        }
	else
	{
	    (void) printf("No data file: use -i or -w option\n");
	    exit(0);
	}
	data->new_data = NO;

	printf("Available operations are\n");
	printf("\tAdd data (a)\n");
	printf("\tBuy shares (b)\n");
	printf("\tCompare data (c)\n");
	printf("\tExperiment (e)\n");
	printf("\tFragment trade (f)\n");
	printf("\tInitiation (i)\n");
	printf("\tModify data (m)\n");
	printf("\tPrint open trade (o)\n");
	printf("\tRank data (r)\n");
	printf("\tSimulate data (s)\n");
	printf("\tTrade (t)\n");
	printf("\tXgraph data (x)\n");
	printf("\tDo nothing (n)\n");
	CreateJVM();

	while (YES)
	{
	    printf("Enter option: ");
	    scanf("%s",string);
	    switch(string[0])
	    {
	    case 'a':
	    	AddData(data,NO);
		WriteMarketData(data);
    		PromptForDataMap(account);
		break;
	    case 'b':
	    	InvestShares(account);
		WriteAccountData(account);
		break;
	    case 'c':
	    	DataCompare(data);
		closing_out(account);
	    case 'e':
	    	if (!ExperimentTrade(account))
		    closing_out(account);
		WriteAccountData(account);
		break;
	    case 'f':
		FragmentTrade(account);
		WriteAccountData(account);
		closing_out(account);
	    case 'i':
		InitTrade(data);
		closing_out(account);
	    case 'm':
	    	ModifyMarketData(data);
    		PromptForDataMap(account);
		WriteMarketData(data);
		WriteAccountData(account);
		break;
	    case 'o':
	    	TradeInfo(account);
		PrintDataStates(account);
		break;
	    case 'r':
		RankData(data,account->data_map);
		closing_out(account);
	    case 's':
	    	InvestSimulation(data);
		closing_out(account);
	    case 't':
	    	if (!TradeShares(account))
		    closing_out(account);
		WriteAccountData(account);
		break;
	    case 'n':
		closing_out(account);
	    case 'x':
    		XgraphData(data,account->data_map);
		break;
	    default:
		printf("Unknown option\n");
		break;
	    }
	    printf("Type yes to continue: ");
	    scanf("%s",string);
	    if (string[0] != 'y' && string[0] != 'Y')
		break;
	}
	exit(0);
}	/* end main */

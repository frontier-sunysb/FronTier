#include "swap.h"

int main(int argc, char **argv)
{
	DATA_SET *data;
	boolean input_data = NO;
	boolean get_web_data = NO;
	char string[100];

	FT_Init(argc,argv,NULL);
	FT_ScalarMemoryAlloc((POINTER*)&data,sizeof(DATA_SET));

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
                zero_scalar(data->account_name,256);
                strcpy(data->account_name,argv[1]);
                argc -= 2;
                argv += 2;
                break;
	    }
	}

	if (input_data)
	{
	    ReadMarketData(data);
	    ReadAccountData(data);
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
	printf("\tCompare data (c)\n");
	printf("\tSimulate data (s)\n");
	printf("\tModify data (m)\n");
	printf("\tBuy shares (b)\n");
	printf("\tTrade (t)\n");
	printf("\tInitiation (i)\n");
	printf("\tRank data (r)\n");
	printf("\tPrint open trade (o)\n");
	printf("\tExperiment (e)\n");
	printf("\tDo nothing (n)\n");
	CreateJVM();

	while (YES)
	{
	    printf("Enter option: ");
	    scanf("%s",string);
	    switch(string[0])
	    {
	    case 'a':
	    	AddData(data);
		WriteMarketData(data);
    		PromptForDataMap(data);
    		XgraphData(data);
    		DataTrend(data);
		break;
	    case 'm':
	    	ModifyData(data);
		WriteMarketData(data);
		WriteAccountData(data);
    		XgraphData(data);
    		DataTrend(data);
		break;
	    case 'b':
	    	InvestShares(data);
		WriteAccountData(data);
		break;
	    case 't':
	    	if (!TradeShares(data))
		    closing_out(data);
		WriteAccountData(data);
		break;
	    case 'c':
	    	DataCompare(data);
		closing_out(data);
	    case 'i':
		InitTrade(data);
		closing_out(data);
	    case 'r':
		RankData(data);
		closing_out(data);
	    case 'o':
	    	TradeInfo(data);
		closing_out(data);
	    case 's':
	    	InvestSimulation(data);
		closing_out(data);
	    case 'e':
	    	if (!ExperimentTrade(data))
		    closing_out(data);
		WriteAccountData(data);
	    case 'n':
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
	closing_out(data);
}	/* end main */

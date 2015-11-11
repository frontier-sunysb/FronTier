#include "swap.h"

int main(int argc, char **argv)
{
	MARKET_DATA *data;
	boolean input_market_data = NO;
	boolean get_web_data = NO;
	char string[100];

	FT_Init(argc,argv,NULL);
	FT_ScalarMemoryAlloc((POINTER*)&data,sizeof(MARKET_DATA));

	argc--;
	argv++;

	while (argc >= 1)
	{
	    if (argv[0][0] != '-')
	    {
	    	printf("Usage: swap [-i data-file] [-w]\n");
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
		input_market_data = YES;
                zero_scalar(data->data_name,256);
                strcpy(data->data_name,argv[1]);
                argc -= 2;
                argv += 2;
                break;
	    }
	}

	if (input_market_data)
	{
	    ReadMarketData(data);
	}
	else if (get_web_data)
        {
            GetStockFile(data);
            SaveStockFile(data);
	    exit(0);
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
	printf("\tInitiation (i)\n");
	printf("\tModify data (m)\n");
	printf("\tRank data (r)\n");
	printf("\tSimulate data (s)\n");
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
		break;
	    case 'c':
	    	DataCompare(data);
		break;
	    case 'i':
		InitTrade(data);
		break;
	    case 'm':
	    	ModifyMarketData(data);
		WriteMarketData(data);
		break;
	    case 'r':
		RankData(data,NULL);
		break;
	    case 'n':
		break;
	    case 's':
		InvestSimulation(data);
		break;
	    case 'x':
    		XgraphData(data,NULL);
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

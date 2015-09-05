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
		input_data = YES;
                zero_scalar(data->data_name,256);
                strcpy(data->data_name,argv[1]);
                argc -= 2;
                argv += 2;
                break;
	    }
	}

	if (input_data)
	    ReadData(data);
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

	printf("Available operations are\n");
	printf("\tAdd data (a)\n");
	printf("\tCompare data (c)\n");
	printf("\tSimulate data (s)\n");
	printf("\tModify data (m)\n");
	printf("\tBuy shares (b)\n");
	printf("\tTrade (t)\n");
	printf("\tInformation (i)\n");
	printf("\tDo nothing (n)\n");

	while (YES)
	{
	    printf("Enter option: ");
	    scanf("%s",string);
	    switch(string[0])
	    {
	    case 'a':
	    	AddData(data);
		break;
	    case 'c':
	    	DataCompare(data);
		break;
	    case 's':
	    	InvestSimulation(data);
		break;
	    case 'm':
	    	ModifyData(data);
		break;
	    case 'b':
	    	InvestShares(data);
		break;
	    case 't':
	    	TradeShares(data);
		break;
	    case 'i':
	    	PromptForDataMap(data);
	    	XgraphData(data);
	    	DataInfo(data);
	    	DataTrend(data);
		break;
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
	WriteData(data);
	exit(0);
}	/* end main */

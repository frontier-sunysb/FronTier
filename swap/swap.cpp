#include "swap.h"

int main(int argc, char **argv)
{
	DATA_SET *data;
	boolean modify_data = NO;
	boolean compare_data = NO;
	boolean simulate_data = NO;
	boolean print_data = NO;
	boolean trade_data = NO;
	boolean get_web_data = NO;
	boolean input_data = NO;

	FT_Init(argc,argv,NULL);
	FT_ScalarMemoryAlloc((POINTER*)&data,sizeof(DATA_SET));

	argc--;
	argv++;

	while (argc >= 1)
	{
	    if (argv[0][0] != '-')
	    {
	    	printf("Usage: swap -i data-file [-m] [-c] [-s] [-t]\n");
		printf("For: modify, compare, simulate, trade data\n");
            	exit(1);
	    }
	    switch (argv[0][1]) {
	    case 'w':
                get_web_data = YES;
		argc -= 1;
                argv += 1;
                break;
	    case 'p':
		print_data = YES;
		argc -= 1;
                argv += 1;
		break;
	    case 'm':
		modify_data = YES;
		argc -= 1;
                argv += 1;
		break;
	    case 'c':
		compare_data = YES;
		argc -= 1;
                argv += 1;
		break;
	    case 's':
		simulate_data = YES;
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
	    case 't':
		trade_data = YES;
		argc -= 1;
                argv += 1;
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

	if (modify_data)
	{
	    ModifyData(data);
	}
	if (print_data)
	{
	    AddData(data);
	    PromptForDataMap(data);
	    XgraphData(data);
	    DataInfo(data);
	    DataTrend(data);
	}
	if (compare_data)
	    DataCompare(data);
	if (trade_data)
	{
	    InvestShares(data);
	    TradeShares(data);
	}

	if (simulate_data)
	    InvestSimulation(data);

	WriteData(data);
	exit(0);
}

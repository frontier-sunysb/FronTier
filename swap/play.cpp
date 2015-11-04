#include "swap.h"

int main(int argc, char **argv)
{
	MARKET_DATA *data;
	char string[256];
	int i,n;

	FT_Init(argc,argv,NULL);
	FT_ScalarMemoryAlloc((POINTER*)&data,sizeof(MARKET_DATA));

	argc--;
	argv++;

	while (argc >= 1)
	{
	    if (argv[0][0] != '-')
	    {
	    	printf("Usage: status [-i data-file]\n");
		printf("For: existing input\n");
            	exit(1);
	    }
	    switch (argv[0][1]) {
	    case 'i':
                zero_scalar(data->data_name,256);
                strcpy(data->data_name,argv[1]);
                argc -= 2;
                argv += 2;
                break;
	    }
	}

	ReadMarketData(data);
	InvestSimulation(data);
}	/* end main */

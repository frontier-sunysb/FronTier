#include "swap.h"

int main(int argc, char **argv)
{
	MARKET_DATA *data;
	PORTFOLIO **accounts;
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
	    	printf("Usage: status [-i data-file] [-a n accounts-file]\n");
		printf("For: existing input or getting data from web\n");
            	exit(1);
	    }
	    switch (argv[0][1]) {
	    case 'i':
                zero_scalar(data->data_name,256);
                strcpy(data->data_name,argv[1]);
                argc -= 2;
                argv += 2;
                break;
	    case 'a':
                zero_scalar(string,256);
                strcpy(string,argv[1]);
		n = atoi(string);
		FT_VectorMemoryAlloc((POINTER*)&accounts,n,sizeof(PORTFOLIO*));
                argc -= 2;
                argv += 2;
		for (i = 0; i < n; ++i)
		{
		    FT_ScalarMemoryAlloc((POINTER*)&accounts[i],
					sizeof(PORTFOLIO));
                    zero_scalar(accounts[i]->account_name,256);
                    strcpy(accounts[i]->account_name,argv[i]);
		}
                argc -= n;
                argv += n;
                break;
	    }
	}

	ReadMarketData(data);
	for (i = 0; i < n; ++i)
	{
	    accounts[i]->data = data;
	    ReadAccountData(accounts[i]);
	}
	CreateJVM();

	ReportDataStates(accounts,n);
}	/* end main */

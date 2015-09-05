#include <FronTier.h>
#include <iostream>
#include <string>
#include <sstream>
#include <dirent.h>
#include <sys/types.h>
#include <fstream>
#include "boost/asio.hpp"

using namespace std;

#define         MAX_TRACE               212
#define         MAX_NUM_TRADE           200

const char Xcolor[][100] ={"red","orange","blue",
			   "green","violet","navy",
			   "cyan","yellow","fuchsia",
			   "gold","pink","aqua"};

enum _TRADE_STATUS {
	OPEN			= 0,
	CLOSED,
	PARTIAL_CLOSED
};
typedef enum _TRADE_STATUS TRADE_STATUS;

struct _ASSET
{
	char asset_name[100];
	char color[20];
	double *value;
	double *norm_value;
	double base_value;
	double A;		/* fit norm_value to: y = A*exp(lambda*t) */
	double lambda;
};

struct _DATA_SET
{
	char data_name[258];
	int num_days;
	int num_assets;
	int num_backtrace;
	int istar;
	struct _ASSET *assets;
	char **date;
	boolean new_data;
	boolean *data_map;
	int *shares;
	int num_trades;
	struct _TRADE *trades;
};

#define		MAX_NUM_STAGES		20
struct _TRADE
{
	int num_stages;
	TRADE_STATUS status;
	int index_buy[MAX_NUM_STAGES];
	int index_sell[MAX_NUM_STAGES];
	int num_shares_buy[MAX_NUM_STAGES];
	int num_shares_sell[MAX_NUM_STAGES];
	double price_buy[MAX_NUM_STAGES];
	double price_sell[MAX_NUM_STAGES];
};

typedef struct _ASSET ASSET;
typedef struct _DATA_SET DATA_SET;
typedef struct _TRADE TRADE;

/* sim.cpp */
extern void InvestSimulation(DATA_SET*);

/* sub.cpp */
extern double GetLeastSquare(DATA_SET*,int,int,double*,double*);
extern double LeastSquareQuadr(double*,double*,int,double*,double*,double*);
extern double LeastSquareLinear(double*,double*,int,double*,double*);
extern void ContinueBaseAdjust(DATA_SET*);
extern void AdjustBase(DATA_SET*,int,int);
extern void ModifyData(DATA_SET*);
extern void ReadData(DATA_SET*);
extern void WriteData(DATA_SET*);
extern void AddData(DATA_SET*);
extern void XgraphData(DATA_SET*);
extern void DataInfo(DATA_SET*);
extern void DataCompare(DATA_SET*);
extern void DataTrend(DATA_SET*);
extern void PromptForDataMap(DATA_SET*);
extern DATA_SET *CopyData(DATA_SET*);

/* web.cpp */
extern string GetWebData(const std::string& server,const std::string& file);
extern void GetStockFile(DATA_SET*);
extern void SaveStockFile(DATA_SET*);

/* trade.cpp */
extern void InvestShares(DATA_SET*);
extern void TradeShares(DATA_SET*);
extern void PrintOpenTrade(DATA_SET*);
extern void PrintClosedTrade(FILE*,TRADE,DATA_SET*);
extern void PrintClosedTradeLoop(FILE*,TRADE,DATA_SET*);
extern void SaveDeleteClosedTradeLoop(DATA_SET*);
extern void SortTradeOrder(DATA_SET*);

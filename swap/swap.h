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
	double *price;
	double *norm_price;
	//double base;
};

struct _MARKET_DATA
{
	char data_name[258];
	int num_segs;
	int num_assets;
	int num_backtrace;
	struct _ASSET *assets;
	double *base;
	boolean new_data;
};

struct _PORTFOLIO
{
	struct _MARKET_DATA *data;
	char account_name[256];
	boolean *data_map;
	int *shares;
	int num_trades;
	struct _TRADE *trades;
	double polar_ratio;
	int eindex;		// Index of share equivalence
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

struct _STATE_INFO {
	char stock_max[20];	// name of stick with maximum normal price
	char stock_min[20];	// name of stick with minimum normal price
	double dnp_max;		// np diff between max and next to max
	double dnp_min;		// np diff between min and next to min
        double svalue;         	// max surplus value
        double dvalue;         	// max deficit value
        double dnp;            	// np diff between surplue and deficit
        char sname[20];         // name of stock with max surplus
        char dname[20];         // name of stock with max deficit
	int S[5];		// normalized share of five states
	int C;			// current normalized shares
};

typedef struct _ASSET ASSET;
typedef struct _DATA_SET DATA_SET;
typedef struct _TRADE TRADE;
typedef struct _PORTFOLIO PORTFOLIO;
typedef struct _MARKET_DATA MARKET_DATA;
typedef struct _STATE_INFO STATE_INFO;


/* sim.cpp */
extern void InvestSimulation(MARKET_DATA*);
extern void CompareToFiveStates(PORTFOLIO*);
extern void PrintCurrentLinearProfile(PORTFOLIO*);
extern void PrintDataStates(PORTFOLIO*);
extern void ReportDataStates(PORTFOLIO**,int);

/* sub.cpp */
extern double GetLeastSquare(MARKET_DATA*,int,int,double*,double*);
extern void XgraphData(MARKET_DATA*,boolean*);
extern void DataTrend(MARKET_DATA*,boolean*);
extern void ReadMarketData(MARKET_DATA*);
extern void WriteMarketData(MARKET_DATA*);
extern void AddData(MARKET_DATA*,boolean);
extern void ModifyMarketData(MARKET_DATA*);
extern void DataCompare(MARKET_DATA*);
extern void RankData(MARKET_DATA*,boolean*);
extern void InitTrade(MARKET_DATA*);
extern void PrintAssetList(MARKET_DATA*);
extern void PrintEquityIndex(MARKET_DATA*);
extern void TimelyRecordMarketData(MARKET_DATA*);
extern void ComputeAllNormPrice(MARKET_DATA*,int);
extern void FreeMarketData(MARKET_DATA*);
extern MARKET_DATA *CopyMarketData(MARKET_DATA*);
extern void LeastSquareQuadr(double*,double*,int,double*,double*,double*);
extern void LeastSquareLinear(double*,double*,int,double*,double*);
extern void LeastSquareConstant(double*,double*,int,double*);
extern void TradeInfo(PORTFOLIO*);
extern void PromptForDataMap(PORTFOLIO*);
extern void ReadAccountData(PORTFOLIO*);
extern void WriteAccountData(PORTFOLIO*);
extern void PrintAccountValue(PORTFOLIO*);

/* web.cpp */
extern string GetWebData(const std::string& server,const std::string& file);
extern void GetStockFile(MARKET_DATA*);
extern void SaveStockFile(MARKET_DATA*);
extern void GetCurrentPrice(char**,double*,int);
extern void CreateJVM();
extern void DestroyJVM();

/* trade.cpp */
extern void InvestShares(PORTFOLIO*);
extern void PrintOpenTrade(PORTFOLIO*);
extern void PrintOpenTradeForIndex(PORTFOLIO*,int);
extern void SaveDeleteClosedTradeLoop(PORTFOLIO*);
extern void SortTradeOrder(PORTFOLIO*);
extern void WrapPartialTradeLoop(PORTFOLIO*);
extern void ReOrganizeTrade(PORTFOLIO*);
extern boolean TradeShares(PORTFOLIO*);
extern boolean ExperimentTrade(PORTFOLIO*);
extern double PrintClosedTradeLoop(FILE*,TRADE,PORTFOLIO*);

#define		closing_out(data)				\
		{						\
			PrintDataStates((data));		\
			PrintAccountValue((data));		\
			DestroyJVM();				\
			exit(0);				\
		}						

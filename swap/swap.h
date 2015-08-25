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

const char Xcolor[][100] ={"red","orange","blue",
			   "green","violet","navy",
			   "cyan","yellow","fuchsia",
			   "light-grey","dark-gray","aqua",
			   "gold","pink"};

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
	struct _ASSET *assets;
	char **date;
	boolean new_data;
	boolean *data_map;
};

typedef struct _ASSET ASSET;
typedef struct _DATA_SET DATA_SET;

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

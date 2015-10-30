#include "swap.h"
#include <jni.h>  

const int NUM_STOCKS = 100;
int stock_days = 0;
vector < vector <float> > stock_data;
vector < string > filename;
char **stock_file_name;

static boolean IsCash(char*);

extern string GetWebData(
	const std::string& server, 
	const std::string& file)
{
	try
	{
	    boost::asio::ip::tcp::iostream s(server, "http");
		
	    s.expires_from_now(boost::posix_time::seconds(60));

	    if (!s){ throw "Unable to connect: " + s.error().message(); }

	    // ask for the file
	    s << "GET " << file << " HTTP/1.0\r\n";
	    s << "Host: " << server << "\r\n";
	    s << "Accept: */*\r\n";
	    s << "Connection: close\r\n\r\n";

	    // Check that response is OK.
	    string http_version;
	    s >> http_version;
	    unsigned int status_code;
	    s >> status_code;
	    string status_message;
	    getline(s, status_message);
	    if (!s && http_version.substr(0, 5) != "HTTP/")
	    { 
		throw "Invalid response\n"; 
	    }
	    if (status_code != 200){throw "Response returned with status code "
				 + status_code; }

	    // Process headers, which are terminated by a blank line.
	    string header;
	    while (getline(s, header) && header != "\r")
	    {
		;
	    }

	    // Write the remaining data to output.
	    stringstream ss;
	    ss << s.rdbuf();
	    return ss.str();
	}
	catch(exception& e)
	{
	    return e.what();
	}
}	/* end GetWebData */

extern void GetStockFile(
	DATA_SET *data)
{
	string stocks[NUM_STOCKS];
	string server_addr = "ichart.finance.yahoo.com";
	string file_addr = "";
	string file_name = "";
	int    s_month;
	string s_day,s_year;
	int    e_month;
        string e_day,e_year;
	int i;
	int  n_stock;
	char str_smonth[20];
	char str_emonth[20];

	string file_folder = "stock_data";
	system("exec rm -rf ./stock_data");
	system("exec mkdir ./stock_data");
	printf("Enter number of stocks: ");
        scanf("%d", &n_stock);
	printf("Enter starting year: ");
	cin >> s_year;
	printf("Enter starting month: ");
        scanf("%d", &s_month);
        printf("Enter starting day: ");
        cin >> s_day;
	printf("Enter end year: ");
        cin >> e_year;
        printf("Enter end month: ");
        scanf("%d", &e_month);
        printf("Enter end day: ");
        cin >> e_day;
	s_month -=1;
	e_month -=1;
	sprintf(str_smonth,"%d",s_month);
	sprintf(str_emonth,"%d",e_month);
	FT_MatrixMemoryAlloc((POINTER*)&stock_file_name,n_stock,100,
				sizeof(char));

	i = 0;
	while(i < n_stock)
	{
	    printf("Enter name of stock %d: ",i);
	    cin >> stocks[i];
	    filename.push_back(stocks[i]);
	    file_addr = "/table.csv?s=" + stocks[i] + "&a=" + 
			str_smonth + "&b=" + s_day + "&c=" + s_year + 
			"&d=" + str_emonth + "&e=" + e_day + "&f=" + 
			e_year + "&g=d&y=0&z=30000";
	    file_name = file_folder + "/" + stocks[i];
	    string result = GetWebData(server_addr.c_str(), 
			file_addr.c_str());
	    ofstream of(file_name.c_str(), std::ios::binary);
	    sprintf(stock_file_name[i],"%s",file_name.c_str());
	    of << result;
	    of.flush();
	    of.close();
	    cout << "Get the data file" << "\t" <<  stocks[i] << endl;
	    i++;
	}
}	/* end GetStockFile */

extern void SaveStockFile(
	DATA_SET *data)
{	
	int i = 0,j;
 	char p; //place holder, stores ","
	float date,open,high,low,close,volume, adjclose;// stocks attribute

	string line;// each line of the file
	vector <float> line_data; // only stores the meaningful daily data

	int file_num = 0;
	int back_trace_days = 0;

	char *dataset_name = new char[20];
	char *color_name =  new char[20];

	
	char *mean_name =  new char[20]; //mean name, basevalue, and 
				  	 //initialization data
	float mean_basevalue;
	float mean_iti_value = 1;
	ifstream f;	// list file in the directory

	FILE *outfile;

	file_num = filename.size(); // get the number of stocks

	for (i = 0; i < file_num; i++)
	{
	    ifstream f;
	    f.open(stock_file_name[i]);
	    while (getline(f, line)) 
	    {
		// Tell whether this line is meaningful, since it 
		// has some strings data, For example, 
		// "DIVIDEND, 20150211,0.688000"
		if (line[0] <= '9' && line[0] >= '0') 
		{
		    stringstream ss;
            	    ss << line;
		    line_data.clear();
		    ss >> date >> p >> open >> p >> high >> p >> 
			low >> p >> close >> p >> volume >> p >> adjclose;
		    line_data.push_back(date);
		    line_data.push_back(open);
		    line_data.push_back(high);
		    line_data.push_back(low);
		    line_data.push_back(close);
		    line_data.push_back(volume);
		    line_data.push_back(adjclose);
		    stock_data.push_back(line_data);
		}
	    }
	    f.close();
	}
	stock_days = stock_data.size()/file_num;
			// get how many days of stock data

	printf("There are total of %d stocks with %d days of data: \n",
			file_num, stock_days);
	printf("Enter your dataset name: ");
	scanf("%s", dataset_name);
	memcpy(data->data_name,dataset_name,sizeof(dataset_name)+1);
	outfile = fopen(dataset_name,"w");
	fprintf(outfile,"%s\n",dataset_name);	
	fprintf(outfile,"%d\n",file_num);
	fprintf(outfile,"%d\n",stock_days);
	printf("Enter number of days for back trace: ");
	scanf("%d",&back_trace_days);
			// input your back_trace_days, which is 
			// less than the "stock_days"
	fprintf(outfile,"%d\n",back_trace_days);
	for(i = 0; i < file_num; i++)
	{
	    fprintf(outfile,"%s\n",filename[i].c_str());
	    fprintf(outfile,"%s\n",Xcolor[i%12]);
	    fprintf(outfile,"%f\n",stock_data[(i + 1) * stock_days-1][4]);
			// try to get each stocks first data, ie, 
			// base value (Inverse order)
	}
	
	fprintf(outfile,"Mean\n");
	fprintf(outfile,"aqua\n");
	fprintf(outfile,"1.0\n");
	
	for(i = stock_days; i > 0; i--)
	{
	    for(j = 0; j < file_num; j++)
	    {
		fprintf(outfile, "%f ", stock_data[stock_days*(j+1) -
			(stock_days - i + 1)][4]);
			//Inverse order to get all the data
	    }
	    fprintf(outfile, "%f\n",mean_iti_value);
			// the initialization of the last data in 
			// each line, which is the "mean base value"
	}
	printf("All the downloaded data will go to file: %s\n",
			data->data_name);
	fclose(outfile);
	ReadMarketData(data);
}	/* end SaveStockFile */

	JNIEnv* env;
	JavaVM *jvm;

extern void GetCurrentPrice(
	char **stock_names,
	double *prices,
	int num_stocks)
{  
	int i;
	JavaVMOption options[1];
	JavaVMInitArgs vm_args;  
	char String[256];
   
	long status;  
	jclass cls;  
	jmethodID mid;  
   
	options->optionString = String;
	sprintf(options[0].optionString,"-Djava.class.path=.:jsoup-1.8.3.jar");
	memset(&vm_args, 0, sizeof(vm_args));
	vm_args.version = JNI_VERSION_1_6;
	vm_args.nOptions = 1;
	vm_args.options = options;
	if (status != JNI_ERR)
	{
	    cls = env->FindClass("DG"); //find class (java .class file) 
	    if (cls != 0)
	    {
		mid = env->GetStaticMethodID(cls,"getRealTimePrice",
				"(Ljava/lang/String;)D"); 
					// the last para is method Sign Name 
		if (mid != 0)
		{	
		    for (i = 0; i < num_stocks; i++)
		    {
			
			int count = 0;
			prices[i] = IsCash(stock_names[i]) ? 1.0 : 0.0;
		    	jstring str_arg = env->NewStringUTF(stock_names[i]);
			while (prices[i] == 0.0)
			{
			    prices[i] = env->CallStaticDoubleMethod(cls,mid,
						str_arg);
			    if (count++ > 5)
			    {
				printf("Cannot get data 5 times\n");
				clean_up(ERROR);
			    }
			}

		    }
		}
	    }
	}
	else
	{
	    printf("JVM Created failed!\n");
	}
}	/* GetCurrentPrice */

extern void DestroyJVM()
{
	jvm->DestroyJavaVM();
}	/* end DestroyJVM */

extern void CreateJVM()
{
        JavaVMOption options[1];
        JavaVMInitArgs vm_args;
	char java_string[256];
	sprintf(java_string,"-Djava.class.path=.:jsoup-1.8.3.jar");
	options[0].optionString = java_string;
        memset(&vm_args, 0, sizeof(vm_args));
        vm_args.version = JNI_VERSION_1_6;
        vm_args.nOptions = 1;
        vm_args.options = options;
        // create java virtual machine
        JNI_CreateJavaVM(&jvm,(void**)&env,&vm_args);
}	/* end CreateJVM */

static boolean IsCash(char *name)
{
        if (strcmp(name,"cash") == 0 ||
            strcmp(name,"CASH") == 0 ||
            strcmp(name,"Cash") == 0)
            return YES;
        else
            return NO;
}       /* end IsCash */

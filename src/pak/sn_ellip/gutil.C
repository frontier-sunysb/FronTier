#include "gutil.h"
#undef vector

#include <vector>
#include <string>
#include <iostream>


/* the one node version of Debug_Print */
void Debug_Print(const char *fmt, ...)
{
	va_list ap;
	va_start(ap, fmt);
	
	FILE *hdebug;
	char debug_name[100]; // "debug?"

	sprintf(debug_name, "%s%d%s", "debug", 0, ".txt");
	hdebug=fopen(debug_name,"a+");

	if(hdebug==NULL)
	{
		printf("\n PID %d: can't open %s for writting", 0, debug_name);
		exit(1);
	}

	(void) vfprintf(hdebug, fmt, ap);
	vprintf(fmt, ap);

	fclose(hdebug);
	
	va_end(ap);
}

/* the one node version of Debug_Print, it gives the output to filename */
void Debug_Printf(char *filename, const char *fmt, ...)
{
	va_list ap;
	va_start(ap, fmt);
	
	FILE *hdebug;
	char debug_name[100]; // "debug?"

	sprintf(debug_name, "%s", filename);
	hdebug=fopen(debug_name,"a+");

	if(hdebug==NULL)
	{
		printf("\n PID %d: can't open %s for writting", 0, debug_name);
		exit(1);
	}

	(void) vfprintf(hdebug, fmt, ap);
	vprintf(fmt, ap);

	fclose(hdebug);
	
	va_end(ap);
}

/* the parallel version of Debug_Print */
void Debug_Print(int my_rank, const char *fmt, ...)
{
	va_list ap;
	va_start(ap, fmt);
	
	FILE *hdebug;
	char debug_name[100]; // "debug?"

	sprintf(debug_name, "%s%d%s", "debug", my_rank, ".txt");
	hdebug=fopen(debug_name,"a+");

	if(hdebug==NULL)
	{
		printf("\n PID %d: can't open %s for writting", my_rank, debug_name);
		exit(1);
	}

	(void) vfprintf(hdebug, fmt, ap);
	vprintf(fmt, ap);

	fclose(hdebug);
	
	va_end(ap);
}

// given n integer P[], sort the numbers from small to large 
void Sort(int n, int *P)
{
	int i, j, pr, tmp;
	
	// sort
	for(i=0; i<n; i++)
	{
		pr = i;
		for(j=i; j<n; j++)
			if(P[j]<P[pr])
				pr = j;
		tmp   = P[i];
		P[i]  = P[pr];
		P[pr] = tmp;
	}
}
void Sort(int n, unsigned long *P)
{
	int i, j, pr;
	unsigned long tmp;
	
	// sort
	for(i=0; i<n; i++)
	{
		pr = i;
		for(j=i; j<n; j++)
			if(P[j]<P[pr])
				pr = j;
		tmp   = P[i];
		P[i]  = P[pr];
		P[pr] = tmp;
	}

}

void Sort(int n, double *P)
{
	int i, j, pr;
	double tmp;
	
	// sort
	for(i=0; i<n; i++)
	{
		pr = i;
		for(j=i; j<n; j++)
			if(P[j]<P[pr])
				pr = j;
		tmp   = P[i];
		P[i]  = P[pr];
		P[pr] = tmp;
	}
}

// the same as Sort() except that the input has another array index[] with the same dimention as
// P[]. index[] is associated with P[] and Index[] is sorted with P[].
void Sort(int n, double *P, int *Index)
{
	int i, j, pr, ptmp;;
	double tmp;
	
	// sort
	for(i=0; i<n; i++)
	{
		pr = i;
		for(j=i; j<n; j++)
			if(P[j]<P[pr])
				pr = j;
		tmp   = P[i];
		P[i]  = P[pr];
		P[pr] = tmp;
		ptmp  	  = Index[i];
		Index[i]  = Index[pr];
		Index[pr] = ptmp;		
	}
}

void Sort(int n, int *P, int *Index)
{
	int i, j, pr, ptmp;;
	int tmp;
	
	// sort
	for(i=0; i<n; i++)
	{
		pr = i;
		for(j=i; j<n; j++)
			if(P[j]<P[pr])
				pr = j;
		tmp   = P[i];
		P[i]  = P[pr];
		P[pr] = tmp;
		ptmp  	  = Index[i];
		Index[i]  = Index[pr];
		Index[pr] = ptmp;		
	}
}

void Sort(double *P, std::vector<int> &Index)
{
	int i, j, pr, ptmp, n = Index.size();
	double tmp;
	
	// sort
	for(i=0; i<n; i++)
	{
		pr = i;
		for(j=i; j<n; j++)
			if(P[j]<P[pr])
				pr = j;
		tmp   = P[i];
		P[i]  = P[pr];
		P[pr] = tmp;
		ptmp  	  = Index[i];
		Index[i]  = Index[pr];
		Index[pr] = ptmp;		
	}
}

/**************************************************************************
 *	SortDelete()
 * given n integer P[], sort the numbers from small to large and then delete
 * redundant element, return the number of elements left.
 **************************************************************************/
int SortDelete(int n, int *P)
{
	int i, j, pr;
	
	// sort
	Sort(n, P);
	// delete redundant points
	pr = 0;
	for(i=0; i<n; i++)
	{
		if(P[i]==P[pr])
			; // do nothing
		else
		{
			pr++;
			P[pr] = P[i];
		}
	}
	return pr+1;
}
int SortDelete(int n, double *P, int *Index)
{
	int i, j, pr; 
	
	// sort
	Sort(n, P, Index);
	// delete redundant points
	pr = 0;
	for(i=0; i<n; i++)
	{
		if(P[i]==P[pr])
			; // do nothing
		else
		{
			pr++;
			P[pr] 	  = P[i];
			Index[pr] = Index[i];
		}
	}
	return pr+1;
}

void Swap(int &a, int &b)
{
	int tmp;
	tmp = a;
	a = b;
	b = tmp;
}


void Swap(double &a, double &b)
{
	double tmp;
	tmp = a;
	a = b;
	b = tmp;
}
/*************************************************************************
 *	SaveAsGRD()	
 * Output the matrix A[][] to a txt file with file extion of .GRD which 
 * could be used by Surfer32.exe.
 *************************************************************************/
void SaveAsGRD(char *filename, int m, int n, double **A)
{
	FILE *hfile = fopen(filename,"w");
	if(hfile==NULL)
	{
		Debug_Print("\n can't open %s in SaveAsGRD()", filename);
		exit(0);
	}
	int i, j;
	double min, max;
	min = A[0][0];	max = A[0][0];
	// find the min/max value
	for(i=0; i<m; i++)
		for(j=0; j<n; j++)
		{
			if(A[i][j]<min)
				min = A[i][j];
			if(A[i][j]>max)
				max = A[i][j];
		}
	fprintf(hfile, "DSAA\n");
	fprintf(hfile, "%d %d\n", n, m);	// should be n m
	fprintf(hfile, "%d %d\n", 1, n);
	fprintf(hfile, "%d %d\n", 1, m);
	fprintf(hfile, "%.4f %.4f\n", min, max);

	for(i=0; i<m; i++)
		for(j=0; j<n; j++)
			fprintf(hfile," %.4f\n", A[i][j]);
	
	fclose(hfile);
}

/*************************************************************************
 *	SaveAsTecplot()	
 * Output the matrix A[][] to a txt file with file extion of .plt which 
 * could be used by Tecplot.exe.
 *************************************************************************/
void SaveAsTecplot(char *filename, int m, int n, double **A)
{
	FILE *hfile = fopen(filename,"w");
	if(hfile==NULL)
	{
		Debug_Print("\n can't open %s in SaveAsTecplot()", filename);
		exit(0);
	}
	int i, j;
	
	// header
	fprintf(hfile, "TITLE = %s	\n", filename);
	fprintf(hfile, "VARIABLES = \"X\" \"Y\" \"Z\" \n");
	fprintf(hfile, "ZONE I=%d J=%d F=POINT	\n", n, m);		// n, m
	// now the data

	for(i=0; i<m; i++)
		for(j=0; j<n; j++)
			fprintf(hfile,"%6d %6d %f\n", i, j, A[i][j]);
	
	fclose(hfile);
}

void Tecplot_Point(char *filename, int n, double **P)
{
	int i;
	FILE *hfile = fopen(filename, "w");
	if(hfile==NULL)
	{
		printf("\n can't open %s in Tecplot_Point().", filename);
		exit(0);
	}
		
	fprintf(hfile, "TITLE = %s	\n", filename);	
	fprintf(hfile, "VARIABLES = \"X\" \"Y\" \n");	
	for(i=0; i<n; i++)
	{
		fprintf(hfile, "ZONE I=%d J=%d F=POINT	\n", 1, 1);	
		fprintf(hfile, " %.8f %.8f \n", P[i][0], P[i][1]);
	}
	fclose(hfile);

}

FILE *FOpen(char *filename)
{
	FILE *hfile = fopen(filename, "w");
	if(hfile==NULL)
	{
		printf("\n FOpen: can't open %s for writing.", filename);
		exit(0);
	}
	return hfile;
}


double TimeDifference(struct timeb &t1, struct timeb &t2)
{
	int tm;
	unsigned short millitm;
	if(t2.millitm<t1.millitm)
	{
		millitm = t2.millitm + 1000 - t1.millitm;
		tm	= t2.time - t1.time - 1;
	}
	else
	{
		millitm = t2.millitm - t1.millitm;
		tm	= t2.time - t1.time;
	}
	
	//double dm = tm + millitm/1000.0;
	return tm + millitm/1000.0;;
}
double Min(double p, double q)
{
	if(p<=q)
		return p;
	else
		return q;
}

double Min(double p, double q, double r)
{
	double a = p;
	if(a<q)
		a = q;
	if(a<r)
		a = r;
	return a;
}

double Min(double p, double q, double r, double s)
{
	double a = p;
	if(a<q)
		a = q;
	if(a<r)
		a = r;
	if(a<s)
		a = s;
	return a;
}
int Max(int p, int q)
{
	if(p<q)
		return q;
	return p;
}
double Max(double p, double q)
{
	if(p<q)
		return q;
	return p;
}
double Max(double p, double q, double r)
{
	double a = p;
	if(a>q)
		a = q;
	if(a>r)
		a = r;
	return a;
}
double Max(double p, double q, double r, double s)
{
	double a = p;
	if(a>q)
		a = q;
	if(a>r)
		a = r;
	if(a>s)
		a = s;
	return a;
}
// find the index for the largest value
int MaxIndex(int n, double *value)
{
	int index = 0;
	for(int i=0; i<n; i++)
		if(value[index]<value[i])
			index = i;
	return index;	
}
int MaxIndex(int n, int *value)
{
	int index = 0;
	for(int i=0; i<n; i++)
		if(value[index]<value[i])
			index = i;
	return index;	
}
int MinIndex(int n, int *value)
{
	int index = 0;
	for(int i=0; i<n; i++)
		if(value[index]>value[i])
			index = i;
	return index;	
}

// convert a integer (0-10000) to a string with size 5 and the first few characters filled with 0;
void IntegerToString(char *str, int index)
{
	if(index<0)
		sprintf(str, "%d", index);
	else if(index<   10)
		sprintf(str, "0000%d", index);
	else if(index<  100)
		sprintf(str, "000%d", index);
	else if(index< 1000)
		sprintf(str, "00%d", index);
	else if(index<10000)
		sprintf(str, "0%d", index);
	else
		sprintf(str, "%d", index);		
}	


/******************************************************************************************
 * if str contains "-p a b", this function reads from str and assign row = a and col = b
 * if str contains "-p a " only, this function lets row = a;
 ******************************************************************************************/
void GetProcessors(int argc, char **argv, int &row, int &col)
{
	using namespace std;
	string str;
	row = 0;
	col = 0;
	
	cout << endl;
	for(int i=0; i<argc; i++)
		str = str+ " " + argv[i];
	cout << "str = " << str	<< endl;
	string::size_type pos = str.find("-p");
	if(pos==string::npos)		// did not find the substring;
		return;	
	string str2 = str.substr(pos+3, str.size()-(pos+3));
	cout << "str2 = " << str2 << endl;
	cout << "not finished yet." << endl;
}


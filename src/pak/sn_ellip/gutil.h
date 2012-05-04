/********************************************************************************************
 *      gutil.c/gutil.h
 * This file contains some useful functions 
 * 
 * Debug_Print(...):          the same as printf, but the output is given in debug0.txt.
 * Debug_Print(my_rank,...):  the same as Debug_Print(), but the output is given in debug?.txt,
 *                            where ? is the 'my_rank' for the current node (MPI).
 *********************************************************************************************/
#ifndef  _GUTIL_H
#define  _GUTIL_H


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdarg.h>
#include <time.h>
#include <sys/timeb.h>

	/* gutil.c */
void Debug_Print(const char *fmt, ...);
void Debug_Printf(char *filename, const char *fmt, ...);
void Debug_Print(int my_rank, const char *fmt, ...);
void Sort(int n, int *P);
void Sort(int n, unsigned long *P);
void Sort(int n, double *P);
void Sort(int n, double *P, int *Index);
void Sort(int n, int *P, int *Index);
/*void Sort(double *P, std::uni_array<int> &Index); */

int SortDelete(int n, int *P);
int SortDelete(int n, double *P, int *index);

void Swap(int &a, int &b);
void Swap(double &a, double &b);

void SaveAsGRD(char *filename, int m, int n, double **A);
void SaveAsTecplot(char *filename, int m, int n, double **A);
void Tecplot_Point(char *filename, int n, double **P);

FILE *FOpen(char *filename);

double TimeDifference(struct timeb &t1, struct timeb &t2);	
double Min(double p, double q);
double Min(double p, double q, double r);
double Min(double p, double q, double r, double s);
int Max(int p, int q);
double Max(double p, double q);
double Max(double p, double q, double r);
double Max(double p, double q, double r, double s);
int MaxIndex(int n, double *value);
int MaxIndex(int n, int *value);
int MinIndex(int n, int *value);

void IntegerToString(char *str, int index);	/* convert a integer (0-10000) to a string with size 5 and the first few characters filled with 0; */


	/* elliptic_solution.c */
/*void elliptic_solution(INIT_DATA*,struct _CHART*,const IO_TYPE*); */

void GetProcessors(int argc, char**argv, int &row, int &col);


#endif  /* #ifndef _GUTIL_H */


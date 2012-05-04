/********************************************************************************************
 *		pgutil.h
 * utility for paralell computing.
 ********************************************************************************************/
#ifndef PGUTIL_H
#define PGUTIL_H

#include <mpi.h>
#include <string.h>
#include <stdio.h>

void ParaSumGather(double Data, int des, double *recvbuf);
void ParaSumGather(int nProcessor, int *Processor, double Data, int des, double *recvbuf);
boolean ParaMax(double data, MPI_Comm comm);
boolean ParaMin(double data, MPI_Comm comm);

boolean ParaPositive(int data, MPI_Comm comm);	/* check whether max(data[])>0 */
void ParaGlobalIndexing(int &startindex, int data, MPI_Comm comm);

#endif

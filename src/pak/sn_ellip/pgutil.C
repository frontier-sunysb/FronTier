/*****************************************************************************************
 *		pgutil.c
 * See pgutil.h
 *****************************************************************************************/
#include <pgutil.h>
#include <string.h>
/**************************************************************************
 *		ParaSumGather()
 *	Similar with MPI_Gather().
 *
 *	SumGather data to processor with id: des
 **************************************************************************/
void ParaSumGather(double Data, int des, double *recvbuf)
{

	MPI_Comm  comm = MPI_COMM_WORLD;
	int comm_size, my_rank;
	MPI_Comm_rank(comm, &my_rank);
        MPI_Comm_size(comm, &comm_size);
	
			
	int TAG=0;
	MPI_Status status;
	double data=Data;
	*recvbuf=0;
	if(my_rank==des)
	{
		int i;
		*recvbuf=data;
		for(i=0; i<comm_size; i++)
		{
			if(i==des)	continue;	// *recvbuf=data;
			MPI_Recv(&data, 1, MPI_DOUBLE, i, TAG, comm, &status);
			*recvbuf += data;
		}
	}
	else
	{
		MPI_Send(&data, 1, MPI_DOUBLE, des, TAG, comm);
	}

}
// same as above, but the processors are those stored in Processor[];
void ParaSumGather(int nProcessor, int *Processor, double Data, int des, double *recvbuf)
{
	MPI_Comm  comm = MPI_COMM_WORLD;
	int comm_size, my_rank;
	MPI_Comm_rank(comm, &my_rank);
//      MPI_Comm_size(comm, &comm_size);
	comm_size = nProcessor;
	
			
	int TAG=0;
	MPI_Status status;
	double data=Data;
	*recvbuf=0;
	if(my_rank==des)
	{
		int i;
		*recvbuf=data;
		for(i=0; i<comm_size; i++)
		{
			if(Processor[i]==des)	continue;	// *recvbuf=data;
			MPI_Recv(&data, 1, MPI_DOUBLE, Processor[i], TAG, comm, &status);
			*recvbuf += data;
		}
	}
	else
	{
		MPI_Send(&data, 1, MPI_DOUBLE, des, TAG, comm);
	}

}

bool ParaMax(double data, MPI_Comm comm)
{
	//MPI_Comm  comm = MPI_COMM_WORLD;
	int size, my_rank;
	MPI_Comm_rank(comm, &my_rank);
	MPI_Comm_size(comm, &size);
	double *Data = new double[size];
	// gather the data
	MPI_Gather(&data, 1, MPI_DOUBLE, Data, 1, MPI_DOUBLE, 0, comm);
	// find the Max;
	if(my_rank==0)
	{
		int index = 0;	
		for(int i=0; i<size; i++)
			if(Data[i]>Data[index])
				index = i;
		for(int i=0; i<size; i++)
			if(i==index)
				Data[i] = 1;
			else
				Data[i] = -1;			
	}
	MPI_Scatter(Data, 1, MPI_DOUBLE, &data, 1, MPI_DOUBLE, 0, comm);
	delete [] Data;
	if(data>0)
		return true;
	else 
		return false;
		
}

bool ParaMin(double data, MPI_Comm comm)
{
	//MPI_Comm  comm = MPI_COMM_WORLD;
	int size, my_rank;
	MPI_Comm_rank(comm, &my_rank);
	MPI_Comm_size(comm, &size);
	double *Data = new double[size];
	// gather the data
	MPI_Gather(&data, 1, MPI_DOUBLE, Data, 1, MPI_DOUBLE, 0, comm);
	// find the Max;
	if(my_rank==0)
	{
		int index = 0;	
		for(int i=0; i<size; i++)
			if(Data[i]<Data[index])
				index = i;
		for(int i=0; i<size; i++)
			if(i==index)
				Data[i] = 1;
			else
				Data[i] = -1;			
	}
	MPI_Scatter(Data, 1, MPI_DOUBLE, &data, 1, MPI_DOUBLE, 0, comm);
	delete [] Data;
	if(data>0)
		return true;
	else 
		return false;
		
}
// check whether max(data[])>0
// see also ParaMax()
bool ParaPositive(int data, MPI_Comm comm)
{
	//MPI_Comm  comm = MPI_COMM_WORLD;
	int size, my_rank;
	MPI_Comm_rank(comm, &my_rank);
	MPI_Comm_size(comm, &size);
	int *Data;
	if(my_rank==0)
		Data = new int[size];
	else
		Data = NULL;
	// gather the data
	MPI_Gather(&data, 1, MPI_INT, Data, 1, MPI_INT, 0, comm);
	// find the Max;
	if(my_rank==0)
	{
		int index = 0;	
		for(int i=0; i<size; i++)
			if(Data[i]>Data[index])
				index = i;
			
		if(Data[index]>0)	// max(data[])			
			for(int i=0; i<size; i++)
				Data[i] = 1;
		else
			for(int i=0; i<size; i++)
				Data[i] = -1;
	}
	MPI_Scatter(Data, 1, MPI_INT, &data, 1, MPI_INT, 0, comm);
	
	if(my_rank==0)	
		delete [] Data;
		
	if(data>0)
		return true;
	else 
		return false;

}

/*****************************************************************************************
 *		ParaGlobalIndexing()
 * each node has a number: number_vertices; let's denote them as n0, n1, n2, n3, n4 for 5 nodes.
 * then computer the following numbers: 
 *			s0 = 0 		 for node 0
 *			s1 = s0 + n0 	 for node 1	
 * 			s2 = s1 + n1	 for node 2
 *			s3 = s2 + n2
 *			s4 = s3 + n3
 *****************************************************************************************/
void ParaGlobalIndexing(int &startindex, int data, MPI_Comm comm)
{
	//MPI_Comm  comm = MPI_COMM_WORLD;
	int size, my_rank;
	MPI_Comm_rank(comm, &my_rank);
	MPI_Comm_size(comm, &size);
	int *Data;
	if(my_rank==0)
		Data = new int[size];
	else
		Data = NULL;
	// gather the data
	MPI_Gather(&data, 1, MPI_INT, Data, 1, MPI_INT, 0, comm);
	// find the Max;
	if(my_rank==0)
	{
		int tmp = 0, sum = 0;	
		//for(int i=0; i<size; i++)
		//	printf("\n Data[%d] = %d.", i, Data[i]);
		for(int i=0; i<size; i++)
		{
			tmp 	= Data[i]; 
			Data[i]	= sum;
			sum    += tmp;
		}		
		//for(int i=0; i<size; i++)
		//	printf("\n Data[%d] = %d.", i, Data[i]);
	}
	MPI_Scatter(Data, 1, MPI_INT, &data, 1, MPI_INT, 0, comm);
	//printf("\n PID %d: data = %d", my_rank, data);
	
	if(my_rank==0)	
		delete [] Data;
		
	startindex = data;
}




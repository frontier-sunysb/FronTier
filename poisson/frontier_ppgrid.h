/*******************************************************************************
 *                      frontier_ppgrid.h
 * This file is used to send/receive 4 set of array of int/double between different 
 * nodes. It is simply a wrapper of the MPI_Sendrecv() function.
  * 
 * In fact, there is no need to define a class here! It is just an easy way of 
 * organizing the fucntions. 
 *******************************************************************************/
#ifndef _FronTierPPGrid_H
#define _FronTierPPGrid_H
#include <FronTier.h>

class FronTierPPGrid {
public:
     static void sendrecv(Front *fr, 
			  void **sendbuff, void **recvbuff, 
			  int *count, MPI_Datatype datatype);

     static void sendrecv(int dir, Front *fr,
			  void **sendbuff, void **recvbuff,
			  int count, MPI_Datatype datatype);

     static void test(Front *fr);
};

#endif

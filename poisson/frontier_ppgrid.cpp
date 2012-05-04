/*******************************************************************************
 *                  frontier_ppgrid.c
 ******************************************************************************/
#include <frontier_ppgrid.h>

/*
 * The dimension is assumed in fr->rect_grid->dim. Depending on the dimension, 
 * the sendbuff/recvbuff are assumed as the following:
 *        for 2D, sendbuff[4][sendcount], recvbuff[4][count] (4 neighbors).
 *        for 3D, sendbuff[6][sendcount], recvbuff[6][count] (6 neighbors).
 * The ordering of the neighbors are 
 *    1) x coordinate, WEST, EAST;
 *    2) y coordinate, SOUTH, NORTH;
 *    3) z coordinate, Down, Up. 
 *  
 * see also:
 *    FT_ParallelExchCellIndex();  
 */
void FronTierPPGrid::sendrecv(Front *fr, 
			      void **sendbuff, void **recvbuff, 
			      int *count, MPI_Datatype datatype)
{     
     int me[3], him[3], dim, i, j, k, index;
     int dest, source;

     INTERFACE *intfc = fr->interf;
     PP_GRID *pp_grid = fr->pp_grid;
     RECT_GRID *gr = fr->rect_grid;
     int *G = pp_grid->gmax;
     int tag = 1;        // not used 
     dim = gr->dim;

     boolean requested[6];
     MPI_Request send_request[6], recv_request[6];
     MPI_Request request[12];

     for(i=0; i<dim*2; i++)
	  requested[i] = NO;

     find_Cartesian_coordinates(pp_mynode(), pp_grid, me);
     for(i=0; i<dim; i++)
     	  for(j=0; j<2; j++)
	  {
	       for(k=0; k<dim; k++)
		    him[k] = me[k];
	       
	       index = i*2 + j;       // index for the buff.

	       if(rect_boundary_type(intfc,i,j)==SUBDOMAIN_BOUNDARY)
	       {
		    him[i] = (me[i] + 2*j - 1 + G[i]) % G[i];
		    source = dest = domain_id(him, G, dim);
		    MPI_Isend(sendbuff[index], count[index], datatype, dest, tag, 
			      MPI_COMM_WORLD, &send_request[index]);
		    MPI_Irecv(recvbuff[index], count[index], datatype, dest, tag,
			      MPI_COMM_WORLD, &recv_request[index]);
		    requested[index] = YES;
	       }
	  }

     int request_count = 0;
     MPI_Status status[12];

     for(i=0; i<dim*2; i++)
	  if(requested[i])
	  {
	       request[request_count]   = send_request[i];
	       request[request_count+1] = recv_request[i];
	       request_count += 2;
	  }
     MPI_Waitall(request_count, request, status); 
}


/*
 * sendrecv works only in a specific direction: 0 for X, 1 for Y, 2 for Z.
 *     
 * The dimension is assumed in fr->rect_grid->dim. Depending on the dimension, 
 * the sendbuff/recvbuff are assumed as the following:
 *        for 2D, sendbuff[2][count], recvbuff[2][count] (4 neighbors).
 *        for 3D, sendbuff[2][count], recvbuff[2][count] (6 neighbors).
 * The ordering of the neighbors are 
 *    1) x coordinate, WEST, EAST;
 *    2) y coordinate, SOUTH, NORTH;
 *    3) z coordinate, Down, Up. 
 *  
 * see also:
 *    FT_ParallelExchCellIndex();  
 */
void FronTierPPGrid::sendrecv(int dir, Front *fr, 
			      void **sendbuff, void **recvbuff, 
			      int count, MPI_Datatype datatype)
{     
     int me[3], him[3], dim, i, j, k, index;
     int dest, source;

     INTERFACE *intfc = fr->interf;
     PP_GRID *pp_grid = fr->pp_grid;
     RECT_GRID *gr = fr->rect_grid;
     int *G = pp_grid->gmax;
     int tag = 1;        // not used 
     dim = gr->dim;

     boolean requested[2] = {NO, NO};
     MPI_Request send_request[2], recv_request[2];
     MPI_Request request[4];

     find_Cartesian_coordinates(pp_mynode(), pp_grid, me);
    
     i = dir;
     for(j=0; j<2; j++)
     {
	  for(k=0; k<dim; k++)
	       him[k] = me[k];
	  
	  index = j;       // index for the buff.
	  
	  if(rect_boundary_type(intfc,i,j)==SUBDOMAIN_BOUNDARY)
	  {
	       him[i] = (me[i] + 2*j - 1 + G[i]) % G[i];
	       source = dest = domain_id(him, G, dim);
	       MPI_Isend(sendbuff[index], count, datatype, dest, tag, 
			 MPI_COMM_WORLD, &send_request[index]);
	       MPI_Irecv(recvbuff[index], count, datatype, dest, tag,
			 MPI_COMM_WORLD, &recv_request[index]);
	       requested[index] = YES;
	  }
     }

     int request_count = 0;
     MPI_Status status[4];

     for(i=0; i<2; i++)
	  if(requested[i])
	  {
	       request[request_count]   = send_request[i];
	       request[request_count+1] = recv_request[i];
	       request_count += 2;
	  }
     MPI_Waitall(request_count, request, status); 
}

void FronTierPPGrid::test(Front *fr)
{
     int dim = fr->rect_grid->dim;
     int N, M;
     int count[4];
     N = dim*2;
     M = 5;

     int **send, **recv;
     int i, j, mynode = pp_mynode();
     
     FT_MatrixMemoryAlloc((POINTER*)&send, N, M, sizeof(int));
     FT_MatrixMemoryAlloc((POINTER*)&recv, N, M, sizeof(int));

    

     printf("ID %d: \n", mynode);
     for(i=0; i<N; i++)
     	  for(j=0; j<M; j++)
	  {
	       send[i][j] = mynode + i*100 + j*1000;
	       recv[i][j] = -1;
	  }

     printf("ID %d: before sendrecv \n", mynode);
     for(i=0; i<N; i++)
     	  for(j=0; j<M; j++)
	  {
	       printf("ID %d: recv[%d][%d] = %4d \n", mynode, i, j, recv[i][j]);
	  }
     //for(i=0; i<N; i++)
     //  count[i] = M;
     count[0] = M;
     count[1] = M;
     sendrecv(0, fr, (void**)send, (void**)recv, count[0], MPI_INT);
     sendrecv(1, fr, (void**)send, (void**)recv, count[1], MPI_INT);
     
     printf("ID %d: after sendrecv \n", mynode);
     for(i=0; i<N; i++)
     	  for(j=0; j<M; j++)
	  {
	       printf("ID %d: recv[%d][%d] = %4d \n", mynode, i, j, recv[i][j]);
	  }
     
}

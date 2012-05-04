/**********************************************************************************************************
 *			test_SN_ELLIP.c
 **********************************************************************************************************/
#include <stdio.h>
#include <SN_ELLIP.h>
#include <VectorMatrix.h>


// copied from test_QUADTREE()
int main(int argc, char**argv)
{
	MPI_Comm comm = MPI_COMM_WORLD;
	int comm_size;
	int comm_rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(comm, &comm_rank);
	MPI_Comm_size(comm, &comm_size);
	
	printf("\n PID %d: there are %d processors in total.", comm_rank, comm_size);		
	
	int comm_row, comm_col;
	comm_row = 2;
	comm_col = 2;

	if(comm_size!=comm_row*comm_col)
	{
		printf("\n comm_size(%d)!=comm_row(%d)*comm_col(%d). ", comm_size, comm_row, comm_col);
		exit(0);
	}
	
	int **Ranks;
	NewMatrix(comm_row, comm_col, &Ranks);
	for(int i=0; i<comm_row; i++)
		for(int j=0; j<comm_col; j++)
			Ranks[i][j] = i*comm_col + j; 
	printf("\n PID %d: comm_row = %d, comm_col = %d", comm_rank, comm_row, comm_col);
	
	SN_ELLIP	ellip;

	int NRank, WRank, SRank, ERank, MRank;
	int pi, pj;
	pi = comm_rank/comm_col;
	pj = comm_rank - pi*comm_col;
	
	
	int nx, ny;
	nx = 10;
	ny = 10;
	double dx1, dx2;
	double domain_x1[2], domain_x2[2], node_x1[2], node_x2[2];
	domain_x1[0] = 0;
	domain_x1[1] = 0;
	domain_x2[0] = 1;
	domain_x2[1] = 1;
	
	dx1 = (domain_x2[0] - domain_x1[0])/comm_col;	// x
	dx2 = (domain_x2[1] - domain_x1[1])/comm_row;	// y
	
	node_x1[0] = pj*dx1 + domain_x1[0];		node_x2[0] = node_x1[0] + dx1;
	node_x1[1] = pi*dx2 + domain_x1[1];		node_x2[1] = node_x1[1] + dx2;
	
	
	
		
	if(pi==(comm_row-1))
		NRank = -1;
	else
		NRank = Ranks[pi+1][pj];
	if(pj==0)
		WRank = -1;
	else
		WRank = Ranks[pi][pj-1];
	if(pi==0)
		SRank = -1;
	else
		SRank = Ranks[pi-1][pj];
	if(pj==(comm_col-1))
		ERank = -1;
	else
		ERank = Ranks[pi][pj+1];
	
	MRank = comm_rank;
	
	Debug_Print(MRank, "\n pi=%d, pj=%d ", pi, pj);
	
	
	double **density, **phi, **fx, **fy;
	int nfluxbuffer = 3;
	
	NewMatrix(nx, ny, &density);
	NewMatrix(nx, ny, &phi);
	NewMatrix(nx+2*nfluxbuffer, ny+2*nfluxbuffer, &fx);
	NewMatrix(nx+2*nfluxbuffer, ny+2*nfluxbuffer, &fy);
	
	double x, y, dx, dy;
	dx = (node_x2[0] - node_x1[0])/nx;
	dy = (node_x2[1] - node_x1[1])/ny;
	
	for(int i=0; i<nx; i++)
		for(int j=0; j<ny; j++)
		{
			x = (i+.5)*dx + node_x1[0];
			y = (i+.5)*dy + node_x2[0];
			density[i][j] = 1;
		}
	
	
	ellip.setParameters(domain_x1, domain_x2, node_x1, node_x2, nx, ny, 1, density);
	ellip.setParaParameters(comm, comm_row, comm_col, NRank, WRank, SRank, ERank);
	ellip.solve();
	ellip.m_solver.Print_b("b");
	ellip.setFluxBuffer(nfluxbuffer);
	ellip.getFlux(fx, fy);
	ellip.printFlux(fx, fy);
	
	char filename[100];
	sprintf(filename, "solution_%d.plt", MRank);	
	ellip.print(filename);	
	
	
	
	Debug_Print(MRank, "\n test_SN_ELLIP().");
	MPI_Finalize();
	
	return 0;
}



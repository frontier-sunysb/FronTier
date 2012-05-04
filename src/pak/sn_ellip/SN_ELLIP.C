/**************************************************************************************************************
 *		SN_ELLIP.c
 **************************************************************************************************************/
 
#include <SN_ELLIP.h>

// return which nodes are first to get information from their neighbors
RECEIVER_SENDER SN_ELLIP::whoami(void)		
{
	// even numbered nodes are first Receiver.
	//int pi, pj;
	m_pj = m_comm_rank/m_comm_col;
	m_pi = m_comm_rank - m_pj*m_comm_col;
	return ((m_pi+m_pj)%2==0)?RECEIVER:SENDER;	
}

// i for x, j for y
//		
//		NORTH
//	0,m_ny		m_nx,m_ny
//
//   West			  East
//
//   	0,0		m_nx,0
//		SOUTH
//
int  SN_ELLIP::getIndex(int i, int j)
{
	return j*m_nx + i;				
}


void SN_ELLIP::setGlobalIndex(void)
{
	//Debug_Print(m_MRank, "\n SN_ELLIP::setGlobalIndex. ");
	// set local index
	int i, j;
	int number_vertices = m_nx*m_ny;
	ParaGlobalIndexing(m_starting_index, number_vertices, m_comm);
	m_ending_index = m_starting_index + number_vertices -1;
	
	//Debug_Print(m_MRank, "\n SN_ELLIP::setGlobalIndex: m_starting_index = %d, m_ending_index = %d. ", m_starting_index, m_ending_index);
	
	
	NewVector(m_ny*m_nx, &m_GlobalIndex);
	for(i=0; i<m_nx; i++)
		for(j=0; j<m_ny; j++)
			m_GlobalIndex[getIndex(i,j)] = m_starting_index + getIndex(i,j);

	// now update global index;
		
	m_whoami = whoami();
	setGlobalIndex_Boundary();
	
	if(m_whoami==RECEIVER)		
		m_whoami = SENDER;
	else
		m_whoami = RECEIVER;	
	setGlobalIndex_Boundary();	
	
	// check 
	if(m_NRank>=0)
	{
		for(i=0; i<m_nx; i++)
			if(m_BoundaryVertexIndex_North[i]<0)
				Debug_Print(m_MRank, "\n PID %d: SN_ELLIP::setGlobalIndex: negative index.", m_MRank);
	}
	if(m_WRank>=0)
	{
		for(i=0; i<m_ny; i++)
			if(m_BoundaryVertexIndex_West[i]<0)
				Debug_Print(m_MRank, "\n PID %d: SN_ELLIP::setGlobalIndex: negative index.", m_MRank);
	}
	if(m_SRank>=0)
	{
		for(i=0; i<m_nx; i++)
			if(m_BoundaryVertexIndex_South[i]<0)
				Debug_Print(m_MRank, "\n PID %d: SN_ELLIP::setGlobalIndex: negative index.", m_MRank);
	}
	if(m_ERank>=0)
	{
		for(i=0; i<m_ny; i++)
			if(m_BoundaryVertexIndex_East[i]<0)
				Debug_Print(m_MRank, "\n PID %d: SN_ELLIP::setGlobalIndex: negative index.", m_MRank);
	}
	
	
}

void SN_ELLIP::setGlobalIndex_Boundary(void)	// GB_GlobalVertexIndexing_UpdateBoundary
{
	//Debug_Print(m_MRank, "\n SN_ELLIP::setGlobalIndex_Boundary. m_whoami = %d", m_whoami);
	int i;
	int TAG=0;
        MPI_Status status[4];
        MPI_Request request[4];	
	int *Buffer[4], nBuffer[4];
	
	NewVector(m_nx, &Buffer[0]);
	NewVector(m_ny, &Buffer[1]);
	NewVector(m_nx, &Buffer[2]);
	NewVector(m_ny, &Buffer[3]);
		
	if(m_whoami==SENDER)
	{
		if(m_NRank>=0)
		{
			for(i=0; i<m_nx; i++)
				Buffer[0][i] = m_GlobalIndex[getIndex(i, m_ny-1)];
			MPI_Isend(Buffer[0], m_nx, MPI_INT, m_NRank, TAG, m_comm, &request[0]);		
		}

		if(m_WRank>=0)
		{
			for(i=0; i<m_ny; i++)
				Buffer[1][i] = m_GlobalIndex[getIndex(0, i)];
			MPI_Isend(Buffer[1], m_ny, MPI_INT, m_WRank, TAG, m_comm, &request[1]);			
		}

		if(m_SRank>=0)
		{
			for(i=0; i<m_nx; i++)
				Buffer[2][i] = m_GlobalIndex[getIndex(i, 0)];
			MPI_Isend(Buffer[2], m_nx, MPI_INT, m_SRank, TAG, m_comm, &request[2]);
		}

		if(m_ERank>=0)
		{
			for(i=0; i<m_ny; i++)
				Buffer[3][i] = m_GlobalIndex[getIndex(m_nx-1, i)];
			MPI_Isend(Buffer[3], m_ny, MPI_INT, m_ERank, TAG, m_comm, &request[3]);
		}

	}
	else	// I am a RECEIVER
	{	
		if(m_NRank>=0)
			MPI_Irecv(Buffer[0], m_nx, MPI_INT, m_NRank, TAG, m_comm, &request[0]);	
	
		if(m_WRank>=0)
			MPI_Irecv(Buffer[1], m_ny, MPI_INT, m_WRank, TAG, m_comm, &request[1]);	
	
		if(m_SRank>=0)
			MPI_Irecv(Buffer[2], m_nx, MPI_INT, m_SRank, TAG, m_comm, &request[2]);	
			
		if(m_ERank>=0)
			MPI_Irecv(Buffer[3], m_ny, MPI_INT, m_ERank, TAG, m_comm, &request[3]);	
	
	} 

	MPI_Barrier(m_comm);	// so that all processors stop here.	
	//Debug_Print(m_MRank, "\n SN_ELLIP::setGlobalIndex_Boundary: after MPI_Barrier().");		
			
	if(m_whoami==SENDER)
	{	
		// nothing;
	}
	else	// I am a RECEIVER
	{	
		if(m_NRank>=0)
			for(i=0; i<m_nx; i++)
			{
				m_BoundaryVertexIndex_North[i] = Buffer[0][i];
				//Debug_Print(m_MRank, "\n m_BoundaryVertexIndex_North[%d] = %d", i, m_BoundaryVertexIndex_North[i]);
			}

		if(m_WRank>=0)
			for(i=0; i<m_ny; i++)
			{
				m_BoundaryVertexIndex_West[i] = Buffer[1][i];
				//Debug_Print(m_MRank, "\n m_BoundaryVertexIndex_West[%d] = %d", i, m_BoundaryVertexIndex_West[i]);
			}

		if(m_SRank>=0)
			for(i=0; i<m_nx; i++)
			{
				m_BoundaryVertexIndex_South[i] = Buffer[2][i];
				//Debug_Print(m_MRank, "\n m_BoundaryVertexIndex_South[%d] = %d", i, m_BoundaryVertexIndex_South[i]);
			}
		
		if(m_ERank>=0)
			for(i=0; i<m_ny; i++)
			{
				m_BoundaryVertexIndex_East[i] = Buffer[3][i];
				//Debug_Print(m_MRank, "\n m_BoundaryVertexIndex_East[%d] = %d", i, m_BoundaryVertexIndex_East[i]);
			}
	} 
	

	DeleteVector(m_nx, &Buffer[0]);
	DeleteVector(m_ny, &Buffer[1]);
	DeleteVector(m_nx, &Buffer[2]);
	DeleteVector(m_ny, &Buffer[3]);
}

void SN_ELLIP::getGlobalSolution(void)
{
	//Debug_Print(m_MRank, "\n SN_ELLIP::getGlobalSolution. ");
		
	m_whoami = whoami();	
	getGlobalSolution_Boundary();
	
	if(m_whoami==RECEIVER)		
		m_whoami = SENDER;
	else
		m_whoami = RECEIVER;	
	getGlobalSolution_Boundary();		
}
void SN_ELLIP::getGlobalSolution_Boundary(void)
{
	//Debug_Print(m_MRank, "\n SN_ELLIP::getGlobalSolution_Boundary. m_whoami = %d", m_whoami);
	int i;
	int TAG=0;
        MPI_Status status[4];
        MPI_Request request[4];	
	double *Buffer[4], nBuffer[4];
	
	NewVector(m_nx, &Buffer[0]);
	NewVector(m_ny, &Buffer[1]);
	NewVector(m_nx, &Buffer[2]);
	NewVector(m_ny, &Buffer[3]);
		
	if(m_whoami==SENDER)
	{
		if(m_NRank>=0)
		{
			for(i=0; i<m_nx; i++)
				Buffer[0][i] = m_unknown[getIndex(i, m_ny-1)];
			MPI_Isend(Buffer[0], m_nx, MPI_DOUBLE, m_NRank, TAG, m_comm, &request[0]);		
		}

		if(m_WRank>=0)
		{
			for(i=0; i<m_ny; i++)
				Buffer[1][i] = m_unknown[getIndex(0, i)];
			MPI_Isend(Buffer[1], m_ny, MPI_DOUBLE, m_WRank, TAG, m_comm, &request[1]);			
		}

		if(m_SRank>=0)
		{
			for(i=0; i<m_nx; i++)
				Buffer[2][i] = m_unknown[getIndex(i, 0)];
			MPI_Isend(Buffer[2], m_nx, MPI_DOUBLE, m_SRank, TAG, m_comm, &request[2]);
		}

		if(m_ERank>=0)
		{
			for(i=0; i<m_ny; i++)
				Buffer[3][i] = m_unknown[getIndex(m_nx-1, i)];
			MPI_Isend(Buffer[3], m_ny, MPI_DOUBLE, m_ERank, TAG, m_comm, &request[3]);
		}

	}
	else	// I am a RECEIVER
	{	
		if(m_NRank>=0)
			MPI_Irecv(Buffer[0], m_nx, MPI_DOUBLE, m_NRank, TAG, m_comm, &request[0]);	
	
		if(m_WRank>=0)
			MPI_Irecv(Buffer[1], m_ny, MPI_DOUBLE, m_WRank, TAG, m_comm, &request[1]);	
	
		if(m_SRank>=0)
			MPI_Irecv(Buffer[2], m_nx, MPI_DOUBLE, m_SRank, TAG, m_comm, &request[2]);	
			
		if(m_ERank>=0)
			MPI_Irecv(Buffer[3], m_ny, MPI_DOUBLE, m_ERank, TAG, m_comm, &request[3]);	
	
	} 

	MPI_Barrier(m_comm);	// so that all processors stop here.	
	//Debug_Print(m_MRank, "\n SN_ELLIP::setGlobalIndex_Boundary: after MPI_Barrier().");		
			
	if(m_whoami==SENDER)
	{	
		// nothing;
	}
	else	// I am a RECEIVER
	{	
		if(m_NRank>=0)
			for(i=0; i<m_nx; i++)
			{
				m_unknown_North[i] = Buffer[0][i];
			}

		if(m_WRank>=0)
			for(i=0; i<m_ny; i++)
			{
				m_unknown_West[i] = Buffer[1][i];
			}

		if(m_SRank>=0)
			for(i=0; i<m_nx; i++)
			{
				m_unknown_South[i] = Buffer[2][i];
			}
		
		if(m_ERank>=0)
			for(i=0; i<m_ny; i++)
			{
				m_unknown_East[i] = Buffer[3][i];
			}
	} 
	

	DeleteVector(m_nx, &Buffer[0]);
	DeleteVector(m_ny, &Buffer[1]);
	DeleteVector(m_nx, &Buffer[2]);
	DeleteVector(m_ny, &Buffer[3]);
}

double SN_ELLIP::computeIntegralTheta(double r, double z, double r1, double z1)
{
	double theta, ntheta = 20, dtheta, xtheta, ytheta, sum;
	dtheta = 2*PI/ntheta;
	
	sum = 0;
	for(theta = 0; theta<2*PI; theta += dtheta)
	{
		xtheta = r1*cos(theta);
		ytheta = r1*sin(theta);
		sum += dtheta/sqrt((r-xtheta)*(r-xtheta)+(0-ytheta)*(0-ytheta)+(z-z1)*(z-z1)); 
	}
	return sum;
}
// compute the boundary condition on the east/north boundary
void SN_ELLIP::computeBoundaryPotential(void)
{
		
	int i, j, k, l;
	double r, z, r1, z1;
	
	for(k=0; k<m_comm_col*m_nx; k++)
		m_BoundaryPotential_North[k] = 0;
		
	for(k=0; k<m_comm_row*m_ny; k++)
		m_BoundaryPotential_East[k] = 0;
		
	//return;
	for(i=0; i<m_nx; i++)
		for(j=0; j<m_ny; j++)
		{
			r1 = (i+1.0/2)*m_dx + m_x1[0];
			z1 = (j+1.0/2)*m_dy + m_x1[1];
			// the north boundary			
			z = m_X2[1];				
			for(k=0; k<m_comm_col*m_nx; k++)
			{
				r = (k+1.0/2)*m_dx + m_X1[0];
				//m_BoundaryPotential_North[k] += m_G*(m_dx*m_dy*r1*m_Density[i][j])*computeIntegralTheta(r,z,r1, z1);	
				//m_BoundaryPotential_North[k] += m_G*(m_dx*m_dy*r1*m_Density[i][j])*computeIntegralTheta(r,z,r1,-z1);
				m_BoundaryPotential_North[k] += -1.0/(4*PI)*r1*func_f(i,j)*m_dx*m_dy*computeIntegralTheta(r,z,r1, z1);	
				m_BoundaryPotential_North[k] += -1.0/(4*PI)*r1*func_f(i,j)*m_dx*m_dy*computeIntegralTheta(r,z,r1, -z1);	
			}
			// the east boundary	
			r = m_X2[0];
			for(k=0; k<m_comm_row*m_ny; k++)
			{
				z = (k+1.0/2)*m_dy + m_X1[1];
				m_BoundaryPotential_East[k] += -1.0/(4*PI)*r1*func_f(i,j)*m_dx*m_dy*computeIntegralTheta(r,z,r1, z1);	
				m_BoundaryPotential_East[k] += -1.0/(4*PI)*r1*func_f(i,j)*m_dx*m_dy*computeIntegralTheta(r,z,r1, -z1);	
			}			
		}
	
	
	
	// now update using MPI		
	//Debug_Print(m_MRank, "\n m_comm_size = %d, m_X2 = {%f,%f}", m_comm_size, m_X2[0], m_X2[1]);
	if(m_comm_size>1)
	{
		for(i=0; i<m_comm_col*m_nx; i++)
			m_buffer[i] = m_BoundaryPotential_North[i];		
		MPI_Allreduce(m_buffer, m_BoundaryPotential_North, m_comm_col*m_nx, MPI_DOUBLE, MPI_SUM, m_comm);
	
		for(i=0; i<m_comm_row*m_ny; i++)
			m_buffer[i] = m_BoundaryPotential_East[i];		
		MPI_Allreduce(m_buffer, m_BoundaryPotential_East, m_comm_row*m_ny, MPI_DOUBLE, MPI_SUM, m_comm);	
	}	
	//for(i=0; i<m_comm_col*m_nx; i++)
	//	Debug_Print(m_MRank, "\n m_BoundaryPotential_North[%d] = %f", i, m_BoundaryPotential_North[i]);		

	//for(i=0; i<m_comm_row*m_ny; i++)
	//	Debug_Print(m_MRank, "\n m_BoundaryPotential_East[%d] = %f", i, m_BoundaryPotential_East[i]);			
}

//
//	  d
//	a e c	= f
//	  b 
void SN_ELLIP::setMatrix(void)
{
//	Debug_Print(m_MRank, "\n SN_ELLIP::setMatrix: (%d,%d)", m_starting_index, m_ending_index);
	m_solver.Create(m_starting_index, m_ending_index);
	int i, j, pi;
	double r, z;
	double a, b, c, d, e, f;	// variables used to facilitate the programming;
	
	
	for(i=1; i<m_nx-1; i++)
		for(j=1; j<m_ny-1; j++)
		{
			r = m_x1[0] + (i+1.0/2)*m_dx;	// the center 
			z = m_x1[1] + (j+1.0/2)*m_dy;
			pi = getIndex(i,j);
			a = r/(m_dx*m_dx) - 1.0/(2*m_dx);	b = r/(m_dy*m_dy);	c = r/(m_dx*m_dx) + 1.0/(2*m_dx);	d = r/(m_dy*m_dy);			
			e = -2*r/(m_dx*m_dx) -2*r/(m_dy*m_dy);	f = func_f(i,j)*r;
			m_solver.Set_A(m_GlobalIndex[pi], m_GlobalIndex[pi], 			e);
			m_solver.Set_A(m_GlobalIndex[pi], m_GlobalIndex[getIndex(i-1,  j)],	a);
			m_solver.Set_A(m_GlobalIndex[pi], m_GlobalIndex[getIndex(i+1,  j)],	c);
			m_solver.Set_A(m_GlobalIndex[pi], m_GlobalIndex[getIndex(  i,j-1)],	b);
			m_solver.Set_A(m_GlobalIndex[pi], m_GlobalIndex[getIndex(  i,j+1)],	d);
			
			m_solver.Set_b(m_GlobalIndex[pi], f);		
		}
	// the boundary edges
	// north, 		Neumann
	j = m_ny-1;
	for(i=1; i<m_nx-1; i++)
	{
		r = m_x1[0] + (i+1.0/2)*m_dx;	// the center 
		z = m_x1[1] + (j+1.0/2)*m_dy;
		pi = getIndex(i,j);
		a = r/(m_dx*m_dx) - 1.0/(2*m_dx);	b = r/(m_dy*m_dy);	c = r/(m_dx*m_dx) + 1.0/(2*m_dx);	d = r/(m_dy*m_dy);			
		e = -2*r/(m_dx*m_dx) -2*r/(m_dy*m_dy);	f = func_f(i,j)*r;
		
		
		m_solver.Set_A(m_GlobalIndex[pi], m_GlobalIndex[getIndex(i-1,  j)],	a);
		m_solver.Set_A(m_GlobalIndex[pi], m_GlobalIndex[getIndex(  i,j-1)],	b);
		m_solver.Set_A(m_GlobalIndex[pi], m_GlobalIndex[getIndex(i+1,  j)],	c);		
		m_solver.Set_A(m_GlobalIndex[pi], m_GlobalIndex[pi], 			e);
		m_solver.Set_b(m_GlobalIndex[pi], f);		
		if(m_NRank>=0)
			m_solver.Set_A(m_GlobalIndex[pi], m_BoundaryVertexIndex_North[i],	d);				
		else
		{
			m_solver.Add_A(m_GlobalIndex[pi], m_GlobalIndex[pi], 			-d);			
			//m_solver.Add_b(m_GlobalIndex[pi], - 2*d*func_u(r,z+m_dy/2));		
			m_solver.Add_b(m_GlobalIndex[pi], - 2*d*m_BoundaryPotential_North[i + m_pi*m_nx]);		
		}		
	}
	// west			Neumann
	i = 0;
	for(j=1; j<m_ny-1; j++)
	{
		r = m_x1[0] + (i+1.0/2)*m_dx;	// the center 
		z = m_x1[1] + (j+1.0/2)*m_dy;
		pi = getIndex(i,j);
		a = r/(m_dx*m_dx) - 1.0/(2*m_dx);	b = r/(m_dy*m_dy);	c = r/(m_dx*m_dx) + 1.0/(2*m_dx);	d = r/(m_dy*m_dy);			
		e = -2*r/(m_dx*m_dx) -2*r/(m_dy*m_dy);	f = func_f(i,j)*r;
		
		m_solver.Set_A(m_GlobalIndex[pi], m_GlobalIndex[getIndex(  i,j-1)],	b);
		m_solver.Set_A(m_GlobalIndex[pi], m_GlobalIndex[getIndex(i+1,  j)],	c);
		m_solver.Set_A(m_GlobalIndex[pi], m_GlobalIndex[getIndex(  i,j+1)],	d);
		m_solver.Set_A(m_GlobalIndex[pi], m_GlobalIndex[pi], 			e);			
		m_solver.Set_b(m_GlobalIndex[pi], f);		
		if(m_WRank>=0)
			m_solver.Set_A(m_GlobalIndex[pi], m_BoundaryVertexIndex_West[j],	a);
		else
			m_solver.Add_A(m_GlobalIndex[pi], m_GlobalIndex[pi], 			a);

	}
	
	// south		Neumann
	j = 0;
	for(i=1; i<m_nx-1; i++)
	{
		r = m_x1[0] + (i+1.0/2)*m_dx;	// the center 
		z = m_x1[1] + (j+1.0/2)*m_dy;
		pi = getIndex(i,j);
		a = r/(m_dx*m_dx) - 1.0/(2*m_dx);	b = r/(m_dy*m_dy);	c = r/(m_dx*m_dx) + 1.0/(2*m_dx);	d = r/(m_dy*m_dy);			
		e = -2*r/(m_dx*m_dx) -2*r/(m_dy*m_dy);	f = func_f(i,j)*r;
		
		
		m_solver.Set_A(m_GlobalIndex[pi], m_GlobalIndex[getIndex(i-1,  j)],	a);
		m_solver.Set_A(m_GlobalIndex[pi], m_GlobalIndex[getIndex(i+1,  j)],	c);		
		m_solver.Set_A(m_GlobalIndex[pi], m_GlobalIndex[getIndex(  i,j+1)],	d);		
		m_solver.Add_A(m_GlobalIndex[pi], m_GlobalIndex[pi], 			e);			
		m_solver.Set_b(m_GlobalIndex[pi], f);		
		if(m_SRank>=0)
			m_solver.Set_A(m_GlobalIndex[pi], m_BoundaryVertexIndex_South[i],	b);
		else
			m_solver.Add_A(m_GlobalIndex[pi], m_GlobalIndex[pi], 			b);
		
	}
	// east,		Diriclet boundary
	i = m_nx-1;
	for(j=1; j<m_ny-1; j++)
	{
		r = m_x1[0] + (i+1.0/2)*m_dx;	// the center 
		z = m_x1[1] + (j+1.0/2)*m_dy;
		pi = getIndex(i,j);
		a = r/(m_dx*m_dx) - 1.0/(2*m_dx);	b = r/(m_dy*m_dy);	c = r/(m_dx*m_dx) + 1.0/(2*m_dx);	d = r/(m_dy*m_dy);			
		e = -2*r/(m_dx*m_dx) -2*r/(m_dy*m_dy);	f = func_f(i,j)*r;
		
		m_solver.Set_A(m_GlobalIndex[pi], m_GlobalIndex[getIndex(i-1, j)],	a);
		m_solver.Set_A(m_GlobalIndex[pi], m_GlobalIndex[getIndex(  i,j-1)],	b);
		m_solver.Set_A(m_GlobalIndex[pi], m_GlobalIndex[getIndex(  i,j+1)],	d);		
		m_solver.Set_A(m_GlobalIndex[pi], m_GlobalIndex[pi], 			e);		
		m_solver.Set_b(m_GlobalIndex[pi], f);		
		if(m_ERank>=0)
			m_solver.Set_A(m_GlobalIndex[pi], m_BoundaryVertexIndex_East[j],	c);
		else
		{
			m_solver.Add_A(m_GlobalIndex[pi], m_GlobalIndex[pi], 			-c);			
			//m_solver.Add_b(m_GlobalIndex[pi], - 2*c*func_u(r+m_dx/2,z));	
			m_solver.Add_b(m_GlobalIndex[pi], - 2*c*m_BoundaryPotential_East[j+m_pj*m_ny]);				
			//m_solver.Add_b(m_GlobalIndex[pi], 0);	
		}
	
	}
	
	// north-west
	i = 0;
	j = m_ny-1;
	r = m_x1[0] + (i+1.0/2)*m_dx;	// the center 
	z = m_x1[1] + (j+1.0/2)*m_dy;
	pi = getIndex(i,j);
	a = r/(m_dx*m_dx) - 1.0/(2*m_dx);	b = r/(m_dy*m_dy);	c = r/(m_dx*m_dx) + 1.0/(2*m_dx);	d = r/(m_dy*m_dy);			
	e = -2*r/(m_dx*m_dx) -2*r/(m_dy*m_dy);	f = func_f(i,j)*r;
	
	
	m_solver.Set_A(m_GlobalIndex[pi], m_GlobalIndex[pi], 			e);
	m_solver.Set_b(m_GlobalIndex[pi], f);		
	
	if(m_WRank>=0)
		m_solver.Set_A(m_GlobalIndex[pi], m_BoundaryVertexIndex_West[j],	a);
	else
		m_solver.Add_A(m_GlobalIndex[pi], m_GlobalIndex[pi], 			a);
	m_solver.Set_A(m_GlobalIndex[pi], m_GlobalIndex[getIndex(  i,j-1)],	b);
	m_solver.Set_A(m_GlobalIndex[pi], m_GlobalIndex[getIndex(i+1,  j)],	c);	
	if(m_NRank>=0)
		m_solver.Set_A(m_GlobalIndex[pi], m_BoundaryVertexIndex_North[i],	d);
	else
	{
		m_solver.Add_A(m_GlobalIndex[pi], m_GlobalIndex[pi], 			-d);
		//m_solver.Add_b(m_GlobalIndex[pi], - 2*d*func_u(r,z+m_dy/2));				
		m_solver.Add_b(m_GlobalIndex[pi], - 2*d*m_BoundaryPotential_North[i + m_pi*m_nx]);				
	}			
	
	
	// south-west
	i = 0;
	j = 0;
	r = m_x1[0] + (i+1.0/2)*m_dx;	// the center 
	z = m_x1[1] + (j+1.0/2)*m_dy;
	pi = getIndex(i,j);
	a = r/(m_dx*m_dx) - 1.0/(2*m_dx);	b = r/(m_dy*m_dy);	c = r/(m_dx*m_dx) + 1.0/(2*m_dx);	d = r/(m_dy*m_dy);			
	e = -2*r/(m_dx*m_dx) -2*r/(m_dy*m_dy);	f = func_f(i,j)*r;
	
	m_solver.Set_A(m_GlobalIndex[pi], m_GlobalIndex[pi], 			e);	
	m_solver.Set_A(m_GlobalIndex[pi], m_GlobalIndex[getIndex(i+1,  j)],	c);
	m_solver.Set_A(m_GlobalIndex[pi], m_GlobalIndex[getIndex(  i,j+1)],	d);
	
	if(m_SRank>=0)
		m_solver.Set_A(m_GlobalIndex[pi], m_BoundaryVertexIndex_South[i],	b);	
	else
		m_solver.Add_A(m_GlobalIndex[pi], m_GlobalIndex[pi],			b);	
	

	if(m_WRank>=0)
		m_solver.Set_A(m_GlobalIndex[pi], m_BoundaryVertexIndex_West[j],	a);
	else
		m_solver.Add_A(m_GlobalIndex[pi], m_GlobalIndex[pi],			a);

	m_solver.Set_b(m_GlobalIndex[pi], f);		
	
	
	// north-east	
	i = m_nx-1;
	j = m_ny-1;
	r = m_x1[0] + (i+1.0/2)*m_dx;	// the center 
	z = m_x1[1] + (j+1.0/2)*m_dy;
	pi = getIndex(i,j);
	a = r/(m_dx*m_dx) - 1.0/(2*m_dx);	b = r/(m_dy*m_dy);	c = r/(m_dx*m_dx) + 1.0/(2*m_dx);	d = r/(m_dy*m_dy);			
	e = -2*r/(m_dx*m_dx) -2*r/(m_dy*m_dy);	f = func_f(i,j)*r;
	
	m_solver.Set_A(m_GlobalIndex[pi], m_GlobalIndex[pi], 			e);	
	m_solver.Set_A(m_GlobalIndex[pi], m_GlobalIndex[getIndex(i-1,  j)],	a);
	m_solver.Set_A(m_GlobalIndex[pi], m_GlobalIndex[getIndex(  i,j-1)],	b);
	m_solver.Set_b(m_GlobalIndex[pi], f);		
	
	if(m_NRank>=0)
		m_solver.Set_A(m_GlobalIndex[pi], m_BoundaryVertexIndex_North[i],	d);
	else
	{
		m_solver.Add_A(m_GlobalIndex[pi], m_GlobalIndex[pi], 			-d);
		//m_solver.Add_b(m_GlobalIndex[pi], - 2*d*func_u(r,z+m_dy/2));				
		m_solver.Add_b(m_GlobalIndex[pi], - 2*d*m_BoundaryPotential_North[i + m_pi*m_nx]);				
	}
	
	if(m_ERank>=0)
		m_solver.Set_A(m_GlobalIndex[pi], m_BoundaryVertexIndex_East[j],	c);
	else
	{
		m_solver.Add_A(m_GlobalIndex[pi], m_GlobalIndex[pi], 			-c);		
		//m_solver.Add_b(m_GlobalIndex[pi], - 2*c*func_u(r+m_dx/2,z));	
		m_solver.Add_b(m_GlobalIndex[pi], - 2*c*m_BoundaryPotential_East[j + m_pj*m_ny]);	
	}
	
	
	// south-east
	i = m_nx-1;
	j = 0;
	r = m_x1[0] + (i+1.0/2)*m_dx;	// the center 
	z = m_x1[1] + (j+1.0/2)*m_dy;
	pi = getIndex(i,j);
	a = r/(m_dx*m_dx) - 1.0/(2*m_dx);	b = r/(m_dy*m_dy);	c = r/(m_dx*m_dx) + 1.0/(2*m_dx);	d = r/(m_dy*m_dy);			
	e = -2*r/(m_dx*m_dx) -2*r/(m_dy*m_dy);	f = func_f(i,j)*r;
	
	m_solver.Set_A(m_GlobalIndex[pi], m_GlobalIndex[pi], 			e);	
	m_solver.Set_A(m_GlobalIndex[pi], m_GlobalIndex[getIndex(i-1,  j)],	a);
	m_solver.Set_A(m_GlobalIndex[pi], m_GlobalIndex[getIndex(  i,j+1)],	d);
	m_solver.Set_b(m_GlobalIndex[pi], f);		
	
	if(m_SRank>=0)
		m_solver.Set_A(m_GlobalIndex[pi], m_BoundaryVertexIndex_South[i],	b);
	else
		m_solver.Add_A(m_GlobalIndex[pi], m_GlobalIndex[pi], 			b);
	
	if(m_ERank>=0)
		m_solver.Set_A(m_GlobalIndex[pi], m_BoundaryVertexIndex_East[j],	c);
	else
	{
		m_solver.Add_A(m_GlobalIndex[pi], m_GlobalIndex[pi], 			-c);		
		//m_solver.Add_b(m_GlobalIndex[pi], - 2*c*func_u(r+m_dx/2,z));	
		m_solver.Add_b(m_GlobalIndex[pi], - 2*c*m_BoundaryPotential_East[j + m_pj*m_ny]);	
	}
}

// see also SN_ELLIP::getGlobalSolution_Boundary()
void SN_ELLIP::getFlux_Boundary(double **fx, double **fy)
{
	//Debug_Print(m_MRank, "\n SN_ELLIP::getFlux_Boundary");
	int i, j, index;
	int TAG=0;
        MPI_Status status[4];
        MPI_Request request[4];	
	double *Buffer[4]; int nBuffer[4];
	
	nBuffer[0] = (m_nx+m_nfluxbuffer*2)*m_nfluxbuffer*2;
	NewVector(nBuffer[0], &Buffer[0]);		// north
	nBuffer[1] = (m_ny+m_nfluxbuffer*2)*m_nfluxbuffer*2;
	NewVector(nBuffer[1], &Buffer[1]);		// west
	nBuffer[2] = (m_nx+m_nfluxbuffer*2)*m_nfluxbuffer*2;
	NewVector(nBuffer[2], &Buffer[2]);		// south
	nBuffer[3] = (m_ny+m_nfluxbuffer*2)*m_nfluxbuffer*2;
	NewVector(nBuffer[3], &Buffer[3]);		// east
		
	if(m_whoami==SENDER)
	{
		if(m_NRank>=0)
		{
			index = 0;
			for(i=0; i<m_nx+m_nfluxbuffer*2; i++)
				for(j=m_ny; j<m_ny+m_nfluxbuffer; j++)
					Buffer[0][index++] = fx[i][j];			
			for(i=0; i<m_nx+m_nfluxbuffer*2; i++)
				for(j=m_ny; j<m_ny+m_nfluxbuffer; j++)
					Buffer[0][index++] = fy[i][j];	
			MPI_Isend(Buffer[0], index, MPI_DOUBLE, m_NRank, TAG, m_comm, &request[0]);		
		}

		if(m_WRank>=0)
		{
			index = 0;
			for(i=m_nfluxbuffer; i<m_nfluxbuffer+m_nfluxbuffer; i++)
				for(j=0; j<m_ny+m_nfluxbuffer*2; j++)
					Buffer[1][index++] = fx[i][j];			
			for(i=m_nfluxbuffer; i<m_nfluxbuffer+m_nfluxbuffer; i++)
				for(j=0; j<m_ny+m_nfluxbuffer*2; j++)
					Buffer[1][index++] = fy[i][j];	
			MPI_Isend(Buffer[1], index, MPI_DOUBLE, m_WRank, TAG, m_comm, &request[1]);			
		}

		if(m_SRank>=0)
		{
			index = 0;
			for(i=0; i<m_nx+m_nfluxbuffer*2; i++)
				for(j=m_nfluxbuffer; j<m_nfluxbuffer+m_nfluxbuffer; j++)
					Buffer[2][index++] = fx[i][j];			
			for(i=0; i<m_nx+m_nfluxbuffer*2; i++)
				for(j=m_nfluxbuffer; j<m_nfluxbuffer+m_nfluxbuffer; j++)
					Buffer[2][index++] = fy[i][j];	
			MPI_Isend(Buffer[2], index, MPI_DOUBLE, m_SRank, TAG, m_comm, &request[2]);
		}

		if(m_ERank>=0)
		{
			index = 0;
			for(i=m_nx; i<m_nx+m_nfluxbuffer; i++)
				for(j=0; j<m_ny+m_nfluxbuffer*2; j++)
					Buffer[3][index++] = fx[i][j];			
			for(i=m_nx; i<m_nx+m_nfluxbuffer; i++)
				for(j=0; j<m_ny+m_nfluxbuffer*2; j++)
					Buffer[3][index++] = fy[i][j];	
			MPI_Isend(Buffer[3], index, MPI_DOUBLE, m_ERank, TAG, m_comm, &request[3]);
		}

	}
	else	// I am a RECEIVER
	{	
		if(m_NRank>=0)
			MPI_Irecv(Buffer[0], nBuffer[0], MPI_DOUBLE, m_NRank, TAG, m_comm, &request[0]);	
	
		if(m_WRank>=0)
			MPI_Irecv(Buffer[1], nBuffer[1], MPI_DOUBLE, m_WRank, TAG, m_comm, &request[1]);	
	
		if(m_SRank>=0)
			MPI_Irecv(Buffer[2], nBuffer[2], MPI_DOUBLE, m_SRank, TAG, m_comm, &request[2]);	
			
		if(m_ERank>=0)
			MPI_Irecv(Buffer[3], nBuffer[3], MPI_DOUBLE, m_ERank, TAG, m_comm, &request[3]);	
	
	} 

	MPI_Barrier(m_comm);	// so that all processors stop here.	
	//Debug_Print(m_MRank, "\n SN_ELLIP::setGlobalIndex_Boundary: after MPI_Barrier().");		
			
	if(m_whoami==SENDER)
	{	
		// nothing;
	}
	else	// I am a RECEIVER
	{	
		if(m_NRank>=0)
		{
			index = 0;
			for(i=0; i<m_nx+m_nfluxbuffer*2; i++)
				for(j=m_ny+m_nfluxbuffer; j<m_ny+m_nfluxbuffer*2; j++)
					fx[i][j] = Buffer[0][index++];			
			for(i=0; i<m_nx+m_nfluxbuffer*2; i++)
				for(j=m_ny+m_nfluxbuffer; j<m_ny+m_nfluxbuffer*2; j++)
					fy[i][j] = Buffer[0][index++];	
		}
		else
		{
			index = 0;
			for(i=0; i<m_nx+m_nfluxbuffer*2; i++)
				for(j=m_ny+m_nfluxbuffer; j<m_ny+m_nfluxbuffer*2; j++)
					fx[i][j] = 0;			
			for(i=0; i<m_nx+m_nfluxbuffer*2; i++)
				for(j=m_ny+m_nfluxbuffer; j<m_ny+m_nfluxbuffer*2; j++)
					fy[i][j] = 0;	
		}

		if(m_WRank>=0)
		{
			index = 0;
			for(i=0; i<m_nfluxbuffer; i++)
				for(j=0; j<m_ny+m_nfluxbuffer*2; j++)
					fx[i][j] = Buffer[1][index++];			
			for(i=0; i<m_nfluxbuffer; i++)
				for(j=0; j<m_ny+m_nfluxbuffer*2; j++)
					fy[i][j] = Buffer[1][index++];				
		}
		else
		{
			index = 0;
			for(i=0; i<m_nfluxbuffer; i++)
				for(j=0; j<m_ny+m_nfluxbuffer*2; j++)
					fx[i][j] = - fx[2*m_nfluxbuffer-1-i][j];			
			for(i=0; i<m_nfluxbuffer; i++)
				for(j=0; j<m_ny+m_nfluxbuffer*2; j++)
					fy[i][j] = fy[2*m_nfluxbuffer-1-i][j];				
		}

		if(m_SRank>=0)
		{
			index = 0;
			for(i=0; i<m_nx+m_nfluxbuffer*2; i++)
				for(j=0; j<m_nfluxbuffer; j++)
					fx[i][j] = Buffer[2][index++];			
			for(i=0; i<m_nx+m_nfluxbuffer*2; i++)
				for(j=0; j<m_nfluxbuffer; j++)
					fy[i][j] = Buffer[2][index++];	
		}
		else
		{
			index = 0;
			for(i=0; i<m_nx+m_nfluxbuffer*2; i++)
				for(j=0; j<m_nfluxbuffer; j++)
					fx[i][j] = fx[i][2*m_nfluxbuffer-1-j];			
			for(i=0; i<m_nx+m_nfluxbuffer*2; i++)
				for(j=0; j<m_nfluxbuffer; j++)
					fy[i][j] = - fy[i][2*m_nfluxbuffer-1-j];
		}
		
		if(m_ERank>=0)
		{
			index = 0;
			for(i=m_nx+m_nfluxbuffer; i<m_nx+m_nfluxbuffer*2; i++)
				for(j=0; j<m_ny+m_nfluxbuffer*2; j++)
					fx[i][j] = Buffer[3][index++];			
			for(i=m_nx+m_nfluxbuffer; i<m_nx+m_nfluxbuffer*2; i++)
				for(j=0; j<m_ny+m_nfluxbuffer*2; j++)
					fy[i][j] = Buffer[3][index++];				
		}
		else
		{
			index = 0;
			for(i=m_nx+m_nfluxbuffer; i<m_nx+m_nfluxbuffer*2; i++)
				for(j=0; j<m_ny+m_nfluxbuffer*2; j++)
					fx[i][j] = 0;			
			for(i=m_nx+m_nfluxbuffer; i<m_nx+m_nfluxbuffer*2; i++)
				for(j=0; j<m_ny+m_nfluxbuffer*2; j++)
					fy[i][j] = 0;				
		}
	} 
	

	DeleteVector(m_nx, &Buffer[0]);
	DeleteVector(m_ny, &Buffer[1]);
	DeleteVector(m_nx, &Buffer[2]);
	DeleteVector(m_ny, &Buffer[3]);
}


SN_ELLIP::SN_ELLIP()
{
	m_nx 		= -1;
	m_ny 		= -1;
	m_comm		= MPI_COMM_WORLD;
	m_comm_col	= 1;
	m_comm_row	= 1;
	m_comm_size	= 1;
	
	m_NRank		= -1;
	m_WRank		= -1;
	m_SRank		= -1;
	m_ERank		= -1;
	
	m_BoundaryVertexIndex_North 	= NULL;
	m_BoundaryVertexIndex_West 	= NULL;
	m_BoundaryVertexIndex_South 	= NULL;
	m_BoundaryVertexIndex_East 	= NULL;
	m_BoundaryPotential_North	= NULL;
	m_BoundaryPotential_East	= NULL;
	
	m_GlobalIndex	= NULL;	
	m_unknown 	= NULL;
	
	m_unknown_North = NULL;
	m_unknown_West  = NULL;
	m_unknown_South = NULL;
	m_unknown_East  = NULL;
	
	m_nfluxbuffer 	= 0;
	
}
SN_ELLIP::~SN_ELLIP()
{
	if(m_NRank>=0)
	{
		DeleteVector(m_nx, &m_BoundaryVertexIndex_North);
		DeleteVector(m_nx, &m_unknown_North);
	}
	if(m_WRank>=0)
	{
		DeleteVector(m_ny, &m_BoundaryVertexIndex_West);
		DeleteVector(m_ny, &m_unknown_West);		
	}
	if(m_SRank>=0)
	{
		DeleteVector(m_ny, &m_BoundaryVertexIndex_South);
		DeleteVector(m_nx, &m_unknown_South);
	}
	if(m_ERank>=0)
	{
		DeleteVector(m_ny, &m_BoundaryVertexIndex_East);
		DeleteVector(m_ny, &m_unknown_East);
	}
		
	DeleteVector(m_ny*m_nx, &m_GlobalIndex);
	DeleteVector(m_ny*m_nx, &m_unknown);
	DeleteVector(m_comm_col*m_nx, &m_BoundaryPotential_North);
	DeleteVector(m_comm_row*m_ny, &m_BoundaryPotential_East);		
	
	delete [] m_buffer;
}
	


void SN_ELLIP::setParameters(double X1[2], double X2[2], double x1[2], double x2[2], int nx, int ny, double g, double **density)
{
	m_X1[0] = X1[0];
	m_X1[1] = X1[1];
	m_X2[0] = X2[0];
	m_X2[1] = X2[1];
	
	m_x1[0] = x1[0];
	m_x1[1] = x1[1];
	m_x2[0] = x2[0];
	m_x2[1] = x2[1];
	m_nx = nx;
	m_ny = ny;
	m_dx = (m_x2[0] - m_x1[0])/m_nx;
	m_dy = (m_x2[1] - m_x1[1])/m_ny;
	
	m_G = g;
	m_Density = density;
	
	NewVector(m_ny*m_nx, &m_GlobalIndex);
	NewVector(m_ny*m_nx, &m_unknown);
}
void SN_ELLIP::setParaParameters(MPI_Comm comm, int row, int col, int nrank, int wrank, int srank, int erank)
{
	m_comm 		= comm;
	m_comm_col	= col;
	m_comm_row	= row;
	MPI_Comm_size(comm, &m_comm_size);
	MPI_Comm_rank(comm, &m_comm_rank);
	m_NRank 	= nrank;
	m_WRank		= wrank;
	m_SRank		= srank;
	m_ERank		= erank;
	MPI_Comm_rank(comm, &m_MRank);
	
	if(m_NRank>=0)	
	{
		NewVector(m_nx, &m_BoundaryVertexIndex_North);
		NewVector(m_nx, &m_unknown_North);		
	}	
	if(m_WRank>=0)	
	{
		NewVector(m_nx, &m_BoundaryVertexIndex_West);
		NewVector(m_nx, &m_unknown_West);		
	}
	if(m_SRank>=0)	
	{
		NewVector(m_nx, &m_BoundaryVertexIndex_South);
		NewVector(m_nx, &m_unknown_South);		
	}
	if(m_ERank>=0)	
	{
		NewVector(m_nx, &m_BoundaryVertexIndex_East);
		NewVector(m_nx, &m_unknown_East);		
	}
	
	NewVector(m_comm_col*m_nx, &m_BoundaryPotential_North);
	NewVector(m_comm_row*m_ny, &m_BoundaryPotential_East);
	
	m_buffer = new double[Max(m_comm_col*m_nx, m_comm_row*m_ny)];
	
	Debug_Print(m_MRank, "\n NRank=%d, WRank=%d, SRank=%d, ERank=%d", m_NRank, m_WRank, m_SRank, m_ERank);
	Debug_Print(m_MRank, "\n m_x1={%f,%f}, m_x2={%f,%f}", m_x1[0], m_x1[1], m_x2[0], m_x2[1]);
	Debug_Print(m_MRank, "\n m_nx=%d, m_ny=%d", m_nx, m_ny);
	Debug_Print(m_MRank, "\n m_dx=%f, m_dy=%f", m_dx, m_dy);
	Debug_Print(m_MRank, "\n whoami = %d", whoami());
	
	whoami();	// set m_pi, m_pj
}


double SN_ELLIP::func_f(int i, int j)
{	
	return 4*PI*m_G*m_Density[i][j];
}
 
void SN_ELLIP::solve(void)
{
	//Debug_Print(m_MRank, "\n SN_ELLIP::solve. ");
	setGlobalIndex();
	computeBoundaryPotential();
	setMatrix();
	
	//m_solver.Print_A("A");
	//m_solver.Print_b("b");
	//Debug_Print(m_MRank, "\n after setMatrix(). ");
	MPI_Barrier(m_comm);	// so that all processors stop here.	
	m_solver.Solve();
	
	//Debug_Print(m_MRank, "\n. m_solver.Solve(). ");
	NewVector(m_nx*m_ny, &m_unknown);		// DeleteVector() will be called in FEM2D::~FEM2D();
	m_solver.Get_x(m_unknown);		
	getGlobalSolution();				// get unknowns from neighbours;
}

void SN_ELLIP::getPotential(double **phi)
{
	int i, j;
	for(i=0; i<m_nx; i++)
		for(j=0; j<m_ny; j++)	
			phi[i][j] = m_unknown[getIndex(i,j)];
}

void SN_ELLIP::setFluxBuffer(int n)
{
	m_nfluxbuffer = n;
}

// fx[m_nx+2m_nfluxbuffer][m_ny+2m_nfluxbuffer], fy[m_nx+2m_nfluxbuffer][m_ny+2m_nfluxbuffer]
void SN_ELLIP::getFlux(double **fx, double **fy)		
{
	//Debug_Print(m_MRank, "\n SN_ELLIP::getFlux");
	int i, j, bi, bj;
	for(bi=m_nfluxbuffer; bi<m_nx+m_nfluxbuffer; bi++)
		for(bj=m_nfluxbuffer; bj<m_ny+m_nfluxbuffer; bj++)
		{
			i = bi - m_nfluxbuffer;
			j = bj - m_nfluxbuffer;
			//fx[bi][bj] = m_MRank+1;		// debug
			
			// fx
			if(bi==m_nfluxbuffer)
			{
				if(m_WRank>=0)
					fx[bi][bj] = (m_unknown[getIndex(i+1,j)] - m_unknown_West[j]) / (2*m_dx);
				else
					fx[bi][bj] = (m_unknown[getIndex(i+1,j)] - m_unknown[getIndex(i,j)]) / (2*m_dx);	// need to change 
			}
			else if(bi==(m_nx+m_nfluxbuffer-1))
			{
				if(m_ERank>=0)
					fx[bi][bj] = (m_unknown_East[j] - m_unknown[getIndex(i-1,j)]) / (2*m_dx);
				else
					fx[bi][bj] = (m_unknown[getIndex(i,j)] - m_unknown[getIndex(i-1,j)]) / (m_dx);	// need to change 
			}
			else
				fx[bi][bj] = (m_unknown[getIndex(i+1,j)] - m_unknown[getIndex(i-1,j)]) / (2*m_dx);
				
			// fy
			if(bj==m_nfluxbuffer)
			{
				if(m_SRank>=0)
					fy[bi][bj] = (m_unknown[getIndex(i,j+1)] - m_unknown_South[i]) / (2*m_dy);
				else
				//	fy[i][j] = (m_unknown[getIndex(i,j+1)] - m_unknown[getIndex(i,j)]) / (m_dy);	// need to change 
					fy[bi][bj] = (m_unknown[getIndex(i,j+1)] - m_unknown[getIndex(i,j)]) / (2*m_dy);	// need to change 
			}
			else if(bj==(m_ny+m_nfluxbuffer-1))
			{
				if(m_NRank>=0)
					fy[bi][bj] = (m_unknown_North[i] - m_unknown[getIndex(i,j-1)]) / (2*m_dy);
				else
					fy[bi][bj] = (m_unknown[getIndex(i,j)] - m_unknown[getIndex(i,j-1)]) / (m_dy);	// need to change 
			}
			else
				fy[bi][bj] = (m_unknown[getIndex(i,j+1)] - m_unknown[getIndex(i,j-1)]) / (2*m_dy);
		}	
	
	// now update the buffer regions;
	m_whoami = whoami();	
	getFlux_Boundary(fx, fy);
	
	if(m_whoami==RECEIVER)		
		m_whoami = SENDER;
	else
		m_whoami = RECEIVER;	
	getFlux_Boundary(fx, fy);
	
	if(m_whoami==RECEIVER)		
		m_whoami = SENDER;
	else
		m_whoami = RECEIVER;	
	getFlux_Boundary(fx, fy);
}

void SN_ELLIP::print(char *filename)
{
	FILE *hfile = fopen(filename,"w");
	if(hfile==NULL)
	{
		printf("\n can't open %s in SN_ELLIP::print", filename);
		exit(0);
	}
	int i, j;
	
	// header
	fprintf(hfile, "TITLE = %s	\n", filename);
	fprintf(hfile, "VARIABLES = \"X\" \"Y\" \"Z\" \n");
	fprintf(hfile, "ZONE I=%d J=%d F=POINT	\n", m_ny, m_nx);		// n, m
	// now the data
	int index;
	for(i=0; i<m_nx; i++)
	{
		for(j=0; j<m_ny; j++)
		{
			index = getIndex(i,j);
			fprintf(hfile,"%6d %6d %f\n", i, j, m_unknown[index]);
		}
		//fprintf(hfile, "%6d %6d %f\n", i, m_ny, m_BoundaryPotential_North[i]);
	}
		
	//for(j=0; j<m_ny; j++)
	//	fprintf(hfile, "%6d %6d %f\n", m_nx, j, m_BoundaryPotential_East[j]);	
	//fprintf(hfile, "%6d %6d %f\n", m_nx, m_ny, (m_BoundaryPotential_East[m_ny-1]+m_BoundaryPotential_North[m_nx-1])/2);
	fclose(hfile);
}

void SN_ELLIP::printFlux(double **fx, double **fy)
{
//	Debug_Print(m_MRank, "\n SN_ELLIP::printFlux");
	char filename[100];
	int i, j;
	FILE *hfile;
	
	// fx
	sprintf(filename, "fx_%d.plt", m_MRank);
	hfile = fopen(filename,"w");
	if(hfile==NULL)
	{
		printf("\n can't open %s in SN_ELLIP::printFlux", filename);
		exit(0);
	}
	
	
	// header
	fprintf(hfile, "TITLE = %s	\n", filename);
	fprintf(hfile, "VARIABLES = \"X\" \"Y\" \"Z\" \n");
	fprintf(hfile, "ZONE I=%d J=%d F=POINT	\n", m_ny+m_nfluxbuffer*2, m_nx+m_nfluxbuffer*2);		// n, m
	// now the data
	for(i=0; i<m_nx+m_nfluxbuffer*2; i++)
		for(j=0; j<m_ny+m_nfluxbuffer*2; j++)
			fprintf(hfile,"%6d %6d %f\n", i, j, fx[i][j]);
		
	fclose(hfile);
	
	sprintf(filename, "fy_%d.plt", m_MRank);
	hfile = fopen(filename,"w");
	if(hfile==NULL)
	{
		printf("\n can't open %s in SN_ELLIP::printFlux", filename);
		exit(0);
	}
	
	
	// header
	fprintf(hfile, "TITLE = %s	\n", filename);
	fprintf(hfile, "VARIABLES = \"X\" \"Y\" \"Z\" \n");
	fprintf(hfile, "ZONE I=%d J=%d F=POINT	\n", m_ny+m_nfluxbuffer*2, m_nx+m_nfluxbuffer*2);		// n, m
	// now the data
	for(i=0; i<m_nx+m_nfluxbuffer*2; i++)
		for(j=0; j<m_ny+m_nfluxbuffer*2; j++)
			fprintf(hfile,"%6d %6d %f\n", i, j, fy[i][j]);
		
	fclose(hfile);
}

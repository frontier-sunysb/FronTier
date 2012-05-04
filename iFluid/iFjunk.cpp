/*******************************************************************
 * 			iFjunk.cpp	
 *******************************************************************/
#include "iFluid.h"
#include "solver.h"

//--------------------------------------------------------------------------
// 		   Incompress_Solver_Smooth_3D_Cartesian	
//--------------------------------------------------------------------------

void Incompress_Solver_Smooth_3D_Cartesian::
	compDiffWithSmoothProperty_1st_coupled(void)
{
        COMPONENT comp;
        int index,index_nb[18],size;
        int I,I_nb[18];
        double coords[MAXD],crx_coords[MAXD];
	double coeff[18],mu[6],mu_edge[6],mu0,rho,rhs,U0_nb[18],U1_nb[18],
	       U2_nb[18];
	int flag[6]; //denote whether this is dirichlet or neumann boundary
        L_STATE state;
        int i,j,k,l,nb,icoords[MAXD];
        double speed;
	double *x;
	GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
	POINTER intfc_state;
	HYPER_SURF *hs;
	PetscInt num_iter;
	double rel_residual;

	max_speed = 0.0;
	setIndexMap();

	size = iupper - ilower;
	FT_VectorMemoryAlloc((POINTER*)&x,3*size,sizeof(double));

        PETSc solver;
        solver.Create(3*ilower, 3*iupper-1, 15, 15);
	// 7u + 4v + 4w for the first equation
	// 7v + 4u + 4w for the second equation
	// 7w + 4u + 4v for the third equation

	solver.Reset_A();
	solver.Reset_b();
	solver.Reset_x();

	for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I  = ijk_to_I[i][j][k];
            if (I == -1) continue;

            index  = d_index3d(i,j,k,top_gmax);	
	//6 neighbours of the center cell
            index_nb[0] = d_index3d(i-1,j,k,top_gmax);
            index_nb[1] = d_index3d(i+1,j,k,top_gmax);
            index_nb[2] = d_index3d(i,j-1,k,top_gmax);
            index_nb[3] = d_index3d(i,j+1,k,top_gmax);
	    index_nb[4] = d_index3d(i,j,k-1,top_gmax);
	    index_nb[5] = d_index3d(i,j,k+1,top_gmax);

	//xy cut neighbours
	    index_nb[6] = d_index3d(i-1,j-1,k,top_gmax);
	    index_nb[7] = d_index3d(i+1,j-1,k,top_gmax);
	    index_nb[8] = d_index3d(i+1,j+1,k,top_gmax);
	    index_nb[9] = d_index3d(i-1,j+1,k,top_gmax);
	//yz cut neighbours
	    index_nb[10] = d_index3d(i,j-1,k-1,top_gmax);
	    index_nb[11] = d_index3d(i,j+1,k-1,top_gmax);
	    index_nb[12] = d_index3d(i,j+1,k+1,top_gmax);
	    index_nb[13] = d_index3d(i,j-1,k+1,top_gmax);
	//xz cut neighbours
	    index_nb[14] = d_index3d(i-1,j,k-1,top_gmax);
	    index_nb[15] = d_index3d(i+1,j,k-1,top_gmax);
	    index_nb[16] = d_index3d(i+1,j,k+1,top_gmax);
	    index_nb[17] = d_index3d(i-1,j,k+1,top_gmax);

        
	//6 neighbours of the center cell
            I_nb[0] = ijk_to_I[i-1][j][k];
            I_nb[1] = ijk_to_I[i+1][j][k];
            I_nb[2] = ijk_to_I[i][j-1][k];
            I_nb[3] = ijk_to_I[i][j+1][k];
            I_nb[4] = ijk_to_I[i][j][k-1];
            I_nb[5] = ijk_to_I[i][j][k+1];
	
	//xy cut neighbours
            I_nb[6] = ijk_to_I[i-1][j-1][k];
	    I_nb[7] = ijk_to_I[i+1][j-1][k];
	    I_nb[8] = ijk_to_I[i+1][j+1][k];
	    I_nb[9] = ijk_to_I[i-1][j+1][k];
	//yz cut neighbours
	    I_nb[10] = ijk_to_I[i][j-1][k-1];
	    I_nb[11] = ijk_to_I[i][j+1][k-1];
	    I_nb[12] = ijk_to_I[i][j+1][k+1];
	    I_nb[13] = ijk_to_I[i][j-1][k+1];
	//xz cut neighbours
	    I_nb[14] = ijk_to_I[i-1][j][k-1];
	    I_nb[15] = ijk_to_I[i+1][j][k-1];
	    I_nb[16] = ijk_to_I[i+1][j][k+1];
	    I_nb[17] = ijk_to_I[i-1][j][k+1];

	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;
	    comp = top_comp[index];

	    mu0 = cell_center[index].m_state.m_mu;
	    rho = cell_center[index].m_state.m_rho;
            for (nb = 0; nb < 6; nb++)
            {
                if (FT_StateStructAtGridCrossing(front,icoords,dir[nb],
                            comp,&intfc_state,&hs,crx_coords) &&
                            wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	       	{
		    flag[nb] = 1;
		    U0_nb[nb] = getStateVel[0](intfc_state);
		    U1_nb[nb] = getStateVel[1](intfc_state);
		    U2_nb[nb] = getStateVel[2](intfc_state);

		    if (wave_type(hs) == DIRICHLET_BOUNDARY || 
			wave_type(hs) == NEUMANN_BOUNDARY)
		    {
			mu[nb] = mu0;
			mu_edge[nb] = mu0;
		    }
			
		    else
		    {
			mu[nb] = 0.5*(mu0 + 
				cell_center[index_nb[nb]].m_state.m_mu);
			mu_edge[nb] = cell_center[index_nb[nb]].m_state.m_mu;
		    }
		}
                    else
		    {
			flag[nb] = 0;
                    	U0_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[0];
			U1_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[1];
			U2_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[2];
			mu_edge[nb] = 0.5*(mu0 + 
				cell_center[index_nb[nb]].m_state.m_mu);
			mu[nb] = cell_center[index_nb[nb]].m_state.m_mu;
		    }
	    }

	    //xy cut neighbour cells
	    
	    if (flag[0] == 1 && flag[2] == 1) //neighbourcell 6
	    {
		U0_nb[6] = 0.5*(U0_nb[0] + U0_nb[2]);
		U1_nb[6] = 0.5*(U1_nb[0] + U1_nb[2]);
		U2_nb[6] = 0.5*(U2_nb[0] + U2_nb[2]);
	    }
	    else
	    {
		U0_nb[6] = cell_center[index_nb[6]].m_state.m_U[0];
		U1_nb[6] = cell_center[index_nb[6]].m_state.m_U[1];
		U2_nb[6] = cell_center[index_nb[6]].m_state.m_U[2];
	    }
	    if (flag[1] == 1 && flag[2] == 1) //neighbourcell 7 
	    {
		U0_nb[7] = 0.5*(U0_nb[1] + U0_nb[2]);
		U1_nb[7] = 0.5*(U1_nb[1] + U1_nb[2]);
		U2_nb[7] = 0.5*(U2_nb[1] + U2_nb[2]);
	    }
	    else
	    {
		U0_nb[7] = cell_center[index_nb[7]].m_state.m_U[0];
		U1_nb[7] = cell_center[index_nb[7]].m_state.m_U[1];
		U2_nb[7] = cell_center[index_nb[7]].m_state.m_U[2];
	    }
	    if (flag[1] == 1 && flag[3] == 1) //neighbourcell 8 
	    {
		U0_nb[8] = 0.5*(U0_nb[1] + U0_nb[3]);
		U1_nb[8] = 0.5*(U1_nb[1] + U1_nb[3]);
		U2_nb[8] = 0.5*(U2_nb[1] + U2_nb[3]);
	    }
	    else
	    {
		U0_nb[8] = cell_center[index_nb[8]].m_state.m_U[0];
		U1_nb[8] = cell_center[index_nb[8]].m_state.m_U[1];
		U2_nb[8] = cell_center[index_nb[8]].m_state.m_U[2];
	    }
	    if (flag[0] == 1 && flag[3] == 1) //neighbourcell 9 
	    {
		U0_nb[9] = 0.5*(U0_nb[0] + U0_nb[3]);
		U1_nb[9] = 0.5*(U1_nb[0] + U1_nb[3]);
		U2_nb[9] = 0.5*(U2_nb[0] + U2_nb[3]);
	    }
	    else
	    {
		U0_nb[9] = cell_center[index_nb[9]].m_state.m_U[0];
		U1_nb[9] = cell_center[index_nb[9]].m_state.m_U[1];
		U2_nb[9] = cell_center[index_nb[9]].m_state.m_U[2];
	    }
	    //yz cut neighbours
	    if (flag[2] == 1 && flag[4] == 1) //neighbourcell 10 
	    {
		U0_nb[10] = 0.5*(U0_nb[2] + U0_nb[4]);
		U1_nb[10] = 0.5*(U1_nb[2] + U1_nb[4]);
		U2_nb[10] = 0.5*(U2_nb[2] + U2_nb[4]);
	    }
	    else
	    {
		U0_nb[10] = cell_center[index_nb[10]].m_state.m_U[0];
		U1_nb[10] = cell_center[index_nb[10]].m_state.m_U[1];
		U2_nb[10] = cell_center[index_nb[10]].m_state.m_U[2];
	    }
	    if (flag[3] == 1 && flag[4] == 1) //neighbourcell 11 
	    {
		U0_nb[11] = 0.5*(U0_nb[3] + U0_nb[4]);
		U1_nb[11] = 0.5*(U1_nb[3] + U1_nb[4]);
		U2_nb[11] = 0.5*(U2_nb[3] + U2_nb[4]);
	    }
	    else
	    {
		U0_nb[11] = cell_center[index_nb[11]].m_state.m_U[0];
		U1_nb[11] = cell_center[index_nb[11]].m_state.m_U[1];
		U2_nb[11] = cell_center[index_nb[11]].m_state.m_U[2];
	    }
	    if (flag[3] == 1 && flag[5] == 1) //neighbourcell 12 
	    {
		U0_nb[12] = 0.5*(U0_nb[3] + U0_nb[5]);
		U1_nb[12] = 0.5*(U1_nb[3] + U1_nb[5]);
		U2_nb[12] = 0.5*(U2_nb[3] + U2_nb[5]);
	    }
	    else
	    {
		U0_nb[12] = cell_center[index_nb[12]].m_state.m_U[0];
		U1_nb[12] = cell_center[index_nb[12]].m_state.m_U[1];
		U2_nb[12] = cell_center[index_nb[12]].m_state.m_U[2];
	    }
	    if (flag[2] == 1 && flag[5] == 1) //neighbourcell 13
	    {
		U0_nb[13] = 0.5*(U0_nb[2] + U0_nb[5]);
		U1_nb[13] = 0.5*(U1_nb[2] + U1_nb[5]);
		U2_nb[13] = 0.5*(U2_nb[2] + U2_nb[5]);
	    }
	    else
	    {
		U0_nb[13] = cell_center[index_nb[13]].m_state.m_U[0];
		U1_nb[13] = cell_center[index_nb[13]].m_state.m_U[1];
		U2_nb[13] = cell_center[index_nb[13]].m_state.m_U[2];
	    }
	    //xz cut neighbours
	    if (flag[0] == 1 && flag[4] == 1) //neighbourcell 14
	    {
		U0_nb[14] = 0.5*(U0_nb[0] + U0_nb[4]);
		U1_nb[14] = 0.5*(U1_nb[0] + U1_nb[4]);
		U2_nb[14] = 0.5*(U2_nb[0] + U2_nb[4]);
	    }
	    else
	    {
		U0_nb[14] = cell_center[index_nb[14]].m_state.m_U[0];
		U1_nb[14] = cell_center[index_nb[14]].m_state.m_U[1];
		U2_nb[14] = cell_center[index_nb[14]].m_state.m_U[2];
	    }
	    if (flag[1] == 1 && flag[4] == 1) //neighbourcell 15
	    {
		U0_nb[15] = 0.5*(U0_nb[1] + U0_nb[4]);
		U1_nb[15] = 0.5*(U1_nb[1] + U1_nb[4]);
		U2_nb[15] = 0.5*(U2_nb[1] + U2_nb[4]);
	    }
	    else
	    {
		U0_nb[15] = cell_center[index_nb[15]].m_state.m_U[0];
		U1_nb[15] = cell_center[index_nb[15]].m_state.m_U[1];
		U2_nb[15] = cell_center[index_nb[15]].m_state.m_U[2];
	    }
	    if (flag[1] == 1 && flag[5] == 1) //neighbourcell 16
	    {
		U0_nb[16] = 0.5*(U0_nb[1] + U0_nb[5]);
		U1_nb[16] = 0.5*(U1_nb[1] + U1_nb[5]);
		U2_nb[16] = 0.5*(U2_nb[1] + U2_nb[5]);
	    }
	    else
	    {
		U0_nb[16] = cell_center[index_nb[16]].m_state.m_U[0];
		U1_nb[16] = cell_center[index_nb[16]].m_state.m_U[1];
		U2_nb[16] = cell_center[index_nb[16]].m_state.m_U[2];
	    }
	    if (flag[0] == 1 && flag[5] == 1) //neighbourcell 17
	    {
		U0_nb[17] = 0.5*(U0_nb[0] + U0_nb[5]);
		U1_nb[17] = 0.5*(U1_nb[0] + U1_nb[5]);
		U2_nb[17] = 0.5*(U2_nb[0] + U2_nb[5]);
	    }
	    else
	    {
		U0_nb[17] = cell_center[index_nb[17]].m_state.m_U[0];
		U1_nb[17] = cell_center[index_nb[17]].m_state.m_U[1];
		U2_nb[17] = cell_center[index_nb[17]].m_state.m_U[2];
	    }



            getRectangleCenter(index, coords);
            computeSourceTerm(coords, state);


    //Setting the coeffecients for the first equation


	    coeff[0] = m_dt/rho * mu_edge[0]/(top_h[0]*top_h[0]);
	    coeff[1] = m_dt/rho * mu_edge[1]/(top_h[0]*top_h[0]);
	    coeff[2] = 0.5*m_dt/rho * mu_edge[2]/(top_h[1]*top_h[1]);
	    coeff[3] = 0.5*m_dt/rho * mu_edge[3]/(top_h[1]*top_h[1]);
	    coeff[4] = 0.5*m_dt/rho * mu_edge[4]/(top_h[2]*top_h[2]);
	    coeff[5] = 0.5*m_dt/rho * mu_edge[5]/(top_h[2]*top_h[2]);

	    coeff[6] =  0.5*m_dt/rho * mu[2]/(4*top_h[0]*top_h[1]);
	    coeff[7] = -0.5*m_dt/rho * mu[2]/(4*top_h[0]*top_h[1]);
	    coeff[8] =  0.5*m_dt/rho * mu[3]/(4*top_h[0]*top_h[1]);
	    coeff[9] = -0.5*m_dt/rho * mu[3]/(4*top_h[0]*top_h[1]);

	    coeff[14] =  0.5*m_dt/rho * mu[4]/(4*top_h[0]*top_h[2]);
	    coeff[15] = -0.5*m_dt/rho * mu[4]/(4*top_h[0]*top_h[2]);
	    coeff[16] =  0.5*m_dt/rho * mu[5]/(4*top_h[0]*top_h[2]);
	    coeff[17] = -0.5*m_dt/rho * mu[5]/(4*top_h[0]*top_h[2]);

	    solver.Set_A(I*3,I*3,1.0+coeff[0]+coeff[1]+coeff[2]+coeff[3]
				+coeff[4]+coeff[5]);
	    rhs = (1.0-coeff[0]-coeff[1]-coeff[2]-coeff[3]-coeff[4]-
				coeff[5])*cell_center[index].m_state.m_U[0];

	    for (nb = 0; nb < 6; nb++)
	    {
		if (I_nb[nb] != -1)
		{
		    solver.Set_A(I*3,I_nb[nb]*3,-coeff[nb]);
		    rhs += coeff[nb]*U0_nb[nb];
		}
		else
		    rhs += 2.0*coeff[nb]*U0_nb[nb];
	    }
	    for (nb = 6; nb < 10; nb++)
	    {
		if (I_nb[nb] != -1)
		{
		    solver.Set_A(I*3,I_nb[nb]*3+1,-coeff[nb]);
		    rhs += coeff[nb]*U1_nb[nb];
		}
		else
		    rhs += 2.0*coeff[nb]*U1_nb[nb];
	    }
	    for (nb = 14; nb < 18; nb++)
	    {
		if (I_nb[nb] != -1)
		{
		    solver.Set_A(I*3,I_nb[nb]*3+2,-coeff[nb]);
		    rhs += coeff[nb]*U2_nb[nb];
		}
		else
		    rhs += 2.0*coeff[nb]*U2_nb[nb];
	    }

	    rhs += m_dt*state.m_U[0];
	    rhs += m_dt*cell_center[index].m_state.f_surf[0];
	    rhs -= m_dt*cell_center[index].m_state.grad_q[0]/rho;

	    solver.Set_b(I*3, rhs);

    //Setting the coeffecients for the second equation

	    coeff[0] = 0.5*m_dt/rho * mu_edge[0]/(top_h[0]*top_h[0]);
	    coeff[1] = 0.5*m_dt/rho * mu_edge[1]/(top_h[0]*top_h[0]);
	    coeff[2] = m_dt/rho * mu_edge[2]/(top_h[1]*top_h[1]);
	    coeff[3] = m_dt/rho * mu_edge[3]/(top_h[1]*top_h[1]);
	    coeff[4] = 0.5*m_dt/rho * mu_edge[4]/(top_h[2]*top_h[2]);
	    coeff[5] = 0.5*m_dt/rho * mu_edge[5]/(top_h[2]*top_h[2]);

	    coeff[6] =  0.5*m_dt/rho * mu[0]/(4*top_h[0]*top_h[1]);
	    coeff[7] = -0.5*m_dt/rho * mu[1]/(4*top_h[0]*top_h[1]);
	    coeff[8] =  0.5*m_dt/rho * mu[1]/(4*top_h[0]*top_h[1]);
	    coeff[9] = -0.5*m_dt/rho * mu[0]/(4*top_h[0]*top_h[1]);

	    coeff[10] =  0.5*m_dt/rho * mu[4]/(4*top_h[1]*top_h[2]);
	    coeff[11] = -0.5*m_dt/rho * mu[4]/(4*top_h[1]*top_h[2]);
	    coeff[12] =  0.5*m_dt/rho * mu[5]/(4*top_h[1]*top_h[2]);
	    coeff[13] = -0.5*m_dt/rho * mu[5]/(4*top_h[1]*top_h[2]);

	    solver.Set_A(I*3+1,I*3+1,1.0+coeff[0]+coeff[1]+coeff[2]+coeff[3]
				+coeff[4]+coeff[5]);
	    rhs = (1.0-coeff[0]-coeff[1]-coeff[2]-coeff[3]-coeff[4]-coeff[5])
				*cell_center[index].m_state.m_U[1];

	    for (nb = 0; nb < 6; nb++)
	    {
		if (I_nb[nb] != -1)
		{
		    solver.Set_A(I*3+1,I_nb[nb]*3+1,-coeff[nb]);
		    rhs += coeff[nb]*U1_nb[nb];
		}
		else
		    rhs += 2.0*coeff[nb]*U1_nb[nb];
	    }
	    for (nb = 6; nb < 10; nb++)
	    {
		if (I_nb[nb] != -1)
		{
		    solver.Set_A(I*3+1,I_nb[nb]*3,-coeff[nb]);
		    rhs += coeff[nb]*U0_nb[nb];
		}
		else
		    rhs += 2.0*coeff[nb]*U0_nb[nb];
	    }
	    for (nb = 10; nb < 14; nb++)
	    {
		if (I_nb[nb] != -1)
		{
		    solver.Set_A(I*3+1,I_nb[nb]*3+2,-coeff[nb]);
		    rhs += coeff[nb]*U2_nb[nb];
		}
		else
		    rhs += 2.0*coeff[nb]*U2_nb[nb];
	    }

	    rhs += m_dt*state.m_U[1];
	    rhs += m_dt*cell_center[index].m_state.f_surf[1];
	    rhs -= m_dt*cell_center[index].m_state.grad_q[1]/rho;

	    solver.Set_b(I*3+1, rhs);

    //Setting the coeffecients for the third equation

	    coeff[0] = 0.5*m_dt/rho * mu_edge[0]/(top_h[0]*top_h[0]);
	    coeff[1] = 0.5*m_dt/rho * mu_edge[1]/(top_h[0]*top_h[0]);
	    coeff[2] = 0.5*m_dt/rho * mu_edge[2]/(top_h[1]*top_h[1]);
	    coeff[3] = 0.5*m_dt/rho * mu_edge[3]/(top_h[1]*top_h[1]);
	    coeff[4] = m_dt/rho * mu_edge[4]/(top_h[2]*top_h[2]);
	    coeff[5] = m_dt/rho * mu_edge[5]/(top_h[2]*top_h[2]);

	    coeff[14] =  0.5*m_dt/rho * mu[0]/(4*top_h[0]*top_h[2]);
	    coeff[15] = -0.5*m_dt/rho * mu[1]/(4*top_h[0]*top_h[2]);
	    coeff[16] =  0.5*m_dt/rho * mu[1]/(4*top_h[0]*top_h[2]);
	    coeff[17] = -0.5*m_dt/rho * mu[0]/(4*top_h[0]*top_h[2]);

	    coeff[10] =  0.5*m_dt/rho * mu[2]/(4*top_h[1]*top_h[2]);
	    coeff[11] = -0.5*m_dt/rho * mu[3]/(4*top_h[1]*top_h[2]);
	    coeff[12] =  0.5*m_dt/rho * mu[3]/(4*top_h[1]*top_h[2]);
	    coeff[13] = -0.5*m_dt/rho * mu[2]/(4*top_h[1]*top_h[2]);

	    solver.Set_A(I*3+2,I*3+2,1.0+coeff[0]+coeff[1]+coeff[2]+coeff[3]
				+coeff[4]+coeff[5]);
	    rhs = (1.0-coeff[0]-coeff[1]-coeff[2]-coeff[3]-coeff[4]
				-coeff[5])*cell_center[index].m_state.m_U[2];

	    for (nb = 0; nb < 6; nb++)
	    {
		if (I_nb[nb] != -1)
		{
		    solver.Set_A(I*3+2,I_nb[nb]*3+2,-coeff[nb]);
		    rhs += coeff[nb]*U2_nb[nb];
		}
		else
		    rhs += 2.0*coeff[nb]*U2_nb[nb];
	    }
	    for (nb = 14; nb < 18; nb++)
	    {
		if (I_nb[nb] != -1)
		{
		    solver.Set_A(I*3+2,I_nb[nb]*3,-coeff[nb]);
		    rhs += coeff[nb]*U0_nb[nb];
		}
		else
		    rhs += 2.0*coeff[nb]*U0_nb[nb];
	    }
	    for (nb = 10; nb < 14; nb++)
	    {
		if (I_nb[nb] != -1)
		{
		    solver.Set_A(I*3+2,I_nb[nb]*3+1,-coeff[nb]);
		    rhs += coeff[nb]*U1_nb[nb];
		}
		else
		    rhs += 2.0*coeff[nb]*U1_nb[nb];
	    }

	    rhs += m_dt*state.m_U[2];
	    rhs += m_dt*cell_center[index].m_state.f_surf[2];
	    rhs -= m_dt*cell_center[index].m_state.grad_q[2]/rho;

	    solver.Set_b(I*3+2, rhs);
        }

        solver.SetMaxIter(40000);
        solver.SetTol(1e-14);

	start_clock("Before Petsc Solve");
        solver.Solve_GMRES();
	solver.GetNumIterations(&num_iter);
	solver.GetFinalRelativeResidualNorm(&rel_residual);

	stop_clock("After Petsc Solve");

        // get back the solution
        solver.Get_x(x);

        if (debugging("PETSc"))
            (void) printf("Incompress_Solver_Smooth_3D_Cartesian::"
			  "compDiffWithSmoothProperty: "
                          "num_iter = %d, rel_residual = %g. \n",
                          num_iter,rel_residual);


	for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I = ijk_to_I[i][j][k];
            index = d_index3d(i,j,k,top_gmax);
            if (I >= 0)
            {
                cell_center[index].m_state.m_U[0] = x[I*3-ilower*3];
                cell_center[index].m_state.m_U[1] = x[I*3-ilower*3+1];
		cell_center[index].m_state.m_U[2] = x[I*3-ilower*3+2];
                speed = fabs(cell_center[index].m_state.m_U[0]) +
                        fabs(cell_center[index].m_state.m_U[1]) +
			fabs(cell_center[index].m_state.m_U[2]);
                if (speed > max_speed)
                    max_speed = speed;
            }
            else
            {
                cell_center[index].m_state.m_U[0] = 0.0;
                cell_center[index].m_state.m_U[1] = 0.0;
		cell_center[index].m_state.m_U[2] = 0.0;
            }
        }
        for (l = 0; l < 3; ++l)
        {
	    for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)

            {
                index  = d_index3d(i,j,k,top_gmax);
                array[index] = cell_center[index].m_state.m_U[l];
            }
            scatMeshArray();
	    for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)

            {
                index  = d_index3d(i,j,k,top_gmax);
                cell_center[index].m_state.m_U[l] = array[index];
            }
        }
        pp_global_max(&max_speed,1);

        FT_FreeThese(1,x);
}       /* end compDiffWithSmoothProperty3d */

double Incompress_Solver_Smooth_3D_Cartesian::getVorticityX(int i, int j, int k)
{
	int index0,index00,index01,index10,index11;
	double v00,v01,v10,v11;
	double dy,dz;
	double vorticity;

	dy = top_h[1];
	dz = top_h[2];
	index0 = d_index3d(i,j,k,top_gmax);
	index00 = d_index3d(i,j-1,k,top_gmax);
	index01 = d_index3d(i,j+1,k,top_gmax);
	index10 = d_index3d(i,j,k-1,top_gmax);
	index11 = d_index3d(i,j,k+1,top_gmax);
	v00 = -cell_center[index00].m_state.m_U[2];
	v01 =  cell_center[index01].m_state.m_U[2];
	v10 =  cell_center[index10].m_state.m_U[1];
	v11 = -cell_center[index11].m_state.m_U[1];

	vorticity = (v00 + v01)/2.0/dz + (v10 + v11)/2.0/dy;
	return vorticity;
}	/* end getVorticityX */

double Incompress_Solver_Smooth_3D_Cartesian::getVorticityY(int i, int j, int k)
{
	int index0,index00,index01,index10,index11;
	double v00,v01,v10,v11;
	double dx,dz;
	double vorticity;

	dx = top_h[0];
	dz = top_h[2];
	index0 = d_index3d(i,j,k,top_gmax);
	index00 = d_index3d(i,j,k-1,top_gmax);
	index01 = d_index3d(i,j,k+1,top_gmax);
	index10 = d_index3d(i-1,j,k,top_gmax);
	index11 = d_index3d(i+1,j,k,top_gmax);
	v00 = -cell_center[index00].m_state.m_U[0];
	v01 =  cell_center[index01].m_state.m_U[0];
	v10 =  cell_center[index10].m_state.m_U[2];
	v11 = -cell_center[index11].m_state.m_U[2];

	vorticity = (v00 + v01)/2.0/dx + (v10 + v11)/2.0/dz;
	return vorticity;
}	/* end getVorticityY */

double Incompress_Solver_Smooth_3D_Cartesian::getVorticityZ(int i, int j, int k)
{
	int index0,index00,index01,index10,index11;
	double v00,v01,v10,v11;
	double dx,dy;
	double vorticity;

	dx = top_h[0];
	dy = top_h[1];
	index0 = d_index3d(i,j,k,top_gmax);
	index00 = d_index3d(i-1,j,k,top_gmax);
	index01 = d_index3d(i+1,j,k,top_gmax);
	index10 = d_index3d(i,j-1,k,top_gmax);
	index11 = d_index3d(i,j+1,k,top_gmax);
	v00 = -cell_center[index00].m_state.m_U[1];
	v01 =  cell_center[index01].m_state.m_U[1];
	v10 =  cell_center[index10].m_state.m_U[0];
	v11 = -cell_center[index11].m_state.m_U[0];

	vorticity = (v00 + v01)/2.0/dy + (v10 + v11)/2.0/dx;
	return vorticity;
}	/* end getVorticityZ */


void Incompress_Solver_Smooth_3D_Cartesian::
	compDiffWithSmoothProperty_2nd_decoupled(void)
{
    COMPONENT comp;
    int index,index_nb[6],size;
    int I,I_nb[6];
    int i,j,k,l,nb,icoords[MAXD];
    L_STATE source_term;
    double coords[MAXD], crx_coords[MAXD];
    double mu[6],mu0,rho,rhs;

    // U_nb contains states at neighbor cell or states on the boundary.
    double U_nb[6], U_nb_new[6], U_center;

    double speed;
    double *x;
    GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
    POINTER intfc_state;
    HYPER_SURF *hs;

    setIndexMap();


    max_speed = 0.0;

    size = iupper - ilower;
    FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));

    for (l = 0; l < dim; ++l)
    {
	PETSc solver;

	solver.Create(ilower, iupper-1, 7, 7);
	solver.Reset_A();
	solver.Reset_b();
	solver.Reset_x();

	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
	for (i = imin; i <= imax; i++)
	{
		    I  = ijk_to_I[i][j][k];
		    if (I == -1) continue;

		    index  = d_index3d(i,j,k,top_gmax);
		    index_nb[0] = d_index3d(i-1,j,k,top_gmax);
		    index_nb[1] = d_index3d(i+1,j,k,top_gmax);
		    index_nb[2] = d_index3d(i,j-1,k,top_gmax);
		    index_nb[3] = d_index3d(i,j+1,k,top_gmax);
		    index_nb[4] = d_index3d(i,j,k-1,top_gmax);
		    index_nb[5] = d_index3d(i,j,k+1,top_gmax);

		    icoords[0] = i;
		    icoords[1] = j;
		    icoords[2] = k;
		    comp = top_comp[index];

		    I_nb[0] = ijk_to_I[i-1][j][k]; //west
		    I_nb[1] = ijk_to_I[i+1][j][k]; //east
		    I_nb[2] = ijk_to_I[i][j-1][k]; //south
		    I_nb[3] = ijk_to_I[i][j+1][k]; //north
		    I_nb[4] = ijk_to_I[i][j][k-1]; //lower
		    I_nb[5] = ijk_to_I[i][j][k+1]; //upper


		    mu0   = cell_center[index].m_state.m_mu;
		    rho   = cell_center[index].m_state.m_rho;
		    U_center =  cell_center[index].m_state.m_U[l];

		    for (nb = 0; nb < 6; nb++)
		    {
			if (FT_StateStructAtGridCrossing(front,icoords,dir[nb],
				comp,&intfc_state,&hs,crx_coords) &&
				wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
			{
			    // old boundary condition
			    U_nb[nb] = getStateVel[l](intfc_state);
			    // new boundary condition
			    FT_StateStructAtGridCrossing(front,icoords,dir[nb],
				    comp,&intfc_state,&hs,crx_coords);
			    U_nb_new[nb] = getStateVel[l](intfc_state);

			    if (wave_type(hs) == DIRICHLET_BOUNDARY ||
				    wave_type(hs) == NEUMANN_BOUNDARY)
				mu[nb] = mu0;
			    else
				mu[nb] = 1.0/2*(mu0 +
					cell_center[index_nb[nb]].m_state.m_mu);
			}
			else
			{
			    U_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[l];
			    mu[nb] = 1.0/2*(mu0 +
				    cell_center[index_nb[nb]].m_state.m_mu);
			}
		    }


		    getRectangleCenter(index, coords);
		    computeSourceTerm(coords, source_term);


		    rhs = 0;
		    double dh[3] = {top_h[0], top_h[1], top_h[2]};
		    for(nb = 0; nb<6; nb++)
		    {
			// use dh[1]*dh[2] as the face area
			if(nb<2)
			{
			    dh[0] = top_h[0];
			    dh[1] = top_h[1];
			    dh[2] = top_h[2];
			}
			else if(nb<4 && nb>=2)
			{
			    dh[0] = top_h[1];
			    dh[1] = top_h[2];
			    dh[2] = top_h[0];
			}
			else if(nb<6 && nb>=4)
			{
			    dh[0] = top_h[2];
			    dh[1] = top_h[0];
			    dh[2] = top_h[1];
			}

			if(I_nb[nb]>=0)
			{
			    // u^{*}
			    solver.Add_A(I, I,(-1) * -0.5*m_dt/rho*mu[nb] /(dh[0]*dh[0]));
			    solver.Add_A(I, I_nb[nb], (-1) *  0.5*m_dt/rho*mu[nb] /(dh[0]*dh[0]));
			    // u^{n}
			    rhs += 0.5*m_dt/rho*mu[nb] *
				    (U_nb[nb]-U_center) /(dh[0]*dh[0]);
			}
			else		// boundary
			{
			    // u^{*}
			    solver.Add_A(I, I,       (-1) * -0.5*m_dt/rho*mu[nb] /(dh[0]*dh[0]) * 2);
			    rhs += 0.5*m_dt/rho*mu[nb] /(dh[0]*dh[0]) * 2 * U_nb_new[nb];
			    // u^{n}
			    rhs += 0.5*m_dt/rho*mu[nb] *
				    (U_nb[nb]-U_center) /(dh[0]*dh[0]) * 2;
			}
		    }

		    // interior
		    rhs += m_dt*source_term.m_U[l];
		    rhs += m_dt*cell_center[index].m_state.f_surf[l];
		    rhs -= m_dt*cell_center[index].m_state.grad_q[l]/rho;
		    rhs -= m_dt * cell_center[index].m_state.m_adv[l]; 
						//advection source term

		    solver.Add_A(I, I, 1.0);
		    rhs += U_center;

		    solver.Add_b(I, rhs);
	}

	solver.SetMaxIter(40000);
	solver.SetTol(1e-15);
	solver.Solve_GMRES();

	// get back the solution
	solver.Get_x(x);

	PetscInt num_iter;
	double rel_residual;
	solver.GetNumIterations(&num_iter);
	solver.GetFinalRelativeResidualNorm(&rel_residual);

	if (debugging("PETSc"))
	    (void) printf("L_CARTESIAN::"
		    "compDiffWithSmoothProperty_2nd_decoupled: "
		    "num_iter = %d, rel_residual = %g. \n",
		    num_iter,rel_residual);

	for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
		for (i = imin; i <= imax; i++)
		{
		    I = ijk_to_I[i][j][k];
		    index = d_index3d(i,j,k,top_gmax);
		    if (I >= 0)
		    {
			cell_center[index].m_state.m_U[l] = x[I-ilower];
		    }
		    else
		    {
			cell_center[index].m_state.m_U[l] = 0.0;
		    }
		}

	for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
		for (i = imin; i <= imax; i++)
		{
		    index  = d_index3d(i,j,k,top_gmax);
		    array[index] = cell_center[index].m_state.m_U[l];
		}
	scatMeshArray();
	for (k = 0; k <= top_gmax[2]; k++)
	    for (j = 0; j <= top_gmax[1]; j++)
		for (i = 0; i <= top_gmax[0]; i++)
		{
		    index  = d_index3d(i,j,k,top_gmax);
		    cell_center[index].m_state.m_U[l] = array[index];
		}
    }
    for (k = kmin; k <= kmax; k++)
    for (j = jmin; j <= jmax; j++)
    for (i = imin; i <= imax; i++)
    {
	index = d_index3d(i,j,k,top_gmax);
	speed = fabs(cell_center[index].m_state.m_U[0]) +
		fabs(cell_center[index].m_state.m_U[1]) +
		fabs(cell_center[index].m_state.m_U[2]);
	if (speed > max_speed)
	    max_speed = speed;
    }
    pp_global_max(&max_speed,1);
    free_these(1,x);

}

void Incompress_Solver_Smooth_3D_Cartesian::
		compDiffWithSmoothProperty_2nd_coupled(void)
{
        COMPONENT comp;
        int index,index_nb[18],size;
        int I,I_nb[18];
        double coords[MAXD],crx_coords[MAXD];
	double coeff0[6],coeff1[6],coeff2[6],coeff_temp;
	double mu[6],mu_edge[6],mu0,rho,rhs;
	double U0_nb[18],U1_nb[18],U2_nb[18],U0_center,U1_center,U2_center;
	int flag[6]; //denote whether this is dirichlet or neumann boundary
        L_STATE state;
        int i,j,k,l,nb,icoords[MAXD];
        double speed;
	double *x;
	GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
	POINTER intfc_state;
	HYPER_SURF *hs;
	PetscInt num_iter;
	double rel_residual;

	max_speed = 0.0;
	setIndexMap();

	size = iupper - ilower;
	FT_VectorMemoryAlloc((POINTER*)&x,3*size,sizeof(double));

        PETSc solver;
        solver.Create(3*ilower, 3*iupper-1, 25, 25);
	// 7u + 9v + 9w for the first equation
	// 7v + 9u + 9w for the second equation
	// 7w + 9u + 9v for the third equation

	solver.Reset_A();
	solver.Reset_b();
	solver.Reset_x();

	for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I  = ijk_to_I[i][j][k];
            if (I == -1) continue;

            index  = d_index3d(i,j,k,top_gmax);	
	//6 neighbours of the center cell
            index_nb[0] = d_index3d(i-1,j,k,top_gmax);
            index_nb[1] = d_index3d(i+1,j,k,top_gmax);
            index_nb[2] = d_index3d(i,j-1,k,top_gmax);
            index_nb[3] = d_index3d(i,j+1,k,top_gmax);
	    index_nb[4] = d_index3d(i,j,k-1,top_gmax);
	    index_nb[5] = d_index3d(i,j,k+1,top_gmax);

	//xy cut neighbours
	    index_nb[6] = d_index3d(i-1,j-1,k,top_gmax);
	    index_nb[7] = d_index3d(i+1,j-1,k,top_gmax);
	    index_nb[8] = d_index3d(i+1,j+1,k,top_gmax);
	    index_nb[9] = d_index3d(i-1,j+1,k,top_gmax);
	//yz cut neighbours
	    index_nb[10] = d_index3d(i,j-1,k-1,top_gmax);
	    index_nb[11] = d_index3d(i,j+1,k-1,top_gmax);
	    index_nb[12] = d_index3d(i,j+1,k+1,top_gmax);
	    index_nb[13] = d_index3d(i,j-1,k+1,top_gmax);
	//xz cut neighbours
	    index_nb[14] = d_index3d(i-1,j,k-1,top_gmax);
	    index_nb[15] = d_index3d(i+1,j,k-1,top_gmax);
	    index_nb[16] = d_index3d(i+1,j,k+1,top_gmax);
	    index_nb[17] = d_index3d(i-1,j,k+1,top_gmax);

        
	//6 neighbours of the center cell
            I_nb[0] = ijk_to_I[i-1][j][k];
            I_nb[1] = ijk_to_I[i+1][j][k];
            I_nb[2] = ijk_to_I[i][j-1][k];
            I_nb[3] = ijk_to_I[i][j+1][k];
            I_nb[4] = ijk_to_I[i][j][k-1];
            I_nb[5] = ijk_to_I[i][j][k+1];
	
	//xy cut neighbours
            I_nb[6] = ijk_to_I[i-1][j-1][k];
	    I_nb[7] = ijk_to_I[i+1][j-1][k];
	    I_nb[8] = ijk_to_I[i+1][j+1][k];
	    I_nb[9] = ijk_to_I[i-1][j+1][k];
	//yz cut neighbours
	    I_nb[10] = ijk_to_I[i][j-1][k-1];
	    I_nb[11] = ijk_to_I[i][j+1][k-1];
	    I_nb[12] = ijk_to_I[i][j+1][k+1];
	    I_nb[13] = ijk_to_I[i][j-1][k+1];
	//xz cut neighbours
	    I_nb[14] = ijk_to_I[i-1][j][k-1];
	    I_nb[15] = ijk_to_I[i+1][j][k-1];
	    I_nb[16] = ijk_to_I[i+1][j][k+1];
	    I_nb[17] = ijk_to_I[i-1][j][k+1];

	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;
	    comp = top_comp[index];

	    mu0 = cell_center[index].m_state.m_mu;
	    rho = cell_center[index].m_state.m_rho;
	    U0_center = cell_center[index].m_state.m_U[0];
	    U1_center = cell_center[index].m_state.m_U[1];
	    U2_center = cell_center[index].m_state.m_U[2];

            for (nb = 0; nb < 6; nb++)
            {
                if (FT_StateStructAtGridCrossing(front,icoords,dir[nb],
                            comp,&intfc_state,&hs,crx_coords) &&
                            wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	       	{
		    flag[nb] = 1;
		    U0_nb[nb] = getStateVel[0](intfc_state);
		    U1_nb[nb] = getStateVel[1](intfc_state);
		    U2_nb[nb] = getStateVel[2](intfc_state);

		    if (wave_type(hs) == DIRICHLET_BOUNDARY || 
			wave_type(hs) == NEUMANN_BOUNDARY)
		    {
			mu[nb] = mu0;
			mu_edge[nb] = mu0;
		    }
		    else
		    {
			mu[nb] = 0.5*(mu0 + 
				cell_center[index_nb[nb]].m_state.m_mu);
			mu_edge[nb] = cell_center[index_nb[nb]].m_state.m_mu;
		    }
		}
		else
		{
		    flag[nb] = 0;
                    U0_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[0];
		    U1_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[1];
		    U2_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[2];
		    mu_edge[nb] = 0.5*(mu0 + 
				cell_center[index_nb[nb]].m_state.m_mu);
		    mu[nb] = cell_center[index_nb[nb]].m_state.m_mu;
		}
	    }

	    for (nb = 6; nb < 18; nb++) //cut corner values, interior
	    {
		if (I_nb[nb] != -1)
		{
		    U0_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[0];
		    U1_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[1];
		    U2_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[2];
		}
	    }


	    // source term
            getRectangleCenter(index, coords);
            computeSourceTerm(coords, state);

	    //Setting the coeffecients for U0 in the first equation

	    coeff0[0] = m_dt/rho * mu_edge[0]/(top_h[0]*top_h[0]);
	    coeff0[1] = m_dt/rho * mu_edge[1]/(top_h[0]*top_h[0]);
	    coeff0[2] = 0.5*m_dt/rho * mu_edge[2]/(top_h[1]*top_h[1]);
	    coeff0[3] = 0.5*m_dt/rho * mu_edge[3]/(top_h[1]*top_h[1]);
	    coeff0[4] = 0.5*m_dt/rho * mu_edge[4]/(top_h[2]*top_h[2]);
	    coeff0[5] = 0.5*m_dt/rho * mu_edge[5]/(top_h[2]*top_h[2]);

	    solver.Set_A(I*3,I*3,1.0);
	    rhs = U0_center;

	    for (nb = 0; nb < 6; nb++)
	    {
		if (I_nb[nb] != -1)
		{
		    solver.Set_A(I*3,I_nb[nb]*3,-coeff0[nb]);
		    rhs += coeff0[nb]*U0_nb[nb];
		}
		else
		{
		    coeff0[nb] = 2.0*coeff0[nb];
		    rhs += 2.0*coeff0[nb]*U0_nb[nb];
		}
		solver.Add_A(I*3,I*3,coeff0[nb]);
		rhs -= coeff0[nb]*U0_center;
	    }

	    //set the coefficients for U1 in the first equation (mu*v_x)_y
	    //traverse the four corners
	    //(i-1/2,j-1/2,k),(i+1/2,j-1/2,k),(i+1/2,j+1/2,k),(i-1/2,j+1/2,k)
	    
	    //corner (i-1/2,j-1/2,k)

	    coeff_temp = mu_edge[2]*m_dt/(top_h[0]*top_h[1])/rho;

	    if (flag[0] == 1 && flag[2] == 0)
		rhs += coeff_temp*U1_nb[0];
	    else if(flag[0] == 0 && flag[2] == 1)
		rhs += coeff_temp*U1_nb[2];
	    else if(flag[0] == 1 && flag[2] == 1)
		rhs += coeff_temp*U1_nb[0];
	    else {
		rhs += coeff_temp*(U1_nb[0]+U1_nb[2]+U1_nb[6]+U1_center)/8.0;

		solver.Add_A(I*3,I*3+1,       -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[0]*3+1, -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[2]*3+1, -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[6]*3+1, -coeff_temp/8.0);
	    }
		    
	    //corner (i+1/2,j-1/2,k)

	    coeff_temp = -mu_edge[2]*m_dt/(top_h[0]*top_h[1])/rho;

	    if (flag[2] == 1 && flag[1] == 0)
		rhs += coeff_temp*U1_nb[2];
	    else if(flag[2] == 0 && flag[1] == 1)
		rhs += coeff_temp*U1_nb[1];
	    else if(flag[2] == 1 && flag[1] == 1)
		rhs += coeff_temp*U1_nb[2];
	    else {
		rhs += coeff_temp*(U1_nb[2]+U1_nb[1]+U1_nb[7]+U1_center)/8.0;

		solver.Add_A(I*3,I*3+1,       -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[2]*3+1, -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[1]*3+1, -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[7]*3+1, -coeff_temp/8.0);
	    }

	    
	    //corner (i+1/2,j+1/2,k)

	    coeff_temp = mu_edge[3]*m_dt/(top_h[0]*top_h[1])/rho;

	    if (flag[1] == 1 && flag[3] == 0)
		rhs += coeff_temp*U1_nb[1];
	    else if(flag[1] == 0 && flag[3] == 1)
		rhs += coeff_temp*U1_nb[3];
	    else if(flag[1] == 1 && flag[3] == 1)
		rhs += coeff_temp*U1_nb[1];
	    else {
		rhs += coeff_temp*(U1_nb[1]+U1_nb[3]+U1_nb[8]+U1_center)/8.0;

		solver.Add_A(I*3,I*3+1,       -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[1]*3+1, -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[3]*3+1, -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[8]*3+1, -coeff_temp/8.0);
	    }

	    //corner (i-1/2,j+1/2,k)

	    coeff_temp = -mu_edge[3]*m_dt/(top_h[0]*top_h[1])/rho;

	    if (flag[3] == 1 && flag[0] == 0)
		rhs += coeff_temp*U1_nb[3];
	    else if(flag[3] == 0 && flag[0] == 1)
		rhs += coeff_temp*U1_nb[0];
	    else if(flag[3] == 1 && flag[0] == 1)
		rhs += coeff_temp*U1_nb[3];
	    else {
		rhs += coeff_temp*(U1_nb[3]+U1_nb[0]+U1_nb[9]+U1_center)/8.0;

		solver.Add_A(I*3,I*3+1,       -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[3]*3+1, -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[0]*3+1, -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[9]*3+1, -coeff_temp/8.0);
	    }

	    //set the coefficients for U2 in the first equation (mu*w_x)_z
	    //traverse the four corners
	    //(i-1/2,j,k-1/2),(i+1/2,j,k-1/2),(i+1/2,j,k+1/2),(i-1/2,j,k+1/2)
	    
	    //corner (i-1/2,j,k-1/2)

	    coeff_temp = mu_edge[4]*m_dt/(top_h[0]*top_h[2])/rho;

	    if (flag[0] == 1 && flag[4] == 0)
		rhs += coeff_temp*U2_nb[0];
	    else if(flag[0] == 0 && flag[4] == 1)
		rhs += coeff_temp*U2_nb[4];
	    else if(flag[0] == 1 && flag[4] == 1)
		rhs += coeff_temp*U2_nb[0];
	    else {
		rhs += coeff_temp*(U2_nb[0]+U2_nb[4]+U2_nb[14]+U2_center)/8.0;

		solver.Add_A(I*3,I*3+2,       -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[0]*3+2, -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[4]*3+2, -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[14]*3+2, -coeff_temp/8.0);
	    }
		    
	    //corner (i+1/2,j,k-1/2)

	    coeff_temp = -mu_edge[4]*m_dt/(top_h[0]*top_h[2])/rho;

	    if (flag[1] == 1 && flag[4] == 0)
		rhs += coeff_temp*U2_nb[1];
	    else if(flag[1] == 0 && flag[4] == 1)
		rhs += coeff_temp*U2_nb[4];
	    else if(flag[1] == 1 && flag[4] == 1)
		rhs += coeff_temp*U2_nb[1];
	    else {
		rhs += coeff_temp*(U2_nb[1]+U2_nb[4]+U2_nb[15]+U2_center)/8.0;

		solver.Add_A(I*3,I*3+2,       -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[1]*3+2, -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[4]*3+2, -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[15]*3+2, -coeff_temp/8.0);
	    }

	    
	    //corner (i+1/2,j,k+1/2)

	    coeff_temp = mu_edge[5]*m_dt/(top_h[0]*top_h[2])/rho;

	    if (flag[1] == 1 && flag[5] == 0)
		rhs += coeff_temp*U2_nb[1];
	    else if(flag[1] == 0 && flag[5] == 1)
		rhs += coeff_temp*U2_nb[5];
	    else if(flag[1] == 1 && flag[5] == 1)
		rhs += coeff_temp*U2_nb[1];
	    else {
		rhs += coeff_temp*(U2_nb[1]+U2_nb[5]+U2_nb[16]+U2_center)/8.0;

		solver.Add_A(I*3,I*3+2,       -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[1]*3+2, -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[5]*3+2, -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[16]*3+2, -coeff_temp/8.0);
	    }

	    //corner (i-1/2,j,k+1/2)

	    coeff_temp = -mu_edge[5]*m_dt/(top_h[0]*top_h[2])/rho;

	    if (flag[5] == 1 && flag[0] == 0)
		rhs += coeff_temp*U2_nb[5];
	    else if(flag[5] == 0 && flag[0] == 1)
		rhs += coeff_temp*U2_nb[0];
	    else if(flag[5] == 1 && flag[0] == 1)
		rhs += coeff_temp*U2_nb[5];
	    else {
		rhs += coeff_temp*(U2_nb[5]+U2_nb[0]+U2_nb[17]+U2_center)/8.0;

		solver.Add_A(I*3,I*3+2,       -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[5]*3+2, -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[0]*3+2, -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[17]*3+2, -coeff_temp/8.0);
	    }

	    rhs += m_dt*state.m_U[0];
	    rhs += m_dt*cell_center[index].m_state.f_surf[0];
	    rhs -= m_dt*cell_center[index].m_state.grad_q[0]/rho;
	    rhs -= m_dt*cell_center[index].m_state.m_adv[0];

	    solver.Set_b(I*3, rhs);

	    //Setting the coeffecients for U1 in the second equation

	    coeff1[0] = 0.5*m_dt/rho * mu_edge[0]/(top_h[0]*top_h[0]);
	    coeff1[1] = 0.5*m_dt/rho * mu_edge[1]/(top_h[0]*top_h[0]);
	    coeff1[2] = m_dt/rho * mu_edge[2]/(top_h[1]*top_h[1]);
	    coeff1[3] = m_dt/rho * mu_edge[3]/(top_h[1]*top_h[1]);
	    coeff1[4] = 0.5*m_dt/rho * mu_edge[4]/(top_h[2]*top_h[2]);
	    coeff1[5] = 0.5*m_dt/rho * mu_edge[5]/(top_h[2]*top_h[2]);

	    solver.Set_A(I*3+1,I*3+1,1.0);
	    rhs = U1_center;

	    for (nb = 0; nb < 6; nb++)
	    {
		if (I_nb[nb] != -1)
		{
		    solver.Set_A(I*3+1,I_nb[nb]*3+1,-coeff1[nb]);
		    rhs += coeff1[nb]*U1_nb[nb];
		}
		else
		{
		    coeff1[nb] = 2.0*coeff1[nb];
		    rhs += 2.0*coeff1[nb]*U1_nb[nb];
		}
		solver.Add_A(I*3+1,I*3+1,coeff1[nb]);
		rhs -= coeff1[nb]*U1_center;
	    }

	    //set the coefficients for U0 in the second equation (mu*u_y)_x
	    //traverse the four corners
	    //(i-1/2,j-1/2,k),(i+1/2,j-1/2,k),(i+1/2,j+1/2,k),(i-1/2,j+1/2,k)

	    //corner (i-1/2,j-1/2,k)
	    
	    coeff_temp = mu_edge[0]*m_dt/(top_h[0]*top_h[1])/rho;

	    if (flag[0] == 1 && flag[2] == 0)
		rhs += coeff_temp*U0_nb[0];
	    else if (flag[0] == 0 && flag[2] == 1)
		rhs += coeff_temp*U0_nb[2];
	    else if (flag[0] == 1 && flag[2] == 1)
		rhs += coeff_temp*U0_nb[0];
	    else {
		rhs += coeff_temp*(U0_nb[0]+U0_nb[2]+U0_nb[6]+U0_center)/8.0;

		solver.Add_A(I*3+1,I*3,        -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[0]*3,  -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[2]*3,  -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[6]*3,  -coeff_temp/8.0);
	    }

	    //corner (i+1/2,j-1/2,k)
	    
	    coeff_temp = -mu_edge[1]*m_dt/(top_h[0]*top_h[1])/rho;

	    if (flag[1] == 1 && flag[2] == 0)
		rhs += coeff_temp*U0_nb[1];
	    else if (flag[1] == 0 && flag[2] == 1)
		rhs += coeff_temp*U0_nb[2];
	    else if (flag[1] == 1 && flag[2] == 1)
		rhs += coeff_temp*U0_nb[1];
	    else {
		rhs += coeff_temp*(U0_nb[1]+U0_nb[2]+U0_nb[7]+U0_center)/8.0;

		solver.Add_A(I*3+1,I*3,        -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[1]*3,  -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[2]*3,  -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[7]*3,  -coeff_temp/8.0);
	    }
	    //corner (i+1/2,j+1/2,k)
	    
	    coeff_temp = mu_edge[1]*m_dt/(top_h[0]*top_h[1])/rho;

	    if (flag[1] == 1 && flag[3] == 0)
		rhs += coeff_temp*U0_nb[1];
	    else if (flag[1] == 0 && flag[3] == 1)
		rhs += coeff_temp*U0_nb[3];
	    else if (flag[1] == 1 && flag[3] == 1)
		rhs += coeff_temp*U0_nb[1];
	    else {
		rhs += coeff_temp*(U0_nb[1]+U0_nb[3]+U0_nb[8]+U0_center)/8.0;

		solver.Add_A(I*3+1,I*3,        -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[1]*3,  -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[3]*3,  -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[8]*3,  -coeff_temp/8.0);
	    }
	    //corner (i-1/2,j+1/2,k)
	    
	    coeff_temp = -mu_edge[0]*m_dt/(top_h[0]*top_h[1])/rho;

	    if (flag[0] == 1 && flag[3] == 0)
		rhs += coeff_temp*U0_nb[0];
	    else if (flag[0] == 0 && flag[3] == 1)
		rhs += coeff_temp*U0_nb[3];
	    else if (flag[0] == 1 && flag[3] == 1)
		rhs += coeff_temp*U0_nb[0];
	    else {
		rhs += coeff_temp*(U0_nb[0]+U0_nb[3]+U0_nb[9]+U0_center)/8.0;

		solver.Add_A(I*3+1,I*3,        -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[0]*3,  -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[3]*3,  -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[9]*3,  -coeff_temp/8.0);
	    }

	    //set the coefficients for U2 in the second equation (mu*w_y)_z
	    //traverse the four corners
	    //(i,j-1/2,k-1/2),(i,j+1/2,k-1/2),(i,j+1/2,k+1/2),(i,j-1/2,k+1/2)

	    //corner (i,j-1/2,k-1/2)
	    
	    coeff_temp = mu_edge[4]*m_dt/(top_h[1]*top_h[2])/rho;

	    if (flag[2] == 1 && flag[4] == 0)
		rhs += coeff_temp*U2_nb[2];
	    else if (flag[2] == 0 && flag[4] == 1)
		rhs += coeff_temp*U2_nb[4];
	    else if (flag[2] == 1 && flag[4] == 1)
		rhs += coeff_temp*U2_nb[2];
	    else {
		rhs += coeff_temp*(U2_nb[2]+U2_nb[4]+U2_nb[10]+U2_center)/8.0;

		solver.Add_A(I*3+1,I*3+2,        -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[2]*3+2,  -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[4]*3+2,  -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[10]*3+2,  -coeff_temp/8.0);
	    }
	
	    //corner (i,j+1/2,k-1/2)
	    
	    coeff_temp = -mu_edge[4]*m_dt/(top_h[1]*top_h[2])/rho;

	    if (flag[3] == 1 && flag[4] == 0)
		rhs += coeff_temp*U2_nb[3];
	    else if (flag[3] == 0 && flag[4] == 1)
		rhs += coeff_temp*U2_nb[4];
	    else if (flag[3] == 1 && flag[4] == 1)
		rhs += coeff_temp*U2_nb[3];
	    else {
		rhs += coeff_temp*(U2_nb[3]+U2_nb[4]+U2_nb[11]+U2_center)/8.0;

		solver.Add_A(I*3+1,I*3+2,        -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[3]*3+2,  -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[4]*3+2,  -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[11]*3+2,  -coeff_temp/8.0);
	    }

	    
	    //corner (i,j+1/2,k+1/2)
	    
	    coeff_temp = mu_edge[5]*m_dt/(top_h[1]*top_h[2])/rho;

	    if (flag[3] == 1 && flag[5] == 0)
		rhs += coeff_temp*U2_nb[3];
	    else if (flag[3] == 0 && flag[5] == 1)
		rhs += coeff_temp*U2_nb[5];
	    else if (flag[3] == 1 && flag[5] == 1)
		rhs += coeff_temp*U2_nb[3];
	    else {
		rhs += coeff_temp*(U2_nb[3]+U2_nb[5]+U2_nb[12]+U2_center)/8.0;

		solver.Add_A(I*3+1,I*3+2,        -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[3]*3+2,  -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[5]*3+2,  -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[12]*3+2,  -coeff_temp/8.0);
	    }

	    //corner (i,j-1/2,k+1/2)
	    
	    coeff_temp = -mu_edge[5]*m_dt/(top_h[1]*top_h[2])/rho;

	    if (flag[2] == 1 && flag[5] == 0)
		rhs += coeff_temp*U2_nb[2];
	    else if (flag[2] == 0 && flag[5] == 1)
		rhs += coeff_temp*U2_nb[5];
	    else if (flag[2] == 1 && flag[5] == 1)
		rhs += coeff_temp*U2_nb[2];
	    else {
		rhs += coeff_temp*(U2_nb[2]+U2_nb[5]+U2_nb[13]+U2_center)/8.0;

		solver.Add_A(I*3+1,I*3+2,        -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[2]*3+2,  -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[5]*3+2,  -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[13]*3+2,  -coeff_temp/8.0);
	    }

	    rhs += m_dt*state.m_U[1];
	    rhs += m_dt*cell_center[index].m_state.f_surf[1];
	    rhs -= m_dt*cell_center[index].m_state.grad_q[1]/rho;
	    rhs -= m_dt*cell_center[index].m_state.m_adv[1];

	    solver.Set_b(I*3+1, rhs);

	    //Setting the coeffecients of U2 for the third equation

	    coeff2[0] = 0.5*m_dt/rho * mu_edge[0]/(top_h[0]*top_h[0]);
	    coeff2[1] = 0.5*m_dt/rho * mu_edge[1]/(top_h[0]*top_h[0]);
	    coeff2[2] = 0.5*m_dt/rho * mu_edge[2]/(top_h[1]*top_h[1]);
	    coeff2[3] = 0.5*m_dt/rho * mu_edge[3]/(top_h[1]*top_h[1]);
	    coeff2[4] = m_dt/rho * mu_edge[4]/(top_h[2]*top_h[2]);
	    coeff2[5] = m_dt/rho * mu_edge[5]/(top_h[2]*top_h[2]);

	    solver.Set_A(I*3+2,I*3+2,1.0);
	    rhs = U2_center;

	    for (nb = 0; nb < 6; nb++)
	    {
		if (I_nb[nb] != -1)
		{
		    solver.Set_A(I*3+2,I_nb[nb]*3+2,-coeff2[nb]);
		    rhs += coeff2[nb]*U2_nb[nb];
		}
		else
		{
		    coeff2[nb] = 2.0*coeff2[nb];
		    rhs += 2.0*coeff2[nb]*U2_nb[nb];
		}
		solver.Add_A(I*3+2,I*3+2,coeff2[nb]);
		rhs -= coeff2[nb]*U2_center;
	    }

	    //set the coefficients for U0 in the thrid equation (mu*u_z)_x
	    //traverse the four corners
	    //(i-1/2,j,k-1/2),(i+1/2,j,k-1/2),(i+1/2,j,k+1/2),(i-1/2,j,k+1/2)

	    //corner (i-1/2,j,k-1/2)

	    coeff_temp = mu_edge[0]*m_dt/(top_h[0]*top_h[2])/rho;

	    if (flag[0] == 1 && flag[4] == 0)
		rhs += coeff_temp*U0_nb[0];
	    else if (flag[0] == 0 && flag[4] == 1)
		rhs += coeff_temp*U0_nb[4];
	    else if (flag[0] == 1 && flag[4] == 1)
		rhs += coeff_temp*U0_nb[0];
	    else {
		rhs += coeff_temp*(U0_nb[0]+U0_nb[4]+U0_nb[14]+U0_center)/8.0;

		solver.Add_A(I*3+2,I*3,        -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[0]*3,  -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[4]*3,  -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[14]*3, -coeff_temp/8.0);
	    }

	    //corner (i+1/2,j,k-1/2)

	    coeff_temp = -mu_edge[1]*m_dt/(top_h[0]*top_h[2])/rho;

	    if (flag[1] == 1 && flag[4] == 0)
		rhs += coeff_temp*U0_nb[1];
	    else if (flag[1] == 0 && flag[4] == 1)
		rhs += coeff_temp*U0_nb[4];
	    else if (flag[1] == 1 && flag[4] == 1)
		rhs += coeff_temp*U0_nb[1];
	    else {
		rhs += coeff_temp*(U0_nb[1]+U0_nb[4]+U0_nb[15]+U0_center)/8.0;

		solver.Add_A(I*3+2,I*3,        -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[1]*3,  -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[4]*3,  -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[15]*3, -coeff_temp/8.0);
	    }

	    //corner (i+1/2,j,k+1/2)

	    coeff_temp = mu_edge[1]*m_dt/(top_h[0]*top_h[2])/rho;

	    if (flag[1] == 1 && flag[5] == 0)
		rhs += coeff_temp*U0_nb[1];
	    else if (flag[1] == 0 && flag[5] == 1)
		rhs += coeff_temp*U0_nb[5];
	    else if (flag[1] == 1 && flag[5] == 1)
		rhs += coeff_temp*U0_nb[1];
	    else {
		rhs += coeff_temp*(U0_nb[1]+U0_nb[5]+U0_nb[16]+U0_center)/8.0;

		solver.Add_A(I*3+2,I*3,        -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[1]*3,  -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[5]*3,  -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[16]*3, -coeff_temp/8.0);
	    }

	    //corner (i-1/2,j,k+1/2)

	    coeff_temp = -mu_edge[0]*m_dt/(top_h[0]*top_h[2])/rho;

	    if (flag[0] == 1 && flag[5] == 0)
		rhs += coeff_temp*U0_nb[0];
	    else if (flag[0] == 0 && flag[5] == 1)
		rhs += coeff_temp*U0_nb[5];
	    else if (flag[0] == 1 && flag[5] == 1)
		rhs += coeff_temp*U0_nb[0];
	    else {
		rhs += coeff_temp*(U0_nb[0]+U0_nb[5]+U0_nb[17]+U0_center)/8.0;

		solver.Add_A(I*3+2,I*3,        -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[0]*3,  -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[5]*3,  -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[17]*3, -coeff_temp/8.0);
	    }

	    //set the coefficients for U1 in the thrid equation (mu*v_z)_y
	    //traverse the four corners
	    //(i,j-1/2,k-1/2),(i,j+1/2,k-1/2),(i,j+1/2,k+1/2),(i,j-1/2,k+1/2)

	    //corner (i,j-1/2,k-1/2)

	    coeff_temp = mu_edge[2]*m_dt/(top_h[1]*top_h[2])/rho;

	    if (flag[2] == 1 && flag[4] == 0)
		rhs += coeff_temp*U1_nb[2];
	    else if (flag[2] == 0 && flag[4] == 1)
		rhs += coeff_temp*U1_nb[4];
	    else if (flag[2] == 1 && flag[4] == 1)
		rhs += coeff_temp*U1_nb[2];
	    else {
		rhs += coeff_temp*(U1_nb[2]+U1_nb[4]+U1_nb[10]+U1_center)/8.0;

		solver.Add_A(I*3+2,I*3+1,        -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[2]*3+1,  -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[4]*3+1,  -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[10]*3+1, -coeff_temp/8.0);
	    }

	    //corner (i,j+1/2,k-1/2)

	    coeff_temp = -mu_edge[3]*m_dt/(top_h[1]*top_h[2])/rho;

	    if (flag[3] == 1 && flag[4] == 0)
		rhs += coeff_temp*U1_nb[3];
	    else if (flag[3] == 0 && flag[4] == 1)
		rhs += coeff_temp*U1_nb[4];
	    else if (flag[3] == 1 && flag[4] == 1)
		rhs += coeff_temp*U1_nb[3];
	    else {
		rhs += coeff_temp*(U1_nb[3]+U1_nb[4]+U1_nb[11]+U1_center)/8.0;

		solver.Add_A(I*3+2,I*3+1,        -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[3]*3+1,  -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[4]*3+1,  -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[11]*3+1, -coeff_temp/8.0);
	    }

	    //corner (i,j+1/2,k+1/2)

	    coeff_temp = mu_edge[3]*m_dt/(top_h[1]*top_h[2])/rho;

	    if (flag[3] == 1 && flag[5] == 0)
		rhs += coeff_temp*U1_nb[3];
	    else if (flag[3] == 0 && flag[5] == 1)
		rhs += coeff_temp*U1_nb[5];
	    else if (flag[3] == 1 && flag[5] == 1)
		rhs += coeff_temp*U1_nb[3];
	    else {
		rhs += coeff_temp*(U1_nb[3]+U1_nb[5]+U1_nb[12]+U1_center)/8.0;

		solver.Add_A(I*3+2,I*3+1,        -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[3]*3+1,  -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[5]*3+1,  -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[12]*3+1, -coeff_temp/8.0);
	    }

	    //corner (i,j-1/2,k+1/2)

	    coeff_temp = -mu_edge[2]*m_dt/(top_h[1]*top_h[2])/rho;

	    if (flag[2] == 1 && flag[5] == 0)
		rhs += coeff_temp*U1_nb[2];
	    else if (flag[2] == 0 && flag[5] == 1)
		rhs += coeff_temp*U1_nb[5];
	    else if (flag[2] == 1 && flag[5] == 1)
		rhs += coeff_temp*U1_nb[2];
	    else {
		rhs += coeff_temp*(U1_nb[2]+U1_nb[5]+U1_nb[13]+U1_center)/8.0;

		solver.Add_A(I*3+2,I*3+1,        -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[2]*3+1,  -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[5]*3+1,  -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[13]*3+1, -coeff_temp/8.0);
	    }

	    rhs += m_dt*state.m_U[2];
	    rhs += m_dt*cell_center[index].m_state.f_surf[2];
	    rhs -= m_dt*cell_center[index].m_state.grad_q[2]/rho;
	    rhs -= m_dt*cell_center[index].m_state.m_adv[2];

	    solver.Set_b(I*3+2, rhs);
        }

        solver.SetMaxIter(40000);
        solver.SetTol(1e-14);

	start_clock("Before Petsc Solve");
        solver.Solve_GMRES();
	solver.GetNumIterations(&num_iter);
	solver.GetFinalRelativeResidualNorm(&rel_residual);

	stop_clock("After Petsc Solve");

        // get back the solution
        solver.Get_x(x);

        if (debugging("PETSc"))
            (void) printf("Incompress_Solver_Smooth_3D_Cartesian::"
			"compDiffWithSmoothProperty_2nd_coupled: "
                        "num_iter = %d, rel_residual = %g. \n",
                        num_iter,rel_residual);


	for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I = ijk_to_I[i][j][k];
            index = d_index3d(i,j,k,top_gmax);
            if (I >= 0)
            {
                cell_center[index].m_state.m_U[0] = x[I*3-ilower*3];
                cell_center[index].m_state.m_U[1] = x[I*3-ilower*3+1];
		cell_center[index].m_state.m_U[2] = x[I*3-ilower*3+2];
                speed = fabs(cell_center[index].m_state.m_U[0]) +
                        fabs(cell_center[index].m_state.m_U[1]) +
			fabs(cell_center[index].m_state.m_U[2]);
                if (speed > max_speed)
                    max_speed = speed;
            }
            else
            {
                cell_center[index].m_state.m_U[0] = 0.0;
                cell_center[index].m_state.m_U[1] = 0.0;
		cell_center[index].m_state.m_U[2] = 0.0;
            }
        }
        for (l = 0; l < 3; ++l)
        {
	    for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)

            {
                index  = d_index3d(i,j,k,top_gmax);
                array[index] = cell_center[index].m_state.m_U[l];
            }
            scatMeshArray();
	    for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)

            {
                index  = d_index3d(i,j,k,top_gmax);
                cell_center[index].m_state.m_U[l] = array[index];
            }
        }
        pp_global_max(&max_speed,1);

        free_these(1,x);
}       /* end compDiffWithSmoothProperty3d_2nd_coupled */

/*******************************************************************
 * 	               iFjunk.cpp	for 2D
 *******************************************************************/

void Incompress_Solver_Smooth_2D_Cartesian::
	compDiffWithSmoothProperty_1st_coupled(void)
{
        COMPONENT comp;
	int index,index_nb[8],size;
	int I,I_nb[8];
	double coords[MAXD],crx_coords[MAXD];
	double coeff[8],mu[4],mu_edge[4],mu0,rho,rhs,U0_nb[8],U1_nb[8];
	int flag[4]; //denote whether this is dirichlet or neumann boundary
	L_STATE state;
	int i,j,k,nb,icoords[MAXD];
	double speed;
	double *x;
	GRID_DIRECTION dir[4] = {WEST,EAST,SOUTH,NORTH};
	POINTER intfc_state;
	HYPER_SURF *hs;
	PetscInt num_iter;
	double rel_residual;

	max_speed = 0.0;
	setIndexMap();

	size = iupper - ilower;
	FT_VectorMemoryAlloc((POINTER*)&x,2*size,sizeof(double));

	PETSc solver;

	solver.Create(2*ilower, 2*iupper-1, 9, 9);
		// two buffer zones, 5 u and 4 v for the first  equation
	        //                   5 v and 4 u for the second equation

	solver.Reset_A();
	solver.Reset_b();
	solver.Reset_x();

	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    I  = ij_to_I[i][j];
	    if (I == -1) continue;

	    index  = d_index2d(i,j,top_gmax);
	    index_nb[0] = d_index2d(i-1,j,top_gmax);
	    index_nb[1] = d_index2d(i+1,j,top_gmax);
	    index_nb[2] = d_index2d(i,j-1,top_gmax);
	    index_nb[3] = d_index2d(i,j+1,top_gmax);
	    index_nb[4] = d_index2d(i-1,j-1,top_gmax);
	    index_nb[5] = d_index2d(i+1,j-1,top_gmax);
	    index_nb[6] = d_index2d(i+1,j+1,top_gmax);
	    index_nb[7] = d_index2d(i-1,j+1,top_gmax);

	    /*
	     *       |		|
	     *   7   |    3	|   6
	     *-------|----------|-------
	     *	     |		|
	     * 	 0   | 		|   1
	     *       |		|
	     *-------|----------|--------
	     *   4   |    2	|   5
	     *       |		|
	     */

	    I_nb[0] = ij_to_I[i-1][j];
	    I_nb[1] = ij_to_I[i+1][j];
	    I_nb[2] = ij_to_I[i][j-1];
	    I_nb[3] = ij_to_I[i][j+1];
	    I_nb[4] = ij_to_I[i-1][j-1];
	    I_nb[5] = ij_to_I[i+1][j-1];
	    I_nb[6] = ij_to_I[i+1][j+1];
	    I_nb[7] = ij_to_I[i-1][j+1];

	    icoords[0] = i;
	    icoords[1] = j;
	    comp = top_comp[index];

	    mu0 = cell_center[index].m_state.m_mu;
	    rho = cell_center[index].m_state.m_rho;

	    for (nb = 0; nb < 4; nb++)
	    {
		if (FT_StateStructAtGridCrossing(front,icoords,dir[nb],
			    comp,&intfc_state,&hs,crx_coords) &&
		            wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
		{
		    flag[nb] = 1;
		    U0_nb[nb] = getStateVel[0](intfc_state);
		    U1_nb[nb] = getStateVel[1](intfc_state);
		    
		    if (wave_type(hs) == DIRICHLET_BOUNDARY ||
		        wave_type(hs) == NEUMANN_BOUNDARY)
		    {
			mu_edge[nb] = mu0;
			mu[nb] = mu0;
		    }
		    else
		    {
			mu_edge[nb] = 0.5*(mu0 + 
				cell_center[index_nb[nb]].m_state.m_mu);
			mu[nb] = cell_center[index_nb[nb]].m_state.m_mu;
		    }
		}
		else
		{
		    flag[nb] = 0;
		    U0_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[0];
		    U1_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[1];
		    mu_edge[nb] = 0.5*(mu0 + 
				cell_center[index_nb[nb]].m_state.m_mu);
		    mu[nb] = cell_center[index_nb[nb]].m_state.m_mu;
		}
	    }

	    // Neighbour cell 4
	    if (flag[0] == 1 && flag[2] == 1)
	    {
		U0_nb[4] = 0.5*(U0_nb[0] + U0_nb[2]);
		U1_nb[4] = 0.5*(U1_nb[0] + U1_nb[2]);
	    }
	    else
	    {
		U0_nb[4] = cell_center[index_nb[4]].m_state.m_U[0];
		U1_nb[4] = cell_center[index_nb[4]].m_state.m_U[1];
	    }
	    // Neighbour cell 5
	    if (flag[1] == 1 && flag[2] == 1)
	    {
		U0_nb[5] = 0.5*(U0_nb[2] + U0_nb[1]);
		U1_nb[5] = 0.5*(U1_nb[2] + U1_nb[1]);
	    }
	    else
	    {
		U0_nb[5] = cell_center[index_nb[5]].m_state.m_U[0];
		U1_nb[5] = cell_center[index_nb[5]].m_state.m_U[1];
	    }
	    // Neighbour cell 6 
	    if (flag[1] == 1 && flag[3] == 1)
	    {
		U0_nb[6] = 0.5*(U0_nb[1] + U0_nb[3]);
		U1_nb[6] = 0.5*(U1_nb[1] + U1_nb[3]);
	    }
	    else
	    {
		U0_nb[6] = cell_center[index_nb[6]].m_state.m_U[0];
		U1_nb[6] = cell_center[index_nb[6]].m_state.m_U[1];
	    }
	    // Neighbour cell 7
	    if (flag[3] == 1 && flag[0] == 1)
	    {
		U0_nb[7] = 0.5*(U0_nb[3] + U0_nb[0]);
		U1_nb[7] = 0.5*(U1_nb[3] + U1_nb[0]);
	    }
	    else
	    {
		U0_nb[7] = cell_center[index_nb[7]].m_state.m_U[0];
		U1_nb[7] = cell_center[index_nb[7]].m_state.m_U[1];
	    }


	    // source term
	    getRectangleCenter(index, coords);
	    computeSourceTerm(coords, state);

	    //Setting the coefficients for the first equation
	    coeff[0] = m_dt/rho * mu_edge[0]/(top_h[0]*top_h[0]);	// down
	    coeff[1] = m_dt/rho * mu_edge[1]/(top_h[0]*top_h[0]);	// right
	    coeff[2] = 0.5*m_dt/rho * mu_edge[2]/(top_h[1]*top_h[1]);	// up	
	    coeff[3] = 0.5*m_dt/rho * mu_edge[3]/(top_h[1]*top_h[1]);	// left

	    coeff[4] = 0.5*m_dt/rho * mu[2]/(4*top_h[0]*top_h[1]);
	    coeff[5] = -0.5*m_dt/rho * mu[2]/(4*top_h[0]*top_h[1]);
	    coeff[6] = 0.5*m_dt/rho * mu[3]/(4*top_h[0]*top_h[1]);
	    coeff[7] = -0.5*m_dt/rho * mu[3]/(4*top_h[0]*top_h[1]);

	    solver.Set_A(I*2,I*2,1+coeff[0]+coeff[1]+coeff[2]+coeff[3]);
	    rhs = (1-coeff[0]-coeff[1]-coeff[2]-coeff[3])*
				cell_center[index].m_state.m_U[0];
	    

	    for (nb = 0; nb < 4; nb++)
	    {
		if (I_nb[nb] != -1)
		{
		    solver.Set_A(I*2,I_nb[nb]*2,-coeff[nb]);
		    rhs += coeff[nb]*U0_nb[nb];
		}
		else
		    rhs += 2.0*coeff[nb]*U0_nb[nb];
	    }
	    for (nb = 4; nb < 8; nb++)
	    {
		if (I_nb[nb] != -1)
		{
		    solver.Set_A(I*2,I_nb[nb]*2+1,-coeff[nb]);
		    rhs += coeff[nb]*U1_nb[nb];
		}
		else
		    rhs += 2.0*coeff[nb]*U1_nb[nb];
	    }

	    rhs += m_dt*state.m_U[0];
	    rhs += m_dt*cell_center[index].m_state.f_surf[0];
	    rhs -= m_dt*cell_center[index].m_state.grad_q[0]/rho;

	    solver.Set_b(I*2, rhs);


	    /****************************************************************/
	    //Setting the coefficients for the second equation

	    coeff[0] = 0.5*m_dt/rho * mu_edge[0]/(top_h[0]*top_h[0]);	// down
	    coeff[1] = 0.5*m_dt/rho * mu_edge[1]/(top_h[0]*top_h[0]);	// right
	    coeff[2] = m_dt/rho * mu_edge[2]/(top_h[1]*top_h[1]);	// up	
	    coeff[3] = m_dt/rho * mu_edge[3]/(top_h[1]*top_h[1]);	// left

	    coeff[4] = 0.5*m_dt/rho * mu[0]/(4*top_h[0]*top_h[1]);
	    coeff[5] = -0.5*m_dt/rho * mu[1]/(4*top_h[0]*top_h[1]);
	    coeff[6] = 0.5*m_dt/rho * mu[1]/(4*top_h[0]*top_h[1]);
	    coeff[7] = -0.5*m_dt/rho * mu[0]/(4*top_h[0]*top_h[1]);


	    solver.Set_A(I*2+1,I*2+1,1+coeff[0]+coeff[1]+coeff[2]+coeff[3]);
	    rhs = (1-coeff[0]-coeff[1]-coeff[2]-coeff[3])*
				cell_center[index].m_state.m_U[1];

	    for (nb = 0; nb < 4; nb++)
	    {
		if (I_nb[nb] != -1)
		{
		    solver.Set_A(I*2+1,I_nb[nb]*2+1,-coeff[nb]);
		    rhs += coeff[nb]*U1_nb[nb];
		}
		else
		    rhs += 2.0*coeff[nb]*U1_nb[nb];
	    }
	    for (nb = 4; nb < 8; nb++)
	    {
		if (I_nb[nb] != -1)
		{
		    solver.Set_A(I*2+1,I_nb[nb]*2,-coeff[nb]);
		    rhs += coeff[nb]*U0_nb[nb];
		}
		else
		    rhs += 2.0*coeff[nb]*U0_nb[nb];
	    }

	    rhs += m_dt*state.m_U[1];
	    rhs += m_dt*cell_center[index].m_state.f_surf[1];
	    rhs -= m_dt*cell_center[index].m_state.grad_q[1]/rho;

	    solver.Set_b(I*2+1, rhs);
	}

	solver.SetMaxIter(40000);
	solver.SetTol(1e-14);

	start_clock("Before Petsc solver");
	
	solver.Solve_GMRES();
	solver.GetNumIterations(&num_iter);
	solver.GetFinalRelativeResidualNorm(&rel_residual);

	//if (rel_residual > 1)
	//{
	//    solver.Reset_x();
	//    solver.Solve_GMRES();
	//    solver.GetNumIterations(&num_iter);
	//    solver.GetFinalRelativeResidualNorm(&rel_residual);
	//}

	stop_clock("After Petsc solver");

	// get back the solution
	solver.Get_x(x);

	if (debugging("PETSc"))
	{
	    printf("\nIncompress_Solver_Smooth_2D_Cartesian::"
			"compDiffWithSmoothProperty: "
	       		"num_iter = %d, rel_residual = %g. \n", 
			num_iter,rel_residual); 
	}

	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    I = ij_to_I[i][j];
	    index = d_index2d(i,j,top_gmax);
	    if (I >= 0)
	    {
		cell_center[index].m_state.m_U[0] = x[I*2-ilower*2];
		cell_center[index].m_state.m_U[1] = x[I*2+1-ilower*2];
		speed = fabs(cell_center[index].m_state.m_U[0]) +
		    	fabs(cell_center[index].m_state.m_U[1]);
		if (speed > max_speed)
		    max_speed = speed;
	    }
	    else
	    {
		cell_center[index].m_state.m_U[0] = 0.0;
		cell_center[index].m_state.m_U[1] = 0.0;
		speed = fabs(cell_center[index].m_state.m_U[0]) +
		    	fabs(cell_center[index].m_state.m_U[1]);
		if (speed > max_speed)
		    max_speed = speed;
	    }
	}
	for (k = 0; k < 2; ++k) //scatter the velocity
	{
	    for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
	    {	
	    	index  = d_index2d(i,j,top_gmax);
	    	array[index] = cell_center[index].m_state.m_U[k];
	    }
	    scatMeshArray();
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {	
	    	index  = d_index2d(i,j,top_gmax);
	    	cell_center[index].m_state.m_U[k] = array[index];
	    }
	}
	pp_global_max(&max_speed,1);
	FT_FreeThese(1,x);
}	/* end compDiffWithSmoothProperty2D */


void Incompress_Solver_Smooth_2D_Cartesian::
	compDiffWithSmoothProperty_2nd_coupled(void)
{
        COMPONENT comp;
	int index,index_nb[8],size;
	int I,I_nb[8];
	double coords[MAXD],crx_coords[MAXD];
	double coeff0[4],coeff1[4],coeff_temp;
	double mu[4],mu_edge[4],mu0,rho,rhs;
	double U0_nb[8],U1_nb[8],U0_center, U1_center;
	int flag[4]; //denote whether this is dirichlet or neumann boundary
	L_STATE state;
	int i,j,k,nb,icoords[MAXD];
	double speed;
	double *x;
	GRID_DIRECTION dir[4] = {WEST,EAST,SOUTH,NORTH};
	POINTER intfc_state;
	HYPER_SURF *hs;
	PetscInt num_iter;
	double rel_residual;

	max_speed = 0.0;
	setIndexMap();

	size = iupper - ilower;
	FT_VectorMemoryAlloc((POINTER*)&x,2*size,sizeof(double));

	PETSc solver;

	solver.Create(2*ilower, 2*iupper-1, 14, 14);
		// two buffer zones, 5 u and 9 v for the first  equation
	        //                   5 v and 9 u for the second equation

	solver.Reset_A();
	solver.Reset_b();
	solver.Reset_x();

	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    I  = ij_to_I[i][j];
	    if (I == -1) continue;

	    index  = d_index2d(i,j,top_gmax);
	    index_nb[0] = d_index2d(i-1,j,top_gmax);
	    index_nb[1] = d_index2d(i+1,j,top_gmax);
	    index_nb[2] = d_index2d(i,j-1,top_gmax);
	    index_nb[3] = d_index2d(i,j+1,top_gmax);
	    index_nb[4] = d_index2d(i-1,j-1,top_gmax);
	    index_nb[5] = d_index2d(i+1,j-1,top_gmax);
	    index_nb[6] = d_index2d(i+1,j+1,top_gmax);
	    index_nb[7] = d_index2d(i-1,j+1,top_gmax);

	    /*
	     *       |		|
	     *   7   |    3	|   6
	     *-------|----------|-------
	     *	     |		|
	     * 	 0   | 		|   1
	     *       |		|
	     *-------|----------|--------
	     *   4   |    2	|   5
	     *       |		|
	     */

	    I_nb[0] = ij_to_I[i-1][j];
	    I_nb[1] = ij_to_I[i+1][j];
	    I_nb[2] = ij_to_I[i][j-1];
	    I_nb[3] = ij_to_I[i][j+1];
	    I_nb[4] = ij_to_I[i-1][j-1];
	    I_nb[5] = ij_to_I[i+1][j-1];
	    I_nb[6] = ij_to_I[i+1][j+1];
	    I_nb[7] = ij_to_I[i-1][j+1];

	    icoords[0] = i;
	    icoords[1] = j;
	    comp = top_comp[index];

	    mu0 = cell_center[index].m_state.m_mu;
	    rho = cell_center[index].m_state.m_rho;
	    U0_center = cell_center[index].m_state.m_U[0];
	    U1_center = cell_center[index].m_state.m_U[1];

	    for (nb = 0; nb < 4; nb++)
	    {
		if (FT_StateStructAtGridCrossing(front,icoords,dir[nb],
			    comp,&intfc_state,&hs,crx_coords) &&
		            wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
		{
		    flag[nb] = 1;
		    U0_nb[nb] = getStateVel[0](intfc_state);
		    U1_nb[nb] = getStateVel[1](intfc_state);
		    
		    if (wave_type(hs) == DIRICHLET_BOUNDARY ||
		        wave_type(hs) == NEUMANN_BOUNDARY)
		    {
			mu_edge[nb] = mu0;
			mu[nb] = mu0;
		    }
		    else
		    {
			mu_edge[nb] = 0.5*(mu0 + 
				cell_center[index_nb[nb]].m_state.m_mu);
			mu[nb] = cell_center[index_nb[nb]].m_state.m_mu;
		    }
		}
		else
		{
		    flag[nb] = 0;
		    U0_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[0];
		    U1_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[1];
		    mu_edge[nb] = 0.5*(mu0 + 
				cell_center[index_nb[nb]].m_state.m_mu);
		    mu[nb] = cell_center[index_nb[nb]].m_state.m_mu;
		}
	    }

	    for (nb = 4; nb < 8; nb++) // corner cell value, interior
	    {
		if (I_nb[nb] != -1)
		{
		    U0_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[0];
		    U1_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[1];
		}
	    }

	    // source term
	    getRectangleCenter(index, coords);
	    computeSourceTerm(coords, state);

	    //Setting the coefficients and matrix for the U0 in first equation
	    coeff0[0] = m_dt/rho * mu_edge[0]/(top_h[0]*top_h[0]);	// down
	    coeff0[1] = m_dt/rho * mu_edge[1]/(top_h[0]*top_h[0]);	// right
	    coeff0[2] = 0.5*m_dt/rho * mu_edge[2]/(top_h[1]*top_h[1]);	// up	
	    coeff0[3] = 0.5*m_dt/rho * mu_edge[3]/(top_h[1]*top_h[1]);	// left


	    solver.Set_A(I*2,I*2,1.0);
	    rhs = U0_center;

	    for (nb = 0; nb < 4; nb++)
	    {
		if (I_nb[nb] != -1)
		{
		    solver.Set_A(I*2,I_nb[nb]*2,-coeff0[nb]);
		    rhs += coeff0[nb]*U0_nb[nb];
		}
		else
		{
		    coeff0[nb] = 2.0*coeff0[nb];
		    rhs += 2.0*coeff0[nb]*U0_nb[nb];
		}
		
		solver.Add_A(I*2,I*2,coeff0[nb]);
		rhs -= coeff0[nb]*U0_center;
	    }

	    //set the coefficients and matrix for U1 in first equation
	    //traverse the four corners

	    //corner 4

	    coeff_temp = mu_edge[2]*m_dt/(top_h[0]*top_h[1])/rho;

	    if (flag[0] == 1 && flag[2] == 0)
		rhs += coeff_temp*U1_nb[0];
	    else if(flag[0] == 0 && flag[2] ==1)
		rhs += coeff_temp*U1_nb[2];
	    else if(flag[0] ==1 && flag[2] == 1)
		rhs += coeff_temp*U1_nb[0];
	    else {
		rhs += coeff_temp*(U1_nb[0]+U1_nb[2]+U1_nb[4]+U1_center)/8.0;

		solver.Add_A(I*2,I*2+1,       -coeff_temp/8.0);
		solver.Add_A(I*2,I_nb[0]*2+1, -coeff_temp/8.0);
		solver.Add_A(I*2,I_nb[2]*2+1, -coeff_temp/8.0);
		solver.Add_A(I*2,I_nb[4]*2+1, -coeff_temp/8.0);
	    }

	    //corner 5

	    coeff_temp = -mu_edge[2]*m_dt/(top_h[0]*top_h[1])/rho;

	    if (flag[2] == 1 && flag[1] == 0)
		rhs += coeff_temp*U1_nb[2];
	    else if(flag[2] == 0 && flag[1] == 1)
		rhs += coeff_temp*U1_nb[1];
	    else if(flag[2] == 1 && flag[1] == 1)
		rhs += coeff_temp*U1_nb[2];
	    else {
		rhs += coeff_temp*(U1_nb[2]+U1_nb[1]+U1_nb[5]+U1_center)/8.0;

		solver.Add_A(I*2,I*2+1,        -coeff_temp/8.0);
		solver.Add_A(I*2,I_nb[2]*2+1,  -coeff_temp/8.0);
		solver.Add_A(I*2,I_nb[1]*2+1,  -coeff_temp/8.0);
		solver.Add_A(I*2,I_nb[5]*2+1,  -coeff_temp/8.0);
	    }

	    //corner 6

	    coeff_temp = mu_edge[3]*m_dt/(top_h[0]*top_h[1])/rho;

	    if (flag[1] == 1 && flag[3] == 0)
		rhs += coeff_temp*U1_nb[1];
	    else if(flag[1] == 0 && flag[3] == 1)
		rhs += coeff_temp*U1_nb[3];
	    else if(flag[1] == 1 && flag[3] == 1)
		rhs += coeff_temp*U1_nb[1];
	    else {
		rhs += coeff_temp*(U1_nb[1]+U1_nb[3]+U1_nb[6]+U1_center)/8.0;

		solver.Add_A(I*2,I*2+1,        -coeff_temp/8.0);
		solver.Add_A(I*2,I_nb[1]*2+1,  -coeff_temp/8.0);
		solver.Add_A(I*2,I_nb[3]*2+1,  -coeff_temp/8.0);
		solver.Add_A(I*2,I_nb[6]*2+1,  -coeff_temp/8.0);
	    }

	    //corner 7

	    coeff_temp = -mu_edge[3]*m_dt/(top_h[0]*top_h[1])/rho;

	    if (flag[3] == 1 && flag[0] == 0)
		rhs += coeff_temp*U1_nb[3];
	    else if(flag[3] == 0 && flag[0] == 1)
		rhs += coeff_temp*U1_nb[0];
	    else if(flag[3] == 1 && flag[0] == 1)
		rhs += coeff_temp*U1_nb[3];
	    else {
		rhs += coeff_temp*(U1_nb[3]+U1_nb[0]+U1_nb[7]+U1_center)/8.0;

		solver.Add_A(I*2,I*2+1,        -coeff_temp/8.0);
		solver.Add_A(I*2,I_nb[3]*2+1,  -coeff_temp/8.0);
		solver.Add_A(I*2,I_nb[0]*2+1,  -coeff_temp/8.0);
		solver.Add_A(I*2,I_nb[7]*2+1,  -coeff_temp/8.0);
	    }

	    rhs += m_dt*state.m_U[0];
	    rhs += m_dt*cell_center[index].m_state.f_surf[0];
	    rhs -= m_dt*cell_center[index].m_state.grad_q[0]/rho;
	    rhs -= m_dt*cell_center[index].m_state.m_adv[0];

	    solver.Set_b(I*2, rhs);


	    /****************************************************************/
	    //Setting the coefficients for U1 in the second equation

	    coeff1[0] = 0.5*m_dt/rho * mu_edge[0]/(top_h[0]*top_h[0]);	// down
	    coeff1[1] = 0.5*m_dt/rho * mu_edge[1]/(top_h[0]*top_h[0]);	// right
	    coeff1[2] = m_dt/rho * mu_edge[2]/(top_h[1]*top_h[1]);	// up	
	    coeff1[3] = m_dt/rho * mu_edge[3]/(top_h[1]*top_h[1]);	// left


	    solver.Set_A(I*2+1,I*2+1,1.0);
	    rhs = U1_center;


	    for (nb = 0; nb < 4; nb++)
	    {
		if (I_nb[nb] != -1)
		{
		    solver.Set_A(I*2+1,I_nb[nb]*2+1,-coeff1[nb]);
		    rhs += coeff1[nb]*U1_nb[nb];
		}
		else
		{
		    coeff1[nb] = 2.0*coeff1[nb];
		    rhs += 2.0*coeff1[nb]*U1_nb[nb];
		}

		solver.Add_A(I*2+1,I*2+1,coeff1[nb]);
		rhs -= coeff1[nb]*U1_center;
	    }

	    //set the coefficients and matrix for U0 in second equation
	    //traverse the four corners

	    //corner 4

	    coeff_temp = mu_edge[0]*m_dt/(top_h[0]*top_h[1])/rho;

	    if (flag[0] == 1 && flag[2] == 0)
		rhs += coeff_temp*U0_nb[0];
	    else if(flag[0] == 0 && flag[2] ==1)
		rhs += coeff_temp*U0_nb[2];
	    else if(flag[0] ==1 && flag[2] == 1)
		rhs += coeff_temp*U0_nb[0];
	    else {
		rhs += coeff_temp*(U0_nb[0]+U0_nb[2]+U0_nb[4]+U0_center)/8.0;

		solver.Add_A(I*2+1,I*2,       -coeff_temp/8.0);
		solver.Add_A(I*2+1,I_nb[0]*2, -coeff_temp/8.0);
		solver.Add_A(I*2+1,I_nb[2]*2, -coeff_temp/8.0);
		solver.Add_A(I*2+1,I_nb[4]*2, -coeff_temp/8.0);
	    }

	    //corner 5

	    coeff_temp = -mu_edge[1]*m_dt/(top_h[0]*top_h[1])/rho;

	    if (flag[2] == 1 && flag[1] == 0)
		rhs += coeff_temp*U0_nb[2];
	    else if(flag[2] == 0 && flag[1] ==1)
		rhs += coeff_temp*U0_nb[1];
	    else if(flag[2] ==1 && flag[1] == 1)
		rhs += coeff_temp*U0_nb[2];
	    else {
		rhs += coeff_temp*(U0_nb[2]+U0_nb[1]+U0_nb[5]+U0_center)/8.0;

		solver.Add_A(I*2+1,I*2,       -coeff_temp/8.0);
		solver.Add_A(I*2+1,I_nb[2]*2, -coeff_temp/8.0);
		solver.Add_A(I*2+1,I_nb[1]*2, -coeff_temp/8.0);
		solver.Add_A(I*2+1,I_nb[5]*2, -coeff_temp/8.0);
	    }
	    //corner 6

	    coeff_temp = mu_edge[1]*m_dt/(top_h[0]*top_h[1])/rho;

	    if (flag[1] == 1 && flag[3] == 0)
		rhs += coeff_temp*U0_nb[1];
	    else if(flag[1] == 0 && flag[3] ==1)
		rhs += coeff_temp*U0_nb[3];
	    else if(flag[1] ==1 && flag[3] == 1)
		rhs += coeff_temp*U0_nb[1];
	    else {
		rhs += coeff_temp*(U0_nb[1]+U0_nb[3]+U0_nb[6]+U0_center)/8.0;

		solver.Add_A(I*2+1,I*2,       -coeff_temp/8.0);
		solver.Add_A(I*2+1,I_nb[1]*2, -coeff_temp/8.0);
		solver.Add_A(I*2+1,I_nb[3]*2, -coeff_temp/8.0);
		solver.Add_A(I*2+1,I_nb[6]*2, -coeff_temp/8.0);
	    }
	    //corner 7

	    coeff_temp = -mu_edge[0]*m_dt/(top_h[0]*top_h[1])/rho;

	    if (flag[3] == 1 && flag[0] == 0)
		rhs += coeff_temp*U0_nb[3];
	    else if(flag[3] == 0 && flag[0] ==1)
		rhs += coeff_temp*U0_nb[0];
	    else if(flag[3] ==1 && flag[0] == 1)
		rhs += coeff_temp*U0_nb[3];
	    else {
		rhs += coeff_temp*(U0_nb[3]+U0_nb[0]+U0_nb[7]+U0_center)/8.0;

		solver.Add_A(I*2+1,I*2,       -coeff_temp/8.0);
		solver.Add_A(I*2+1,I_nb[3]*2, -coeff_temp/8.0);
		solver.Add_A(I*2+1,I_nb[0]*2, -coeff_temp/8.0);
		solver.Add_A(I*2+1,I_nb[7]*2, -coeff_temp/8.0);
	    }
	   
   	    rhs += m_dt*state.m_U[1];
	    rhs += m_dt*cell_center[index].m_state.f_surf[1];
	    rhs -= m_dt*cell_center[index].m_state.grad_q[1]/rho;
	    rhs -= m_dt*cell_center[index].m_state.m_adv[1];

	    solver.Set_b(I*2+1, rhs);
	}

	solver.SetMaxIter(40000);
	solver.SetTol(1e-14);

	start_clock("Before Petsc solver");
	
	solver.Solve_GMRES();
	solver.GetNumIterations(&num_iter);
	solver.GetFinalRelativeResidualNorm(&rel_residual);

	//if (rel_residual > 1)
	//{
	//    solver.Reset_x();
	//    solver.Solve_GMRES();
	//    solver.GetNumIterations(&num_iter);
	//    solver.GetFinalRelativeResidualNorm(&rel_residual);
	//}

	stop_clock("After Petsc solver");

	// get back the solution
	solver.Get_x(x);

	if (debugging("PETSc"))
	{
	    printf("\nIncompress_Solver_Smooth_2D_Cartesian::"
			"compDiffWithSmoothProperty_2nd_coupled: "
	       		"num_iter = %d, rel_residual = %g. \n", 
			num_iter,rel_residual); 
	}

	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    I = ij_to_I[i][j];
	    index = d_index2d(i,j,top_gmax);
	    if (I >= 0)
	    {
		cell_center[index].m_state.m_U[0] = x[I*2-ilower*2];
		cell_center[index].m_state.m_U[1] = x[I*2+1-ilower*2];
		speed = fabs(cell_center[index].m_state.m_U[0]) +
		    	fabs(cell_center[index].m_state.m_U[1]);
		if (speed > max_speed)
		    max_speed = speed;
	    }
	    else
	    {
		cell_center[index].m_state.m_U[0] = 0.0;
		cell_center[index].m_state.m_U[1] = 0.0;
		speed = fabs(cell_center[index].m_state.m_U[0]) +
		    	fabs(cell_center[index].m_state.m_U[1]);
		if (speed > max_speed)
		    max_speed = speed;
	    }
	}
	for (k = 0; k < 2; ++k) //scatter the velocity
	{
	    for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
	    {	
	    	index  = d_index2d(i,j,top_gmax);
	    	array[index] = cell_center[index].m_state.m_U[k];
	    }
	    scatMeshArray();
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {	
	    	index  = d_index2d(i,j,top_gmax);
	    	cell_center[index].m_state.m_U[k] = array[index];
	    }
	}
	pp_global_max(&max_speed,1);
	free_these(1,x);
}	/* end compDiffWithSmoothProperty2D */


void Incompress_Solver_Smooth_2D_Cartesian::
	compDiffWithSmoothProperty_2nd_decoupled(void)
{
    COMPONENT comp;
    int index,index_nb[4],size;
    int I,I_nb[4];
    int i,j,l,nb,icoords[MAXD];
    L_STATE source_term;
    double coords[MAXD],crx_coords[MAXD];
    double mu[4],mu0,rho,rhs;

    // U_nb contains states at neighbor cell or states on the boundary.
    double U_nb[4],U_nb_new[4], U_center;

    double speed;
    double *x;
    GRID_DIRECTION dir[4] = {WEST,EAST,SOUTH,NORTH};
    POINTER intfc_state;
    HYPER_SURF *hs;

    max_speed = 0.0;

    setIndexMap();

    size = iupper - ilower;
    FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));

    for (l = 0; l < dim; ++l)
    {
	PETSc solver;
	solver.Create(ilower, iupper-1, 5, 5);
	solver.Reset_A();
	solver.Reset_b();
	solver.Reset_x();

	for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		I  = ij_to_I[i][j];
		if (I == -1) continue;

		index  = d_index2d(i,j,top_gmax);
		index_nb[0] = d_index2d(i-1,j,top_gmax);
		index_nb[1] = d_index2d(i+1,j,top_gmax);
		index_nb[2] = d_index2d(i,j-1,top_gmax);
		index_nb[3] = d_index2d(i,j+1,top_gmax);
		icoords[0] = i;
		icoords[1] = j;
		comp = top_comp[index];

		I_nb[0] = ij_to_I[i-1][j]; // left or west
		I_nb[1] = ij_to_I[i+1][j]; // right or east
		I_nb[2] = ij_to_I[i][j-1]; // down or south
		I_nb[3] = ij_to_I[i][j+1]; // up or north


		mu0 = cell_center[index].m_state.m_mu;
		rho = cell_center[index].m_state.m_rho;
		U_center =  cell_center[index].m_state.m_U[l];

		for (nb = 0; nb < 4; nb++)
		{
		    if (FT_StateStructAtGridCrossing(front,icoords,dir[nb],
			    comp,&intfc_state,&hs,crx_coords) &&
			    wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
		    {
			// old boundary condition
			U_nb[nb] = getStateVel[l](intfc_state);
			// new boundary condition
			FT_StateStructAtGridCrossing(front,icoords,dir[nb],
				comp,&intfc_state,&hs,crx_coords);
			U_nb_new[nb] = getStateVel[l](intfc_state);

			if (wave_type(hs) == DIRICHLET_BOUNDARY ||
				wave_type(hs) == NEUMANN_BOUNDARY)
			    mu[nb] = mu0;
			else
			    mu[nb] = 1.0/2*(mu0 +
				    cell_center[index_nb[nb]].m_state.m_mu);
		    }
		    else
		    {
			U_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[l];
			mu[nb] = 1.0/2*(mu0 +
				cell_center[index_nb[nb]].m_state.m_mu);
		    }
		}

		getRectangleCenter(index, coords);
		computeSourceTerm(coords, source_term);

		rhs = 0;
		double dh[2] = {top_h[0], top_h[1]};
		for(nb = 0; nb < 4; nb++)
		{
		    // use dh[1] as the edge length
		    if(nb >= 2)
			std::swap(dh[0],dh[1]);

		    if(I_nb[nb] >= 0)
		    {
			// u^{*}
			solver.Add_A(I,I,
				(-1)*-0.5*m_dt/rho*mu[nb]/(dh[0]*dh[0]));
			solver.Add_A(I,I_nb[nb],
				(-1)*0.5*m_dt/rho*mu[nb]/(dh[0]*dh[0]));
			// u^{n}
			rhs += 0.5*m_dt/rho*mu[nb] *
				(U_nb[nb]-U_center) /(dh[0]*dh[0]);
		    }
		    else		// boundary
		    {
			// u^{*}
			solver.Add_A(I,I,
				(-1)*-0.5*m_dt/rho*mu[nb]/(dh[0]*dh[0]) * 2);
			rhs += 0.5*m_dt/rho*mu[nb]/(dh[0]*dh[0])*2*U_nb_new[nb];
			// u^{n}
			rhs += 0.5*m_dt/rho*mu[nb] *
				(U_nb[nb]-U_center) /(dh[0]*dh[0]) * 2;
		    }
		}

		// interior
		rhs += m_dt*source_term.m_U[l];
		rhs += m_dt*cell_center[index].m_state.f_surf[l];

		rhs -= m_dt*cell_center[index].m_state.grad_q[l]/rho;
		rhs -= m_dt * cell_center[index].m_state.m_adv[l];

		solver.Add_A(I, I, 1.0);
		rhs += U_center;
		solver.Add_b(I, rhs);
	    }

	solver.SetMaxIter(40000);
	solver.SetTol(1e-15);
	solver.Solve_GMRES();

	// get back the solution
	solver.Get_x(x);

	int num_iter;
	double rel_residual;
	solver.GetNumIterations(&num_iter);
	solver.GetFinalRelativeResidualNorm(&rel_residual);

	if (debugging("PETSc"))
	    (void) printf("Incompress_Solver_Smooth_2D_Cartesian::"
		    "compDiffWithSmoothProperty_2nd_decoupled: "
		    "num_iter = %d, rel_residual = %g. \n",
		    num_iter,rel_residual);

	for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		I = ij_to_I[i][j];
		index = d_index2d(i,j,top_gmax);
		if (I >= 0)
		{
		    cell_center[index].m_state.m_U[l] = x[I-ilower];
		}
		else
		{
		    cell_center[index].m_state.m_U[l] = 0.0;
		}
	    }

	for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		index  = d_index2d(i,j,top_gmax);
		array[index] = cell_center[index].m_state.m_U[l];
	    }
	scatMeshArray();
	for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {
		index  = d_index2d(i,j,top_gmax);
		cell_center[index].m_state.m_U[l] = array[index];
	    }
    }

    for (j = jmin; j <= jmax; j++)
	for (i = imin; i <= imax; i++)
	{
	    index = d_index2d(i,j,top_gmax);
	    speed = fabs(cell_center[index].m_state.m_U[0]) +
		    fabs(cell_center[index].m_state.m_U[1]);
	    if (speed > max_speed)
		max_speed = speed;
	}
    pp_global_max(&max_speed,1);

    free_these(1,x);
}

/* Old computeAdvection functions */

void Incompress_Solver_Smooth_2D_Cartesian::computeAdvection(void)
{
	int i,j,k,index,index00,index01,index10,index11,size;
	L_STATE state;
	COMPONENT comp;
	double speed;
	double *u, *v;
	double u0,v0,u00,u01,u10,u11,v00,v01,v10,v11;
	double crx_coords[MAXD];
	int icoords[MAXD];
	POINTER intfc_state;
	HYPER_SURF *hs;

	max_speed = 0.0;
	size = (top_gmax[0]+1)*(top_gmax[1]+1);
	FT_VectorMemoryAlloc((POINTER*)&u,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&v,size,sizeof(double));

	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index = d_index2d(i,j,top_gmax);
	    u[index] = cell_center[index].m_state.m_U[0];
	    v[index] = cell_center[index].m_state.m_U[1];
	}

	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{	
	    icoords[0] = i;
	    icoords[1] = j;
	    index  = d_index2d(i,j,top_gmax);
	    comp = top_comp[index];
	    if (comp == SOLID_COMP)
	    {
	    	cell_center[index].m_state.m_U[0] = 0.0; 
	    	cell_center[index].m_state.m_U[1] = 0.0; 
		continue;
	    }
	    u0 = u[index];
	    v0 = v[index];
	    index00 = d_index2d(i-1,j,top_gmax);
	    if ((*findStateAtCrossing)(front,icoords,WEST,comp,
			        &intfc_state,&hs,crx_coords)) 
	    {
	        u00 = getStateXvel(intfc_state);
	        v00 = getStateYvel(intfc_state);
		if (wave_type(hs) == DIRICHLET_BOUNDARY &&
                    boundary_state_function(hs) &&
                    strcmp(boundary_state_function_name(hs),
                    "flowThroughBoundaryState") == 0)
                {
                    u00 = u0;
                    v00 = v0;
                }
	    }
	    else
	    {
		u00 = u[index00];
		v00 = v[index00];
	    }
	    index01 = d_index2d(i+1,j,top_gmax);
	    if ((*findStateAtCrossing)(front,icoords,EAST,comp,
			        &intfc_state,&hs,crx_coords)) 
	    {
	    	u01 = getStateXvel(intfc_state);
	    	v01 = getStateYvel(intfc_state);
		if (wave_type(hs) == DIRICHLET_BOUNDARY &&
                    boundary_state_function(hs) &&
                    strcmp(boundary_state_function_name(hs),
                    "flowThroughBoundaryState") == 0)
                {
                    u01 = u0;
                    v01 = v0;
                }
	    }
	    else
	    {
		u01 = u[index01];
		v01 = v[index01];
	    }
	    index10 = d_index2d(i,j-1,top_gmax);
	    if ((*findStateAtCrossing)(front,icoords,SOUTH,comp,
			        &intfc_state,&hs,crx_coords)) 
	    {
	    	u10 = getStateXvel(intfc_state);
	    	v10 = getStateYvel(intfc_state);
		if (wave_type(hs) == DIRICHLET_BOUNDARY &&
                    boundary_state_function(hs) &&
                    strcmp(boundary_state_function_name(hs),
                    "flowThroughBoundaryState") == 0)
                {
                    u10 = u0;
                    v10 = v0;
                }
	    }
	    else
	    {
		u10 = u[index10];
		v10 = v[index10];
	    }
	    index11 = d_index2d(i,j+1,top_gmax);
	    if ((*findStateAtCrossing)(front,icoords,NORTH,comp,
			        &intfc_state,&hs,crx_coords)) 
	    {
		u11 = getStateXvel(intfc_state);
		v11 = getStateYvel(intfc_state);
		 if (wave_type(hs) == DIRICHLET_BOUNDARY &&
                    boundary_state_function(hs) &&
                    strcmp(boundary_state_function_name(hs),
                    "flowThroughBoundaryState") == 0)
                {
                    u11 = u0;
                    v11 = v0;
                }
	    }
	    else
	    {
		u11 = u[index11];
		v11 = v[index11];
	    }

	    cell_center[index].m_state.m_U[0] += -m_dt*(
				burger_flux(u00,u0,u01)/top_h[0] +
				linear_flux(v0,u10,u0,u11)/top_h[1]);

	    cell_center[index].m_state.m_U[1] += - m_dt*(
				linear_flux(u0,v00,v0,v01)/top_h[0] +
			 	burger_flux(v10,v0,v11)/top_h[1]);

	    speed = fabs(cell_center[index].m_state.m_U[0]) +
		    fabs(cell_center[index].m_state.m_U[1]);
	    if (speed > max_speed)
		max_speed = speed;
	}
	for (k = 0; k < 2; ++k)
	{
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {	
	    	index  = d_index2d(i,j,top_gmax);
	    	array[index] = cell_center[index].m_state.m_U[k];
	    }
	    scatMeshArray();
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {	
	    	index  = d_index2d(i,j,top_gmax);
	    	cell_center[index].m_state.m_U[k] = array[index];
	    }
	}
	pp_global_max(&max_speed,1);
	FT_FreeThese(2,u,v);
}	/* end computeAdvection2d */

void Incompress_Solver_Smooth_3D_Cartesian::computeAdvection(void)
{
	int i,j,k,l;
	int index,index00,index01,index10,index11,index20,index21,size;
	L_STATE state;
	COMPONENT comp;
	double speed;
	double *u, *v, *w;
	double u0,u00,u01,u10,u11,u20,u21;
	double v0,v00,v01,v10,v11,v20,v21;
	double w0,w00,w01,w10,w11,w20,w21;
	double crx_coords[MAXD];
	int icoords[MAXD];
	POINTER intfc_state;
	HYPER_SURF *hs;

	if (debugging("trace"))
	    (void) printf("Entering Incompress_Solver_Smooth_3D_Cartesian::"
			  "computeAdvection()\n");
	max_speed = 0.0;
	size = (top_gmax[0]+1)*(top_gmax[1]+1)*(top_gmax[2]+1);
	FT_VectorMemoryAlloc((POINTER*)&u,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&v,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&w,size,sizeof(double));

	for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    u[index] = cell_center[index].m_state.m_U[0];
	    v[index] = cell_center[index].m_state.m_U[1];
	    w[index] = cell_center[index].m_state.m_U[2];
	}

	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{	
	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;
	    index  = d_index3d(i,j,k,top_gmax);
	    comp = top_comp[index];
	    if (comp == SOLID_COMP)
	    {
	    	cell_center[index].m_state.m_U[0] = 0.0; 
	    	cell_center[index].m_state.m_U[1] = 0.0; 
	    	cell_center[index].m_state.m_U[2] = 0.0; 
		continue;
	    }
	    u0 = u[index];
	    v0 = v[index];
	    w0 = w[index];
	    index00 = d_index3d(i-1,j,k,top_gmax);
	    if ((*findStateAtCrossing)(front,icoords,WEST,comp,
			        &intfc_state,&hs,crx_coords))
	    {
		u00 = getStateXvel(intfc_state);
		v00 = getStateYvel(intfc_state);
		w00 = getStateZvel(intfc_state);
		if (wave_type(hs) == DIRICHLET_BOUNDARY &&
		    boundary_state_function(hs) &&
		    strcmp(boundary_state_function_name(hs),
                    "flowThroughBoundaryState") == 0)
		{
		    u00 = u0;
		    v00 = v0;
		    w00 = w0;
		}
	    }
	    else
	    {
		u00 = u[index00];
		v00 = v[index00];
		w00 = w[index00];
	    }
	    index01 = d_index3d(i+1,j,k,top_gmax);
	    if ((*findStateAtCrossing)(front,icoords,EAST,comp,
			        &intfc_state,&hs,crx_coords)) 
	    {
		u01 = getStateXvel(intfc_state);
		v01 = getStateYvel(intfc_state);
		w01 = getStateZvel(intfc_state);
		if (wave_type(hs) == DIRICHLET_BOUNDARY &&
		    boundary_state_function(hs) &&
		    strcmp(boundary_state_function_name(hs),
                    "flowThroughBoundaryState") == 0)
		{
		    u01 = u0;
		    v01 = v0;
		    w01 = w0;
		}
	    }
	    else
	    {
		u01 = u[index01];
		v01 = v[index01];
		w01 = w[index01];
	    }
	    index10 = d_index3d(i,j-1,k,top_gmax);
	    if ((*findStateAtCrossing)(front,icoords,SOUTH,comp,
			        &intfc_state,&hs,crx_coords)) 
	    {
		u10 = getStateXvel(intfc_state);
		v10 = getStateYvel(intfc_state);
		w10 = getStateZvel(intfc_state);
		if (wave_type(hs) == DIRICHLET_BOUNDARY &&
		    boundary_state_function(hs) &&
		    strcmp(boundary_state_function_name(hs),
                    "flowThroughBoundaryState") == 0)
		{
		    u10 = u0;
		    v10 = v0;
		    w10 = w0;
		}
	    }
	    else
	    {
		u10 = u[index10];
		v10 = v[index10];
		w10 = w[index10];
	    }
	    index11 = d_index3d(i,j+1,k,top_gmax);
	    if ((*findStateAtCrossing)(front,icoords,NORTH,comp,
			        &intfc_state,&hs,crx_coords)) 
	    {
		u11 = getStateXvel(intfc_state);
		v11 = getStateYvel(intfc_state);
		w11 = getStateZvel(intfc_state);
		if (wave_type(hs) == DIRICHLET_BOUNDARY &&
		    boundary_state_function(hs) &&
		    strcmp(boundary_state_function_name(hs),
                    "flowThroughBoundaryState") == 0)
		{
		    u11 = u0;
		    v11 = v0;
		    w11 = w0;
		}
	    }
	    else
	    {
		u11 = u[index11];
		v11 = v[index11];
		w11 = w[index11];
	    }
	    index20 = d_index3d(i,j,k-1,top_gmax);
	    if ((*findStateAtCrossing)(front,icoords,LOWER,comp,
			        &intfc_state,&hs,crx_coords))
	    {
		u20 = getStateXvel(intfc_state);
		v20 = getStateYvel(intfc_state);
		w20 = getStateZvel(intfc_state);
		if (wave_type(hs) == DIRICHLET_BOUNDARY &&
		    boundary_state_function(hs) &&
		    strcmp(boundary_state_function_name(hs),
                    "flowThroughBoundaryState") == 0)
		{
		    u20 = u0;
		    v20 = v0;
		    w20 = w0;
		}
	    }
	    else
	    {
		u20 = u[index20];
		v20 = v[index20];
		w20 = w[index20];
	    }
	    index21 = d_index3d(i,j,k+1,top_gmax);
	    if ((*findStateAtCrossing)(front,icoords,UPPER,comp,
			        &intfc_state,&hs,crx_coords)) 
	    {
		u21 = getStateXvel(intfc_state);
		v21 = getStateYvel(intfc_state);
		w21 = getStateZvel(intfc_state);
		if (wave_type(hs) == DIRICHLET_BOUNDARY &&
		    boundary_state_function(hs) &&
		    strcmp(boundary_state_function_name(hs),
                    "flowThroughBoundaryState") == 0)
		{
		    u21 = u0;
		    v21 = v0;
		    w21 = w0;
		}
	    }
	    else
	    {
		u21 = u[index21];
		v21 = v[index21];
		w21 = w[index21];
	    }

	    cell_center[index].m_state.m_U[0] += -m_dt*(
				burger_flux(u00,u0,u01)/top_h[0] +
				linear_flux(v0,u10,u0,u11)/top_h[1] +
				linear_flux(w0,u20,u0,u21)/top_h[2]);

	    cell_center[index].m_state.m_U[1] += - m_dt*(
				linear_flux(u0,v00,v0,v01)/top_h[0] +
			 	burger_flux(v10,v0,v11)/top_h[1] +
				linear_flux(w0,v20,v0,v21)/top_h[2]);

	    cell_center[index].m_state.m_U[2] += - m_dt*(
				linear_flux(u0,w00,w0,w01)/top_h[0] +
				linear_flux(v0,w10,w0,w11)/top_h[1] +
			 	burger_flux(w20,w0,w21)/top_h[2]);

	    speed = fabs(cell_center[index].m_state.m_U[0]) +
		    fabs(cell_center[index].m_state.m_U[1]) +
		    fabs(cell_center[index].m_state.m_U[2]);
	    if (speed > max_speed)
		max_speed = speed;
	}
	for (l = 0; l < 3; ++l)
	{
	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {	
	    	index  = d_index3d(i,j,k,top_gmax);
	    	array[index] = cell_center[index].m_state.m_U[l];
	    }
	    scatMeshArray();
	    for (k = 0; k <= top_gmax[2]; k++)
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {	
	    	index  = d_index3d(i,j,k,top_gmax);
	    	cell_center[index].m_state.m_U[l] = array[index];
	    }
	}
	pp_global_max(&max_speed,1);
	FT_FreeThese(3,u,v,w);
	if (debugging("trace"))
	    (void) printf("Leaving Incompress_Solver_Smooth_3D_Cartesian::"
			  "computeAdvection()\n");
}	/* end computeAdvection3d */

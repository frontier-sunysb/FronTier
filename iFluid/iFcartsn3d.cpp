/*******************************************************************
 * 			iFcartsn3d.cpp	
 *******************************************************************/
#include "iFluid.h"
#include "solver.h"

//--------------------------------------------------------------------------
// 		   Incompress_Solver_Smooth_3D_Cartesian	
//--------------------------------------------------------------------------

static double (*getStateVel[3])(POINTER) = {getStateXvel,getStateYvel,
                                        getStateZvel};

void Incompress_Solver_Smooth_3D_Cartesian::computeAdvection(void)
{
	int i,j,k,index;
	static HYPERB_SOLVER hyperb_solver(*front);
	double speed;

	static double *rho;
	if (rho == NULL)
	{
	    int size = (top_gmax[0]+1)*(top_gmax[1]+1)*(top_gmax[2]+1);
	    FT_VectorMemoryAlloc((POINTER*)&rho,size,sizeof(double));
	}
	for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    rho[index] = cell_center[index].m_state.m_rho;
	}
	hyperb_solver.rho = rho;

	hyperb_solver.obst_comp = SOLID_COMP;
	switch (iFparams->adv_order)
	{
	case 1:
	    hyperb_solver.order = 1;
	    hyperb_solver.numericalFlux = upwind_flux;
	    break;
	case 4:
	    hyperb_solver.order = 4;
	    hyperb_solver.numericalFlux = weno5_flux;
	    break;
	default:
	    (void) printf("Advection order %d not implemented!\n",
					iFparams->adv_order);
	    clean_up(ERROR);
	}
	hyperb_solver.dt = m_dt;
	hyperb_solver.var = field->vel;
	hyperb_solver.soln = field->vel;
	hyperb_solver.soln_comp1 = LIQUID_COMP1;
	hyperb_solver.soln_comp2 = LIQUID_COMP2;
	hyperb_solver.rho1 = iFparams->rho1;
	hyperb_solver.rho2 = iFparams->rho2;
	hyperb_solver.getStateVel[0] = getStateXvel;
	hyperb_solver.getStateVel[1] = getStateYvel;
	hyperb_solver.getStateVel[2] = getStateZvel;
	hyperb_solver.findStateAtCrossing = findStateAtCrossing;
	hyperb_solver.solveRungeKutta();

	for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    cell_center[index].m_state.m_U[0] = field->vel[0][index];
	    cell_center[index].m_state.m_U[1] = field->vel[1][index];
	    cell_center[index].m_state.m_U[2] = field->vel[2][index];
	    speed = sqrt(sqr(field->vel[0][index]) +
			 sqr(field->vel[1][index]) +
			 sqr(field->vel[2][index]));
	    if (max_speed < speed) max_speed = speed;
	}
}

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

void Incompress_Solver_Smooth_3D_Cartesian::computeNewVelocity(void)
{
	int i, j, k, l, index;
	double grad_phi[3], rho;
	COMPONENT comp;
	double speed;
	int icoords[MAXD];

	max_speed = 0.0;

	for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    array[index] = cell_center[index].m_state.m_phi;
	}
	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    comp = top_comp[index];
	    if (!ifluid_comp(comp))
	    {
		for (l = 0; l < 3; ++l)
		    cell_center[index].m_state.m_U[l] = 0.0;
		continue;
	    }
	    rho = cell_center[index].m_state.m_rho;
	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;

	    computeFieldPointGrad(icoords,array,grad_phi);
	    speed = 0.0;
	    for (l = 0; l < 3; ++l)
	    {
	    	cell_center[index].m_state.m_U[l] -= accum_dt/rho*grad_phi[l];
		speed += fabs(cell_center[index].m_state.m_U[l]);
	    }
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
}	/* end computeNewVelocity3d */


void Incompress_Solver_Smooth_3D_Cartesian::
	computeSourceTerm(double *coords, L_STATE &state) 
{
	int i;
	for (i = 0; i < dim; ++i)
	    state.m_U[i] = iFparams->gravity[i];

	//state.m_P = HUGE_VAL;
}

void Incompress_Solver_Smooth_3D_Cartesian::
	computeSourceTerm(double *coords, double t, L_STATE &state) 
{
	computeSourceTerm(coords, state);
}

// for initial condition: 
// 		setInitialCondition();	
// this function should be called before solve()
// for the source term of the momentum equation: 	
// 		computeSourceTerm();
void Incompress_Solver_Smooth_3D_Cartesian::solve(double dt)
{
	static boolean first = YES;
	if (first)
	{
	    accum_dt = 0.0;
	    first = NO;
	}
	m_dt = dt;
	max_speed = 0.0;

	start_clock("solve");
	setDomain();

	setComponent();
	if (debugging("trace"))
	    printf("Passed setComponent()\n");

	setGlobalIndex();
	if (debugging("trace"))
	    printf("Passed setGlobalIndex()\n");

	start_clock("setSmoothedProperties");
	setSmoothedProperties();
	stop_clock("setSmoothedProperties");
	if (debugging("trace"))
	    printf("Passed setSmoothedProperties()\n");
	
	// 1) solve for intermediate velocity
	start_clock("computeAdvection");
	computeAdvection();
	stop_clock("computeAdvection");
	if (debugging("step_size"))
	    (void) printf("max_speed after computeAdvection(): %20.14f\n",
				max_speed);
	if (debugging("sample_velocity"))
	    sampleVelocity();
	
	start_clock("compDiffWithSmoothProperty");
	compDiffWithSmoothProperty_1st_decoupled();
	stop_clock("compDiffWithSmoothProperty");
	if (debugging("step_size"))
	    (void) printf("max_speed after compDiffusion(): %20.14f\n",
				max_speed);
	if (debugging("sample_velocity"))
	    sampleVelocity();

        start_clock("compSGS");
        //compSGS();	//Subgrid model by Hyunkyun Lim
        stop_clock("compSGS");

	if (debugging("step_size"))
	    printf("max_speed after compDiffWithSmoothProperty(): %20.14f\n",
				max_speed);

	// 2) projection step
	accum_dt += m_dt;
	if (accum_dt >= min_dt)
	{
	    start_clock("computeProjection");
	    computeProjection();
	    stop_clock("computeProjection");

	    start_clock("computePressure");
	    computePressure();
	    stop_clock("computePressure");
	    if (debugging("trace"))
		printf("min_pressure = %f  max_pressure = %f\n",
			min_pressure,max_pressure);

	    start_clock("computeNewVelocity");
	    computeNewVelocity();
	    stop_clock("computeNewVelocity");
	    accum_dt = 0.0;
	}

	if (debugging("step_size"))
	    printf("max_speed after computeNewVelocity(): %20.14f\n",
				max_speed);
	if (debugging("sample_velocity"))
	    sampleVelocity();

	start_clock("copyMeshStates");
	copyMeshStates();
	stop_clock("copyMeshStates");

	setAdvectionDt();
	stop_clock("solve");
}	/* end solve */


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


void Incompress_Solver_Smooth_3D_Cartesian::copyMeshStates()
{
	int i,j,k,d,index;
	double **vel = field->vel;
	double *pres = field->pres;
	double **vort3d = field->vort3d;

	for (i = imin; i <= imax; ++i)
	for (j = jmin; j <= jmax; ++j)
	for (k = kmin; k <= kmax; ++k)
	{
	    index  = d_index3d(i,j,k,top_gmax);
	    if (ifluid_comp(top_comp[index]))
	    {
		pres[index] = cell_center[index].m_state.m_P;
	    	vel[0][index] = cell_center[index].m_state.m_U[0];
	    	vel[1][index] = cell_center[index].m_state.m_U[1];
	    	vel[2][index] = cell_center[index].m_state.m_U[2];
		vort3d[0][index] = getVorticityX(i,j,k);
		vort3d[1][index] = getVorticityY(i,j,k);
		vort3d[2][index] = getVorticityZ(i,j,k);
	    }
	    else
	    {
	    	pres[index] = 0.0;
		for (d = 0; d < 3; ++d)
		{
		    vel[d][index] = 0.0;
		    vort3d[d][index] = 0.0;
		}
	    }
	}
	FT_ParallelExchGridArrayBuffer(pres,front);
	FT_ParallelExchGridArrayBuffer(vel[0],front);
	FT_ParallelExchGridArrayBuffer(vel[1],front);
	FT_ParallelExchGridArrayBuffer(vel[2],front);
	FT_ParallelExchGridArrayBuffer(vort3d[0],front);
	FT_ParallelExchGridArrayBuffer(vort3d[1],front);
	FT_ParallelExchGridArrayBuffer(vort3d[2],front);
}	/* end copyMeshStates */


void Incompress_Solver_Smooth_3D_Cartesian::
	compDiffWithSmoothProperty_1st_decoupled(void)
{
        COMPONENT comp;
        int index,index_nb[6],size;
        int I,I_nb[6];
	int i,j,k,l,nb,icoords[MAXD];
        L_STATE state;
	double coords[MAXD], crx_coords[MAXD];
	double coeff[6],mu[6],mu0,rho,rhs,U_nb[6];
        double speed;
        double *x;
	GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
	POINTER intfc_state;
	HYPER_SURF *hs;
	PetscInt num_iter;
	double rel_residual;
	double aII;

	if (debugging("trace"))
	    (void) printf("Entering Incompress_Solver_Smooth_3D_Cartesian::"
			"compDiffWithSmoothProperty_1st_decoupled()\n");

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

            	for (nb = 0; nb < 6; nb++)
            	{
                    if (FT_StateStructAtGridCrossing(front,icoords,dir[nb],
                                comp,&intfc_state,&hs,crx_coords) &&
                                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
		    {
			if (wave_type(hs) == DIRICHLET_BOUNDARY &&
                    	    boundary_state_function(hs) &&
                    	    strcmp(boundary_state_function_name(hs),
                    	    "flowThroughBoundaryState") == 0)
			{
                    	    U_nb[nb] = cell_center[index].m_state.m_U[l];
			}
			else
			    U_nb[nb] = getStateVel[l](intfc_state);
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

            	coeff[0] = 0.5*m_dt/rho*mu[0]/(top_h[0]*top_h[0]);
            	coeff[1] = 0.5*m_dt/rho*mu[1]/(top_h[0]*top_h[0]);
            	coeff[2] = 0.5*m_dt/rho*mu[2]/(top_h[1]*top_h[1]);
            	coeff[3] = 0.5*m_dt/rho*mu[3]/(top_h[1]*top_h[1]);
            	coeff[4] = 0.5*m_dt/rho*mu[4]/(top_h[2]*top_h[2]);
            	coeff[5] = 0.5*m_dt/rho*mu[5]/(top_h[2]*top_h[2]);

            	getRectangleCenter(index, coords);
            	computeSourceTerm(coords, state);

		aII = 1+coeff[0]+coeff[1]+coeff[2]+coeff[3]+coeff[4]+coeff[5];
		rhs = (1-coeff[0]-coeff[1]-coeff[2]-coeff[3]-coeff[4]-coeff[5])*
		      		cell_center[index].m_state.m_U[l];

		for(nb = 0; nb < 6; nb++)
		{
	            if (!(*findStateAtCrossing)(front,icoords,dir[nb],comp,
			        &intfc_state,&hs,crx_coords))
		    {
			solver.Set_A(I,I_nb[nb],-coeff[nb]);
			rhs += coeff[nb]*U_nb[nb];
		    }
		    else
		    {
			if (wave_type(hs) == DIRICHLET_BOUNDARY &&
                    	    boundary_state_function(hs) &&
                    	    strcmp(boundary_state_function_name(hs),
                    	    "flowThroughBoundaryState") == 0)
			{
			    aII -= coeff[nb];
			    rhs += coeff[nb]*U_nb[nb];
			}
			else
			    rhs += 2.0*coeff[nb]*U_nb[nb];
		    }
		}
		rhs += m_dt*state.m_U[l];
		rhs += m_dt*cell_center[index].m_state.f_surf[l];
		rhs -= m_dt*cell_center[index].m_state.grad_q[l]/rho;
            	solver.Set_A(I,I,aII);

		solver.Set_b(I, rhs);
            }

            solver.SetMaxIter(40000);
            solver.SetTol(1e-14);

	    start_clock("Befor Petsc solve");
            solver.Solve_GMRES();
            solver.GetNumIterations(&num_iter);
            solver.GetFinalRelativeResidualNorm(&rel_residual);

	    stop_clock("After Petsc solve");

            // get back the solution
            solver.Get_x(x);

            if (debugging("PETSc"))
                (void) printf("L_CARTESIAN::"
			"compDiffWithSmoothProperty_decoupled: "
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
        FT_FreeThese(1,x);

	if (debugging("trace"))
	    (void) printf("Leaving Incompress_Solver_Smooth_3D_Cartesian::"
			"compDiffWithSmoothProperty_1st_decoupled()\n");
}       /* end compDiffWithSmoothProperty3d_decoupled */


void Incompress_Solver_Smooth_3D_Cartesian::computePressurePmI(void)
{
        int i,j,k,index;

        for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
	{
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_P += cell_center[index].m_state.m_phi;
	    cell_center[index].m_state.m_q = cell_center[index].m_state.m_P;
	    if (min_pressure > cell_center[index].m_state.m_P)
		min_pressure = cell_center[index].m_state.m_P;
	    if (max_pressure < cell_center[index].m_state.m_P)
		max_pressure = cell_center[index].m_state.m_P;
	}
}        /* end computePressurePmI3d */

void Incompress_Solver_Smooth_3D_Cartesian::computePressurePmII(void)
{
        int i,j,k,index;
        double mu0;

        for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
	{
            index = d_index3d(i,j,k,top_gmax);
            mu0 = 0.5*cell_center[index].m_state.m_mu;
            cell_center[index].m_state.m_P += 
				cell_center[index].m_state.m_phi -
                        	accum_dt*mu0*cell_center[index].m_state.div_U;
	    cell_center[index].m_state.m_q = cell_center[index].m_state.m_P;
	    if (min_pressure > cell_center[index].m_state.m_P)
		min_pressure = cell_center[index].m_state.m_P;
	    if (max_pressure < cell_center[index].m_state.m_P)
		max_pressure = cell_center[index].m_state.m_P;
	}
}        /* end computePressurePmII3d */

void Incompress_Solver_Smooth_3D_Cartesian::computePressurePmIII(void)
{
        int i,j,k,index;
        double mu0;

        for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
	{
            index = d_index3d(i,j,k,top_gmax);
            mu0 = 0.5*cell_center[index].m_state.m_mu;
            cell_center[index].m_state.m_P = 
				cell_center[index].m_state.m_phi -
                        	accum_dt*mu0*cell_center[index].m_state.div_U;
	    cell_center[index].m_state.m_q = 0.0;
	    if (min_pressure > cell_center[index].m_state.m_P)
		min_pressure = cell_center[index].m_state.m_P;
	    if (max_pressure < cell_center[index].m_state.m_P)
		max_pressure = cell_center[index].m_state.m_P;
	}
}        /* end computePressurePmIII3d */

void Incompress_Solver_Smooth_3D_Cartesian::computePressure(void)
{
        int i,j,k,index;
	min_pressure =  HUGE;
	max_pressure = -HUGE;
        for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
	{
            index = d_index3d(i,j,k,top_gmax);
	    if (min_pressure > cell_center[index].m_state.m_P)
		min_pressure = cell_center[index].m_state.m_P;
	    if (max_pressure < cell_center[index].m_state.m_P)
		max_pressure = cell_center[index].m_state.m_P;
	}
	min_pressure =  HUGE;
	max_pressure = -HUGE;
	switch (iFparams->num_scheme.projc_method)
	{
	case BELL_COLELLA:
	    computePressurePmI();
	    break;
	case KIM_MOIN:
	    computePressurePmII();
	    break;
	case SIMPLE:
	case PEROT_BOTELLA:
	    computePressurePmIII();
	    break;
	case ERROR_PROJC_SCHEME:
	default:
	    (void) printf("Unknown computePressure() scheme!\n");
	    clean_up(ERROR);
	}
	computeGradientQ();
}	/* end computePressure */

void Incompress_Solver_Smooth_3D_Cartesian::computeGradientQ(void)
{
	int i,j,k,l,index;
	double *grad_q;
	int icoords[MAXD];

	for (k = 0; k < top_gmax[2]; ++k)
	for (j = 0; j < top_gmax[1]; ++j)
	for (i = 0; i < top_gmax[0]; ++i)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    array[index] = cell_center[index].m_state.m_q;
	}
	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
	for (i = imin; i <= imax; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    grad_q = cell_center[index].m_state.grad_q;
	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;

	    computeFieldPointGrad(icoords,array,grad_q);
	}
	for (l = 0; l < dim; ++l)
	{
	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
	    	index = d_index3d(i,j,k,top_gmax);
		array[index] = cell_center[index].m_state.grad_q[l];
	    }
	    scatMeshArray();
	    for (k = 0; k <= top_gmax[2]; k++)
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {
	    	index = d_index3d(i,j,k,top_gmax);
		cell_center[index].m_state.grad_q[l] = array[index];
	    }
	}
}	/* end computeGradientQ3d */

#define		MAX_TRI_FOR_INTEGRAL		100
void Incompress_Solver_Smooth_3D_Cartesian::surfaceTension(
	double *coords,
	HYPER_SURF_ELEMENT *hse,
	HYPER_SURF *hs,
	double *force,
	double sigma)
{
	int i,j,k,num_tris;
	TRI *tri,*tri_list[MAX_TRI_FOR_INTEGRAL];
	double kappa_tmp,kappa,mag_nor,area,delta;
	double median[MAXD],nor[MAXD];
	POINT *p;

	TriAndFirstRing(hse,hs,&num_tris,tri_list);
	for (i = 0; i < num_tris; ++i)
	{
	    kappa = 0.0;
	    tri = tri_list[i];
	    for (j = 0; j < 3; ++j) median[j] = 0.0;
	    for (j = 0; j < 3; ++j)
	    {
		p = Point_of_tri(tri)[j];
		for (k = 0; k < 3; ++k) 
		    median[k] += Coords(p)[k];
	    	GetFrontCurvature(p,Hyper_surf_element(tri),hs,
				&kappa_tmp,front);
		kappa += kappa_tmp;
		nor[j] = Tri_normal(tri)[j];
	    }
	    kappa /= 3.0;
	    mag_nor = mag_vector(nor,3);
	    area = 0.5*mag_nor;
	    for (j = 0; j < 3; ++j)  
	    {
		nor[j] /= mag_nor;
		median[j] /= 3.0;
	    }
	    delta = smoothedDeltaFunction(coords,median);
	    if (delta == 0.0) continue;
	    for (j = 0; j < dim; ++j) 
	    {
		force[j] += delta*sigma*area*kappa*nor[j];
	    }
	}
}	/* end surfaceTension3d */

void Incompress_Solver_Smooth_3D_Cartesian::setInitialCondition()
{
	int i;
	COMPONENT comp;
	double coords[MAXD];
	int size = (int)cell_center.size();

	FT_MakeGridIntfc(front);
	setDomain();

        m_rho[0] = iFparams->rho1;
        m_rho[1] = iFparams->rho2;
        m_mu[0] = iFparams->mu1;
        m_mu[1] = iFparams->mu2;
	m_comp[0] = iFparams->m_comp1;
	m_comp[1] = iFparams->m_comp2;
	m_smoothing_radius = iFparams->smoothing_radius;
	m_sigma = iFparams->surf_tension;
	mu_min = rho_min = HUGE;
	for (i = 0; i < 2; ++i)
	{
	    if (ifluid_comp(m_comp[i]))
	    {
        	mu_min = std::min(mu_min,m_mu[i]);
        	rho_min = std::min(rho_min,m_rho[i]);
	    }
	}

	// Initialize state at cell_center
        for (i = 0; i < size; i++)
        {
            getRectangleCenter(i, coords);
	    cell_center[i].m_state.setZero();
	    comp = top_comp[i];
	    if (getInitialState != NULL)
	    	(*getInitialState)(comp,coords,cell_center[i].m_state,dim,
						iFparams);
        }
	computeGradientQ();
        copyMeshStates();
	setAdvectionDt();
}       /* end setInitialCondition */

double Incompress_Solver_Smooth_3D_Cartesian::computeFieldPointDiv(
        int *icoords,
        double **field)
{
        int index;
        COMPONENT comp;
        int i, j,k;
	int index_nb[6];
        double div,u_nb[2],v_nb[2],w_nb[2];
        double crx_coords[MAXD];
        POINTER intfc_state;
        HYPER_SURF *hs;
	double h[MAXD];

	i = icoords[0];
	j = icoords[1];
	k = icoords[2];
	index = d_index3d(i,j,k,top_gmax);
        comp = top_comp[index];

	index_nb[0] = d_index3d(i-1,j,k,top_gmax);
	index_nb[1] = d_index3d(i+1,j,k,top_gmax);
	index_nb[2] = d_index3d(i,j-1,k,top_gmax);
	index_nb[3] = d_index3d(i,j+1,k,top_gmax);
	index_nb[4] = d_index3d(i,j,k-1,top_gmax);
	index_nb[5] = d_index3d(i,j,k+1,top_gmax);

	h[0] = 2.0*top_h[0];
	if ((*findStateAtCrossing)(front,icoords,WEST,
		comp,&intfc_state,&hs,crx_coords) &&
	        wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE &&
		wave_type(hs) != ELASTIC_BOUNDARY)
	{
	    if (wave_type(hs) == NEUMANN_BOUNDARY ||
                wave_type(hs) == GROWING_BODY_BOUNDARY ||
                wave_type(hs) == MOVABLE_BODY_BOUNDARY)
            {
                u_nb[0] = getStateXvel(intfc_state);
            }
            else
            {
                u_nb[0] = field[0][index];
                h[0] = top_h[0];
	    }
	}
	else
	    u_nb[0] = field[0][index_nb[0]];

	if ((*findStateAtCrossing)(front,icoords,EAST,
		comp,&intfc_state,&hs,crx_coords) &&
	        wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE &&
		wave_type(hs) != ELASTIC_BOUNDARY)
	{
	    if (wave_type(hs) == NEUMANN_BOUNDARY ||
                wave_type(hs) == GROWING_BODY_BOUNDARY ||
                wave_type(hs) == MOVABLE_BODY_BOUNDARY)
            {
                u_nb[1] = getStateXvel(intfc_state);
            }
            else
            {
                u_nb[1] = field[0][index];
                h[0] = top_h[0];
            }
	}
	else
	    u_nb[1] = field[0][index_nb[1]];

	h[1] = 2.0*top_h[1];
	if ((*findStateAtCrossing)(front,icoords,SOUTH,
		comp,&intfc_state,&hs,crx_coords) &&
	        wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE &&
		wave_type(hs) != ELASTIC_BOUNDARY)
	{
	    if (wave_type(hs) == NEUMANN_BOUNDARY ||
                wave_type(hs) == GROWING_BODY_BOUNDARY ||
                wave_type(hs) == MOVABLE_BODY_BOUNDARY)
            {
                v_nb[0] = getStateYvel(intfc_state);
            }
            else
            {
                v_nb[0] = field[1][index];
                h[1] = top_h[1];
            }
	}
	else
	    v_nb[0] = field[1][index_nb[2]];

	if ((*findStateAtCrossing)(front,icoords,NORTH,
		comp,&intfc_state,&hs,crx_coords) &&
	        wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE &&
		wave_type(hs) != ELASTIC_BOUNDARY)
	{
	    if (wave_type(hs) == NEUMANN_BOUNDARY ||
                wave_type(hs) == GROWING_BODY_BOUNDARY ||
                wave_type(hs) == MOVABLE_BODY_BOUNDARY)
            {
                v_nb[1] = getStateYvel(intfc_state);
            }
            else
            {
                v_nb[1] = field[1][index];
                h[1] = top_h[1];
            }
	}
	else
	    v_nb[1] = field[1][index_nb[3]];

	h[2] = 2.0*top_h[2];
	if ((*findStateAtCrossing)(front,icoords,LOWER,
		comp,&intfc_state,&hs,crx_coords) && 
	        wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE &&
		wave_type(hs) != ELASTIC_BOUNDARY)
	{
	    if (wave_type(hs) == NEUMANN_BOUNDARY ||
                wave_type(hs) == GROWING_BODY_BOUNDARY ||
                wave_type(hs) == MOVABLE_BODY_BOUNDARY)
            {
                w_nb[0] = getStateZvel(intfc_state);
            }
            else
            {
                w_nb[0] = field[2][index];
                h[2] = top_h[2];
            }
	}
	else
	    w_nb[0] = field[2][index_nb[4]];

	if ((*findStateAtCrossing)(front,icoords,UPPER,
		comp,&intfc_state,&hs,crx_coords) &&
	        wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE &&
		wave_type(hs) != ELASTIC_BOUNDARY)
	{
	    if (wave_type(hs) == NEUMANN_BOUNDARY ||
                wave_type(hs) == GROWING_BODY_BOUNDARY ||
                wave_type(hs) == MOVABLE_BODY_BOUNDARY)
            {
                w_nb[1] = getStateZvel(intfc_state);
            }
            else
            {
                w_nb[1] = field[2][index];
                h[2] = top_h[2];
            }
	}
	else
	    w_nb[1] = field[2][index_nb[5]];

        div = (u_nb[1] - u_nb[0])/h[0] + (v_nb[1] - v_nb[0])/h[1] + 
			(w_nb[1] - w_nb[0])/h[2];
        return div;
}       /* end computeFieldPointDiv */

void Incompress_Solver_Smooth_3D_Cartesian::computeFieldPointGrad(
        int *icoords,
        double *field,
        double *grad_field)
{
        int index;
        COMPONENT comp;
        int i,j,k,nb;
	int index_nb[6];
        double p_nbedge[6],p0;  //the p values on the cell edges and cell center
        double crx_coords[MAXD];
        POINTER intfc_state;
        HYPER_SURF *hs;
        GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
	double h[MAXD];

	i = icoords[0];
	j = icoords[1];
	k = icoords[2];
	
	index = d_index3d(i,j,k,top_gmax);
        comp = top_comp[index];
	p0 = field[index];

	index_nb[0] = d_index3d(i-1,j,k,top_gmax);
	index_nb[1] = d_index3d(i+1,j,k,top_gmax);
	index_nb[2] = d_index3d(i,j-1,k,top_gmax);
	index_nb[3] = d_index3d(i,j+1,k,top_gmax);
	index_nb[4] = d_index3d(i,j,k-1,top_gmax);
	index_nb[5] = d_index3d(i,j,k+1,top_gmax);
	for (i = 0; i < 3; ++i)
	    h[i] = 2*top_h[i];

	for (nb = 0; nb < 6; nb++)
	{
	    if(FT_StateStructAtGridCrossing(front,icoords,dir[nb],
			comp,&intfc_state,&hs,crx_coords) &&
		        wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    {
		p_nbedge[nb] = p0; // what is the better way to do this?
		/*
		if (wave_type(hs) == DIRICHLET_BOUNDARY &&
                    boundary_state_function(hs) &&
                    strcmp(boundary_state_function_name(hs),
                    "flowThroughBoundaryState") == 0)
		*/
		if (wave_type(hs) == DIRICHLET_BOUNDARY)
		    h[nb/2] = top_h[nb/2];
	    }
	    else
	    {
		p_nbedge[nb] = field[index_nb[nb]];
	    }
	}
	grad_field[0] = (p_nbedge[1] - p_nbedge[0])/h[0];
	grad_field[1] = (p_nbedge[3] - p_nbedge[2])/h[1];
	grad_field[2] = (p_nbedge[5] - p_nbedge[4])/h[2];
}

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

void Incompress_Solver_Smooth_3D_Cartesian::computeProjectionCim(void)
{
	(void) printf("Not implemented yet!\n");
	clean_up(ERROR);
}	/* end computeProjectionCim */

void Incompress_Solver_Smooth_3D_Cartesian::computeProjection(void)
{
        switch (iFparams->num_scheme.ellip_method)
        {
        case SIMPLE_ELLIP:
            computeProjectionSimple();
            return;
        case CIM_ELLIP:
            computeProjectionCim();
            return;
        }
}       /* end computeProjection */

void Incompress_Solver_Smooth_3D_Cartesian::computeProjectionSimple(void)
{
	static ELLIPTIC_SOLVER elliptic_solver(*front);
        int index;
        int i,j,k,l,icoords[MAXD];
        double **vel = iFparams->field->vel;
        double sum_div;
        double value;
	double min_phi,max_phi;

        for (l = 0; l < dim; ++l)
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index  = d_index3d(i,j,k,top_gmax);
            vel[l][index] = cell_center[index].m_state.m_U[l];
        }

        /* Compute velocity divergence */
        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            icoords[0] = i;
            icoords[1] = j;
            icoords[2] = k;
            index  = d_index(icoords,top_gmax,dim);
            source[index] = computeFieldPointDiv(icoords,vel);
            diff_coeff[index] = 1.0/cell_center[index].m_state.m_rho;
        }
	FT_ParallelExchGridArrayBuffer(source,front);
        FT_ParallelExchGridArrayBuffer(diff_coeff,front);
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index  = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.div_U = source[index];
            source[index] /= accum_dt;
            array[index] = cell_center[index].m_state.m_phi;
        }

        if(debugging("step_size"))
        {
            sum_div = 0.0;
            min_value =  HUGE;
            max_value = -HUGE;
	    min_phi =  HUGE;
	    max_phi = -HUGE;
            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                value = fabs(cell_center[index].m_state.div_U);
                sum_div = sum_div + cell_center[index].m_state.div_U;
            	if (min_phi > cell_center[index].m_state.m_phi)
		    min_phi = cell_center[index].m_state.m_phi;
            	if (max_phi < cell_center[index].m_state.m_phi)
		    max_phi = cell_center[index].m_state.m_phi;
	    	if (min_value > cell_center[index].m_state.div_U) 
		    min_value = cell_center[index].m_state.div_U;
	    	if (max_value < cell_center[index].m_state.div_U) 
		    max_value = cell_center[index].m_state.div_U;
            }
            pp_global_sum(&sum_div,1);
            pp_global_min(&min_value,1);
            pp_global_max(&max_value,1);
	    (void) printf("Before computeProjection:\n");
            (void) printf("\nThe summation of divergence of U is %.16g\n",
                                        sum_div);
            (void) printf("\nThe min value of divergence of U is %.16g\n",
                                        min_value);
            (void) printf("\nThe max value of divergence of U is %.16g\n",
                                        max_value);
	    (void) printf("min_phi = %f  max_phi = %f\n",min_phi,max_phi);
        }
        elliptic_solver.D = diff_coeff;
        elliptic_solver.source = source;
        elliptic_solver.ilower = ilower;
        elliptic_solver.iupper = iupper;
        elliptic_solver.soln = array;
        elliptic_solver.ijk_to_I = ijk_to_I;
        elliptic_solver.set_solver_domain();
        elliptic_solver.getStateVar = getStatePres;
        elliptic_solver.findStateAtCrossing = findStateAtCrossing;
        elliptic_solver.solve(array);

	min_phi =  HUGE;
	max_phi = -HUGE;
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index  = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_phi = array[index];
            if (min_phi > cell_center[index].m_state.m_phi)
		min_phi = cell_center[index].m_state.m_phi;
            if (max_phi < cell_center[index].m_state.m_phi)
		max_phi = cell_center[index].m_state.m_phi;
        }
        if(debugging("step_size"))
	{
	    (void) printf("After computeProjection:\n");
	    (void) printf("min_phi = %f  max_phi = %f\n",min_phi,max_phi);
	}
}	/* end computeProjectionSimple */

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
	    if (max_speed < speed) 
	    {
		max_speed = speed;
		icrds_max[0] = i;
		icrds_max[1] = j;
		icrds_max[2] = k;
	    }
	}
}

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
	    {
		icrds_max[0] = i;
		icrds_max[1] = j;
		icrds_max[2] = k;
		max_speed = speed;
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
}	/* end computeNewVelocity3d */


void Incompress_Solver_Smooth_3D_Cartesian::
	computeSourceTerm(double *coords, L_STATE &state) 
{
	int i;
	for (i = 0; i < dim; ++i)
	    state.m_U[i] = iFparams->gravity[i];
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
	if (debugging("sample_velocity"))
	    sampleVelocity();
	
	
	// 1) solve for intermediate velocity
	start_clock("computeAdvection");
	computeAdvection();
	stop_clock("computeAdvection");
	if (debugging("step_size"))
	{
	    (void) printf("max_speed after computeAdvection(): %20.14f\n",
				max_speed);
	    (void) printf("max speed occured at (%d %d %d)\n",icrds_max[0],
				icrds_max[1],icrds_max[2]);
	}
	if (debugging("sample_velocity"))
	    sampleVelocity();
	
	start_clock("computeDiffusion");
	computeDiffusion();
	stop_clock("computeDiffusion");
	if (debugging("step_size"))
	{
	    (void) printf("max_speed after computeDiffusion(): %20.14f\n",
				max_speed);
	    (void) printf("max speed occured at (%d %d %d)\n",icrds_max[0],
				icrds_max[1],icrds_max[2]);
	}
	if (debugging("sample_velocity"))
	    sampleVelocity();

        start_clock("compSGS");
        //compSGS();	//Subgrid model by Hyunkyun Lim
        stop_clock("compSGS");

	if (debugging("step_size"))
	{
	    (void) printf("max_speed after computeDiffusion(): %20.14f\n",
				max_speed);
	    (void) printf("max speed occured at (%d %d %d)\n",icrds_max[0],
				icrds_max[1],icrds_max[2]);
	}

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
	{
	    printf("max_speed after computeNewVelocity(): %20.14f\n",
				max_speed);
	    (void) printf("max speed occured at (%d %d %d)\n",icrds_max[0],
				icrds_max[1],icrds_max[2]);
	}
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
	computeDiffusion(void)
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
			"computeDiffusion()\n");

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
		//rhs -= m_dt*cell_center[index].m_state.grad_q[l]/rho;
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
			"computeDiffusion: "
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
	    {
		icrds_max[0] = i;
		icrds_max[1] = j;
		icrds_max[2] = k;
                max_speed = speed;
	    }
	}
        pp_global_max(&max_speed,1);
        FT_FreeThese(1,x);

	if (debugging("trace"))
	    (void) printf("Leaving Incompress_Solver_Smooth_3D_Cartesian::"
			"computeDiffusion()\n");
}       /* end computeDiffusion */


void Incompress_Solver_Smooth_3D_Cartesian::computePressurePmI(void)
{
        int i,j,k,index;

	if (debugging("trace"))
	    (void) printf("Entering computePressurePmI()\n");
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
	if (debugging("trace"))
	    (void) printf("Leaving computePressurePmI()\n");
}        /* end computePressurePmI3d */

void Incompress_Solver_Smooth_3D_Cartesian::computePressurePmII(void)
{
        int i,j,k,index;
        double mu0;

	if (debugging("trace"))
	    (void) printf("Entering computePressurePmII()\n");
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
	if (debugging("trace"))
	    (void) printf("Leaving computePressurePmII()\n");
}        /* end computePressurePmII3d */

void Incompress_Solver_Smooth_3D_Cartesian::computePressurePmIII(void)
{
        int i,j,k,index;
        double mu0;

	if (debugging("trace"))
	    (void) printf("Entering computePressurePmIII()\n");
        for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
	{
            index = d_index3d(i,j,k,top_gmax);
            mu0 = 0.5*cell_center[index].m_state.m_mu;
            cell_center[index].m_state.m_P = 
				cell_center[index].m_state.m_phi -
                        	accum_dt*mu0*cell_center[index].m_state.div_U;
	    cell_center[index].m_state.m_q = cell_center[index].m_state.m_P;
	    if (min_pressure > cell_center[index].m_state.m_P)
		min_pressure = cell_center[index].m_state.m_P;
	    if (max_pressure < cell_center[index].m_state.m_P)
		max_pressure = cell_center[index].m_state.m_P;
	}
	if (debugging("trace"))
	    (void) printf("Leaving computePressurePmIII()\n");
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
	    computePressurePmI(); /* computePressurePmIII() not working */
	    //computePressurePmIII();
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
	    {
	    	(*getInitialState)(comp,coords,cell_center[i].m_state,dim,
						iFparams);
		cell_center[i].m_state.m_P = getPressure(front,coords,NULL);
                cell_center[i].m_state.m_phi = getPhiFromPres(front,
                                cell_center[i].m_state.m_P);
	    }
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
                wave_type(hs) == MOVABLE_BODY_BOUNDARY ||
		(wave_type(hs) == DIRICHLET_BOUNDARY && boundary_state(hs)))
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
                wave_type(hs) == MOVABLE_BODY_BOUNDARY ||
		(wave_type(hs) == DIRICHLET_BOUNDARY && boundary_state(hs)))
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
                wave_type(hs) == MOVABLE_BODY_BOUNDARY ||
		(wave_type(hs) == DIRICHLET_BOUNDARY && boundary_state(hs)))
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
                wave_type(hs) == MOVABLE_BODY_BOUNDARY ||
		(wave_type(hs) == DIRICHLET_BOUNDARY && boundary_state(hs)))
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
                wave_type(hs) == MOVABLE_BODY_BOUNDARY ||
		(wave_type(hs) == DIRICHLET_BOUNDARY && boundary_state(hs)))
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
                wave_type(hs) == MOVABLE_BODY_BOUNDARY ||
		(wave_type(hs) == DIRICHLET_BOUNDARY && boundary_state(hs)))
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

	if (debugging("trace"))
	    (void) printf("Entering computeProjectionSimple()\n");
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
        elliptic_solver.getStateVar = getStatePhi;
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
	if (debugging("trace"))
	    (void) printf("Leaving computeProjectionSimple()\n");
}	/* end computeProjectionSimple */

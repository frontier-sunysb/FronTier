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
	    rho[index] = field->rho[index];
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
	double grad_phi[MAXD], rho;
	COMPONENT comp;
	double speed;
	int icoords[MAXD];
	double **vel = field->vel;
	double *phi = field->phi;
	double mag_grad_phi,max_grad_phi,ave_grad_phi;
	int icrds_max[MAXD];
	double vmin[MAXD],vmax[MAXD];

	max_speed = max_grad_phi = ave_grad_phi = 0.0;
	for (i = 0; i < dim; ++i)
        {
            vmin[i] = HUGE;
            vmax[i] = -HUGE;
        }

	for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    array[index] = phi[index];
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
		    vel[l][index] = 0.0;
		continue;
	    }
	    rho = field->rho[index];
	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;

	    computeFieldPointGrad(icoords,array,grad_phi);
	    speed = 0.0;
	    for (l = 0; l < 3; ++l)
	    {
	    	vel[l][index] -= accum_dt/rho*grad_phi[l];
		speed += fabs(vel[l][index]);
	    }
	    mag_grad_phi = Mag3d(grad_phi);
	    ave_grad_phi += mag_grad_phi;
	    if (mag_grad_phi > max_grad_phi)
	    {
		max_grad_phi = mag_grad_phi;	
		icrds_max[0] = i;
		icrds_max[1] = j;
		icrds_max[2] = k;
	    }
	    if (speed > iFparams->ub_speed)
	    {
	    	for (l = 0; l < 3; ++l)
		    vel[l][index] *= iFparams->ub_speed/speed;
		speed = iFparams->ub_speed;
	    }
	    if (speed > max_speed)
	    {
		icrds_max[0] = i;
		icrds_max[1] = j;
		icrds_max[2] = k;
		max_speed = speed;
	    }
	    for (l = 0; l < dim; ++l)
            {
                if (vmin[l] > vel[l][index]) vmin[l] = vel[l][index];
                if (vmax[l] < vel[l][index]) vmax[l] = vel[l][index];
            }
	}
	for (l = 0; l < 3; ++l)
	{
	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
	    {	
	    	index  = d_index3d(i,j,k,top_gmax);
	    	array[index] = vel[l][index];
	    }
	    scatMeshArray();
	    for (k = 0; k <= top_gmax[2]; k++)
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {	
	    	index  = d_index3d(i,j,k,top_gmax);
	    	vel[l][index] = array[index];
	    }
	}
	if (debugging("step_size"))
	{
	    (void) printf("Max gradient phi = %f  occuring at: %d %d %d\n",
			max_grad_phi,icrds_max[0],icrds_max[1],icrds_max[2]);
	    (void) printf("Ave gradient phi = %f\n",ave_grad_phi/(imax-imin+1)
			/(jmax-jmin+1)/(kmax-kmin+1));
	}
	if (debugging("check_div"))
	{
	    checkVelocityDiv("After computeNewVelocity3d()",vmin,vmax);
	}
	pp_global_max(&max_speed,1);
}	/* end computeNewVelocity3d */


void Incompress_Solver_Smooth_3D_Cartesian::
	computeSourceTerm(double *coords, double *source) 
{
	int i;
	for (i = 0; i < dim; ++i)
	    source[i] = iFparams->gravity[i];
} 	/* computeSourceTerm */

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

	int i,j,k,index;
	start_clock("solve");
	setDomain();

	setComponent();
	if (debugging("trace"))
	    printf("Passed setComponent()\n");

	paintAllGridPoint(TO_SOLVE);
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


void Incompress_Solver_Smooth_3D_Cartesian::copyMeshStates()
{
	int i,j,k,d,index;
	double **vel = field->vel;
	double *pres = field->pres;

	for (i = imin; i <= imax; ++i)
	for (j = jmin; j <= jmax; ++j)
	for (k = kmin; k <= kmax; ++k)
	{
	    index  = d_index3d(i,j,k,top_gmax);
	    if (!ifluid_comp(top_comp[index]))
	    {
	    	pres[index] = 0.0;
		for (d = 0; d < 3; ++d)
		{
		    vel[d][index] = 0.0;
		}
	    }
	}
	FT_ParallelExchGridArrayBuffer(pres,front);
	FT_ParallelExchGridArrayBuffer(vel[0],front);
	FT_ParallelExchGridArrayBuffer(vel[1],front);
	FT_ParallelExchGridArrayBuffer(vel[2],front);
}	/* end copyMeshStates */

void Incompress_Solver_Smooth_3D_Cartesian::
	computeDiffusion(void)
{
	//return computeDiffusionCN();
	return computeDiffusionImplicit();
}	/* end computeDiffusion */

void Incompress_Solver_Smooth_3D_Cartesian::
	computeDiffusionCN(void)
{
        COMPONENT comp;
        int index,index_nb[6],size;
        int I,I_nb[6];
	int i,j,k,l,nb,icoords[MAXD];
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
	double source[MAXD];
	double **vel = field->vel;
	double **f_surf = field->f_surf;
	INTERFACE *grid_intfc = front->grid_intfc;
	double vel_min = HUGE;

	if (debugging("trace"))
	    (void) printf("Entering Incompress_Solver_Smooth_3D_Cartesian::"
			"computeDiffusionCN()\n");

        setIndexMap();
	max_speed = 0.0;

        size = iupper - ilower;
        FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));

	for (l = 0; l < dim; ++l)
        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index  = d_index3d(i,j,k,top_gmax);
            if (vel_min > vel[l][index]) vel_min = vel[l][index];
        }
	if (vel_min < 0.0)
	    vel_min *= 5.0;

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


            	mu0   = field->mu[index];
            	rho   = field->rho[index];

            	for (nb = 0; nb < 6; nb++)
            	{
                    if (FT_StateStructAtGridCrossing(front,grid_intfc,icoords,
				dir[nb],comp,&intfc_state,&hs,crx_coords) &&
                                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
		    {
			if (wave_type(hs) == DIRICHLET_BOUNDARY &&
                    	    boundary_state_function(hs) &&
                    	    strcmp(boundary_state_function_name(hs),
                    	    "flowThroughBoundaryState") == 0)
                    	    U_nb[nb] = vel[l][index] - vel_min;
			else
			    U_nb[nb] = getStateVel[l](intfc_state) - vel_min;
			if (wave_type(hs) == DIRICHLET_BOUNDARY || 
			    neumann_type_bdry(wave_type(hs)))
			    mu[nb] = mu0;
			else
			    mu[nb] = 1.0/2*(mu0 + field->mu[index_nb[nb]]);
		    }
                    else
		    {
                    	U_nb[nb] = vel[l][index_nb[nb]] - vel_min;
			mu[nb] = 1.0/2*(mu0 + field->mu[index_nb[nb]]);
		    }
            	}

            	coeff[0] = 0.5*m_dt/rho*mu[0]/(top_h[0]*top_h[0]);
            	coeff[1] = 0.5*m_dt/rho*mu[1]/(top_h[0]*top_h[0]);
            	coeff[2] = 0.5*m_dt/rho*mu[2]/(top_h[1]*top_h[1]);
            	coeff[3] = 0.5*m_dt/rho*mu[3]/(top_h[1]*top_h[1]);
            	coeff[4] = 0.5*m_dt/rho*mu[4]/(top_h[2]*top_h[2]);
            	coeff[5] = 0.5*m_dt/rho*mu[5]/(top_h[2]*top_h[2]);

            	getRectangleCenter(index, coords);
            	computeSourceTerm(coords, source);

		aII = 1+coeff[0]+coeff[1]+coeff[2]+coeff[3]+coeff[4]+coeff[5];
		rhs = (1-coeff[0]-coeff[1]-coeff[2]-coeff[3]-coeff[4]-coeff[5])*
		      		(vel[l][index] - vel_min);

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
		rhs += m_dt*source[l];
		rhs += m_dt*f_surf[l][index];
		//rhs -= m_dt*grad_q[l][index]/rho;
            	solver.Set_A(I,I,aII);
		solver.Set_b(I, rhs);
            }

            solver.SetMaxIter(40000);
            solver.SetTol(1e-14);

	    start_clock("Befor Petsc solve");
            solver.Solve();
            solver.GetNumIterations(&num_iter);
            solver.GetFinalRelativeResidualNorm(&rel_residual);

	    stop_clock("After Petsc solve");

            // get back the solution
            solver.Get_x(x);

            if (debugging("PETSc"))
                (void) printf("L_CARTESIAN::"
			"computeDiffusionCN: "
                        "num_iter = %d, rel_residual = %g. \n",
                        num_iter,rel_residual);

	    for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                I = ijk_to_I[i][j][k];
                index = d_index3d(i,j,k,top_gmax);
                if (I >= 0)
                    array[index] = x[I-ilower] + vel_min;
                else
                    array[index] = 0.0;
            }

            scatMeshArray();
	    for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index3d(i,j,k,top_gmax);
                vel[l][index] = array[index];
            }
        }
	for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
            speed = fabs(vel[0][index]) + fabs(vel[1][index]) +
                    fabs(vel[2][index]);
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
			"computeDiffusionCN()\n");
}       /* end computeDiffusion */

void Incompress_Solver_Smooth_3D_Cartesian::
	computeDiffusionExplicit(void)
{
        COMPONENT comp;
        int index,index_nb[6],size;
        int I,I_nb[6];
	int i,j,k,l,nb,icoords[MAXD];
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
	double source[MAXD];
	double **vel = field->vel;
	double **f_surf = field->f_surf;
	INTERFACE *grid_intfc = front->grid_intfc;

	if (debugging("trace"))
	    (void) printf("Entering Incompress_Solver_Smooth_3D_Cartesian::"
			"computeDiffusionExplicit()\n");

        setIndexMap();
	max_speed = 0.0;

        size = iupper - ilower;
        FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));

	for (l = 0; l < dim; ++l)
	{
	    printf("Before solving direction %d\n",l);
            for (k = kmin; k <= kmax; k++)
	    {
            	index = d_index3d(13,13,k,top_gmax);
		printf("vel[%d][%d] = %20.14f\n",l,k,vel[l][index]);
	    }
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

            	mu0   = field->mu[index];
            	rho   = field->rho[index];

            	for (nb = 0; nb < 6; nb++)
            	{
                    if (FT_StateStructAtGridCrossing(front,grid_intfc,icoords,
				dir[nb],comp,&intfc_state,&hs,crx_coords) &&
                                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
		    {
			if (wave_type(hs) == DIRICHLET_BOUNDARY &&
                    	    boundary_state_function(hs) &&
                    	    strcmp(boundary_state_function_name(hs),
                    	    "flowThroughBoundaryState") == 0)
			{
                    	    U_nb[nb] = vel[l][index];
			}
			else
			{
			    U_nb[nb] = getStateVel[l](intfc_state);
			}
			if (wave_type(hs) == DIRICHLET_BOUNDARY || 
			    neumann_type_bdry(wave_type(hs)))
			    mu[nb] = mu0;
			else
			    mu[nb] = 1.0/2*(mu0 + field->mu[index_nb[nb]]);
		    }
                    else
		    {
                    	U_nb[nb] = vel[l][index_nb[nb]];
			mu[nb] = 1.0/2*(mu0 + field->mu[index_nb[nb]]);
		    }
            	}

            	coeff[0] = m_dt/rho*mu[0]/(top_h[0]*top_h[0]);
            	coeff[1] = m_dt/rho*mu[1]/(top_h[0]*top_h[0]);
            	coeff[2] = m_dt/rho*mu[2]/(top_h[1]*top_h[1]);
            	coeff[3] = m_dt/rho*mu[3]/(top_h[1]*top_h[1]);
            	coeff[4] = m_dt/rho*mu[4]/(top_h[2]*top_h[2]);
            	coeff[5] = m_dt/rho*mu[5]/(top_h[2]*top_h[2]);

            	getRectangleCenter(index, coords);
            	computeSourceTerm(coords, source);

		rhs = (-coeff[0]-coeff[1]-coeff[2]-coeff[3]-coeff[4]-coeff[5])*
		      		vel[l][index];
		if (i == 10 &&  j == 10 && k == 5)
		{
		    printf("step 1: vel[%d][%d] = %20.14f\n",l,k,vel[l][index]);
		    printf("coeffs = %20.14f %20.14f %20.14f\n",
				coeff[0],coeff[1],coeff[2]);
		}

		int num_nb = 0;
		for(nb = 0; nb < 6; nb++)
		{
		    rhs += coeff[nb]*U_nb[nb];
		}
		rhs += m_dt*source[l];
		rhs += m_dt*f_surf[l][index];
                x[I-ilower] = vel[l][index] + rhs;
            }

	    for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                I = ijk_to_I[i][j][k];
                index = d_index3d(i,j,k,top_gmax);
                if (I >= 0)
                    vel[l][index] = x[I-ilower];
                else
                    vel[l][index] = 0.0;
            	speed = fabs(vel[0][index]) + fabs(vel[1][index]) +
                    			fabs(vel[2][index]);
            	if (speed > max_speed)
	    	{
		    icrds_max[0] = i;
		    icrds_max[1] = j;
		    icrds_max[2] = k;
                    max_speed = speed;
	    	}
            }
	    printf("After solving direction %d\n",l);
            for (k = kmin; k <= kmax; k++)
	    {
            	index = d_index3d(13,13,k,top_gmax);
		printf("vel[%d][%d] = %20.14f\n",l,k,vel[l][index]);
	    }
	    index = d_index3d(13,14,8,top_gmax);
	    printf("At 13,14,8: vel = %f %f %f\n",vel[0][index],
				vel[1][index],vel[2][index]);
	    printf("m_dt = %20.14f\n",m_dt);

	    for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index  = d_index3d(i,j,k,top_gmax);
                array[index] = vel[l][index];
            }
            scatMeshArray();
	    for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index3d(i,j,k,top_gmax);
                vel[l][index] = array[index];
            }
        }
	for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
            speed = fabs(vel[0][index]) + fabs(vel[1][index]) +
                    fabs(vel[2][index]);
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
	if (max_speed > 100000.0)
	{
	    printf("In computeDiffusion() max_speed too large: %f\n",
			max_speed);
	    printf("icrds_max = %d %d %d\n",icrds_max[0],icrds_max[1],
			icrds_max[2]);
	    clean_up(0);
	}

	if (debugging("trace"))
	    (void) printf("Leaving Incompress_Solver_Smooth_3D_Cartesian::"
			"computeDiffusionExplicit()\n");
}       /* end computeDiffusionExplicit */

void Incompress_Solver_Smooth_3D_Cartesian::
	computeDiffusionImplicit(void)
{
        COMPONENT comp;
        int index,index_nb[6],size;
        int I,I_nb[6];
	int i,j,k,l,nb,icoords[MAXD];
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
	double source[MAXD];
	double **vel = field->vel;
	double **f_surf = field->f_surf;
	INTERFACE *grid_intfc = front->grid_intfc;
	double vel_min = HUGE;
	double vel_max = -HUGE;

	if (debugging("trace"))
	    (void) printf("Entering Incompress_Solver_Smooth_3D_Cartesian::"
			"computeDiffusionImplicit()\n");

        setIndexMap();
	max_speed = 0.0;

        size = iupper - ilower;
        FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));

	for (l = 0; l < dim; ++l)
        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
            index  = d_index3d(i,j,k,top_gmax);
	    if (vel_min > vel[l][index]) vel_min = vel[l][index];
	    if (vel_max < vel[l][index]) vel_max = vel[l][index];
	}
	if (vel_max - vel_min < 0.1)
	    vel_min -= 0.01;
	else
            vel_min -= 0.1*(vel_max - vel_min);

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


            	mu0   = field->mu[index];
            	rho   = field->rho[index];

            	for (nb = 0; nb < 6; nb++)
            	{
                    if (FT_StateStructAtGridCrossing(front,grid_intfc,icoords,
				dir[nb],comp,&intfc_state,&hs,crx_coords) &&
                                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
		    {
			if (wave_type(hs) == DIRICHLET_BOUNDARY &&
                    	    boundary_state_function(hs) &&
                    	    strcmp(boundary_state_function_name(hs),
                    	    "flowThroughBoundaryState") == 0)
			{
                    	    U_nb[nb] = vel[l][index] - vel_min;
			}
			else
			{
			    U_nb[nb] = getStateVel[l](intfc_state) - vel_min;
			}
			if (wave_type(hs) == DIRICHLET_BOUNDARY || 
			    neumann_type_bdry(wave_type(hs)))
			    mu[nb] = mu0;
			else
			    mu[nb] = 1.0/2*(mu0 + field->mu[index_nb[nb]]);
		    }
                    else
		    {
                    	U_nb[nb] = vel[l][index_nb[nb]] - vel_min;
			mu[nb] = 1.0/2*(mu0 + field->mu[index_nb[nb]]);
		    }
            	}

            	coeff[0] = m_dt/rho*mu[0]/(top_h[0]*top_h[0]);
            	coeff[1] = m_dt/rho*mu[1]/(top_h[0]*top_h[0]);
            	coeff[2] = m_dt/rho*mu[2]/(top_h[1]*top_h[1]);
            	coeff[3] = m_dt/rho*mu[3]/(top_h[1]*top_h[1]);
            	coeff[4] = m_dt/rho*mu[4]/(top_h[2]*top_h[2]);
            	coeff[5] = m_dt/rho*mu[5]/(top_h[2]*top_h[2]);

            	getRectangleCenter(index, coords);
            	computeSourceTerm(coords, source);

		aII = 1+coeff[0]+coeff[1]+coeff[2]+coeff[3]+coeff[4]+coeff[5];
		rhs = vel[l][index] - vel_min;

		for(nb = 0; nb < 6; nb++)
		{
	            if (!(*findStateAtCrossing)(front,icoords,dir[nb],comp,
			        &intfc_state,&hs,crx_coords))
			solver.Set_A(I,I_nb[nb],-coeff[nb]);
		    else
		    {
			if (wave_type(hs) == DIRICHLET_BOUNDARY &&
                    	    boundary_state_function(hs) &&
                    	    strcmp(boundary_state_function_name(hs),
                    	    "flowThroughBoundaryState") == 0)
			    aII -= coeff[nb];
			else
			    rhs += coeff[nb]*U_nb[nb];
		    }
		}
		rhs += m_dt*source[l];
		rhs += m_dt*f_surf[l][index];
		//rhs -= m_dt*grad_q[l][index]/rho;
            	solver.Set_A(I,I,aII);

		solver.Set_b(I, rhs);
            }

            solver.SetMaxIter(40000);
            solver.SetTol(1e-14);

	    start_clock("Befor Petsc solve");
            solver.Solve();
            solver.GetNumIterations(&num_iter);
            solver.GetFinalRelativeResidualNorm(&rel_residual);

	    stop_clock("After Petsc solve");

            // get back the solution
            solver.Get_x(x);

            if (debugging("PETSc"))
                (void) printf("L_CARTESIAN::"
			"computeDiffusionImplicit: "
                        "num_iter = %d, rel_residual = %g. \n",
                        num_iter,rel_residual);

	    for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                I = ijk_to_I[i][j][k];
                index = d_index3d(i,j,k,top_gmax);
                if (I >= 0)
                    vel[l][index] = x[I-ilower] + vel_min;
                else
                    vel[l][index] = 0.0;
            	speed = fabs(vel[0][index]) + fabs(vel[1][index]) +
                    			fabs(vel[2][index]);
            	if (speed > max_speed)
	    	{
		    icrds_max[0] = i;
		    icrds_max[1] = j;
		    icrds_max[2] = k;
                    max_speed = speed;
	    	}
            }
	    printf("After solving direction %d\n",l);
            for (k = kmin; k <= kmax; k++)
	    {
            	index = d_index3d(13,13,k,top_gmax);
		printf("vel[%d][%d] = %20.14f\n",l,k,vel[l][index]);
	    }
	    index = d_index3d(13,14,8,top_gmax);

	    for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index  = d_index3d(i,j,k,top_gmax);
                array[index] = vel[l][index];
            }
            scatMeshArray();
	    for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index3d(i,j,k,top_gmax);
                vel[l][index] = array[index];
            }
        }
	for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
            speed = fabs(vel[0][index]) + fabs(vel[1][index]) +
                    fabs(vel[2][index]);
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
			"computeDiffusionImplicit()\n");
}       /* end computeDiffusionImplicit */

void Incompress_Solver_Smooth_3D_Cartesian::computePressurePmI(void)
{
        int i,j,k,index;
	double *pres = field->pres;
	double *phi = field->phi;
	double *q = field->q;
	int icrds_max[MAXD],icrds_min[MAXD];

	if (debugging("trace"))
	    (void) printf("Entering computePressurePmI()\n");
        for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
	{
            index = d_index3d(i,j,k,top_gmax);
            //pres[index] += phi[index];
            pres[index] = phi[index];
	    q[index] = pres[index];
	    if (i < imin || i > imax || j < jmin || j > jmax || 
		k < kmin || k > kmax)
		continue;
	    if (min_pressure > pres[index])
	    {
		min_pressure = pres[index];
		icrds_min[0] = i;
		icrds_min[1] = j;
		icrds_min[2] = k;
	    }
	    if (max_pressure < pres[index])
	    {
		max_pressure = pres[index];
		icrds_max[0] = i;
		icrds_max[1] = j;
		icrds_max[2] = k;
	    }
	}
	if (debugging("step_size"))
	{
	    (void) printf(" Max pressure = %f  occuring at: %d %d %d\n",
			max_pressure,icrds_max[0],icrds_max[1],icrds_max[2]);
	    (void) printf(" Min pressure = %f  occuring at: %d %d %d\n",
			min_pressure,icrds_min[0],icrds_min[1],icrds_min[2]);
	    (void) printf("Diff pressure = %f\n",max_pressure-min_pressure);
	}
	if (debugging("trace"))
	    (void) printf("Leaving computePressurePmI()\n");
}        /* end computePressurePmI3d */

void Incompress_Solver_Smooth_3D_Cartesian::computePressurePmII(void)
{
        int i,j,k,index;
        double mu0;
	double *pres = field->pres;
	double *phi = field->phi;
	double *q = field->q;
	double *div_U = field->div_U;

	if (debugging("trace"))
	    (void) printf("Entering computePressurePmII()\n");
        for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
	{
            index = d_index3d(i,j,k,top_gmax);
            mu0 = 0.5*field->mu[index];
            pres[index] += phi[index] - accum_dt*mu0*div_U[index];
	    q[index] = pres[index];
	    if (min_pressure > pres[index])
		min_pressure = pres[index];
	    if (max_pressure < pres[index])
		max_pressure = pres[index];
	}
	if (debugging("trace"))
	    (void) printf("Leaving computePressurePmII()\n");
}        /* end computePressurePmII3d */

void Incompress_Solver_Smooth_3D_Cartesian::computePressurePmIII(void)
{
        int i,j,k,index;
        double mu0;
	double *pres = field->pres;
	double *phi = field->phi;
	double *q = field->q;
	double *div_U = field->div_U;

	if (debugging("trace"))
	    (void) printf("Entering computePressurePmIII()\n");
        for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
	{
            index = d_index3d(i,j,k,top_gmax);
            mu0 = 0.5*field->mu[index];
            pres[index] = phi[index] - accum_dt*mu0*div_U[index];
	    q[index] = pres[index];
	    if (min_pressure > pres[index])
		min_pressure = pres[index];
	    if (max_pressure < pres[index])
		max_pressure = pres[index];
	}
	if (debugging("trace"))
	    (void) printf("Leaving computePressurePmIII()\n");
}        /* end computePressurePmIII3d */

void Incompress_Solver_Smooth_3D_Cartesian::computePressure(void)
{
        int i,j,k,index;
	double *pres = field->pres;
	min_pressure =  HUGE;
	max_pressure = -HUGE;
        for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
	{
            index = d_index3d(i,j,k,top_gmax);
	    if (min_pressure > pres[index])
		min_pressure = pres[index];
	    if (max_pressure < pres[index])
		max_pressure = pres[index];
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
	double **grad_q = field->grad_q;
	int icoords[MAXD];
	double *q = field->q;
	double point_grad_q[MAXD];

	for (k = 0; k < top_gmax[2]; ++k)
	for (j = 0; j < top_gmax[1]; ++j)
	for (i = 0; i < top_gmax[0]; ++i)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    array[index] = q[index];
	}
	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
	for (i = imin; i <= imax; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;

	    computeFieldPointGrad(icoords,array,point_grad_q);
	    for (l = 0; l < 3; ++l)
		grad_q[l][index] = point_grad_q[l];
	}
	for (l = 0; l < dim; ++l)
	{
	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
	    	index = d_index3d(i,j,k,top_gmax);
		array[index] = grad_q[l][index];
	    }
	    scatMeshArray();
	    for (k = 0; k <= top_gmax[2]; k++)
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {
	    	index = d_index3d(i,j,k,top_gmax);
		grad_q[l][index] = array[index];
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
	double *pres = field->pres;
	double *phi = field->phi;

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
	    //cell_center[i].m_state.setZero();
	    comp = top_comp[i];
	    if (getInitialState != NULL)
	    {
	    	(*getInitialState)(comp,coords,field,i,dim,iFparams);
		pres[i] = getPressure(front,coords,NULL);
                phi[i] = getPhiFromPres(front,pres[i]);
	    }
        }
	computeGradientQ();
        copyMeshStates();
	setAdvectionDt();
}       /* end setInitialCondition */

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
        case DOUBLE_ELLIP:
            computeProjectionDouble();
            return;
        case DUAL_ELLIP:
            computeProjectionDual();
            return;
        case CIM_ELLIP:
            computeProjectionCim();
            return;
        }
}       /* end computeProjection */

void Incompress_Solver_Smooth_3D_Cartesian::computeProjectionDouble(void)
{
}	/* end computeProjectionDouble */

void Incompress_Solver_Smooth_3D_Cartesian::computeProjectionDual(void)
{
}	/* end computeProjectionDual */

void Incompress_Solver_Smooth_3D_Cartesian::computeProjectionSimple(void)
{
	static ELLIPTIC_SOLVER elliptic_solver(*front);
        int index;
        int i,j,k,l,icoords[MAXD];
        double **vel = field->vel;
        double *phi = field->phi;
        double *div_U = field->div_U;
        double sum_div,L1_div;
        double value;
	double min_phi,max_phi;
	int icrds_max[MAXD],icrds_min[MAXD];
	double vmin[MAXD],vmax[MAXD];

	if (debugging("trace"))
	    (void) printf("Entering computeProjectionSimple()\n");

	for (l = 0; l < dim; ++l)
        {
            vmin[l] = HUGE;
            vmax[l] = -HUGE;
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
            diff_coeff[index] = 1.0/field->rho[index];
	    if (debugging("check_div"))
            {
                for (l = 0; l < dim; ++l)
                {
                    if (vmin[l] > field->vel[l][index])
                        vmin[l] = field->vel[l][index];
                    if (vmax[l] < field->vel[l][index])
                        vmax[l] = field->vel[l][index];
                }
            }
        }
	FT_ParallelExchGridArrayBuffer(source,front);
        FT_ParallelExchGridArrayBuffer(diff_coeff,front);
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index  = d_index3d(i,j,k,top_gmax);
            div_U[index] = source[index];
            source[index] /= accum_dt;
            array[index] = phi[index];
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
                value = fabs(div_U[index]);
                sum_div = sum_div + fabs(div_U[index]);
	    	if (min_value > div_U[index]) 
		{
		    min_value = div_U[index];
		    icrds_min[0] = i;
		    icrds_min[1] = j;
		    icrds_min[2] = k;
		}
	    	if (max_value < div_U[index]) 
		{
		    max_value = div_U[index];
		    icrds_max[0] = i;
		    icrds_max[1] = j;
		    icrds_max[2] = k;
		}
            }
            pp_global_min(&min_value,1);
            pp_global_max(&max_value,1);
	    L1_div = sum_div/(imax-imin+1)/(jmax-jmin+1)/(kmax-kmin+1);
	    (void) printf("Before computeProjection:\n");
            (void) printf("Sum div(U) = %f\n",sum_div);
            (void) printf("Min div(U) = %f  ",min_value);
	    (void) printf("occuring at: %d %d %d\n",icrds_min[0],icrds_min[1],
					icrds_min[2]);
            (void) printf("Max div(U) = %f  ",max_value);
	    (void) printf("occuring at: %d %d %d\n",icrds_max[0],icrds_max[1],
					icrds_max[2]);
            (void) printf("L1  div(U) = %f\n",L1_div);
        }
	if (debugging("check_div"))
        {
            checkVelocityDiv("Before computeProjection()",vmin,vmax);
        }
        elliptic_solver.D = diff_coeff;
        elliptic_solver.source = source;
        elliptic_solver.soln = array;
        elliptic_solver.set_solver_domain();
        elliptic_solver.getStateVar = getStatePhi;
        elliptic_solver.findStateAtCrossing = findStateAtCrossing;
	paintAllGridPoint(NOT_SOLVED);
	while (paintToSolveGridPoint())
	{
	    setGlobalIndex();
	    setIndexMap();
            elliptic_solver.ijk_to_I = ijk_to_I;
            elliptic_solver.ilower = ilower;
            elliptic_solver.iupper = iupper;
	    if (iFparams->total_div_cancellation)
            	elliptic_solver.dsolve(array);
	    else
            	elliptic_solver.solve(array);
	    paintSolvedGridPoint();
	}


	min_phi =  HUGE;
	max_phi = -HUGE;
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index  = d_index3d(i,j,k,top_gmax);
            phi[index] = array[index];
            if (min_phi > phi[index])
	    {
		min_phi = phi[index];
		icrds_min[0] = i;
		icrds_min[1] = j;
		icrds_min[2] = k;
	    }
            if (max_phi < phi[index])
	    {
		max_phi = phi[index];
		icrds_max[0] = i;
		icrds_max[1] = j;
		icrds_max[2] = k;
	    }
        }
        if(debugging("step_size"))
	{
	    (void) printf("After computeProjection:\n");
	    (void) printf("min_phi = %f  occuring at: %d %d %d\n",
			min_phi,icrds_min[0],icrds_min[1],icrds_min[2]);
	    (void) printf("max_phi = %f  occuring at: %d %d %d\n",
			max_phi,icrds_max[0],icrds_max[1],icrds_max[2]);
	    (void) printf("diff_phi = %f\n",max_phi-min_phi);
	}
	if (debugging("trace"))
	    (void) printf("Leaving computeProjectionSimple()\n");
}	/* end computeProjectionSimple */

void Incompress_Solver_Smooth_3D_Cartesian::solveTest(const char *msg)
{
	// This function is reserved for various debugging tests.
	int k,index;
	int icoords[MAXD];
	double **vel = field->vel;
	(void) printf("%s\n",msg);
        for (k = kmin; k <= kmax; k++)
        {
            icoords[0] = 30;
            icoords[1] = 30;
            icoords[2] = k;
            index  = d_index(icoords,top_gmax,dim);
            source[index] = computeFieldPointDiv(icoords,vel);
            source[index] /= m_dt;
	    printf("div[%d]/dt = %14.8f\n",k,source[index]);
        }
}	/* end solveTest */

static int parab_find_state_at_crossing(Front*,int*,GRID_DIRECTION,int,
                                POINTER*,HYPER_SURF**,double*);

void Incompress_Solver_Smooth_3D_Cartesian::
	computeDiffusionParab(void)
{
	static PARABOLIC_SOLVER parab_solver(*front);

        COMPONENT comp;
        int index;
	int i,j,k,l,icoords[MAXD];
	double source[MAXD];
	double **vel = field->vel;
	double **f_surf = field->f_surf;
	static double *nu;

	if (debugging("trace"))
	    (void) printf("Entering Incompress_Solver_Smooth_3D_Cartesian::"
			"computeDiffusionParab()\n");

	if (nu == NULL)
	{
	    int size = 1;
	    for (i = 0; i < dim; ++i)  size *= (top_gmax[i] + 1);
	    FT_VectorMemoryAlloc((POINTER*)&nu,size,sizeof(double));
	}
        setIndexMap();
        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index  = d_index3d(i,j,k,top_gmax);
	    nu[index] = field->mu[index]/field->rho[index];
	}
	FT_ParallelExchGridArrayBuffer(nu,front);

	parab_solver.soln_comp = LIQUID_COMP2;
	parab_solver.obst_comp = SOLID_COMP;
	parab_solver.ilower = ilower;
        parab_solver.iupper = iupper;
	parab_solver.dt = m_dt;
	parab_solver.order = 2;
	parab_solver.a = NULL;
	parab_solver.findStateAtCrossing = parab_find_state_at_crossing;
	parab_solver.set_solver_domain();
	switch(dim)
        {
        case 2:
            parab_solver.ij_to_I = ij_to_I;
            break;
        case 3:
            parab_solver.ijk_to_I = ijk_to_I;
            break;
        }
	
	for (l = 0; l < dim; ++l)
	{
	    parab_solver.var = vel[l];
	    parab_solver.soln = vel[l];
	    parab_solver.getStateVarFunc = getStateVel[l];
	    parab_solver.runge_kutta();
	}

	if (debugging("trace"))
	    (void) printf("Leaving Incompress_Solver_Smooth_3D_Cartesian::"
			"computeDiffusionParab()\n");
}       /* end computeDiffusionParab */

static int parab_find_state_at_crossing(
	Front *front,
        int *icoords,
        GRID_DIRECTION dir,
        int comp,
        POINTER *state,
        HYPER_SURF **hs,
        double *crx_coords)
{
	boolean status;
        INTERFACE *grid_intfc = front->grid_intfc;

        status = FT_StateStructAtGridCrossing(front,grid_intfc,icoords,dir,
                                comp,state,hs,crx_coords);
        if (status == NO) return NO_PDE_BOUNDARY;

	switch (wave_type(*hs))
	{
	case FIRST_PHYSICS_WAVE_TYPE:
	    return NO_PDE_BOUNDARY;
	case DIRICHLET_BOUNDARY:
	case GROWING_BODY_BOUNDARY:
	case MOVABLE_BODY_BOUNDARY:
	case NEUMANN_BOUNDARY:
	case ICE_PARTICLE_BOUNDARY:
            return DIRICHLET_PDE_BOUNDARY;
	default:
	    (void) printf("In parab_find_state_at_crossing()\n");
	    (void) printf("Unknown wave type %s\n",
				f_wave_type_as_string(wave_type(*hs)));
	    clean_up(ERROR);
        }
}	/* end parab_find_state_at_crossing */

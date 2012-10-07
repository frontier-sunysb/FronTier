/*******************************************************************
 * 	               iFcartsn2d.cpp	
 *******************************************************************/
#include "iFluid.h"
#include "solver.h"

static double (*getStateVel[3])(POINTER) = {getStateXvel,getStateYvel,
                                        getStateZvel};


void Incompress_Solver_Smooth_2D_Cartesian::computeAdvection(void)
{
	int i,j,l,index;
	static HYPERB_SOLVER hyperb_solver(*front);

	static double *rho;
	if (rho == NULL)
	{
	    int size = (top_gmax[0]+1)*(top_gmax[1]+1);
	    FT_VectorMemoryAlloc((POINTER*)&rho,size,sizeof(double));
	}
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index = d_index2d(i,j,top_gmax);
	    rho[index] = cell_center[index].m_state.m_rho;
	    for (l = 0; l < dim; ++l)
                field->vel[l][index] = cell_center[index].m_state.m_U[l];
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
	hyperb_solver.getStateVel[0] = getStateXvel;
	hyperb_solver.getStateVel[1] = getStateYvel;
	hyperb_solver.rho1 = iFparams->rho1;
	hyperb_solver.rho2 = iFparams->rho2;
	hyperb_solver.findStateAtCrossing = ifluid_find_state_at_crossing;
	hyperb_solver.solveRungeKutta();

	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index = d_index2d(i,j,top_gmax);
	    cell_center[index].m_state.m_U[0] = field->vel[0][index];
	    cell_center[index].m_state.m_U[1] = field->vel[1][index];
	}
}

void Incompress_Solver_Smooth_2D_Cartesian::computeProjection(void)
{
	setIndexMap();
	switch (iFparams->num_scheme.ellip_method)
	{
	case SIMPLE_ELLIP:
	    computeProjectionSimple();
	    return;
	case CIM_ELLIP:
	    computeProjectionCim();
	    return;
	}
}	/* end computeProjection */

void Incompress_Solver_Smooth_2D_Cartesian::computeProjectionCim(void)
{
	static CIM_ELLIPTIC_SOLVER elliptic_solver(*front);
	int index;
	int i,j,l,icoords[MAXD];
	double **vel = iFparams->field->vel;
	double sum_div;
	double value;

	sum_div = 0.0;
	max_value = 0.0;

	for (l = 0; l < dim; ++l)
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index  = d_index2d(i,j,top_gmax);
	    vel[l][index] = cell_center[index].m_state.m_U[l];
	}

	/* Compute velocity divergence */
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    icoords[0] = i;
	    icoords[1] = j;
	    index  = d_index(icoords,top_gmax,dim);
	    source[index] = computeFieldPointDiv(icoords,vel);
	    diff_coeff[index] = 1.0/cell_center[index].m_state.m_rho;
	}
	FT_ParallelExchGridArrayBuffer(source,front);
	FT_ParallelExchGridArrayBuffer(diff_coeff,front);
	for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
	{
	    index  = d_index2d(i,j,top_gmax);
	    cell_center[index].m_state.div_U = source[index];
	    source[index] /= -accum_dt;
	    array[index] = cell_center[index].m_state.m_phi;
	}

	if(debugging("step_size"))
	{
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		index = d_index2d(i,j,top_gmax);
	        value = fabs(cell_center[index].m_state.div_U);
		sum_div = sum_div + cell_center[index].m_state.div_U;
	    }
	    pp_global_sum(&sum_div,1);
	    (void) printf("\nThe summation of divergence of U is %.16g\n",
					sum_div);
	    pp_global_max(&max_value,1);
	    (void) printf("\nThe max value of divergence of U is %.16g\n",
					max_value);
	    max_value = 0.0;
	}
	elliptic_solver.w_type = FIRST_PHYSICS_WAVE_TYPE;
	elliptic_solver.neg_comp = LIQUID_COMP2;
	elliptic_solver.pos_comp = LIQUID_COMP1;
	elliptic_solver.solutionJump = p_jump;
        elliptic_solver.gradJumpDotN = grad_p_jump_n;
        elliptic_solver.gradJumpDotT = grad_p_jump_t;

        elliptic_solver.diff_coeff[0] = 1.0/m_rho[0];
        elliptic_solver.diff_coeff[1] = 1.0/m_rho[1];
        elliptic_solver.kappa = diff_coeff;
        elliptic_solver.source = source;
        elliptic_solver.ilower = ilower;
        elliptic_solver.iupper = iupper;
        elliptic_solver.soln = array;
        elliptic_solver.ij_to_I = ij_to_I;
        elliptic_solver.I_to_ij = I_to_ij;
	elliptic_solver.getStateVar = getStatePhi;
	elliptic_solver.findStateAtCrossing = ifluid_find_state_at_crossing;
	elliptic_solver.size = iupper - ilower;
	elliptic_solver.set_solver_domain();
	elliptic_solver.solve(array);
	//viewTopVariable(front,array,NO,0.0,0.0,(char*)"test-cim",
	//			(char*)"array");

	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index  = d_index2d(i,j,top_gmax);
	    cell_center[index].m_state.m_phi = array[index];
	}
}	/* end computeProjectionCim */

void Incompress_Solver_Smooth_2D_Cartesian::computeProjectionSimple(void)
{
	static ELLIPTIC_SOLVER elliptic_solver(*front);
	int index;
	int i,j,l,icoords[MAXD];
	double **vel = iFparams->field->vel;
	double sum_div;
	double value;

	sum_div = 0.0;
	max_value = 0.0;

	for (l = 0; l < dim; ++l)
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index  = d_index2d(i,j,top_gmax);
	    if (!ifluid_comp(top_comp[index]))
		continue;
	    vel[l][index] = cell_center[index].m_state.m_U[l];
	}

	/* Compute velocity divergence */
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    icoords[0] = i;
	    icoords[1] = j;
	    index  = d_index(icoords,top_gmax,dim);
	    if (!ifluid_comp(top_comp[index]))
		continue;
	    source[index] = computeFieldPointDiv(icoords,vel);
	    diff_coeff[index] = 1.0/cell_center[index].m_state.m_rho;
	}
	FT_ParallelExchGridArrayBuffer(source,front);
	FT_ParallelExchGridArrayBuffer(diff_coeff,front);
	for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
	{
	    index  = d_index2d(i,j,top_gmax);
	    if (!ifluid_comp(top_comp[index]))
		continue;
	    cell_center[index].m_state.div_U = source[index];
	    source[index] /= accum_dt;
	    array[index] = cell_center[index].m_state.m_phi;
	}

	if(debugging("step_size"))
	{
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		index = d_index2d(i,j,top_gmax);
	        if (!ifluid_comp(top_comp[index]))
		    continue;
	        value = fabs(cell_center[index].m_state.div_U);
		sum_div = sum_div + cell_center[index].m_state.div_U;
	    }
	    pp_global_sum(&sum_div,1);
	    (void) printf("\nThe summation of divergence of U is %.16g\n",
					sum_div);
	    pp_global_max(&max_value,1);
	    (void) printf("\nThe max value of divergence of U is %.16g\n",
					max_value);
	    max_value = 0.0;
	}
        elliptic_solver.D = diff_coeff;
        elliptic_solver.source = source;
        elliptic_solver.ilower = ilower;
        elliptic_solver.iupper = iupper;
        elliptic_solver.soln = array;
        elliptic_solver.ij_to_I = ij_to_I;
	elliptic_solver.set_solver_domain();
	elliptic_solver.getStateVar = getStatePres;
	elliptic_solver.findStateAtCrossing = ifluid_find_projection_crossing;
	elliptic_solver.solve(array);
	//viewTopVariable(front,array,NO,0.0,0.0,(char*)"test-simple",
	//			(char*)"array");

	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index  = d_index2d(i,j,top_gmax);
	    cell_center[index].m_state.m_phi = array[index];
	}
}	/* end computeProjectionSimple */

void Incompress_Solver_Smooth_2D_Cartesian::computeNewVelocity(void)
{
	int i, j, k, index;
	double grad_phi[2], rho;
	COMPONENT comp;
	double speed;
	int icoords[MAXD];

	max_speed = 0.0;

	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index  = d_index2d(i,j,top_gmax);
	    array[index] = cell_center[index].m_state.m_phi;
	}
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    index = d_index2d(i,j,top_gmax);
	    comp = top_comp[index];
	    if (!ifluid_comp(comp))
	    {
		cell_center[index].m_state.m_U[0] = 0.0;
		cell_center[index].m_state.m_U[1] = 0.0;
		continue;
	    }
	    rho = cell_center[index].m_state.m_rho;
	    icoords[0] = i;
	    icoords[1] = j;
	    computeFieldPointGrad(icoords,array,grad_phi);
	    cell_center[index].m_state.m_U[0] -= accum_dt/rho*grad_phi[0];
	    cell_center[index].m_state.m_U[1] -= accum_dt/rho*grad_phi[1];
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
}	/* end computeNewVelocity2d */


void Incompress_Solver_Smooth_2D_Cartesian::computeSourceTerm(double *coords, L_STATE &state) 
{
	int i;
	for (i = 0; i < dim; ++i)
	    state.m_U[i] = iFparams->gravity[i];

	state.m_P = HUGE_VAL;
}
void Incompress_Solver_Smooth_2D_Cartesian::computeSourceTerm(double *coords, double t, L_STATE &state) 
{
	computeSourceTerm(coords, state);
}


// for initial condition: 
// 		setInitialCondition();	
// this function should be called before solve()
// for the source term of the momentum equation: 	
// 		computeSourceTerm();
void Incompress_Solver_Smooth_2D_Cartesian::solve(double dt)
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
	    printf("max_speed after computeAdvection(): %20.14f\n",max_speed);
	if (debugging("sample_velocity"))
	    sampleVelocity();
	
	start_clock("computeDiffusion");
	computeDiffusion();
	stop_clock("computeDiffusion");
	if (debugging("sample_velocity"))
	    sampleVelocity();

        start_clock("compSGS");
        //compSGS();	//Subgrid model by Hyunkyun Lim
        stop_clock("compSGS");

	if (debugging("step_size"))
	    printf("max_speed after computeDiffusion(): %20.14f\n",
				max_speed);

	// 2) projection step
	accum_dt += m_dt;
	printf("min_dt = %f  accum_dt = %f\n",min_dt,accum_dt);
	if (accum_dt >= min_dt)
	{
	    start_clock("computeProjection");
	    computeProjection();
	    stop_clock("computeProjection");

	    start_clock("computePressure");
	    computePressure();
	    stop_clock("computePressure");

	    start_clock("computeNewVelocity");
	    computeNewVelocity();
	    stop_clock("computeNewVelocity");
	    accum_dt = 0.0;
	}

	if (debugging("sample_velocity"))
	    sampleVelocity();

	if (debugging("step_size"))
	    printf("max_speed after computeNewVelocity(): %20.14f\n",
				max_speed);

	start_clock("copyMeshStates");
	copyMeshStates();
	stop_clock("copyMeshStates");

	setAdvectionDt();
	stop_clock("solve");
}	/* end solve */


double Incompress_Solver_Smooth_2D_Cartesian::getVorticity(int i, int j)
{
	int index0,index00,index01,index10,index11;
	double v00,v01,v10,v11;
	double dx,dy;
	double vorticity;

	dx = top_h[0];
	dy = top_h[1];
	index0 = d_index2d(i,j,top_gmax);
	index00 = d_index2d(i-1,j,top_gmax);
	index01 = d_index2d(i+1,j,top_gmax);
	index10 = d_index2d(i,j-1,top_gmax);
	index11 = d_index2d(i,j+1,top_gmax);
	v00 = -cell_center[index00].m_state.m_U[1];
	v01 =  cell_center[index01].m_state.m_U[1];
	v10 =  cell_center[index10].m_state.m_U[0];
	v11 = -cell_center[index11].m_state.m_U[0];

	vorticity = (v00 + v01)/2.0/dy + (v10 + v11)/2.0/dx;
	return vorticity;
}	/* end getVorticity */

void Incompress_Solver_Smooth_2D_Cartesian::copyMeshStates(void)
{
	int i,j,index;
	double **vel = field->vel;
	double *pres = field->pres;
	double *vort = field->vort;
	
	for (i = imin; i <= imax; ++i)
	for (j = jmin; j <= jmax; ++j)
	{
	    index  = d_index2d(i,j,top_gmax);
	    if (ifluid_comp(top_comp[index]))
	    {
		pres[index] = cell_center[index].m_state.m_P;
		//phi[index] = cell_center[index].m_state.m_phi;
	    	vel[0][index] = cell_center[index].m_state.m_U[0];
	    	vel[1][index] = cell_center[index].m_state.m_U[1];
		vort[index] = getVorticity(i,j);
	    }
	    else
	    {
	    	pres[index] = 0.0;
	    	vel[0][index] = 0.0;
	    	vel[1][index] = 0.0;
		vort[index] = 0.0;
	    }
	}
	FT_ParallelExchGridArrayBuffer(pres,front);
	FT_ParallelExchGridArrayBuffer(vort,front);
	FT_ParallelExchGridArrayBuffer(vel[0],front);
	FT_ParallelExchGridArrayBuffer(vel[1],front);
}	/* end copyMeshStates */



void Incompress_Solver_Smooth_2D_Cartesian::
	computeDiffusion(void)
{
	COMPONENT comp;
        int index,index_nb[4],size;
        int I,I_nb[4];
        int i,j,l,nb,icoords[MAXD];
        L_STATE state;
        double coords[MAXD],crx_coords[MAXD];
	double coeff[4],mu[4],mu0,rho,rhs,U_nb[4];
        double speed;
        double *x;
	GRID_DIRECTION dir[4] = {WEST,EAST,SOUTH,NORTH};
	POINTER intfc_state;
	HYPER_SURF *hs;
	PetscInt num_iter;
	double rel_residual;
	double aII;
	int status;

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

            	for (nb = 0; nb < 4; nb++)
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

            	getRectangleCenter(index, coords);
            	computeSourceTerm(coords, state);

        	//first equation  decoupled, some terms may be lost
		aII = 1+coeff[0]+coeff[1]+coeff[2]+coeff[3];
            	rhs = (1-coeff[0]-coeff[1]-coeff[2]-coeff[3])
				*cell_center[index].m_state.m_U[l];

            	for (nb = 0; nb < 4; nb++)
            	{
		    status = (*findStateAtCrossing)(front,icoords,dir[nb],comp,
                                &intfc_state,&hs,crx_coords);
	    	    if (status == NO_PDE_BOUNDARY)
		    {
            	    	solver.Set_A(I,I_nb[nb],-coeff[nb]);
            	    	rhs += coeff[nb]*U_nb[nb];
		    }
		    else
		    {
                        if (status == DIRICHLET_PDE_BOUNDARY &&
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

	    //if(rel_residual > 1)
	    //{
	    //	solver.Reset_x();
	    //	solver.Solve_GMRES();
	    //	solver.GetNumIterations(&num_iter);
            //	solver.GetFinalRelativeResidualNorm(&rel_residual);
	    //}
	    stop_clock("After Petsc solve");


            // get back the solution
            solver.Get_x(x);

            if (debugging("PETSc"))
                (void) printf("Incompress_Solver_Smooth_2D_Cartesian::"
			"computeDiffusion: "
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

        FT_FreeThese(1,x);
}       /* end computeDiffusion */


void Incompress_Solver_Smooth_2D_Cartesian::computePressurePmI(void)
{
        int i,j,index;

	for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
	{
            index = d_index2d(i,j,top_gmax);
            cell_center[index].m_state.m_P += cell_center[index].m_state.m_phi;
	    cell_center[index].m_state.m_q = cell_center[index].m_state.m_P;
	}
}        /* end computePressurePmI2d */


void Incompress_Solver_Smooth_2D_Cartesian::computePressurePmII(void)
{
        int i,j,index;
        double mu0;

	for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
	{
            index = d_index2d(i,j,top_gmax);
            mu0 = 0.5*cell_center[index].m_state.m_mu;
	    cell_center[index].m_state.m_P +=
			cell_center[index].m_state.m_phi -
                        accum_dt*mu0*cell_center[index].m_state.div_U;
	    cell_center[index].m_state.m_q = cell_center[index].m_state.m_P;
	}
}        /* end computePressurePmII2d */

void Incompress_Solver_Smooth_2D_Cartesian::computePressurePmIII(void)
{
        int i,j,index;
        double mu0;

	for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
	{
            index = d_index2d(i,j,top_gmax);
            mu0 = 0.5*cell_center[index].m_state.m_mu;
            cell_center[index].m_state.m_P = 
				cell_center[index].m_state.m_phi -
                        	accum_dt*mu0*cell_center[index].m_state.div_U;
	    cell_center[index].m_state.m_q = 0.0;
	}
}        /* end computePressurePmIII2d */


void Incompress_Solver_Smooth_2D_Cartesian::computePressure(void)
{
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
	    (void) printf("Unknown computePressure scheme!\n");
	    clean_up(ERROR);
	}
	computeGradientQ();
}	/* end computePressure */


void Incompress_Solver_Smooth_2D_Cartesian::computeGradientQ(void)
{
	int i,j,l,index;
	double *grad_q;
	int icoords[MAXD];

	for (j = 0; j < top_gmax[1]; ++j)
	for (i = 0; i < top_gmax[0]; ++i)
	{
	    index = d_index2d(i,j,top_gmax);
	    array[index] = cell_center[index].m_state.m_q;
	}
	for (j = jmin; j <= jmax; j++)
	for (i = imin; i <= imax; i++)
	{
	    index = d_index2d(i,j,top_gmax);
	    grad_q = cell_center[index].m_state.grad_q;
	    icoords[0] = i;
	    icoords[1] = j;
	    computeFieldPointGrad(icoords,array,grad_q);
	}
	for (l = 0; l < dim; ++l)
	{
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
	    	index = d_index2d(i,j,top_gmax);
		array[index] = cell_center[index].m_state.grad_q[l];
	    }
	    scatMeshArray();
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {
	    	index = d_index2d(i,j,top_gmax);
		cell_center[index].m_state.grad_q[l] = array[index];
	    }
	}
}	/* end computeGradientQ2d */


void Incompress_Solver_Smooth_2D_Cartesian::surfaceTension(
	double *coords,
	HYPER_SURF_ELEMENT *hse,
	HYPER_SURF *hs,
	double *force,
	double sigma)
{
	int i,k,nb;
	BOND *pb,*bonds[100];
	double kappa,nor[MAXD];
	double kappa0,nor0[MAXD];
	double kappa1,nor1[MAXD];
	double len,delta;
	HYPER_SURF_ELEMENT *phse;
	double p[MAXD];
	

	BondAndNeighbors(hse,hs,&nb,bonds,3);

	for (i = 0; i < dim; ++i) force[i] = 0.0;
	for (k = 0; k < nb; ++k)
	{
	    pb = bonds[k];
	    for (i = 0; i < dim; ++i) 
		p[i] = 0.5*(Coords(pb->start)[i] + Coords(pb->end)[i]);
	    delta = smoothedDeltaFunction(coords,p);
	    if (delta == 0.0) continue;

	    len = bond_length(pb);
	    phse = Hyper_surf_element(pb);
	    GetFrontNormal(pb->start,phse,hs,nor0,front);
	    GetFrontNormal(pb->end,phse,hs,nor1,front);
	    for (i = 0; i < dim; ++i) nor[i] = 0.5*(nor0[i] + nor1[i]);
	    GetFrontCurvature(pb->start,phse,hs,&kappa0,front);
	    GetFrontCurvature(pb->end,phse,hs,&kappa1,front);
	    kappa = 0.5*(kappa0 + kappa1);
	    for (i = 0; i < dim; ++i) 
	    {
		force[i] += delta*sigma*len*kappa*nor[i];
	    }
	}
}	/* end surfaceTension2d */

void Incompress_Solver_Smooth_2D_Cartesian::setInitialCondition()
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
	printf("Test setInitialCondition() step 5\n");
}       /* end setInitialCondition */

double Incompress_Solver_Smooth_2D_Cartesian::computeFieldPointDiv(
        int *icoords,
        double **field)
{
        int index;
        COMPONENT comp;
        int i, j;
	int index_nb[4];
        double div,u_nb[2],v_nb[2];
        double crx_coords[MAXD];
        POINTER intfc_state;
        HYPER_SURF *hs;
	double h[MAXD];

	i = icoords[0];
	j = icoords[1];
	index = d_index2d(i,j,top_gmax);
        comp = top_comp[index];

	index_nb[0] = d_index2d(i-1,j,top_gmax);
	index_nb[1] = d_index2d(i+1,j,top_gmax);
	index_nb[2] = d_index2d(i,j-1,top_gmax);
	index_nb[3] = d_index2d(i,j+1,top_gmax);


	h[0] = 2.0*top_h[0];
	if ((*findStateAtCrossing)(front,icoords,WEST,
		comp,&intfc_state,&hs,crx_coords) != NO_PDE_BOUNDARY)
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
		comp,&intfc_state,&hs,crx_coords) != NO_PDE_BOUNDARY)
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
		comp,&intfc_state,&hs,crx_coords) != NO_PDE_BOUNDARY)
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
		comp,&intfc_state,&hs,crx_coords) != NO_PDE_BOUNDARY)
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


        div = (u_nb[1] - u_nb[0])/h[0] + (v_nb[1] - v_nb[0])/h[1];
        return div;
}       /* end computeFieldPointDiv */

void Incompress_Solver_Smooth_2D_Cartesian::computeFieldPointGrad(
        int *icoords,
        double *field,
        double *grad_field)
{
        int index;
        COMPONENT comp;
        int i,j,nb;
	int index_nb[4];
        double p_nbedge[4],p0;  //the p values on the cell edges and cell center
        double crx_coords[MAXD];
        POINTER intfc_state;
        HYPER_SURF *hs;
        GRID_DIRECTION dir[4] = {WEST,EAST,SOUTH,NORTH};      
	double h[MAXD];

	i = icoords[0];
	j = icoords[1];
	index = d_index2d(i,j,top_gmax);
        comp = top_comp[index];
	p0 = field[index];

	index_nb[0] = d_index2d(i-1,j,top_gmax);
	index_nb[1] = d_index2d(i+1,j,top_gmax);
	index_nb[2] = d_index2d(i,j-1,top_gmax);
	index_nb[3] = d_index2d(i,j+1,top_gmax);
	for (i = 0; i < 2; ++i)
            h[i] = 2*top_h[i];

	for (nb = 0; nb < 4; nb++)
	{
	    /*
	    if(FT_StateStructAtGridCrossing(front,icoords,dir[nb],
			comp,&intfc_state,&hs,crx_coords) &&
		        wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    {
		p_nbedge[nb] = p0;
		if (wave_type(hs) == DIRICHLET_BOUNDARY)
                    h[nb/2] = top_h[nb/2];
	    }
	    */
	    if(FT_StateStructAtGridCrossing(front,icoords,dir[nb],
			comp,&intfc_state,&hs,crx_coords)) 
	    {
		p_nbedge[nb] = p0;
		if (wave_type(hs) == DIRICHLET_BOUNDARY ||
		    wave_type(hs) == FIRST_PHYSICS_WAVE_TYPE)
                    h[nb/2] = top_h[nb/2];
	    }
	    else
		p_nbedge[nb] = field[index_nb[nb]];
	}
	grad_field[0] = (p_nbedge[1] - p_nbedge[0])/h[0];
	grad_field[1] = (p_nbedge[3] - p_nbedge[2])/h[1];
}

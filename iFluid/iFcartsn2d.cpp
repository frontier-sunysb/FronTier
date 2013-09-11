/*******************************************************************
 * 	               iFcartsn2d.cpp	
 *******************************************************************/
#include "iFluid.h"
#include "solver.h"

static double (*getStateVel[3])(POINTER) = {getStateXvel,getStateYvel,
                                        getStateZvel};


void Incompress_Solver_Smooth_2D_Cartesian::computeAdvection(void)
{
	int i,j,index;
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
	hyperb_solver.getStateVel[0] = getStateXvel;
	hyperb_solver.getStateVel[1] = getStateYvel;
	hyperb_solver.rho1 = iFparams->rho1;
	hyperb_solver.rho2 = iFparams->rho2;
	hyperb_solver.findStateAtCrossing = ifluid_find_state_at_crossing;
	hyperb_solver.solveRungeKutta();
}

void Incompress_Solver_Smooth_2D_Cartesian::computeProjection(void)
{
	setIndexMap();
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
}	/* end computeProjection */

void Incompress_Solver_Smooth_2D_Cartesian::computeProjectionCim(void)
{
	static CIM_ELLIPTIC_SOLVER elliptic_solver(*front);
	int index;
	int i,j,l,icoords[MAXD];
	double **vel = field->vel;
	double *phi = field->phi;
	double *div_U = field->div_U;
	double sum_div;
	double value;

	sum_div = 0.0;
	max_value = 0.0;

	/* Compute velocity divergence */
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    icoords[0] = i;
	    icoords[1] = j;
	    index  = d_index(icoords,top_gmax,dim);
	    source[index] = computeFieldPointDiv(icoords,vel);
	    diff_coeff[index] = 1.0/field->rho[index];
	}
	FT_ParallelExchGridArrayBuffer(source,front);
	FT_ParallelExchGridArrayBuffer(diff_coeff,front);
	for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
	{
	    index  = d_index2d(i,j,top_gmax);
	    div_U[index] = source[index];
	    source[index] /= -accum_dt;
	    array[index] = phi[index];
	}

	if(debugging("step_size"))
	{
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		index = d_index2d(i,j,top_gmax);
	        value = fabs(div_U[index]);
		sum_div = sum_div + div_U[index];
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
	    phi[index] = array[index];
	}
}	/* end computeProjectionCim */

void Incompress_Solver_Smooth_2D_Cartesian::computeProjectionDouble(void)
{
	iFparams->total_div_cancellation = YES;
	computeProjectionSimple(); 
	return;
}	/* end computeProjectionDouble */

void Incompress_Solver_Smooth_2D_Cartesian::computeProjectionSimple(void)
{
	static ELLIPTIC_SOLVER elliptic_solver(*front);
	int index;
	int i,j,l,icoords[MAXD];
	double **vel = field->vel;
	double *phi = field->phi;
	double *div_U = field->div_U;
	double sum_div;
	double value;
	double vmin[MAXD],vmax[MAXD];

	sum_div = 0.0;
	max_value = 0.0;
	for (l = 0; l < dim; ++l)
	{
	    vmin[l] = HUGE;
	    vmax[l] = -HUGE;
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
	for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
	{
	    index  = d_index2d(i,j,top_gmax);
	    if (!ifluid_comp(top_comp[index]))
		continue;
	    div_U[index] = source[index];
	    source[index] /= accum_dt;
	    array[index] = phi[index];
	}

	if(debugging("step_size"))
	{
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		index = d_index2d(i,j,top_gmax);
	        if (!ifluid_comp(top_comp[index]))
		    continue;
	        value = fabs(div_U[index]);
		sum_div = sum_div + div_U[index];
		if (max_value < value)
		    max_value = value;
	    }
	    pp_global_sum(&sum_div,1);
	    (void) printf("\nThe summation of divergence of U is %.16g\n",
					sum_div);
	    pp_global_max(&max_value,1);
	    (void) printf("\nThe max value of divergence of U is %.16g\n",
					max_value);
	    max_value = 0.0;
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
            elliptic_solver.ij_to_I = ij_to_I;
            elliptic_solver.ilower = ilower;
            elliptic_solver.iupper = iupper;
	    if (iFparams->total_div_cancellation)
	    	elliptic_solver.dsolve(array);
	    else
	    	elliptic_solver.solve(array);
	    paintSolvedGridPoint();
	}
	//viewTopVariable(front,array,NO,0.0,0.0,(char*)"test-simple",
	//			(char*)"array");

	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index  = d_index2d(i,j,top_gmax);
	    phi[index] = array[index];
	}
}	/* end computeProjectionSimple */

void Incompress_Solver_Smooth_2D_Cartesian::computeNewVelocity(void)
{
	int i,j,k,l,index;
	double grad_phi[2],rho;
	COMPONENT comp;
	double speed;
	int icoords[MAXD];
	double **vel = field->vel;
	double *phi = field->phi;
	double vmin[MAXD],vmax[MAXD];

	max_speed = 0.0;
	for (i = 0; i < dim; ++i)
        {
            vmin[i] = HUGE;
            vmax[i] = -HUGE;
        }

	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index  = d_index2d(i,j,top_gmax);
	    array[index] = phi[index];
	}
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    index = d_index2d(i,j,top_gmax);
	    comp = top_comp[index];
	    if (!ifluid_comp(comp))
	    {
		vel[0][index] = 0.0;
		vel[1][index] = 0.0;
		continue;
	    }
	    rho = field->rho[index];
	    icoords[0] = i;
	    icoords[1] = j;
	    computeFieldPointGrad(icoords,array,grad_phi);
	    vel[0][index] -= accum_dt/rho*grad_phi[0];
	    vel[1][index] -= accum_dt/rho*grad_phi[1];
	    speed = fabs(vel[0][index]) + fabs(vel[1][index]);
	    if (speed > max_speed)
		max_speed = speed;
	    for (l = 0; l < dim; ++l)
	    {
		if (vmin[l] > vel[l][index]) vmin[l] = vel[l][index];
		if (vmax[l] < vel[l][index]) vmax[l] = vel[l][index];
	    }
	}
	for (k = 0; k < 2; ++k)
	{
	    for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
	    {	
	    	index  = d_index2d(i,j,top_gmax);
	    	array[index] = vel[k][index];
	    }
	    scatMeshArray();
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {	
	    	index  = d_index2d(i,j,top_gmax);
	    	vel[k][index] = array[index];
	    }
	}
	if (debugging("check_div"))
        {
	    checkVelocityDiv("After computeNewVelocity2d()",vmin,vmax);
        }
	pp_global_max(&max_speed,1);
}	/* end computeNewVelocity2d */


void Incompress_Solver_Smooth_2D_Cartesian::computeSourceTerm(
	double *coords, 
	double *source) 
{
	int i;
	for (i = 0; i < dim; ++i)
	    source[i] = iFparams->gravity[i];
}	/* end computeSourceTerm */

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
	double **vel = field->vel;

	dx = top_h[0];
	dy = top_h[1];
	index0 = d_index2d(i,j,top_gmax);
	index00 = d_index2d(i-1,j,top_gmax);
	index01 = d_index2d(i+1,j,top_gmax);
	index10 = d_index2d(i,j-1,top_gmax);
	index11 = d_index2d(i,j+1,top_gmax);
	v00 = -vel[1][index00];
	v01 =  vel[1][index01];
	v10 =  vel[0][index10];
	v11 = -vel[0][index11];

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
	return computeDiffusionCN();
}

void Incompress_Solver_Smooth_2D_Cartesian::
	computeDiffusionCN(void)
{
        COMPONENT comp;
        int index,index_nb[4],size;
        int I,I_nb[4];
        int i,j,k,l,nb,icoords[MAXD];
        double coords[MAXD], crx_coords[MAXD];
        double coeff[4],mu[4],mu0,rho,rhs,U_nb[4];
        double speed;
        double *x;
        GRID_DIRECTION dir[4] = {WEST,EAST,SOUTH,NORTH};
        POINTER intfc_state;
        HYPER_SURF *hs;
        PetscInt num_iter;
        double rel_residual;
        double aII;
        double source[MAXD];
        double **vel = field->vel;
        double **f_surf = field->f_surf;
        INTERFACE *grid_intfc = front->grid_intfc;
	int status;
        double vel_min = HUGE;

        if (debugging("trace"))
            (void) printf("Entering Incompress_Solver_Smooth_3D_Cartesian::"
                        "computeDiffusionCN()\n");
	start_clock("computeDiffusionCN");

        setIndexMap();
        max_speed = 0.0;

        size = iupper - ilower;
        FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));

	for (l = 0; l < dim; ++l)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index  = d_index2d(i,j,top_gmax);
            if (vel_min > vel[l][index]) vel_min = vel[l][index];
        }
	if (vel_min < 0.0)
	    vel_min *= 5.0;

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


		mu0 = field->mu[index];
		rho = field->rho[index];

            	for (nb = 0; nb < 4; nb++)
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

            	getRectangleCenter(index, coords);
            	computeSourceTerm(coords, source);

        	//first equation  decoupled, some terms may be lost
		aII = 1+coeff[0]+coeff[1]+coeff[2]+coeff[3];
            	rhs = (1-coeff[0]-coeff[1]-coeff[2]-coeff[3])
				*(vel[l][index] - vel_min);

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
                        if (status == CONST_P_PDE_BOUNDARY)
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
                    array[index] = x[I-ilower] + vel_min;
                else
                    array[index] = 0.0;
            }
            scatMeshArray();
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index2d(i,j,top_gmax);
                vel[l][index] = array[index];
            }
        }
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    index = d_index2d(i,j,top_gmax);
            speed = fabs(vel[0][index]) + fabs(vel[1][index]);
            if (speed > max_speed)
                max_speed = speed;
	}
        pp_global_max(&max_speed,1);
	stop_clock("computeDiffusionCN");

        FT_FreeThese(1,x);
        if (debugging("trace"))
            (void) printf("Leaving Incompress_Solver_Smooth_3D_Cartesian::"
                        "computeDiffusionCN()\n");
}	/* end computeDiffusionCN */

void Incompress_Solver_Smooth_2D_Cartesian::computePressurePmI(void)
{
        int i,j,index;
	double *pres = field->pres;
	double *phi = field->phi;
	double *q = field->q;

	for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
	{
            index = d_index2d(i,j,top_gmax);
            pres[index] += phi[index];
	    q[index] = pres[index];
	}
}        /* end computePressurePmI2d */


void Incompress_Solver_Smooth_2D_Cartesian::computePressurePmII(void)
{
        int i,j,index;
        double mu0;
	double *pres = field->pres;
	double *phi = field->phi;
	double *q = field->q;
	double *div_U = field->div_U;

	for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
	{
            index = d_index2d(i,j,top_gmax);
            mu0 = 0.5*field->mu[index];
	    pres[index] += phi[index] - accum_dt*mu0*div_U[index];
	    q[index] = pres[index];
	}
}        /* end computePressurePmII2d */

void Incompress_Solver_Smooth_2D_Cartesian::computePressurePmIII(void)
{
        int i,j,index;
        double mu0;
	double *pres = field->pres;
	double *phi = field->phi;
	double *q = field->q;
	double *div_U = field->div_U;

	for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
	{
            index = d_index2d(i,j,top_gmax);
            mu0 = 0.5*field->mu[index];
            pres[index] = phi[index] -
                        	accum_dt*mu0*div_U[index];
	    q[index] = 0.0;
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
	double **grad_q = field->grad_q;
	int icoords[MAXD];
	double *q = field->q;
	double point_grad_q[MAXD];

	for (j = 0; j < top_gmax[1]; ++j)
	for (i = 0; i < top_gmax[0]; ++i)
	{
	    index = d_index2d(i,j,top_gmax);
	    array[index] = q[index];
	}
	for (j = jmin; j <= jmax; j++)
	for (i = imin; i <= imax; i++)
	{
	    index = d_index2d(i,j,top_gmax);
	    icoords[0] = i;
	    icoords[1] = j;
	    computeFieldPointGrad(icoords,array,point_grad_q);
	}
	for (l = 0; l < dim; ++l)
	{
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
	    	index = d_index2d(i,j,top_gmax);
		array[index] = grad_q[l][index];
	    }
	    scatMeshArray();
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {
	    	index = d_index2d(i,j,top_gmax);
		grad_q[l][index] = array[index];
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
	    //cell_center[i].m_state.setZero();
	    comp = top_comp[i];
	    if (getInitialState != NULL)
	    	(*getInitialState)(comp,coords,field,i,dim,iFparams);
        }
	computeGradientQ();
        copyMeshStates();
	setAdvectionDt();
}       /* end setInitialCondition */

void Incompress_Solver_Smooth_2D_Cartesian::computeProjectionDual(void)
{
	static DUAL_ELLIPTIC_SOLVER dual_elliptic_solver(*front);
	int index;
	int i,j,l,icoords[MAXD];
	double **vel = field->vel;
	double *phi = field->phi;
	double *div_U = field->div_U;
	double sum_div;
	double value;
	printf("Entering computeProjectionDual()\n");
	FT_MakeCompGridIntfc(front);
	setDualDomain();
	setDualGlobalIndex();
	setDualIndexMap();
	printf(" top_gmax = %d %d\n",top_gmax[0],top_gmax[1]); 
	printf("ctop_gmax = %d %d\n",ctop_gmax[0],ctop_gmax[1]); 
	printf(" top_L = %f %f\n",top_L[0],top_L[1]); 
	printf("ctop_L = %f %f\n",ctop_L[0],ctop_L[1]); 
	printf(" top_U = %f %f\n",top_U[0],top_U[1]); 
	printf("ctop_U = %f %f\n",ctop_U[0],ctop_U[1]); 
	printf(" imin = %d %d   imax = %d %d\n",imin,jmin,imax,jmax);
	printf("cimin = %d %d  cimax = %d %d\n",cimin,cjmin,cimax,cjmax);
	printf(" ilower = %d   iupper = %d\n",ilower,iupper);
	printf("cilower = %d  ciupper = %d\n",cilower,ciupper);

	sum_div = 0.0;
	max_value = 0.0;

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
	    diff_coeff[index] = 1.0/field->rho[index];
	}
	FT_ParallelExchGridArrayBuffer(source,front);
	FT_ParallelExchGridArrayBuffer(diff_coeff,front);
	for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
	{
	    index  = d_index2d(i,j,ctop_gmax);
	    if (!ifluid_comp(ctop_comp[index]))
		continue;
	    div_U[index] = source[index];
	    source[index] /= accum_dt;
	    array[index] = phi[index];
	}

	if(debugging("step_size"))
	{
	    for (j = cjmin; j <= cjmax; j++)
	    for (i = cimin; i <= cimax; i++)
	    {
		index = d_index2d(i,j,ctop_gmax);
	        if (!ifluid_comp(ctop_comp[index]))
		    continue;
	        value = fabs(div_U[index]);
		sum_div = sum_div + div_U[index];
		if (max_value < value)
		    max_value = value;
	    }
	    pp_global_sum(&sum_div,1);
	    (void) printf("\nThe summation of divergence of U is %.16g\n",
					sum_div);
	    pp_global_max(&max_value,1);
	    (void) printf("\nThe max value of divergence of U is %.16g\n",
					max_value);
	    max_value = 0.0;
	}
        dual_elliptic_solver.D = diff_coeff;
        dual_elliptic_solver.source = source;
        dual_elliptic_solver.ilower = cilower;
        dual_elliptic_solver.iupper = ciupper;
        dual_elliptic_solver.soln = array;
        dual_elliptic_solver.obst_comp = SOLID_COMP;
        dual_elliptic_solver.ij_to_I = cij_to_I;
	dual_elliptic_solver.set_solver_domain();
	dual_elliptic_solver.getStateVar = getStatePhi;
	dual_elliptic_solver.findStateAtCrossing = 
		ifluid_find_state_at_dual_crossing;
	dual_elliptic_solver.solve(array);

	for (j = 0; j <= ctop_gmax[1]; j++)
	for (i = 0; i <= ctop_gmax[0]; i++)
	{
	    index  = d_index2d(i,j,ctop_gmax);
	    phi[index] = array[index];
	}

	FT_FreeCompGridIntfc(front);
	clean_up(0);
}	/* end computeProjectionDual */

void Incompress_Solver_Smooth_2D_Cartesian::
        computeDiffusionImplicit(void)
{
        COMPONENT comp;
        int index,index_nb[4],size;
        int I,I_nb[4];
        int i,j,k,l,nb,icoords[MAXD];
        double coords[MAXD], crx_coords[MAXD];
        double coeff[4],mu[4],mu0,rho,rhs,U_nb[4];
        double speed;
        double *x;
        GRID_DIRECTION dir[4] = {WEST,EAST,SOUTH,NORTH};
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
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index  = d_index2d(i,j,top_gmax);
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

                I_nb[0] = ij_to_I[i-1][j]; //west
                I_nb[1] = ij_to_I[i+1][j]; //east
                I_nb[2] = ij_to_I[i][j-1]; //south
                I_nb[3] = ij_to_I[i][j+1]; //north


                mu0   = field->mu[index];
                rho   = field->rho[index];

                for (nb = 0; nb < 4; nb++)
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

                getRectangleCenter(index, coords);
                computeSourceTerm(coords, source);

                aII = 1+coeff[0]+coeff[1]+coeff[2]+coeff[3];
                rhs = vel[l][index] - vel_min;

                for(nb = 0; nb < 4; nb++)
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
            //solver.Solve_GMRES();
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

            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                I = ij_to_I[i][j];
                index = d_index2d(i,j,top_gmax);
                if (I >= 0)
                    array[index] = x[I-ilower] + vel_min;
                else
                    array[index] = 0.0;
            }
            scatMeshArray();
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index2d(i,j,top_gmax);
                vel[l][index] = array[index];
            }
        }
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index2d(i,j,top_gmax);
            speed = fabs(vel[0][index]) + fabs(vel[1][index]);
            if (speed > max_speed)
            {
                icrds_max[0] = i;
                icrds_max[1] = j;
                max_speed = speed;
            }
        }
        pp_global_max(&max_speed,1);
        FT_FreeThese(1,x);

        if (debugging("trace"))
            (void) printf("Leaving Incompress_Solver_Smooth_3D_Cartesian::"
                        "computeDiffusionImplicit()\n");
}       /* end computeDiffusionImplicit */

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
}       /* end parab_find_state_at_crossing */

void Incompress_Solver_Smooth_2D_Cartesian::
        computeDiffusionParab(void)
{
        static PARABOLIC_SOLVER parab_solver(*front);

        COMPONENT comp;
        int index;
        int i,j,l,icoords[MAXD];
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
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index  = d_index2d(i,j,top_gmax);
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

void Incompress_Solver_Smooth_2D_Cartesian::
	computeDiffusionExplicit(void)
{
        COMPONENT comp;
        int index,index_nb[4],size;
        int I,I_nb[4];
	int i,j,l,nb,icoords[MAXD];
	double coords[MAXD], crx_coords[MAXD];
	double coeff[4],mu[4],mu0,rho,rhs,U_nb[4];
        double speed;
        double *x;
	GRID_DIRECTION dir[4] = {WEST,EAST,SOUTH,NORTH};
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
	    (void) printf("Entering Incompress_Solver_Smooth_2D_Cartesian::"
			"computeDiffusionExplicit()\n");

        setIndexMap();
	max_speed = 0.0;

        size = iupper - ilower;
        FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));

	for (l = 0; l < dim; ++l)
	{
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

            	mu0   = field->mu[index];
            	rho   = field->rho[index];

            	for (nb = 0; nb < 4; nb++)
            	{
                    if (FT_StateStructAtGridCrossing(front,grid_intfc,icoords,
				dir[nb],comp,&intfc_state,&hs,crx_coords) &&
                                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
		    {
			if (wave_type(hs) == DIRICHLET_BOUNDARY &&
                    	    boundary_state_function(hs) &&
                    	    strcmp(boundary_state_function_name(hs),
                    	    "flowThroughBoundaryState") == 0)
                    	    U_nb[nb] = vel[l][index];
			else
			    U_nb[nb] = getStateVel[l](intfc_state);
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

            	getRectangleCenter(index, coords);
            	computeSourceTerm(coords, source);

		rhs = (-coeff[0]-coeff[1]-coeff[2]-coeff[3])*vel[l][index];

		int num_nb = 0;
		for(nb = 0; nb < 4; nb++)
		{
		    rhs += coeff[nb]*U_nb[nb];
		}
		rhs += m_dt*source[l];
		rhs += m_dt*f_surf[l][index];
                x[I-ilower] = vel[l][index] + rhs;
            }

            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                I = ij_to_I[i][j];
                index = d_index2d(i,j,top_gmax);
                if (I >= 0)
                    array[index] = x[I-ilower];
                else
                    array[index] = 0.0;
            }
            scatMeshArray();
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index2d(i,j,top_gmax);
                vel[l][index] = array[index];
            }
        }
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    index = d_index2d(i,j,top_gmax);
            speed = fabs(vel[0][index]) + fabs(vel[1][index]);
            if (speed > max_speed)
	    {
		icrds_max[0] = i;
		icrds_max[1] = j;
                max_speed = speed;
	    }
	}
        pp_global_max(&max_speed,1);
        FT_FreeThese(1,x);

	if (debugging("trace"))
	    (void) printf("Leaving Incompress_Solver_Smooth_2D_Cartesian::"
			"computeDiffusionExplicit()\n");
}       /* end computeDiffusionExplicit */

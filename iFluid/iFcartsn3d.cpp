/***************************************************************
FronTier is a set of libraries that implements differnt types of 
Front Traking algorithms. Front Tracking is a numerical method for 
the solution of partial differential equations whose solutions have 
discontinuities.  

Copyright (C) 1999 by The University at Stony Brook. 

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
****************************************************************/

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
	int i,j,k,l,index;
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
	    for (l = 0; l < dim; ++l)
	    {
		if (vmin[l] > field->vel[l][index])
		    vmin[l] = field->vel[l][index];
		if (vmax[l] < field->vel[l][index])
		    vmax[l] = field->vel[l][index];
	    }
	}
	pp_global_max(&max_speed,1);
	pp_global_min(vmin,dim);
	pp_global_max(vmax,dim);
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

	if (iFparams->num_scheme.ellip_method == DUAL_ELLIP)
        {
            computeNewVelocityDual();
            return;
        }
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
	FT_ParallelExchGridVectorArrayBuffer(vel,front);
	if (debugging("step_size"))
	{
	    (void) printf("Max gradient phi = %f  occuring at: %d %d %d\n",
			max_grad_phi,icrds_max[0],icrds_max[1],icrds_max[2]);
	    (void) printf("Ave gradient phi = %f\n",ave_grad_phi/(imax-imin+1)
			/(jmax-jmin+1)/(kmax-kmin+1));
	}
	if (debugging("check_div"))
	{
	    checkVelocityDiv("After computeNewVelocity3d()");
	}
	pp_global_max(&max_speed,1);
	pp_global_min(vmin,dim);
	pp_global_max(vmax,dim);
}	/* end computeNewVelocity3d */

boolean Incompress_Solver_Smooth_3D_Cartesian::InsideSolid(int* icoords)
{
	int icu[MAXD],icl[MAXD];
	for (int m = 0; m < dim; ++m)
	    icl[m] = icoords[m] - offset[m];
	for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
        for (int k = 0; k < 2; ++k)
        {
            icu[0] = icl[0] + i;
            icu[1] = icl[1] + j;
            icu[2] = icl[2] + k;
            int index = d_index(icu,ctop_gmax,dim);
	    if (ifluid_comp(ctop_comp[index]))
		return NO;
	}
	return YES;
}

void Incompress_Solver_Smooth_3D_Cartesian::computeNewVelocityDual(void)
{
	int i,j,k,l,m,n,mm,nn,ii;
	int index,index_tmp,index_u,index_l,dual_compu,dual_compl;
	double grad_phi[MAXD],rho;
        COMPONENT comp;
        double v_ave,speed,denom;
        int icoords[MAXD];
	int ic[MAXD],icu[MAXD],icl[MAXD],ictmp[MAXD];
        double **vel = field->vel;
        double *d_phi = field->d_phi;
	static double **V;
	int symmetry[MAXD];

	computeProjectionDual();

	if (V == NULL)
            FT_MatrixMemoryAlloc((POINTER*)&V,3,
			(top_gmax[0]+1)*(top_gmax[1]+1)*(top_gmax[2]+1), 
                	sizeof(double));

	max_speed = 0.0;
        for (i = 0; i < dim; ++i)
        {
            vmin[i] = HUGE;
            vmax[i] = -HUGE;
        }
        for (k = 0; k <= top_gmax[2]; ++k)
        for (j = 0; j <= top_gmax[1]; ++j)
        for (i = 0; i <= top_gmax[0]; ++i)
	{
            index = d_index3d(i,j,k,top_gmax);
	    for (l = 0; l < dim; ++l)
		V[l][index] = vel[l][index];
	}

        for (k = kmin; k <= kmax; ++k)
        for (j = jmin; j <= jmax; ++j)
        for (i = imin; i <= imax; ++i)
        {
	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;
            index = d_index(icoords,top_gmax,dim);
	    if (InsideSolid(icoords))
	    {
		for (l = 0; l < dim; ++l)
                    vel[l][index] = 0.0;
		continue;
	    }
	    rho = field->rho[index];
	    computeDualFieldPointGrad(icoords,d_phi,grad_phi);
	    speed = 0.0;
	    for (l = 0; l < dim; ++l)
	    {
		v_ave = 0.0; 
		denom = 0.0;
		icl[l] = icoords[l] - offset[l];
		icu[l] = icoords[l] - offset[l] + 1;
		for (m = 0; m < 2; ++m)
                for (n = 0; n < 2; ++n)
                {
                    icl[(l+1)%dim] = icu[(l+1)%dim] = icoords[(l+1)%dim] -
                        offset[(l+1)%dim] + m;
                    icl[(l+2)%dim] = icu[(l+2)%dim] = icoords[(l+2)%dim] -
                        offset[(l+2)%dim] + n;
                    index_u = d_index(icu,ctop_gmax,dim);
                    index_l = d_index(icl,ctop_gmax,dim);
                    dual_compu = ctop_comp[index_u];
                    dual_compl = ctop_comp[index_l];
		    if (ifluid_comp(dual_compu)||ifluid_comp(dual_compl))
		    {
			for (ii = 0; ii < dim; ++ii)
			    ictmp[ii] = icu[ii] - 1 + offset[ii];
			for (mm = 0; mm < 2; ++mm)
			for (nn = 0; nn < 2; ++nn)
			{
			    ic[l] = ictmp[l];
			    ic[(l+1)%dim] = ictmp[(l+1)%dim] + mm;
                            ic[(l+2)%dim] = ictmp[(l+2)%dim] + nn;
			    index_tmp = d_index(ic,top_gmax,dim);
			    v_ave += V[l][index_tmp];
                            denom += 1.0;
			}
		    }
		}
		v_ave /= denom;
		vel[l][index] = v_ave - accum_dt/rho*grad_phi[l];
                if (vmin[l] > vel[l][index]) vmin[l] = vel[l][index];
                if (vmax[l] < vel[l][index]) vmax[l] = vel[l][index];
		speed += sqr(vel[l][index]);
	    }
            speed = sqrt(speed);
            if (speed > max_speed) max_speed = speed;
	}

	FT_ParallelExchGridVectorArrayBuffer(vel,front);

        if (debugging("check_div"))
	    checkVelocityDiv("Before extractFlowThroughVelocity()");
	extractFlowThroughVelocity();
        if (debugging("check_div"))
            checkVelocityDiv("After extractFlowThroughVelocity()");

        pp_global_max(&max_speed,1);
        pp_global_min(vmin,dim);
        pp_global_max(vmax,dim);
}	/* end computeNewVelocityDual */

void Incompress_Solver_Smooth_3D_Cartesian::
	computeSourceTerm(double *coords, double *source) 
{
        int i;
	
        if(iFparams->if_buoyancy)
        {
	    int ic[MAXD],index;
            rect_in_which(coords,ic,top_grid);
            index = d_index(ic,top_gmax,dim);
            for (i = 0; i < dim; ++i)
                source[i] = field->ext_accel[i][index];
        }
        else
        {
            for (i = 0; i < dim; ++i)
                source[i] = iFparams->gravity[i];
        }
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
	if (iFparams->num_scheme.ellip_method == DUAL_ELLIP)
	    updateComponent();
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

        if(iFparams->if_ref_pres == YES)
            setReferencePressure();

	setAdvectionDt();
	stop_clock("solve");
}	/* end solve */


void Incompress_Solver_Smooth_3D_Cartesian::copyMeshStates()
{
	int i,j,k,d,index;
	double *pres = field->pres;

	for (i = imin; i <= imax; ++i)
	for (j = jmin; j <= jmax; ++j)
	for (k = kmin; k <= kmax; ++k)
	{
	    index  = d_index3d(i,j,k,top_gmax);
	    if (!ifluid_comp(top_comp[index]))
	    	pres[index] = 0.0;
	}
	FT_ParallelExchGridArrayBuffer(pres,front,NULL);
}	/* end copyMeshStates */

void Incompress_Solver_Smooth_3D_Cartesian::
	computeDiffusion(void)
{
	return computeDiffusionCN();
	//return computeDiffusionImplicit();
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

	if (debugging("trace"))
	    (void) printf("Entering Incompress_Solver_Smooth_3D_Cartesian::"
			"computeDiffusionCN()\n");

        setIndexMap();
	max_speed = 0.0;
	for (l = 0; l < dim; ++l)
        {
            vmin[l] = HUGE;
            vmax[l] = -HUGE;
        }

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
		      		(vel[l][index]);

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
                    vel[l][index] = x[I-ilower];
                else
                    vel[l][index] = 0.0;
            }
        }
	FT_ParallelExchGridVectorArrayBuffer(vel,front);

	max_speed = 0.0;
	for (i = 0; i < 3; ++i)
	{
	    vmin[i] = HUGE;	vmax[i] = -HUGE;
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
	    for (l = 0; l < dim; ++l)
	    {
		if (vmin[l] > vel[l][index]) vmin[l] = vel[l][index];
		if (vmax[l] < vel[l][index]) vmax[l] = vel[l][index];
	    }
	}
        pp_global_max(&max_speed,1);
        pp_global_min(vmin,dim);
        pp_global_max(vmax,dim);
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
            for (k = kmin; k <= kmax; k++)
	    {
            	index = d_index3d(13,13,k,top_gmax);
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
	}
	FT_ParallelExchGridVectorArrayBuffer(vel,front);
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

	if (debugging("trace"))
	    (void) printf("Entering Incompress_Solver_Smooth_3D_Cartesian::"
			"computeDiffusionImplicit()\n");

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

		aII = 1+coeff[0]+coeff[1]+coeff[2]+coeff[3]+coeff[4]+coeff[5];
		rhs = vel[l][index];

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
	}
	FT_ParallelExchGridVectorArrayBuffer(vel,front);
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
	    for (l = 0; l < dim; ++l)
		grad_q[l][index] = point_grad_q[l];
	}
	FT_ParallelExchGridVectorArrayBuffer(grad_q,front);
}	/* end computeGradientQ */

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
	if (iFparams->num_scheme.ellip_method == DUAL_ELLIP)
            FT_MakeCompGridIntfc(front);
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
        case DUAL_ELLIP:
            computeProjectionSimple();
            return;
        case DOUBLE_ELLIP:
            computeProjectionDouble();
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
	static DUAL_ELLIPTIC_SOLVER dual_elliptic_solver(*front);
	int index;
	int i,j,k,l,icoords[MAXD];
	double **vel = field->vel;
	double *phi = field->phi;
	double *div_U = field->div_U;
	double *d_phi = field->d_phi;
	double sum_div;
	double value;
	setDualGlobalIndex();
	setDualIndexMap();

	updateComponent();

	sum_div = 0.0;
	max_value = 0.0;

	/* Compute velocity divergence */
	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    index  = d_index3d(i,j,k,top_gmax);
	    diff_coeff[index] = 1.0/field->rho[index];
	}
	FT_ParallelExchGridArrayBuffer(diff_coeff,front,NULL);
	for (k = ckmin; k <= ckmax; k++)
	for (j = cjmin; j <= cjmax; j++)
        for (i = cimin; i <= cimax; i++)
	{
	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;
	    index  = d_index(icoords,ctop_gmax,dim);
	    if (!ifluid_comp(ctop_comp[index]))
		continue;
	    source[index] = computeDualFieldPointDiv(icoords,vel);
	}

	FT_ParallelExchCompGridArrayBuffer(source,front,NULL);

	for (k = 0; k <= ctop_gmax[2]; k++)
	for (j = 0; j <= ctop_gmax[1]; j++)
        for (i = 0; i <= ctop_gmax[0]; i++)
	{
	    index  = d_index3d(i,j,k,ctop_gmax);
	    if (!ifluid_comp(ctop_comp[index]))
		continue;
	    div_U[index] = source[index];
	    source[index] /= accum_dt;
	    array[index] = d_phi[index];
	}

	if(debugging("step_size"))
	{
	    for (k = ckmin; k <= ckmax; k++)
	    for (j = cjmin; j <= cjmax; j++)
	    for (i = cimin; i <= cimax; i++)
	    {
		index = d_index3d(i,j,k,ctop_gmax);
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
        dual_elliptic_solver.ijk_to_I = cijk_to_I;
	dual_elliptic_solver.set_solver_domain();
	dual_elliptic_solver.getStateVar = getStatePhi;
	dual_elliptic_solver.findStateAtCrossing = 
		ifluid_find_state_at_dual_crossing;
	dual_elliptic_solver.solve(array);

	for (k = 0; k <= ctop_gmax[2]; k++)
	for (j = 0; j <= ctop_gmax[1]; j++)
	for (i = 0; i <= ctop_gmax[0]; i++)
	{
	    index  = d_index3d(i,j,k,ctop_gmax);
	    d_phi[index] = array[index];
	}
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
	FT_ParallelExchGridArrayBuffer(source,front,NULL);
        FT_ParallelExchGridArrayBuffer(diff_coeff,front,NULL);
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
            checkVelocityDiv("Before computeProjection()");
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
	static boolean first = YES;

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
	FT_ParallelExchGridArrayBuffer(nu,front,NULL);

	parab_solver.soln_comp = LIQUID_COMP2;
	parab_solver.obst_comp = SOLID_COMP;
	parab_solver.ilower = ilower;
        parab_solver.iupper = iupper;
	parab_solver.dt = m_dt;
	parab_solver.order = 2;
	parab_solver.a = NULL;
	parab_solver.findStateAtCrossing = parab_find_state_at_crossing;
	parab_solver.first = first;
	parab_solver.set_solver_domain();
	first = NO;
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

void Incompress_Solver_Smooth_3D_Cartesian::vtk_plot_scalar(
        char *outname, const char* varname)
{
        std::vector<int> ph_index;
        int i,j,k,index;
        char dirname[256],filename[256];
        FILE *outfile;
        double coord_x,coord_y,coord_z,xmin,ymin,zmin;
        COMPONENT comp;
        int pointsx,pointsy,pointsz,num_points,num_cells,num_cell_list;
        int icoords[3],p_gmax[3];

        int ii,jj,kk;
        double ih,jh,kh;

        sprintf(filename, "%s/vtk/vtk.ts%s",outname,
                right_flush(front->step,7));
        if (pp_numnodes() > 1)
            sprintf(filename,"%s-nd%s",filename,right_flush(pp_mynode(),4));
        //cell-based liquid phase
        ph_index.clear();
        if (!create_directory(filename,NO))
        {
            printf("Cannot create directory %s\n",filename);
            clean_up(ERROR);
        }
        sprintf(filename,"%s/%s.vtk",filename,varname);
        outfile = fopen(filename,"w");
        fprintf(outfile,"# vtk DataFile Version 3.0\n");
        fprintf(outfile,"%s\n",varname);
        fprintf(outfile,"ASCII\n");
        fprintf(outfile,"DATASET UNSTRUCTURED_GRID\n");

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            if (ifluid_comp(cell_center[index].comp))
                ph_index.push_back(index);
        }

        pointsx = top_gmax[0] + 2;
        pointsy = top_gmax[1] + 2;
        pointsz = top_gmax[2] + 2;
        num_points = pointsx*pointsy*pointsz;

        num_cells = (int)ph_index.size();
        num_cell_list = 9*num_cells;
        fprintf(outfile,"POINTS %d double\n", num_points);

        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            coord_x = cell_center[index].m_coords[0] - top_h[0]/2.0;
            coord_y = cell_center[index].m_coords[1] - top_h[1]/2.0;
            coord_z = cell_center[index].m_coords[2] - top_h[2]/2.0;
            fprintf(outfile,"%f %f %f\n",coord_x,coord_y,coord_z);
        }
        for (i = 0; i <= top_gmax[0]; i++)
        for (j = 0; j <= top_gmax[1]; j++)
        {
            k = top_gmax[2];
            index = d_index3d(i,j,k,top_gmax);
            coord_x = cell_center[index].m_coords[0] - top_h[0]/2.0;
            coord_y = cell_center[index].m_coords[1] - top_h[1]/2.0;
            coord_z = cell_center[index].m_coords[2] + top_h[2]/2.0;
            fprintf(outfile,"%f %f %f\n",coord_x,coord_y,coord_z);
        }

        for (i = 0; i <= top_gmax[0]; i++)
        for (k = 0; k <= top_gmax[2]; k++)
        {
            j = top_gmax[1];
            index = d_index3d(i,j,k,top_gmax);
            coord_x = cell_center[index].m_coords[0] - top_h[0]/2.0;
            coord_y = cell_center[index].m_coords[1] + top_h[1]/2.0;
            coord_z = cell_center[index].m_coords[2] - top_h[2]/2.0;
            fprintf(outfile,"%f %f %f\n",coord_x,coord_y,coord_z);
        }

        for (j = 0; j <= top_gmax[1]; j++)
        for (k = 0; k <= top_gmax[2]; k++)
        {
            i = top_gmax[0];
            index = d_index3d(i,j,k,top_gmax);
            coord_x = cell_center[index].m_coords[0] + top_h[0]/2.0;
            coord_y = cell_center[index].m_coords[1] - top_h[1]/2.0;
            coord_z = cell_center[index].m_coords[2] - top_h[2]/2.0;
            fprintf(outfile,"%f %f %f\n",coord_x,coord_y,coord_z);
        }

        for (i = 0; i <= top_gmax[0]; i++)
        {
            j = top_gmax[1];
            k = top_gmax[2];
            index = d_index3d(i,j,k,top_gmax);
            coord_x = cell_center[index].m_coords[0] - top_h[0]/2.0;
            coord_y = cell_center[index].m_coords[1] + top_h[1]/2.0;
            coord_z = cell_center[index].m_coords[2] + top_h[2]/2.0;
            fprintf(outfile,"%f %f %f\n",coord_x,coord_y,coord_z);
        }

        for (j = 0; j <= top_gmax[1]; j++)
        {
            i = top_gmax[0];
            k = top_gmax[2];
            index = d_index3d(i,j,k,top_gmax);
            coord_x = cell_center[index].m_coords[0] + top_h[0]/2.0;
            coord_y = cell_center[index].m_coords[1] - top_h[1]/2.0;
            coord_z = cell_center[index].m_coords[2] + top_h[2]/2.0;
            fprintf(outfile,"%f %f %f\n",coord_x,coord_y,coord_z);
        }
        for (k = 0; k <= top_gmax[2]; k++)
        {
            i = top_gmax[0];
            j = top_gmax[1];
            index = d_index3d(i,j,k,top_gmax);
            coord_x = cell_center[index].m_coords[0] + top_h[0]/2.0;
            coord_y = cell_center[index].m_coords[1] + top_h[1]/2.0;
            coord_z = cell_center[index].m_coords[2] - top_h[2]/2.0;
            fprintf(outfile,"%f %f %f\n",coord_x,coord_y,coord_z);
        }

        i = top_gmax[0];
        j = top_gmax[1];
        k = top_gmax[2];
        index = d_index3d(i,j,k,top_gmax);
        coord_x = cell_center[index].m_coords[0] + top_h[0]/2.0;
        coord_y = cell_center[index].m_coords[1] + top_h[1]/2.0;
        coord_z = cell_center[index].m_coords[2] + top_h[2]/2.0;
        fprintf(outfile,"%f %f %f\n",coord_x,coord_y,coord_z);


        fprintf(outfile,"CELLS %i %i\n", num_cells,num_cell_list);
        for (i = 0; i < num_cells; i++)
        {
            int index0,index1,index2,index3,index4,index5,index6,index7;
            index = ph_index[i];
            icoords[0] = cell_center[index].icoords[0];
            icoords[1] = cell_center[index].icoords[1];
            icoords[2] = cell_center[index].icoords[2];
            index0 = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
            index1 =
                d_index3d(icoords[0]+1,icoords[1],icoords[2],top_gmax);
            index2 =
                d_index3d(icoords[0],icoords[1]+1,icoords[2],top_gmax);
            index3 =
                d_index3d(icoords[0]+1,icoords[1]+1,icoords[2],top_gmax);
            index4 =
                d_index3d(icoords[0],icoords[1],icoords[2]+1,top_gmax);
            index5 =
                d_index3d(icoords[0]+1,icoords[1],icoords[2]+1,top_gmax);
            index6 =
                d_index3d(icoords[0],icoords[1]+1,icoords[2]+1,top_gmax);
            index7 =
                d_index3d(icoords[0]+1,icoords[1]+1,icoords[2]+1,top_gmax);

            fprintf(outfile,"8 %i %i %i %i %i %i %i %i\n",
                index0,index1,index2,index3,index4,index5,index6,index7);
        }

        fprintf(outfile, "CELL_TYPES %i\n", num_cells);
        for (i = 0; i < num_cells; i++)
            fprintf(outfile,"11\n");

        fprintf(outfile, "CELL_DATA %i\n", num_cells);
        fprintf(outfile, "SCALARS %s double\n",varname);
        fprintf(outfile, "LOOKUP_TABLE default\n");
        for (i = 0; i < num_cells; i++)
        {
            index = ph_index[i];
            if(strcmp(varname,"pres") == 0)
              fprintf(outfile,"%f\n",field->pres[index]);
        }
        fclose(outfile);
}

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

void Incompress_Solver_Smooth_3D_Cartesian::setParallelVelocity(void)
{
        FILE *infile;
        int i,j,id,k,l,index,G_index;
        char fname[100];
        COMPONENT comp;
        double coords[MAXD];
        int size = (int)cell_center.size();
        int myid = pp_mynode();
        int numprocs = pp_numnodes();

        int G_icoords[MAXD],pp_icoords[MAXD],icoords[MAXD];
        int local_gmax[MAXD], global_gmax[MAXD];
        int G_size, L_size;
        PP_GRID *pp_grid = front->pp_grid;
        double *local_L = pp_grid->Zoom_grid.L;
        double *local_U = pp_grid->Zoom_grid.U;
        double *GU_buff,*GV_buff, *GW_buff, *U_buff, *V_buff, *W_buff;

        for (i = 0; i < dim; i++)
        {
            global_gmax[i] = pp_grid->Global_grid.gmax[i]-1;
            local_gmax[i] = pp_grid->Zoom_grid.gmax[i]-1;
        }
        FT_MakeGridIntfc(front);
        setDomain();
        G_size = 1;
        L_size = 1;
        for (i = 0; i < dim; i++)
        {
            G_size = G_size * (global_gmax[i]+1);
            L_size = L_size * (top_gmax[i]+1);
        }
        uni_array(&U_buff,L_size,sizeof(double));
        uni_array(&V_buff,L_size,sizeof(double));
        uni_array(&W_buff,L_size,sizeof(double));
        if (myid == 0)
        {
            uni_array(&GU_buff,G_size,sizeof(double));
            uni_array(&GV_buff,G_size,sizeof(double));
            uni_array(&GW_buff,G_size,sizeof(double));

	    if (getInitialState != NULL)
                (*setInitialVelocity)(comp,pp_grid->Global_grid.gmax,
				   GU_buff,GV_buff,GW_buff,
				   dim,iFparams);
            for (id = 0; id < numprocs; id++)
            {
                find_Cartesian_coordinates(id,pp_grid,pp_icoords);
                for (k = kmin; k <= kmax; ++k)
                for (j = jmin; j <= jmax; ++j)
                for (i = imin; i <= imax; ++i)
                {
                    icoords[0] = i;
                    icoords[1] = j;
		    icoords[2] = k;
                    G_icoords[0] = pp_icoords[0]*(local_gmax[0]+1)+icoords[0]-imin;
                    G_icoords[1] = pp_icoords[1]*(local_gmax[1]+1)+icoords[1]-jmin;
                    G_icoords[2] = pp_icoords[2]*(local_gmax[2]+1)+icoords[2]-kmin;
                    G_index = d_index(G_icoords,global_gmax,dim);
                    index = d_index(icoords,top_gmax,dim);
                    U_buff[index] = GU_buff[G_index];
                    V_buff[index] = GV_buff[G_index];
                    W_buff[index] = GW_buff[G_index];
                }
                if (id == 0)
                {
                    for (i = 0; i < L_size; i++)
                    {
                        field->vel[0][i] = U_buff[i];
                        field->vel[1][i] = V_buff[i];
                        field->vel[2][i] = W_buff[i];
                    }
                }
                else
                {
                    pp_send(1,(POINTER)(U_buff),sizeof(double)*L_size,id);
                    pp_send(2,(POINTER)(V_buff),sizeof(double)*L_size,id);
                    pp_send(3,(POINTER)(W_buff),sizeof(double)*L_size,id);
                }
            }
            FT_FreeThese(3,GU_buff,GV_buff,GW_buff);
        }
        else
        {
            pp_recv(1,0,(POINTER)(U_buff),sizeof(double)*L_size);
            pp_recv(2,0,(POINTER)(V_buff),sizeof(double)*L_size);
            pp_recv(3,0,(POINTER)(W_buff),sizeof(double)*L_size);
            for (i = 0; i < L_size; i++)
            {
                field->vel[0][i] = U_buff[i];
                field->vel[1][i] = V_buff[i];
                field->vel[2][i] = W_buff[i];
            }
        }


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
        FT_FreeThese(3,U_buff,V_buff,W_buff);

        computeGradientQ();
        copyMeshStates();
        setAdvectionDt();
}

void Incompress_Solver_Smooth_3D_Cartesian::extractFlowThroughVelocity()
{
	int index,index_nb,index_op;
	int i,j,k,l,nb,icoords[MAXD],icoords_nb[MAXD],icoords_op[MAXD];
	GRID_DIRECTION dir[3][2] = {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
	POINTER intfc_state;
        HYPER_SURF *hs;
	double crx_coords[MAXD];
	double vel_save,**vel = field->vel;
	double div;
	int status;

	/* Extract velocity for zero interior divergence */
	for (k = 0; k <= top_gmax[2]; ++k)
        for (j = 0; j <= top_gmax[1]; ++j)
        for (i = 0; i <= top_gmax[0]; ++i)
        {
            if (i >= imin && i <= imax && j >= jmin &&
                j <= jmax && k >= kmin && k <= kmax)
                continue;
            icoords[0] = i;
            icoords[1] = j;
            icoords[2] = k;
            index = d_index(icoords,top_gmax,dim);
            if (ifluid_comp(top_comp[index])) continue;
            for (l = 0; l < dim; ++l)
                icoords_nb[l] = icoords_op[l] = icoords[l];
            for (l = 0; l < dim; ++l)
            for (nb = 0; nb < 2; ++nb)
            {
                status = (*findStateAtCrossing)(front,icoords,dir[l][nb],
                                top_comp[index],&intfc_state,&hs,crx_coords);
                if (status == NO_PDE_BOUNDARY) continue;
                icoords_nb[l] = (nb == 0) ? icoords[l] - 1 : icoords[l] + 1;
                icoords_op[l] = (nb == 0) ? icoords[l] - 2 : icoords[l] + 2;
                index_nb = d_index(icoords_nb,top_gmax,dim);
                index_op = d_index(icoords_op,top_gmax,dim);
                if (!ifluid_comp(top_comp[index_nb]))
                    continue;
                vel[l][index] = 0.0;
                div = computeFieldPointDiv(icoords_nb,vel);
                vel[l][index] = (nb == 0) ? vel[l][index] - 2.0*top_h[l]*div
                                : vel[l][index] + 2.0*top_h[l]*div;
            }
        }

}	/* end extractFlowThroughVelocity */

void Incompress_Solver_Smooth_3D_Cartesian::updateComponent(void)
{
        int i,j,k,l,icoords[MAXD];
        int index;

	/*update the component of pressure on dual grid*/
        for (k = kmin; k <= kmax; ++k)
        for (j = jmin; j <= jmax; ++j)
        for (i = imin; i <= imax; ++i)
        {
            icoords[0] = i;
            icoords[1] = j;
            icoords[2] = k;
            index = d_index(icoords,top_gmax,dim);
            if (ifluid_comp(top_comp[index]))
            {
                int cl[MAXD], cu[MAXD];
                for (l = 0; l < dim; ++l)
                    cl[l] = icoords[l] - offset[l];
                for (int m = 0; m < 2; ++m)
                for (int n = 0; n < 2; ++n)
                for (int r = 0; r < 2; ++r)
                {
                    cu[0] = cl[0] + m;
                    cu[1] = cl[1] + n;
                    cu[2] = cl[2] + r;
		    if (cu[0]<cimin || cu[0]>cimax || cu[1]<cjmin 
			|| cu[1]>cjmax || cu[2]<ckmin || cu[2]>ckmax)
			continue;
                    int index_tmp = d_index(cu,ctop_gmax,dim);
		    if (!ifluid_comp(ctop_comp[index_tmp]))
		    	ctop_comp[index_tmp] = top_comp[index];
                }
            }
        }

        /*Set rho for boundary layer on computational grid*/
	if (field->rho == NULL)
	    return;
        for (k = kmin; k <= kmax; ++k)
        for (j = jmin; j <= jmax; ++j)
        for (i = imin; i <= imax; ++i)
        {
            icoords[0] = i;
            icoords[1] = j;
            icoords[2] = k;
            index = d_index(icoords,top_gmax,dim);
            if (!ifluid_comp(top_comp[index])&&!InsideSolid(icoords))
            {
                int cl[MAXD],cu[MAXD],indexl,indexu;
		boolean VelSet = NO;
		for (l = 0; l < dim && !VelSet; ++l)
		{
		    cl[l] = icoords[l]-1;
		    cu[l] = icoords[l]+1;
		    for (int n = -1; n < 2 && !VelSet; ++n)
                    for (int r = -1; r < 2 && !VelSet; ++r)
		    {
			cl[(l+1)%dim] = cu[(l+1)%dim] = icoords[(l+1)%dim]+n;
			cl[(l+2)%dim] = cu[(l+2)%dim] = icoords[(l+2)%dim]+r;
			indexl = d_index(cl,top_gmax,dim);
			if (ifluid_comp(top_comp[indexl]))
			{
			    field->rho[index] = field->rho[indexl];
			    VelSet = YES;
			    continue;
			}
			indexu = d_index(cu,top_gmax,dim);
			if (ifluid_comp(top_comp[indexu]))
                        {
                            field->rho[index] = field->rho[indexu];
                            VelSet = YES;
                            continue;
                        }
		    }
		}
            }
        }
}	/* end updateComponent */

void Incompress_Solver_Smooth_3D_Cartesian::computeVelDivergence()
{
	double *div_U = field->div_U;
	double **vel = field->vel;
	int i,j,k,index,icoords[MAXD];
	double Lnorm[3];

	/* Compute velocity divergence */
	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;
	    index  = d_index(icoords,top_gmax,dim);
	    if (!ifluid_comp(top_comp[index]))
		div_U[index] = 0.0;
	    div_U[index] = computeFieldPointDiv(icoords,vel);
	}
}	/* end computeVelDivergence */

void Incompress_Solver_Smooth_3D_Cartesian::computeVarIncrement(
	double *var_old,
	double *var_new,
	boolean use_dual_grid)
{
	int i,j,k,index,size;
	double mag;
	double Lnorm[3];

	if (use_dual_grid)
	{
	    ;	// To add dual grid computation.
	}
	else
	{
	    Lnorm[0] = Lnorm[1] = Lnorm[2] = 0.0;
	    size = (imax - imin + 1)*(jmax - jmin + 1)*(kmax - kmin + 1);
	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
		mag = fabs(var_new[index] - var_old[index]);
		Lnorm[0] += mag;
		Lnorm[1] += sqr(mag);
		if (Lnorm[2] < mag)
		{
		    Lnorm[2] = mag;
		}
            }
	    Lnorm[0] /= size;
	    Lnorm[1] = sqrt(Lnorm[1]/size);
	    (void) printf("L-1 norm = %20.14f  L-1/dt = %20.14f\n",
					Lnorm[0],Lnorm[0]/m_dt);
	    (void) printf("L-2 norm = %20.14f  L-2/dt = %20.14f\n",
					Lnorm[1],Lnorm[1]/m_dt);
	    (void) printf("L-I norm = %20.14f  L-I/dt = %20.14f\n",
					Lnorm[2],Lnorm[2]/m_dt);
	    Lnorm[0] = 0.0;
	    for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
	    {
		mag = fabs(var_new[index]);
		Lnorm[0] += mag;
	    }
	    Lnorm[0] /= size;
	    (void) printf("L-1 norm of old variable = %20.14f\n",Lnorm[0]);
	}
}	/* end computeVarIncrement */	

/***************************************************************
FronTier is a set of libraries that implements different types of 
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

#include "solver.h"

void PARABOLIC_SOLVER::solve(
	double *var_in,
	double *var_out)
{
	switch (dim)
	{
	case 1:
	    return solve1d(var_in,var_out);
	case 2:
	    return solve2d(var_in,var_out);
	case 3:
	    return solve3d(var_in,var_out);
	}
}	/* end solve */

void PARABOLIC_SOLVER::solve2d(
	double *var_in,
	double *var_out)
{
	int i,j,k,l,m,ic,icn,I,I_nb,icoords[MAXD];
	int gmin[MAXD],ipn[MAXD];
	double crx_coords[MAXD];
	double C_nb,lambda,eta,coeff,coeff_nb,rhs;
	COMPONENT comp;
	PETSc solver;
	double *x;
        int num_iter = 0;
	int size;
        double rel_residual = 0;
        int fr_crx_grid_seg;
	const GRID_DIRECTION dir[3][2] =
                {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
	HYPER_SURF *hs;
	POINTER state;

        start_clock("solve2d");

        for (i = 0; i < dim; ++i) gmin[i] = 0;

        size = iupper - ilower;

	solver.Create(ilower, iupper-1, 5, 5);
	for (i = imin; i <= imax; ++i)
	for (j = jmin; j <= jmax; ++j)
	{
	    icoords[0] = i;
	    icoords[1] = j;
	    ic = d_index(icoords,top_gmax,dim);
	    comp = top_comp[ic];
	    I = ij_to_I[i][j];
	    if (comp != soln_comp) 
	    	continue;
	    rhs = var_in[ic];
	    coeff = 1.0;
	    for (l = 0; l < dim; ++l)
	    {
                for (m = 0; m < 2; ++m)
                {
		    if (nu == NULL)
                        lambda = D*sub_dt/sqr(top_h[l]);
                    else
                        lambda = nu[ic]*sub_dt/sqr(top_h[l]);
                    if (next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax))
		    {
                    	icn = d_index(ipn,top_gmax,dim);
		    	I_nb = ij_to_I[ipn[0]][ipn[1]];
			eta = (a == NULL) ? 0.0 :
				0.25*(a[l][ic] + a[l][icn])*sub_dt/top_h[l];
		    }
		    else
		    {
		    	eta = 0.0;
			I_nb = -1;
		    }
		    fr_crx_grid_seg =(*findStateAtCrossing)(front,icoords,
				dir[l][m],comp,&state,&hs,crx_coords);

                    if (fr_crx_grid_seg == NO_PDE_BOUNDARY)
                    {
			if (nu != NULL)
                        {
                            lambda += nu[icn]*sub_dt/sqr(top_h[l]);
                            lambda *= 0.5;
                        }
			coeff_nb = -lambda;
                        coeff += lambda;
		        if (m == 0) coeff_nb -= eta;
                        else coeff_nb += eta;
			solver.Add_A(I,I_nb,coeff_nb);
                    }
		    else if (fr_crx_grid_seg == DIRICHLET_PDE_BOUNDARY ||
		    	     fr_crx_grid_seg == MOVING_BOUNDARY)
		    {
			coeff_nb = -lambda;
			if (wave_type(hs) == DIRICHLET_BOUNDARY &&
			    boundary_state_function(hs) &&
                            strcmp(boundary_state_function_name(hs),
                            "flowThroughBoundaryState") == 0)
			{
			    C_nb = var_in[ic];
			    if (m == 0) coeff -= eta;
                            else coeff += eta;
			}
			else
			{
                            coeff += lambda;
			    C_nb = (*getStateVarFunc)(state);
			    rhs -= coeff_nb*C_nb;
                            if (m == 0) rhs += eta*C_nb;
                            else rhs -= eta*C_nb;
			}
		    }
		    else if (fr_crx_grid_seg != NEUMANN_PDE_BOUNDARY)
			coeff += lambda;
	
                }
	    }
	    solver.Add_A(I,I,coeff);
	    solver.Add_b(I,rhs);
	}
	start_clock("petsc_solve");
	solver.SetMaxIter(500);   
	solver.SetTol(1e-10);   
	solver.Solve();

	if (debugging("PETSc"))
	{
	    (void) printf("C_CARTESIAN::computeAdvectionImplicit: "
	       		"num_iter = %d, rel_residual = %g \n", 
			num_iter, rel_residual);
	}

	FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));
	solver.Get_x(x);
	stop_clock("petsc_solve");

        for (i = imin; i <= imax; i++)
        for (j = jmin; j <= jmax; j++)
        {
	    icoords[0] = i;
	    icoords[1] = j;
	    I = ij_to_I[i][j];
	    ic = d_index(icoords,top_gmax,dim);
	    comp = top_comp[ic];
	    if (comp == soln_comp)
	    	array[ic] = x[I-ilower];
	    else if (comp == obst_comp)
                array[ic] = var_obst;
	}
	FT_ParallelExchGridArrayBuffer(array,front,NULL);
        for (i = 0; i <= top_gmax[0]; ++i)
        for (j = 0; j <= top_gmax[1]; ++j)
        {
	    ic = d_index2d(i,j,top_gmax);
            var_out[ic] = array[ic];
	}
	FT_FreeThese(1,x);
        stop_clock("solve2d");
}	/* end solve2d */

void PARABOLIC_SOLVER::set_solver_domain(void)
{
	RECT_GRID *rgr = &topological_grid(front->grid_intfc);
        struct Table *T = table_of_interface(front->grid_intfc);
        int *lbuf = front->rect_grid->lbuf;
        int *ubuf = front->rect_grid->ubuf;
	int i;
	dim = Dimension(front->interf);
        top_comp = T->components;
        top_gmax = rgr->gmax;
	top_h = rgr->h;
	top_L = rgr->L;

	if (first)
	{
	    switch(dim)
	    {
	    case 1:
		imin = (lbuf[0] == 0) ? 1 : lbuf[0];
		imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
		FT_VectorMemoryAlloc((POINTER*)&array,top_gmax[0]+1,FLOAT);
		break;
	    case 2:
		imin = (lbuf[0] == 0) ? 1 : lbuf[0];
	    	jmin = (lbuf[1] == 0) ? 1 : lbuf[1];
		imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
	    	jmax = (ubuf[1] == 0) ? top_gmax[1] - 1 : top_gmax[1] - ubuf[1];
		FT_VectorMemoryAlloc((POINTER*)&array,
                                (top_gmax[0]+1)*(top_gmax[1]+1),FLOAT);
		break;
	    case 3:
		imin = (lbuf[0] == 0) ? 1 : lbuf[0];
	    	jmin = (lbuf[1] == 0) ? 1 : lbuf[1];
		kmin = (lbuf[2] == 0) ? 1 : lbuf[2];
		imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
	    	jmax = (ubuf[1] == 0) ? top_gmax[1] - 1 : top_gmax[1] - ubuf[1];
		kmax = (ubuf[2] == 0) ? top_gmax[2] - 1 : top_gmax[2] - ubuf[2];
		FT_VectorMemoryAlloc((POINTER*)&array,
                        (top_gmax[0]+1)*(top_gmax[1]+1)*(top_gmax[2]+1),FLOAT);
		break;
	    }
	    array_size = 1;
	    for (i = 0; i < dim; ++i)
	    	array_size *= (top_gmax[i] + 1);
	}
}	/* end set_solver_domain */

PARABOLIC_SOLVER::PARABOLIC_SOLVER(Front &front):front(&front)
{
}

void PARABOLIC_SOLVER::solve1d(
	double *var_in,
	double *var_out)
{
	int i,l,m,ic,icn,I,I_nb,icoords[MAXD];
	int gmin[MAXD],ipn[MAXD];
	double crx_coords[MAXD];
	double C_nb,lambda,coeff,coeff_nb,rhs;
	COMPONENT comp;
	PETSc solver;
	double *x;
        int num_iter = 0;
	int size;
        double rel_residual = 0;
        int fr_crx_grid_seg;
	const GRID_DIRECTION dir[3][2] =
                {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
	HYPER_SURF *hs;
        POINTER state;

        start_clock("solve1d");

        for (i = 0; i < dim; ++i) gmin[i] = 0;

        size = iupper - ilower;

	solver.Create(ilower, iupper-1, 3, 3);
	for (i = imin; i <= imax; ++i)
	{
	    icoords[0] = i;
	    ic = d_index(icoords,top_gmax,dim);
	    comp = top_comp[ic];
	    I = i_to_I[i];
	    if (comp != soln_comp) 
	    	continue;
	    rhs = var_in[ic];
	    coeff = 1.0;
	    for (l = 0; l < dim; ++l)
	    {
                for (m = 0; m < 2; ++m)
                {
		    if (nu == NULL)
                        lambda = D*sub_dt/sqr(top_h[l]);
                    else
                        lambda = nu[ic]*sub_dt/sqr(top_h[l]);
                    next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax);
                    icn = d_index(ipn,top_gmax,dim);
		    I_nb = i_to_I[ipn[0]];
                    fr_crx_grid_seg = (*findStateAtCrossing)(front,icoords,
				dir[l][m],comp,&state,&hs,crx_coords); 

                    if (fr_crx_grid_seg == NO_PDE_BOUNDARY)
                    {
			if (nu != NULL)
                        {
                            lambda += nu[icn]*sub_dt/sqr(top_h[l]);
                            lambda *= 0.5;
                        }
                        coeff_nb = -lambda;
                        coeff += lambda;
			solver.Add_A(I,I_nb,coeff_nb);
                    }
		    else if (fr_crx_grid_seg == DIRICHLET_PDE_BOUNDARY ||
		    	     fr_crx_grid_seg == MOVING_BOUNDARY)
		    {
			coeff_nb = -lambda;
			if (wave_type(hs) == DIRICHLET_BOUNDARY &&
                            boundary_state_function(hs) &&
                            strcmp(boundary_state_function_name(hs),
                            "flowThroughBoundaryState") == 0)
                            C_nb = var_in[ic];
                        else
			{
                            coeff += lambda;
                            C_nb = (*getStateVarFunc)(state);
                            rhs -= coeff_nb*C_nb;
			}
		    }
		    else if (fr_crx_grid_seg != NEUMANN_PDE_BOUNDARY)
                        coeff += lambda;

                }
	    }
	    solver.Add_A(I,I,coeff);
	    solver.Add_b(I,rhs);
	}
	start_clock("petsc_solve");
	solver.SetMaxIter(500);   
	solver.SetTol(1e-10);   
	solver.Solve();

	if (debugging("PETSc"))
	{
	    (void) printf("C_CARTESIAN::computeAdvectionImplicit: "
	       		"num_iter = %d, rel_residual = %g \n", 
			num_iter, rel_residual);
	}

	FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));
	solver.Get_x(x);
	stop_clock("petsc_solve");

        for (i = imin; i <= imax; i++)
        {
	    icoords[0] = i;
	    I = i_to_I[i];
	    ic = d_index(icoords,top_gmax,dim);
	    comp = top_comp[ic];
	    if (comp == soln_comp)
	    	array[ic] = x[I-ilower];
	    else if (comp == obst_comp)
                array[ic] = var_obst;
	}
	FT_ParallelExchGridArrayBuffer(array,front,NULL);
        for (i = 0; i <= top_gmax[0]; ++i)
        {
	    icoords[0] = i;
	    ic = d_index(icoords,top_gmax,dim);
            var_out[ic] = array[ic];
	}
	FT_FreeThese(1,x);
        stop_clock("solve1d");
}	/* end solve1d */

void PARABOLIC_SOLVER::solve3d(
	double *var_in,
	double *var_out)
{
	int i,j,k,l,m,ic,icn,I,I_nb,icoords[MAXD];
	int gmin[MAXD],ipn[MAXD];
	double crx_coords[MAXD];
	double C_nb,lambda,eta,coeff,coeff_nb,rhs;
	COMPONENT comp;
	PETSc solver;
	double *x;
        int num_iter = 0;
	int size;
        double rel_residual = 0;
        int fr_crx_grid_seg;
	const GRID_DIRECTION dir[3][2] =
                {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
	HYPER_SURF *hs;
	POINTER state;

        start_clock("solve3d");

        for (i = 0; i < dim; ++i) gmin[i] = 0;

        size = iupper - ilower;

	solver.Create(ilower, iupper-1, 7, 7);
	for (i = imin; i <= imax; ++i)
	for (j = jmin; j <= jmax; ++j)
	for (k = kmin; k <= kmax; ++k)
	{
	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;
	    ic = d_index(icoords,top_gmax,dim);
	    comp = top_comp[ic];
	    I = ijk_to_I[i][j][k];
	    if (comp != soln_comp) 
	    	continue;
	    rhs = var_in[ic];
	    coeff = 1.0;
	    for (l = 0; l < dim; ++l)
	    {
                for (m = 0; m < 2; ++m)
                {
		    if (nu == NULL)
                        lambda = D*sub_dt/sqr(top_h[l]);
                    else
                        lambda = nu[ic]*sub_dt/sqr(top_h[l]);
                    if (next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax))
		    {
                    	icn = d_index(ipn,top_gmax,dim);
		    	I_nb = ij_to_I[ipn[0]][ipn[1]];
			eta = (a == NULL) ? 0.0 :
				0.25*(a[l][ic] + a[l][icn])*sub_dt/top_h[l];
		    }
		    else
		    {
		    	eta = 0.0;
			I_nb = -1;
		    }
                    fr_crx_grid_seg =(*findStateAtCrossing)(front,icoords,
				dir[l][m],comp,&state,&hs,crx_coords);
		    if (m == 0) coeff_nb -= eta;
		    else coeff_nb += eta;

                    if (fr_crx_grid_seg == NO_PDE_BOUNDARY)
                    {
			if (nu != NULL)
			{
			    lambda += nu[icn]*sub_dt/sqr(top_h[l]);
			    lambda *= 0.5;
			}
			coeff_nb = -lambda;
			coeff += lambda;
			if (m == 0) coeff_nb -= eta;
                        else coeff_nb += eta;
			solver.Add_A(I,I_nb,coeff_nb);
                    }
		    else if (fr_crx_grid_seg == DIRICHLET_PDE_BOUNDARY ||
		    	     fr_crx_grid_seg == MOVING_BOUNDARY)
		    {
			coeff_nb = -lambda;
			if (wave_type(hs) == DIRICHLET_BOUNDARY &&
                            boundary_state_function(hs) &&
                            strcmp(boundary_state_function_name(hs),
                            "flowThroughBoundaryState") == 0)
			{
                            C_nb = var_in[ic];
			    if (m == 0) coeff -= eta;
                            else coeff += eta;
			}
                        else
			{
			    coeff += lambda;
                            C_nb = (*getStateVarFunc)(state);
			    rhs -= coeff_nb*C_nb;
			    if (m == 0) rhs += eta*C_nb;
                            else rhs -= eta*C_nb;
			}
		    }
		    else if (fr_crx_grid_seg != NEUMANN_PDE_BOUNDARY)
			coeff += lambda;
                }
	    }
	    solver.Add_A(I,I,coeff);
	    solver.Add_b(I,rhs);
	}
	start_clock("petsc_solve");
	solver.SetMaxIter(500);   
	solver.SetTol(1e-10);   
	solver.Solve();

	if (debugging("PETSc"))
	{
	    (void) printf("C_CARTESIAN::computeAdvectionImplicit: "
	       		"num_iter = %d, rel_residual = %g \n", 
			num_iter, rel_residual);
	}

	FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));
	solver.Get_x(x);
	stop_clock("petsc_solve");

        for (i = imin; i <= imax; i++)
        for (j = jmin; j <= jmax; j++)
        for (k = kmin; k <= kmax; k++)
        {
	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;
	    I = ijk_to_I[i][j][k];
	    ic = d_index(icoords,top_gmax,dim);
	    comp = top_comp[ic];
	    if (comp == soln_comp)
	    	array[ic] = x[I-ilower];
	    else if (comp == obst_comp)
                array[ic] = var_obst;
	}
	FT_ParallelExchGridArrayBuffer(array,front,NULL);
        for (i = 0; i <= top_gmax[0]; ++i)
        for (j = 0; j <= top_gmax[1]; ++j)
        for (k = 0; k <= top_gmax[2]; ++k)
        {
	    ic = d_index3d(i,j,k,top_gmax);
	    var_out[ic] = array[ic];
	}
	FT_FreeThese(1,x);
        stop_clock("solve3d");
}	/* end solve3d */

void PARABOLIC_SOLVER::solveIM()
{
	static double *var_mid;

	bdry_influx = 0.0;
	if (var_mid == NULL)
	{
	    if (order > 1)
	    	FT_VectorMemoryAlloc((POINTER*)&var_mid,array_size,
				sizeof(double));
	}
	switch (order)
	{
	case 1:
	    sub_dt = dt;
	    solve(var,soln);
	    bdry_influx = compBdryFlux(soln,sub_dt);
	    break;
	case 2:
	    sub_dt = dt;
	    solve(var,var_mid);
	    bdry_influx = compBdryFlux(var_mid,0.5*dt);
	    addMeshState(var_mid,0.5,var,0.5,var_mid);
	    sub_dt = 0.5*dt;
	    solve(var_mid,soln);
	    bdry_influx += compBdryFlux(soln,0.5*dt);
	    break;
	case 4:
	default:
	    (void)printf("ERROR: %d-th order RK method not implemented\n",
					order);
	    clean_up(ERROR);
	}
}	/* end solveIM */

void PARABOLIC_SOLVER::addMeshState(
	double *ans,
	double C1,
	double *var1,
	double C2,
	double *var2)
{
	int i;
	for (i = 0; i < array_size; ++i)
	    ans[i] = C1*var1[i] + C2*var2[i];
}	/* end addMeshState */

double PARABOLIC_SOLVER::checkSolver(
	int *icoords,
	boolean print_details,
	double *var_in)
{
	int i,l,m;
	int comp;
	double w[2];
	int id0,index_nb;
	double dw[2],coefs[2],lhs,rhs;
	GRID_DIRECTION dir[3][2] = {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
	HYPER_SURF *hs;
        double crx_coords[MAXD];
        int status;
        POINTER intfc_state;
	int icnb[MAXD];
	double denom = 0.0;

	if (print_details)
	{
	    (void) printf("\nEntering checkSolver()\n");
	    (void) printf("icoords = ");
	    for (i = 0; i < dim; ++i)
	    	(void) printf("%d ",icoords[i]);
	    (void) printf("\n");
	    (void) printf("nu = %p  D = %20.14f\n",(void*)nu,D);
	}

	id0 = d_index(icoords,top_gmax,dim);
	comp = top_comp[id0];
	lhs = 0;
	for (l = 0; l < dim; ++l)
	{
	    if (print_details)
	    	(void) printf("Direction %d:\n",l);
	    for (i = 0; i < dim; ++i)
		icnb[i] = icoords[i];
	    for (m = 0; m < 2; ++m)
	    {
		status = (*findStateAtCrossing)(front,icoords,dir[l][m],comp,
                                &intfc_state,&hs,crx_coords);
		icnb[l] = (m == 0) ? icoords[l] - 1 : icoords[l] + 1;
		index_nb = d_index(icnb,top_gmax,dim);
		switch(status)
		{
		case NO_PDE_BOUNDARY:
	    	    if (print_details)
		    	(void) printf("Side %d NO_PDE_BOUNDARY\n",m);
		    if (nu == NULL)
		    	coefs[m] = D*sub_dt;
		    else
		    	coefs[m] = 0.5*(nu[id0] + nu[index_nb])*sub_dt;
		    w[0] = array[id0];
		    w[1] = array[index_nb];
		    break;
		case DIRICHLET_PDE_BOUNDARY:
	    	    if (print_details)
		    	(void) printf("Side %d DIRICHLET_PDE_BOUNDARY\n",m);
		    if (nu == NULL)
		    	coefs[m] = D*sub_dt;
		    else
		    	coefs[m] = nu[id0]*sub_dt;
		    w[0] = array[id0];
		    w[1] = (*getStateVarFunc)(intfc_state);
	    	    if (print_details)
		    	printf("intfc_state = %20.14f\n",w[1]);
		    break;
		default:
		    (void) printf("Side %d Unknown BOUNDARY\n",m);
		    clean_up(ERROR);
		}
		dw[m] = (w[1] - w[0])/top_h[l];
		if (denom < fabs(coefs[m]*dw[m]/2.0/top_h[l]))
		    denom = fabs(coefs[m]*dw[m]/2.0/top_h[l]);
	    }
	    if (print_details)
	    {
	    	(void) printf("Coefs: %20.14f %20.14f\n",coefs[0],coefs[1]);
	    	(void) printf("C*dw: %20.14f %20.14f\n",coefs[0]*dw[0],
					coefs[1]*dw[1]);
	    }
	    lhs += (coefs[1]*dw[1] + coefs[0]*dw[0])/top_h[l];
	}
	if (source != NULL)
	    rhs = source[id0];
	rhs = array[id0] - (var_in[id0]);
	if (print_details)
	{
	    double error;
	    if (fabs(lhs-rhs) == 0.0 && denom == 0.0)
		error = 0.0;
	    else
		error = fabs(lhs-rhs)/denom;
	    (void) printf("Solution = %20.14f\n",array[id0]);
	    (void) printf("LHS = %20.14f  RHS = %20.14f\n",lhs,rhs);
	    (void) printf("LHS - RHS = %20.14f\n",lhs-rhs);
	    (void) printf("Relative error = %20.14g\n",error);
	    (void) printf("Leaving checkSolver()\n\n");
	}
	return fabs(lhs-rhs)/denom;
}	/* end checkSolver */

double PARABOLIC_SOLVER::compBdryFlux(
	double *var,
	double dt)
{
	int i,j,k,ic,icoords[MAXD];
	double S0,Snb;
	double influx;
	int comp;
	double crx_coords[MAXD];
	int fr_crx_grid_seg;
	HYPER_SURF *hs;
        POINTER state;

	influx = 0.0;
	switch (dim)
	{
	case 1:
	    icoords[0] = imin;
	    ic = d_index(icoords,top_gmax,dim);
	    comp = top_comp[ic];
	    if (comp == soln_comp)
	    {
		S0 = var[ic];
		fr_crx_grid_seg = (*findStateAtCrossing)(front,
                                        icoords,WEST,comp,&state,&hs,
                                        crx_coords);
		if (fr_crx_grid_seg)
		{
		    Snb = (*getStateVarFunc)(state);
		    influx += D*dt*(Snb - S0)/top_h[0];
		}
	    }
	    icoords[0] = imax;
	    ic = d_index(icoords,top_gmax,dim);
	    comp = top_comp[ic];
	    if (comp == soln_comp)
	    {
		S0 = var[ic];
		fr_crx_grid_seg = (*findStateAtCrossing)(front,
                                        icoords,EAST,comp,&state,&hs,
                                        crx_coords);
		if (fr_crx_grid_seg)
		{
		    Snb = (*getStateVarFunc)(state);
		    influx += D*dt*(Snb - S0)/top_h[0];
		}
	    }
	    break;
	case 2:
	    if (FT_RectBoundaryType(front,1,0) == DIRICHLET_BOUNDARY)
	    {
	      for (i = imin; i <= imax; ++i)
	      {
	    	icoords[0] = i;
	    	icoords[1] = jmin;
	    	ic = d_index(icoords,top_gmax,dim);
	    	comp = top_comp[ic];
	    	if (comp == soln_comp)
	    	{
		    S0 = var[ic];
		    fr_crx_grid_seg = (*findStateAtCrossing)(front,
                                        icoords,SOUTH,comp,&state,&hs,
                                        crx_coords);
		    if (fr_crx_grid_seg)
		    {
			Snb = (*getStateVarFunc)(state);
		    	influx += D*dt*(Snb - S0)/top_h[1]*top_h[0];
		    }
	    	}
	      }
	    }
	    if (FT_RectBoundaryType(front,1,1) == DIRICHLET_BOUNDARY)
	    {
	      for (i = imin; i <= imax; ++i)
	      {
	    	icoords[0] = i;
	    	icoords[1] = jmax;
	    	ic = d_index(icoords,top_gmax,dim);
	    	comp = top_comp[ic];
	    	if (comp == soln_comp)
	    	{
		    S0 = var[ic];
		    fr_crx_grid_seg = (*findStateAtCrossing)(front,
                                        icoords,NORTH,comp,&state,&hs,
                                        crx_coords);
		    if (fr_crx_grid_seg)
		    {
			Snb = (*getStateVarFunc)(state);
		    	influx += D*dt*(Snb - S0)/top_h[1]*top_h[0];
		    }
	    	}
	      }
	    }
	    if (FT_RectBoundaryType(front,0,0) == DIRICHLET_BOUNDARY)
	    {
	      for (j = jmin; j <= jmax; ++j)
	      {
	    	icoords[1] = j;
	    	icoords[0] = imin;
	    	ic = d_index(icoords,top_gmax,dim);
	    	comp = top_comp[ic];
	    	if (comp == soln_comp)
	    	{
		    S0 = var[ic];
		    fr_crx_grid_seg = (*findStateAtCrossing)(front,
                                        icoords,WEST,comp,&state,&hs,
                                        crx_coords);
		    if (fr_crx_grid_seg)
		    {
			Snb = (*getStateVarFunc)(state);
		    	influx += D*dt*(Snb - S0)/top_h[0]*top_h[1];
		    }
	    	}
	      }
	    }
	    if (FT_RectBoundaryType(front,0,1) == DIRICHLET_BOUNDARY)
	    {
	      for (j = jmin; j <= jmax; ++j)
	      {
	    	icoords[1] = j;
	    	icoords[0] = imax;
	    	ic = d_index(icoords,top_gmax,dim);
	    	comp = top_comp[ic];
	    	if (comp == soln_comp)
	    	{
		    S0 = var[ic];
		    fr_crx_grid_seg = (*findStateAtCrossing)(front,
                                        icoords,EAST,comp,&state,&hs,
                                        crx_coords);
		    if (fr_crx_grid_seg)
		    {
			Snb = (*getStateVarFunc)(state);
		    	influx += D*dt*(Snb - S0)/top_h[0]*top_h[1];
		    }
	    	}
	      }
	    }
	    break;
	case 3:
	    for (i = imin; i <= imax; ++i)
	    for (j = jmin; j <= jmax; ++j)
	    {
	    	icoords[0] = i;
	    	icoords[1] = j;
	    	icoords[2] = kmin;
	    	ic = d_index(icoords,top_gmax,dim);
	    	comp = top_comp[ic];
	    	if (comp == soln_comp)
	    	{
		    S0 = var[ic];
		    fr_crx_grid_seg = FT_StateVarAtGridCrossing(front,
                                        icoords,LOWER,comp,getStateVarFunc,
                                        &Snb,crx_coords);
		    if (fr_crx_grid_seg)
		    {
		    	influx += D*dt*(Snb - S0)/top_h[2]*top_h[0]*top_h[1];
		    }
	    	}
	    	icoords[2] = kmax;
	    	ic = d_index(icoords,top_gmax,dim);
	    	comp = top_comp[ic];
	    	if (comp == soln_comp)
	    	{
		    S0 = var[ic];
		    fr_crx_grid_seg = FT_StateVarAtGridCrossing(front,
                                        icoords,UPPER,comp,getStateVarFunc,
                                        &Snb,crx_coords);
		    if (fr_crx_grid_seg)
		    {
		    	influx += D*dt*(Snb - S0)/top_h[2]*top_h[0]*top_h[1];
		    }
	    	}
	    }
	    for (j = jmin; j <= jmax; ++j)
	    for (k = kmin; k <= kmax; ++k)
	    {
	    	icoords[1] = j;
	    	icoords[2] = k;
	    	icoords[0] = imin;
	    	ic = d_index(icoords,top_gmax,dim);
	    	comp = top_comp[ic];
	    	if (comp == soln_comp)
	    	{
		    S0 = var[ic];
		    fr_crx_grid_seg = FT_StateVarAtGridCrossing(front,
                                        icoords,WEST,comp,getStateVarFunc,
                                        &Snb,crx_coords);
		    if (fr_crx_grid_seg)
		    {
		    	influx += D*dt*(Snb - S0)/top_h[2]*top_h[0]*top_h[1];
		    }
	    	}
	    	icoords[0] = imax;
	    	ic = d_index(icoords,top_gmax,dim);
	    	comp = top_comp[ic];
	    	if (comp == soln_comp)
	    	{
		    S0 = var[ic];
		    fr_crx_grid_seg = FT_StateVarAtGridCrossing(front,
                                        icoords,EAST,comp,getStateVarFunc,
                                        &Snb,crx_coords);
		    if (fr_crx_grid_seg)
		    {
		    	influx += D*dt*(Snb - S0)/top_h[2]*top_h[0]*top_h[1];
		    }
	    	}
	    }
	    for (i = imin; i <= imax; ++i)
	    for (k = kmin; k <= kmax; ++k)
	    {
	    	icoords[0] = i;
	    	icoords[2] = k;
	    	icoords[1] = jmin;
	    	ic = d_index(icoords,top_gmax,dim);
	    	comp = top_comp[ic];
	    	if (comp == soln_comp)
	    	{
		    S0 = var[ic];
		    fr_crx_grid_seg = FT_StateVarAtGridCrossing(front,
                                        icoords,SOUTH,comp,getStateVarFunc,
                                        &Snb,crx_coords);
		    if (fr_crx_grid_seg)
		    {
		    	influx += D*dt*(Snb - S0)/top_h[2]*top_h[0]*top_h[1];
		    }
	    	}
	    	icoords[1] = jmax;
	    	ic = d_index(icoords,top_gmax,dim);
	    	comp = top_comp[ic];
	    	if (comp == soln_comp)
	    	{
		    S0 = var[ic];
		    fr_crx_grid_seg = FT_StateVarAtGridCrossing(front,
                                        icoords,NORTH,comp,getStateVarFunc,
                                        &Snb,crx_coords);
		    if (fr_crx_grid_seg)
		    {
		    	influx += D*dt*(Snb - S0)/top_h[2]*top_h[0]*top_h[1];
		    }
	    	}
	    }
	    break;
	default:
	    (void) printf("Unknown dimension %d\n",dim);
	    clean_up(ERROR);
	}
	return influx;
}	/* end compBdryFlux */

void PARABOLIC_SOLVER::solveEX()
{
	int i,j,k,l,m,ic,icn,icoords[MAXD];
	int gmin[MAXD],ipn[MAXD];
	double crx_coords[MAXD];
	double crx_coords_old[MAXD];
	double solute,solute_nb[2],dgrad[MAXD],grad_plus[MAXD],grad_minus[MAXD];
	double coef;
	COMPONENT comp;
	HYPER_SURF *hs;
	POINTER state;
	int fr_crx_grid_seg;
	const GRID_DIRECTION dir[3][2] = 
		{{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
	double coords[MAXD];
	double v[MAXD],v_plus[MAXD],v_minus[MAXD];
	double eta_new,eta_old; // cell irregularity ratio
	double kappa[3][2];	// side regularity number:
				// 1: for regular; 0: for irregular
	double xl_new,xr_new;		// TMP for 1d test.
	double xl_old,xr_old;		// TMP for 1d test.
	double v_new,v_old;		// TMP for 1d test.
	double sink;

	if (debugging("trace")) 
	    printf("Entering solveEX()\n");
	start_clock("solveEX");

	bdry_influx = compBdryFlux(var,dt);

	coef = D*dt;
	for (i = 0; i < dim; ++i) gmin[i] = 0;

	switch (dim)
	{
        case 1:
            for (i = imin; i <= imax; ++i)
            {
	    	icoords[0] = i;
	    	ic = d_index1d(i,top_gmax);
	    	comp = top_comp[ic];
	    	if (comp != soln_comp) 
	    	    continue;
                array[ic] = solute = var[ic];
		for (l = 0; l < dim; ++l)
		{
		    v[l] = 0.0;
		    v_plus[l] = 0.0;
		    v_minus[l] = 0.0;
		}
		if (a != NULL)
		{
		    for (l = 0; l < dim; ++l)
		    {
			v[l] = a[l][ic];
			v_plus[l] = std::max(0.0,v[l]);
			v_minus[l] = std::min(0.0,v[l]);
		    }
		}
		for (l = 0; l < dim; ++l)
		{
	            dgrad[l] = 0.0;
	 	    grad_plus[l] = 0.0;
		    grad_minus[l] = 0.0;

                    for (m = 0; m < 2; ++m)
                    {
                        fr_crx_grid_seg = FT_StateVarAtGridCrossing(front,
                                        icoords,dir[l][m],comp,getStateVarFunc,
                                        &solute_nb[m],crx_coords);
                        if (!fr_crx_grid_seg)
                        {
                            next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax);
                            icn = d_index1d(ipn[0],top_gmax);
			    solute_nb[m] = var[icn];
                        }
                        dgrad[l] += (solute_nb[m] - solute)/top_h[l];
                    }
		    grad_plus[l] = (solute_nb[1] - solute)/top_h[l];
		    grad_minus[l] = (solute - solute_nb[0])/top_h[l];

		    array[ic] += coef*dgrad[l]/top_h[l] - 
			dt*(v_plus[l]*grad_minus[l]+v_minus[l]*grad_plus[l]);
		}
            }
            break;
	case 2:
	    for (i = imin; i <= imax; ++i)
	    for (j = jmin; j <= jmax; ++j)
	    {
	    	icoords[0] = i;
	    	icoords[1] = j;
	    	ic = d_index2d(i,j,top_gmax);
	    	comp = top_comp[ic];
	    	if (comp != soln_comp) 
	    	    continue;
                array[ic] = solute = var[ic];
		for (l = 0; l < dim; ++l)
		{
		    v[l] = 0.0;
		    v_plus[l] = 0.0;
		    v_minus[l] = 0.0;
		}
		if (a != NULL)
		{
		    for (l = 0; l < dim; ++l)
		    {
			v[l] = a[l][ic];
			v_plus[l] = std::max(0.0,v[l]);
			v_minus[l] = std::min(0.0,v[l]);
		    }
		}
		for (l = 0; l < dim; ++l)
		{
	            dgrad[l] = 0.0;
	 	    grad_plus[l] = 0.0;
		    grad_minus[l] = 0.0;

                    for (m = 0; m < 2; ++m)
                    {
                        fr_crx_grid_seg = FT_StateVarAtGridCrossing(front,
                                        icoords,dir[l][m],comp,getStateVarFunc,
                                        &solute_nb[m],crx_coords);
                        if (!fr_crx_grid_seg)
                        {
                            next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax);
                            icn = d_index2d(ipn[0],ipn[1],top_gmax);
			    solute_nb[m] = var[icn];
                        }
                        dgrad[l] += (solute_nb[m] - solute)/top_h[l];
                    }
		    grad_plus[l] = (solute_nb[1] - solute)/top_h[l];
		    grad_minus[l] = (solute - solute_nb[0])/top_h[l];

		    array[ic] += coef*dgrad[l]/top_h[l] - 
			dt*(v_plus[l]*grad_minus[l]+v_minus[l]*grad_plus[l]);
		}
	    }
	    break;
	case 3:
	    for (i = imin; i <= imax; ++i)
	    for (j = jmin; j <= jmax; ++j)
	    for (k = kmin; k <= kmax; ++k)
	    {
	    	icoords[0] = i;
	    	icoords[1] = j;
	    	icoords[2] = k;
	    	ic = d_index3d(i,j,k,top_gmax);
	    	comp = top_comp[ic];
	    	if (comp != soln_comp) 
	    	    continue;
                array[ic] = solute = var[ic];
		for (l = 0; l < dim; ++l)
		{
		    v[l] = 0.0;
		    v_plus[l] = 0.0;
		    v_minus[l] = 0.0;
		}
		if (a != NULL)
		{
		    for ( l = 0; l < dim; ++l)
		    {	
			v[l] = a[l][ic];
			v_plus[l] = std::max(0.0,v[l]);
			v_minus[l] = std::min(0.0,v[l]);
		    }
		}
		for (l = 0; l < dim; ++l)
		{
	            dgrad[l] = 0.0;
	            grad_plus[l] = 0.0;
		    grad_minus[l] = 0.0;

                    for (m = 0; m < 2; ++m)
                    {
                        fr_crx_grid_seg = FT_StateVarAtGridCrossing(front,
                                        icoords,dir[l][m],comp,getStateVarFunc,
                                        &solute_nb[m],crx_coords);
                        if (!fr_crx_grid_seg)
                        {
                            next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax);
                            icn = d_index3d(ipn[0],ipn[1],ipn[2],top_gmax);
                	    solute_nb[m] = var[icn];
                        }
                        dgrad[l] += (solute_nb[m] - solute)/top_h[l];
                    }

		    grad_plus[l] = (solute_nb[1] - solute)/top_h[l];
		    grad_minus[l] = (solute - solute_nb[0])/top_h[l];

		    array[ic] += coef*dgrad[l]/top_h[l] - 
			dt*(v_plus[l]*grad_minus[l]+v_minus[l]*grad_plus[l]);
		}
	    }
	    break;
	}
	FT_ParallelExchGridArrayBuffer(array,front,NULL);
	switch (dim)
	{
        case 1:
            for (i = 0; i <= top_gmax[0]; ++i)
            {
		ic = d_index1d(i,top_gmax);
		if (top_comp[ic] == soln_comp)
		    soln[ic] = array[ic];
		else if (top_comp[ic] == obst_comp)
		    soln[ic] = var_obst;
	    }
	    break;
        case 2:
            for (i = 0; i <= top_gmax[0]; ++i)
            for (j = 0; j <= top_gmax[1]; ++j)
            {
		ic = d_index2d(i,j,top_gmax);
		if (top_comp[ic] == soln_comp)
		    soln[ic] = array[ic];
		else if (top_comp[ic] == obst_comp)
		    soln[ic] = var_obst;
	    }
	    break;
        case 3:
            for (i = 0; i <= top_gmax[0]; ++i)
            for (j = 0; j <= top_gmax[1]; ++j)
            for (k = 0; k <= top_gmax[2]; ++k)
            {
		ic = d_index3d(i,j,k,top_gmax);
		if (top_comp[ic] == soln_comp)
		    soln[ic] = array[ic];
		else if (top_comp[ic] == obst_comp)
		    soln[ic] = var_obst;
	    }
	    break;
	}

	stop_clock("solveEX");
	if (debugging("trace")) 
	    printf("Leaving solveEX()\n");
}	/* solveEX */

void PARABOLIC_SOLVER::solveCEX()
{
	int i,j,k,l,m,ic,icn,icoords[MAXD];
	int gmin[MAXD],ipn[MAXD];
	double crx_coords[MAXD];
	double solute,solute_nb[2],dgrad[MAXD],grad_plus[MAXD],grad_minus[MAXD];
	double coef;
	COMPONENT comp;
	HYPER_SURF *hs;
	POINTER state;
	int fr_crx_grid_seg;
	const GRID_DIRECTION dir[3][2] = 
		{{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
	double v[MAXD],v_plus[MAXD],v_minus[MAXD];
	double v_new,v_old,v0;			// cell volumes.
	double eta_new,eta_old; // cell irregularity ratio
	double kappa[3][2];	// side regularity number:
				// 1: for regular; 0: for irregular
	double sink;

	if (debugging("trace")) 
	    printf("Entering solveCEX()\n");
	start_clock("solveCEX");

	bdry_influx = compBdryFlux(var,dt);

	coef = D*dt;
	v0 = 1.0;
	for (i = 0; i < dim; ++i) 
	{
	    gmin[i] = 0;
	    v0 *= top_h[i];
	}

	switch (dim)
	{
        case 1:
            for (i = imin; i <= imax; ++i)
            {
                icoords[0] = i;
                ic = d_index1d(i,top_gmax);
                comp = top_comp[ic];
                if (comp != soln_comp)
                    continue;
                array[ic] = solute = var[ic];
                for (l = 0; l < dim; ++l)
                {
		    v[l] = 0.0;
		    v_plus[l] = 0.0;
		    v_minus[l] = 0.0;
		}
		if (a != NULL)
		{
		    for (l = 0; l < dim; ++l)
		    {
			v[l] = a[l][ic];
			v_plus[l] = std::max(0.0,v[l]);
			v_minus[l] = std::min(0.0,v[l]);
		    }

		}
		for (l = 0; l < dim; ++l)
		{
                    dgrad[l] = 0.0;
		    grad_plus[l] = 0.0;
		    grad_minus[l] = 0.0;

                    for (m = 0; m < 2; ++m)
                    {
			kappa[l][m] = 1.0;
			fr_crx_grid_seg = (*findStateAtCrossing)(front,icoords,
					dir[l][m],comp,&state,&hs,crx_coords);
                        if (fr_crx_grid_seg == NO_PDE_BOUNDARY)
                        {
                             next_ip_in_dir(icoords,dir[l][m],ipn,gmin,
							top_gmax);
                             icn = d_index1d(ipn[0],top_gmax);
			     solute_nb[m] = var[icn];
                        }
			else if (fr_crx_grid_seg == DIRICHLET_PDE_BOUNDARY)
			{
			    solute_nb[m] = (*getStateVarFunc)(state);
			}
			else if (fr_crx_grid_seg == MOVING_BOUNDARY)
			{
			    kappa[l][m] = 0.0;
			}
                        dgrad[l] += kappa[l][m]*(solute_nb[m] - solute)
						/top_h[l];
                    }
		    grad_plus[l] = (solute_nb[1] - solute)/top_h[l];
		    grad_minus[l] = (solute - solute_nb[0])/top_h[l];
                }
		v_new = cell_part[ic].vol_new[1];
		v_old = cell_part[ic].vol_old[1];
		eta_new = v_new/v0;
		eta_old = v_old/v0;
		sink = source[ic]/v0;
		array[ic] *= eta_old;
		for (l = 0; l < dim; ++l)
		    array[ic] += coef*dgrad[l]/top_h[l]-
			dt*(v_plus[l]*grad_minus[l]+v_minus[l]*grad_plus[l]);
		array[ic] += sink;
		array[ic] /= eta_new;
            }
            break;
	case 2:
	    for (i = imin; i <= imax; ++i)
	    for (j = jmin; j <= jmax; ++j)
	    {
		boolean irregular = NO;
	    	icoords[0] = i;
	    	icoords[1] = j;
	    	ic = d_index2d(i,j,top_gmax);
	    	comp = top_comp[ic];
	    	if (comp != soln_comp) 
	    	    continue;
                array[ic] = solute = var[ic];
		for (l = 0; l < dim; ++l)
		{
		    v[l] = 0.0;
		    v_plus[l] = 0.0;
		    v_minus[l] = 0.0;
		}
		if (a != NULL)
		{
		    for (l = 0; l < dim; ++l)
		    {
			v[l] = a[l][ic];
			v_plus[l] = std::max(0.0,v[l]);
			v_minus[l] = std::min(0.0,v[l]);
		    }
		}
		for (l = 0; l < dim; ++l)
		{
	            dgrad[l] = 0.0;
	 	    grad_plus[l] = 0.0;
		    grad_minus[l] = 0.0;

                    for (m = 0; m < 2; ++m)
                    {
			kappa[l][m] = 1.0;
			fr_crx_grid_seg = (*findStateAtCrossing)(front,icoords,
					dir[l][m],comp,&state,&hs,crx_coords);
                        if (fr_crx_grid_seg == NO_PDE_BOUNDARY)
                        {
                            next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax);
                            icn = d_index2d(ipn[0],ipn[1],top_gmax);
			    solute_nb[m] = var[icn];
                        }
			else if (fr_crx_grid_seg == DIRICHLET_PDE_BOUNDARY)
			{
			    solute_nb[m] = (*getStateVarFunc)(state);
			}
			else if (fr_crx_grid_seg == MOVING_BOUNDARY)
			{
			    irregular = YES;
			    kappa[l][m] = 0.0;
			}
                        dgrad[l] += kappa[l][m]*(solute_nb[m] - solute)
						/top_h[l];
                    }
		    grad_plus[l] = (solute_nb[1] - solute)/top_h[l];
		    grad_minus[l] = (solute - solute_nb[0])/top_h[l];
		}
	    	v_new = cell_part[ic].vol_new[1];
	    	v_old = cell_part[ic].vol_old[1];
		eta_new = v_new/v0;
                eta_old = v_old/v0;
                sink = source[ic]/v0;
                array[ic] *= eta_old;
		for (l = 0; l < dim; ++l)
		    array[ic] += coef*dgrad[l]/top_h[l] - 
			dt*(v_plus[l]*grad_minus[l]+v_minus[l]*grad_plus[l]);
	    	array[ic] += sink;
		array[ic] /= eta_new;
	    }
	    break;
	case 3:
	    for (i = imin; i <= imax; ++i)
	    for (j = jmin; j <= jmax; ++j)
	    for (k = kmin; k <= kmax; ++k)
	    {
	    	icoords[0] = i;
	    	icoords[1] = j;
	    	icoords[2] = k;
	    	ic = d_index3d(i,j,k,top_gmax);
	    	comp = top_comp[ic];
	    	if (comp != soln_comp) 
	    	    continue;
                array[ic] = solute = var[ic];
		for (l = 0; l < dim; ++l)
		{
		    v[l] = 0.0;
		    v_plus[l] = 0.0;
		    v_minus[l] = 0.0;
		}
		if (a != NULL)
		{
		    for ( l = 0; l < dim; ++l)
		    {	
			v[l] = a[l][ic];
			v_plus[l] = std::max(0.0,v[l]);
			v_minus[l] = std::min(0.0,v[l]);
		    }
		}
		for (l = 0; l < dim; ++l)
		{
	            dgrad[l] = 0.0;
	            grad_plus[l] = 0.0;
		    grad_minus[l] = 0.0;

                    for (m = 0; m < 2; ++m)
                    {
                        fr_crx_grid_seg = FT_StateVarAtGridCrossing(front,
                                        icoords,dir[l][m],comp,getStateVarFunc,
                                        &solute_nb[m],crx_coords);
                        if (!fr_crx_grid_seg)
                        {
                            next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax);
                            icn = d_index3d(ipn[0],ipn[1],ipn[2],top_gmax);
                	    solute_nb[m] = var[icn];
                        }
                        dgrad[l] += (solute_nb[m] - solute)/top_h[l];
                    }

		    grad_plus[l] = (solute_nb[1] - solute)/top_h[l];
		    grad_minus[l] = (solute - solute_nb[0])/top_h[l];

		    array[ic] += coef*dgrad[l]/top_h[l] - 
			dt*(v_plus[l]*grad_minus[l]+v_minus[l]*grad_plus[l]);
		}
	    }
	    break;
	}
	FT_ParallelExchGridArrayBuffer(array,front,NULL);
	switch (dim)
	{
        case 1:
            for (i = 0; i <= top_gmax[0]; ++i)
            {
		ic = d_index1d(i,top_gmax);
		if (top_comp[ic] == soln_comp)
		    soln[ic] = array[ic];
		else if (top_comp[ic] == obst_comp)
		    soln[ic] = var_obst;
	    }
	    break;
        case 2:
            for (i = 0; i <= top_gmax[0]; ++i)
            for (j = 0; j <= top_gmax[1]; ++j)
            {
		ic = d_index2d(i,j,top_gmax);
		if (top_comp[ic] == soln_comp)
		    soln[ic] = array[ic];
		else if (top_comp[ic] == obst_comp)
		    soln[ic] = var_obst;
	    }
	    break;
        case 3:
            for (i = 0; i <= top_gmax[0]; ++i)
            for (j = 0; j <= top_gmax[1]; ++j)
            for (k = 0; k <= top_gmax[2]; ++k)
            {
		ic = d_index3d(i,j,k,top_gmax);
		if (top_comp[ic] == soln_comp)
		    soln[ic] = array[ic];
		else if (top_comp[ic] == obst_comp)
		    soln[ic] = var_obst;
	    }
	    break;
	}

	stop_clock("solveCEX");
	if (debugging("trace")) 
	    printf("Leaving solveCEX()\n");
}	/* solveCEX */

void PARABOLIC_SOLVER::solveCN(void)
{
	int i,j,k,l,m,ic,icn,I,I_nb,icoords[MAXD];
	int gmin[MAXD],ipn[MAXD];
	double crx_coords[MAXD];
	double C0,C_nb,lambda,coeff,coeff_nb,rhs;
	COMPONENT comp;
	PETSc solver;
	double *x;
        int num_iter = 0;
	int size;
        double rel_residual = 0;
        boolean fr_crx_grid_seg;
        const GRID_DIRECTION dir[3][2] =
                {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};

	if (debugging("trace")) printf("Entering solveCN()\n");
	start_clock("solveCN");

	bdry_influx = compBdryFlux(var,0.5*dt);

	for (i = 0; i < dim; ++i) gmin[i] = 0;

	size = iupper - ilower;
        
	start_clock("set_coefficients");
	switch(dim)
	{
	case 1:
	    solver.Create(ilower, iupper-1, 3, 3);
	    for (i = imin; i <= imax; ++i)
	    {
	    	icoords[0] = i;
	    	ic = d_index1d(i,top_gmax);
	    	comp = top_comp[ic];
	    	if (comp != soln_comp) 
	    	    continue;
		I = i_to_I[i];
                C0 = var[ic];
		rhs = C0;
		coeff = 1.0;
		for (l = 0; l < dim; ++l)
		{
		    lambda = D*dt/sqr(top_h[l]);
                    for (m = 0; m < 2; ++m)
                    {
                        next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax);
                        icn = d_index1d(ipn[0],top_gmax);
			I_nb = i_to_I[ipn[0]];
			coeff_nb = -0.5*lambda;
                        fr_crx_grid_seg = FT_StateVarAtGridCrossing(front,
                                        icoords,dir[l][m],comp,getStateVarFunc,
                                        &C_nb,crx_coords);
                        if (!fr_crx_grid_seg)
                        {
			    solver.Add_A(I,I_nb,coeff_nb);
			    C_nb = var[icn];
			    rhs -= coeff_nb*(C_nb - C0);
		    	    coeff += 0.5*lambda;
                        }
			else
			{
			    rhs -= coeff_nb*(2.0*C_nb - C0);
		    	    coeff += 0.5*lambda;
			}
                    }
		}
		solver.Add_A(I,I,coeff);
		solver.Add_b(I,rhs);
	    }
	    break;
	case 2:
	    solver.Create(ilower, iupper-1, 5, 5);
	    for (i = imin; i <= imax; ++i)
	    for (j = jmin; j <= jmax; ++j)
	    {
	    	icoords[0] = i;
	    	icoords[1] = j;
	    	ic = d_index2d(i,j,top_gmax);
	    	comp = top_comp[ic];
		I = ij_to_I[i][j];
	    	if (comp != soln_comp) 
	    	    continue;
                C0 = var[ic];
		rhs = C0;
		coeff = 1.0;
		for (l = 0; l < dim; ++l)
		{
		    lambda = D*dt/sqr(top_h[l]);
                    for (m = 0; m < 2; ++m)
                    {
                        next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax);
                        icn = d_index2d(ipn[0],ipn[1],top_gmax);
			I_nb = ij_to_I[ipn[0]][ipn[1]];
			coeff_nb = -0.5*lambda;
                        fr_crx_grid_seg = FT_StateVarAtGridCrossing(front,
                                        icoords,dir[l][m],comp,getStateVarFunc,
                                        &C_nb,crx_coords);
                        if (!fr_crx_grid_seg)
                        {
			    solver.Add_A(I,I_nb,coeff_nb);
			    C_nb = var[icn];
			    rhs -= coeff_nb*(C_nb - C0);
			    coeff += 0.5*lambda;
                        }
			else
			{
			    rhs -= coeff_nb*(2.0*C_nb - C0);
			    coeff += 0.5*lambda;
			}
                    }
		}
		solver.Add_A(I,I,coeff);
		solver.Add_b(I,rhs);
	    }
	    break;
	case 3:
	    solver.Create(ilower, iupper-1, 7, 7);
	    for (i = imin; i <= imax; ++i)
	    for (j = jmin; j <= jmax; ++j)
	    for (k = kmin; k <= kmax; ++k)
	    {
	    	icoords[0] = i;
	    	icoords[1] = j;
	    	icoords[2] = k;
	    	ic = d_index3d(i,j,k,top_gmax);
	    	comp = top_comp[ic];
		I = ijk_to_I[i][j][k];
	    	if (comp != soln_comp) 
	    	    continue;
                C0 = var[ic];
		rhs = C0;
		coeff = 1.0;
		for (l = 0; l < dim; ++l)
		{
		    lambda = D*dt/sqr(top_h[l]);
                    for (m = 0; m < 2; ++m)
                    {
                        next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax);
                        icn = d_index3d(ipn[0],ipn[1],ipn[2],top_gmax);
			I_nb = ijk_to_I[ipn[0]][ipn[1]][ipn[2]];
			coeff_nb = -0.5*lambda;
                        fr_crx_grid_seg = FT_StateVarAtGridCrossing(front,
                                        icoords,dir[l][m],comp,getStateVarFunc,
                                        &C_nb,crx_coords);
                        if (!fr_crx_grid_seg)
                        {
			    solver.Add_A(I,I_nb,coeff_nb);
			    C_nb = var[icn];
			    rhs -= coeff_nb*(C_nb - C0);
			    coeff += 0.5*lambda;
                        }
			else
			{
			    rhs -= coeff_nb*(2.0*C_nb - C0);
			    coeff += 0.5*lambda;
			}
                    }
		}
		solver.Add_A(I,I,coeff);
		solver.Add_b(I,rhs);
	    }
	    break;
	}
	stop_clock("set_coefficients");
	start_clock("petsc_solve");
	solver.SetMaxIter(500);   
	solver.SetTol(1e-8);   
	solver.Solve();

	if (debugging("PETSc"))
	{
	    (void) printf("C_CARTESIAN::computeAdvectionCN: "
	       		"num_iter = %d, rel_residual = %g \n", 
			num_iter, rel_residual);
	}

	FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));
	solver.Get_x(x);
	stop_clock("petsc_solve");

	start_clock("scatter_data");
	switch (dim)
	{
        case 1:
            for (i = imin; i <= imax; i++)
            {
	    	I = i_to_I[i];
	    	ic = d_index1d(i,top_gmax);
	    	comp = top_comp[ic];
	    	if (comp == soln_comp)
	    	    array[ic] = x[I-ilower];
	    	else if (comp == obst_comp)
	    	    array[ic] = var_obst;
	    }
	    break;
        case 2:
            for (i = imin; i <= imax; i++)
            for (j = jmin; j <= jmax; j++)
            {
	    	I = ij_to_I[i][j];
	    	ic = d_index2d(i,j,top_gmax);
	    	comp = top_comp[ic];
	    	if (comp == soln_comp)
	    	    array[ic] = x[I-ilower];
	    	else if (comp == obst_comp)
	    	    array[ic] = var_obst;
	    }
	    break;
        case 3:
            for (i = imin; i <= imax; i++)
            for (j = jmin; j <= jmax; j++)
            for (k = kmin; k <= kmax; k++)
            {
	    	I = ijk_to_I[i][j][k];
	    	ic = d_index3d(i,j,k,top_gmax);
	    	comp = top_comp[ic];
	    	if (comp == soln_comp)
	    	    array[ic] = x[I-ilower];
	    	else if (comp == obst_comp)
	    	    array[ic] = var_obst;
	    }
	    break;
	}
	FT_ParallelExchGridArrayBuffer(array,front,NULL);
	switch (dim)
	{
        case 1:
            for (i = 0; i <= top_gmax[0]; ++i)
            {
		ic = d_index1d(i,top_gmax);
		soln[ic] = array[ic];
	    }
	    break;
        case 2:
            for (i = 0; i <= top_gmax[0]; ++i)
            for (j = 0; j <= top_gmax[1]; ++j)
            {
		ic = d_index2d(i,j,top_gmax);
		soln[ic] = array[ic];
	    }
	    break;
        case 3:
            for (i = 0; i <= top_gmax[0]; ++i)
            for (j = 0; j <= top_gmax[1]; ++j)
            for (k = 0; k <= top_gmax[2]; ++k)
            {
		ic = d_index3d(i,j,k,top_gmax);
		soln[ic] = array[ic];
	    }
	    break;
	}
	stop_clock("scatter_data");
	FT_FreeThese(1,x);

	bdry_influx += compBdryFlux(soln,0.5*dt);

	stop_clock("solveCN");
	if (debugging("trace")) printf("Leaving solveCN()\n");
}	/* end solveCN */

void PARABOLIC_SOLVER::solveCCN(void)
{
	int i,j,k,l,m,ic,icn,I,I_nb,icoords[MAXD];
	int gmin[MAXD],ipn[MAXD];
	double crx_coords[MAXD];
        double crx_coords_old[MAXD];
	double C0,C_nb,lambda,coeff,coeff_nb,rhs;
	COMPONENT comp;
	PETSc solver;
	double *x;
        int num_iter = 0;
	int size;
        double rel_residual = 0;
        int fr_crx_grid_seg;
        const GRID_DIRECTION dir[3][2] =
                {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
	HYPER_SURF *hs;
	POINTER state;
        double v_new,v_old,v0;                  // cell volumes.
        double eta_new,eta_old; // cell irregularity ratio
        double kappa[3][2];     // side regularity number:
                                // 1: for regular; 0: for irregular
	double sink;

	start_clock("solveCCN");

	bdry_influx = compBdryFlux(var,0.5*dt);

	v0 = 1.0;
	for (i = 0; i < dim; ++i)
        {
            gmin[i] = 0;
            v0 *= top_h[i];
        }

	size = iupper - ilower;
        
	start_clock("set_coefficients");
	switch(dim)
	{
	case 1:
	    solver.Create(ilower, iupper-1, 3, 3);
	    for (i = imin; i <= imax; ++i)
	    {
	    	icoords[0] = i;
	    	ic = d_index1d(i,top_gmax);
	    	comp = top_comp[ic];
	    	if (comp != soln_comp) 
	    	    continue;
		I = i_to_I[i];
                C0 = var[ic];
		rhs = 0.0;
		coeff = 0.0;
		for (l = 0; l < dim; ++l)
		{
		    lambda = D*dt/sqr(top_h[l]);
                    for (m = 0; m < 2; ++m)
                    {
                        next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax);
                        icn = d_index1d(ipn[0],top_gmax);
			I_nb = i_to_I[ipn[0]];
			coeff_nb = -0.5*lambda;
			fr_crx_grid_seg = (*findStateAtCrossing)(front,icoords,
                                        dir[l][m],comp,&state,&hs,crx_coords);
                        if (!fr_crx_grid_seg)
                        {
			    solver.Add_A(I,I_nb,coeff_nb);
			    C_nb = var[icn];
			    rhs -= coeff_nb*(C_nb - C0);
		    	    coeff += 0.5*lambda;
                        }
			else if (fr_crx_grid_seg == DIRICHLET_PDE_BOUNDARY)
			{
			    C_nb = (*getStateVarFunc)(state);
			    rhs -= coeff_nb*(2.0*C_nb - C0);
		    	    coeff += 0.5*lambda;
			}
                    }
		}
		v_new = cell_part[ic].vol_new[1];
		v_old = cell_part[ic].vol_old[1];
                eta_new = v_new/v0;
                eta_old = v_old/v0;
		rhs += eta_old*C0;
		coeff += eta_new;
                sink = source[ic]/v0;
		rhs += sink;
		solver.Add_A(I,I,coeff);
		solver.Add_b(I,rhs);
	    }
	    break;
	case 2:
	    solver.Create(ilower, iupper-1, 5, 5);
	    for (i = imin; i <= imax; ++i)
	    for (j = jmin; j <= jmax; ++j)
	    {
	    	icoords[0] = i;
	    	icoords[1] = j;
	    	ic = d_index2d(i,j,top_gmax);
	    	comp = top_comp[ic];
		I = ij_to_I[i][j];
	    	if (comp != soln_comp) 
	    	    continue;
                C0 = var[ic];
		rhs = 0.0;
		coeff = 0.0;
		for (l = 0; l < dim; ++l)
		{
		    lambda = D*dt/sqr(top_h[l]);
                    for (m = 0; m < 2; ++m)
                    {
                        next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax);
                        icn = d_index2d(ipn[0],ipn[1],top_gmax);
			I_nb = ij_to_I[ipn[0]][ipn[1]];
			coeff_nb = -0.5*lambda;
			fr_crx_grid_seg = (*findStateAtCrossing)(front,icoords,
                                        dir[l][m],comp,&state,&hs,crx_coords);
                        if (!fr_crx_grid_seg)
                        {
			    solver.Add_A(I,I_nb,coeff_nb);
			    C_nb = var[icn];
			    rhs -= coeff_nb*(C_nb - C0);
			    coeff += 0.5*lambda;
                        }
			else if (fr_crx_grid_seg == DIRICHLET_PDE_BOUNDARY)
			{
			    C_nb = (*getStateVarFunc)(state);
			    rhs -= coeff_nb*(2.0*C_nb - C0);
			    coeff += 0.5*lambda;
			}
                    }
		}
		v_new = cell_part[ic].vol_new[1];
		v_old = cell_part[ic].vol_old[1];
                eta_new = v_new/v0;
                eta_old = v_old/v0;
		rhs += eta_old*C0;
		coeff += eta_new;
		sink = source[ic]/v0;
		rhs += sink;
		solver.Add_A(I,I,coeff);
		solver.Add_b(I,rhs);
	    }
	    break;
	case 3:
	    solver.Create(ilower, iupper-1, 7, 7);
	    for (i = imin; i <= imax; ++i)
	    for (j = jmin; j <= jmax; ++j)
	    for (k = kmin; k <= kmax; ++k)
	    {
	    	icoords[0] = i;
	    	icoords[1] = j;
	    	icoords[2] = k;
	    	ic = d_index3d(i,j,k,top_gmax);
	    	comp = top_comp[ic];
		I = ijk_to_I[i][j][k];
	    	if (comp != soln_comp) 
	    	    continue;
                C0 = var[ic];
		rhs = C0;
		coeff = 1.0;
		for (l = 0; l < dim; ++l)
		{
		    lambda = D*dt/sqr(top_h[l]);
                    for (m = 0; m < 2; ++m)
                    {
                        next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax);
                        icn = d_index3d(ipn[0],ipn[1],ipn[2],top_gmax);
			I_nb = ijk_to_I[ipn[0]][ipn[1]][ipn[2]];
			coeff_nb = -0.5*lambda;
                        fr_crx_grid_seg = FT_StateVarAtGridCrossing(front,
                                        icoords,dir[l][m],comp,getStateVarFunc,
                                        &C_nb,crx_coords);
                        if (!fr_crx_grid_seg)
                        {
			    solver.Add_A(I,I_nb,coeff_nb);
			    C_nb = var[icn];
			    rhs -= coeff_nb*(C_nb - C0);
			    coeff += 0.5*lambda;
                        }
			else
			{
			    rhs -= coeff_nb*(2.0*C_nb - C0);
			    coeff += 0.5*lambda;
			}
                    }
		}
		solver.Add_A(I,I,coeff);
		solver.Add_b(I,rhs);
	    }
	    break;
	}
	stop_clock("set_coefficients");
	start_clock("petsc_solve");
	solver.SetMaxIter(500);   
	solver.SetTol(1e-8);   
	solver.Solve();

	if (debugging("PETSc"))
	{
	    (void) printf("C_CARTESIAN::computeAdvectionCN: "
	       		"num_iter = %d, rel_residual = %g \n", 
			num_iter, rel_residual);
	}

	FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));
	solver.Get_x(x);
	stop_clock("petsc_solve");

	start_clock("scatter_data");
	switch (dim)
	{
        case 1:
            for (i = imin; i <= imax; i++)
            {
	    	I = i_to_I[i];
	    	ic = d_index1d(i,top_gmax);
	    	comp = top_comp[ic];
	    	if (comp == soln_comp)
	    	    array[ic] = x[I-ilower];
	    	else if (comp == obst_comp)
	    	    array[ic] = var_obst;
	    }
	    break;
        case 2:
            for (i = imin; i <= imax; i++)
            for (j = jmin; j <= jmax; j++)
            {
	    	I = ij_to_I[i][j];
	    	ic = d_index2d(i,j,top_gmax);
	    	comp = top_comp[ic];
	    	if (comp == soln_comp)
	    	    array[ic] = x[I-ilower];
	    	else if (comp == obst_comp)
	    	    array[ic] = var_obst;
	    }
	    break;
        case 3:
            for (i = imin; i <= imax; i++)
            for (j = jmin; j <= jmax; j++)
            for (k = kmin; k <= kmax; k++)
            {
	    	I = ijk_to_I[i][j][k];
	    	ic = d_index3d(i,j,k,top_gmax);
	    	comp = top_comp[ic];
	    	if (comp == soln_comp)
	    	    array[ic] = x[I-ilower];
	    	else if (comp == obst_comp)
	    	    array[ic] = var_obst;
	    }
	    break;
	}
	FT_ParallelExchGridArrayBuffer(array,front,NULL);
	switch (dim)
	{
        case 1:
            for (i = 0; i <= top_gmax[0]; ++i)
            {
		ic = d_index1d(i,top_gmax);
		soln[ic] = array[ic];
	    }
	    break;
        case 2:
            for (i = 0; i <= top_gmax[0]; ++i)
            for (j = 0; j <= top_gmax[1]; ++j)
            {
		ic = d_index2d(i,j,top_gmax);
		soln[ic] = array[ic];
	    }
	    break;
        case 3:
            for (i = 0; i <= top_gmax[0]; ++i)
            for (j = 0; j <= top_gmax[1]; ++j)
            for (k = 0; k <= top_gmax[2]; ++k)
            {
		ic = d_index3d(i,j,k,top_gmax);
		soln[ic] = array[ic];
	    }
	    break;
	}
	stop_clock("scatter_data");
	FT_FreeThese(1,x);

	bdry_influx += compBdryFlux(soln,0.5*dt);

	stop_clock("solveCCN");
}	/* end solveCCN */

void PARABOLIC_SOLVER::solveCIM()
{
	static double *var_mid;

	bdry_influx = 0.0;
	if (var_mid == NULL)
	{
	    if (order > 1)
	    	FT_VectorMemoryAlloc((POINTER*)&var_mid,array_size,
				sizeof(double));
	}
	switch (order)
	{
	case 1:
	    sub_dt = dt;
	    Csolve(var,soln);
	    bdry_influx = compBdryFlux(soln,sub_dt);
	    break;
	case 2:
	    sub_dt = dt;
	    Csolve(var,var_mid);
	    bdry_influx = compBdryFlux(var_mid,0.5*dt);
	    addMeshState(var_mid,0.5,var,0.5,var_mid);
	    sub_dt = 0.5*dt;
	    Csolve(var_mid,soln);
	    bdry_influx += compBdryFlux(soln,0.5*dt);
	    break;
	case 4:
	default:
	    (void)printf("ERROR: %d-th order RK method not implemented\n",
					order);
	    clean_up(ERROR);
	}

}	/* end solveCIM */

void PARABOLIC_SOLVER::Csolve(
	double *var_in,
	double *var_out)
{
	switch (dim)
	{
	case 1:
	    return Csolve1d(var_in,var_out);
	case 2:
	    return Csolve2d(var_in,var_out);
	case 3:
	    return Csolve3d(var_in,var_out);
	}
}	/* end Csolve */


void PARABOLIC_SOLVER::Csolve1d(
	double *var_in,
	double *var_out)
{
	int i,l,m,ic,icn,I,I_nb,icoords[MAXD];
	int gmin[MAXD],ipn[MAXD];
	double crx_coords[MAXD];
	double crx_coords_old[MAXD];
	double C_nb,lambda,coeff,coeff_nb,rhs;
	COMPONENT comp;
	PETSc solver;
	double *x;
        int num_iter = 0;
	int size;
        double rel_residual = 0;
        int fr_crx_grid_seg;
	const GRID_DIRECTION dir[3][2] =
                {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
	HYPER_SURF *hs;
        POINTER state;
        double v_new,v_old,v0;                  // cell volumes.
        double eta_new,eta_old; // cell irregularity ratio
        double kappa[3][2];     // side regularity number:
                                // 1: for regular; 0: for irregular
	double sink;

        start_clock("Csolve1d");

	v0 = 1.0;
	for (i = 0; i < dim; ++i)
	{
	    gmin[i] = 0;
	    v0 *= top_h[i];
	}

        size = iupper - ilower;

	solver.Create(ilower, iupper-1, 3, 3);
	for (i = imin; i <= imax; ++i)
	{
	    icoords[0] = i;
	    ic = d_index(icoords,top_gmax,dim);
	    comp = top_comp[ic];
	    I = i_to_I[i];
	    if (comp != soln_comp) 
	    	continue;
	    rhs = 0.0;
	    coeff = 0.0;
	    for (l = 0; l < dim; ++l)
	    {
                for (m = 0; m < 2; ++m)
                {
		    if (nu == NULL)
                        lambda = D*sub_dt/sqr(top_h[l]);
                    else
                        lambda = nu[ic]*sub_dt/sqr(top_h[l]);
                    next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax);
                    icn = d_index(ipn,top_gmax,dim);
		    I_nb = i_to_I[ipn[0]];
                    fr_crx_grid_seg = (*findStateAtCrossing)(front,icoords,
				dir[l][m],comp,&state,&hs,crx_coords); 

                    if (fr_crx_grid_seg == NO_PDE_BOUNDARY)
                    {
			if (nu != NULL)
                        {
                            lambda += nu[icn]*sub_dt/sqr(top_h[l]);
                            lambda *= 0.5;
                        }
                        coeff_nb = -lambda;
                        coeff += lambda;
			solver.Add_A(I,I_nb,coeff_nb);
                    }
		    else if (fr_crx_grid_seg == DIRICHLET_PDE_BOUNDARY)
		    {
			coeff_nb = -lambda;
			if (wave_type(hs) == DIRICHLET_BOUNDARY &&
                            boundary_state_function(hs) &&
                            strcmp(boundary_state_function_name(hs),
                            "flowThroughBoundaryState") == 0)
                            C_nb = var_in[ic];
                        else
			{
                            coeff += lambda;
                            C_nb = (*getStateVarFunc)(state);
                            rhs -= coeff_nb*C_nb;
			}
		    }
		    else if (fr_crx_grid_seg != NEUMANN_PDE_BOUNDARY &&
		    	     fr_crx_grid_seg != MOVING_BOUNDARY)
                        coeff += lambda;

                }
	    }
	    v_new = cell_part[ic].vol_new[1];
	    v_old = cell_part[ic].vol_old[1];
	    eta_new = v_new/v0;
	    eta_old = v_old/v0;
	    rhs += eta_old*var_in[ic];
	    coeff += eta_new;
	    sink = source[ic]/v0;
	    rhs += sink;
	    solver.Add_A(I,I,coeff);
	    solver.Add_b(I,rhs);
	}
	start_clock("petsc_solve");
	solver.SetMaxIter(500);   
	solver.SetTol(1e-10);   
	solver.Solve();

	if (debugging("PETSc"))
	{
	    (void) printf("C_CARTESIAN::computeAdvectionImplicit: "
	       		"num_iter = %d, rel_residual = %g \n", 
			num_iter, rel_residual);
	}

	FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));
	solver.Get_x(x);
	stop_clock("petsc_solve");

        for (i = imin; i <= imax; i++)
        {
	    icoords[0] = i;
	    I = i_to_I[i];
	    ic = d_index(icoords,top_gmax,dim);
	    comp = top_comp[ic];
	    if (comp == soln_comp)
	    	array[ic] = x[I-ilower];
	    else if (comp == obst_comp)
                array[ic] = var_obst;
	}
	FT_ParallelExchGridArrayBuffer(array,front,NULL);
        for (i = 0; i <= top_gmax[0]; ++i)
        {
	    icoords[0] = i;
	    ic = d_index(icoords,top_gmax,dim);
            var_out[ic] = array[ic];
	}
	FT_FreeThese(1,x);
        stop_clock("Csolve1d");
}	/* end Csolve1d */

void PARABOLIC_SOLVER::Csolve3d(
	double *var_in,
	double *var_out)
{
}

void PARABOLIC_SOLVER::Csolve2d(
	double *var_in,
	double *var_out)
{
	int i,j,k,l,m,ic,icn,I,I_nb,icoords[MAXD];
	int gmin[MAXD],ipn[MAXD];
	double crx_coords[MAXD];
	double C_nb,lambda,eta,coeff,coeff_nb,rhs;
	COMPONENT comp;
	PETSc solver;
	double *x;
        int num_iter = 0;
	int size;
        double rel_residual = 0;
        int fr_crx_grid_seg;
	const GRID_DIRECTION dir[3][2] =
                {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
	HYPER_SURF *hs;
	POINTER state;
        double v_new,v_old,v0;                  // cell volumes.
        double eta_new,eta_old; // cell irregularity ratio
        double kappa[3][2];     // side regularity number:
                                // 1: for regular; 0: for irregular
	double sink;

        start_clock("Csolve2d");

	v0 = 1.0;
        for (i = 0; i < dim; ++i)
        {
            gmin[i] = 0;
            v0 *= top_h[i];
        }

        size = iupper - ilower;

	solver.Create(ilower, iupper-1, 5, 5);
	for (i = imin; i <= imax; ++i)
	for (j = jmin; j <= jmax; ++j)
	{
	    icoords[0] = i;
	    icoords[1] = j;
	    ic = d_index(icoords,top_gmax,dim);
	    comp = top_comp[ic];
	    I = ij_to_I[i][j];
	    if (comp != soln_comp) 
	    	continue;
	    rhs = 0.0;
	    coeff = 0.0;
	    for (l = 0; l < dim; ++l)
	    {
                for (m = 0; m < 2; ++m)
                {
		    if (nu == NULL)
                        lambda = D*sub_dt/sqr(top_h[l]);
                    else
                        lambda = nu[ic]*sub_dt/sqr(top_h[l]);
                    if (next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax))
		    {
                    	icn = d_index(ipn,top_gmax,dim);
		    	I_nb = ij_to_I[ipn[0]][ipn[1]];
			eta = (a == NULL) ? 0.0 :
				0.25*(a[l][ic] + a[l][icn])*sub_dt/top_h[l];
		    }
		    else
		    {
		    	eta = 0.0;
			I_nb = -1;
		    }
		    fr_crx_grid_seg =(*findStateAtCrossing)(front,icoords,
				dir[l][m],comp,&state,&hs,crx_coords);

                    if (fr_crx_grid_seg == NO_PDE_BOUNDARY)
                    {
			if (nu != NULL)
                        {
                            lambda += nu[icn]*sub_dt/sqr(top_h[l]);
                            lambda *= 0.5;
                        }
			coeff_nb = -lambda;
                        coeff += lambda;
		        if (m == 0) coeff_nb -= eta;
                        else coeff_nb += eta;
			solver.Add_A(I,I_nb,coeff_nb);
                    }
		    else if (fr_crx_grid_seg == DIRICHLET_PDE_BOUNDARY)
		    {
			coeff_nb = -lambda;
			if (wave_type(hs) == DIRICHLET_BOUNDARY &&
			    boundary_state_function(hs) &&
                            strcmp(boundary_state_function_name(hs),
                            "flowThroughBoundaryState") == 0)
			{
			    C_nb = var_in[ic];
			    if (m == 0) coeff -= eta;
                            else coeff += eta;
			}
			else
			{
                            coeff += lambda;
			    C_nb = (*getStateVarFunc)(state);
			    rhs -= coeff_nb*C_nb;
                            if (m == 0) rhs += eta*C_nb;
                            else rhs -= eta*C_nb;
			}
		    }
		    else if (fr_crx_grid_seg != NEUMANN_PDE_BOUNDARY &&
		    	     fr_crx_grid_seg != MOVING_BOUNDARY)
			coeff += lambda;
	
                }
	    }
	    v_new = cell_part[ic].vol_new[1];
	    v_old = cell_part[ic].vol_old[1];
            eta_new = v_new/v0;
            eta_old = v_old/v0;
	    if (eta_new < 0.49)
	    {
	    	printf("Less than half size cell (%d, %d)\n",i,j);
		printf("eta_old = %20.14f  eta_new = %20.14f\n",
				eta_old,eta_new);
		clean_up(ERROR);
	    }
            rhs += eta_old*var_in[ic];
            coeff += eta_new;
            sink = source[ic]/v0;
            rhs += sink;
	    solver.Add_A(I,I,coeff);
	    solver.Add_b(I,rhs);
	}
	start_clock("petsc_solve");
	solver.SetMaxIter(500);   
	solver.SetTol(1e-10);   
	solver.Solve();

	if (debugging("PETSc"))
	{
	    (void) printf("C_CARTESIAN::computeAdvectionImplicit: "
	       		"num_iter = %d, rel_residual = %g \n", 
			num_iter, rel_residual);
	}

	FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));
	solver.Get_x(x);
	stop_clock("petsc_solve");

        for (i = imin; i <= imax; i++)
        for (j = jmin; j <= jmax; j++)
        {
	    icoords[0] = i;
	    icoords[1] = j;
	    I = ij_to_I[i][j];
	    ic = d_index(icoords,top_gmax,dim);
	    comp = top_comp[ic];
	    if (comp == soln_comp)
	    	array[ic] = x[I-ilower];
	    else if (comp == obst_comp)
                array[ic] = var_obst;
	}
	FT_ParallelExchGridArrayBuffer(array,front,NULL);
        for (i = 0; i <= top_gmax[0]; ++i)
        for (j = 0; j <= top_gmax[1]; ++j)
        {
	    ic = d_index2d(i,j,top_gmax);
            var_out[ic] = array[ic];
	}
	FT_FreeThese(1,x);
        stop_clock("Csolve2d");
}	/* end Csolve2d */

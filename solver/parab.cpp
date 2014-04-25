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
	double v[MAXD];
	HYPER_SURF *hs;
	POINTER state;

        start_clock("solve2d");

        for (i = 0; i < dim; ++i) gmin[i] = 0;

        size = iupper - ilower;

	var_min = HUGE;
	var_max = -HUGE;
	for (i = imin; i <= imax; ++i)
	for (j = jmin; j <= jmax; ++j)
	{
	    ic = d_index2d(i,j,top_gmax);
	    comp = top_comp[ic];
	    if (comp != soln_comp) continue;
	    if (var_in[ic] < var_min) var_min = var_in[ic];
	    if (var_in[ic] > var_max) var_max = var_in[ic];
	}
	if (var_max - var_min < 0.1) 
	    var_min -= 0.01;
	else
	    var_min -= 0.1*(var_max - var_min);

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
	    rhs = var_in[ic] - var_min;
	    coeff = 1.0;
	    for (l = 0; l < dim; ++l) v[l] = 0.0;
	    if (a != NULL)
	    {
		for (l = 0; l < dim; ++l)
		    v[l] = a[l][ic];
	    }
	    for (l = 0; l < dim; ++l)
	    {
		eta = 0.5*v[l]*sub_dt/top_h[l];
                for (m = 0; m < 2; ++m)
                {
		    if (nu == NULL)
                        lambda = D*sub_dt/sqr(top_h[l]);
                    else
                        lambda = nu[ic]*sub_dt/sqr(top_h[l]);
                    next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax);
                    icn = d_index(ipn,top_gmax,dim);
		    I_nb = ij_to_I[ipn[0]][ipn[1]];
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
			    C_nb = var_in[ic] - var_min;
			    if (m == 0) coeff -= eta;
                            else coeff += eta;
			}
			else
			{
                            coeff += lambda;
			    C_nb = (*getStateVarFunc)(state) - var_min;
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
	    	array[ic] = x[I-ilower] + var_min;
	    else if (comp == obst_comp)
                array[ic] = var_obst;
	}
	FT_ParallelExchGridArrayBuffer(array,front,NULL);
        for (i = 0; i <= top_gmax[0]; ++i)
        for (j = 0; j <= top_gmax[1]; ++j)
        {
	    icoords[0] = i;
	    icoords[1] = j;
	    ic = d_index(icoords,top_gmax,dim);
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

	var_min = HUGE;
	var_max = -HUGE;
	for (i = imin; i <= imax; ++i)
	{
	    ic = d_index1d(i,top_gmax);
	    comp = top_comp[ic];
	    if (comp != soln_comp) continue;
	    if (var_in[ic] < var_min) var_min = var_in[ic];
	    if (var_in[ic] > var_max) var_max = var_in[ic];
	}
	if (var_max - var_min < 0.1)
            var_min -= 0.01;
        else
            var_min -= 0.1*(var_max - var_min);

	solver.Create(ilower, iupper-1, 3, 3);
	for (i = imin; i <= imax; ++i)
	{
	    icoords[0] = i;
	    ic = d_index(icoords,top_gmax,dim);
	    comp = top_comp[ic];
	    I = i_to_I[i];
	    if (comp != soln_comp) 
	    	continue;
	    rhs = var_in[ic] - var_min;
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
		    else if (fr_crx_grid_seg == DIRICHLET_PDE_BOUNDARY)
		    {
			coeff_nb = -lambda;
			if (wave_type(hs) == DIRICHLET_BOUNDARY &&
                            boundary_state_function(hs) &&
                            strcmp(boundary_state_function_name(hs),
                            "flowThroughBoundaryState") == 0)
                            C_nb = var_in[ic] - var_min;
                        else
			{
                            coeff += lambda;
                            C_nb = (*getStateVarFunc)(state) - var_min;
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
	    	array[ic] = x[I-ilower] + var_min;
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
	double v[MAXD];
	HYPER_SURF *hs;
	POINTER state;

        start_clock("solve3d");

        for (i = 0; i < dim; ++i) gmin[i] = 0;

        size = iupper - ilower;

	var_min = HUGE;
	var_max = -HUGE;
	for (i = imin; i <= imax; ++i)
	for (j = jmin; j <= jmax; ++j)
	for (k = kmin; k <= kmax; ++k)
	{
	    ic = d_index3d(i,j,k,top_gmax);
	    comp = top_comp[ic];
	    if (comp != soln_comp) continue;
	    if (var_in[ic] < var_min) var_min = var_in[ic];
	    if (var_in[ic] > var_max) var_max = var_in[ic];
	}
	if (var_max - var_min < 0.1)
            var_min -= 0.01;
        else
            var_min -= 0.1*(var_max - var_min);

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
	    rhs = var_in[ic] - var_min;
	    coeff = 1.0;
	    for (l = 0; l < dim; ++l) v[l] = 0.0;
	    if (a != NULL)
	    {
		for (l = 0; l < dim; ++l)
		    v[l] = a[l][ic];
	    }
	    for (l = 0; l < dim; ++l)
	    {
		eta = 0.5*v[l]*sub_dt/top_h[l];
                for (m = 0; m < 2; ++m)
                {
		    if (nu == NULL)
                        lambda = D*sub_dt/sqr(top_h[l]);
                    else
                        lambda = nu[ic]*sub_dt/sqr(top_h[l]);
                    next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax);
                    icn = d_index(ipn,top_gmax,dim);
		    I_nb = ijk_to_I[ipn[0]][ipn[1]][ipn[2]];
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
		    else if (fr_crx_grid_seg == DIRICHLET_PDE_BOUNDARY)
		    {
			coeff_nb = -lambda;
			if (wave_type(hs) == DIRICHLET_BOUNDARY &&
                            boundary_state_function(hs) &&
                            strcmp(boundary_state_function_name(hs),
                            "flowThroughBoundaryState") == 0)
			{
                            C_nb = var_in[ic] - var_min;
			    if (m == 0) coeff -= eta;
                            else coeff += eta;
			}
                        else
			{
			    coeff += lambda;
                            C_nb = (*getStateVarFunc)(state) - var_min;
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
	    	array[ic] = x[I-ilower] + var_min;
	    else if (comp == obst_comp)
                array[ic] = var_obst;
	}
	FT_ParallelExchGridArrayBuffer(array,front,NULL);
        for (i = 0; i <= top_gmax[0]; ++i)
        for (j = 0; j <= top_gmax[1]; ++j)
        for (k = 0; k <= top_gmax[2]; ++k)
        {
	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;
	    ic = d_index(icoords,top_gmax,dim);
	    var_out[ic] = array[ic];
	}
	FT_FreeThese(1,x);
        stop_clock("solve3d");
}	/* end solve3d */

void PARABOLIC_SOLVER::runge_kutta()
{
	static double *var_mid;

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
	    break;
	case 2:
	    sub_dt = dt;
	    solve(var,var_mid);
	    addMeshState(var_mid,0.5,var,0.5,var_mid);
	    sub_dt = 0.5*dt;
	    solve(var_mid,soln);
	    break;
	case 4:
	default:
	    (void)printf("ERROR: %d-th order RK method not implemented\n",
					order);
	    clean_up(ERROR);
	    }
}	/* end runge_kutta */

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
}	/* end runge_kutta */

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
		    w[1] = (*getStateVarFunc)(intfc_state) - var_min;
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
	rhs = array[id0] - (var_in[id0] - var_min);
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

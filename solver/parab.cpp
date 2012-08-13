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
	    {
		if (comp == obst_comp)
		    array[ic] = var_obst;
	    	continue;
	    }
	    rhs = var_in[ic];
	    coeff = 1.0;
	    for (l = 0; l < dim; ++l) v[l] = 0.0;
	    if (a != NULL)
	    {
		for (l = 0; l < dim; ++l)
		    v[l] = a[l][ic];
	    }
	    for (l = 0; l < dim; ++l)
	    {
		lambda = D*sub_dt/sqr(top_h[l]);
		eta = 0.5*v[l]*sub_dt/top_h[l];
		coeff += 2.0*lambda;
                for (m = 0; m < 2; ++m)
                {
                    next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax);
                    icn = d_index(ipn,top_gmax,dim);
		    I_nb = ij_to_I[ipn[0]][ipn[1]];
		    coeff_nb = -lambda;
		    fr_crx_grid_seg =(*findStateAtCrossing)(front,icoords,
				dir[l][m],comp,&state,&hs,crx_coords);
		    if (m == 0) coeff_nb -= eta;
		    else coeff_nb += eta;
                    if (fr_crx_grid_seg == NO_PDE_BOUNDARY)
                    {
			solver.Add_A(I,I_nb,coeff_nb);
                    }
		    else if (fr_crx_grid_seg == DIRICHLET_PDE_BOUNDARY)
		    {
			if (wave_type(hs) == DIRICHLET_BOUNDARY &&
			    boundary_state_function(hs) &&
                            strcmp(boundary_state_function_name(hs),
                            "flowThroughBoundaryState") == 0)
			    C_nb = var_in[ic];
			else
			    C_nb = (*getStateVarFunc)(state);
			rhs -= coeff_nb*C_nb;
			if (m == 0) rhs += eta*C_nb;
			else rhs -= eta*C_nb;
		    }
		    else if (fr_crx_grid_seg == NEUMANN_PDE_BOUNDARY)
			coeff -= lambda;
	
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
	FT_ParallelExchGridArrayBuffer(array,front);
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
	static boolean first = YES;
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
	    first = NO;
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
	double C_nb,lambda,eta,coeff,coeff_nb,rhs;
	COMPONENT comp;
	PETSc solver;
	double *x;
        int num_iter = 0;
	int size;
        double rel_residual = 0;
        boolean fr_crx_grid_seg;
	const GRID_DIRECTION dir[3][2] =
                {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
	double v[MAXD];

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
	    {
	    	if (comp == obst_comp) 
		    array[ic] = var_obst;
	    	continue;
	    }
	    rhs = var_in[ic];
	    coeff = 1.0;
	    for (l = 0; l < dim; ++l) v[l] = 0.0;
	    if (a != NULL)
	    {
		for (l = 0; l < dim; ++l)
		    v[l] = a[l][ic];
	    }
	    for (l = 0; l < dim; ++l)
	    {
		lambda = D*sub_dt/sqr(top_h[l]);
		eta = 0.5*v[l]*sub_dt/top_h[l];
		coeff += 2.0*lambda;
                for (m = 0; m < 2; ++m)
                {
                    next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax);
                    icn = d_index(ipn,top_gmax,dim);
		    I_nb = i_to_I[ipn[0]];
		    coeff_nb = -lambda;
                    fr_crx_grid_seg = FT_StateVarAtGridCrossing(front,
                                    icoords,dir[l][m],comp,getStateVarFunc,
                                    &C_nb,crx_coords);
		    if (m == 0) coeff_nb -= eta;
		    else coeff_nb += eta;
                    if (!fr_crx_grid_seg)
                    {
			solver.Add_A(I,I_nb,coeff_nb);
                    }
		    else
		    {
			rhs -= coeff_nb*C_nb;
			if (m == 0) rhs += eta*C_nb;
			else rhs -= eta*C_nb;
		    }
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
	FT_ParallelExchGridArrayBuffer(array,front);
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
        boolean fr_crx_grid_seg;
	const GRID_DIRECTION dir[3][2] =
                {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
	double v[MAXD];

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
	    {
	    	if (comp == obst_comp) 
		    array[ic] = var_obst;
	    	continue;
	    }
	    rhs = var_in[ic];
	    coeff = 1.0;
	    for (l = 0; l < dim; ++l) v[l] = 0.0;
	    if (a != NULL)
	    {
		for (l = 0; l < dim; ++l)
		    v[l] = a[l][ic];
	    }
	    for (l = 0; l < dim; ++l)
	    {
		lambda = D*sub_dt/sqr(top_h[l]);
		eta = 0.5*v[l]*sub_dt/top_h[l];
		coeff += 2.0*lambda;
                for (m = 0; m < 2; ++m)
                {
                    next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax);
                    icn = d_index(ipn,top_gmax,dim);
		    I_nb = ijk_to_I[ipn[0]][ipn[1]][ipn[2]];
		    coeff_nb = -lambda;
                    fr_crx_grid_seg = FT_StateVarAtGridCrossing(front,
                                    icoords,dir[l][m],comp,getStateVarFunc,
                                    &C_nb,crx_coords);
		    if (m == 0) coeff_nb -= eta;
		    else coeff_nb += eta;
                    if (!fr_crx_grid_seg)
                    {
			solver.Add_A(I,I_nb,coeff_nb);
                    }
		    else
		    {
			rhs -= coeff_nb*C_nb;
			if (m == 0) rhs += eta*C_nb;
			else rhs -= eta*C_nb;
		    }
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
	FT_ParallelExchGridArrayBuffer(array,front);
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
	    addMeshState(var,0.5,var,0.5,var_mid);
	    sub_dt = 0.5*dt;
	    solve(var,soln);
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

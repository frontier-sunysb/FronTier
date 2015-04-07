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

ELLIPTIC_SOLVER::ELLIPTIC_SOLVER(Front &front):front(&front)
{
	porosity = 0.0;
}

void ELLIPTIC_SOLVER::set_solver_domain(void)
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
	top_L = rgr->L;
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

void ELLIPTIC_SOLVER::solve(double *soln)
{
	switch (dim)
	{
	case 1:
	    return solve1d(soln);
	case 2:
	    return solve2d(soln);
	case 3:
	    return solve3d(soln);
	}
}	/* end solve */

void ELLIPTIC_SOLVER::dsolve(double *soln)
{
	switch (dim)
	{
	case 1:
	    (void) printf("dsolve1d() not implemented!\n");
	    clean_up(ERROR);
	case 2:
	    return dsolve2d(soln);
	case 3:
	    return dsolve3d(soln);
	}
}	/* end dsolve */

void ELLIPTIC_SOLVER::solve1d(double *soln)
{
	int index,index_nb[2],size;
	double k0,k_nb[2];
	double rhs,coeff[2];
	int I,I_nb[2];
	int i,l,icoords[MAXD];
	COMPONENT comp;
	double aII;
	int num_nb;
	GRID_DIRECTION dir[2] = {WEST,EAST};
	boolean use_neumann_solver = YES;
	PetscInt num_iter = 0;
	double rel_residual = 0.0;
	HYPER_SURF *hs;
	double crx_coords[MAXD];
        int status;
	POINTER intfc_state;

	PETSc solver;
	solver.Create(ilower, iupper-1, 3, 3);
	solver.Reset_A();
	solver.Reset_b();
	solver.Reset_x();
	size = iupper - ilower;
	max_soln = -HUGE;
	min_soln = HUGE;

        for (i = imin; i <= imax; i++)
	{
	    index  = d_index1d(i,top_gmax);
	    comp = top_comp[index];
	    I = i_to_I[i];
	    if (I == -1) continue;

	    index_nb[0] = d_index1d(i-1,top_gmax);
	    index_nb[1] = d_index1d(i+1,top_gmax);
	    I_nb[0] = i_to_I[i-1];
	    I_nb[1] = i_to_I[i+1];
	    icoords[0] = i;
	
	    k0 = D[index];
	    num_nb = 0;
	    for (l = 0; l < 2; ++l)
	    {
		status = (*findStateAtCrossing)(front,icoords,dir[l],comp,
                                &intfc_state,&hs,crx_coords);
		if (status != CONST_V_PDE_BOUNDARY)
		    num_nb++;
                if (status == CONST_V_PDE_BOUNDARY ||
		    status == CONST_P_PDE_BOUNDARY)
		    index_nb[l] = index;
		k_nb[l] = 0.5*(k0 + D[index_nb[l]]);
	    	coeff[l] = k_nb[l]/(top_h[l/2]*top_h[l/2]); 
	    }

	    rhs = source[index];

	    aII = 0.0;
	    for (l = 0; l < 2; ++l)
	    {
		if (num_nb == 0) break;
		status = (*findStateAtCrossing)(front,icoords,dir[l],comp,
                                &intfc_state,&hs,crx_coords);
                if (status == NO_PDE_BOUNDARY)
                {
                    solver.Set_A(I,I_nb[l],coeff[l]);
                    aII += -coeff[l];
                }
                else if (status == CONST_P_PDE_BOUNDARY)
                {
                    aII += -coeff[l];
		    rhs += -coeff[l]*getStateVar(intfc_state);
		    use_neumann_solver = NO;
                }
	    }
	    /*
	     * This change reflects the need to treat point with only one
	     * interior neighbor (a convex point). Not sure why PETSc cannot
	     * handle such case. If we have better understanding, this should
	     * be changed back.
	     */
	    if(num_nb > 0)
	    {
                solver.Set_A(I,I,aII);
	    }
            else
            {
	
		(void) printf("WARNING: isolated value!\n");
                solver.Set_A(I,I,1.0);
		rhs = soln[index];
            }
            solver.Set_b(I,rhs);
	}
	use_neumann_solver = pp_min_status(use_neumann_solver);
	
	solver.SetMaxIter(40000);
	solver.SetTol(1e-10);

	start_clock("Petsc Solver");
	if (use_neumann_solver)
	{
	    (void) printf("\nUsing Neumann Solver!\n");
	    if (size < 4)
	    {
	    	(void) printf("Isolated small region for solve1d()\n");
		stop_clock("Petsc Solver");
	    }
	    solver.Solve_withPureNeumann();
	    solver.GetNumIterations(&num_iter);
	    solver.GetFinalRelativeResidualNorm(&rel_residual);
	    if(rel_residual > 1)
	    {
		(void) printf("\n The solution diverges! The residual "
		       "is %g. Solve again using GMRES!\n",rel_residual);
		solver.Reset_x();
		solver.Solve_withPureNeumann_GMRES();
		solver.GetNumIterations(&num_iter);
		solver.GetFinalRelativeResidualNorm(&rel_residual);
	    }

	}
	else
	{
	    (void) printf("\nUsing non-Neumann Solver!\n");
	    solver.Solve();
	    solver.GetNumIterations(&num_iter);
	    solver.GetFinalRelativeResidualNorm(&rel_residual);

	    if(rel_residual > 1)
	    {
		(void) printf("\n The solution diverges! The residual "
		       "is %g. Solve again using GMRES!\n",rel_residual);
		solver.Reset_x();
		solver.Solve_GMRES();
		solver.GetNumIterations(&num_iter);
		solver.GetFinalRelativeResidualNorm(&rel_residual);
	    }

	}
	stop_clock("Petsc Solver");

	double *x;
	FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));
	solver.Get_x(x);

	if (debugging("PETSc"))
	    (void) printf("In poisson_solver(): "
	       		"num_iter = %d, rel_residual = %g \n", 
			num_iter, rel_residual);
	
        for (i = imin; i <= imax; i++)
	{
	    index = d_index1d(i,top_gmax);
	    I = i_to_I[i];
	    if (I == -1) continue;
	    soln[index] = x[I-ilower];
	    if (max_soln < soln[index]) max_soln = soln[index];
	    if (min_soln > soln[index]) min_soln = soln[index];
	}
	pp_global_max(&max_soln,1);
	pp_global_min(&min_soln,1);

	if(debugging("step_size"))
	{
	    printf("\nThe max solution value is %.16g\n",max_soln);
	    printf("\nThe min solution value is %.16g\n",min_soln);
	}

	FT_FreeThese(1,x);
}	/* end solve1d */

void ELLIPTIC_SOLVER::solve2d(double *soln)
{
	int index,index_nb[4],size;
	double k0,k_nb[4];
	double rhs,coeff[4];
	int I,I_nb[4];
	int i,j,l,icoords[MAXD];
	COMPONENT comp;
	double aII;
	int num_nb;
	GRID_DIRECTION dir[4] = {WEST,EAST,SOUTH,NORTH};
	boolean use_neumann_solver = YES;
	PetscInt num_iter = 0;
	double rel_residual = 0.0;
	HYPER_SURF *hs;
	double crx_coords[MAXD];
        int status;
	POINTER intfc_state;
	int icrds_max[MAXD],icrds_min[MAXD];

	PETSc solver;
	solver.Create(ilower, iupper-1, 5, 5);
	solver.Reset_A();
	solver.Reset_b();
	solver.Reset_x();
	size = iupper - ilower;
	max_soln = -HUGE;
	min_soln = HUGE;

	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    index  = d_index2d(i,j,top_gmax);
	    comp = top_comp[index];
	    I = ij_to_I[i][j];
	    if (I == -1) continue;

	    index_nb[0] = d_index2d(i-1,j,top_gmax);
	    index_nb[1] = d_index2d(i+1,j,top_gmax);
	    index_nb[2] = d_index2d(i,j-1,top_gmax);
	    index_nb[3] = d_index2d(i,j+1,top_gmax);
	    I_nb[0] = ij_to_I[i-1][j];
	    I_nb[1] = ij_to_I[i+1][j];
	    I_nb[2] = ij_to_I[i][j-1];
	    I_nb[3] = ij_to_I[i][j+1];
	    icoords[0] = i;
	    icoords[1] = j;
	
	    k0 = D[index];
	    num_nb = 0;
	    for (l = 0; l < 4; ++l)
	    {
		status = (*findStateAtCrossing)(front,icoords,dir[l],comp,
                                &intfc_state,&hs,crx_coords);
		if (status != CONST_V_PDE_BOUNDARY)
		    num_nb++;
                if (status == CONST_V_PDE_BOUNDARY ||
		    status == CONST_P_PDE_BOUNDARY)
		    index_nb[l] = index;
		k_nb[l] = 0.5*(k0 + D[index_nb[l]]);
	    	coeff[l] = k_nb[l]/(top_h[l/2]*top_h[l/2]); 
	    }

	    rhs = source[index];

	    aII = 0.0;
	    for (l = 0; l < 4; ++l)
	    {
		if (num_nb == 0) break;
		status = (*findStateAtCrossing)(front,icoords,dir[l],comp,
                                &intfc_state,&hs,crx_coords);
		if (status == NO_PDE_BOUNDARY)
                {
                    solver.Set_A(I,I_nb[l],coeff[l]);
                    aII += -coeff[l];
                }
                else if (status == CONST_P_PDE_BOUNDARY)
                {
		    rhs += -coeff[l]*getStateVar(intfc_state);
                    aII += -coeff[l];
		    use_neumann_solver = NO;
                }
		else if (status == CONST_V_PDE_BOUNDARY &&
                        wave_type(hs) == NEUMANN_BOUNDARY)
                {
                    solver.Set_A(I,I_nb[l],porosity*coeff[l]);
                    aII += -porosity*coeff[l];
                }
	    }
	    /*
	     * This change reflects the need to treat point with only one
	     * interior neighbor (a convex point). Not sure why PETSc cannot
	     * handle such case. If we have better understanding, this should
	     * be changed back.
	     */
	    if(num_nb > 0)
	    {
                solver.Set_A(I,I,aII);
	    }
            else
            {
	    	if (debugging("linear_solver"))
		    (void) printf("WARNING: isolated value!\n");
                solver.Set_A(I,I,1.0);
		rhs = soln[index];
            }
            solver.Set_b(I,rhs);
	}
	use_neumann_solver = pp_min_status(use_neumann_solver);
	
	solver.SetMaxIter(40000);
	solver.SetTol(1e-10);

	start_clock("Petsc Solver");
	if (use_neumann_solver)
	{
	    if (debugging("linear_solver"))
	    	(void) printf("\nUsing Neumann Solver!\n");
	    if (size < 20)
	    {
	    	if (debugging("linear_solver"))
	    	    (void) printf("Isolated small region for solve2d()\n");
		stop_clock("Petsc Solver");
		return;
	    }
	    solver.Solve_withPureNeumann();
	    solver.GetNumIterations(&num_iter);
	    solver.GetFinalRelativeResidualNorm(&rel_residual);
	    if(rel_residual > 1)
	    {
		(void) printf("\n The solution diverges! The residual "
		       "is %g. Solve again using GMRES!\n",rel_residual);
		solver.Reset_x();
		solver.Solve_withPureNeumann_GMRES();
		solver.GetNumIterations(&num_iter);
		solver.GetFinalRelativeResidualNorm(&rel_residual);
	    }

	}
	else
	{
	    if (debugging("linear_solver"))
	        (void) printf("\nUsing non-Neumann Solver!\n");
	    solver.Solve();
	    solver.GetNumIterations(&num_iter);
	    solver.GetFinalRelativeResidualNorm(&rel_residual);

	    if(rel_residual > 1)
	    {
		(void) printf("\n The solution diverges! The residual "
		       "is %g. Solve again using GMRES!\n",rel_residual);
		solver.Reset_x();
		solver.Solve_GMRES();
		solver.GetNumIterations(&num_iter);
		solver.GetFinalRelativeResidualNorm(&rel_residual);
	    }

	}
	stop_clock("Petsc Solver");

	double *x;
	FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));
	solver.Get_x(x);

	if (debugging("PETSc"))
	    (void) printf("In poisson_solver(): "
	       		"num_iter = %d, rel_residual = %g \n", 
			num_iter, rel_residual);
	
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    index = d_index2d(i,j,top_gmax);
	    I = ij_to_I[i][j];
	    if (I == -1) continue;
	    else soln[index] = x[I-ilower];
	    if (max_soln < soln[index]) 
	    {
		icrds_max[0] = i;
		icrds_max[1] = j;
		max_soln = soln[index];
	    }
	    if (min_soln > soln[index]) 
	    {
		icrds_min[0] = i;
		icrds_min[1] = j;
		min_soln = soln[index];
	    }
	}
	pp_global_max(&max_soln,1);
	pp_global_min(&min_soln,1);

	if (debugging("step_size"))
	{
	    (void) printf("Max solution = %20.14f occuring at: %d %d\n",
                        max_soln,icrds_max[0],icrds_max[1]);
            checkSolver(icrds_max,YES);
            (void) printf("Min solution = %20.14f occuring at: %d %d\n",
                        min_soln,icrds_min[0],icrds_min[1]);
            checkSolver(icrds_min,YES);
	}
	if (debugging("elliptic_error"))
        {
            double error,max_error = 0.0;
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                icoords[0] = i;
                icoords[1] = j;
                if (ij_to_I[i][j] == -1) continue;
                error = checkSolver(icoords,NO);
                if (error > max_error)
                {
                    max_error = error;
                    icrds_max[0] = i;
                    icrds_max[1] = j;
                }
            }
            (void) printf("In dual elliptic solver:\n");
            (void) printf("Max relative elliptic error: %20.14f\n",max_error);
            (void) printf("Occuring at (%d %d)\n",icrds_max[0],icrds_max[1]);
            error = checkSolver(icrds_max,YES);
	}

	FT_FreeThese(1,x);
}	/* end solve2d */

void ELLIPTIC_SOLVER::solve3d(double *soln)
{
	int index,index_nb[6],size;
	double k0,k_nb[6];
	double rhs,coeff[6];
	int I,I_nb[6];
	int i,j,k,l,icoords[MAXD],icrds_max[MAXD],icrds_min[MAXD];
	COMPONENT comp;
	double aII;
	int num_nb;
	GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
	boolean use_neumann_solver = YES;
	PetscInt num_iter = 0;
	double rel_residual = 0.0;
	HYPER_SURF *hs;
	double crx_coords[MAXD];
	int status;
	POINTER intfc_state;

	PETSc solver;
	solver.Create(ilower, iupper-1, 7, 7);
	solver.Reset_A();
	solver.Reset_b();
	solver.Reset_x();
	size = iupper - ilower;
	max_soln = -HUGE;
	min_soln = HUGE;

	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    index  = d_index3d(i,j,k,top_gmax);
	    comp = top_comp[index];
	    I = ijk_to_I[i][j][k];
	    if (I == -1) continue;

	    index_nb[0] = d_index3d(i-1,j,k,top_gmax);
	    index_nb[1] = d_index3d(i+1,j,k,top_gmax);
	    index_nb[2] = d_index3d(i,j-1,k,top_gmax);
	    index_nb[3] = d_index3d(i,j+1,k,top_gmax);
	    index_nb[4] = d_index3d(i,j,k-1,top_gmax);
	    index_nb[5] = d_index3d(i,j,k+1,top_gmax);
	    I_nb[0] = ijk_to_I[i-1][j][k];
	    I_nb[1] = ijk_to_I[i+1][j][k];
	    I_nb[2] = ijk_to_I[i][j-1][k];
	    I_nb[3] = ijk_to_I[i][j+1][k];
	    I_nb[4] = ijk_to_I[i][j][k-1];
	    I_nb[5] = ijk_to_I[i][j][k+1];
	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;
	
	    k0 = D[index];
	    num_nb = 0;
	    for (l = 0; l < 6; ++l)
	    {
		status = (*findStateAtCrossing)(front,icoords,dir[l],comp,
                                &intfc_state,&hs,crx_coords);
		if (status != CONST_V_PDE_BOUNDARY)
		    num_nb++;
                if (status == CONST_V_PDE_BOUNDARY ||
		    status == CONST_P_PDE_BOUNDARY)
		    index_nb[l] = index;
		k_nb[l] = 0.5*(k0 + D[index_nb[l]]);
	    	coeff[l] = k_nb[l]/(top_h[l/2]*top_h[l/2]); 
	    }

	    rhs = source[index];

	    aII = 0.0;
	    for (l = 0; l < 6; ++l)
	    {
		if (num_nb == 0) break;
		status = (*findStateAtCrossing)(front,icoords,dir[l],comp,
                                &intfc_state,&hs,crx_coords);
                if (status == NO_PDE_BOUNDARY)
		{
		    solver.Set_A(I,I_nb[l],coeff[l]);
                    aII += -coeff[l];
		}
		else if (status == CONST_P_PDE_BOUNDARY)
		{
                    aII += -coeff[l];
		    rhs += -coeff[l]*getStateVar(intfc_state);
		    use_neumann_solver = NO;
		}
	    }
	    /*
	     * This change reflects the need to treat point with only one
	     * interior neighbor (a convex point). Not sure why PETSc cannot
	     * handle such case. If we have better understanding, this should
	     * be changed back.
	     */
	    if(num_nb > 0)
	    {
                solver.Set_A(I,I,aII);
	    }
            else
            {
		(void) printf("WARNING: isolated value!\n");
                solver.Set_A(I,I,1.0);
		rhs = soln[index];
            }
	    if (num_nb == 1 && rhs < 0.0) rhs *= -1.0;
            solver.Set_b(I,rhs);
	}
	use_neumann_solver = pp_min_status(use_neumann_solver);
	
	solver.SetMaxIter(40000);
	solver.SetTol(1e-10);

	start_clock("Petsc Solver");
	if (use_neumann_solver)
	{
	    if (skip_neumann_solver)
	    {
	    	(void) printf("Skip isolated small region for solve3d()\n");
		stop_clock("Petsc Solver");
		return;
	    }
	    printf("\nUsing Neumann Solver!\n");
	    solver.Solve_withPureNeumann();
	    solver.GetNumIterations(&num_iter);
	    solver.GetFinalRelativeResidualNorm(&rel_residual);
	    if(rel_residual > 1)
	    {
		(void) printf("\n The solution diverges! The residual "
		       "is %g. Solve again using GMRES!\n",rel_residual);
		solver.Reset_x();
		solver.Solve_withPureNeumann_GMRES();
		solver.GetNumIterations(&num_iter);
		solver.GetFinalRelativeResidualNorm(&rel_residual);
	    }

	}
	else
	{
	    (void) printf("\nUsing non-Neumann Solver!\n");
	    solver.Solve();
	    solver.GetNumIterations(&num_iter);
	    solver.GetFinalRelativeResidualNorm(&rel_residual);

	    if(rel_residual > 1)
	    {
		(void) printf("\n The solution diverges! The residual "
		       "is %g. Solve again using GMRES!\n",rel_residual);
		solver.Reset_x();
		solver.Solve_GMRES();
		solver.GetNumIterations(&num_iter);
		solver.GetFinalRelativeResidualNorm(&rel_residual);
	    }

	}
	stop_clock("Petsc Solver");

	double *x;
	FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));
	solver.Get_x(x);

	if (debugging("PETSc"))
	    (void) printf("In poisson_solver(): "
	       		"num_iter = %d, rel_residual = %g \n", 
			num_iter, rel_residual);
	
	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    I = ijk_to_I[i][j][k];
	    if (I == -1) continue;
	    soln[index] = x[I-ilower];
	    if (max_soln < soln[index]) 
	    {
		max_soln = soln[index];
		icrds_max[0] = i;
		icrds_max[1] = j;
		icrds_max[2] = k;
	    }
	    if (min_soln > soln[index]) 
	    {
		min_soln = soln[index];
		icrds_min[0] = i;
		icrds_min[1] = j;
		icrds_min[2] = k;
	    }
	}
	pp_global_max(&max_soln,1);
	pp_global_min(&min_soln,1);

	if (debugging("step_size"))
	{
	    printf("Max solution = %20.14f occuring at: %d %d %d\n",
                        max_soln,icrds_max[0],icrds_max[1],icrds_max[2]);
            checkSolver(icrds_max,YES);
            printf("Min solution = %20.14f occuring at: %d %d %d\n",
                        min_soln,icrds_min[0],icrds_min[1],icrds_min[2]);
            checkSolver(icrds_min,YES);
	}
	if (debugging("elliptic_error"))
        {
            double error,max_error = 0.0;
            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                icoords[0] = i;
                icoords[1] = j;
                icoords[2] = k;
                if (ijk_to_I[i][j][k] == -1) continue;
                error = checkSolver(icoords,NO);
                if (error > max_error)
                {
                    max_error = error;
                    icrds_max[0] = i;
                    icrds_max[1] = j;
                    icrds_max[2] = k;
                }
            }
            (void) printf("In dual elliptic solver:\n");
            (void) printf("Max relative elliptic error: %20.14f\n",max_error);
            (void) printf("Occuring at (%d %d %d)\n",icrds_max[0],
                                icrds_max[1],icrds_max[2]);
            error = checkSolver(icrds_max,YES);
        }

	FT_FreeThese(1,x);
}	/* end solve3d */

void ELLIPTIC_SOLVER::dsolve2d(double *soln)
{
	int index,index_nb,size;
	double rhs,coeff[2][2];
	int I,I_nb;
	int i,j,l,icoords[MAXD],icnb[MAXD];
	int icrds_max[MAXD],icrds_min[MAXD];
	COMPONENT comp;
	double aII;
	GRID_DIRECTION dir[2][2] = {{WEST,EAST},{SOUTH,NORTH}};
	boolean use_neumann_solver = YES;
	PetscInt num_iter = 0;
	double rel_residual = 0.0;
	HYPER_SURF *hs;
	double crx_coords[MAXD];
	int status;
	POINTER intfc_state;
	int idir,nb;
	double h2[MAXD];
	double *x;

	PETSc solver;
	solver.Create(ilower, iupper-1, 9, 9);
	solver.Reset_A();
	solver.Reset_b();
	solver.Reset_x();
	size = iupper - ilower;
	max_soln = -HUGE;
	min_soln = HUGE;

	for (i = 0; i < dim; ++i)
	    h2[i] = 4.0*sqr(top_h[i]);

	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    icoords[0] = i;
	    icoords[1] = j;
	    I = ij_to_I[i][j];

	    index  = d_index(icoords,top_gmax,dim);
	    comp = top_comp[index];
	    if (I == -1) continue;
	
	    rhs = source[index];
	    aII = 0.0;

	    for (idir = 0; idir < dim; ++idir)
	    for (nb = 0; nb < 2; ++nb)
	    {
	    	for (l = 0; l < dim; ++l)
	    	    icnb[l] = icoords[l];
		status = (*findStateAtCrossing)(front,icoords,dir[idir][nb],
				comp,&intfc_state,&hs,crx_coords);
                if (status == NO_PDE_BOUNDARY)
		{
		    icnb[idir] = (nb == 0) ? icoords[idir] - 1 : 
					icoords[idir] + 1;
		    index_nb = d_index(icnb,top_gmax,dim);
		    status = (*findStateAtCrossing)(front,icnb,dir[idir][nb],
				    comp,&intfc_state,&hs,crx_coords);
                    if (status == NO_PDE_BOUNDARY)
		    {
		    	coeff[idir][nb] = D[index_nb]/h2[idir];
			icnb[idir] = (nb == 0) ? icoords[idir] - 2:
                                        icoords[idir] + 2;
			I_nb = ij_to_I[icnb[0]][icnb[1]];
		    	solver.Set_A(I,I_nb,coeff[idir][nb]);
                    	aII += -coeff[idir][nb];
		    }
		    else if (status == CONST_P_PDE_BOUNDARY)
		    {
		    	coeff[idir][nb] = 0.5*(D[index] + D[index_nb])/h2[idir];
                    	aII += -coeff[idir][nb];
			rhs += -coeff[idir][nb]*getStateVar(intfc_state);
			use_neumann_solver = NO;
		    }
		    else if (status == CONST_V_PDE_BOUNDARY)
		    {
		    	coeff[idir][nb] = 0.5*(D[index] + D[index_nb])/h2[idir];
			I_nb = ij_to_I[icnb[0]][icnb[1]];
		    	solver.Set_A(I,I_nb,coeff[idir][nb]);
                    	aII += -coeff[idir][nb];
		    }
		}
		else if (status == CONST_P_PDE_BOUNDARY)
		{
		    coeff[idir][nb] = D[index]/h2[idir];
		    icnb[idir] = (nb == 0) ? icoords[idir] + 1:
                                        icoords[idir] - 1;
		    I_nb = ij_to_I[icnb[0]][icnb[1]];
		    solver.Set_A(I,I_nb,-coeff[idir][nb]);
		    rhs += -coeff[idir][nb]*getStateVar(intfc_state);
		    use_neumann_solver = NO;
		}
	    }
            solver.Set_A(I,I,aII);
            solver.Set_b(I,rhs);
	}
	use_neumann_solver = pp_min_status(use_neumann_solver);
	
	solver.SetMaxIter(40000);
	solver.SetTol(1e-10);

	start_clock("Petsc Solver");
	if (use_neumann_solver)
	{
	    (void) printf("\nUsing Neumann Solver!\n");
	    if (size < 6)
	    {
	    	(void) printf("Isolated small region for dsolve2d()\n");
		stop_clock("Petsc Solver");
		return;
	    }
	    (void) printf("Neumann solver is not yet working for dsolve2d()\n");
	    clean_up(ERROR);
	    solver.Solve_withPureNeumann();
	    solver.GetNumIterations(&num_iter);
	    solver.GetFinalRelativeResidualNorm(&rel_residual);
	    if(rel_residual > 1)
	    {
		(void) printf("\n The solution diverges! The residual "
		       "is %g. Solve again using GMRES!\n",rel_residual);
		solver.Reset_x();
		solver.Solve_withPureNeumann_GMRES();
		solver.GetNumIterations(&num_iter);
		solver.GetFinalRelativeResidualNorm(&rel_residual);
	    }

	}
	else
	{
	    (void) printf("\nUsing non-Neumann Solver!\n");
	    solver.Solve();
	    solver.GetNumIterations(&num_iter);
	    solver.GetFinalRelativeResidualNorm(&rel_residual);

	    if(rel_residual > 1)
	    {
		(void) printf("\n The solution diverges! The residual "
		       "is %g. Solve again using GMRES!\n",rel_residual);
		solver.Reset_x();
		solver.Solve_GMRES();
		solver.GetNumIterations(&num_iter);
		solver.GetFinalRelativeResidualNorm(&rel_residual);
	    }

	}
	stop_clock("Petsc Solver");

	FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));
	solver.Get_x(x);

	if (debugging("PETSc"))
	    (void) printf("In poisson_solver(): "
	       		"num_iter = %d, rel_residual = %g \n", 
			num_iter, rel_residual);
	
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    index = d_index2d(i,j,top_gmax);
	    I = ij_to_I[i][j];
	    if (I == -1) continue;
	    soln[index] = x[I-ilower];
	    if (max_soln < soln[index]) 
	    {
		max_soln = soln[index];
		icrds_max[0] = i;
		icrds_max[1] = j;
	    }
	    if (min_soln > soln[index]) 
	    {
		min_soln = soln[index];
		icrds_min[0] = i;
		icrds_min[1] = j;
	    }
	}
	pp_global_max(&max_soln,1);
	pp_global_min(&min_soln,1);

	if (debugging("step_size"))
	{
	    (void) printf("Max solution = %20.14f occuring at: %d %d\n",
                        max_soln,icrds_max[0],icrds_max[1]);
            dcheckSolver(icrds_max,YES);
            (void) printf("Min solution = %20.14f occuring at: %d %d\n",
                        min_soln,icrds_min[0],icrds_min[1]);
            dcheckSolver(icrds_min,YES);
	}
	if (debugging("elliptic_error"))
        {
            double error,max_error = 0.0;
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                icoords[0] = i;
                icoords[1] = j;
                if (ij_to_I[i][j] == -1) continue;
                error = dcheckSolver(icoords,NO);
                if (error > max_error)
                {
                    max_error = error;
                    icrds_max[0] = i;
                    icrds_max[1] = j;
                }
            }
            (void) printf("In dual elliptic solver:\n");
            (void) printf("Max relative elliptic error: %20.14f\n",max_error);
            (void) printf("Occuring at (%d %d)\n",icrds_max[0],icrds_max[1]);
            error = dcheckSolver(icrds_max,YES);
        }

	FT_FreeThese(1,x);
}	/* end dsolve2d */

void ELLIPTIC_SOLVER::dsolve3d(double *soln)
{
	int index,index_nb,size;
	double rhs,coeff[3][2];
	int I,I_nb;
	int i,j,k,l,icoords[MAXD],icnb[MAXD];
	int icrds_max[MAXD],icrds_min[MAXD];
	COMPONENT comp;
	double aII;
	GRID_DIRECTION dir[3][2] = {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
	boolean use_neumann_solver = YES;
	PetscInt num_iter = 0;
	double rel_residual = 0.0;
	HYPER_SURF *hs;
	double crx_coords[MAXD];
	int status;
	POINTER intfc_state;
	int idir,nb;
	double h2[MAXD];
	double *x;

	PETSc solver;
	solver.Create(ilower, iupper-1, 13, 13);
	solver.Reset_A();
	solver.Reset_b();
	solver.Reset_x();
	size = iupper - ilower;
	max_soln = -HUGE;
	min_soln = HUGE;

	for (i = 0; i < dim; ++i)
	    h2[i] = 4.0*sqr(top_h[i]);

	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;
	    I = ijk_to_I[i][j][k];

	    index  = d_index(icoords,top_gmax,dim);
	    comp = top_comp[index];
	    if (I == -1) continue;
	
	    rhs = source[index];
	    aII = 0.0;

	    for (idir = 0; idir < dim; ++idir)
	    for (nb = 0; nb < 2; ++nb)
	    {
	    	for (l = 0; l < dim; ++l)
	    	    icnb[l] = icoords[l];
		status = (*findStateAtCrossing)(front,icoords,dir[idir][nb],
				comp,&intfc_state,&hs,crx_coords);
                if (status == NO_PDE_BOUNDARY)
		{
		    icnb[idir] = (nb == 0) ? icoords[idir] - 1 : 
					icoords[idir] + 1;
		    index_nb = d_index(icnb,top_gmax,dim);
		    status = (*findStateAtCrossing)(front,icnb,dir[idir][nb],
				    comp,&intfc_state,&hs,crx_coords);
                    if (status == NO_PDE_BOUNDARY)
		    {
		    	coeff[idir][nb] = D[index_nb]/h2[idir];
			icnb[idir] = (nb == 0) ? icoords[idir] - 2:
                                        icoords[idir] + 2;
			I_nb = ijk_to_I[icnb[0]][icnb[1]][icnb[2]];
		    	solver.Set_A(I,I_nb,coeff[idir][nb]);
                    	aII += -coeff[idir][nb];
		    }
		    else if (status == CONST_P_PDE_BOUNDARY)
		    {
		    	coeff[idir][nb] = 0.5*(D[index] + D[index_nb])/h2[idir];
                    	aII += -coeff[idir][nb];
			rhs += -coeff[idir][nb]*getStateVar(intfc_state);
			use_neumann_solver = NO;
		    }
		    else if (status == CONST_V_PDE_BOUNDARY)
		    {
		    	coeff[idir][nb] = 0.5*(D[index] + D[index_nb])/h2[idir];
			I_nb = ijk_to_I[icnb[0]][icnb[1]][icnb[2]];
		    	solver.Set_A(I,I_nb,coeff[idir][nb]);
                    	aII += -coeff[idir][nb];
		    }
		}
		else if (status == CONST_P_PDE_BOUNDARY)
		{
		    coeff[idir][nb] = D[index]/h2[idir];
		    icnb[idir] = (nb == 0) ? icoords[idir] + 1:
                                        icoords[idir] - 1;
		    I_nb = ijk_to_I[icnb[0]][icnb[1]][icnb[2]];
		    solver.Set_A(I,I_nb,-coeff[idir][nb]);
		    rhs += -coeff[idir][nb]*getStateVar(intfc_state);
		    use_neumann_solver = NO;
		}
	    }
	    /*
	     * This change reflects the need to treat point with only one
	     * interior neighbor (a convex point). Not sure why PETSc cannot
	     * handle such case. If we have better understanding, this should
	     * be changed back.
	     */
            solver.Set_A(I,I,aII);
            solver.Set_b(I,rhs);
	}
	use_neumann_solver = pp_min_status(use_neumann_solver);
	
	solver.SetMaxIter(40000);
	solver.SetTol(1e-10);

	start_clock("Petsc Solver");
	if (use_neumann_solver)
	{
	    (void) printf("\nUsing Neumann Solver!\n");
	    if (size < 8)
	    {
	    	(void) printf("Isolated small region for dsolve3d()\n");
		stop_clock("Petsc Solver");
		return;
	    }
	    (void) printf("Neumann solver not working for dsolve3d()\n");
	    clean_up(ERROR);
	    solver.Solve_withPureNeumann();
	    solver.GetNumIterations(&num_iter);
	    solver.GetFinalRelativeResidualNorm(&rel_residual);
	    if(rel_residual > 1)
	    {
		(void) printf("\n The solution diverges! The residual "
		       "is %g. Solve again using GMRES!\n",rel_residual);
		solver.Reset_x();
		solver.Solve_withPureNeumann_GMRES();
		solver.GetNumIterations(&num_iter);
		solver.GetFinalRelativeResidualNorm(&rel_residual);
	    }

	}
	else
	{
	    (void) printf("\nUsing non-Neumann Solver!\n");
	    solver.Solve();
	    solver.GetNumIterations(&num_iter);
	    solver.GetFinalRelativeResidualNorm(&rel_residual);

	    if(rel_residual > 1)
	    {
		(void) printf("\n The solution diverges! The residual "
		       "is %g. Solve again using GMRES!\n",rel_residual);
		solver.Reset_x();
		solver.Solve_GMRES();
		solver.GetNumIterations(&num_iter);
		solver.GetFinalRelativeResidualNorm(&rel_residual);
	    }

	}
	stop_clock("Petsc Solver");

	FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));
	solver.Get_x(x);

	if (debugging("PETSc"))
	    (void) printf("In poisson_solver(): "
	       		"num_iter = %d, rel_residual = %g \n", 
			num_iter, rel_residual);
	
	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    I = ijk_to_I[i][j][k];
	    if (I == -1) continue;
	    soln[index] = x[I-ilower];
	    if (max_soln < soln[index]) 
	    {
		max_soln = soln[index];
		icrds_max[0] = i;
		icrds_max[1] = j;
		icrds_max[2] = k;
	    }
	    if (min_soln > soln[index]) 
	    {
		min_soln = soln[index];
		icrds_min[0] = i;
		icrds_min[1] = j;
		icrds_min[2] = k;
	    }
	}
	pp_global_max(&max_soln,1);
	pp_global_min(&min_soln,1);

	if (debugging("step_size"))
	{
	    (void) printf("Max solution = %20.14f occuring at: %d %d %d\n",
                        max_soln,icrds_max[0],icrds_max[1],icrds_max[2]);
            dcheckSolver(icrds_max,YES);
            (void) printf("Min solution = %20.14f occuring at: %d %d %d\n",
                        min_soln,icrds_min[0],icrds_min[1],icrds_min[2]);
            dcheckSolver(icrds_min,YES);
	}
	if (debugging("elliptic_error"))
        {
            double error,max_error = 0.0;
            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                icoords[0] = i;
                icoords[1] = j;
                icoords[2] = k;
                if (ijk_to_I[i][j][k] == -1) continue;
                error = dcheckSolver(icoords,NO);
                if (error > max_error)
                {
                    max_error = error;
                    icrds_max[0] = i;
                    icrds_max[1] = j;
                    icrds_max[2] = k;
                }
            }
            (void) printf("In dual elliptic solver:\n");
            (void) printf("Max relative elliptic error: %20.14f\n",max_error);
            (void) printf("Occuring at (%d %d %d)\n",icrds_max[0],
                                icrds_max[1],icrds_max[2]);
            error = dcheckSolver(icrds_max,YES);
        }

	FT_FreeThese(1,x);
}	/* end dsolve3d */

double ELLIPTIC_SOLVER::checkSolver(
	int *icoords,
	boolean print_details)
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
		    coefs[m] = 0.5*(D[id0] + D[index_nb]);
		    w[0] = soln[id0];
		    w[1] = soln[index_nb];
		    break;
		case CONST_V_PDE_BOUNDARY:
	    	    if (print_details)
		    	(void) printf("Side %d CONST_V_PDE_BOUNDARY\n",m);
		    coefs[m] = D[id0];
		    w[0] = soln[id0];
		    w[1] = soln[id0];
		    break;
		case CONST_P_PDE_BOUNDARY:
	    	    if (print_details)
		    	(void) printf("Side %d CONST_P_PDE_BOUNDARY\n",m);
		    coefs[m] = D[id0];
		    w[0] = soln[id0];
		    w[1] = getStateVar(intfc_state);
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
	    	(void) printf("Coefs: %f %f\n",coefs[0],coefs[1]);
	    	(void) printf("C*dw: %f %f\n",coefs[0]*dw[0],coefs[1]*dw[1]);
	    }
	    lhs += (coefs[1]*dw[1] + coefs[0]*dw[0])/top_h[l];
	}
	rhs = source[id0];
	if (print_details)
	{
	    (void) printf("Solution = %20.14f\n",soln[id0]);
	    (void) printf("LHS = %20.14f  RHS = %20.14f\n",lhs,rhs);
	    (void) printf("LHS - RHS = %20.14f\n",lhs-rhs);
	    (void) printf("Relative error = %20.14g\n",fabs(lhs-rhs)/denom);
	    (void) printf("Leaving checkSolver()\n\n");
	}
	return fabs(lhs-rhs)/denom;
}	/* end checkSolver */

double ELLIPTIC_SOLVER::dcheckSolver(
	int *icoords,
	boolean print_details)
{
	int i,j,l,m;
	int comp;
	double w[2];
	int id0,index_nb;
	double dw[2],coefs[2],lhs,rhs;
	GRID_DIRECTION dir[3][2] = {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
	HYPER_SURF *hs;
        double crx_coords[MAXD],w_crx;
        int status;
        POINTER intfc_state;
	int icnb[MAXD];
	double denom = 0.0;

	if (print_details)
	{
	    (void) printf("\nEntering dcheckSolver()\n");
	    (void) printf("icoords = ");
	    for (i = 0; i < dim; ++i)
	    	(void) printf("%d ",icoords[i]);
	    (void) printf("\n");
	}

	id0 = d_index(icoords,top_gmax,dim);
	comp = top_comp[id0];
	lhs = 0;
	for (l = 0; l < dim; ++l)
	{
	    if (print_details)
	    	printf("Direction %d:\n",l);
	    for (i = 0; i < dim; ++i)
	    	icnb[i] = icoords[i];
	    for (m = 0; m < 2; ++m)
	    {
		status = (*findStateAtCrossing)(front,icoords,dir[l][m],comp,
                                &intfc_state,&hs,crx_coords);
		if (status == NO_PDE_BOUNDARY)
		{
		    icnb[l] = (m == 0) ? icoords[l] - 1 : icoords[l] + 1;
		    index_nb = d_index(icnb,top_gmax,dim);
		    status = (*findStateAtCrossing)(front,icnb,dir[l][m],comp,
                                &intfc_state,&hs,crx_coords);
		    if (status == NO_PDE_BOUNDARY)
		    {
		    	if (print_details)
		    	    (void) printf("Side %d NO_PDE_BOUNDARY\n",m);
		    	coefs[m] = D[index_nb];
		    	icnb[l] = (m == 0) ? icoords[l] - 2 : icoords[l] + 2;
		    	index_nb = d_index(icnb,top_gmax,dim);
			w[0] = soln[id0];
			w[1] = soln[index_nb];
		    }
		    else if (status == CONST_V_PDE_BOUNDARY)
		    {
		    	if (print_details)
		    	    (void) printf("Side %d-1 CONST_V_PDE_BOUNDARY\n",m);
		    	coefs[m] = 0.5*(D[id0] + D[index_nb]);
			w[0] = soln[id0];
			w[1] = soln[index_nb];
		    }
		    else if (status == CONST_P_PDE_BOUNDARY)
		    {
		    	if (print_details)
		    	    (void) printf("Side %d-1 CONST_P_PDE_BOUNDARY\n",m);
		    	coefs[m] = 0.5*(D[id0] + D[index_nb]);
			w[0] = soln[id0];
		    	w[1] = getStateVar(intfc_state);
		    }
		}
		else if (status == CONST_V_PDE_BOUNDARY)
		{
		    if (print_details)
		    	(void) printf("Side %d-0 CONST_V_PDE_BOUNDARY\n",m);
		    coefs[m] = D[id0];
		    w[0] = soln[id0];
		    w[1] = soln[id0];
		}
		else if (status == CONST_P_PDE_BOUNDARY)
		{
		    if (print_details)
		    	(void) printf("Side %d-0 CONST_P_PDE_BOUNDARY\n",m);
		    icnb[l] = (m == 0) ? icoords[l] + 1 : icoords[l] - 1;
		    index_nb = d_index(icnb,top_gmax,dim);
		    coefs[m] = D[id0];
		    w[0] = soln[index_nb];
		    w[1] = getStateVar(intfc_state);
		}
		dw[m] = (w[1] - w[0])/2.0/top_h[l];
		if (denom < fabs(coefs[m]*dw[m]/2.0/top_h[l]))
		    denom = fabs(coefs[m]*dw[m]/2.0/top_h[l]);
	    }
	    if (print_details)
            {
	    	(void) printf("Coefs: %f %f\n",coefs[0],coefs[1]);
	    	(void) printf("C*dw: %f %f\n",coefs[0]*dw[0],coefs[1]*dw[1]);
	    }
	    lhs += (coefs[1]*dw[1] + coefs[0]*dw[0])/2.0/top_h[l];
	}
	rhs = source[id0];
	if (print_details)
        {
	    (void) printf("Solution = %20.14f\n",soln[id0]);
	    (void) printf("LHS = %20.14f  RHS = %20.14f\n",lhs,rhs);
	    (void) printf("LHS - RHS = %20.14f\n",lhs-rhs);
	    (void) printf("Relative error = %20.14g\n",fabs(lhs-rhs)/denom);
	    (void) printf("Leaving dcheckSolver()\n\n");
	}
	return fabs(lhs-rhs)/denom;
}	/* end dcheckSolver */


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

DUAL_ELLIPTIC_SOLVER::DUAL_ELLIPTIC_SOLVER(Front &front):front(&front)
{
}

void DUAL_ELLIPTIC_SOLVER::set_solver_domain(void)
{
	static boolean first = YES;
	RECT_GRID *rgr = &topological_grid(front->grid_intfc);
	RECT_GRID *crgr = &topological_grid(front->comp_grid_intfc);
	struct Table *T = table_of_interface(front->grid_intfc);
        struct Table *cT = table_of_interface(front->comp_grid_intfc);
        int *lbuf = front->rect_grid->lbuf;
        int *ubuf = front->rect_grid->ubuf;
	int i;

	dim = Dimension(front->interf);
	top_comp = T->components;
	top_gmax = rgr->gmax;
        top_h = rgr->h;
        ctop_comp = cT->components;
        ctop_gmax = crgr->gmax;
	ctop_L = crgr->L;

	if (first)
	{
	    first = NO;
	    cimin = (lbuf[0] == 0) ? 1 : lbuf[0] + 1;
	    cjmin = (lbuf[1] == 0) ? 1 : lbuf[1] + 1;
	    ckmin = (lbuf[2] == 0) ? 1 : lbuf[2] + 1;
	    cimax = (ubuf[0] == 0) ? ctop_gmax[0] - 1 : ctop_gmax[0] - ubuf[0];
            cjmax = (ubuf[1] == 0) ? ctop_gmax[1] - 1 : ctop_gmax[1] - ubuf[1];
            ckmax = (ubuf[2] == 0) ? ctop_gmax[2] - 1 : ctop_gmax[2] - ubuf[2];
	    array_size = 1;
	    for (i = 0; i < dim; ++i)
	    {
		array_size *= (ctop_gmax[i] + 1);
		offset[i] = (lbuf[i] == 0) ? 1 : 0;
	    }
	    FT_VectorMemoryAlloc((POINTER*)&array,array_size,FLOAT);
	}
}	/* end set_solver_domain */

void DUAL_ELLIPTIC_SOLVER::solve(double *soln)
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

void DUAL_ELLIPTIC_SOLVER::solve1d(double *soln)
{
}	/* end solve1d */

void DUAL_ELLIPTIC_SOLVER::solve2d(double *soln)
{
	int index,index_nb[4],size;
	double k_nb[4];
	double rhs,coeff[4];
	int I,I_nb[4];
	int i,j,ii,jj,l,icoords[MAXD];
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

	for (j = cjmin; j <= cjmax; j++)
        for (i = cimin; i <= cimax; i++)
	{
	    index  = d_index2d(i,j,ctop_gmax);
	    comp = ctop_comp[index];
	    I = ij_to_I[i][j];
	    if (I == -1) continue;
	    I_nb[0] = ij_to_I[i-1][j];
	    I_nb[1] = ij_to_I[i+1][j];
	    I_nb[2] = ij_to_I[i][j-1];
	    I_nb[3] = ij_to_I[i][j+1];
	    icoords[0] = i;
	    icoords[1] = j;

	    get_dual_D(icoords,k_nb);

	    num_nb = 0;
	    for (l = 0; l < 4; ++l)
	    {
		status = (*findStateAtCrossing)(front,icoords,dir[l],comp,
                                &intfc_state,&hs,crx_coords);
		if (status != CONST_V_PDE_BOUNDARY)
		    num_nb++;
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
	    if (size < 6)
	    {
	    	(void) printf("Isolated small region for solve2d()\n");
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
	
	for (j = cjmin; j <= cjmax; j++)
        for (i = cimin; i <= cimax; i++)
	{
	    index = d_index2d(i,j,ctop_gmax);
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
	FT_ParallelExchCompGridArrayBuffer(soln,front);

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
	    for (j = cjmin; j <= cjmax; j++)
            for (i = cimin; i <= cimax; i++)
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

void DUAL_ELLIPTIC_SOLVER::solve3d(double *soln)
{
}	/* end solve3d */

void DUAL_ELLIPTIC_SOLVER::get_dual_D(
	int *icoords,
	double *D_nb)
{
	int ii,jj,kk;
	/* The following dual indexing [x][y] are:
	x: 0-west,  1-east;  
	y: 0-south, 1-north; 
	*/
	static int ***ic_dual2d;	
	static int **index_dual2d;
	static COMPONENT **comp_dual2d;
	/* The following dual indexing [x][y][z] are:
	x: 0-west,  1-east;  
	y: 0-south, 1-north; 
	z: 0-lower, 1-upper;
	*/
	static int ****ic_dual3d;	
	static int ***index_dual3d;
	static COMPONENT ***comp_dual3d;

	switch (dim)
	{
	case 2:
	    if (ic_dual2d == NULL)
	    {
	    	FT_TriArrayMemoryAlloc((POINTER*)&ic_dual2d,2,2,MAXD,
					sizeof(INT));
	    	FT_MatrixMemoryAlloc((POINTER*)&index_dual2d,2,2,
					sizeof(INT));
	    	FT_MatrixMemoryAlloc((POINTER*)&comp_dual2d,2,2,
					sizeof(COMPONENT));
	    }
	    for (ii = 0; ii < 2; ++ii)
	    for (jj = 0; jj < 2; ++jj)
	    {
	    	ic_dual2d[ii][jj][0] = icoords[0]+ii+offset[0]-1;
	    	ic_dual2d[ii][jj][1] = icoords[1]+jj+offset[1]-1;
		index_dual2d[ii][jj] = d_index(ic_dual2d[ii][jj],top_gmax,dim);
		comp_dual2d[ii][jj] = top_comp[index_dual2d[ii][jj]];
	    }
	    D_nb[0] = dual_average_D_2d(0,0,index_dual2d,comp_dual2d); // WEST
	    D_nb[1] = dual_average_D_2d(0,1,index_dual2d,comp_dual2d); // EAST
	    D_nb[2] = dual_average_D_2d(1,0,index_dual2d,comp_dual2d); // SOUTH
	    D_nb[3] = dual_average_D_2d(1,1,index_dual2d,comp_dual2d); // NORTH
	    break;
	case 3:
	    if (ic_dual3d == NULL)
	    {
	    	FT_QuadArrayMemoryAlloc((POINTER*)&ic_dual3d,2,2,2,MAXD,
					sizeof(INT));
	    	FT_TriArrayMemoryAlloc((POINTER*)&index_dual3d,2,2,2,
					sizeof(INT));
	    	FT_TriArrayMemoryAlloc((POINTER*)&comp_dual3d,2,2,2,
					sizeof(COMPONENT));
	    }
	    for (ii = 0; ii < 2; ++ii)
	    for (jj = 0; jj < 2; ++jj)
	    for (kk = 0; kk < 2; ++kk)
	    {
	    	ic_dual3d[ii][jj][kk][0] = icoords[0]+ii+offset[0]-1;
	    	ic_dual3d[ii][jj][kk][1] = icoords[1]+jj+offset[1]-1;
	    	ic_dual3d[ii][jj][kk][2] = icoords[2]+kk+offset[2]-1;
		index_dual3d[ii][jj][kk] = d_index(ic_dual3d[ii][jj][kk],
					top_gmax,dim);
		comp_dual3d[ii][jj][kk] = top_comp[index_dual3d[ii][jj][kk]];
	    }
	    D_nb[0] = dual_average_D_3d(0,0,index_dual3d,comp_dual3d); // WEST
	    D_nb[1] = dual_average_D_3d(0,1,index_dual3d,comp_dual3d); // EAST
	    D_nb[2] = dual_average_D_3d(1,0,index_dual3d,comp_dual3d); // SOUTH
	    D_nb[3] = dual_average_D_3d(1,1,index_dual3d,comp_dual3d); // NORTH
	    D_nb[4] = dual_average_D_3d(2,0,index_dual3d,comp_dual3d); // LOWER
	    D_nb[5] = dual_average_D_3d(2,1,index_dual3d,comp_dual3d); // UPPER
	    break;
	}
}	/* end get_dual_D */

double DUAL_ELLIPTIC_SOLVER::dual_average_D_2d(
	int dir,
	int nb,
	int **dual_index,
	COMPONENT **dual_comp)
{
	int i,n = 0;
	double coeff = 0.0;
	switch (dir)
	{
	case 0:
	    for (i = 0; i < 2; ++i)
	    {
		if (dual_comp[nb][i] != obst_comp)
		{
		    coeff += D[dual_index[nb][i]];
		    n++;
		}
	    }
	    if (n != 0) return coeff/n;
	    for (i = 0; i < 2; ++i)
	    {
		if (dual_comp[(nb+1)%2][i] != obst_comp)
		{
		    coeff += D[dual_index[(nb+1)%2][i]];
		    n++;
		}
	    }
	    if (n != 0) return coeff/n;
	    else
	    {
		(void) printf("In dual_average_D(), "
			      "all dual cells are obst_comp\n");
		clean_up(ERROR);
	    }
	    break;
	case 1:
	    for (i = 0; i < 2; ++i)
	    {
		if (dual_comp[i][nb] != obst_comp)
		{
		    coeff += D[dual_index[i][nb]];
		    n++;
		}
	    }
	    if (n != 0) return coeff/n;
	    for (i = 0; i < 2; ++i)
	    {
		if (dual_comp[i][(nb+1)%2] != obst_comp)
		{
		    coeff += D[dual_index[i][(nb+1)%2]];
		    n++;
		}
	    }
	    if (n != 0) return coeff/n;
	    else
	    {
		(void) printf("In dual_average_D(), "
			      "all dual cells are obst_comp\n");
		clean_up(ERROR);
	    }
	    break;
	}
}	/* end dual_average_D_2d */

double DUAL_ELLIPTIC_SOLVER::dual_average_D_3d(
	int dir,
	int nb,
	int ***dual_index,
	COMPONENT ***dual_comp)
{
}	/* end dual_average_D_3d */

double DUAL_ELLIPTIC_SOLVER::checkSolver(
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
	double D_nb[6];

	if (print_details)
	{
	    (void) printf("\nEntering checkSolver()\n");
	    (void) printf("icoords = ");
	    for (i = 0; i < dim; ++i)
	    	(void) printf("%d ",icoords[i]);
	    (void) printf("\n");
	}

	id0 = d_index(icoords,ctop_gmax,dim);
	comp = ctop_comp[id0];
	lhs = 0;
	get_dual_D(icoords,D_nb);
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
		index_nb = d_index(icnb,ctop_gmax,dim);
		coefs[m] = D_nb[l*2+m];

		switch(status)
		{
		case NO_PDE_BOUNDARY:
	    	    if (print_details)
		    	(void) printf("Side %d NO_PDE_BOUNDARY\n",m);
		    w[0] = soln[id0];
		    w[1] = soln[index_nb];
		    break;
		case CONST_V_PDE_BOUNDARY:
	    	    if (print_details)
		    	(void) printf("Side %d CONST_V_PDE_BOUNDARY\n",m);
		    w[0] = soln[id0];
		    w[1] = soln[id0];
		    break;
		case CONST_P_PDE_BOUNDARY:
	    	    if (print_details)
		    	(void) printf("Side %d CONST_P_PDE_BOUNDARY\n",m);
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


/*******************************************************************
 * 		C_CARTESIAN.c
 *******************************************************************/

#include "cim.h"
#include "solver.h"

#define MYZERO 1e-13

//--------------------------------------------------------------------------
// 		C_CARTESIAN
//--------------------------------------------------------------------------

C_CARTESIAN::~C_CARTESIAN()
{
}

void C_CARTESIAN::setIndexMap(void)
{
	static boolean first = YES;
	int i,j,k,ic,index;
	int llbuf[MAXD],uubuf[MAXD];

	if (debugging("trace")) printf("Entering setIndexMap()\n");

	size = iupper - ilower;
	if (first)
	{
	    first = NO;
	    switch (dim)
	    {
	    case 1:
	    	FT_VectorMemoryAlloc((POINTER*)&i_to_I,top_gmax[0]+1,INT);
	    	FT_MatrixMemoryAlloc((POINTER*)&I_to_i,size,1,INT);
	    	break;
	    case 2:
	    	FT_MatrixMemoryAlloc((POINTER*)&ij_to_I,top_gmax[0]+1,
					top_gmax[1]+1,INT);
	    	FT_MatrixMemoryAlloc((POINTER*)&I_to_ij,size,2,INT);
	    	break;
	    case 3:
	    	FT_TriArrayMemoryAlloc((POINTER*)&ijk_to_I,top_gmax[0]+1,
					top_gmax[1]+1,top_gmax[2]+1,INT);
	    	FT_MatrixMemoryAlloc((POINTER*)&I_to_ijk,size,3,INT);
	    	break;
	    }
	}

	index = 0;
	for (i = 0; i < dim; ++i)
	{
	    llbuf[i] = lbuf[i] != 0 ? lbuf[i] : 1;
	    uubuf[i] = ubuf[i] != 0 ? ubuf[i] : 1;
	}
	switch (dim)
	{
	case 1:
	    for (i = 0; i <= top_gmax[0]; i++)
		i_to_I[i] = -1;

	    for (i = imin; i <= imax; i++)
	    {
		ic = d_index1d(i,top_gmax);
	    	i_to_I[i] = index + ilower;
	    	index++;
	    }
	    FT_ParallelExchCellIndex(front,llbuf,uubuf,(POINTER)i_to_I);
	    break;
	case 2:
	    for (i = 0; i <= top_gmax[0]; i++)
	    for (j = 0; j <= top_gmax[1]; j++)
		ij_to_I[i][j] = -1;

	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		ic = d_index2d(i,j,top_gmax);
	    	ij_to_I[i][j] = index + ilower;
		I_to_ij[index + ilower][0] = i;
		I_to_ij[index + ilower][1] = j;
	    	index++;
	    }
	    FT_ParallelExchCellIndex(front,llbuf,uubuf,(POINTER)ij_to_I);
	    break;
	case 3:
	    for (i = 0; i <= top_gmax[0]; i++)
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (k = 0; k <= top_gmax[2]; k++)
		ijk_to_I[i][j][k] = -1;

	    for (i = imin; i <= imax; i++)
	    for (j = jmin; j <= jmax; j++)
	    for (k = kmin; k <= kmax; k++)
	    {
		ic = d_index3d(i,j,k,top_gmax);
	    	ijk_to_I[i][j][k] = index + ilower;
		I_to_ijk[index + ilower][0] = i;
		I_to_ijk[index + ilower][1] = j;
		I_to_ijk[index + ilower][2] = k;
	    	index++;
	    }
	    FT_ParallelExchCellIndex(front,llbuf,uubuf,(POINTER)ijk_to_I);
	    break;
	}
	if (debugging("trace")) printf("Leaving setIndexMap()\n");

}	/* end setIndexMap */

C_CARTESIAN::C_CARTESIAN(Front &front):front(&front)
{
}

void C_CARTESIAN::setDomain()
{
	static boolean first = YES;
	INTERFACE *grid_intfc;
	Table *T;
	int i;

	grid_intfc = front->grid_intfc;
	lbuf = front->rect_grid->lbuf;
	ubuf = front->rect_grid->ubuf;

	top_grid = &topological_grid(grid_intfc);
	top_gmax = top_grid->gmax;
	top_L = top_grid->L;
	top_U = top_grid->U;
	top_h = top_grid->h;

	dim = grid_intfc->dim;
	T = table_of_interface(grid_intfc);
	top_comp = T->components;
	hmin = top_h[0];
	for (i = 1; i < dim; ++i)
	    if (hmin > top_h[i]) hmin = top_h[i];

	switch (dim)
	{
	case 1:
            if (first)
            {
                FT_VectorMemoryAlloc((POINTER*)&array,top_gmax[0]+1,FLOAT);
                FT_VectorMemoryAlloc((POINTER*)&source,top_gmax[0]+1,FLOAT);
                FT_VectorMemoryAlloc((POINTER*)&diff_coeff,top_gmax[0]+1,FLOAT);
                first = NO;
            }	
	    imin = (lbuf[0] == 0) ? 1 : lbuf[0];
	    imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
	    break;
	case 2:
	    if (first)
	    {
	    	FT_VectorMemoryAlloc((POINTER*)&array,
				(top_gmax[0]+1)*(top_gmax[1]+1),FLOAT);
	    	FT_VectorMemoryAlloc((POINTER*)&source,
				(top_gmax[0]+1)*(top_gmax[1]+1),FLOAT);
	    	FT_VectorMemoryAlloc((POINTER*)&diff_coeff,
				(top_gmax[0]+1)*(top_gmax[1]+1),FLOAT);
	    	first = NO;
	    }
	    imin = (lbuf[0] == 0) ? 1 : lbuf[0];
	    jmin = (lbuf[1] == 0) ? 1 : lbuf[1];
	    imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
	    jmax = (ubuf[1] == 0) ? top_gmax[1] - 1 : top_gmax[1] - ubuf[1];
	    break;
	case 3:
	    if (first)
	    {
	    	FT_VectorMemoryAlloc((POINTER*)&array,
			(top_gmax[0]+1)*(top_gmax[1]+1)*(top_gmax[2]+1),FLOAT);
	    	FT_VectorMemoryAlloc((POINTER*)&source,
			(top_gmax[0]+1)*(top_gmax[1]+1)*(top_gmax[2]+1),FLOAT);
	    	FT_VectorMemoryAlloc((POINTER*)&diff_coeff,
			(top_gmax[0]+1)*(top_gmax[1]+1)*(top_gmax[2]+1),FLOAT);
	    	first = NO;
	    }
	    imin = (lbuf[0] == 0) ? 1 : lbuf[0];
	    jmin = (lbuf[1] == 0) ? 1 : lbuf[1];
	    kmin = (lbuf[2] == 0) ? 1 : lbuf[2];
	    imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
	    jmax = (ubuf[1] == 0) ? top_gmax[1] - 1 : top_gmax[1] - ubuf[1];
	    kmax = (ubuf[2] == 0) ? top_gmax[2] - 1 : top_gmax[2] - ubuf[2];
	    break;
	}
}	/* end setDomain */

void C_CARTESIAN::scatMeshArray()
{
	FT_ParallelExchGridArrayBuffer(array,front);
}

void C_CARTESIAN::setGlobalIndex()
{
	int i,j,k,ic;
	int num_nodes = pp_numnodes();
	int myid = pp_mynode();
	static boolean first = YES;

	if (first)
	{
	    first = NO;
	    FT_VectorMemoryAlloc((POINTER*)&n_dist,num_nodes,sizeof(int));
	}
	NLblocks = 0;
	switch (dim)
	{
	case 1:
	    for (i = imin; i <= imax; i++)
	    {
		ic = d_index1d(i,top_gmax);
		NLblocks++;
	    }
	    break;
	case 2:
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		ic = d_index2d(i,j,top_gmax);
		NLblocks++;
	    }
	    break;
	case 3:
	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		ic = d_index3d(i,j,k,top_gmax);
		NLblocks++;
	    }
	    break;
	}

	for (i = 0; i < num_nodes; ++i) n_dist[i] = 0;
	n_dist[myid] = NLblocks;
	pp_global_imax(n_dist,num_nodes);
	ilower = 0;
        iupper = n_dist[0];

        for (i = 1; i <= myid; i++)
        {
            ilower += n_dist[i-1];
            iupper += n_dist[i];
        }	
}

void C_CARTESIAN::solve()
{
        if (debugging("trace")) printf("Entering solve()\n");
        start_clock("solve");

	FT_MakeGridIntfc(front);
        setDomain();
        if (debugging("trace")) printf("Passing setDomain()\n");

        setGlobalIndex();
        if (debugging("trace")) printf("Passing setGlobalIndex()\n");

        setIndexMap();
        if (debugging("trace")) printf("Passing setIndexMap()\n");

	//old_solve();
	cim_solve();

	stop_clock("solve");
}

void C_CARTESIAN::cim_solve()
{
	switch (dim)
	{
	case 2:
	    cim_solve2d();
	    return;
	case 3:
	    cim_solve3d();
	    return;
	}
}	/* end cim_solve */

void C_CARTESIAN::cim_solve2d()
{
	int i,j,k,l,index;
	double P[MAXD],L1,L2,L_inf,err;
	int ic[MAXD];
	static CIM_ELLIPTIC_SOLVER elliptic_solver(*front);
	double *uex;
	CIM_PARAMS *jparams = (CIM_PARAMS*)front->extra1;;

	jparams->dim = dim;
	FT_VectorMemoryAlloc((POINTER*)&uex,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&exact_soln,
		(top_gmax[0]+1)*(top_gmax[1]+1),sizeof(double));
	for (i = imin; i <= imax; ++i)
	for (j = jmin; j <= jmax; ++j)
	{
	    int index,DD;
	    k = ij_to_I[i][j];
    	    for (l = 0; l < dim; ++l) 
	    {
	    	ic[l] = I_to_ij[k][l];
            	P[l] = cell_edge(ic[l],l,top_grid);
            }
	    index = d_index(ic,top_gmax,dim);
	    DD = (top_comp[index] == 2) ? 1 : -1;
            uex[k] = exact_solution((POINTER)jparams,DD,P);
	    source[index] = exact_source((POINTER)jparams,DD,P);
	    diff_coeff[index] = exact_eps(DD);
	    exact_soln[index] = uex[k];
	}
	FT_ParallelExchGridArrayBuffer(exact_soln,front);

	elliptic_solver.w_type = w_type;
	elliptic_solver.neg_comp = neg_comp;
	elliptic_solver.pos_comp = pos_comp;
	elliptic_solver.solutionJump = exact_jump_u;
	elliptic_solver.gradJumpDotN = exact_jump_eps_gradu_dot_n;
	elliptic_solver.gradJumpDotT = exact_jump_gradu_dot_t;
        elliptic_solver.jparams = (POINTER)jparams;
	elliptic_solver.diff_coeff[0] = exact_eps(1);
	elliptic_solver.diff_coeff[1] = exact_eps(-1);
	elliptic_solver.findStateAtCrossing = cim_find_state_at_crossing;
	elliptic_solver.getStateVar = getStateU;
        elliptic_solver.source = source;
        elliptic_solver.ilower = ilower;
        elliptic_solver.iupper = iupper;
        elliptic_solver.soln = array;
	elliptic_solver.ij_to_I = ij_to_I;
	elliptic_solver.I_to_ij = I_to_ij;
	elliptic_solver.size = size;
        elliptic_solver.set_solver_domain();
	elliptic_solver.solve(array);

        L1 = L2 = L_inf = 0.0;

	double u_max,u_min;
	u_max = -HUGE;
	u_min = HUGE;
	printf("Comparison of solutions:\n");
	for (i = imin; i <= imax; ++i)
	for (j = jmin; j <= jmax; ++j)
	{
	    k = ij_to_I[i][j];
	    index = d_index2d(i,j,top_gmax);
	    if (u_max < array[index]) u_max = array[index];
	    if (u_min > array[index]) u_min = array[index];
            err = fabs(array[index] - uex[k]);
	    L1 += err;
	    L2 += sqr(err);
            if (err > L_inf) 
		L_inf = err;
        }
	L1 /= size;
	L2 /= size;
	L2 = sqrt(L2);
        printf("L1 = %18.16f  L2 = %18.16f  L_inf = %18.16f\n",
				L1,L2,L_inf);
	printf("u_max = %f  u_min = %f\n",u_max,u_min);
	numer_soln = array;
        
        FT_FreeThese(1,uex);

}	/* end cim_solve2d */

void C_CARTESIAN::cim_solve3d()
{
	int ii,jj,kk,k,l,index;
	double P[MAXD],L1,L2,L_inf,err;
	int ic[MAXD];
	static CIM_ELLIPTIC_SOLVER elliptic_solver(*front);
	double *uex;
	CIM_PARAMS *jparams = (CIM_PARAMS*)front->extra1;;

	jparams->dim = dim;
	FT_VectorMemoryAlloc((POINTER*)&uex,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&exact_soln,
		(top_gmax[0]+1)*(top_gmax[1]+1)*(top_gmax[2]+1),
		sizeof(double));
	for (ii = imin; ii <= imax; ++ii)
	for (jj = jmin; jj <= jmax; ++jj)
	for (kk = kmin; kk <= kmax; ++kk)
	{
	    int index,DD;
	    k = ijk_to_I[ii][jj][kk];
    	    for (l = 0; l < dim; ++l) 
	    {
	    	ic[l] = I_to_ijk[k][l];
            	P[l] = cell_edge(ic[l],l,top_grid);
            }
	    index = d_index(ic,top_gmax,dim);
	    DD = (top_comp[index] == 2) ? 1 : -1;
            uex[k] = exact_solution((POINTER)jparams,DD,P);
	    source[index] = exact_source((POINTER)jparams,DD,P);
	    diff_coeff[index] = exact_eps(DD);
	    exact_soln[index] = uex[k];
	}
	FT_ParallelExchGridArrayBuffer(exact_soln,front);

	elliptic_solver.w_type = w_type;
	elliptic_solver.neg_comp = neg_comp;
	elliptic_solver.pos_comp = pos_comp;
	elliptic_solver.solutionJump = exact_jump_u;
	elliptic_solver.gradJumpDotN = exact_jump_eps_gradu_dot_n;
	elliptic_solver.gradJumpDotT = exact_jump_gradu_dot_t;
        elliptic_solver.jparams = (POINTER)jparams;
	elliptic_solver.diff_coeff[0] = exact_eps(1);
	elliptic_solver.diff_coeff[1] = exact_eps(-1);
	elliptic_solver.findStateAtCrossing = cim_find_state_at_crossing;
	elliptic_solver.getStateVar = getStateU;
        elliptic_solver.source = source;
        elliptic_solver.ilower = ilower;
        elliptic_solver.iupper = iupper;
        elliptic_solver.soln = array;
	elliptic_solver.ijk_to_I = ijk_to_I;
	elliptic_solver.I_to_ijk = I_to_ijk;
	elliptic_solver.size = size;
        elliptic_solver.set_solver_domain();
	elliptic_solver.solve(array);

        L1 = L2 = L_inf = 0.0;

	double u_max,u_min;
	u_max = -HUGE;
	u_min = HUGE;
	printf("Comparison of solutions:\n");
	for (ii = imin; ii <= imax; ++ii)
	for (jj = jmin; jj <= jmax; ++jj)
	for (kk = kmin; kk <= kmax; ++kk)
	{
	    k = ijk_to_I[ii][jj][kk];
	    index = d_index3d(ii,jj,kk,top_gmax);
	    if (u_max < array[index]) u_max = array[index];
	    if (u_min > array[index]) u_min = array[index];
	    if (ii == 20 && jj == 20)
	    	printf("Solns[%d] = %f %f\n",kk,uex[k],array[index]);
            err = fabs(array[index] - uex[k]);
	    L1 += err;
	    L2 += sqr(err);
            if (err > L_inf) 
		L_inf = err;
        }
	L1 /= size;
	L2 /= size;
	L2 = sqrt(L2);
        printf("L1 = %18.16f  L2 = %18.16f  L_inf = %18.16f\n",
				L1,L2,L_inf);
	printf("u_max = %f  u_min = %f\n",u_max,u_min);
	numer_soln = array;
        
        FT_FreeThese(1,uex);

}	/* end cim_solve3d */

void C_CARTESIAN::old_solve()
{
	int i,j,k,l,index;
	double P[MAXD],L1,L2,L_inf,err;
	int ic[MAXD];
	static ELLIPTIC_SOLVER elliptic_solver(*front);
	double *uex;
	CIM_PARAMS *jparams = (CIM_PARAMS*)front->extra1;

	jparams->dim = dim;
	FT_VectorMemoryAlloc((POINTER*)&uex,size,sizeof(double));
	for (i = imin; i <= imax; ++i)
	for (j = jmin; j <= jmax; ++j)
	{
	    int index,DD;
	    k = ij_to_I[i][j];
    	    for (l = 0; l < dim; ++l) 
	    {
	    	ic[l] = I_to_ij[k][l];
            	P[l] = cell_edge(ic[l],l,top_grid);
            }
	    index = d_index(ic,top_gmax,dim);
	    DD = (top_comp[index] == 2) ? 1 : -1;
            uex[k] = exact_solution((POINTER)jparams,DD,P);
	    source[index] = exact_source((POINTER)jparams,DD,P);
	    diff_coeff[index] = exact_eps(DD);
	}

	elliptic_solver.findStateAtCrossing = cim_find_state_at_crossing;
	elliptic_solver.getStateVar = getStateU;
        elliptic_solver.source = source;
        elliptic_solver.D = diff_coeff;
        elliptic_solver.ilower = ilower;
        elliptic_solver.iupper = iupper;
        elliptic_solver.soln = array;
	elliptic_solver.ij_to_I = ij_to_I;
        elliptic_solver.set_solver_domain();
	elliptic_solver.solve(array);
	viewTopVariable(front,array,NO,0.0,0.0,(char*)"test-cim",
                                (char*)"array");

        L1 = L2 = L_inf = 0.0;

	double u_max,u_min;
	u_max = -HUGE;
	u_min = HUGE;
	printf("Comparison of solutions:\n");
	for (i = imin; i <= imax; ++i)
	for (j = jmin; j <= jmax; ++j)
	{
	    k = ij_to_I[i][j];
	    index = d_index2d(i,j,top_gmax);
	    if (u_max < array[index]) u_max = array[index];
	    if (u_min > array[index]) u_min = array[index];
            err = fabs(array[index] - uex[k]);
	    L1 += err;
	    L2 += sqr(err);
            if (err > L_inf) 
		L_inf = err;
        }
	L1 /= size;
	L2 /= size;
	L2 = sqrt(L2);
        printf("L1 = %18.16f  L2 = %18.16f  L_inf = %18.16f\n",
				L1,L2,L_inf);
	printf("u_max = %f  u_min = %f\n",u_max,u_min);
        
        FT_FreeThese(1,uex);

}	/* end cim_solve */

static double getStateNumSoln(POINTER);
static double getStateExcSoln(POINTER);

void C_CARTESIAN::initMovieVariables()
{
	static HDF_MOVIE_VAR *hdf_movie_var;
	int n;

	if (dim == 3) return;
	if (hdf_movie_var == NULL)
	{
	    FT_ScalarMemoryAlloc((POINTER*)&hdf_movie_var,
				sizeof(HDF_MOVIE_VAR));
	    hdf_movie_var->num_var = n = 0;
	    FT_MatrixMemoryAlloc((POINTER*)&hdf_movie_var->var_name,2,100,
				sizeof(char));
	    FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->top_var,2,
				sizeof(double*));
	    FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->obstacle_comp,2,
                                        sizeof(COMPONENT));
	    FT_VectorMemoryAlloc(
				(POINTER*)&hdf_movie_var->preset_bound,2,
					sizeof(boolean));
	    FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->var_min,2,
				sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->var_max,2,
				sizeof(double));
	    sprintf(hdf_movie_var->var_name[n],"numerical_solution");
	    hdf_movie_var->get_state_var[n] = getStateNumSoln;
	    hdf_movie_var->top_var[n] = numer_soln;
	    hdf_movie_var->preset_bound[n] = NO;
	    hdf_movie_var->num_var = ++n;
	    sprintf(hdf_movie_var->var_name[n],"exact_solution");
	    hdf_movie_var->get_state_var[n] = getStateExcSoln;
	    hdf_movie_var->top_var[n] = exact_soln;
	    hdf_movie_var->preset_bound[n] = NO;
	    hdf_movie_var->num_var = ++n;
	}
	hdf_movie_var->preset_bound[0] = YES;
	hdf_movie_var->var_min[0] = 0.0;
	hdf_movie_var->var_max[0] = 1.0;
	hdf_movie_var->preset_bound[1] = YES;
	hdf_movie_var->var_min[1] = 0.0;
	hdf_movie_var->var_max[1] = 1.0;
	front->hdf_movie_var = hdf_movie_var;
}	/* end initMovieVariables */

static double getStateNumSoln(POINTER state)
{
        STATE *fstate = (STATE*)state;
        return fstate->u;
}       /* end getStateNumSoln */

static double getStateExcSoln(POINTER state)
{
        STATE *fstate = (STATE*)state;
        return fstate->uex;
}       /* end getStateNumSoln */

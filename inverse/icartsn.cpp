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

/*******************************************************************
 * 		C_CARTESIAN.c
 *******************************************************************/

#include "inverse.h"
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
	params = (CIM_PARAMS*)front->extra1;
	hmin = top_h[0];
	for (i = 1; i < dim; ++i)
	    if (hmin > top_h[i]) hmin = top_h[i];

	switch (dim)
	{
	case 1:
            if (first)
            {
		FT_ScalarMemoryAlloc((POINTER*)&field,sizeof(FIELD));
                FT_VectorMemoryAlloc((POINTER*)&array,top_gmax[0]+1,FLOAT);
                FT_VectorMemoryAlloc((POINTER*)&source,top_gmax[0]+1,FLOAT);
                FT_VectorMemoryAlloc((POINTER*)&diff_coeff,top_gmax[0]+1,FLOAT);
                FT_VectorMemoryAlloc((POINTER*)&field->u,top_gmax[0]+1,FLOAT);
                FT_VectorMemoryAlloc((POINTER*)&field->uex,top_gmax[0]+1,FLOAT);
                first = NO;
            }	
	    imin = (lbuf[0] == 0) ? 1 : lbuf[0];
	    imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
	    break;
	case 2:
	    if (first)
	    {
		FT_ScalarMemoryAlloc((POINTER*)&field,sizeof(FIELD));
	    	FT_VectorMemoryAlloc((POINTER*)&array,
				(top_gmax[0]+1)*(top_gmax[1]+1),FLOAT);
	    	FT_VectorMemoryAlloc((POINTER*)&source,
				(top_gmax[0]+1)*(top_gmax[1]+1),FLOAT);
	    	FT_VectorMemoryAlloc((POINTER*)&diff_coeff,
				(top_gmax[0]+1)*(top_gmax[1]+1),FLOAT);
	    	FT_VectorMemoryAlloc((POINTER*)&field->u,
				(top_gmax[0]+1)*(top_gmax[1]+1),FLOAT);
	    	FT_VectorMemoryAlloc((POINTER*)&field->uex,
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
		FT_ScalarMemoryAlloc((POINTER*)&field,sizeof(FIELD));
	    	FT_VectorMemoryAlloc((POINTER*)&array,
			(top_gmax[0]+1)*(top_gmax[1]+1)*(top_gmax[2]+1),FLOAT);
	    	FT_VectorMemoryAlloc((POINTER*)&source,
			(top_gmax[0]+1)*(top_gmax[1]+1)*(top_gmax[2]+1),FLOAT);
	    	FT_VectorMemoryAlloc((POINTER*)&diff_coeff,
			(top_gmax[0]+1)*(top_gmax[1]+1)*(top_gmax[2]+1),FLOAT);
	    	FT_VectorMemoryAlloc((POINTER*)&field->u,
			(top_gmax[0]+1)*(top_gmax[1]+1)*(top_gmax[2]+1),FLOAT);
	    	FT_VectorMemoryAlloc((POINTER*)&field->uex,
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
	params->field = field;
}	/* end setDomain */

void C_CARTESIAN::scatMeshArray()
{
	FT_ParallelExchGridArrayBuffer(array,front,NULL);
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

	cim_solve();

	stop_clock("solve");
}

void C_CARTESIAN::cim_solve()
{
	int ii,jj,kk,k,l,index,DD;
	double P[MAXD];
	static CIM_ELLIPTIC_SOLVER elliptic_solver(*front);
	CIM_PARAMS *params = (CIM_PARAMS*)front->extra1;;
	INTERFACE *intfc = front->interf;
        POINTER sl,sr;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
        double U;

	params->dim = dim;

	switch (dim)
	{
	case 2:
	    for (ii = imin; ii <= imax; ++ii)
	    for (jj = jmin; jj <= jmax; ++jj)
	    {
            	P[0] = cell_edge(ii,0,top_grid);
            	P[1] = cell_edge(jj,1,top_grid);
	    	index = d_index2d(ii,jj,top_gmax);
	    	DD = (top_comp[index] == 2) ? 1 : -1;
	    	source[index] = sourceFunc((POINTER)params,DD,P);
	    	diff_coeff[index] = exact_eps(DD);
	    }
	    elliptic_solver.ij_to_I = ij_to_I;
	    elliptic_solver.I_to_ij = I_to_ij;
	    break;
	case 3:
	    for (ii = imin; ii <= imax; ++ii)
	    for (jj = jmin; jj <= jmax; ++jj)
	    for (kk = kmin; kk <= kmax; ++kk)
	    {
            	P[0] = cell_edge(ii,0,top_grid);
            	P[1] = cell_edge(jj,1,top_grid);
            	P[2] = cell_edge(kk,2,top_grid);
	    	index = d_index3d(ii,jj,kk,top_gmax);
	    	DD = (top_comp[index] == 2) ? 1 : -1;
	    	source[index] = sourceFunc((POINTER)params,DD,P);
	    	diff_coeff[index] = exact_eps(DD);
	    }
	    elliptic_solver.ijk_to_I = ijk_to_I;
	    elliptic_solver.I_to_ijk = I_to_ijk;
	    break;
	}

	elliptic_solver.w_type = w_type;
	elliptic_solver.neg_comp = neg_comp;
	elliptic_solver.pos_comp = pos_comp;
	elliptic_solver.solutionJump = jumpU;
	elliptic_solver.gradJumpDotN = jumpEpsGradDotNorm;
	elliptic_solver.gradJumpDotT = jumpGradDotTan;
	elliptic_solver.exactSolution = exact_solution;
        elliptic_solver.jparams = (POINTER)params;
	elliptic_solver.diff_coeff[0] = exact_eps(1);
	elliptic_solver.diff_coeff[1] = exact_eps(-1);
	elliptic_solver.findStateAtCrossing = cim_find_state_at_crossing;
	elliptic_solver.getStateVar = getStateU;
	elliptic_solver.assignStateVar = assignStateU;
        elliptic_solver.source = source;
        elliptic_solver.ilower = ilower;
        elliptic_solver.iupper = iupper;
        elliptic_solver.soln = array;
	elliptic_solver.size = size;
	elliptic_solver.solve_front_state = YES;
        elliptic_solver.set_solver_domain();
	elliptic_solver.solve(array);

	switch (dim)
	{
	case 2:
	    for (ii = 0; ii <= top_gmax[0]; ++ii)
	    for (jj = 0; jj <= top_gmax[1]; ++jj)
	    {
	    	index = d_index2d(ii,jj,top_gmax);
	    	field->u[index] = array[index];
	    }
	    break;
	case 3:
	    for (ii = 0; ii <= top_gmax[0]; ++ii)
	    for (jj = 0; jj <= top_gmax[1]; ++jj)
	    for (kk = 0; kk <= top_gmax[2]; ++kk)
	    {
	    	index = d_index3d(ii,jj,kk,top_gmax);
	    	field->u[index] = array[index];
	    }
	    break;
	}
}	/* end cim_solve2d */

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
	    hdf_movie_var->top_var[n] = field->u;
	    hdf_movie_var->preset_bound[n] = NO;
	    hdf_movie_var->num_var = ++n;
	    sprintf(hdf_movie_var->var_name[n],"exact_solution");
	    hdf_movie_var->get_state_var[n] = getStateExcSoln;
	    hdf_movie_var->top_var[n] = field->uex;
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

void C_CARTESIAN::readBaseFront(
	CIM_PARAMS *cim_params,
	int i)
{
	char *dir_name = cim_params->base_dir_name;
	int RestartStep = cim_params->steps[i];
	F_BASIC_DATA *f_basic = cim_params->f_basic;
	char restart_state_name[200];

        FT_ScalarMemoryAlloc((POINTER*)&base_front,sizeof(Front));
        f_basic->RestartRun = YES;
	f_basic->size_of_intfc_state = sizeof(STATE);

        sprintf(f_basic->restart_name,"%s/intfc-ts%s",dir_name,
                        right_flush(RestartStep,7));

        FT_StartUp(base_front,f_basic);
	sprintf(restart_state_name,"%s/state.ts%s",dir_name,
                        right_flush(RestartStep,7));
	printf("restart_state_name = %s\n",restart_state_name);
	readBaseStates(restart_state_name);
}	/* end readBaseFront */


void C_CARTESIAN::printFrontInteriorStates(char *out_name)
{
	int i,j,k,index;
	char filename[100];
	FILE *outfile;
	double *u = field->u;

	sprintf(filename,"%s/state.ts%s",out_name,
			right_flush(front->step,7));
#if defined(__MPI__)
	if (pp_numnodes() > 1)
            sprintf(filename,"%s-nd%s",filename,right_flush(pp_mynode(),4));
#endif /* defined(__MPI__) */
	sprintf(filename,"%s-u",filename);
	outfile = fopen(filename,"w");
	
	print_front_states(outfile,front);

	fprintf(outfile,"\nInterior u states:\n");
	switch (dim)
	{
	case 1:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    {
		index = d_index1d(i,top_gmax);
	        fprintf(outfile,"%24.18g\n",u[index]);
	    }
	    break;
	case 2:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    {
		index = d_index2d(i,j,top_gmax);
	        fprintf(outfile,"%24.18g\n",u[index]);
	    }
	    break;
	case 3:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    for (k = 0; k <= top_gmax[2]; ++k)
	    {
		index = d_index3d(i,j,k,top_gmax);
	        fprintf(outfile,"%24.18g\n",u[index]);
	    }
	    break;
	}
	fclose(outfile);
}

void C_CARTESIAN::readBaseStates(
	char *restart_name)
{
	FILE *infile;
	int i,j,k,index;
	char fname[100];
	double *u;
	int *base_gmax;
	RECT_GRID *base_grid;
	static CIM_PARAMS params;

	sprintf(fname,"%s-u",restart_name);
        infile = fopen(fname,"r");

        /* Initialize states in the interior regions */

	read_front_states(infile,base_front);

	FT_MakeGridIntfc(base_front);
	base_grid = &topological_grid(base_front->grid_intfc);
	base_gmax = base_grid->gmax;
	FT_ScalarMemoryAlloc((POINTER*)&params.field,sizeof(FIELD));

	next_output_line_containing_string(infile,"Interior u states:");

	switch (dim)
	{
	case 1:
	    FT_VectorMemoryAlloc((POINTER*)&u,base_gmax[0]+1,FLOAT);
	    for (i = 0; i <= base_gmax[0]; ++i)
	    {
		index = d_index1d(i,base_gmax);
	    	fscanf(infile,"%lf",&u[index]);
	    }
	    break;
	case 2:
	    FT_VectorMemoryAlloc((POINTER*)&u,
                        (base_gmax[0]+1)*(base_gmax[1]+1),FLOAT);
	    for (i = 0; i <= base_gmax[0]; ++i)
	    for (j = 0; j <= base_gmax[1]; ++j)
	    {
		index = d_index2d(i,j,base_gmax);
	    	fscanf(infile,"%lf",&u[index]);
	    }
	    break;
	case 3:
	    FT_VectorMemoryAlloc((POINTER*)&u,(base_gmax[0]+1)*
			(base_gmax[1]+1)*(base_gmax[2]+1),FLOAT);
	    for (i = 0; i <= base_gmax[0]; ++i)
	    for (j = 0; j <= base_gmax[1]; ++j)
	    for (k = 0; k <= base_gmax[2]; ++k)
	    {
		index = d_index3d(i,j,k,base_gmax);
	    	fscanf(infile,"%lf",&u[index]);
	    }
	    break;
	}
	params.field->u = u;
	base_front->extra1 = (POINTER)&params;
	fclose(infile);
}	/* end readBaseStates */

void C_CARTESIAN::compareWithExacySoln()
{
	int ii,jj,kk,index,k,l,DD;
	double P[MAXD],L1,L2,L_inf,err;
	double *uex,*u;
	double u_max,u_min;
	CIM_PARAMS *params = (CIM_PARAMS*)front->extra1;

        L1 = L2 = L_inf = 0.0;

	u_max = -HUGE;
	u_min = HUGE;
	u   = field->u;
	uex = field->uex;

	printf("Comparison of solutions:\n");
	switch (dim)
	{
	case 2:
	    for (ii = imin; ii <= imax; ++ii)
	    for (jj = jmin; jj <= jmax; ++jj)
	    {
            	P[0] = cell_edge(ii,0,top_grid);
            	P[1] = cell_edge(jj,1,top_grid);
	    	index = d_index2d(ii,jj,top_gmax);
	    	DD = (top_comp[index] == 2) ? 1 : -1;
            	uex[index] = exact_solution((POINTER)params,DD,P);
	    	if (u_max < u[index]) u_max = u[index];
	    	if (u_min > u[index]) u_min = u[index];
            	err = fabs(u[index] - uex[index]);
	    	L1 += err;
	    	L2 += sqr(err);
            	if (err > L_inf) 
		    L_inf = err;
            }
	    break;
	case 3:
	    for (ii = imin; ii <= imax; ++ii)
	    for (jj = jmin; jj <= jmax; ++jj)
	    for (kk = kmin; kk <= kmax; ++kk)
	    {
            	P[0] = cell_edge(ii,0,top_grid);
            	P[1] = cell_edge(jj,1,top_grid);
            	P[2] = cell_edge(kk,2,top_grid);
	    	index = d_index3d(ii,jj,kk,top_gmax);
	    	DD = (top_comp[index] == 2) ? 1 : -1;
            	uex[index] = exact_solution((POINTER)params,DD,P);
	    	if (u_max < u[index]) u_max = u[index];
	    	if (u_min > u[index]) u_min = u[index];
            	err = fabs(u[index] - uex[index]);
	    	L1 += err;
	    	L2 += sqr(err);
            	if (err > L_inf) 
		    L_inf = err;
            }
	    break;
	}
	L1 /= size;
	L2 /= size;
	L2 = sqrt(L2);
        printf("L1 = %18.16f  L2 = %18.16f  L_inf = %18.16f\n",
				L1,L2,L_inf);
	printf("u_max = %f  u_min = %f\n",u_max,u_min);
}	/* end compareWithExacySoln */

void C_CARTESIAN::compareWithBaseSoln()
{
	int ii,jj,kk,index,k,l,DD;
	double P[MAXD],L1,L2,L_inf,err;
	double *uex,*u;
	double u_max,u_min;
	CIM_PARAMS *params = (CIM_PARAMS*)front->extra1;
	CIM_PARAMS *base_params;
	FIELD *base_field;
	RECT_GRID *base_grid;
	int *base_gmax;
	int *base_comp;
	int base_ic[MAXD];
	Table *base_T;
	int N,base_index;

        L1 = L2 = L_inf = 0.0;

	u_max = -HUGE;
	u_min = HUGE;
	u   = field->u;
	uex = field->uex;

	readBaseFront(params,0);
	base_params = (CIM_PARAMS*)base_front->extra1;
	base_field = base_params->field;
	base_grid = &topological_grid(base_front->grid_intfc);
	base_gmax = base_grid->gmax;
	base_T = table_of_interface(base_front->grid_intfc);
	base_comp = base_T->components;

	N = 0;
	switch (dim)
	{
	case 2:
	    for (ii = imin; ii <= imax; ++ii)
	    for (jj = jmin; jj <= jmax; ++jj)
	    {
		double R;
            	P[0] = cell_edge(ii,0,top_grid);
            	P[1] = cell_edge(jj,1,top_grid);
	    	index = d_index2d(ii,jj,top_gmax);
		rect_in_which(P,base_ic,base_grid);
		base_index = d_index(base_ic,base_gmax,2);
		if (top_comp[index] != base_comp[base_index])
		    continue;

	    	DD = (top_comp[index] == 2) ? 1 : -1;
		R = sqrt(sqr(P[0]) + sqr(P[1]));
		FT_IntrpStateVarAtCoords(base_front,top_comp[index],
			P,base_field->u,getStateU,&uex[index],NULL);
	    	if (u_max < u[index]) u_max = u[index];
	    	if (u_min > u[index]) u_min = u[index];
            	err = fabs(u[index] - uex[index]);
	    	L1 += err;
	    	L2 += sqr(err);
            	if (err > L_inf) 
		{
		    L_inf = err;
		}
		N++;
            }
	    break;
	case 3:
	    for (ii = imin; ii <= imax; ++ii)
	    for (jj = jmin; jj <= jmax; ++jj)
	    for (kk = kmin; kk <= kmax; ++kk)
	    {
            	P[0] = cell_edge(ii,0,top_grid);
            	P[1] = cell_edge(jj,1,top_grid);
            	P[2] = cell_edge(kk,2,top_grid);
		rect_in_which(P,base_ic,base_grid);
		base_index = d_index(base_ic,base_gmax,3);
		if (top_comp[index] != base_comp[base_index])
		    continue;

	    	index = d_index3d(ii,jj,kk,top_gmax);
	    	DD = (top_comp[index] == 2) ? 1 : -1;
		FT_IntrpStateVarAtCoords(base_front,top_comp[index],
			P,base_field->u,getStateU,&uex[index],NULL);
	    	if (u_max < u[index]) u_max = u[index];
	    	if (u_min > u[index]) u_min = u[index];
            	err = fabs(u[index] - uex[index]);
	    	L1 += err;
	    	L2 += sqr(err);
            	if (err > L_inf) 
		    L_inf = err;
		N++;
            }
	    break;
	}
	L1 /= N;
	L2 /= N;
	L2 = sqrt(L2);
        printf("L1 = %18.16f  L2 = %18.16f  L_inf = %18.16f\n",
				L1,L2,L_inf);
	printf("u_max = %f  u_min = %f\n",u_max,u_min);
}	/* end compareWithBaseSoln */

void C_CARTESIAN::initMesh(void)
{
	int i,j,k, index;
	double crds[MAXD];
	int num_cells;

	// init vertices,edges & cell_center
	if (debugging("trace")) printf("Entering initMesh()\n");
	FT_MakeGridIntfc(front);
	if (debugging("trace")) printf("Passed FT_MakeGridIntfc()\n");
	setDomain();
	if (debugging("trace")) printf("Passed setDomain()\n");
	FT_FreeGridIntfc(front);
	if (debugging("trace")) printf("Leaving initMesh()\n");
}

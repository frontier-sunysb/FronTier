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
 * 		C_CARTESIAN.c
 *******************************************************************/
/*
      To-Do: Line 392, the values of ilower and iupper 
      
*/

#include "crystal.h"
#include "crystal_basic.h"
#include<solver.h>

static int find_state_at_crossing(Front*,int*,GRID_DIRECTION,int,
                                POINTER*,HYPER_SURF**,double*);


//----------------------------------------------------------------
//		C_RECTANGLE
//----------------------------------------------------------------

C_RECTANGLE::C_RECTANGLE(): index(-1), comp(-1)
{
}

void C_RECTANGLE::setCoords(
	double *crds,
	int dim)
{
	int i;
	for (i = 0; i < dim; ++i)
	    coords[i] = crds[i];
}
//--------------------------------------------------------------------------
// 		C_CARTESIAN
//--------------------------------------------------------------------------

C_CARTESIAN::~C_CARTESIAN()
{
}

//---------------------------------------------------------------
//	initMesh
// include the following parts
// 1) setup edges:		
// 2) setup cell_center
//---------------------------------------------------------------
void C_CARTESIAN::initMesh(void)
{
	int i,j,k, index;
	double crds[MAXD];
	int num_cells;
	C_RECTANGLE       rectangle;

	// init vertices,edges & cell_center
	if (debugging("trace")) printf("Entering initMesh()\n");
	FT_MakeGridIntfc(front);
	if (debugging("trace")) printf("Passed FT_MakeGridIntfc()\n");
	setDomain();
	if (debugging("trace")) printf("Passed setDomain()\n");

	num_cells = 1;
	for (i = 0; i < dim; ++i)
	{
	    num_cells *= (top_gmax[i] + 1);
	}
	cell_center.insert(cell_center.end(),num_cells,rectangle);
	
	// setup vertices
	// left to right, down to up
	switch (dim)
	{
	case 1:
	    for (i = 0; i <= top_gmax[0]; i++)
	    {
	    	crds[0] = top_L[0] + top_h[0]*i;
		index = d_index1d(i,top_gmax);
	    	cell_center[index].setCoords(crds,dim);
	    	cell_center[index].icoords[0] = i;
	    }
	case 2:
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {
	    	crds[0] = top_L[0] + top_h[0]*i;
	    	crds[1] = top_L[1] + top_h[1]*j;
		index = d_index2d(i,j,top_gmax);
	    	cell_center[index].setCoords(crds,dim);
	    	cell_center[index].icoords[0] = i;
	    	cell_center[index].icoords[1] = j;
	    }
	    break;
	case 3:
	    for (k = 0; k <= top_gmax[2]; k++)
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {
	    	crds[0] = top_L[0] + top_h[0]*i;
	    	crds[1] = top_L[1] + top_h[1]*j;
	    	crds[2] = top_L[2] + top_h[2]*k;
		index = d_index3d(i,j,k,top_gmax);
	    	cell_center[index].setCoords(crds,dim);
	    	cell_center[index].icoords[0] = i;
	    	cell_center[index].icoords[1] = j;
	    	cell_center[index].icoords[2] = k;
	    }
	}
	setComponent();
	FT_FreeGridIntfc(front);
	if (debugging("trace")) printf("Leaving initMesh()\n");
}

void C_CARTESIAN::setComponent(void)
{
	int i;
	static STATE    *state = NULL;
	double *coords;
	int size = (int)cell_center.size();
	
	// cell center components
	if(state == NULL)
            FT_ScalarMemoryAlloc((POINTER*)&state,sizeof(STATE));

	for (i = 0; i < size; i++)
	{
            coords = cell_center[i].coords;
	    if (cell_center[i].comp != -1 &&
                cell_center[i].comp != top_comp[i])
            {
		if (!FrontNearestIntfcState(front,coords,top_comp[i],
                                (POINTER)state))
                {
                    (void) printf("In setComponent()\n");
                    (void) printf("FrontNearestIntfcState() failed\n");
                    (void) printf("old_comp = %d new_comp = %d\n",
                                        cell_center[i].comp,top_comp[i]);
                    clean_up(ERROR);
                }
	    	field->solute[i] = state->solute;
	    }
	    cell_center[i].comp = top_comp[i];
	}
}

void C_CARTESIAN::setInitialCondition(void)
{
	int i,j;
	COMPONENT c;
	double rho_s;
	double coords[MAXD],distance;
	int size = (int)cell_center.size();
	double **point;
	double ans[MAXD],t[MAXD];
	HYPER_SURF_ELEMENT *hse;
	HYPER_SURF *hs;
	RECT_GRID *gr = front->rect_grid;
	double h;
	int range;

	if (debugging("trace")) 
	    (void) printf("Entering setInitialCondition() for crystal\n");
	/* Initialize states at the interface */
	FT_MakeGridIntfc(front);
	setDomain();
	if (cRparams->reaction_type == DEPOSITION_ONLY)
	{
	    h = gr->h[0];
	    for (j = 0; j < dim; j++)
	    	h = std::min(gr->h[j],h);
	    range = (int)(2*cRparams->gap/h+1);
	} // Needed only for deposition
	for (i = 0; i < size; i++)
	{
	    c = top_comp[i];

	    if (c == SOLUTE_COMP)
	    {
	    	getRectangleCenter(i,coords);
		if (cRparams->reaction_type == DEPOSITION_ONLY)
		{
		    distance = HUGE;
		    if (dim == 1) 
		    {
		    	POINT **p,*pc;
		    	for (p = front->interf->points; p && *p; ++p)
			    if (wave_type(*p) == GROWING_BODY_BOUNDARY)
			    {
			   	pc = *p;
			   	break;
			    }

		    	distance = coords[0] - Coords(pc)[0];
		    }
		    else
		    {
		    	if (FT_FindNearestIntfcPointInRange(front,c,coords,
				ans,t,&hse,&hs,range))
		    	{
			    distance = 0.0;
			    for (j = 0; j < dim; j++)
			    	distance += sqr(coords[j] - ans[j]);
		    	    distance = sqrt(distance);
		    	}
		    }
		    if (distance < cRparams->gap)
	    	    	field->solute[i] = cRparams->C_eq;
	            else
	    	    	field->solute[i] = cRparams->C_0;
		} // No need for dissolution, need for precipitation		
		else
	    	    field->solute[i] = cRparams->C_0;

	    }
	    else if (c == CRYSTAL_COMP)
	    {
	    	getRectangleCenter(i,coords);
		if (cRparams->crystal_dens_func != NULL)
		    rho_s = (*cRparams->crystal_dens_func)(
			    	cRparams->crystal_dens_params,coords);
		else
		    rho_s = cRparams->rho_s;
		field->solute[i] = rho_s;
	    }
	    else
	    	field->solute[i] = 0.0;
	}
	if (debugging("trace")) 
	    (void) printf("Leaving setInitialCondition() for crystal\n");
}	/* end setInitialCondition */

void C_CARTESIAN::setIndexMap(void)
{
	static boolean first = YES;
	int i,j,k,ic,index;
	int llbuf[MAXD],uubuf[MAXD];

	if (debugging("trace")) printf("Entering setIndexMap()\n");
	if (first)
	{
	    first = NO;
	    switch (dim)
	    {
	    case 1:
	    	FT_VectorMemoryAlloc((POINTER*)&i_to_I,top_gmax[0]+1,INT);
	    	break;
	    case 2:
	    	FT_MatrixMemoryAlloc((POINTER*)&ij_to_I,top_gmax[0]+1,
					top_gmax[1]+1,INT);
	    	break;
	    case 3:
	    	FT_TriArrayMemoryAlloc((POINTER*)&ijk_to_I,top_gmax[0]+1,
					top_gmax[1]+1,top_gmax[2]+1,INT);
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
	    for (i = imin; i <= imax; i++)
	    {
		ic = d_index1d(i,top_gmax);
		if (cell_center[ic].comp == SOLUTE_COMP)
		{
	    	    i_to_I[i] = index + ilower;
	    	    index++;
		}
		else
		    i_to_I[i] = -1;
	    }
	    FT_ParallelExchCellIndex(front,llbuf,uubuf,(POINTER)i_to_I);
	    break;
	case 2:
	    for (i = imin; i <= imax; i++)
	    for (j = jmin; j <= jmax; j++)
	    {
		ic = d_index2d(i,j,top_gmax);
		if (cell_center[ic].comp == SOLUTE_COMP)
		{
	    	    ij_to_I[i][j] = index + ilower;
	    	    index++;
		}
		else
		    ij_to_I[i][j] = -1;
	    }
	    FT_ParallelExchCellIndex(front,llbuf,uubuf,(POINTER)ij_to_I);
	    break;
	case 3:
	    for (i = imin; i <= imax; i++)
	    for (j = jmin; j <= jmax; j++)
	    for (k = kmin; k <= kmax; k++)
	    {
		ic = d_index3d(i,j,k,top_gmax);
		if (cell_center[ic].comp == SOLUTE_COMP)
		{
	    	    ijk_to_I[i][j][k] = index + ilower;
	    	    index++;
		}
		else
		    ijk_to_I[i][j][k] = -1;
	    }
	    FT_ParallelExchCellIndex(front,llbuf,uubuf,(POINTER)ijk_to_I);
	    break;
	}
	if (debugging("trace")) printf("Leaving setIndexMap()\n");

}	/* end setIndexMap */

void C_CARTESIAN::computeAdvection(void)
{
	if (debugging("trace")) printf("Entering computeAdvection()\n");
	switch (cRparams->num_scheme)
	{
	case UNSPLIT_EXPLICIT:
	    return computeAdvectionExplicit();
	case UNSPLIT_IMPLICIT:
	    return computeAdvectionImplicit();
	case CRANK_NICOLSON:
	    return computeAdvectionCN();
	}
}
    
void C_CARTESIAN::computeAdvectionCN(void)
{
	int i,j,k,l,m,ic,icn,I,I_nb,icoords[MAXD];
	int gmin[MAXD],ipn[MAXD];
	double crx_coords[MAXD];
	double C0,C_nb,D,rho_s,lambda,coeff,coeff_nb,rhs;
	COMPONENT comp;
	PETSc solver;
	double *x;
        int num_iter = 0;
	int size;
        double rel_residual = 0;
        boolean fr_crx_grid_seg;
        const GRID_DIRECTION dir[3][2] =
                {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
	double coords[MAXD];

	D = cRparams->D;
	rho_s = 0.0;

	if (m_dt < 0.1*sqr(hmin)/D/(double)dim)
	    return computeAdvectionExplicit();

	start_clock("computeAdvectionCN");
	if (debugging("trace")) printf("Entering computeAdvectionCN()\n");

	for (i = 0; i < dim; ++i) gmin[i] = 0;

	setIndexMap();
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
	    	if (comp != SOLUTE_COMP) 
	    	{
		    getRectangleCenter(ic,coords);
		    if (cRparams->crystal_dens_func != NULL)
			rho_s = (*cRparams->crystal_dens_func)(
				    cRparams->crystal_dens_params,coords);
		    else
			rho_s = cRparams->rho_s;
		    array[ic] = rho_s;
	    	    continue;
	        }
		I = i_to_I[i];
                C0 = field->solute[ic];
		rhs = C0;
		coeff = 1.0;
		for (l = 0; l < dim; ++l)
		{
		    lambda = D*m_dt/sqr(top_h[l]);
                    for (m = 0; m < 2; ++m)
                    {
                        next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax);
                        icn = d_index1d(ipn[0],top_gmax);
			I_nb = i_to_I[ipn[0]];
			coeff_nb = -0.5*lambda;
                        fr_crx_grid_seg = FT_StateVarAtGridCrossing(front,
                                        icoords,dir[l][m],comp,getStateSolute,
                                        &C_nb,crx_coords);
                        if (!fr_crx_grid_seg)
                        {
			    solver.Add_A(I,I_nb,coeff_nb);
			    C_nb = field->solute[icn];
			    rhs -= coeff_nb*(C_nb - C0);
		    	    coeff += 0.5*lambda;
                        }
			else
			{
			    rhs -= coeff_nb*2.0*C_nb;
		    	    coeff += lambda;
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
	    	if (comp != SOLUTE_COMP) 
	    	{
		    getRectangleCenter(ic,coords);
		    if (cRparams->crystal_dens_func != NULL)
			rho_s = (*cRparams->crystal_dens_func)(
				    cRparams->crystal_dens_params,coords);
		    else
			rho_s = cRparams->rho_s;
		    array[ic] = rho_s;
	    	    continue;
	        }
                C0 = field->solute[ic];
		rhs = C0;
		coeff = 1.0;
		for (l = 0; l < dim; ++l)
		{
		    lambda = D*m_dt/sqr(top_h[l]);
                    for (m = 0; m < 2; ++m)
                    {
                        next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax);
                        icn = d_index2d(ipn[0],ipn[1],top_gmax);
			I_nb = ij_to_I[ipn[0]][ipn[1]];
			coeff_nb = -0.5*lambda;
                        fr_crx_grid_seg = FT_StateVarAtGridCrossing(front,
                                        icoords,dir[l][m],comp,getStateSolute,
                                        &C_nb,crx_coords);
                        if (!fr_crx_grid_seg)
                        {
			    solver.Add_A(I,I_nb,coeff_nb);
			    C_nb = field->solute[icn];
			    rhs -= coeff_nb*(C_nb - C0);
			    coeff += 0.5*lambda;
                        }
			else
			{
			    rhs -= coeff_nb*2.0*C_nb;
			    coeff += lambda;
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
	    	if (comp != SOLUTE_COMP) 
	    	{
		    getRectangleCenter(ic,coords);
		    if (cRparams->crystal_dens_func != NULL)
			rho_s = (*cRparams->crystal_dens_func)(
				    cRparams->crystal_dens_params,coords);
		    else
			rho_s = cRparams->rho_s;
		    array[ic] = rho_s;
	    	    continue;
	        }
                C0 = field->solute[ic];
		rhs = C0;
		coeff = 1.0;
		for (l = 0; l < dim; ++l)
		{
		    lambda = D*m_dt/sqr(top_h[l]);
                    for (m = 0; m < 2; ++m)
                    {
                        next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax);
                        icn = d_index3d(ipn[0],ipn[1],ipn[2],top_gmax);
			I_nb = ijk_to_I[ipn[0]][ipn[1]][ipn[2]];
			coeff_nb = -0.5*lambda;
                        fr_crx_grid_seg = FT_StateVarAtGridCrossing(front,
                                        icoords,dir[l][m],comp,getStateSolute,
                                        &C_nb,crx_coords);
                        if (!fr_crx_grid_seg)
                        {
			    solver.Add_A(I,I_nb,coeff_nb);
			    C_nb = field->solute[icn];
			    rhs -= coeff_nb*(C_nb - C0);
			    coeff += 0.5*lambda;
                        }
			else
			{
			    rhs -= coeff_nb*2.0*C_nb;
			    coeff += lambda;
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
	    	comp = cell_center[ic].comp;
	    	if (comp == SOLUTE_COMP)
	    	    array[ic] = x[I-ilower];
	    	else
		{
		    getRectangleCenter(ic,coords);
		    if (cRparams->crystal_dens_func != NULL)
			rho_s = (*cRparams->crystal_dens_func)(
				    cRparams->crystal_dens_params,coords);
		    else
			rho_s = cRparams->rho_s;
	    	    array[ic] = rho_s;
		}
	    }
	    break;
        case 2:
            for (i = imin; i <= imax; i++)
            for (j = jmin; j <= jmax; j++)
            {
	    	I = ij_to_I[i][j];
	    	ic = d_index2d(i,j,top_gmax);
	    	comp = cell_center[ic].comp;
	    	if (comp == SOLUTE_COMP)
	    	    array[ic] = x[I-ilower];
	    	else
		{
		    getRectangleCenter(ic,coords);
	    	    if (cRparams->crystal_dens_func != NULL)
			rho_s = (*cRparams->crystal_dens_func)(
				    cRparams->crystal_dens_params,coords);
		    else
			rho_s = cRparams->rho_s;
	    	    array[ic] = rho_s;
		}
	    }
	    break;
        case 3:
            for (i = imin; i <= imax; i++)
            for (j = jmin; j <= jmax; j++)
            for (k = kmin; k <= kmax; k++)
            {
	    	I = ijk_to_I[i][j][k];
	    	ic = d_index3d(i,j,k,top_gmax);
	    	comp = cell_center[ic].comp;
	    	if (comp == SOLUTE_COMP)
	    	    array[ic] = x[I-ilower];
	    	else
		{
		    getRectangleCenter(ic,coords);
	    	    if (cRparams->crystal_dens_func != NULL)
			rho_s = (*cRparams->crystal_dens_func)(
				    cRparams->crystal_dens_params,coords);
		    else
			rho_s = cRparams->rho_s;
	    	    array[ic] = rho_s;
		}
	    }
	    break;
	}
	scatMeshArray();
	switch (dim)
	{
        case 1:
            for (i = 0; i <= top_gmax[0]; ++i)
            {
		ic = d_index1d(i,top_gmax);
		field->solute[ic] = array[ic];
	    }
	    break;
        case 2:
            for (i = 0; i <= top_gmax[0]; ++i)
            for (j = 0; j <= top_gmax[1]; ++j)
            {
		ic = d_index2d(i,j,top_gmax);
		field->solute[ic] = array[ic];
	    }
	    break;
        case 3:
            for (i = 0; i <= top_gmax[0]; ++i)
            for (j = 0; j <= top_gmax[1]; ++j)
            for (k = 0; k <= top_gmax[2]; ++k)
            {
		ic = d_index3d(i,j,k,top_gmax);
		field->solute[ic] = array[ic];
	    }
	    break;
	}
	stop_clock("scatter_data");
	FT_FreeThese(1,x);

	if (debugging("trace")) printf("Leaving computeAdvectionCN()\n");
	stop_clock("computeAdvectionCN");
}	/* end computeAdvectionCN */

void C_CARTESIAN::computeAdvectionImplicit(void)
{
	static PARABOLIC_SOLVER parab_solver(*front);
	static boolean first = YES;

        if (m_dt < 0.1*sqr(hmin)/cRparams->D/(double)dim)
            return computeAdvectionExplicit();

	start_clock("computeAdvectionIM");
        setIndexMap();

	parab_solver.soln_comp = SOLUTE_COMP;
	parab_solver.obst_comp = CRYSTAL_COMP;
	parab_solver.var = field->solute;
	parab_solver.soln = field->solute;
	parab_solver.a = cRparams->field->vel;
	parab_solver.getStateVarFunc = getStateSolute;
	parab_solver.findStateAtCrossing = find_state_at_crossing;
	parab_solver.source = NULL;
	parab_solver.D = cRparams->D;
	parab_solver.nu = NULL;
	parab_solver.order = cRparams->pde_order;
	parab_solver.var_obst = cRparams->rho_s;
	parab_solver.ilower = ilower;
	parab_solver.iupper = iupper;
	parab_solver.dt = m_dt;
	parab_solver.first = first;
	parab_solver.set_solver_domain();
	first = NO;

	switch(dim)
	{
	case 1:
	    parab_solver.i_to_I = i_to_I;
	    break;
	case 2:
	    parab_solver.ij_to_I = ij_to_I;
	    break;
	case 3:
	    parab_solver.ijk_to_I = ijk_to_I;
	    break;
	}
	parab_solver.runge_kutta();
	stop_clock("computeAdvectionIM");

}	/* end computeAdvectionImplicit */

void C_CARTESIAN::solve(double dt)
{

	if (debugging("trace")) printf("Entering crystal solve()\n");
	m_dt = dt;
	start_clock("solve");

	setDomain();
	if (debugging("trace")) printf("Passing setDomain()\n");

	setComponent();
	if (debugging("trace")) printf("Passing setComponent()\n");

	setGlobalIndex();
	if (debugging("trace")) printf("Passing setGlobalIndex()\n");

	if (debugging("sample_solute"))
            sampleSolute();
	computeAdvection();
	if (debugging("trace")) printf("Passing computeAdvection()\n");
	if (debugging("sample_solute"))
            sampleSolute();

	setAdvectionDt();
	if (debugging("trace")) printf("Passing setAdvectionDt()\n");

	setBoundary();
	if (debugging("trace")) printf("Passing setBoundary()\n");

	stop_clock("solve");
	if (debugging("trace")) printf("Leaving crystal solve()\n");
}


void C_CARTESIAN::setAdvectionDt()
{
	double D,k;

	cRparams = (CRT_PARAMS*)front->extra2;
	D = cRparams->D;
	k = cRparams->k;

	if (cRparams->num_scheme == UNSPLIT_EXPLICIT)
	    max_dt = 0.5*sqr(hmin)/D/(double)dim;
	else
	    max_dt = 0.5*hmin/D/(double)dim;

	if (cRparams->point_prop_scheme == EXPLICIT_EULER)
            max_dt = std::min(max_dt,0.25*hmin/k);
	pp_global_min(&max_dt,1);

	if (debugging("trace"))
	{
	    printf("In setAdvectionDt: max_dt = %24.18g\n",max_dt);
	}
}	/* end setAdvectionDt */

void C_CARTESIAN::getRectangleIndex(int index, int &i, int &j)
{
	i = cell_center[index].icoords[0];
	j = cell_center[index].icoords[1];
}

void C_CARTESIAN::getRectangleIndex(int index, int &i, int &j, int &k)
{
	i = cell_center[index].icoords[0];
	j = cell_center[index].icoords[1];
	k = cell_center[index].icoords[2];
}


void C_CARTESIAN::getRectangleCenter(
	int index, 
	double *coords)
{
	int i;
	for (i = 0; i < dim; ++i)
	    coords[i] = cell_center[index].coords[i];
}

void C_CARTESIAN::getRectangleCenter(
	int index0, 
	int index1, 
	double *coords)
{
	int i;
	for (i = 0; i < dim; ++i)
	{
	    coords[i] = 0.5*(cell_center[index0].coords[i] +
	    		     cell_center[index1].coords[i]);
	}
}

void C_CARTESIAN::save(char *filename)
{
	INTERFACE *intfc    = front->interf;
		
	FILE *hfile = fopen(filename, "w");
	if(hfile==NULL)
	{
		printf("\n can't open %s in "
		       "SaveAsTecplot_rect_grid_and_interface().", filename);
		exit(0);
	}
	
	// secondly print out the interface
		
	if(exists_interface(intfc))
	{
	    CURVE		**cur;
	    CURVE		*curve;
	    BOND		*bond;
			
	    for(cur=intfc->curves; cur && *cur; cur++)
	    {
		curve = *cur;
		fprintf(hfile, "ZONE I=%d J=%d F=POINT \n", 
				curve->num_points, 1);
		bond=curve->first;
		fprintf(hfile, "%.4f %.4f \n",bond->start->_coords[0], 
				bond->start->_coords[1]);
		for(bond=curve->first; bond!=NULL; bond=bond->next)
		    fprintf(hfile, "%.4f %.4f \n",bond->end->_coords[0], 
		    		bond->end->_coords[1]);
		}					
	}		
	fclose(hfile);
}

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
	top_grid = &topological_grid(grid_intfc);
	lbuf = front->rect_grid->lbuf;
	ubuf = front->rect_grid->ubuf;
	top_gmax = top_grid->gmax;
	top_L = top_grid->L;
	top_U = top_grid->U;
	top_h = top_grid->h;
	dim = grid_intfc->dim;
	T = table_of_interface(grid_intfc);
	top_comp = T->components;
	cRparams = (CRT_PARAMS*)front->extra2;
	hmin = top_h[0];
	for (i = 1; i < dim; ++i)
	    if (hmin > top_h[i]) hmin = top_h[i];

	switch (dim)
	{
	case 1:
            if (first)
            {
		FT_ScalarMemoryAlloc((POINTER*)&field,sizeof(CRT_FIELD));
                FT_VectorMemoryAlloc((POINTER*)&array,top_gmax[0]+1,FLOAT);
                FT_VectorMemoryAlloc((POINTER*)&field->solute,top_gmax[0]+1,
						FLOAT);
                first = NO;
            }	
	    imin = (lbuf[0] == 0) ? 1 : lbuf[0];
	    imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
	    cRparams->field = field;
	    break;
	case 2:
	    if (first)
	    {
		FT_ScalarMemoryAlloc((POINTER*)&field,sizeof(CRT_FIELD));
	    	FT_VectorMemoryAlloc((POINTER*)&array,
				(top_gmax[0]+1)*(top_gmax[1]+1),FLOAT);
	    	FT_VectorMemoryAlloc((POINTER*)&field->solute,
				(top_gmax[0]+1)*(top_gmax[1]+1),FLOAT);
	    	first = NO;
	    }
	    imin = (lbuf[0] == 0) ? 1 : lbuf[0];
	    jmin = (lbuf[1] == 0) ? 1 : lbuf[1];
	    imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
	    jmax = (ubuf[1] == 0) ? top_gmax[1] - 1 : top_gmax[1] - ubuf[1];
	    cRparams->field = field;
	    break;
	case 3:
	    if (first)
	    {
		FT_ScalarMemoryAlloc((POINTER*)&field,sizeof(CRT_FIELD));
	    	FT_VectorMemoryAlloc((POINTER*)&array,
			(top_gmax[0]+1)*(top_gmax[1]+1)*(top_gmax[2]+1),FLOAT);
	    	FT_VectorMemoryAlloc((POINTER*)&field->solute,
			(top_gmax[0]+1)*(top_gmax[1]+1)*(top_gmax[2]+1),FLOAT);
	    	first = NO;
	    }
	    imin = (lbuf[0] == 0) ? 1 : lbuf[0];
	    jmin = (lbuf[1] == 0) ? 1 : lbuf[1];
	    kmin = (lbuf[2] == 0) ? 1 : lbuf[2];
	    imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
	    jmax = (ubuf[1] == 0) ? top_gmax[1] - 1 : top_gmax[1] - ubuf[1];
	    kmax = (ubuf[2] == 0) ? top_gmax[2] - 1 : top_gmax[2] - ubuf[2];
	    cRparams->field = field;
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
		if (cell_center[ic].comp != SOLUTE_COMP) continue;
		NLblocks++;
	    }
	    break;
	case 2:
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		ic = d_index2d(i,j,top_gmax);
		if (cell_center[ic].comp != SOLUTE_COMP) continue;
		NLblocks++;
	    }
	    break;
	case 3:
	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		ic = d_index3d(i,j,k,top_gmax);
		if (cell_center[ic].comp != SOLUTE_COMP) continue;
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


void C_CARTESIAN::initMovieVariables()
{
	char string[100];
	boolean set_bound = NO;
	double var_min,var_max;
	FILE *infile = fopen(InName(front),"r");

	if (CursorAfterStringOpt(infile,"Type y to set movie bounds:"))
        {
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'Y' || string[0] == 'y')
                set_bound = YES;
        }
	switch (dim)
	{
	case 1:
	    CursorAfterString(infile,"Type y to make movie of solute:");
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'Y' || string[0] == 'y')
	    {
		if (set_bound)
		{
		    CursorAfterString(infile,"Enter min and max solute:");
                    fscanf(infile,"%lf %lf",&var_min,&var_max);
                    (void) printf("%f %f\n",var_min,var_max);
		}
		FT_AddHdfMovieVariable(front,set_bound,YES,CRYSTAL_COMP,
                                "solute",0,field->solute,getStateSolute,
                                var_max,var_min);
	    }
	    break;
	case 2:
	    CursorAfterString(infile,"Type y to make movie of solute:");
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'Y' || string[0] == 'y')
	    {
		if (set_bound)
		{
		    CursorAfterString(infile,"Enter min and max solute:");
                    fscanf(infile,"%lf %lf",&var_min,&var_max);
                    (void) printf("%f %f\n",var_min,var_max);
		}
		FT_AddHdfMovieVariable(front,set_bound,YES,CRYSTAL_COMP,
                                "solute",0,field->solute,getStateSolute,
                                var_max,var_min);
		FT_AddVtkScalarMovieVariable(front,"solute",field->solute);
	    }
	    break;
	case 3:
	    CursorAfterString(infile,"Type y to make yz cross section movie:");
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'Y' || string[0] == 'y')
	    {
		if (set_bound)
		{
		    CursorAfterString(infile,"Enter min and max solute:");
                    fscanf(infile,"%lf %lf",&var_min,&var_max);
                    (void) printf("%f %f\n",var_min,var_max);
		}
		FT_AddHdfMovieVariable(front,set_bound,YES,CRYSTAL_COMP,
                                "solute-yz",0,field->solute,getStateSolute,
                                var_max,var_min);
		FT_AddVtkScalarMovieVariable(front,"solute",field->solute);
	    }
	    CursorAfterString(infile,"Type y to make xz cross section movie:");
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'Y' || string[0] == 'y')
	    {
		if (set_bound)
		{
		    CursorAfterString(infile,"Enter min and max solute:");
                    fscanf(infile,"%lf %lf",&var_min,&var_max);
                    (void) printf("%f %f\n",var_min,var_max);
		}
		FT_AddHdfMovieVariable(front,set_bound,YES,CRYSTAL_COMP,
                                "solute-xz",0,field->solute,getStateSolute,
                                var_max,var_min);
		FT_AddVtkScalarMovieVariable(front,"solute",field->solute);
	    }
	    CursorAfterString(infile,"Type y to make xy cross section movie:");
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'Y' || string[0] == 'y')
	    {
		if (set_bound)
		{
		    CursorAfterString(infile,"Enter min and max solute:");
                    fscanf(infile,"%lf %lf",&var_min,&var_max);
                    (void) printf("%f %f\n",var_min,var_max);
		}
		FT_AddHdfMovieVariable(front,set_bound,YES,CRYSTAL_COMP,
                                "solute-xy",0,field->solute,getStateSolute,
                                var_max,var_min);
		FT_AddVtkScalarMovieVariable(front,"solute",field->solute);
	    }
	}
	fclose(infile);
}	/* end initMovieVariables */

void C_CARTESIAN::printFrontInteriorStates(char *out_name)
{
	int i,j,k,index;
	char filename[100];
	FILE *outfile;
	double *solute = field->solute;

	sprintf(filename,"%s/state.ts%s",out_name,
			right_flush(front->step,7));
#if defined(__MPI__)
	if (pp_numnodes() > 1)
            sprintf(filename,"%s-nd%s",filename,right_flush(pp_mynode(),4));
#endif /* defined(__MPI__) */
	sprintf(filename,"%s-solute",filename);
	outfile = fopen(filename,"w");
	
	solute_print_front_states(outfile,front);

	fprintf(outfile,"\nInterior solute states:\n");
	switch (dim)
	{
	case 1:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    {
		index = d_index1d(i,top_gmax);
	        fprintf(outfile,"%24.18g\n",solute[index]);
	    }
	    break;
	case 2:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    {
		index = d_index2d(i,j,top_gmax);
	        fprintf(outfile,"%24.18g\n",solute[index]);
	    }
	    break;
	case 3:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    for (k = 0; k <= top_gmax[2]; ++k)
	    {
		index = d_index3d(i,j,k,top_gmax);
	        fprintf(outfile,"%24.18g\n",solute[index]);
	    }
	    break;
	}
	fclose(outfile);
}

void C_CARTESIAN::readFrontInteriorStates(char *restart_name)
{
	FILE *infile;
	int i,j,k,index;
	char fname[100];
	double *solute = field->solute;

	sprintf(fname,"%s-solute",restart_name);
        infile = fopen(fname,"r");

        /* Initialize states in the interior regions */

	solute_read_front_states(infile,front);

	FT_MakeGridIntfc(front);
	setDomain();

	next_output_line_containing_string(infile,"Interior solute states:");

	switch (dim)
	{
	case 1:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    {
		index = d_index1d(i,top_gmax);
	    	fscanf(infile,"%lf",&solute[index]);
	    }
	    break;
	case 2:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    {
		index = d_index2d(i,j,top_gmax);
	    	fscanf(infile,"%lf",&solute[index]);
	    }
	    break;
	case 3:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    for (k = 0; k <= top_gmax[2]; ++k)
	    {
		index = d_index3d(i,j,k,top_gmax);
	    	fscanf(infile,"%lf",&solute[index]);
	    }
	    break;
	}
	fclose(infile);
}	/* end readInteriorStates */

void C_CARTESIAN::setBoundary()
{
	int i,j,k,index0,index1;
	INTERFACE *intfc = front->interf;
	double *solute = field->solute;

	switch (dim)
	{
	case 1:
	    if (rect_boundary_type(intfc,0,0) == NEUMANN_BOUNDARY)
	    {
		index0 = d_index1d(0,top_gmax);
		index1 = d_index1d(1,top_gmax);
	    	solute[index0] = solute[index1];
	    }
	    if (rect_boundary_type(intfc,0,1) == NEUMANN_BOUNDARY)
	    {
		index0 = d_index1d(top_gmax[0],top_gmax);
		index1 = d_index1d(top_gmax[0]-1,top_gmax);
	    	solute[index0] = solute[index1];
	    }
	    break;
	case 2:
	    if (rect_boundary_type(intfc,0,0) == NEUMANN_BOUNDARY)
	    {
		for (j = 0; j <= top_gmax[1]; ++j)
		{
		    index0 = d_index2d(0,j,top_gmax);
		    index1 = d_index2d(1,j,top_gmax);
	    	    solute[index0] = solute[index1];
		}
	    }
	    if (rect_boundary_type(intfc,0,1) == NEUMANN_BOUNDARY)
	    {
		for (j = 0; j <= top_gmax[1]; ++j)
		{
		    index0 = d_index2d(top_gmax[0],j,top_gmax);
		    index1 = d_index2d(top_gmax[0]-1,j,top_gmax);
	    	    solute[index0] = solute[index1];
		}
	    }
	    if (rect_boundary_type(intfc,1,0) == NEUMANN_BOUNDARY)
	    {
		for (i = 0; i <= top_gmax[0]; ++i)
		{
		    index0 = d_index2d(i,0,top_gmax);
		    index1 = d_index2d(i,1,top_gmax);
	    	    solute[index0] = solute[index1];
		}
	    }
	    if (rect_boundary_type(intfc,1,1) == NEUMANN_BOUNDARY)
	    {
		for (i = 0; i <= top_gmax[0]; ++i)
		{
		    index0 = d_index2d(i,top_gmax[1],top_gmax);
		    index1 = d_index2d(i,top_gmax[1]-1,top_gmax);
	    	    solute[index0] = solute[index1];
		}
	    }
	    break;
	case 3:
	    if (rect_boundary_type(intfc,0,0) == NEUMANN_BOUNDARY)
	    {
		for (j = 0; j <= top_gmax[1]; ++j)
		for (k = 0; k <= top_gmax[2]; ++k)
		{
		    index0 = d_index3d(0,j,k,top_gmax);
		    index1 = d_index3d(1,j,k,top_gmax);
	    	    solute[index0] = solute[index1];
		}
	    }
	    if (rect_boundary_type(intfc,0,1) == NEUMANN_BOUNDARY)
	    {
		for (j = 0; j <= top_gmax[1]; ++j)
		for (k = 0; k <= top_gmax[2]; ++k)
		{
		    index0 = d_index3d(top_gmax[0],j,k,top_gmax);
		    index1 = d_index3d(top_gmax[0]-1,j,k,top_gmax);
	    	    solute[index0] = solute[index1];
		}
	    }
	    if (rect_boundary_type(intfc,1,0) == NEUMANN_BOUNDARY)
	    {
		for (i = 0; i <= top_gmax[0]; ++i)
		for (k = 0; k <= top_gmax[2]; ++k)
		{
		    index0 = d_index3d(i,0,k,top_gmax);
		    index1 = d_index3d(i,1,k,top_gmax);
	    	    solute[index0] = solute[index1];
		}
	    }
	    if (rect_boundary_type(intfc,1,1) == NEUMANN_BOUNDARY)
	    {
		for (i = 0; i <= top_gmax[0]; ++i)
		for (k = 0; k <= top_gmax[2]; ++k)
		{
		    index0 = d_index3d(i,top_gmax[1],k,top_gmax);
		    index1 = d_index3d(i,top_gmax[1]-1,k,top_gmax);
	    	    solute[index0] = solute[index1];
		}
	    }
	    if (rect_boundary_type(intfc,2,0) == NEUMANN_BOUNDARY)
	    {
		for (i = 0; i <= top_gmax[0]; ++i)
		for (j = 0; j <= top_gmax[1]; ++j)
		{
		    index0 = d_index3d(i,j,0,top_gmax);
		    index1 = d_index3d(i,j,1,top_gmax);
	    	    solute[index0] = solute[index1];
		}
	    }
	    if (rect_boundary_type(intfc,2,1) == NEUMANN_BOUNDARY)
	    {
		for (i = 0; i <= top_gmax[0]; ++i)
		for (j = 0; j <= top_gmax[1]; ++j)
		{
		    index0 = d_index3d(i,j,top_gmax[2],top_gmax);
		    index1 = d_index3d(i,j,top_gmax[2]-1,top_gmax);
	    	    solute[index0] = solute[index1];
		}
	    }
	    break;
	}
}	/* end setBoundary */

void C_CARTESIAN::computeAdvectionExplicit(void)
{
	int i,j,k,l,m,ic,icn,icoords[MAXD];
	int gmin[MAXD],ipn[MAXD];
	double crx_coords[MAXD];
	double solute,solute_nb[2],dgrad[MAXD],grad_plus[MAXD],grad_minus[MAXD];
	double rho_s,coef;
	COMPONENT comp;
	boolean fr_crx_grid_seg;
	const GRID_DIRECTION dir[3][2] = 
		{{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
	double coords[MAXD];
	double v[MAXD],**vel,v_plus[MAXD],v_minus[MAXD];

	start_clock("computeAdvectionExplicit");

	coef = cRparams->D*m_dt;
	rho_s = 0.0;
	vel = cRparams->field->vel;

	for (i = 0; i < dim; ++i) gmin[i] = 0;

	switch (dim)
	{
        case 1:
            for (i = imin; i <= imax; ++i)
            {
                icoords[0] = i;
                ic = d_index1d(i,top_gmax);
                comp = top_comp[ic];
                if (comp != SOLUTE_COMP)
                {
		    getRectangleCenter(ic,coords);
		    if (cRparams->crystal_dens_func != NULL)
			rho_s = (*cRparams->crystal_dens_func)(
				    cRparams->crystal_dens_params,coords);
		    else
			rho_s = cRparams->rho_s;
                    array[ic] = rho_s;
                    continue;
                }
                array[ic] = solute = field->solute[ic];
                for (l = 0; l < dim; ++l)
                {
		    v[l] = 0.0;
		    v_plus[l] = 0.0;
		    v_minus[l] = 0.0;
		}
		if (vel != NULL)
		{
		    for (l = 0; l < dim; ++l)
		    {
			v[l] = vel[l][ic];
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
                                        icoords,dir[l][m],comp,getStateSolute,
                                        &solute_nb[m],crx_coords);
                        if (!fr_crx_grid_seg)
                        {
                             next_ip_in_dir(icoords,dir[l][m],ipn,gmin,
							top_gmax);
                             icn = d_index1d(ipn[0],top_gmax);
			     solute_nb[m] = field->solute[icn];
                        }
                        dgrad[l] += (solute_nb[m] - solute)/top_h[l];
                    }
			
		    grad_plus[l] = (solute_nb[1] - solute)/top_h[l];
		    grad_minus[l] = (solute - solute_nb[0])/top_h[l];
                    
		    array[ic] += coef*dgrad[l]/top_h[l]-m_dt*(v_plus[l]*grad_minus[l]+v_minus[l]*grad_plus[l]);
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
	    	if (comp != SOLUTE_COMP) 
	    	{
		    getRectangleCenter(ic,coords);
		    if (cRparams->crystal_dens_func != NULL)
			rho_s = (*cRparams->crystal_dens_func)(
				    cRparams->crystal_dens_params,coords);
		    else
			rho_s = cRparams->rho_s;
		    array[ic] = rho_s;
	    	    continue;
	        }
                array[ic] = solute = field->solute[ic];
		for (l = 0; l < dim; ++l)
		{
		    v[l] = 0.0;
		    v_plus[l] = 0.0;
		    v_minus[l] = 0.0;
		}
		if (vel != NULL)
		{
		    for (l = 0; l < dim; ++l)
		    {
			v[l] = vel[l][ic];
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
                                        icoords,dir[l][m],comp,getStateSolute,
                                        &solute_nb[m],crx_coords);
                        if (!fr_crx_grid_seg)
                        {
                            next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax);
                            icn = d_index2d(ipn[0],ipn[1],top_gmax);
			    solute_nb[m] = field->solute[icn];
                        }
                        dgrad[l] += (solute_nb[m] - solute)/top_h[l];
                    }
		    grad_plus[l] = (solute_nb[1] - solute)/top_h[l];
		    grad_minus[l] = (solute - solute_nb[0])/top_h[l];

		    array[ic] += coef*dgrad[l]/top_h[l] - m_dt*(v_plus[l]*grad_minus[l]+v_minus[l]*grad_plus[l]);
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
	    	if (comp != SOLUTE_COMP) 
	    	{
		    getRectangleCenter(ic,coords);
		    if (cRparams->crystal_dens_func != NULL)
			rho_s = (*cRparams->crystal_dens_func)(
				    cRparams->crystal_dens_params,coords);
		    else
			rho_s = cRparams->rho_s;
		    array[ic] = rho_s;
	    	    continue;
	        }
                array[ic] = solute = field->solute[ic];
		for (l = 0; l < dim; ++l)
		{
		    v[l] = 0.0;
		    v_plus[l] = 0.0;
		    v_minus[l] = 0.0;
		}
		if (vel != NULL)
		{
		    for ( l = 0; l < dim; ++l)
		    {	
			v[l] = vel[l][ic];
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
                                        icoords,dir[l][m],comp,getStateSolute,
                                        &solute_nb[m],crx_coords);
                        if (!fr_crx_grid_seg)
                        {
                            next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax);
                            icn = d_index3d(ipn[0],ipn[1],ipn[2],top_gmax);
                	    solute_nb[m] = field->solute[icn];
                        }
                        dgrad[l] += (solute_nb[m] - solute)/top_h[l];
                    }

		    grad_plus[l] = (solute_nb[1] - solute)/top_h[l];
		    grad_minus[l] = (solute - solute_nb[0])/top_h[l];

		    array[ic] += coef*dgrad[l]/top_h[l] - m_dt*(v_plus[l]*grad_minus[l]+v_minus[l]*grad_plus[l]);
		}
	    }
	    break;
	}
	scatMeshArray();
	switch (dim)
	{
        case 1:
            for (i = 0; i <= top_gmax[0]; ++i)
            {
		ic = d_index1d(i,top_gmax);
		field->solute[ic] = array[ic];
	    }
	    break;
        case 2:
            for (i = 0; i <= top_gmax[0]; ++i)
            for (j = 0; j <= top_gmax[1]; ++j)
            {
		ic = d_index2d(i,j,top_gmax);
		field->solute[ic] = array[ic];
	    }
	    break;
        case 3:
            for (i = 0; i <= top_gmax[0]; ++i)
            for (j = 0; j <= top_gmax[1]; ++j)
            for (k = 0; k <= top_gmax[2]; ++k)
            {
		ic = d_index3d(i,j,k,top_gmax);
		field->solute[ic] = array[ic];
	    }
	    break;
	}

	stop_clock("computeAdvectionExplicit");
}	/* computeAdvectionExplicit */

void C_CARTESIAN::oneDimPlot(char *outname)
{
	xgraphOneDimPlot(outname);
}	/* end solutePlot */

void C_CARTESIAN::xgraphOneDimPlot(char *outname)
{
	int i,index;
	char filename[100];
	FILE *outfile;

	if (debugging("trace"))
	    printf("Entering xgraphSolute1()\n");
        sprintf(filename,"%s/solute-xg.ts%s",
			outname,right_flush(front->step,7));
#if defined(__MPI__)
        sprintf(filename,"%s-nd%s",filename,right_flush(pp_mynode(),4));
#endif /* defined(__MPI__) */
        outfile = fopen(filename,"w");
	fprintf(outfile,"\"OP at %6.3f\"\n",front->time);
	for (i = 0; i <= top_gmax[0]; ++i)
	{
	    index = d_index1d(i,top_gmax);
	    fprintf(outfile,"%24.18g  %24.18g\n",
	    		cell_center[index].coords[0],
	    		field->solute[index]);
	}
	fclose(outfile);
	if (debugging("trace"))
	    printf("Leaving xgraphSolute1()\n");
}	/* end xgraphOneDimPlot */

void C_CARTESIAN::vtk_plot_concentration2d(
	char *outname)
{

    	std::vector<int> ph_index;
	int i,j,k,index;
	char filename[256];
	FILE *outfile;
	double coord_x,coord_y,coord_z,rho_s,xmin,ymin;
	COMPONENT comp;
	int pointsx,pointsy,num_points,num_cells,num_cell_list;
	int icoords[2],p_gmax[2];
	int count = 0;

	sprintf(filename,"%s/vtk.ts%s",outname,
		right_flush(front->step,7));
	if (pp_numnodes() > 1)
	    sprintf(filename,"%s-nd%s",filename,right_flush(pp_mynode(),4));
	if (create_directory(filename,YES) == FUNCTION_FAILED)
        {
            screen("ERROR in vtk_plot_concentration2d() directory %s ",
                           "doesn't exist ",
                           "and can't be made\n",filename);
            clean_up(ERROR);
        }

	//cell-based liquid phase
	ph_index.clear();
	sprintf(filename,"%s/liquid.vtk",filename);
	outfile = fopen(filename,"w");
	fprintf(outfile,"# vtk DataFile Version 3.0\n");
	fprintf(outfile,"Solute concentration\n");        
	fprintf(outfile,"ASCII\n");
	fprintf(outfile,"DATASET UNSTRUCTURED_GRID\n");

	// No buffer considered to plot
	for (j = jmin; j <= jmax; j++)
	for (i = imin; i <= imax; i++)
	{
	    index = d_index2d(i,j,top_gmax);
	    comp = cell_center[index].comp;
	    if (comp == SOLUTE_COMP)
		ph_index.push_back(index);
	}

	pointsx = imax - imin + 2;
	pointsy = jmax - jmin + 2;
	num_points = pointsx*pointsy;

	num_cells = (int) ph_index.size();
	num_cell_list = 5*num_cells;

	p_gmax[0] = pointsx - 1;
	p_gmax[1] = pointsy - 1;

	index = d_index2d(imin,jmin,top_gmax);
	xmin = cell_center[index].coords[0] - top_h[0]/2.0;
	ymin = cell_center[index].coords[1] - top_h[1]/2.0;

	fprintf(outfile, "POINTS %d double\n", num_points);
	for (j = 0; j < pointsy; j++)
	for (i = 0; i < pointsx; i++)
	{
	    coord_x = xmin + i*top_h[0];
	    coord_y = ymin + j*top_h[1];
	    coord_z = 0.0;
	    fprintf(outfile,"%f %f %f\n",coord_x,coord_y,coord_z);
	}
	
	fprintf(outfile, "CELLS %i %i\n", num_cells, num_cell_list);
	for (i = 0; i < (int)ph_index.size(); i++)
	{
	    int index0,index1,index2,index3;
	    index = ph_index[i];
	    icoords[0] = cell_center[index].icoords[0];
	    icoords[1] = cell_center[index].icoords[1];
	    index0 = d_index2d(icoords[0]-4,icoords[1]-1,p_gmax);
	    index1 = d_index2d(icoords[0]-4+1,icoords[1]-1,p_gmax);
	    index2 = d_index2d(icoords[0]-4,icoords[1]-1+1,p_gmax);
	    index3 = d_index2d(icoords[0]-4+1,icoords[1]-1+1,p_gmax);

	    fprintf(outfile, "4 %i %i %i %i\n",
	    	index0,index1,index2,index3);
	}
	
	fprintf(outfile, "CELL_TYPES %i\n", num_cells);
	for (i = 0; i < num_cells; i++)
	    fprintf(outfile, "8\n");

	fprintf(outfile, "CELL_DATA %i\n", num_cells);
	fprintf(outfile,"SCALARS concentration double\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	for (i = 0; i < num_cells; i++)
	{
	    index = ph_index[i];
	    fprintf(outfile,"%f\n",cRparams->field->solute[index]/cRparams->C_eq);
	}

	fclose(outfile);

	//cell-based solid phase
	sprintf(filename,"%s/vtk.ts%s",outname,
		right_flush(front->step,7));
	if (pp_numnodes() > 1)
	    sprintf(filename,"%s-nd%s",filename,right_flush(pp_mynode(),4));
	
	ph_index.clear();
	sprintf(filename,"%s/solid.vtk",filename);
	outfile = fopen(filename,"w");
	fprintf(outfile,"# vtk DataFile Version 3.0\n");
	fprintf(outfile,"Solid density\n");
	fprintf(outfile,"ASCII\n");
	fprintf(outfile,"DATASET UNSTRUCTURED_GRID\n");

	for (j = jmin; j <= jmax; j++)
	for (i = imin; i <= imax; i++)
	{
	    index = d_index2d(i,j,top_gmax);
	    if (cell_center[index].comp == CRYSTAL_COMP)
		ph_index.push_back(index);
	}

	pointsx = imax - imin + 2;
	pointsy = jmax - jmin + 2;
	num_points = pointsx*pointsy;

	num_cells = (int) ph_index.size();
	num_cell_list = 5*num_cells;

	p_gmax[0] = pointsx - 1;
	p_gmax[1] = pointsy - 1;

	index = d_index2d(imin,jmin,top_gmax);
	xmin = cell_center[index].coords[0] - top_h[0]/2.0;
	ymin = cell_center[index].coords[1] - top_h[1]/2.0;

	fprintf(outfile,"POINTS %d double\n", num_points);
	for (j = 0; j < pointsy; j++)
	for (i = 0; i < pointsx; i++)
	{
	    coord_x = xmin + i*top_h[0];
	    coord_y = ymin + j*top_h[1];
	    coord_z = 0.0;
	    fprintf(outfile,"%f %f %f\n",coord_x,coord_y,coord_z);
	}

	fprintf(outfile,"CELLS %i %i\n", num_cells,num_cell_list);
	for (i = 0; i < num_cells; i++)
	{
	    int index0,index1,index2,index3;
	    index = ph_index[i];
	    icoords[0] = cell_center[index].icoords[0];
	    icoords[1] = cell_center[index].icoords[1];
	    index0 = d_index2d(icoords[0]-4,icoords[1]-1,p_gmax); 
	    index1 = d_index2d(icoords[0]-4+1,icoords[1]-1,p_gmax);
	    index2 = d_index2d(icoords[0]-4,icoords[1]-1+1,p_gmax);
	    index3 = d_index2d(icoords[0]-4+1,icoords[1]-1+1,p_gmax);
		
	    fprintf(outfile,"4 %i %i %i %i\n",
		    index0,index1,index2,index3);
	}

	fprintf(outfile, "CELL_TYPES %i\n", num_cells);
	for (i = 0; i < num_cells; i++)
	    fprintf(outfile, "8\n");
	
	fprintf(outfile, "CELL_DATA %i\n", num_cells);
	fprintf(outfile, "SCALARS density double\n");
	fprintf(outfile, "LOOKUP_TABLE default\n");
	for (i = 0; i < num_cells; i++)
	{
	    index = ph_index[i];
	    if (cRparams->crystal_dens_func != NULL)
		rho_s = (*cRparams->crystal_dens_func)(
			cRparams->crystal_dens_params,
			cell_center[index].coords);
	    else
		rho_s = cRparams->rho_s;
	    fprintf(outfile,"%f\n",rho_s);
	    
	    if (rho_s <= 0.6)
		count++;
	}

	fclose(outfile);
}       /* end vtk_plot_concentration2d */   

static int find_state_at_crossing(
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

        if (wave_type(*hs) == FIRST_PHYSICS_WAVE_TYPE) 
	    return NO_PDE_BOUNDARY;
	else if (wave_type(*hs) == DIRICHLET_BOUNDARY)
	{
            return DIRICHLET_PDE_BOUNDARY;
	}
	else if (wave_type(*hs) == GROWING_BODY_BOUNDARY)
	{
            return DIRICHLET_PDE_BOUNDARY;
	}
	else if (wave_type(*hs) == NEUMANN_BOUNDARY)
            return NEUMANN_PDE_BOUNDARY;
}       /* find_state_at_crossing */

void C_CARTESIAN::initSampleSolute(char *in_name)
{
        FILE *infile;
        static SAMPLE *sample;
        char *sample_type;
        double *sample_line;

        infile = fopen(in_name,"r");
        FT_ScalarMemoryAlloc((POINTER*)&sample,sizeof(SAMPLE));
        sample_type = sample->sample_type;
        sample_line = sample->sample_coords;

        if (dim == 2)
        {
            CursorAfterString(infile,"Enter the sample line type:");
            fscanf(infile,"%s",sample_type);
            (void) printf(" %s\n",sample_type);
            CursorAfterString(infile,"Enter the sample line coordinate:");
            fscanf(infile,"%lf",sample_line);
            (void) printf(" %f\n",sample_line[0]);
        }
        else if (dim == 3)
        {
            CursorAfterString(infile,"Enter the sample line type:");
            fscanf(infile,"%s",sample_type);
            (void) printf(" %s\n",sample_type);
            CursorAfterString(infile,"Enter the sample line coordinate:");
            fscanf(infile,"%lf %lf",sample_line,sample_line+1);
            (void) printf(" %f %f\n",sample_line[0],sample_line[1]);
        }
        CursorAfterString(infile,"Enter the start step for sample: ");
        fscanf(infile,"%d",&sample->start_step);
        (void) printf("%d\n",sample->start_step);
        CursorAfterString(infile,"Enter the end step for sample: ");
        fscanf(infile,"%d",&sample->end_step);
        (void) printf("%d\n",sample->end_step);
        CursorAfterString(infile,"Enter the step interval for sample: ");
        fscanf(infile,"%d",&sample->step_interval);
        (void) printf("%d\n",sample->step_interval);
        front->sample = sample;
        fclose(infile);
}       /* end initSampleSolute */

void C_CARTESIAN::sampleSolute()
{
        if (dim == 2)
            sampleSolute2d();
        else if (dim == 3)
            sampleSolute3d();
}       /* end sampleConcentration */

void C_CARTESIAN::sampleSolute2d()
{
	int i,j,index;
        SAMPLE *sample = front->sample;
        char *sample_type = sample->sample_type;
        double *line = sample->sample_coords;
        char *out_name = front->out_name;
        double coords[MAXD];
        double var1,var2,var;
        FILE *sfile;
        char sname[100];
        static int count = 0;
        static int step = 0;
        static int l = -1;
        static double lambda;
        char dirname[256];
        static char **sample_color;

	if (sample_color == NULL)
        {
            FT_MatrixMemoryAlloc((POINTER*)&sample_color,10,20,sizeof(char));
            sprintf(sample_color[0],"red");
            sprintf(sample_color[1],"blue");
            sprintf(sample_color[2],"green");
            sprintf(sample_color[3],"violet");
            sprintf(sample_color[4],"orange");
            sprintf(sample_color[5],"yellow");
            sprintf(sample_color[6],"pink");
            sprintf(sample_color[7],"cyan");
            sprintf(sample_color[8],"light-gray");
            sprintf(sample_color[9],"dark-gray");
        }
	if (pp_numnodes() > 1)
            return;
        if (front->step < sample->start_step || front->step > sample->end_step)
            return;
        if ((front->step - sample->start_step)%sample->step_interval)
            return;
        if (step != front->step)
        {
            step = front->step;
            count = 0;
        }
        sprintf(dirname, "%s/sample-%d", out_name,step);
        if (!create_directory(dirname,NO))
        {
            screen("Cannot create directory %s\n",dirname);
            clean_up(ERROR);
        }
	switch (sample_type[0])
        {
        case 'x':
            if (l == -1)
            {
                double x1,x2;
                do
                {
                    ++l;
                    index = d_index2d(l,0,top_gmax);
                    getRectangleCenter(index, coords);
                } while(line[0] >= coords[0]);
                --l;
                index = d_index2d(l,0,top_gmax);
                getRectangleCenter(index,coords);
                x1 = coords[0];
                index = d_index2d(l+1,0,top_gmax);
                getRectangleCenter(index,coords);
                x2 = coords[0];
                lambda = (line[0] - x1) / (x2 - line[0]);
            }
            i = l;
            sprintf(sname, "%s/solute-%d.xg",dirname,count);
            sfile = fopen(sname,"w");
            fprintf(sfile,"Next\n");
            fprintf(sfile,"color=%s\n",sample_color[count]);
            fprintf(sfile,"thickness=1.5\n");
            for (j = jmin; j <= jmax; ++j)
            {
                index = d_index2d(i,j,top_gmax);
                var1 = field->solute[index];
                index = d_index2d(i+1,j,top_gmax);
                var2 = field->solute[index];
                var = (var1 + lambda*var2) / (1.0 + lambda);
                getRectangleCenter(index,coords);
                fprintf(sfile,"%20.14f %20.14f\n",coords[1],var);
            }
            fclose(sfile);
            break;
	case 'y':
            if (l == -1)
            {
                double y1,y2;
                do
                {
                    ++l;
                    index = d_index2d(0,l,top_gmax);
                    getRectangleCenter(index, coords);
                } while (line[0] >= coords[1]);
                --l;
                index = d_index2d(0,l,top_gmax);
                getRectangleCenter(index,coords);
                y1 = coords[1];
                index = d_index2d(0,l+1,top_gmax);
                getRectangleCenter(index,coords);
                y2 = coords[1];
                lambda = (line[0] - y1) / (y2 - line[0]);
            }
            j = l;
            sprintf(sname, "%s/solute-%d.xg",dirname,count);
            sfile = fopen(sname,"w");
            fprintf(sfile,"Next\n");
            fprintf(sfile,"color=%s\n",sample_color[count]);
            fprintf(sfile,"thickness=1.5\n");
            for (i = imin; i <= imax; ++i)
            {
                index = d_index2d(i,j,top_gmax);
                var1 = field->solute[index];
                index = d_index2d(i,j+1,top_gmax);
                var2 = field->solute[index];
                var = (var1 + lambda*var2) / (1.0 + lambda);
                getRectangleCenter(index,coords);
                fprintf(sfile,"%20.14f   %20.14f\n",coords[0],var);
            }
            fclose(sfile);
            break;
        default:
            printf("Incorrect input for sample solute!\n");
            break;
        }
        count++;
}       /* sampleSolute2d */

void C_CARTESIAN::sampleSolute3d()
{
	int i,j,k,index;
        SAMPLE *sample = front->sample;
        char *sample_type = sample->sample_type;
        double *sample_line = sample->sample_coords;
        char *out_name = front->out_name;
        double coords[MAXD];
        double var1,var2,var_tmp1,var_tmp2,var;
        FILE *sfile;
        char sname[100];
        static int count = 0;
        static int step = 0;
        static int l = -1, m = -1;
        static double lambda1,lambda2;
        char dirname[256];
        static char **sample_color;

	if (sample_color == NULL)
        {
            FT_MatrixMemoryAlloc((POINTER*)&sample_color,10,20,sizeof(char));
            sprintf(sample_color[0],"red");
            sprintf(sample_color[1],"blue");
            sprintf(sample_color[2],"green");
            sprintf(sample_color[3],"violet");
            sprintf(sample_color[4],"orange");
            sprintf(sample_color[5],"yellow");
            sprintf(sample_color[6],"pink");
            sprintf(sample_color[7],"cyan");
            sprintf(sample_color[8],"light-gray");
            sprintf(sample_color[9],"dark-gray");
        }
        if (pp_numnodes() > 1)
            return;
        if (front->step < sample->start_step || front->step > sample->end_step)
            return;
        if ((front->step - sample->start_step)%sample->step_interval)
            return;
        if (step != front->step)
        {
            step = front->step;
            count = 0;
        }
        sprintf(dirname, "%s/sample-%d", out_name,step);
        if (!create_directory(dirname,NO))
        {
            screen("Cannot create directory %s\n",dirname);
            clean_up(ERROR);
        }

	switch (sample_type[0])
        {
        case 'x':
            if (l == -1)
            {
                double x1,x2;
                do
                {
                    ++l;
                    index = d_index3d(l,0,0,top_gmax);
                    getRectangleCenter(index,coords);
                }while(sample_line[0]>=coords[0]);
                --l;
                index = d_index3d(l,0,0,top_gmax);
                getRectangleCenter(index,coords);
                x1 = coords[0];
                index = d_index3d(l+1,0,0,top_gmax);
                getRectangleCenter(index,coords);
                x2 = coords[0];
                lambda1 = (sample_line[0] - x1) / (x2 - sample_line[0]);
            }
	    switch (sample_type[1])
            {
                case 'y':
                    if (m == -1)
                    {
                        double y1,y2;
                        do
                        {
                            ++m;
                            index = d_index3d(0,m,0,top_gmax);
                            getRectangleCenter(index,coords);
                        }while(sample_line[1]>=coords[1]);
                        --m;
                        index = d_index3d(0,m,0,top_gmax);
                        getRectangleCenter(index,coords);
                        y1 = coords[1];
                        index = d_index3d(0,m+1,0,top_gmax);
                        getRectangleCenter(index,coords);
                        y2 = coords[1];
                        lambda2 = (sample_line[1] - y1)/(y2 - sample_line[1]);
                    }
                    i = l;
                    j = m;
		    sprintf(sname, "%s/solute-%d.xg",dirname,count);
                    sfile = fopen(sname,"w");
                    fprintf(sfile,"Next\n");
                    fprintf(sfile,"color=%s\n",sample_color[count]);
                    fprintf(sfile,"thickness=1.5\n");
                    for (k = kmin; k <= kmax; ++k)
                    {
                        index = d_index3d(i,j,k,top_gmax);
                        var1 = field->solute[index];
                        index = d_index3d(i+1,j,k,top_gmax);
                        var2 = field->solute[index];
                        var_tmp1 = (var1 + lambda1*var2)/(1.0 + lambda1);

                        index = d_index3d(i,j+1,k,top_gmax);
                        var1 = field->solute[index];
                        index = d_index3d(i+1,j+1,k,top_gmax);
                        var2 = field->solute[index];
                        var_tmp2 = (var1 + lambda1*var2)/(1.0 + lambda1);

                        var = (var_tmp1 + lambda2*var_tmp2)/(1.0 + lambda2);
                        getRectangleCenter(index,coords);
                        fprintf(sfile,"%20.14f   %20.14f\n",coords[2],var);
                    }
                    fclose(sfile);
                    break;
		case 'z':
                    if (m == -1)
                    {
                        double z1,z2;
                        do
                        {
                            ++m;
                            index = d_index3d(0,0,m,top_gmax);
                            getRectangleCenter(index,coords);
                        }while(sample_line[1]>=coords[2]);
                        --m;
                        index = d_index3d(0,0,m,top_gmax);
                        getRectangleCenter(index,coords);
                        z1 = coords[2];
                        index = d_index3d(0,0,m+1,top_gmax);
                        getRectangleCenter(index,coords);
                        z2 = coords[2];
                        lambda2 = (sample_line[1] - z1)/(z2 - sample_line[1]);
                    }
                    i = l;
                    k = m;
		    sprintf(sname, "%s/solute-%d.xg",dirname,count);
                    sfile = fopen(sname,"w");
                    fprintf(sfile,"Next\n");
                    fprintf(sfile,"color=%s\n",sample_color[count]);
                    fprintf(sfile,"thickness=1.5\n");
                    for (j = jmin; j <= jmax; ++j)
                    {
                        index = d_index3d(i,j,k,top_gmax);
                        var1 = field->solute[index];
                        index = d_index3d(i+1,j,k,top_gmax);
                        var2 = field->solute[index];
                        var_tmp1 = (var1 + lambda1*var2)/(1.0 + lambda1);

                        index = d_index3d(i,j,k+1,top_gmax);
                        var1 = field->solute[index];
                        index = d_index3d(i+1,j,k+1,top_gmax);
                        var2 = field->solute[index];
                        var_tmp2 = (var1 + lambda1*var2)/(1.0 + lambda1);

                        var = (var_tmp1 + lambda2*var_tmp2)/(1.0 + lambda2);
                        getRectangleCenter(index,coords);
                        fprintf(sfile,"%20.14f   %20.14f\n",coords[2],var);
                    }
                    fclose(sfile);
                    break;
                default:
                    printf("Incorrect input for sample solute!\n");
                    break;
            }
	    break;
	case 'y':
            if (l == -1)
            {
                double y1,y2;
                do
                {
                    ++l;
                    index = d_index3d(0,l,0,top_gmax);
                    getRectangleCenter(index, coords);
                }while(sample_line[0]>=coords[1]);
                --l;
                index = d_index3d(0,l,0,top_gmax);
                getRectangleCenter(index,coords);
                y1 = coords[1];
                index = d_index3d(0,l+1,0,top_gmax);
                getRectangleCenter(index,coords);
                y2 = coords[1];
                lambda1 = (sample_line[0] - y1)/(y2 - sample_line[0]);
            }
	    switch (sample_type[1])
            {
                case 'z':
                    if (m == -1)
                    {
                        double z1,z2;
                        do
                        {
                            ++m;
                            index = d_index3d(0,0,m,top_gmax);
                            getRectangleCenter(index,coords);
                        }while(sample_line[1]>=coords[2]);
                        --m;
                        index = d_index3d(0,0,m,top_gmax);
                        getRectangleCenter(index,coords);
                        z1 = coords[2];
                        index = d_index3d(0,0,m+1,top_gmax);
                        getRectangleCenter(index,coords);
                        z2 = coords[2];
                        lambda2 = (sample_line[1] - z1)/(z2 - sample_line[1]);
                    }
                    j = l;
                    k = m;
		    sprintf(sname, "%s/solute-%d.xg",dirname,count);
                    sfile = fopen(sname,"w");
                    fprintf(sfile,"Next\n");
                    fprintf(sfile,"color=%s\n",sample_color[count]);
                    fprintf(sfile,"thickness=1.5\n");
                    for (i = imin; i <= imax; ++i)
                    {
                        index = d_index3d(i,j,k,top_gmax);
                        var1 = field->solute[index];
                        index = d_index3d(i,j+1,k,top_gmax);
                        var2 = field->solute[index];
                        var_tmp1 = (var1 + lambda1*var2)/(1.0 + lambda1);

                        index = d_index3d(i,j,k+1,top_gmax);
                        var1 = field->solute[index];
                        index = d_index3d(i,j+1,k+1,top_gmax);
                        var2 = field->solute[index];
                        var_tmp2 = (var1 + lambda1*var2)/(1.0 + lambda1);

                        var = (var_tmp1 + lambda2*var_tmp2)/(1.0 + lambda2);
                        getRectangleCenter(index,coords);
                        fprintf(sfile,"%20.14f   %20.14f\n",coords[0],var);
                    }
                    fclose(sfile);
                    break;
		default:
                    printf("Incorrect input for sample solute!\n");
                    break;
            }
        default:
            printf("Incorrect input for sample solute!\n");
            break;
        }
        count++;
}       /* end sampleSolute3d */

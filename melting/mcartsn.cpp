/*******************************************************************
 * 		CARTESIAN.c
 *******************************************************************/
/*
      To-Do: Line 392, the values of ilower and iupper 
      
*/

#include <iFluid.h>
#include "solver.h"
#include "melting.h"

static int find_state_at_crossing(Front*,int*,GRID_DIRECTION,int,
                                POINTER*,HYPER_SURF**,double*);

//----------------------------------------------------------------
//		RECTANGLE
//----------------------------------------------------------------

//RECTANGLE::RECTANGLE()
RECTANGLE::RECTANGLE(): index(-1), comp(-1)
{
}

void RECTANGLE::setCoords(
	double *crds,
	int dim)
{
	int i;
	for (i = 0; i < dim; ++i)
	    coords[i] = crds[i];
}
//--------------------------------------------------------------------------
// 		CARTESIAN
//--------------------------------------------------------------------------

CARTESIAN::~CARTESIAN()
{
}

//---------------------------------------------------------------
//	initMesh
// include the following parts
// 1) setup edges:		
// 2) setup cell_center
//---------------------------------------------------------------
void CARTESIAN::initMesh(void)
{
	int i,j,k, index;
	double crds[MAXD];
	int icoords[MAXD];
	int num_cells;
	int cell_index;
	RECTANGLE       rectangle;

	// init vertices,edges & cell_center
	if (debugging("trace")) printf("Entering initMesh()\n");
	FT_MakeGridIntfc(front);
	setDomain();

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

void CARTESIAN::setComponent(void)
{
	int i;

        static STATE *state = NULL;
        double *coords;
        int *icoords;
        int size = (int)cell_center.size();
        HYPER_SURF_ELEMENT *hse;
        HYPER_SURF *hs;
        double t[MAXD],point[MAXD];
        int n;

        if (state == NULL)
            FT_ScalarMemoryAlloc((POINTER*)&state, sizeof(STATE));

        for (i = 0; i < size; i++)
        {
            coords = cell_center[i].coords;
            if (cell_center[i].comp != -1 &&
                cell_center[i].comp != top_comp[i])
            {
                if (FT_FindNearestIntfcPointInRange(front,top_comp[i],coords,
                                point,t,&hse,&hs,2))
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
                    field->temperature[i] = state->temperature;
                }
		else
                {
                    double temp_nb = 0.0;
                    int ii,jj,ic[MAXD],index;
                    icoords = cell_center[i].icoords;
                    n = 0;
                    for (ii = 0; ii < dim; ++ii)
                    {
                        for (jj = 0; jj < dim; ++jj) ic[jj] = icoords[jj];
                        ic[ii] = (icoords[ii] == 0) ? 0 : icoords[ii] - 1;
                        index = d_index(ic,top_gmax,dim);
                        if (cell_center[index].comp != -1 &&
                            cell_center[index].comp == top_comp[index])
                        {
                            temp_nb += field->temperature[index];
                            n++;
                        }
                        ic[ii] = (icoords[ii] == top_gmax[ii]) ? top_gmax[ii]
                                        : icoords[ii] + 1;
                        index = d_index(ic,top_gmax,dim);
                        if (cell_center[index].comp != -1 &&
                            cell_center[index].comp == top_comp[index])
                        {
                            temp_nb += field->temperature[index];
                            n++;
                        }
                    }
                    field->temperature[i] = temp_nb/n;
                }
	    }
	    cell_center[i].comp = top_comp[i];
        }
}	/* end setComponent */

void CARTESIAN::setInitialCondition(void)
{
	int i;
	double coords[MAXD];
	INTERFACE *intfc = front->interf;
	POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
	STATE *sl,*sr;
	int c;

	FT_MakeGridIntfc(front);
        setDomain();

	/* Initialize states at the interface */
        next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
	    FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
            sl->temperature = sr->temperature =
                        0.5*(eqn_params->T0[0] + eqn_params->T0[1]);
        }

	// cell_center
	for (i = 0; i < cell_center.size(); i++)
	{
	    c = top_comp[i];
	    if (c == LIQUID_COMP)
	    	field->temperature[i] = eqn_params->T0[1];
	    else if (c == SOLID_COMP)
	    	field->temperature[i] = eqn_params->T0[0];
	    else
	    	field->temperature[i] = 0.0;
	}
}	/* end setInitialCondition */

void CARTESIAN::setIndexMap(COMPONENT sub_comp)
{
	static boolean first = YES;
	int i,j,k,ic,index;
	int llbuf[MAXD],uubuf[MAXD];
	int count;

	if (debugging("trace")) printf("Entering setIndexMap()\n");
	if (first)
	{
	    first = NO;
	    switch (dim)
	    {
	    case 1:
		count = (imax - imin + 1);
	    	FT_VectorMemoryAlloc((POINTER*)&i_to_I,top_gmax[0]+1,INT);
	    	FT_MatrixMemoryAlloc((POINTER*)&I_to_i,count,1,INT);
	    	break;
	    case 2:
		count = (imax - imin + 1)*(jmax - jmin + 1);
	    	FT_MatrixMemoryAlloc((POINTER*)&ij_to_I,top_gmax[0]+1,
					top_gmax[1]+1,INT);
	    	FT_MatrixMemoryAlloc((POINTER*)&I_to_ij,count,2,INT);
	    	break;
	    case 3:
		count = (imax - imin + 1)*(jmax - jmin + 1)*(kmax - kmin + 1);
	    	FT_TriArrayMemoryAlloc((POINTER*)&ijk_to_I,top_gmax[0]+1,
					top_gmax[1]+1,top_gmax[2]+1,INT);
	    	FT_MatrixMemoryAlloc((POINTER*)&I_to_ijk,count,3,INT);
	    	break;
	    }
	}

	index = 0;
	for (i = 0; i < dim; ++i)
	{
	    llbuf[i] = lbuf[i] != 0 ? lbuf[0] : 1;
	    uubuf[i] = ubuf[i] != 0 ? ubuf[0] : 1;
	}
	switch (dim)
	{
	case 1:
	    for (i = imin; i <= imax; i++)
	    {
		ic = d_index1d(i,top_gmax);
		if (sub_comp == NO_COMP || cell_center[ic].comp == sub_comp)
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
		if (sub_comp == NO_COMP || cell_center[ic].comp == sub_comp)
		{
	    	    ij_to_I[i][j] = index + ilower;
		    I_to_ij[index + ilower][0] = i;
                    I_to_ij[index + ilower][1] = j;
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
		if (sub_comp == NO_COMP || cell_center[ic].comp == sub_comp)
		{
	    	    ijk_to_I[i][j][k] = index + ilower;
		    I_to_ijk[index + ilower][0] = i;
                    I_to_ijk[index + ilower][1] = j;
                    I_to_ijk[index + ilower][2] = k;
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

void CARTESIAN::computeAdvection()
{
	int i;
	COMPONENT sub_comp[2];

	if (eqn_params->num_scheme == UNSPLIT_IMPLICIT_CIM)
	    return computeAdvectionCim();

	sub_comp[0] = SOLID_COMP;
	sub_comp[1] = LIQUID_COMP;

	for (i = 0; i < 2; ++i)
	{
	    setGlobalIndex(sub_comp[i]);
	    eqn_params->D = eqn_params->k[i]/eqn_params->rho[i]
				/eqn_params->Cp[i];
	    if (eqn_params->num_scheme == UNSPLIT_EXPLICIT)
	    	computeAdvectionExplicit(sub_comp[i]);
	    else if (eqn_params->num_scheme == UNSPLIT_IMPLICIT)
	    	computeAdvectionImplicit(sub_comp[i]);
	    else if (eqn_params->num_scheme == CRANK_NICOLSON)
	    	computeAdvectionCN(sub_comp[i]);
	}
}
    
void CARTESIAN::computeAdvectionCN(COMPONENT sub_comp)
{
	int i,j,k,l,m,ic,icn,I,I_nb,icoords[MAXD];
	int gmin[MAXD],ipn[MAXD];
	double crx_coords[MAXD];
	double T0,T_nb,D,lambda,coeff,coeff_nb,rhs;
	COMPONENT comp;
	PETSc solver;
	double *x;
        int num_iter = 0;
        double rel_residual = 0;
        boolean fr_crx_grid_seg;
        const GRID_DIRECTION dir[3][2] =
                {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
	double *Temp = field->temperature;

	start_clock("computeAdvectionCN");
	if (debugging("trace")) printf("Entering computeAdvectionCN()\n");

	D = eqn_params->D;
	
	for (i = 0; i < dim; ++i) gmin[i] = 0;

	setIndexMap(sub_comp);
	if (debugging("trace")) 
	{
	    int domain_size = 1;
	    printf("ilower = %d  iupper = %d\n",ilower,iupper);
	    for (i = 0; i < dim; ++i)
		domain_size *= (imax-imin+1);
	    printf("domain_size = %d\n",domain_size);
	}
        
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
	    	if (comp != sub_comp) 
	    	    continue;
		I = i_to_I[i];
                T0 = Temp[ic];
		rhs = T0;
		coeff = 1.0;
		for (l = 0; l < dim; ++l)
		{
		    lambda = D*m_dt/sqr(top_h[l]);
		    coeff += lambda;
                    for (m = 0; m < 2; ++m)
                    {
                        next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax);
                        icn = d_index1d(ipn[0],top_gmax);
			I_nb = i_to_I[ipn[0]];
			coeff_nb = -0.5*lambda;
                        fr_crx_grid_seg = FT_StateVarAtGridCrossing(front,
                                icoords,dir[l][m],comp,getStateTemperature,
                                &T_nb,crx_coords);
                        if (!fr_crx_grid_seg)
                        {
			    solver.Add_A(I,I_nb,coeff_nb);
			    T_nb = Temp[icn];
			    rhs -= coeff_nb*(T_nb - T0);
                        }
			else
			{
			    rhs -= coeff_nb*(2.0*T_nb - T0);
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
	    	if (comp != sub_comp) 
	    	    continue;
                T0 = Temp[ic];
		rhs = T0;
		coeff = 1.0;
		for (l = 0; l < dim; ++l)
		{
		    lambda = D*m_dt/sqr(top_h[l]);
		    coeff += lambda;
                    for (m = 0; m < 2; ++m)
                    {
                        next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax);
                        icn = d_index2d(ipn[0],ipn[1],top_gmax);
			I_nb = ij_to_I[ipn[0]][ipn[1]];
			coeff_nb = -0.5*lambda;
                        fr_crx_grid_seg = FT_StateVarAtGridCrossing(front,
                                icoords,dir[l][m],comp,getStateTemperature,
                                &T_nb,crx_coords);
                        if (!fr_crx_grid_seg)
                        {
			    solver.Add_A(I,I_nb,coeff_nb);
			    T_nb = Temp[icn];
			    rhs -= coeff_nb*(T_nb - T0);
                        }
			else
			{
			    rhs -= coeff_nb*(2.0*T_nb - T0);
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
	    	if (comp != sub_comp) 
	    	    continue;
                T0 = Temp[ic];
		rhs = T0;
		coeff = 1.0;
		for (l = 0; l < dim; ++l)
		{
		    lambda = D*m_dt/sqr(top_h[l]);
		    coeff += lambda;
                    for (m = 0; m < 2; ++m)
                    {
                        next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax);
                        icn = d_index3d(ipn[0],ipn[1],ipn[2],top_gmax);
			I_nb = ijk_to_I[ipn[0]][ipn[1]][ipn[2]];
			coeff_nb = -0.5*lambda;
                        fr_crx_grid_seg = FT_StateVarAtGridCrossing(front,
                                icoords,dir[l][m],comp,getStateTemperature,
                                &T_nb,crx_coords);
                        if (!fr_crx_grid_seg)
                        {
			    solver.Add_A(I,I_nb,coeff_nb);
			    T_nb = Temp[icn];
			    rhs -= coeff_nb*(T_nb - T0);
                        }
			else
			{
			    rhs -= coeff_nb*(2.0*T_nb - T0);
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
	    (void) printf("CARTESIAN::computeAdvectionCN: "
	       		"num_iter = %d, rel_residual = %g \n", 
			num_iter, rel_residual);
	}

	FT_VectorMemoryAlloc((POINTER*)&x,cell_center.size(),sizeof(double));
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
	    	if (comp == sub_comp)
	    	    array[ic] = x[I-ilower];
	    	else
	    	    array[ic] = 0.0;
	    }
	    break;
        case 2:
            for (i = imin; i <= imax; i++)
            for (j = jmin; j <= jmax; j++)
            {
	    	I = ij_to_I[i][j];
	    	ic = d_index2d(i,j,top_gmax);
	    	comp = cell_center[ic].comp;
	    	if (comp == sub_comp)
	    	    array[ic] = x[I-ilower];
	    	else
	    	    array[ic] = 0.0;
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
	    	if (comp == sub_comp)
	    	    array[ic] = x[I-ilower];
	    	else
	    	    array[ic] = 0.0;
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
	    	comp = cell_center[ic].comp;
	    	if (comp == sub_comp)
		    Temp[ic] = array[ic];
	    }
	    break;
        case 2:
            for (i = 0; i <= top_gmax[0]; ++i)
            for (j = 0; j <= top_gmax[1]; ++j)
            {
		ic = d_index2d(i,j,top_gmax);
	    	comp = cell_center[ic].comp;
	    	if (comp == sub_comp)
		    Temp[ic] = array[ic];
	    }
	    break;
        case 3:
            for (i = 0; i <= top_gmax[0]; ++i)
            for (j = 0; j <= top_gmax[1]; ++j)
            for (k = 0; k <= top_gmax[2]; ++k)
            {
		ic = d_index3d(i,j,k,top_gmax);
	    	comp = cell_center[ic].comp;
	    	if (comp == sub_comp)
		    Temp[ic] = array[ic];
	    }
	    break;
	}
	stop_clock("scatter_data");
	FT_FreeThese(1,x);

	if (debugging("trace")) printf("Leaving computeAdvectionCN()\n");
	stop_clock("computeAdvectionCN");
}	/* end computeAdvectionCN */

void CARTESIAN::computeAdvectionCim()
{
	//static CIM_PARAB_SOLVER parab_solver(*front);
	static CIM_PARAB_SOLVER parab_solver(*front);
	int i,j,k,index;

	printf("Entering computeAdvectionCim()\n");
	printf("m_dt = %f  min_dt = %f\n",m_dt,min_dt);
	/*
	if (m_dt < min_dt)
	    return computeAdvectionExplicit(sub_comp);
	*/

	if (debugging("trace")) printf("Entering computeAdvectionCim()\n");
	start_clock("computeAdvectionCim");
	setGlobalIndex(NO_COMP);
	setIndexMap(NO_COMP);

	switch(dim)
        {
        case 1:
            parab_solver.i_to_I = i_to_I;
            parab_solver.I_to_i = I_to_i;
	    for (i = imin; i <= imax; ++i)
	    {
		index = d_index1d(i,top_gmax);
		source[index] = -field->temperature[index]/m_dt;
	    }
            break;
        case 2:
            parab_solver.ij_to_I = ij_to_I;
            parab_solver.I_to_ij = I_to_ij;
	    for (i = imin; i <= imax; ++i)
	    for (j = jmin; j <= jmax; ++j)
	    {
		index = d_index2d(i,j,top_gmax);
		source[index] = -field->temperature[index]/m_dt;
	    }
            break;
        case 3:
            parab_solver.ijk_to_I = ijk_to_I;
            parab_solver.I_to_ijk = I_to_ijk;
	    for (i = imin; i <= imax; ++i)
	    for (j = jmin; j <= jmax; ++j)
	    for (k = kmin; k <= kmax; ++k)
	    {
		index = d_index3d(i,j,k,top_gmax);
		source[index] = -field->temperature[index]/m_dt;
	    }
            break;
        }
	printf("Assigning class variables\n");
	parab_solver.w_type = GROWING_BODY_BOUNDARY;
	parab_solver.neg_comp = SOLID_COMP;
        parab_solver.pos_comp = LIQUID_COMP;
        parab_solver.source = source;
	parab_solver.solutionJump = jumpT;
	parab_solver.gradJumpDotN = jumpEpsGradDotNorm;
        parab_solver.gradJumpDotT = jumpGradDotTan;
        parab_solver.diff_coeff[0] = eqn_params->k[0]/eqn_params->rho[0]
                                /eqn_params->Cp[0];
        parab_solver.diff_coeff[1] = eqn_params->k[1]/eqn_params->rho[1]
                                /eqn_params->Cp[1];
	parab_solver.jparams = (POINTER)eqn_params;
	parab_solver.findStateAtCrossing = find_state_at_crossing;
	parab_solver.getStateVar = getStateTemperature;
        parab_solver.assignStateVar = assignStateTemperature;
	parab_solver.ilower = ilower;
        parab_solver.iupper = iupper;
        parab_solver.soln = field->temperature;
        parab_solver.size = iupper - ilower;
	parab_solver.solve_front_state = YES;
        parab_solver.set_solver_domain();
        parab_solver.solve(field->temperature);

	stop_clock("computeAdvectionCim");
	if (debugging("trace")) printf("Leaving computeAdvectionCim()\n");
	clean_up(0);
	return;
}	/* end computeAdvectionImplicit */

void CARTESIAN::computeAdvectionImplicit(COMPONENT sub_comp)
{
	static PARABOLIC_SOLVER parab_solver(*front);
	static double *soln;
	int i,j,k,index;

	if (soln == NULL)
            FT_VectorMemoryAlloc((POINTER*)&soln,comp_size,FLOAT);

	if (debugging("trace")) printf("Entering computeAdvectionImplicit()\n");
	start_clock("computeAdvectionImplicit");
	setIndexMap(sub_comp);
	parab_solver.soln_comp = sub_comp;
        parab_solver.obst_comp = ERROR_COMP;
        parab_solver.var = field->temperature;
        parab_solver.soln = soln;
        parab_solver.getStateVarFunc = getStateTemperature;
        parab_solver.findStateAtCrossing = find_state_at_crossing;
        parab_solver.source = NULL;
        parab_solver.D = eqn_params->D;
        parab_solver.order = eqn_params->pde_order;
        parab_solver.ilower = ilower;
        parab_solver.iupper = iupper;
        parab_solver.dt = m_dt;
        parab_solver.set_solver_domain();

	if (sub_comp == LIQUID_COMP)
	    parab_solver.a = eqn_params->field->vel;
	else
	    parab_solver.a = NULL;

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
	switch(dim)
        {
        case 1:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    {
	    	index = d_index1d(i,top_gmax);
	    	if (top_comp[index] == sub_comp)
		    field->temperature[index] = soln[index];
	    }
            break;
        case 2:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    {
	    	index = d_index2d(i,j,top_gmax);
	    	if (top_comp[index] == sub_comp)
		    field->temperature[index] = soln[index];
	    }
            break;
        case 3:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    for (k = 0; k <= top_gmax[2]; ++k)
	    {
	    	index = d_index3d(i,j,k,top_gmax);
	    	if (top_comp[index] == sub_comp)
		    field->temperature[index] = soln[index];
	    }
            break;
        }
	stop_clock("computeAdvectionImplicit");
	if (debugging("trace")) printf("Leaving computeAdvectionImplicit()\n");
	return;
}	/* end computeAdvectionImplicit */


// for initial condition: 
// 		setInitialCondition();	
// this function should be called before solve()
// for the source term of the momentum equation: 	
// 		computeSourceTerm();
void CARTESIAN::solve(double dt)
{

	if (debugging("trace")) printf("Entering solve()\n");
	m_dt = dt;
	start_clock("solve");

	setDomain();
        if (debugging("trace")) printf("Passing setDomain()\n");

	if (debugging("sample_temperature"))
            sampleTemperature();

	setComponent();
	if (debugging("trace")) printf("Passing setComponent()\n");

	if (debugging("sample_temperature"))
            sampleTemperature();

	computeAdvection();
	if (debugging("trace")) printf("Passing liquid computeAdvection()\n");

	if (debugging("sample_temperature"))
            sampleTemperature();

	setAdvectionDt();
	if (debugging("trace")) printf("Passing setAdvectionDt()\n");

	setBoundary();
	if (debugging("trace")) printf("Passing setBoundary()\n");

	stop_clock("solve");
	if (debugging("trace")) printf("Leaving solve()\n");
}


void CARTESIAN::setAdvectionDt()
{
	double D,Dl,Ds;
	static double m_dt_expl,m_dt_impl;  // explicit and implicit time step
	static boolean first = YES;

	if (first)
	{
	    first = NO;
	    eqn_params = (PARAMS*)front->extra2;
	    Dl = eqn_params->k[1]/eqn_params->rho[1]/eqn_params->Cp[1];
	    Ds = eqn_params->k[0]/eqn_params->rho[0]/eqn_params->Cp[0];
	    D = std::max(Dl,Ds);
	    m_dt_expl = 0.5*sqr(hmin)/D/(double)dim;
	    m_dt_impl = 0.5*hmin/D/(double)dim;
	    min_dt = 0.1*sqr(hmin)/D/(double)dim;
	}

	if (eqn_params->num_scheme == UNSPLIT_EXPLICIT)
	    m_dt = m_dt_expl;
	else
	{
	    // For smooth transition to implicit step
	    double tstep = (double)front->step;
	    double smooth_factor; 
	    smooth_factor = 1.0/(1.0 + sqr(tstep/20.0));
	    m_dt = m_dt_impl - (m_dt_impl - m_dt_expl)*smooth_factor;
	}
	if (debugging("trace"))
	{
	    printf("In setAdvectionDt: m_dt = %24.18g min_dt = %f\n",
				m_dt,min_dt);
	}
}	/* end setAdvectionDt */

void CARTESIAN::getVelocity(double *p, double *U)
{
	// locate the point
	int icoords[MAXD];
	int i,j,k;
	double c1[MAXD], c2[MAXD], denominator;

	if (!rect_in_which(p,icoords,top_grid))
	{
	    for (i=0; i<2; i++)
	    {
	    	U[i] = 0.0;
	    }
	    return;
	}

	switch (dim)
	{
	case 2:
	    break;
	}
}

void CARTESIAN::getRectangleIndex(int index, int &i, int &j)
{
	i = cell_center[index].icoords[0];
	j = cell_center[index].icoords[1];
}

void CARTESIAN::getRectangleIndex(int index, int &i, int &j, int &k)
{
	i = cell_center[index].icoords[0];
	j = cell_center[index].icoords[1];
	k = cell_center[index].icoords[2];
}


int CARTESIAN::getRectangleComponent(int index)
{	
	return getComponent(cell_center[index].icoords);
}

void CARTESIAN::getRectangleCenter(
	int index, 
	double *coords)
{
	int i;
	for (i = 0; i < dim; ++i)
	    coords[i] = cell_center[index].coords[i];
}

void CARTESIAN::getRectangleCenter(
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

int CARTESIAN::getComponent(
	double *p)
{
	return component(p,front->interf);
}

int CARTESIAN::getComponent(
	int *icoords)
{
	int index;
	switch (dim)
	{
	case 1:
	    index = d_index1d(icoords[0],top_gmax);
	    return top_comp[index];
	case 2:
	    index = d_index2d(icoords[0],icoords[1],top_gmax);
	    return top_comp[index];
	case 3:
	    index = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
	    return top_comp[index];
	}
}

void CARTESIAN::save(char *filename)
{
	
	RECT_GRID *rect_grid = front->rect_grid;
	INTERFACE *intfc    = front->interf;
		
	int i, j;
	int xmax = rect_grid->gmax[0];
	int ymax = rect_grid->gmax[1];
	double x, y;
	
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

CARTESIAN::CARTESIAN(Front &front):front(&front)
{
}

void CARTESIAN::makeGridIntfc()
{
	static boolean first = YES;
	INTERFACE *grid_intfc;
	Table *T;
	static double *temperature;
	int i;

	FT_MakeGridIntfc(front);
	if (debugging("trace")) printf("Passed FT_MakeGridIntfc()\n");

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
	eqn_params = (PARAMS*)front->extra2;
	hmin = top_h[0];
	for (i = 1; i < dim; ++i)
	    if (hmin > top_h[i]) hmin = top_h[i];

	switch (dim)
	{
	case 1:
            if (first)
            {
                FT_VectorMemoryAlloc((POINTER*)&array,top_gmax[0]+1,FLOAT);
                FT_VectorMemoryAlloc((POINTER*)&temperature,top_gmax[0]+1,FLOAT);
                first = NO;
            }	
	    imin = (lbuf[0] == 0) ? 1 : lbuf[0];
	    imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
	    eqn_params->field->temperature = temperature;
	    eqn_params = (PARAMS*)front->extra2;
	    break;
	case 2:
	    if (first)
	    {
	    	FT_VectorMemoryAlloc((POINTER*)&array,(top_gmax[0]+1)*(top_gmax[1]+1),FLOAT);
	    	FT_VectorMemoryAlloc((POINTER*)&temperature,(top_gmax[0]+1)*(top_gmax[1]+1),FLOAT);
	    	first = NO;
	    }
	    imin = (lbuf[0] == 0) ? 1 : lbuf[0];
	    jmin = (lbuf[1] == 0) ? 1 : lbuf[1];
	    imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
	    jmax = (ubuf[1] == 0) ? top_gmax[1] - 1 : top_gmax[1] - ubuf[1];
	    eqn_params->field->temperature = temperature;
	    break;
	case 3:
	    if (first)
	    {
	    	FT_VectorMemoryAlloc((POINTER*)&array,
			(top_gmax[0]+1)*(top_gmax[1]+1)*(top_gmax[2]+1),FLOAT);
	    	FT_VectorMemoryAlloc((POINTER*)&temperature,
			(top_gmax[0]+1)*(top_gmax[1]+1)*(top_gmax[2]+1),FLOAT);
	    	first = NO;
	    }
	    imin = (lbuf[0] == 0) ? 1 : lbuf[0];
	    jmin = (lbuf[1] == 0) ? 1 : lbuf[1];
	    kmin = (lbuf[2] == 0) ? 1 : lbuf[2];
	    imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
	    jmax = (ubuf[1] == 0) ? top_gmax[1] - 1 : top_gmax[1] - ubuf[1];
	    kmax = (ubuf[2] == 0) ? top_gmax[2] - 1 : top_gmax[2] - ubuf[2];
	    eqn_params->field->temperature = temperature;
	    break;
	}
}	/* end makeGridIntfc */

void CARTESIAN::deleteGridIntfc()
{
	FT_FreeGridIntfc(front);
}

void CARTESIAN::scatMeshArray()
{
	FT_ParallelExchGridArrayBuffer(array,front);
}

void CARTESIAN::setGlobalIndex(COMPONENT sub_comp)
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
		if (sub_comp != NO_COMP && cell_center[ic].comp != sub_comp) 
		    continue;
		NLblocks++;
	    }
	    break;
	case 2:
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		ic = d_index2d(i,j,top_gmax);
		if (sub_comp != NO_COMP && cell_center[ic].comp != sub_comp) 
		    continue;
		NLblocks++;
	    }
	    break;
	case 3:
	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		ic = d_index3d(i,j,k,top_gmax);
		if (sub_comp != NO_COMP && cell_center[ic].comp != sub_comp) 
		    continue;
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

void CARTESIAN::printFrontInteriorState(char *out_name)
{
	int i,j,k,l,index;
	char filename[100];
	FILE *outfile;
	INTERFACE *intfc = front->interf;
	STATE *sl,*sr;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
	double *Temp = field->temperature;

	sprintf(filename,"%s/state.ts%s",out_name,right_flush(front->step,7));
#if defined(__MPI__)
	if (pp_numnodes() > 1)
            sprintf(filename,"%s-nd%s",filename,right_flush(pp_mynode(),4));
#endif /* defined(__MPI__) */
	outfile = fopen(filename,"w");
	
        /* Initialize states at the interface */
        fprintf(outfile,"Interface states:\n");
        int count = 0;
        next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
            count++;
            FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
            fprintf(outfile,"%24.18g %24.18g\n",sl->temperature,
				sr->temperature);
        }

	fprintf(outfile,"\nInterior states:\n");
	switch (dim)
	{
	case 1:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    {
		index = d_index1d(i,top_gmax);
	        fprintf(outfile,"%24.18g\n",Temp[index]);
	    }
	    break;
	case 2:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    {
		index = d_index2d(i,j,top_gmax);
	        fprintf(outfile,"%24.18g\n",Temp[index]);
	    }
	    break;
	case 3:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    for (k = 0; k <= top_gmax[2]; ++k)
	    {
		index = d_index3d(i,j,k,top_gmax);
	        fprintf(outfile,"%24.18g\n",Temp[index]);
	    }
	    break;
	}
	fclose(outfile);
}

void CARTESIAN::readFrontInteriorState(char *restart_name)
{
	FILE *infile;
	int i,j,k,l,index;
	INTERFACE *intfc = front->interf;
	STATE *sl,*sr;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
	double x;
	double *Temp = field->temperature;

	infile = fopen(restart_name,"r");

        /* Initialize states at the interface */
        next_output_line_containing_string(infile,"Interface states:");
        next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
            FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
            fscanf(infile,"%lf",&x);
            sl->temperature = x;
            fscanf(infile,"%lf",&x);
            sr->temperature = x;
        }

	FT_MakeGridIntfc(front);
        setDomain();

        /* Initialize states in the interior regions */

	next_output_line_containing_string(infile,"Interior states:");

	switch (dim)
	{
	case 1:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    {
		index = d_index1d(i,top_gmax);
	    	fscanf(infile,"%lf",&Temp[index]);
	    }
	    break;
	case 2:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    {
		index = d_index2d(i,j,top_gmax);
	    	fscanf(infile,"%lf",&Temp[index]);
	    }
	    break;
	case 3:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    for (k = 0; k <= top_gmax[2]; ++k)
	    {
		index = d_index3d(i,j,k,top_gmax);
	    	fscanf(infile,"%lf",&Temp[index]);
	    }
	    break;
	}
	fclose(infile);
}

void CARTESIAN::setBoundary()
{
	int i,j,k,index0,index1;
	INTERFACE *intfc = front->interf;
	double *Temp = field->temperature;

	switch (dim)
	{
	case 1:
	    if (rect_boundary_type(intfc,0,0) == NEUMANN_BOUNDARY)
	    {
		index0 = d_index1d(0,top_gmax);
		index1 = d_index1d(1,top_gmax);
	    	Temp[index0]  = Temp[index1];
	    }
	    if (rect_boundary_type(intfc,0,1) == NEUMANN_BOUNDARY)
	    {
		index0 = d_index1d(top_gmax[0],top_gmax);
		index1 = d_index1d(top_gmax[0]-1,top_gmax);
	    	Temp[index0]  = Temp[index1];
	    }
	    break;
	case 2:
	    if (rect_boundary_type(intfc,0,0) == NEUMANN_BOUNDARY)
	    {
		for (j = 0; j <= top_gmax[1]; ++j)
		{
		    index0 = d_index2d(0,j,top_gmax);
		    index1 = d_index2d(1,j,top_gmax);
	    	    Temp[index0]  = Temp[index1];
		}
	    }
	    if (rect_boundary_type(intfc,0,1) == NEUMANN_BOUNDARY)
	    {
		for (j = 0; j <= top_gmax[1]; ++j)
		{
		    index0 = d_index2d(top_gmax[0],j,top_gmax);
		    index1 = d_index2d(top_gmax[0]-1,j,top_gmax);
	    	    Temp[index0]  = Temp[index1];
		}
	    }
	    if (rect_boundary_type(intfc,1,0) == NEUMANN_BOUNDARY)
	    {
		for (i = 0; i <= top_gmax[0]; ++i)
		{
		    index0 = d_index2d(i,0,top_gmax);
		    index1 = d_index2d(i,1,top_gmax);
	    	    Temp[index0]  = Temp[index1];
		}
	    }
	    if (rect_boundary_type(intfc,1,1) == NEUMANN_BOUNDARY)
	    {
		for (i = 0; i <= top_gmax[0]; ++i)
		{
		    index0 = d_index2d(i,top_gmax[1],top_gmax);
		    index1 = d_index2d(i,top_gmax[1]-1,top_gmax);
	    	    Temp[index0]  = Temp[index1];
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
	    	    Temp[index0]  = Temp[index1];
		}
	    }
	    if (rect_boundary_type(intfc,0,1) == NEUMANN_BOUNDARY)
	    {
		for (j = 0; j <= top_gmax[1]; ++j)
		for (k = 0; k <= top_gmax[2]; ++k)
		{
		    index0 = d_index3d(top_gmax[0],j,k,top_gmax);
		    index1 = d_index3d(top_gmax[0]-1,j,k,top_gmax);
	    	    Temp[index0]  = Temp[index1];
		}
	    }
	    if (rect_boundary_type(intfc,1,0) == NEUMANN_BOUNDARY)
	    {
		for (i = 0; i <= top_gmax[0]; ++i)
		for (k = 0; k <= top_gmax[2]; ++k)
		{
		    index0 = d_index3d(i,0,k,top_gmax);
		    index1 = d_index3d(i,1,k,top_gmax);
	    	    Temp[index0]  = Temp[index1];
		}
	    }
	    if (rect_boundary_type(intfc,1,1) == NEUMANN_BOUNDARY)
	    {
		for (i = 0; i <= top_gmax[0]; ++i)
		for (k = 0; k <= top_gmax[2]; ++k)
		{
		    index0 = d_index3d(i,top_gmax[1],k,top_gmax);
		    index1 = d_index3d(i,top_gmax[1]-1,k,top_gmax);
	    	    Temp[index0]  = Temp[index1];
		}
	    }
	    if (rect_boundary_type(intfc,2,0) == NEUMANN_BOUNDARY)
	    {
		for (i = 0; i <= top_gmax[0]; ++i)
		for (j = 0; j <= top_gmax[1]; ++j)
		{
		    index0 = d_index3d(i,j,0,top_gmax);
		    index1 = d_index3d(i,j,1,top_gmax);
	    	    Temp[index0]  = Temp[index1];
		}
	    }
	    if (rect_boundary_type(intfc,2,1) == NEUMANN_BOUNDARY)
	    {
		for (i = 0; i <= top_gmax[0]; ++i)
		for (j = 0; j <= top_gmax[1]; ++j)
		{
		    index0 = d_index3d(i,j,top_gmax[2],top_gmax);
		    index1 = d_index3d(i,j,top_gmax[2]-1,top_gmax);
	    	    Temp[index0]  = Temp[index1];
		}
	    }
	    break;
	}
}	/* end setBoundary */

void CARTESIAN::oneDimPlot(char *outname)
{
	xgraphOneDimPlot(outname);
}	/* end temperaturePlot */

void CARTESIAN::xgraphOneDimPlot(char *outname)
{
	int i,index;
	char filename[100];
	FILE *outfile;
	double *Temp = field->temperature;

	if (debugging("trace"))
	    printf("Entering xgraphTemp1()\n");
        sprintf(filename,"%s/tmp-xg.ts%s",outname,right_flush(front->step,7));
#if defined(__MPI__)
        sprintf(filename,"%s-nd%s",filename,right_flush(pp_mynode(),4));
#endif /* defined(__MPI__) */
        outfile = fopen(filename,"w");
	fprintf(outfile,"\"Solid temp at %6.3f\"\n",front->time);
	for (i = 0; i <= top_gmax[0]; ++i)
	{
	    index = d_index1d(i,top_gmax);
	    if (cell_center[index].comp == SOLID_COMP)
	    	fprintf(outfile,"%24.18g  %24.18g\n",
	    		cell_center[index].coords[0],
	    		Temp[index]);
	}
	fprintf(outfile,"\n\n");
	fprintf(outfile,"\"Liquid temp at %6.3f\"\n",front->time);
	for (i = 0; i <= top_gmax[0]; ++i)
	{
	    index = d_index1d(i,top_gmax);
	    if (cell_center[index].comp == LIQUID_COMP)
	    	fprintf(outfile,"%24.18g  %24.18g\n",
	    		cell_center[index].coords[0],
	    		Temp[index]);
	}
	fclose(outfile);
	if (debugging("trace"))
	    printf("Leaving xgraphTemp1()\n");
}	/* end xgraphOneDimPlot */

void CARTESIAN::vtk_plot_temperature2d(
	char *outname)
{
	std::vector<int>ph_index;
        int i,j,k,index;
        char dirname[256],filename[256];
        FILE *outfile;
        double coord_x,coord_y,coord_z,xmin,ymin;
        COMPONENT comp;
        int pointsx,pointsy,num_points,num_cells,num_cell_list;
        int icoords[2],p_gmax[2];

        sprintf(filename, "%s/vtk.ts%s",outname,
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
        sprintf(filename,"%s/liquid.vtk",filename);
        outfile = fopen(filename,"w");
        fprintf(outfile,"# vtk DataFile Version 3.0\n");
        fprintf(outfile,"liquid temperature\n");
        fprintf(outfile,"ASCII\n");
        fprintf(outfile,"DATASET UNSTRUCTURED_GRID\n");

        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index2d(i,j,top_gmax);
            if (cell_center[index].comp == LIQUID_COMP)
                ph_index.push_back(index);
        }

        pointsx = top_gmax[0] + 2;
        pointsy = top_gmax[1] + 2;

        num_points = pointsx*pointsy;

        num_cells = (int) ph_index.size();
        num_cell_list = 5*num_cells;


        index = d_index2d(imin,jmin,top_gmax);
        xmin = cell_center[index].coords[0] - top_h[0]/2.0;
        ymin = cell_center[index].coords[1] - top_h[1]/2.0;

        fprintf(outfile,"POINTS %d double\n", num_points);
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index2d(i,j,top_gmax);
            coord_x = cell_center[index].coords[0] - top_h[0]/2.0;
            coord_y = cell_center[index].coords[1] - top_h[1]/2.0;
            coord_z = 0.0;
            fprintf(outfile,"%f %f %f\n",coord_x,coord_y,coord_z);
        }

	for (i = 0; i <= top_gmax[0]; i++)
        {
            j = top_gmax[1];
            index = d_index2d(i,j,top_gmax);
            coord_x = cell_center[index].coords[0] - top_h[0]/2.0;
            coord_y = cell_center[index].coords[1] + top_h[1]/2.0;
            coord_z = 0.0;
            fprintf(outfile,"%f %f %f\n",coord_x,coord_y,coord_z);
        }

        for (j = 0; j <= top_gmax[1]; j++)
        {
            i = top_gmax[0];
            index = d_index2d(i,j,top_gmax);
            coord_x = cell_center[index].coords[0] + top_h[0]/2.0;
            coord_y = cell_center[index].coords[1] - top_h[1]/2.0;
            coord_z = 0.0;
            fprintf(outfile,"%f %f %f\n",coord_x,coord_y,coord_z);
        }

        i = top_gmax[0];
        j = top_gmax[1];
        index = d_index2d(i,j,top_gmax);
        coord_x = cell_center[index].coords[0] + top_h[0]/2.0;
        coord_y = cell_center[index].coords[1] + top_h[1]/2.0;
        coord_z = 0.0;
        fprintf(outfile,"%f %f %f\n",coord_x,coord_y,coord_z);

	fprintf(outfile,"CELLS %i %i\n", num_cells,num_cell_list);
        for (i = 0; i < num_cells; i++)
        {
            int index0,index1,index2,index3;
            index = ph_index[i];
            icoords[0] = cell_center[index].icoords[0];
            icoords[1] = cell_center[index].icoords[1];

            index0 = d_index2d(icoords[0],icoords[1],top_gmax);
            index1 = d_index2d(icoords[0]+1,icoords[1],top_gmax);
            index2 = d_index2d(icoords[0],icoords[1]+1,top_gmax);
            index3 = d_index2d(icoords[0]+1,icoords[1]+1,top_gmax);

            fprintf(outfile,"4 %i %i %i %i\n",
                    index0,index1,index2,index3);
        }

        fprintf(outfile, "CELL_TYPES %i\n", num_cells);
        for (i = 0; i < num_cells; i++)
            fprintf(outfile, "8\n");

        fprintf(outfile, "CELL_DATA %i\n", num_cells);
        fprintf(outfile, "SCALARS temperature double\n");
        fprintf(outfile, "LOOKUP_TABLE default\n");
        for (i = 0; i < num_cells; i++)
        {
            index = ph_index[i];
            fprintf(outfile,"%f\n", eqn_params->field->temperature[index]);
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
        fprintf(outfile,"solid temperature\n");
        fprintf(outfile,"ASCII\n");
        fprintf(outfile,"DATASET UNSTRUCTURED_GRID\n");

        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index2d(i,j,top_gmax);
            if (cell_center[index].comp == SOLID_COMP)
                ph_index.push_back(index);
        }

        pointsx = top_gmax[0] + 2;
        pointsy = top_gmax[1] + 2;
        num_points = pointsx*pointsy;

        num_cells = (int) ph_index.size();
        num_cell_list = 5*num_cells;

        index = d_index2d(imin,jmin,top_gmax);
        xmin = cell_center[index].coords[0] - top_h[0]/2.0;
        ymin = cell_center[index].coords[1] - top_h[1]/2.0;

        fprintf(outfile,"POINTS %d double\n", num_points);
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index2d(i,j,top_gmax);
            coord_x = cell_center[index].coords[0] - top_h[0]/2.0;
            coord_y = cell_center[index].coords[1] - top_h[1]/2.0;
            coord_z = 0.0;
            fprintf(outfile,"%f %f %f\n",coord_x,coord_y,coord_z);
        }

	for (i = 0; i <= top_gmax[0]; i++)
        {
            j = top_gmax[1];
            index = d_index2d(i,j,top_gmax);
            coord_x = cell_center[index].coords[0] - top_h[0]/2.0;
            coord_y = cell_center[index].coords[1] + top_h[1]/2.0;
            coord_z = 0.0;
            fprintf(outfile,"%f %f %f\n",coord_x,coord_y,coord_z);
        }

        for (j = 0; j <= top_gmax[1]; j++)
        {
            i = top_gmax[0];
            index = d_index2d(i,j,top_gmax);
            coord_x = cell_center[index].coords[0] + top_h[0]/2.0;
            coord_y = cell_center[index].coords[1] - top_h[1]/2.0;
            coord_z = 0.0;
            fprintf(outfile,"%f %f %f\n",coord_x,coord_y,coord_z);
        }

        i = top_gmax[0];
        j = top_gmax[1];
        index = d_index2d(i,j,top_gmax);
        coord_x = cell_center[index].coords[0] + top_h[0]/2.0;
        coord_y = cell_center[index].coords[1] + top_h[1]/2.0;
        coord_z = 0.0;
        fprintf(outfile,"%f %f %f\n",coord_x,coord_y,coord_z);



        fprintf(outfile,"CELLS %i %i\n", num_cells,num_cell_list);
        for (i = 0; i < num_cells; i++)
        {
            int index0,index1,index2,index3;
            index = ph_index[i];
            icoords[0] = cell_center[index].icoords[0];
            icoords[1] = cell_center[index].icoords[1];

            index0 = d_index2d(icoords[0],icoords[1],top_gmax);
            index1 = d_index2d(icoords[0]+1,icoords[1],top_gmax);
            index2 = d_index2d(icoords[0],icoords[1]+1,top_gmax);
            index3 = d_index2d(icoords[0]+1,icoords[1]+1,top_gmax);

            fprintf(outfile,"4 %i %i %i %i\n",
                    index0,index1,index2,index3);
        }

	fprintf(outfile, "CELL_TYPES %i\n", num_cells);
        for (i = 0; i < num_cells; i++)
            fprintf(outfile, "8\n");

        fprintf(outfile, "CELL_DATA %i\n", num_cells);
        fprintf(outfile, "SCALARS temperature double\n");
        fprintf(outfile, "LOOKUP_TABLE default\n");
        for (i = 0; i < num_cells; i++)
        {
            index = ph_index[i];
            fprintf(outfile,"%f\n",eqn_params->field->temperature[index]);
        }

        fclose(outfile);
}       /* end vtk_plot_temperature2d */

void CARTESIAN::vtk_plot_temperature3d(
        char *outname)
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

        sprintf(filename, "%s/vtk.ts%s",outname,
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
        sprintf(filename,"%s/liquid.vtk",filename);
        outfile = fopen(filename,"w");
        fprintf(outfile,"# vtk DataFile Version 3.0\n");
        fprintf(outfile,"liquid temperature\n");
        fprintf(outfile,"ASCII\n");
        fprintf(outfile,"DATASET UNSTRUCTURED_GRID\n");

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            if (cell_center[index].comp == LIQUID_COMP)
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
            coord_x = cell_center[index].coords[0] - top_h[0]/2.0;
            coord_y = cell_center[index].coords[1] - top_h[1]/2.0;
            coord_z = cell_center[index].coords[2] - top_h[2]/2.0;
            fprintf(outfile,"%f %f %f\n",coord_x,coord_y,coord_z);
        }

	for (i = 0; i <= top_gmax[0]; i++)
        for (j = 0; j <= top_gmax[1]; j++)
        {
            k = top_gmax[2];
            index = d_index3d(i,j,k,top_gmax);
            coord_x = cell_center[index].coords[0] - top_h[0]/2.0;
            coord_y = cell_center[index].coords[1] - top_h[1]/2.0;
            coord_z = cell_center[index].coords[2] + top_h[2]/2.0;
            fprintf(outfile,"%f %f %f\n",coord_x,coord_y,coord_z);
        }

        for (i = 0; i <= top_gmax[0]; i++)
        for (k = 0; k <= top_gmax[2]; k++)
        {
            j = top_gmax[1];
            index = d_index3d(i,j,k,top_gmax);
            coord_x = cell_center[index].coords[0] - top_h[0]/2.0;
            coord_y = cell_center[index].coords[1] + top_h[1]/2.0;
            coord_z = cell_center[index].coords[2] - top_h[2]/2.0;
            fprintf(outfile,"%f %f %f\n",coord_x,coord_y,coord_z);
        }

        for (j = 0; j <= top_gmax[1]; j++)
        for (k = 0; k <= top_gmax[2]; k++)
        {
            i = top_gmax[0];
            index = d_index3d(i,j,k,top_gmax);
            coord_x = cell_center[index].coords[0] + top_h[0]/2.0;
            coord_y = cell_center[index].coords[1] - top_h[1]/2.0;
            coord_z = cell_center[index].coords[2] - top_h[2]/2.0;
            fprintf(outfile,"%f %f %f\n",coord_x,coord_y,coord_z);
        }

	for (i = 0; i <= top_gmax[0]; i++)
        {
            j = top_gmax[1];
            k = top_gmax[2];
            index = d_index3d(i,j,k,top_gmax);
            coord_x = cell_center[index].coords[0] - top_h[0]/2.0;
            coord_y = cell_center[index].coords[1] + top_h[1]/2.0;
            coord_z = cell_center[index].coords[2] + top_h[2]/2.0;
            fprintf(outfile,"%f %f %f\n",coord_x,coord_y,coord_z);
        }

        for (j = 0; j <= top_gmax[1]; j++)
        {
            i = top_gmax[0];
            k = top_gmax[2];
            index = d_index3d(i,j,k,top_gmax);
            coord_x = cell_center[index].coords[0] + top_h[0]/2.0;
            coord_y = cell_center[index].coords[1] - top_h[1]/2.0;
            coord_z = cell_center[index].coords[2] + top_h[2]/2.0;
            fprintf(outfile,"%f %f %f\n",coord_x,coord_y,coord_z);
        }

        for (k = 0; k <= top_gmax[2]; k++)
        {
            i = top_gmax[0];
            j = top_gmax[1];
            index = d_index3d(i,j,k,top_gmax);
            coord_x = cell_center[index].coords[0] + top_h[0]/2.0;
            coord_y = cell_center[index].coords[1] + top_h[1]/2.0;
            coord_z = cell_center[index].coords[2] - top_h[2]/2.0;
            fprintf(outfile,"%f %f %f\n",coord_x,coord_y,coord_z);
        }

	i = top_gmax[0];
        j = top_gmax[1];
        k = top_gmax[2];
        index = d_index3d(i,j,k,top_gmax);
        coord_x = cell_center[index].coords[0] + top_h[0]/2.0;
        coord_y = cell_center[index].coords[1] + top_h[1]/2.0;
        coord_z = cell_center[index].coords[2] + top_h[2]/2.0;
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
        fprintf(outfile, "SCALARS temperature double\n");
        fprintf(outfile, "LOOKUP_TABLE default\n");
        for (i = 0; i < num_cells; i++)
        {
            index = ph_index[i];
            fprintf(outfile,"%f\n",eqn_params->field->temperature[index]);
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
        fprintf(outfile,"solid temperature\n");
        fprintf(outfile,"ASCII\n");
        fprintf(outfile,"DATASET UNSTRUCTURED_GRID\n");

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            if (cell_center[index].comp == SOLID_COMP)
                ph_index.push_back(index);
        }


        pointsx = top_gmax[0] + 2;
        pointsy = top_gmax[1] + 2;
        pointsz = top_gmax[2] + 2;
        num_points = pointsx*pointsy*pointsz;

        num_cells = (int)ph_index.size();
        num_cell_list = 9*num_cells;

        p_gmax[0] = pointsx - 1;
        p_gmax[1] = pointsy - 1;
        p_gmax[2] = pointsz - 1;

        fprintf(outfile,"POINTS %d double\n", num_points);

	for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            coord_x = cell_center[index].coords[0] - top_h[0]/2.0;
            coord_y = cell_center[index].coords[1] - top_h[1]/2.0;
            coord_z = cell_center[index].coords[2] - top_h[2]/2.0;
            fprintf(outfile,"%f %f %f\n",coord_x,coord_y,coord_z);
        }

        for (i = 0; i <= top_gmax[0]; i++)
        for (j = 0; j <= top_gmax[1]; j++)
        {
            k = top_gmax[2];
            index = d_index3d(i,j,k,top_gmax);
            coord_x = cell_center[index].coords[0] - top_h[0]/2.0;
            coord_y = cell_center[index].coords[1] - top_h[1]/2.0;
            coord_z = cell_center[index].coords[2] + top_h[2]/2.0;
            fprintf(outfile,"%f %f %f\n",coord_x,coord_y,coord_z);
        }

        for (i = 0; i <= top_gmax[0]; i++)
        for (k = 0; k <= top_gmax[2]; k++)
        {
            j = top_gmax[1];
            index = d_index3d(i,j,k,top_gmax);
            coord_x = cell_center[index].coords[0] - top_h[0]/2.0;
            coord_y = cell_center[index].coords[1] + top_h[1]/2.0;
            coord_z = cell_center[index].coords[2] - top_h[2]/2.0;
            fprintf(outfile,"%f %f %f\n",coord_x,coord_y,coord_z);
        }

        for (j = 0; j <= top_gmax[1]; j++)
        for (k = 0; k <= top_gmax[2]; k++)
        {
            i = top_gmax[0];
            index = d_index3d(i,j,k,top_gmax);
            coord_x = cell_center[index].coords[0] + top_h[0]/2.0;
            coord_y = cell_center[index].coords[1] - top_h[1]/2.0;
            coord_z = cell_center[index].coords[2] - top_h[2]/2.0;
            fprintf(outfile,"%f %f %f\n",coord_x,coord_y,coord_z);
        }

	for (i = 0; i <= top_gmax[0]; i++)
        {
            j = top_gmax[1];
            k = top_gmax[2];
            index = d_index3d(i,j,k,top_gmax);
            coord_x = cell_center[index].coords[0] - top_h[0]/2.0;
            coord_y = cell_center[index].coords[1] + top_h[1]/2.0;
            coord_z = cell_center[index].coords[2] + top_h[2]/2.0;
            fprintf(outfile,"%f %f %f\n",coord_x,coord_y,coord_z);
        }

        for (j = 0; j <= top_gmax[1]; j++)
        {
            i = top_gmax[0];
            k = top_gmax[2];
            index = d_index3d(i,j,k,top_gmax);
            coord_x = cell_center[index].coords[0] + top_h[0]/2.0;
            coord_y = cell_center[index].coords[1] - top_h[1]/2.0;
            coord_z = cell_center[index].coords[2] + top_h[2]/2.0;
            fprintf(outfile,"%f %f %f\n",coord_x,coord_y,coord_z);
        }

        for (k = 0; k <= top_gmax[2]; k++)
        {
            i = top_gmax[0];
            j = top_gmax[1];
            index = d_index3d(i,j,k,top_gmax);
            coord_x = cell_center[index].coords[0] + top_h[0]/2.0;
            coord_y = cell_center[index].coords[1] + top_h[1]/2.0;
            coord_z = cell_center[index].coords[2] - top_h[2]/2.0;
            fprintf(outfile,"%f %f %f\n",coord_x,coord_y,coord_z);
        }

	i = top_gmax[0];
        j = top_gmax[1];
        k = top_gmax[2];
        index = d_index3d(i,j,k,top_gmax);
        coord_x = cell_center[index].coords[0] + top_h[0]/2.0;
        coord_y = cell_center[index].coords[1] + top_h[1]/2.0;
        coord_z = cell_center[index].coords[2] + top_h[2]/2.0;
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
        fprintf(outfile, "SCALARS temperature double\n");
        fprintf(outfile, "LOOKUP_TABLE default\n");
        for (i = 0; i < num_cells; i++)
        {
            index = ph_index[i];
            fprintf(outfile,"%f\n",eqn_params->field->temperature[index]);
        }
        fclose(outfile);
}       /* end vtk_plot_temperature3d */

void CARTESIAN::computeAdvectionExplicit(COMPONENT sub_comp)
{
	int i,j,k,l,m,ic,icn,icoords[MAXD],nc;
	int gmin[MAXD],ipn[MAXD];
	int index0;
	double coords[MAXD],crx_coords[MAXD];
	double temperature,temperature_nb[2],dgrad[MAXD],grad_plus[MAXD],
	       grad_minus[MAXD];
	double coef;
	COMPONENT comp;
	boolean fr_crx_grid_seg;
	const GRID_DIRECTION dir[3][2] = 
		{{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
	double *Temp = field->temperature;
	double v[MAXD],**vel,v_plus[MAXD],v_minus[MAXD];

	start_clock("computeAdvectionExplicit");

	coef = eqn_params->D*m_dt;
	vel = eqn_params->field->vel;

	for (i = 0; i < dim; ++i) gmin[i] = 0;

	switch (dim)
	{
        case 1:
            for (i = imin; i <= imax; ++i)
            {
                icoords[0] = i;
                ic = d_index1d(i,top_gmax);
                comp = top_comp[ic];
                if (comp != sub_comp)
                     continue;
                array[ic] = temperature = Temp[ic];
		for (l = 0; l < dim; ++l)
		{
		    v[l] = 0.0;
		    v_plus[l] = 0.0;
		    v_minus[l] = 0.0;
		}
		if (sub_comp == LIQUID_COMP && vel != NULL)
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
                                icoords,dir[l][m],comp,getStateTemperature,
                                &temperature_nb[m],crx_coords);
                        if (!fr_crx_grid_seg)
                        {
                             next_ip_in_dir(icoords,dir[l][m],ipn,gmin,
							top_gmax);
                             icn = d_index1d(ipn[0],top_gmax);
			     temperature_nb[m] = Temp[icn];
                        }
                        dgrad[l] += (temperature_nb[m] - temperature)/top_h[l];
                    }
		    grad_plus[l] = (temperature_nb[1] - temperature)/top_h[l];
		    grad_minus[l] = (temperature - temperature_nb[0])/top_h[l];
                    array[ic] += coef*dgrad[l]/top_h[l]-m_dt*(v_plus[l]*grad_minus[l]+
					v_minus[l]*grad_plus[l]);
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
	    	if (comp != sub_comp) 
	    	    continue;
                array[ic] = temperature = Temp[ic];
		for (l = 0; l < dim; ++l)
                {
                    v[l] = 0.0;
                    v_plus[l] = 0.0;
                    v_minus[l] = 0.0;
                }
                if (sub_comp == LIQUID_COMP && vel != NULL)
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
                                icoords,dir[l][m],comp,getStateTemperature,
                                &temperature_nb[m],crx_coords);
                        if (!fr_crx_grid_seg)
                        {
                            next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax);
                            icn = d_index2d(ipn[0],ipn[1],top_gmax);
			    temperature_nb[m] = Temp[icn];
                        }
                        dgrad[l] += (temperature_nb[m] - temperature)/top_h[l];
                    }
		    grad_plus[l] = (temperature_nb[1] - temperature)/top_h[l];
		    grad_minus[l] = (temperature - temperature_nb[0])/top_h[l];
		    array[ic] += coef*dgrad[l]/top_h[l] - m_dt*(v_plus[l]*grad_minus[l]+
					v_minus[l]*grad_plus[l]);
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
	    	if (comp != sub_comp) 
	    	    continue;
                array[ic] = temperature = Temp[ic];
		for (l = 0; l < dim; ++l)
                {
                    v[l] = 0.0;
                    v_plus[l] = 0.0;
                    v_minus[l] = 0.0;
                }
                if (sub_comp == LIQUID_COMP && vel != NULL)
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
                                icoords,dir[l][m],comp,getStateTemperature,
                                &temperature_nb[m],crx_coords);
                        if (!fr_crx_grid_seg)
                        {
                            next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax);
                            icn = d_index3d(ipn[0],ipn[1],ipn[2],top_gmax);
                	    temperature_nb[m] = Temp[icn];
                        }
                        dgrad[l] += (temperature_nb[m] - temperature)/top_h[l];
                    }
		    grad_plus[l] = (temperature_nb[1] - temperature)/top_h[l];
		    grad_minus[l] = (temperature - temperature_nb[0])/top_h[l];
		    array[ic] += coef*dgrad[l]/top_h[l] - m_dt*(v_plus[l]*grad_minus[l]+
					v_minus[l]*grad_plus[l]);
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
	    	comp = top_comp[ic];
		if (comp == sub_comp)
		    Temp[ic] = array[ic];
	    }
	    break;
        case 2:
            for (i = 0; i <= top_gmax[0]; ++i)
            for (j = 0; j <= top_gmax[1]; ++j)
            {
		ic = d_index2d(i,j,top_gmax);
	    	comp = top_comp[ic];
		if (comp == sub_comp)
		    Temp[ic] = array[ic];
	    }
	    break;
        case 3:
            for (i = 0; i <= top_gmax[0]; ++i)
            for (j = 0; j <= top_gmax[1]; ++j)
            for (k = 0; k <= top_gmax[2]; ++k)
            {
		ic = d_index3d(i,j,k,top_gmax);
	    	comp = top_comp[ic];
		if (comp == sub_comp)
		    Temp[ic] = array[ic];
	    }
	    break;
	}

	stop_clock("computeAdvectionExplicit");
}	/* computeAdvectionExplicit */

void CARTESIAN::setDomain()
{
	static boolean first = YES;
	INTERFACE *grid_intfc;
	Table *T;
	int i;

	if (debugging("trace")) printf("Passed FT_MakeGridIntfc()\n");

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
	eqn_params = (PARAMS*)front->extra2;
	hmin = top_h[0];
	for (i = 1; i < dim; ++i)
	    if (hmin > top_h[i]) hmin = top_h[i];

	switch (dim)
	{
	case 1:
            if (first)
            {
		comp_size = top_gmax[0]+1;
		FT_ScalarMemoryAlloc((POINTER*)&field,sizeof(PHASE_FIELD));
                FT_VectorMemoryAlloc((POINTER*)&array,comp_size,FLOAT);
                FT_VectorMemoryAlloc((POINTER*)&source,comp_size,FLOAT);
                FT_VectorMemoryAlloc((POINTER*)&field->temperature,comp_size,
			FLOAT);
                first = NO;
            }	
	    imin = (lbuf[0] == 0) ? 1 : lbuf[0];
	    imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
	    eqn_params->field = field;
	    break;
	case 2:
	    if (first)
	    {
		comp_size = (top_gmax[0]+1)*(top_gmax[1]+1);
		FT_ScalarMemoryAlloc((POINTER*)&field,sizeof(PHASE_FIELD));
	    	FT_VectorMemoryAlloc((POINTER*)&array,comp_size,FLOAT);
                FT_VectorMemoryAlloc((POINTER*)&source,comp_size,FLOAT);
	    	FT_VectorMemoryAlloc((POINTER*)&field->temperature,
			comp_size,FLOAT);
	    	first = NO;
	    }
	    imin = (lbuf[0] == 0) ? 1 : lbuf[0];
	    jmin = (lbuf[1] == 0) ? 1 : lbuf[1];
	    imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
	    jmax = (ubuf[1] == 0) ? top_gmax[1] - 1 : top_gmax[1] - ubuf[1];
	    eqn_params->field = field;
	    break;
	case 3:
	    if (first)
	    {
		comp_size = (top_gmax[0]+1)*(top_gmax[1]+1)*(top_gmax[2]+1);
		FT_ScalarMemoryAlloc((POINTER*)&field,sizeof(PHASE_FIELD));
	    	FT_VectorMemoryAlloc((POINTER*)&array,comp_size,FLOAT);
                FT_VectorMemoryAlloc((POINTER*)&source,comp_size,FLOAT);
	    	FT_VectorMemoryAlloc((POINTER*)&field->temperature,
			comp_size,FLOAT);
	    	first = NO;
	    }
	    imin = (lbuf[0] == 0) ? 1 : lbuf[0];
	    jmin = (lbuf[1] == 0) ? 1 : lbuf[1];
	    kmin = (lbuf[2] == 0) ? 1 : lbuf[2];
	    imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
	    jmax = (ubuf[1] == 0) ? top_gmax[1] - 1 : top_gmax[1] - ubuf[1];
	    kmax = (ubuf[2] == 0) ? top_gmax[2] - 1 : top_gmax[2] - ubuf[2];
	    eqn_params->field = field;
	    break;
	}
}	/* end setDomain */

void CARTESIAN::initMovieVariables()
{
	int n;
	static HDF_MOVIE_VAR *hdf_movie_var;
	MOVIE_OPTION *movie_option = eqn_params->movie_option;

        if (debugging("trace"))
            printf("Entering initMovieVariables()\n");
	if (hdf_movie_var == NULL)
	{
	    FT_ScalarMemoryAlloc((POINTER*)&hdf_movie_var,
				sizeof(HDF_MOVIE_VAR));
	    switch (dim)
	    {
	    case 1:
		hdf_movie_var->num_var = 1;
                FT_MatrixMemoryAlloc((POINTER*)&hdf_movie_var->var_name,1,
                                        100,sizeof(char));
                FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->top_var,1,
                                        sizeof(double*));
                FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->obstacle_comp,1,
                                        sizeof(COMPONENT));
                FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->preset_bound,1,
                                        sizeof(boolean));
                FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->var_min,1,
                                        sizeof(double));
                FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->var_max,1,
                                        sizeof(double));
		if (movie_option->plot_temperature)
                {
                    sprintf(hdf_movie_var->var_name[0],"temperature");
                    hdf_movie_var->get_state_var[0] = getStateTemperature;
                    hdf_movie_var->top_var[0] = 
				eqn_params->field->temperature;
                    hdf_movie_var->obstacle_comp[0] = ERROR_COMP;
                }
		break;
	    case 2:
	    	hdf_movie_var->num_var = 1;
	    	FT_MatrixMemoryAlloc((POINTER*)&hdf_movie_var->var_name,1,
					100,sizeof(char));
	    	FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->top_var,1,
					sizeof(double*));
	    	FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->obstacle_comp,1,
					sizeof(COMPONENT));
		FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->preset_bound,1,
                                        sizeof(boolean));
                FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->var_min,1,
                                        sizeof(double));
                FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->var_max,1,
                                        sizeof(double));
		if (movie_option->plot_temperature)
		{
	    	    sprintf(hdf_movie_var->var_name[0],"temperature");
	    	    hdf_movie_var->get_state_var[0] = getStateTemperature;
		    hdf_movie_var->top_var[0] = 
				eqn_params->field->temperature;
		    hdf_movie_var->obstacle_comp[0] = ERROR_COMP;
		}
		break;
	    case 3:
	    	hdf_movie_var->num_var = n = 0;
	    	FT_MatrixMemoryAlloc((POINTER*)&hdf_movie_var->var_name,3,100,
					sizeof(char));
	    	FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->top_var,3,
					sizeof(double*));
		FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->idir,3,
					sizeof(int));
	    	FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->obstacle_comp,
					3,sizeof(COMPONENT));
		FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->preset_bound,1,
                                        sizeof(boolean));
                FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->var_min,1,
                                        sizeof(double));
                FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->var_max,1,
                                        sizeof(double));
		if (movie_option->plot_temperature)
		{
		    if (movie_option->plot_cross_section[0])
		    {
	    	    	sprintf(hdf_movie_var->var_name[n],"temperature-yz");
	    	    	hdf_movie_var->get_state_var[n] = getStateTemperature;
		    	hdf_movie_var->top_var[n] = 
				eqn_params->field->temperature;
			hdf_movie_var->idir[n] = 0;
		    	hdf_movie_var->obstacle_comp[n] = ERROR_COMP;
			hdf_movie_var->num_var = ++n;
		    }
		    if (movie_option->plot_cross_section[1])
		    {
			sprintf(hdf_movie_var->var_name[n],"temperature-xz");
			hdf_movie_var->get_state_var[n] = getStateTemperature;
		    	hdf_movie_var->top_var[n] = 
				eqn_params->field->temperature;
			hdf_movie_var->idir[n] = 1;
		    	hdf_movie_var->obstacle_comp[n] = ERROR_COMP;
			hdf_movie_var->num_var = ++n;
		    }
		    if (movie_option->plot_cross_section[2])
		    {
	    	    	sprintf(hdf_movie_var->var_name[n],"temperature-xy");
	    	    	hdf_movie_var->get_state_var[n] = getStateTemperature;
		    	hdf_movie_var->top_var[n] = 
				eqn_params->field->temperature;
			hdf_movie_var->idir[n] = 2;
		    	hdf_movie_var->obstacle_comp[n] = ERROR_COMP;
			hdf_movie_var->num_var = ++n;
		    }
		}
	    }
	    front->hdf_movie_var = hdf_movie_var;
	}
        if (debugging("trace"))
            printf("Leaving initMovieVariables()\n");
}	/* end initMovieVariables */


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

void CARTESIAN::initSampleTemperature(char *in_name)
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
}       /* end initSampleTemperature */

void CARTESIAN::sampleTemperature()
{
        if (dim == 2)
            sampleTemperature2d();
        else if (dim == 3)
            sampleTemperature3d();
}       /* end sampleTemperature */

void CARTESIAN::sampleTemperature2d()
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
	    sprintf(sname, "%s/t-%d.xg",dirname,count);
	    sfile = fopen(sname,"w");
	    fprintf(sfile,"Next\n");
	    fprintf(sfile,"color=%s\n",sample_color[count]);
	    fprintf(sfile,"thickness=1.5\n");
	    for (j = jmin; j <= jmax; ++j)
	    {
		index = d_index2d(i,j,top_gmax);
		var1 = field->temperature[index];
		index = d_index2d(i+1,j,top_gmax);
		var2 = field->temperature[index];
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
	    sprintf(sname, "%s/t-%d.xg",dirname,count);
            sfile = fopen(sname,"w");
            fprintf(sfile,"Next\n");
            fprintf(sfile,"color=%s\n",sample_color[count]);
            fprintf(sfile,"thickness=1.5\n");
            for (i = imin; i <= imax; ++i)
            {
                index = d_index2d(i,j,top_gmax);
                var1 = field->temperature[index];
                index = d_index2d(i,j+1,top_gmax);
                var2 = field->temperature[index];
                var = (var1 + lambda*var2) / (1.0 + lambda);
                getRectangleCenter(index,coords);
                fprintf(sfile,"%20.14f   %20.14f\n",coords[0],var);
            }
            fclose(sfile);
            break;
	default:
            printf("Incorrect input for sample temperature!\n");
            break;
	}
	count++;
}	/* end sampleTemperature2d */

void CARTESIAN::sampleTemperature3d()
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
                    sprintf(sname, "%s/t-%d.xg",dirname,count);
                    sfile = fopen(sname,"w");
                    fprintf(sfile,"Next\n");
                    fprintf(sfile,"color=%s\n",sample_color[count]);
                    fprintf(sfile,"thickness=1.5\n");
                    for (k = kmin; k <= kmax; ++k)
                    {
                        index = d_index3d(i,j,k,top_gmax);
                        var1 = field->temperature[index];
                        index = d_index3d(i+1,j,k,top_gmax);
                        var2 = field->temperature[index];
                        var_tmp1 = (var1 + lambda1*var2)/(1.0 + lambda1);

                        index = d_index3d(i,j+1,k,top_gmax);
                        var1 = field->temperature[index];
                        index = d_index3d(i+1,j+1,k,top_gmax);
                        var2 = field->temperature[index];
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
                    sprintf(sname, "%s/t-%d.xg",dirname,count);
                    sfile = fopen(sname,"w");
                    fprintf(sfile,"Next\n");
                    fprintf(sfile,"color=%s\n",sample_color[count]);
                    fprintf(sfile,"thickness=1.5\n");
                    for (j = jmin; j <= jmax; ++j)
                    {
                        index = d_index3d(i,j,k,top_gmax);
                        var1 = field->temperature[index];
                        index = d_index3d(i+1,j,k,top_gmax);
                        var2 = field->temperature[index];
                        var_tmp1 = (var1 + lambda1*var2)/(1.0 + lambda1);

                        index = d_index3d(i,j,k+1,top_gmax);
                        var1 = field->temperature[index];
                        index = d_index3d(i+1,j,k+1,top_gmax);
                        var2 = field->temperature[index];
                        var_tmp2 = (var1 + lambda1*var2)/(1.0 + lambda1);

                        var = (var_tmp1 + lambda2*var_tmp2)/(1.0 + lambda2);
                        getRectangleCenter(index,coords);
                        fprintf(sfile,"%20.14f   %20.14f\n",coords[2],var);
                    }
                    fclose(sfile);
                    break;
                default:
                    printf("Incorrect input for sample temperature!\n");
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
                    sprintf(sname, "%s/t-%d.xg",dirname,count);
                    sfile = fopen(sname,"w");
                    fprintf(sfile,"Next\n");
                    fprintf(sfile,"color=%s\n",sample_color[count]);
                    fprintf(sfile,"thickness=1.5\n");
                    for (i = imin; i <= imax; ++i)
                    {
                        index = d_index3d(i,j,k,top_gmax);
                        var1 = field->temperature[index];
                        index = d_index3d(i,j+1,k,top_gmax);
                        var2 = field->temperature[index];
                        var_tmp1 = (var1 + lambda1*var2)/(1.0 + lambda1);

                        index = d_index3d(i,j,k+1,top_gmax);
                        var1 = field->temperature[index];
                        index = d_index3d(i,j+1,k+1,top_gmax);
                        var2 = field->temperature[index];
                        var_tmp2 = (var1 + lambda1*var2)/(1.0 + lambda1);

                        var = (var_tmp1 + lambda2*var_tmp2)/(1.0 + lambda2);
                        getRectangleCenter(index,coords);
                        fprintf(sfile,"%20.14f   %20.14f\n",coords[0],var);
                    }
                    fclose(sfile);
                    break;
                default:
                    printf("Incorrect input for sample temperature!\n");
                    break;
            }
        default:
            printf("Incorrect input for sample temperature!\n");
            break;
        }
        count++;
}       /* end sampleTemperature3d */

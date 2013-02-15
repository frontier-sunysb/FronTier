/*******************************************************************
 * 		CARTESIAN.c
 *******************************************************************/
/*
      To-Do: Line 392, the values of ilower and iupper 
      
*/

#include "solver.h"
#include "finance.h"

//-------------------------------------------------------------------
//		INC_STATE
//-------------------------------------------------------------------

void INC_STATE::setZero(void)
{
	P = 0.0;
	sigma[0] = 0.0;
	sigma[1] = 0.0;
	sigma[2] = 0.0;
	r = 0.0;
}


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
	int num_cells;

	// init vertices,edges & cell_center
	if (debugging("trace")) printf("Entering initMesh()\n");
	RECTANGLE       rectangle;

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
	    for (j = 0; j <= top_gmax[1]; j++)
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
	int size = (int)cell_center.size();
	
	// cell center components
	for (i = 0; i < size; i++)
	{
	    cell_center[i].comp = 
	    		getComponent(cell_center[i].icoords);
	}
}

void CARTESIAN::setInitialCondition(void)
{
	int i;
	double coords[MAXD];
	INTERFACE *intfc = front->interf;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
        STATE *sl,*sr;
	int size = (int)cell_center.size();

	FT_MakeGridIntfc(front);
	setDomain();

        next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
            FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
	    getInitialState(Coords(p),sl);
	    getInitialState(Coords(p),sr);
        }

	// cell_center
	for (i = 0; i < size; i++)
	{
	    getRectangleCenter(i, coords);
	    getInitialState(coords,cell_center[i].state);
	}
	copyOptionPrice();
}	/* end setInitialCondition */

void CARTESIAN::setIndexMap(void)
{
	static boolean first = YES;
	int i,j,k,ic,index;
	int llbuf[MAXD],uubuf[MAXD];

	if (first)
	{
	    first = NO;
	    switch (dim)
	    {
	    case 1:
	    	FT_VectorMemoryAlloc((POINTER*)&i_to_I,top_gmax[0]+1,INT);
	    	break;
	    case 2:
	    	FT_MatrixMemoryAlloc((POINTER*)&ij_to_I,top_gmax[0]+1,top_gmax[1]+1,INT);
	    	break;
	    case 3:
	    	FT_TriArrayMemoryAlloc((POINTER*)&ijk_to_I,top_gmax[0]+1,top_gmax[1]+1,top_gmax[2]+1,
					INT);
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
	    for (i = 0; i <= top_gmax[0]; i++)
		i_to_I[i] = -1;
	    for (i = imin; i <= imax; i++)
	    {
		ic = d_index1d(i,top_gmax);
		if (cell_center[ic].comp == BLACK_SCHOLES_COMP)
		{
	    	    i_to_I[i] = index + ilower;
	    	    index++;
		}
	    }
	    FT_ParallelExchCellIndex(front,llbuf,uubuf,(POINTER)i_to_I);
	    break;
	case 2:
	    for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
                    ij_to_I[i][j] = -1;
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		ic = d_index2d(i,j,top_gmax);
		if (cell_center[ic].comp == BLACK_SCHOLES_COMP)
		{
	    	    ij_to_I[i][j] = index + ilower;
	    	    index++;
		}
	    }
	    FT_ParallelExchCellIndex(front,llbuf,uubuf,(POINTER)ij_to_I);
	    break;
	case 3:
	    for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
                    ijk_to_I[i][j][k] = -1;
	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		ic = d_index3d(i,j,k,top_gmax);
		if (cell_center[ic].comp == BLACK_SCHOLES_COMP)
		{
	    	    ijk_to_I[i][j][k] = index + ilower;
	    	    index++;
		}
	    }
	    FT_ParallelExchCellIndex(front,llbuf,uubuf,(POINTER)ijk_to_I);
	    break;
	}
}	/* end setIndexMap */

void CARTESIAN::computeAdvection(void)
{
	if (debugging("trace")) printf("Entering computeAdvection()\n");
	return computeAdvectionImplicit();
}
    
void CARTESIAN::computeAdvectionExplicit(void)
{
	switch (dim)
	{
	case 1:
	    computeAdvectionExplicit1d();
	    return;
	case 2:
	    computeAdvectionExplicit2d();
	    return;
	case 3:
	    computeAdvectionExplicit3d();
	    return;
	}
}	/* end computeAdvectionExplicit */

void CARTESIAN::computeAdvectionExplicit1d(void)
{
	int i;
	double *sigma,r,D,E;
	int index0,index1,index2;
	double P0,P1,P2,crx_P;
	double S1;
        double coords[MAXD],crx_coords[MAXD];
	int icoords[MAXD];
	boolean fr_crx_grid_seg;
	COMPONENT comp;
	STATE *state;
	HYPER_SURF *hs;
	INTERFACE *grid_intfc = front->grid_intfc;

	if (debugging("trace")) printf("Entering computeAdvectionExplicit()\n");

	eqn_params = (PARAMS*)front->extra1;
	sigma = eqn_params->sigma;
	r = eqn_params->r;
	D = eqn_params->D;
	E = eqn_params->E;
	for (i = imin; i <= imax; i++)
	{
	    index1 = d_index1d(i,top_gmax);
	    comp = top_comp[index1];
	    if (comp != BLACK_SCHOLES_COMP) 
	    {
		coords[0] = top_L[0] + i*top_h[0];
                if (eqn_params->f_type == AMRI_CALL_OPTION)
                    array[index1] = coords[0] - E;
                else if (eqn_params->f_type == AMRI_PUT_OPTION)
                    array[index1] = E - coords[0];
		continue;
	    }

	    icoords[0] = i;
	    coords[0] = top_L[0] + i*top_h[0];
	    P1 = cell_center[index1].state.P;
	    index0 = d_index1d(i-1,top_gmax);
	    P0 = cell_center[index0].state.P;
	    S1 = cell_center[index1].coords[0]/top_h[0];
	    index2 = d_index1d(i+1,top_gmax);
	    P2 = cell_center[index2].state.P;

	    fr_crx_grid_seg = FT_StateStructAtGridCrossing(front,grid_intfc,
			icoords,WEST,comp,(POINTER*)&state,&hs,crx_coords);
            if (fr_crx_grid_seg)
            {
                if (wave_type(hs) == DIRICHLET_BOUNDARY)
                {
                    P0 = getStatePrice(state);
                }
                else
                {
                    if (eqn_params->f_type == AMRI_PUT_OPTION)
                    {
                        P0 = extend_from_put_exc(P1,coords,crx_coords,top_h,E);
                        P0 = E - coords[0] + top_h[0];
                    }
                }
            }
	    fr_crx_grid_seg = FT_StateStructAtGridCrossing(front,grid_intfc,
			icoords,EAST,comp,(POINTER*)&state,&hs,crx_coords);
            if (fr_crx_grid_seg)
            {
                if (wave_type(hs) == DIRICHLET_BOUNDARY)
                {
                    P2 = getStatePrice(state);
                }
                else
                {
                    if (eqn_params->f_type == AMRI_CALL_OPTION)
                    {
                        P2 = extend_from_call_exc(P1,coords,crx_coords,top_h,E);
                        //P2 = coords[0] + top_h[0] - E;
                    }
                }
                printf("P2 = %f\n",P2);
	    }
            

	    /*
	    fr_crx_grid_seg = FT_StateVarAtGridCrossing(front,icoords,
                              WEST,comp,price_of_state,&crx_P,crx_coords);
            if (fr_crx_grid_seg)
            {
                if (eqn_params->f_type == AMRI_PUT_OPTION)
                {
                    P0 = extend_from_put_exc(P1,coords,crx_coords,top_h,E);
                    P0 = E - coords[0] + top_h[0];
                }
            }

	    fr_crx_grid_seg = FT_StateVarAtGridCrossing(front,icoords,
                              EAST,comp,price_of_state,&crx_P,crx_coords);
	    if (fr_crx_grid_seg)
	    {
		if (eqn_params->f_type == AMRI_CALL_OPTION)
		{
		    double a,b,c;
		    double S;
		    a = eqn_params->a;
		    b = eqn_params->b;
		    c = eqn_params->c;
		    S = top_L[0] + (i+1)*top_h[0];
		    P2 = a*sqr(S) + b*S + c;
		    P2 = S - E;
		}
	    }
	    */

	    array[index1] = P1 + m_dt*(
			0.5*sqr(sigma[0]*S1)*(P0 - 2.0*P1 + P2) 
			+ (r-D)*S1*(P2 - P1)
			- r*P1);
	}
	scatMeshArray();
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index1 = d_index1d(i,top_gmax);
	    cell_center[index1].state.P = array[index1];
	}

	if (debugging("trace")) printf("Leaving computeAdvectionExplicit()\n");
}	/* end computeAdvectionExplicit1d */

void CARTESIAN::computeAdvectionExplicit2d(void)
{
	(void) printf("computeAdvectionExplicit2d() not implemented!\n");
	clean_up(ERROR);
}	/* end computeAdvectionExplicit2d */

void CARTESIAN::computeAdvectionExplicit3d(void)
{
	(void) printf("computeAdvectionExplicit3d() not implemented!\n");
	clean_up(ERROR);
}	/* end computeAdvectionExplicit3d */

void CARTESIAN::computeAdvectionImplicit(void)
{
	switch (dim)
	{
	case 1:
	    computeAdvectionImplicit1d();
	    return;
	case 2:
	    computeAdvectionImplicit2d();
	    return;
	case 3:
	    computeAdvectionImplicit3d();
	    return;
	}
}	/* end computeAdvectionImplicit */

void CARTESIAN::computeAdvectionImplicit1d(void)
{
	double *sigma,r,E,D;
	double P0,P1,P2,P_crx;
        int i,index[3],I[3],size;
	double coeff[3],rhs;
	PETSc solver;
	double *x;
	int num_iter = 0;
	double rel_residual = 0;
	double S;
	COMPONENT comp;
	double coords[MAXD],crx_coords[MAXD],ext_coords[MAXD];
        int icoords[MAXD];
	boolean fr_crx_grid_seg;
	HYPER_SURF *hs;
	STATE *state;
	double time = front->time;
	INTERFACE *grid_intfc = front->grid_intfc;

	if (debugging("trace")) 
	    printf("Entering computeAdvectionImplicit()\n");
	
	eqn_params = (PARAMS*)front->extra1;
	sigma = eqn_params->sigma;
	r = eqn_params->r;
	E = eqn_params->E;
        D = eqn_params->D;

	solver.Create(ilower, iupper-1, 3, 3);
	size = iupper - ilower;
	setIndexMap();
        
        for (i = imin; i <= imax; i++)
	{
	    index[0] = d_index1d(i-1,top_gmax);
	    index[1] = d_index1d(i,top_gmax);
	    index[2] = d_index1d(i+1,top_gmax);
	    comp = top_comp[index[1]];

	    if (comp != BLACK_SCHOLES_COMP) continue;
	    icoords[0] = i;

	    I[0] = i_to_I[i-1];
	    I[1] = i_to_I[i];
	    I[2] = i_to_I[i+1];
	    coords[0] = top_L[0] + i*top_h[0];
            S = coords[0]/top_h[0];
             
	    coeff[0] = 0.5*((r-D)*S - sqr(sigma[0]*S))*m_dt;
	    coeff[1] = 1 + (sqr(sigma[0]*S)+r)*m_dt;
	    coeff[2] = 0.5*(-(r-D)*S - sqr(sigma[0]*S))*m_dt;
	    P0 = cell_center[index[0]].state.P;
	    P1 = cell_center[index[1]].state.P;
	    P2 = cell_center[index[2]].state.P;
	    rhs = P1;

	    fr_crx_grid_seg = FT_StateStructAtGridCrossing(front,grid_intfc,
			icoords,WEST,comp,(POINTER*)&state,&hs,crx_coords);
            if (fr_crx_grid_seg)
            {
		P_crx = getStatePrice(state);
		ext_coords[0] = coords[0] - top_h[0];
                P0 = linear_extension(P1,coords[0],P_crx,crx_coords[0],
				ext_coords[0],top_h[0]);
                rhs += -coeff[0]*P0;
            }
	    fr_crx_grid_seg = FT_StateStructAtGridCrossing(front,grid_intfc,
			icoords,EAST,comp,(POINTER*)&state,&hs,crx_coords);
            if (fr_crx_grid_seg)
            {
		P_crx = getStatePrice(state);
		ext_coords[0] = coords[0] + top_h[0];
                P2 = linear_extension(P1,coords[0],P_crx,crx_coords[0],
				ext_coords[0],top_h[0]);
                rhs += -coeff[2]*P2;
            }

            if (I[0] != -1) solver.Add_A(I[1],I[0],coeff[0]);
            solver.Add_A(I[1],I[1],coeff[1]);
            if (I[2] != -1) solver.Add_A(I[1],I[2],coeff[2]);
            solver.Add_b(I[1],rhs);
	}
	solver.SetMaxIter(40000);
        solver.SetTol(1e-14);
        solver.Solve_GMRES();
        solver.GetNumIterations(&num_iter);
        solver.GetFinalRelativeResidualNorm(&rel_residual);

	if (debugging("PETSc"))
	{
	    (void) printf("CARTESIAN::computeAdvectionImplicit: "
	       		"num_iter = %d, rel_residual = %g \n", 
			num_iter, rel_residual);
	}

	FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));
	solver.Get_x(x);

        for (i = imin; i <= imax; i++)
	{
	    I[0] = i_to_I[i];
	    index[0]  = d_index1d(i,top_gmax);
	    comp = top_comp[index[0]];
            if (comp == BLACK_SCHOLES_COMP)
            {
                array[index[0]] = x[I[0]-ilower];
            }
            else
            {
                coords[0] = top_L[0] + i*top_h[0];
                if (eqn_params->f_type == AMRI_CALL_OPTION)
                    array[index[0]] = coords[0] - E;
                else if (eqn_params->f_type == AMRI_PUT_OPTION)
                    array[index[0]] = E - coords[0];
            }
	}
	scatMeshArray();
	
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index[0]  = d_index1d(i,top_gmax);
	    cell_center[index[0]].state.P = array[index[0]];
	}
	
	FT_FreeThese(1,x);
	if (debugging("trace")) printf("Leaving computeAdvectionImplicit()\n");
}	/* end computeAdvectionImplicit1d */

void CARTESIAN::computeAdvectionImplicit2d(void)
{
	(void) printf("computeAdvectionImplicit2d() not implemented!\n");
	clean_up(ERROR);
}	/* end computeAdvectionImplicit2d */

void CARTESIAN::computeAdvectionImplicit3d(void)
{
	(void) printf("computeAdvectionImplicit3d() not implemented!\n");
	clean_up(ERROR);
}	/* end computeAdvectionImplicit3d */

void CARTESIAN::computeSourceTerm(double *coords, double t, INC_STATE &state) 
{
	//computeSourceTerm(coords, state);
}

void CARTESIAN::getInitialState(double *coords, STATE *state) 
{
	switch (dim)
	{
	case 1:
	    if (eqn_params->f_type == EURO_CALL_OPTION ||
		eqn_params->f_type == AMRI_CALL_OPTION)
	    	*state = (STATE)std::max(0.0,coords[0]-eqn_params->E);
	    else if (eqn_params->f_type == EURO_PUT_OPTION ||
		eqn_params->f_type == AMRI_PUT_OPTION)
	    	*state = (STATE)std::max(0.0,eqn_params->E-coords[0]);
	    break;
	case 2:
	    /*TODO*/
	    break;
	case 3:
	    /*TODO*/
	    break;
	}
}

void CARTESIAN::getInitialState(double *coords, INC_STATE &state) 
{
	switch (dim)
	{
	case 1:
	    if (eqn_params->f_type == EURO_CALL_OPTION ||
		eqn_params->f_type == AMRI_CALL_OPTION)
	    	state.P = std::max(0.0,coords[0]-eqn_params->E);
	    else if (eqn_params->f_type == EURO_PUT_OPTION ||
		eqn_params->f_type == AMRI_PUT_OPTION)
	    	state.P = std::max(0.0,eqn_params->E-coords[0]);
	    break;
	case 2:
	    /*TODO*/
	    break;
	case 3:
	    /*TODO*/
	    break;
	}
}

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

	setComponent();
	if (debugging("trace")) printf("Passing setComponent()\n");

	setGlobalIndex();
	if (debugging("trace")) printf("Passing setGlobalIndex()\n");

	copyBackOptionPrice();
	computeAdvection();
	if (debugging("trace")) printf("Passing computeAdvection()\n");

	copyOptionPrice();
	if (debugging("trace")) printf("Passing copyOptionPrice()\n");
	
	setAdvectionDt();
	if (debugging("trace")) printf("Passing setAdvectionDt()\n");

	stop_clock("solve");
	if (debugging("trace")) printf("Leaving solve()\n");
}

void CARTESIAN::setAdvectionDt()
{
	double *sigma,r;
	double max_D_coeff,max_drift;
	double dt_parab,dt_hyper;

	eqn_params = (PARAMS*)front->extra1;
	sigma = eqn_params->sigma;
	r = eqn_params->r;

	max_drift = r*top_U[0];
	dt_hyper = top_h[0]/max_drift;
	m_dt = dt_hyper;
}	/* end setAdvectionDt */

void CARTESIAN::getRectangleCenter(
	int index, 
	double *coords)
{
	int i;
	for (i = 0; i < dim; ++i)
	    coords[i] = cell_center[index].coords[i];
}

int CARTESIAN::getComponent(
	int *icoords)
{
	int index;
	index = d_index(icoords,top_gmax,dim);
	return top_comp[index];
}

void CARTESIAN::save(char *filename)
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

CARTESIAN::CARTESIAN(Front &front):front(&front)
{
}

void CARTESIAN::setDomain()
{
	static boolean first = YES;
	INTERFACE *grid_intfc;
	Table *T;
	static double *option_price;
	static double *rate_of_change;

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
	eqn_params = (PARAMS*)front->extra1;

	switch (dim)
	{
	case 1:
            if (first)
            {
		FT_ScalarMemoryAlloc((POINTER*)&field,sizeof(FIELD));
                FT_VectorMemoryAlloc((POINTER*)&array,top_gmax[0]+1,FLOAT);
                FT_VectorMemoryAlloc((POINTER*)&option_price,
					top_gmax[0]+1,FLOAT);
                FT_VectorMemoryAlloc((POINTER*)&rate_of_change,
					top_gmax[0]+1,FLOAT);
                first = NO;
            }	
	    imin = (lbuf[0] == 0) ? 1 : lbuf[0];
	    imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
	    eqn_params->option_price = option_price;
	    eqn_params->rate_of_change = rate_of_change;
	    eqn_params = (PARAMS*)front->extra1;
	    break;
	case 2:
	    if (first)
	    {
		FT_ScalarMemoryAlloc((POINTER*)&field,sizeof(FIELD));
	    	FT_VectorMemoryAlloc((POINTER*)&array,
					(top_gmax[0]+1)*(top_gmax[1]+1),FLOAT);
	    	FT_VectorMemoryAlloc((POINTER*)&option_price,
					(top_gmax[0]+1)*(top_gmax[1]+1),FLOAT);
	    	first = NO;
	    }
	    imin = (lbuf[0] == 0) ? 1 : lbuf[0];
	    jmin = (lbuf[1] == 0) ? 1 : lbuf[1];
	    imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
	    jmax = (ubuf[1] == 0) ? top_gmax[1] - 1 : top_gmax[1] - ubuf[1];
	    eqn_params->option_price = option_price;
	    eqn_params->rate_of_change = rate_of_change;
	    break;
	case 3:
	    if (first)
	    {
		FT_ScalarMemoryAlloc((POINTER*)&field,sizeof(FIELD));
	    	FT_VectorMemoryAlloc((POINTER*)&array,
			(top_gmax[0]+1)*(top_gmax[1]+1)*(top_gmax[2]+1),FLOAT);
	    	FT_VectorMemoryAlloc((POINTER*)&option_price,
			(top_gmax[0]+1)*(top_gmax[1]+1)*(top_gmax[2]+1),FLOAT);
	    	first = NO;
	    }
	    imin = (lbuf[0] == 0) ? 1 : lbuf[0];
	    jmin = (lbuf[1] == 0) ? 1 : lbuf[1];
	    kmin = (lbuf[2] == 0) ? 1 : lbuf[2];
	    imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
	    jmax = (ubuf[1] == 0) ? top_gmax[1] - 1 : top_gmax[1] - ubuf[1];
	    kmax = (ubuf[2] == 0) ? top_gmax[2] - 1 : top_gmax[2] - ubuf[2];
	    eqn_params->option_price = option_price;
	    eqn_params->rate_of_change = rate_of_change;
	    break;
	}
	field->option_price = option_price;
	eqn_params->field = field;
}	/* end setDomain */

void CARTESIAN::scatMeshArray()
{
	FT_ParallelExchGridArrayBuffer(array,front);
}

void CARTESIAN::setGlobalIndex()
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
		if (cell_center[ic].comp != BLACK_SCHOLES_COMP) continue;
		NLblocks++;
	    }
	    break;
	case 2:
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		ic = d_index2d(i,j,top_gmax);
		if (cell_center[ic].comp != BLACK_SCHOLES_COMP) continue;
		NLblocks++;
	    }
	    break;
	case 3:
	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		ic = d_index3d(i,j,k,top_gmax);
		if (cell_center[ic].comp != BLACK_SCHOLES_COMP) continue;
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


void CARTESIAN::copyOptionPrice()
{
	int i,j,k,index;
	double *option_price = eqn_params->option_price;
	double *rate_of_change = eqn_params->rate_of_change;

	switch (dim)
	{
	case 1:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    {
	    	index  = d_index1d(i,top_gmax);
	    	rate_of_change[index] = (cell_center[index].state.P -
				option_price[index])/m_dt;
	    	option_price[index] = cell_center[index].state.P;
	    }
	    break;
	case 2:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    {
	    	index  = d_index2d(i,j,top_gmax);
	    	rate_of_change[index] = (cell_center[index].state.P -
				option_price[index])/m_dt;
	    	option_price[index] = cell_center[index].state.P;
	    }
	    break;
	case 3:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    for (k = 0; k <= top_gmax[2]; ++k)
	    {
	    	index  = d_index3d(i,j,k,top_gmax);
	    	rate_of_change[index] = (cell_center[index].state.P -
				option_price[index])/m_dt;
	    	option_price[index] = cell_center[index].state.P;
	    }
	    break;
	}
}	/* end copyOptionPrice */

void CARTESIAN::copyBackOptionPrice()
{
	int i,j,k,index;
	double *option_price = eqn_params->option_price;

	switch (dim)
	{
	case 1:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    {
	    	index  = d_index1d(i,top_gmax);
	    	cell_center[index].state.P = option_price[index];
	    }
	    break;
	case 2:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    {
	    	index  = d_index2d(i,j,top_gmax);
	    	cell_center[index].state.P = option_price[index];
	    }
	    break;
	case 3:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    for (k = 0; k <= top_gmax[2]; ++k)
	    {
	    	index  = d_index3d(i,j,k,top_gmax);
	    	cell_center[index].state.P = option_price[index];
	    }
	    break;
	}
}	/* end copyBackOptionPrice */

void CARTESIAN::printOptionPrice(char *out_name)
{
	int i,j,k,index;
	char filename[100];
	FILE *outfile;

	sprintf(filename,"%s/state.ts%s",out_name,right_flush(front->step,7));
#if defined(__MPI__)
        sprintf(filename,"%s-nd%s",filename,right_flush(pp_mynode(),4));
#endif /* defined(__MPI__) */
	outfile = fopen(filename,"w");
	
	fprintf(outfile,"\nInterior states:\n");
	switch (dim)
	{
	case 1:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    {
		index = d_index1d(i,top_gmax);
	        fprintf(outfile,"%20.16f\n",cell_center[index].state.P);
	    }
	    break;
	case 2:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    {
		index = d_index2d(i,j,top_gmax);
	        fprintf(outfile,"%20.16f\n",cell_center[index].state.P);
	    }
	    break;
	case 3:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    for (k = 0; k <= top_gmax[2]; ++k)
	    {
		index = d_index3d(i,j,k,top_gmax);
	        fprintf(outfile,"%20.16f\n",
			cell_center[index].state.P);
	    }
	    break;
	}
	fclose(outfile);
}

void CARTESIAN::readOptionPrice(char *restart_name)
{
	FILE *infile;
	int i,j,k,index;

	FT_MakeGridIntfc(front);
        setDomain();

	infile = fopen(restart_name,"r");

	next_output_line_containing_string(infile,"Interior states:");

	switch (dim)
	{
	case 1:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    {
		index = d_index1d(i,top_gmax);
	    	fscanf(infile,"%lf",&cell_center[index].state.P);
	    }
	    break;
	case 2:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    {
		index = d_index2d(i,j,top_gmax);
	    	fscanf(infile,"%lf",&cell_center[index].state.P);
	    }
	    break;
	case 3:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    for (k = 0; k <= top_gmax[2]; ++k)
	    {
		index = d_index3d(i,j,k,top_gmax);
	    	fscanf(infile,"%lf",&cell_center[index].state.P);
	    }
	    break;
	}
	fclose(infile);
	copyOptionPrice();
}

void CARTESIAN::initMovieVariables()
{
	int n;
	static HDF_MOVIE_VAR *hdf_movie_var;
	F_MOVIE_OPTION *movie_option = eqn_params->movie_option;

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
                if (movie_option->plot_price)
                {
                    sprintf(hdf_movie_var->var_name[0],"price");
                    hdf_movie_var->get_state_var[0] = getStatePrice;
                    hdf_movie_var->top_var[0] = eqn_params->field->option_price;
                    hdf_movie_var->obstacle_comp[0] = EXCERCISE_COMP;
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
		if (movie_option->plot_price)
		{
	    	    sprintf(hdf_movie_var->var_name[0],"price");
	    	    hdf_movie_var->get_state_var[0] = getStatePrice;
		    hdf_movie_var->top_var[0] = eqn_params->field->option_price;
		    hdf_movie_var->obstacle_comp[0] = EXCERCISE_COMP;
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
		FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->preset_bound,3,
                                        sizeof(boolean));
                FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->var_min,3,
                                        sizeof(double));
                FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->var_max,3,
                                        sizeof(double));
		if (movie_option->plot_price)
		{
		    if (movie_option->plot_cross_section[0])
		    {
	    	    	sprintf(hdf_movie_var->var_name[n],"price-yz");
	    	    	hdf_movie_var->get_state_var[n] = getStatePrice;
		    	hdf_movie_var->top_var[n] = 
					eqn_params->field->option_price;
			hdf_movie_var->idir[n] = 0;
		    	hdf_movie_var->obstacle_comp[n] = EXCERCISE_COMP;
			hdf_movie_var->num_var = ++n;
		    }
		    if (movie_option->plot_cross_section[1])
		    {
			sprintf(hdf_movie_var->var_name[n],"price-xz");
			hdf_movie_var->get_state_var[n] = getStatePrice;
		    	hdf_movie_var->top_var[n] = 
					eqn_params->field->option_price;
			hdf_movie_var->idir[n] = 1;
		    	hdf_movie_var->obstacle_comp[n] = EXCERCISE_COMP;
			hdf_movie_var->num_var = ++n;
		    }
		    if (movie_option->plot_cross_section[2])
		    {
	    	    	sprintf(hdf_movie_var->var_name[n],"price-xy");
	    	    	hdf_movie_var->get_state_var[n] = getStatePrice;
		    	hdf_movie_var->top_var[n] = 
					eqn_params->field->option_price;
			hdf_movie_var->idir[n] = 2;
		    	hdf_movie_var->obstacle_comp[n] = EXCERCISE_COMP;
			hdf_movie_var->num_var = ++n;
		    }
		}
	    }
	    front->hdf_movie_var = hdf_movie_var;
	}
}	/* end initMovieVariables */

void CARTESIAN::printTimeData(char *out_name)
{
	char filename[100];
	static FILE *exS_file;
	static FILE *exV_file;
	static FILE *op_file;
	static FILE *v_file;
	double time = front->time;
	double Sf,Pf,OP,sample_S;
	double p[MAXD],sample_P;

	if (op_file == NULL)
	{
	    sprintf(filename,"%s/sample_op.xg",out_name);
	    op_file = fopen(filename,"w");
	    fprintf(op_file,"\"Option Price\"\n");

	    sprintf(filename,"%s/volatility.xg",out_name);
	    v_file = fopen(filename,"w");
	    fprintf(v_file,"\"Volatility\"\n");

	    if (eqn_params->f_type == AMRI_PUT_OPTION ||
		eqn_params->f_type == AMRI_CALL_OPTION)
	    {
	    	sprintf(filename,"%s/exc_S.xg",out_name);
	    	exS_file = fopen(filename,"w");
	    	fprintf(exS_file,"\"Excercise Point\"\n");

	    	sprintf(filename,"%s/exc_V.xg",out_name);
	    	exV_file = fopen(filename,"w");
	    	fprintf(exV_file,"\"Excercise Value\"\n");
	    }
	}

	p[0] = sample_S = eqn_params->sample_S;
	FT_IntrpStateVarAtCoords(front,BLACK_SCHOLES_COMP,p,
		eqn_params->field->option_price,getStatePrice,
		&sample_P,NULL);

	fprintf(op_file,"%f %f\n",time,sample_P);
	fflush(op_file);

	fprintf(v_file,"%f %f\n",time,eqn_params->sigma[0]);
	fflush(v_file);

	if (eqn_params->f_type == AMRI_PUT_OPTION ||
            eqn_params->f_type == AMRI_CALL_OPTION)
	{
	    excercise_point(front->interf,&Sf,&Pf);
	    fprintf(exS_file,"%f %f\n",time,Sf);
	    fflush(exS_file);
	    fprintf(exV_file,"%f %f\n",time,Pf);
	    fflush(exV_file);
	}
}


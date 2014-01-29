/*******************************************************************
 * 		VCARTESIAN.c
 *******************************************************************/
/*
      To-Do: Line 392, the values of ilower and iupper 
      
*/

#include <iFluid.h>
#include "solver.h"
#include "climate.h"

static int find_state_at_crossing(Front*,int*,GRID_DIRECTION,int,
                                POINTER*,HYPER_SURF**,double*);

// 		VCARTESIAN
//--------------------------------------------------------------------------

VCARTESIAN::~VCARTESIAN()
{
}

//---------------------------------------------------------------
//	initMesh
// include the following parts
// 1) setup edges:		
// 2) setup cell_center
//---------------------------------------------------------------
void VCARTESIAN::initMesh(void)
{
	int i,j,k,index;
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

void VCARTESIAN::setComponent(void)
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
                    field->vapor[i] = state->vapor;
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
                            temp_nb += field->vapor[index];
                            n++;
                        }
                        ic[ii] = (icoords[ii] == top_gmax[ii]) ? top_gmax[ii]
                                        : icoords[ii] + 1;
                        index = d_index(ic,top_gmax,dim);
                        if (cell_center[index].comp != -1 &&
                            cell_center[index].comp == top_comp[index])
                        {
                            temp_nb += field->vapor[index];
                            n++;
                        }
                    }
                    field->vapor[i] = temp_nb/n;
                }
	    }
	    cell_center[i].comp = top_comp[i];
        }
}	/* end setComponent */

void VCARTESIAN::setInitialCondition(void)
{
	int i;
	double coords[MAXD];
	INTERFACE *intfc = front->interf;
	POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
	STATE *sl,*sr;
	int c;
	
	short unsigned int seed[3] = {2,72,7172};
        GAUSS_PARAMS gauss_params;
        gauss_params.mu = 3.2;
        gauss_params.sigma = 0.2;

	FT_MakeGridIntfc(front);
        setDomain();

	/* Initialize states at the interface */
        next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
	    FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
            sl->vapor = sr->vapor = 0;
        }

	// cell_center
	for (i = 0; i < cell_center.size(); i++)
	{
	    field->vapor[i] = 0;
	    field->supersat[i] = 0;
	}
	for (i = 0; i < cell_center.size(); i++)
	{
	    c = top_comp[i];
	    if (c == LIQUID_COMP2)
	    	field->vapor[i] = gauss_center_limit((POINTER)&gauss_params,seed);
		//field->vapor[i] = 3.2;
	    else if (c == SOLID_COMP)
	    	field->vapor[i] = 0;
	    else
	    	field->vapor[i] = 3.2;
	}
}	/* end setInitialCondition */

void VCARTESIAN::setInitialVapor(void)
{
	int i,j,k,index;
	double coords[MAXD];
	INTERFACE *intfc = front->interf;
	POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
	STATE *sl,*sr;
	int c;
	char fname[100];
	FILE *infile;
	
	short unsigned int seed[3] = {2,72,7172};
        GAUSS_PARAMS gauss_params;
        gauss_params.mu = 3.2;
        gauss_params.sigma = 0.2;

	FT_MakeGridIntfc(front);
        setDomain();

	/* Initialize states at the interface */
        next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
	    FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
            sl->vapor = sr->vapor = 0;
        }

	// cell_center
	for (i = 0; i < cell_center.size(); i++)
	{
	    field->vapor[i] = 0;
	    field->supersat[i] = 0;
	}
	switch(dim)
	{
	    case 2: 
        	sprintf(fname,"./vapor2d");
        	infile = fopen(fname,"r");
        	for (j = jmin; j <= jmax; ++j)
        	for (i = imin; i <= imax; ++i)
        	{
            	    index = d_index2d(i,j,top_gmax);
            	    fscanf(infile,"%lf",&field->vapor[index]);
        	}
		break;
	    case 3: 
		sprintf(fname,"./vapor3d");
                infile = fopen(fname,"r");
		for (k = kmin; k <= kmax; ++k)
                for (j = jmin; j <= jmax; ++j)
                for (i = imin; i <= imax; ++i)
                {
                    index = d_index3d(i,j,k,top_gmax);
                    fscanf(infile,"%lf",&field->vapor[index]);
                }
                break;
	    default:
		printf("Unknown dim = %d\n",dim);
		clean_up(ERROR);
	}
}	/* end setInitialVapor */

void VCARTESIAN::setIndexMap(COMPONENT sub_comp)
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
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		ic = d_index2d(i,j,top_gmax);
		if (sub_comp == NO_COMP || cell_center[ic].comp == sub_comp)
		{
	    	    ij_to_I[i][j] = index + ilower;
		    I_to_ij[index][0] = i;
                    I_to_ij[index][1] = j;
	    	    index++;
		}
		else
		    ij_to_I[i][j] = -1;
	    }
	    FT_ParallelExchCellIndex(front,llbuf,uubuf,(POINTER)ij_to_I);
	    break;
	case 3:
	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		ic = d_index3d(i,j,k,top_gmax);
		if (sub_comp == NO_COMP || cell_center[ic].comp == sub_comp)
		{
	    	    ijk_to_I[i][j][k] = index + ilower;
		    I_to_ijk[index][0] = i;
                    I_to_ijk[index][1] = j;
                    I_to_ijk[index][2] = k;
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

void VCARTESIAN::computeAdvection()
{
	int i;
	COMPONENT sub_comp[2];

	if (eqn_params->num_scheme == UNSPLIT_IMPLICIT_CIM)
	    return computeAdvectionCim();

	sub_comp[0] = SOLID_COMP;
	sub_comp[1] = LIQUID_COMP2;
	
	for (i = 0; i < 2; ++i)
	{
		if(sub_comp[i] == SOLID_COMP)
		    continue;
	        setGlobalIndex(sub_comp[i]);
	        if (eqn_params->num_scheme == UNSPLIT_EXPLICIT)
	    	    computeAdvectionExplicit(sub_comp[i]);
	        else if (eqn_params->num_scheme == UNSPLIT_EXPLICIT_CIM)
	    	    computeAdvectionExplicitCim(sub_comp[i]);
	        else if (eqn_params->num_scheme == UNSPLIT_IMPLICIT)
	    	    computeAdvectionImplicit(sub_comp[i]);
	        else if (eqn_params->num_scheme == CRANK_NICOLSON)
	    	    computeAdvectionCN(sub_comp[i]);
	}
	printf("End solving advection\n\n");
}

void VCARTESIAN::computeAdvectionCN(COMPONENT sub_comp)
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
	double *Temp = field->vapor;
	double v[MAXD];
	double eta;

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
                                icoords,dir[l][m],comp,getStateVapor,
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
	    for (j = jmin; j <= jmax; ++j)
	    for (i = imin; i <= imax; ++i)
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
		if (source != NULL)
		    rhs += m_dt*source[ic];
		coeff = 1.0;
	 	for (l = 0; l < dim; ++l) v[l] = 0.0;
                if (field->vel != NULL)
                {
                    for (l = 0; l < dim; ++l)
                        v[l] = field->vel[l][ic];
                }
                for (l = 0; l < dim; ++l)
                {
                    lambda = D*m_dt/sqr(top_h[l]);
                    eta = v[l]*m_dt/top_h[l];
		    coeff += lambda;
		    
                    for (m = 0; m < 2; ++m)
                    {
                        next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax);
                        icn = d_index2d(ipn[0],ipn[1],top_gmax);
                        I_nb = ij_to_I[ipn[0]][ipn[1]];
                        coeff_nb = -0.5*lambda;
                        fr_crx_grid_seg = FT_StateVarAtGridCrossing(front,
                                icoords,dir[l][m],comp,getStateVapor,
                                &T_nb,crx_coords);
                        if (!fr_crx_grid_seg)
                        {
                            solver.Add_A(I,I_nb,coeff_nb);
                            T_nb = Temp[icn];
                            rhs -= -0.5*lambda*(T_nb - T0);
			    if(v[l] > 0 && m == 0)
				rhs += eta * (T_nb - T0);
			    if(v[l] < 0 && m == 1)
				rhs += eta * (T0 - T_nb);
                        }
                        else
                        {
                            rhs -= -0.5*lambda*(2.0*T_nb - T0);
                        }
                    }
                }
                solver.Add_A(I,I,coeff);
                solver.Add_b(I,rhs);
	    }
	    break;
	case 3:
	    solver.Create(ilower, iupper-1, 7, 7);
	    for (k = kmin; k <= kmax; ++k)
	    for (j = jmin; j <= jmax; ++j)
	    for (i = imin; i <= imax; ++i)
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
                if (source != NULL)
                    rhs += m_dt*source[ic];
		coeff = 1.0;
                for (l = 0; l < dim; ++l) v[l] = 0.0;
                if (field->vel != NULL)
                {
                    for (l = 0; l < dim; ++l)
                        v[l] = field->vel[l][ic];
                }
		for (l = 0; l < dim; ++l)
		{
                    eta = v[l]*m_dt/top_h[l];
		    lambda = D*m_dt/sqr(top_h[l]);
		    coeff += lambda;
                    for (m = 0; m < 2; ++m)
                    {
                        next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax);
                        icn = d_index3d(ipn[0],ipn[1],ipn[2],top_gmax);
			I_nb = ijk_to_I[ipn[0]][ipn[1]][ipn[2]];
			coeff_nb = -0.5*lambda;
                        fr_crx_grid_seg = FT_StateVarAtGridCrossing(front,
                                icoords,dir[l][m],comp,getStateVapor,
                                &T_nb,crx_coords);
                        if (!fr_crx_grid_seg)
                        {
			    solver.Add_A(I,I_nb,coeff_nb);
                            T_nb = Temp[icn];
                            rhs -= -0.5*lambda*(T_nb - T0);
                            if(v[l] > 0 && m == 0)
                                rhs += eta * (T_nb - T0);
                            if(v[l] < 0 && m == 1)
                                rhs += eta * (T0 - T_nb);
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
	    (void) printf("VCARTESIAN::computeAdvectionCN: "
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
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
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
            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
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
            for (j = 0; j <= top_gmax[1]; ++j)
            for (i = 0; i <= top_gmax[0]; ++i)
            {
		ic = d_index2d(i,j,top_gmax);
	    	comp = cell_center[ic].comp;
	    	if (comp == sub_comp)
		    Temp[ic] = array[ic];
	    }
	    break;
        case 3:
            for (k = 0; k <= top_gmax[2]; ++k)
            for (j = 0; j <= top_gmax[1]; ++j)
            for (i = 0; i <= top_gmax[0]; ++i)
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

void VCARTESIAN::computeAdvectionCim()
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
		source[index] = -field->vapor[index]/m_dt;
	    }
            break;
        case 2:
            parab_solver.ij_to_I = ij_to_I;
            parab_solver.I_to_ij = I_to_ij;
	    for (j = jmin; j <= jmax; ++j)
	    for (i = imin; i <= imax; ++i)
	    {
		index = d_index2d(i,j,top_gmax);
		source[index] = -field->vapor[index]/m_dt;
	    }
            break;
        case 3:
            parab_solver.ijk_to_I = ijk_to_I;
            parab_solver.I_to_ijk = I_to_ijk;
	    for (k = kmin; k <= kmax; ++k)
	    for (j = jmin; j <= jmax; ++j)
	    for (i = imin; i <= imax; ++i)
	    {
		index = d_index3d(i,j,k,top_gmax);
		source[index] = -field->vapor[index]/m_dt;
	    }
            break;
        }
	printf("Assigning class variables\n");
	parab_solver.w_type = ICE_PARTICLE_BOUNDARY;
	parab_solver.neg_comp = SOLID_COMP;
        parab_solver.pos_comp = LIQUID_COMP2;
        parab_solver.source = source;
	parab_solver.solutionJump = jumpT;
	parab_solver.gradJumpDotN = jumpEpsGradDotNorm;
        parab_solver.gradJumpDotT = jumpGradDotTan;
        parab_solver.diff_coeff[0] = eqn_params->D;
        parab_solver.diff_coeff[1] = eqn_params->D;
	parab_solver.jparams = (POINTER)eqn_params;
	parab_solver.findStateAtCrossing = find_state_at_crossing;
	parab_solver.getStateVar = getStateVapor;
        parab_solver.assignStateVar = assignStateTemperature;
	parab_solver.ilower = ilower;
        parab_solver.iupper = iupper;
        parab_solver.soln = field->vapor;
        parab_solver.size = iupper - ilower;
	parab_solver.solve_front_state = YES;
        parab_solver.set_solver_domain();
        parab_solver.solve(field->vapor);

	stop_clock("computeAdvectionCim");
	if (debugging("trace")) printf("Leaving computeAdvectionCim()\n");
	clean_up(0);
	return;
}	/* end computeAdvectionImplicit */

void VCARTESIAN::computeAdvectionImplicit(COMPONENT sub_comp)
{
	static PARABOLIC_SOLVER parab_solver(*front);
	static double *soln;
	int i,j,k,index;
	static boolean first = YES;
        
	if (soln == NULL)
            FT_VectorMemoryAlloc((POINTER*)&soln,comp_size,FLOAT);

	printf("Entering computeAdvectionImplicit()\n");
	if (debugging("trace")) printf("Entering computeAdvectionImplicit()\n");
	start_clock("computeAdvectionImplicit");
	setIndexMap(sub_comp);
	parab_solver.soln_comp = sub_comp;
        parab_solver.obst_comp = SOLID_COMP;
        parab_solver.var = field->vapor;
        parab_solver.soln = soln;
        parab_solver.getStateVarFunc = getStateVapor;
        parab_solver.findStateAtCrossing = find_state_at_crossing;
        parab_solver.source = source;
        parab_solver.D = eqn_params->D;
        parab_solver.order = eqn_params->pde_order;
        parab_solver.ilower = ilower;
        parab_solver.iupper = iupper;
        parab_solver.dt = m_dt;
	parab_solver.first = first;
        parab_solver.set_solver_domain();
	first = NO;
	
	if (sub_comp == LIQUID_COMP2)
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
		    field->vapor[index] = soln[index];
	    }
            break;
        case 2:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    {
	    	index = d_index2d(i,j,top_gmax);
	    	if (top_comp[index] == sub_comp)
		    field->vapor[index] = soln[index];
                /*for debugging*/
	    }
            break;
        case 3:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    for (k = 0; k <= top_gmax[2]; ++k)
	    {
	    	index = d_index3d(i,j,k,top_gmax);
	    	if (top_comp[index] == sub_comp)
		    field->vapor[index] = soln[index];
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
void VCARTESIAN::solve(double dt)
{
	printf("solving adv-diff equation\n");

	if (debugging("trace")) printf("Entering solve()\n");
	m_dt = dt;
	start_clock("solve");

	setDomain();
        if (debugging("trace")) printf("Passing setDomain()\n");

	setComponent();
	if (debugging("trace")) printf("Passing setComponent()\n");

	if(eqn_params->prob_type == PARTICLE_TRACKING)
	    computeSource();
	else
	    source = NULL;
	if (debugging("trace")) printf("Passing computeSource()\n");

	computeAdvection();
	if (debugging("trace")) printf("Passing computeAdvection()\n");
	
	computeSupersat();
	if (debugging("trace")) printf("Passing computeSupersat()\n");

	setAdvectionDt();
	if (debugging("trace")) printf("Passing setAdvectionDt()\n");

	setBoundary();
	if (debugging("trace")) printf("Passing setBoundary()\n");

	stop_clock("solve");
	if (debugging("trace")) printf("Leaving solve()\n");
}

void VCARTESIAN::computeSupersat()
{	
	int i,j,k, index;
	double *temp  = eqn_params->field->temperature;
	double *vapor  = eqn_params->field->vapor;
	double sat_vap_pre, sat_vap_rat;
	double *super = eqn_params->field->supersat;
	switch(dim)
	{
	    case 2:
		for (i = 0; i <= top_gmax[0]; ++i)
            	for (j = 0; j <= top_gmax[1]; ++j)
		{
		    index = d_index2d(i,j,top_gmax);
		    sat_vap_pre = 611.2*exp(17.67*(temp[index]-273.15)
					           /(temp[index]-29.65));
		    sat_vap_rat = 621.97 * sat_vap_pre
					/(eqn_params->field->pres[index]
					  -sat_vap_pre);
                    if(top_comp[index] == SOLID_COMP)
                        super[index] = 0;
                    else
                        super[index] = vapor[index]/sat_vap_rat - 1;
		    if(eqn_params->field->pres[index] == 0)
			super[index] = 0;
		}	
		break;
	    case 3:
                for (i = 0; i <= top_gmax[0]; ++i)
                for (j = 0; j <= top_gmax[1]; ++j)
		for (k = 0; k <= top_gmax[2]; ++k)
                {
                    index = d_index3d(i,j,k,top_gmax);
                    sat_vap_pre = 611.2*exp(17.67*(temp[index]-273.15)
						/(temp[index]-29.65));
                    sat_vap_rat = 621.97 * sat_vap_pre
					/(eqn_params->field->pres[index]);
		    if(top_comp[index] == SOLID_COMP)
			super[index] = 0;
		    else
                        super[index] = vapor[index]/sat_vap_rat - 1;
		    if(eqn_params->field->pres[index] == 0)
			super[index] = 0;
                } 
		break;
	}
}


void VCARTESIAN::setAdvectionDt()
{
	double D,Dl,Ds;
	static double m_dt_expl,m_dt_impl;  // explicit and implicit time step
	static boolean first = YES;

	if (first)
	{
	    first = NO;
	    eqn_params = (PARAMS*)front->extra2;
	    D = eqn_params->D;
	    m_dt_expl = 0.5*sqr(hmin)/D/(double)dim;
	    m_dt_impl = 0.5*hmin/D/(double)dim;
	    min_dt = 0.1*sqr(hmin)/D/(double)dim;
	}

	if (eqn_params->num_scheme == UNSPLIT_EXPLICIT ||
	    eqn_params->num_scheme == UNSPLIT_EXPLICIT_CIM)
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

void VCARTESIAN::getVelocity(double *p, double *U)
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

void VCARTESIAN::getRectangleIndex(int index, int &i, int &j)
{
	i = cell_center[index].icoords[0];
	j = cell_center[index].icoords[1];
}

void VCARTESIAN::getRectangleIndex(int index, int &i, int &j, int &k)
{
	i = cell_center[index].icoords[0];
	j = cell_center[index].icoords[1];
	k = cell_center[index].icoords[2];
}


int VCARTESIAN::getRectangleComponent(int index)
{	
	return getComponent(cell_center[index].icoords);
}

void VCARTESIAN::getRectangleCenter(
	int index, 
	double *coords)
{
	int i;
	for (i = 0; i < dim; ++i)
	    coords[i] = cell_center[index].coords[i];
}

void VCARTESIAN::getRectangleCenter(
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

int VCARTESIAN::getComponent(
	double *p)
{
	return component(p,front->interf);
}

int VCARTESIAN::getComponent(
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

void VCARTESIAN::save(char *filename)
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

VCARTESIAN::VCARTESIAN(Front &front):front(&front)
{
}

void VCARTESIAN::makeGridIntfc()
{
	static boolean first = YES;
	INTERFACE *grid_intfc;
	Table *T;
	static double *vapor;
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
                FT_VectorMemoryAlloc((POINTER*)&vapor,top_gmax[0]+1,FLOAT);
                first = NO;
            }	
	    imin = (lbuf[0] == 0) ? 1 : lbuf[0];
	    imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
	    eqn_params->field->vapor = vapor;
	    eqn_params = (PARAMS*)front->extra2;
	    break;
	case 2:
	    if (first)
	    {
	    	FT_VectorMemoryAlloc((POINTER*)&array,(top_gmax[0]+1)*(top_gmax[1]+1),FLOAT);
	    	FT_VectorMemoryAlloc((POINTER*)&vapor,(top_gmax[0]+1)*(top_gmax[1]+1),FLOAT);
	    	first = NO;
	    }
	    imin = (lbuf[0] == 0) ? 1 : lbuf[0];
	    jmin = (lbuf[1] == 0) ? 1 : lbuf[1];
	    imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
	    jmax = (ubuf[1] == 0) ? top_gmax[1] - 1 : top_gmax[1] - ubuf[1];
	    eqn_params->field->vapor = vapor;
	    break;
	case 3:
	    if (first)
	    {
	    	FT_VectorMemoryAlloc((POINTER*)&array,
			(top_gmax[0]+1)*(top_gmax[1]+1)*(top_gmax[2]+1),FLOAT);
	    	FT_VectorMemoryAlloc((POINTER*)&vapor,
			(top_gmax[0]+1)*(top_gmax[1]+1)*(top_gmax[2]+1),FLOAT);
	    	first = NO;
	    }
	    imin = (lbuf[0] == 0) ? 1 : lbuf[0];
	    jmin = (lbuf[1] == 0) ? 1 : lbuf[1];
	    kmin = (lbuf[2] == 0) ? 1 : lbuf[2];
	    imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
	    jmax = (ubuf[1] == 0) ? top_gmax[1] - 1 : top_gmax[1] - ubuf[1];
	    kmax = (ubuf[2] == 0) ? top_gmax[2] - 1 : top_gmax[2] - ubuf[2];
	    eqn_params->field->vapor = vapor;
	    break;
	}
}	/* end makeGridIntfc */

void VCARTESIAN::deleteGridIntfc()
{
	FT_FreeGridIntfc(front);
}

void VCARTESIAN::scatMeshArray()
{
	FT_ParallelExchGridArrayBuffer(array,front);
}

void VCARTESIAN::setGlobalIndex(COMPONENT sub_comp)
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


void VCARTESIAN::printFrontInteriorState(char *out_name)
{
        int i,j,k,l,index;
        char filename[100];
        FILE *outfile;
        INTERFACE *intfc = front->interf;
        STATE *sl,*sr;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
        double *Temp = field->vapor;
        double *Supersat = field->supersat;

        sprintf(filename,"%s/state.ts%s-vapor",out_name,
                        right_flush(front->step,7));
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
            fprintf(outfile,"%24.18g %24.18g\n",sl->vapor,
                                sr->vapor);
            fprintf(outfile,"%24.18g %24.18g\n",sl->supersat,
                                sr->supersat);
        }
        fprintf(outfile,"\nInterior states:\n");
        switch (dim)
        {
        case 1:
            for (i = 0; i <= top_gmax[0]; ++i)
            {
                index = d_index1d(i,top_gmax);
                fprintf(outfile,"%24.18g\n",Temp[index]);
                fprintf(outfile,"%24.18g\n",Supersat[index]);
            }
            break;
        case 2:
            for (i = 0; i <= top_gmax[0]; ++i)
            for (j = 0; j <= top_gmax[1]; ++j)
            {
                index = d_index2d(i,j,top_gmax);
                fprintf(outfile,"%24.18g\n",Temp[index]);
                fprintf(outfile,"%24.18g\n",Supersat[index]);
            }
            break;
        case 3:
            for (i = 0; i <= top_gmax[0]; ++i)
            for (j = 0; j <= top_gmax[1]; ++j)
            for (k = 0; k <= top_gmax[2]; ++k)
            {
                index = d_index3d(i,j,k,top_gmax);
                fprintf(outfile,"%24.18g\n",Temp[index]);
                fprintf(outfile,"%24.18g\n",Supersat[index]);
            }
            break;
        }
        fclose(outfile);
}

void VCARTESIAN::readFrontInteriorState(char *restart_name)
{
        FILE *infile;
        int i,j,k,l,index;
        INTERFACE *intfc = front->interf;
        STATE *sl,*sr;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
        double x;
        double *Temp = field->vapor;
        double *Supersat = field->supersat;

        char fname[100];
        sprintf(fname,"%s-vapor",restart_name);
        infile = fopen(fname,"r");

        /* Initialize states at the interface */
        next_output_line_containing_string(infile,"Interface states:");
        next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
            FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
            fscanf(infile,"%lf",&x);
            sl->vapor = x;
            fscanf(infile,"%lf",&x);
            sr->vapor = x;
            fscanf(infile,"%lf",&x);
            sl->supersat = x;
            fscanf(infile,"%lf",&x);
            sr->supersat = x;
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
                fscanf(infile,"%lf",&Supersat[index]);
            }
            break;
        case 2:
            for (i = 0; i <= top_gmax[0]; ++i)
            for (j = 0; j <= top_gmax[1]; ++j)
            {
                index = d_index2d(i,j,top_gmax);
                fscanf(infile,"%lf",&Temp[index]);
                fscanf(infile,"%lf",&Supersat[index]);
            }
            break;
        case 3:
            for (i = 0; i <= top_gmax[0]; ++i)
            for (j = 0; j <= top_gmax[1]; ++j)
            for (k = 0; k <= top_gmax[2]; ++k)
            {
                index = d_index3d(i,j,k,top_gmax);
                fscanf(infile,"%lf",&Temp[index]);
                fscanf(infile,"%lf",&Supersat[index]);
            }
            break;
        }
        fclose(infile);
}

void VCARTESIAN::setBoundary()
{
	int i,j,k,index0,index1;
	INTERFACE *intfc = front->interf;
	double *Temp = field->vapor;

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

void VCARTESIAN::vtk_plot3d(
        const char *varname,double *var)
{
	std::vector<int> ph_index;
	char *outname = front->out_name;
        int i,j,k,index;
        char dirname[256],filename[256];
        FILE *outfile;
        double coord_x,coord_y,coord_z,xmin,ymin,zmin;
        COMPONENT comp;
        int pointsx,pointsy,pointsz,num_points,num_cells,num_cell_list;
        int icoords[3],p_gmax[3];

        int ii,jj,kk;
        double ih,jh,kh;

        sprintf(filename, "%s/vtk/vtk.ts%s",outname,
                right_flush(front->step,7));
        //cell-based liquid phase
        ph_index.clear();
        if (!create_directory(filename,NO))
        {
            printf("Cannot create directory %s\n",filename);
            clean_up(ERROR);
        }
#if defined(__MPI__)
        if (pp_numnodes() > 1)
            sprintf(filename,"%s-nd%s",filename,right_flush(pp_mynode(),4));
#endif /* defined(__MPI__) */

        sprintf(filename,"%s/%s.vtk",filename,varname);
        outfile = fopen(filename,"w");
        fprintf(outfile,"# vtk DataFile Version 3.0\n");
        fprintf(outfile,"%s\n",varname);
        fprintf(outfile,"ASCII\n");
        fprintf(outfile,"DATASET UNSTRUCTURED_GRID\n");

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            if (cell_center[index].comp == LIQUID_COMP2)
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
        fprintf(outfile, "SCALARS %s double\n",varname);
        fprintf(outfile, "LOOKUP_TABLE default\n");
	if(var != NULL)
	{
	    for (i = 0; i < num_cells; i++)
            {
                index = ph_index[i];
                fprintf(outfile,"%f\n",var[index]);
            }
	}
	else
	{
	 	printf("The var is NULL!\n");
		clean_up(ERROR);
	}
        fclose(outfile);
}       /* end vtk_plot_vapor3d */

void VCARTESIAN::recordMixingLine()
{
	PARAMS *eqn_params = (PARAMS*)front->extra2;
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	PARTICLE *particle_array = eqn_params->particle_array;
	int num_drops = eqn_params->num_drops;
	static boolean first = YES;
	static double prev_RN[2] = {0,0};
	double mu = iFparams->mu2;
	double **vel = field->vel;
	double slope;
	double t_evap, t_mix;
	double DspRat,Omega,rv,r0; /*dissipation rate*/
	double Vortex[3]; /*vortex field*/
	double supersat;
	int i, j, k, l, m, m1, m2, ic, I, icoords[MAXD],size, count, nzeros;
	int gmin[MAXD], ipn[MAXD];
	int I_nb[3][2];
	const GRID_DIRECTION dir[3][2] =
                {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
	FILE *file;
	char fname[100];
	char *out_name = front->out_name;
	double L = top_U[0] - top_L[0];

	if (3 != dim)
	    return;
	double a3 = 1;
	for (i = 0; i < dim; i++) 
	{
	    gmin[i] = 0;
	    a3 *= top_h[i]; 
	}
	/*compute mean energy dissipation rate*/
	DspRat = 0;
	supersat = 0;
	count = 0;
	size = 0;
        for (k = kmin; k < kmax; k++)
        for (j = jmin; j < jmax; j++)
        for (i = imin; i < imax; i++)
        {
	    size++;
            ic = d_index3d(i,j,k,top_gmax);
	    if (field->drops[ic] == 0)
	    {
   	        supersat += field->supersat[ic];
		count ++;
	    }
            icoords[0] = i;
            icoords[1] = j;
            icoords[2] = k;
            I = ijk_to_I[i][j][k];
            for (m = 0; m < dim; ++m)
	    for (l = 0; l < 2; ++l)
            {
		if(next_ip_in_dir(icoords,dir[m][l],ipn,gmin,top_gmax))
                    I_nb[m][l] = ijk_to_I[ipn[0]][ipn[1]][ipn[2]];
		else
		{
		    printf("In recordMixingLine(), cannot find next ip\n");
	        printf("gmin=[%d %d %d]\n",gmin[0],gmin[1],gmin[2]);
		printf("gmax=[%d %d %d]\n",top_gmax[0],top_gmax[1],top_gmax[2]);
		printf("icoords=[%d %d %d]\n",icoords[0],icoords[1],icoords[2]);
		    clean_up(ERROR);
		}
            }
	    /*compute vortex: Nabla X U*/
	    for (l = 0; l < dim; l++)
	    {
		m1 = (l+1)%dim;
		m2 = (l+2)%dim;
	        Vortex[l]  =  (vel[m2][I_nb[m1][1]]-vel[m2][I_nb[m1][0]])/(2*top_h[m1]);
	        Vortex[l] -=  (vel[m1][I_nb[m2][1]]-vel[m1][I_nb[m2][0]])/(2*top_h[m2]);
	    }
	    Omega = 0; /*compute enstrophy: Omega=|Vortex|^2*/
	    for (l = 0; l < dim; ++l)
	    {
		Omega += sqr(Vortex[l]);
	    }
	    /*compute dissipation rate using: 2*mu*Omega*/ 
	    DspRat += 2*mu*Omega;
	}
	DspRat   /= size;
	supersat /= count;

	t_mix = pow(L*L/DspRat,1.0/3.0);
	
	rv = 0;
	r0 = 0;
	nzeros = 0;
	for (i = 0; i < num_drops; i++)
	{   
	    if (particle_array[i].radius != 0)
		nzeros ++;
	    rv += pow(particle_array[i].radius,3.0);
	    r0 += pow(particle_array[i].R0,3.0);
	}
	if (nzeros == 0)
	{
	   printf("No dropletes contained\n");
	   clean_up(0);
	}

	rv /= nzeros;
	r0 /= num_drops;

	if (supersat < 0)
	    t_evap = pow(rv,2.0/3.0)/(-eqn_params->K*supersat);
	else
	    t_evap = HUGE;

	slope = (double(nzeros)/double(num_drops)-prev_RN[1])
		/((rv/r0)-prev_RN[0]);
	prev_RN[0] = rv/r0;
	prev_RN[1] = double(nzeros)/double(num_drops);
	if (first)
	{
	    sprintf(fname,"%s/mixing",out_name);
	    file = fopen(fname,"w");
	    fprintf(file,"%%Damkoehler number VS slope of R-N diagram\n");
	    fprintf(file,"%%t_mix    t_evap    Damkoehler    slope\n");
	    fclose(file);

            sprintf(fname,"%s/RN",out_name);
            file = fopen(fname,"w");
            fprintf(file,"%%Rv    supersat    Rv^3/R0^3    N/N0\n");
            fclose(file);

            sprintf(fname,"%s/dissipation",out_name);
            file = fopen(fname,"w");
            fprintf(file,"%%Dissipation Rate\n");
            fclose(file);

	    first = NO;
	    return;
	}
	sprintf(fname,"%s/mixing",out_name);
	file = fopen(fname,"a");
	fprintf(file,"%15.14f  %15.14f  %15.14f  %15.14f\n",
		t_mix,t_evap,t_mix/t_evap,slope);
	fclose(file);
	/*plot mixing line: N/N0 VS. (Rv/Rv0)^3*/
        sprintf(fname,"%s/RN",out_name);
        file = fopen(fname,"a");
	fprintf(file,"%15.14e  %15.14f  %15.14f  %15.14f\n", 
		pow(rv,1.0/3.0),supersat,rv/r0, double(nzeros)/double(num_drops));
	fclose(file);
        sprintf(fname,"%s/dissipation",out_name);
        file = fopen(fname,"a");
	fprintf(file,"%20.19f\n",DspRat);
	fclose(file);
}

void VCARTESIAN::pointExplicitCimSolver(
	int *icoords,
	COMPONENT sub_comp)
{
	int l,m,ic,icn;
	int gmin[MAXD],ipn[MAXD],ipn2[MAXD];
	double coords[MAXD],crx_coords[MAXD];
	double vapor,vapor_nb[2],dgrad[MAXD],grad_plus[MAXD],
	       grad_minus[MAXD];
	double coef;
	COMPONENT comp;
	boolean fr_crx_grid_seg;
	const GRID_DIRECTION dir[3][2] = 
		{{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
	double *Temp = field->vapor;
	double **vel = field->vel;
	double v[MAXD],v_plus[MAXD],v_minus[MAXD];
	STATE *state;
	HYPER_SURF *hs;

	coef = eqn_params->D*m_dt;
        ic = d_index(icoords,top_gmax,dim);
        comp = top_comp[ic];
        if (comp != sub_comp) return;
        array[ic] = vapor = Temp[ic];
	for (l = 0; l < dim; ++l) gmin[l] = 0;

	for (l = 0; l < dim; ++l)
	{
	    v[l] = 0.0;
	    v_plus[l] = 0.0;
	    v_minus[l] = 0.0;
	}
	if (sub_comp == LIQUID_COMP2 && vel != NULL)
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
		fr_crx_grid_seg = FT_StateStructAtGridCrossing(front,
				front->grid_intfc,icoords,dir[l][m],comp,
				(POINTER*)&state,&hs,crx_coords);
                if (!fr_crx_grid_seg)
                {
                    next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax);
                    icn = d_index(ipn,top_gmax,dim);
		    vapor_nb[m] = Temp[icn];
                }
		else
		    vapor_nb[m] = getStateVapor((POINTER)state);

		if (!fr_crx_grid_seg) 
		{
                    dgrad[l] += (vapor_nb[m] - vapor)/top_h[l];
		} 
		else if (wave_type(hs) == NEUMANN_BOUNDARY)
		{
		    ;
		}
		else if (wave_type(hs) == DIRICHLET_BOUNDARY)
		{
                    dgrad[l] += (vapor_nb[m] - vapor)/top_h[l];
		}
		else
		{
		    double a, v1, v2, P[MAXD];
		    getRectangleCenter(ic, P);
		    a = fabs(crx_coords[l] - P[l])/top_h[l];
                    next_ip_in_dir(icoords,dir[l][1-m],ipn,gmin,top_gmax);
                    icn = d_index(ipn,top_gmax,dim);
		    v1 = Temp[icn]-vapor_nb[m];
                    next_ip_in_dir(ipn,dir[l][1-m],ipn2,gmin,top_gmax);
                    icn = d_index(ipn2,top_gmax,dim);
		    v2 = Temp[icn]-vapor_nb[m];
		    dgrad[l] += ((1-a)*v2+2*(a*a+a-1)/(1+a)*v1-(1+a)*
				(vapor-vapor_nb[m]))/top_h[l] - 
				(v1-vapor)/top_h[l];
		}
            }
	    grad_plus[l] = (vapor_nb[1] - vapor)/top_h[l];
	    grad_minus[l] = (vapor - vapor_nb[0])/top_h[l];
            array[ic] += coef*dgrad[l]/top_h[l]-m_dt*(v_plus[l]*
				grad_minus[l]+v_minus[l]*grad_plus[l]);
        }
}	/* end pointExplicitCimSolver */

void VCARTESIAN::computeAdvectionExplicitCim(COMPONENT sub_comp)
{
	int i,j,k,l,m,ic,icn,icoords[MAXD],nc;
	int gmin[MAXD],ipn[MAXD],ipn2[MAXD];
	int index0;
	double coords[MAXD],crx_coords[MAXD];
	double vapor,vapor_nb[2],dgrad[MAXD],grad_plus[MAXD],
	       grad_minus[MAXD];
	double coef;
	COMPONENT comp;
	boolean fr_crx_grid_seg;
	const GRID_DIRECTION dir[3][2] = 
		{{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
	double *Temp = field->vapor;
	double v[MAXD],**vel,v_plus[MAXD],v_minus[MAXD];

	start_clock("computeAdvectionExplicitCim");

	coef = eqn_params->D*m_dt;
	vel = eqn_params->field->vel;

	for (i = 0; i < dim; ++i) gmin[i] = 0;

	switch (dim)
	{
        case 1:
            for (i = imin; i <= imax; ++i)
            {
                icoords[0] = i;
		pointExplicitCimSolver(icoords,sub_comp);
            }
            break;
	case 2:
	    for (j = jmin; j <= jmax; ++j)
	    for (i = imin; i <= imax; ++i)
	    {
	    	icoords[0] = i;
	    	icoords[1] = j;
		pointExplicitCimSolver(icoords,sub_comp);
	    }
	    break;
	case 3:
	    for (k = kmin; k <= kmax; ++k)
	    for (j = jmin; j <= jmax; ++j)
	    for (i = imin; i <= imax; ++i)
	    {
	    	icoords[0] = i;
	    	icoords[1] = j;
	    	icoords[2] = k;
		pointExplicitCimSolver(icoords,sub_comp);
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

	stop_clock("computeAdvectionExplicitCim");
}	/* computeAdvectionExplicit */
void VCARTESIAN::computeAdvectionExplicit(COMPONENT sub_comp)
{
	int i,j,k,l,m,ic,icn,icoords[MAXD],nc;
	int gmin[MAXD],ipn[MAXD];
	int index0;
	double coords[MAXD],crx_coords[MAXD];
	double vapor,vapor_nb[2],dgrad[MAXD],grad_plus[MAXD],
	       grad_minus[MAXD];
	double coef;
	COMPONENT comp;
	boolean fr_crx_grid_seg;
	const GRID_DIRECTION dir[3][2] = 
		{{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
	double *Temp = field->vapor;
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
                array[ic] = vapor = Temp[ic];
		for (l = 0; l < dim; ++l)
		{
		    v[l] = 0.0;
		    v_plus[l] = 0.0;
		    v_minus[l] = 0.0;
		}
		if (sub_comp == LIQUID_COMP2 && vel != NULL)
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
                                icoords,dir[l][m],comp,getStateVapor,
                                &vapor_nb[m],crx_coords);
                        if (!fr_crx_grid_seg)
                        {
                             next_ip_in_dir(icoords,dir[l][m],ipn,gmin,
							top_gmax);
                             icn = d_index1d(ipn[0],top_gmax);
			     vapor_nb[m] = Temp[icn];
                        }
                        dgrad[l] += (vapor_nb[m] - vapor)/top_h[l];
                    }
		    grad_plus[l] = (vapor_nb[1] - vapor)/top_h[l];
		    grad_minus[l] = (vapor - vapor_nb[0])/top_h[l];
                    array[ic] += coef*dgrad[l]/top_h[l]-m_dt*
			(v_plus[l]*grad_minus[l]+ v_minus[l]*grad_plus[l]);
                }
            }
            break;
	case 2:
	    for (j = jmin; j <= jmax; ++j)
	    for (i = imin; i <= imax; ++i)
	    {
	    	icoords[0] = i;
	    	icoords[1] = j;
	    	ic = d_index2d(i,j,top_gmax);
	    	comp = top_comp[ic];
	    	if (comp != sub_comp) 
	    	    continue;
                array[ic] = vapor = Temp[ic];
		for (l = 0; l < dim; ++l)
                {
                    v[l] = 0.0;
                    v_plus[l] = 0.0;
                    v_minus[l] = 0.0;
                }
                if (sub_comp == LIQUID_COMP2 && vel != NULL)
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
                                icoords,dir[l][m],comp,getStateVapor,
                                &vapor_nb[m],crx_coords);
                        if (!fr_crx_grid_seg)
                        {
                            next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax);
                            icn = d_index2d(ipn[0],ipn[1],top_gmax);
			    vapor_nb[m] = Temp[icn];
                        }
                        dgrad[l] += (vapor_nb[m] - vapor)/top_h[l];
                    }
		    grad_plus[l] = (vapor_nb[1] - vapor)/top_h[l];
		    grad_minus[l] = (vapor - vapor_nb[0])/top_h[l];
		    array[ic] += coef*dgrad[l]/top_h[l] - m_dt*
			(v_plus[l]*grad_minus[l]+ v_minus[l]*grad_plus[l]);
		}
	    }
	    break;
	case 3:
	    for (k = kmin; k <= kmax; ++k)
	    for (j = jmin; j <= jmax; ++j)
	    for (i = imin; i <= imax; ++i)
	    {
	    	icoords[0] = i;
	    	icoords[1] = j;
	    	icoords[2] = k;
	    	ic = d_index3d(i,j,k,top_gmax);
	    	comp = top_comp[ic];
	    	if (comp != sub_comp) 
	    	    continue;
                array[ic] = vapor = Temp[ic];
		for (l = 0; l < dim; ++l)
                {
                    v[l] = 0.0;
                    v_plus[l] = 0.0;
                    v_minus[l] = 0.0;
                }
                if (sub_comp == LIQUID_COMP2 && vel != NULL)
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
                                icoords,dir[l][m],comp,getStateVapor,
                                &vapor_nb[m],crx_coords);
                        if (!fr_crx_grid_seg)
                        {
                            next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax);
                            icn = d_index3d(ipn[0],ipn[1],ipn[2],top_gmax);
                	    vapor_nb[m] = Temp[icn];
                        }
                        dgrad[l] += (vapor_nb[m] - vapor)/top_h[l];
                    }
		    grad_plus[l] = (vapor_nb[1] - vapor)/top_h[l];
		    grad_minus[l] = (vapor - vapor_nb[0])/top_h[l];
		    array[ic] += coef*dgrad[l]/top_h[l] - m_dt*
			(v_plus[l]*grad_minus[l]+ v_minus[l]*grad_plus[l]);
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

void VCARTESIAN::checkField()
{
	int i,j,k,index,count;
	printf("Enter checkField()\n");
	double prev = 0;

	count = 0;
	for(j = 0; j < comp_size; j++)
        {
                if(field->temperature[j] != prev)
		{
                    printf("%s[%d] = %20.14f\n",
                        "temperature",j,
                        field->temperature[j]);
		    prev = field->temperature[j];
		    count ++;
		}
		if(count > 20)
		{
		    printf("......\n");
		    break;
		}
        }

	count = 0;
        for(j = 0; j < comp_size; j++)
        {
		if(field->vapor[j] != prev)
                { 
		    printf("%s[%d] = %20.14f\n",
                        "vapor",j,
                        field->vapor[j]);
		    prev = field->vapor[j];
		    count ++;
		}
		if(count > 20)
                {
                    printf("......\n");
                    break;
                }
        }

	count = 0;
        for(j = 0; j < comp_size; j++)
        {
		if(field->supersat[j] != prev)
                {
		    printf("%s[%d] = %20.14f\n",
                        "supersat",j,
                        field->supersat[j]);
		    prev = field->supersat[j];
		    count++;
		}
		if(count > 20)
                {
                    printf("......\n");
                    break;
                }
        }
	
	int s_count = 0, l_count = 0, o_count = 0;
	if (dim == 2)
	{
	    for(j = jmin; j <= jmax; j++)
            {
                int i = (int)(imax/2);
		index = d_index2d(i,j,top_gmax);
                    printf("%s[%d][%d] = (%f,%f)\n",
                        "velo",i,j,
                        field->vel[0][index],field->vel[1][index]);
            }
            for(i = imin; i <= imax; i++)
            {
		j = (int)(jmax/2);
		index = d_index2d(i,j,top_gmax);
                    printf("%s[%d][%d] = %f\n",
                        "pres",i,j,
                        field->pres[index]);
            }
	    index = d_index2d(imin,jmin,top_gmax);
	    printf("pres[%d][%d] = %f\n",imin,jmin,field->pres[index]);

            for (j = jmin; j <= jmax; ++j)
            for (i = imin; i <= imax; ++i)
            {
                index = d_index2d(i,j,top_gmax);
	        if(top_comp[index] == SOLID_COMP)
	    	    s_count ++; 
	        if(top_comp[index] == LIQUID_COMP2)
		    l_count ++;
	        if(field->vapor[index] == 0)
		    o_count ++;
	    }
	}
	else if (dim == 3)
	{
	    for(k = kmin; k <= kmax; ++k)
            {
                i = (int)(imax/2);
		j = (int)(jmax/2); 
		index = d_index3d(i,j,k,top_gmax);
                    printf("%s[%d][%d][%d] = (%f,%f,%f)\n",
                        "velo",i,j,k,
                        field->vel[0][index],field->vel[1][index],field->vel[2][index]);
            }
            for(i = imin; i <= imax; i++)
            {
		j = (int)(jmax/2);
		k = (int)(kmax/2);
		index = d_index3d(i,j,k,top_gmax);
                    printf("%s[%d][%d][%d] = %f\n",
                        "pres",i,j,k,
                        field->pres[index]);
            }
	    index = d_index3d(imin,jmin,kmin,top_gmax);
	    printf("pres[%d][%d][%d] = %f\n",imin,jmin,kmin,field->pres[index]);


	    for (k = kmin; k <= kmax; ++k)
            for (j = jmin; j <= jmax; ++j)
            for (i = imin; i <= imax; ++i)
            {
                index = d_index3d(i,j,k,top_gmax);
                if(top_comp[index] == SOLID_COMP)
                    s_count ++;
                if(top_comp[index] == LIQUID_COMP2)
                    l_count ++;
                if(field->vapor[index] == 0)
                    o_count ++;
            }
	}
	printf("s_num=%d, l_num=%d, 0_num=%d\n",s_count,l_count,o_count);
}

void VCARTESIAN::recordField(char *outname, const char *varname)
{
	int i,j,k,index;
	FILE *outfile;
	char filename[256];
        sprintf(filename, "%s/record-%s",outname,varname);

        if (!create_directory(filename,NO))
        {
            printf("Cannot create directory %s\n",filename);
            clean_up(ERROR);
        }
        sprintf(filename,"%s/%s-%4.2f",filename,varname,front->time);
        if (pp_numnodes() > 1)
           return;
        outfile = fopen(filename,"w");
	switch(dim)
	{
	    case 2:
	        for (j = jmin; j <= jmax; ++j)
		for (i = imin; i <= imax; ++i)
        	{
		    index = d_index2d(i,j,top_gmax);
	            if(strcmp(varname,"temperature") == 0)
        	        fprintf(outfile,"%20.14f\n",field->temperature[index]);
	            else if(strcmp(varname,"vapor") == 0)
        	        fprintf(outfile,"%f\n",field->vapor[index]);
	            else if(strcmp(varname,"supersat") == 0)
        	        fprintf(outfile,"%f\n",field->supersat[index]);
		}
	        if(strcmp(varname,"velocity") == 0)
	  	{
	            for (j = jmin; j <= jmax; ++j)
		    for (i = imin; i <= imax; ++i)
		    {
			index = d_index2d(i,j,top_gmax);
			fprintf(outfile,"%f  %f\n",
				field->vel[0][index],field->vel[1][index]);
		    }
		}
		break;
	    case 3:
		for (k = kmin; k <= kmax; ++k)
                for (j = jmin; j <= jmax; ++j)
                for (i = imin; i <= imax; ++i)
		{
		    index = d_index3d(i,j,k,top_gmax);
	            if(strcmp(varname,"temperature") == 0)
        	        fprintf(outfile,"%20.14f\n",field->temperature[index]);
	            else if(strcmp(varname,"vapor") == 0)
        	        fprintf(outfile,"%f\n",field->vapor[index]);
	            else if(strcmp(varname,"supersat") == 0)
        	        fprintf(outfile,"%f\n",field->supersat[index]);
		}
	        if(strcmp(varname,"velocity") == 0)
	  	{
		    for (k = kmin; k <= kmax; ++k)
                    for (j = jmin; j <= jmax; ++j)
		    for (i = imin; i <= imax; ++i)
		    {
			index = d_index3d(i,j,k,top_gmax);
			fprintf(outfile,"%f  %f  %f\n",
				field->vel[0][index],field->vel[1][index],
				field->vel[2][index]);
		    }
		}
		break;
	    default:
		printf("Unknown dim = %d\n",dim);
		clean_up(ERROR);
	}
	fclose(outfile);
}

void VCARTESIAN::recordTKE()
{
        double **vel = field->vel;
        int i,j,k,count,index;
        static boolean first = YES;
	double E;
        char fname[256];
        FILE *file;
        sprintf(fname,"%s/TKE.m",front->out_name);
        if (pp_mynode() != 0)
           return;
        file = fopen(fname,"a");
	if (first)
	{
	    first = NO;
            fprintf(file,"%%turbulence kinetic energy vs. time\n");
            fprintf(file,"%%time  TKE\n");
	}
	
	E = 0.0;
	count = 0;
	switch(dim)
	{
	    case 2:
                for (j = jmin; j <= jmax; ++j)
                for (i = imin; i <= imax; ++i)
                {
                    index = d_index2d(i,j,top_gmax);
		    E += 0.5*vel[0][index]*vel[0][index];
		    E += 0.5*vel[1][index]*vel[1][index];
		    count++;	
		}
		E /= count;
		break;
	    case 3:
                for (k = kmin; k <= kmax; ++k)
                for (j = jmin; j <= jmax; ++j)
		for (i = imin; i <= imax; ++i)
                {
                    index = d_index3d(i,j,k,top_gmax);
		    E += 0.5*vel[0][index]*vel[0][index];
                    E += 0.5*vel[1][index]*vel[1][index];
		    E += 0.5*vel[2][index]*vel[2][index];
		    count++;
		}
		E /= count;
		break;
	    default:
		printf("Unknown dim = %d\n",dim);
                clean_up(ERROR);
	}
	fprintf(file,"%f  %20.19f\n",front->time,E);
	fclose(file);
}


void VCARTESIAN::recordWaterBalance()
{
        double *vapor = field->vapor;
	double vapor_mass,liquid_mass;
	double rho,r,a3 = 1;
        int i,j,k,index;
        static boolean first = YES;
	double water_balance[4];
	IF_PARAMS* iFparams = (IF_PARAMS*)front->extra1;
        PARAMS* eqn_params = (PARAMS*)front->extra2;
        PARTICLE* particle_array = eqn_params->particle_array;
        int num_drops = eqn_params->num_drops;
	double temp;
        char fname[256];
        FILE *file;

        sprintf(fname,"%s/water_balance.m",front->out_name);
        file = fopen(fname,"a");

	for(i = 0; i < dim; i++)
            a3 *= top_h[i];

        if (first == YES)
        {
            first = NO;
	    fprintf(file,"%%water balance vs. time\n");
	    fprintf(file,"%%time  vapor  liquid  total\n");
        }
	/*compute water mass in vapor*/
        if (dim == 2)
        {
	    vapor_mass = 0.0;
	    for(j = jmin; j <= jmax; j++)
	    for(i = imin; i <= imax; i++)
	    {
                index = d_index2d(i,j,top_gmax);
		rho = iFparams->rho2;
		vapor_mass += rho * a3 *vapor[index]; 
	    }
        }
        else if (dim == 3)
        {
	    vapor_mass = 0.0;
	    int count = 0;
	    for(k = kmin; k <= kmax; k++)
            for(j = jmin; j <= jmax; j++)
            for(i = imin; i <= imax; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
		rho = iFparams->rho2;
                vapor_mass += rho * a3 *vapor[index];
            }
        }
	water_balance[0] = front->time;
        water_balance[1] = 0.001 * vapor_mass;
	/*compute water mass in droplets*/
	liquid_mass = 0.0;
	if (eqn_params->prob_type == PARTICLE_TRACKING)
	{
            for (i = 0; i < num_drops; i++)
            {
	        r = particle_array[i].radius;
	        liquid_mass += 4.0/3.0*PI*r*r*r*particle_array[i].rho;
	    }
	} 
	water_balance[2] = liquid_mass;
	water_balance[3] = liquid_mass + 0.001 * vapor_mass;

#if defined(__MPI__)
	double local_water[3];
	for (i = 1; i < 4; i++)
	    local_water[i-1] = water_balance[i];

	if (pp_mynode() != 0)
	{
	    pp_send(1, (POINTER)local_water,sizeof(double)*4,0);
	}
	else
	{
	    for (i = 1; i < pp_numnodes(); i++)
	    {
		pp_recv(1, i, (POINTER)local_water,sizeof(double)*4);
		for (i = 1; i < 4; i++)
		    water_balance[i] += local_water[i-1];
	    }
	}
#endif /* defined(__MPI__) */
	for (i = 0 ; i < 4; i++)
	    fprintf(file,"%20.19f  ",water_balance[i]);
	fprintf(file,"\n");	
	fclose(file);
}

void VCARTESIAN::recordSampleRadius()
{
	int i,j;
        char fname[256];
        FILE *file;
	static boolean first = YES;

        sprintf(fname,"%s/sample_radius.m",front->out_name);
        file = fopen(fname,"a");
	if (first)
	{
            fprintf(file,"%%sample_radius vs. time\n");
            fprintf(file,"%%time radius\n");
	    first = NO;
	}
	PARAMS* eqn_params = (PARAMS*)front->extra2;
	PARTICLE* particle_array = eqn_params->particle_array;
        fprintf(file,"%f  %20.19f\n",front->time,
                         particle_array[0].radius);
	fclose(file);
}

void VCARTESIAN::setDomain()
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
            if (first == YES)
            {
		comp_size = top_gmax[0]+1;
                FT_VectorMemoryAlloc((POINTER*)&array,comp_size,FLOAT);
                FT_VectorMemoryAlloc((POINTER*)&source,comp_size,FLOAT);
                FT_VectorMemoryAlloc((POINTER*)&field->vapor,comp_size,
			FLOAT);
		FT_VectorMemoryAlloc((POINTER*)&field->supersat,comp_size,
			FLOAT);
                first = NO;
            }	
	    imin = (lbuf[0] == 0) ? 1 : lbuf[0];
	    imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
	    eqn_params->field = field;
	    break;
	case 2:
	    if (first == YES)
	    {
		comp_size = (top_gmax[0]+1)*(top_gmax[1]+1);
	    	FT_VectorMemoryAlloc((POINTER*)&array,comp_size,FLOAT);
                FT_VectorMemoryAlloc((POINTER*)&source,comp_size,FLOAT);
	    	FT_VectorMemoryAlloc((POINTER*)&field->vapor,
			comp_size,FLOAT);
                FT_VectorMemoryAlloc((POINTER*)&field->supersat,comp_size,
                        FLOAT);
                FT_VectorMemoryAlloc((POINTER*)&field->drops,comp_size,
                        FLOAT);
	    	first = NO;
	    }
	    imin = (lbuf[0] == 0) ? 1 : lbuf[0];
	    jmin = (lbuf[1] == 0) ? 1 : lbuf[1];
	    imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
	    jmax = (ubuf[1] == 0) ? top_gmax[1] - 1 : top_gmax[1] - ubuf[1];
	    eqn_params->field = field;
	    break;
	case 3:
	    if (first == YES)
	    {
		comp_size = (top_gmax[0]+1)*(top_gmax[1]+1)*(top_gmax[2]+1);
	    	FT_VectorMemoryAlloc((POINTER*)&array,comp_size,FLOAT);
                FT_VectorMemoryAlloc((POINTER*)&source,comp_size,FLOAT);
	    	FT_VectorMemoryAlloc((POINTER*)&field->vapor,
			comp_size,FLOAT);
                FT_VectorMemoryAlloc((POINTER*)&field->supersat,comp_size,
                        FLOAT);
                FT_VectorMemoryAlloc((POINTER*)&field->drops,comp_size,
                        FLOAT);
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

void VCARTESIAN::initMovieVariables()
{
	if (dim == 2)
	{	
 	    FT_AddHdfMovieVariable(front,NO,YES,SOLID_COMP,
                                "vapor",0,field->vapor,
				getStateVapor,0,0);
	    FT_AddHdfMovieVariable(front,NO,YES,SOLID_COMP,
                                "supersat",0,field->supersat,
				getStateSuper,0,0);
	}
        if (debugging("trace"))
            printf("Leaving initMovieVariables()\n");
}	/* end initMovieVariables */

void VCARTESIAN::augmentMovieVariables(const char* varname)
{
        int i;
        static HDF_MOVIE_VAR *hdf_movie_var;
        int n,offset,num_var;
        MOVIE_OPTION *movie_option = eqn_params->movie_option;

        offset = front->hdf_movie_var->num_var;
        num_var = offset+1;
	
        FT_ScalarMemoryAlloc((POINTER*)&hdf_movie_var,
                                sizeof(HDF_MOVIE_VAR));
        FT_MatrixMemoryAlloc((POINTER*)&hdf_movie_var->var_name,num_var,
                                        100,sizeof(char));
        FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->top_var,num_var,
                                        sizeof(double*));
        FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->idir,num_var,
                                        sizeof(int));
        FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->obstacle_comp,
                                        num_var,sizeof(COMPONENT));
        FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->preset_bound,num_var,
                                        sizeof(boolean));
        FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->var_min,num_var,
                                        sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->var_max,num_var,
                                        sizeof(double));
        for (i = 0; i < front->hdf_movie_var->num_var; ++i)
        {
            strcpy(hdf_movie_var->var_name[i],
                        front->hdf_movie_var->var_name[i]);
            hdf_movie_var->get_state_var[i] =
                        front->hdf_movie_var->get_state_var[i];
            hdf_movie_var->top_var[i] =
                        front->hdf_movie_var->top_var[i];
            hdf_movie_var->obstacle_comp[i] =
                        front->hdf_movie_var->obstacle_comp[i];
        }
        hdf_movie_var->num_var = n = offset;
        sprintf(hdf_movie_var->var_name[n],"%s",varname);
        hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
        if(strcmp(varname,"vapor") == 0)
	{
            hdf_movie_var->get_state_var[n] = getStateVapor;
	    hdf_movie_var->top_var[n] = eqn_params->field->vapor;
	}
	if(strcmp(varname,"supersat") == 0)
	{ 
	    hdf_movie_var->get_state_var[n] = getStateSuper;
	    hdf_movie_var->top_var[n] = eqn_params->field->supersat;
	}

        hdf_movie_var->num_var = ++n;	
	front->hdf_movie_var = hdf_movie_var;
}


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

void VCARTESIAN::recordRadius(char *out_name)
{
	INTERFACE *intfc = front->interf;
	PARAMS* eqn_params = (PARAMS*)front->extra2;
	char fname[100];
	FILE *file;
	double radius;
	CURVE **c;
	SURFACE **s;
	int i,index;
	
	sprintf(fname,"%s/record-radius",out_name);
        if (!create_directory(fname,NO))
        {
            printf("Cannot create directory %s\n",fname);
            clean_up(ERROR);
        }
        sprintf(fname,"%s/radius-%4.2f",fname,front->time);
#if defined(__MPI__)
        if (pp_numnodes() > 1)
            sprintf(fname,"%s-nd%s",fname,right_flush(pp_mynode(),4));
#endif /* defined(__MPI__) */
	
        file = fopen(fname,"w");

	if (eqn_params->prob_type == PARTICLE_TRACKING)
	{
	    for (i = 0; i < eqn_params->num_drops; i++)
                fprintf(file,"%16.15f\n",eqn_params->particle_array[i].radius);
	}
	fclose(file);
}

void VCARTESIAN::computeSource()
{
	/*compute condensation rate*/
	/*this is a source term for vapor mixing ratio equation*/
        int i, index, size, num_drops;
        int ic[MAXD];
        PARAMS* eqn_params = (PARAMS*)front->extra2;
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
        PARTICLE* particle_array = eqn_params->particle_array;
        num_drops = eqn_params->num_drops;
	double *supersat = eqn_params->field->supersat;
        static boolean first = YES;
        double *coords;
	/*for computing the coefficient*/
	double rho_0 = iFparams->rho2;
	double a3 = 1.0;
	double coeff;

	for(i = 0; i < dim; i++)
	    a3 *= top_h[i];

        for (i = 0; i < comp_size; i++)
        {
	    source[i] = 0.0;
	    field->drops[i] = 0;
	}
	
        for (i = 0; i < num_drops; i++)
        {
            coords = particle_array[i].center;
            rect_in_which(coords,ic,top_grid);
            index = d_index(ic,top_gmax,dim);
	    coeff = 4.0*PI*particle_array[i].rho * eqn_params->K
					       / (rho_0 * a3);
	    source[index] += -1000.0 * coeff * supersat[index]
				     * particle_array[i].radius;
	    field->drops[index] += 1;
        }
}


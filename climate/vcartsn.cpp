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
                                INCLUDE_BOUNDARIES,point,t,&hse,&hs,2))
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
	PARAMS* eqn_params = (PARAMS*)front->extra2;
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
            sl->vapor = sr->vapor = eqn_params->qv0;
        }

	// cell_center
	for (i = 0; i < cell_center.size(); i++)
	{
	    field->vapor[i] = 0;
	    field->supersat[i] = 0;
	}
	for (i = 0; i < cell_center.size(); i++)
	{
	    getRectangleCenter(i,coords);
	    c = top_comp[i];
	    if (c == LIQUID_COMP2)
	    {
		if (eqn_params->init_vapor_state == RAND_STATE)
	    	    field->vapor[i] = gauss_center_limit(
				      (POINTER)&gauss_params,seed);
		else if (eqn_params->init_vapor_state == CONST_STATE)
		    field->vapor[i] = eqn_params->qv0;
		else
		    getInitialState(&c,coords,field,i,dim,eqn_params);
	    }
	    else if (c == SOLID_COMP)
	    	field->vapor[i] = 0;
	    else
	    	field->vapor[i] = eqn_params->qv0;
	}
	FT_ParallelExchGridArrayBuffer(field->vapor,front,NULL);
}	/* end setInitialCondition */


void VCARTESIAN::setParallelVapor(void)
{
        FILE *infile;
        int i,j,id,k,l,index,G_index;
        char fname[100];
        COMPONENT comp;
        double coords[MAXD];
        int size = (int)cell_center.size();
        int myid = pp_mynode();
        int numprocs = pp_numnodes();

        int G_icoords[MAXD],pp_icoords[MAXD],icoords[MAXD];
        int local_gmax[MAXD], global_gmax[MAXD];
        int G_size, L_size;
        PP_GRID *pp_grid = front->pp_grid;
        double *local_L = pp_grid->Zoom_grid.L;
        double *local_U = pp_grid->Zoom_grid.U;
        double *GV_buff, *V_buff;
	
        for (i = 0; i < dim; i++)
        {
            global_gmax[i] = pp_grid->Global_grid.gmax[i]-1;
            local_gmax[i] = pp_grid->Zoom_grid.gmax[i]-1;
        }
        FT_MakeGridIntfc(front);
        setDomain();
        G_size = 1;
        L_size = 1;
        for (i = 0; i < dim; i++)
        {
            G_size = G_size * (global_gmax[i]+1);
            L_size = L_size * (top_gmax[i]+1);
        }
        uni_array(&V_buff,L_size,sizeof(double));
        if (myid == 0)
        {
            uni_array(&GV_buff,G_size,sizeof(double));
	    /*setInitialVapor(front,LIQUID_COMP,GV_buff);*/
	    /*Please set GV_buff before sending out*/
            for (id = 0; id < numprocs; id++)
            {
                find_Cartesian_coordinates(id,pp_grid,pp_icoords);
		switch (dim)
		{
		    case 2:
                	for (j = jmin; j <= jmax; ++j)
                	for (i = imin; i <= imax; ++i)
                	{
                    	    icoords[0] = i;
                    	    icoords[1] = j;
                    	    G_icoords[0] = pp_icoords[0]*(local_gmax[0]+1)+icoords[0]-imin;
                    	    G_icoords[1] = pp_icoords[1]*(local_gmax[1]+1)+icoords[1]-jmin;
                    	    G_index = d_index(G_icoords,global_gmax,dim);
                    	    index = d_index(icoords,top_gmax,dim);
                    	    V_buff[index] = GV_buff[G_index];
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
                    	    G_icoords[0] = pp_icoords[0]*(local_gmax[0]+1)+icoords[0]-imin;
                    	    G_icoords[1] = pp_icoords[1]*(local_gmax[1]+1)+icoords[1]-jmin;
                    	    G_icoords[2] = pp_icoords[2]*(local_gmax[2]+1)+icoords[2]-kmin;
                    	    G_index = d_index(G_icoords,global_gmax,dim);
                    	    index = d_index(icoords,top_gmax,dim);
                    	    V_buff[index] = GV_buff[G_index];
                	}
			break;
            	    Default:
                  	printf("Unknown dim = %d\n",dim);
                	clean_up(ERROR);
		}
                if (id == 0)
                {
                    for (i = 0; i < L_size; i++)
                        field->vapor[i] = V_buff[i];
                }
                else
                {
                    pp_send(1,(POINTER)(V_buff),sizeof(double)*L_size,id);
                }
            }
            FT_FreeThese(1,GV_buff);
        }
        else
        {
            pp_recv(1,0,(POINTER)(V_buff),sizeof(double)*L_size);
            for (i = 0; i < L_size; i++)
            {
                field->vapor[i] = V_buff[i];
            }
        }

        FT_FreeThese(1,V_buff);
        setAdvectionDt();
}

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
	printf("computeAdvectionCim() to be implemented\n");
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
	int i,j,k,size,index;
	double *temp  = eqn_params->field->temperature;
	double *vapor  = eqn_params->field->vapor;
	double sat_vap_pre, sat_vap_rat;
	double *super = eqn_params->field->supersat;

	/*compute coeffecient for condensation*/
        double Lh, Rv, Rd, rhoL, Kc, es, D, Cp, ksi;
        double T,p;
        double A;
	/*see instructions in climate.h*/
        Lh = eqn_params->Lh;
        Rv = eqn_params->Rv;
        Rd = eqn_params->Rd;
        rhoL = eqn_params->rho_l;
        Kc = eqn_params->Kc;
        D = eqn_params->D;
        ksi = Rd/Rv;

        T = 0;
        p = 0;
	for (index = 0; index < comp_size; index++)
        {
            T += field->temperature[index];
            p += field->pres[index];
        }
	size = comp_size;
#if defined(__MPI)
	pp_gsync();
        pp_global_sum(&T,1);
        pp_global_sum(&p,1);
        pp_global_isum(&size,1);
#endif
        T /= size;
        p /= size;

        es = 611.2*exp(17.67*(T-273.15)/(T-29.65));
        eqn_params->K = 1/((Lh/(Rv*T)-1)*Lh*rhoL/(Kc*T)+rhoL*Rv*T/(D*es));
	printf("Condensation coeffecient = %e\n",eqn_params->K);

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
		printf("sat_vap_rat = %f\n",sat_vap_rat); 
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
	FT_ParallelExchGridArrayBuffer(array,front,NULL);
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

        sprintf(dirname,"%s/vtk",OutName(front));
        if (pp_mynode() == 0)
        {
            if (!create_directory(dirname,YES))
            {
                screen("Cannot create directory %s\n",dirname);
                clean_up(ERROR);
            }
        }
        pp_gsync();
        sprintf(dirname,"%s/vtk.ts%s",dirname,right_flush(front->step,7));
        if (pp_numnodes() > 1)
            sprintf(dirname,"%s-nd%s",dirname,right_flush(pp_mynode(),4));
        if (!create_directory(dirname,YES))
        {
            screen("Cannot create directory %s\n",dirname);
            clean_up(ERROR);
        }

        //cell-based liquid phase
        ph_index.clear();

        sprintf(filename,"%s/%s.vtk",dirname,varname);
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

void VCARTESIAN::recordClusteringIndex()
{
	double Mean, Var, CL, sigma, CSL;
	static boolean first = YES;
	FILE *file;
	char fname[100];
	char *out_name = front->out_name;
	double *array = field->drops;
	int size = comp_size;
	sprintf(fname,"%s/cluster",out_name);
	if (first)
	{
	    if(pp_mynode() == 0)
	    {
	        file = fopen(fname,"w");
	        fprintf(file,"%%Clustering Index\n");
	        fprintf(file,"%%CL	sigma	CSL = CL/sigma\n");
	        fclose(file);
	    }
	    first = NO;
	}
	Deviation(array,size,Mean,Var);
	CL = Var/Mean - 1.0;
	printf("Mean = %f Variance = %f CL = %f\n",Mean,Var,CL);
	sigma = sqrt(2.0/(size));
	CSL = CL/sigma;
	if(pp_mynode() == 0)
	{
	    file = fopen(fname,"a");
	    fprintf(file,"%f    %e    %f\n",CL,sigma,CSL);
	    fclose(file);
	}
	return;
}

double VCARTESIAN::computeReactTimeScale(double R0, double S0, int N)
{
	PARAMS *eqn_params = (PARAMS*)front->extra2;
        IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	
	int i,j,k,size,index;
	double R2, S, t, dt = 0.001;
	double Lh, Rv, Rd, rhoL, K, es, D, Cp, ksi;
	double T, p;  
	double A, B;
	Lh = eqn_params->Lh;
	Rv = eqn_params->Rv;
	Rd = eqn_params->Rd;
	rhoL = eqn_params->rho_l;
	K = eqn_params->Kc;
	D = eqn_params->D;
	Cp = eqn_params->Cp;
	ksi = Rd/Rv;

	T = 0;
	p = 0;
	size = 0;
        for (k = kmin; k < kmax; k++)
        for (j = jmin; j < jmax; j++)
        for (i = imin; i < imax; i++)
        {
	    size ++;
	    index = d_index3d(i,j,k,top_gmax);
	    T += field->temperature[index];
	    p += field->pres[index];
	}
#if defined(__MPI)
	pp_gsync();
	pp_global_sum(&T,1);
	pp_global_sum(&p,1);
	pp_global_isum(&size,1);
#endif
	T /= size;
	p /= size;

	es = 611.2*exp(17.67*(T-273.15)/(T-29.65));
	A = 1/((Lh/(Rv*T)-1)*Lh*rhoL/(K*T)+rhoL*Rv*T/(D*es));
	B = 4*PI*N*rhoL*(Rd*T/(ksi*es)+ksi*sqr(Lh)/(p*T*Cp))
	    /((Lh/(Rv*T)-1)*Lh*rhoL/(K*T)+rhoL*Rv*T/(D*es));
	t = 0;
	R2 = sqr(R0);
	S = S0;
	printf("A = %e, B = %e\n",A,B);
	while(R2 > 0 && S < -0.005)
	{
	    R2 += 2*A*S*dt;
	    S  /= 1+B*dt*sqrt(R2);
	    t  += dt;
	}
	return t;
}


void VCARTESIAN::recordMixingLine()
{
	PARAMS *eqn_params = (PARAMS*)front->extra2;
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	PARTICLE *particle_array = eqn_params->particle_array;
	int num_drops = eqn_params->num_drops;
	static boolean first = YES;
	static double prev_RN[2] = {0,0};
	double mu = iFparams->mu2, eta;
	double **vel = field->vel;
	double slope;
	double t_evap, t_mix, t_react;
	double Dev; /*Deviation of radius distribution*/
	double alpha; /**/
	double DspRat,rv,r0,rm; /*dissipation rate*/
	double Sij[MAXD][MAXD];
	double NL; /*transition scale number*/
	double avg_supersat;  /*averaged supersaturation in clear air*/
	int i, j, k, l, m, m1, m2, ic, I, icoords[MAXD],size, count, nzeros;
	int gmin[MAXD], ipn[MAXD];
	int I_nb[3][2];
	const GRID_DIRECTION dir[3][2] =
                {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
	FILE *file;
	char fname[100];
	char *out_name = front->out_name;
	int  myid = pp_mynode();
	double L = front->pp_grid->Global_grid.U[0] 
		 - front->pp_grid->Global_grid.L[0];
	double ReduceBuff[8];

	static int max_array_size = 0;
	static double* radius_array;

	if (debugging("trace"))
	    printf("Entering record mixing line\n");
        if (eqn_params->num_drops > max_array_size)
        {
            max_array_size = eqn_params->num_drops;
            free_these(1,radius_array);
            FT_VectorMemoryAlloc((POINTER*)&radius_array,max_array_size,FLOAT);
        }

	if (3 != dim)
	    return;
	double a3 = 1.0;
	for (i = 0; i < dim; i++) 
	{
	    gmin[i] = 0;
	    a3 *= top_h[i]; 
	}
	/*compute mean energy dissipation rate*/
	DspRat = 0.0;
	avg_supersat = 0.0;
	count = 0;
	size = 0;
        for (k = kmin; k < kmax; k++)
        for (j = jmin; j < jmax; j++)
        for (i = imin; i < imax; i++)
        {
	    size++;
            ic = d_index3d(i,j,k,top_gmax);
	    if (field->supersat[ic] < 0)
	    {
   	        avg_supersat += field->supersat[ic];
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
                    I_nb[m][l] = d_index3d(ipn[0],ipn[1],ipn[2],top_gmax);
		else
		{
		    printf("In recordMixingLine(), cannot find next ip\n");
	            printf("gmin=[%d %d %d]\n",gmin[0],gmin[1],gmin[2]);
		    printf("gmax=[%d %d %d]\n",top_gmax[0],top_gmax[1],top_gmax[2]);
		    printf("icoords=[%d %d %d]\n",icoords[0],icoords[1],icoords[2]);
		    clean_up(ERROR);
		}
            }
	    /*compute strain rate tensor Sij*/
	    for(l = 0; l < dim; ++l)
	    for(m = 0; m < dim; ++m)
	    {
		Sij[l][m]  = 0.5*(vel[l][I_nb[m][1]]-vel[l][I_nb[m][0]])/top_h[m];
		Sij[l][m] += 0.5*(vel[m][I_nb[l][1]]-vel[m][I_nb[l][0]])/top_h[l];
	    }
	    /*compute dissipation rate using: 0.5*mu*Sij^2*/ 
            for(l = 0; l < dim; ++l)
            for(m = 0; m < dim; ++m)
	    {
		DspRat += 0.5*mu*sqr(Sij[l][m]);
	    }
	}
	
	/*compute Rv, t_evap*/
	rv = 0;
	r0 = 0;
	nzeros = 0;
	for (i = 0; i < eqn_params->num_drops; i++)
	{   
	    radius_array[i] = particle_array[i].radius;
	    if (particle_array[i].radius != 0)
		nzeros ++;
	    rv += pow(particle_array[i].radius,3.0);
	    r0 += pow(particle_array[i].R0,3.0);
	}

#if defined(__MPI__)
	ReduceBuff[0] = rv;
	ReduceBuff[1] = r0;
	ReduceBuff[2] = nzeros;
	ReduceBuff[3] = num_drops;
	ReduceBuff[4] = DspRat;
	ReduceBuff[5] = avg_supersat;
	ReduceBuff[6] = size;
	ReduceBuff[7] = count;
	pp_gsync();
	pp_global_sum(ReduceBuff,8);
	rv = ReduceBuff[0];
	r0 = ReduceBuff[1];
	nzeros    = (int)ReduceBuff[2];
	num_drops = (int)ReduceBuff[3];
	DspRat       = ReduceBuff[4];
	avg_supersat = ReduceBuff[5];
	size         = (int)ReduceBuff[6];
	count        = (int)ReduceBuff[7];
#endif 
	rv /= nzeros;
	r0 /= num_drops;
	DspRat   /= size;
	avg_supersat /= count;
	/*compute t_mix*/
	t_mix = pow(L*L/DspRat,1.0/3.0);
	printf("avg_supersat = %f, count = %d\n",avg_supersat,count);

	/*compute mean value and deviation of radius distribution*/
	Deviation(radius_array,eqn_params->num_drops,rm,Dev);
	if (avg_supersat < 0)
	{
	    t_evap  = sqr(rm)/(-eqn_params->K*avg_supersat);
	    t_react = computeReactTimeScale(rm,avg_supersat,nzeros); 
	}
	else
	{
	    t_evap = HUGE;
	    t_react = 0;
	}
	if (myid != 0) /*Post analysis is only done by the master processor*/
	    return; /*return if processor id is slaver processor*/
	slope = (double(nzeros)/double(num_drops)-prev_RN[1])
		/((rv/r0)-prev_RN[0]);
	prev_RN[0] = rv/r0;
	prev_RN[1] = double(nzeros)/double(num_drops);

	if (first)
	{
	    sprintf(fname,"%s/mixing",out_name);
	    file = fopen(fname,"w");
	    fprintf(file,"%%Damkoehler number VS slope of R-N diagram\n");
	    fprintf(file,"%%t_mix    t_evap    supersat    Damkoehler    slope\n");
	    fclose(file);

            sprintf(fname,"%s/RN",out_name);
            file = fopen(fname,"w");
            fprintf(file,"%%Rv    Rm    stdDev    Rv^3/R0^3    N/N0\n");
            fclose(file);

            sprintf(fname,"%s/transition",out_name);
            file = fopen(fname,"w");
            fprintf(file,"%%time    epsilon    eta    t_react    NL\n");
            fclose(file);
	    first = NO;
	    return;
	}
	sprintf(fname,"%s/mixing",out_name);
	file = fopen(fname,"a");
	fprintf(file,"%15.14f  %15.14f  %15.14f  %15.14f  %15.14f\n",
		t_mix,t_evap,avg_supersat,t_mix/t_evap,slope);
	fclose(file);
	/*plot mixing line: N/N0 VS. (Rv/Rv0)^3*/
        sprintf(fname,"%s/RN",out_name);
        file = fopen(fname,"a");
	fprintf(file,"%15.14f  %15.14f  %15.14f  %15.14f  %15.14f\n", 
		pow(rv,1.0/3.0),rm,sqrt(Dev),rv/r0, double(nzeros)/double(num_drops));
	fclose(file);
        sprintf(fname,"%s/transition",out_name);
        file = fopen(fname,"a");
	eta = pow(pow(mu,3.0)/DspRat,0.25);
	NL = pow(DspRat,0.5)*pow(t_react,1.5)/eta;
	fprintf(file,"%f  %15.14f  %15.14f  %15.14f  %15.14f\n",
		front->time,DspRat,eta,t_react,NL);
	fclose(file);
	if (debugging("trace"))
	    printf("Leaving record mixing line\n");
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
	
	count = 0;
        for(j = 0; j < comp_size; j++)
        {
		if(field->drops[j] != prev)
                {
		    printf("%s[%d] = %20.14f\n",
                        "drops",j,
                        field->drops[j]);
		    prev = field->drops[j];
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
}


void VCARTESIAN::recordField(char *outname, const char *varname)
{
	int i, j, k, index;
	FILE* outfile;
	char filename[256];
	double **vel = field->vel;
	
        sprintf(filename, "%s/record-%s",outname,varname);
        if (!create_directory(filename,NO))
        {
            printf("Cannot create directory %s\n",filename);
            clean_up(ERROR);
        }
	pp_gsync();

	if (pp_numnodes() > 1)
	{
	    sprintf(filename, "%s/record-%s-nd%03d",filename,varname,pp_mynode());
	    create_directory(filename,YES);
	}
        sprintf(filename,"%s/%s-%4.2f",filename,varname,front->time);
        outfile = fopen(filename,"w");
        switch (dim)
        {
            case 2:
                for (j = jmin; j <= jmax; j++)
                for (i = imin; i <= imax; i++)
                {
                    index = d_index2d(i,j,top_gmax);
		    if(strcmp(varname,"vapor") == 0)
                        fprintf(outfile,"%f\n",field->vapor[index]);
		    else if(strcmp(varname,"supersat") == 0)
			fprintf(outfile,"%f\n",field->supersat[index]);
		    else if(strcmp(varname,"temperature") == 0)
			fprintf(outfile,"%f\n",field->temperature[index]);
		    else if(strcmp(varname,"mrad") == 0)
			fprintf(outfile,"%15.14f\n",field->mrad[index]);
		    else if(strcmp(varname,"velocity") == 0)
			fprintf(outfile,"%f  %f\n",vel[0][index],vel[1][index]);
		    else
			printf("WARNING: Unknown field: %s\n",varname);
			
                }
                break;

                case 3:
                for (k = kmin; k <= kmax; k++)
                for (j = jmin; j <= jmax; j++)
                for (i = imin; i <= imax; i++)
                {
		    index = d_index3d(i,j,k,top_gmax);
                    if(strcmp(varname,"vapor") == 0)
                        fprintf(outfile,"%f\n",field->vapor[index]);
                    else if(strcmp(varname,"supersat") == 0)
                        fprintf(outfile,"%f\n",field->supersat[index]);
                    else if(strcmp(varname,"temperature") == 0)
                        fprintf(outfile,"%f\n",field->temperature[index]);
                    else if(strcmp(varname,"mrad") == 0)
                        fprintf(outfile,"%15.14f\n",field->mrad[index]);
                    else if(strcmp(varname,"velocity") == 0)
                        fprintf(outfile,"%f  %f  %f\n",
                            vel[0][index],vel[1][index],vel[2][index]);
                    else
                        printf("WARNING: Unknown field: %s\n",varname);
                }
                break;
        }
	fclose(outfile);
	return;
}

void VCARTESIAN::recordPDF(char *outname, const char *varname)
{
	int i,j,k,index;
	FILE *outfile;
	FILE *superfile;
	char supername[256];
	char filename[256];
	double *PDF;
	double bin_size,mid_bin;
	double var_min,var_max;
	double mean_super,var_super;
	double **vel = field->vel;
	int num_bins = 100;

	if(strcmp(varname,"all") == 0)
	{
	    recordPDF(outname,"vapor");
	    recordPDF(outname,"supersat");
	    recordPDF(outname,"temperature");
	    recordPDF(outname,"velocity");
	    recordPDF(outname,"numdensity");
	    recordPDF(outname,"cloud");
	    return;
	}

	if (strcmp(varname,"velocity") == 0)
	{
	    recordPDF(outname,"xvel");
	    recordPDF(outname,"yvel");
	    if (dim == 3)
		recordPDF(outname,"zvel");	
	    return;
	}
	if (pp_mynode() == 0)
	{
            sprintf(filename, "%s/PDF-%s",outname,varname);

            if (!create_directory(filename,NO))
            {
                printf("Cannot create directory %s\n",filename);
                clean_up(ERROR);
            }
            sprintf(filename,"%s/%s-%4.2f",filename,varname,front->time);
            outfile = fopen(filename,"w");
	}
        if(strcmp(varname,"temperature") == 0)
	{
	    PDF = ComputePDF(field->temperature,comp_size,bin_size,num_bins,var_min,var_max);
	}
        else if(strcmp(varname,"vapor") == 0)
	{
	    PDF = ComputePDF(field->vapor,comp_size,bin_size,num_bins,var_min,var_max);
	}
        else if(strcmp(varname,"supersat") == 0)
	{
	    PDF = ComputePDF(field->supersat,comp_size,bin_size,num_bins,var_min,var_max);
	    Deviation(field->supersat,comp_size,mean_super,var_super);
	    if (pp_mynode() == 0)
	    {
	        sprintf(supername, "%s/supersat",outname);
	        superfile = fopen(supername,"a");
		fprintf(superfile,"%f  %15.14f  %15.14f\n",front->time,mean_super,var_super);
		fclose(superfile);
	    }
	}
	else if (strcmp(varname,"numdensity") == 0)
	{
	    recordClusteringIndex();
	    PDF = ComputePDF(field->drops,comp_size,bin_size,num_bins,var_min,var_max);
	}
	 else if(strcmp(varname,"cloud") == 0)
	{
	    PDF = ComputePDF(field->cloud,comp_size,bin_size,num_bins,var_min,var_max);
	}
        else if(strcmp(varname,"xvel") == 0)
	{
	    PDF = ComputePDF(field->vel[0],comp_size,bin_size,num_bins,var_min,var_max);
	}
        else if(strcmp(varname,"yvel") == 0)
	{
	    PDF = ComputePDF(field->vel[1],comp_size,bin_size,num_bins,var_min,var_max);
	}
        else if(strcmp(varname,"zvel") == 0)
	{
	    PDF = ComputePDF(field->vel[2],comp_size,bin_size,num_bins,var_min,var_max);
	}
	else
	{
	   printf("WARNING: Unknown field: %s\n",varname);
	   return; 
	}
	if (pp_mynode() == 0)
	{
	    for (i = 0; i < num_bins;i++)
	    {
	        mid_bin = var_min + (0.5+i)*bin_size;
	        fprintf(outfile,"%f  %f\n",mid_bin,PDF[i]);
	    }
	    fclose(outfile);
	}
}

void VCARTESIAN::recordTKE()
{
        double **vel = field->vel;
        int i,j,k,count,index;
        static boolean first = YES;
	double E;
        char fname[256];
        FILE *file;
        sprintf(fname,"%s/TKE",front->out_name);
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
	double vapor_mass,liquid_mass,alpha;
	double rho,r,a3 = 1.0;
        int i,j,k,index,nzeros;
        static boolean first = YES;
	static double liquid_mass0 = 0.0;
	double water_balance[4];
	IF_PARAMS* iFparams = (IF_PARAMS*)front->extra1;
        PARAMS* eqn_params = (PARAMS*)front->extra2;
        PARTICLE* particle_array = eqn_params->particle_array;
        int num_drops = eqn_params->num_drops;
	double temp;
        char fname[256];
        FILE *file;
	double ReduceBuff[7];

	pp_gsync();
	if (pp_mynode() == 0)
	{
            sprintf(fname,"%s/water_balance",front->out_name);
            file = fopen(fname,"a");
	}

	for(i = 0; i < dim; i++)
            a3 *= top_h[i];

        if (first == YES && pp_mynode() == 0)
        {
	    fprintf(file,"%%water balance vs. time\n");
	    fprintf(file,"%%time  vapor  liquid  total  alpha\n");
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
	nzeros = 0;
	if (eqn_params->prob_type == PARTICLE_TRACKING)
	{
            for (i = 0; i < num_drops; i++)
            {
	        r = particle_array[i].radius;
		if (0 != r)
		    nzeros ++;
	        liquid_mass += 4.0/3.0*PI*r*r*r*particle_array[i].rho;
	    }
	    if (first) 
	    {
		liquid_mass0 = liquid_mass;
		first = NO; /*IMPORTANT: preserve initial LWC*/
		if (pp_mynode() == 0)
		    fclose(file);
		return;
	    }
	} 
	water_balance[2] = liquid_mass;
	water_balance[3] = liquid_mass + 0.001 * vapor_mass;

#if defined(__MPI__)
	for (i = 0; i < 3; i++)	
	    ReduceBuff[i] = water_balance[i+1];
	ReduceBuff[3] = nzeros;
	ReduceBuff[4] = num_drops;
	ReduceBuff[5] = liquid_mass;
	ReduceBuff[6] = liquid_mass0;
	pp_gsync();
	pp_global_sum(ReduceBuff,7);
	for (i = 0; i < 3; i++)
	    water_balance[i+1] = ReduceBuff[i];
	nzeros = (int)ReduceBuff[3];
	num_drops = (int)ReduceBuff[4];
	liquid_mass = ReduceBuff[5];
	liquid_mass0 = ReduceBuff[6];
#endif 

	if (!eqn_params->no_droplets && nzeros == 0)
	{
	     printf("No droplets included\n");
	     front->time_limit_reached = YES;
	}
	if(liquid_mass != liquid_mass0)
	{
	    alpha = (log(double(nzeros)/double(num_drops)))
		    /(log(liquid_mass/liquid_mass0));
	}
	else
	{
	    alpha = 0.0;
	}
	if (pp_mynode() == 0)
	{
	    for (i = 0 ; i < 4; i++)
	        fprintf(file,"%20.19f  ",water_balance[i]);
	
	    fprintf(file,"%20.19f\n",alpha);	
	    fclose(file);
	}
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

	if (field == NULL)
	    FT_ScalarMemoryAlloc((POINTER*)&field,sizeof(PHASE_FIELD));

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
		FT_VectorMemoryAlloc((POINTER*)&field->mrad,comp_size,
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
	    	FT_VectorMemoryAlloc((POINTER*)&field->cloud,
			comp_size,FLOAT);
                FT_VectorMemoryAlloc((POINTER*)&field->supersat,comp_size,
                        FLOAT);
                FT_VectorMemoryAlloc((POINTER*)&field->mrad,comp_size,
                        FLOAT);
                FT_VectorMemoryAlloc((POINTER*)&field->drops,comp_size,
                        FLOAT);
                FT_MatrixMemoryAlloc((POINTER*)&field->ext_accel,2,comp_size,
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
	    	FT_VectorMemoryAlloc((POINTER*)&field->cloud,
			comp_size,FLOAT);
                FT_VectorMemoryAlloc((POINTER*)&field->supersat,comp_size,
                        FLOAT);
                FT_VectorMemoryAlloc((POINTER*)&field->mrad,comp_size,
                        FLOAT);
                FT_VectorMemoryAlloc((POINTER*)&field->drops,comp_size,
                        FLOAT);
                FT_MatrixMemoryAlloc((POINTER*)&field->ext_accel,3,comp_size,
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
	PARAMS *params = (PARAMS*)front->extra2;
	if (dim == 2)
	{	
 	    FT_AddHdfMovieVariable(front,NO,YES,SOLID_COMP,
                                "vapor",0,field->vapor,
				getStateVapor,0,0);
	    FT_AddHdfMovieVariable(front,NO,YES,SOLID_COMP,
                                "supersat",0,field->supersat,
				getStateSuper,0,0);
	    FT_AddHdfMovieVariable(front,NO,YES,SOLID_COMP,
                                "cloud",0,field->cloud,
				getStateSuper,0,0);
            if (params->movie_option->plot_particles == YES)
                FT_AddVtkScalarMovieVariable(front,"Cloud",field->cloud);

	}
	else
	{
	    /* Added for vtk movie of scalar field */
            if (params->movie_option->plot_vapor == YES)
            	FT_AddVtkScalarMovieVariable(front,"Vapor",field->vapor);
            if (params->movie_option->plot_particles == YES)
                FT_AddVtkScalarMovieVariable(front,"Cloud",field->cloud);
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

void VCARTESIAN::recordRadius(char *out_name)
{
	INTERFACE *intfc = front->interf;
	PARAMS* eqn_params = (PARAMS*)front->extra2;
	char fname[100];
	FILE *file;
	double *PDF;
	double min_radius,max_radius;
	double bin_size, bin_mid;
	int i,bin_num;
	static double *radius_array;
	static int max_array_size = 0;
	
 	if (debugging("trace"))	
	    printf("Entering record radius\n");
	if (pp_mynode() == 0)
	{
	    sprintf(fname,"%s/record-radius",out_name);
            if (!create_directory(fname,NO))
            {
                printf("Cannot create directory %s\n",fname);
                clean_up(ERROR);
            }
            sprintf(fname,"%s/radius-%4.2f",fname,front->time);
            file = fopen(fname,"w");
	}
	if (eqn_params->num_drops > max_array_size)
	{
	    max_array_size = eqn_params->num_drops;
	    free_these(1,radius_array);
	    FT_VectorMemoryAlloc((POINTER*)&radius_array,max_array_size,FLOAT);
	}
	for (i = 0; i < eqn_params->num_drops; i++)
	    radius_array[i] = eqn_params->particle_array[i].radius;
	bin_num = 100;
	PDF = ComputePDF(radius_array,eqn_params->num_drops,bin_size,bin_num,min_radius,max_radius);
	printf("max_radius = %e, min_radius = %e, %d drops contained, %d bins used\n",
		max_radius, min_radius, eqn_params->num_drops,bin_num);
	if (pp_mynode() == 0)
	{
	    for (i = 0; i < bin_num; i ++)
	    {
	        bin_mid = min_radius+(0.5+i)*bin_size;
	        fprintf(file,"%15.14f  %15.14f\n",bin_mid,PDF[i]);
	    }
	    fclose(file);
	}
	if (debugging("trace"))
	    printf("Leaving record radius\n");
}

void VCARTESIAN::recordSampleParticles()
{
	int i, j, index;
	int ic[MAXD];
	PARAMS* eqn_params = (PARAMS*)front->extra2;
        PARTICLE* particle_array = eqn_params->particle_array;
	double *coords, Supersat, W, Radius;
	static boolean first = YES;
	FILE *file;
	char fname[100];

	/*For random choice samples*/
	UNIFORM_PARAMS uniform_params;
	unsigned short int xsubi[3];
	int SampleNum = 2000, num_drops = eqn_params->num_drops;
	static int *SampleIndex;
	double x;

	if (pp_numnodes() > 1)
	    return;
	
	uniform_params.a = 0.0;
	uniform_params.b = num_drops;
	xsubi[0] = 10;
	xsubi[1] = 100;
	xsubi[2] = 1000;
	if (first)
	{
	    first = NO;
	    FT_VectorMemoryAlloc((POINTER*)&SampleIndex,SampleNum,sizeof(int));
	    for (i = 0; i < SampleNum; i++)
	        SampleIndex[i] = 0;
	    if (pp_mynode() == 0)
	    for (i = 0; i < SampleNum; i++)
	    {
		x = dist_uniform((POINTER)&uniform_params,xsubi);
		SampleIndex[i] = round(x);
	    }
	    pp_global_imax(SampleIndex,pp_numnodes());
	}
	if (debugging("trace"))
            printf("Entering record sample radius\n");
	
	if (pp_mynode() == 0)
        {
            sprintf(fname,"%s/record-sample",OutName(front));
            if (!create_directory(fname,NO))
            {
                printf("Cannot create directory %s\n",fname);
                clean_up(ERROR);
            }
            sprintf(fname,"%s/particle-%f",fname,front->time);
        }
        file = fopen(fname,"w");
	for (j = 0; j < SampleNum; j++)
        {
	    i = SampleIndex[j];
            coords = particle_array[i].center;
            rect_in_which(coords,ic,top_grid);
            index = d_index(ic,top_gmax,dim);
	    W = field->vel[dim-1][index];
	    Supersat = field->supersat[index];
	    Radius = particle_array[i].radius;
	    fprintf(file,"%d  %f  %15.14f  %f  %f\n",
		    particle_array[i].Gindex,front->time,Radius,Supersat,W);
	}
	fclose(file);
}

void VCARTESIAN::initPresetParticles()
{
	int i,j,k,l,id;
	int index,count,G_count,num_drops;
	GAUSS_PARAMS gauss_params;
	UNIFORM_PARAMS uniform_params;
	double r_bar, sigma, x;
	unsigned short int xsubi[3];
	double cL[MAXD]; /*bound of a grid cell*/
	char msg[200];
	char *inname = InName(front);
	FILE *infile = fopen(inname,"r");
	double *supersat = field->supersat;
	PARAMS *eqn_params  = (PARAMS*)front->extra2;
	PARTICLE* particle_array = eqn_params->particle_array;
	if (particle_array != NULL)
	    FT_FreeThese(1,particle_array);

	CursorAfterString(infile,"Enter number of water drops:");
        fscanf(infile,"%d",&num_drops);
        (void) printf("%d\n",num_drops);
        eqn_params->num_drops = num_drops;

        sprintf(msg,"Enter mean radius of water drop:");
        CursorAfterString(infile,msg);
        fscanf(infile,"%lf",&r_bar);
        (void) printf("%f\n",r_bar);
        sprintf(msg,"Enter standard deviation of radius:");
        CursorAfterString(infile,msg);
        fscanf(infile,"%lf",&sigma);
        (void) printf("%f\n",sigma);	
	fclose(infile);
       	xsubi[0] = 10;
       	xsubi[1] = 100;
       	xsubi[2] = 1000;

       	gauss_params.mu = r_bar;
       	gauss_params.sigma = sigma;
       	uniform_params.a = 0.0;
       	uniform_params.b = 1.0;

	for (i = 0; i < comp_size; i++)
		field->drops[i] = 0;
	count = 0;
	switch (dim)
	{
	    case 2:
            	for (j = jmin; j <= jmax; j++)
		for (i = imin; i <= imax; i++)
		{
		    index = d_index2d(i,j,top_gmax); 
		    if(supersat[index] >= 0)
			count ++;	
		}
		G_count = count;

#if defined(__MPI__)
		pp_gsync();
		pp_global_isum(&G_count,1);
#endif
		printf("Cloud number = %d\n",G_count);
		if (G_count != 0)
		    num_drops = eqn_params->num_drops/G_count;/*num_drops in a cell*/
		else
		{
		    printf("No droplets in the domain!\n");
		    clean_up(ERROR);
		}
		eqn_params->num_drops = num_drops*count;
		FT_VectorMemoryAlloc((POINTER*)&particle_array,
                                      eqn_params->num_drops,sizeof(PARTICLE));
		count = 0;
                for (j = jmin; j <= jmax; j++)
                for (i = imin; i <= imax; i++)
                {
                    index = d_index2d(i,j,top_gmax);
                    if (supersat[index] < 0)
                        continue;
                    cL[0] = top_L[0] + i*top_h[0];
                    cL[1] = top_L[1] + j*top_h[1];
		    for (l = 0; l < num_drops; l++)
		    {
			for (id = 0; id < dim; ++id)
			{
			    x = dist_uniform((POINTER)&uniform_params,xsubi);
			    particle_array[l+count].center[id] = cL[id]+x*top_h[id];
			    particle_array[l+count].vel[id] = 0.0;
			}
			particle_array[l+count].radius = gauss_center_limit((POINTER)&gauss_params,xsubi);
			particle_array[l+count].R0 = particle_array[l].radius;
			particle_array[l+count].flag = YES;
			particle_array[l+count].rho = eqn_params->rho_l;
		    }
		    count += num_drops;
		    field->drops[index] = num_drops;
                }
	    	break;
	    case 3:
	    	for (k = kmin; k <= kmax; k++)
            	for (j = jmin; j <= jmax; j++)
		for (i = imin; i <= imax; i++)
		{
		    index = d_index3d(i,j,k,top_gmax);
		    if (supersat[index] >= 0)
			count ++;
		}
		G_count = count;
#if defined(__MPI__)
		pp_gsync();
		pp_global_isum(&G_count,1);
#endif
		printf("Cloud number = %d\n",G_count);
		num_drops = eqn_params->num_drops/G_count;/*num_drops in a cell*/
		if(eqn_params->num_drops < G_count)
			num_drops = 1;
		eqn_params->num_drops = num_drops*count;
		FT_VectorMemoryAlloc((POINTER*)&particle_array,
                                      eqn_params->num_drops,sizeof(PARTICLE));
		count = 0;
                for (k = kmin; k <= kmax; k++)
                for (j = jmin; j <= jmax; j++)
                for (i = imin; i <= imax; i++)
                {
                    index = d_index3d(i,j,k,top_gmax);
                    if (supersat[index] < 0)
                        continue;
		    cL[0] = top_L[0] + i*top_h[0];
		    cL[1] = top_L[1] + j*top_h[1];
		    cL[2] = top_L[2] + k*top_h[2]; 
		    for (l = 0; l < num_drops; l++)
		    {
			for (id = 0; id < dim; ++id)
			{
			    x = dist_uniform((POINTER)&uniform_params,xsubi);
			    particle_array[l+count].center[id] = cL[id]+x*top_h[id];
			    particle_array[l+count].vel[id] = 0.0;
			}
			particle_array[l+count].radius = gauss_center_limit((POINTER)&gauss_params,xsubi);
			particle_array[l+count].R0 = particle_array[l].radius;
			particle_array[l+count].flag = YES;
			particle_array[l+count].rho = eqn_params->rho_l;
		    }
		    count += num_drops;
		    field->drops[index] = num_drops;
                }
	
		break;
	    default:
		printf("Unknow dim = %d\n",dim);
		clean_up(ERROR);
	}
	setParticleGlobalIndex(particle_array,eqn_params->num_drops);
	eqn_params->particle_array = particle_array;
	computeSource();
}

static void computeFluctuation(Front* front, double **ext_accel, int size, int dim)
{
	int i, j;
	FILE* outfile;
	char filename[256]; 
	double mean_buoyancy[MAXD] = {0.0, 0.0, 0.0};
	for (i = 0; i < size; i++)
	for (j = 0; j < dim; j++)
	{
	    mean_buoyancy[j] += ext_accel[j][i];
	}
	
	/*find mean buoyancy*/
	pp_gsync();
	pp_global_sum(mean_buoyancy,3);
	for (j = 0; j < dim; j++)
	    mean_buoyancy[j] /= (pp_numnodes()*size);
	if (pp_mynode() == 0)
        {
            sprintf(filename,"%s/buoyancy",OutName(front));
            outfile = fopen(filename,"a");
            fprintf(outfile,"%f  %15.14f\n",
		    front->time,mean_buoyancy[dim-1]);
            fclose(outfile);
        }
	/*remove mean buoyancy to obtain neutral buoyancy*/
	for (i = 0; i < size; i++)
        {
            ext_accel[0][i] = 0.0;
            ext_accel[dim-2][i] = 0.0;
            ext_accel[dim-1][i] -= mean_buoyancy[dim-1];
        }
}

void VCARTESIAN::computeSource()
{
	/*compute condensation rate*/
	/*this is a source term for vapor mixing ratio equation*/
        int i, j, index, size, num_drops;
        int ic[MAXD];
        PARAMS* eqn_params = (PARAMS*)front->extra2;
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
        PARTICLE* particle_array = eqn_params->particle_array;
        num_drops = eqn_params->num_drops;
	double *supersat = eqn_params->field->supersat;
	double *T = eqn_params->field->temperature;
	double *qv = eqn_params->field->vapor;
	double *qc = eqn_params->field->cloud;
	double T0 = eqn_params->T0;
	double q0 = eqn_params->qv0;
        static boolean first = YES;
        double *coords;
	/*for computing the coefficient*/
	double rho_0 = iFparams->rho2;
	double a3 = 1.0;
	double coeff;

	static double maxsource = -HUGE, minsource = HUGE;

	for(i = 0; i < dim; i++)
	    a3 *= top_h[i];

        for (i = 0; i < comp_size; i++)
        {
	    source[i] = 0.0;
	    field->drops[i] = 0;
	    field->cloud[i] = 0.0;
	    field->mrad[i] = 0.0;
	}
	/*caculate num_drops in each cell: drops[index]*/
	/*compute source term for vapor equation: source[index]*/
	/*compute cloud water mixing ratio: qc[index]*/
        for (i = 0; i < num_drops; i++)
        {
            coords = particle_array[i].center;
            rect_in_which(coords,ic,top_grid);
            index = d_index(ic,top_gmax,dim);
	    coeff = 4.0*PI*particle_array[i].rho * eqn_params->K
					       / (rho_0 * a3);
	    source[index] += -1000.0 * coeff * supersat[index]
				     * particle_array[i].radius;
	    if (particle_array[i].radius != 0)
	    {
	        field->drops[index] += 1;
		field->mrad[index] += particle_array[i].radius;
		qc[index] += (4.0/3.0)*PI
				    *   pow(particle_array[i].radius,3)
				    *   particle_array[i].rho
				    /   (a3 * rho_0);
	    }
        }
	/*compute mean radius in a cell*/
	for (index = 0; index < comp_size; index++)
	{
		if (field->drops[index] != 0)
		    field->mrad[index] /= field->drops[index];
	}
	/*compute source for Navier Stokes equation:ext_accel[dim][index]*/
	for (index = 0; index < comp_size; index++)
	for (j = 0; j < dim; j++)
	{
	    field->ext_accel[j][index] = -iFparams->gravity[j]
	     *((T[index]-T0)/T0 + 0.608 * 0.001 * (qv[index] - q0) - qc[index]);
   	}
	/*remove mean value to obtain neutral buoyancy*/
	computeFluctuation(front,field->ext_accel,comp_size,dim);
}

/*******************************************************************
 * 		G_CARTESIAN.c
 *******************************************************************/
#include "cFluid.h"

class ToFill{
public:
int icoords[3];
};

EXPORT  void    tecplot_interface_states(const char*, INTERFACE	*);

static double (*getStateMom[MAXD])(Locstate) =
               {getStateXmom,getStateYmom,getStateZmom};

//----------------------------------------------------------------
//		L_RECTANGLE
//----------------------------------------------------------------

L_RECTANGLE::L_RECTANGLE(): m_index(-1), comp(-1)
{
}

void L_RECTANGLE::setCoords(
	double *coords,
	int dim)
{
	int i;
	for (i = 0; i < dim; ++i)
	    m_coords[i] = coords[i];
}
//--------------------------------------------------------------------------
// 		G_CARTESIAN
//--------------------------------------------------------------------------

G_CARTESIAN::~G_CARTESIAN()
{
}

//---------------------------------------------------------------
//	initMesh
// include the following parts
// 1) setup cell_center
//---------------------------------------------------------------
void G_CARTESIAN::initMesh(void)
{
	int i,j,k, index;
	double coords[2];
	int num_cells;

	// init cell_center
	L_RECTANGLE       rectangle;

	if (debugging("trace"))
	    (void) printf("Entering g_cartesian.initMesh()\n");
	/*TMP*/
	min_dens = 0.0001;
	min_pres = 0.0001;
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
	    	coords[0] = top_L[0] + top_h[0]*i;
		index = d_index1d(i,top_gmax);
	    	cell_center[index].setCoords(coords,dim);
	    	cell_center[index].icoords[0] = i;
	    }
	    break;
	case 2:
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {
	    	coords[0] = top_L[0] + top_h[0]*i;
	    	coords[1] = top_L[1] + top_h[1]*j;
		index = d_index2d(i,j,top_gmax);
	    	cell_center[index].setCoords(coords,dim);
	    	cell_center[index].icoords[0] = i;
	    	cell_center[index].icoords[1] = j;
	    }
	    break;
	case 3:
	    for (k = 0; k <= top_gmax[2]; k++)
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {
	    	coords[0] = top_L[0] + top_h[0]*i;
	    	coords[1] = top_L[1] + top_h[1]*j;
	    	coords[2] = top_L[2] + top_h[2]*k;
		index = d_index3d(i,j,k,top_gmax);
	    	cell_center[index].setCoords(coords,dim);
	    	cell_center[index].icoords[0] = i;
	    	cell_center[index].icoords[1] = j;
	    	cell_center[index].icoords[2] = k;
	    }
	}
	
	setComponent();
	FT_FreeGridIntfc(front);
	if (debugging("trace"))
	    (void) printf("Leaving g_cartesian.initMesh()\n");
}

void G_CARTESIAN::setComponent(void)
{
	int		i,j, ind;
	double 		*coords;
	int 		*icoords;
	COMPONENT 	old_comp,new_comp;
	double		***Gvel = eqn_params->Gvel;
	double		**Gdens = eqn_params->Gdens;
	double		**Gpres = eqn_params->Gpres;
	static STATE 	*state = NULL;
	double		*dens = field.dens;
	double		*engy = field.engy;
	double		*pres = field.pres;
	double		**momn = field.momn;
	int		size = (int)cell_center.size();
	
	// cell center components
	if(state == NULL)
	    FT_ScalarMemoryAlloc((POINTER*)&state,sizeof(STATE));

	for (i = 0; i < size; i++)
	{
	    icoords = cell_center[i].icoords;
	    coords = cell_center[i].m_coords;
	    old_comp = cell_center[i].comp;
	    new_comp = top_comp[i];
	    if (cell_center[i].comp != -1 &&
		cell_center[i].comp != top_comp[i] && gas_comp(new_comp))
	    {
		if (!FrontNearestIntfcState(front,coords,new_comp,
				(POINTER)state))
		{
		    (void) printf("In setComponent()\n");
		    (void) printf("FrontNearestIntfcState() failed\n");
		    (void) printf("old_comp = %d new_comp = %d\n",
					old_comp,new_comp);
		    clean_up(ERROR);
		}

		//GFM
		state->dim = dim;
		state->eos = &eqn_params->eos[new_comp];
		if (gas_comp(old_comp) && gas_comp(new_comp))
		{
		    if(new_comp == GAS_COMP1)
			ind = 0;
		    else
			ind = 1;

		    state->dens = Gdens[ind][i];
		    state->pres = Gpres[ind][i];
		    for(j = 0; j < dim; ++j)
			state->momn[j] = Gvel[ind][j][i]*Gdens[ind][i];
		    state->engy = EosEnergy(state);
		}

		dens[i] = state->dens;
		pres[i] = state->pres;
		engy[i] = state->engy;
		for (j = 0; j < dim; ++j)
		    momn[j][i] = state->momn[j];
	    }
	    cell_center[i].comp = top_comp[i];
	}
}	/* end setComponent() */

void G_CARTESIAN::computeAdvection(void)
{
	int order;
	switch (eqn_params->num_scheme)
	{
	case TVD_FIRST_ORDER:
	case WENO_FIRST_ORDER:
	    nrad = 3;
	    order = 1;
	    break;
	case TVD_SECOND_ORDER:
	case WENO_SECOND_ORDER:
	    nrad = 3;
	    order = 2;
	    break;
	case TVD_FOURTH_ORDER:
	case WENO_FOURTH_ORDER:
	    nrad = 3;
	    order = 4;
	    break;
	default:
	    order = -1;
	}
	solveRungeKutta(order);
}	/* end computeAdvection */


void G_CARTESIAN::solveRungeKutta(int order)
{
	static SWEEP *st_field,st_tmp;
	static FSWEEP *st_flux;
	static double **a,*b;
	double delta_t;
	int i,j;

	/* Allocate memory for Runge-Kutta of order */
	start_clock("solveRungeKutta");
	if (st_flux == NULL)
	{
	    FT_VectorMemoryAlloc((POINTER*)&b,order,sizeof(double));
	    FT_MatrixMemoryAlloc((POINTER*)&a,order,order,sizeof(double));

	    FT_VectorMemoryAlloc((POINTER*)&st_field,order,sizeof(SWEEP));
	    FT_VectorMemoryAlloc((POINTER*)&st_flux,order,sizeof(FSWEEP));
	    for (i = 0; i < order; ++i)
	    {
	    	allocMeshVst(&st_tmp);
	    	allocMeshVst(&st_field[i]);
	    	allocMeshFlux(&st_flux[i]);
	    }
	    /* Set coefficient a, b, c for different order of RK method */
	    switch (order)
	    {
	    case 1:
		b[0] = 1.0;
	    	break;
	    case 2:
	    	a[0][0] = 1.0;
	    	b[0] = 0.5;  b[1] = 0.5;
	    	break;
	    case 4:
	    	a[0][0] = 0.5;
	    	a[1][0] = 0.0;  a[1][1] = 0.5;
	    	a[2][0] = 0.0;  a[2][1] = 0.0;  a[2][2] = 1.0;
	    	b[0] = 1.0/6.0;  b[1] = 1.0/3.0;
	    	b[2] = 1.0/3.0;  b[3] = 1.0/6.0;
	    	break;
	    default:
	    	(void)printf("ERROR: %d-th order RK method not implemented\n",
					order);
	    	clean_up(ERROR);
	    }
	}
	delta_t = m_dt;

	/* Compute flux and advance field */

	copyToMeshVst(&st_field[0]);
	computeMeshFlux(st_field[0],&st_flux[0],delta_t);
	
	for (i = 0; i < order-1; ++i)
	{
	    copyMeshVst(st_field[0],&st_field[i+1]);
	    for (j = 0; j <= i; ++j)
	    {
		if (a[i][j] != 0.0)
		    addMeshFluxToVst(&st_field[i+1],st_flux[j],a[i][j]);
	    }
	    computeMeshFlux(st_field[i+1],&st_flux[i+1],delta_t);
	}
	for (i = 0; i < order; ++i)
	{
	    if (b[i] != 0.0)
		addMeshFluxToVst(&st_field[0],st_flux[i],b[i]);
	}
	copyFromMeshVst(st_field[0]);
	stop_clock("solveRungeKutta");
}	/* end solveRungeKutta */

void G_CARTESIAN::computeMeshFlux(
	SWEEP m_vst,
	FSWEEP *m_flux,
	double delta_t)
{
	int dir;

	if(eqn_params->tracked)
	{
	    start_clock("get_ghost_state");
	    get_ghost_state(m_vst, 2, 0);
	    get_ghost_state(m_vst, 3, 1);
	    scatMeshGhost();
	    stop_clock("get_ghost_state");
	    start_clock("solve_exp_value");
	    solve_exp_value();
	    stop_clock("solve_exp_value");
	}

	resetFlux(m_flux);
	for (dir = 0; dir < dim; ++dir)
	{
	    addFluxInDirection(dir,&m_vst,m_flux,delta_t);
	}
	addSourceTerm(&m_vst,m_flux,delta_t);
}	/* end computeMeshFlux */

void G_CARTESIAN::resetFlux(FSWEEP *m_flux)
{
	int i,j;
	int size = (int)cell_center.size();
	for (i = 0; i < size; i++)
	{
	    m_flux->dens_flux[i] = 0.0;
	    m_flux->engy_flux[i] = 0.0;
	    for (j = 0; j < MAXD; ++j)
	    	m_flux->momn_flux[j][i] = 0.0;
	}
}	/* resetFlux */

void G_CARTESIAN::addFluxInDirection(
	int dir,
	SWEEP *m_vst,
	FSWEEP *m_flux,
	double delta_t)
{
	switch (dim)
	{
	case 1:
	    return addFluxInDirection1d(dir,m_vst,m_flux,delta_t);
	case 2:
	    return addFluxInDirection2d(dir,m_vst,m_flux,delta_t);
	case 3:
	    return addFluxInDirection3d(dir,m_vst,m_flux,delta_t);
	}
}	/* end addFluxInDirection */

void G_CARTESIAN::addFluxInDirection1d(
	int dir,
	SWEEP *m_vst,
	FSWEEP *m_flux,
	double delta_t)
{
	int i,n,index;
	SCHEME_PARAMS scheme_params;
	EOS_PARAMS	*eos;
	static SWEEP vst;
	static FSWEEP vflux;
	static boolean first = YES;
	COMPONENT comp;
	int seg_min,seg_max;
	static int icoords[MAXD];
	
	start_clock("addFluxInDirection1d");
	if (first)
	{
	    first = NO;
	    allocDirVstFlux(&vst,&vflux);
	}

	scheme_params.lambda = delta_t/top_h[dir];
	scheme_params.beta = 0.0;

	seg_min = imin[0];
	while (seg_min <= imax[0])
	{
	    for (; seg_min <= imax[0]; ++seg_min)
	    {
		i = seg_min;
	    	index = d_index1d(i,top_gmax);
	    	comp = top_comp[index];
	    	if (gas_comp(comp)) break;
	    }
	    if (seg_min > imax[0]) break;
	    for (i = 0; i <= top_gmax[0]; ++i)
	    {
	    	vst.dens[i] = 0.0; 
	    	vst.pres[i] = 0.0; 
	    	vst.engy[i] = 0.0; 
	    	vst.momn[0][i] = vst.momn[1][i] = vst.momn[2][i] = 0.0;
	    }
	    i = seg_min;
	    index = d_index1d(i,top_gmax);
	    comp = top_comp[index];
	    n = 0;
	    vst.dens[n+nrad] = m_vst->dens[index];
                   vst.engy[n+nrad] = m_vst->engy[index];
                   vst.pres[n+nrad] = m_vst->pres[index];
                   vst.momn[0][n+nrad] = m_vst->momn[0][index];
                   vst.momn[1][n+nrad] = 0.0;
                   vst.momn[2][n+nrad] = 0.0;
	    seg_max = i;
	    n++;
	    for (i = seg_min+1; i <= imax[0]; i++)
	    {
		index = d_index1d(i,top_gmax);
		if (needBufferFromIntfc(comp,top_comp[index]))
		    break;
		else
		{
	    	    vst.dens[n+nrad] = m_vst->dens[index];
	    	    vst.engy[n+nrad] = m_vst->engy[index];
	    	    vst.pres[n+nrad] = m_vst->pres[index];
	    	    vst.momn[0][n+nrad] = m_vst->momn[0][index];
	    	    vst.momn[1][n+nrad] = 0.0;
	    	    vst.momn[2][n+nrad] = 0.0;
		    n++;
		}
		seg_max = i;
	    }
	    icoords[0] = seg_min;
	    appendGhostBuffer(&vst,m_vst,n,icoords,0,0);
	    icoords[0] = seg_max;
	    appendGhostBuffer(&vst,m_vst,n,icoords,0,1);
	    
	    eos = &(eqn_params->eos[comp]);
	    EosSetTVDParams(&scheme_params, eos);
	    numericalFlux((POINTER)&scheme_params,&vst,&vflux,n);
	    
	    n = 0;
	    for (i = seg_min; i <= seg_max; ++i)
	    {
	    	index = d_index1d(i,top_gmax);
	    	m_flux->dens_flux[index] += vflux.dens_flux[n+nrad];
	    	m_flux->engy_flux[index] += vflux.engy_flux[n+nrad];
	    	m_flux->momn_flux[0][index] += 
				vflux.momn_flux[0][n+nrad];
	    	m_flux->momn_flux[1][index] = 0.0;
	    	m_flux->momn_flux[2][index] = 0.0;
		n++;
	    }
	    seg_min = seg_max + 1;
	}
	stop_clock("addFluxInDirection1d");
}	/* end addFluxInDirection1d */

void G_CARTESIAN::addFluxInDirection2d(
	int dir,
	SWEEP *m_vst,
	FSWEEP *m_flux,
	double delta_t)
{
	int i,j,n,index;
	SCHEME_PARAMS scheme_params;
	EOS_PARAMS	*eos;
	static SWEEP vst;
	static FSWEEP vflux;
	static boolean first = YES;
	COMPONENT comp;
	int seg_min,seg_max;
	static int icoords[MAXD];
	
	start_clock("addFluxInDirection2d");
	if (first)
	{
	    first = NO;
	    allocDirVstFlux(&vst,&vflux);
	}

	scheme_params.lambda = delta_t/top_h[dir];
	scheme_params.beta = 0.0;
	switch (dir)
	{
	case 0:
	    for (j = imin[1]; j <= imax[1]; j++)
	    {
		seg_min = imin[0];
		while (seg_min <= imax[0])
		{
		    for (; seg_min <= imax[0]; ++seg_min)
		    {
			i = seg_min;
		    	index = d_index2d(i,j,top_gmax);
		    	comp = top_comp[index];
		    	if (gas_comp(comp)) break;
		    }
		    if (seg_min > imax[0]) break;
		    for (i = 0; i <= top_gmax[0]; ++i)
		    {
		    	vst.dens[i] = 0.0; 
		    	vst.pres[i] = 0.0; 
		    	vst.engy[i] = 0.0; 
		    	vst.momn[0][i] = vst.momn[1][i] = vst.momn[2][i] = 0.0;
		    }
		    i = seg_min;
		    index = d_index2d(i,j,top_gmax);
		    comp = top_comp[index];
		    n = 0;
		    vst.dens[n+nrad] = m_vst->dens[index];
                    vst.engy[n+nrad] = m_vst->engy[index];
                    vst.pres[n+nrad] = m_vst->pres[index];
                    vst.momn[0][n+nrad] = m_vst->momn[0][index];
                    vst.momn[1][n+nrad] = m_vst->momn[1][index];
                    vst.momn[2][n+nrad] = 0.0;
		    seg_max = i;
		    n++;
		    for (i = seg_min+1; i <= imax[0]; i++)
		    {
			index = d_index2d(i,j,top_gmax);
			if (needBufferFromIntfc(comp,top_comp[index]))
			    break;
			else
			{
		    	    vst.dens[n+nrad] = m_vst->dens[index];
		    	    vst.engy[n+nrad] = m_vst->engy[index];
		    	    vst.pres[n+nrad] = m_vst->pres[index];
		    	    vst.momn[0][n+nrad] = m_vst->momn[0][index];
		    	    vst.momn[1][n+nrad] = m_vst->momn[1][index];
		    	    vst.momn[2][n+nrad] = 0.0;
			    n++;
			}
			seg_max = i;
		    }
		    icoords[1] = j;
		    icoords[0] = seg_min;
		    appendGhostBuffer(&vst,m_vst,n,icoords,0,0);
		    icoords[0] = seg_max;
		    appendGhostBuffer(&vst,m_vst,n,icoords,0,1);
		    
		    eos = &(eqn_params->eos[comp]);
		    EosSetTVDParams(&scheme_params, eos);
		    numericalFlux((POINTER)&scheme_params,&vst,&vflux,n);
		    
		    n = 0;
		    for (i = seg_min; i <= seg_max; ++i)
		    {
		    	index = d_index2d(i,j,top_gmax);
		    	m_flux->dens_flux[index] += vflux.dens_flux[n+nrad];
		    	m_flux->engy_flux[index] += vflux.engy_flux[n+nrad];
		    	m_flux->momn_flux[0][index] += 
					vflux.momn_flux[0][n+nrad];
		    	m_flux->momn_flux[1][index] += 
					vflux.momn_flux[1][n+nrad];
		    	m_flux->momn_flux[2][index] = 0.0;
			n++;
		    }
		    seg_min = seg_max + 1;
		}
	    }
	    break;
	case 1:
	    for (i = imin[0]; i <= imax[0]; i++)
	    {
		seg_min = imin[1];
		while (seg_min <= imax[1])
		{
		    for (; seg_min <= imax[1]; ++seg_min)
		    {
			j = seg_min;
		    	index = d_index2d(i,j,top_gmax);
		    	comp = top_comp[index];
		    	if (gas_comp(comp)) break;
		    }
		    if (seg_min > imax[1]) break;
		    for (j = 0; j <= top_gmax[1]; ++j)
		    {
		    	vst.dens[j] = 0.0; 
		    	vst.pres[j] = 0.0; 
		    	vst.engy[j] = 0.0; 
		    	vst.momn[0][j] = vst.momn[1][j] = vst.momn[2][j] = 0.0;
		    }
		    j = seg_min;
		    index = d_index2d(i,j,top_gmax);
		    comp = top_comp[index];
		    n = 0;
		    vst.dens[n+nrad] = m_vst->dens[index];
                    vst.engy[n+nrad] = m_vst->engy[index];
                    vst.pres[n+nrad] = m_vst->pres[index];
                    vst.momn[0][n+nrad] = m_vst->momn[1][index];
                    vst.momn[1][n+nrad] = m_vst->momn[0][index];
                    vst.momn[2][n+nrad] = 0.0;
		    seg_max = j;
		    n++;
		    for (j = seg_min+1; j <= imax[1]; j++)
		    {
			index = d_index2d(i,j,top_gmax);
			if (needBufferFromIntfc(comp,top_comp[index]))
			    break;
			else
			{
		    	    vst.dens[n+nrad] = m_vst->dens[index];
		    	    vst.engy[n+nrad] = m_vst->engy[index];
		    	    vst.pres[n+nrad] = m_vst->pres[index];
		    	    vst.momn[0][n+nrad] = m_vst->momn[1][index];
		    	    vst.momn[1][n+nrad] = m_vst->momn[0][index];
		    	    vst.momn[2][n+nrad] = 0.0;
			    n++;
			}
			seg_max = j;
		    }
		    icoords[0] = i;
		    icoords[1] = seg_min;
		    appendGhostBuffer(&vst,m_vst,n,icoords,1,0);
		    icoords[1] = seg_max;
		    appendGhostBuffer(&vst,m_vst,n,icoords,1,1);
		    
		    eos = &(eqn_params->eos[comp]);
		    EosSetTVDParams(&scheme_params, eos);
		    numericalFlux((POINTER)&scheme_params,&vst,&vflux,n);
		    
		    n = 0;
		    for (j = seg_min; j <= seg_max; ++j)
		    {
		    	index = d_index2d(i,j,top_gmax);
		    	m_flux->dens_flux[index] += vflux.dens_flux[n+nrad];
		    	m_flux->engy_flux[index] += vflux.engy_flux[n+nrad];
		    	m_flux->momn_flux[1][index] += 
					vflux.momn_flux[0][n+nrad];
		    	m_flux->momn_flux[0][index] += 
					vflux.momn_flux[1][n+nrad];
		    	m_flux->momn_flux[2][index] = 0.0;
			n++;
		    }
		    seg_min = seg_max + 1;
		}
	    }
	    break;
	}
	stop_clock("addFluxInDirection2d");
}	/* end addFluxInDirection2d */

void G_CARTESIAN::addFluxInDirection3d(
	int dir,
	SWEEP *m_vst,
	FSWEEP *m_flux,
	double delta_t)
{
	int		i,j,k,n,index;
	SCHEME_PARAMS	scheme_params;
	EOS_PARAMS	*eos;
	static SWEEP 	vst;
	static FSWEEP 	vflux;
	static boolean 	first = YES;
	COMPONENT 	comp;
	int 		seg_min,seg_max;
	int 		icoords[3];
	
	start_clock("addFluxInDirection3d");
	if (first)
	{
	    first = NO;
	    allocDirVstFlux(&vst,&vflux);
	}
	
	scheme_params.lambda = delta_t/top_h[dir];
	scheme_params.beta = 0.0;
	
	switch (dir)
	{
	case 0:
	    for (k = imin[2]; k <= imax[2]; k++)
	    for (j = imin[1]; j <= imax[1]; j++)
	    {
		seg_min = imin[0];
		while (seg_min <= imax[0])
		{
		    for (; seg_min <= imax[0]; ++seg_min)
                    {
                        i = seg_min;
                        index = d_index3d(i,j,k,top_gmax);
                        comp = top_comp[index];
                        if (gas_comp(comp)) break;
                    }
                    if (seg_min > imax[0]) break;
		    for (i = 0; i <= top_gmax[1]; ++i)
		    {
		    	vst.dens[i] = 0.0; 
		    	vst.pres[i] = 0.0; 
		    	vst.engy[i] = 0.0; 
		    	vst.momn[0][i] = vst.momn[1][i] = vst.momn[2][i] = 0.0;
		    }
		    i = seg_min;
		    index = d_index3d(i,j,k,top_gmax);
		    comp = top_comp[index];
		    n = 0;
		    vst.dens[n+nrad] = m_vst->dens[index];
                    vst.engy[n+nrad] = m_vst->engy[index];
                    vst.pres[n+nrad] = m_vst->pres[index];
                    vst.momn[0][n+nrad] = m_vst->momn[0][index];
                    vst.momn[1][n+nrad] = m_vst->momn[1][index];
                    vst.momn[2][n+nrad] = m_vst->momn[2][index];
		    seg_max = i;
		    n++;
		    for (i = seg_min+1; i <= imax[0]; i++)
		    {
			index = d_index3d(i,j,k,top_gmax);
			if (needBufferFromIntfc(comp,top_comp[index]))
                            break;
			else
			{
		    	    vst.dens[n+nrad] = m_vst->dens[index];
		    	    vst.engy[n+nrad] = m_vst->engy[index];
		    	    vst.pres[n+nrad] = m_vst->pres[index];
		    	    vst.momn[0][n+nrad] = m_vst->momn[0][index];
		    	    vst.momn[1][n+nrad] = m_vst->momn[1][index];
		    	    vst.momn[2][n+nrad] = m_vst->momn[2][index];
			    n++;
			}
			seg_max = i;
		    }
		    
		    icoords[1] = j;
		    icoords[2] = k;
		    icoords[0] = seg_min;
		    appendGhostBuffer(&vst,m_vst,n,icoords,0,0);
		    icoords[0] = seg_max;
		    appendGhostBuffer(&vst,m_vst,n,icoords,0,1);
		    
		    eos = &(eqn_params->eos[comp]);
		    EosSetTVDParams(&scheme_params, eos);
		    numericalFlux((POINTER)&scheme_params,&vst,&vflux,n);
		    
		    n = 0;
		    for (i = seg_min; i <= seg_max; ++i)
		    {
		    	index = d_index3d(i,j,k,top_gmax);
		    	m_flux->dens_flux[index] += vflux.dens_flux[n+nrad];
		    	m_flux->engy_flux[index] += vflux.engy_flux[n+nrad];
		    	m_flux->momn_flux[0][index] += 
					vflux.momn_flux[0][n+nrad];
		    	m_flux->momn_flux[1][index] += 
					vflux.momn_flux[1][n+nrad];
		    	m_flux->momn_flux[2][index] +=
					vflux.momn_flux[2][n+nrad];
			n++;
		    }

		    seg_min = seg_max + 1;
		}
	    }
	    break;
	case 1:
	    for (k = imin[2]; k <= imax[2]; k++)
	    for (i = imin[0]; i <= imax[0]; i++)
	    {
		seg_min = imin[1];
		while (seg_min <= imax[1])
		{
		    for (; seg_min <= imax[1]; ++seg_min)
                    {
                        j = seg_min;
                        index = d_index3d(i,j,k,top_gmax);
                        comp = top_comp[index];
                        if (gas_comp(comp)) break;
                    }
                    if (seg_min > imax[1]) break;
		    for (j = 0; j <= top_gmax[1]; ++j)
		    {
		    	vst.dens[j] = 0.0; 
		    	vst.pres[j] = 0.0; 
		    	vst.engy[j] = 0.0; 
		    	vst.momn[0][j] = vst.momn[1][j] = vst.momn[2][j] = 0.0;
		    }
		    j = seg_min;
		    index = d_index3d(i,j,k,top_gmax);
		    comp = top_comp[index];
		    n = 0;
		    vst.dens[n+nrad] = m_vst->dens[index];
                    vst.engy[n+nrad] = m_vst->engy[index];
                    vst.pres[n+nrad] = m_vst->pres[index];
                    vst.momn[0][n+nrad] = m_vst->momn[1][index];
                    vst.momn[1][n+nrad] = m_vst->momn[2][index];
                    vst.momn[2][n+nrad] = m_vst->momn[0][index];
		    seg_max = j;
		    n++;
		    
		    for (j = seg_min+1; j <= imax[1]; j++)
		    {
			index = d_index3d(i,j,k,top_gmax);
			if (needBufferFromIntfc(comp,top_comp[index]))
			    break;
			else
			{
		    	    vst.dens[n+nrad] = m_vst->dens[index];
		    	    vst.engy[n+nrad] = m_vst->engy[index];
		    	    vst.pres[n+nrad] = m_vst->pres[index];
		    	    vst.momn[0][n+nrad] = m_vst->momn[1][index];
		    	    vst.momn[1][n+nrad] = m_vst->momn[2][index];
		    	    vst.momn[2][n+nrad] = m_vst->momn[0][index];
			    n++;
			}
			seg_max = j;
		    }
		    icoords[0] = i;
		    icoords[2] = k;
		    icoords[1] = seg_min;
		    appendGhostBuffer(&vst,m_vst,n,icoords,1,0);
		    icoords[1] = seg_max;
		    appendGhostBuffer(&vst,m_vst,n,icoords,1,1);
		    
		    eos = &(eqn_params->eos[comp]);
		    EosSetTVDParams(&scheme_params, eos);
		    numericalFlux((POINTER)&scheme_params,&vst,&vflux,n);
		    
		    n = 0;
		    for (j = seg_min; j <= seg_max; ++j)
		    {
		    	index = d_index3d(i,j,k,top_gmax);
		    	m_flux->dens_flux[index] += vflux.dens_flux[n+nrad];
		    	m_flux->engy_flux[index] += vflux.engy_flux[n+nrad];
		    	m_flux->momn_flux[1][index] += 
					vflux.momn_flux[0][n+nrad];
		    	m_flux->momn_flux[0][index] += 
					vflux.momn_flux[2][n+nrad];
		    	m_flux->momn_flux[2][index] += 
					vflux.momn_flux[1][n+nrad];
			n++;
		    }
		    seg_min = seg_max + 1;
		}
	    }
	    break;
	case 2:
	    for (j = imin[1]; j <= imax[1]; j++)
	    for (i = imin[0]; i <= imax[0]; i++)
	    {
		seg_min = imin[2];
		while (seg_min <= imax[2])
		{
		    for (; seg_min <= imax[2]; ++seg_min)
                    {
                        k = seg_min;
                        index = d_index3d(i,j,k,top_gmax);
                        comp = top_comp[index];
                        if (gas_comp(comp)) break;
                    }
                    if (seg_min > imax[2]) break;
		    for (k = 0; k <= top_gmax[2]; ++k)
		    {
		    	vst.dens[k] = 0.0; 
		    	vst.pres[k] = 0.0; 
		    	vst.engy[k] = 0.0; 
		    	vst.momn[0][k] = vst.momn[1][k] = vst.momn[2][k] = 0.0;
		    }
		    k = seg_min;
		    index = d_index3d(i,j,k,top_gmax);
		    comp = top_comp[index];
		    n = 0;
		    vst.dens[n+nrad] = m_vst->dens[index];
                    vst.engy[n+nrad] = m_vst->engy[index];
                    vst.pres[n+nrad] = m_vst->pres[index];
                    vst.momn[0][n+nrad] = m_vst->momn[2][index];
                    vst.momn[1][n+nrad] = m_vst->momn[0][index];
                    vst.momn[2][n+nrad] = m_vst->momn[1][index];
		    seg_max = k;
		    n++;
		    
		    for (k = seg_min+1; k <= imax[2]; k++)
		    {
			index = d_index3d(i,j,k,top_gmax);
			if (needBufferFromIntfc(comp,top_comp[index]))
			    break;
			else
			{
		    	    vst.dens[n+nrad] = m_vst->dens[index];
		    	    vst.engy[n+nrad] = m_vst->engy[index];
		    	    vst.pres[n+nrad] = m_vst->pres[index];
		    	    vst.momn[0][n+nrad] = m_vst->momn[2][index];
		    	    vst.momn[1][n+nrad] = m_vst->momn[0][index];
		    	    vst.momn[2][n+nrad] = m_vst->momn[1][index];
			    n++;
			}
			seg_max = k;
		    }
		    icoords[0] = i;
		    icoords[1] = j;
		    icoords[2] = seg_min;
		    appendGhostBuffer(&vst,m_vst,n,icoords,2,0);
		    icoords[2] = seg_max;
		    appendGhostBuffer(&vst,m_vst,n,icoords,2,1);
		    
		    eos = &(eqn_params->eos[comp]);
		    EosSetTVDParams(&scheme_params, eos);
		    numericalFlux((POINTER)&scheme_params,&vst,&vflux,n);
		    
		    n = 0;
		    for (k = seg_min; k <= seg_max; ++k)
		    {
		    	index = d_index3d(i,j,k,top_gmax);
		    	m_flux->dens_flux[index] += vflux.dens_flux[n+nrad];
		    	m_flux->engy_flux[index] += vflux.engy_flux[n+nrad];
		    	m_flux->momn_flux[2][index] += 
					vflux.momn_flux[0][n+nrad];
		    	m_flux->momn_flux[0][index] += 
					vflux.momn_flux[1][n+nrad];
		    	m_flux->momn_flux[1][index] += 
					vflux.momn_flux[2][n+nrad];
			n++;
		    }
		    seg_min = seg_max + 1;
		}
	    }
	    break;
	}
	stop_clock("addFluxInDirection3d");
}

void G_CARTESIAN::scatMeshFlux(FSWEEP *m_flux)
{
	int i,j,k,l,index;

	switch (dim)
	{
	case 1:
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index1d(i,top_gmax);
		array[index] = m_flux->dens_flux[index];
	    }
	    scatMeshArray();
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index1d(i,top_gmax);
		m_flux->dens_flux[index] = array[index];
	    }
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index1d(i,top_gmax);
		array[index] = m_flux->engy_flux[index];
	    }
	    scatMeshArray();
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index1d(i,top_gmax);
		m_flux->engy_flux[index] = array[index];
	    }
	    for (l = 0; l < dim; ++l)
	    {
	    	for (i = imin[0]; i <= imax[0]; ++i)
	    	{
		    index = d_index1d(i,top_gmax);
		    array[index] = m_flux->momn_flux[l][index];
	    	}
	    	scatMeshArray();
            	for (i = 0; i <= top_gmax[0]; i++)
            	{
                    index = d_index1d(i,op_gmax);
		    m_flux->momn_flux[l][index] = array[index];
	    	}
	    }
	    break;
	case 2:
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index2d(i,j,top_gmax);
		array[index] = m_flux->dens_flux[index];
	    }
	    scatMeshArray();
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index2d(i,j,top_gmax);
		m_flux->dens_flux[index] = array[index];
	    }
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index2d(i,j,top_gmax);
		array[index] = m_flux->engy_flux[index];
	    }
	    scatMeshArray();
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index2d(i,j,top_gmax);
		m_flux->engy_flux[index] = array[index];
	    }
	    for (l = 0; l < dim; ++l)
	    {
	    	for (j = imin[1]; j <= imax[1]; ++j)
	    	for (i = imin[0]; i <= imax[0]; ++i)
	    	{
		    index = d_index2d(i,j,top_gmax);
		    array[index] = m_flux->momn_flux[l][index];
	    	}
	    	scatMeshArray();
            	for (j = 0; j <= top_gmax[1]; j++)
            	for (i = 0; i <= top_gmax[0]; i++)
            	{
                    index = d_index2d(i,j,top_gmax);
		    m_flux->momn_flux[l][index] = array[index];
	    	}
	    }
	    break;
	case 3:
	    for (k = imin[2]; k <= imax[2]; ++k)
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index3d(i,j,k,top_gmax);
		array[index] = m_flux->dens_flux[index];
	    }
	    scatMeshArray();
            for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
		index = d_index3d(i,j,k,top_gmax);
		m_flux->dens_flux[index] = array[index];
	    }
	    for (k = imin[2]; k <= imax[2]; ++k)
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index3d(i,j,k,top_gmax);
		array[index] = m_flux->engy_flux[index];
	    }
	    scatMeshArray();
            for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index2d(i,j,top_gmax);
		m_flux->engy_flux[index] = array[index];
	    }
	    for (l = 0; l < dim; ++l)
	    {
	    	for (k = imin[2]; k <= imax[2]; ++k)
	    	for (j = imin[1]; j <= imax[1]; ++j)
	    	for (i = imin[0]; i <= imax[0]; ++i)
	    	{
		    index = d_index3d(i,j,k,top_gmax);
		    array[index] = m_flux->momn_flux[l][index];
	    	}
	    	scatMeshArray();
            	for (k = 0; k <= top_gmax[2]; k++)
            	for (j = 0; j <= top_gmax[1]; j++)
            	for (i = 0; i <= top_gmax[0]; i++)
            	{
		    index = d_index3d(i,j,k,top_gmax);
		    m_flux->momn_flux[l][index] = array[index];
	    	}
	    }
	}
}	/* end scatMeshFlux */

void G_CARTESIAN::addSourceTerm(
	SWEEP *m_vst,
	FSWEEP *m_flux,
	double delta_t)
{
	int i,j,k,l,index;
	double *gravity = eqn_params->gravity;

	switch (dim)
	{
	case 1:
            for (i = imin[0]; i <= imax[0]; i++)
            {
		index = d_index1d(i,top_gmax);
		if (!gas_comp(top_comp[index]))
		{
		    for (l = 0; l < dim; ++l)
		    {
		    	m_flux->momn_flux[l][index] = 0.0; 
		    	m_flux->engy_flux[index] = 0.0; 
		    }
		}
		else
		{
		    for (l = 0; l < dim; ++l)
		    {
		    	m_flux->momn_flux[l][index] += 
				delta_t*gravity[l]*m_vst->dens[index];
		    	m_flux->engy_flux[index] += 
				delta_t*gravity[l]*m_vst->momn[l][index];
		    }
		}
	    }
	    break;
	case 2:
            for (j = imin[1]; j <= imax[1]; j++)
            for (i = imin[0]; i <= imax[0]; i++)
            {
		index = d_index2d(i,j,top_gmax);
		if (!gas_comp(top_comp[index]))
		{
		    for (l = 0; l < dim; ++l)
		    {
		    	m_flux->momn_flux[l][index] = 0.0; 
		    	m_flux->engy_flux[index] = 0.0; 
		    }
		}
		else
		{
		    for (l = 0; l < dim; ++l)
		    {
		    	m_flux->momn_flux[l][index] += 
				delta_t*gravity[l]*m_vst->dens[index];
		    	m_flux->engy_flux[index] += 
				delta_t*gravity[l]*m_vst->momn[l][index];
		    }
		}
	    }
	    break;
	case 3:
            for (k = imin[2]; k <= imax[2]; k++)
            for (j = imin[1]; j <= imax[1]; j++)
            for (i = imin[0]; i <= imax[0]; i++)
            {
		index = d_index3d(i,j,k,top_gmax);
		if (!gas_comp(top_comp[index]))
		{
		    for (l = 0; l < dim; ++l)
		    {
		    	m_flux->momn_flux[l][index] = 0.0; 
		    	m_flux->engy_flux[index] = 0.0; 
		    }
		}
		else
		{
		    for (l = 0; l < dim; ++l)
		    {
		    	m_flux->momn_flux[l][index] += 
				delta_t*gravity[l]*m_vst->dens[index];
		    	m_flux->engy_flux[index] += 
				delta_t*gravity[l]*m_vst->momn[l][index];
		    }
		}
	    }
	}
	
}	/* end addSourceTerm */

// for initial condition: 
// 		setInitialCondition();	
// this function should be called before solve()
// for the source term of the momentum equation: 	
// 		computeSourceTerm();
void G_CARTESIAN::solve(double dt)
{
	m_dt = dt;
	max_speed = 0.0;

	if (debugging("trace"))
	    printf("Entering solve()\n");
	start_clock("solve");
	setDomain();

	setComponent();
	
	if (debugging("trace"))
	    printf("Passed setComponent()\n");

	// 1) solve for intermediate velocity
	start_clock("computeAdvection");
	computeAdvection();
	if (debugging("trace"))
	    printf("max_speed after computeAdvection(): %20.14f\n",max_speed);
	stop_clock("computeAdvection");
	
	if (debugging("sample_velocity"))
	{
	    sampleVelocity();
	}

	start_clock("copyMeshStates");
	copyMeshStates();
	stop_clock("copyMeshStates");

	setAdvectionDt();
	stop_clock("solve");
	if (debugging("trace"))
	    printf("Leaving solve()\n");
}	/* end solve */


// check http://en.wikipedia.org/wiki/Bilinear_interpolation
void G_CARTESIAN::getVelocity(double *p, double *U)
{
        double **vel = eqn_params->vel;

        FT_IntrpStateVarAtCoords(front,NO_COMP,p,vel[0],getStateXvel,&U[0],
					NULL);
        if (dim > 1)
            FT_IntrpStateVarAtCoords(front,NO_COMP,p,vel[1],getStateYvel,&U[1],
					NULL);
        if (dim > 2)
            FT_IntrpStateVarAtCoords(front,NO_COMP,p,vel[2],getStateZvel,&U[2],
					NULL);
}

void G_CARTESIAN::getRectangleIndex(int index, int &i, int &j)
{
	i = cell_center[index].icoords[0];
	j = cell_center[index].icoords[1];
}

void G_CARTESIAN::getRectangleIndex(int index, int &i, int &j, int &k)
{
	i = cell_center[index].icoords[0];
	j = cell_center[index].icoords[1];
	k = cell_center[index].icoords[2];
}


int G_CARTESIAN::getRectangleComponent(int index)
{	
	return getComponent(cell_center[index].icoords);
}

void G_CARTESIAN::getRectangleCenter(
	int index, 
	double *coords)
{
	int i;
	for (i = 0; i < dim; ++i)
	    coords[i] = cell_center[index].m_coords[i];
}

void G_CARTESIAN::getRectangleCenter(
	int index0, 
	int index1, 
	double *coords)
{
	int i;
	for (i = 0; i < dim; ++i)
	{
	    coords[i] = 0.5*(cell_center[index0].m_coords[i] +
	    		     cell_center[index1].m_coords[i]);
	}
}


double G_CARTESIAN::getDistance(double *c0, double *c1)
{
	return sqrt( (c0[0]-c1[0])*(c0[0]-c1[0])
		    +(c0[1]-c1[1])*(c0[1]-c1[1]) );
}


// input : p[]
// output: q[]

void G_CARTESIAN::getNearestInterfacePoint(
	double *p, 
	double *q)
{
	INTERFACE *intfc = front->interf;
	double t;
	HYPER_SURF_ELEMENT *phse;
	HYPER_SURF *phs;
	nearest_interface_point(p,getComponent(p),intfc,NO_BOUNDARIES,
				NULL,q,&t,&phse,&phs);
}

int G_CARTESIAN::getComponent(
	double *p)
{
	return component(p,front->interf);
}

int G_CARTESIAN::getComponent(
	int *icoords)
{
	int index;
	switch (dim)
	{
	case 2:
	    index = d_index2d(icoords[0],icoords[1],top_gmax);
	    return top_comp[index];
	case 3:
	    index = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
	    return top_comp[index];
	default:
	    return NO_COMP;
	}
}

void G_CARTESIAN::save(char *filename)
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

G_CARTESIAN::G_CARTESIAN(Front &front):front(&front)
{
}

void G_CARTESIAN::setDomain()
{
	static boolean first = YES;
	INTERFACE *grid_intfc;
	Table *T;
	int i,size;

	grid_intfc = front->grid_intfc;
	top_grid = &topological_grid(grid_intfc);
	T = table_of_interface(grid_intfc);
	top_comp = T->components;
	eqn_params = (EQN_PARAMS*)front->extra1;
	
	if (first)
	{
	    first = NO;
	    dim = grid_intfc->dim;

	    hmin = HUGE;
	    size = 1;
	    
            for (i = 0; i < 3; ++i)
	    	top_gmax[i] = 0;

            for (i = 0; i < dim; ++i)
	    {
	    	lbuf[i] = front->rect_grid->lbuf[i];
	    	ubuf[i] = front->rect_grid->ubuf[i];
	    	top_gmax[i] = top_grid->gmax[i];
	    	top_L[i] = top_grid->L[i];
	    	top_U[i] = top_grid->U[i];
	    	top_h[i] = top_grid->h[i];

                if (hmin > top_h[i]) hmin = top_h[i];
	        size *= (top_gmax[i]+1);
	    	imin[i] = (lbuf[i] == 0) ? 1 : lbuf[i];
	    	imax[i] = (ubuf[i] == 0) ? top_gmax[i] - 1 : 
				top_gmax[i] - ubuf[i];
	    }

	    FT_VectorMemoryAlloc((POINTER*)&eqn_params->dens,size,
					sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&eqn_params->pres,size,
					sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&eqn_params->engy,size,
					sizeof(double));
	    FT_MatrixMemoryAlloc((POINTER*)&eqn_params->vel,dim,size,
					sizeof(double));
	    FT_MatrixMemoryAlloc((POINTER*)&eqn_params->mom,dim,size,
					sizeof(double));
	    //GFM
	    FT_MatrixMemoryAlloc((POINTER*)&eqn_params->gnor,dim,size,
					sizeof(double));
	    FT_MatrixMemoryAlloc((POINTER*)&eqn_params->Gdens,2,size,
					sizeof(double));
	    FT_MatrixMemoryAlloc((POINTER*)&eqn_params->Gpres,2,size,
					sizeof(double));
	    FT_TriArrayMemoryAlloc((POINTER*)&eqn_params->Gvel,2,dim,size,
					sizeof(double));

	    FT_VectorMemoryAlloc((POINTER*)&array,size,sizeof(double));
	    if (dim == 2)
	    	FT_VectorMemoryAlloc((POINTER*)&eqn_params->vort,size,
					sizeof(double));
	    field.dens = eqn_params->dens;
	    field.engy = eqn_params->engy;
	    field.pres = eqn_params->pres;
	    field.momn = eqn_params->mom;
	    field.vel = eqn_params->vel;
	}
}

void G_CARTESIAN::allocMeshVst(
	SWEEP *vst)
{
	int i,size;

	size = 1;
        for (i = 0; i < dim; ++i)
	    size *= (top_gmax[i]+1);

	FT_VectorMemoryAlloc((POINTER*)&vst->dens,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&vst->engy,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&vst->pres,size,sizeof(double));
	FT_MatrixMemoryAlloc((POINTER*)&vst->momn,MAXD,size,sizeof(double));
}	/* end allocMeshVstFlux */

void G_CARTESIAN::allocMeshFlux(
	FSWEEP *flux)
{
	int i,size;

	size = 1;
        for (i = 0; i < dim; ++i)
	    size *= (top_gmax[i]+1);

	FT_VectorMemoryAlloc((POINTER*)&flux->dens_flux,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&flux->engy_flux,size,sizeof(double));
	FT_MatrixMemoryAlloc((POINTER*)&flux->momn_flux,MAXD,size,sizeof(double));
}	/* end allocMeshVstFlux */

void G_CARTESIAN::allocDirVstFlux(
        SWEEP *vst,
        FSWEEP *flux)
{
	int i,size;

	size = 1;
        for (i = 0; i < dim; ++i)
	    if (size < top_gmax[i]+7) 
		size = top_gmax[i]+7;
	FT_VectorMemoryAlloc((POINTER*)&vst->dens,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&vst->engy,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&vst->pres,size,sizeof(double));
	FT_MatrixMemoryAlloc((POINTER*)&vst->momn,MAXD,size,sizeof(double));

	FT_VectorMemoryAlloc((POINTER*)&flux->dens_flux,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&flux->engy_flux,size,sizeof(double));
	FT_MatrixMemoryAlloc((POINTER*)&flux->momn_flux,MAXD,size,sizeof(double));
}	/* end allocDirMeshVstFlux */

void G_CARTESIAN::scatMeshArray()
{
	FT_ParallelExchGridArrayBuffer(array,front);
}

void G_CARTESIAN::checkVst(SWEEP *vst)
{
	int i,j,index;
	for (j = imin[1]; j < imax[1]; j++)
	for (i = imin[0]; i < imax[0]; i++)
	{	
	    index  = d_index2d(i,j,top_gmax);
	    if (isnan(vst->dens[index]))
		printf("At %d %d: dens is nan\n",i,j);
	    if (vst->dens[index] < 0.0)
		printf("At %d %d: dens is negative\n",i,j);
	}
}

void G_CARTESIAN::checkFlux(FSWEEP *flux)
{
	int i,j,index;
	//for (j = imin[1]; j < imax[1]; j++)
	j = 140;
	for (i = imin[0]; i <= imax[0]; i++)
	{	
	    index  = d_index2d(i,j,top_gmax);
	    printf("%d %f  %f\n",i,flux->momn_flux[1][index],
				flux->engy_flux[index]);
	}
}

void G_CARTESIAN::printFrontInteriorStates(char *out_name)
{
	int i,j,k,l,index;
	char filename[100];
	FILE *outfile;
	INTERFACE *intfc = front->interf;
        STATE *sl,*sr;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
	double *dens = field.dens;
	double *engy = field.engy;
	double **momn = field.momn;

	sprintf(filename,"%s/state.ts%s",out_name,
			right_flush(front->step,7));
	if (pp_numnodes() > 1)
            sprintf(filename,"%s-nd%s",filename,right_flush(pp_mynode(),4));
	sprintf(filename,"%s-gas",filename);
	outfile = fopen(filename,"w");

        /* Initialize states at the interface */
        fprintf(outfile,"Interface gas states:\n");
        next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
            FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
            fprintf(outfile,"%24.18g %24.18g\n",getStateDens(sl),
				getStateDens(sr));
            fprintf(outfile,"%24.18g %24.18g\n",getStateEngy(sl),
				getStateEngy(sr));
	    for (i = 0; i < dim; ++i)
            	fprintf(outfile,"%24.18g %24.18g\n",getStateMom[i](sl),
				getStateMom[i](sr));
        }
	
	fprintf(outfile,"\nInterior gas states:\n");
	switch (dim)
	{
	case 1:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    {
		index = d_index1d(i,top_gmax);
	        fprintf(outfile,"%24.18g\n",dens[index]);
	        fprintf(outfile,"%24.18g\n",engy[index]);
	    	for (l = 0; l < dim; ++l)
	            fprintf(outfile,"%24.18g\n",momn[l][index]);
	    }
	    break;
	case 2:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    {
		index = d_index2d(i,j,top_gmax);
	        fprintf(outfile,"%24.18g\n",dens[index]);
	        fprintf(outfile,"%24.18g\n",engy[index]);
	    	for (l = 0; l < dim; ++l)
	            fprintf(outfile,"%24.18g\n",momn[l][index]);
	    }
	    break;
	case 3:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    for (k = 0; k <= top_gmax[2]; ++k)
	    {
		index = d_index3d(i,j,k,top_gmax);
	        fprintf(outfile,"%24.18g\n",dens[index]);
	        fprintf(outfile,"%24.18g\n",engy[index]);
	    	for (l = 0; l < dim; ++l)
	            fprintf(outfile,"%24.18g\n",momn[l][index]);
	    }
	}
	fclose(outfile);
}

void G_CARTESIAN::readInteriorStates(char *restart_name)
{
	FILE *infile;
	int i,j,k,l,index;
	STATE st_tmp;
	char fname[100];
	int		comp;
	EOS_PARAMS	*eos = eqn_params->eos;
	double *dens = field.dens;
	double *engy = field.engy;
	double *pres = field.pres;
	double **momn = field.momn;

	setDomain();
	m_dens[0] = eqn_params->rho1;		
	m_dens[1] = eqn_params->rho2;		
	m_mu[0] = eqn_params->mu1;		
	m_mu[1] = eqn_params->mu2;		
	if (eqn_params->prob_type == FLUID_SOLID_CIRCLE ||
	    eqn_params->prob_type == FLUID_RIGID_BODY ||
	    eqn_params->prob_type == FLUID_CRYSTAL)
	    m_comp[0] = SOLID_COMP;
	else
	    m_comp[0] = GAS_COMP1;
	m_comp[1] = GAS_COMP2;
	m_smoothing_radius = top_h[0] < top_h[1] ? top_h[1] : top_h[0];
	m_smoothing_radius *= 2.0;
	
	st_tmp.dim = eqn_params->dim;

	sprintf(fname,"%s-gas",restart_name);
	infile = fopen(fname,"r");
	
	next_output_line_containing_string(infile,"Interior gas states:");

	switch (dim)
	{
	case 2:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    {
		index = d_index2d(i,j,top_gmax);
		comp = top_comp[index];
		st_tmp.eos = &(eos[comp]);
	    	
		fscanf(infile,"%lf",&dens[index]);
	    	fscanf(infile,"%lf",&engy[index]);
		st_tmp.dens = dens[index];
		st_tmp.engy = engy[index];
		for (l = 0; l < dim; ++l)
		{
	    	    fscanf(infile,"%lf",&momn[l][index]);
		    st_tmp.momn[l] = momn[l][index];
		}
		pres[index] = EosPressure(&st_tmp);
	    }
	    break;
	case 3:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    for (k = 0; k <= top_gmax[2]; ++k)
	    {
		index = d_index3d(i,j,k,top_gmax);
		comp = top_comp[index];
		st_tmp.eos = &(eos[comp]);

	    	fscanf(infile,"%lf",&dens[index]);
	    	fscanf(infile,"%lf",&engy[index]);
		st_tmp.dens = dens[index];
		st_tmp.engy = engy[index];
		for (l = 0; l < dim; ++l)
		{
	    	    fscanf(infile,"%lf",&momn[l][index]);
		    st_tmp.momn[l] = momn[l][index];
		}
		pres[index] = EosPressure(&st_tmp);
	    }
	}
	fclose(infile);
	scatMeshStates();
	copyMeshStates();
}


void G_CARTESIAN::setAdvectionDt()
{
	double d = (double)dim;
	pp_global_max(&max_speed,1);
	if (max_speed != 0.0)
	    max_dt = hmin/max_speed/d;
	else
	    max_dt = 0.0;
	if (debugging("trace"))
	    printf("In setAdvectionDt: max_dt = %24.18g\n",max_dt);
}	/* end setAdvectionDt */


void G_CARTESIAN::augmentMovieVariables()
{
	int i;
	static HDF_MOVIE_VAR *hdf_movie_var;
	int offset,num_var;

	hdf_movie_var = front->hdf_movie_var;
	offset = front->hdf_movie_var->num_var;
	if (hdf_movie_var == NULL)
	    return initMovieVariables();
	else
	{
	    num_var = offset + dim + 3;
	    FT_ScalarMemoryAlloc((POINTER*)&hdf_movie_var,
				sizeof(HDF_MOVIE_VAR));
	    hdf_movie_var->num_var = num_var;
	    FT_MatrixMemoryAlloc((POINTER*)&hdf_movie_var->var_name,
				num_var,100,sizeof(char));
	    FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->top_var,
				num_var,sizeof(double*));
	    FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->obstacle_comp,
				num_var,sizeof(COMPONENT));
	    FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->preset_bound,
				num_var,sizeof(boolean));
	    FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->var_min,
				num_var,sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->var_max,
				num_var,sizeof(double));
	    for (i = 0; i < front->hdf_movie_var->num_var; ++i)
	    {
		strcpy(hdf_movie_var->var_name[i],
			front->hdf_movie_var->var_name[i]);
		hdf_movie_var->get_state_var[i] =
			front->hdf_movie_var->get_state_var[i];
		hdf_movie_var->top_var[i] = 
			front->hdf_movie_var->top_var[i];
	    }
	    sprintf(hdf_movie_var->var_name[offset+0],"dens");
	    sprintf(hdf_movie_var->var_name[offset+1],"pres");
	    sprintf(hdf_movie_var->var_name[offset+2],"vort");
	    sprintf(hdf_movie_var->var_name[offset+3],"xvel");
	    sprintf(hdf_movie_var->var_name[offset+4],"yvel");
	    hdf_movie_var->get_state_var[offset+0] = getStateDens;
	    hdf_movie_var->get_state_var[offset+1] = getStatePres;
	    hdf_movie_var->get_state_var[offset+2] = getStateVort;
	    hdf_movie_var->get_state_var[offset+3] = getStateXvel;
	    hdf_movie_var->get_state_var[offset+4] = getStateYvel;
	    if (dim == 3)
	    {
	    	sprintf(hdf_movie_var->var_name[offset+5],"zvel");
	    	hdf_movie_var->get_state_var[offset+5] = getStateZvel;
	    }
	}
	hdf_movie_var->top_var[offset+0] = eqn_params->dens;
	hdf_movie_var->top_var[offset+1] = eqn_params->pres;
	hdf_movie_var->top_var[offset+2] = eqn_params->vort;
	hdf_movie_var->top_var[offset+3] = eqn_params->vel[0];
	hdf_movie_var->top_var[offset+4] = eqn_params->vel[1];
	if (dim == 3)
	    hdf_movie_var->top_var[offset+5] = eqn_params->vel[2];
	FT_FreeThese(2,front->hdf_movie_var->var_name,
			front->hdf_movie_var->top_var);
	FT_FreeThese(1,front->hdf_movie_var);
	front->hdf_movie_var = hdf_movie_var;
	front->hdf_movie_var->num_var = num_var;

}	/* end augmentMovieVariables */

void G_CARTESIAN::initMovieVariables()
{
	static HDF_MOVIE_VAR *hdf_movie_var;
	int n;
	MOVIE_OPTION *movie_option = eqn_params->movie_option;

        if (debugging("trace"))
            (void) printf("Entering initMovieVariable()\n");
	if (hdf_movie_var == NULL)
	{
	    FT_ScalarMemoryAlloc((POINTER*)&hdf_movie_var,
				sizeof(HDF_MOVIE_VAR));
	    if (eqn_params->tracked == NO)
		hdf_movie_var->untracked = YES;
	    switch (dim)
	    {
	    case 1:
		hdf_movie_var->num_var = n = 0;
	    	FT_MatrixMemoryAlloc((POINTER*)&hdf_movie_var->var_name,3,100,
				sizeof(char));
	    	FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->top_var,3,
				sizeof(double*));
	    	FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->obstacle_comp,3,
				sizeof(COMPONENT));
	    	FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->preset_bound,3,
				sizeof(boolean));
	    	FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->var_min,3,
				sizeof(double));
	    	FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->var_max,3,
				sizeof(double));
		if (movie_option->plot_dens)
		{
	    	    sprintf(hdf_movie_var->var_name[n],"density");
	    	    hdf_movie_var->get_state_var[n] = getStateDens;
	    	    hdf_movie_var->top_var[n] = eqn_params->dens;
		    if (movie_option->set_bounds)
		    {
			hdf_movie_var->preset_bound[n] = YES;
			hdf_movie_var->var_min[n] = movie_option->min_dens;
			hdf_movie_var->var_max[n] = movie_option->max_dens;
		    }
		    else hdf_movie_var->preset_bound[n] = NO;
		    hdf_movie_var->num_var = ++n;
		}
		if (movie_option->plot_pres)
		{
	    	    sprintf(hdf_movie_var->var_name[n],"pressure");
	    	    hdf_movie_var->get_state_var[n] = getStatePres;
	    	    hdf_movie_var->top_var[n] = eqn_params->pres;
		    if (movie_option->set_bounds)
		    {
			hdf_movie_var->preset_bound[n] = YES;
			hdf_movie_var->var_min[n] = movie_option->min_pres;
			hdf_movie_var->var_max[n] = movie_option->max_pres;
		    }
		    else hdf_movie_var->preset_bound[n] = NO;
		    hdf_movie_var->num_var = ++n;
		}
		if (movie_option->plot_velo)
		{
	    	    sprintf(hdf_movie_var->var_name[n],"velocity");
	    	    hdf_movie_var->get_state_var[n] = getStateXvel;
	    	    hdf_movie_var->top_var[n] = eqn_params->vel[0];
		    if (movie_option->set_bounds)
		    {
			hdf_movie_var->preset_bound[n] = YES;
			hdf_movie_var->var_min[n] = movie_option->min_velo;
			hdf_movie_var->var_max[n] = movie_option->max_velo;
		    }
		    else hdf_movie_var->preset_bound[n] = NO;
		    hdf_movie_var->num_var = ++n;
		}
		break;
	    case 2:
		hdf_movie_var->num_var = n = 0;
	    	FT_MatrixMemoryAlloc((POINTER*)&hdf_movie_var->var_name,5,100,
					sizeof(char));
	    	FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->top_var,5,
					sizeof(double*));
	    	FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->obstacle_comp,5,
					sizeof(COMPONENT));
	    	FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->preset_bound,5,
				sizeof(boolean));
	    	FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->var_min,5,
				sizeof(double));
	    	FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->var_max,5,
				sizeof(double));
		if (movie_option->plot_dens)
		{
	    	    sprintf(hdf_movie_var->var_name[n],"dens");
	    	    hdf_movie_var->get_state_var[n] = getStateDens;
	    	    hdf_movie_var->top_var[n] = eqn_params->dens;
	    	    hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    if (movie_option->set_bounds)
		    {
			hdf_movie_var->preset_bound[n] = YES;
			hdf_movie_var->var_min[n] = movie_option->min_dens;
			hdf_movie_var->var_max[n] = movie_option->max_dens;
		    }
		    else hdf_movie_var->preset_bound[n] = NO;
		    hdf_movie_var->num_var = ++n;
		}
		if (movie_option->plot_pres)
		{
	    	    sprintf(hdf_movie_var->var_name[n],"pres");
	    	    hdf_movie_var->get_state_var[n] = getStatePres;
	    	    hdf_movie_var->top_var[n] = eqn_params->pres;
	    	    hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    if (movie_option->set_bounds)
		    {
			hdf_movie_var->preset_bound[n] = YES;
			hdf_movie_var->var_min[n] = movie_option->min_pres;
			hdf_movie_var->var_max[n] = movie_option->max_pres;
		    }
		    else hdf_movie_var->preset_bound[n] = NO;
		    hdf_movie_var->num_var = ++n;
		}
		if (movie_option->plot_vort)
		{
	    	    sprintf(hdf_movie_var->var_name[n],"vort");
	    	    hdf_movie_var->get_state_var[n] = getStateVort;
	    	    hdf_movie_var->top_var[n] = eqn_params->vort;
	    	    hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    hdf_movie_var->num_var = ++n;
		}
		if (movie_option->plot_velo)
		{
	    	    sprintf(hdf_movie_var->var_name[n],"xvel");
	    	    hdf_movie_var->get_state_var[n] = getStateXvel;
	    	    hdf_movie_var->top_var[n] = eqn_params->vel[0];
	    	    hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    if (movie_option->set_bounds)
		    {
			hdf_movie_var->preset_bound[n] = YES;
			hdf_movie_var->var_min[n] = movie_option->min_velo;
			hdf_movie_var->var_max[n] = movie_option->max_velo;
		    }
		    else hdf_movie_var->preset_bound[n] = NO;
		    hdf_movie_var->num_var = ++n;
	    	    sprintf(hdf_movie_var->var_name[n],"yvel");
	    	    hdf_movie_var->get_state_var[n] = getStateYvel;
	    	    hdf_movie_var->top_var[n] = eqn_params->vel[1];
	    	    hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    if (movie_option->set_bounds)
		    {
			hdf_movie_var->preset_bound[n] = YES;
			hdf_movie_var->var_min[n] = movie_option->min_velo;
			hdf_movie_var->var_max[n] = movie_option->max_velo;
		    }
		    else hdf_movie_var->preset_bound[n] = NO;
		    hdf_movie_var->num_var = ++n;
		}
		break;
	    case 3:
		hdf_movie_var->num_var = n = 0;
	    	FT_MatrixMemoryAlloc((POINTER*)&hdf_movie_var->var_name,15,100,
					sizeof(char));
	    	FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->top_var,15,
					sizeof(double*));
	    	FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->idir,15,
					sizeof(int));
	    	FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->obstacle_comp,15,
					sizeof(COMPONENT));
	    	FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->preset_bound,15,
				sizeof(boolean));
	    	FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->var_min,15,
				sizeof(double));
	    	FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->var_max,15,
				sizeof(double));
		if (movie_option->plot_cross_section[0])
		{
		    if (movie_option->plot_dens)
		    {
	    	    	sprintf(hdf_movie_var->var_name[n],"dens-yz");
	    	    	hdf_movie_var->get_state_var[n] = getStateDens;
	    		hdf_movie_var->top_var[n] = eqn_params->dens;
	    		hdf_movie_var->idir[n] = 0;
	    	    	hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    	hdf_movie_var->num_var = ++n;
		    }
		    if (movie_option->plot_pres)
		    {
	    	    	sprintf(hdf_movie_var->var_name[n],"pres-yz");
	    	    	hdf_movie_var->get_state_var[n] = getStatePres;
	    		hdf_movie_var->top_var[n] = eqn_params->pres;
	    		hdf_movie_var->idir[n] = 0;
	    	    	hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    	hdf_movie_var->num_var = ++n;
		    }
		    if (movie_option->plot_velo)
		    {
	    	    	sprintf(hdf_movie_var->var_name[n],"velo-yz-y");
	    	    	hdf_movie_var->get_state_var[n] = getStateYvel;
	    		hdf_movie_var->top_var[n] = eqn_params->vel[1];
	    		hdf_movie_var->idir[n] = 0;
	    	    	hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    	hdf_movie_var->num_var = ++n;
	    	    	sprintf(hdf_movie_var->var_name[n],"velo-yz-z");
	    	    	hdf_movie_var->get_state_var[n] = getStateZvel;
	    		hdf_movie_var->top_var[n] = eqn_params->vel[2];
	    		hdf_movie_var->idir[n] = 0;
	    	    	hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    	hdf_movie_var->num_var = ++n;
		    }
		}
		if (movie_option->plot_cross_section[1])
		{
		    if (movie_option->plot_dens)
		    {
	    	    	sprintf(hdf_movie_var->var_name[n],"dens-xz");
	    	    	hdf_movie_var->get_state_var[n] = getStateDens;
	    		hdf_movie_var->top_var[n] = eqn_params->dens;
	    		hdf_movie_var->idir[n] = 0;
	    	    	hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    	hdf_movie_var->num_var = ++n;
		    }
		    if (movie_option->plot_pres)
		    {
	    	    	sprintf(hdf_movie_var->var_name[n],"pres-xz");
	    	    	hdf_movie_var->get_state_var[n] = getStatePres;
	    		hdf_movie_var->top_var[n] = eqn_params->pres;
	    		hdf_movie_var->idir[n] = 1;
	    	    	hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    	hdf_movie_var->num_var = ++n;
		    }
		    if (movie_option->plot_velo)
		    {
	    	    	sprintf(hdf_movie_var->var_name[n],"velo-xz-x");
	    	    	hdf_movie_var->get_state_var[n] = getStateXvel;
	    		hdf_movie_var->top_var[n] = eqn_params->vel[0];
	    		hdf_movie_var->idir[n] = 1;
	    	    	hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    	hdf_movie_var->num_var = ++n;
	    	    	sprintf(hdf_movie_var->var_name[n],"velo-xz-z");
	    	    	hdf_movie_var->get_state_var[n] = getStateZvel;
	    		hdf_movie_var->top_var[n] = eqn_params->vel[2];
	    		hdf_movie_var->idir[n] = 1;
	    	    	hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    	hdf_movie_var->num_var = ++n;
		    }
		}
		if (movie_option->plot_cross_section[2])
		{
		    if (movie_option->plot_dens)
		    {
	    	    	sprintf(hdf_movie_var->var_name[n],"dens-xy");
	    	    	hdf_movie_var->get_state_var[n] = getStateDens;
	    		hdf_movie_var->top_var[n] = eqn_params->dens;
	    		hdf_movie_var->idir[n] = 0;
	    	    	hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    	hdf_movie_var->num_var = ++n;
		    }
		    if (movie_option->plot_pres)
		    {
	    	    	sprintf(hdf_movie_var->var_name[n],"pres-xy");
	    	    	hdf_movie_var->get_state_var[n] = getStatePres;
	    		hdf_movie_var->top_var[n] = eqn_params->pres;
	    		hdf_movie_var->idir[n] = 2;
	    	    	hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    	hdf_movie_var->num_var = ++n;
		    }
		    if (movie_option->plot_velo)
		    {
	    	    	sprintf(hdf_movie_var->var_name[n],"velo-xy-x");
	    	    	hdf_movie_var->get_state_var[n] = getStateXvel;
	    		hdf_movie_var->top_var[n] = eqn_params->vel[0];
	    		hdf_movie_var->idir[n] = 2;
	    	    	hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    	hdf_movie_var->num_var = ++n;
	    	    	sprintf(hdf_movie_var->var_name[n],"velo-xy-y");
	    	    	hdf_movie_var->get_state_var[n] = getStateYvel;
	    		hdf_movie_var->top_var[n] = eqn_params->vel[1];
	    		hdf_movie_var->idir[n] = 2;
	    	    	hdf_movie_var->obstacle_comp[n] = SOLID_COMP;
		    	hdf_movie_var->num_var = ++n;
		    }
		}
	    }
	}
	front->hdf_movie_var = hdf_movie_var;
        if (debugging("trace"))
            (void) printf("Leaving initMovieVariable()\n");
}	/* end initMovieVariables */

double G_CARTESIAN::getVorticityX(int i, int j, int k)
{
	int index0,index00,index01,index10,index11;
	double v00,v01,v10,v11;
	double dy,dz;
	double vorticity;
	double *dens = field.dens;
	double **momn = field.momn;

	dy = top_h[1];
	dz = top_h[2];
	index0 = d_index3d(i,j,k,top_gmax);
	index00 = d_index3d(i,j-1,k,top_gmax);
	index01 = d_index3d(i,j+1,k,top_gmax);
	index10 = d_index3d(i,j,k-1,top_gmax);
	index11 = d_index3d(i,j,k+1,top_gmax);
	v00 = -momn[2][index00]/dens[index00];
	v01 =  momn[2][index01]/dens[index01];
	v10 =  momn[1][index10]/dens[index10];
	v11 = -momn[1][index11]/dens[index11];

	vorticity = (v00 + v01)/2.0/dz + (v10 + v11)/2.0/dy;
	return vorticity;
}	/* end getVorticityX */

double G_CARTESIAN::getVorticityY(int i, int j, int k)
{
	int index0,index00,index01,index10,index11;
	double v00,v01,v10,v11;
	double dx,dz;
	double vorticity;
	double *dens = field.dens;
	double **momn = field.momn;

	dx = top_h[0];
	dz = top_h[2];
	index0 = d_index3d(i,j,k,top_gmax);
	index00 = d_index3d(i,j,k-1,top_gmax);
	index01 = d_index3d(i,j,k+1,top_gmax);
	index10 = d_index3d(i-1,j,k,top_gmax);
	index11 = d_index3d(i+1,j,k,top_gmax);
	v00 = -momn[0][index00]/dens[index00];
	v01 =  momn[0][index01]/dens[index01];
	v10 =  momn[2][index10]/dens[index10];
	v11 = -momn[2][index11]/dens[index11];

	vorticity = (v00 + v01)/2.0/dx + (v10 + v11)/2.0/dz;
	return vorticity;
}	/* end getVorticityY */

double G_CARTESIAN::getVorticityZ(int i, int j, int k)
{
	int index0,index00,index01,index10,index11;
	double v00,v01,v10,v11;
	double dx,dy;
	double vorticity;
	double *dens = field.dens;
	double **momn = field.momn;

	dx = top_h[0];
	dy = top_h[1];
	index0 = d_index3d(i,j,k,top_gmax);
	index00 = d_index3d(i-1,j,k,top_gmax);
	index01 = d_index3d(i+1,j,k,top_gmax);
	index10 = d_index3d(i,j-1,k,top_gmax);
	index11 = d_index3d(i,j+1,k,top_gmax);
	v00 = -momn[1][index00]/dens[index00];
	v01 =  momn[1][index01]/dens[index01];
	v10 =  momn[0][index10]/dens[index10];
	v11 = -momn[0][index11]/dens[index11];

	vorticity = (v00 + v01)/2.0/dy + (v10 + v11)/2.0/dx;
	return vorticity;
}	/* end getVorticityZ */

double G_CARTESIAN::getVorticity(int i, int j)
{
	int index0,index00,index01,index10,index11;
	double v00,v01,v10,v11;
	double dx,dy;
	double vorticity;
	double *dens = field.dens;
	double **momn = field.momn;

	dx = top_h[0];
	dy = top_h[1];
	index0 = d_index2d(i,j,top_gmax);
	index00 = d_index2d(i-1,j,top_gmax);
	index01 = d_index2d(i+1,j,top_gmax);
	index10 = d_index2d(i,j-1,top_gmax);
	index11 = d_index2d(i,j+1,top_gmax);
	v00 = -momn[1][index00]/dens[index00];
	v01 =  momn[1][index01]/dens[index01];
	v10 =  momn[0][index10]/dens[index10];
	v11 = -momn[0][index11]/dens[index11];

	vorticity = (v00 + v01)/2.0/dy + (v10 + v11)/2.0/dx;
	return vorticity;
}	/* end getVorticity */

void G_CARTESIAN::copyMeshStates()
{
	int i,j,k,l,index;
	double **vel = eqn_params->vel;
	double **mom = eqn_params->mom;
	double *dens = eqn_params->dens;
	double *pres = eqn_params->pres;
	double *engy = eqn_params->engy;
	double *vort = eqn_params->vort;

	switch (dim)
	{
	case 1:
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index1d(i,top_gmax);
		for (l = 0; l < dim; ++l)
		    vel[l][index] = mom[l][index]/dens[index];	
	    }
	    FT_ParallelExchGridArrayBuffer(dens,front);
	    FT_ParallelExchGridArrayBuffer(pres,front);
	    FT_ParallelExchGridArrayBuffer(engy,front);
	    FT_ParallelExchGridArrayBuffer(mom[0],front);
	    FT_ParallelExchGridArrayBuffer(vel[0],front);
	    break;
	case 2:
	    for (i = imin[0]; i <= imax[0]; ++i)
	    for (j = imin[1]; j <= imax[1]; ++j)
	    {
		index = d_index2d(i,j,top_gmax);
		for (l = 0; l < dim; ++l)
		    vel[l][index] = mom[l][index]/dens[index];	
	    }
	    FT_ParallelExchGridArrayBuffer(dens,front);
	    FT_ParallelExchGridArrayBuffer(pres,front);
	    FT_ParallelExchGridArrayBuffer(engy,front);
	    FT_ParallelExchGridArrayBuffer(vort,front);
	    for (l = 0; l < dim; ++l)
	    {
	    	FT_ParallelExchGridArrayBuffer(mom[l],front);
	    	FT_ParallelExchGridArrayBuffer(vel[l],front);
	    }
	    break;
	case 3:
	    for (i = imin[0]; i <= imax[0]; ++i)
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (k = imin[2]; k <= imax[2]; ++k)
	    {
		index = d_index3d(i,j,k,top_gmax);
		for (l = 0; l < dim; ++l)
		    vel[l][index] = mom[l][index]/dens[index];	
	    }
	    FT_ParallelExchGridArrayBuffer(dens,front);
	    FT_ParallelExchGridArrayBuffer(pres,front);
	    FT_ParallelExchGridArrayBuffer(engy,front);
	    for (l = 0; l < dim; ++l)
	    {
	    	FT_ParallelExchGridArrayBuffer(mom[l],front);
	    	FT_ParallelExchGridArrayBuffer(vel[l],front);
	    }
	    break;
	}
}	/* end copyMeshStates */

void G_CARTESIAN::compSGS(void)
{
        int i,j,k,index,index0,index1,index2,index3,index4,size;  
        double *u, *v;
        double ulx,urx,vlx,vrx;
        double uly,ury,vly,vry;
        double ux,uy,vx,vy;
        double s, *s11, *s12, *s22;
        double *ss11, *ss12, *ss22;
        double *tau00, *tau01, *tau10, *tau11;
        double *vel_u, *vel_v, *vel_uu, *vel_uv, *vel_vv;
        double sum_vel_u,sum_vel_v,sum_vel_uu,sum_vel_uv,sum_vel_vv;
        double sum_s11,sum_s12,sum_s22,sum_ss11,sum_ss12,sum_ss22,sum_s;
        double *ma11, *ma12, *la11, *la12, *la22;
        double *cs, *cs_ave, *deno, *nume, *co_coords_y;
        double coords[2];
        int    *r, num_r;
        int    ii,jj,iii,jjj;
        const int nn = pp_numnodes();
        num_r = (int)(((top_U[1]-top_L[1])/top_h[1])+1);
	double **momn = field.momn;

        size = (top_gmax[0]+1)*(top_gmax[1]+1);
        FT_VectorMemoryAlloc((POINTER*)&u,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&v,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&s11,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&s12,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&s22,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&ss11,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&ss12,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&ss22,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&tau00,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&tau01,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&tau10,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&tau11,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&vel_u,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&vel_v,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&vel_uu,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&vel_uv,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&vel_vv,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&co_coords_y,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&ma11,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&ma12,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&la11,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&la12,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&la22,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&r,size,sizeof(int));
        FT_VectorMemoryAlloc((POINTER*)&cs,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&cs_ave,num_r,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&deno,num_r,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&nume,num_r,sizeof(double));

        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index2d(i,j,top_gmax);
            u[index] = momn[0][index];
            v[index] = momn[1][index];
            getRectangleCenter(index, coords);
            co_coords_y[index] = coords[1] + (top_h[1]/2.0);
        }

        for (j = imin[1]; j <= imax[1]; j++)
        for (i = imin[0]-2; i <= imax[0]+2; i++)
        {
            index0  = d_index2d(i,j,top_gmax);
            index1  = d_index2d(i-1,j,top_gmax);
            ulx = u[index1];
            uly = v[index1];
            index2  = d_index2d(i+1,j,top_gmax);
            urx = u[index2];
            ury = v[index2];
            index3  = d_index2d(i,j-1,top_gmax);
            vlx = u[index3];
            vly = v[index3];
            index4  = d_index2d(i,j+1,top_gmax);
            vrx = u[index4];
            vry = v[index4];

            ux = (urx - ulx) / (2.0*top_h[0]);
            uy = (ury - uly) / (2.0*top_h[1]);
            vx = (vrx - vlx) / (2.0*top_h[0]);
            vy = (vry - vly) / (2.0*top_h[1]);
            s11[index0] = ux;
            s12[index0] = (uy + vx)/2;
            s22[index0] = vy;
            s = sqrt(2*( (s11[index0]*s11[index0]) 
                         + (2*(s12[index0]*s12[index0]))
                         + (s22[index0]*s22[index0])));
            ss11[index0] = s*s11[index0];
            ss12[index0] = s*s12[index0];
            ss22[index0] = s*s22[index0];
            vel_u[index0] = u[index0];
            vel_v[index0] = v[index0];  
            vel_uu[index0] = u[index0]*u[index0]; 
            vel_uv[index0] = u[index0]*v[index0];  
            vel_vv[index0] = v[index0]*v[index0];      
        }

        for (j = imin[1]; j <= (imax[1]/2); j++)
        {
            jj = (2*j)-1;
            for (i = imin[0]-1; i <= (imax[0]/2)+1; i++)
            {
                ii = (2*i)-4;
                index = d_index2d(ii,jj,top_gmax);
                sum_vel_u = sum_vel_v = 0.0;
                sum_vel_uu = sum_vel_uv = sum_vel_vv = 0.0;
                sum_s11 = sum_s12 = sum_s22 = 0.0;
                sum_ss11 = sum_ss12 = sum_ss22 = sum_s = 0.0;
                for(jjj = jj; jjj < jj+2; jjj++)
                for(iii = ii; iii < ii+2; iii++)
                {
                    index0  = d_index2d(iii,jjj,top_gmax);
                    sum_vel_u += vel_u[index0];
                    sum_vel_v += vel_v[index0];
                    sum_vel_uu += vel_uu[index0];
                    sum_vel_uv += vel_uv[index0];
                    sum_vel_vv += vel_vv[index0];
                    sum_s11 += s11[index0];
                    sum_s12 += s12[index0];
                    sum_s22 += s22[index0];
                    sum_ss11 += ss11[index0];
                    sum_ss12 += ss12[index0];
                    sum_ss22 += ss22[index0];
                    sum_s += sqrt(2*( (s11[index0]*s11[index0]) 
                                  + (2*(s12[index0]*s12[index0])) 
                                  + (s22[index0]*s22[index0])));
                } 
                ma11[index] = (2.0*top_h[1]*top_h[1]*(sum_ss11/4.0))
                        - (2.0*4*top_h[1]*top_h[1]*(sum_s/4.0)*(sum_s11/4.0));
                ma12[index] = (2.0*top_h[1]*top_h[1]*(sum_ss12/4.0))
                        - (2.0*4*top_h[1]*top_h[1]*(sum_s/4.0)*(sum_s12/4.0));
                la11[index] = (sum_vel_uu/4.0)-((sum_vel_u/4.0)*
			(sum_vel_u/4.0));
                la12[index] = (sum_vel_uv/4.0)-((sum_vel_u/4.0)*
			(sum_vel_v/4.0));
                la12[index] = (sum_vel_vv/4.0)-((sum_vel_v/4.0)*
			(sum_vel_v/4.0));
                r[index] = (int)(co_coords_y[index]/(2*top_h[1]));
            }
        }

        for (k = 0; k < num_r; k++)
        {
            deno[k] = 0.0;
            nume[k] = 0.0;
        }

        for (k = 0; k < num_r; k++)
        for (j = imin[1]; j <= (imax[1]/2); j++)
        {
            jj = (2*j)-1;
            for (i = imin[0]-1; i <= (imax[0]/2)+1; i++)
            {
                ii = (2*i)-4;
                index0 = d_index2d(ii,jj,top_gmax);
                if(k == r[index0])
                {
                    deno[k] += (ma11[index0]*ma11[index0]) + 
				(ma12[index0]*ma12[index0]);
                    nume[k] += (((la11[index0]/2.0)-(la22[index0]/2.0))*
				ma11[index0]) + (la12[index0]*ma12[index0]);
                }
            }
        }

        pp_gsync();
        
        if (nn > 1)
        {
           for (k = 0; k < num_r; k++)
           {
              pp_global_sum(&deno[k],1L);
              pp_global_sum(&nume[k],1L);
           }
        }

        for (k = 0; k < num_r; k++)
        {
            if(deno[k] < 10e-16)
                cs_ave[k] = 0.0;
            else
                cs_ave[k] = nume[k]/deno[k];
        }

        for (j = imin[1]; j <= (imax[1]/2); j++)
        {
            jj = (2*j)-1;
            for (i = imin[0]-1; i <= (imax[0]/2)+1; i++)
            {
                ii = (2*i)-4;
                index = d_index2d(ii,jj,top_gmax);
                for(jjj = jj; jjj < jj+2; jjj++)
                for(iii = ii; iii < ii+2; iii++)
                {
                    index0 = d_index2d(iii,jjj,top_gmax);
                    cs[index0] = cs_ave[r[index]];
                }
            }
        }

        for (j = imin[1]; j <= imax[1]; j++)
        for (i = imin[0]-1; i <= imax[0]+1; i++)
        {
            index0  = d_index2d(i,j,top_gmax);
            s = sqrt(2*( (s11[index0]*s11[index0]) 
                          + (2*(s12[index0]*s12[index0]))
                          + (s22[index0]*s22[index0])));
            tau00[index0] = - 2.0*cs[index0]*top_h[0]*top_h[0]*
                              s*((s11[index0]/2.0)-(s22[index0]/2.0));
            tau01[index0] = - 2.0*cs[index0]*top_h[0]*top_h[0]*s*(s12[index0]);
            tau10[index0] = - 2.0*cs[index0]*top_h[0]*top_h[0]*s*(s12[index0]);
            tau11[index0] = - 2.0*cs[index0]*top_h[0]*top_h[0]*
                              s*((s22[index0]/2.0)-(s11[index0]/2.0));
        }

        for (j = imin[1]; j <= imax[1]; j++)
        for (i = imin[0]; i <= imax[0]; i++)
        {
            index0 = d_index2d(i,j,top_gmax);
            index1 = d_index2d(i-1,j,top_gmax);
            index2 = d_index2d(i+1,j,top_gmax);
            index3 = d_index2d(i,j-1,top_gmax);
            index4 = d_index2d(i,j+1,top_gmax);

            momn[0][index0] += -m_dt*(
                              ((tau00[index2]-tau00[index1])/(2.0*top_h[0])) + 
                                ((tau01[index4]-tau01[index3])/(2.0*top_h[1])));
            momn[1][index0] += -m_dt*(
                              ((tau10[index2]-tau10[index1])/(2.0*top_h[0])) + 
                              ((tau11[index4]-tau11[index3])/(2.0*top_h[1])));
        }
        FT_FreeThese(2,u,v);
        FT_FreeThese(4,tau00,tau01,tau10,tau11);
        FT_FreeThese(6,s11,s12,s22,ss11,ss12,ss22);
        FT_FreeThese(5,vel_u,vel_v,vel_uu,vel_uv,vel_vv);
        FT_FreeThese(11,co_coords_y,ma11,ma12,la11,la12,la22,r,cs,cs_ave,
					deno,nume);
}       /* end compSGS */

void G_CARTESIAN::sampleVelocity()
{
	switch (dim)
	{
	case 2:
	    return sampleVelocity2d();
	case 3:
	    return sampleVelocity3d();
	}
}	/* end sampleVelocity */

void G_CARTESIAN::sampleVelocity3d()
{
        int i,j,k,index;
        double coords[MAXD];
        double velo1,velo2,velo_tmp1,velo_tmp2,velo;
        FILE *sfile;
        char sname[100];
        static int count = 0;
        static int step = 0;
        static int l=-1,m=-1;
        static double lambda1,lambda2;
	SAMPLE *sample = front->sample;
	char *sample_type = sample->sample_type;
	double *sample_line = sample->sample_coords;
	char *out_name = front-> out_name;
	double dens;

	if (front->step < sample->start_step || front->step > sample->end_step)
	    return;
	if ((front->step - sample->start_step)%sample->step_interval)
	    return;
        if (step != front->step)
        {
            step = front->step;
            count = 0;
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
                    getRectangleCenter(index, coords);
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
                    sprintf(sname, "%s/x-%d-%d.xg",out_name,step,count);
                    sfile = fopen(sname,"w");
                    for (k = imin[2]; k <= imax[2]; ++k)
                    {
                        index = d_index3d(i,j,k,top_gmax);
                        dens = field.dens[index];
                        velo1 = field.momn[0][index]/dens;
                        index = d_index3d(i+1,j,k,top_gmax);
                        dens = field.dens[index];
                        velo2 = field.momn[0][index]/dens;
                        velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        index = d_index3d(i,j+1,k,top_gmax);
                        dens = field.dens[index];
                        velo1 = field.momn[0][index]/dens;
                        index = d_index3d(i+1,j+1,k,top_gmax);
                        dens = field.dens[index];
                        velo2 = field.momn[0][index]/dens;
                        velo_tmp2 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        velo = (velo_tmp1 + lambda2*velo_tmp2)/(1.0 + lambda2);
                        getRectangleCenter(index,coords);
                        fprintf(sfile,"%20.14f   %20.14f\n",coords[2],velo);
                    }
                    fclose(sfile);

                    sprintf(sname,"%s/y-%d-%d.xg",out_name,step,count);
                    sfile = fopen(sname,"w");
                    for (k = imin[2]; k <= imax[2]; ++k)
                    {
                        index = d_index3d(i,j,k,top_gmax);
                        dens = field.dens[index];
                        velo1 = field.momn[1][index]/dens;
                        index = d_index3d(i+1,j,k,top_gmax);
                        dens = field.dens[index];
                        velo2 = field.momn[1][index]/dens;
                        velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        index = d_index3d(i,j+1,k,top_gmax);
                        dens = field.dens[index];
                        velo1 = field.momn[1][index]/dens;
                        index = d_index3d(i+1,j+1,k,top_gmax);
                        dens = field.dens[index];
                        velo2 = field.momn[1][index]/dens;
                        velo_tmp2 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        velo = (velo_tmp1 + lambda2*velo_tmp2)/(1.0 + lambda2);
                        getRectangleCenter(index,coords);
                        fprintf(sfile,"%20.14f   %20.14f\n",coords[2],velo);
                    }
                    fclose(sfile);

                    sprintf(sname,"%s/z-%d-%d.xg",out_name,step,count++);
                    sfile = fopen(sname,"w");
                    for (k = imin[2]; k <= imax[2]; ++k)
                    {
                        index = d_index3d(i,j,k,top_gmax);
                        dens = field.dens[index];
                        velo1 = field.momn[2][index]/dens;
                        index = d_index3d(i+1,j,k,top_gmax);
                        dens = field.dens[index];
                        velo2 = field.momn[2][index]/dens;
                        velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        index = d_index3d(i,j+1,k,top_gmax);
                        dens = field.dens[index];
                        velo1 = field.momn[2][index]/dens;
                        index = d_index3d(i+1,j+1,k,top_gmax);
                        dens = field.dens[index];
                        velo2 = field.momn[2][index]/dens;
                        velo_tmp2 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        velo = (velo_tmp1 + lambda2*velo_tmp2)/(1.0 + lambda2);
                        getRectangleCenter(index,coords);
                        fprintf(sfile,"%20.14f   %20.14f\n",coords[2],velo);
                    }
                    fclose(sfile);

                    printf("sample line: x = %20.14f, y = %20.14f\n",coords[0],
                        coords[1]);

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
                    sprintf(sname, "%s/x-%d-%d.xg",out_name,step,count);
                    sfile = fopen(sname,"w");
                    for (j = imin[1]; j <= imax[1]; ++j)
                    {
                        index = d_index3d(i,j,k,top_gmax);
                        dens = field.dens[index];
                        velo1 = field.momn[0][index]/dens;
                        index = d_index3d(i+1,j,k,top_gmax);
                        dens = field.dens[index];
                        velo2 = field.momn[0][index]/dens;
                        velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        index = d_index3d(i,j,k+1,top_gmax);
                        dens = field.dens[index];
                        velo1 = field.momn[0][index]/dens;
                        index = d_index3d(i+1,j,k+1,top_gmax);
                        dens = field.dens[index];
                        velo2 = field.momn[0][index]/dens;
                        velo_tmp2 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        velo = (velo_tmp1 + lambda2*velo_tmp2)/(1.0 + lambda2);
                        getRectangleCenter(index,coords);
                        fprintf(sfile,"%20.14f   %20.14f\n",coords[1],velo);
                    }
                    fclose(sfile);

                    sprintf(sname,"%s/y-%d-%d.xg",out_name,step,count);
                    sfile = fopen(sname,"w");
                    for (j = imin[1]; j <= imax[1]; ++j)
                    {
                        index = d_index3d(i,j,k,top_gmax);
                        dens = field.dens[index];
                        velo1 = field.momn[1][index]/dens;
                        index = d_index3d(i+1,j,k,top_gmax);
                        dens = field.dens[index];
                        velo2 = field.momn[1][index]/dens;
                        velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        index = d_index3d(i,j,k+1,top_gmax);
                        dens = field.dens[index];
                        velo1 = field.momn[1][index]/dens;
                        index = d_index3d(i+1,j,k+1,top_gmax);
                        dens = field.dens[index];
                        velo2 = field.momn[1][index]/dens;
                        velo_tmp2 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        velo = (velo_tmp1 + lambda2*velo_tmp2)/(1.0 + lambda2);
                        getRectangleCenter(index,coords);
                        fprintf(sfile,"%20.14f   %20.14f\n",coords[1],velo);
                    }
                    fclose(sfile);

                    sprintf(sname,"%s/z-%d-%d.xg",out_name,step,count++);
                    sfile = fopen(sname,"w");
                    for (j = imin[1]; j <= imax[1]; ++j)
                    {
                        index = d_index3d(i,j,k,top_gmax);
                        dens = field.dens[index];
                        velo1 = field.momn[2][index]/dens;
                        index = d_index3d(i+1,j,k,top_gmax);
                        dens = field.dens[index];
                        velo2 = field.momn[2][index]/dens;
                        velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        index = d_index3d(i,j,k+1,top_gmax);
                        dens = field.dens[index];
                        velo1 = field.momn[2][index]/dens;
                        index = d_index3d(i+1,j,k+1,top_gmax);
                        dens = field.dens[index];
                        velo2 = field.momn[2][index]/dens;
                        velo_tmp2 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        velo = (velo_tmp1 + lambda2*velo_tmp2)/(1.0 + lambda2);
                        getRectangleCenter(index,coords);
                        fprintf(sfile,"%20.14f   %20.14f\n",coords[1],velo);
                    }
                    fclose(sfile);

                    printf("sample line: x = %20.14f, z = %20.14f\n",coords[0],
                        coords[2]);

                    break;

                    default:
                        printf("Incorrect input for sample velocity!\n");
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
                    sprintf(sname, "%s/x-%d-%d.xg",out_name,step,count);
                    sfile = fopen(sname,"w");
                    for (i = imin[0]; i <= imax[0]; ++i)
                    {
                        index = d_index3d(i,j,k,top_gmax);
                        dens = field.dens[index];
                        velo1 = field.momn[0][index]/dens;
                        index = d_index3d(i,j+1,k,top_gmax);
                        dens = field.dens[index];
                        velo2 = field.momn[0][index]/dens;
                        velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        index = d_index3d(i,j,k+1,top_gmax);
                        dens = field.dens[index];
                        velo1 = field.momn[0][index]/dens;
                        index = d_index3d(i,j+1,k+1,top_gmax);
                        dens = field.dens[index];
                        velo2 = field.momn[0][index]/dens;
                        velo_tmp2 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        velo = (velo_tmp1 + lambda2*velo_tmp2)/(1.0 + lambda2);
                        getRectangleCenter(index,coords);
                        fprintf(sfile,"%20.14f   %20.14f\n",coords[0],velo);
                    }
                    fclose(sfile);

                    sprintf(sname, "%s/y-%d-%d.xg",out_name,step,count);
                    sfile = fopen(sname,"w");
                    for (i = imin[0]; i <= imax[0]; ++i)
                    {
                        index = d_index3d(i,j,k,top_gmax);
                        dens = field.dens[index];
                        velo1 = field.momn[1][index]/dens;
                        index = d_index3d(i,j+1,k,top_gmax);
                        dens = field.dens[index];
                        velo2 = field.momn[1][index]/dens;
                        velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        index = d_index3d(i,j,k+1,top_gmax);
                        dens = field.dens[index];
                        velo1 = field.momn[1][index]/dens;
                        index = d_index3d(i,j+1,k+1,top_gmax);
                        dens = field.dens[index];
                        velo2 = field.momn[1][index]/dens;
                        velo_tmp2 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        velo = (velo_tmp1 + lambda2*velo_tmp2)/(1.0 + lambda2);
                        getRectangleCenter(index,coords);
                        fprintf(sfile,"%20.14f   %20.14f\n",coords[0],velo);
                    }
                    fclose(sfile);

                    sprintf(sname, "%s/z-%d-%d.xg",out_name,step,count++);
                    sfile = fopen(sname,"w");
                    for (i = imin[0]; i <= imax[0]; ++i)
                    {
                        index = d_index3d(i,j,k,top_gmax);
                        dens = field.dens[index];
                        velo1 = field.momn[2][index]/dens;
                        index = d_index3d(i,j+1,k,top_gmax);
                        dens = field.dens[index];
                        velo2 = field.momn[2][index]/dens;
                        velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        index = d_index3d(i,j,k+1,top_gmax);
                        dens = field.dens[index];
                        velo1 = field.momn[2][index]/dens;
                        index = d_index3d(i,j+1,k+1,top_gmax);
                        dens = field.dens[index];
                        velo2 = field.momn[2][index]/dens;
                        velo_tmp2 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        velo = (velo_tmp1 + lambda2*velo_tmp2)/(1.0 + lambda2);
                        getRectangleCenter(index,coords);
                        fprintf(sfile,"%20.14f   %20.14f\n",coords[0],velo);
                    }
                    fclose(sfile);

                    printf("sample line: y = %20.14f, z = %20.14f\n",coords[1],
                        coords[2]);

                    break;

                default:
                    printf("Incorrect input for sample velocity!\n");
                    break;
            }
        default:
            printf("Incorrect input for sample velocity!\n");
            break;
        }
}	/* end sampleVelocity3d */

void G_CARTESIAN::sampleVelocity2d()
{
	int i,j,index;
	SAMPLE *sample = front->sample;
        char *sample_type = sample->sample_type;
        double *line = sample->sample_coords;
        char *out_name = front->out_name;
        double coords[MAXD];
        double velo1,velo2,velo;
        FILE *sfile;
        char sname[100];
        static int count = 0;
        static int step = 0;
        static int l = -1;
        static double lambda;
	double dens;

	if (front->step < sample->start_step || front->step > sample->end_step)
            return;
        if ((front->step - sample->start_step)%sample->step_interval)
            return;
        if (step != front->step)
            step = front->step;
	
        switch (sample_type[0])
        {
        case 'x':
            sprintf(sname, "%s/vertical-x-%d-%d.xg",out_name,step,count);
            sfile = fopen(sname,"w");
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
            for (j = imin[1]; j <= imax[1]; ++j)
            {
                index = d_index2d(i,j,top_gmax);
		dens = field.dens[index];
                velo1 = field.momn[0][index]/dens;
                index = d_index2d(i+1,j,top_gmax);
		dens = field.dens[index];
                velo2 = field.momn[0][index]/dens;
                velo = (velo1 + lambda*velo2) / (1.0 + lambda);
                getRectangleCenter(index,coords);
                fprintf(sfile,"%20.14f   %20.14f\n",coords[1],velo);
            }
            fclose(sfile);
            sprintf(sname,"%s/vertical-y-%d-%d.xg",out_name,step,count++);
            sfile = fopen(sname,"w");
            for (j = imin[1]; j <= imax[1]; ++j)
            {
                index = d_index2d(i,j,top_gmax);
                dens = field.dens[index];
                velo1 = field.momn[1][index]/dens;
                index = d_index2d(i+1,j,top_gmax);
                dens = field.dens[index];
                velo2 = field.momn[1][index]/dens;
                velo = (velo1 + lambda*velo2) / (1.0 + lambda);
                getRectangleCenter(index,coords);
                fprintf(sfile,"%20.14f   %20.14f\n",coords[1],velo);
            }
            fclose(sfile);
            break;
        case 'y':
            sprintf(sname, "%s/horizontal-x-%d-%d.xg",out_name,step,count);
            sfile = fopen(sname,"w");
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
            for (i = imin[0]; i <= imax[0]; ++i)
            {
                index = d_index2d(i,j,top_gmax);
                dens = field.dens[index];
                velo1 = field.momn[0][index]/dens;
                index = d_index2d(i,j+1,top_gmax);
                dens = field.dens[index];
                velo2 = field.momn[0][index]/dens;
                velo = (velo1 + lambda*velo2) / (1.0 + lambda);
                getRectangleCenter(index,coords);
                fprintf(sfile,"%20.14f   %20.14f\n",coords[0],velo);
            }
            fclose(sfile);
            sprintf(sname,"%s/horizontal-y-%d-%d.xg",out_name,step,count++);
            sfile = fopen(sname,"w");
            for (i = imin[0]; i <= imax[0]; ++i)
            {
                index = d_index2d(i,j,top_gmax);
                dens = field.dens[index];
                velo1 = field.momn[1][index]/dens;
                index = d_index2d(i,j+1,top_gmax);
                dens = field.dens[index];
                velo2 = field.momn[1][index]/dens;
                velo = (velo1 + lambda*velo2) / (1.0 + lambda);
                getRectangleCenter(index,coords);
                fprintf(sfile,"%20.14f   %20.14f\n",coords[0],velo);
            }
            fclose(sfile);
            break;
        }
}	/* end sampleVelocity2d */

void G_CARTESIAN::numericalFlux(
	POINTER scheme_params,
	SWEEP *sweep,
	FSWEEP *fsweep,
	int n)
{
	switch (eqn_params->num_scheme)
	{
	case TVD_FIRST_ORDER:
	case TVD_SECOND_ORDER:
	case TVD_FOURTH_ORDER:
	    TVD_flux(scheme_params,sweep,fsweep,n);
	    break;
	case WENO_FIRST_ORDER:
	case WENO_SECOND_ORDER:
	case WENO_FOURTH_ORDER:
	    WENO_flux(scheme_params,sweep,fsweep,n);
	    break;
	default:
	    (void) printf("Unknow numerical scheme\n");
	    clean_up(ERROR);
	}
}	/* numericalFlux */


void G_CARTESIAN::scatMeshVst(SWEEP *m_vst)
{
	int i,j,k,l,index;

	switch (dim)
	{
	case 1:
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index1d(i,top_gmax);
		array[index] = m_vst->dens[index];
	    }
	    scatMeshArray();
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index1d(i,top_gmax);
		m_vst->dens[index] = array[index];
	    }
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index1d(i,top_gmax);
		array[index] = m_vst->engy[index];
	    }
	    scatMeshArray();
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index1d(i,top_gmax);
		m_vst->engy[index] = array[index];
	    }
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index1d(i,top_gmax);
		array[index] = m_vst->pres[index];
	    }
	    scatMeshArray();
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index1d(i,top_gmax);
		m_vst->pres[index] = array[index];
	    }
	    for (l = 0; l < dim; ++l)
	    {
	    	for (i = imin[0]; i <= imax[0]; ++i)
	    	{
		    index = d_index1d(i,top_gmax);
		    array[index] = m_vst->momn[l][index];
	    	}
	    	scatMeshArray();
            	for (i = 0; i <= top_gmax[0]; i++)
            	{
                    index = d_index1d(i,top_gmax);
		    m_vst->momn[l][index] = array[index];
	    	}
	    }
	    break;
	case 2:
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index2d(i,j,top_gmax);
		array[index] = m_vst->dens[index];
	    }
	    scatMeshArray();
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index2d(i,j,top_gmax);
		m_vst->dens[index] = array[index];
	    }
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index2d(i,j,top_gmax);
		array[index] = m_vst->engy[index];
	    }
	    scatMeshArray();
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index2d(i,j,top_gmax);
		m_vst->engy[index] = array[index];
	    }
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index2d(i,j,top_gmax);
		array[index] = m_vst->pres[index];
	    }
	    scatMeshArray();
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index2d(i,j,top_gmax);
		m_vst->pres[index] = array[index];
	    }
	    for (l = 0; l < dim; ++l)
	    {
	    	for (j = imin[1]; j <= imax[1]; ++j)
	    	for (i = imin[0]; i <= imax[0]; ++i)
	    	{
		    index = d_index2d(i,j,top_gmax);
		    array[index] = m_vst->momn[l][index];
	    	}
	    	scatMeshArray();
            	for (j = 0; j <= top_gmax[1]; j++)
            	for (i = 0; i <= top_gmax[0]; i++)
            	{
                    index = d_index2d(i,j,top_gmax);
		    m_vst->momn[l][index] = array[index];
	    	}
	    }
	    break;
	case 3:
	    for (k = imin[2]; k <= imax[2]; ++k)
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index3d(i,j,k,top_gmax);
                array[index] = m_vst->dens[index];
	    }
	    scatMeshArray();
            for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
		index = d_index3d(i,j,k,top_gmax);
                m_vst->dens[index] = array[index];
	    }
	    for (k = imin[2]; k <= imax[2]; ++k)
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index3d(i,j,k,top_gmax);
                array[index] = m_vst->engy[index];
	    }
	    scatMeshArray();
            for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index3d(i,j,k,top_gmax);
                m_vst->engy[index] = array[index];
	    }
	    for (k = imin[2]; k <= imax[2]; ++k)
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index3d(i,j,k,top_gmax);
                array[index] = m_vst->pres[index];
	    }
	    scatMeshArray();
            for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index3d(i,j,k,top_gmax);
                m_vst->pres[index] = array[index];
	    }
	    for (l = 0; l < dim; ++l)
	    {
	    	for (k = imin[2]; k <= imax[2]; ++k)
	    	for (j = imin[1]; j <= imax[1]; ++j)
	    	for (i = imin[0]; i <= imax[0]; ++i)
	    	{
		    index = d_index3d(i,j,k,top_gmax);
                    array[index] = m_vst->momn[l][index];
	    	}
	    	scatMeshArray();
            	for (k = 0; k <= top_gmax[2]; k++)
            	for (j = 0; j <= top_gmax[1]; j++)
            	for (i = 0; i <= top_gmax[0]; i++)
            	{
		    index = d_index3d(i,j,k,top_gmax);
                    m_vst->momn[l][index] = array[index];
	    	}
	    }
	}
}	/* end scatMeshStates */

void G_CARTESIAN::copyMeshVst(
	SWEEP m_vst_orig,
	SWEEP *m_vst)
{
	int i,j,k,l,index;
	switch (dim)
	{
	case 1:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    {
		index = d_index1d(i,top_gmax);
		m_vst->dens[index] = m_vst_orig.dens[index];
		m_vst->engy[index] = m_vst_orig.engy[index];
		m_vst->pres[index] = m_vst_orig.pres[index];
		for (l = 0; l < dim; ++l)
		    m_vst->momn[l][index] = m_vst_orig.momn[l][index];
	    }
	    break;
	case 2:
	    for (j = 0; j <= top_gmax[1]; ++j)
	    for (i = 0; i <= top_gmax[0]; ++i)
	    {
		index = d_index2d(i,j,top_gmax);
		m_vst->dens[index] = m_vst_orig.dens[index];
		m_vst->engy[index] = m_vst_orig.engy[index];
		m_vst->pres[index] = m_vst_orig.pres[index];
		for (l = 0; l < dim; ++l)
		    m_vst->momn[l][index] = m_vst_orig.momn[l][index];
	    }
	    break;
	case 3:
	    for (k = 0; k <= top_gmax[2]; ++k)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    for (i = 0; i <= top_gmax[0]; ++i)
	    {
		index = d_index3d(i,j,k,top_gmax);
		m_vst->dens[index] = m_vst_orig.dens[index];
		m_vst->engy[index] = m_vst_orig.engy[index];
		m_vst->pres[index] = m_vst_orig.pres[index];
		for (l = 0; l < dim; ++l)
		    m_vst->momn[l][index] = m_vst_orig.momn[l][index];
	    }
	}
}	/* end copyMeshVst */

void G_CARTESIAN::copyToMeshVst(
	SWEEP *m_vst)
{
	int i,j,k,l,index;
	double *dens = field.dens;
	double *engy = field.engy;
	double *pres = field.pres;
	double **momn = field.momn;
	switch (dim)
	{
	case 1:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    {
		index = d_index1d(i,top_gmax);
		m_vst->dens[index] = dens[index];
		m_vst->engy[index] = engy[index];
		m_vst->pres[index] = pres[index];
		for (l = 0; l < dim; ++l)
		    m_vst->momn[l][index] = momn[l][index];
	    }
	    break;
	case 2:
	    for (j = 0; j <= top_gmax[1]; ++j)
	    for (i = 0; i <= top_gmax[0]; ++i)
	    {
		index = d_index2d(i,j,top_gmax);
		m_vst->dens[index] = dens[index];
		m_vst->engy[index] = engy[index];
		m_vst->pres[index] = pres[index];
		for (l = 0; l < dim; ++l)
		    m_vst->momn[l][index] = momn[l][index];
	    }
	    break;
	case 3:
	    for (k = 0; k <= top_gmax[2]; ++k)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    for (i = 0; i <= top_gmax[0]; ++i)
	    {
		index = d_index3d(i,j,k,top_gmax);
		m_vst->dens[index] = dens[index];
		m_vst->engy[index] = engy[index];
		m_vst->pres[index] = pres[index];
		for (l = 0; l < dim; ++l)
		    m_vst->momn[l][index] = momn[l][index];
	    }
	}
}	/* end copyToMeshVst */

void G_CARTESIAN::copyFromMeshVst(
	SWEEP m_vst)
{
	int i,j,k,l,index;
	STATE state;
	COMPONENT comp;
	double *dens = field.dens;
	double *engy = field.engy;
	double *pres = field.pres;
	double **momn = field.momn;
	
	//GFM
	if(eqn_params->tracked)
	{
	    get_ghost_state(m_vst, 2, 0);
	    get_ghost_state(m_vst, 3, 1);
	    scatMeshGhost();
	}

	state.dim = dim;
	switch (dim)
	{
	case 1:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    {
		index = d_index1d(i,top_gmax);
		comp = top_comp[index];
		state.dens = m_vst.dens[index];
		state.engy = m_vst.engy[index];
		state.pres = m_vst.pres[index];
		for (l = 0; l < dim; ++l)
		    state.momn[l] = m_vst.momn[l][index];
		if (gas_comp(top_comp[index]))
		{
		    state.eos = &(eqn_params->eos[comp]);
		    checkCorrectForTolerance(&state);
		}
		dens[index] = state.dens;
		engy[index] = state.engy;
		pres[index] = state.pres;
		for (l = 0; l < dim; ++l)
		    momn[l][index] = state.momn[l];
	    }
	    break;
	case 2:
	    for (j = 0; j <= top_gmax[1]; ++j)
	    for (i = 0; i <= top_gmax[0]; ++i)
	    {
		index = d_index2d(i,j,top_gmax);
		comp = top_comp[index];
		state.dens = m_vst.dens[index];
		state.engy = m_vst.engy[index];
		state.pres = m_vst.pres[index];
		for (l = 0; l < dim; ++l)
		    state.momn[l] = m_vst.momn[l][index];
		if (gas_comp(top_comp[index]))
		{
		    state.eos = &(eqn_params->eos[comp]);
		    checkCorrectForTolerance(&state);
		}
		dens[index] = state.dens;
		engy[index] = state.engy;
		pres[index] = state.pres;
		for (l = 0; l < dim; ++l)
		    momn[l][index] = state.momn[l];
	    }
	    break;
	case 3:
	    for (k = 0; k <= top_gmax[2]; ++k)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    for (i = 0; i <= top_gmax[0]; ++i)
	    {
		index = d_index3d(i,j,k,top_gmax);
		comp = top_comp[index];
		state.dens = m_vst.dens[index];
		state.engy = m_vst.engy[index];
		state.pres = m_vst.pres[index];
		for (l = 0; l < dim; ++l)
		    state.momn[l] = m_vst.momn[l][index];
		if (gas_comp(top_comp[index]))
		{
		    state.eos = &(eqn_params->eos[comp]);
		    checkCorrectForTolerance(&state);
		}
		dens[index] = state.dens;
		engy[index] = state.engy;
		pres[index] = state.pres;
		for (l = 0; l < dim; ++l)
		    momn[l][index] = state.momn[l];
	    }
	}
}	/* end copyFromMeshVst */

void G_CARTESIAN::appendStencilBuffer2d(
	SWEEP *vst,
	SWEEP *m_vst,
	int i,
	int dir)
{
	int		i1,i2,k,offset,index0,index;
	INTERFACE 	*intfc = front->interf;
	GRID_DIRECTION 	ldir[3] = {WEST,SOUTH,LOWER};
	GRID_DIRECTION 	rdir[3] = {EAST,NORTH,UPPER};
	HYPER_SURF 	*hs;
	double 		crx_coords[MAXD];
	STATE 		*state;
	int		comp, icoords[3];
	INTERFACE	*grid_intfc = front->grid_intfc;

	switch (dir)
	{
	case 0:
	    i2 = i;
	    if (rect_boundary_type(intfc,dir,0) == NEUMANN_BOUNDARY)
	    {
		i1 = imin[0];
		index0 = d_index2d(i1,i2,top_gmax);
		for (k = 1; k <= 3; ++k)
		{
		    i1 = imin[0] + k - 1;
		    index = d_index2d(i1,i2,top_gmax);
		    vst->dens[3-k] = m_vst->dens[index];
		    vst->engy[3-k] = m_vst->engy[index];
		    vst->pres[3-k] = m_vst->pres[index];
		    vst->momn[0][3-k] = -m_vst->momn[1][index];
		    vst->momn[1][3-k] = m_vst->momn[0][index];
		    vst->momn[2][3-k] = m_vst->momn[2][index];
		}
	    }
	    else if (rect_boundary_type(intfc,dir,0) == SUBDOMAIN_BOUNDARY)
	    {
		for (k = 1; k <= 3; ++k)
		{
		    i1 = imin[0] - k;
		    index = d_index2d(i1,i2,top_gmax);
		    vst->dens[3-k] = m_vst->dens[index];
		    vst->engy[3-k] = m_vst->engy[index];
		    vst->pres[3-k] = m_vst->pres[index];
		    vst->momn[0][3-k] = m_vst->momn[0][index];
		    vst->momn[1][3-k] = m_vst->momn[1][index];
		    vst->momn[2][3-k] = m_vst->momn[2][index];
		}
	    }
	    else
	    {
		i1 = imin[0];
		index = d_index2d(i1,i2,top_gmax);
		comp = top_comp[index];
		icoords[0] = i1;
		icoords[1] = i2;
		if (!FT_StateStructAtGridCrossing(front,grid_intfc,icoords,
			ldir[dir],comp,(POINTER*)&state,&hs,crx_coords))
		{
		    printf("In appendStencilBuffer2d()\n");
		    printf("ERROR: No crossing found!\n");
		    print_int_vector("icoords=", icoords, 2, "\n");
		    printf("direction: %s side %d\n",
		           grid_direction_name(ldir[dir]), 0);
		    clean_up(ERROR);
		}
		switch (wave_type(hs))
		{
		case DIRICHLET_BOUNDARY:
		    setDirichletStates(state,vst,m_vst,hs,icoords,dir,0,0,1);
		    break;
		default: 
		    printf("ERROR appendStencilBuffer2d: "
		    	   "unknown boundary type %d\n", wave_type(hs));
		    clean_up(ERROR);
		}
	    }

	    if (rect_boundary_type(intfc,dir,1) == NEUMANN_BOUNDARY)
	    {
		i1 = imax[0];
		index0 = d_index2d(i1,i2,top_gmax);
		for (k = 1; k <= 3; ++k)
		{
		    i1 = imax[0] - k + 1;
		    offset = imax[0] - imin[0] + 3;
		    index = d_index2d(i1,i2,top_gmax);
		    vst->dens[offset+k] = m_vst->dens[index];
		    vst->engy[offset+k] = m_vst->engy[index];
		    vst->pres[offset+k] = m_vst->pres[index];
		    vst->momn[0][offset+k] = -m_vst->momn[1][index];
		    vst->momn[1][offset+k] = m_vst->momn[0][index];
		    vst->momn[2][offset+k] = m_vst->momn[2][index];
		}
	    }
	    else if (rect_boundary_type(intfc,dir,1) == SUBDOMAIN_BOUNDARY)
	    {
		for (k = 1; k <= 3; ++k)
		{
		    i1 = imax[0] + k;
		    index = d_index2d(i1,i2,top_gmax);
		    offset = imax[0] - imin[0] + 3;
		    vst->dens[offset+k] = m_vst->dens[index];
		    vst->engy[offset+k] = m_vst->engy[index];
		    vst->pres[offset+k] = m_vst->pres[index];
		    vst->momn[0][offset+k] = m_vst->momn[0][index];
		    vst->momn[1][offset+k] = m_vst->momn[1][index];
		    vst->momn[2][offset+k] = m_vst->momn[2][index];
		}
	    }
	    else
	    {
		i1 = imax[0];
		index = d_index2d(i1,i2,top_gmax);
		comp = top_comp[index];
		icoords[0] = i1;
		icoords[1] = i2;
		if (!FT_StateStructAtGridCrossing(front,grid_intfc,icoords,
			rdir[dir],comp,(POINTER*)&state,&hs,crx_coords))
		{
		    printf("In appendStencilBuffer2d()\n");
		    printf("ERROR: No crossing found!\n");
		    print_int_vector("icoords=", icoords, 2, "\n");
		    printf("direction: %s side %d\n",
		           grid_direction_name(ldir[0]), 1);
		    clean_up(ERROR);
		}
		switch (wave_type(hs))
		{
		case DIRICHLET_BOUNDARY:
		    offset = imax[dir] - imin[dir] + nrad;
		    setDirichletStates(state,vst,m_vst,hs,icoords,dir,1,
					offset,1);
		    break;
		default: 
		    printf("ERROR appendStencilBuffer2d"
		    	   "unknown boundary type %d\n", wave_type(hs));
		    clean_up(ERROR);
		}
	    }
	    break;
	case 1:
	    i1 = i;
	    if (rect_boundary_type(intfc,dir,0) == NEUMANN_BOUNDARY)
	    {
		i2 = imin[1];
		index0 = d_index2d(i1,i2,top_gmax);
		for (k = 1; k <= 3; ++k)
		{
		    i2 = imin[1] + k - 1;
		    index = d_index2d(i1,i2,top_gmax);
		    vst->dens[3-k] = m_vst->dens[index];
		    vst->engy[3-k] = m_vst->engy[index];
		    vst->pres[3-k] = m_vst->pres[index];
		    vst->momn[0][3-k] = -m_vst->momn[1][index];
		    vst->momn[1][3-k] = m_vst->momn[0][index];
		    vst->momn[2][3-k] = m_vst->momn[2][index];
		}
	    }
	    else if (rect_boundary_type(intfc,dir,0) == SUBDOMAIN_BOUNDARY)
	    {
		i2 = imin[1];
		index0 = d_index2d(i1,i2,top_gmax);
		for (k = 1; k <= 3; ++k)
		{
		    i2 = imin[1] - k;
		    index = d_index2d(i1,i2,top_gmax);
		    vst->dens[3-k] = m_vst->dens[index];
		    vst->engy[3-k] = m_vst->engy[index];
		    vst->pres[3-k] = m_vst->pres[index];
		    vst->momn[0][3-k] = m_vst->momn[1][index];
		    vst->momn[1][3-k] = m_vst->momn[0][index];
		    vst->momn[2][3-k] = 0.0;
		}
	    }
	    else
	    {
		i2 = imin[1];
		index = d_index2d(i1,i2,top_gmax);
		comp = top_comp[index];
		icoords[0] = i1;
		icoords[1] = i2;
		if (!FT_StateStructAtGridCrossing(front,grid_intfc,icoords,
			ldir[dir],comp,(POINTER*)&state,&hs,crx_coords))
		{
		    printf("In appendStencilBuffer2d()\n");
		    printf("ERROR: No crossing found!\n");
		    print_int_vector("icoords=", icoords, 2, "\n");
		    printf("direction: %s side %d\n",
		           grid_direction_name(ldir[dir]), 0);
		    clean_up(ERROR);
		}
		switch (wave_type(hs))
		{
		case DIRICHLET_BOUNDARY:
		    setDirichletStates(state,vst,m_vst,hs,icoords,dir,0,0,1);
		    break;
		default: 
		    printf("ERROR appendStencilBuffer2d"
		    	   "unknown boundary type %d\n", wave_type(hs));
		    clean_up(ERROR);
		}
	    }

	    if (rect_boundary_type(intfc,dir,1) == NEUMANN_BOUNDARY)
	    {
		i2 = imax[1];
		index0 = d_index2d(i1,i2,top_gmax);
		for (k = 1; k <= 3; ++k)
		{
		    i2 = imax[1] - k + 1;
		    offset = imax[1] - imin[1] + 3;
		    index = d_index2d(i1,i2,top_gmax);
		    vst->dens[offset+k] = m_vst->dens[index];
		    vst->engy[offset+k] = m_vst->engy[index];
		    vst->pres[offset+k] = m_vst->pres[index];
		    vst->momn[0][offset+k] = -m_vst->momn[1][index];
		    vst->momn[1][offset+k] = m_vst->momn[0][index];
		    vst->momn[2][offset+k] = m_vst->momn[2][index];
		}
	    }
	    else if (rect_boundary_type(intfc,dir,1) == SUBDOMAIN_BOUNDARY)
	    {
		i2 = imax[1];
		index0 = d_index2d(i1,i2,top_gmax);
		for (k = 1; k <= 3; ++k)
		{
		    i2 = imax[1] + k;
		    offset = imax[1] - imin[1] + 3;
		    index = d_index2d(i1,i2,top_gmax);
		    vst->dens[offset+k] = m_vst->dens[index];
		    vst->engy[offset+k] = m_vst->engy[index];
		    vst->pres[offset+k] = m_vst->pres[index];
		    vst->momn[0][offset+k] = m_vst->momn[1][index];
		    vst->momn[1][offset+k] = m_vst->momn[0][index];
		    vst->momn[2][offset+k] = 0.0;
		}
	    }
	    else
	    {
		i2 = imax[1];
		index = d_index2d(i1,i2,top_gmax);
		comp = top_comp[index];
		icoords[0] = i1;
		icoords[1] = i2;
		if (!FT_StateStructAtGridCrossing(front,grid_intfc,icoords,
			rdir[dir],comp,(POINTER*)&state,&hs,crx_coords))
		{
		    printf("In appendStencilBuffer2d()\n");
		    printf("ERROR: No crossing found!\n");
		    print_int_vector("icoords=", icoords, 2, "\n");
		    printf("direction: %s side %d\n",
		           grid_direction_name(ldir[dir]), 0);
		    clean_up(ERROR);
		}
		switch (wave_type(hs))
		{
		case DIRICHLET_BOUNDARY:
		    offset = imax[dir] - imin[dir] + nrad;
		    setDirichletStates(state,vst,m_vst,hs,icoords,dir,1,
					offset,1);
		    break;
		default: 
		    printf("ERROR appendStencilBuffer2d"
		    	   "unknown boundary type %d\n", wave_type(hs));
		    clean_up(ERROR);
		}
	    }
	}

}	/* end appendStencilBuffer2d */

void G_CARTESIAN::appendStencilBuffer3d(
	SWEEP *vst,
	SWEEP *m_vst,
	int i1,
	int i2,
	int dir)
{
	int i,j,k,l,offset,index;
	INTERFACE *intfc = front->interf;

	switch (dir)
	{
	case 0:
	    j = i1;	k = i2;
	    if (rect_boundary_type(intfc,dir,0) == SUBDOMAIN_BOUNDARY)
	    {
		for (l = 1; l <= 3; ++l)
		{
		    i = imin[0] - l;
		    index = d_index3d(i,j,k,top_gmax);
		    vst->dens[3-l] = m_vst->dens[index];
		    vst->engy[3-l] = m_vst->engy[index];
		    vst->pres[3-l] = m_vst->pres[index];
		    vst->momn[0][3-l] = m_vst->momn[0][index];
		    vst->momn[1][3-l] = m_vst->momn[1][index];
		    vst->momn[2][3-l] = m_vst->momn[2][index];
		}
	    }
	    if (rect_boundary_type(intfc,dir,1) == SUBDOMAIN_BOUNDARY)
	    {
		for (l = 1; l <= 3; ++l)
		{
		    i = imax[0] + l;
		    index = d_index3d(i,j,k,top_gmax);
		    offset = imax[0] - imin[0] + 3;
		    vst->dens[offset+l] = m_vst->dens[index];
		    vst->engy[offset+l] = m_vst->engy[index];
		    vst->pres[offset+l] = m_vst->pres[index];
		    vst->momn[0][offset+l] = m_vst->momn[0][index];
		    vst->momn[1][offset+l] = m_vst->momn[1][index];
		    vst->momn[2][offset+l] = m_vst->momn[2][index];
		}
	    }
	    break;
	case 1:
	    k = i1;	i = i2;
	    if (rect_boundary_type(intfc,dir,0) == SUBDOMAIN_BOUNDARY)
	    {
		for (l = 1; l <= 3; ++l)
		{
		    j = imin[1] - l;
		    index = d_index3d(i,j,k,top_gmax);
		    vst->dens[3-l] = m_vst->dens[index];
		    vst->engy[3-l] = m_vst->engy[index];
		    vst->pres[3-l] = m_vst->pres[index];
		    vst->momn[0][3-l] = m_vst->momn[1][index];
		    vst->momn[1][3-l] = m_vst->momn[2][index];
		    vst->momn[2][3-l] = m_vst->momn[0][index];
		}
	    }
	    else if (rect_boundary_type(intfc,dir,0) == NEUMANN_BOUNDARY)
	    {
		for (l = 1; l <= 3; ++l)
		{
		    j = imin[1] + l - 1;
		    index = d_index3d(i,j,k,top_gmax);
		    vst->dens[3-l] = m_vst->dens[index];
		    vst->engy[3-l] = m_vst->engy[index];
		    vst->pres[3-l] = m_vst->pres[index];
		    vst->momn[0][3-l] = -m_vst->momn[1][index];
		    vst->momn[1][3-l] = m_vst->momn[2][index];
		    vst->momn[2][3-l] = m_vst->momn[0][index];
		}
	    }
	    if (rect_boundary_type(intfc,dir,1) == SUBDOMAIN_BOUNDARY)
	    {
		for (l = 1; l <= 3; ++l)
		{
		    j = imax[1] + l;
		    index = d_index3d(i,j,k,top_gmax);
		    offset = imax[1] - imin[1] + 3;
		    vst->dens[offset+l] = m_vst->dens[index];
		    vst->engy[offset+l] = m_vst->engy[index];
		    vst->pres[offset+l] = m_vst->pres[index];
		    vst->momn[0][offset+l] = m_vst->momn[1][index];
		    vst->momn[1][offset+l] = m_vst->momn[2][index];
		    vst->momn[2][offset+l] = m_vst->momn[0][index];
		}
	    }
	    else if (rect_boundary_type(intfc,dir,1) == NEUMANN_BOUNDARY)
	    {
		for (l = 1; l <= 3; ++l)
		{
		    j = imax[1] - l + 1;
		    offset = imax[1] - imin[1] + 3;
		    index = d_index3d(i,j,k,top_gmax);
		    vst->dens[offset+l] = m_vst->dens[index];
		    vst->engy[offset+l] = m_vst->engy[index];
		    vst->pres[offset+l] = m_vst->pres[index];
		    vst->momn[0][offset+l] = -m_vst->momn[1][index];
		    vst->momn[1][offset+l] = m_vst->momn[2][index];
		    vst->momn[2][offset+l] = m_vst->momn[0][index];
		}
	    }
	    break;
	case 2:
	    i = i1;	j = i2;
	    if (rect_boundary_type(intfc,dir,0) == NEUMANN_BOUNDARY)
	    {
		for (l = 1; l <= 3; ++l)
		{
		    k = imin[2] + l - 1;
		    index = d_index3d(i,j,k,top_gmax);
		    vst->dens[3-l] = m_vst->dens[index];
		    vst->engy[3-l] = m_vst->engy[index];
		    vst->pres[3-l] = m_vst->pres[index];
		    vst->momn[0][3-l] = -m_vst->momn[2][index];
		    vst->momn[1][3-l] = m_vst->momn[0][index];
		    vst->momn[2][3-l] = m_vst->momn[1][index];
		}
	    }
	    if (rect_boundary_type(intfc,dir,1) == NEUMANN_BOUNDARY)
	    {
		for (l = 1; l <= 3; ++l)
		{
		    k = imax[2] - l + 1;
		    offset = imax[2] - imin[2] + 3;
		    index = d_index3d(i,j,k,top_gmax);
		    vst->dens[offset+l] = m_vst->dens[index];
		    vst->engy[offset+l] = m_vst->engy[index];
		    vst->pres[offset+l] = m_vst->pres[index];
		    vst->momn[0][offset+l] = -m_vst->momn[2][index];
		    vst->momn[1][offset+l] = m_vst->momn[0][index];
		    vst->momn[2][offset+l] = m_vst->momn[1][index];
		}
	    }
	    
	    if (rect_boundary_type(intfc,dir,0) == SUBDOMAIN_BOUNDARY)
	    {
		for (l = 1; l <= 3; ++l)
		{
		    k = imin[2] - l; 
		    index = d_index3d(i,j,k,top_gmax);
		    vst->dens[3-l] = m_vst->dens[index];
		    vst->engy[3-l] = m_vst->engy[index];
		    vst->pres[3-l] = m_vst->pres[index];
		    vst->momn[0][3-l] = m_vst->momn[2][index];
		    vst->momn[1][3-l] = m_vst->momn[0][index];
		    vst->momn[2][3-l] = m_vst->momn[1][index];
		}
	    }
	    if (rect_boundary_type(intfc,dir,1) == SUBDOMAIN_BOUNDARY)
	    {
		for (l = 1; l <= 3; ++l)
		{
		    k = imax[2] + l;
		    offset = imax[2] - imin[2] + 3;
		    index = d_index3d(i,j,k,top_gmax);
		    vst->dens[offset+l] = m_vst->dens[index];
		    vst->engy[offset+l] = m_vst->engy[index];
		    vst->pres[offset+l] = m_vst->pres[index];
		    vst->momn[0][offset+l] = m_vst->momn[2][index];
		    vst->momn[1][offset+l] = m_vst->momn[0][index];
		    vst->momn[2][offset+l] = m_vst->momn[1][index];
		}
	    }
	}
}	/* end appendStencilBuffer3d */

void G_CARTESIAN::scatMeshStates()
{
	SWEEP vst;
	allocMeshVst(&vst);
	copyToMeshVst(&vst);
	scatMeshVst(&vst);
	copyFromMeshVst(vst);
	freeVst(&vst);
}	/* end scatMeshStates */

void G_CARTESIAN::freeVst(
	SWEEP *vst)
{
	FT_FreeThese(4,vst->dens,vst->engy,vst->pres,vst->momn);
}	/* end freeVstFlux */

void G_CARTESIAN::freeFlux(
	FSWEEP *flux)
{
	FT_FreeThese(3,flux->dens_flux,flux->engy_flux,flux->momn_flux);
}

void G_CARTESIAN::addMeshFluxToVst(
	SWEEP *m_vst,
	FSWEEP m_flux,
	double chi)
{
	int 		i,j,k,l,index;
	double		ke,c,u;
	EOS_PARAMS	*eos;
	STATE		st;
	int		comp;
	double		temp;

	switch (dim)
	{
	case 1:
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index1d(i,top_gmax);
		comp = top_comp[index];
		if (!gas_comp(comp))
		{
		    m_vst->dens[index] = 0.0;
		    m_vst->engy[index] = 0.0;
		    for (l = 0; l < dim; ++l)
		    	m_vst->momn[l][index] = 0.0; 
		    continue;
		}
		eos = &(eqn_params->eos[comp]);

		m_vst->dens[index] += chi*m_flux.dens_flux[index];
		m_vst->engy[index] += chi*m_flux.engy_flux[index];
		ke = u = 0.0;
		for (l = 0; l < dim; ++l)
		{
		    m_vst->momn[l][index] += 
			chi*m_flux.momn_flux[l][index];
		    ke += sqr(m_vst->momn[l][index]);
		    u += sqr(m_vst->momn[l][index]);
		}
		
		CovertVstToState(&st, m_vst, eos, index, dim);
		checkCorrectForTolerance(&st);
		m_vst->dens[index] = st.dens;
		m_vst->pres[index] = st.pres;
		m_vst->engy[index] = st.engy;
		u = sqrt(u)/m_vst->dens[index];
		c = EosSoundSpeed(&st);
		temp = std::max((std::max(u,fabs(u-c))),(fabs(u+c)));
                if (max_speed < temp)
                    max_speed = temp;
	    }
	    scatMeshVst(m_vst);
	    break;
	case 2:
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index2d(i,j,top_gmax);
		comp = top_comp[index];
		if (!gas_comp(comp))
		{
		    m_vst->dens[index] = 0.0;
		    m_vst->engy[index] = 0.0;
		    for (l = 0; l < dim; ++l)
		    	m_vst->momn[l][index] = 0.0; 
		    continue;
		}
		eos = &(eqn_params->eos[comp]);

		m_vst->dens[index] += chi*m_flux.dens_flux[index];
		m_vst->engy[index] += chi*m_flux.engy_flux[index];
		ke = u = 0.0;
		for (l = 0; l < dim; ++l)
		{
		    m_vst->momn[l][index] += 
			chi*m_flux.momn_flux[l][index];
		    ke += sqr(m_vst->momn[l][index]);
		    u += sqr(m_vst->momn[l][index]);
		}
		
		CovertVstToState(&st, m_vst, eos, index, dim);
		checkCorrectForTolerance(&st);
		m_vst->dens[index] = st.dens;
		m_vst->pres[index] = st.pres;
		m_vst->engy[index] = st.engy;
		u = sqrt(u)/m_vst->dens[index];
		c = EosSoundSpeed(&st);
		temp = std::max((std::max(u,fabs(u-c))),(fabs(u+c)));
                if (max_speed < temp)
                    max_speed = temp;
	    }
	    scatMeshVst(m_vst);
	    break;
	case 3:
	    for (k = imin[2]; k <= imax[2]; ++k)
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index3d(i,j,k,top_gmax);
		comp = top_comp[index];
		eos = &(eqn_params->eos[comp]);

		m_vst->dens[index] += chi*m_flux.dens_flux[index];
		m_vst->engy[index] += chi*m_flux.engy_flux[index];
		ke = u = 0.0;
		for (l = 0; l < dim; ++l)
		{
		    m_vst->momn[l][index] += 
				chi*m_flux.momn_flux[l][index];
		    ke += sqr(m_vst->momn[l][index]);
		    u += sqr(m_vst->momn[l][index]);
		}
		
		CovertVstToState(&st, m_vst, eos, index, dim);
		checkCorrectForTolerance(&st);
		m_vst->dens[index] = st.dens;
		m_vst->pres[index] = st.pres;
		m_vst->engy[index] = st.engy;
		u = sqrt(u)/m_vst->dens[index];
		c = EosSoundSpeed(&st);
		temp = std::max((std::max(u,fabs(u-c))),(fabs(u+c)));
                if (max_speed < temp)
                    max_speed = temp;
	    }
	    scatMeshVst(m_vst);
	}
}	/* end addMeshFluxToVst */


void G_CARTESIAN::appendGhostBuffer(
	SWEEP *vst,
	SWEEP *m_vst,
	int n,
	int *icoords,
	int idir,
	int nb)
{
	int		i,j,k,index,ic[MAXD];
	GRID_DIRECTION 	ldir[3] = {WEST,SOUTH,LOWER};
	GRID_DIRECTION 	rdir[3] = {EAST,NORTH,UPPER};
	HYPER_SURF 	*hs;
	COMPONENT 	comp;
	double 		crx_coords[MAXD];
	STATE 		*state,ghost_st;
	int		ind2[2][2] = {{0,1},{1,0}};
	int		ind3[3][3] = {{0,1,2},{1,2,0},{2,0,1}};
	int 		ic_next[MAXD];
	INTERFACE	*grid_intfc = front->grid_intfc;

	if (debugging("append_buffer"))
		printf("Entering appendGhostBuffer()\n");
	for (i = 0; i < dim; ++i) ic[i] = icoords[i];
	
	index = d_index(ic,top_gmax,dim);
	comp = cell_center[index].comp;
	
	switch(nb)
	{
	case 0:
	    for (i = 1; i <= nrad; ++i)
	    {
		ic[idir] = icoords[idir] - i;
		index = d_index(ic,top_gmax,dim);
		    
		if (!needBufferFromIntfc(comp,cell_center[index].comp))
		{
		    vst->dens[nrad-i] = m_vst->dens[index];
		    vst->engy[nrad-i] = m_vst->engy[index];
		    vst->pres[nrad-i] = m_vst->pres[index];

		    for (j = 0; j < 3; j++)
			vst->momn[j][nrad-i] = 0.0;
		    if (dim == 1)
			vst->momn[0][nrad-i] = m_vst->momn[0][index];
		    else if (dim == 2)
			for(j = 0; j < 2; j++)
			    vst->momn[j][nrad-i] = 
			    	 	m_vst->momn[ind2[idir][j]][index];
		    else if (dim == 3)
			for (j = 0; j < 3; j++)
			    vst->momn[j][nrad-i] = 
			    	 	m_vst->momn[ind3[idir][j]][index];
		}
		else
		{
		    for (k = 0; k < dim; ++k)
			ic_next[k] = ic[k];
		    ic_next[idir]++;
		    if (!FT_StateStructAtGridCrossing(front,grid_intfc,ic_next,
			ldir[idir],comp,(POINTER*)&state,&hs,crx_coords))
		    {
		    	(void) printf("In appendGhostBuffer()\n");
		    	(void) printf("ERROR: No crossing found!\n");
		    	(void) print_int_vector("icoords=", ic_next,3,"\n");
		    	(void) printf("direction: %s side %d\n",
		           		grid_direction_name(rdir[idir]), nb);
		    	clean_up(ERROR);
		    }
		    switch (wave_type(hs))
		    {
		    case NEUMANN_BOUNDARY:
		    case MOVABLE_BODY_BOUNDARY:
		    	setNeumannStates(vst,m_vst,hs,state,ic_next,idir,
						nb,0,i,comp);
		    	break;
		    case DIRICHLET_BOUNDARY:
		    	setDirichletStates(state,vst,m_vst,hs,ic_next,
					idir,nb,0,i);
		    	break;
		    case FIRST_PHYSICS_WAVE_TYPE:
		    	//GFM
		    	GFMGhostState(ic,comp,&ghost_st);
		    	for (k = i; k <= nrad; ++k)
		    	{
		    	    vst->dens[nrad-k] = ghost_st.dens;
		    	    vst->engy[nrad-k] = ghost_st.engy;
		    	    vst->pres[nrad-k] = ghost_st.pres;
			
			    for (j=0; j < 3; j++)
			    	    vst->momn[j][nrad-k] = 0.0;
			    if (dim == 1)
				vst->momn[0][nrad-k] = ghost_st.momn[0];
			    else if (dim == 2)
			    	for (j=0; j < 2; j++)
				    vst->momn[j][nrad-k] = 
				     	    ghost_st.momn[ind2[idir][j]];
			    else if (dim == 3)
			    	for (j = 0; j < 3; j++)
				    vst->momn[j][nrad-k] = 
				     	    ghost_st.momn[ind3[idir][j]];
		    	}
		    	break;
		    default:
		    	(void) printf("In appendGhostBuffer(): ");
		    	(void) print_wave_type("Unknown wave type ",
					wave_type(hs),"\n",front->interf);
		    	(void) print_int_vector("icoords=", icoords,3,"\n");
		    	clean_up(ERROR);
		    }
		    break;
		}
	    }
	    break;
	case 1:
	    for (i = 1; i <= nrad; ++i)
	    {
		ic[idir] = icoords[idir] + i;
		index = d_index(ic,top_gmax,dim);
		if (!needBufferFromIntfc(comp,cell_center[index].comp))
		{
		    vst->dens[n+nrad+i-1] = m_vst->dens[index];
		    vst->engy[n+nrad+i-1] = m_vst->engy[index];
		    vst->pres[n+nrad+i-1] = m_vst->pres[index];
		    
		    for (j = 0; j < 3; j++)
			vst->momn[j][n+nrad+i-1] = 0.0;
		    if (dim == 1)
			vst->momn[0][n+nrad+i-1] = 
			         	m_vst->momn[0][index];
		    else if (dim == 2)
			for(j = 0; j < 2; j++)
			    	vst->momn[j][n+nrad+i-1] = 
			         	m_vst->momn[ind2[idir][j]][index];
		    else if (dim == 3)
			for (j = 0; j < 3; j++)
			    vst->momn[j][n+nrad+i-1] = 
			         	m_vst->momn[ind3[idir][j]][index];
		}
		else
		{
		    for (k = 0; k < dim; ++k)
			ic_next[k] = ic[k];
		    ic_next[idir]--;
		    if (!FT_StateStructAtGridCrossing(front,grid_intfc,ic_next,
			rdir[idir],comp,(POINTER*)&state,&hs,crx_coords))
		    {
		    	(void) printf("In appendGhostBuffer()\n");
		    	(void) printf("ERROR: No crossing found!\n");
		    	(void) print_int_vector("icoords=",ic_next,3,"\n");
		    	(void) printf("direction: %s side %d\n",
		            	grid_direction_name(rdir[idir]), nb);
		    	clean_up(ERROR);
		    }
		    switch (wave_type(hs))
		    {
		    case NEUMANN_BOUNDARY:
		    case MOVABLE_BODY_BOUNDARY:
		    	setNeumannStates(vst,m_vst,hs,state,ic_next,idir,
						nb,n,i,comp);
		    	break;
		    case DIRICHLET_BOUNDARY:
		    	setDirichletStates(state,vst,m_vst,hs,ic_next,idir,nb,
						n,i);
		    	break;
		    case FIRST_PHYSICS_WAVE_TYPE:
		    	//GFM
		    	GFMGhostState(ic,comp,&ghost_st);

		    	for (k = i; k <= nrad; ++k)
		    	{
		    	    vst->dens[n+nrad+k-1] = ghost_st.dens;
		    	    vst->engy[n+nrad+k-1] = ghost_st.engy;
		    	    vst->pres[n+nrad+k-1] = ghost_st.pres;
			
			    for(j=0; j<3; j++)
			    	vst->momn[j][n+nrad+k-1] = 0.0;
			    if (dim == 1)
				vst->momn[0][n+nrad+k-1] = ghost_st.momn[0];
			    else if (dim == 2)
			    	for(j = 0; j < 2; j++)
				    vst->momn[j][n+nrad+k-1] = 
				     	    ghost_st.momn[ind2[idir][j]];
			    else if (dim == 3)
			    	for(j = 0; j < 3; j++)
				    vst->momn[j][n+nrad+k-1] = 
				     	    ghost_st.momn[ind3[idir][j]];
		    	}
		    	break;
		    default:
		    	(void) printf("In appendGhostBuffer(): ");
		    	(void) print_wave_type("Unknown wave type ",
				wave_type(hs),"\n",front->interf);
		    	(void) print_int_vector("icoords=",icoords,3,"\n");
		    	(void) printf("nb = %d\n",nb);
		    	clean_up(ERROR);
		    }
		    break;
		}
	    }
	}
}	/* end appendGhostBuffer */

//ghost fluid method.

void G_CARTESIAN::solve_exp_value()
{
	int		i, j, k, n;
	int		index;
	double		**gnor = eqn_params->gnor;

	fflush(NULL);

	get_normal_from_front();

	if (dim == 1)
	{
	    for(k=0; k<dim; k++)
	    {
		for (i = imin[0]; i <= imax[0]; ++i)
		{
	 	    index = d_index1d(i,top_gmax);
		    array[index] = gnor[k][index];
		}
		scatMeshArray();
        	for (i = 0; i <= top_gmax[0]; i++)
        	{
	    	    index  = d_index1d(i,top_gmax);
	    	    gnor[k][index] = array[index];
		}
	    }
	}
	else if (dim == 2)
	{
	    for(k=0; k<dim; k++)
	    {
		for (j = imin[1]; j <= imax[1]; ++j)
		for (i = imin[0]; i <= imax[0]; ++i)
		{
	 	    index = d_index2d(i,j,top_gmax);
		    array[index] = gnor[k][index];
		}
		scatMeshArray();
        	for (j = 0; j <= top_gmax[1]; j++)
        	for (i = 0; i <= top_gmax[0]; i++)
        	{
	    	    index  = d_index2d(i,j,top_gmax);
	    	    gnor[k][index] = array[index];
		}
	    }
	}
	else
	{
	    for(k=0; k<dim; k++)
	    {
		for (n = imin[2]; n <= imax[2]; ++n)
		for (j = imin[1]; j <= imax[1]; ++j)
		for (i = imin[0]; i <= imax[0]; ++i)
		{
	    	    index = d_index3d(i,j,n,top_gmax);
	    	    array[index] = gnor[k][index];
		}
		scatMeshArray();
        	for (n = 0; n <= top_gmax[2]; n++)
        	for (j = 0; j <= top_gmax[1]; j++)
        	for (i = 0; i <= top_gmax[0]; i++)
        	{
	    	    index  = d_index3d(i,j,n,top_gmax);
	    	    gnor[k][index] = array[index];
		}
	    }
	}
}

void G_CARTESIAN::scatMeshGhost()
{
	int		i, j, k, n, index;
	double		***Gvel = eqn_params->Gvel;
	double		**Gdens = eqn_params->Gdens;
	double		**Gpres = eqn_params->Gpres;

	if(dim == 2)
	{
	for(k=0; k<2; k++)
	{
	for (j = imin[1]; j <= imax[1]; ++j)
	for (i = imin[0]; i <= imax[0]; ++i)
	{
	    index = d_index2d(i,j,top_gmax);
	    array[index] = Gdens[k][index];
	}
	scatMeshArray();
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
	    index  = d_index2d(i,j,top_gmax);
	    Gdens[k][index] = array[index];
	}
	
	for (j = imin[1]; j <= imax[1]; ++j)
	for (i = imin[0]; i <= imax[0]; ++i)
	{
	    index = d_index2d(i,j,top_gmax);
	    array[index] = Gpres[k][index];
	}
	scatMeshArray();
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
	    index  = d_index2d(i,j,top_gmax);
	    Gpres[k][index] = array[index];
	}

	for (j = imin[1]; j <= imax[1]; ++j)
	for (i = imin[0]; i <= imax[0]; ++i)
	{
	    index = d_index2d(i,j,top_gmax);
	    array[index] = Gvel[k][0][index];
	}
	scatMeshArray();
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
	    index  = d_index2d(i,j,top_gmax);
	    Gvel[k][0][index] = array[index];
	}

	for (j = imin[1]; j <= imax[1]; ++j)
	for (i = imin[0]; i <= imax[0]; ++i)
	{
	    index = d_index2d(i,j,top_gmax);
	    array[index] = Gvel[k][1][index];
	}
	scatMeshArray();
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
	    index  = d_index2d(i,j,top_gmax);
	    Gvel[k][1][index] = array[index];
	}
	}    //for k
	}    //if dim == 2
	else if(dim == 3)
	{
	for(k=0; k<2; k++)
	{
	for (n = imin[2]; n <= imax[2]; ++n)
	for (j = imin[1]; j <= imax[1]; ++j)
	for (i = imin[0]; i <= imax[0]; ++i)
	{
	    index = d_index3d(i,j,n,top_gmax);
	    array[index] = Gdens[k][index];
	}
	scatMeshArray();
        for (n = 0; n <= top_gmax[2]; n++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
	    index  = d_index3d(i,j,n,top_gmax);
	    Gdens[k][index] = array[index];
	}

	for (n = imin[2]; n <= imax[2]; ++n)
	for (j = imin[1]; j <= imax[1]; ++j)
	for (i = imin[0]; i <= imax[0]; ++i)
	{
	    index = d_index3d(i,j,n,top_gmax);
	    array[index] = Gpres[k][index];
	}
	scatMeshArray();
        for (n = 0; n <= top_gmax[2]; n++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
	    index  = d_index3d(i,j,n,top_gmax);
	    Gpres[k][index] = array[index];
	}
	
	for (n = imin[2]; n <= imax[2]; ++n)
	for (j = imin[1]; j <= imax[1]; ++j)
	for (i = imin[0]; i <= imax[0]; ++i)
	{
	    index = d_index3d(i,j,n,top_gmax);
	    array[index] = Gvel[k][0][index];
	}
	scatMeshArray();
        for (n = 0; n <= top_gmax[2]; n++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
	    index  = d_index3d(i,j,n,top_gmax);
	    Gvel[k][0][index] = array[index];
	}

	for (n = imin[2]; n <= imax[2]; ++n)
	for (j = imin[1]; j <= imax[1]; ++j)
	for (i = imin[0]; i <= imax[0]; ++i)
	{
	    index = d_index3d(i,j,n,top_gmax);
	    array[index] = Gvel[k][1][index];
	}
	scatMeshArray();
        for (n = 0; n <= top_gmax[2]; n++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
	    index  = d_index3d(i,j,n,top_gmax);
	    Gvel[k][1][index] = array[index];
	}
	
	for (n = imin[2]; n <= imax[2]; ++n)
	for (j = imin[1]; j <= imax[1]; ++j)
	for (i = imin[0]; i <= imax[0]; ++i)
	{
	    index = d_index3d(i,j,n,top_gmax);
	    array[index] = Gvel[k][2][index];
	}
	scatMeshArray();
        for (n = 0; n <= top_gmax[2]; n++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
	    index  = d_index3d(i,j,n,top_gmax);
	    Gvel[k][2][index] = array[index];
	}
	}    //for k
	}    //for dim == 3
}

#define	corner_index(p,i,gr)	irint(floor(((p)-(gr)->L[i])/(gr)->h[i]-0.5))

boolean	find_block(
	double		*f,
	int		*icrds,
	double		*p,
	RECT_GRID	*gr)
{
	int	i;
	int	dim=gr->dim;

	for(i=0; i<dim; i++)
	{
	    icrds[i] = corner_index(p[i],i,gr);
	    if(icrds[i] < -gr->lbuf[i] || icrds[i] >= gr->gmax[i]+gr->ubuf[i]-1)
		return  NO;
	    f[i] = p[i] - (gr->L[i]+(0.5+icrds[i])*gr->h[i]);
	    f[i] /= gr->h[i];
	}
	return  YES;
}

boolean G_CARTESIAN::get_ave_normal(
	int		*ic,
	int		***norset)
{
	double		f;
	int		i, j, k, n, ic1[3], ic2[3], dir, num;
	boolean		found;
	int		index0, index;
	double		**gnor = eqn_params->gnor;
	
	found = NO;

	for(i=0; i<dim; i++)
	for(j=0; j<2; j++)
	{
		dir = j == 0 ? -1 : 1;
		ft_assign(ic1, ic, 3*INT);
		ic1[i] = ic[i] + dir;

		if(ic1[i] < 0 || ic1[i] > top_gmax[i])
		    continue;
		if(norset[ic1[0]][ic1[1]][ic1[2]] == 1)
		    found = YES;
	}
	if(!found)
	    return NO;
	
	index0  = d_index(ic,top_gmax,dim);

	gnor[0][index0] = 0.0;
	if (dim > 1)
	    gnor[1][index0] = 0.0;
	if (dim > 2)
	    gnor[2][index0] = 0.0;

	num = 0;
	for(i=ic[0]-1; i<=ic[0]+1; i++)
	for(j=ic[1]-1; j<=ic[1]+1; j++)
	for(k=ic[2]-1; k<=ic[2]+1; k++)
	{
	    if(i < 0 || i > top_gmax[0] || 
	       j < 0 || j > top_gmax[1] || 
	       k < 0 || k > top_gmax[2]) 
		continue;
	    if(norset[i][j][k] != 1)
		continue;

	    ic2[0] = i;
	    ic2[1] = j;
	    ic2[2] = k;
	    index  = d_index(ic2,top_gmax,dim);
		    
	    //do not use length weighted normal direction
	    gnor[0][index0] += gnor[0][index];
	    if(dim > 1)
	    	gnor[1][index0] += gnor[1][index];
	    if(dim > 2)
		gnor[2][index0] += gnor[2][index];
	    num++;
	}
	
	f = 0.0;
	for(n=0; n<dim; n++)
	    f += sqr(gnor[n][index0]);
	f = sqrt(f);

	if(f < 1.0e-6)
	{
	    gnor[0][index0] = 0.0;
	    if (dim > 1) gnor[1][index0] = 0.0;
	    if (dim > 2) gnor[2][index0] = 0.0;
	}
	else
	{
	    gnor[0][index0] /= f;
	    if (dim > 1) gnor[1][index0] /= f;
	    if (dim > 2) gnor[2][index0] /= f;
	}

	return YES;
}

boolean	find_block(double*,int*,double*,RECT_GRID*);

void	get_normal_from_front();

//it will fill gnor field by interface normals
void G_CARTESIAN::get_normal_from_front()
{
	INTERFACE               *intfc = front->interf;
	RECT_GRID		*rgr = front->rect_grid;
	HYPER_SURF              *hs;
	HYPER_SURF_ELEMENT      *hse;
	POINT                   *p;
	int			i,j,k,n, num;
	int			ic[3];
	double			curv,nor[3],d[3],f,d1,d2,d3,*pt,tol;
	int			ix, iy, iz, index;
	boolean			found;
	double			**gnor = eqn_params->gnor;
	static	int		***norset;
	int ict[3];
	double ptt[3];
	boolean status;

	if (norset == NULL)
	    FT_TriArrayMemoryAlloc((POINTER*)&norset,top_gmax[0]+1,
				top_gmax[1]+1,top_gmax[2]+1,INT);

	tol = hmin*1.0e-6;

	for (i = 0; i <= top_gmax[0]; i++)
	for (j = 0; j <= top_gmax[1]; j++)
	for (k = 0; k <= top_gmax[2]; k++)
	{
	    ic[0] = i;
	    ic[1] = j;
	    ic[2] = k;
	    index = d_index(ic,top_gmax,dim);
	    gnor[0][index] = 0.0;
	    if (dim > 1)
	    	gnor[1][index] = 0.0;
	    if (dim > 2)
		gnor[2][index] = 0.0;
	    norset[i][j][k] = 0;
	}

	(void) next_point(intfc,NULL,NULL,NULL);
	while (next_point(intfc,&p,&hse,&hs))
	{
	    if (Boundary_point(p))
	    {
		p->_nor[0] = 0.0;
		p->_nor[1] = 0.0;
		p->_nor[2] = 0.0;
		continue;
	    }

	    normal(p,hse,hs,nor,front);
	    curv = p->curvature;

	    pt = Coords(p);
	   
	    status = rect_in_which(pt,ict,top_grid);
	    if (!status) continue;
	    for(i = 0; i < dim; i++)
		ptt[i] = top_grid->L[i] + ict[i]*top_grid->h[i];

	    for(i = 0; i < dim; i++)
	    {
		d[i] = fabs(pt[i]-ptt[i])/rgr->h[i];
	        if(d[i] < -tol || d[i] > 1.0 + tol)
		{
		    status = NO;
		}
	    }
	    if (status == NO) continue;

	    if (dim == 1)
	    {
	    	for(i = 0; i < 2; i++)
		{
		    ix = ict[0] + i;
		    d1 = (i == 0) ? fabs(1.0-d[0]) : d[0];
		    f = d1;

		    index = d_index1d(ix,top_gmax);
		    gnor[0][index] += nor[0]*f;
		    norset[ix][0][0] = 1;
		}
	    }
	    else if (dim == 2)
	    {
	    	for (i = 0; i < 2; i++)
		for (j = 0; j < 2; j++)
		{
		    ix = ict[0] + i;
		    iy = ict[1] + j;
		    d1 = (i == 0) ? fabs(1.0-d[0]) : d[0];
		    d2 = (j == 0) ? fabs(1.0-d[1]) : d[1];
		    f = d1*d2;
		    index = d_index2d(ix,iy,top_gmax);
		    gnor[0][index] += nor[0]*f;
		    gnor[1][index] += nor[1]*f;
		    norset[ix][iy][0] = 1;
		}
	    }
	    else if (dim == 3)
	    {
	    	for (i = 0; i < 2; i++)
		for (j = 0; j < 2; j++)
		for (k = 0; k < 2; k++)
		{
		    ix = ict[0] + i;
		    iy = ict[1] + j;
		    iz = ict[2] + k;
		    d1 = (i == 0) ? fabs(1.0-d[0]) : d[0];
		    d2 = (j == 0) ? fabs(1.0-d[1]) : d[1];
		    d3 = (k == 0) ? fabs(1.0-d[2]) : d[2];
		    f = d1*d2*d3;
		    index = d_index3d(ix,iy,iz,top_gmax);
		    gnor[0][index] += nor[0]*f;
		    gnor[1][index] += nor[1]*f;
		    gnor[2][index] += nor[2]*f;
		    norset[ix][iy][iz] = 1;
		}
	    }
	}

	for (i = 0; i <= top_gmax[0]; i++)
	for (j = 0; j <= top_gmax[1]; j++)
	for (k = 0; k <= top_gmax[2]; k++)
	{
	    //make sure Vel(st) is assigned 
	    if (norset[i][j][k] != 1)
		continue;

	    ic[0] = i;
	    ic[1] = j;
	    ic[2] = k;
	    index = d_index(ic,top_gmax,dim);
	    f = 0.0;
	    for(n=0; n<dim; n++)
		f += sqr(gnor[n][index]);
	    f = sqrt(f);

	    if (f < 1.0e-10)
	    {
		gnor[0][index] = 0.0;
		if (dim > 1) gnor[1][index] = 0.0;
		if (dim > 2) gnor[2][index] = 0.0;
	    }
	    else
	    {
		gnor[0][index] /= f;
		if (dim > 1) gnor[1][index] /= f;
		if (dim > 2) gnor[2][index] /= f;
	    }
	}

	found = YES;
	num = 1;
	while (found && num > 0)
	{
	    found = NO;
	    num = 0;

	    for (i = 0; i <= top_gmax[0]; i++)
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (k = 0; k <= top_gmax[2]; k++)
	    {
		if(norset[i][j][k] != 0)
		    continue;

		found = YES;
		ic[0] = i;
		ic[1] = j;
		ic[2] = k;

		if(get_ave_normal(ic,norset))
		{
		    num++;
		    norset[i][j][k] = 2;
		}
	    }

	    for (i = 0; i <= top_gmax[0]; i++)
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (k = 0; k <= top_gmax[2]; k++)
		if(norset[i][j][k] == 2)
		    norset[i][j][k] = 1;
	}
}

void G_CARTESIAN::tecplot_interior_states(
			char	*bname)
{
	char		s[1000];
	double		coords[3];
	int		ix, iy, iz, comp, i, imin[3], imax[3];
	double		**vel = eqn_params->vel;
	double		*dens = eqn_params->dens;
	double		*pres = eqn_params->pres;
	double		**gnor = eqn_params->gnor;
	double		***Gvel = eqn_params->Gvel;
	double		**Gdens = eqn_params->Gdens;
	double		**Gpres = eqn_params->Gpres;
	FILE		*fp;
	int		index;

	sprintf(s,"%s-%d.plt", bname,pp_mynode());
	printf("tecplot_interior_states  file name %s \n",s);

	fp = fopen(s, "w");
	if(fp == NULL)
	{
	    printf("WARNING tecplot_interior_states, can not open file %s\n", s);
	    return; 
	}
	
	for(i=0; i<3; i++)
	{
	    imin[i] = 0;
	    imax[i] = top_gmax[i];
	}

	fprintf(fp, "TITLE = \"inner states\" ");
	if(dim == 2)
	{
	    fprintf(fp, "VARIABLES = \"x\", \"y\", \"comp\",  ");
	    fprintf(fp, "\"dens\", \"press\", \"u\", \"v\",  " );
	    fprintf(fp, "\"nx\", \"ny\",  " );
	    fprintf(fp, "\"dens1\", \"press1\", \"u1\", \"v1\",  " );
	    fprintf(fp, "\"dens2\", \"press2\", \"u2\", \"v2\"  \n" );
	}
	else
	{
	    fprintf(fp, "VARIABLES = \"x\", \"y\", \"z\", \"comp\",  ");
	    fprintf(fp, "\"dens\", \"press\", \"u\", \"v\", \"w\", " );
	    fprintf(fp, "\"nx\", \"ny\", \"nz\", " );
	    fprintf(fp, "\"dens1\", \"press1\", \"u1\", \"v1\", \"w1\" " );
	    fprintf(fp, "\"dens2\", \"press2\", \"u2\", \"v2\"  \"w2\"\n" );
	}

	if(dim == 2)
	    fprintf(fp, "ZONE i=%d, j=%d \n", imax[0]-imin[0]+1, imax[1]-imin[1]+1);
	else
	    fprintf(fp, "ZONE i=%d, j=%d, k=%d \n", imax[0]-imin[0]+1, imax[1]-imin[1]+1, imax[2]-imin[2]+1);

	if(dim == 2)
	{
	    for(iy=imin[1]; iy <= imax[1]; iy++)
		  for(ix=imin[0]; ix <= imax[0]; ix++)
		  {
			index = d_index2d(ix,iy,top_gmax);
			
			getRectangleCenter(index, coords);
			comp = cell_center[index].comp;

			fprintf(fp, "%f ", coords[0]);
			fprintf(fp, "%f ", coords[1]);
			fprintf(fp, "%d ", comp);
		
			if(!gas_comp(comp))
			{
			    fprintf(fp, "%12.5e %12.5e   %12.5e %12.5e   %12.5e %12.5e  ", 
					0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
			    fprintf(fp, "%12.5e %12.5e   %12.5e %12.5e  ", 
			    		0.0, 0.0, 0.0, 0.0);
			    fprintf(fp, "%12.5e %12.5e   %12.5e %12.5e  ", 
			    		0.0, 0.0, 0.0, 0.0);
			}
			else
			{
			    fprintf(fp, "%12.5e %12.5e   %12.5e %12.5e   %12.5e %12.5e  ", 
					dens[index], pres[index], 
					vel[0][index],
					vel[1][index],
					gnor[0][index],
					gnor[1][index]);
			    fprintf(fp, "%12.5e %12.5e   %12.5e %12.5e  ",
			    		Gdens[0][index], Gpres[0][index],
					Gvel[0][0][index], Gvel[0][1][index]);
			    fprintf(fp, "%12.5e %12.5e   %12.5e %12.5e  ",
			    		Gdens[1][index], Gpres[1][index],
					Gvel[1][0][index], Gvel[1][1][index]);
			}

			fprintf(fp, "\n");
		  }
	}
	else
	{
	for(iz=imin[2]; iz <= imax[2]; iz++)
	    for(iy=imin[1]; iy <= imax[1]; iy++)
		  for(ix=imin[0]; ix <= imax[0]; ix++)
		  {
			index = d_index3d(ix,iy,iz,top_gmax);
			
			getRectangleCenter(index, coords);
			comp = cell_center[index].comp;

			fprintf(fp, "%f %f %f ", coords[0], coords[1], coords[2]);
			fprintf(fp, "%d ", comp);
		
			if(!gas_comp(comp))
			{
			    fprintf(fp, "%12.5e %12.5e  %12.5e %12.5e %12.5e  %12.5e %12.5e %12.5e ", 
					0.0,0.0, 0.0,0.0,0.0,  0.0,0.0,0.0);
			    fprintf(fp, "%12.5e %12.5e   %12.5e %12.5e %12.5e ", 
			    		0.0, 0.0, 0.0, 0.0, 0.0);
			    fprintf(fp, "%12.5e %12.5e   %12.5e %12.5e %12.5e ", 
			    		0.0, 0.0, 0.0, 0.0, 0.0);
			}
			else
			{
			    fprintf(fp, "%12.5e %12.5e   %12.5e %12.5e %12.5e  %12.5e %12.5e %12.5e ", 
					dens[index], pres[index], 
					vel[0][index],
					vel[1][index],
					vel[2][index],
					gnor[0][index],
					gnor[1][index],
					gnor[2][index]);
			    fprintf(fp, "%12.5e %12.5e   %12.5e %12.5e %12.5e ",
			    		Gdens[0][index], Gpres[0][index],
					Gvel[0][0][index], Gvel[0][1][index], Gvel[0][2][index]);
			    fprintf(fp, "%12.5e %12.5e   %12.5e %12.5e %12.5e ",
			    		Gdens[1][index], Gpres[1][index],
					Gvel[1][0][index], Gvel[1][1][index], Gvel[1][2][index]);
			}
			fprintf(fp, "\n");
		  }
	}

	fclose(fp);

}

EXPORT  void    tecplot_surface_states(
	const char	*bname,
	FILE		*file,
	SURFACE		*s)
{
	TRI	*tri;
	POINT	*p;
	int	i,npts,ntri;
	Locstate  sl, sr;

	if (bname != NULL)//direct call
	{
	    if ((file = fopen(bname,"w")) == NULL)
	    {
		screen("WARNING in tecplot_surface(), "
		       "can't open %s\n",bname);
		return;
	    }
	    (void) fprintf(file,"TITLE = \"tecplot surface\"\n"
		   	    "VARIABLES = \"x\", \"y\", \"z\", \"PL\", \"PR\", \"DL\", \"DR\" "
			    "\"u\", \"v\", \"w\", \"u1\", \"v1\", \"w1\" \n");
	}
	
	//called from tecplot_interface
	if (file == NULL)
	{
	    screen("ERROR, in tecplot_surface, file is NULL\n");
	    clean_up(ERROR);
	}
	if (!(first_tri(s)))
	{
	    screen("WARNING, first bond of the curve is NULL\n");
	    return;
	}

	//count number of points(npts) and number of tris(ntri)
	for (tri=first_tri(s),ntri=0; !at_end_of_tri_list(tri,s); tri=tri->next,ntri++)
	{
	    for (i = 0; i < 3; i++)
	    {
		Index_of_point(Point_of_tri(tri)[i]) = -1;
	    }
	}
	for (tri=first_tri(s),npts=0; !at_end_of_tri_list(tri,s); tri=tri->next)
	{
	    for (i = 0; i < 3; i++)
	    {
		p = Point_of_tri(tri)[i];
		if (Index_of_point(p) == -1)
		{
		    Index_of_point(p) = ++npts;
		}
	    }
	}
	//end counting
	
	fprint_wave_type(file, "ZONE T=\"", wave_type(s), "\"", s->interface);
    	fprintf(file, " N=%d E=%d\nF=FEPOINT, ET=TRIANGLE\n",npts,ntri);

	for (tri=first_tri(s); !at_end_of_tri_list(tri,s); tri=tri->next)
	{
	    for (i = 0; i < 3; i++)
	    {
		Index_of_point(Point_of_tri(tri)[i]) = -1;
	    }
	}
	for (tri=first_tri(s),npts=0; !at_end_of_tri_list(tri,s); tri=tri->next)
	{
	    for (i = 0; i < 3; i++)
	    {
		p = Point_of_tri(tri)[i];
		if (Index_of_point(p) == -1)
		{
		    Index_of_point(p) = ++npts;
		    
		    FT_GetStatesAtPoint(p,Hyper_surf_element(tri),Hyper_surf(s),&sl,&sr);
		    //sl = is_obstacle_state(sl) ? sr : sl;
		    //sr = is_obstacle_state(sr) ? sl : sr;
	     	   
		    //if(NO && wave_type(s) != PASSIVE_BOUNDARY)
		    //  fprintf(file,"%15.8e %15.8e %15.8e  %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e\n",Coords(p)[0],
		    //	   Coords(p)[1],Coords(p)[2], pressure(sl), pressure(sr), density(sl), density(sr), 
			//   vel(0,sl), vel(1,sl), vel(2,sl), vel(0,sr), vel(1,sr), vel(2,sr));
		        //fprintf(file,"%15.8e %15.8e %15.8e   %15.8e %15.8e %15.8e\n",Coords(p)[0],
		    	//   Coords(p)[1],Coords(p)[2], vel(0,sr), vel(1,sr), vel(2,sr));
		    //else
		        fprintf(file,"%15.8e %15.8e %15.8e  %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e\n",Coords(p)[0],
		    	 Coords(p)[1],Coords(p)[2], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
		}
	    }
	}
	for (tri=first_tri(s); !at_end_of_tri_list(tri,s); tri=tri->next)
	{
	    for (i = 0; i < 3; i++)
	    {
		fprintf(file,"%d ",Index_of_point(Point_of_tri(tri)[i]));
	    }
	    fprintf(file,"\n");
	}

	if (ntri != s->num_tri)
	{
	    printf("WARNING, num of tri in surface is wrong\n"); 
	}
	if (bname != NULL)
	    fclose(file);

}	/* end tecplot_surface_states */


EXPORT  void    tecplot_interface_states(const char*, INTERFACE	*);

EXPORT  void    tecplot_interface_states(
	const char	*bname,
	INTERFACE	*intfc)
{
	SURFACE	**s;
	char    bname1[200];
	FILE	*file;

	sprintf(bname1, "%s.plt", bname);
	if ((file = fopen(bname1,"w")) == NULL)
	{
	    screen("WARNING in tecplot_interface_states(), "
	           "can't open %s\n",bname1);
	    return;
	}
	(void) fprintf(file,"TITLE = \"tecplot interface\"\n"
		   	    "VARIABLES = \"x\", \"y\", \"z\", \"PL\", \"PR\", \"DL\", \"DR\" "
			    "\"u\", \"v\", \"w\", \"u1\", \"v1\", \"w1\" \n");

	for (s = intfc->surfaces; s && *s; ++s)
	{
	    tecplot_surface_states(NULL,file,*s);
	}
	fclose(file);
}	/* end tecplot_interface */

boolean G_CARTESIAN::get_ave_state(
	SWEEP 		m_vst,
	int		*ic,
	int		***norset,
	int		comp,
	int		ind)
{
	int		i, j, k, l, num, ic1[3], dir;
	float		gd, gp, gvel[3];
	boolean		found;
	double		**momn = m_vst.momn;
	double		*dens = m_vst.dens;
	double		*pres = m_vst.pres;
	double		***Gvel = eqn_params->Gvel;
	double		**Gdens = eqn_params->Gdens;
	double		**Gpres = eqn_params->Gpres;
	int		index, index0;
	int		icoords[MAXD];

	found = NO;

	for(i=0; i<dim; i++)
	    for(j=0; j<2; j++)
	    {
		dir = j == 0 ? -1 : 1;
		ft_assign(ic1, ic, 3*INT);
		ic1[i] = ic[i] + dir;

		if(ic1[i] < 0 || ic1[i] > top_gmax[i])
		    continue;
		if(norset[ic1[0]][ic1[1]][ic1[2]] == 1)
		    found = YES;
	    }

	if(!found)
	    return NO;

	index0 = d_index(ic,top_gmax,dim);

	num = 0;
	gd = 0.0;
	gp = 0.0;
	gvel[0] = 0.0;
	gvel[1] = 0.0;
	gvel[2] = 0.0;

	for (i = ic[0]-1; i <= ic[0]+1; i++)
	for (j = ic[1]-1; j <= ic[1]+1; j++)
	for (k = ic[2]-1; k <= ic[2]+1; k++)
	{
	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;
	    if(i < 0 || i > top_gmax[0] || 
	       j < 0 || j > top_gmax[1] ||
	       k < 0 || k > top_gmax[2])
		continue;
	    if(norset[i][j][k] != 1)
		continue;

	    index = d_index(icoords,top_gmax,dim);

	    if(cell_center[index].comp == comp)
	    {
		gd += dens[index];
		gp += pres[index];
		for (l = 0; l < dim; ++l)
		    gvel[l] += momn[l][index]/dens[index];
	    }
	    else
	    {
		gd += Gdens[ind][index];
		gp += Gpres[ind][index];
		for (l = 0; l < dim; ++l)
		    gvel[l] += Gvel[ind][l][index];
	    }
	    num++;
	}

	Gdens[ind][index0] = gd/num;
	Gpres[ind][index0] = gp/num;
	for (l = 0; l < dim; ++l)
	    Gvel[ind][l][index0] = gvel[l]/num;

	return YES;
}

/*
void G_CARTESIAN::get_ghost_state(
	SWEEP 		m_vst,
	int		comp,
	int		ind)
{
	int			i,j,k,l;
	int			ic[3],index;
	int			c, num;
	boolean			found;
	double			**momn = m_vst.momn;
	double			*dens = m_vst.dens;
	double			*pres = m_vst.pres;
	double			***Gvel = eqn_params->Gvel;
	double			**Gdens = eqn_params->Gdens;
	double			**Gpres = eqn_params->Gpres;
	static	int		***norset;


	if (norset == NULL)
	    FT_TriArrayMemoryAlloc((POINTER*)&norset,top_gmax[0]+1,
				top_gmax[1]+1,top_gmax[2]+1, INT);
	
	for (i = 0; i <= top_gmax[0]; i++)
	for (j = 0; j <= top_gmax[1]; j++)
	for (k = 0; k <= top_gmax[2]; k++)
	{
	    ic[0] = i;
	    ic[1] = j;
	    ic[2] = k;
	    index = d_index(ic,top_gmax,dim);
	    c = cell_center[index].comp;
		    
	    if(c == comp)
	    {
		norset[i][j][k] = 1;
		Gdens[ind][index] = dens[index];
		Gpres[ind][index] = pres[index];
		for (l = 0; l < dim; ++l)
		    Gvel[ind][l][index] = momn[l][index]/dens[index];
	    }
	    else
		norset[i][j][k] = 0;
	}

	found = YES;
	num = 1;
	while(found && num > 0)
	{
	    found = NO;

	    num = 0;
	    for (i = 0; i <= top_gmax[0]; i++)
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (k = 0; k <= top_gmax[2]; k++)
	    {
		if (norset[i][j][k] != 0)
		    continue;

		found = YES;
		ic[0] = i;
		ic[1] = j;
		ic[2] = k;

		if (get_ave_state(m_vst,ic,norset,comp,ind))
		{
		    num++;
		    norset[i][j][k] = 2;
		}
	    }

	    for (i = 0; i <= top_gmax[0]; i++)
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (k = 0; k <= top_gmax[2]; k++)
		if(norset[i][j][k] == 2)
		    norset[i][j][k] = 1;
	}
}
*/

void G_CARTESIAN::get_ghost_state(
	SWEEP 		m_vst,
	int		comp,
	int		ind)
{
	int			i,j,k;
	int			ic[3],index;
	int			c, num;
	boolean			found;
	double			**momn = m_vst.momn;
	double			*dens = m_vst.dens;
	double			*pres = m_vst.pres;
	double			***Gvel = eqn_params->Gvel;
	double			**Gdens = eqn_params->Gdens;
	double			**Gpres = eqn_params->Gpres;
	static	int		***norset;
	static 	int 		loop_count = 0;
	std::list<ToFill> resetThese;
	std::list<ToFill> fillThese;


	if (norset == NULL)
	{
	    int ft_vec_size = 0;
	    for (i = 0; i < dim; ++i)
	    {
		if (top_gmax[i]+8 > ft_vec_size)
		    ft_vec_size = top_gmax[i]+8;
	    }
	    if(dim == 1)
	    	FT_TriArrayMemoryAlloc((POINTER*)&norset,ft_vec_size,
				   1,1,INT);

	    if(dim == 2)
	    	FT_TriArrayMemoryAlloc((POINTER*)&norset,ft_vec_size,
				   ft_vec_size,1,INT);
	    if(dim == 3)
	    	FT_TriArrayMemoryAlloc((POINTER*)&norset,ft_vec_size,
				   ft_vec_size,ft_vec_size,INT);
	}

	ToFill aghst;
	
	for (i=0; i<=top_gmax[0]; i++)
	for (j=0; j<=top_gmax[1]; j++)
	for (k=0; k<=top_gmax[2]; k++)
	{
	    if (dim == 1)
	    	index = d_index1d(i,top_gmax);
	    else if (dim == 2)
	    	index = d_index2d(i,j,top_gmax);
	    else if (dim == 3)
	    	index = d_index3d(i,j,k,top_gmax);
	    c = cell_center[index].comp;

	    // for each cell that has component "comp" we
	    // set G* values and mark norset  for that cell to 1
	    if(c == comp)
	    {
		norset[i][j][k] = 1;
		Gdens[ind][index] = dens[index];
		Gpres[ind][index] = pres[index];
		Gvel[ind][0][index] = momn[0][index]/dens[index];
		if(dim > 1)
		    Gvel[ind][1][index] = momn[1][index]/dens[index];
		if(dim > 2)
		    Gvel[ind][2][index] = momn[2][index]/dens[index];
	    }
	    else
	    {
		aghst.icoords[0] = i;
		aghst.icoords[1] = j;
		aghst.icoords[2] = k;
		// hardcoded the stencil size to 4.... 
		if(withinStencilLen(aghst.icoords, 3) )
		{
		    fillThese.push_back(aghst);
		}
		norset[i][j][k] = 0;
	    }
	}

	found = YES;
	num = 1;
	while(found && (num > 0))
	{
	    std::list<ToFill>::iterator it;

	    found = NO;
	    loop_count++;
	    num = 0;

	    resetThese.clear();	
	    for (it=fillThese.begin() ; it != fillThese.end(); )
	    {
		found = YES;
		ic[0] = it->icoords[0]; 
		ic[1] = it->icoords[1]; 
		ic[2] = it->icoords[2]; 

		// if no neighbors are 1, return 0.
		if(get_ave_state(m_vst, ic,norset,comp,ind))
		{
		    num++;
		    norset[ ic[0] ][ ic[1] ][ ic[2] ] = 2;
 		    aghst.icoords[0] = ic[0];
		    aghst.icoords[1] = ic[1];
		    aghst.icoords[2] = ic[2]; 
		    resetThese.push_back(aghst);
		    it=fillThese.erase(it);// erase returns the next valid entery after the one we just erased.
		}
		else
		{
		    ++it;
		}
	    }
	    for (it=resetThese.begin(); it != resetThese.end(); it++)
		 norset[it->icoords[0]][it->icoords[1]][it->icoords[2]] = 1;
	}
	fillThese.clear();
	resetThese.clear();	
	loop_count = 0;
}

void G_CARTESIAN::GFMGhostState(
	int	*ic,
	int	comp,
	STATE	*ghost_st)
{
	int		i, index;
	double		ncor;
	double		***Gvel = eqn_params->Gvel;
	double		**Gdens = eqn_params->Gdens;
	double		**Gpres = eqn_params->Gpres;
	double		**Gnor = eqn_params->gnor;
	EOS_PARAMS	*eos = eqn_params->eos;

	index = d_index(ic,top_gmax,dim);
	ghost_st->eos = &(eos[comp]);
	ghost_st->dim = dim;

	ncor = 0.0;
	for(i=0; i<dim; i++)
	    ncor += (Gvel[1][i][index] - Gvel[0][i][index])*Gnor[i][index];
		    
	if(comp == 2)
	{
	    ghost_st->pres = Gpres[1][index];
	    ghost_st->dens = Gdens[0][index];
	    for(i=0; i<dim; i++)
		ghost_st->vel[i] = Gvel[0][i][index] + ncor*Gnor[i][index];
	}
	else
	{
	    ghost_st->pres = Gpres[0][index];
	    ghost_st->dens = Gdens[1][index];
	    for(i=0; i<dim; i++)
		ghost_st->vel[i] = Gvel[1][i][index] - ncor*Gnor[i][index];
	}
	for(i=0; i<dim; i++)
	    ghost_st->momn[i] = ghost_st->dens*ghost_st->vel[i];
	
	ghost_st->engy = EosEnergy(ghost_st);
}

void G_CARTESIAN::setNeumannStates(
	SWEEP		*vst,
	SWEEP		*m_vst,
	HYPER_SURF 	*hs,
	STATE		*state,
	int		*icoords,
	int		idir,
	int		nb,
	int		n,
	int		istart,
	COMPONENT	comp)
{
	int 		i,j,index;
	int             ind2[2][2] = {{0,1},{1,0}};
        int             ind3[3][3] = {{0,1,2},{1,2,0},{2,0,1}};
	int 		ic[MAXD];
	double		*vel_ref = state->vel;
	double		coords[MAXD],coords_ref[MAXD],crx_coords[MAXD];
	double		nor[MAXD],vn,v[MAXD];
	GRID_DIRECTION 	ldir[3] = {WEST,SOUTH,LOWER};
	GRID_DIRECTION 	rdir[3] = {EAST,NORTH,UPPER};
	GRID_DIRECTION  dir;
	STATE		st_tmp;

	st_tmp.eos = state->eos;
	st_tmp.dim = dim;
	index = d_index(icoords,top_gmax,dim);
	for (i = 0; i < dim; ++i)
	{
	    coords[i] = top_L[i] + icoords[i]*top_h[i];
	    ic[i] = icoords[i];
	}
	dir = (nb == 0) ? ldir[idir] : rdir[idir];
	FT_NormalAtGridCrossing(front,icoords,dir,comp,nor,&hs,crx_coords);

	if (debugging("neumann_buffer"))
	{
	    (void) printf("Entering setNeumannStates()\n");
	    (void) printf("comp = %d\n",comp);
	    (void) printf("icoords = %d %d %d\n",icoords[0],icoords[1],
				icoords[2]);
	    (void) printf("idir = %d nb = %d\n",idir,nb);
	    (void) printf("istart = %d nrad = %d n = %d\n",istart,nrad,n);
	    (void) print_general_vector("coords = ",coords,dim,"\n");
	    (void) print_general_vector("crx_coords = ",crx_coords,dim,"\n");
	    (void) print_general_vector("nor = ",nor,dim,"\n");
	    (void) print_general_vector("vel_ref = ",vel_ref,dim,"\n");
	}

	for (i = istart; i <= nrad; ++i)
	{
	    /* Find ghost point */
	    ic[idir] = (nb == 0) ? icoords[idir] - i : icoords[idir] + i;
	    for (j = 0; j < dim; ++j)
		coords_ref[j] = top_L[j] + ic[j]*top_h[j];

	    /* Reflect ghost point through intfc-mirror at crossing */
	    coords_ref[idir] = 2.0*crx_coords[idir] - coords_ref[idir];
	    vn = 0.0;
	    for (j = 0; j < dim; ++j)
	    {
		v[j] = coords_ref[j] - crx_coords[j];
		vn += v[j]*nor[j];
	    }
	    for (j = 0; j < dim; ++j)
		v[j] = 2.0*vn*nor[j] - v[j];
	    for (j = 0; j < dim; ++j)
		coords_ref[j] = crx_coords[j] + v[j];
			
	    /* Interpolate the state at the reflected point */
	    FT_IntrpStateVarAtCoords(front,comp,coords_ref,
		m_vst->dens,getStateDens,&st_tmp.dens,&m_vst->dens[index]);
	    FT_IntrpStateVarAtCoords(front,comp,coords_ref,
		m_vst->engy,getStateEngy,&st_tmp.engy,&m_vst->engy[index]);
	    FT_IntrpStateVarAtCoords(front,comp,coords_ref,
		m_vst->pres,getStatePres,&st_tmp.pres,&m_vst->pres[index]);
	    FT_IntrpStateVarAtCoords(front,comp,coords_ref,
			m_vst->momn[0],getStateXmom,&st_tmp.momn[0],
			&m_vst->momn[0][index]);
	    if (dim > 1)
		FT_IntrpStateVarAtCoords(front,comp,coords_ref,
			m_vst->momn[1],getStateYmom,&st_tmp.momn[1],
			&m_vst->momn[1][index]);
	    if (dim > 2)
		FT_IntrpStateVarAtCoords(front,comp,coords_ref,
			m_vst->momn[2],getStateZmom,&st_tmp.momn[2],
			&m_vst->momn[2][index]);
		/* Galileo Transformation */
	    vn = 0.0;
	    for (j = 0; j < dim; j++)
	    {
		v[j] = st_tmp.momn[j]/st_tmp.dens - vel_ref[j];
		vn += v[j]*nor[j];
	    }
	    for (j = 0; j < dim; j++)
	    {
		v[j] += vel_ref[j] - 2.0*vn*nor[j];
		st_tmp.momn[j] = v[j]*st_tmp.dens;
	    }
	    st_tmp.pres = EosPressure(&st_tmp);
	    if (st_tmp.pres < min_pres) st_tmp.pres = min_pres;
	    st_tmp.engy = EosEnergy(&st_tmp);

	    if (nb == 0)
	    {
		vst->dens[nrad-i] = st_tmp.dens;
		vst->engy[nrad-i] = st_tmp.engy;
		vst->pres[nrad-i] = st_tmp.pres;
	    	for (j = 0; j < 3; j++)
		    vst->momn[j][nrad-i] = 0.0;
		if (dim == 1)
		   vst->momn[0][nrad-i] = st_tmp.momn[0];
	    	else if (dim == 2)
		    for (j = 0; j < 2; j++)
		    	vst->momn[j][nrad-i] = 
				st_tmp.momn[ind2[idir][j]];
	    	else if (dim == 3)
		    for (j = 0; j < 3; j++)
		    	vst->momn[j][nrad-i] = 
				st_tmp.momn[ind3[idir][j]];
	    }
	    else
	    {
		/* Debug selectively!
		if (debugging("crx_reflection"))
		{
	            sprintf(fname,"intfc-%d-%d",count,i);
	            xgraph_2d_reflection(fname,front->grid_intfc,coords,
				crx_coords,coords_ref,nor);
		}
		*/
		vst->dens[n+nrad+i-1] = st_tmp.dens;
		vst->engy[n+nrad+i-1] = st_tmp.engy;
		vst->pres[n+nrad+i-1] = st_tmp.pres;
	    	for (j = 0; j < 3; j++)
		    vst->momn[j][n+nrad+i-1] = 0.0;
		if (dim == 1)
		   vst->momn[0][n+nrad+i-1] = st_tmp.momn[0];
	    	else if (dim == 2)
		    for (j = 0; j < 2; j++)
		    	vst->momn[j][n+nrad+i-1] = 
				st_tmp.momn[ind2[idir][j]];
	    	else if (dim == 3)
		    for (j = 0; j < 3; j++)
		    	vst->momn[j][n+nrad+i-1] = 
				st_tmp.momn[ind3[idir][j]];
	    }
	}
	if (debugging("neumann_buffer"))
	    (void) printf("Leaving setNeumannStates()\n");
}	/* end setNeumannStates */

void G_CARTESIAN::setDirichletStates(
	STATE		*crx_state,
	SWEEP		*vst,
	SWEEP		*m_vst,
	HYPER_SURF 	*hs,
	int		*icoords,
	int		dir,
	int		nb,
	int		n,
	int		istart)
{
	int		j, k, index;
	STATE 		*state;
	int		ind2[2][2] = {{0,1},{1,0}};
	int		ind3[3][3] = {{0,1,2},{1,2,0},{2,0,1}};

	if (nb == 0)
	{
	  if (boundary_state(hs) != NULL)
	  {
	    //preset state bdry
	    state = (STATE*)boundary_state(hs);
	    for (k = istart; k <= nrad; ++k)
	    {
		vst->dens[nrad-k] = state->dens;
		vst->engy[nrad-k] = state->engy;
		vst->pres[nrad-k] = state->pres;
		
		for (j = 0; j < 3; j++)
                    vst->momn[j][nrad-k] = 0.0;
		if (dim == 1)
		    vst->momn[0][nrad-k] = state->momn[0];
		else if (dim == 2)
		  for (j = 0; j < 2; j++)
		    vst->momn[j][nrad-k] = state->momn[ind2[dir][j]];
		else if (dim == 3)
		  for (j = 0; j < 3; j++)
		    vst->momn[j][nrad-k] = state->momn[ind3[dir][j]];
	    }
	  }
	  else if (boundary_state_function(hs) &&
              strcmp(boundary_state_function_name(hs),
	      "cF_flowThroughBoundaryState") == 0)
	  {
	    //flow through bdry
	    for (k = istart; k <= nrad; ++k)
	    {
		index = d_index(icoords,top_gmax, dim);
		vst->dens[nrad-k] = m_vst->dens[index];
		vst->engy[nrad-k] = m_vst->engy[index];
		vst->pres[nrad-k] = m_vst->pres[index];
		
		for (j = 0; j < 3; j++)
                    vst->momn[j][nrad-k] = 0.0;
		if (dim == 1)
		    vst->momn[0][nrad-k] = m_vst->momn[0][index];
		else if (dim == 2)
		  for (j = 0; j < 2; j++)
		    vst->momn[j][nrad-k] = m_vst->momn[ind2[dir][j]][index];
		else if (dim == 3)
		  for (j = 0; j < 3; j++)
		    vst->momn[j][nrad-k] = m_vst->momn[ind3[dir][j]][index];
	    }
	  }
	  else
	  {
	    (void) printf("Unimplemented Dirichlet boundary type!\n");
	    clean_up(ERROR);
	  }
	}
	else
	{
	  if (boundary_state(hs) != NULL)
	  {
	    state = (STATE*)boundary_state(hs);
	    for (k = istart; k <= nrad; ++k)
	    {
		vst->dens[n+nrad+k-1] = state->dens;
		vst->engy[n+nrad+k-1] = state->engy;
		vst->pres[n+nrad+k-1] = state->pres;
		
		for (j = 0; j < 3; j++)
                    vst->momn[j][n+nrad+k-1] = 0.0;
		if (dim == 1)
		    vst->momn[0][n+nrad+k-1] = state->momn[0];
		else if (dim == 2)
		  for (j = 0; j < 2; j++)
		    vst->momn[j][n+nrad+k-1] = state->momn[ind2[dir][j]];
		else if (dim == 3)
		  for (j = 0; j < 3; j++)
		    vst->momn[j][n+nrad+k-1] = state->momn[ind3[dir][j]];
	    }
	  }
	  else if (boundary_state_function(hs) &&
              strcmp(boundary_state_function_name(hs),
	      "cF_flowThroughBoundaryState") == 0)
	  {
	    for (k = istart; k <= nrad; ++k)
	    {
		index = d_index(icoords,top_gmax, dim);
		vst->dens[n+nrad+k-1] = m_vst->dens[index];
		vst->engy[n+nrad+k-1] = m_vst->engy[index];
		vst->pres[n+nrad+k-1] = m_vst->pres[index];
		
		for (j = 0; j < 3; j++)
                    vst->momn[j][n+nrad+k-1] = 0.0;
		if (dim == 1)
		    vst->momn[0][n+nrad+k-1] = m_vst->momn[0][index];
		else if (dim == 2)
		  for (j = 0; j < 2; j++)
		    vst->momn[j][n+nrad+k-1] = 
					m_vst->momn[ind2[dir][j]][index];
		else if (dim == 3)
		  for (j = 0; j < 3; j++)
		    vst->momn[j][n+nrad+k-1] = 
					m_vst->momn[ind3[dir][j]][index];
	    }
	  }
	  else
	  {
	    (void) printf("Unimplemented Dirichlet boundary type!\n");
	    clean_up(ERROR);
	  }
	}
}

void G_CARTESIAN::initSampleVelocity(char *in_name)
{
        FILE *infile;
	static SAMPLE *sample;
	char *sample_type;
	double *sample_line;

	infile = fopen(in_name,"r");
	FT_ScalarMemoryAlloc((POINTER*)&sample,sizeof(SAMPLE));
	sample_type = sample->sample_type;
	sample_line = sample->sample_coords;
	dim = front->rect_grid->dim;

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
}	/* end initSampleVelocity */

void G_CARTESIAN::checkCorrectForTolerance(STATE *state)
{
	if (state->dens < min_dens)
	    state->dens = min_dens;
	if (state->pres < min_pres)
	    state->pres = min_pres;
	state->engy = EosEnergy(state);
}	/* end checkCorrectForTolerance */

boolean G_CARTESIAN::needBufferFromIntfc(
	COMPONENT domain_comp,
	COMPONENT comp)
{
	if (eqn_params->tracked)
	    return (domain_comp != comp) ? YES : NO;
	else
	    return (gas_comp(comp)) ? NO : YES;
}	/* needBufferFromIntfc */


bool G_CARTESIAN::withinStencilLen( int *icrds, int stencil )
{
        int istart = std::max(0, icrds[0] - stencil);
        int jstart = std::max(0, icrds[1] - stencil);
        int kstart = std::max(0, icrds[2] - stencil);
        int iend = std::min(top_gmax[0],icrds[0]+stencil);
        int jend = std::min(top_gmax[1],icrds[1]+stencil);
        int kend = std::min(top_gmax[2],icrds[2]+stencil);

        int index  =  d_index(icrds,top_gmax,dim);
        int mycomp =  top_comp[index];

        int i,j,k;
        if (dim == 2)
        {
            for (i = istart; i <= iend; i++)
            for (j = jstart; j <= jend; j++)
            {
                int ic[3];
                ic[0] = i; ic[1] = j; ic[2] = 0;
                index  =  d_index(ic,top_gmax,dim);
                if(mycomp != top_comp[index])
                    return true;
            }
            return NO;
        }
        else if (dim == 3)
	{
            for (i = istart; i <= iend; i++)
            for (j = jstart; j <= jend; j++)
            for (k = kstart; k <= kend; k++)
	    {
                int ic[3];
                ic[0] = i; ic[1] = j; ic[2] = k;
                index  =  d_index(ic,top_gmax,dim);
                if(mycomp != top_comp[index])
                    return true;
	    }
            return NO;
	}
	return YES;
}

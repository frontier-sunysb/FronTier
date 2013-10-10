/*******************************************************************
 * 		         iFbasic.cpp
 *******************************************************************/
#include "iFluid.h"
#include "solver.h"

//----------------------------------------------------------------
//		L_RECTANGLE
//----------------------------------------------------------------

static double (*getStateVel[3])(POINTER) = {getStateXvel,getStateYvel,
                                        getStateZvel};

L_RECTANGLE::L_RECTANGLE(): comp(-1)
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
//               Incompress_Solver_Basis
//               Pure virtual class	
//--------------------------------------------------------------------------



//-------------------------------------------------------------------------------
//               Incompress_Solver_Smooth_Basis
//------------------------------------------------------------------------------
Incompress_Solver_Smooth_Basis::Incompress_Solver_Smooth_Basis(Front &front):front(&front)
{
}

//---------------------------------------------------------------
//	initMesh
// include the following parts
// 1) setup cell_center
//---------------------------------------------------------------
void Incompress_Solver_Smooth_Basis::initMesh(void)
{
	int i,j,k, index;
	double coords[MAXD];
	int num_cells;

	// init cell_center
	L_RECTANGLE       rectangle;

	if (debugging("trace"))
            (void) printf("Entering initMesh()\n");
	FT_MakeGridIntfc(front);
	setDomain();

	num_cells = 1;
	for (i = 0; i < dim; ++i)
	{
	    num_cells *= (top_gmax[i] + 1);
	}
	cell_center.clear();
	cell_center.insert(cell_center.end(),num_cells,rectangle);
	
	// setup vertices
	// left to right, down to up
	switch (dim)
	{
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
            (void) printf("Leaving initMesh()\n");
}

void Incompress_Solver_Smooth_Basis::setComponent(void)
{
	int i;
	static POINTER state;
        double coords[MAXD];
	int size = (int)cell_center.size();
	double **vel = field->vel;
	double *pres = field->pres;
	
	for (i = 0; i < size; i++)
	{
            cell_center[i].comp =
                        getComponent(cell_center[i].icoords);
        }
	if(state == NULL)
            FT_ScalarMemoryAlloc((POINTER*)&state,front->sizest);

	for (i = 0; i < size; i++)
	{
            if (cell_center[i].comp != -1 &&
                cell_center[i].comp != top_comp[i])
            {
		getRectangleCenter(i, coords);
                if (!FrontNearestIntfcState(front,coords,top_comp[i],
                                (POINTER)state))
                {
                    (void) printf("In setComponent()\n");
                    (void) printf("FrontNearestIntfcState() failed\n");
                    (void) printf("old_comp = %d new_comp = %d\n",
                                        cell_center[i].comp,top_comp[i]);
                    clean_up(ERROR);
                }
		vel[0][i] = getStateXvel(state);
		vel[1][i] = getStateYvel(state);
		if (dim == 3)
		    vel[3][i] = getStateZvel(state);
		pres[i] = getStatePres(state);
            }
            cell_center[i].comp = top_comp[i];
	}
}

void Incompress_Solver_Smooth_Basis::setIndexMap(void)
{
	static boolean first = YES;
	int i,j,k,ic,index;
	int llbuf[MAXD],uubuf[MAXD];
	int size = iupper - ilower;
	static int old_size;

	if (first)
	{
	    first = NO;
	    switch (dim)
	    {
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
	    old_size = size;
	}
	else if (old_size < size)
	{
	    switch (dim)
	    {
	    case 2:
		FT_FreeThese(1,I_to_ij);
		FT_MatrixMemoryAlloc((POINTER*)&I_to_ij,size,2,INT);
	    	break;
	    case 3:
		FT_FreeThese(1,I_to_ijk);
		FT_MatrixMemoryAlloc((POINTER*)&I_to_ijk,size,3,INT);
	    	break;
	    }
	    old_size = size;
	}

	index = 0;
	for (i = 0; i < dim; ++i)
	{
	    llbuf[i] = lbuf[i] != 0 ? lbuf[i] : 1;
	    uubuf[i] = ubuf[i] != 0 ? ubuf[i] : 1;
	}
	switch (dim)
	{
	case 2:
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
		    ij_to_I[i][j] = -1;
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		ic = d_index2d(i,j,top_gmax);
		if (domain_status[ic] != TO_SOLVE)
		    continue;
                if (cell_center[ic].comp != SOLID_COMP)
                {
                    ij_to_I[i][j] = index + ilower;
		    I_to_ij[index][0] = i;
                    I_to_ij[index][1] = j;
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
		if (domain_status[ic] != TO_SOLVE)
		    continue;
		if (cell_center[ic].comp != SOLID_COMP)
		{
                    ijk_to_I[i][j][k] = index + ilower;
		    I_to_ijk[index][0] = i;
                    I_to_ijk[index][1] = j;
                    I_to_ijk[index][2] = k;
                    index++;
                }
	    }
	    FT_ParallelExchCellIndex(front,llbuf,uubuf,(POINTER)ijk_to_I);
	    break;
	}
}	/* end setIndexMap */

// for initial condition: 
// 		setInitialCondition();	
// this function should be called before solve()
// for the source term of the momentum equation: 	
// 		computeSourceTerm();

void Incompress_Solver_Smooth_Basis::getVelocity(double *p, double *U)
{
	int i;
	double **vel = field->vel;

	if (!FT_IntrpStateVarAtCoords(front,NO_COMP,p,vel[0],getStateXvel,
					&U[0],NULL))
	{
	    for (i = 0; i < dim; ++i) U[i] = 0.0;
	    return;
	}
	FT_IntrpStateVarAtCoords(front,NO_COMP,p,vel[1],getStateYvel,&U[1],
					NULL);
	if (dim == 3)
	    FT_IntrpStateVarAtCoords(front,NO_COMP,p,vel[2],getStateZvel,
					&U[2],NULL);
}

void Incompress_Solver_Smooth_Basis::getRectangleIndex(int index, int &i, int &j)
{
	i = cell_center[index].icoords[0];
	j = cell_center[index].icoords[1];
}

void Incompress_Solver_Smooth_Basis::getRectangleIndex(int index, int &i, int &j, int &k)
{
	i = cell_center[index].icoords[0];
	j = cell_center[index].icoords[1];
	k = cell_center[index].icoords[2];
}


int Incompress_Solver_Smooth_Basis::getRectangleComponent(int index)
{	
	return getComponent(cell_center[index].icoords);
}

void Incompress_Solver_Smooth_Basis::getRectangleCenter(
	int index, 
	double *coords)
{
	int i;
	for (i = 0; i < dim; ++i)
	    coords[i] = cell_center[index].m_coords[i];
}

double Incompress_Solver_Smooth_Basis::getDistance(double *c0, double *c1)
{
	int i;
	double dist;
	dist = 0.0;
	for (i = 0; i < dim; ++i)
	    dist += sqr(c0[i] - c1[i]);
	dist = sqrt(dist);
	return dist;
}


// input : p[]
// output: q[]

void Incompress_Solver_Smooth_Basis::getNearestInterfacePoint(
	COMPONENT comp,
	double *p, 
	double *q,
	double *nor,		// Normal vector
	double *kappa)		// curvature
{
	INTERFACE *intfc = front->interf;
	int i,j,dim = front->rect_grid->dim;
	double t[MAXD],mag_nor;
	HYPER_SURF_ELEMENT *hse;
	HYPER_SURF *hs;
	POINT *pts[MAXD];
	BOND *bond;
	TRI *tri;
	static double **pts_nor,*pts_kappa;

	if (pts_kappa == NULL)
	{
	    FT_VectorMemoryAlloc((POINTER*)&pts_kappa,MAXD,FLOAT);
	    FT_MatrixMemoryAlloc((POINTER*)&pts_nor,MAXD,MAXD,FLOAT);
	}
	nearest_interface_point(p,comp,intfc,NO_BOUNDARIES,NULL,q,t,&hse,&hs);
	if (hse != NULL)
	{
	    switch (dim)
	    {
	    case 2:
		bond = Bond_of_hse(hse);
		pts[0] = bond->start;
		pts[1] = bond->end;
		for (i = 0; i < dim; ++i)
		{
		    GetFrontCurvature(pts[i],hse,hs,&pts_kappa[i],front);
		    GetFrontNormal(pts[i],hse,hs,pts_nor[i],front);
		}
		*kappa = (1.0 - t[0])*pts_kappa[0] + t[0]*pts_kappa[1];
		for (i = 0; i < dim; ++i)
		    nor[i] = (1.0 - t[0])*pts_nor[0][i] + t[0]*pts_nor[1][i];
		mag_nor = mag_vector(nor,dim);
		for (i = 0; i < dim; ++i)
		    nor[i] /= mag_nor;
		break;
	    case 3:
		tri = Tri_of_hse(hse);	
		for (i = 0; i < dim; ++i)
		{
		    pts[i] = Point_of_tri(tri)[i];
		    GetFrontCurvature(pts[i],hse,hs,&pts_kappa[i],front);
		    GetFrontNormal(pts[i],hse,hs,pts_nor[i],front);
		}
		for (i = 0; i < dim; ++i)
		{
		    *kappa = t[i]*pts_kappa[i];
		    for (j = 0; j < dim; ++j)
			nor[j] = t[i]*pts_nor[i][j];
		}
		mag_nor = mag_vector(nor,dim);
		for (i = 0; i < dim; ++i)
		    nor[i] /= mag_nor;
	    }
	}
	else
	{
	    *kappa = 0.0;
	    nor[0] = 1.0;
	    for (i = 1; i < dim; ++i) nor[i] = 0.0;
	}
}	/* end getNearestInterfacePoint */

int Incompress_Solver_Smooth_Basis::getComponent(
	double *p)
{
	return component(p,front->interf);
}

int Incompress_Solver_Smooth_Basis::getComponent(
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
	}
	return NO_COMP;
}

void Incompress_Solver_Smooth_Basis::save(char *filename)
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

void Incompress_Solver_Smooth_Basis::setDomain()
{
	static int current_size = 0;
	INTERFACE *grid_intfc;
	Table *T;
	int i,size;

	iFparams = (IF_PARAMS*)front->extra1;

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
	
	hmin = top_h[0];
	size = top_gmax[0]+1;
        for (i = 1; i < dim; ++i)
	{
            if (hmin > top_h[i]) hmin = top_h[i];
	    size *= (top_gmax[i]+1);
	}

	switch (dim)
	{
	case 2:
	    if (size > current_size)
	    {
		if (field != NULL)
		{
		    FT_FreeThese(15,field,array,source,diff_coeff,field->mu,
				field->rho,field->pres,field->phi,field->q,
				field->div_U,field->vort,field->vel,
				field->grad_q,field->f_surf,domain_status);
		}
		FT_ScalarMemoryAlloc((POINTER*)&field,sizeof(IF_FIELD));
		iFparams->field = field;
	    	FT_VectorMemoryAlloc((POINTER*)&array,size,sizeof(double));
	    	FT_VectorMemoryAlloc((POINTER*)&source,size,sizeof(double));
	    	FT_VectorMemoryAlloc((POINTER*)&diff_coeff,size,sizeof(double));
	    	FT_VectorMemoryAlloc((POINTER*)&field->mu,size,
					sizeof(double));
	    	FT_VectorMemoryAlloc((POINTER*)&field->rho,size,
					sizeof(double));
	    	FT_VectorMemoryAlloc((POINTER*)&field->pres,size,
					sizeof(double));
	    	FT_VectorMemoryAlloc((POINTER*)&field->phi,size,
					sizeof(double));
	    	FT_VectorMemoryAlloc((POINTER*)&field->q,size,
					sizeof(double));
	    	FT_VectorMemoryAlloc((POINTER*)&field->div_U,size,
					sizeof(double));
	    	FT_VectorMemoryAlloc((POINTER*)&field->vort,size,
					sizeof(double));
	    	FT_MatrixMemoryAlloc((POINTER*)&field->vel,2,size,
					sizeof(double));
	    	FT_MatrixMemoryAlloc((POINTER*)&field->grad_q,2,size,
					sizeof(double));
	    	FT_MatrixMemoryAlloc((POINTER*)&field->f_surf,2,size,
					sizeof(double));
	    	FT_VectorMemoryAlloc((POINTER*)&domain_status,size,INT);
		current_size = size;
	    }
	    imin = (lbuf[0] == 0) ? 1 : lbuf[0];
	    jmin = (lbuf[1] == 0) ? 1 : lbuf[1];
	    imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
	    jmax = (ubuf[1] == 0) ? top_gmax[1] - 1 : top_gmax[1] - ubuf[1];
	    break;
	case 3:
	    if (size > current_size)
	    {
		if (field != NULL)
		{
		    FT_FreeThese(14,field,array,source,diff_coeff,field->mu,
				field->rho,field->pres,field->phi,field->q,
				field->div_U,field->vel,
				field->grad_q,field->f_surf,domain_status);
		}
		FT_ScalarMemoryAlloc((POINTER*)&field,sizeof(IF_FIELD));
		iFparams->field = field;
	    	FT_VectorMemoryAlloc((POINTER*)&array,size,sizeof(double));
	    	FT_VectorMemoryAlloc((POINTER*)&source,size,sizeof(double));
		FT_VectorMemoryAlloc((POINTER*)&diff_coeff,size,sizeof(double));
	    	FT_VectorMemoryAlloc((POINTER*)&field->mu,size,
					sizeof(double));
	    	FT_VectorMemoryAlloc((POINTER*)&field->rho,size,
					sizeof(double));
	    	FT_VectorMemoryAlloc((POINTER*)&field->pres,size,
					sizeof(double));
	    	FT_VectorMemoryAlloc((POINTER*)&field->phi,size,
					sizeof(double));
	    	FT_VectorMemoryAlloc((POINTER*)&field->q,size,
					sizeof(double));
	    	FT_VectorMemoryAlloc((POINTER*)&field->div_U,size,
					sizeof(double));
	    	FT_MatrixMemoryAlloc((POINTER*)&field->vel,3,size,
					sizeof(double));
	    	FT_MatrixMemoryAlloc((POINTER*)&field->grad_q,3,size,
					sizeof(double));
	    	FT_MatrixMemoryAlloc((POINTER*)&field->f_surf,3,size,
					sizeof(double));
	    	FT_VectorMemoryAlloc((POINTER*)&domain_status,size,INT);
		current_size = size;
	    }
	    imin = (lbuf[0] == 0) ? 1 : lbuf[0];
	    jmin = (lbuf[1] == 0) ? 1 : lbuf[1];
	    kmin = (lbuf[2] == 0) ? 1 : lbuf[2];
	    imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
	    jmax = (ubuf[1] == 0) ? top_gmax[1] - 1 : top_gmax[1] - ubuf[1];
	    kmax = (ubuf[2] == 0) ? top_gmax[2] - 1 : top_gmax[2] - ubuf[2];
	    break;
	}
}

void Incompress_Solver_Smooth_Basis::scatMeshArray()
{
	FT_ParallelExchGridArrayBuffer(array,front);
}

void Incompress_Solver_Smooth_Basis::setGlobalIndex()
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
		if (cell_center[ic].comp == SOLID_COMP) continue;
		NLblocks++;
	    }
	    break;
	case 2:
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		ic = d_index2d(i,j,top_gmax);
		if (cell_center[ic].comp == SOLID_COMP) continue;
		if (domain_status[ic] != TO_SOLVE) continue;
		NLblocks++;
	    }
	    break;
	case 3:
	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		ic = d_index3d(i,j,k,top_gmax);
		if (cell_center[ic].comp == SOLID_COMP) continue;
		if (domain_status[ic] != TO_SOLVE) continue;
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

void Incompress_Solver_Smooth_Basis::printFrontInteriorStates(char *out_name)
{
	int i,j,k,l,index;
	char filename[100];
	FILE *outfile;
	double **vel = field->vel;

	sprintf(filename,"%s/state.ts%s",out_name,
			right_flush(front->step,7));
#if defined(__MPI__)
	if (pp_numnodes() > 1)
            sprintf(filename,"%s-nd%s",filename,right_flush(pp_mynode(),4));
#endif /* defined(__MPI__) */
	sprintf(filename,"%s-ifluid",filename);
	outfile = fopen(filename,"w");

	fluid_print_front_states(outfile,front);
	
	fprintf(outfile,"\nInterior ifluid states:\n");
	switch (dim)
	{
	case 2:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    {
		index = d_index2d(i,j,top_gmax);
	        fprintf(outfile,"%24.18g\n",field->rho[index]);
	        fprintf(outfile,"%24.18g\n",field->pres[index]);
	        fprintf(outfile,"%24.18g\n",field->phi[index]);
	        fprintf(outfile,"%24.18g\n",field->mu[index]);
	    	for (l = 0; l < dim; ++l)
	            fprintf(outfile,"%24.18g\n",vel[l][index]);
	    }
	    break;
	case 3:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    for (k = 0; k <= top_gmax[2]; ++k)
	    {
		index = d_index3d(i,j,k,top_gmax);
	        fprintf(outfile,"%24.18g\n",field->rho[index]);
	        fprintf(outfile,"%24.18g\n",field->pres[index]);
	        fprintf(outfile,"%24.18g\n",field->phi[index]);
	        fprintf(outfile,"%24.18g\n",field->mu[index]);
	    	for (l = 0; l < dim; ++l)
	            fprintf(outfile,"%24.18g\n",vel[l][index]);
	    }
	}
	fclose(outfile);
}

void Incompress_Solver_Smooth_Basis::readFrontInteriorStates(char *restart_name)
{
	FILE *infile;
	int i,j,k,l,index;
	char fname[100];
	double **vel = field->vel;


	m_rho[0] = iFparams->rho1;		
	m_rho[1] = iFparams->rho2;		
	m_mu[0] = iFparams->mu1;		
	m_mu[1] = iFparams->mu2;		
	m_sigma = iFparams->surf_tension;
	m_comp[0] = iFparams->m_comp1;
	m_comp[1] = iFparams->m_comp2;
	m_smoothing_radius = iFparams->smoothing_radius;
	mu_min = HUGE;
	for(i = 0; i < 2; i++)
	{
	    if(ifluid_comp(m_comp[i]))
	    {
		mu_min = std::min(mu_min,m_mu[i]);
		rho_min = std::min(rho_min,m_rho[i]);
	    }
	}

	sprintf(fname,"%s-ifluid",restart_name);
	infile = fopen(fname,"r");
	
	/* Initialize states at the interface */
	fluid_read_front_states(infile,front);

	FT_MakeGridIntfc(front);
	setDomain();

	next_output_line_containing_string(infile,"Interior ifluid states:");

	switch (dim)
	{
	case 2:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    {
		index = d_index2d(i,j,top_gmax);
	    	fscanf(infile,"%lf",&field->rho[index]);
	    	fscanf(infile,"%lf",&field->pres[index]);
	    	fscanf(infile,"%lf",&field->phi[index]);
	    	fscanf(infile,"%lf",&field->mu[index]);
		for (l = 0; l < dim; ++l)
	    	    fscanf(infile,"%lf",&vel[l][index]);
	    }
	    break;
	case 3:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    for (k = 0; k <= top_gmax[2]; ++k)
	    {
		index = d_index3d(i,j,k,top_gmax);
	    	fscanf(infile,"%lf",&field->rho[index]);
	    	fscanf(infile,"%lf",&field->pres[index]);
	    	fscanf(infile,"%lf",&field->phi[index]);
	    	fscanf(infile,"%lf",&field->mu[index]);
		for (l = 0; l < dim; ++l)
	    	    fscanf(infile,"%lf",&vel[l][index]);
	    }
	}
	fclose(infile);
	computeGradientQ();
	copyMeshStates();
	setAdvectionDt();
}


void Incompress_Solver_Smooth_Basis::setAdvectionDt()
{
	pp_global_max(&max_speed,1);
	if (max_speed != 0.0)
	    max_dt = hmin/max_speed;
	else
	    max_dt = HUGE;
	min_dt = 0.0000001*sqr(hmin)/mu_min;
	if (debugging("trace"))
	{
	    if (max_dt == HUGE)
	    	(void) printf("In setAdvectionDt: \n"
			"max_dt = HUGE min_dt = %24.18g\n",min_dt);
	    else
	    	(void) printf("In setAdvectionDt:\n"
			"max_dt = %24.18g min_dt = %24.18g\n",max_dt,min_dt);
	}
}	/* end setAdvectionDt */


void Incompress_Solver_Smooth_Basis::initMovieVariables()
{
	boolean set_bound = NO;
	FILE *infile = fopen(InName(front),"r");
	char string[100];
	double var_max,var_min;

	if (CursorAfterStringOpt(infile,"Type y to set movie bounds:"))
        {
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'Y' || string[0] == 'y')
                set_bound = YES;
        }

	    /* Begin hdf movies */
	switch (dim)
	{
	case 2:
	    CursorAfterString(infile,"Type y to make movie of pressure:");
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'Y' || string[0] == 'y')
	    {
		if (set_bound)
		{
		    CursorAfterString(infile,"Enter min and max pressure:");
                    fscanf(infile,"%lf %lf",&var_min,&var_min);
                    (void) printf("%f %f\n",var_min,var_max);
		}
		FT_AddHdfMovieVariable(front,set_bound,YES,SOLID_COMP,
				"pres",0,field->pres,getStatePres,
				var_max,var_min);
	    }
	    CursorAfterString(infile,"Type y to make movie of vorticity:");
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'Y' || string[0] == 'y')
	    {
		if (set_bound)
		{
		    CursorAfterString(infile,"Enter min and max vorticity:");
                    fscanf(infile,"%lf %lf",&var_min,&var_min);
                    (void) printf("%f %f\n",var_min,var_max);
		}
		FT_AddHdfMovieVariable(front,set_bound,YES,SOLID_COMP,
				"vort",0,field->vort,getStateVort,
				var_max,var_min);
	    }
	    CursorAfterString(infile,"Type y to make movie of velocity:");
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'Y' || string[0] == 'y')
	    {
		if (set_bound)
		{
		    CursorAfterString(infile,"Enter min and max velocity:");
                    fscanf(infile,"%lf %lf",&var_min,&var_min);
                    (void) printf("%f %f\n",var_min,var_max);
		}
		FT_AddHdfMovieVariable(front,set_bound,YES,SOLID_COMP,
				"xvel",0,field->vel[0],getStateXvel,
				var_max,var_min);
		FT_AddHdfMovieVariable(front,set_bound,YES,SOLID_COMP,
				"yvel",0,field->vel[1],getStateYvel,
				var_max,var_min);
	    }
	    break;
	case 3:
	    CursorAfterString(infile,"Type y to make yz cross section movie:");
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'Y' || string[0] == 'y')
	    {
	    	CursorAfterString(infile,"Type y to make movie of pressure:");
            	fscanf(infile,"%s",string);
            	(void) printf("%s\n",string);
            	if (string[0] == 'Y' || string[0] == 'y')
		    FT_AddHdfMovieVariable(front,set_bound,YES,SOLID_COMP,
				"pres-yz",0,field->pres,getStatePres,
				0.0,0.0);
	    	CursorAfterString(infile,"Type y to make movie of velocity:");
            	fscanf(infile,"%s",string);
            	(void) printf("%s\n",string);
            	if (string[0] == 'Y' || string[0] == 'y')
		{
		    FT_AddHdfMovieVariable(front,set_bound,YES,SOLID_COMP,
				"velo-yz-y",0,field->vel[1],getStateYvel,
				0.0,0.0);
		    FT_AddHdfMovieVariable(front,set_bound,YES,SOLID_COMP,
				"velo-yz-z",0,field->vel[2],getStateZvel,
				0.0,0.0);
		}
	    }
	    CursorAfterString(infile,"Type y to make xz cross section movie:");
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'Y' || string[0] == 'y')
	    {
	    	CursorAfterString(infile,"Type y to make movie of pressure:");
            	fscanf(infile,"%s",string);
            	(void) printf("%s\n",string);
            	if (string[0] == 'Y' || string[0] == 'y')
		    FT_AddHdfMovieVariable(front,set_bound,YES,SOLID_COMP,
				"pres-xz",1,field->pres,getStatePres,
				0.0,0.0);
	    	CursorAfterString(infile,"Type y to make movie of velocity:");
            	fscanf(infile,"%s",string);
            	(void) printf("%s\n",string);
            	if (string[0] == 'Y' || string[0] == 'y')
		{
		    FT_AddHdfMovieVariable(front,set_bound,YES,SOLID_COMP,
				"velo-xz-x",1,field->vel[0],getStateXvel,
				0.0,0.0);
		    FT_AddHdfMovieVariable(front,set_bound,YES,SOLID_COMP,
				"velo-xz-z",1,field->vel[2],getStateZvel,
				0.0,0.0);
		}
	    }
	    CursorAfterString(infile,"Type y to make xy cross section movie:");
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'Y' || string[0] == 'y')
	    {
	    	CursorAfterString(infile,"Type y to make movie of pressure:");
            	fscanf(infile,"%s",string);
            	(void) printf("%s\n",string);
            	if (string[0] == 'Y' || string[0] == 'y')
		    FT_AddHdfMovieVariable(front,set_bound,YES,SOLID_COMP,
				"pres-xy",2,field->pres,getStatePres,
				0.0,0.0);
	    	CursorAfterString(infile,"Type y to make movie of velocity:");
            	fscanf(infile,"%s",string);
            	(void) printf("%s\n",string);
            	if (string[0] == 'Y' || string[0] == 'y')
		{
		    FT_AddHdfMovieVariable(front,set_bound,YES,SOLID_COMP,
				"velo-xy-x",2,field->vel[0],getStateXvel,
				0.0,0.0);
		    FT_AddHdfMovieVariable(front,set_bound,YES,SOLID_COMP,
				"velo-xy-y",2,field->vel[1],getStateYvel,
				0.0,0.0);
		}
	    }
	}
	/* Added for vtk movie of vector field */
	CursorAfterString(infile,"Type y to make vector velocity field movie:");
        fscanf(infile,"%s",string);
        (void) printf("%s\n",string);
        if (string[0] == 'Y' || string[0] == 'y')
	    FT_AddVtkVectorMovieVariable(front,"VELOCITY",field->vel);

	fclose(infile);
}	/* end initMovieVariables */


void Incompress_Solver_Smooth_Basis::computeSubgridModel(void)
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
	double **vel = field->vel;

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
            u[index] = vel[0][index];
            v[index] = vel[1][index];
            getRectangleCenter(index, coords);
            co_coords_y[index] = coords[1] + (top_h[1]/2.0);
        }

        for (j = jmin; j <= jmax; j++)
        for (i = imin-2; i <= imax+2; i++)
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

        for (j = jmin; j <= (jmax/2); j++)
        {
            jj = (2*j)-1;
            for (i = imin-1; i <= (imax/2)+1; i++)
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
        for (j = jmin; j <= (jmax/2); j++)
        {
            jj = (2*j)-1;
            for (i = imin-1; i <= (imax/2)+1; i++)
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

        for (j = jmin; j <= (jmax/2); j++)
        {
            jj = (2*j)-1;
            for (i = imin-1; i <= (imax/2)+1; i++)
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

        for (j = jmin; j <= jmax; j++)
        for (i = imin-1; i <= imax+1; i++)
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

        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index0 = d_index2d(i,j,top_gmax);
            index1 = d_index2d(i-1,j,top_gmax);
            index2 = d_index2d(i+1,j,top_gmax);
            index3 = d_index2d(i,j-1,top_gmax);
            index4 = d_index2d(i,j+1,top_gmax);

            vel[0][index0] += -m_dt*(
                              ((tau00[index2]-tau00[index1])/(2.0*top_h[0])) + 
                                ((tau01[index4]-tau01[index3])/(2.0*top_h[1])));
            vel[1][index0] += -m_dt*(
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

//-------------------------------------------------------------------------------
//               Incompress_Solver_Smooth_2D_Basis
//------------------------------------------------------------------------------

double Incompress_Solver_Smooth_2D_Basis::getSmoothingFunction(double phi)
{
	// Heaviside function [1]
	if (phi < -m_smoothing_radius)	
	    return 0;
	else if (phi > m_smoothing_radius)
	    return 1;
	else
	    return 1.0/2 + phi/(2*m_smoothing_radius) + 
		   1/(2*PI)*sin(PI*phi/m_smoothing_radius);
}

double Incompress_Solver_Smooth_2D_Basis::getSmoothingFunctionD(double *center, double *point)
{
        if (fabs(center[0]-point[0]) < 2*top_h[0] && 
		fabs(center[1]-point[1]) < 2*top_h[1])
	    return ((1.0 + cos((PI*(center[0]-point[0]))/(2.0*top_h[0]))) 
		    *(1.0 + cos((PI*(center[1]-point[1]))/(2.0*top_h[1])))) 
		    /(16.0*top_h[0]*top_h[1]);
	else
            return 0.0;
}	/* end getSmoothingFunctionD in 2D */

double Incompress_Solver_Smooth_2D_Basis::smoothedDeltaFunction(double *p, double *center)
{
	int i;
	double len,d[MAXD],r,delta;

	for (i = 0; i < dim; ++i) d[i] = p[i] - center[i];
	len = mag_vector(d,dim);
	r = len/hmin/m_smoothing_radius;	
	delta = (fabs(r) > 1.0) ? 0.0 : 0.5*(1.0 + cos(PI*r));
	return delta;
}	/* end smoothedDeltaFunction */

double Incompress_Solver_Smooth_2D_Basis::smoothedStepFunction(double *p, double *center, int sign)
{
	int i;
	double dist,dp[MAXD],x,H;

	for (i = 0; i < dim; ++i) dp[i] = p[i] - center[i];
	dist = mag_vector(dp,dim);
	x = sign*dist/hmin/m_smoothing_radius;	
	if (x < -1.0) H = 0.0;
	else if (x > 1.0) H = 1.0;
	else H = 0.5*(1.0 + x) + 0.5/PI*sin(PI*x);
	return H;
}	/* end smoothedStepFunction */

void Incompress_Solver_Smooth_2D_Basis::sampleVelocity()
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
	double **vel = field->vel;
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
	sprintf(dirname, "%s/samples/sample-%d",out_name,step);
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
	    sprintf(sname, "%s/x-%d.xg",dirname,count);
            sfile = fopen(sname,"w");
	    fprintf(sfile,"Next\n");
	    fprintf(sfile,"color=%s\n",sample_color[count]);
	    fprintf(sfile,"thickness=1.5\n");
            for (j = jmin; j <= jmax; ++j)
            {
                index = d_index2d(i,j,top_gmax);
                var1 = vel[0][index];
                index = d_index2d(i+1,j,top_gmax);
                var2 = vel[0][index];
                var = (var1 + lambda*var2) / (1.0 + lambda);
                getRectangleCenter(index,coords);
                fprintf(sfile,"%20.14f   %20.14f\n",coords[1],var);
            }
            fclose(sfile);
	    sprintf(sname, "%s/y-%d.xg",dirname,count);
            sfile = fopen(sname,"w");
	    fprintf(sfile,"Next\n");
	    fprintf(sfile,"color=%s\n",sample_color[count]);
	    fprintf(sfile,"thickness=1.5\n");
            for (j = jmin; j <= jmax; ++j)
            {
                index = d_index2d(i,j,top_gmax);
                var1 = vel[1][index];
                index = d_index2d(i+1,j,top_gmax);
                var2 = vel[1][index];
                var = (var1 + lambda*var2) / (1.0 + lambda);
                getRectangleCenter(index,coords);
                fprintf(sfile,"%20.14f   %20.14f\n",coords[1],var);
            }
            fclose(sfile);
	    sprintf(sname, "%s/p-%d.xg",dirname,count);
            sfile = fopen(sname,"w");
	    fprintf(sfile,"Next\n");
	    fprintf(sfile,"color=%s\n",sample_color[count]);
	    fprintf(sfile,"thickness=1.5\n");
            for (j = jmin; j <= jmax; ++j)
            {
                index = d_index2d(i,j,top_gmax);
                var1 = field->pres[index];
                index = d_index2d(i+1,j,top_gmax);
                var2 = field->pres[index];
                var = (var1 + lambda*var2) / (1.0 + lambda);
                getRectangleCenter(index,coords);
                fprintf(sfile,"%20.14f   %20.14f\n",coords[1],var);
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
	    sprintf(sname, "%s/x-%d.xg",dirname,count);
            sfile = fopen(sname,"w");
	    fprintf(sfile,"Next\n");
	    fprintf(sfile,"color=%s\n",sample_color[count]);
	    fprintf(sfile,"thickness=1.5\n");
            for (i = imin; i <= imax; ++i)
            {
                index = d_index2d(i,j,top_gmax);
                var1 = vel[0][index];
                index = d_index2d(i,j+1,top_gmax);
                var2 = vel[0][index];
                var = (var1 + lambda*var2) / (1.0 + lambda);
                getRectangleCenter(index,coords);
                fprintf(sfile,"%20.14f   %20.14f\n",coords[0],var);
            }
            fclose(sfile);
	    sprintf(sname, "%s/y-%d.xg",dirname,count);
            sfile = fopen(sname,"w");
	    fprintf(sfile,"Next\n");
	    fprintf(sfile,"color=%s\n",sample_color[count]);
	    fprintf(sfile,"thickness=1.5\n");
            for (i = imin; i <= imax; ++i)
            {
                index = d_index2d(i,j,top_gmax);
                var1 = vel[1][index];
                index = d_index2d(i,j+1,top_gmax);
                var2 = vel[1][index];
                var = (var1 + lambda*var2) / (1.0 + lambda);
                getRectangleCenter(index,coords);
                fprintf(sfile,"%20.14f   %20.14f\n",coords[0],var);
            }
            fclose(sfile);
	    sprintf(sname, "%s/p-%d.xg",dirname,count);
            sfile = fopen(sname,"w");
	    fprintf(sfile,"Next\n");
	    fprintf(sfile,"color=%s\n",sample_color[count]);
	    fprintf(sfile,"thickness=1.5\n");
            for (i = imin; i <= imax; ++i)
            {
                index = d_index2d(i,j,top_gmax);
                var1 = field->pres[index];
                index = d_index2d(i,j+1,top_gmax);
                var2 = field->pres[index];
                var = (var1 + lambda*var2) / (1.0 + lambda);
                getRectangleCenter(index,coords);
                fprintf(sfile,"%20.14f   %20.14f\n",coords[0],var);
            }
            fclose(sfile);
            break;
        }
	count++;
}	/* end sampleVelocity2d */

void Incompress_Solver_Smooth_2D_Basis::setSmoothedProperties(void)
{
	boolean status;
	int i,j,l,index,sign;
	COMPONENT comp;
	double t[MAXD],force[MAXD];
	double center[MAXD],point[MAXD],H,D;
	HYPER_SURF_ELEMENT *hse;
	HYPER_SURF *hs;
	int range = (int)(m_smoothing_radius+1);
	double **f_surf = field->f_surf;
	double *mu = field->mu;
	double *rho = field->rho;

	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    index  = d_index2d(i,j,top_gmax);			
	    comp  = cell_center[index].comp;
	    if (!ifluid_comp(comp)) continue;

	    getRectangleCenter(index, center);
	    status = FT_FindNearestIntfcPointInRange(front,comp,center,point,
				t,&hse,&hs,range);

	    for (l = 0; l < dim; ++l) force[l] = 0.0;

	    if (status == YES && 
		ifluid_comp(positive_component(hs)) &&
		ifluid_comp(negative_component(hs)) &&
		positive_component(hs) != negative_component(hs))
	    {
		sign = (comp == m_comp[0]) ? -1 : 1;
		D = smoothedDeltaFunction(center,point);
		H = smoothedStepFunction(center,point,sign);
		mu[index] = m_mu[0] + (m_mu[1]-m_mu[0])*H;
		rho[index] = m_rho[0] + (m_rho[1]-m_rho[0])*H; 
		
		if (m_sigma != 0.0 && D != 0.0)
		{
		    surfaceTension(center,hse,hs,force,m_sigma);
		    for (l = 0; l < dim; ++l)
		    {
			force[l] /= -rho[index];
			f_surf[l][index] = force[l];
		    }
		}
	    }
	    else
	    {
		switch (comp)
		{
		case LIQUID_COMP1:
		    mu[index] = m_mu[0];
		    rho[index] = m_rho[0];
		    break;
		case LIQUID_COMP2:
		    mu[index] = m_mu[1];
		    rho[index] = m_rho[1];
		    break;
		}
	    }
	}
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{	
	    index  = d_index2d(i,j,top_gmax);
	    array[index] = mu[index];
	}
	scatMeshArray();
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{	
	    index  = d_index2d(i,j,top_gmax);
	    mu[index] = array[index];
	}
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    index  = d_index2d(i,j,top_gmax);
	    array[index] = rho[index];
	}
	scatMeshArray();
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index  = d_index2d(i,j,top_gmax);
	    rho[index] = array[index];
	}
	for (l = 0; l < dim; ++l)
	{
	    for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
	    {
	    	index  = d_index2d(i,j,top_gmax);
	    	array[index] = f_surf[l][index];
	    }
	    scatMeshArray();
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {
	    	index  = d_index2d(i,j,top_gmax);
	    	f_surf[l][index] = array[index];
	    }
	}
}	/* end setSmoothedProperties2d */

//-------------------------------------------------------------------------------
//               Incompress_Solver_Smooth_3D_Basis
//------------------------------------------------------------------------------

double Incompress_Solver_Smooth_3D_Basis::getSmoothingFunction(double phi)
{
	// Heaviside function [1]
	if (phi < -m_smoothing_radius)	
	    return 0;
	else if (phi > m_smoothing_radius)
	    return 1;
	else
	    return 1.0/2 + phi/(2*m_smoothing_radius) + 
		   1/(2*PI)*sin(PI*phi/m_smoothing_radius);
}

double Incompress_Solver_Smooth_3D_Basis::getSmoothingFunctionD(double *center, double *point)
{
        if (fabs(center[0]-point[0]) < 2*top_h[0] && 
		fabs(center[1]-point[1]) < 2*top_h[1] && 
		fabs(center[2]-point[2]) < 2*top_h[2])
	    return ((1.0 + cos((PI*(center[0]-point[0]))/(2.0*top_h[0]))) 
		    *(1.0 + cos((PI*(center[1]-point[1]))/(2.0*top_h[1]))) 
		    *(1.0 + cos((PI*(center[2]-point[2]))/(2.0*top_h[2])))) 
		/(64.0*top_h[0]*top_h[1]*top_h[2]);
	else
            return 0.0;
}	/* end getSmoothingFunctionD in 3D */

double Incompress_Solver_Smooth_3D_Basis::smoothedDeltaFunction(double *p, double *center)
{
	int i;
	double len,d[MAXD],r,delta;

	for (i = 0; i < dim; ++i) d[i] = p[i] - center[i];
	len = mag_vector(d,dim);
	r = len/hmin/m_smoothing_radius;	
	delta = (fabs(r) > 1.0) ? 0.0 : 0.5*(1.0 + cos(PI*r));
	return delta;
}	/* end smoothedDeltaFunction */

double Incompress_Solver_Smooth_3D_Basis::smoothedStepFunction(double *p, double *center, int sign)
{
	int i;
	double dist,dp[MAXD],x,H;

	for (i = 0; i < dim; ++i) dp[i] = p[i] - center[i];
	dist = mag_vector(dp,dim);
	x = sign*dist/hmin/m_smoothing_radius;	
	if (x < -1.0) H = 0.0;
	else if (x > 1.0) H = 1.0;
	else H = 0.5*(1.0 + x) + 0.5/PI*sin(PI*x);
	return H;
}	/* end smoothedStepFunction */

void Incompress_Solver_Smooth_3D_Basis::sampleVelocity()
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
	char dirname[256];
	double **vel = field->vel;
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

        sprintf(dirname, "%s/samples/sample-%d",out_name,step);
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
                    sprintf(sname, "%s/x-%d.xg",dirname,count);
                    sfile = fopen(sname,"w");
		    fprintf(sfile,"Next\n");
		    fprintf(sfile,"color=%s\n",sample_color[count]);
		    fprintf(sfile,"thickness=1.5\n");
                    for (k = kmin; k <= kmax; ++k)
                    {
                        index = d_index3d(i,j,k,top_gmax);
                        velo1 = vel[0][index];
                        index = d_index3d(i+1,j,k,top_gmax);
                        velo2 = vel[0][index];
                        velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        index = d_index3d(i,j+1,k,top_gmax);
                        velo1 = vel[0][index];
                        index = d_index3d(i+1,j+1,k,top_gmax);
                        velo2 = vel[0][index];
                        velo_tmp2 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        velo = (velo_tmp1 + lambda2*velo_tmp2)/(1.0 + lambda2);
                        getRectangleCenter(index,coords);
                        fprintf(sfile,"%20.14f   %20.14f\n",coords[2],velo);
                    }
                    fclose(sfile);

                    sprintf(sname,"%s/y-%d.xg",dirname,count);
                    sfile = fopen(sname,"w");
		    fprintf(sfile,"Next\n");
		    fprintf(sfile,"color=%s\n",sample_color[count]);
		    fprintf(sfile,"thickness=1.5\n");
                    for (k = kmin; k <= kmax; ++k)
                    {
                        index = d_index3d(i,j,k,top_gmax);
                        velo1 = vel[1][index];
                        index = d_index3d(i+1,j,k,top_gmax);
                        velo2 = vel[1][index];
                        velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        index = d_index3d(i,j+1,k,top_gmax);
                        velo1 = vel[1][index];
                        index = d_index3d(i+1,j+1,k,top_gmax);
                        velo2 = vel[1][index];
                        velo_tmp2 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        velo = (velo_tmp1 + lambda2*velo_tmp2)/(1.0 + lambda2);
                        getRectangleCenter(index,coords);
                        fprintf(sfile,"%20.14f   %20.14f\n",coords[2],velo);
                    }
                    fclose(sfile);

                    sprintf(sname,"%s/z-%d.xg",dirname,count);
                    sfile = fopen(sname,"w");
		    fprintf(sfile,"Next\n");
		    fprintf(sfile,"color=%s\n",sample_color[count]);
		    fprintf(sfile,"thickness=1.5\n");
                    for (k = kmin; k <= kmax; ++k)
                    {
                        index = d_index3d(i,j,k,top_gmax);
                        velo1 = vel[2][index];
                        index = d_index3d(i+1,j,k,top_gmax);
                        velo2 = vel[2][index];
                        velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        index = d_index3d(i,j+1,k,top_gmax);
                        velo1 = vel[2][index];
                        index = d_index3d(i+1,j+1,k,top_gmax);
                        velo2 = vel[2][index];
                        velo_tmp2 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        velo = (velo_tmp1 + lambda2*velo_tmp2)/(1.0 + lambda2);
                        getRectangleCenter(index,coords);
                        fprintf(sfile,"%20.14f   %20.14f\n",coords[2],velo);
                    }
                    fclose(sfile);

                    sprintf(sname,"%s/p-%d.xg",dirname,count);
                    sfile = fopen(sname,"w");
		    fprintf(sfile,"Next\n");
		    fprintf(sfile,"color=%s\n",sample_color[count]);
		    fprintf(sfile,"thickness=1.5\n");
                    for (k = kmin; k <= kmax; ++k)
                    {
                        index = d_index3d(i,j,k,top_gmax);
                        velo1 = field->pres[index];
                        index = d_index3d(i+1,j,k,top_gmax);
                        velo2 = field->pres[index];
                        velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        index = d_index3d(i,j+1,k,top_gmax);
                        velo1 = field->pres[index];
                        index = d_index3d(i+1,j+1,k,top_gmax);
                        velo2 = field->pres[index];
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
                    sprintf(sname, "%s/x-%d.xg",dirname,count);
                    sfile = fopen(sname,"w");
		    fprintf(sfile,"Next\n");
		    fprintf(sfile,"color=%s\n",sample_color[count]);
		    fprintf(sfile,"thickness=1.5\n");
                    for (j = jmin; j <= jmax; ++j)
                    {
                        index = d_index3d(i,j,k,top_gmax);
                        velo1 = vel[0][index];
                        index = d_index3d(i+1,j,k,top_gmax);
                        velo2 = vel[0][index];
                        velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        index = d_index3d(i,j,k+1,top_gmax);
                        velo1 = vel[0][index];
                        index = d_index3d(i+1,j,k+1,top_gmax);
                        velo2 = vel[0][index];
                        velo_tmp2 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        velo = (velo_tmp1 + lambda2*velo_tmp2)/(1.0 + lambda2);
                        getRectangleCenter(index,coords);
                        fprintf(sfile,"%20.14f   %20.14f\n",coords[1],velo);
                    }
                    fclose(sfile);

                    sprintf(sname, "%s/y-%d.xg",dirname,count);
                    sfile = fopen(sname,"w");
		    fprintf(sfile,"Next\n");
		    fprintf(sfile,"color=%s\n",sample_color[count]);
		    fprintf(sfile,"thickness=1.5\n");
                    for (j = jmin; j <= jmax; ++j)
                    {
                        index = d_index3d(i,j,k,top_gmax);
                        velo1 = vel[1][index];
                        index = d_index3d(i+1,j,k,top_gmax);
                        velo2 = vel[1][index];
                        velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        index = d_index3d(i,j,k+1,top_gmax);
                        velo1 = vel[1][index];
                        index = d_index3d(i+1,j,k+1,top_gmax);
                        velo2 = vel[1][index];
                        velo_tmp2 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        velo = (velo_tmp1 + lambda2*velo_tmp2)/(1.0 + lambda2);
                        getRectangleCenter(index,coords);
                        fprintf(sfile,"%20.14f   %20.14f\n",coords[1],velo);
                    }
                    fclose(sfile);

                    sprintf(sname, "%s/z-%d.xg",dirname,count);
                    sfile = fopen(sname,"w");
		    fprintf(sfile,"Next\n");
		    fprintf(sfile,"color=%s\n",sample_color[count]);
		    fprintf(sfile,"thickness=1.5\n");
                    for (j = jmin; j <= jmax; ++j)
                    {
                        index = d_index3d(i,j,k,top_gmax);
                        velo1 = vel[2][index];
                        index = d_index3d(i+1,j,k,top_gmax);
                        velo2 = vel[2][index];
                        velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        index = d_index3d(i,j,k+1,top_gmax);
                        velo1 = vel[2][index];
                        index = d_index3d(i+1,j,k+1,top_gmax);
                        velo2 = vel[2][index];
                        velo_tmp2 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        velo = (velo_tmp1 + lambda2*velo_tmp2)/(1.0 + lambda2);
                        getRectangleCenter(index,coords);
                        fprintf(sfile,"%20.14f   %20.14f\n",coords[1],velo);
                    }
                    fclose(sfile);

                    sprintf(sname, "%s/p-%d.xg",dirname,count);
                    sfile = fopen(sname,"w");
		    fprintf(sfile,"Next\n");
		    fprintf(sfile,"color=%s\n",sample_color[count]);
		    fprintf(sfile,"thickness=1.5\n");
                    for (j = jmin; j <= jmax; ++j)
                    {
                        index = d_index3d(i,j,k,top_gmax);
                        velo1 = field->pres[index];
                        index = d_index3d(i+1,j,k,top_gmax);
                        velo2 = field->pres[index];
                        velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        index = d_index3d(i,j,k+1,top_gmax);
                        velo1 = field->pres[index];
                        index = d_index3d(i+1,j,k+1,top_gmax);
                        velo2 = field->pres[index];
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
                    sprintf(sname, "%s/x-%d.xg",dirname,count);
                    sfile = fopen(sname,"w");
		    fprintf(sfile,"Next\n");
		    fprintf(sfile,"color=%s\n",sample_color[count]);
		    fprintf(sfile,"thickness=1.5\n");
                    for (i = imin; i <= imax; ++i)
                    {
                        index = d_index3d(i,j,k,top_gmax);
                        velo1 = vel[0][index];
                        index = d_index3d(i,j+1,k,top_gmax);
                        velo2 = vel[0][index];
                        velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        index = d_index3d(i,j,k+1,top_gmax);
                        velo1 = vel[0][index];
                        index = d_index3d(i,j+1,k+1,top_gmax);
                        velo2 = vel[0][index];
                        velo_tmp2 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        velo = (velo_tmp1 + lambda2*velo_tmp2)/(1.0 + lambda2);
                        getRectangleCenter(index,coords);
                        fprintf(sfile,"%20.14f   %20.14f\n",coords[0],velo);
                    }
                    fclose(sfile);

                    sprintf(sname, "%s/y-%d.xg",dirname,count);
                    sfile = fopen(sname,"w");
		    fprintf(sfile,"Next\n");
		    fprintf(sfile,"color=%s\n",sample_color[count]);
		    fprintf(sfile,"thickness=1.5\n");
                    for (i = imin; i <= imax; ++i)
                    {
                        index = d_index3d(i,j,k,top_gmax);
                        velo1 = vel[1][index];
                        index = d_index3d(i,j+1,k,top_gmax);
                        velo2 = vel[1][index];
                        velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        index = d_index3d(i,j,k+1,top_gmax);
                        velo1 = vel[1][index];
                        index = d_index3d(i,j+1,k+1,top_gmax);
                        velo2 = vel[1][index];
                        velo_tmp2 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        velo = (velo_tmp1 + lambda2*velo_tmp2)/(1.0 + lambda2);
                        getRectangleCenter(index,coords);
                        fprintf(sfile,"%20.14f   %20.14f\n",coords[0],velo);
                    }
                    fclose(sfile);

                    sprintf(sname, "%s/z-%d.xg",dirname,count);
                    sfile = fopen(sname,"w");
		    fprintf(sfile,"Next\n");
		    fprintf(sfile,"color=%s\n",sample_color[count]);
		    fprintf(sfile,"thickness=1.5\n");
                    for (i = imin; i <= imax; ++i)
                    {
                        index = d_index3d(i,j,k,top_gmax);
                        velo1 = vel[2][index];
                        index = d_index3d(i,j+1,k,top_gmax);
                        velo2 = vel[2][index];
                        velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        index = d_index3d(i,j,k+1,top_gmax);
                        velo1 = vel[2][index];
                        index = d_index3d(i,j+1,k+1,top_gmax);
                        velo2 = vel[2][index];
                        velo_tmp2 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        velo = (velo_tmp1 + lambda2*velo_tmp2)/(1.0 + lambda2);
                        getRectangleCenter(index,coords);
                        fprintf(sfile,"%20.14f   %20.14f\n",coords[0],velo);
                    }
                    fclose(sfile);

                    sprintf(sname, "%s/p-%d.xg",dirname,count);
                    sfile = fopen(sname,"w");
		    fprintf(sfile,"Next\n");
		    fprintf(sfile,"color=%s\n",sample_color[count]);
		    fprintf(sfile,"thickness=1.5\n");
                    for (i = imin; i <= imax; ++i)
                    {
                        index = d_index3d(i,j,k,top_gmax);
                        velo1 = field->pres[index];
                        index = d_index3d(i,j+1,k,top_gmax);
                        velo2 = field->pres[index];
                        velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        index = d_index3d(i,j,k+1,top_gmax);
                        velo1 = field->pres[index];
                        index = d_index3d(i,j+1,k+1,top_gmax);
                        velo2 = field->pres[index];
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
	count++;
}	/* end sampleVelocity in 3D */

void Incompress_Solver_Smooth_Basis::initSampleVelocity(char *in_name)
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
}	/* end initSampleVelocity */

void Incompress_Solver_Smooth_3D_Basis::setSmoothedProperties(void)
{
	boolean status;
	int i,j,k,l,index,sign; 
	COMPONENT comp;
        double t[MAXD],force[MAXD];
	double center[MAXD],point[MAXD],H,D;
	HYPER_SURF_ELEMENT *hse;
        HYPER_SURF *hs;
	double **f_surf = field->f_surf;
	double *mu = field->mu;
	double *rho = field->rho;

	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    index  = d_index3d(i,j,k,top_gmax);			
	    comp  = cell_center[index].comp;
	    if (!ifluid_comp(comp)) continue;

	    getRectangleCenter(index, center);
            status = FT_FindNearestIntfcPointInRange(front,comp,center,point,
				t,&hse,&hs,(int)m_smoothing_radius);
            for (l = 0; l < dim; ++l) force[l] = 0.0;

	    if (status  == YES && 
		ifluid_comp(positive_component(hs)) &&
                ifluid_comp(negative_component(hs)) && 
		positive_component(hs) != negative_component(hs))
	    {
		sign = (comp == m_comp[0]) ? -1 : 1;
                D = smoothedDeltaFunction(center,point);
                H = smoothedStepFunction(center,point,sign);
                mu[index] = m_mu[0] + (m_mu[1]-m_mu[0])*H;
                rho[index] = m_rho[0] + (m_rho[1]-m_rho[0])*H;

                if (m_sigma != 0.0 && D != 0.0)
                {
                    surfaceTension(center,hse,hs,force,m_sigma);
                    for (l = 0; l < dim; ++l)
                    {
                        force[l] /= -rho[index];
			f_surf[l][index] = force[l];
                    }
                }
	    }
	    else
	    {
		switch (comp)
		{
		case LIQUID_COMP1:
		    mu[index] = m_mu[0];
		    rho[index] = m_rho[0];
		    break;
		case LIQUID_COMP2:
		    mu[index] = m_mu[1];
		    rho[index] = m_rho[1];
		    break;
		}
	    }
	}
	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{	
	    index  = d_index3d(i,j,k,top_gmax);
	    array[index] = mu[index];
	}
	scatMeshArray();
	for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{	
	    index  = d_index3d(i,j,k,top_gmax);
	    mu[index] = array[index];
	}
	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    index  = d_index3d(i,j,k,top_gmax);
	    array[index] = rho[index];
	}
	scatMeshArray();
	for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index  = d_index3d(i,j,k,top_gmax);
	    rho[index] = array[index];
	}
	for (l = 0; l < dim; ++l)
	{
	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
	    {
	    	index  = d_index3d(i,j,k,top_gmax);
	    	array[index] = f_surf[l][index];
	    }
	    scatMeshArray();
	    for (k = 0; k <= top_gmax[2]; k++)
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {
	    	index  = d_index3d(i,j,k,top_gmax);
	    	f_surf[l][index] = array[index];
	    }
	}
}	/* end setSmoothedProperties in 3D */

// Flux of Riemann solution of Burgers equation u_t + uu_x = 0

double burger_flux(	
	double ul,
	double um,
	double ur)
{
	double u_Rl,u_Rr;
	if (ul < um)
	{
	    if (ul > 0.0) u_Rl = ul;
	    else if (um < 0.0) u_Rl = um;
	    else u_Rl = 0.0;
	}
	else
	{
	    if (ul + um > 0.0) u_Rl = ul;
	    else u_Rl = um;
	}

	if (um < ur)
	{
	    if (um > 0.0) u_Rr = um;
	    else if (ur < 0.0) u_Rr = ur;
	    else u_Rr = 0.0;
	}
	else
	{
	    if (um + ur > 0.0) u_Rr = um;
	    else u_Rr = ur;
	}
	return 0.5*(u_Rr*u_Rr - u_Rl*u_Rl);
}	/* end flux */

// Flux of Riemann solution of linear equation u_t + au_x = 0

double linear_flux(	
	double a,
	double ul,
	double um,
	double ur)
{
	if (a > 0.0)
	    return a*(um - ul);
	else
	    return a*(ur - um);
}	/* end net_uwind_flux */


void Incompress_Solver_Smooth_Basis::paintAllGridPoint(int status)
{
	int i,j,k,ic;	
	switch(dim)
	{
	case 2:
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		ic = d_index2d(i,j,top_gmax);
		domain_status[ic] = status;
	    }
	    break;
	case 3:
	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		ic = d_index3d(i,j,k,top_gmax);
		domain_status[ic] = status;
	    }
	    break;
	}
}	/* end paintAllGridPoint */

void Incompress_Solver_Smooth_Basis::paintSolvedGridPoint()
{
	int i,j,k,ic;	
	switch(dim)
	{
	case 2:
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		ic = d_index2d(i,j,top_gmax);
		if (domain_status[ic] == TO_SOLVE)
		    domain_status[ic] = SOLVED;
	    }
	    break;
	case 3:
	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		ic = d_index3d(i,j,k,top_gmax);
		if (domain_status[ic] == TO_SOLVE)
		    domain_status[ic] = SOLVED;
	    }
	    break;
	}
}	/* end paintSolvedGridPoin */

boolean Incompress_Solver_Smooth_Basis::paintToSolveGridPoint()
{
	GRID_DIRECTION  dir[6] = {EAST,NORTH,UPPER,WEST,SOUTH,LOWER};
	int idir,i,j,k,n,ic,*ip,ipn[MAXD];
	int ip_seed[MAXD];
	static int **ips,ip_size;
	boolean seed_found = NO;
	int paint_method;

	paint_method = BY_CROSSING;

	start_clock("paint_color");
	smin[0] = imin;	smin[1] = jmin;	smin[2] = kmin;
	smax[0] = imax;	smax[1] = jmax;	smax[2] = kmax;
	if (ips == NULL)
	{
	    ip_size = 1;
	    for (i = 0; i < dim; ++i)
	    	ip_size *= (smax[i] - smin[i] + 1); 
	    FT_MatrixMemoryAlloc((POINTER*)&ips,ip_size,dim,sizeof(int));
	}

	seed_found = NO;
	switch (dim)
	{
	case 2:
	    for (j = jmin; j <= jmax; j++)
            {
                for (i = imin; i <= imax; i++)
                {
                    ic = d_index2d(i,j,top_gmax);
                    if (domain_status[ic] == NOT_SOLVED)
                    {
                        ip_seed[0] = i;
                        ip_seed[1] = j;
                        domain_status[ic] = TO_SOLVE;
                        seed_found = YES;
                        break;
                    }
                }
                if (seed_found) break;
            }
	    break;
	case 3:
	    for (k = kmin; k <= kmax; k++)
	    {
	    	for (j = jmin; j <= jmax; j++)
	    	{
		    for (i = imin; i <= imax; i++)
		    {
		    	ic = d_index3d(i,j,k,top_gmax);	
		    	if (domain_status[ic] == NOT_SOLVED)
		    	{
			    ip_seed[0] = i;
			    ip_seed[1] = j;
			    ip_seed[2] = k;
			    domain_status[ic] = TO_SOLVE;
			    seed_found = YES;
			    break;
		    	}
		    }
		    if (seed_found) break;
	    	}
	    	if (seed_found) break;
	    }
	    break;
	}
	if (!seed_found)
	{
	    stop_clock("paint_color");
	    return NO;
	}

	/* Start traversing through connected neighbors */

	n = 0;
	for (i = 0; i < dim; ++i)
	    ips[n][i] = ip_seed[i];
	ip = ips[n++];

	while (n != 0)
	{
start_loop:
	    for (idir = 0; idir < dim*2; ++idir)
	    {
	    	if (nextConnectedPoint(ip,dir[idir],ipn,paint_method,smin,smax))
	    	{
		    for (i = 0; i < dim; ++i)
		    {
		    	ips[n][i] = ipn[i];
		    }
		    ic = d_index(ipn,top_gmax,dim);
		    domain_status[ic] = TO_SOLVE;
		    ip = ips[n];
		    n++;
		    if (n >= ip_size)
		    {
			(void) printf("Dead loop! Check code!\n");
			clean_up(ERROR);
		    }
		    goto start_loop;
	    	}
	    }
	    n--;
	    ip = ips[n];
	}
	stop_clock("paint_color");
	return YES;
}	/* end paintToSolveGridPoint */

boolean Incompress_Solver_Smooth_Basis::nextConnectedPoint(
	int *ip,
	GRID_DIRECTION dir,
	int *ipn,
	int paint_method,
	int *smin,
	int *smax)
{
	int id,idn;
	double crx_coords[MAXD];
	HYPER_SURF *hs;
	POINTER intfc_state;

	if (!next_ip_in_dir(ip,dir,ipn,smin,smax))
	    return NO;
	idn = d_index(ipn,top_gmax,dim);
	if (domain_status[idn] != NOT_SOLVED)
	    return NO;

	switch (paint_method)
	{
	case BY_COMPONENT:
	    id = d_index(ip,top_gmax,dim); 
	    idn = d_index(ipn,top_gmax,dim); 
	    if (top_comp[id] == top_comp[idn]) 
		return YES;
	    else 
		return NO;
	case BY_CROSSING:
	    id = d_index(ip,top_gmax,dim); 
	    if ((*findStateAtCrossing)(front,ip,dir,top_comp[id],
                                &intfc_state,&hs,crx_coords))
		return NO;
	    else
		return YES;
	case BY_WALL:
	    idn = d_index(ipn,top_gmax,dim); 
	    if (top_comp[idn] == SOLID_COMP) 
		return NO;
	    else 
		return YES;
	}
}	/* end connectedPoint */

void Incompress_Solver_Smooth_Basis::checkVelocityDiv(
	const char *mesg,
	double *vmin,
	double *vmax)
{
	double **vel = field->vel;
	double div_tmp,denom;
	double div_max = -HUGE;
	double div_min =  HUGE;
	int i,j,k,l,index;
	int icoords[MAXD];
	
	switch (dim)
	{
	case 2:
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
	    {
		icoords[0] = i;
                icoords[1] = j;
		index = d_index2d(i,j,top_gmax);
		div_tmp = computeFieldPointDiv(icoords,vel);
		if (div_max < div_tmp) div_max = div_tmp;
                if (div_min > div_tmp) div_min = div_tmp;
	    }
	    break;
	case 3:
	    for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
	    {
		icoords[0] = i;
                icoords[1] = j;
                icoords[2] = k;
		index = d_index3d(i,j,k,top_gmax);
		div_tmp = computeFieldPointDiv(icoords,vel);
		if (div_max < div_tmp) div_max = div_tmp;
                if (div_min > div_tmp) div_min = div_tmp;
	    }
	    break;
	}
	denom = 0.0;
	for (l = 0; l < dim; ++l)
	    denom += fabs((vmax - vmin)/(top_U[l] - top_L[l]));
	(void) printf("Check divergence at %s\n",mesg);
	(void) printf("Absolute divergence: max = %20.14f  div_min = %20.14f\n",
			div_max,div_min);
	(void) printf("Relative divergence: max = %20.14f  div_min = %20.14f\n",
			div_max/denom,div_min/denom);
}	/* end checkVelocityDiv */

void Incompress_Solver_Smooth_Basis::setDualDomain()
{
	static boolean first = YES;
	INTERFACE *comp_grid_intfc;
	Table *T;
	int i,size;

	comp_grid_intfc = front->comp_grid_intfc;
	ctop_grid = &topological_grid(comp_grid_intfc);
	ctop_gmax = ctop_grid->gmax;
	ctop_L = ctop_grid->L;
	ctop_U = ctop_grid->U;
	T = table_of_interface(comp_grid_intfc);
	ctop_comp = T->components;
	cimin = (lbuf[0] == 0) ? 1 : lbuf[0];
	cjmin = (lbuf[1] == 0) ? 1 : lbuf[1];
	ckmin = (lbuf[2] == 0) ? 1 : lbuf[2];
	cimax = (ubuf[0] == 0) ? ctop_gmax[0] - 1 : ctop_gmax[0] - ubuf[0] - 1;
	cjmax = (ubuf[1] == 0) ? ctop_gmax[1] - 1 : ctop_gmax[1] - ubuf[1] - 1;
	ckmax = (ubuf[2] == 0) ? ctop_gmax[2] - 1 : ctop_gmax[2] - ubuf[2] - 1;

	printf("In setDualDomain()\n");
	printf("Lower bound: cimin = %d  BL = %f\n",cimin,
			ctop_L[0]+cimin*top_h[0]);
	printf("Upper bound: cimax = %d  BL = %f\n",cimax,
			ctop_L[0]+cimax*top_h[0]);
	printf("Lower bound: cjmin = %d  BL = %f\n",cjmin,
			ctop_L[1]+cjmin*top_h[1]);
	printf("Upper bound: cjmax = %d  BL = %f\n",cjmax,
			ctop_L[1]+cjmax*top_h[1]);

	size = ctop_gmax[0]+1;
        for (i = 1; i < dim; ++i)
	    size *= (ctop_gmax[i]+1);
	printf("size = %d\n",size);
	switch (dim)
	{
	case 2:
	    if (first)
	    {
		iFparams->field = field;
	    	FT_VectorMemoryAlloc((POINTER*)&field->d_pres,size,
					sizeof(double));
	    	FT_VectorMemoryAlloc((POINTER*)&field->d_phi,size,
					sizeof(double));
	    	first = NO;
	    }
	    break;
	case 3:
	    if (first)
	    {
		iFparams->field = field;
	    	FT_VectorMemoryAlloc((POINTER*)&field->d_pres,size,
					sizeof(double));
	    	FT_VectorMemoryAlloc((POINTER*)&field->d_phi,size,
					sizeof(double));
	    	FT_VectorMemoryAlloc((POINTER*)&domain_status,size,INT);
	    	first = NO;
	    }
	}
}	/* end setDualDomain */

void Incompress_Solver_Smooth_Basis::setDualGlobalIndex()
{
	int i,j,k,ic;
	int num_nodes = pp_numnodes();
	int myid = pp_mynode();
	static boolean first = YES;

	if (first)
	{
	    first = NO;
	    FT_VectorMemoryAlloc((POINTER*)&cn_dist,num_nodes,sizeof(int));
	}
	cNLblocks = 0;
	switch (dim)
	{
	case 1:
	    for (i = cimin; i <= cimax; i++)
	    {
		ic = d_index1d(i,ctop_gmax);
		if (ctop_comp[ic] == SOLID_COMP) continue;
		cNLblocks++;
	    }
	    break;
	case 2:
	    for (j = cjmin; j <= cjmax; j++)
	    for (i = cimin; i <= cimax; i++)
	    {
		ic = d_index2d(i,j,ctop_gmax);
		if (ctop_comp[ic] == SOLID_COMP) continue;
		cNLblocks++;
	    }
	    break;
	case 3:
	    for (k = ckmin; k <= ckmax; k++)
	    for (j = cjmin; j <= cjmax; j++)
	    for (i = cimin; i <= cimax; i++)
	    {
		ic = d_index3d(i,j,k,ctop_gmax);
		if (ctop_comp[ic] == SOLID_COMP) continue;
		if (domain_status[ic] != TO_SOLVE) continue;
		cNLblocks++;
	    }
	    break;
	}

	for (i = 0; i < num_nodes; ++i) n_dist[i] = 0;
	cn_dist[myid] = cNLblocks;
	pp_global_imax(cn_dist,num_nodes);
	cilower = 0;
        ciupper = cn_dist[0];

        for (i = 1; i <= myid; i++)
        {
            cilower += cn_dist[i-1];
            ciupper += cn_dist[i];
        }	
}	/* setDualGlobalIndex */

void Incompress_Solver_Smooth_Basis::setDualIndexMap(void)
{
	static boolean first = YES;
	int i,j,k,ic,index;
	int llbuf[MAXD],uubuf[MAXD];
	int size = ciupper - cilower;
	static int old_size;

	if (first)
	{
	    first = NO;
	    switch (dim)
	    {
	    case 2:
	    	FT_MatrixMemoryAlloc((POINTER*)&cij_to_I,ctop_gmax[0]+1,
				ctop_gmax[1]+1,INT);
	    	break;
	    case 3:
	    	FT_TriArrayMemoryAlloc((POINTER*)&cijk_to_I,ctop_gmax[0]+1,
				ctop_gmax[1]+1,ctop_gmax[2]+1,INT);
	    	break;
	    }
	    old_size = size;
	}

	index = 0;
	for (i = 0; i < dim; ++i)
	{
	    llbuf[i] = lbuf[i] != 0 ? lbuf[i] : 1;
	    uubuf[i] = ubuf[i] != 0 ? ubuf[i] : 1;
	}
	switch (dim)
	{
	case 2:
	    for (j = 0; j <= ctop_gmax[1]; j++)
	    for (i = 0; i <= ctop_gmax[0]; i++)
		    cij_to_I[i][j] = -1;
	    for (j = cjmin; j <= cjmax; j++)
	    for (i = cimin; i <= cimax; i++)
	    {
		ic = d_index2d(i,j,ctop_gmax);
                if (ctop_comp[ic] != SOLID_COMP)
                {
                    cij_to_I[i][j] = index + cilower;
                    index++;
                }
	    }
	    FT_ParallelExchCellIndex(front,llbuf,uubuf,(POINTER)cij_to_I);
	    break;
	case 3:
	    for (k = 0; k <= ctop_gmax[2]; k++)
	    for (j = 0; j <= ctop_gmax[1]; j++)
	    for (i = 0; i <= ctop_gmax[0]; i++)
		    cijk_to_I[i][j][k] = -1;
	    for (k = ckmin; k <= ckmax; k++)
	    for (j = cjmin; j <= cjmax; j++)
	    for (i = cimin; i <= cimax; i++)
	    {
		ic = d_index3d(i,j,k,ctop_gmax);
		if (domain_status[ic] != TO_SOLVE)
		    continue;
                if (ctop_comp[ic] != SOLID_COMP)
		{
                    cijk_to_I[i][j][k] = index + cilower;
                    index++;
                }
	    }
	    FT_ParallelExchCellIndex(front,llbuf,uubuf,(POINTER)cijk_to_I);
	    break;
	}
}	/* end setDualIndexMap */

double Incompress_Solver_Smooth_Basis::computeFieldPointDiv(
        int *icoords,
        double **field)
{
	int icnb[MAXD];
        int i,j,index,index_nb;
        COMPONENT comp;
	double div,u_edge[3][2];
        double crx_coords[MAXD];
        POINTER intfc_state;
        HYPER_SURF *hs;
	int status;
        GRID_DIRECTION dir[3][2] = {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
	int idir,nb;

	index = d_index(icoords,top_gmax,dim);
        comp = top_comp[index];

	for (idir = 0; idir < dim; idir++)
	{
	    for (j = 0; j < dim; ++j)
	    	icnb[j] = icoords[j];
	    for (nb = 0; nb < 2; nb++)
	    {
	    	icnb[idir] = (nb == 0) ? icoords[idir] - 1 : icoords[idir] + 1;
	    	index_nb = d_index(icnb,top_gmax,dim);
	    	status = (*findStateAtCrossing)(front,icoords,dir[idir][nb],
			comp,&intfc_state,&hs,crx_coords);
	    	if (status == NO_PDE_BOUNDARY)
	    	    u_edge[idir][nb] = field[idir][index_nb];
	    	else if (status == CONST_V_PDE_BOUNDARY)
	    	    u_edge[idir][nb] = getStateVel[idir](intfc_state);
	    	else if (status == CONST_P_PDE_BOUNDARY)
	    	    u_edge[idir][nb] = field[idir][index];
	    }
	}

	div = 0.0;
	for (i = 0; i < dim; ++i)
	    div += 0.5*(u_edge[i][1] - u_edge[i][0])/top_h[i];
        return div;
}       /* end computeFieldPointDiv */

void Incompress_Solver_Smooth_Basis::computeFieldPointGrad(
        int *icoords,
        double *field,
        double *grad_field)
{
        int index,index_nb,icnb[MAXD];
        COMPONENT comp;
        int i,j,idir,nb;
	double p_edge[3][2],p0;
        double crx_coords[MAXD];
        POINTER intfc_state;
        HYPER_SURF *hs;
        GRID_DIRECTION dir[3][2] = {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
	int status;

	index = d_index(icoords,top_gmax,dim);
        comp = top_comp[index];
	p0 = field[index];

	for (idir = 0; idir < dim; idir++)
	{
	    for (j = 0; j < dim; ++j)
	    	icnb[j] = icoords[j];
	    for (nb = 0; nb < 2; nb++)
	    {
	    	icnb[idir] = (nb == 0) ? icoords[idir] - 1 : icoords[idir] + 1;
	    	index_nb = d_index(icnb,top_gmax,dim);
	    	status = (*findStateAtCrossing)(front,icoords,dir[idir][nb],
				comp,&intfc_state,&hs,crx_coords);
	    	if (status == NO_PDE_BOUNDARY)
		    p_edge[idir][nb] = field[index_nb];
	    	else if (status ==CONST_P_PDE_BOUNDARY)
		    p_edge[idir][nb] = getStatePhi(intfc_state);
	    	else if (status ==CONST_V_PDE_BOUNDARY)
		    p_edge[idir][nb] = p0;
	    }
	}
	for (i = 0; i < dim; ++i)
	    grad_field[i] = 0.5*(p_edge[i][1] - p_edge[i][0])/top_h[i];
}

void Incompress_Solver_Smooth_Basis::setReferencePressure()
{
	int i,j,k,index;
	double *pres = iFparams->field->pres;
        double diff_pres;

	switch(dim)
        {
	case 2:
	    index = d_index2d(imin,jmin,top_gmax);
            diff_pres = pres[index] - iFparams->ref_pres;
	    for (j = 0; j <= top_gmax[1]; ++j)
            for (i = 0; i <= top_gmax[0]; ++i)
	    {
		index = d_index2d(i,j,top_gmax);
		if (!ifluid_comp(top_comp[index]))
		    continue;
		pres[index] -= diff_pres;
	    }
	    break;
	case 3:
            index = d_index3d(imin,jmin,kmin,top_gmax);
            diff_pres = pres[index] - iFparams->ref_pres;
	    for (k = 0; k <= top_gmax[2]; ++k)
	    for (j = 0; j <= top_gmax[1]; ++j)
            for (i = 0; i <= top_gmax[0]; ++i)
	    {
		index = d_index3d(i,j,k,top_gmax);
		if (!ifluid_comp(top_comp[index]))
		    continue;
		pres[index] -= diff_pres;
	    }
	    break;
	}
}	/* end computeFieldPointGrad */

extern int ifluid_find_state_at_crossing(
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
	if (status == NO) 
	    return NO_PDE_BOUNDARY;
	if (wave_type(*hs) == FIRST_PHYSICS_WAVE_TYPE) 
	    return NO_PDE_BOUNDARY;
	if (wave_type(*hs) == NEUMANN_BOUNDARY) 
	    return CONST_V_PDE_BOUNDARY;
	if (wave_type(*hs) == GROWING_BODY_BOUNDARY) 
	    return CONST_V_PDE_BOUNDARY;
	if (wave_type(*hs) == MOVABLE_BODY_BOUNDARY) 
	    return CONST_V_PDE_BOUNDARY;
	if (wave_type(*hs) == ICE_PARTICLE_BOUNDARY) 
	    return CONST_V_PDE_BOUNDARY;
	if (wave_type(*hs) == DIRICHLET_BOUNDARY) 
	{
	    if (boundary_state(*hs))
	    	return CONST_V_PDE_BOUNDARY;
	    else
	    	return CONST_P_PDE_BOUNDARY;
	}
}	/* ifluid_find_state_at_crossing */

extern boolean neumann_type_bdry(int w_type)
{
	switch (w_type)
	{
	case NEUMANN_BOUNDARY:
	case MOVABLE_BODY_BOUNDARY:
	case GROWING_BODY_BOUNDARY:
	case ICE_PARTICLE_BOUNDARY:
	    return YES;
	default:
	    return NO;
	}
}	/* end neumann_type_bdry */


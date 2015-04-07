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

#include "solver.h"

#define		rho_of_comp(comp)				\
		(((comp) == soln_comp1) ? rho1 : rho2)

#define		soln_comp(comp)					\
		((comp) == soln_comp1 || (comp) == soln_comp2)

HYPERB_SOLVER::HYPERB_SOLVER(Front &front):front(&front)
{
	porosity = 0.0;
}

void HYPERB_SOLVER::solveRungeKutta()
{
	static boolean first = YES;
	int i,j;

	/* Allocate memory for Runge-Kutta of order */
	start_clock("solveRungeKutta");
	setSolverDomain();

	if (first)
	{
	    first = NO;
	    FT_VectorMemoryAlloc((POINTER*)&b,order,sizeof(double));
	    FT_MatrixMemoryAlloc((POINTER*)&a,order,order,sizeof(double));

	    FT_VectorMemoryAlloc((POINTER*)&st_field,order,sizeof(SWEEP));
	    FT_VectorMemoryAlloc((POINTER*)&st_flux,order,sizeof(FSWEEP));

	    FT_MatrixMemoryAlloc((POINTER*)&st_tmp.vel,dim,size,sizeof(double));
	    for (i = 0; i < order; ++i)
	    {
	    	FT_MatrixMemoryAlloc((POINTER*)&st_field[i].vel,dim,size,
				sizeof(double));
	    	FT_MatrixMemoryAlloc((POINTER*)&st_flux[i].vel_flux,dim,size,
				sizeof(double));
	    	FT_VectorMemoryAlloc((POINTER*)&st_field[i].rho,size,
				sizeof(double));
	    }
	    /* Set coefficient a, b, c for different order of RK method */
	    switch (order)
	    {
	    case 1:
		b[0] = 1.0;
		nrad = 1;
	    	break;
	    case 2:
	    	a[0][0] = 1.0;
	    	b[0] = 0.5;  b[1] = 0.5;
		nrad = 2;
	    	break;
	    case 4:
	    	a[0][0] = 0.5;
	    	a[1][0] = 0.0;  a[1][1] = 0.5;
	    	a[2][0] = 0.0;  a[2][1] = 0.0;  a[2][2] = 1.0;
	    	b[0] = 1.0/6.0;  b[1] = 1.0/3.0;
	    	b[2] = 1.0/3.0;  b[3] = 1.0/6.0;
		nrad = 3;
	    	break;
	    case 5:
		nrad = 3;
	    	break;
	    default:
	    	(void)printf("ERROR: %d-th order RK method not implemented\n",
					order);
	    	clean_up(ERROR);
	    }
	}

	/* Compute flux and advance field */

	copyToMeshVst(&st_field[0]);
	computeMeshFlux(st_field[0],&st_flux[0]);
	
	for (i = 0; i < order-1; ++i)
	{
	    copyMeshVst(st_field[0],&st_field[i+1]);
	    for (j = 0; j <= i; ++j)
	    {
		if (a[i][j] != 0.0)
		{
		    addMeshFluxToVst(&st_field[i+1],st_flux[j],a[i][j]);
		}
	    }
	    computeMeshFlux(st_field[i+1],&st_flux[i+1]);
	}
	for (i = 0; i < order; ++i)
	{
	    if (b[i] != 0.0)
	    {
		addMeshFluxToVst(&st_field[0],st_flux[i],b[i]);
	    }
	}
	copyFromMeshVst(st_field[0]);
	stop_clock("solveRungeKutta");
}	/* end solveRungeKutta */

void HYPERB_SOLVER::computeMeshFlux(
	SWEEP m_vst,
	FSWEEP *m_flux)
{
	int dir;

	resetFlux(m_flux);
	for (dir = 0; dir < dim; ++dir)
	{
	    addFluxInDirection(dir,&m_vst,m_flux);
	}
	addSourceTerm(&m_vst,m_flux);
}	/* end computeMeshFlux */

void HYPERB_SOLVER::resetFlux(FSWEEP *m_flux)
{
	int i,j;
	for (i = 0; i < size; i++)
	{
	    for (j = 0; j < dim; ++j)
	    	m_flux->vel_flux[j][i] = 0.0;
	}
}	/* resetFlux */

void HYPERB_SOLVER::addFluxInDirection(
	int dir,
	SWEEP *m_vst,
	FSWEEP *m_flux)
{
	switch (dim)
	{
	case 1:
	    return addFluxInDirection1d(dir,m_vst,m_flux);
	case 2:
	    return addFluxInDirection2d(dir,m_vst,m_flux);
	case 3:
	    return addFluxInDirection3d(dir,m_vst,m_flux);
	}
}	/* end addFluxInDirection */

void HYPERB_SOLVER::allocDirectionalVstFlux(
	SWEEP *vst,
	FSWEEP *flux)
{
	int i,size;
	size = 1;
	for (i = 0; i < dim; ++i)
            if (size < top_gmax[i]+2*nrad)
                size = top_gmax[i]+2*nrad;
	FT_MatrixMemoryAlloc((POINTER*)&vst->vel,dim,size,sizeof(double));
	FT_MatrixMemoryAlloc((POINTER*)&flux->vel_flux,dim,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&vst->rho,size,sizeof(double));
}	/* end allocDirectionalVstFlux */

void HYPERB_SOLVER::addFluxInDirection1d(
	int dir,
	SWEEP *m_vst,
	FSWEEP *m_flux)
{
}	/* end addFluxInDirection1d */

void HYPERB_SOLVER::addFluxInDirection2d(
	int dir,
	SWEEP *m_vst,
	FSWEEP *m_flux)
{
	static SWEEP vst;
	static FSWEEP vflux;
	static boolean first = YES;
	POINTER intfc_state;
	double crx_coords[MAXD];
	HYPER_SURF *hs;
	int i,j,l,n,seg_min,seg_max,index;
	int icoords[MAXD];
	int comp;
	double lambda = -dt/top_h[dir];

	if (first)
        {
            first = NO;
            allocDirectionalVstFlux(&vst,&vflux);
        }
	switch (dir)
	{
	case 0:
	    for (j = jmin; j <= jmax; j++)
	    {
		icoords[1] = j;
		seg_min = imin;
	    	for (i = 0; i <= top_gmax[0]; i++)
		    vst.rho[i] = 0.0;
		while (seg_min <= imax)
		{
		    for (; seg_min <= imax; ++seg_min)
		    {
			i = seg_min;
			icoords[0] = i;
		    	index = d_index2d(i,j,top_gmax);
		    	comp = top_comp[index];
			if (soln_comp(comp))
			    break;
		    }
		    if (seg_min > imax) break;
		    i = seg_min;
		    index = d_index2d(i,j,top_gmax);
		    comp = top_comp[index];
		    n = 0;
		    for (l = 0; l < dim; ++l)
                    	vst.vel[l][n+nrad] = m_vst->vel[l][index];
		    vst.rho[n+nrad] = rho_of_comp(top_comp[index]);
		    seg_max = i;
		    n++;
		    for (i = seg_min+1; i <= imax; i++)
            	    {
                	index = d_index2d(i,j,top_gmax);
		    	for (l = 0; l < dim; ++l)
                    	    vst.vel[l][n+nrad] = m_vst->vel[l][index];
		    	vst.rho[n+nrad] = rho_of_comp(top_comp[index]);
                    	n++;
                	seg_max = i;
			icoords[0] = i;
			comp = top_comp[index];
			if ((*findStateAtCrossing)(front,icoords,EAST,comp,
                                &intfc_state,&hs,crx_coords) != 
				NO_PDE_BOUNDARY)
			    break;
            	    }
		    icoords[1] = j;
		    icoords[0] = seg_min;
		    appendGhostBuffer(&vst,m_vst,n,icoords,0,0);
		    icoords[0] = seg_max;
		    appendGhostBuffer(&vst,m_vst,n,icoords,0,1);
		    
		    (*numericalFlux)(&vst,&vflux,lambda,n,dir,dim);
		    
		    n = 0;
		    for (i = seg_min; i <= seg_max; ++i)
		    {
		    	index = d_index2d(i,j,top_gmax);
		    	for (l = 0; l < dim; ++l)
		    	    m_flux->vel_flux[l][index] += vflux.vel_flux[l][n];
			n++;
		    }
		    seg_min = seg_max + 1;
		}
	    }
	    break;
	case 1:
	    for (i = imin; i <= imax; i++)
	    {
		icoords[0] = i;
		seg_min = jmin;
	    	for (j = 0; j <= top_gmax[1]; j++)
		    vst.rho[j] = 0.0;
		while (seg_min <= jmax)
		{
		    for (; seg_min <= jmax; ++seg_min)
		    {
			j = seg_min;
		    	index = d_index2d(i,j,top_gmax);
		    	comp = top_comp[index];
			if (soln_comp(comp))
			    break;
		    }
		    if (seg_min > jmax) break;
		    j = seg_min;
		    index = d_index2d(i,j,top_gmax);
		    comp = top_comp[index];
		    n = 0;
		    for (l = 0; l < dim; ++l)
                    	vst.vel[l][n+nrad] = m_vst->vel[l][index];
		    vst.rho[n+nrad] = rho_of_comp(top_comp[index]);
		    seg_max = j;
		    n++;
		    for (j = seg_min+1; j <= jmax; j++)
            	    {
                	index = d_index2d(i,j,top_gmax);
		    	for (l = 0; l < dim; ++l)
                    	    vst.vel[l][n+nrad] = m_vst->vel[l][index];
		    	vst.rho[n+nrad] = rho_of_comp(top_comp[index]);
                    	n++;
                	seg_max = j;
			icoords[1] = j;
			comp = top_comp[index];
			if ((*findStateAtCrossing)(front,icoords,NORTH,comp,
                                &intfc_state,&hs,crx_coords) != 
				NO_PDE_BOUNDARY)
			    break;
            	    }
		    icoords[0] = i;
		    icoords[1] = seg_min;
		    appendGhostBuffer(&vst,m_vst,n,icoords,1,0);
		    icoords[1] = seg_max;
		    appendGhostBuffer(&vst,m_vst,n,icoords,1,1);
		    
		    (*numericalFlux)(&vst,&vflux,lambda,n,dir,dim);
		    
		    n = 0;
		    for (j = seg_min; j <= seg_max; ++j)
		    {
		    	index = d_index2d(i,j,top_gmax);
		    	for (l = 0; l < dim; ++l)
		    	    m_flux->vel_flux[l][index] += vflux.vel_flux[l][n];
			n++;
		    }
		    seg_min = seg_max + 1;
		}
	    }
	    break;
	}
}	/* end addFluxInDirection2d */

void HYPERB_SOLVER::addFluxInDirection3d(
	int dir,
	SWEEP *m_vst,
	FSWEEP *m_flux)
{
	static SWEEP vst;
	static FSWEEP vflux;
	static boolean first = YES;
	POINTER intfc_state;
	double crx_coords[MAXD];
	HYPER_SURF *hs;
	int i,j,k,l,n,seg_min,seg_max,index;
	int icoords[MAXD];
	int comp;
	double lambda = -dt/top_h[dir];

	if (first)
        {
            first = NO;
            allocDirectionalVstFlux(&vst,&vflux);
        }
	switch (dir)
	{
	case 0:
	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    {
		icoords[2] = k;
		icoords[1] = j;
		seg_min = imin;
	    	for (i = 0; i <= top_gmax[0]; i++)
		    vst.rho[i] = 0.0;
		while (seg_min <= imax)
		{
		    for (; seg_min <= imax; ++seg_min)
		    {
			i = seg_min;
		    	index = d_index3d(i,j,k,top_gmax);
		    	comp = top_comp[index];
			if (soln_comp(comp))
			    break;
		    }
		    if (seg_min > imax) break;
		    i = seg_min;
		    index = d_index3d(i,j,k,top_gmax);
		    comp = top_comp[index];
		    n = 0;
		    for (l = 0; l < dim; ++l)
                    	vst.vel[l][n+nrad] = m_vst->vel[l][index];
		    vst.rho[n+nrad] = rho_of_comp(top_comp[index]);
		    seg_max = i;
		    n++;
		    for (i = seg_min+1; i <= imax; i++)
            	    {
                	index = d_index3d(i,j,k,top_gmax);
		    	for (l = 0; l < dim; ++l)
                    	    vst.vel[l][n+nrad] = m_vst->vel[l][index];
		    	vst.rho[n+nrad] = rho_of_comp(top_comp[index]);
                    	n++;
                	seg_max = i;
			icoords[0] = i;
			comp = top_comp[index];
			if ((*findStateAtCrossing)(front,icoords,EAST,comp,
                                &intfc_state,&hs,crx_coords) != 
				NO_PDE_BOUNDARY)
                            break;
            	    }
		    icoords[2] = k;
		    icoords[1] = j;
		    icoords[0] = seg_min;
		    appendGhostBuffer(&vst,m_vst,n,icoords,0,0);
		    icoords[0] = seg_max;
		    appendGhostBuffer(&vst,m_vst,n,icoords,0,1);
		    
		    (*numericalFlux)(&vst,&vflux,lambda,n,dir,dim);
		    
		    n = 0;
		    for (i = seg_min; i <= seg_max; ++i)
		    {
		    	index = d_index3d(i,j,k,top_gmax);
		    	for (l = 0; l < dim; ++l)
		    	    m_flux->vel_flux[l][index] += vflux.vel_flux[l][n];
			n++;
		    }
		    seg_min = seg_max + 1;
		}
	    }
	    break;
	case 1:
	    for (i = imin; i <= imax; i++)
	    for (k = kmin; k <= kmax; k++)
	    {
		icoords[0] = i;
		icoords[2] = k;
		seg_min = jmin;
	    	for (j = 0; j <= top_gmax[1]; j++)
		    vst.rho[j] = 0.0;
		while (seg_min <= jmax)
		{
		    for (; seg_min <= jmax; ++seg_min)
		    {
			j = seg_min;
		    	index = d_index3d(i,j,k,top_gmax);
		    	comp = top_comp[index];
			if (soln_comp(comp))
			    break;
		    }
		    if (seg_min > jmax) break;
		    j = seg_min;
		    index = d_index3d(i,j,k,top_gmax);
		    comp = top_comp[index];
		    n = 0;
		    for (l = 0; l < dim; ++l)
                    	vst.vel[l][n+nrad] = m_vst->vel[l][index];
		    vst.rho[n+nrad] = rho_of_comp(top_comp[index]);
		    seg_max = j;
		    n++;
		    for (j = seg_min+1; j <= jmax; j++)
            	    {
		    	index = d_index3d(i,j,k,top_gmax);
		    	for (l = 0; l < dim; ++l)
                    	    vst.vel[l][n+nrad] = m_vst->vel[l][index];
		    	vst.rho[n+nrad] = rho_of_comp(top_comp[index]);
                    	n++;
                	seg_max = j;
			icoords[1] = j;
			comp = top_comp[index];
                        if ((*findStateAtCrossing)(front,icoords,NORTH,comp,
                                &intfc_state,&hs,crx_coords) != 
				NO_PDE_BOUNDARY)
                            break;
            	    }
		    icoords[0] = i;
		    icoords[2] = k;
		    icoords[1] = seg_min;
		    appendGhostBuffer(&vst,m_vst,n,icoords,1,0);
		    icoords[1] = seg_max;
		    appendGhostBuffer(&vst,m_vst,n,icoords,1,1);
		    
		    (*numericalFlux)(&vst,&vflux,lambda,n,dir,dim);
		    
		    n = 0;
		    for (j = seg_min; j <= seg_max; ++j)
		    {
		    	index = d_index3d(i,j,k,top_gmax);
		    	for (l = 0; l < dim; ++l)
		    	    m_flux->vel_flux[l][index] += vflux.vel_flux[l][n];
			n++;
		    }
		    seg_min = seg_max + 1;
		}
	    }
	    break;
	case 2:
	    for (i = imin; i <= imax; i++)
	    for (j = jmin; j <= jmax; j++)
	    {
		icoords[0] = i;
		icoords[1] = j;
		seg_min = kmin;
	    	for (k = 0; k <= top_gmax[2]; k++)
		    vst.rho[k] = 0.0;
		while (seg_min <= kmax)
		{
		    for (; seg_min <= kmax; ++seg_min)
		    {
			k = seg_min;
		    	index = d_index3d(i,j,k,top_gmax);
		    	comp = top_comp[index];
			if (soln_comp(comp))
			    break;
		    }
		    if (seg_min > kmax) break;
		    k = seg_min;
		    index = d_index3d(i,j,k,top_gmax);
		    comp = top_comp[index];
		    n = 0;
		    for (l = 0; l < dim; ++l)
                    	vst.vel[l][n+nrad] = m_vst->vel[l][index];
		    vst.rho[n+nrad] = rho_of_comp(top_comp[index]);
		    seg_max = k;
		    n++;
		    for (k = seg_min+1; k <= kmax; k++)
            	    {
		    	index = d_index3d(i,j,k,top_gmax);
		    	for (l = 0; l < dim; ++l)
                    	    vst.vel[l][n+nrad] = m_vst->vel[l][index];
		    	vst.rho[n+nrad] = rho_of_comp(top_comp[index]);
                    	n++;
                	seg_max = k;
			icoords[2] = k;
			comp = top_comp[index];
			if ((*findStateAtCrossing)(front,icoords,UPPER,comp,
                                &intfc_state,&hs,crx_coords) != 
				NO_PDE_BOUNDARY)
                            break;
            	    }
		    icoords[0] = i;
		    icoords[1] = j;
		    icoords[2] = seg_min;
		    appendGhostBuffer(&vst,m_vst,n,icoords,2,0);
		    icoords[2] = seg_max;
		    appendGhostBuffer(&vst,m_vst,n,icoords,2,1);
		    
		    (*numericalFlux)(&vst,&vflux,lambda,n,dir,dim);
		    
		    n = 0;
		    for (k = seg_min; k <= seg_max; ++k)
		    {
		    	index = d_index3d(i,j,k,top_gmax);
		    	for (l = 0; l < dim; ++l)
		    	    m_flux->vel_flux[l][index] += vflux.vel_flux[l][n];
			n++;
		    }
		    seg_min = seg_max + 1;
		}
	    }
	    break;
	}
}	/* end addFluxInDirection3d */

void HYPERB_SOLVER::setSolverDomain(void)
{
	static boolean first = YES;
	RECT_GRID *rgr = &topological_grid(front->grid_intfc);
        struct Table *T = table_of_interface(front->grid_intfc);
        int *lbuf = front->rect_grid->lbuf;
        int *ubuf = front->rect_grid->ubuf;
	int i;

	dim = Dimension(front->interf);
        top_comp = T->components;
        top_gmax = rgr->gmax;
	top_h = rgr->h;
	top_L = rgr->L;
	if (first)
	{
	    first = NO;
	    size = 1;
	    for (i = 0; i < dim; ++i)
	    	size *= (top_gmax[i] + 1);
	    switch(dim)
	    {
	    case 1:
		imin = (lbuf[0] == 0) ? 1 : lbuf[0];
		imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
		FT_VectorMemoryAlloc((POINTER*)&array,size,sizeof(double));
		break;
	    case 2:
		imin = (lbuf[0] == 0) ? 1 : lbuf[0];
	    	jmin = (lbuf[1] == 0) ? 1 : lbuf[1];
		imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
	    	jmax = (ubuf[1] == 0) ? top_gmax[1] - 1 : top_gmax[1] - ubuf[1];
		FT_VectorMemoryAlloc((POINTER*)&array,size,sizeof(double));
		break;
	    case 3:
		imin = (lbuf[0] == 0) ? 1 : lbuf[0];
	    	jmin = (lbuf[1] == 0) ? 1 : lbuf[1];
		kmin = (lbuf[2] == 0) ? 1 : lbuf[2];
		imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
	    	jmax = (ubuf[1] == 0) ? top_gmax[1] - 1 : top_gmax[1] - ubuf[1];
		kmax = (ubuf[2] == 0) ? top_gmax[2] - 1 : top_gmax[2] - ubuf[2];
		FT_VectorMemoryAlloc((POINTER*)&array,size,sizeof(double));
		break;
	    }
	}
}	/* end setSolverDomain */

void HYPERB_SOLVER::copyToMeshVst(SWEEP *state)
{
	int i,j,k,l,index;
        switch (dim)
        {
        case 1:
            for (i = 0; i <= top_gmax[0]; ++i)
            {
                index = d_index1d(i,top_gmax);
		for (l = 0; l < dim; ++l)
		    state->vel[l][index] = var[l][index];
            }
            break;
	case 2:
            for (j = 0; j <= top_gmax[1]; ++j)
            for (i = 0; i <= top_gmax[0]; ++i)
            {
                index = d_index2d(i,j,top_gmax);
		for (l = 0; l < dim; ++l)
		    state->vel[l][index] = var[l][index];
		state->rho[index] = rho[index];
            }
            break;
	case 3:
            for (k = 0; k <= top_gmax[2]; ++k)
            for (j = 0; j <= top_gmax[1]; ++j)
            for (i = 0; i <= top_gmax[0]; ++i)
            {
                index = d_index3d(i,j,k,top_gmax);
		for (l = 0; l < dim; ++l)
		    state->vel[l][index] = var[l][index];
            }
        }
}

void HYPERB_SOLVER::copyMeshVst(
	SWEEP st1, 
	SWEEP *st2)
{
	int i,j,k,l,index;
        switch (dim)
        {
        case 1:
            for (i = 0; i <= top_gmax[0]; ++i)
            {
                index = d_index1d(i,top_gmax);
		for (l = 0; l < dim; ++l)
		    st2->vel[l][index] = st1.vel[l][index];
            }
            break;
	case 2:
            for (j = 0; j <= top_gmax[1]; ++j)
            for (i = 0; i <= top_gmax[0]; ++i)
            {
                index = d_index2d(i,j,top_gmax);
		for (l = 0; l < dim; ++l)
		    st2->vel[l][index] = st1.vel[l][index];
            }
            break;
	case 3:
            for (k = 0; k <= top_gmax[2]; ++k)
            for (j = 0; j <= top_gmax[1]; ++j)
            for (i = 0; i <= top_gmax[0]; ++i)
            {
                index = d_index3d(i,j,k,top_gmax);
		for (l = 0; l < dim; ++l)
		    st2->vel[l][index] = st1.vel[l][index];
            }
        }
}

void HYPERB_SOLVER::addMeshFluxToVst(
	SWEEP *state, 
	FSWEEP flux, 
	double chi)
{
        int             i,j,k,l,index;

	max_speed = -HUGE;
        switch (dim)
        {
        case 1:
            for (l = 0; l < dim; ++l)
            {
            	for (i = imin; i <= imax; ++i)
                {
		    index = d_index1d(i,top_gmax);
		    array[index] = state->vel[l][index] + 
				chi*flux.vel_flux[l][index];
                    if (max_speed < fabs(array[index]))
                    	max_speed = array[index];
                }
		FT_ParallelExchGridArrayBuffer(array,front,NULL);
            	for (i = 0; i <= top_gmax[0]; ++i)
		{
		    index = d_index1d(i,top_gmax);
		    state->vel[l][index] = array[index];
		}
            }
            break;
        case 2:
            for (l = 0; l < dim; ++l)
            {
                for (j = jmin; j <= jmax; ++j)
                for (i = imin; i <= imax; ++i)
                {
		    index = d_index2d(i,j,top_gmax);
		    array[index] = state->vel[l][index] + 
				chi*flux.vel_flux[l][index];
                    if (max_speed < fabs(array[index]))
                    	max_speed = array[index];
                }
		FT_ParallelExchGridArrayBuffer(array,front,NULL);
            	for (j = 0; j <= top_gmax[1]; ++j)
            	for (i = 0; i <= top_gmax[0]; ++i)
		{
		    index = d_index2d(i,j,top_gmax);
		    state->vel[l][index] = array[index];
		}
            }
            break;
        case 3:
            for (l = 0; l < dim; ++l)
            {
            	for (k = kmin; k <= kmax; ++k)
            	for (j = jmin; j <= jmax; ++j)
            	for (i = imin; i <= imax; ++i)
                {
		    index = d_index3d(i,j,k,top_gmax);
		    array[index] = state->vel[l][index] + 
				chi*flux.vel_flux[l][index];
                    if (max_speed < fabs(array[index]))
                    	max_speed = array[index];
                }
		FT_ParallelExchGridArrayBuffer(array,front,NULL);
            	for (k = 0; k <= top_gmax[2]; ++k)
            	for (j = 0; j <= top_gmax[1]; ++j)
            	for (i = 0; i <= top_gmax[0]; ++i)
		{
		    index = d_index3d(i,j,k,top_gmax);
		    state->vel[l][index] = array[index];
		}
            }
	}
}

void HYPERB_SOLVER::copyFromMeshVst(SWEEP state)
{
	int i,j,k,l,index;
	int c;
        switch (dim)
        {
        case 1:
            for (i = 0; i <= top_gmax[0]; ++i)
            {
                index = d_index1d(i,top_gmax);
		if (!soln_comp(top_comp[index])) continue;
		for (l = 0; l < dim; ++l)
		    soln[l][index] = state.vel[l][index];
            }
            break;
	case 2:
            for (j = 0; j <= top_gmax[1]; ++j)
            for (i = 0; i <= top_gmax[0]; ++i)
            {
                index = d_index2d(i,j,top_gmax);
		if (!soln_comp(top_comp[index])) continue;
		for (l = 0; l < dim; ++l)
		    soln[l][index] = state.vel[l][index];
            }
            break;
	case 3:
            for (k = 0; k <= top_gmax[2]; ++k)
            for (j = 0; j <= top_gmax[1]; ++j)
            for (i = 0; i <= top_gmax[0]; ++i)
            {
                index = d_index3d(i,j,k,top_gmax);
		if (!soln_comp(top_comp[index])) continue;
		for (l = 0; l < dim; ++l)
		    soln[l][index] = state.vel[l][index];
            }
        }
}

void HYPERB_SOLVER::appendGhostBuffer(
	SWEEP *vst,
	SWEEP *m_vst,
	int n,
	int *icoords,
	int idir,
	int nb)
{
	int             i,j,k,index,ic[MAXD],ic_next[MAXD];
        GRID_DIRECTION  ldir[3] = {WEST,SOUTH,LOWER};
        GRID_DIRECTION  rdir[3] = {EAST,NORTH,UPPER};
        HYPER_SURF      *hs;
        COMPONENT       comp;
        double          crx_coords[MAXD];
        POINTER         state;
	int		is_crxing;

	if (debugging("append_buffer"))
                printf("Entering appendGhostBuffer()\n");
	for (i = 0; i < dim; ++i) ic[i] = icoords[i];
	index = d_index(ic,top_gmax,dim);
        comp = top_comp[index];

	switch(nb)
	{
	case 0:
	    for (i = 1; i <= nrad; ++i)
	    {
		ic[idir] = icoords[idir] - i + 1;
		is_crxing = (*findStateAtCrossing)(front,ic,ldir[idir],comp,
                                	&state,&hs,crx_coords);
		    
		if (is_crxing == NO_PDE_BOUNDARY)
		{
		    ic[idir]--;
		    index = d_index(ic,top_gmax,dim);
		    for (j = 0; j < dim; j++)
		    {
			vst->vel[j][nrad-i] = m_vst->vel[j][index];
		    }
		    vst->rho[nrad-i] = rho_of_comp(top_comp[index]);
		}
		else
		{
		    switch (wave_type(hs))
		    {
		    case NEUMANN_BOUNDARY:
		    case MOVABLE_BODY_BOUNDARY:
		    case GROWING_BODY_BOUNDARY:
		    case ICE_PARTICLE_BOUNDARY:
		    	setNeumannStates(vst,m_vst,hs,state,ic,idir,
					nb,0,i,comp);
		    	break;
		    case DIRICHLET_BOUNDARY:
		    	setDirichletStates(vst,m_vst,hs,state,ic,idir,
					nb,0,i,comp);
		    	break;
		    case ELASTIC_BOUNDARY:
		    	setElasticStates(vst,m_vst,hs,state,ic,idir,
					nb,0,i,comp);
		    	break;
		    default:
		    	(void) printf("In appendGhostBuffer(): ");
		    	(void) print_wave_type("Unknown wave type ",
					wave_type(hs),"\n",front->interf);
		    	(void) print_int_vector("icoords=", icoords,dim,"\n");
			print_curve(Curve_of_hs(hs));
		    	clean_up(ERROR);
		    }
		    break;
		}
	    }
	    break;
	case 1:
	    for (i = 1; i <= nrad; ++i)
	    {
		ic[idir] = icoords[idir] + i - 1;
		is_crxing = (*findStateAtCrossing)(front,ic,rdir[idir],comp,
                                	&state,&hs,crx_coords);
		if (is_crxing == NO_PDE_BOUNDARY)
		{
		    ic[idir]++;
		    index = d_index(ic,top_gmax,dim);
		    for (j = 0; j < dim; j++)
			vst->vel[j][n+nrad+i-1] = m_vst->vel[j][index];
		    vst->rho[n+nrad+i-1] = rho_of_comp(top_comp[index]);
		}
		else
		{
		    switch (wave_type(hs))
		    {
		    case NEUMANN_BOUNDARY:
		    case MOVABLE_BODY_BOUNDARY:
		    case GROWING_BODY_BOUNDARY:
		    case ICE_PARTICLE_BOUNDARY:
		    	setNeumannStates(vst,m_vst,hs,state,ic,idir,
					nb,n,i,comp);
		    	break;
		    case DIRICHLET_BOUNDARY:
		    	setDirichletStates(vst,m_vst,hs,state,ic,idir,
					nb,n,i,comp);
		    	break;
		    case ELASTIC_BOUNDARY:
		    	setElasticStates(vst,m_vst,hs,state,ic,idir,
					nb,n,i,comp);
		    	break;
		    default:
		    	(void) printf("In appendGhostBuffer(): ");
		    	(void) print_wave_type("Unknown wave type ",
				wave_type(hs),"\n",front->interf);
		    	(void) print_int_vector("icoords=",icoords,dim,"\n");
		    	(void) printf("nb = %d\n",nb);
		    	clean_up(ERROR);
		    }
		    break;
		}
	    }
	}

	if (debugging("append_buffer"))
	    (void) printf("Leaving appendGhostBuffer()\n");
}	/* end appendGhostBuffer */

void HYPERB_SOLVER::setNeumannStates(
	SWEEP *vst, 
	SWEEP *m_vst,
	HYPER_SURF *hs,
	POINTER state,
	int *icoords,
	int idir,
	int nb,
	int n,
	int istart,
	int comp)
{
	int 		i,j,index;
	int 		ic[MAXD];
	double		coords[MAXD],coords_ref[MAXD],crx_coords[MAXD];
	double		nor[MAXD],vn,v[MAXD];
	GRID_DIRECTION 	ldir[3] = {WEST,SOUTH,LOWER};
	GRID_DIRECTION 	rdir[3] = {EAST,NORTH,UPPER};
	GRID_DIRECTION  dir;
	double	vel_ref[MAXD],v_tmp[MAXD];

	index = d_index(icoords,top_gmax,dim);
	for (i = 0; i < dim; ++i)
	{
	    vel_ref[i] = (*getStateVel[i])(state);
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
	    for (j = 0; j < dim; ++j)
	    	FT_IntrpStateVarAtCoords(front,comp,coords_ref,m_vst->vel[j],
			getStateVel[j],&v_tmp[j],&m_vst->vel[j][index]);

		/* Galileo Transformation */
	    vn = 0.0;
	    for (j = 0; j < dim; j++)
	    {
		/* Relative velocity of reflected point to boundary */
		v[j] = v_tmp[j] - vel_ref[j];
		/* Normal component of the relative velocity */
		vn += v[j]*nor[j];
	    }
	    /*
	    for (j = 0; j < dim; j++)
	    {
		v[j] += vel_ref[j] - 2.0*vn*nor[j];
		v_tmp[j] = v[j];
	    }
	    */
	    /* Only normal component is reflected, 
	       relative tangent velocity is zero */
	    for (j = 0; j < dim; j++)
		v_tmp[j] = vel_ref[j] - 2.0*vn*nor[j];

	    if (nb == 0)
	    {
	    	for (j = 0; j < dim; j++)
		    vst->vel[j][nrad-i] = (1.0 - porosity)*v_tmp[j] +
                        porosity*m_vst->vel[j][d_index(ic,top_gmax,dim)];
		vst->rho[nrad-i] = rho_of_comp(comp);
	    }
	    else
	    {
	    	for (j = 0; j < dim; j++)
		    vst->vel[j][n+nrad+i-1] = (1.0 - porosity)*v_tmp[j] +
                        porosity*m_vst->vel[j][d_index(ic,top_gmax,dim)];
		vst->rho[n+nrad+i-1] = rho_of_comp(comp);
	    }
	}
	if (debugging("neumann_buffer"))
	    (void) printf("Leaving setNeumannStates()\n");
}

void HYPERB_SOLVER::setDirichletStates(
	SWEEP		*vst,
	SWEEP		*m_vst,
	HYPER_SURF 	*hs,
	POINTER 	state,
	int		*icoords,
	int		dir,
	int		nb,
	int		n,
	int		istart,
	int 		comp)
{
	int		j, k, index;

	if (nb == 0)
	{
	  if (boundary_state(hs) != NULL)
	  {
	    //preset state bdry
	    state = (POINTER)boundary_state(hs);
	    for (k = istart; k <= nrad; ++k)
	    {
		for (j = 0; j < dim; j++)
		    vst->vel[j][nrad-k] = (*getStateVel[j])(state);
		vst->rho[nrad-k] = rho_of_comp(comp);
	    }
	  }
	  else if (boundary_state_function(hs) &&
              strcmp(boundary_state_function_name(hs),
              "iF_splitBoundaryState") == 0)
	  {
	    for (k = istart; k <= nrad; ++k)
	    {
		for (j = 0; j < dim; j++)
		    vst->vel[j][nrad-k] = (*getStateVel[j])(state);
		vst->rho[nrad-k] = rho_of_comp(comp);
	    }
	  }
	  else if (boundary_state_function(hs) &&
              strcmp(boundary_state_function_name(hs),
	      "flowThroughBoundaryState") == 0)
	  {
	    //flow through bdry
	    for (k = istart; k <= nrad; ++k)
	    {
		index = d_index(icoords,top_gmax, dim);
		
		for (j = 0; j < dim; j++)
		    vst->vel[j][nrad-k] = m_vst->vel[j][index];
		vst->rho[nrad-k] = rho_of_comp(comp);
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
	    state = (POINTER)boundary_state(hs);
	    for (k = istart; k <= nrad; ++k)
	    {
		for (j = 0; j < dim; j++)
                    vst->vel[j][n+nrad+k-1] = (*getStateVel[j])(state);
		vst->rho[n+nrad+k-1] = rho_of_comp(comp);
	    }
	  }
	  else if (boundary_state_function(hs) &&
              strcmp(boundary_state_function_name(hs),
	      "flowThroughBoundaryState") == 0)
	  {
	    for (k = istart; k <= nrad; ++k)
	    {
		index = d_index(icoords,top_gmax, dim);
		
		for (j = 0; j < dim; j++)
		    vst->vel[j][n+nrad+k-1] = m_vst->vel[j][index];
		vst->rho[n+nrad+k-1] = rho_of_comp(comp);
	    }
	  }
	  else
	  {
	    (void) printf("Unimplemented Dirichlet boundary type!\n");
	    clean_up(ERROR);
	  }
	}
}

void HYPERB_SOLVER::setElasticStates(
	SWEEP		*vst,
	SWEEP		*m_vst,
	HYPER_SURF 	*hs,
	POINTER 	state,
	int		*icoords,
	int		dir,
	int		nb,
	int		n,
	int		istart,
	int 		comp)
{
	int		j, k, index;

	if (nb == 0)
	{
	    for (k = istart; k <= nrad; ++k)
	    {
		for (j = 0; j < dim; j++)
		    vst->vel[j][nrad-k] = (*getStateVel[j])(state);
		vst->rho[nrad-k] = rho_of_comp(comp);
	    }
	}
	else
	{
	    for (k = istart; k <= nrad; ++k)
	    {
		for (j = 0; j < dim; j++)
                    vst->vel[j][n+nrad+k-1] = (*getStateVel[j])(state);
		vst->rho[n+nrad+k-1] = rho_of_comp(comp);
	    }
	}
}

void HYPERB_SOLVER::addSourceTerm(SWEEP *state, FSWEEP * flux)
{
}

void upwind_flux(
	SWEEP *vst,
	FSWEEP *vflux,
	double lambda,
	int n,
	int dir,
	int dim)
{
	static int max_n = 0;
	static double **q;
	int i,j;
	double **vel = vst->vel;
	double *rho = vst->rho;
	double **vel_flux = vflux->vel_flux;
	double ul,ur,a;
	
	if (debugging("upwind"))
	{
	    (void) printf("Entering upwind_flux()\n");
	    (void) printf("dir = %d  n = %d\n",dir,n);
	}
	if (max_n < n+2)
	{
	    if (max_n != 0)
		FT_FreeThese(1,q);
	    FT_MatrixMemoryAlloc((POINTER*)&q,dim,n+1,sizeof(double));
	    max_n = n+2;
	}
	for (i = 0; i < dim; ++i)
	{
	    for (j = 0; j < n+1; ++j)
	    {
		if (i == dir)
		{
		    ul = vel[i][j];
		    ur = vel[i][j+1];
		    if (ul < ur) // Rarefaction
		    {
			if (ul > 0.0) q[i][j] = 0.5*sqr(ul);
			else if (ur < 0.0) q[i][j] = 0.5*sqr(ur);
			else q[i][j] = 0.0;
		    }
		    else	// Shock
		    {
			q[i][j] = (ul+ur > 0.0) ? 0.5*sqr(ul) : 0.5*sqr(ur);
		    }
		}
	   }
	}
	
	for (i = 0; i < dim; ++i)
	{
	    for (j = 0; j < n; ++j)
	    {
		vel_flux[i][j] = lambda*(q[i][j+1] - q[i][j]);
	    }
	}
	for (i = 0; i < dim; ++i)
	{
	    	if (i == dir) continue;
	    	for (j = 0; j < n; ++j)
	    	{
		    double um;

		    a = vel[dir][j+1];
		    ul = vel[i][j];
		    um = vel[i][j+1];
		    ur = vel[i][j+2];
		    if (a > 0.0)
		    	vel_flux[i][j] = lambda*a*(um - ul);	
		    else
		    	vel_flux[i][j] = lambda*a*(ur - um);	
	    	}
	}
}	/* end upwind_flux */


double (*flux_func)(double,double);
double (*dF)(double,double);

static double linear_flux(double,double);
static double burger_flux(double,double);
static double linear_dF(double,double);
static double burger_dF(double,double);

static double linear_flux(
	double wave_speed,
	double u)
{
        return u;
}

static double linear_dF(
	double wave_speed,
	double u)
{
	return 1.0;
        //return wave_speed;
}

static double burger_flux(
	double wave_speed,
	double u)
{
        return 0.5*u*u;
}

static double burger_dF(
	double wave_speed,
	double u)
{
        return u;
}

void Weno5_Get_Flux(
	double *u, 
	double *v, 
	double *flux, 
	double lambda,
	int n)
{
	int nrad = 3; /* 5th-order weno */
	int i,j,k;
	int extend_size = n + 2*nrad;
	const double eps = 1.e-8;
	const int p = 2;
	static double *f; /* f(u(x)) */
	static double *fp, *fm, *fm_burg;
	static int max_n = 0;
	double max_df, norm;
	double aa;

	/* coefficients for 2rd-order ENO interpolation stencil */
	double a[3][3] = {{1.0/3.0, -7.0/6.0, 11.0/6.0}, 
			  {-1.0/6.0, 5.0/6.0, 1.0/3.0}, 
			  {1.0/3.0, 5.0/6.0, -1.0/6.0}};

	/* Optimal weights C_k */
	double c[3] = {0.1,0.6,0.3};

	double is[3]; /* a smoothness measurement of the flux function */
	double alpha[3];
	double omega[3]; /* weights for stencils */
	double q[3]; /* ENO approximation of flux */
	double sum;

	if (max_n < n)
        {
            max_n = n;
            if (f != NULL)
            {
                FT_FreeThese(3,f,fp,fm);
            }
            FT_VectorMemoryAlloc((POINTER*)&f,extend_size,sizeof(double));
            FT_VectorMemoryAlloc((POINTER*)&fp,n+1,sizeof(double));
            FT_VectorMemoryAlloc((POINTER*)&fm,n+1,sizeof(double));
            FT_VectorMemoryAlloc((POINTER*)&fm_burg,n+1,sizeof(double));
        }

	/* Find the maximum of fabs(df(u))*/
	max_df = 0;
	for (i = 0; i < extend_size; ++i) 
	{
	    norm = dF(v[i],u[i]) > 0 ? dF(v[i],u[i]) : -dF(v[i],u[i]);
	    max_df = max_df > norm ? max_df : norm;
	}

	/* f[i] = 0.5 * (f_{i - nrad} + max_df * u[i])
	 * ensure that df[i] > 0*/
	for (i = 0; i < extend_size; ++i)
	    f[i] = 0.5 * (flux_func(v[i],u[i]) + max_df * u[i]);

	for (j = 0; j < n + 1; ++j)
	/* To approximate flux at x_{j+1/2}, we use stencil
	 S_k = {x_{j+k-2}, x_{j+k-1}, x_{j+k}} = {x[j+k-2+nrad],
	 x[j+k-1+nrad], x[j+k+nrad]} */
	{
		/* compute the weights omega[] */
	    is[0] = 13.0/12.0*(f[j] - 2.0*f[j+1] + f[j+2])*(f[j] - 2.0
		    *f[j+1] + f[j+2]) + 0.25*(f[j] - 4.0*f[j+1] + 3.0
		    *f[j+2])*(f[j] - 4.0*f[j+1] + 3.0 * f[j+2]);
	    is[1] = 13.0/12.0*(f[j+1] - 2.0*f[j+2] + f[j+3])
		    *(f[j+1] - 2.0*f[j+2] + f[j+3]) + 0.25*(f[j+1]
		    -f[j+3])*(f[j+1] - f[j+3]);
	    is[2] = 13.0/12.0*(f[j+2] - 2.0*f[j+3] + f[j+4])
		    *(f[j+2] - 2.0*f[j+3] + f[j+4]) + 0.25*(3.0*f[j+2] 
		    - 4.0*f[j+3] + f[j+4])*(3.0*f[j+2] - 4.0*f[j+3] + f[j+4]);

	    for (i = 0; i < nrad; ++i)
		alpha[i] = c[i]/pow(eps + is[i],p);

	    sum = alpha[0] + alpha[1] + alpha[2];
	    for (i = 0; i < nrad; ++i)
		omega[i] = alpha[i]/sum;

		/* compute ENO approximation of the flux */
	    for (i = 0; i < nrad; ++i) 
	    {
		q[i] = 0.0;
		for (k = 0; k < nrad; ++k)
		    q[i] += a[i][k] * f[j + i + k];
	    }

	    /* compute the linear combination of the r candidate
	     stencils to get a higher order
	     approximation of the flux */
	    fp[j] = 0.0;
	    for (i = 0; i < nrad; ++i)
		fp[j] += omega[i]*q[i];
	}

	/* f[i] = 0.5 * (f_{i - nrad} - max_df * u[i])
	 * ensure that df[i] < 0*/
	for (i = 0; i < extend_size; ++i)
	    f[i] = 0.5*(flux_func(v[i],u[i]) + max_df*u[i]);
	for (j = 0; j < n + 1; ++j)
	/* To approximate flux at x_{j+1/2}, we use stencil S_k =
	 {x_{j+1-k+2}, x_{j+1-k+1}, x_{j+1-k}} = {x[j+1-k+2+nrad],
	 x[j+1-k+1+nrad], x[j+1-k+nrad]} */

	{
	    /* compute the weights omega[] */
	    is[0] = 13.0/12.0*(f[j+5] - 2.0*f[j+4] + f[j+3])
		    *(f[j+5] - 2.0*f[j+4] + f[j+3]) + 0.25*(f[j+5]
		    - 4.0*f[j+4] + 3.0*f[j+3])*(f[j+5] - 4.0*f[j+4]
		    + 3.0*f[j+3]);
	    is[1] = 13.0/12.0*(f[j+4] - 2.0*f[j+3] + f[j+2])
		    *(f[j+4] - 2.0*f[j+3] + f[j+2]) + 0.25*(f[j+4]
		    - f[j+2])*(f[j+4] - f[j+2]);
	    is[2] = 13.0/12.0*(f[j+3] - 2.0*f[j+2] + f[j+1])
		    *(f[j+3] - 2.0*f[j+2] + f[j+1]) + 0.25*(3.0*f[j+3] 
		    - 4.0*f[j+2] + f[j+1])*(3.0*f[j+3] - 4.0*f[j+2] + f[j+1]);

	    for (i = 0; i < nrad; ++i)
		alpha[i] = c[i]/pow(eps + is[i], p);

	    sum = alpha[0] + alpha[1] + alpha[2];
	    for (i = 0; i < nrad; ++i)
		omega[i] = alpha[i]/sum;

		/* compute ENO approximation of the flux */
	    for (i = 0; i < nrad; ++i) 
	    {
		q[i] = 0.0;
		for (k = 0; k < nrad; ++k)
		    q[i] += a[i][k]*f[j+5-i-k];
	    }

		/* compute the linear combination of the r candidate stencils
	     to get a higher order approximation of the flux */
	    fm[j] = 0.0;
	    for (i = 0; i < nrad; ++i)
		fm[j] += omega[i]*q[i];
	}
	
/*only for burger's equations*/
	for (i = 0; i < extend_size; ++i)
	    f[i] = 0.5*(flux_func(v[i],u[i]) - max_df*u[i]);
	for (j = 0; j < n + 1; ++j)
	/* To approximate flux at x_{j+1/2}, we use stencil S_k =
	 {x_{j+1-k+2}, x_{j+1-k+1}, x_{j+1-k}} = {x[j+1-k+2+nrad],
	 x[j+1-k+1+nrad], x[j+1-k+nrad]} */

	{
	    /* compute the weights omega[] */
	    is[0] = 13.0/12.0*(f[j+5] - 2.0*f[j+4] + f[j+3])
		    *(f[j+5] - 2.0*f[j+4] + f[j+3]) + 0.25*(f[j+5]
		    - 4.0*f[j+4] + 3.0*f[j+3])*(f[j+5] - 4.0*f[j+4]
		    + 3.0*f[j+3]);
	    is[1] = 13.0/12.0*(f[j+4] - 2.0*f[j+3] + f[j+2])
		    *(f[j+4] - 2.0*f[j+3] + f[j+2]) + 0.25*(f[j+4]
		    - f[j+2])*(f[j+4] - f[j+2]);
	    is[2] = 13.0/12.0*(f[j+3] - 2.0*f[j+2] + f[j+1])
		    *(f[j+3] - 2.0*f[j+2] + f[j+1]) + 0.25*(3.0*f[j+3] 
		    - 4.0*f[j+2] + f[j+1])*(3.0*f[j+3] - 4.0*f[j+2] + f[j+1]);

	    for (i = 0; i < nrad; ++i)
		alpha[i] = c[i]/pow(eps + is[i], p);

	    sum = alpha[0] + alpha[1] + alpha[2];
	    for (i = 0; i < nrad; ++i)
		omega[i] = alpha[i]/sum;

		/* compute ENO approximation of the flux */
	    for (i = 0; i < nrad; ++i) 
	    {
		q[i] = 0.0;
		for (k = 0; k < nrad; ++k)
		    q[i] += a[i][k]*f[j+5-i-k];
	    }

		/* compute the linear combination of the r candidate stencils
		 to get a higher order approximation of the flux */
	    fm_burg[j] = 0.0;
	    for (i = 0; i < nrad; ++i)
		fm_burg[j] += omega[i]*q[i];
	}
/*end fm_burg*/

	/*upwinding strategy*/
	for (j = 0; j < n; ++j)
	{
	    double fL, fR;
	    aa = 0.5*(v[j+nrad-1]+v[j+nrad]);
	    if (aa >= 0)
		fL = fp[j];
	    else
		fL = fm[j];

	    aa = 0.5*(v[j+nrad]+v[j+nrad+1]);
	    if (aa >= 0)
		fR = fp[j+1];
	    else
		fR = fm[j+1];
	    if (u[j] == v[j]) /*using different strategy for burger and linear*/
		flux[j] = lambda*(fp[j+1]+fm_burg[j+1] - fp[j]-fm_burg[j]);
	    else
	        flux[j] = lambda*(fR - fL);
	
	}

}

void weno5_flux(
	SWEEP *vst,
	FSWEEP *vflux,
	double lambda,
	int n,
	int dir,
	int dim)
{
	double **vel = vst->vel;
	double **vel_flux = vflux->vel_flux;
	int i, j;
	int nrad = 3;
	double a;

	if (debugging("weno"))
	{
	    (void) printf("Entering weno5_flux()\n");
	    (void) printf("dir = %d  n = %d\n",dir,n);
	}
	for (i = 0; i < dim; ++i)
	{
	    if (i == dir)
	    {
		dF = burger_dF;
		flux_func = burger_flux;
	    }
	    else
	    {
		dF = linear_dF;
		flux_func = linear_flux;
	    }
	    Weno5_Get_Flux(vel[i],vel[dir],vel_flux[i],lambda,n);
	    if (i != dir)
	    for (j = 0; j < n; j++)
	    {
	            vel_flux[i][j] *= vel[dir][j+nrad];
	    }
	}
}	/* end weno5_flux */


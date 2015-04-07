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

static int CIM1(CIM_COEF *CC,CIM_STRUCT CS,double ep0,double ep1,int dim);
static int CIM2(CIM_COEF *CC,CIM_STRUCT CS,double ep0,double ep1,int dim);
static boolean Linear_Solver(MTX,double*,double*);
static int tangential_direction(int i,double *T,double *N,int dim);
static boolean use_neumann;

CIM_ELLIPTIC_SOLVER::CIM_ELLIPTIC_SOLVER(Front &front):front(&front)
{
}

void CIM_ELLIPTIC_SOLVER::set_solver_domain(void)
{
	static boolean first = YES;
	top_grid = &topological_grid(front->grid_intfc);
        struct Table *T = table_of_interface(front->grid_intfc);
        int *lbuf = front->rect_grid->lbuf;
        int *ubuf = front->rect_grid->ubuf;

	dim = Dimension(front->interf);
        top_comp = T->components;
        top_gmax = top_grid->gmax;
	top_h = top_grid->h;
	top_L = top_grid->L;
	if (first)
	{
	    first = NO;
	    switch(dim)
	    {
	    case 1:
		imin = (lbuf[0] == 0) ? 1 : lbuf[0];
		imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
		break;
	    case 2:
		imin = (lbuf[0] == 0) ? 1 : lbuf[0];
	    	jmin = (lbuf[1] == 0) ? 1 : lbuf[1];
		imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
	    	jmax = (ubuf[1] == 0) ? top_gmax[1] - 1 : top_gmax[1] - ubuf[1];
		break;
	    case 3:
		imin = (lbuf[0] == 0) ? 1 : lbuf[0];
	    	jmin = (lbuf[1] == 0) ? 1 : lbuf[1];
		kmin = (lbuf[2] == 0) ? 1 : lbuf[2];
		imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
	    	jmax = (ubuf[1] == 0) ? top_gmax[1] - 1 : top_gmax[1] - ubuf[1];
		kmax = (ubuf[2] == 0) ? top_gmax[2] - 1 : top_gmax[2] - ubuf[2];
		break;
	    }
	}
}	/* end set_solver_domain */

void CIM_ELLIPTIC_SOLVER::solve(double *soln)
{
	switch (dim)
	{
	case 1:
	    return solve1d(soln);
	case 2:
	    return solve2d(soln);
	case 3:
	    return solve3d(soln);
	}
}	/* end solve */

void CIM_ELLIPTIC_SOLVER::solve1d(double *soln)
{
}	/* end solve1d */

void CIM_ELLIPTIC_SOLVER::solve2d(double *soln)
{
	int i,j,k,index,NzMax;

	NEWCIM = 1;
	NzMax = 12;

        u.N[0]  = imax - imin + 1;
        u.N[1]  = jmax - jmin + 1;
	u.NA = size;

	FT_VectorMemoryAlloc((POINTER*)&Order,size,sizeof(char));
	FT_VectorMemoryAlloc((POINTER*)&b,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&u.u,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&V,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&D,size,sizeof(char));
        for (i = 0; i < 2*dim; ++i) 
	{
	    FT_VectorMemoryAlloc((POINTER*)&NB[i],size,sizeof(int));
	    FT_VectorMemoryAlloc((POINTER*)&S[i],size,sizeof(char));
    	}


    	// Generate Domain indicator D, Neighborhood NB
	Generate_Domain_and_NB();
	for (i = imin; i <= imax; ++i)
	for (j = jmin; j <= jmax; ++j)
	{
	    k = ij_to_I[i][j];
	    index = d_index2d(i,j,top_gmax);
            b[k] = source[index];    	
        }
        Set_Intersection();

        Check_Type();

        A.N = size;
        A.K = NzMax*size;
	FT_VectorMemoryAlloc((POINTER*)&A.i,A.K,sizeof(int));
	FT_VectorMemoryAlloc((POINTER*)&A.j,A.K,sizeof(int));
	FT_VectorMemoryAlloc((POINTER*)&A.a,A.K,sizeof(double));
    	A.K = HCIM_Matrix_Generation();

        if (Linear_Solver(A,u.u,b) == YES) 
	{
	    double u_max,u_min;
	    u_max = -HUGE;
	    u_min = HUGE;
	    for (i = imin; i <= imax; ++i)
	    for (j = jmin; j <= jmax; ++j)
	    {
		index = d_index2d(i,j,top_gmax);
		k = ij_to_I[i][j];
		soln[index] = u.u[k];
	    }
        } 
	else 
	{
            (void) printf("Linear solver failed\n");
	    clean_up(ERROR);
        }
	FT_ParallelExchGridArrayBuffer(soln,front,NULL);
	if (solve_front_state)
	    cimSolveFrontState();
        
        FT_FreeThese(3,A.a,A.j,A.i);

        for (i = 0; i < 2*dim; ++i) 
	{
    	    FT_FreeThese(2,S[i],NB[i]);
    	}

        FT_FreeThese(5,b,Order,D,u.u,V);
}	/* end solve2d */

void CIM_ELLIPTIC_SOLVER::solve3d(double *soln)
{
	int ii,jj,kk,i,k,index,NzMax;

	NEWCIM = 1;
	NzMax = 12;

        u.N[0]  = imax - imin + 1;
        u.N[1]  = jmax - jmin + 1;
        u.N[2]  = kmax - kmin + 1;
	u.NA = size;

	FT_VectorMemoryAlloc((POINTER*)&Order,size,sizeof(char));
	FT_VectorMemoryAlloc((POINTER*)&b,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&u.u,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&V,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&D,size,sizeof(char));
        for (i = 0; i < 2*dim; ++i) 
	{
	    FT_VectorMemoryAlloc((POINTER*)&NB[i],size,sizeof(int));
	    FT_VectorMemoryAlloc((POINTER*)&S[i],size,sizeof(char));
    	}


    	// Generate Domain indicator D, Neighborhood NB
	Generate_Domain_and_NB();
	for (ii = imin; ii <= imax; ++ii)
	for (jj = jmin; jj <= jmax; ++jj)
	for (kk = kmin; kk <= kmax; ++kk)
	{
	    k = ijk_to_I[ii][jj][kk];
	    index = d_index3d(ii,jj,kk,top_gmax);
            b[k] = source[index];    	
        }
        Set_Intersection();

        Check_Type();

        A.N = size;
        A.K = NzMax*size;
	FT_VectorMemoryAlloc((POINTER*)&A.i,A.K,sizeof(int));
	FT_VectorMemoryAlloc((POINTER*)&A.j,A.K,sizeof(int));
	FT_VectorMemoryAlloc((POINTER*)&A.a,A.K,sizeof(double));
    	A.K = HCIM_Matrix_Generation();

        if (Linear_Solver(A,u.u,b) == YES) 
	{
	    double u_max,u_min;
	    u_max = -HUGE;
	    u_min = HUGE;
	    for (ii = imin; ii <= imax; ++ii)
	    for (jj = jmin; jj <= jmax; ++jj)
	    for (kk = kmin; kk <= kmax; ++kk)
	    {
		index = d_index3d(ii,jj,kk,top_gmax);
		k = ijk_to_I[ii][jj][kk];
		soln[index] = u.u[k];
	    }
        } 
	else 
	{
            (void) printf("Linear solver failed\n");
	    clean_up(ERROR);
        }
	FT_ParallelExchGridArrayBuffer(soln,front,NULL);
        
        FT_FreeThese(3,A.a,A.j,A.i);

        for (i = 0; i < 2*dim; ++i) 
	{
    	    FT_FreeThese(2,S[i],NB[i]);
    	}

        FT_FreeThese(5,b,Order,D,u.u,V);
}	/* end solve3d */


int CIM_ELLIPTIC_SOLVER::Generate_Domain_and_NB()
{
    	int ii,jj,kk,i,j,k,index,SD[3],ic[MAXD],ic2[MAXD];

    	SD[0] = SD[1] = SD[2] = 0;
	if (dim == 2)
	{
	    for (ii = imin; ii <= imax; ++ii)
	    for (jj = jmin; jj <= jmax; ++jj)
	    {
	    	ic[0] = ii;
	    	ic[1] = jj;
	    	k = ij_to_I[ii][jj];
	    	index = d_index(ic,top_gmax,dim);
	    	D[k] = (top_comp[index] == pos_comp) ? 1 : -1;

            	SD[D[k]+1]++;
            	for (i = 0; i < dim; ++i) 
	    	{
            	    for (j = 0; j < dim; ++j) 
		    	ic2[j] = ic[j];
            	    ic2[i] -= 1;
            	    NB[2*i][k] = ij_to_I[ic2[0]][ic2[1]];
            	    ic2[i] += 2;
            	    NB[2*i+1][k] = ij_to_I[ic2[0]][ic2[1]];
            	}
    	    }
	}
	else if (dim == 3)
	{
	    for (ii = imin; ii <= imax; ++ii)
	    for (jj = jmin; jj <= jmax; ++jj)
	    for (kk = kmin; kk <= kmax; ++kk)
	    {
	    	ic[0] = ii;
	    	ic[1] = jj;
	    	ic[2] = kk;
	    	k = ijk_to_I[ii][jj][kk];
	    	index = d_index(ic,top_gmax,dim);
	    	D[k] = (top_comp[index] == pos_comp) ? 1 : -1;

            	SD[D[k]+1]++;
            	for (i = 0; i < dim; ++i) 
	    	{
            	    for (j = 0; j < dim; ++j) 
		    	ic2[j] = ic[j];
            	    ic2[i] -= 1;
            	    NB[2*i][k] = ijk_to_I[ic2[0]][ic2[1]][ic2[2]];
            	    ic2[i] += 2;
            	    NB[2*i+1][k] = ijk_to_I[ic2[0]][ic2[1]][ic2[2]];
            	}
    	    }
	}
    	return YES;
}	/* end Generate_Domain_and_NB */

void CIM_ELLIPTIC_SOLVER::Set_Intersection()
{
    
    	int ii,jj,kk,k,l;

	if (dim == 2)
	{
	    for (ii = imin; ii <= imax; ++ii)
	    for (jj = jmin; jj <= jmax; ++jj)
	    {
	    	k = ij_to_I[ii][jj];
	    	for (l = 0; l < 2*dim; ++l) 
	    	{            
		    if (NB[l][k] >= 0) 
		    {
		    	S[l][k] = (D[k] == D[NB[l][k]]) ? 0 : 1;
		    } 
		    else 
		    {
		    	S[l][k] = 0;
		    }
	    	}
    	    }
	}
	else if (dim == 3)
	{
	    for (ii = imin; ii <= imax; ++ii)
	    for (jj = jmin; jj <= jmax; ++jj)
	    for (kk = kmin; kk <= kmax; ++kk)
	    {
	    	k = ijk_to_I[ii][jj][kk];
	    	for (l = 0; l < 2*dim; ++l) 
	    	{            
		    if (NB[l][k] >= 0) 
		    {
		    	S[l][k] = (D[k] == D[NB[l][k]]) ? 0 : 1;
		    } 
		    else 
		    {
		    	S[l][k] = 0;
		    }
	    	}
    	    }
	}
}	/* end Set_Intersection */

int CIM_ELLIPTIC_SOLVER::Check_Type_k(int k)
{
    	int Max, Min, SD[MAXD];
    	int i,j,s,decide,find;

    	// Check what we can do at a point P
    	//  0: max = 0: interior points
    	// 20: max = 1: CIM2 applicable
    	// 21: max = 1: CIM2 is not applicable at d-b-d approach (Solved by GVM)
    	//     CIM1 in CIM2007
    	// 22: max = 1: CIM2 is not applicable at cross approach in original CIM
    	//     CIM1 in CIM2007 (maybe 23, 24, 25 in CIM2011)
    	// 23: max = 1: CIM2 is applicable at cross approach in new CIM
    	//     CIM1 in CIM2007, Solved in CIM2011
    	// 24: max = 1: CIM2 is applicable at cross approach with principal 
	//     derivatives
    	//     CIM1 in CIM2007, Solved in CIM2011 but not coded
    	// 25: max = 1: CIM2 is applicable at cross approach with ghost values
    	//     CIM1 in CIM2007, Solved in CIM2011 but not coded    
    	// 30: max = 2, min = 1, GVM applicable
    	//     CIM1 in CIM2007, Solved in CIM2011
    	// 31: max = 2, min = 0, GVM applicable
    	//     CIM1 in CIM2007, Solved in CIM2011
    	// 10: max = min = 2: surround by interface (Solved by GVM with Local 
	//     Poisson Solver)
    	//     CIM1 in CIM2007, Solved in CIM2011 but not coded
    	// 11: max = 2, min = 0, GVM is not applicable (Solved by refining mesh)
    	//     CIM1 in CIM2007, NOT Solved in CIM2011

    	Ord[0] = Ord[1] = Ord[2] = Ord[3] = Ord[4] = 0;

	{
            Max = 0; Min = 2; s = 0; decide = -1;
            for (i = 0; i < dim; ++i) 
	    {
            	SD[i] = S[2*i][k]+S[2*i+1][k];
            	if (SD[i] > Max) Max = SD[i];
            	if (SD[i] < Min) Min = SD[i];
            }
            if(Max == 0) 
	    {
            	// Interior point
            	decide = 0;
            } 
	    else if (Max == 2 && Min == 2) 
	    {
            	// Type I exceptional
            	decide = 10;
            } 
	    else if (Max == 2 && Min == 1) 
	    {
            	decide = 30;
            } 
	    else if (Max == 2 && Min == 0) 
	    {
            	for (i = 0; i < dim; ++i) 
		{
                    if (SD[i] == 0) 
		    {
                    	if (NB[2*i][k] < 0 || NB[2*i+1][k] < 0) 
			{
                            printf("CIM2 at boundary! Refine mesh\n");
                    	} 
			else 
			{
                            if (S[2*i][NB[2*i][k]] == 0 || 
				S[2*i+1][NB[2*i+1][k]] == 0) 
			    {
                            	decide = 11;
                            } 
			    else 
			    {
                            	decide = 31;
                            }
                    	}
                    }
            	}
            } 
	    else if (Max == 1) 
	    {
            	// Check CIM2 applicability
            	for (i = 0; i < 2*dim; ++i) 
		{
                    if (NB[i][k] < 0) 
		    {
                    	if (S[i][k] == 1) 
			{
                            printf("CIM2 at boundary! Refine mesh!\n");
                    	} 
			else 
			{
                            printf("CIM2 at boundary! Be careful\n");
                    	}
                    } 
		    else 
		    {
                    	if (S[i][k]+S[i][NB[i][k]] == 2) 
			{
                            decide = 21;
                            break;
                    	}
                    }
            	}
            	if (decide == -1) 
		{
                    if (NEWCIM == 0) 
		    {
                    	// CIM 2007
                    	find = 1;
                    	for (i = 0; i < 2*dim; ++i) 
			{
                            if (S[i][k] == 1) 
			    {
                            	for (j = 0; j < 2*dim; ++j) 
				{
                                    if ((i/2) != (j/2)) 
				    {
                                    	if (S[j][k] == 0 && 
					    NB[i+1-2*(i%2)][k] >= 0)
					    if(S[j][NB[i+1-2*(i%2)][k]] != 0) 
						find = 0;
                                    }
                            	}
                            }
                    	}
                    	if (find == 0) 
			{
			    // Type II exceptional
                            decide = 22;
                    	} 
			else 
			{
                            decide = 20;
                    	}
                    } 
		    else 
		    {
                    	// CIM 2011
                    	for (i = 0; i < dim; ++i) 
			{
                            for (j = i+1; j < dim; ++j) 
			    {
                            	find = 0;
                            	if (S[2*i][k] == 0 && S[2*j][k] == 0 && 
				    S[2*i][NB[2*j][k]] == 0 && 
				    S[2*j][NB[2*i][k]] == 0) 
				    find = 1;
                            	if (S[2*i+1][k] == 0 && S[2*j][k] == 0 && 
				    S[2*i+1][NB[2*j][k]] == 0 && 
				    S[2*j][NB[2*i+1][k]] == 0) 
				    find = 1;
                            	if (S[2*i][k] == 0 && S[2*j+1][k] == 0 && 
				    S[2*i][NB[2*j+1][k]] == 0 && 
				    S[2*j+1][NB[2*i][k]] == 0) 
				    find = 1;
                            	if (S[2*i+1][k] == 0 && S[2*j+1][k] == 0 && 
				    S[2*i+1][NB[2*j+1][k]] == 0 && 
				    S[2*j+1][NB[2*i+1][k]] == 0) 
				    find = 1;
                            	if (find == 0) 
				{                                
                                    if (S[2*i][k] == 0 && S[2*j][k] == 0 && 
					S[2*j+1][NB[2*i][k]] == 0) 
					find = 1;
                                    if (S[2*i][k] == 0 && S[2*j][k] == 0 && 
					S[2*i+1][NB[2*j][k]] == 0) 
					find = 1;
                                    if (S[2*i+1][k] == 0 && S[2*j][k] == 0 && 
					S[2*j+1][NB[2*i+1][k]] == 0) 
					find = 1;
                                    if (S[2*i+1][k] == 0 && S[2*j][k] == 0 && 
					S[2*i][NB[2*j][k]] == 0) 
					find = 1;
                                    if (S[2*i][k] == 0 && S[2*j+1][k] == 0 && 
					S[2*j][NB[2*i][k]] == 0) 
					find = 1;
                                    if (S[2*i][k] == 0 && S[2*j+1][k] == 0 && 
					S[2*i+1][NB[2*j+1][k]] == 0) 
					find = 1;
                                    if (S[2*i+1][k] == 0 && S[2*j+1][k] == 0 &&
					S[2*j][NB[2*i+1][k]] == 0) 
					find = 1;
                                    if (S[2*i+1][k] == 0 && S[2*j+1][k] == 0 &&
					S[2*i][NB[2*j+1][k]] == 0) 
					find = 1;                              
                                    if (find == 0) 
				    {
                                    	decide = 25;
                                    } 
				    else 
				    {
                                    	decide = 24;
                                    }
                            	} 
				else 
				{
                                    decide = 23;
                            	}
                            }
                    	}
                    }                
            	}
            }
            if (decide == -1) 
	    {
            	printf("NOT EXPECTED!! Check code!!\n"); 
            	return 0;
            }
            if (NEWCIM == 0) 
	    {
            	switch(decide) 
		{
                case 0:
                    Order[k] = 0;
                    break;
                case 20:
                    Order[k] = 2;
                    break;
                default:
                    Order[k] = 1;
                    break;
        	}
            } 
	    else 
	    {
            	switch(decide) 
		{
                case 0:
                    Order[k] = 0;
                    break;
                case 10:
                case 24:
                case 25:
                    Order[k] = 1;
                    break;
                case 11:
		    Order[k] = 1;
		    Ord[4]++;
		    break;
		case 21:
		    if (NEWCIM == 1) 
		    {
			Order[k] = 1;
		    } 
		    else 
		    {
			Order[k] = 1;
		    }
		    break;
                case 30:
                case 31:
                    Order[k] = 3;
                    break;
                default:
                    Order[k] = 2;
                    break;
            	}
            }
            Ord[(int)Order[k]]++;
	}
}	/* end Check_Type */

int CIM_ELLIPTIC_SOLVER::Check_Type()
{
    	int k,ii,jj,kk;

	switch (dim)
	{
	case 2:
	    for (ii = imin; ii <= imax; ++ii)
	    for (jj = jmin; jj <= jmax; ++jj)
	    {
	    	k = ij_to_I[ii][jj];
	    	Check_Type_k(k);
    	    }    
	    break;
	case 3:
	    for (ii = imin; ii <= imax; ++ii)
	    for (jj = jmin; jj <= jmax; ++jj)
	    for (kk = kmin; kk <= kmax; ++kk)
	    {
	    	k = ijk_to_I[ii][jj][kk];
	    	Check_Type_k(k);
    	    }    
	    break;
	}
    	return 1;
}	/* enc Check_Type */

int CIM_ELLIPTIC_SOLVER::Search_Ghost_Value(int k)
{
    	int i, j, k0, k1;
    	D[k] = -D[k];
    	for (i = 0; i < 2*dim; ++i) 
	{
            S[i][k] = 1-S[i][k];
            if (NB[i][k] >= 0)
            	S[i+1-2*(i%2)][NB[i][k]] = 1-S[i+1-2*(i%2)][NB[i][k]];
            else 
	    {
            	printf("Too close to boundary!!\n");
            }
    	}
    	for (i = 0; i < dim; ++i) 
	{
            if (S[2*i][k] + S[2*i+1][k] == 2) 
	    {
            	k0 = NB[2*i][k];
            	k1 = NB[2*i+1][k];
            	if (k0 >= 0 && Order[k0] != 7 && Order[k0] != 1 &&
               	    k1 >= 0 && Order[k1] != 7 && Order[k1] != 1) 
		{
                    printf("(1)Ghost value at %d->%d\n", k, k0);
                    Order[k0] = 7;
                    Search_Ghost_Value(k0);
                    printf("(1)Ghost value at %d->%d\n", k, k1);
                    Order[k1] = 7;
                    Search_Ghost_Value(k1);
            	} 
		else if (k0 < 0 || k1 < 0) 
		{
                    printf("Ghost value at boundary\n");
            	}
            }
    	}
    	for (i = 0; i < 2*dim; ++i) 
	{
            k0 = NB[i][k];
            for (j = 0; j < dim; ++j) 
	    {
            	if (S[2*j][k0]+S[2*j+1][k0] == 2 && Order[k0] != 7) 
		{
                    if (Order[k0] != 3) 	
			printf("(2)Ghost value at %d->%d\n", k, k0);
                    Order[k0] = 7;
                    Search_Ghost_Value(k0);
            	}
            }
    	}
    	return 1;
}

int CIM_ELLIPTIC_SOLVER::HCIM_Matrix_Generation()
{
	int ii,jj,kk,k,n;

	n = 0;
	switch (dim)
	{
	case 2:
	    for (ii = imin; ii <= imax; ++ii)
	    for (jj = jmin; jj <= jmax; ++jj)
	    {
	        k = ij_to_I[ii][jj];
		HCIM_Matrix_Generation_k(k,&n);
	    }
	    break;
	case 3:
	    for (ii = imin; ii <= imax; ++ii)
	    for (jj = jmin; jj <= jmax; ++jj)
	    for (kk = kmin; kk <= kmax; ++kk)
	    {
	        k = ijk_to_I[ii][jj][kk];
		HCIM_Matrix_Generation_k(k,&n);
	    }
	    break;
	}
    	return n;
}	/* end HCIM_Matrix_Generation */

int CIM_ELLIPTIC_SOLVER::HCIM_Matrix_Generation_k(
	int k,
	int *n)
{
    	int i,j,index,s,N,Ri;
    	int *row, *col;
    	int icoords[MAXD], ic[MAXD];
    	double P[MAXD], Q[MAXD], R[MAXD];
    	double epi, epo;
    	double J_u[2*MAXD], J_eps_un[2*MAXD], J_ut[2*MAXD];
	double *h = top_h;
    	CIM_STRUCT CS;
    	CIM_COEF CC;
    	double v, *val;
	HYPER_SURF *hs;
	double crx_coords[MAXD];
	COMPONENT comp;
	GRID_DIRECTION dir[3][2] = {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
	int status;
	int icrds[MAXD];
	int idir,nb;
	double kp;
	int index_nb;
	POINTER state;
    
    	row = A.i;
    	col = A.j;
    	val = A.a;
    
	Ri = 0;
	use_neumann = YES;

	{
	    for (i = 0; i < dim; ++i) 
	    {
		if (dim == 2)
		    icoords[i] = I_to_ij[k-ilower][i];
		else if (dim == 3)
		    icoords[i] = I_to_ijk[k-ilower][i];
            	P[i] = cell_edge(icoords[i],i,top_grid);
            }
	    index = d_index(icoords,top_gmax,dim);
	    comp = top_comp[index];
            epi = (D[k] ==  1) ? diff_coeff[0] : diff_coeff[1];
            epo = (D[k] == -1) ? diff_coeff[0] : diff_coeff[1];
            if (Order[k] == 0) 
	    {
            	// Standard Finite Difference 
		v = 0;
		boolean neumann_nb = NO;
		for (i = 0; i < 2*dim; ++i) 
		{
		    idir = i/2;
		    nb = i%2;
		    status = (*findStateAtCrossing)(front,icoords,
				dir[idir][nb],comp,&state,&hs,crx_coords);
		    if (NB[i][k] >= 0 && NB[2*idir+((i+1)%2)][k]>=0) 
		    {
			row[*n] = k;
			col[*n] = NB[i][k];
			val[*n] = -1.0/sqr(h[idir])*epi;
			v += 1.0/sqr(h[idir])*epi;
			(*n)++;
		    } 
		    else if (status == NEUMANN_PDE_BOUNDARY)
			continue;
		    else if(NB[i][k] < 0 && (NB[2*idir+((i+1)%2)][k]>=0))
		    {
			b[k] += 8.0/3/sqr(h[idir])*epi*getStateVar(state);
		    	v += 8.0/3/sqr(h[idir])*epi;
			use_neumann = NO;
		    } 
		    else
		    {
			row[*n] = k;
			col[*n] = NB[i][k];
			val[*n] = -4.0/3/sqr(h[idir])*epi;
		    	v += 4.0/3/sqr(h[idir])*epi;
			(*n)++;
		    }
		}
		row[*n] = col[*n] = k;
		val[*n] = v;
		(*n)++;
	    } 
	    else if (Order[k] == 2) 
	    {
            	// CIM2
		for (i = 0; i < dim; ++i) 
		{
            	    CS.s[i] = S[2*i+1][k]-S[2*i][k];
            	    CS.a[i] = 1.0;
            	    for (j = 0; j < dim; ++j) 
		    {
			CS.n[i][j] = (i == j) ? 1 : 0;
		    }
		}
		for (i = 0; i < dim; ++i) 
		{
                    J_u[i] = J_eps_un[i] = J_ut[i] = 0;
            	    if (CS.s[i] != 0) 
		    {
			for (j = 0; j < dim; ++j) 
			    icrds[j] = icoords[j];
			if (CS.s[i] == -1)
			{
			    status = FT_NormalAtGridCrossing(front,icrds,
					dir[i][0],comp,R,&hs,crx_coords);
			    if (comp == negative_component(hs))
			    {
				for (j = 0; j < dim; ++j) 
				    R[i] *= -1;
			    }
			}
			else
			{
			    status = FT_NormalAtGridCrossing(front,icrds,
					dir[i][1],comp,R,&hs,crx_coords);
			    if (comp == negative_component(hs))
			    {
				for (j = 0; j < dim; ++j) 
				    R[i] *= -1;
			    }
			}
			for (j = 0; j < dim; ++j) 
			    Q[j] = crx_coords[j];	
			CS.a[i] = (crx_coords[i] - P[i])/top_h[i]/CS.s[i];

			J_u[i] = (*solutionJump)(jparams,D[k],Q);
			J_eps_un[i] = (*gradJumpDotN)(jparams,D[k],R,Q);
			J_ut[i] = (*gradJumpDotT)(jparams,D[k],i,R,Q);
			for (j = 0; j < dim; ++j) 
			{
			    CS.n[i][j] = R[j];
			    CS.c[i][j] = CS.s[j];
			}
		    }
	        }
	        N = CIM2(&CC,CS,epi,epo,dim);
	        for (j = 0; j < N; ++j) 
		{
		    row[*n] = k;
		    for (i = 0; i < dim; ++i) 
		    {
			ic[i] = icoords[i] + CC.indx[i][j];
		    }
		    if (dim == 2)
		    	s = ij_to_I[ic[0]][ic[1]];
		    else if (dim == 3)
		    	s = ijk_to_I[ic[0]][ic[1]][ic[2]];
		    if (s >= 0) 
		    {
			col[*n] = s;
			v = 0;
			for (i = 0; i < dim; ++i) 
			{
			    v += CC.coef[i][j]/h[i]/h[i]*epi;
			}
			val[*n] = -v;
			(*n)++;
		    } 
		    else 
		    {
			// boundary treatment
			printf("CIM2 at boundary!!\n");
		    }
		}
		for (i = 0; i < dim; ++i) 
		{
                    for (j = 0; j < dim; ++j) 
		    {
                    	b[k] += (CC.J[i][3*j]*J_eps_un[j]+CC.J[i][3*j+1]*
				J_u[j]/h[i]+CC.J[i][3*j+2]*J_ut[j])
				/h[i]*epi;
                    }
            	}
            } 
	    else 
	    { 
		for (i = 0; i < 2*dim; ++i) 
		{
            	    CS.s[i] = S[i][k];
            	    CS.a[i] = 1.0;
            	    for (j = 0; j < dim; ++j) 
		    {
			CS.n[i][j] = (i == j) ? 1 : 0;
		    }
		}
		for (i = 0; i < 2*dim; ++i) 
		{
		    idir = i/2;
		    nb = (2*(i%2)-1);
                    J_u[i] = J_eps_un[i] = J_ut[i] = 0;
            	    if (CS.s[i] != 0) 
		    {
			for (j = 0; j < dim; ++j) 
			    icrds[j] = icoords[j];
			if (nb == -1)
			{
			    status = FT_NormalAtGridCrossing(front,icrds,
					dir[idir][0],comp,R,&hs,crx_coords);
			    if (comp == negative_component(hs))
			    {
				for (j = 0; j < dim; ++j) 
				    R[i] *= -1;
			    }
			}
			else
			{
			    status = FT_NormalAtGridCrossing(front,icrds,
					dir[idir][1],comp,R,&hs,crx_coords);
			    if (comp == negative_component(hs))
			    {
				for (j = 0; j < dim; ++j) 
				    R[i] *= -1;
			    }
			}
			for (j = 0; j < dim; ++j) 
			    Q[j] = crx_coords[j];	
			CS.a[i] = (crx_coords[i] - P[i])/top_h[idir]/CS.s[i];
					
			J_u[i] = (*solutionJump)(jparams,D[k],Q);
			J_eps_un[i] = (*gradJumpDotN)(jparams,D[k],R,Q);
			J_ut[i] = (*gradJumpDotT)(jparams,D[k],i/2,R,Q);
			for (j = 0; j < dim; ++j) 
			{
			    CS.n[i][j] = R[j];
			}
		    }
	        }
	        N = CIM1(&CC,CS,epi,epo,dim); // CIM1
	        for (j = 0; j < N; ++j) 
		{
		    row[*n] = k;
		    for (i = 0; i < dim; ++i) 
		    {
			ic[i] = icoords[i] + CC.indx[i][j];
		    }
		    v = 0;
		    for (i = 0; i < dim; ++i) 
		    {
			v += CC.coef[i][j];
		    }
		    if (dim == 2)
		    	s = ij_to_I[ic[0]][ic[1]];
		    else if (dim == 3)
		    	s = ijk_to_I[ic[0]][ic[1]][ic[2]];
		    if (s >= 0) 
		    {
			col[*n] = s;
			v = 0;
			for  (i = 0; i < dim; ++i) 
			{
			    v += CC.coef[i][j]/h[i]/h[i]*epi;
			}
			val[*n] = -v;
			(*n)++;
		    } 
		    else 
		    {
			// boundary treatment
			printf("CIM1 at boundary!!\n");
		    }
		}
		for (i = 0; i < dim; ++i) 
		{
                    for (j = 0; j < 2*dim; ++j) 
		    {
                    	b[k] += (CC.J[i][3*j]*J_eps_un[j]+CC.J[i][3*j+1]*
					J_u[j]/h[i]+CC.J[i][3*j+2]*
					J_ut[j])/h[i]*epi;
                    }
            	}
	    }
    	}
}	/* end HCIM_Matrix_Generation_k */

#define EPS    1e-13

static int Gauss_Elimination_of_CIM(double **M,double *x,int N);

static int CIM1(
	CIM_COEF *CC, 
	CIM_STRUCT CS, 
	double ep0, 
	double ep1,
	int dim)
{
    	// CC.indx[i][j]
    	// CC.coef[i][j]
    	// j-th coefficient, i: u_{x_i x_i}
    	// total:  2*dim+1 coefficients, 
    	//         6*dim jump conditions
    
	int i, j, s;
	double ep[2*MAXD], coef[2*MAXD][2*MAXD+1], J[2*MAXD][6*MAXD];
	double **M, **MI;

	FT_VectorMemoryAlloc((POINTER*)&M,2*dim,sizeof(double*));
	FT_VectorMemoryAlloc((POINTER*)&M[0],2*dim*2*dim,sizeof(double));
    	for (i = 1; i < 2*dim; ++i) 
	    M[i] = M[i-1]+2*dim;
	FT_VectorMemoryAlloc((POINTER*)&MI,2*dim,sizeof(double*));
	FT_VectorMemoryAlloc((POINTER*)&MI[0],2*dim*2*dim,sizeof(double));
    	for (i = 1; i < 2*dim; ++i) 
	    MI[i] = MI[i-1]+2*dim;
    	for (i = 0; i < 2*dim; ++i) 
	{
            for (j = 0; j < 2*dim; ++j) 
	    {
            	M[i][j] = MI[i][j] = 0;
            }
            M[i][i] = MI[i][i] = 1;
    	}

	for (i = 0; i < 2*dim; ++i) 
	{
	    for (j = 0; j < 6*dim; ++j) 
	    {
		J[i][j] = 0;
	    }
	    for (j = 0; j < 2*dim+1; j++) 
	    {
		coef[i][j] = 0;
	    }
	}

	for (i = 0; i < dim; ++i) 
	{
	    for (j = 0; j < 6*dim; ++j) 
	    {
		CC->J[i][j] = 0;
	    }
	    for (j = 0; j < 2*dim+1; j++) 
	    {
		CC->indx[i][j] = 0;
		CC->coef[i][j] = 0;
	    }
	}

	for (i = 0; i < 2*dim; ++i) 
	{
	    if (CS.s[i] == 0) 
	    {
		coef[i][0]   = 1-2*(i%2);
		coef[i][i+1] = -coef[i][0];
		ep[i] = ep0;
	    } 
	    else 
	    {
		ep[i] = (1-CS.a[i])*ep0+CS.a[i]*ep1;
		coef[i][0]   = (1-2*(i%2))*ep1/ep[i];
		coef[i][i+1] = -coef[i][0];
		M[i][i] += (1-CS.a[i])/ep[i]*(ep1-ep0)*
				(1-CS.n[i][i/2]*CS.n[i][i/2]);
		J[i][3*i]   = -(1-CS.a[i])/ep[i];
		J[i][3*i+1] =  (1-2*(i%2))*ep1/ep[i];
	    }
	}
	for (i = 0; i < 2*dim; ++i) 
	{
	    if (CS.s[i] == 1) 
	    {
		for (j = 0; j < dim; ++j) 
		{
		    if (j != i/2) 
		    {
			if (CS.s[2*j + (i%2)] == 0) 
			{
			    M[i][2*j + (i%2)] -= (1-CS.a[i])/ep[i]*(ep1-
					ep0)*(CS.n[i][i/2]*CS.n[i][j]);
			} 
			else if (CS.s[2*j+1-(i%2)] == 0) 
			{
			    M[i][2*j+1-(i%2)] -= (1-CS.a[i])/ep[i]*(ep1-
					ep0)*(CS.n[i][i/2]*CS.n[i][j]);
			} 
			else if (fabs(M[2*j + (i%2)][2*j + (i%2)]) >= 
				 fabs( M[2*j+1-(i%2)][2*j+1-(i%2)])) 
			{
			    M[i][2*j + (i%2)] -= (1-CS.a[i])/ep[i]*(ep1-
					ep0)*(CS.n[i][i/2]*CS.n[i][j]);
			} 
			else 
			{
			    M[i][2*j+1-(i%2)] -= (1-CS.a[i])/ep[i]*(ep1-
					ep0)*(CS.n[i][i/2]*CS.n[i][j]);
			}
		    }
		}
	    }
	}
	
	for(i=0;i<2*dim;++i) {
		Gauss_Elimination_of_CIM(M, MI[i], 2*dim);
	}
	for (j = 0; j < 2*dim+1; ++j) 
	{
	    for (s = 0; s < 2*dim; ++s) 
	    {
		M[0][s] = 0;
		for (i = 0; i < 2*dim; ++i) 
		{
		    M[0][s] += MI[i][s]*coef[i][j];
		}
	    }
	    for (i = 0; i < 2*dim; ++i) 
	    {
		coef[i][j] = M[0][i];
	    }
	}
	// put index
	for (i = 0; i < 2*dim; ++i) 
	{
            CC->indx[i/2][i+1] = 2*(i % 2) - 1;
    	}
    
    	// put coefficients
	for (i = 0; i < dim; ++i) 
	{
            for (j = 0; j < 2*dim+1; ++j) 
	    {
            	CC->coef[i][j] = coef[2*i+1][j]-coef[2*i][j];
            }
    	}
	for (i = 0; i < 2*dim; ++i) 
	{
	    J[i][3*i+2]  = J[i][3*i]*ep1*sqrt(1-CS.n[i][i/2]*CS.n[i][i/2]);
	    J[i][3*i]   *= CS.n[i][i/2];
	}
	// here
	for (j = 0; j < 6*dim; ++j) 
	{
	    for (s = 0; s < 2*dim; ++s) 
	    {
		M[0][s] = 0;
		for (i = 0; i < 2*dim; ++i) 
		{
		    M[0][s] += MI[i][s]*J[i][j];
		}
	    }
	    for (i = 0; i < 2*dim; ++i) 
	    {
		J[i][j] = M[0][i];
	    }
	}
	for (i = 0; i < dim; ++i) 
	{
	    for (j = 0; j < 6*dim; ++j) 
	    {
		CC->J[i][j] = J[2*i+1][j] - J[2*i][j];
	    }
	}
	FT_FreeThese(4,M[0],M,MI[0],MI);
    	return 2*dim+1;
}

static int CIM2(
	CIM_COEF *CC, 
	CIM_STRUCT CS, 
	double ep0, 
	double ep1,
	int dim)
{
    	// CC.indx[i][j]
    	// CC.coef[i][j]
    	// j-th coefficient, i: u_{x_i x_i}
    	// total:  2d: 8 coefficients, 3d: 11-13 coefficients
    	//         3*dim jump conditions
    
	int i, j, k, s, p, N, newid, indx[MAXD][20*MAXD]; 
	double a, b, ep, r0, r1, L, coef[MAXD][20*MAXD], J[MAXD][3*MAXD]; 
	double **M, **MI, t[MAXD][MAXD];

    	// Initial setup	
	for (i = 0; i < dim; ++i) 
	{
	    for (j = 0; j < 3*dim; ++j) 
	    {
		J[i][j] = 0;
	    }
	    for (j = 0; j < 20*dim; j++) 
	    {
		indx[i][j] = 0;
		coef[i][j] = 0;
	    }
	    tangential_direction(i,t[i],CS.n[i],dim);
	}
	k = 0;
	// Dimension-by-Dimension approach
	for (i = 0; i < dim; ++i) 
	{
	    s = CS.s[i];
	    if (s == 0) 
	    {
		indx[i][k] = 0; indx[i][k+1] = -1; indx[i][k+2] =  1;
		coef[i][k] =-2; coef[i][k+1] =  1; coef[i][k+2] =  1;
		J[i][3*i+1] = 0;
		J[i][3*i] = 0;
		k += 3;
	    } 
	    else 
	    {
		a = CS.a[i];
		b = 1-a;
		ep = (0.5+a)*(b+b*b)*ep0+(0.5+b)*(a+a*a)*ep1;
		if (fabs(ep) < CIM_TOL && debugging("cim_solver")) 
		{
                    printf("determinant in 1D too small:%e, %f, %f, %f\n",
					ep,a,ep0,ep1);
            	}
		r0 = ep0/ep; r1 = ep1/ep;
		indx[i][k]   = 0; 
		coef[i][k]   = -(b+b*b)*r0-(1+a)*(1+2*b)*r1;
		indx[i][k+1] = -s; 
		coef[i][k+1] = (b+b*b)*r0+a*(1+2*b)*r1;
		indx[i][k+2] = s; 
		coef[i][k+2] = (1+b)*(1+b)*r1;
		indx[i][k+3] = 2*s;
		coef[i][k+3] = -b*b*r1;	
		J[i][i*3+1] = -(1+2*b)*r1;
		J[i][i*3] = -s*(b+b*b)/ep;
		k += 4;
	    }
	}

	if (dim == 1) 
	{
	    for (j = 0; j < N; ++j) 
	    {
		CC->indx[0][j] = indx[0][j];
		CC->coef[0][j] = coef[0][j];
	    }
	    CC->J[0][0] = J[0][0];
	    CC->J[0][1] = J[0][1];
	    return k;
	}
	for (i = 0; i < dim; ++i) 
	{
	    s = CS.s[i];
	    if (s != 0) 
	    {
    		for (j = 0; j < dim; ++j) 
		{
    		    if(i == j) 
		    {
    			indx[i][k]   =  0;
    			indx[i][k+1] = -s;
    			coef[i][k] =  s*J[i][i*3]*(ep1-ep0)*t[i][i]*
						t[i][i];
    			coef[i][k+1] = -coef[i][k];
    			k += 2;
    		    } 
		    else 
		    {                                        
    			indx[i][k]   = 0;
    			indx[j][k]   = -tr(-CS.c[i][j]);
    			indx[i][k+1] = 0;
    			indx[j][k+1] = tr(CS.c[i][j]);
    			indx[i][k+2] = -CS.s[i];
    			indx[j][k+2] = tr(CS.c[i][j]);
    			indx[i][k+3] = -CS.s[i];
    			indx[j][k+3] = -tr(-CS.c[i][j]); 
    			coef[i][k]   = -(1+CS.a[i])*J[i][i*3]*(ep1-ep0)*
				(t[i][i])*(t[i][j])/(2-fabs(1.0*CS.c[i][j]));
    			coef[i][k+1] = -coef[i][k];
    			coef[i][k+2] = -(  CS.a[i])*J[i][i*3]*(ep1-ep0)*
				(t[i][i])*(t[i][j])/(2-fabs(1.0*CS.c[i][j]));
    			coef[i][k+3] = -coef[i][k+2];
    			k += 4;
    		    }
    		}
            }
	}
	FT_VectorMemoryAlloc((POINTER*)&M,dim,sizeof(double*));
	FT_VectorMemoryAlloc((POINTER*)&M[0],dim*dim,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&MI,dim,sizeof(double*));
	FT_VectorMemoryAlloc((POINTER*)&MI[0],dim*dim,sizeof(double));

	for (i = 1; i < dim; ++i) 
	{
	    M[i]  = M[i-1]+dim;
	    MI[i] = MI[i-1]+dim;
	}
	for (i = 0; i < dim; ++i) 
	{
	    for (j = 0; j < dim; ++j) 
	    {
		M[i][j] = 0;
		L = (i==j) ? (0.5+CS.a[i])*CS.s[i] : 0.5*CS.c[i][j];
		M[i][j] -= J[i][i*3]*(ep1-ep0)*t[i][i]*t[i][j]*L;
	    }
	    M[i][i] += 1;
	}

	for (i = 0; i < dim; ++i) 
	{
	    for (j = 0; j < dim; ++j) 
	    {
		MI[i][j] = (i==j) ? 1 : 0;
	    }
	    Gauss_Elimination_of_CIM(M, MI[i], dim);
	}

	for (j = 0; j < k; ++j) 
	{
	    for (s = 0; s < dim; ++s) 
	    {
		M[0][s] = 0;
		for (i = 0; i < dim; ++i) 
		{
		    M[0][s] += MI[i][s]*coef[i][j];
		}
	    }
	    for (i = 0; i < dim; ++i) 
	    {
		coef[i][j] = M[0][i];
	    }
	}
	N = 0;
	for (j = 0; j < k; ++j) 
	{
	    newid = 1;
	    for (s = 0; s < N; ++s) 
	    {
		p = 0;
		for (i = 0; i < dim; ++i) 
		{
		    if (CC->indx[i][s] != indx[i][j]) 
		    {
			p = 1;
			break;
		    }
		}
		if (p == 0) 
		{
		    newid = 0;
		    break;
		}
	    }
	    if (newid == 1) 
	    {
		for (i = 0; i < dim; ++i) 
		{
		    CC->indx[i][N] = indx[i][j];
		    CC->coef[i][N] = 0.0;
		}
		N += 1;
	    }
	}
    	if (N > 14) 
	{
            for (i = 0; i < dim; ++i) 
	    {
            	for (j = 0; j < dim; ++j) 
		{
                    printf("%d ", CS.c[i][j]);
            	}
            	printf("\n");
            }
    	}
	for (j = 0; j < k; ++j) 
	{
	    for (s = 0; s < N; ++s) 
	    {
		p = 0;
		for (i = 0; i < dim; ++i) 
		{
		    if (CC->indx[i][s] != indx[i][j]) 
		    {
			p = 1;
			break;
		    }
		}
		if (p == 0) 
		{
		    break;
		}
	    }
	    for (i = 0; i < dim; ++i) 
	    {
		CC->coef[i][s] += coef[i][j];
	    }
	}
	
	for (i = 0; i < dim; ++i) 
	{
	    J[i][i*3+2] = J[i][i*3]*ep1*t[i][i];
	    J[i][i*3]  *= CS.n[i][i];
	}
	for (j = 0; j < 3*dim; ++j) 
	{
	    for (s = 0; s < dim; ++s) 
	    {
		M[0][s] = 0;
		for (i = 0; i < dim; ++i) 
		{
		    M[0][s] += MI[i][s]*J[i][j];
		}
	    }
	    for (i = 0; i < dim; ++i) 
	    {
		CC->J[i][j] = M[0][i];
	    }
	}
	FT_FreeThese(4,M[0],M,MI[0],MI);
	return N;
}

static int Gauss_Elimination_of_CIM(
	double **M, 
	double *x, 
	int N)
{
    	int i, j, k, row_max, result;
    	double *Swap_A, a_max, v;
    	double **A, *A0;

	FT_VectorMemoryAlloc((POINTER*)&A,N,sizeof(double*));
	FT_VectorMemoryAlloc((POINTER*)&A0,N*N,sizeof(double));
    	A[0] = A0;
    	for (i = 1; i < N; ++i) 
	    A[i] = A[i-1] + N;   
    	for (i = 0; i < N; ++i) 
	{
            for (j = 0; j < N; ++j) 
	    {
            	A[i][j] = M[i][j];
            }
    	} 
    
    	for (i = 0; i < N-1; ++i) 
	{
            a_max = 0.0;
            row_max = i;
            for (j = i; j < N; ++j) 
	    {
            	if ((v = fabs(A[j][i])) > a_max) 
		{
                    a_max = v;
                    row_max = j;
            	}  
            }
            if (a_max < CIM_TOL) 
	    {
            	printf("No pivoting element! Matrix is close to singular.\n");
            	result = 0;
            }
            j = row_max;
            // swap pivoting element
            if (j != i) 
	    {
            	Swap_A = A[i]; A[i] = A[j]; A[j] = Swap_A;
            	v = x[i]; x[i] = x[j]; x[j] = v;
            }
            for (j = i+1; j < N; ++j) 
	    {
            	v = A[j][i] / A[i][i];
            	x[j] -= v*x[i];
            	for (k = i+1; k < N; ++k) 
		{
                    A[j][k] -= v*A[i][k];
            	}
            }
    	}
    	v = 1.0;
    	for (i = 0; i < N; ++i) 
	    v *= A[i][i];
    	for (i = N-1; i >= 0; i--) 
	{
            for (j = i+1; j < N; ++j) 
	    {
            	x[i] -= A[i][j]*x[j];
            }
            x[i] /= A[i][i];
    	}
    	result = 1;
    	FT_FreeThese(2,A0,A);
    	return result;
}


static boolean Linear_Solver(
	MTX A,
        double *x,
        double *b)
{
	int i,j,k;
	PETSc solver;

	solver.Create(0,A.N-1,9,9);
	for (k = 0; k < A.K; ++k)
	{
	    i = A.i[k];
	    j = A.j[k];
	    solver.Add_A(i,j,A.a[k]);
	}
	for (k = 0; k < A.N; ++k)
	    solver.Add_b(k,b[k]);

	start_clock("petsc_solve");
        solver.SetMaxIter(500);
        solver.SetTol(1e-10);
	if (use_neumann)
	{
	    printf("Using Neumann Solver\n");
            solver.Solve_withPureNeumann();
	}
	else
	{
	    printf("Using Dirichlet Solver\n");
            solver.Solve();
	}
	solver.Get_x(x);
        stop_clock("petsc_solve");

	return YES;
}	/* Linear_Solver */

static int tangential_direction(
	int i, 
	double *T, 
	double *N,
	int dim)
{
    	int j;
    	double L;
    
    	for (j = 0; j < dim; ++j) 
	    T[j] = 0.0;
	L = sqrt(1.0-N[i]*N[i]);
	if ( L > CIM_TOL ) 
	{
	    for(j=0;j<dim;++j) 
	    {
		T[j] = -N[i]*N[j];
		T[j] += (i==j) ? 1 : 0;
		T[j] /= L;
	    }
	    j = 1;
	} 
	else 
	{
    	    j = 0;
    	}
    	return j;
}

void CIM_ELLIPTIC_SOLVER::cimSolveFrontState()	
{
        INTERFACE *intfc = front->interf;
        POINTER sl,sr;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
	double ul,ur;

	if (debugging("trace"))
	    (void) printf("Entering cimSolveFrontState()\n");
	next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
            if (Boundary_hs(hs))
                continue;
            FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
	    cimIntfcPointState(Coords(p),negative_component(hs),&ul,&ur);
            assignStateVar(ul,sl);
            assignStateVar(ur,sr);
        }
}	/* end cimSolveFrontState */

void CIM_ELLIPTIC_SOLVER::cimIntfcPointState(
	double *coords,
	int comp,
	double *ul,
	double *ur)
{
	// ul = u at coords in omega_-
	// ur = u at coords in omega_+
	int icoords[MAXD],ic[MAXD];
	int i,j,k,index;
	GRID_DIRECTION forward_dir[3] = {EAST,NORTH,UPPER};
	GRID_DIRECTION backward_dir[3] = {WEST,SOUTH,LOWER};
	double nor[MAXD],crx_coords[MAXD];
	boolean status;
	HYPER_SURF *hs;
	
	int N;
	double u, uJ;
   	double epi, epo;
   	double J_u[2*MAXD], J_eps_un[2*MAXD], J_ut[2*MAXD];
	double *h = top_h;
   	CIM_STRUCT CS;
   	CIM_COEF CC;
	GRID_DIRECTION dir[3][2] = {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
	static double max_diff = 0.0;
	double P[MAXD],R[MAXD],Q[MAXD];
	double a[MAXD];
	double uex;
	double u_xx[MAXD], u_b[MAXD], u_xy[MAXD][MAXD], uc, ut, us, ud, ue;
	double u_at_itfc;
	int ii;

	rect_in_which(coords,icoords,top_grid);
	index = d_index(icoords,top_gmax,dim);
	k  = ij_to_I[icoords[0]][icoords[1]];
	uJ = (*solutionJump)(jparams,D[k],coords);
	u = soln[index];
	for (i = 0; i < dim; ++i) 
	{
	    ic[i] = icoords[i];
	    P[i]  = top_L[i] + icoords[i]*top_h[i];
	    a[i]  = (coords[i]-P[i])/top_h[i];
	}


	epi = (D[k] ==  1) ? diff_coeff[0] : diff_coeff[1];
    	epo = (D[k] == -1) ? diff_coeff[0] : diff_coeff[1];

	for (i = 0;i < dim; ++i) ic[i] = icoords[i];

	for (i = 0; i < dim; ++i) 
	{
            CS.s[i] = S[2*i+1][k]-S[2*i][k];
            CS.a[i] = 1.0;
            for (j = 0; j < dim; ++j) 
	    {
		CS.n[i][j] = (i == j) ? 1 : 0;
	    }
	}
	for (i = 0; i < dim; ++i) 
	{
            J_u[i] = J_eps_un[i] = J_ut[i] = 0;
            if (CS.s[i] != 0) 
	    {
		for (j = 0; j < dim; ++j) 
		    ic[j] = icoords[j];
		if (CS.s[i] == -1)
		{
		    status = FT_NormalAtGridCrossing(front,ic,
					dir[i][0],comp,R,&hs,Q);
		    if (comp == negative_component(hs))
		    {
			for (j = 0; j < dim; ++j) 
		    	    R[i] *= -1;
		    }
		}
		else
		{
		    status = FT_NormalAtGridCrossing(front,ic,
					dir[i][1],comp,R,&hs,Q);
		    if (comp == negative_component(hs))
		    {
			for (j = 0; j < dim; ++j) 
			    R[i] *= -1;
		    }
		}
		CS.a[i] = (Q[i] - P[i])/top_h[i]/CS.s[i];
		J_u[i] = (*solutionJump)(jparams,D[k],Q);
		J_eps_un[i] = (*gradJumpDotN)(jparams,D[k],R,Q);
		J_ut[i] = (*gradJumpDotT)(jparams,D[k],i,R,Q);
		for (j = 0; j < dim; ++j) 
		{
		    CS.n[i][j] = R[j];
		    CS.c[i][j] = CS.s[j];
		}
	    }
    	}
	N = CIM2(&CC,CS,epi,epo,dim);


	uc = soln[index];
	u_at_itfc = -uc;
	for(i = 0; i < dim; ++i) 
	{	
	    u_xx[i] = u_b[i] = 0.0;
	    for (j = 0; j < N; ++j)
    	    {
       		for (ii = 0; ii < dim; ++ii)
           	    ic[ii] = icoords[ii]+CC.indx[ii][j];
	  	index = d_index(ic,top_gmax,dim);
        	u_xx[i] += CC.coef[i][j]*soln[index];
    	    }
            for (j = 0; j < dim; ++j)
            {
       	    	u_xx[i] += (CC.J[i][3*j]*J_eps_un[j]*h[i]
                       +CC.J[i][3*j+1]*J_u[j]
                       +CC.J[i][3*j+2]*J_ut[j]*h[i]);
            }
	    if(CS.s[i] != 0) 
	    {
		for (ii = 0; ii < dim; ++ii) 
		    ic[ii] = icoords[ii];
		ic[i] -= CS.s[i];
		ut = soln[d_index(ic,top_gmax,dim)];
		u_b[i] = uc + a[i]*CS.s[i]*(uc-ut)
			+(1+a[i]*CS.s[i])*a[i]*CS.s[i]*0.5*u_xx[i];
	    }
	    else
	    {
		for (ii = 0; ii < dim; ++ii) 
		    ic[ii] = icoords[ii];
		ic[i] += 1;
		ut = soln[d_index(ic,top_gmax,dim)];
		ic[i] -= 2;
		us = soln[d_index(ic,top_gmax,dim)];
		u_b[i] = uc + a[i]*0.5*(ut-us) + 0.5*a[i]*a[i]*(ut+us-2*uc);
	    }
	    u_at_itfc += u_b[i];
	}
	
    	for (i = 0; i < dim; ++i) 
        for (j = i+1; j < dim; ++j) 
	{
            if(CS.s[i] != 0 && CS.s[j] != 0) 
	    {
		for (ii = 0; ii < dim; ++ii) 
		    ic[ii] = icoords[ii];
		ic[i] -= CS.s[i];
		ut = soln[d_index(ic,top_gmax,dim)];
		ic[i] += CS.s[i];
		ic[j] -= CS.s[j];
		us = soln[d_index(ic,top_gmax,dim)];
		ic[i] -= CS.s[i];
		ud = soln[d_index(ic,top_gmax,dim)];
                u_xy[i][j] = (uc - ut - us + ud)*CS.s[i]*CS.s[j];
            } 
	    else if (CS.s[i] !=0 && CS.s[j] == 0) 
	    {
		for (ii = 0; ii < dim; ++ii) 
		    ic[ii] = icoords[ii];
		ic[j]++; 
		ud = soln[d_index(ic,top_gmax,dim)];
		ic[j] -= 2;
		ut = soln[d_index(ic,top_gmax,dim)];
		ic[i] -= CS.s[i];
		ue = soln[d_index(ic,top_gmax,dim)];
		ic[j] += 2;
		us = soln[d_index(ic,top_gmax,dim)];
                u_xy[i][j] = (ud - ut - us + ue)*CS.s[i]*0.5;
            } 
	    else if(CS.s[i] ==0 && CS.s[j] != 0) 
	    {
		for (ii = 0; ii < dim; ++ii) 
		    ic[ii] = icoords[ii];
		ic[i]++; 
		ud = soln[d_index(ic,top_gmax,dim)];
		ic[i] -= 2;
		ut = soln[d_index(ic,top_gmax,dim)];
		ic[j] -= CS.s[j];
		ue = soln[d_index(ic,top_gmax,dim)];
		ic[i] += 2;
		us = soln[d_index(ic,top_gmax,dim)];
                u_xy[i][j] = (ud - ut - us + ue)*CS.s[j]*0.5;
            } 
	    else 
	    {
		for ( ii = 0; ii < dim; ++ii) 
		    ic[ii] = icoords[ii];
		ic[i] -= (a[i] > 0)? 1 : -1;
		ut = soln[d_index(ic,top_gmax,dim)];
		ic[i] += (a[i] > 0)? 1 : -1;
		ic[j] -= (a[j] > 0)? 1 : -1;
		us = soln[d_index(ic,top_gmax,dim)];
		ic[i] -= (a[i] > 0)? 1 : -1;
		ud = soln[d_index(ic,top_gmax,dim)];
                u_xy[i][j] = (uc - ut - us + ud);
		u_xy[i][j] *= (a[i] > 0)? 1 : -1;
		u_xy[i][j] *= (a[j] > 0)? 1 : -1;
            }
            u_at_itfc += a[i]*a[j]*u_xy[i][j];
    	}
	if (D[k] < 0) 
	{
	    (*ul) = u_at_itfc;
	    (*ur) = u_at_itfc+uJ;
	}
	else
	{
	    (*ur) = u_at_itfc;
	    (*ul) = u_at_itfc+uJ;
	}
}	/* end cimIntfcPointState */

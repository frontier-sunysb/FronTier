int CIM_ELLIPTIC_SOLVER::HCIM_Matrix_Generation()
{
    	int i,j,n,k,ii,jj,index,s,N,Ri;
    	int *row, *col;
    	int icoords[MAXD], ic[MAXD];
    	double P[MAXD], Q[MAXD], R[MAXD];
    	double epi, epo;
    	double J_u[2*MAXD], J_eps_un[2*MAXD], J_ut[2*MAXD];
    	double v, *val;
	double *h = top_h;
	HYPER_SURF *hs;
	double crx_coords[MAXD];
	COMPONENT comp;
	GRID_DIRECTION dir[2][2] = {{WEST,EAST},{SOUTH,NORTH}};
	int status;
	int icrds[MAXD];
	int idir,nb;
	double kp;
	int index_nb;
    
    	row = A.i;
    	col = A.j;
    	val = A.a;
    
    	CIM_STRUCT CS;
    	CIM_COEF CC;

	n = 0; Ri = 0;
	use_neumann = YES;
	for (ii = imin; ii <= imax; ++ii)
	for (jj = jmin; jj <= jmax; ++jj)
	{
	    k = ij_to_I[ii][jj];
	    for (i = 0; i < dim; ++i) 
	    {
		icoords[i] = I_to_ij[k][i];
            	P[i] = cell_edge(icoords[i],i,top_grid);
            }
	    index = d_index(icoords,top_gmax,dim);
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
		    if (NB[i][k] >= 0) 
		    {
			row[n] = k;
			col[n] = NB[i][k];
			val[n] = -1.0/sqr(h[idir])*epi;
			n++;
		    } 
		    else 
		    {
			// boundary treatment
			HYPER_SURF *hs;
			GRID_DIRECTION dir[2][2] ={{WEST,EAST},{SOUTH,NORTH}};
			POINTER state;
			double crx_coords[MAXD];

			comp = top_comp[index];
			status = (*findStateAtCrossing)(front,icoords,
				dir[idir][nb],comp,&state,&hs,crx_coords);
			if (status == DIRICHLET_PDE_BOUNDARY)
			{
			    b[k] += 1.0/sqr(h[idir])*epi*getStateVar(state);
			    use_neumann = NO;
			}
			else if (status == NEUMANN_PDE_BOUNDARY)
			{
			    neumann_nb = YES;
		    	    v -= 1.0/sqr(h[idir])*epi;
			}
		    }
		    v += 1.0/sqr(h[idir])*epi;
		}
		row[n] = col[n] = k;
		val[n] = v;
		n++;
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
			    status = FT_NormalAtGridCrossing(front,icrds,
					dir[i][0],comp,R,&hs,crx_coords);
			else
			    status = FT_NormalAtGridCrossing(front,icrds,
					dir[i][1],comp,R,&hs,crx_coords);
			for (j = 0; j < dim; ++j) 
			    Q[j] = crx_coords[j];	
			CS.a[i] = (crx_coords[i] - P[i])/top_h[i]/CS.s[i];

			J_u[i] = (*solutionJump)(D[k],Q);
			J_eps_un[i] = (*gradJumpDotN)(D[k],R,Q,dim);
			J_ut[i] = (*gradJumpDotT)(D[k],i,R,Q,dim);
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
		    row[n] = k;
		    for (i = 0; i < dim; ++i) 
		    {
			ic[i] = icoords[i] + CC.indx[i][j];
		    }
		    s = ij_to_I[ic[0]][ic[1]];
		    if (s >= 0) 
		    {
			col[n] = s;
			v = 0;
			for (i = 0; i < dim; ++i) 
			{
			    v += CC.coef[i][j]/h[i]/h[i]*epi;
			}
			val[n] = -v;
			n++;
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
		printf("Order[%d] = %d\n",k,Order[k]);
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
			printf("Position 2\n");
			if (nb == -1)
			    status = FT_NormalAtGridCrossing(front,icrds,
					dir[idir][0],comp,R,&hs,crx_coords);
			else
			    status = FT_NormalAtGridCrossing(front,icrds,
					dir[idir][1],comp,R,&hs,crx_coords);
			printf("status = %d\n",status);
			printf("crx_coords = %f %f R = %f %f\n",
					crx_coords[0],crx_coords[1],R[0],R[1]);
			for (j = 0; j < dim; ++j) 
			    Q[j] = crx_coords[j];	
			CS.a[i] = (crx_coords[i] - P[i])/top_h[i]/CS.s[i];
					
			J_u[i] = (*solutionJump)(D[k],Q);
			J_eps_un[i] = (*gradJumpDotN)(D[k],R,Q,dim);
			J_ut[i] = (*gradJumpDotT)(D[k],i/2,R,Q,dim);
			for (j = 0; j < dim; ++j) 
			{
			    CS.n[i][j] = R[j];
			}
		    }
	        }
	        N = CIM1(&CC,CS,epi,epo,dim); // CIM1
	        for (j = 0; j < N; ++j) 
		{
		    row[n] = k;
		    for (i = 0; i < dim; ++i) 
		    {
			ic[i] = icoords[i] + CC.indx[i][j];
		    }
		    v = 0;
		    for (i = 0; i < dim; ++i) 
		    {
			v += CC.coef[i][j];
		    }
		    s = ij_to_I[ic[0]][ic[1]];
		    if (s >= 0) 
		    {
			col[n] = s;
			v = 0;
			for  (i = 0; i < dim; ++i) 
			{
			    v += CC.coef[i][j]/h[i]/h[i]*epi;
			}
			val[n] = -v;
			n++;
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
    	return n;
}	/* end HCIM_Matrix_Generation */

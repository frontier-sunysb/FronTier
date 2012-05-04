/*
*
*	Embedded Boundary method for parabolic/elliptic solver 
*       
*/

#include <FronTier.h>
#include "poisson.h"

static  int   count_local_grid_blocks(Front*,BLK_EDGE**, int**);
static	void  create_global_element_list(Front*,int,int,int*, int**, int**);
static  void  set_right_hand_side(Vec,P_CELL2D**,int,int,int**);
static  void  finite_difference_stencil(POINTER,P_CELL2D*,double*,double*,
			int,int);
static  void  init_stencil_element2D(BLK_EDGE**,P_CELL2D**,Front*,POINTER,
			RECT_GRID*);
static  void  set_matrix_row(Mat,POINTER,Front*,BLK_EDGE**,P_CELL2D**,int,int,
			int**,int**);
static  void  stencil_linear_interpolator(POINTER,P_CELL2D*,EDGE,double*,
			double*,double*,int);
static 	void  scatter_index(Front*,int**,int**);
static 	void  pack_index_in_dir(int**,int*,int**,int,int);
static 	void  unpack_index_in_dir(int**,int*,int**,int,int);
static 	void  reflect_index_in_dir(int**,int*,int,int);
static 	void  comp_ic_to_top_ic(int*,int*,Front*);
static  void  extract_solution(Vec,POINTER,Front*,int,
                        int,int**,int**);
static 	double quadric_fit_slope(double,double,double,double,double,double);

static  double     tol = 1e-3;
static  int       num_nodes;

static	double     K[2] = {0.5,0.5};

/*
 *   We solve the elliptic problem,
 *
 *   - \Delta phi = div(F),
 *
*/

extern  void            embed_bdry_poison(
        POINTER         wv,
        Front           *fr)
{
        int             myid,i,j;
        int             NLblocks,ilower,iupper,I;
        int             *n_dist;
        int             Nx = fr->rect_grid->gmax[0];
        int             Ny = fr->rect_grid->gmax[1];
        
        int             ix,iy;
        int             nlb = 0;
        int             icoord[2], ICOORDS[2];
        int             **I_to_ijk;   /* local  i, j, k, local I */
        int             **IJK_to_I;  /* global I, J, K, global I */
        double           hx = fr->rect_grid->h[0];
        double           hy = fr->rect_grid->h[1];
        double           J_comp[2];
        double           norm,min,max,x,y,*coords;
        P_CELL2D        **p_cell;
        POINTER        state;

            /*PETSC Variables*/
        KSP             ksp;
        PC              pc;
        Mat             A;
        Vec             phi, F;
        PetscScalar     *vphi;
        MatNullSpace    nullsp;

	/* Variables for debugging (extracting values) */
	int             indices[1];
        int             its;
        FILE            *file1;
        POINT           *pt;
        BOND            *bond;
        CURVE           *c,**cv;
	BLK_EDGE	**blk_edge;
	double		***grad_phi;
	INTERFACE	*emb_grid_intfc = fr->emb_grid_intfc;
	struct Table	*T;
	LAPLACE		*laplace = (LAPLACE*)wv;
	double		*solute = laplace->solute;
	int		comm_size;

	printf("Entering embed_bdry_poison()\n");
	num_nodes = pp_numnodes();

	MPI_Comm_size(PETSC_COMM_WORLD,&comm_size);
        FT_MatrixMemoryAlloc((POINTER*)&p_cell,Nx,Ny,sizeof(P_CELL2D));

	T = table_of_interface(emb_grid_intfc);
	blk_edge = T->blk_edge;
	printf("Calling init_stencil_element2D()\n");
        init_stencil_element2D(blk_edge,p_cell,fr,wv,fr->rect_grid);

        myid = pp_mynode();
	printf("Calling find_Cartesian_coordinates()\n");
	find_Cartesian_coordinates(myid,fr->pp_grid,ICOORDS);

	FT_VectorMemoryAlloc((POINTER*)&n_dist,num_nodes,sizeof(int));
	for (i = 0; i < num_nodes; ++i) n_dist[i] = 0;
	FT_MatrixMemoryAlloc((POINTER*)&I_to_ijk,Nx*Ny,2,sizeof( int));
	FT_MatrixMemoryAlloc((POINTER*)&IJK_to_I,(Nx+4),(Ny+4),sizeof( int));
	FT_TriArrayMemoryAlloc((POINTER*)&grad_phi,Nx,Ny,2,sizeof(double));
	
	printf("Calling count_local_grid_blocks()\n");
        NLblocks = count_local_grid_blocks(fr,blk_edge,I_to_ijk);
	n_dist[myid] = NLblocks;
	pp_global_imax(n_dist,num_nodes);

        ilower = 0;
        iupper = n_dist[0]-1;

        for (i = 1; i <= myid; i++)
        {
            ilower += n_dist[i-1];
            iupper += n_dist[i];
        }
	printf("ilower = %d  iupper = %d\n",ilower,iupper);
	printf("Calling create_global_element_list()\n");
        create_global_element_list(fr,ilower,iupper,ICOORDS,I_to_ijk,IJK_to_I);

        MatCreateMPIAIJ(PETSC_COMM_WORLD,iupper-ilower+1,iupper-ilower+1,
              PETSC_DECIDE,PETSC_DECIDE,7,PETSC_NULL,7,PETSC_NULL,&A);
 
        VecCreate(PETSC_COMM_WORLD,&F);
        VecSetSizes(F,(iupper-ilower + 1),PETSC_DECIDE);
        VecSetFromOptions(F);
        VecDuplicate(F,&phi); 
	
	for (I = ilower; I <= iupper; I++)
	{
            set_matrix_row(A,wv,fr,blk_edge,p_cell,I,ilower,I_to_ijk,IJK_to_I);
            set_right_hand_side(F,p_cell,I,ilower,I_to_ijk);
        }
        
        MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);

        VecAssemblyBegin(F);
        VecAssemblyEnd(F);

        MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,PETSC_NULL,&nullsp);
        KSPCreate(PETSC_COMM_WORLD,&ksp);
        KSPSetOperators(ksp,A,A,DIFFERENT_NONZERO_PATTERN);

        KSPGetPC(ksp,&pc);
        PCSetType(pc,PCBJACOBI);

        KSPSetNullSpace(ksp,nullsp);
        KSPSetType(ksp,KSPGMRES);
        KSPSetTolerances(ksp,1.e-15,PETSC_DEFAULT,PETSC_DEFAULT,30);
        KSPSetComputeSingularValues(ksp, PETSC_TRUE);
      
        KSPSetFromOptions(ksp);
        KSPSetUp(ksp);

        KSPSolve(ksp,F,phi);
        KSPGetIterationNumber(ksp,&its);
        KSPGetResidualNorm(ksp,&norm);

        KSPComputeExtremeSingularValues(ksp, &max, &min);
        printf("\nNumber of iteration is %d\n",its);
	printf("Nesidual norm = %6.5e\n",norm);
	printf("Extreme test SV is [%7.6e ; %7.6e]\n",min,max);

        VecGetArray(phi,&vphi);
        if (its == 0 || norm > 0.001)
        {
            printf("\n WARNING,No convergent happens!");
            VecDestroy(F);
            VecDestroy(phi);
            MatDestroy(A);
            KSPDestroy(ksp);
            MatNullSpaceDestroy(nullsp);

            FT_FreeThese(3,n_dist,I_to_ijk,IJK_to_I);
            FT_FreeThese(1,p_cell);
            return;
        }

	/* Extract solution from the linear solution */

	printf("Calling extract_solution()\n");
	extract_solution(phi,wv,fr,ilower,iupper,I_to_ijk,IJK_to_I);

        VecDestroy(F);
        VecDestroy(phi);
        MatDestroy(A);
        KSPDestroy(ksp);
        MatNullSpaceDestroy(nullsp);
               
        FT_FreeThese(3,n_dist,I_to_ijk,IJK_to_I);
        FT_FreeThese(1,p_cell);
        return;
}

static  void            set_matrix_row(
	Mat             A,
	POINTER         wave,
	Front           *fr,
        BLK_EDGE        **blk_edge,
        P_CELL2D        **p_cell,
        int             I,
	int             ilower,
        int             **I_to_ijk,
        int             **IJK_to_I)
{
        double          hx = fr->rect_grid->h[0];
	double          hy = fr->rect_grid->h[1];
	double		*L = fr->rect_grid->L;
	double          hxs = hx*hx;
	double          hys = hy*hy;
        double          beta;
	int             i,ix,iy;
        int             icoords[2],ic_top[2],idxm[1];
        static POINTER  state[4];
	CRXING		**crxings;
	INTERFACE	*grid_intfc = fr->grid_intfc;
	RECT_GRID       *gr = &topological_grid(grid_intfc);
	GRID_DIRECTION	dir[4] = {WEST,EAST,SOUTH,NORTH};
	LAPLACE         *laplace = (LAPLACE*)wave;
        POINTER         sparams = laplace->sparams;
        double          D = laplace->D;
        double          (*sfunc)(POINTER,double*) = laplace->sfunc;

        idxm[0] = I;

	FT_VectorMemoryAlloc((POINTER*)&crxings,10,sizeof(CRXING*));
        icoords[0] = ix = I_to_ijk[I-ilower][0];
        icoords[1] = iy = I_to_ijk[I-ilower][1];
	comp_ic_to_top_ic(icoords,ic_top,fr);
 
 	if (state[0] == NULL)
	{
            for(i = 0; i < 4; ++i)
              FT_ScalarMemoryAlloc((POINTER*)&state[i],size_of_state(fr->interf));
	}

        if(blk_edge[ix][iy].ctype == INTERNAL)
        {
	    /* Interior point */
	    int             cols[5];
	    double          coords[MAXD],sigma[4];
            double          values[5];
            double          x,y;

	    coords[0] = gr->L[0] + ic_top[0]*gr->h[0];
            coords[1] = gr->L[1] + ic_top[1]*gr->h[1];
	    p_cell[ix][iy].r_side = -(*sfunc)(sparams,coords);

	    cols[0] = I; 
	    cols[1] = IJK_to_I[ix+1][iy+2];  
	    cols[2] = IJK_to_I[ix+3][iy+2];  
	    cols[3] = IJK_to_I[ix+2][iy+1];
	    cols[4] = IJK_to_I[ix+2][iy+3];
           
	    for (i = 0; i < 4; ++i)
	    	sigma[i] = laplace->D;
	    values[0] = (sigma[0]+sigma[1])/hxs + (sigma[2]+sigma[3])/hys;
	    values[1] = -sigma[0]/hxs;
	    values[2] = -sigma[1]/hxs;
	    values[3] = -sigma[2]/hys;
	    values[4] = -sigma[3]/hys;
           
            MatSetValues(A,1,idxm,5,cols,values,INSERT_VALUES);
        }
        else
        {
            int  sten_index[3][3];   
	    		//set up the matrix element for partial cell here
            int  i,j;

            for (i = 0; i < 3; ++i)
	    {
              	for (j = 0; j < 3; ++j)
		{
                    sten_index[i][j] = IJK_to_I[ix+i+1][iy+j+1];
		}
	    }

            for(i = 0; i < 3; ++i)
	    {
               	for(j = 0; j < 3; ++j)
               	{
                    if(p_cell[ix][iy].set[i][j] == YES && sten_index[i][j] < 0)
                    	printf("\n WARNING!, negative index for interpolated"
			       "stencil!");

                   if(sten_index[i][j] >= 0 && p_cell[ix][iy].set[i][j] == YES)
                   	MatSetValues(A,1,idxm,1,&sten_index[i][j],
		   		&(p_cell[ix][iy].stencil[i][j]),INSERT_VALUES); 
               	}
	    }
        }
        return;
}

static   void	        set_right_hand_side(
	Vec             F,
        P_CELL2D        **p_cell,
	int             I,
	int             ilower,
        int             **I_to_ijk)
{
        int             ix,iy,indices[1];
	double           values[1];

        ix = I_to_ijk[I-ilower][0];
        iy = I_to_ijk[I-ilower][1];
        
	indices[0]=I;       
	
        values[0] = p_cell[ix][iy].r_side;
        VecSetValues(F,1,indices,values,INSERT_VALUES); 
        
        return;
}

static   int             count_local_grid_blocks(
	Front           *fr,
        BLK_EDGE        **blk_edge,
        int             **I_to_ijk)
{        
        int             Nx = fr->rect_grid->gmax[0];
	int             Ny = fr->rect_grid->gmax[1];
	int             ix, iy;
	int             nlb = 0;

	for (iy = 0; iy < Ny; iy++)
	{
	    for (ix = 0; ix < Nx; ix++)
	    {
                if(blk_edge[ix][iy].ctype)
		{
                      I_to_ijk[nlb][0] = ix;
		      I_to_ijk[nlb][1] = iy;
		      nlb++;
                }
	    }
	}
	printf("nlb = %d\n",nlb);
        return nlb;
      
} /* end count_local_grid_blocks() */

static  void            create_global_element_list(
	Front           *fr,					   
        int             ilower,
	int             iupper,
	int             *ICOORDS,
        int             **I_to_ijk,
        int             **IJK_to_I)
{
        int             **bfs,**bfr;
        int             Nx = fr->rect_grid->gmax[0];
	int             Ny = fr->rect_grid->gmax[1];
        int             NNx = fr->pp_grid->gmax[0];
      	int             NNy = fr->pp_grid->gmax[1];
	int             I, ix, iy;
	int             icoords[2];
	int             send,recv;
        int             IIX = ICOORDS[0];
        int             IIY = ICOORDS[1];  /* subdomain icoords */ 
        int             his_id, tag;

        MPI_Status      status; 

        for (iy = 0; iy < Ny + 4; iy++)
	{
            for (ix = 0; ix < Nx + 4; ix++)
            {
                IJK_to_I[ix][iy] = -1;
            }
	}
        for (I = ilower; I <= iupper; I++)
        {
            ix = I_to_ijk[I-ilower][0];
            iy = I_to_ijk[I-ilower][1];
            IJK_to_I[ix+2][iy+2] = I;
        }
	//scatter_index(fr,I_to_ijk,IJK_to_I);
	int lbuf[2] = {1,2};
	int ubuf[2] = {1,2};
	FT_ParallelExchCellIndex(fr,lbuf,ubuf,(POINTER)IJK_to_I);

	return;
}

static void   init_stencil_element2D(
         BLK_EDGE      **blk_edge,
         P_CELL2D      **p_cell,
         Front         *fr,
         POINTER       wv,
         RECT_GRID     *gr)
{
       int     	xmax,ymax;
       int     	ix,iy;
       int     	i,j,k,dir;
       int     	icoords[MAXD];
       double   coords[MAXD],coord[MAXD];      
       	       	//coords is the cell center and coord is the cell edge center
       double   u[3],B[3],uxB[3];
       double   sigma,beta; 
       static 	POINTER   state = NULL;
       double   *h = gr->h;

       if(!state)
           FT_ScalarMemoryAlloc((POINTER*)&state,size_of_state(fr->interf));
 
       	xmax = gr->gmax[0];
       	ymax = gr->gmax[1];

       	u[2] = 0.0;
       	for (iy = 0;  iy < ymax;  ++iy)           
       	{
            for (ix = 0;  ix < xmax;  ++ix)
            {
              	icoords[0] = ix; icoords[1] = iy; 
	      	coords[0] = gr->L[0] + (ix + 0.5)*gr->h[0];
	      	coords[1] = gr->L[1] + (iy + 0.5)*gr->h[1];

              	for(i = 0; i < 3; ++i)
		{
                    for(j = 0; j < 3; ++j)
                    {
                    	p_cell[ix][iy].stencil[i][j] = 0.0;
                    	p_cell[ix][iy].set[i][j] = NO;
                    }
		}

              	p_cell[ix][iy].r_side = 0.0;

              	if (blk_edge[ix][iy].ctype == PARTIAL)  //if it is partial cell
              	{
                    for (i = 0; i < 2; ++i)
		    {
                     	for (j = 0; j < 2; ++j)
                     	{
			    if (blk_edge[ix][iy].edge[i][j].etype 
			    		== PARTIAL_EDGE)
			    {
			    	coord[0] = blk_edge[ix][iy].edge[i][j].cen[0];
				coord[1] = blk_edge[ix][iy].edge[i][j].cen[1];
			    }
			    else if (blk_edge[ix][iy].edge[i][j].etype
			    		== FULL_EDGE)
			    {
                            	dir = (j == 0 ? -1 : 1);
				coord[i] = coords[i]+ dir*h[i]/2;
				coord[(i+1)%2] = coords[(i+1)%2];
			    }
			    else
			    {
			    	continue;
			    }

                            if (blk_edge[ix][iy].edge[i][j].etype == 
			    			PARTIAL_EDGE)
			    {
                              	stencil_linear_interpolator(wv,
					&(p_cell[ix][iy]),
                                	blk_edge[ix][iy].edge[i][j],
					coord,coords,h,i);    
			    }
                            else if(blk_edge[ix][iy].edge[i][j].etype == 
			    			FULL_EDGE)   
			    {
				finite_difference_stencil(wv,&(p_cell[ix][iy]),
						coord,h,i,dir);
			    }
                    	}
		    }
              	}
            }
       	}
}	/* end init_stencil_element2D */

LOCAL  void  stencil_linear_interpolator(
	POINTER	  wave,
       	P_CELL2D  *p_cell,
        EDGE      edge,
        double     *coord,
        double     *coords,
        double     *h,
        int       i)
{
        int       k,m;
        int       e[2][2];
        double     ratio,crds_new[2];    //new coords of the edge center
	LAPLACE   *laplace = (LAPLACE*)wave;
        POINTER   sparams = laplace->sparams;
        void      (*flux_func)(POINTER,double*,double*) = laplace->flux_func;
        double    (*sfunc)(POINTER,double*) = laplace->sfunc;
        double    flux[2];
        double    s1,s2,old_rhs,new_rhs;
 
        for (k = 0; k < 2; ++k)
            for(m = 0; m < 2; ++m)
                e[k][m] = 0;
 
        for (k = 0; k < 2; ++k)
        {
            crds_new[k] = edge.cen[k]-coords[k];

            if (k == i)
                 continue;
            e[1][k] = (crds_new[k] >= 0? 1 : -1);
            ratio = fabs(crds_new[k]/h[k]);
        }   
        if (crds_new[i] >= 0)
            e[0][i] = 1;
        else
            e[0][i] = -1;    
        
	(*flux_func)(sparams,coord,flux);
        p_cell->r_side += -flux[i]*e[0][i]*edge.length;

        p_cell->stencil[1][1] += (1-ratio)*edge.length/h[i];
        p_cell->set[1][1] = YES;

        p_cell->stencil[1+e[0][0]][1+e[0][1]] += (ratio-1)*edge.length/h[i];
        p_cell->set[1+e[0][0]][1+e[0][1]] = YES;

        p_cell->stencil[1+e[1][0]][1+e[1][1]] += ratio*edge.length/h[i];
        p_cell->set[1+e[1][0]][1+e[1][1]] = YES;

        p_cell->stencil[1+e[0][0]+e[1][0]][1+e[0][1]+e[1][1]] 
                                    += -ratio*edge.length/h[i];
        p_cell->set[1+e[0][0]+e[1][0]][1+e[0][1]+e[1][1]] = YES;
}	/* end stencil_linear_interpolator */

LOCAL  void  finite_difference_stencil(
	POINTER	  wave,
        P_CELL2D  *p_cell,
        double     *coord,
        double     *h,
        int       i,
        int       dir)
{
	LAPLACE   *laplace = (LAPLACE*)wave;
        POINTER   sparams = laplace->sparams;
       	void      (*flux_func)(POINTER,double*,double*) = laplace->flux_func;
        double    (*sfunc)(POINTER,double*) = laplace->sfunc;
        double    flux[2];
        double    s1,s2,old_rhs,new_rhs;

	(*flux_func)(sparams,coord,flux);
        p_cell->r_side += -flux[i]*dir*h[(i+1)%2];

        p_cell->stencil[1][1] += h[(i+1)%2]/h[i];
        p_cell->set[1][1] = YES;

        if (i == 0)
        {
            p_cell->set[1+dir][1] = YES;
            p_cell->stencil[1+dir][1] += -h[(i+1)%2]/h[i];
        }
        else if(i == 1)
        {
            p_cell->set[1][1+dir] = YES;
            p_cell->stencil[1][1+dir] += -h[(i+1)%2]/h[i];
        }
}	/* end finite_difference_stencil */	

static double  quadric_fit_slope(
        double        x0,
	double        x1,
	double        x2,
	double        y0,
	double        y1,
	double        y2)
	   
{           
	double   K1,K2;
        double   slope;
	    
	K1 = (y1-y0)/(x1-x0);
	K2 = (y2-y1)/(x2-x1);
	K2 = (K2-K1)/(x2-x0);

	slope = 2.0*K2*x0+(K1-K2*(x0+x1));
	return slope;
}

static void scatter_index(
	Front *fr,
	int **I_to_ijk,
	int **IJK_to_I)
{
	INTERFACE *intfc = fr->interf;
	int       me[MAXD], him[MAXD];
	int       myid, dst_id;
	PP_GRID   *pp_grid = fr->pp_grid;
	RECT_GRID *gr = fr->rect_grid;
	int       *G = pp_grid->gmax;
	int       i, j, k;
        int       dim = gr->dim;
	int       *gmax = gr->gmax;
	int       **bfs,**bfr;
	int	  lbuf,ubuf;
	int	  size,max_size;
	int	  index_tag = 8;

	lbuf = 2;
	ubuf = 2;

	max_size = 0;
	for (i = 0; i < dim; ++i) 
	    if (max_size < gmax[i] + 4) max_size = gmax[i] + 4;
	FT_MatrixMemoryAlloc((POINTER*)&bfs,2,max_size,sizeof(int));
        FT_MatrixMemoryAlloc((POINTER*)&bfr,2,max_size,sizeof(int));

	find_Cartesian_coordinates(pp_mynode(),pp_grid,me);
	for (i = 0; i < dim; ++i)
	{
	    for (j = 0; j < 2; ++j)
	    {
		for (k = 0; k < dim; ++k)
		    him[k] = me[k];

		size = (gmax[(i+1)%2] + 4)*2;
	    	if (rect_boundary_type(intfc,i,j) == SUBDOMAIN_BOUNDARY)
		{
		    him[i] = (me[i] + 2*j - 1 + G[i])%G[i];
		    dst_id = domain_id(him,G,dim);
		    pack_index_in_dir(IJK_to_I,gmax,bfs,i,j);
		    pp_send(index_tag,*bfs,size*INT,dst_id);
		}
		else if (rect_boundary_type(intfc,i,j) == REFLECTION_BOUNDARY)
		{
		    reflect_index_in_dir(IJK_to_I,gmax,i,j);
		}

		if (rect_boundary_type(intfc,i,(j+1)%2) == SUBDOMAIN_BOUNDARY)
		{
		    him[i] = (me[i] - 2*j + 1 + G[i])%G[i];
		    dst_id = domain_id(him,G,dim);
		    pp_recv(index_tag,dst_id,*bfr,size*INT);
		    unpack_index_in_dir(IJK_to_I,gmax,bfr,i,(j+1)%2);
		}
	    }
	}
}	/* end scatter_index */

static void pack_index_in_dir(
	int **IJK_to_I,
	int *gmax,
	int **bfs,
	int dir,
	int nb)
{
	int i;
	if (dir == 0)
	{
	    switch(nb)
	    {
	    case 0:
	    	for (i = 0; i < gmax[1] + 4; ++i)
		{
		    bfs[0][i] = IJK_to_I[2][i];
		    bfs[1][i] = IJK_to_I[3][i];
		}
		break;
	    case 1:
	    	for (i = 0; i < gmax[1] + 4; ++i)
		{
		    bfs[0][i] = IJK_to_I[gmax[0]][i];
		    bfs[1][i] = IJK_to_I[gmax[0]+1][i];
		}
	    }
	}
	else if (dir == 1)
	{
	    switch(nb)
	    {
	    case 0:
	    	for (i = 0; i < gmax[0] + 4; ++i)
		{
		    bfs[0][i] = IJK_to_I[i][2];
		    bfs[1][i] = IJK_to_I[i][3];
		}
		break;
	    case 1:
	    	for (i = 0; i < gmax[0] + 4; ++i)
		{
		    bfs[0][i] = IJK_to_I[i][gmax[1]];
		    bfs[1][i] = IJK_to_I[i][gmax[1]+1];
		}
	    }
	}
}	/* end pack_index_in_dir */
	
static void unpack_index_in_dir(
	int **IJK_to_I,
	int *gmax,
	int **bfr,
	int dir,
	int nb)
{
	int i;
	if (dir == 0)
	{
	    switch(nb)
	    {
	    case 0:
	    	for (i = 0; i < gmax[1] + 4; ++i)
		{
		    IJK_to_I[0][i] = bfr[0][i];
		    IJK_to_I[1][i] = bfr[1][i];
		}
		break;
	    case 1:
	    	for (i = 0; i < gmax[1] + 4; ++i)
		{
		    IJK_to_I[gmax[0]+2][i] = bfr[0][i];
		    IJK_to_I[gmax[0]+3][i] = bfr[1][i];
		}
	    }
	}
	else if (dir == 1)
	{
	    switch(nb)
	    {
	    case 0:
	    	for (i = 0; i < gmax[0] + 4; ++i)
		{
		    IJK_to_I[i][0] = bfr[0][i];
		    IJK_to_I[i][1] = bfr[1][i];
		}
		break;
	    case 1:
	    	for (i = 0; i < gmax[0] + 4; ++i)
		{
		    IJK_to_I[i][gmax[1]+2] = bfr[0][i];
		    IJK_to_I[i][gmax[1]+3] = bfr[1][i];
		}
	    }
	}
}	/* end unpack_index_in_dir */
	
static void reflect_index_in_dir(
	int **IJK_to_I,
	int *gmax,
	int dir,
	int nb)
{
	int i;
	if (dir == 0)
	{
	    switch(nb)
	    {
	    case 0:
	    	for (i = 0; i < gmax[1] + 4; ++i)
		{
		    IJK_to_I[0][i] = IJK_to_I[3][i];
		    IJK_to_I[1][i] = IJK_to_I[2][i];
		}
		break;
	    case 1:
	    	for (i = 0; i < gmax[1] + 4; ++i)
		{
		    IJK_to_I[gmax[0]+2][i] = IJK_to_I[gmax[0]+1][i];
		    IJK_to_I[gmax[0]+3][i] = IJK_to_I[gmax[0]][i];
		}
	    }
	}
	else if (dir == 1)
	{
	    switch(nb)
	    {
	    case 0:
	    	for (i = 0; i < gmax[0] + 4; ++i)
		{
		    IJK_to_I[i][0] = IJK_to_I[i][3];
		    IJK_to_I[i][1] = IJK_to_I[i][2];
		}
		break;
	    case 1:
	    	for (i = 0; i < gmax[0] + 4; ++i)
		{
		    IJK_to_I[i][gmax[1]+2] = IJK_to_I[i][gmax[1]+1];
		    IJK_to_I[i][gmax[1]+3] = IJK_to_I[i][gmax[1]];
		}
	    }
	}
}	/* end reflect_index_in_dir */
	

static void comp_ic_to_top_ic(
	int *icoords,
	int *ic_top,
	Front *front)
{
	INTERFACE *grid_intfc = front->grid_intfc;
	RECT_GRID *top_gr = &topological_grid(grid_intfc);
	RECT_GRID *comp_gr = front->rect_grid;
	int i;
	for (i = 0; i < 2; ++i)
	{
	    if (comp_gr->lbuf[i] != 0)
	    	ic_top[i] = icoords[i] + comp_gr->lbuf[i];
	    else
	    	ic_top[i] = icoords[i] + 1;
	}
}

static  void	        extract_solution(
	Vec             phi,
	POINTER         wave,
	Front           *fr,
	int             ilower,
	int             iupper,
        int             **I_to_ijk,
        int             **IJK_to_I)
{
        int             Nx = fr->rect_grid->gmax[0];
        int             Ny = fr->rect_grid->gmax[1];
        int             indices,ix,iy;
        PetscScalar     *values;
	RECT_GRID	*grid;
	LAPLACE		*laplace = (LAPLACE*)wave;
	double		*solute = laplace->solute;
	int		i,j,*gmax;
        double          (*sfunc)(POINTER,double*) = laplace->sfunc;

	grid = &topological_grid(fr->grid_intfc);
	gmax = grid->gmax;

	for (iy = 0; iy < Ny; iy++)
	{
            for (ix = 0; ix < Nx; ix++)
            {
		indices = IJK_to_I[ix+2][iy+2];
		i = ix + 1;	  j = iy + 1;
		if (fr->rect_grid->lbuf[0] != 0)
		    i = ix + fr->rect_grid->lbuf[0];
		if (fr->rect_grid->lbuf[1] != 0)
		    j = iy + fr->rect_grid->lbuf[1];

		double	coords[2];
		coords[0] = grid->L[0]+i*grid->h[0];
		coords[1] = grid->L[1]+j*grid->h[1];
		if (indices < 0)
                {
		    solute[d_index2d(i,j,gmax)] = 0.0;
                }
                else
                { 
                    VecGetArray(phi,&values);
		    solute[d_index2d(i,j,gmax)] = values[indices-ilower];
                    VecRestoreArray(phi,&values);
                }
            }
	}
}	/* extract_potential */

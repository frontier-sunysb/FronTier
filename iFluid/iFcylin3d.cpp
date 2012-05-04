/*******************************************************************
 * 			iFcartsn3d.cpp	
 *******************************************************************/
#include "iFluid.h"
#include "solver.h"

//--------------------------------------------------------------------------
// 		   Incompress_Solver_Smooth_3D_Cylindrical
//--------------------------------------------------------------------------
//
//
//
//
//-----------------------------------------------------------------------------------------------------
//  utility functions related to the computation of advection source term
//-------------------------------------------------------------------------------------------------------

static double (*getStateVel[3])(POINTER) = {getStateXvel,getStateYvel,
                                        getStateZvel};

void Incompress_Solver_Smooth_3D_Cylindrical::compAdvectionTerm_decoupled(void)
{
        int index;
    	int i,j,k,icoords[MAXD];
    	int I;

	setIndexMap();

    	for (k = kmin; k <= kmax; k++)
    	for (j = jmin; j <= jmax; j++)
	for (i = imin; i <= imax; i++)
	{
	     I  = ijk_to_I[i][j][k];
	     if (I == -1) 
	     {
		 continue;
	     }

	     index = d_index3d(i,j,k,top_gmax);

	     icoords[0] = i;
	     icoords[1] = j;
	     icoords[2] = k;

	     double convectionTerm[3];
	     getAdvectionTerm_decoupled(icoords, convectionTerm);

	     cell_center[index].m_state.m_adv[0] = convectionTerm[0];
	     cell_center[index].m_state.m_adv[1] = convectionTerm[1];
	     cell_center[index].m_state.m_adv[2] = convectionTerm[2];
	}

}

void Incompress_Solver_Smooth_3D_Cylindrical::compAdvectionTerm_coupled(void)
{
        int index;
    	int i,j,k,icoords[MAXD];
    	int I;

	setIndexMap();

    	for (k = kmin; k <= kmax; k++)
    	for (j = jmin; j <= jmax; j++)
	for (i = imin; i <= imax; i++)
	{
	     I  = ijk_to_I[i][j][k];
	     if (I == -1) 
	     {
		 continue;
	     }

	     index = d_index3d(i,j,k,top_gmax);

	     icoords[0] = i;
	     icoords[1] = j;
	     icoords[2] = k;

	     double convectionTerm[3];
	     getAdvectionTerm_coupled(icoords, convectionTerm);

	     cell_center[index].m_state.m_adv[0] = convectionTerm[0];
	     cell_center[index].m_state.m_adv[1] = convectionTerm[1];
	     cell_center[index].m_state.m_adv[2] = convectionTerm[2];
	}

}


void Incompress_Solver_Smooth_3D_Cylindrical::getAdvectionTerm_decoupled(
	int *icoords,
	double convectionTerm[3])
{
    bool bNoBoundary;
    int ICoords[3];
    L_STATE sl, sr, state_west, state_east, state_south, state_north, state_lower, state_upper;

    double dtheta = top_h[0];
    double dz 	  = top_h[1];
    double dr     = top_h[2];
    int index = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
    double r = cell_center[index].m_coords[2];

    convectionTerm[0] = 0;
    convectionTerm[1] = 0;
    convectionTerm[2] = 0;

    // WEST
    bNoBoundary = getNeighborOrBoundaryState(icoords,WEST,sl,m_t_int);
    if(bNoBoundary)
    {
	ICoords[0] = icoords[0] - 1;
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep(ICoords,EAST,sl);
	getFaceVelocity_middleStep(icoords,WEST,sr);
	getRiemannSolution(COORD_X,sl,sr,state_west);
    }
    else
	state_west = sl;

    // EAST
    bNoBoundary = getNeighborOrBoundaryState(icoords,EAST,sr,m_t_int);
    if(bNoBoundary)
    {
	ICoords[0] = icoords[0] + 1;
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep(ICoords,WEST,sr);
    	getFaceVelocity_middleStep(icoords,EAST,sl);
	getRiemannSolution(COORD_X,sl,sr,state_east);
    }
    else
	state_east = sr;

    // SOUTH
    bNoBoundary = getNeighborOrBoundaryState(icoords,SOUTH,sl,m_t_int);
    if(bNoBoundary)
    {
	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1] - 1;
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep(ICoords,NORTH,sl);
	getFaceVelocity_middleStep(icoords,SOUTH,sr);
	getRiemannSolution(COORD_Y,sl,sr,state_south);
    }
    else
	state_south = sl;

    // NORTH
    bNoBoundary = getNeighborOrBoundaryState(icoords,NORTH,sr,m_t_int);
    if(bNoBoundary)
    {
	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1] + 1;
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep(ICoords,SOUTH,sr);
	getFaceVelocity_middleStep(icoords,NORTH,sl);
	getRiemannSolution(COORD_Y,sl,sr,state_north);
    }
    else
	state_north = sr;

    // LOWER
    bNoBoundary = getNeighborOrBoundaryState(icoords,LOWER,sl,m_t_int);
    if(bNoBoundary)
    {
	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2] - 1;
	getFaceVelocity_middleStep(ICoords,UPPER,sl);
	getFaceVelocity_middleStep(icoords,LOWER,sr);
	getRiemannSolution(COORD_Z,sl,sr,state_lower);
    }
    else
	state_lower = sl;

    // UPPER
    bNoBoundary = getNeighborOrBoundaryState(icoords,UPPER,sr,m_t_int);
    if(bNoBoundary)
    {
	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2] + 1;
	getFaceVelocity_middleStep(ICoords,LOWER,sr);
	getFaceVelocity_middleStep(icoords,UPPER,sl);
	getRiemannSolution(COORD_Z,sl,sr,state_upper);
    }
    else
	state_upper = sr;

    convectionTerm[0] =
	    1.0/2*(state_east.m_U[0]+state_west.m_U[0])   /r * (state_east.m_U[0]-state_west.m_U[0])  /dtheta +
	    1.0/2*(state_north.m_U[1]+state_south.m_U[1])    * (state_north.m_U[0]-state_south.m_U[0])/dz +
	    1.0/2*(state_upper.m_U[2]+state_lower.m_U[2])    * (state_upper.m_U[0]-state_lower.m_U[0])/dr +
	    1/r * 1.0/2*(state_east.m_U[0]+state_west.m_U[0]) * 1.0/2*(state_upper.m_U[2]+state_lower.m_U[2]);

    convectionTerm[1] =
	    1.0/2*(state_east.m_U[0]+state_west.m_U[0])   /r * (state_east.m_U[1]-state_west.m_U[1])  /dtheta +
	    1.0/2*(state_north.m_U[1]+state_south.m_U[1])    * (state_north.m_U[1]-state_south.m_U[1])/dz +
	    1.0/2*(state_upper.m_U[2]+state_lower.m_U[2])    * (state_upper.m_U[1]-state_lower.m_U[1])/dr;

    convectionTerm[2] =
	    1.0/2*(state_east.m_U[0]+state_west.m_U[0])   /r * (state_east.m_U[2]-state_west.m_U[2])  /dtheta +
	    1.0/2*(state_north.m_U[1]+state_south.m_U[1])    * (state_north.m_U[2]-state_south.m_U[2])/dz +
	    1.0/2*(state_upper.m_U[2]+state_lower.m_U[2])    * (state_upper.m_U[2]-state_lower.m_U[2])/dr -
	    1/r * sqr(1.0/2*(state_east.m_U[0]+state_west.m_U[0]));
}


void Incompress_Solver_Smooth_3D_Cylindrical::getAdvectionTerm_coupled(
	int *icoords,
	double convectionTerm[3])
{
    bool bNoBoundary;
    int ICoords[3];
    L_STATE sl, sr, state_west, state_east, state_south, state_north, state_lower, state_upper;

    double dtheta = top_h[0];
    double dz 	  = top_h[1];
    double dr     = top_h[2];
    int index = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
    double r = cell_center[index].m_coords[2];

    convectionTerm[0] = 0;
    convectionTerm[1] = 0;
    convectionTerm[2] = 0;

    // WEST
    bNoBoundary = getNeighborOrBoundaryState(icoords,WEST,sl,m_t_int);
    if(bNoBoundary)
    {
	ICoords[0] = icoords[0] - 1;
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep_coupled(ICoords,EAST,sl);
	getFaceVelocity_middleStep_coupled(icoords,WEST,sr);
	getRiemannSolution(COORD_X,sl,sr,state_west);
    }
    else
	state_west = sl;

    // EAST
    bNoBoundary = getNeighborOrBoundaryState(icoords,EAST,sr,m_t_int);
    if(bNoBoundary)
    {
	ICoords[0] = icoords[0] + 1;
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep_coupled(ICoords,WEST,sr);
    	getFaceVelocity_middleStep_coupled(icoords,EAST,sl);
	getRiemannSolution(COORD_X,sl,sr,state_east);
    }
    else
	state_east = sr;

    // SOUTH
    bNoBoundary = getNeighborOrBoundaryState(icoords,SOUTH,sl,m_t_int);
    if(bNoBoundary)
    {
	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1] - 1;
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep_coupled(ICoords,NORTH,sl);
	getFaceVelocity_middleStep_coupled(icoords,SOUTH,sr);
	getRiemannSolution(COORD_Y,sl,sr,state_south);
    }
    else
	state_south = sl;

    // NORTH
    bNoBoundary = getNeighborOrBoundaryState(icoords,NORTH,sr,m_t_int);
    if(bNoBoundary)
    {
	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1] + 1;
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep_coupled(ICoords,SOUTH,sr);
	getFaceVelocity_middleStep_coupled(icoords,NORTH,sl);
	getRiemannSolution(COORD_Y,sl,sr,state_north);
    }
    else
	state_north = sr;

    // LOWER
    bNoBoundary = getNeighborOrBoundaryState(icoords,LOWER,sl,m_t_int);
    if(bNoBoundary)
    {
	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2] - 1;
	getFaceVelocity_middleStep_coupled(ICoords,UPPER,sl);
	getFaceVelocity_middleStep_coupled(icoords,LOWER,sr);
	getRiemannSolution(COORD_Z,sl,sr,state_lower);
    }
    else
	state_lower = sl;

    // UPPER
    bNoBoundary = getNeighborOrBoundaryState(icoords,UPPER,sr,m_t_int);
    if(bNoBoundary)
    {
	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2] + 1;
	getFaceVelocity_middleStep_coupled(ICoords,LOWER,sr);
	getFaceVelocity_middleStep_coupled(icoords,UPPER,sl);
	getRiemannSolution(COORD_Z,sl,sr,state_upper);
    }
    else
	state_upper = sr;

    convectionTerm[0] =
	    1.0/2*(state_east.m_U[0]+state_west.m_U[0])   /r * (state_east.m_U[0]-state_west.m_U[0])  /dtheta +
	    1.0/2*(state_north.m_U[1]+state_south.m_U[1])    * (state_north.m_U[0]-state_south.m_U[0])/dz +
	    1.0/2*(state_upper.m_U[2]+state_lower.m_U[2])    * (state_upper.m_U[0]-state_lower.m_U[0])/dr +
	    1/r * 1.0/2*(state_east.m_U[0]+state_west.m_U[0]) * 1.0/2*(state_upper.m_U[2]+state_lower.m_U[2]);

    convectionTerm[1] =
	    1.0/2*(state_east.m_U[0]+state_west.m_U[0])   /r * (state_east.m_U[1]-state_west.m_U[1])  /dtheta +
	    1.0/2*(state_north.m_U[1]+state_south.m_U[1])    * (state_north.m_U[1]-state_south.m_U[1])/dz +
	    1.0/2*(state_upper.m_U[2]+state_lower.m_U[2])    * (state_upper.m_U[1]-state_lower.m_U[1])/dr;

    convectionTerm[2] =
	    1.0/2*(state_east.m_U[0]+state_west.m_U[0])   /r * (state_east.m_U[2]-state_west.m_U[2])  /dtheta +
	    1.0/2*(state_north.m_U[1]+state_south.m_U[1])    * (state_north.m_U[2]-state_south.m_U[2])/dz +
	    1.0/2*(state_upper.m_U[2]+state_lower.m_U[2])    * (state_upper.m_U[2]-state_lower.m_U[2])/dr -
	    1/r * sqr(1.0/2*(state_east.m_U[0]+state_west.m_U[0]));
}


void Incompress_Solver_Smooth_3D_Cylindrical::getFaceVelocity_middleStep(
	int *icoords,
	GRID_DIRECTION dir,
	L_STATE &state_face)
{
    L_STATE state;
    int index = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
    state = cell_center[index].m_state;
    double r = cell_center[index].m_coords[2];
    double rho = cell_center[index].m_state.m_rho;

    double dx = 0.0, dy = 0.0 ,dz = 0.0, slope_x_limited[3] = {0,0,0}, slope_y_limited[3] = {0,0,0}, slope_z_limited[3] = {0,0,0};

    switch(dir)
    {
    case WEST:
	dx = -top_h[0];
	break;
    case EAST:
	dx =  top_h[0];
	break;
    case SOUTH:
	dy = -top_h[1];
	break;
    case NORTH:
	dy =  top_h[1];
	break;
    case LOWER:
	dz = -top_h[2];
	break;
    case UPPER:
	dz =  top_h[2];
	break;
    default:
	assert(false);
    }

    state_face.m_U[0] = state.m_U[0];
    state_face.m_U[1] = state.m_U[1];
    state_face.m_U[2] = state.m_U[2];

//    return;

    getLimitedSlope(icoords,COORD_X,slope_x_limited);
    getLimitedSlope(icoords,COORD_Y,slope_y_limited);
    getLimitedSlope(icoords,COORD_Z,slope_z_limited);

    // dx/2, dy/2, dz/2
    state_face.m_U[0] += dx/2 * slope_x_limited[0] + dy/2 * slope_y_limited[0] + dz/2 * slope_z_limited[0];
    state_face.m_U[1] += dx/2 * slope_x_limited[1] + dy/2 * slope_y_limited[1] + dz/2 * slope_z_limited[1];
    state_face.m_U[2] += dx/2 * slope_x_limited[2] + dy/2 * slope_y_limited[2] + dz/2 * slope_z_limited[2];

    //    return;
    // dt/2
    double diffusion[3];
    getDifffusion(icoords,diffusion);
    state_face.m_U[0] += m_dt/2 *
	    ( diffusion[0]/rho - (1/r*state.m_U[0]*slope_x_limited[0] + state.m_U[1]*slope_y_limited[0] + state.m_U[2]*slope_z_limited[0] + 1/r*(state.m_U[0]*state.m_U[2]) ) - cell_center[index].m_state.grad_q[0]/rho );
    state_face.m_U[1] += m_dt/2 *
	    ( diffusion[1]/rho - (1/r*state.m_U[0]*slope_x_limited[1] + state.m_U[1]*slope_y_limited[1] + state.m_U[2]*slope_z_limited[1]                                   ) - cell_center[index].m_state.grad_q[1]/rho );
    state_face.m_U[2] += m_dt/2 *
	    ( diffusion[2]/rho - (1/r*state.m_U[0]*slope_x_limited[2] + state.m_U[1]*slope_y_limited[2] + state.m_U[2]*slope_z_limited[2] - 1/r*(state.m_U[0]*state.m_U[0]) ) - cell_center[index].m_state.grad_q[2]/rho );

    //    return;
    // rhs
    double coords[3];
    L_STATE source_term;

    getRectangleCenter(index, coords);
    // the source term does not contain density.
    computeSourceTerm(coords, source_term);

    state_face.m_U[0] += m_dt/2 * source_term.m_U[0];
    state_face.m_U[1] += m_dt/2 * source_term.m_U[1];
    state_face.m_U[2] += m_dt/2 * source_term.m_U[2];
}

void Incompress_Solver_Smooth_3D_Cylindrical::getFaceVelocity_middleStep_coupled(
	int *icoords,
	GRID_DIRECTION dir,
	L_STATE &state_face)
{
    L_STATE state;
    int index = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
    state = cell_center[index].m_state;
    double r = cell_center[index].m_coords[2];
    double rho = cell_center[index].m_state.m_rho;

    double dx = 0.0, dy = 0.0 ,dz = 0.0, slope_x_limited[3] = {0,0,0}, slope_y_limited[3] = {0,0,0}, slope_z_limited[3] = {0,0,0};

    switch(dir)
    {
    case WEST:
	dx = -top_h[0];
	break;
    case EAST:
	dx =  top_h[0];
	break;
    case SOUTH:
	dy = -top_h[1];
	break;
    case NORTH:
	dy =  top_h[1];
	break;
    case LOWER:
	dz = -top_h[2];
	break;
    case UPPER:
	dz =  top_h[2];
	break;
    default:
	assert(false);
    }

    state_face.m_U[0] = state.m_U[0];
    state_face.m_U[1] = state.m_U[1];
    state_face.m_U[2] = state.m_U[2];

//    return;

    getLimitedSlope(icoords,COORD_X,slope_x_limited);
    getLimitedSlope(icoords,COORD_Y,slope_y_limited);
    getLimitedSlope(icoords,COORD_Z,slope_z_limited);

    // dx/2, dy/2, dz/2
    state_face.m_U[0] += dx/2 * slope_x_limited[0] + dy/2 * slope_y_limited[0] + dz/2 * slope_z_limited[0];
    state_face.m_U[1] += dx/2 * slope_x_limited[1] + dy/2 * slope_y_limited[1] + dz/2 * slope_z_limited[1];
    state_face.m_U[2] += dx/2 * slope_x_limited[2] + dy/2 * slope_y_limited[2] + dz/2 * slope_z_limited[2];

    //    return;
    // dt/2
    double diffusion[3];
    getDiffusion_coupled(icoords,diffusion);
    state_face.m_U[0] += m_dt/2 *
	    ( diffusion[0]/rho - (1/r*state.m_U[0]*slope_x_limited[0] + state.m_U[1]*slope_y_limited[0] + state.m_U[2]*slope_z_limited[0] + 1/r*(state.m_U[0]*state.m_U[2]) ) - cell_center[index].m_state.grad_q[0]/rho );
    state_face.m_U[1] += m_dt/2 *
	    ( diffusion[1]/rho - (1/r*state.m_U[0]*slope_x_limited[1] + state.m_U[1]*slope_y_limited[1] + state.m_U[2]*slope_z_limited[1]                                   ) - cell_center[index].m_state.grad_q[1]/rho );
    state_face.m_U[2] += m_dt/2 *
	    ( diffusion[2]/rho - (1/r*state.m_U[0]*slope_x_limited[2] + state.m_U[1]*slope_y_limited[2] + state.m_U[2]*slope_z_limited[2] - 1/r*(state.m_U[0]*state.m_U[0]) ) - cell_center[index].m_state.grad_q[2]/rho );

    //    return;
    // rhs
    double coords[3];
    L_STATE source_term;

    getRectangleCenter(index, coords);
    // the source term does not contain density.
    computeSourceTerm(coords, source_term);

    state_face.m_U[0] += m_dt/2 * source_term.m_U[0];
    state_face.m_U[1] += m_dt/2 * source_term.m_U[1];
    state_face.m_U[2] += m_dt/2 * source_term.m_U[2];
}


/**
* compute mu * (Uxx+Uyy+Uzz)
* @param icoords
* @param diffusion
* @param gradP
*/

void Incompress_Solver_Smooth_3D_Cylindrical::getDifffusion(
	int *icoords,
	double diffusion[3])
{
    double dU2_theta[3], dU2_z[3], dU2_r[3];
    double dU_theta[3], dU_z[3], dU_r[3];
    getDU2(icoords,COORD_X, dU2_theta, dU_theta);
    getDU2(icoords,COORD_Y, dU2_z,     dU_z);
    getDU2(icoords,COORD_Z, dU2_r,     dU_r);

    int index = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
    double mu = cell_center[index].m_state.m_mu;

    double r = cell_center[index].m_coords[2];
    L_STATE U = cell_center[index].m_state;

    diffusion[0] = mu * (1/(r*r)*dU2_theta[0] + 1/r*dU_r[0] + dU2_z[0] + dU2_r[0]);
    diffusion[0]+= mu * ( 2/(r*r)*dU_theta[2] - 1/(r*r)*U.m_U[0]);

    diffusion[1] = mu * (1/(r*r)*dU2_theta[1] + 1/r*dU_r[1] + dU2_z[1] + dU2_r[1]);

    diffusion[2] = mu * (1/(r*r)*dU2_theta[2] + 1/r*dU_r[2] + dU2_z[2] + dU2_r[2]);
    diffusion[2]+= mu * (-2/(r*r)*dU_theta[0] - 1/(r*r)*U.m_U[2]);
}

void Incompress_Solver_Smooth_3D_Cylindrical::getDiffusion_coupled(
	int *icoords,
	double diffusion[3])
{
    int index, index_nb[18];
    double mu[6], mu_edge[6], mu0;
    double r0, rr, r[6], r_edge[6];
    L_STATE Unb;
    double U_nb[3][18], U_center[3], U_face[3][6];
    int nb;
    GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
    bool bNoBoundary[6];
    double dh[3],dh0[3],dh1[3];
    double dtheta,dz,dr;

    dtheta = top_h[0];
    dz = top_h[1];
    dr = top_h[2];

    diffusion[0] = 0.0;
    diffusion[1] = 0.0;
    diffusion[2] = 0.0;

    int i = icoords[0];
    int j = icoords[1];
    int k = icoords[2];

    index = d_index3d(i,j,k,top_gmax);

    //6 neighbours of the center cell
    index_nb[0] = d_index3d(i-1,j,k,top_gmax);
    index_nb[1] = d_index3d(i+1,j,k,top_gmax);
    index_nb[2] = d_index3d(i,j-1,k,top_gmax);
    index_nb[3] = d_index3d(i,j+1,k,top_gmax);
    index_nb[4] = d_index3d(i,j,k-1,top_gmax);
    index_nb[5] = d_index3d(i,j,k+1,top_gmax);
    
    //theta-z cut neighbours
    index_nb[6] = d_index3d(i-1,j-1,k,top_gmax);
    index_nb[7] = d_index3d(i+1,j-1,k,top_gmax);
    index_nb[8] = d_index3d(i+1,j+1,k,top_gmax);
    index_nb[9] = d_index3d(i-1,j+1,k,top_gmax);
	
    //z-r cut neighbours
    index_nb[10] = d_index3d(i,j-1,k-1,top_gmax);
    index_nb[11] = d_index3d(i,j+1,k-1,top_gmax);
    index_nb[12] = d_index3d(i,j+1,k+1,top_gmax);
    index_nb[13] = d_index3d(i,j-1,k+1,top_gmax);

    //theta-r cut neighbours
    index_nb[14] = d_index3d(i-1,j,k-1,top_gmax);
    index_nb[15] = d_index3d(i+1,j,k-1,top_gmax);
    index_nb[16] = d_index3d(i+1,j,k+1,top_gmax);
    index_nb[17] = d_index3d(i-1,j,k+1,top_gmax);

    mu0 = cell_center[index].m_state.m_mu;
    r0 = cell_center[index].m_coords[2];
    rr = r0*r0;

    for (int l = 0; l < 3; l++)
	U_center[l] = cell_center[index].m_state.m_U[l];

    for (nb = 0; nb < 6; nb++)
    {
	bNoBoundary[nb] = getNeighborOrBoundaryState(icoords, dir[nb], Unb, m_t_old);
	U_nb[0][nb] = Unb.m_U[0];
	U_nb[1][nb] = Unb.m_U[1];
	U_nb[2][nb] = Unb.m_U[2];

	if(!bNoBoundary[nb])
	{
	    mu[nb] = mu0;
	    mu_edge[nb] = mu0;
	    U_face[0][nb] = U_nb[0][nb];
	    U_face[1][nb] = U_nb[1][nb];
	    U_face[2][nb] = U_nb[2][nb];
	}
	else
	{
	    mu[nb] = cell_center[index_nb[nb]].m_state.m_mu;
	    mu_edge[nb] = 0.5*(mu[nb] + mu0);
	    U_face[0][nb] = 0.5*(U_nb[0][nb] + U_center[0]);
	    U_face[1][nb] = 0.5*(U_nb[1][nb] + U_center[1]);
	    U_face[2][nb] = 0.5*(U_nb[2][nb] + U_center[2]);
	}
	r[nb] = cell_center[index_nb[nb]].m_coords[2];
	r_edge[nb] = 0.5 * (r[nb] + r0);
    }

    //traverse the corners on 3 cut planes

    //corner (i-1/2,j-1/2,k)

    if (!bNoBoundary[0] && bNoBoundary[2])
    {
	U_nb[0][6] = U_nb[0][0];
	U_nb[1][6] = U_nb[1][0];
	U_nb[2][6] = U_nb[2][0];
    }
    else if(bNoBoundary[0] && !bNoBoundary[2])
    {
	U_nb[0][6] = U_nb[0][2];
	U_nb[1][6] = U_nb[1][2];
	U_nb[2][6] = U_nb[2][2];
    }
    else if(!bNoBoundary[0] && !bNoBoundary[2])
    {
	U_nb[0][6] = U_nb[0][0];
	U_nb[1][6] = U_nb[1][0];
	U_nb[2][6] = U_nb[2][0];
    }
    else
    {
	U_nb[0][6] = (U_nb[0][0]+U_nb[0][2]+cell_center[index_nb[6]].m_state.m_U[0]+U_center[0])/4.0;
	U_nb[1][6] = (U_nb[1][0]+U_nb[1][2]+cell_center[index_nb[6]].m_state.m_U[1]+U_center[1])/4.0;
	U_nb[2][6] = (U_nb[2][0]+U_nb[2][2]+cell_center[index_nb[6]].m_state.m_U[2]+U_center[2])/4.0;
    }

    //corner (i+1/2,j-1/2,k)

    if (!bNoBoundary[1] && bNoBoundary[2])
    {
	U_nb[0][7] = U_nb[0][1];
	U_nb[1][7] = U_nb[1][1];
	U_nb[2][7] = U_nb[2][1];
    }
    else if(bNoBoundary[1] && !bNoBoundary[2])
    {
	U_nb[0][7] = U_nb[0][2];
	U_nb[1][7] = U_nb[1][2];
	U_nb[2][7] = U_nb[2][2];
    }
    else if(!bNoBoundary[1] && !bNoBoundary[2])
    {
	U_nb[0][7] = U_nb[0][1];
	U_nb[1][7] = U_nb[1][1];
	U_nb[2][7] = U_nb[2][1];
    }
    else
    {
	U_nb[0][7] = (U_nb[0][1]+U_nb[0][2]+cell_center[index_nb[7]].m_state.m_U[0]+U_center[0])/4.0;
	U_nb[1][7] = (U_nb[1][1]+U_nb[1][2]+cell_center[index_nb[7]].m_state.m_U[1]+U_center[1])/4.0;
	U_nb[2][7] = (U_nb[2][1]+U_nb[2][2]+cell_center[index_nb[7]].m_state.m_U[2]+U_center[2])/4.0;
    }


    //corner (i+1/2,j+1/2,k)

    if (!bNoBoundary[1] && bNoBoundary[3])
    {
	U_nb[0][8] = U_nb[0][1];
	U_nb[1][8] = U_nb[1][1];
	U_nb[2][8] = U_nb[2][1];
    }
    else if(bNoBoundary[1] && !bNoBoundary[3])
    {
	U_nb[0][8] = U_nb[0][3];
	U_nb[1][8] = U_nb[1][3];
	U_nb[2][8] = U_nb[2][3];
    }
    else if(!bNoBoundary[1] && !bNoBoundary[3])
    {
	U_nb[0][8] = U_nb[0][1];
	U_nb[1][8] = U_nb[1][1];
	U_nb[2][8] = U_nb[2][1];
    }
    else
    {
	U_nb[0][8] = (U_nb[0][1]+U_nb[0][3]+cell_center[index_nb[8]].m_state.m_U[0]+U_center[0])/4.0;
	U_nb[1][8] = (U_nb[1][1]+U_nb[1][3]+cell_center[index_nb[8]].m_state.m_U[1]+U_center[1])/4.0;
	U_nb[2][8] = (U_nb[2][1]+U_nb[2][3]+cell_center[index_nb[8]].m_state.m_U[2]+U_center[2])/4.0;
    }

    //corner (i-1/2,j+1/2,k)

    if (!bNoBoundary[0] && bNoBoundary[3])
    {
	U_nb[0][9] = U_nb[0][0];
	U_nb[1][9] = U_nb[1][0];
	U_nb[2][9] = U_nb[2][0];
    }
    else if(bNoBoundary[0] && !bNoBoundary[3])
    {
	U_nb[0][9] = U_nb[0][3];
	U_nb[1][9] = U_nb[1][3];
	U_nb[2][9] = U_nb[2][3];
    }
    else if(!bNoBoundary[0] && !bNoBoundary[3])
    {
	U_nb[0][9] = U_nb[0][0];
	U_nb[1][9] = U_nb[1][0];
	U_nb[2][9] = U_nb[2][0];
    }
    else
    {
	U_nb[0][9] = (U_nb[0][0]+U_nb[0][3]+cell_center[index_nb[9]].m_state.m_U[0]+U_center[0])/4.0;
	U_nb[1][9] = (U_nb[1][0]+U_nb[1][3]+cell_center[index_nb[9]].m_state.m_U[1]+U_center[1])/4.0;
	U_nb[2][9] = (U_nb[2][0]+U_nb[2][3]+cell_center[index_nb[9]].m_state.m_U[2]+U_center[2])/4.0;
    }

    //corner (i,j-1/2,k-1/2)

    if (!bNoBoundary[2] && bNoBoundary[4])
    {
	U_nb[0][10] = U_nb[0][2];
	U_nb[1][10] = U_nb[1][2];
	U_nb[2][10] = U_nb[2][2];
    }
    else if(bNoBoundary[2] && !bNoBoundary[4])
    {
	U_nb[0][10] = U_nb[0][4];
	U_nb[1][10] = U_nb[1][4];
	U_nb[2][10] = U_nb[2][4];
    }
    else if(!bNoBoundary[2] && !bNoBoundary[4])
    {
	U_nb[0][10] = U_nb[0][2];
	U_nb[1][10] = U_nb[1][2];
	U_nb[2][10] = U_nb[2][2];
    }
    else
    {
	U_nb[0][10] = (U_nb[0][2]+U_nb[0][4]+cell_center[index_nb[10]].m_state.m_U[0]+U_center[0])/4.0;
	U_nb[1][10] = (U_nb[1][2]+U_nb[1][4]+cell_center[index_nb[10]].m_state.m_U[1]+U_center[1])/4.0;
	U_nb[2][10] = (U_nb[2][2]+U_nb[2][4]+cell_center[index_nb[10]].m_state.m_U[2]+U_center[2])/4.0;
    }

    //corner (i,j+1/2,k-1/2)

    if (!bNoBoundary[3] && bNoBoundary[4])
    {
	U_nb[0][11] = U_nb[0][3];
	U_nb[1][11] = U_nb[1][3];
	U_nb[2][11] = U_nb[2][3];
    }
    else if(bNoBoundary[3] && !bNoBoundary[4])
    {
	U_nb[0][11] = U_nb[0][4];
	U_nb[1][11] = U_nb[1][4];
	U_nb[2][11] = U_nb[2][4];
    }
    else if(!bNoBoundary[3] && !bNoBoundary[4])
    {
	U_nb[0][11] = U_nb[0][3];
	U_nb[1][11] = U_nb[1][3];
	U_nb[2][11] = U_nb[2][3];
    }
    else
    {
	U_nb[0][11] = (U_nb[0][3]+U_nb[0][4]+cell_center[index_nb[11]].m_state.m_U[0]+U_center[0])/4.0;
	U_nb[1][11] = (U_nb[1][3]+U_nb[1][4]+cell_center[index_nb[11]].m_state.m_U[1]+U_center[1])/4.0;
	U_nb[2][11] = (U_nb[2][3]+U_nb[2][4]+cell_center[index_nb[11]].m_state.m_U[2]+U_center[2])/4.0;
    }

    //corner (i,j+1/2,k+1/2)

    if (!bNoBoundary[3] && bNoBoundary[5])
    {
	U_nb[0][12] = U_nb[0][3];
	U_nb[1][12] = U_nb[1][3];
	U_nb[2][12] = U_nb[2][3];
    }
    else if(bNoBoundary[3] && !bNoBoundary[5])
    {
	U_nb[0][12] = U_nb[0][5];
	U_nb[1][12] = U_nb[1][5];
	U_nb[2][12] = U_nb[2][5];
    }
    else if(!bNoBoundary[3] && !bNoBoundary[5])
    {
	U_nb[0][12] = U_nb[0][3];
	U_nb[1][12] = U_nb[1][3];
	U_nb[2][12] = U_nb[2][3];
    }
    else
    {
	U_nb[0][12] = (U_nb[0][3]+U_nb[0][5]+cell_center[index_nb[12]].m_state.m_U[0]+U_center[0])/4.0;
	U_nb[1][12] = (U_nb[1][3]+U_nb[1][5]+cell_center[index_nb[12]].m_state.m_U[1]+U_center[1])/4.0;
	U_nb[2][12] = (U_nb[2][3]+U_nb[2][5]+cell_center[index_nb[12]].m_state.m_U[2]+U_center[2])/4.0;
    }

    //corner (i,j-1/2,k+1/2)

    if (!bNoBoundary[2] && bNoBoundary[5])
    {
	U_nb[0][13] = U_nb[0][2];
	U_nb[1][13] = U_nb[1][2];
	U_nb[2][13] = U_nb[2][2];
    }
    else if(bNoBoundary[2] && !bNoBoundary[5])
    {
	U_nb[0][13] = U_nb[0][5];
	U_nb[1][13] = U_nb[1][5];
	U_nb[2][13] = U_nb[2][5];
    }
    else if(!bNoBoundary[2] && !bNoBoundary[5])
    {
	U_nb[0][13] = U_nb[0][2];
	U_nb[1][13] = U_nb[1][2];
	U_nb[2][13] = U_nb[2][2];
    }
    else
    {
	U_nb[0][13] = (U_nb[0][2]+U_nb[0][5]+cell_center[index_nb[13]].m_state.m_U[0]+U_center[0])/4.0;
	U_nb[1][13] = (U_nb[1][2]+U_nb[1][5]+cell_center[index_nb[13]].m_state.m_U[1]+U_center[1])/4.0;
	U_nb[2][13] = (U_nb[2][2]+U_nb[2][5]+cell_center[index_nb[13]].m_state.m_U[2]+U_center[2])/4.0;
    }


    //corner (i-1/2,j,k-1/2)

    if (!bNoBoundary[0] && bNoBoundary[4])
    {
	U_nb[0][14] = U_nb[0][0];
	U_nb[1][14] = U_nb[1][0];
	U_nb[2][14] = U_nb[2][0];
    }
    else if(bNoBoundary[0] && !bNoBoundary[4])
    {
	U_nb[0][14] = U_nb[0][4];
	U_nb[1][14] = U_nb[1][4];
	U_nb[2][14] = U_nb[2][4];
    }
    else if(!bNoBoundary[0] && !bNoBoundary[4])
    {
	U_nb[0][14] = U_nb[0][0];
	U_nb[1][14] = U_nb[1][0];
	U_nb[2][14] = U_nb[2][0];
    }
    else
    {
	U_nb[0][14] = (U_nb[0][0]+U_nb[0][4]+cell_center[index_nb[14]].m_state.m_U[0]+U_center[0])/4.0;
	U_nb[1][14] = (U_nb[1][0]+U_nb[1][4]+cell_center[index_nb[14]].m_state.m_U[1]+U_center[1])/4.0;
	U_nb[2][14] = (U_nb[2][0]+U_nb[2][4]+cell_center[index_nb[14]].m_state.m_U[2]+U_center[2])/4.0;
    }

    //corner (i+1/2,j,k-1/2)

    if (!bNoBoundary[1] && bNoBoundary[4])
    {
	U_nb[0][15] = U_nb[0][1];
	U_nb[1][15] = U_nb[1][1];
	U_nb[2][15] = U_nb[2][1];
    }
    else if(bNoBoundary[1] && !bNoBoundary[4])
    {
	U_nb[0][15] = U_nb[0][4];
	U_nb[1][15] = U_nb[1][4];
	U_nb[2][15] = U_nb[2][4];
    }
    else if(!bNoBoundary[1] && !bNoBoundary[4])
    {
	U_nb[0][15] = U_nb[0][1];
	U_nb[1][15] = U_nb[1][1];
	U_nb[2][15] = U_nb[2][1];
    }
    else
    {
	U_nb[0][15] = (U_nb[0][1]+U_nb[0][4]+cell_center[index_nb[15]].m_state.m_U[0]+U_center[0])/4.0;
	U_nb[1][15] = (U_nb[1][1]+U_nb[1][4]+cell_center[index_nb[15]].m_state.m_U[1]+U_center[1])/4.0;
	U_nb[2][15] = (U_nb[2][1]+U_nb[2][4]+cell_center[index_nb[15]].m_state.m_U[2]+U_center[2])/4.0;
    }

    //corner (i+1/2,j,k+1/2)

    if (!bNoBoundary[1] && bNoBoundary[5])
    {
	U_nb[0][16] = U_nb[0][1];
	U_nb[1][16] = U_nb[1][1];
	U_nb[2][16] = U_nb[2][1];
    }
    else if(bNoBoundary[1] && !bNoBoundary[5])
    {
	U_nb[0][16] = U_nb[0][5];
	U_nb[1][16] = U_nb[1][5];
	U_nb[2][16] = U_nb[2][5];
    }
    else if(!bNoBoundary[1] && !bNoBoundary[5])
    {
	U_nb[0][16] = U_nb[0][1];
	U_nb[1][16] = U_nb[1][1];
	U_nb[2][16] = U_nb[2][1];
    }
    else
    {
	U_nb[0][16] = (U_nb[0][1]+U_nb[0][5]+cell_center[index_nb[16]].m_state.m_U[0]+U_center[0])/4.0;
	U_nb[1][16] = (U_nb[1][1]+U_nb[1][5]+cell_center[index_nb[16]].m_state.m_U[1]+U_center[1])/4.0;
	U_nb[2][16] = (U_nb[2][1]+U_nb[2][5]+cell_center[index_nb[16]].m_state.m_U[2]+U_center[2])/4.0;
    }

    //corner (i-1/2,j,k+1/2)

    if (!bNoBoundary[0] && bNoBoundary[5])
    {
	U_nb[0][17] = U_nb[0][0];
	U_nb[1][17] = U_nb[1][0];
	U_nb[2][17] = U_nb[2][0];
    }
    else if(bNoBoundary[0] && !bNoBoundary[5])
    {
	U_nb[0][17] = U_nb[0][5];
	U_nb[1][17] = U_nb[1][5];
	U_nb[2][17] = U_nb[2][5];
    }
    else if(!bNoBoundary[0] && !bNoBoundary[5])
    {
	U_nb[0][17] = U_nb[0][0];
	U_nb[1][17] = U_nb[1][0];
	U_nb[2][17] = U_nb[2][0];
    }
    else
    {
	U_nb[0][17] = (U_nb[0][0]+U_nb[0][5]+cell_center[index_nb[17]].m_state.m_U[0]+U_center[0])/4.0;
	U_nb[1][17] = (U_nb[1][0]+U_nb[1][5]+cell_center[index_nb[17]].m_state.m_U[1]+U_center[1])/4.0;
	U_nb[2][17] = (U_nb[2][0]+U_nb[2][5]+cell_center[index_nb[17]].m_state.m_U[2]+U_center[2])/4.0;
    }

    for (int l = 0; l < 3; l++)
    {
	dh[l] = top_h[l];

	if (bNoBoundary[l*2])
	    dh0[l] = top_h[l];
	else
	    dh0[l] = top_h[l]/2.0;

	if (bNoBoundary[l*2+1])
	    dh1[l] = top_h[l];
	else
	    dh1[l] = top_h[l]/2.0;
    }


    /***************** Diffusion term for Equation 1  ***************/

    ////////////////   Tensor term 1 //////////////
    //first term

    diffusion[0] += ( 
	             mu_edge[3]*(U_nb[0][3]-U_center[0])/dh1[1] 
	            -mu_edge[2]*(U_center[0]-U_nb[0][2])/dh0[1]
		    ) / dh[1];
    //second term

    diffusion[0] += (mu_edge[2]/r_edge[2]*U_nb[1][6] - mu_edge[2]/r_edge[2]*U_nb[1][7]
	            +mu_edge[3]/r_edge[3]*U_nb[1][8] - mu_edge[3]/r_edge[3]*U_nb[1][9])/(dz*dtheta);

    //////////// Tensor term 2 /////////////
    //first term

    diffusion[0] += ( 
	             mu_edge[5]*(U_nb[0][5]-U_center[0])/dh1[2] 
	            -mu_edge[4]*(U_center[0]-U_nb[0][4])/dh0[2]
		    ) / dh[2];

    //second term

    diffusion[0] += (mu_edge[4]/r_edge[4]*U_face[0][4] - mu_edge[5]/r_edge[5]*U_face[0][5]) / dr;

    //third term

    diffusion[0] += (mu_edge[4]/r_edge[4]*U_nb[2][14] - mu_edge[4]/r_edge[4]*U_nb[2][15]
	            +mu_edge[5]/r_edge[5]*U_nb[2][16] - mu_edge[5]/r_edge[5]*U_nb[2][17])/(dr*dtheta);

    ///////////////  Tensor term 3  /////////////////
    //first term

    diffusion[0] += ( 
	             2.0/r0*mu_edge[1]/r_edge[1]*(U_nb[0][1]-U_center[0])/dh1[0] 
	            -2.0/r0*mu_edge[0]/r_edge[0]*(U_center[0]-U_nb[0][0])/dh0[0]
		    ) / dh[0];

    //second term

    diffusion[0] += (-2.0/r0*mu_edge[0]/r_edge[0]*U_face[2][0] 
	             +2.0/r0*mu_edge[1]/r_edge[1]*U_face[2][1]) / dtheta;

    /////////////////  Tensor term 4     //////////////
    //first term

    diffusion[0] += (-2.0*mu0/r0*U_face[0][4] 
	             +2.0*mu0/r0*U_face[0][5]) / dr;

    //second term

    diffusion[0] += (-2.0*mu0/rr*U_center[0]);

    //third term
    
    diffusion[0] += (-2.0*mu0/rr*U_face[2][0] 
	             +2.0*mu0/rr*U_face[2][1]) / dtheta;

    /******************* Diffusion term for Equation 2 ************/

    //////////////  Tensor term 1  /////////////
    //first term

    diffusion[1] += ( 
	             2.0*mu_edge[3]*(U_nb[1][3]-U_center[1])/dh1[1] 
	            -2.0*mu_edge[2]*(U_center[1]-U_nb[1][2])/dh0[1]
		    ) / dh[1];

    //////////////  Tensor term 2 /////////////////
    //first term

    diffusion[1] += (mu_edge[4]*U_nb[2][10] - mu_edge[4]*U_nb[2][11]
	            +mu_edge[5]*U_nb[2][12] - mu_edge[5]*U_nb[2][13])/(dr*dz);

    //second term

    diffusion[1] += ( 
	             mu_edge[5]*(U_nb[1][5]-U_center[1])/dh1[2] 
	            -mu_edge[4]*(U_center[1]-U_nb[1][4])/dh0[2]
		    ) / dh[2];

    /////////////// Tensor term 3 /////////////////
    //first term

    diffusion[1] += (mu_edge[0]/r0*U_nb[0][6] - mu_edge[1]/r0*U_nb[0][7]
	            +mu_edge[1]/r0*U_nb[0][8] - mu_edge[0]/r0*U_nb[0][9])/(dtheta*dz);

    //second term

    diffusion[1] += ( 
	             mu_edge[1]/r_edge[1]/r0*(U_nb[1][1]-U_center[1])/dh1[0] 
	            -mu_edge[0]/r_edge[0]/r0*(U_center[1]-U_nb[1][0])/dh0[0]
		    ) / dh[0];

    //////////////// Tensor term 4  ///////////////////////
    //first term
 
    diffusion[1] += (-mu0/r0*U_face[2][2] 
	             +mu0/r0*U_face[2][3]) / dz;

    //second term
 
    diffusion[1] += (-mu0/r0*U_face[1][4] 
	             +mu0/r0*U_face[1][5]) / dr;


    /**************** Diffusion term for Equation 3  *****************/

    ///////////// Tensor term 1  ////////////////
    //first term

    diffusion[2] += ( 
	             mu_edge[3]*(U_nb[2][3]-U_center[2])/dh1[1] 
	            -mu_edge[2]*(U_center[2]-U_nb[2][2])/dh0[1]
		    ) / dh[1];

    //second term

    diffusion[2] += (mu_edge[2]*U_nb[1][10] - mu_edge[3]*U_nb[1][11]
	            +mu_edge[3]*U_nb[1][12] - mu_edge[2]*U_nb[1][13])/(dr*dz);

    ////////////// Tensor term 2  ////////////////
    //first term

    diffusion[2] += ( 
	             2.0*mu_edge[5]*(U_nb[2][5]-U_center[2])/dh1[2] 
	            -2.0*mu_edge[4]*(U_center[2]-U_nb[2][4])/dh0[2]
		    ) / dh[2];


    ////////////// Tensor term 3  //////////////////
    //first term

    diffusion[2] += (mu_edge[0]/r0*U_nb[0][14] - mu_edge[1]/r0*U_nb[0][15]
	            +mu_edge[1]/r0*U_nb[0][16] - mu_edge[0]/r0*U_nb[0][17])/(dr*dtheta);

    //second term

    diffusion[2] += (mu_edge[0]/r_edge[0]/r0*U_face[0][0] 
	            -mu_edge[1]/r_edge[1]/r0*U_face[0][1]) / dtheta;

    //third term

    diffusion[2] += ( 
	             mu_edge[1]/r_edge[1]/r0*(U_nb[2][1]-U_center[2])/dh1[0] 
	            -mu_edge[0]/r_edge[0]/r0*(U_center[2]-U_nb[2][0])/dh0[0]
		    ) / dh[0];

    ////////// Tensor term 4 ///////////
    //first term

    diffusion[2] += (-2.0*mu0/r0*U_face[2][4] 
	             +2.0*mu0/r0*U_face[2][5]) / dr;

    ///////////// Tensor term 5  /////////////////
    //first term

    diffusion[2] += (2.0*mu0/rr*U_face[0][0] 
	            -2.0*mu0/rr*U_face[0][1]) / dtheta;

    //second term

    diffusion[2] += (-2.0*mu0/rr*U_center[2]);

}

/**
* TODO: the calculation near the boundary might not be correct.
* @param dir
* @param icoords
* @param dU2
* @param dP
*/
void Incompress_Solver_Smooth_3D_Cylindrical::getDU2(
	int *icoords,
	EBM_COORD xyz,
	double dU2[3],
	double dU[3])
{
    double dh0, dh1, dh;
    L_STATE U0, U1, U2;

    U1 = cell_center[d_index3d(icoords[0],icoords[1],icoords[2],top_gmax)].m_state;

    bool bNoBoundary[2];

    if(xyz==COORD_X)
    {
	bNoBoundary[0] = getNeighborOrBoundaryState(icoords,WEST,U0,m_t_old);
	if(bNoBoundary[0])
	    dh0 = top_h[0];
	else
	    dh0 = top_h[0]/2;

	bNoBoundary[1] = getNeighborOrBoundaryState(icoords,EAST,U2,m_t_old);
	if(bNoBoundary[1])
	    dh1 = top_h[0];
	else
	    dh1 = top_h[0]/2;

	dh = top_h[0];
    }
    else if(xyz==COORD_Y)	//
    {
	bNoBoundary[0] = getNeighborOrBoundaryState(icoords,SOUTH,U0,m_t_old);
	if(bNoBoundary[0])
	    dh0 = top_h[1];
	else
	    dh0 = top_h[1]/2;

	bNoBoundary[1] = getNeighborOrBoundaryState(icoords,NORTH,U2,m_t_old);
	if(bNoBoundary[1])
	    dh1 = top_h[1];
	else
	    dh1 = top_h[1]/2;

	dh = top_h[1];
    }
    else
    {
	bNoBoundary[0] = getNeighborOrBoundaryState(icoords,LOWER,U0,m_t_old);
	if(bNoBoundary[0])
	    dh0 = top_h[2];
	else
	    dh0 = top_h[2]/2;

	bNoBoundary[1] = getNeighborOrBoundaryState(icoords,UPPER,U2,m_t_old);
	if(bNoBoundary[1])
	    dh1 = top_h[2];
	else
	    dh1 = top_h[2]/2;

	dh = top_h[2];
    }

    // second order derivative
    dU2[0] = ((U2.m_U[0] - U1.m_U[0])/dh1 - (U1.m_U[0] - U0.m_U[0])/dh0) / dh;
    dU2[1] = ((U2.m_U[1] - U1.m_U[1])/dh1 - (U1.m_U[1] - U0.m_U[1])/dh0) / dh;
    dU2[2] = ((U2.m_U[2] - U1.m_U[2])/dh1 - (U1.m_U[2] - U0.m_U[2])/dh0) / dh;

    // first order derivative
    dU[0] = ((U2.m_U[0] - U1.m_U[0])/dh1 + (U1.m_U[0] - U0.m_U[0])/dh0) / 2;
    dU[1] = ((U2.m_U[1] - U1.m_U[1])/dh1 + (U1.m_U[1] - U0.m_U[1])/dh0) / 2;
    dU[2] = ((U2.m_U[2] - U1.m_U[2])/dh1 + (U1.m_U[2] - U0.m_U[2])/dh0) / 2;

}


void Incompress_Solver_Smooth_3D_Cylindrical::getLimitedSlope(
	int *icoords,
	EBM_COORD xyz,
	double slope[3])
{
    double dh0, dh1;
    L_STATE U0, U1, U2;

    U1 = cell_center[d_index3d(icoords[0],icoords[1],icoords[2],top_gmax)].m_state;

    bool bNoBoundary[2];

    if(xyz==COORD_X)
    {
	bNoBoundary[0] = getNeighborOrBoundaryState(icoords,WEST,U0,m_t_old);
	if(bNoBoundary[0])
	    dh0 = top_h[0];
	else
	    dh0 = top_h[0]/2;

	bNoBoundary[1] = getNeighborOrBoundaryState(icoords,EAST,U2,m_t_old);
	if(bNoBoundary[1])
	    dh1 = top_h[0];
	else
	    dh1 = top_h[0]/2;
    }
    else if(xyz==COORD_Y)	//
    {
	bNoBoundary[0] = getNeighborOrBoundaryState(icoords,SOUTH,U0,m_t_old);
	if(bNoBoundary[0])
	    dh0 = top_h[1];
	else
	    dh0 = top_h[1]/2;

	bNoBoundary[1] = getNeighborOrBoundaryState(icoords,NORTH,U2,m_t_old);
	if(bNoBoundary[1])
	    dh1 = top_h[1];
	else
	    dh1 = top_h[1]/2;
    }
    else
    {
	bNoBoundary[0] = getNeighborOrBoundaryState(icoords,LOWER,U0,m_t_old);
	if(bNoBoundary[0])
	    dh0 = top_h[1];
	else
	    dh0 = top_h[1]/2;

	bNoBoundary[1] = getNeighborOrBoundaryState(icoords,UPPER,U2,m_t_old);
	if(bNoBoundary[1])
	    dh1 = top_h[1];
	else
	    dh1 = top_h[1]/2;
    }
    slope[0] = EBM_minmod((U1.m_U[0]-U0.m_U[0])/dh0, (U2.m_U[0]-U1.m_U[0])/dh1);
    slope[1] = EBM_minmod((U1.m_U[1]-U0.m_U[1])/dh0, (U2.m_U[1]-U1.m_U[1])/dh1);
    slope[2] = EBM_minmod((U1.m_U[2]-U0.m_U[2])/dh0, (U2.m_U[2]-U1.m_U[2])/dh1);
}


double Incompress_Solver_Smooth_3D_Cylindrical::EBM_minmod(
	double x,
	double y)
{

    double sign = x*y;

    if(sign<0)
	return 0;
    else if(sign>=0)
    {
	if(fabs(x)<fabs(y))
	    return x;
	else
	    return y;
    }
    return 0;
}

/**
* get the state from the neighbor cell or boundary.
* @param icoords
* @param dir
* @param comp
* @param state
* @return true,  valid state from neighbor cell
* 	   false, valid state from boundary
*/
bool Incompress_Solver_Smooth_3D_Cylindrical::getNeighborOrBoundaryState(
	int icoords[3],
	GRID_DIRECTION dir,
	L_STATE &state,
	double t)
{
    double crx_coords[MAXD];
    POINTER intfc_state;
    HYPER_SURF *hs;

    int index = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
    int comp = cell_center[index].comp;

    if (FT_StateStructAtGridCrossing(front,icoords,dir,
	    comp,&intfc_state,&hs,crx_coords) &&
	    wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
    {
	state.m_U[0] = getStateVel[0](intfc_state);
	state.m_U[1] = getStateVel[1](intfc_state);
	state.m_U[2] = getStateVel[2](intfc_state);
	return false;
    }
    else
    {
	int index_nb;
	switch(dir)
	{
	case WEST:
	    index_nb = d_index3d(icoords[0]-1,icoords[1],icoords[2],top_gmax);
	    break;
	case EAST:
	    index_nb = d_index3d(icoords[0]+1,icoords[1],icoords[2],top_gmax);
	    break;
	case SOUTH:
	    index_nb = d_index3d(icoords[0],icoords[1]-1,icoords[2],top_gmax);
	    break;
	case NORTH:
	    index_nb = d_index3d(icoords[0],icoords[1]+1,icoords[2],top_gmax);
	    break;
	case LOWER:
	    index_nb = d_index3d(icoords[0],icoords[1],icoords[2]-1,top_gmax);
	    break;
	case UPPER:
	    index_nb = d_index3d(icoords[0],icoords[1],icoords[2]+1,top_gmax);
	    break;
	default:
	    assert(false);
	}
	state = cell_center[index_nb].m_state;
	return true;
    }
}


/**
*
* @param state_left
* @param state_right
* @param ans
*/
void Incompress_Solver_Smooth_3D_Cylindrical::getRiemannSolution(
	EBM_COORD xyz,
	L_STATE &state_left,
	L_STATE &state_right,
	L_STATE &ans)
{
    L_STATE sl, sr;


    // rotate state
    if(xyz==COORD_X)
    {
	sl.m_U[0] = state_left.m_U[0];
	sl.m_U[1] = state_left.m_U[1];
	sl.m_U[2] = state_left.m_U[2];

	sr.m_U[0] = state_right.m_U[0];
	sr.m_U[1] = state_right.m_U[1];
	sr.m_U[2] = state_right.m_U[2];
    }
    else if(xyz==COORD_Y)
    {
	sl.m_U[0] = state_left.m_U[1];
	sl.m_U[1] = state_left.m_U[0];
	sl.m_U[2] = state_left.m_U[2];

	sr.m_U[0] = state_right.m_U[1];
	sr.m_U[1] = state_right.m_U[0];
	sr.m_U[2] = state_right.m_U[2];
    }
    else
    {
	sl.m_U[0] = state_left.m_U[2];
	sl.m_U[1] = state_left.m_U[1];
	sl.m_U[2] = state_left.m_U[0];

	sr.m_U[0] = state_right.m_U[2];
	sr.m_U[1] = state_right.m_U[1];
	sr.m_U[2] = state_right.m_U[0];
    }

    // calculate the Riemann solution
    double uL = sl.m_U[0];
    double uR = sr.m_U[0];

    // BCG, JCP 85, 257-283 (1989)
    // ut + uux = 0
    if(uL>=0 && (uL+uR)>=0)
	ans.m_U[0] = uL;
    else if(uL<0 && uR>0)
	ans.m_U[0] = 0;
    else
	ans.m_U[0] = uR;

    // vt + uvx = 0
    if(ans.m_U[0]>0)
	ans.m_U[1] = sl.m_U[1];
    else if(ans.m_U[0]<0)
	ans.m_U[1] = sr.m_U[1];
    else
	ans.m_U[1] = 1.0/2*(sl.m_U[1]+sr.m_U[1]);

    if(ans.m_U[0]>0)
	ans.m_U[2] = sl.m_U[2];
    else if(ans.m_U[0]<0)
	ans.m_U[2] = sr.m_U[2];
    else
	ans.m_U[2] = 1.0/2*(sl.m_U[2]+sr.m_U[2]);

    // rotate state
    if(xyz==COORD_X)
	; // do nothing
    else if(xyz==COORD_Y)
	std::swap(ans.m_U[0],ans.m_U[1]);
    else
	std::swap(ans.m_U[0],ans.m_U[2]);
}

void Incompress_Solver_Smooth_3D_Cylindrical::compDiffWithSmoothProperty_2nd_decoupled(void)
{
        COMPONENT comp;
        int index,index_nb[6],size;
        int I,I_nb[6];
        double coords[MAXD],crx_coords[MAXD];
	double coeff[18],mu0,rho,rhs,U0_nb[6],U1_nb[6],U2_nb[6],U0_center,U1_center,U2_center;
        L_STATE state;
        int i,j,k,l,nb,icoords[MAXD];
        double speed;
	double *x;
	GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
	POINTER intfc_state;
	HYPER_SURF *hs;
	int num_iter;
	double rel_residual;
	double r, rr, redge[2];

	max_speed = 0.0;
	setIndexMap();

	size = iupper - ilower;
	FT_VectorMemoryAlloc((POINTER*)&x,3*size,sizeof(double));

        PETSc solver;
        solver.Create(3*ilower, 3*iupper-1, 9, 9);
	// 7theta + 2r  for the first equation
	// 7z for the second equation
	// 7r + 2theta for the third equation

	solver.Reset_A();
	solver.Reset_b();
	solver.Reset_x();

	printf("\nIn diffusion solver ,m_dt = %.16g\n", m_dt);

	for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I  = ijk_to_I[i][j][k];
            if (I == -1) continue;

            index  = d_index3d(i,j,k,top_gmax);	
	//6 neighbours of the center cell
            index_nb[0] = d_index3d(i-1,j,k,top_gmax);
            index_nb[1] = d_index3d(i+1,j,k,top_gmax);
            index_nb[2] = d_index3d(i,j-1,k,top_gmax);
            index_nb[3] = d_index3d(i,j+1,k,top_gmax);
	    index_nb[4] = d_index3d(i,j,k-1,top_gmax);
	    index_nb[5] = d_index3d(i,j,k+1,top_gmax);

        
	//6 neighbours of the center cell
            I_nb[0] = ijk_to_I[i-1][j][k];
            I_nb[1] = ijk_to_I[i+1][j][k];
            I_nb[2] = ijk_to_I[i][j-1][k];
            I_nb[3] = ijk_to_I[i][j+1][k];
            I_nb[4] = ijk_to_I[i][j][k-1];
            I_nb[5] = ijk_to_I[i][j][k+1];
	

	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;
	    comp = top_comp[index];

	    mu0 = cell_center[index].m_state.m_mu;
	    rho = cell_center[index].m_state.m_rho;
	    U0_center = cell_center[index].m_state.m_U[0];
	    U1_center = cell_center[index].m_state.m_U[1];
	    U2_center = cell_center[index].m_state.m_U[2];

	  
            for (nb = 0; nb < 6; nb++)
            {
                if (FT_StateStructAtGridCrossing(front,icoords,dir[nb],
                            comp,&intfc_state,&hs,crx_coords) &&
                            wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	       	{
		    U0_nb[nb] = getStateVel[0](intfc_state);
		    U1_nb[nb] = getStateVel[1](intfc_state);
		    U2_nb[nb] = getStateVel[2](intfc_state);
		}
                else
		{
		    U0_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[0];
		    U1_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[1];
		    U2_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[2];
		}
	    }

            getRectangleCenter(index, coords);
            computeSourceTerm(coords, state);


    //Setting the coeffecients for the first equation

	    rr = cell_center[index].m_coords[2] * cell_center[index].m_coords[2];
	    r = cell_center[index].m_coords[2];
	    redge[0] = 0.5 * (cell_center[index].m_coords[2] + cell_center[index_nb[4]].m_coords[2]);
	    redge[1] = 0.5 * (cell_center[index].m_coords[2] + cell_center[index_nb[5]].m_coords[2]);


	    coeff[0] = 0.5*m_dt/rho * mu0/(rr*top_h[0]*top_h[0]);
	    coeff[1] = 0.5*m_dt/rho * mu0/(rr*top_h[0]*top_h[0]);
	    coeff[2] = 0.5*m_dt/rho * mu0/(top_h[1]*top_h[1]);
	    coeff[3] = 0.5*m_dt/rho * mu0/(top_h[1]*top_h[1]);

	    coeff[4] = 0.5*m_dt/rho * mu0*redge[0]/(r*top_h[2]*top_h[2]);
	    coeff[5] = 0.5*m_dt/rho * mu0*redge[1]/(r*top_h[2]*top_h[2]);

	    for (nb = 0; nb < 6; nb++) //For the second order boundary
	    {
		if (I_nb[nb] == -1)
		{
		    coeff[nb] = 2.0*coeff[nb];
		}
	    }

	    solver.Set_A(I*3,I*3,1.0+coeff[0]+coeff[1]+coeff[2]+coeff[3]+coeff[4]+coeff[5]+0.5*m_dt*mu0/(rho*rr));
	    rhs = (1.0-coeff[0]-coeff[1]-coeff[2]-coeff[3]-coeff[4]-coeff[5]-0.5*m_dt*mu0/(rho*rr))*cell_center[index].m_state.m_U[0];


	    coeff[6] = -0.5*m_dt/rho * mu0/(rr*top_h[0]);
	    coeff[7] =  0.5*m_dt/rho * mu0/(rr*top_h[0]);


	    for (nb = 0; nb < 6; nb++)
	    {
		if (I_nb[nb] != -1)
		{
		    solver.Set_A(I*3,I_nb[nb]*3,-coeff[nb]);
		    rhs += coeff[nb]*U0_nb[nb];
		}
		else
		    rhs += 2.0*coeff[nb]*U0_nb[nb];
	    }

	    if (I_nb[0] != -1)
	    {
		solver.Set_A(I*3,I_nb[0]*3+2,-coeff[6]);
		rhs += coeff[6]*U2_nb[0];
	    }
	    else
	    {
		solver.Set_A(I*3,I*3+2,coeff[6]);
		rhs += -coeff[6]*U2_center;
		coeff[6] = coeff[6] * 2.0;
		rhs += 2.0*coeff[6]*U2_nb[0];
	    }

	    if (I_nb[1] != -1)
	    {
		solver.Set_A(I*3,I_nb[1]*3+2,-coeff[7]);
		rhs += coeff[7]*U2_nb[1];
	    }
	    else
	    {
		solver.Set_A(I*3,I*3+2,coeff[7]);
		rhs += -coeff[7]*U2_center;
		coeff[7] = coeff[7] * 2.0;
		rhs += 2.0*coeff[7]*U2_nb[1];
	    }

	    //rhs -= m_dt*cell_center[index].m_state.m_U[0]*cell_center[index].m_state.m_U[2]/r; //Source term in advection step
	    rhs += m_dt*state.m_U[0];
	    rhs += m_dt*cell_center[index].m_state.f_surf[0];
	    rhs -= m_dt*cell_center[index].m_state.grad_q[0]/rho;
	    rhs -= m_dt*cell_center[index].m_state.m_adv[0];



	    solver.Set_b(I*3, rhs);

	    /************************************************************************/


    //Setting the coeffecients for the second equation


	    coeff[0] = 0.5*m_dt/rho * mu0/(rr*top_h[0]*top_h[0]);
	    coeff[1] = 0.5*m_dt/rho * mu0/(rr*top_h[0]*top_h[0]);
	    coeff[2] = 0.5*m_dt/rho * mu0/(top_h[1]*top_h[1]);
	    coeff[3] = 0.5*m_dt/rho * mu0/(top_h[1]*top_h[1]);

	    coeff[4] = 0.5*m_dt/rho * mu0*redge[0]/(r*top_h[2]*top_h[2]);
	    coeff[5] = 0.5*m_dt/rho * mu0*redge[1]/(r*top_h[2]*top_h[2]);

	    for (nb = 0; nb < 6; nb++)
	    {
		if (I_nb[nb] == -1)
		    coeff[nb] = 2.0*coeff[nb];
	    }

	    solver.Set_A(I*3+1,I*3+1,1.0+coeff[0]+coeff[1]+coeff[2]+coeff[3]+coeff[4]+coeff[5]);
	    rhs = (1.0-coeff[0]-coeff[1]-coeff[2]-coeff[3]-coeff[4]-coeff[5])*cell_center[index].m_state.m_U[1];


	    for (nb = 0; nb < 6; nb++)
	    {
		if (I_nb[nb] != -1)
		{
		    solver.Set_A(I*3+1,I_nb[nb]*3+1,-coeff[nb]);
		    rhs += coeff[nb]*U1_nb[nb];
		}
		else
		    rhs += 2.0*coeff[nb]*U1_nb[nb];
	    }

	    rhs += m_dt*state.m_U[1];
	    rhs += m_dt*cell_center[index].m_state.f_surf[1];
	    rhs -= m_dt*cell_center[index].m_state.grad_q[1]/rho;
	    rhs -= m_dt*cell_center[index].m_state.m_adv[1];


	    solver.Set_b(I*3+1, rhs);

	    /************************************************************************/

    //Setting the coeffecients for the third equation

	    coeff[0] = 0.5*m_dt/rho * mu0/(rr*top_h[0]*top_h[0]);
	    coeff[1] = 0.5*m_dt/rho * mu0/(rr*top_h[0]*top_h[0]);
	    coeff[2] = 0.5*m_dt/rho * mu0/(top_h[1]*top_h[1]);
	    coeff[3] = 0.5*m_dt/rho * mu0/(top_h[1]*top_h[1]);

	    coeff[4] = 0.5*m_dt/rho * mu0*redge[0]/(r*top_h[2]*top_h[2]);
	    coeff[5] = 0.5*m_dt/rho * mu0*redge[1]/(r*top_h[2]*top_h[2]);

	    for (nb = 0; nb < 6; nb++)
	    {
		if (I_nb[nb] == -1)
		    coeff[nb] = 2.0*coeff[nb];
	    }

	    solver.Set_A(I*3+2,I*3+2,1.0+coeff[0]+coeff[1]+coeff[2]+coeff[3]+coeff[4]+coeff[5]+0.5*m_dt*mu0/(rho*rr));
	    rhs = (1.0-coeff[0]-coeff[1]-coeff[2]-coeff[3]-coeff[4]-coeff[5]-0.5*m_dt*mu0/(rho*rr))*cell_center[index].m_state.m_U[2];


	    coeff[6] =  0.5*m_dt/rho * mu0/(rr*top_h[0]);
	    coeff[7] = -0.5*m_dt/rho * mu0/(rr*top_h[0]);


	    for (nb = 0; nb < 6; nb++)
	    {
		if (I_nb[nb] != -1)
		{
		    solver.Set_A(I*3+2,I_nb[nb]*3+2,-coeff[nb]);
		    rhs += coeff[nb]*U2_nb[nb];
		}
		else
		    rhs += 2.0*coeff[nb]*U2_nb[nb];
	    }


	    if (I_nb[0] != -1)
	    {
		solver.Set_A(I*3+2, I_nb[0]*3, -coeff[6]);
		rhs += coeff[6]*U0_nb[0];
	    }
	    else
	    {
		solver.Set_A(I*3+2, I*3, coeff[6]);
		rhs += -coeff[6]*U0_center;
		coeff[6] = coeff[6] * 2.0;
		rhs += 2.0*coeff[6]*U0_nb[0];
	    }
	    
	    if (I_nb[1] != -1)
	    {
		solver.Set_A(I*3+2, I_nb[1]*3, -coeff[7]);
		rhs += coeff[7]*U0_nb[1];
	    }
	    else
	    {
		solver.Set_A(I*3+2, I*3, coeff[7]);
		rhs += -coeff[7]*U0_center;
		coeff[7] = coeff[7] * 2.0;
		rhs += 2.0*coeff[7]*U0_nb[1];
	    }

	    //rhs += m_dt* cell_center[index].m_state.m_U[0] * cell_center[index].m_state.m_U[0] / r; //Source term in advection step
	    rhs += m_dt*state.m_U[2];
	    rhs += m_dt*cell_center[index].m_state.f_surf[2];
	    rhs -= m_dt*cell_center[index].m_state.grad_q[2]/rho;
	    rhs -= m_dt*cell_center[index].m_state.m_adv[2];

	    solver.Set_b(I*3+2, rhs);
        }

        solver.SetMaxIter(40000);
        solver.SetTol(1e-15);

	start_clock("Before Petsc Solve");
        solver.Solve_GMRES();
	solver.GetNumIterations(&num_iter);
	solver.GetFinalRelativeResidualNorm(&rel_residual);
	stop_clock("After Petsc Solve");

        // get back the solution
        solver.Get_x(x);

        if (debugging("PETSc"))
            (void) printf("Incompress_Solver_Smooth_3D_Cylindrical::compDiffWithSmoothProperty_2nd_decoupled: "
                        "num_iter = %d, rel_residual = %g. \n",
                        num_iter,rel_residual);


	for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I = ijk_to_I[i][j][k];
            index = d_index3d(i,j,k,top_gmax);
            if (I >= 0)
            {
                cell_center[index].m_state.m_U[0] = x[I*3-ilower*3];
                cell_center[index].m_state.m_U[1] = x[I*3-ilower*3+1];
		cell_center[index].m_state.m_U[2] = x[I*3-ilower*3+2];
                speed = fabs(cell_center[index].m_state.m_U[0]) +
                        fabs(cell_center[index].m_state.m_U[1]) +
			fabs(cell_center[index].m_state.m_U[2]);
                if (speed > max_speed)
                    max_speed = speed;
            }
            else
            {
                cell_center[index].m_state.m_U[0] = 0.0;
                cell_center[index].m_state.m_U[1] = 0.0;
		cell_center[index].m_state.m_U[2] = 0.0;
            }
        }

        for (l = 0; l < 3; ++l)
        {
	    for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)

            {
                index  = d_index3d(i,j,k,top_gmax);
                array[index] = cell_center[index].m_state.m_U[l];
            }
            scatMeshArray();
	    for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)

            {
                index  = d_index3d(i,j,k,top_gmax);
                cell_center[index].m_state.m_U[l] = array[index];
            }
        }

        pp_global_max(&max_speed,1);

        free_these(1,x);
}       /* end 3D_Cylindrical::compDiffWithSmoothProperty_2nd_decoupled */

void Incompress_Solver_Smooth_3D_Cylindrical::compDiff_CellFace(
	PETSc *pSolver,
	int I,
	int I_nb[18],
	double U_center[3],
	double U_nb[3][18],
	int flag[6],
	int equation_index,
	int vel_comp,
	int face_index,
	double coeff)
{
    double rhs;

    if (flag[face_index] == 1)
	rhs = coeff*U_nb[vel_comp][face_index];
    else 
    {
	rhs = coeff*(U_nb[vel_comp][face_index]+U_center[vel_comp])/4.0;

	pSolver->Add_A(I*3+equation_index,I*3+vel_comp,                -coeff/4.0);
	pSolver->Add_A(I*3+equation_index,I_nb[face_index]*3+vel_comp, -coeff/4.0);
    }
    pSolver->Add_b(I*3+equation_index,rhs);

}

void Incompress_Solver_Smooth_3D_Cylindrical::compDiff_CellCorner(
	PETSc *pSolver,
	int I, 
	int I_nb[18],
	double U_center[3],
	double U_nb[3][18], 
	int flag[6],
	int equation_index,
	int vel_comp, 
	int corner_index, 
	double coeff)
{
    double rhs;

    if (corner_index == 6)
    {
	if      (flag[0] == 1 && flag[2] == 0)
	    rhs = coeff*U_nb[vel_comp][0];
	else if (flag[0] == 0 && flag[2] == 1)
	    rhs = coeff*U_nb[vel_comp][2];
	else if (flag[0] == 1 && flag[2] == 1)
	    rhs = coeff*U_nb[vel_comp][0];
	else {
	    rhs = coeff*(U_nb[vel_comp][0]
		        +U_nb[vel_comp][2]
			+U_nb[vel_comp][6]
			+U_center[vel_comp])/8.0;
		    
	    pSolver->Add_A(I*3+equation_index,I*3+vel_comp,       -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[0]*3+vel_comp, -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[2]*3+vel_comp, -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[6]*3+vel_comp, -coeff/8.0);
	}
	pSolver->Add_b(I*3+equation_index,rhs);
    }
    else if (corner_index == 7)
    {
	if      (flag[1] == 1 && flag[2] == 0)
	    rhs = coeff*U_nb[vel_comp][1];
	else if (flag[1] == 0 && flag[2] == 1)
	    rhs = coeff*U_nb[vel_comp][2];
	else if (flag[1] == 1 && flag[2] == 1)
	    rhs = coeff*U_nb[vel_comp][1];
	else {
	    rhs = coeff*(U_nb[vel_comp][1]
		        +U_nb[vel_comp][2]
			+U_nb[vel_comp][7]
			+U_center[vel_comp])/8.0;
		    
	    pSolver->Add_A(I*3+equation_index,I*3+vel_comp,       -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[1]*3+vel_comp, -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[2]*3+vel_comp, -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[7]*3+vel_comp, -coeff/8.0);
	}
	pSolver->Add_b(I*3+equation_index,rhs);
    }
    else if (corner_index == 8)
    {
	if      (flag[1] == 1 && flag[3] == 0)
	    rhs = coeff*U_nb[vel_comp][1];
	else if (flag[1] == 0 && flag[3] == 1)
	    rhs = coeff*U_nb[vel_comp][3];
	else if (flag[1] == 1 && flag[3] == 1)
	    rhs = coeff*U_nb[vel_comp][1];
	else {
	    rhs = coeff*(U_nb[vel_comp][1]
		        +U_nb[vel_comp][3]
			+U_nb[vel_comp][8]
			+U_center[vel_comp])/8.0;
		    
	    pSolver->Add_A(I*3+equation_index,I*3+vel_comp,       -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[1]*3+vel_comp, -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[3]*3+vel_comp, -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[8]*3+vel_comp, -coeff/8.0);
	}
	pSolver->Add_b(I*3+equation_index,rhs);
    }
    else if (corner_index == 9)
    {
	if      (flag[0] == 1 && flag[3] == 0)
	    rhs = coeff*U_nb[vel_comp][0];
	else if (flag[0] == 0 && flag[3] == 1)
	    rhs = coeff*U_nb[vel_comp][3];
	else if (flag[0] == 1 && flag[3] == 1)
	    rhs = coeff*U_nb[vel_comp][0];
	else {
	    rhs = coeff*(U_nb[vel_comp][0]
		        +U_nb[vel_comp][3]
			+U_nb[vel_comp][9]
			+U_center[vel_comp])/8.0;
		    
	    pSolver->Add_A(I*3+equation_index,I*3+vel_comp,       -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[0]*3+vel_comp, -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[3]*3+vel_comp, -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[9]*3+vel_comp, -coeff/8.0);
	}
	pSolver->Add_b(I*3+equation_index,rhs);
    }
    else if (corner_index == 10)
    {
	if      (flag[2] == 1 && flag[4] == 0)
	    rhs = coeff*U_nb[vel_comp][2];
	else if (flag[2] == 0 && flag[4] == 1)
	    rhs = coeff*U_nb[vel_comp][4];
	else if (flag[2] == 1 && flag[4] == 1)
	    rhs = coeff*U_nb[vel_comp][2];
	else {
	    rhs = coeff*(U_nb[vel_comp][2]
		        +U_nb[vel_comp][4]
			+U_nb[vel_comp][10]
			+U_center[vel_comp])/8.0;
		    
	    pSolver->Add_A(I*3+equation_index,I*3+vel_comp,       -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[2]*3+vel_comp, -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[4]*3+vel_comp, -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[10]*3+vel_comp, -coeff/8.0);
	}
	pSolver->Add_b(I*3+equation_index,rhs);
    }
    else if (corner_index == 11)
    {
	if      (flag[3] == 1 && flag[4] == 0)
	    rhs = coeff*U_nb[vel_comp][3];
	else if (flag[3] == 0 && flag[4] == 1)
	    rhs = coeff*U_nb[vel_comp][4];
	else if (flag[3] == 1 && flag[4] == 1)
	    rhs = coeff*U_nb[vel_comp][3];
	else {
	    rhs = coeff*(U_nb[vel_comp][3]
		        +U_nb[vel_comp][4]
			+U_nb[vel_comp][11]
			+U_center[vel_comp])/8.0;
		    
	    pSolver->Add_A(I*3+equation_index,I*3+vel_comp,       -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[3]*3+vel_comp, -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[4]*3+vel_comp, -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[11]*3+vel_comp, -coeff/8.0);
	}
	pSolver->Add_b(I*3+equation_index,rhs);
    }
    else if (corner_index == 12)
    {
	if      (flag[3] == 1 && flag[5] == 0)
	    rhs = coeff*U_nb[vel_comp][3];
	else if (flag[3] == 0 && flag[5] == 1)
	    rhs = coeff*U_nb[vel_comp][5];
	else if (flag[3] == 1 && flag[5] == 1)
	    rhs = coeff*U_nb[vel_comp][3];
	else {
	    rhs = coeff*(U_nb[vel_comp][3]
		        +U_nb[vel_comp][5]
			+U_nb[vel_comp][12]
			+U_center[vel_comp])/8.0;
		    
	    pSolver->Add_A(I*3+equation_index,I*3+vel_comp,       -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[3]*3+vel_comp, -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[5]*3+vel_comp, -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[12]*3+vel_comp, -coeff/8.0);
	}
	pSolver->Add_b(I*3+equation_index,rhs);
    }
    else if (corner_index == 13)
    {
	if      (flag[2] == 1 && flag[5] == 0)
	    rhs = coeff*U_nb[vel_comp][2];
	else if (flag[2] == 0 && flag[5] == 1)
	    rhs = coeff*U_nb[vel_comp][5];
	else if (flag[2] == 1 && flag[5] == 1)
	    rhs = coeff*U_nb[vel_comp][2];
	else {
	    rhs = coeff*(U_nb[vel_comp][2]
		        +U_nb[vel_comp][5]
			+U_nb[vel_comp][13]
			+U_center[vel_comp])/8.0;
		    
	    pSolver->Add_A(I*3+equation_index,I*3+vel_comp,       -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[2]*3+vel_comp, -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[5]*3+vel_comp, -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[13]*3+vel_comp, -coeff/8.0);
	}
	pSolver->Add_b(I*3+equation_index,rhs);
    }
    else if (corner_index == 14)
    {
	if      (flag[0] == 1 && flag[4] == 0)
	    rhs = coeff*U_nb[vel_comp][0];
	else if (flag[0] == 0 && flag[4] == 1)
	    rhs = coeff*U_nb[vel_comp][4];
	else if (flag[0] == 1 && flag[4] == 1)
	    rhs = coeff*U_nb[vel_comp][0];
	else {
	    rhs = coeff*(U_nb[vel_comp][0]
		        +U_nb[vel_comp][4]
			+U_nb[vel_comp][14]
			+U_center[vel_comp])/8.0;
		    
	    pSolver->Add_A(I*3+equation_index,I*3+vel_comp,       -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[0]*3+vel_comp, -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[4]*3+vel_comp, -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[14]*3+vel_comp, -coeff/8.0);
	}
	pSolver->Add_b(I*3+equation_index,rhs);
    }
    else if (corner_index == 15)
    {
	if      (flag[1] == 1 && flag[4] == 0)
	    rhs = coeff*U_nb[vel_comp][1];
	else if (flag[1] == 0 && flag[4] == 1)
	    rhs = coeff*U_nb[vel_comp][4];
	else if (flag[1] == 1 && flag[4] == 1)
	    rhs = coeff*U_nb[vel_comp][1];
	else {
	    rhs = coeff*(U_nb[vel_comp][1]
		        +U_nb[vel_comp][4]
			+U_nb[vel_comp][15]
			+U_center[vel_comp])/8.0;
		    
	    pSolver->Add_A(I*3+equation_index,I*3+vel_comp,       -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[1]*3+vel_comp, -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[4]*3+vel_comp, -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[15]*3+vel_comp, -coeff/8.0);
	}
	pSolver->Add_b(I*3+equation_index,rhs);
    }
    else if (corner_index == 16)
    {
	if      (flag[1] == 1 && flag[5] == 0)
	    rhs = coeff*U_nb[vel_comp][1];
	else if (flag[1] == 0 && flag[5] == 1)
	    rhs = coeff*U_nb[vel_comp][5];
	else if (flag[1] == 1 && flag[5] == 1)
	    rhs = coeff*U_nb[vel_comp][1];
	else {
	    rhs = coeff*(U_nb[vel_comp][1]
		        +U_nb[vel_comp][5]
			+U_nb[vel_comp][16]
			+U_center[vel_comp])/8.0;
		    
	    pSolver->Add_A(I*3+equation_index,I*3+vel_comp,       -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[1]*3+vel_comp, -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[5]*3+vel_comp, -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[16]*3+vel_comp, -coeff/8.0);
	}
	pSolver->Add_b(I*3+equation_index,rhs);
    }
    else if (corner_index == 17)
    {
	if      (flag[0] == 1 && flag[5] == 0)
	    rhs = coeff*U_nb[vel_comp][0];
	else if (flag[0] == 0 && flag[5] == 1)
	    rhs = coeff*U_nb[vel_comp][5];
	else if (flag[0] == 1 && flag[5] == 1)
	    rhs = coeff*U_nb[vel_comp][0];
	else {
	    rhs = coeff*(U_nb[vel_comp][0]
		        +U_nb[vel_comp][5]
			+U_nb[vel_comp][17]
			+U_center[vel_comp])/8.0;
		    
	    pSolver->Add_A(I*3+equation_index,I*3+vel_comp,       -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[0]*3+vel_comp, -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[5]*3+vel_comp, -coeff/8.0);
	    pSolver->Add_A(I*3+equation_index,I_nb[17]*3+vel_comp, -coeff/8.0);
	}
	pSolver->Add_b(I*3+equation_index,rhs);
    }
}


void Incompress_Solver_Smooth_3D_Cylindrical::compDiffWithSmoothProperty_2nd_coupled(void)
{
        COMPONENT comp;
        int index,index_nb[18],size;
        int I,I_nb[18];
        double coords[MAXD],crx_coords[MAXD];
	double coeff_temp, coeff0, coeff1;
	double mu0,rho,rhs,mu[6],mu_edge[6];
	int flag[6];
        L_STATE state;
        int i,j,k,l,nb,icoords[MAXD];
        double speed;
	double *x;
	GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
	POINTER intfc_state;
	HYPER_SURF *hs;
	PetscInt num_iter;
	double rel_residual;
	double r0, rr, r_edge[6],r[6];
	double dr = top_h[2];
	double dz = top_h[1];
	double dtheta = top_h[0];

	double U_center[3];
	double U_nb[3][18];

	max_speed = 0.0;
	setIndexMap();

	size = iupper - ilower;
	FT_VectorMemoryAlloc((POINTER*)&x,3*size,sizeof(double));

        PETSc solver;
        solver.Create(3*ilower, 3*iupper-1, 25, 25);
	// 7theta + 9z + 9r  for the first equation
	// 7z + 9theta + 9r  for the second equation
	// 7r + 9theta + 9z  for the third equation

	solver.Reset_A();
	solver.Reset_b();
	solver.Reset_x();

	printf("\nIn diffusion solver ,m_dt = %.16g\n", m_dt);

	for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I  = ijk_to_I[i][j][k];
            if (I == -1) continue;

            index  = d_index3d(i,j,k,top_gmax);	
	//6 neighbours of the center cell
            index_nb[0] = d_index3d(i-1,j,k,top_gmax);
            index_nb[1] = d_index3d(i+1,j,k,top_gmax);
            index_nb[2] = d_index3d(i,j-1,k,top_gmax);
            index_nb[3] = d_index3d(i,j+1,k,top_gmax);
	    index_nb[4] = d_index3d(i,j,k-1,top_gmax);
	    index_nb[5] = d_index3d(i,j,k+1,top_gmax);

	//theta-z cut neighbours
	    index_nb[6] = d_index3d(i-1,j-1,k,top_gmax);
	    index_nb[7] = d_index3d(i+1,j-1,k,top_gmax);
	    index_nb[8] = d_index3d(i+1,j+1,k,top_gmax);
	    index_nb[9] = d_index3d(i-1,j+1,k,top_gmax);
	
	//z-r cut neighbours
	    index_nb[10] = d_index3d(i,j-1,k-1,top_gmax);
	    index_nb[11] = d_index3d(i,j+1,k-1,top_gmax);
	    index_nb[12] = d_index3d(i,j+1,k+1,top_gmax);
	    index_nb[13] = d_index3d(i,j-1,k+1,top_gmax);

	//theta-r cut neighbours
	    index_nb[14] = d_index3d(i-1,j,k-1,top_gmax);
	    index_nb[15] = d_index3d(i+1,j,k-1,top_gmax);
	    index_nb[16] = d_index3d(i+1,j,k+1,top_gmax);
	    index_nb[17] = d_index3d(i-1,j,k+1,top_gmax);



        
	//6 neighbours of the center cell
            I_nb[0] = ijk_to_I[i-1][j][k];
            I_nb[1] = ijk_to_I[i+1][j][k];
            I_nb[2] = ijk_to_I[i][j-1][k];
            I_nb[3] = ijk_to_I[i][j+1][k];
            I_nb[4] = ijk_to_I[i][j][k-1];
            I_nb[5] = ijk_to_I[i][j][k+1];
		
	//theta-z cut neighbours
            I_nb[6] = ijk_to_I[i-1][j-1][k];
	    I_nb[7] = ijk_to_I[i+1][j-1][k];
	    I_nb[8] = ijk_to_I[i+1][j+1][k];
	    I_nb[9] = ijk_to_I[i-1][j+1][k];
	//z-r cut neighbours
	    I_nb[10] = ijk_to_I[i][j-1][k-1];
	    I_nb[11] = ijk_to_I[i][j+1][k-1];
	    I_nb[12] = ijk_to_I[i][j+1][k+1];
	    I_nb[13] = ijk_to_I[i][j-1][k+1];
	//theta-r cut neighbours
	    I_nb[14] = ijk_to_I[i-1][j][k-1];
	    I_nb[15] = ijk_to_I[i+1][j][k-1];
	    I_nb[16] = ijk_to_I[i+1][j][k+1];
	    I_nb[17] = ijk_to_I[i-1][j][k+1];



	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;
	    comp = top_comp[index];

	    mu0 = cell_center[index].m_state.m_mu;
	    rho = cell_center[index].m_state.m_rho;
	    U_center[0] = cell_center[index].m_state.m_U[0];
	    U_center[1] = cell_center[index].m_state.m_U[1];
	    U_center[2] = cell_center[index].m_state.m_U[2];


	    r0 = cell_center[index].m_coords[2];
	    rr = r0*r0;

	  
            for (nb = 0; nb < 6; nb++)
            {
                if (FT_StateStructAtGridCrossing(front,icoords,dir[nb],
                            comp,&intfc_state,&hs,crx_coords) &&
                            wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	       	{
		    flag[nb] = 1;
		    U_nb[0][nb] = getStateVel[0](intfc_state);
		    U_nb[1][nb] = getStateVel[1](intfc_state);
		    U_nb[2][nb] = getStateVel[2](intfc_state);

		    mu[nb] = mu0;
		    mu_edge[nb] = mu0;
		}
                else
		{
		    flag[nb] = 0;
		    U_nb[0][nb] = cell_center[index_nb[nb]].m_state.m_U[0];
		    U_nb[1][nb] = cell_center[index_nb[nb]].m_state.m_U[1];
		    U_nb[2][nb] = cell_center[index_nb[nb]].m_state.m_U[2];

		    mu[nb] = cell_center[index_nb[nb]].m_state.m_mu;
		    mu_edge[nb] = 0.5*(mu0 + mu[nb]);
		}

		r[nb] = cell_center[index_nb[nb]].m_coords[2];
		r_edge[nb] = 0.5 * (r[nb] + r0);
	    }

	    for (nb = 6; nb < 18; nb++) //corner values for interior
	    {
		if (I_nb[nb] != -1)
		{
		    U_nb[0][nb] = cell_center[index_nb[nb]].m_state.m_U[0];
		    U_nb[1][nb] = cell_center[index_nb[nb]].m_state.m_U[1];
		    U_nb[2][nb] = cell_center[index_nb[nb]].m_state.m_U[2];
		}
	    }

	    //source term
            getRectangleCenter(index, coords);
            computeSourceTerm(coords, state);


	    /****************** The first equation about theta *************/

	    solver.Add_A(I*3,I*3,1.0);
	    rhs = U_center[0];

	    /////////////////    term d(tao_z0)/dz      //////////////////
	    //first term in tao_z0: mu*du_0/dz

	    if (flag[2] == 1)
	    {
		coeff_temp = m_dt/rho * mu_edge[2]/(dz*dz);
		solver.Add_A(I*3,I*3, coeff_temp);
		rhs += 2.0*coeff_temp*U_nb[0][2];
		rhs -= coeff_temp*U_center[0];
	    }
	    else
	    {
		coeff_temp = 0.5*m_dt/rho * mu_edge[2]/(dz*dz);
		solver.Add_A(I*3,I*3, coeff_temp);
		solver.Add_A(I*3,I_nb[2]*3, -coeff_temp);
		rhs += coeff_temp*U_nb[0][2];
		rhs -= coeff_temp*U_center[0];
	    }

	    if (flag[3] == 1)
	    {
		coeff_temp = m_dt/rho * mu_edge[3]/(dz*dz);
		solver.Add_A(I*3,I*3, coeff_temp);
		rhs += 2.0*coeff_temp*U_nb[0][3];
		rhs -= coeff_temp*U_center[0];
	    }
	    else
	    {
		coeff_temp = 0.5*m_dt/rho *mu_edge[3]/(dz*dz);
		solver.Add_A(I*3,I*3, coeff_temp);
		solver.Add_A(I*3,I_nb[3]*3, -coeff_temp);
		rhs += coeff_temp*U_nb[0][3];
		rhs -= coeff_temp*U_center[0];
	    }

	    //second term mu/r * du_z/d0
	    // 4 corners
	    //
	    coeff0 = m_dt/rho * mu_edge[2]/r_edge[2] / (dtheta*dz);
	    coeff1 = m_dt/rho * mu_edge[3]/r_edge[3] / (dtheta*dz);

	    compDiff_CellCorner(&solver, I, I_nb, U_center, U_nb, flag,
		    0,1,6,  coeff0);
	    compDiff_CellCorner(&solver, I, I_nb, U_center, U_nb, flag,
		    0,1,7, -coeff0);
	    compDiff_CellCorner(&solver, I, I_nb, U_center, U_nb, flag,
		    0,1,8,  coeff1);
	    compDiff_CellCorner(&solver, I, I_nb, U_center, U_nb, flag,
		    0,1,9, -coeff1);

	    ////////////      term d(tao_r0)/dr       /////////////////
	    //first term in tao_r0: mu*du_0/dr

	    if (flag[4] == 1)
	    {
		coeff_temp = m_dt/rho * mu_edge[4]/(dr*dr);
		solver.Add_A(I*3,I*3, coeff_temp);
		rhs += 2.0*coeff_temp*U_nb[0][4];
		rhs -= coeff_temp*U_center[0];
	    }
	    else
	    {
		coeff_temp = 0.5*m_dt/rho * mu_edge[4]/(dr*dr);
		solver.Add_A(I*3,I*3, coeff_temp);
		solver.Add_A(I*3,I_nb[4]*3, -coeff_temp);
		rhs += coeff_temp*U_nb[0][4];
		rhs -= coeff_temp*U_center[0];
	    }

	    if (flag[5] == 1)
	    {
		coeff_temp = m_dt/rho * mu_edge[5]/(dr*dr);
		solver.Add_A(I*3,I*3, coeff_temp);
		rhs += 2.0*coeff_temp*U_nb[0][5];
		rhs -= coeff_temp*U_center[0];
	    }
	    else
	    {
		coeff_temp = 0.5*m_dt/rho *mu_edge[5]/(dz*dz);
		solver.Add_A(I*3,I*3, coeff_temp);
		solver.Add_A(I*3,I_nb[5]*3, -coeff_temp);
		rhs += coeff_temp*U_nb[0][5];
		rhs -= coeff_temp*U_center[0];
	    }

	    //second term -u_0/r
	    //
	    coeff0 = m_dt/rho * mu_edge[4]/r_edge[4]/dr;
	    coeff1 = m_dt/rho * mu_edge[5]/r_edge[5]/dr;

	    compDiff_CellFace(&solver, I, I_nb, U_center, U_nb, flag,
		    0,0,4, coeff0);
	    compDiff_CellFace(&solver, I, I_nb, U_center, U_nb, flag,
		    0,0,5,-coeff1);

	    //third term 
	    coeff0 = m_dt/rho * mu_edge[4]/r_edge[4]/(dtheta*dr);
	    coeff1 = m_dt/rho * mu_edge[5]/r_edge[5]/(dtheta*dr);

	    compDiff_CellCorner(&solver, I, I_nb, U_center, U_nb, flag,
		    0,2,14,   coeff0);
	    compDiff_CellCorner(&solver, I, I_nb, U_center, U_nb, flag,
		    0,2,15,  -coeff0);
	    compDiff_CellCorner(&solver, I, I_nb, U_center, U_nb, flag,
		    0,2,16,   coeff1);
	    compDiff_CellCorner(&solver, I, I_nb, U_center, U_nb, flag,
		    0,2,17,  -coeff1);

	    /////////////////   term 1/r * dtao_00/d0   //////////////
	    //first term 2*mu/r*du_0/d0
	    
	    if (flag[0] == 1)
	    {
		coeff_temp = 2.0*m_dt/rho/r0 * mu_edge[0]/r_edge[0] / (dtheta*dtheta);
		solver.Add_A(I*3,I*3, coeff_temp);
		rhs += 2.0*coeff_temp*U_nb[0][0];
		rhs -= coeff_temp*U_center[0];
	    }
	    else
	    {
		coeff_temp = m_dt/rho/r0 * mu_edge[0]/r_edge[0] / (dtheta*dtheta);
		solver.Add_A(I*3,I*3, coeff_temp);
		solver.Add_A(I*3,I_nb[0]*3, -coeff_temp);
		rhs += coeff_temp*U_nb[0][0];
		rhs -= coeff_temp*U_center[0];
	    }

	    if (flag[1] == 1)
	    {
		coeff_temp = 2.0*m_dt/rho/r0 * mu_edge[1]/r_edge[1] / (dtheta*dtheta);
		solver.Add_A(I*3,I*3, coeff_temp);
		rhs += 2.0*coeff_temp*U_nb[0][1];
		rhs -= coeff_temp*U_center[0];
	    }
	    else
	    {
		coeff_temp = m_dt/rho/r0 * mu_edge[1]/r_edge[1] / (dtheta*dtheta);
		solver.Add_A(I*3,I*3, coeff_temp);
		solver.Add_A(I*3,I_nb[1]*3, -coeff_temp);
		rhs += coeff_temp*U_nb[0][1];
		rhs -= coeff_temp*U_center[0];
	    }

	    //second term
	    coeff0 = 2.0*m_dt/rho/r0 * mu_edge[0]/r_edge[0] / dtheta;
	    coeff1 = 2.0*m_dt/rho/r0 * mu_edge[1]/r_edge[1] / dtheta;

	    compDiff_CellFace(&solver, I, I_nb, U_center, U_nb, flag,
		    0,2,0,  -coeff0);
	    compDiff_CellFace(&solver, I, I_nb, U_center, U_nb, flag,
		    0,2,1,   coeff1);

	    /////////////////////   term 2*tao_r0/r  //////////////
	    //first term
	    coeff0 = 2.0*m_dt/rho/r0 * mu0 / dr;
	    coeff1 = 2.0*m_dt/rho/r0 * mu0 / dr;

	    compDiff_CellFace(&solver, I, I_nb, U_center, U_nb, flag,
		    0,0,4, -coeff0);
	    compDiff_CellFace(&solver, I, I_nb, U_center, U_nb, flag,
		    0,0,5,  coeff1);

	    //second term
	    coeff_temp = -m_dt/rho * mu0/rr;
	    rhs += coeff_temp*U_center[0];
	    solver.Add_A(I*3,I*3, -coeff_temp);

	    //third term
	    coeff0 = 2.0*m_dt/rho/rr * mu0 / dtheta;
	    coeff1 = 2.0*m_dt/rho/rr * mu0 / dtheta;

	    compDiff_CellFace(&solver, I, I_nb, U_center, U_nb, flag,
		    0,2,0, -coeff0);
	    compDiff_CellFace(&solver, I, I_nb, U_center, U_nb, flag,
		    0,2,1,  coeff1);

	    rhs += m_dt*state.m_U[0];
	    rhs += m_dt*cell_center[index].m_state.f_surf[0];
	    rhs -= m_dt*cell_center[index].m_state.grad_q[0]/rho;
	    rhs -= m_dt*cell_center[index].m_state.m_adv[0];

	    solver.Add_b(I*3, rhs);

	    /***************** Second Equation ***************/

	    solver.Add_A(I*3+1,I*3+1,1.0);
	    rhs = U_center[1];

	    ////////////////     term d(tao_zz)/dz   ///////////
	    //first term

	    if (flag[2] == 1)
	    {
		coeff_temp = 2.0*m_dt/rho * mu_edge[2]/(dz*dz);
		solver.Add_A(I*3+1,I*3+1, coeff_temp);
		rhs += 2.0*coeff_temp*U_nb[1][2];
		rhs -= coeff_temp*U_center[1];
	    }
	    else
	    {
		coeff_temp = m_dt/rho * mu_edge[2]/(dz*dz);
		solver.Add_A(I*3+1,I*3+1, coeff_temp);
		solver.Add_A(I*3+1,I_nb[2]*3+1, -coeff_temp);
		rhs += coeff_temp*U_nb[1][2];
		rhs -= coeff_temp*U_center[1];
	    }

	    if (flag[3] == 1)
	    {
		coeff_temp = 2.0*m_dt/rho * mu_edge[3]/(dz*dz);
		solver.Add_A(I*3+1,I*3+1, coeff_temp);
		rhs += 2.0*coeff_temp*U_nb[1][3];
		rhs -= coeff_temp*U_center[1];
	    }
	    else
	    {
		coeff_temp = m_dt/rho * mu_edge[3]/(dz*dz);
		solver.Add_A(I*3+1,I*3+1, coeff_temp);
		solver.Add_A(I*3+1,I_nb[3]*3+1, -coeff_temp);
		rhs += coeff_temp*U_nb[1][3];
		rhs -= coeff_temp*U_center[1];
	    }

	    //////////////   term d(tar_rz)/dr     /////////
	    //first term

	    coeff0 = m_dt/rho * mu_edge[4] / (dz*dr);
	    coeff1 = m_dt/rho * mu_edge[5] / (dz*dr);

	    compDiff_CellCorner(&solver, I, I_nb, U_center, U_nb, flag,
		    1,2,10,  coeff0);
	    compDiff_CellCorner(&solver, I, I_nb, U_center, U_nb, flag,
		    1,2,11, -coeff0);
	    compDiff_CellCorner(&solver, I, I_nb, U_center, U_nb, flag,
		    1,2,12,  coeff1);
	    compDiff_CellCorner(&solver, I, I_nb, U_center, U_nb, flag,
		    1,2,13, -coeff1);

	    //second term

	    if (flag[4] == 1)
	    {
		coeff_temp = m_dt/rho * mu_edge[4]/(dr*dr);
		solver.Add_A(I*3+1,I*3+1, coeff_temp);
		rhs += 2.0*coeff_temp*U_nb[1][4];
		rhs -= coeff_temp*U_center[1];
	    }
	    else
	    {
		coeff_temp = 0.5*m_dt/rho * mu_edge[4]/(dr*dr);
		solver.Add_A(I*3+1,I*3+1, coeff_temp);
		solver.Add_A(I*3+1,I_nb[4]*3+1, -coeff_temp);
		rhs += coeff_temp*U_nb[1][4];
		rhs -= coeff_temp*U_center[1];
	    }
	    
	    if (flag[5] == 1)
	    {
		coeff_temp = m_dt/rho * mu_edge[5]/(dr*dr);
		solver.Add_A(I*3+1,I*3+1, coeff_temp);
		rhs += 2.0*coeff_temp*U_nb[1][5];
		rhs -= coeff_temp*U_center[1];
	    }
	    else
	    {
		coeff_temp = 0.5*m_dt/rho * mu_edge[5]/(dr*dr);
		solver.Add_A(I*3+1,I*3+1, coeff_temp);
		solver.Add_A(I*3+1,I_nb[5]*3+1, -coeff_temp);
		rhs += coeff_temp*U_nb[1][5];
		rhs -= coeff_temp*U_center[1];
	    }

	    ///////////////    term  1/r * d tao_0z/d0    //////////////
	    //first term

	    coeff0 = m_dt/rho/r0 * mu_edge[0] / (dtheta*dz);
	    coeff1 = m_dt/rho/r0 * mu_edge[1] / (dtheta*dz);

	    compDiff_CellCorner(&solver, I, I_nb, U_center, U_nb, flag,
		    1,0,6,  coeff0);
	    compDiff_CellCorner(&solver, I, I_nb, U_center, U_nb, flag,
		    1,0,7, -coeff1);
	    compDiff_CellCorner(&solver, I, I_nb, U_center, U_nb, flag,
		    1,0,8,  coeff1);
	    compDiff_CellCorner(&solver, I, I_nb, U_center, U_nb, flag,
		    1,0,9, -coeff0);

	    //second term

	    if (flag[0] == 1)
	    {
		coeff_temp = m_dt/rho/r0 * mu_edge[0]/r_edge[0] / (dtheta*dtheta);
		solver.Add_A(I*3+1,I*3+1, coeff_temp);
		rhs += 2.0*coeff_temp*U_nb[1][0];
		rhs -= coeff_temp*U_center[1];
	    }
	    else
	    {
		coeff_temp = 0.5*m_dt/rho/r0 * mu_edge[0]/r_edge[0] / (dtheta*dtheta);
		solver.Add_A(I*3+1,I*3+1, coeff_temp);
		solver.Add_A(I*3+1,I_nb[0]*3+1, -coeff_temp);
		rhs += coeff_temp*U_nb[1][0];
		rhs -= coeff_temp*U_center[1];
	    }

	    if (flag[1] == 1)
	    {
		coeff_temp = m_dt/rho/r0 * mu_edge[1]/r_edge[1] / (dtheta*dtheta);
		solver.Add_A(I*3+1,I*3+1, coeff_temp);
		rhs += 2.0*coeff_temp*U_nb[1][1];
		rhs -= coeff_temp*U_center[1];
	    }
	    else
	    {
		coeff_temp = 0.5*m_dt/rho/r0 * mu_edge[1]/r_edge[1] / (dtheta*dtheta);
		solver.Add_A(I*3+1,I*3+1, coeff_temp);
		solver.Add_A(I*3+1,I_nb[1]*3+1, -coeff_temp);
		rhs += coeff_temp*U_nb[1][1];
		rhs -= coeff_temp*U_center[1];
	    }

	    //////////    term  tao_rz/r    ////////////
	    //first term

	    coeff0 = m_dt/rho/r0 * mu0 / dz;
	    coeff1 = m_dt/rho/r0 * mu0 / dz;

	    compDiff_CellFace(&solver, I, I_nb, U_center, U_nb, flag,
		    1,2,2, -coeff0);
	    compDiff_CellFace(&solver, I, I_nb, U_center, U_nb, flag,
		    1,2,3,  coeff1);

	    //second term

	    coeff0 = m_dt/rho/r0 * mu0 / dr;
	    coeff1 = m_dt/rho/r0 * mu0 / dr;

	    compDiff_CellFace(&solver, I, I_nb, U_center, U_nb, flag,
		    1,1,4, -coeff0);
	    compDiff_CellFace(&solver, I, I_nb, U_center, U_nb, flag,
		    1,1,5,  coeff1);

	    rhs += m_dt*state.m_U[1];
	    rhs += m_dt*cell_center[index].m_state.f_surf[1];
	    rhs -= m_dt*cell_center[index].m_state.grad_q[1]/rho;
	    rhs -= m_dt*cell_center[index].m_state.m_adv[1];

	    solver.Add_b(I*3+1, rhs);

	    /****************** Third Equation *****************/

	    solver.Add_A(I*3+2,I*3+2,1.0);
	    rhs = U_center[2];

	    /////////////    term d tao_zr / dz    //////////////////
	    //first term

	    if (flag[2] == 1)
	    {
		coeff_temp = m_dt/rho * mu_edge[2]/(dz*dz);
		solver.Add_A(I*3+2,I*3+2, coeff_temp);
		rhs += 2.0*coeff_temp*U_nb[2][2];
		rhs -= coeff_temp*U_center[2];
	    }
	    else
	    {
		coeff_temp = 0.5*m_dt/rho * mu_edge[2]/(dz*dz);
		solver.Add_A(I*3+2,I*3+2, coeff_temp);
		solver.Add_A(I*3+2,I_nb[2]*3+2, -coeff_temp);
		rhs += coeff_temp*U_nb[2][2];
		rhs -= coeff_temp*U_center[2];
	    }

	    if (flag[3] == 1)
	    {
		coeff_temp = m_dt/rho * mu_edge[3]/(dz*dz);
		solver.Add_A(I*3+2,I*3+2, coeff_temp);
		rhs += 2.0*coeff_temp*U_nb[2][3];
		rhs -= coeff_temp*U_center[2];
	    }
	    else
	    {
		coeff_temp = 0.5*m_dt/rho * mu_edge[3]/(dz*dz);
		solver.Add_A(I*3+2,I*3+2, coeff_temp);
		solver.Add_A(I*3+2,I_nb[3]*3+2, -coeff_temp);
		rhs += coeff_temp*U_nb[2][3];
		rhs -= coeff_temp*U_center[2];
	    }

	    //second term

	    coeff0 = m_dt/rho * mu_edge[2] / (dz*dr);
	    coeff1 = m_dt/rho * mu_edge[3] / (dz*dr);

	    compDiff_CellCorner(&solver, I, I_nb, U_center, U_nb, flag,
		    2,1,10,  coeff0);
	    compDiff_CellCorner(&solver, I, I_nb, U_center, U_nb, flag,
		    2,1,11, -coeff1);
	    compDiff_CellCorner(&solver, I, I_nb, U_center, U_nb, flag,
		    2,1,12,  coeff1);
	    compDiff_CellCorner(&solver, I, I_nb, U_center, U_nb, flag,
		    2,1,13, -coeff0);


	    ///////////////    term d tao_rr / dr   /////////////

	    if (flag[4] == 1)
	    {
		coeff_temp = 2.0*m_dt/rho * mu_edge[4]/(dr*dr);
		solver.Add_A(I*3+2,I*3+2, coeff_temp);
		rhs += 2.0*coeff_temp*U_nb[2][4];
		rhs -= coeff_temp*U_center[2];
	    }
	    else
	    {
		coeff_temp = m_dt/rho * mu_edge[4]/(dr*dr);
		solver.Add_A(I*3+2,I*3+2, coeff_temp);
		solver.Add_A(I*3+2,I_nb[4]*3+2, -coeff_temp);
		rhs += coeff_temp*U_nb[2][4];
		rhs -= coeff_temp*U_center[2];
	    }

	    if (flag[5] == 1)
	    {
		coeff_temp = 2.0*m_dt/rho * mu_edge[5]/(dr*dr);
		solver.Add_A(I*3+2,I*3+2, coeff_temp);
		rhs += 2.0*coeff_temp*U_nb[2][5];
		rhs -= coeff_temp*U_center[2];
	    }
	    else
	    {
		coeff_temp = m_dt/rho * mu_edge[5]/(dz*dz);
		solver.Add_A(I*3+2,I*3+2, coeff_temp);
		solver.Add_A(I*3+2,I_nb[5]*3+2, -coeff_temp);
		rhs += coeff_temp*U_nb[2][5];
		rhs -= coeff_temp*U_center[2];
	    }

	    ///////////     term 1/r * d tao_0r / d0   ///////////
	    //first term

	    coeff0 = m_dt/rho/r0 * mu_edge[0] / (dtheta*dr);
	    coeff1 = m_dt/rho/r0 * mu_edge[1] / (dtheta*dr);

	    compDiff_CellCorner(&solver, I, I_nb, U_center, U_nb, flag,
		    2,0,14,  coeff0);
	    compDiff_CellCorner(&solver, I, I_nb, U_center, U_nb, flag,
		    2,0,15, -coeff1);
	    compDiff_CellCorner(&solver, I, I_nb, U_center, U_nb, flag,
		    2,0,16,  coeff1);
	    compDiff_CellCorner(&solver, I, I_nb, U_center, U_nb, flag,
		    2,0,17, -coeff0);

	    //second term

	    coeff0 = m_dt/rho/r0 * mu_edge[0] / r_edge[0] / dtheta;
	    coeff1 = m_dt/rho/r0 * mu_edge[1] / r_edge[1] / dtheta;

	    compDiff_CellFace(&solver, I, I_nb, U_center, U_nb, flag,
		    2,0,0,  coeff0);
	    compDiff_CellFace(&solver, I, I_nb, U_center, U_nb, flag,
		    2,0,1, -coeff1);

	    //third term
	    
	    if (flag[0] == 1)
	    {
		coeff_temp = m_dt/rho/r0 * mu_edge[0]/r_edge[0] / (dtheta*dtheta);
		solver.Add_A(I*3+2,I*3+2, coeff_temp);
		rhs += 2.0*coeff_temp*U_nb[2][0];
		rhs -= coeff_temp*U_center[2];
	    }
	    else
	    {
		coeff_temp = 0.5*m_dt/rho/r0 * mu_edge[0]/r_edge[0] / (dtheta*dtheta);
		solver.Add_A(I*3+2,I*3+2, coeff_temp);
		solver.Add_A(I*3+2,I_nb[0]*3+2, -coeff_temp);
		rhs += coeff_temp*U_nb[2][0];
		rhs -= coeff_temp*U_center[2];
	    }

	    if (flag[1] == 1)
	    {
		coeff_temp = m_dt/rho/r0 * mu_edge[1]/r_edge[1] / (dtheta*dtheta);
		solver.Add_A(I*3+2,I*3+2, coeff_temp);
		rhs += 2.0*coeff_temp*U_nb[2][1];
		rhs -= coeff_temp*U_center[2];
	    }
	    else
	    {
		coeff_temp = 0.5*m_dt/rho/r0 * mu_edge[1]/r_edge[1] / (dtheta*dtheta);
		solver.Add_A(I*3+2,I*3+2, coeff_temp);
		solver.Add_A(I*3+2,I_nb[1]*3+2, -coeff_temp);
		rhs += coeff_temp*U_nb[2][1];
		rhs -= coeff_temp*U_center[2];
	    }

	    /////////////////    term  tao_rr / r    /////////////////////

	    coeff0 = 2.0*m_dt/rho * mu0/r0 / dr;
	    coeff1 = 2.0*m_dt/rho * mu0/r0 / dr;

	    compDiff_CellFace(&solver, I, I_nb, U_center, U_nb, flag,
		    2,2,4,  -coeff0);
	    compDiff_CellFace(&solver, I, I_nb, U_center, U_nb, flag,
		    2,2,5,   coeff1);


	    //////////////////   term  -tao_00 / r    ///////////////////
	    //first term

	    coeff0 = 2.0*m_dt/rho * mu0/rr / dtheta;
	    coeff1 = 2.0*m_dt/rho * mu0/rr / dtheta;

	    compDiff_CellFace(&solver, I, I_nb, U_center, U_nb, flag,
		    2,0,0,   coeff0);
	    compDiff_CellFace(&solver, I, I_nb, U_center, U_nb, flag,
		    2,0,1,  -coeff1);

	    //second term

	    coeff_temp = -m_dt/rho * mu0/rr;
	    rhs += coeff_temp*U_center[2];
	    solver.Add_A(I*3+2,I*3+2, -coeff_temp);


	    rhs += m_dt*state.m_U[2];
	    rhs += m_dt*cell_center[index].m_state.f_surf[2];
	    rhs -= m_dt*cell_center[index].m_state.grad_q[2]/rho;
	    rhs -= m_dt*cell_center[index].m_state.m_adv[2];

	    solver.Add_b(I*3+2, rhs);
        }

        solver.SetMaxIter(40000);
        solver.SetTol(1e-15);

	start_clock("Before Petsc Solve");
        solver.Solve_GMRES();
	solver.GetNumIterations(&num_iter);
	solver.GetFinalRelativeResidualNorm(&rel_residual);
	stop_clock("After Petsc Solve");

        // get back the solution
        solver.Get_x(x);

        if (debugging("PETSc"))
            (void) printf("Incompress_Solver_Smooth_3D_Cylindrical::compDiffWithSmoothProperty_2nd_coupled: "
                        "num_iter = %d, rel_residual = %g. \n",
                        num_iter,rel_residual);


	for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I = ijk_to_I[i][j][k];
            index = d_index3d(i,j,k,top_gmax);
            if (I >= 0)
            {
                cell_center[index].m_state.m_U[0] = x[I*3-ilower*3];
                cell_center[index].m_state.m_U[1] = x[I*3-ilower*3+1];
		cell_center[index].m_state.m_U[2] = x[I*3-ilower*3+2];
                speed = fabs(cell_center[index].m_state.m_U[0]) +
                        fabs(cell_center[index].m_state.m_U[1]) +
			fabs(cell_center[index].m_state.m_U[2]);
                if (speed > max_speed)
                    max_speed = speed;
            }
            else
            {
                cell_center[index].m_state.m_U[0] = 0.0;
                cell_center[index].m_state.m_U[1] = 0.0;
		cell_center[index].m_state.m_U[2] = 0.0;
            }
        }

        for (l = 0; l < 3; ++l)
        {
	    for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)

            {
                index  = d_index3d(i,j,k,top_gmax);
                array[index] = cell_center[index].m_state.m_U[l];
            }
            scatMeshArray();
	    for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)

            {
                index  = d_index3d(i,j,k,top_gmax);
                cell_center[index].m_state.m_U[l] = array[index];
            }
        }

        pp_global_max(&max_speed,1);

        free_these(1,x);
}       /* end 3D_Cylindrical::compDiffWithSmoothProperty_2nd_coupled */


//-------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

void Incompress_Solver_Smooth_3D_Cylindrical::computeAdvection_test(void)
{
	int i,j,k,l;
	int index,index00,index01,index10,index11,index20,index21,size;
	L_STATE state;
	COMPONENT comp;
	double speed;
	double *u, *v, *w;
	double u0,u00,u01,u10,u11,u20,u21;
	double v0,v00,v01,v10,v11,v20,v21;
	double w0,w00,w01,w10,w11,w20,w21;
	double crx_coords[MAXD];
	int icoords[MAXD];
	POINTER intfc_state;
	HYPER_SURF *hs;

	double rho, mu0;

	max_speed = 0.0;

	size = (top_gmax[0]+1)*(top_gmax[1]+1)*(top_gmax[2]+1);
	FT_VectorMemoryAlloc((POINTER*)&u,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&v,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&w,size,sizeof(double));

	for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    u[index] = cell_center[index].m_state.m_U[0];
	    v[index] = cell_center[index].m_state.m_U[1];
	    w[index] = cell_center[index].m_state.m_U[2];
	}

	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{	
	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;
	    index  = d_index3d(i,j,k,top_gmax);
	    comp = top_comp[index];
	    if (comp == SOLID_COMP)
	    {
	    	cell_center[index].m_state.m_U[0] = 0.0; 
	    	cell_center[index].m_state.m_U[1] = 0.0; 
	    	cell_center[index].m_state.m_U[2] = 0.0; 
		continue;
	    }
	    u0 = u[index];
	    v0 = v[index];
	    w0 = w[index];

	    rho = cell_center[index].m_state.m_rho;
	    mu0 = cell_center[index].m_state.m_mu;


	    // To continue
	    index00 = d_index3d(i-1,j,k,top_gmax);
	    if (FT_StateStructAtGridCrossing(front,icoords,WEST,comp,
			        &intfc_state,&hs,crx_coords) &&
		                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    {
		u00 = getStateXvel(intfc_state);
		v00 = getStateYvel(intfc_state);
		w00 = getStateZvel(intfc_state);
	    }
	    else
	    {
		u00 = u[index00];
		v00 = v[index00];
		w00 = w[index00];
	    }
	    index01 = d_index3d(i+1,j,k,top_gmax);
	    if (FT_StateStructAtGridCrossing(front,icoords,EAST,comp,
			        &intfc_state,&hs,crx_coords) &&
		                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    {
		u01 = getStateXvel(intfc_state);
		v01 = getStateYvel(intfc_state);
		w01 = getStateZvel(intfc_state);
	    }
	    else
	    {
		u01 = u[index01];
		v01 = v[index01];
		w01 = w[index01];
	    }
	    index10 = d_index3d(i,j-1,k,top_gmax);
	    if (FT_StateStructAtGridCrossing(front,icoords,SOUTH,comp,
			        &intfc_state,&hs,crx_coords) &&
		                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    {
		u10 = getStateXvel(intfc_state);
		v10 = getStateYvel(intfc_state);
		w10 = getStateZvel(intfc_state);
	    }
	    else
	    {
		u10 = u[index10];
		v10 = v[index10];
		w10 = w[index10];
	    }
	    index11 = d_index3d(i,j+1,k,top_gmax);
	    if (FT_StateStructAtGridCrossing(front,icoords,NORTH,comp,
			        &intfc_state,&hs,crx_coords) &&
		                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    {
		u11 = getStateXvel(intfc_state);
		v11 = getStateYvel(intfc_state);
		w11 = getStateZvel(intfc_state);
	    }
	    else
	    {
		u11 = u[index11];
		v11 = v[index11];
		w11 = w[index11];
	    }
	    index20 = d_index3d(i,j,k-1,top_gmax);
	    if (FT_StateStructAtGridCrossing(front,icoords,LOWER,comp,
			        &intfc_state,&hs,crx_coords) &&
		                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    {
		u20 = getStateXvel(intfc_state);
		v20 = getStateYvel(intfc_state);
		w20 = getStateZvel(intfc_state);
	    }
	    else
	    {
		u20 = u[index20];
		v20 = v[index20];
		w20 = w[index20];
	    }
	    index21 = d_index3d(i,j,k+1,top_gmax);
	    if (FT_StateStructAtGridCrossing(front,icoords,UPPER,comp,
			        &intfc_state,&hs,crx_coords) &&
		                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    {
		u21 = getStateXvel(intfc_state);
		v21 = getStateYvel(intfc_state);
		w21 = getStateZvel(intfc_state);
	    }
	    else
	    {
		u21 = u[index21];
		v21 = v[index21];
		w21 = w[index21];
	    }


	    /*
	    cell_center[index].m_state.m_U[0] += -m_dt*(
				burger_flux(u00,u0,u01)/top_h[0] +
				linear_flux(v0,u10,u0,u11)/top_h[1] +
				linear_flux(w0,u20,u0,u21)/top_h[2]);

	    cell_center[index].m_state.m_U[1] += - m_dt*(
				linear_flux(u0,v00,v0,v01)/top_h[0] +
			 	burger_flux(v10,v0,v11)/top_h[1] +
				linear_flux(w0,v20,v0,v21)/top_h[2]);

	    cell_center[index].m_state.m_U[2] += - m_dt*(
				linear_flux(u0,w00,w0,w01)/top_h[0] +
				linear_flux(v0,w10,w0,w11)/top_h[1] +
			 	burger_flux(w20,w0,w21)/top_h[2]);
	    */

	    double r = cell_center[index].m_coords[2];

            cell_center[index].m_state.m_U[0] += - m_dt*(
                                burger_flux(u00,u0,u01)/(r*top_h[0]) +
                                linear_flux(v0,u10,u0,u11)/top_h[1] +
                                linear_flux(w0,u20,u0,u21)/top_h[2]);

            cell_center[index].m_state.m_U[1] += - m_dt*(
                                linear_flux(u0,v00,v0,v01)/(r*top_h[0]) +
                                burger_flux(v10,v0,v11)/top_h[1] +
                                linear_flux(w0,v20,v0,v21)/top_h[2]);

            cell_center[index].m_state.m_U[2] += - m_dt*(
                                linear_flux(u0,w00,w0,w01)/(r*top_h[0]) +
                                linear_flux(v0,w10,w0,w11)/top_h[1] +
                                burger_flux(w20,w0,w21)/top_h[2]);
          
            //SOURCE TERM
	    //
	    //How about calculating it in the diffusion solver?
	    
            cell_center[index].m_state.m_U[0] += -(m_dt*u0*w0/r);

            cell_center[index].m_state.m_U[2] += (m_dt*u0*u0/r);

	    // First derivative term in the diffusion solver

	    cell_center[index].m_state.m_U[0] += mu0/rho * 2.0/(r*r) * (w01 - w00)/ (2.0*top_h[0]);
	    cell_center[index].m_state.m_U[2] -= mu0/rho * 2.0/(r*r) * (u01 - u00)/ (2.0*top_h[0]);


	    
	    speed = fabs(cell_center[index].m_state.m_U[0]) +
		    fabs(cell_center[index].m_state.m_U[1]) +
		    fabs(cell_center[index].m_state.m_U[2]);
	    if (speed > max_speed)
		max_speed = speed;
	}
	for (l = 0; l < 3; ++l)
	{
	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {	
	    	index  = d_index3d(i,j,k,top_gmax);
	    	array[index] = cell_center[index].m_state.m_U[l];
	    }
	    scatMeshArray();
	    for (k = 0; k <= top_gmax[2]; k++)
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {	
	    	index  = d_index3d(i,j,k,top_gmax);
	    	cell_center[index].m_state.m_U[l] = array[index];
	    }
	}
	pp_global_max(&max_speed,1);
	free_these(3,u,v,w);
/*
	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
	for (i = imin; i <= imax; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    //printf("After advection, u_theta = %.16g, u_z = %.16g, u_r = %.16g\n",cell_center[index].m_state.m_U[0], cell_center[index].m_state.m_U[1], cell_center[index].m_state.m_U[2]);
	}
	*/
}	/* end computeAdvection3d */



void Incompress_Solver_Smooth_3D_Cylindrical::compDiffWithSmoothProperty_1st_decoupled_test(void)
{
        COMPONENT comp;
        int index,index_nb[6],size;
        int I,I_nb[6];
        double coords[MAXD],crx_coords[MAXD];
	double coeff[18],mu0,rho,rhs,U0_nb[6],U1_nb[6],U2_nb[6],U0_center,U1_center,U2_center;
        L_STATE state;
        int i,j,k,l,nb,icoords[MAXD];
        double speed;
	double *x;
	GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
	POINTER intfc_state;
	HYPER_SURF *hs;
	int num_iter;
	double rel_residual;
	double r, rr, redge[2],rnb[2];

	max_speed = 0.0;
	setIndexMap();

	size = iupper - ilower;
	FT_VectorMemoryAlloc((POINTER*)&x,3*size,sizeof(double));

        PETSc solver;
        solver.Create(3*ilower, 3*iupper-1, 7, 7);
	// 7theta + 2r  for the first equation
	// 7z for the second equation
	// 7r + 2theta for the third equation

	solver.Reset_A();
	solver.Reset_b();
	solver.Reset_x();

	printf("\nIn diffusion solver ,m_dt = %.16g\n", m_dt);

	for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I  = ijk_to_I[i][j][k];
            if (I == -1) continue;

            index  = d_index3d(i,j,k,top_gmax);	
	//6 neighbours of the center cell
            index_nb[0] = d_index3d(i-1,j,k,top_gmax);
            index_nb[1] = d_index3d(i+1,j,k,top_gmax);
            index_nb[2] = d_index3d(i,j-1,k,top_gmax);
            index_nb[3] = d_index3d(i,j+1,k,top_gmax);
	    index_nb[4] = d_index3d(i,j,k-1,top_gmax);
	    index_nb[5] = d_index3d(i,j,k+1,top_gmax);

        
	//6 neighbours of the center cell
            I_nb[0] = ijk_to_I[i-1][j][k];
            I_nb[1] = ijk_to_I[i+1][j][k];
            I_nb[2] = ijk_to_I[i][j-1][k];
            I_nb[3] = ijk_to_I[i][j+1][k];
            I_nb[4] = ijk_to_I[i][j][k-1];
            I_nb[5] = ijk_to_I[i][j][k+1];
	

	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;
	    comp = top_comp[index];

	    mu0 = cell_center[index].m_state.m_mu;
	    rho = cell_center[index].m_state.m_rho;
	    U0_center = cell_center[index].m_state.m_U[0];
	    U1_center = cell_center[index].m_state.m_U[1];
	    U2_center = cell_center[index].m_state.m_U[2];

	  
            for (nb = 0; nb < 6; nb++)
            {
                if (FT_StateStructAtGridCrossing(front,icoords,dir[nb],
                            comp,&intfc_state,&hs,crx_coords) &&
                            wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	       	{
		    U0_nb[nb] = getStateVel[0](intfc_state);
		    U1_nb[nb] = getStateVel[1](intfc_state);
		    U2_nb[nb] = getStateVel[2](intfc_state);
		}
                else
		{
		    U0_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[0];
		    U1_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[1];
		    U2_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[2];
		}
	    }

            getRectangleCenter(index, coords);
            computeSourceTerm(coords, state);


    //Setting the coeffecients for the first equation

	    rr = cell_center[index].m_coords[2] * cell_center[index].m_coords[2];
	    r = cell_center[index].m_coords[2];
	    redge[0] = 0.5 * (cell_center[index].m_coords[2] + cell_center[index_nb[4]].m_coords[2]);
	    redge[1] = 0.5 * (cell_center[index].m_coords[2] + cell_center[index_nb[5]].m_coords[2]);
	    rnb[0] = cell_center[index_nb[4]].m_coords[2];
	    rnb[1] = cell_center[index_nb[5]].m_coords[2];


	    coeff[0] = 0.5*m_dt/rho * mu0/(rr*top_h[0]*top_h[0]);
	    coeff[1] = 0.5*m_dt/rho * mu0/(rr*top_h[0]*top_h[0]);
	    coeff[2] = 0.5*m_dt/rho * mu0/(top_h[1]*top_h[1]);
	    coeff[3] = 0.5*m_dt/rho * mu0/(top_h[1]*top_h[1]);

	    //coeff[4] = 0.5*m_dt/rho * mu0*redge[0]/(r*top_h[2]*top_h[2]);
	    //coeff[5] = 0.5*m_dt/rho * mu0*redge[1]/(r*top_h[2]*top_h[2]);

	    coeff[4] = 0.5*m_dt/rho * mu0 / (redge[0]*top_h[2]*top_h[2]);
	    coeff[5] = 0.5*m_dt/rho * mu0 / (redge[1]*top_h[2]*top_h[2]);

	    for (nb = 0; nb < 6; nb++) //For the second order boundary
	    {
		if (I_nb[nb] == -1)
		{
		    //printf("\n i = %d, j = %d, k = %d, neighbour = %d\n",i,j,k,nb);
		    coeff[nb] = 2.0*coeff[nb];
		}
	    }

	    solver.Set_A(I*3,I*3,1.0+coeff[0]+coeff[1]+coeff[2]+coeff[3]+coeff[4]*r+coeff[5]*r);
	    rhs = (1.0-coeff[0]-coeff[1]-coeff[2]-coeff[3]-coeff[4]*r-coeff[5]*r)*cell_center[index].m_state.m_U[0];

	    if (I_nb[4] == -1)
		coeff[4] = coeff[4] * redge[0];
	    else
		coeff[4] = coeff[4] * rnb[0];

	    if (I_nb[5] == -1)
		coeff[5] = coeff[5] * redge[1];
	    else
		coeff[5] = coeff[5] * rnb[1];

	    for (nb = 0; nb < 6; nb++)
	    {
		if (I_nb[nb] != -1)
		{
		    solver.Set_A(I*3,I_nb[nb]*3,-coeff[nb]);
		    rhs += coeff[nb]*U0_nb[nb];
		}
		else
		    rhs += 2.0*coeff[nb]*U0_nb[nb];
	    }


	    //rhs -= m_dt*cell_center[index].m_state.m_U[0]*cell_center[index].m_state.m_U[2]/r; //Source term in advection step
	    rhs += m_dt*state.m_U[0];
	    rhs += m_dt*cell_center[index].m_state.f_surf[0];
	    rhs -= m_dt*cell_center[index].m_state.grad_q[0]/rho;

	    solver.Set_b(I*3, rhs);

	    /************************************************************************/


    //Setting the coeffecients for the second equation


	    coeff[0] = 0.5*m_dt/rho * mu0/(rr*top_h[0]*top_h[0]);
	    coeff[1] = 0.5*m_dt/rho * mu0/(rr*top_h[0]*top_h[0]);
	    coeff[2] = 0.5*m_dt/rho * mu0/(top_h[1]*top_h[1]);
	    coeff[3] = 0.5*m_dt/rho * mu0/(top_h[1]*top_h[1]);

	    coeff[4] = 0.5*m_dt/rho * mu0*redge[0]/(r*top_h[2]*top_h[2]);
	    coeff[5] = 0.5*m_dt/rho * mu0*redge[1]/(r*top_h[2]*top_h[2]);

	    for (nb = 0; nb < 6; nb++)
	    {
		if (I_nb[nb] == -1)
		    coeff[nb] = 2.0*coeff[nb];
	    }

	    solver.Set_A(I*3+1,I*3+1,1.0+coeff[0]+coeff[1]+coeff[2]+coeff[3]+coeff[4]+coeff[5]);
	    rhs = (1.0-coeff[0]-coeff[1]-coeff[2]-coeff[3]-coeff[4]-coeff[5])*cell_center[index].m_state.m_U[1];


	    for (nb = 0; nb < 6; nb++)
	    {
		if (I_nb[nb] != -1)
		{
		    solver.Set_A(I*3+1,I_nb[nb]*3+1,-coeff[nb]);
		    rhs += coeff[nb]*U1_nb[nb];
		}
		else
		    rhs += 2.0*coeff[nb]*U1_nb[nb];
	    }

	    rhs += m_dt*state.m_U[1];
	    rhs += m_dt*cell_center[index].m_state.f_surf[1];
	    rhs -= m_dt*cell_center[index].m_state.grad_q[1]/rho;

	    solver.Set_b(I*3+1, rhs);

	    /************************************************************************/

    //Setting the coeffecients for the third equation

	    coeff[0] = 0.5*m_dt/rho * mu0/(rr*top_h[0]*top_h[0]);
	    coeff[1] = 0.5*m_dt/rho * mu0/(rr*top_h[0]*top_h[0]);
	    coeff[2] = 0.5*m_dt/rho * mu0/(top_h[1]*top_h[1]);
	    coeff[3] = 0.5*m_dt/rho * mu0/(top_h[1]*top_h[1]);

	    coeff[4] = 0.5*m_dt/rho * mu0/(redge[0]*top_h[2]*top_h[2]);
	    coeff[5] = 0.5*m_dt/rho * mu0/(redge[1]*top_h[2]*top_h[2]);

	    for (nb = 0; nb < 6; nb++)
	    {
		if (I_nb[nb] == -1)
		    coeff[nb] = 2.0*coeff[nb];
	    }

	    solver.Set_A(I*3+2,I*3+2,1.0+coeff[0]+coeff[1]+coeff[2]+coeff[3]+coeff[4]*r+coeff[5]*r);
	    rhs = (1.0-coeff[0]-coeff[1]-coeff[2]-coeff[3]-coeff[4]*r-coeff[5]*r)*cell_center[index].m_state.m_U[2];

	    if (I_nb[4] == -1)
		coeff[4] = coeff[4] * redge[0];
	    else
		coeff[4] = coeff[4] * rnb[0];

	    if (I_nb[5] == -1)
		coeff[5] = coeff[5] * redge[1];
	    else
		coeff[5] = coeff[5] * rnb[1];

	    
	    for (nb = 0; nb < 6; nb++)
	    {
		if (I_nb[nb] != -1)
		{
		    solver.Set_A(I*3+2,I_nb[nb]*3+2,-coeff[nb]);
		    rhs += coeff[nb]*U2_nb[nb];
		}
		else
		    rhs += 2.0*coeff[nb]*U2_nb[nb];
	    }

	    //rhs += m_dt* cell_center[index].m_state.m_U[0] * cell_center[index].m_state.m_U[0] / r; //Source term in advection step
	    rhs += m_dt*state.m_U[2];
	    rhs += m_dt*cell_center[index].m_state.f_surf[2];
	    rhs -= m_dt*cell_center[index].m_state.grad_q[2]/rho;

	    solver.Set_b(I*3+2, rhs);
        }

        solver.SetMaxIter(40000);
        solver.SetTol(1e-12);

	start_clock("Before Petsc Solve");
        solver.Solve_GMRES();
	solver.GetNumIterations(&num_iter);
	solver.GetFinalRelativeResidualNorm(&rel_residual);
	stop_clock("After Petsc Solve");

        // get back the solution
        solver.Get_x(x);

        if (debugging("PETSc"))
            (void) printf("Incompress_Solver_Smooth_3D_Cylindrical::compDiffWithSmoothProperty_1st_decoupled_test: "
                        "num_iter = %d, rel_residual = %g. \n",
                        num_iter,rel_residual);


	for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I = ijk_to_I[i][j][k];
            index = d_index3d(i,j,k,top_gmax);
            if (I >= 0)
            {
                cell_center[index].m_state.m_U[0] = x[I*3-ilower*3];
                cell_center[index].m_state.m_U[1] = x[I*3-ilower*3+1];
		cell_center[index].m_state.m_U[2] = x[I*3-ilower*3+2];
                speed = fabs(cell_center[index].m_state.m_U[0]) +
                        fabs(cell_center[index].m_state.m_U[1]) +
			fabs(cell_center[index].m_state.m_U[2]);
                if (speed > max_speed)
                    max_speed = speed;
            }
            else
            {
                cell_center[index].m_state.m_U[0] = 0.0;
                cell_center[index].m_state.m_U[1] = 0.0;
		cell_center[index].m_state.m_U[2] = 0.0;
            }
        }
        for (l = 0; l < 3; ++l)
        {
	    for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)

            {
                index  = d_index3d(i,j,k,top_gmax);
                array[index] = cell_center[index].m_state.m_U[l];
            }
            scatMeshArray();
	    for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)

            {
                index  = d_index3d(i,j,k,top_gmax);
                cell_center[index].m_state.m_U[l] = array[index];
            }
        }
        pp_global_max(&max_speed,1);

        free_these(1,x);
	/*
	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
	for (i = imin; i <= imax; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    //printf("After diffusion, u_theta = %.16g, u_z = %.16g, u_r = %.16g\n",cell_center[index].m_state.m_U[0], cell_center[index].m_state.m_U[1], cell_center[index].m_state.m_U[2]);
	}
*/
}       /* end compDiffWithSmoothProperty3d */


void Incompress_Solver_Smooth_3D_Cylindrical::computeAdvection(void)
{
	int i,j,k,l;
	int index,index00,index01,index10,index11,index20,index21,size;
	L_STATE state;
	COMPONENT comp;
	double speed;
	double *u, *v, *w;
	double u0,u00,u01,u10,u11,u20,u21;
	double v0,v00,v01,v10,v11,v20,v21;
	double w0,w00,w01,w10,w11,w20,w21;
	double crx_coords[MAXD];
	int icoords[MAXD];
	POINTER intfc_state;
	HYPER_SURF *hs;

	max_speed = 0.0;

	size = (top_gmax[0]+1)*(top_gmax[1]+1)*(top_gmax[2]+1);
	FT_VectorMemoryAlloc((POINTER*)&u,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&v,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&w,size,sizeof(double));

	for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    u[index] = cell_center[index].m_state.m_U[0];
	    v[index] = cell_center[index].m_state.m_U[1];
	    w[index] = cell_center[index].m_state.m_U[2];
	}

	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{	
	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;
	    index  = d_index3d(i,j,k,top_gmax);
	    comp = top_comp[index];
	    if (comp == SOLID_COMP)
	    {
	    	cell_center[index].m_state.m_U[0] = 0.0; 
	    	cell_center[index].m_state.m_U[1] = 0.0; 
	    	cell_center[index].m_state.m_U[2] = 0.0; 
		continue;
	    }
	    u0 = u[index];
	    v0 = v[index];
	    w0 = w[index];
	    // To continue
	    index00 = d_index3d(i-1,j,k,top_gmax);
	    if (FT_StateStructAtGridCrossing(front,icoords,WEST,comp,
			        &intfc_state,&hs,crx_coords) &&
		                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    {
		u00 = getStateXvel(intfc_state);
		v00 = getStateYvel(intfc_state);
		w00 = getStateZvel(intfc_state);
	    }
	    else
	    {
		u00 = u[index00];
		v00 = v[index00];
		w00 = w[index00];
	    }
	    index01 = d_index3d(i+1,j,k,top_gmax);
	    if (FT_StateStructAtGridCrossing(front,icoords,EAST,comp,
			        &intfc_state,&hs,crx_coords) &&
		                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    {
		u01 = getStateXvel(intfc_state);
		v01 = getStateYvel(intfc_state);
		w01 = getStateZvel(intfc_state);
	    }
	    else
	    {
		u01 = u[index01];
		v01 = v[index01];
		w01 = w[index01];
	    }
	    index10 = d_index3d(i,j-1,k,top_gmax);
	    if (FT_StateStructAtGridCrossing(front,icoords,SOUTH,comp,
			        &intfc_state,&hs,crx_coords) &&
		                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    {
		u10 = getStateXvel(intfc_state);
		v10 = getStateYvel(intfc_state);
		w10 = getStateZvel(intfc_state);
	    }
	    else
	    {
		u10 = u[index10];
		v10 = v[index10];
		w10 = w[index10];
	    }
	    index11 = d_index3d(i,j+1,k,top_gmax);
	    if (FT_StateStructAtGridCrossing(front,icoords,NORTH,comp,
			        &intfc_state,&hs,crx_coords) &&
		                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    {
		u11 = getStateXvel(intfc_state);
		v11 = getStateYvel(intfc_state);
		w11 = getStateZvel(intfc_state);
	    }
	    else
	    {
		u11 = u[index11];
		v11 = v[index11];
		w11 = w[index11];
	    }
	    index20 = d_index3d(i,j,k-1,top_gmax);
	    if (FT_StateStructAtGridCrossing(front,icoords,LOWER,comp,
			        &intfc_state,&hs,crx_coords) &&
		                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    {
		u20 = getStateXvel(intfc_state);
		v20 = getStateYvel(intfc_state);
		w20 = getStateZvel(intfc_state);
	    }
	    else
	    {
		u20 = u[index20];
		v20 = v[index20];
		w20 = w[index20];
	    }
	    index21 = d_index3d(i,j,k+1,top_gmax);
	    if (FT_StateStructAtGridCrossing(front,icoords,UPPER,comp,
			        &intfc_state,&hs,crx_coords) &&
		                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    {
		u21 = getStateXvel(intfc_state);
		v21 = getStateYvel(intfc_state);
		w21 = getStateZvel(intfc_state);
	    }
	    else
	    {
		u21 = u[index21];
		v21 = v[index21];
		w21 = w[index21];
	    }


	    /*
	    cell_center[index].m_state.m_U[0] += -m_dt*(
				burger_flux(u00,u0,u01)/top_h[0] +
				linear_flux(v0,u10,u0,u11)/top_h[1] +
				linear_flux(w0,u20,u0,u21)/top_h[2]);

	    cell_center[index].m_state.m_U[1] += - m_dt*(
				linear_flux(u0,v00,v0,v01)/top_h[0] +
			 	burger_flux(v10,v0,v11)/top_h[1] +
				linear_flux(w0,v20,v0,v21)/top_h[2]);

	    cell_center[index].m_state.m_U[2] += - m_dt*(
				linear_flux(u0,w00,w0,w01)/top_h[0] +
				linear_flux(v0,w10,w0,w11)/top_h[1] +
			 	burger_flux(w20,w0,w21)/top_h[2]);
	    */

	    double r = cell_center[index].m_coords[2];

            cell_center[index].m_state.m_U[0] += - m_dt*(
                                burger_flux(u00,u0,u01)/(r*top_h[0]) +
                                linear_flux(v0,u10,u0,u11)/top_h[1] +
                                linear_flux(w0,u20,u0,u21)/top_h[2]);

            cell_center[index].m_state.m_U[1] += - m_dt*(
                                linear_flux(u0,v00,v0,v01)/(r*top_h[0]) +
                                burger_flux(v10,v0,v11)/top_h[1] +
                                linear_flux(w0,v20,v0,v21)/top_h[2]);

            cell_center[index].m_state.m_U[2] += - m_dt*(
                                linear_flux(u0,w00,w0,w01)/(r*top_h[0]) +
                                linear_flux(v0,w10,w0,w11)/top_h[1] +
                                burger_flux(w20,w0,w21)/top_h[2]);
          
            //SOURCE TERM
	    //
	    //How about calculating it in the diffusion solver?
	   /* 
            cell_center[index].m_state.m_U[0] += -(m_dt*u0*w0/r);

            cell_center[index].m_state.m_U[2] += (m_dt*u0*u0/r);
*/
	    
	    speed = fabs(cell_center[index].m_state.m_U[0]) +
		    fabs(cell_center[index].m_state.m_U[1]) +
		    fabs(cell_center[index].m_state.m_U[2]);
	    if (speed > max_speed)
		max_speed = speed;
	}
	for (l = 0; l < 3; ++l)
	{
	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {	
	    	index  = d_index3d(i,j,k,top_gmax);
	    	array[index] = cell_center[index].m_state.m_U[l];
	    }
	    scatMeshArray();
	    for (k = 0; k <= top_gmax[2]; k++)
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {	
	    	index  = d_index3d(i,j,k,top_gmax);
	    	cell_center[index].m_state.m_U[l] = array[index];
	    }
	}
	pp_global_max(&max_speed,1);
	free_these(3,u,v,w);

}	/* end computeAdvection3d */



void Incompress_Solver_Smooth_3D_Cylindrical::compDiffWithSmoothProperty_1st_decoupled(void)
{
        COMPONENT comp;
        int index,index_nb[6],size;
        int I,I_nb[6];
        double coords[MAXD],crx_coords[MAXD];
	double coeff[18],mu0,rho,rhs,U0_nb[6],U1_nb[6],U2_nb[6],U0_center,U1_center,U2_center;
        L_STATE state;
        int i,j,k,l,nb,icoords[MAXD];
        double speed;
	double *x;
	GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
	POINTER intfc_state;
	HYPER_SURF *hs;
	int num_iter;
	double rel_residual;
	double r, rr, redge[2];

	max_speed = 0.0;
	setIndexMap();

	size = iupper - ilower;
	FT_VectorMemoryAlloc((POINTER*)&x,3*size,sizeof(double));

        PETSc solver;
        solver.Create(3*ilower, 3*iupper-1, 9, 9);
	// 7theta + 2r  for the first equation
	// 7z for the second equation
	// 7r + 2theta for the third equation

	solver.Reset_A();
	solver.Reset_b();
	solver.Reset_x();

	printf("\nIn diffusion solver ,m_dt = %.16g\n", m_dt);

	for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I  = ijk_to_I[i][j][k];
            if (I == -1) continue;

            index  = d_index3d(i,j,k,top_gmax);	
	//6 neighbours of the center cell
            index_nb[0] = d_index3d(i-1,j,k,top_gmax);
            index_nb[1] = d_index3d(i+1,j,k,top_gmax);
            index_nb[2] = d_index3d(i,j-1,k,top_gmax);
            index_nb[3] = d_index3d(i,j+1,k,top_gmax);
	    index_nb[4] = d_index3d(i,j,k-1,top_gmax);
	    index_nb[5] = d_index3d(i,j,k+1,top_gmax);

        
	//6 neighbours of the center cell
            I_nb[0] = ijk_to_I[i-1][j][k];
            I_nb[1] = ijk_to_I[i+1][j][k];
            I_nb[2] = ijk_to_I[i][j-1][k];
            I_nb[3] = ijk_to_I[i][j+1][k];
            I_nb[4] = ijk_to_I[i][j][k-1];
            I_nb[5] = ijk_to_I[i][j][k+1];

	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;
	    comp = top_comp[index];

	    mu0 = cell_center[index].m_state.m_mu;
	    rho = cell_center[index].m_state.m_rho;
	    U0_center = cell_center[index].m_state.m_U[0];
	    U1_center = cell_center[index].m_state.m_U[1];
	    U2_center = cell_center[index].m_state.m_U[2];

	  
            for (nb = 0; nb < 6; nb++)
            {
                if (FT_StateStructAtGridCrossing(front,icoords,dir[nb],
                            comp,&intfc_state,&hs,crx_coords) &&
                            wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	       	{
		    U0_nb[nb] = getStateVel[0](intfc_state);
		    U1_nb[nb] = getStateVel[1](intfc_state);
		    U2_nb[nb] = getStateVel[2](intfc_state);
		}
                else
		{
		    U0_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[0];
		    U1_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[1];
		    U2_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[2];
		}
	    }

            getRectangleCenter(index, coords);
            computeSourceTerm(coords, state);


    //Setting the coeffecients for the first equation

	    rr = cell_center[index].m_coords[2] * cell_center[index].m_coords[2];
	    r = cell_center[index].m_coords[2];
	    redge[0] = 0.5 * (cell_center[index].m_coords[2] + cell_center[index_nb[4]].m_coords[2]);
	    redge[1] = 0.5 * (cell_center[index].m_coords[2] + cell_center[index_nb[5]].m_coords[2]);


	    coeff[0] = 0.5*m_dt/rho * mu0/(rr*top_h[0]*top_h[0]);
	    coeff[1] = 0.5*m_dt/rho * mu0/(rr*top_h[0]*top_h[0]);
	    coeff[2] = 0.5*m_dt/rho * mu0/(top_h[1]*top_h[1]);
	    coeff[3] = 0.5*m_dt/rho * mu0/(top_h[1]*top_h[1]);

	    coeff[4] = 0.5*m_dt/rho * mu0*redge[0]/(r*top_h[2]*top_h[2]);
	    coeff[5] = 0.5*m_dt/rho * mu0*redge[1]/(r*top_h[2]*top_h[2]);

	    for (nb = 0; nb < 6; nb++) //For the second order boundary
	    {
		if (I_nb[nb] == -1)
		{
		    coeff[nb] = 2.0*coeff[nb];
		}
	    }

	    solver.Set_A(I*3,I*3,1.0+coeff[0]+coeff[1]+coeff[2]+coeff[3]+coeff[4]+coeff[5]+0.5*m_dt*mu0/(rho*rr));
	    rhs = (1.0-coeff[0]-coeff[1]-coeff[2]-coeff[3]-coeff[4]-coeff[5]-0.5*m_dt*mu0/(rho*rr))*cell_center[index].m_state.m_U[0];


	    coeff[6] = -0.5*m_dt/rho * mu0/(rr*top_h[0]);
	    coeff[7] =  0.5*m_dt/rho * mu0/(rr*top_h[0]);


	    for (nb = 0; nb < 6; nb++)
	    {
		if (I_nb[nb] != -1)
		{
		    solver.Set_A(I*3,I_nb[nb]*3,-coeff[nb]);
		    rhs += coeff[nb]*U0_nb[nb];
		}
		else
		    rhs += 2.0*coeff[nb]*U0_nb[nb];
	    }

	    if (I_nb[0] != -1)
	    {
		solver.Set_A(I*3,I_nb[0]*3+2,-coeff[6]);
		rhs += coeff[6]*U2_nb[0];
	    }
	    else
	    {
		solver.Set_A(I*3,I*3+2,coeff[6]);
		rhs += -coeff[6]*U2_center;
		coeff[6] = coeff[6] * 2.0;
		rhs += 2.0*coeff[6]*U2_nb[0];
	    }

	    if (I_nb[1] != -1)
	    {
		solver.Set_A(I*3,I_nb[1]*3+2,-coeff[7]);
		rhs += coeff[7]*U2_nb[1];
	    }
	    else
	    {
		solver.Set_A(I*3,I*3+2,coeff[7]);
		rhs += -coeff[7]*U2_center;
		coeff[7] = coeff[7] * 2.0;
		rhs += 2.0*coeff[7]*U2_nb[1];
	    }

	    rhs -= m_dt*cell_center[index].m_state.m_U[0]*cell_center[index].m_state.m_U[2]/r; //Source term in advection step
	    rhs += m_dt*state.m_U[0];
	    rhs += m_dt*cell_center[index].m_state.f_surf[0];
	    rhs -= m_dt*cell_center[index].m_state.grad_q[0]/rho;

	    solver.Set_b(I*3, rhs);

	    /************************************************************************/


    //Setting the coeffecients for the second equation


	    coeff[0] = 0.5*m_dt/rho * mu0/(rr*top_h[0]*top_h[0]);
	    coeff[1] = 0.5*m_dt/rho * mu0/(rr*top_h[0]*top_h[0]);
	    coeff[2] = 0.5*m_dt/rho * mu0/(top_h[1]*top_h[1]);
	    coeff[3] = 0.5*m_dt/rho * mu0/(top_h[1]*top_h[1]);

	    coeff[4] = 0.5*m_dt/rho * mu0*redge[0]/(r*top_h[2]*top_h[2]);
	    coeff[5] = 0.5*m_dt/rho * mu0*redge[1]/(r*top_h[2]*top_h[2]);

	    for (nb = 0; nb < 6; nb++)
	    {
		if (I_nb[nb] == -1)
		    coeff[nb] = 2.0*coeff[nb];
	    }

	    solver.Set_A(I*3+1,I*3+1,1.0+coeff[0]+coeff[1]+coeff[2]+coeff[3]+coeff[4]+coeff[5]);
	    rhs = (1.0-coeff[0]-coeff[1]-coeff[2]-coeff[3]-coeff[4]-coeff[5])*cell_center[index].m_state.m_U[1];


	    for (nb = 0; nb < 6; nb++)
	    {
		if (I_nb[nb] != -1)
		{
		    solver.Set_A(I*3+1,I_nb[nb]*3+1,-coeff[nb]);
		    rhs += coeff[nb]*U1_nb[nb];
		}
		else
		    rhs += 2.0*coeff[nb]*U1_nb[nb];
	    }

	    rhs += m_dt*state.m_U[1];
	    rhs += m_dt*cell_center[index].m_state.f_surf[1];
	    rhs -= m_dt*cell_center[index].m_state.grad_q[1]/rho;

	    solver.Set_b(I*3+1, rhs);

	    /************************************************************************/

    //Setting the coeffecients for the third equation

	    coeff[0] = 0.5*m_dt/rho * mu0/(rr*top_h[0]*top_h[0]);
	    coeff[1] = 0.5*m_dt/rho * mu0/(rr*top_h[0]*top_h[0]);
	    coeff[2] = 0.5*m_dt/rho * mu0/(top_h[1]*top_h[1]);
	    coeff[3] = 0.5*m_dt/rho * mu0/(top_h[1]*top_h[1]);

	    coeff[4] = 0.5*m_dt/rho * mu0*redge[0]/(r*top_h[2]*top_h[2]);
	    coeff[5] = 0.5*m_dt/rho * mu0*redge[1]/(r*top_h[2]*top_h[2]);

	    for (nb = 0; nb < 6; nb++)
	    {
		if (I_nb[nb] == -1)
		    coeff[nb] = 2.0*coeff[nb];
	    }

	    solver.Set_A(I*3+2,I*3+2,1.0+coeff[0]+coeff[1]+coeff[2]+coeff[3]+coeff[4]+coeff[5]+0.5*m_dt*mu0/(rho*rr));
	    rhs = (1.0-coeff[0]-coeff[1]-coeff[2]-coeff[3]-coeff[4]-coeff[5]-0.5*m_dt*mu0/(rho*rr))*cell_center[index].m_state.m_U[2];


	    coeff[6] =  0.5*m_dt/rho * mu0/(rr*top_h[0]);
	    coeff[7] = -0.5*m_dt/rho * mu0/(rr*top_h[0]);


	    for (nb = 0; nb < 6; nb++)
	    {
		if (I_nb[nb] != -1)
		{
		    solver.Set_A(I*3+2,I_nb[nb]*3+2,-coeff[nb]);
		    rhs += coeff[nb]*U2_nb[nb];
		}
		else
		    rhs += 2.0*coeff[nb]*U2_nb[nb];
	    }


	    if (I_nb[0] != -1)
	    {
		solver.Set_A(I*3+2, I_nb[0]*3, -coeff[6]);
		rhs += coeff[6]*U0_nb[0];
	    }
	    else
	    {
		solver.Set_A(I*3+2, I*3, coeff[6]);
		rhs += -coeff[6]*U0_center;
		coeff[6] = coeff[6] * 2.0;
		rhs += 2.0*coeff[6]*U0_nb[0];
	    }
	    
	    if (I_nb[1] != -1)
	    {
		solver.Set_A(I*3+2, I_nb[1]*3, -coeff[7]);
		rhs += coeff[7]*U0_nb[1];
	    }
	    else
	    {
		solver.Set_A(I*3+2, I*3, coeff[7]);
		rhs += -coeff[7]*U0_center;
		coeff[7] = coeff[7] * 2.0;
		rhs += 2.0*coeff[7]*U0_nb[1];
	    }

	    rhs += m_dt* cell_center[index].m_state.m_U[0] * cell_center[index].m_state.m_U[0] / r; //Source term in advection step
	    rhs += m_dt*state.m_U[2];
	    rhs += m_dt*cell_center[index].m_state.f_surf[2];
	    rhs -= m_dt*cell_center[index].m_state.grad_q[2]/rho;

	    solver.Set_b(I*3+2, rhs);
        }

        solver.SetMaxIter(40000);
        solver.SetTol(1e-12);

	start_clock("Before Petsc Solve");
        solver.Solve_GMRES();
	solver.GetNumIterations(&num_iter);
	solver.GetFinalRelativeResidualNorm(&rel_residual);
	stop_clock("After Petsc Solve");

        // get back the solution
        solver.Get_x(x);

        if (debugging("PETSc"))
            (void) printf("Incompress_Solver_Smooth_3D_Cartesian::compDiffWithSmoothProperty_1st_decoupled: "
                        "num_iter = %d, rel_residual = %g. \n",
                        num_iter,rel_residual);


	for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I = ijk_to_I[i][j][k];
            index = d_index3d(i,j,k,top_gmax);
            if (I >= 0)
            {
                cell_center[index].m_state.m_U[0] = x[I*3-ilower*3];
                cell_center[index].m_state.m_U[1] = x[I*3-ilower*3+1];
		cell_center[index].m_state.m_U[2] = x[I*3-ilower*3+2];
                speed = fabs(cell_center[index].m_state.m_U[0]) +
                        fabs(cell_center[index].m_state.m_U[1]) +
			fabs(cell_center[index].m_state.m_U[2]);
                if (speed > max_speed)
                    max_speed = speed;
            }
            else
            {
                cell_center[index].m_state.m_U[0] = 0.0;
                cell_center[index].m_state.m_U[1] = 0.0;
		cell_center[index].m_state.m_U[2] = 0.0;
            }
        }
        for (l = 0; l < 3; ++l)
        {
	    for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)

            {
                index  = d_index3d(i,j,k,top_gmax);
                array[index] = cell_center[index].m_state.m_U[l];
            }
            scatMeshArray();
	    for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)

            {
                index  = d_index3d(i,j,k,top_gmax);
                cell_center[index].m_state.m_U[l] = array[index];
            }
        }
        pp_global_max(&max_speed,1);

        free_these(1,x);

}       /* end compDiffWithSmoothProperty3d */

void Incompress_Solver_Smooth_3D_Cylindrical::computeProjectionCim(void)
{
	(void) printf("computeProjectionCim() not implemented\n");
	clean_up(ERROR);
}	/* end computeProjectionCim */

void Incompress_Solver_Smooth_3D_Cylindrical::computeProjection(void)
{
        switch (iFparams->num_scheme.ellip_method)
        {
        case SIMPLE_ELLIP:
            computeProjectionSimple();
            return;
        case CIM_ELLIP:
            computeProjectionCim();
            return;
        }
}       /* end computeProjection */

void Incompress_Solver_Smooth_3D_Cylindrical::computeProjectionSimple(void)
{
	int index, index_nb[6], size;
	double rhs, coeff[6], rho[6], rho0;
	int I,I_nb[6];
	int i,j,k,l,icoords[MAXD];
	double P_max,P_min;
	COMPONENT comp;
	double aII;
	double **vel = iFparams->field->vel;
	
	max_value = 0.0;
	double value;
	double sum_div;
	sum_div = 0.0;
	int num_iter = 0;
	double rel_residual = 0.0;

	PETSc solver;
	solver.Create(ilower, iupper-1, 7, 7);
	solver.Reset_A();
	solver.Reset_b();
	solver.Reset_x();
	size = iupper - ilower;

	setIndexMap();

	for (l = 0; l < dim; ++l)
	for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index  = d_index3d(i,j,k,top_gmax);
	    vel[l][index] = cell_center[index].m_state.m_U[l];
	}

	/* Compute velocity divergence */
	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;
	    index  = d_index3d(i,j,k,top_gmax);
	    array[index] = computeFieldPointDiv(icoords,vel);
	}
	scatMeshArray();
	for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index  = d_index3d(i,j,k,top_gmax);
	    cell_center[index].m_state.div_U = array[index];    
	}


	if(debugging("step_size"))
	{
	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		index = d_index3d(i,j,k,top_gmax);
	        value = fabs(cell_center[index].m_state.div_U);
		sum_div = sum_div + cell_center[index].m_state.div_U * cell_center[index].m_coords[2] * top_h[0]*top_h[1]*top_h[2];
	        if(value > max_value)
		    max_value = value;
	    }
	    pp_global_sum(&sum_div,1);
	    printf("\nThe summation of divergence of U is %.16g\n",sum_div);
	    pp_global_max(&max_value,1);
	    printf("\nThe max value of divergence of U is %.16g\n",max_value);
	    max_value = 0.0;
	}
	

	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    index  = d_index3d(i,j,k,top_gmax);
	    comp = top_comp[index];
	    I = ijk_to_I[i][j][k];
	    if (I == -1) continue;

	    index_nb[0] = d_index3d(i-1,j,k,top_gmax);
	    index_nb[1] = d_index3d(i+1,j,k,top_gmax);
	    index_nb[2] = d_index3d(i,j-1,k,top_gmax);
	    index_nb[3] = d_index3d(i,j+1,k,top_gmax);
	    index_nb[4] = d_index3d(i,j,k-1,top_gmax);
	    index_nb[5] = d_index3d(i,j,k+1,top_gmax);
	    I_nb[0] = ijk_to_I[i-1][j][k];
	    I_nb[1] = ijk_to_I[i+1][j][k];
	    I_nb[2] = ijk_to_I[i][j-1][k];
	    I_nb[3] = ijk_to_I[i][j+1][k];
	    I_nb[4] = ijk_to_I[i][j][k-1];
	    I_nb[5] = ijk_to_I[i][j][k+1];
	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;
	
	    rho0   = cell_center[index].m_state.m_rho;
	    for (l = 0; l < 6; ++l)
	    {
		if (I_nb[l] == -1)
		    index_nb[l] = index;
		rho[l] = 1.0/2*(rho0 + cell_center[index_nb[l]].m_state.m_rho);
		//coeff[l] = 1.0/rho[l]/sqr(top_h[l/2]);
	    }

	    double rr = sqr(cell_center[index].m_coords[2]);
	    double r = cell_center[index].m_coords[2];
	    double redge[2];
	    redge[0] = 0.5 * (cell_center[index].m_coords[2] + cell_center[index_nb[4]].m_coords[2]);
	    redge[1] = 0.5 * (cell_center[index].m_coords[2] + cell_center[index_nb[5]].m_coords[2]);


	    coeff[0] = 1.0/rho[0] * 1.0/(rr*top_h[0]*top_h[0]);
            coeff[1] = 1.0/rho[1] * 1.0/(rr*top_h[0]*top_h[0]);
            coeff[2] = 1.0/rho[2] * 1.0/(top_h[1]*top_h[1]);
            coeff[3] = 1.0/rho[3] * 1.0/(top_h[1]*top_h[1]);
            //coeff[4] = 1.0/rho[4] * ((1.0/(top_h[2]*top_h[2])) - (1.0/(2.0*r*top_h[2])));
            //coeff[5] = 1.0/rho[5] * ((1.0/(top_h[2]*top_h[2])) + (1.0/(2.0*r*top_h[2])));

	    coeff[4] = redge[0]/rho[4] / (r*top_h[2]*top_h[2]);
	    coeff[5] = redge[1]/rho[5] / (r*top_h[2]*top_h[2]);

            rhs = cell_center[index].m_state.div_U/accum_dt;

	    aII = 0.0;
	    for (l = 0; l < 6; ++l)
	    {
	    	if (I_nb[l] != -1)
		{
		    solver.Set_A(I,I_nb[l],coeff[l]);
                    aII += -coeff[l];
		}
	    }
            if (aII != 0.0)
	    {
                solver.Set_A(I,I,aII);
	    }
            else
            {
		printf("\nUsing the original pressure!\n");
                solver.Set_A(I,I,1.0);
                rhs = cell_center[index].m_state.m_P;
            }
            solver.Set_b(I,rhs);
	}
	
	solver.SetMaxIter(40000);
	solver.SetTol(1e-14);

	start_clock("Before Petsc Solver in Projection step");
	solver.Solve_withPureNeumann();
	solver.GetNumIterations(&num_iter);
	solver.GetFinalRelativeResidualNorm(&rel_residual);
	if(rel_residual > 1)
	{
	    printf("\n The solution diverges! The residual is %g. Solve again using GMRES!\n",rel_residual);
	    solver.Reset_x();
	    solver.Solve_withPureNeumann_GMRES();
	    solver.GetNumIterations(&num_iter);
	    solver.GetFinalRelativeResidualNorm(&rel_residual);
	}
	stop_clock("After Petsc Solver in Projection step");

	double *x;
	FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));
	solver.Get_x(x);

	if (debugging("PETSc"))
	    (void) printf("Incompress_Solver_Smooth_3D_Cylindrical::"
			"computeProjection: "
	       		"num_iter = %d, rel_residual = %g \n", 
			num_iter, rel_residual);
	
	P_max = -HUGE;		P_min = HUGE;
	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    I = ijk_to_I[i][j][k];
	    array[index] = x[I-ilower];
	}
	scatMeshArray();
	for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index  = d_index3d(i,j,k,top_gmax);
	    cell_center[index].m_state.m_phi = array[index];
	}

	if(debugging("step_size"))
	{
	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		index = d_index3d(i,j,k,top_gmax);
		value = fabs(cell_center[index].m_state.m_phi);
		if (value > max_value)
		    max_value = value;
	    }
	    pp_global_max(&max_value,1);
	    printf("\nThe max value of phi is %.16g\n",max_value);
	}

	free_these(1,x);
}	/* end computeProjection3dSimple */

void Incompress_Solver_Smooth_3D_Cylindrical::computeNewVelocity(void)
{
	int i, j, k, l, index;
	double grad_phi[3], rho;
	COMPONENT comp;
	double speed;
	int icoords[MAXD];
	double r;

	max_speed = 0.0;

	for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    array[index] = cell_center[index].m_state.m_phi;
	}
	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    comp = top_comp[index];
	    if (!ifluid_comp(comp))
	    {
		for (l = 0; l < 3; ++l)
		    cell_center[index].m_state.m_U[l] = 0.0;
		continue;
	    }
	    rho = cell_center[index].m_state.m_rho;
	    r = cell_center[index].m_coords[2];
	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;

	    computeFieldPointGrad(icoords,array,grad_phi);
	    speed = 0.0;
	    for (l = 0; l < 3; ++l)
	    {
	    	cell_center[index].m_state.m_U[l] -= accum_dt/rho*grad_phi[l];
		speed += fabs(cell_center[index].m_state.m_U[l]);
	    }

	    if (speed > max_speed)
		max_speed = speed;
	}
	for (l = 0; l < 3; ++l)
	{
	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
	    {	
	    	index  = d_index3d(i,j,k,top_gmax);
	    	array[index] = cell_center[index].m_state.m_U[l];
	    }
	    scatMeshArray();
	    for (k = 0; k <= top_gmax[2]; k++)
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {	
	    	index  = d_index3d(i,j,k,top_gmax);
	    	cell_center[index].m_state.m_U[l] = array[index];
	    }
	}
	pp_global_max(&max_speed,1);
}	/* end computeNewVelocity3d */


void Incompress_Solver_Smooth_3D_Cylindrical::computeSourceTerm(double *coords, L_STATE &state) 
{
	int i;
	for (i = 0; i < dim; ++i)
	    state.m_U[i] = iFparams->gravity[i];

	state.m_P = HUGE_VAL;
}
void Incompress_Solver_Smooth_3D_Cylindrical::computeSourceTerm(double *coords, double t, L_STATE &state) 
{
	computeSourceTerm(coords, state);
}

// for initial condition: 
// 		setInitialCondition();	
// this function should be called before solve()
// for the source term of the momentum equation: 	
// 		computeSourceTerm();
void Incompress_Solver_Smooth_3D_Cylindrical::solve(double dt)
{
        printf("\nEntering solve, the dt for this solve is : %.16g\n",dt);

	m_t_old = front->time;
	m_t_int = front->time + dt/2.0;
	m_t_new = front->time + dt;

	static boolean first = YES;
	if (first)
	{
	    accum_dt = 0.0;
	    first = NO;
	}
	m_dt = dt;
	max_speed = 0.0;

	start_clock("solve");
	setDomain();

	setComponent();
	if (debugging("trace"))
	    printf("Passed setComponent()\n");
	setGlobalIndex();
	if (debugging("trace"))
	    printf("Passed setGlobalIndex()\n");
	start_clock("setSmoothedProperties");
	setSmoothedProperties();
	stop_clock("setSmoothedProperties");
	if (debugging("trace"))
	    printf("Passed setSmoothedProperties()\n");
	
	// 1) solve for intermediate velocity
	start_clock("computeAdvection");
	//computeAdvection_test(); //Discretize the equation using the cylindrical Paper, result seems to be similar (First Order)
	//computeAdvection();  //First order discretization
	//compAdvectionTerm_decoupled(); //Second order convection term using slope limiter for first derivative (BELL 1989)
	compAdvectionTerm_coupled();

	stop_clock("computeAdvection");
	if (debugging("trace"))
	    printf("max_speed after computeAdvection(): %20.14f\n",
				max_speed);

	start_clock("compDiffWithSmoothProperty");
	//compDiffWithSmoothProperty_1st_decoupled_source();
	//compDiffWithSmoothProperty_1st_decoupled_test(); //Discretize the equatino using the cylindrical Paper, result seems to be similar
	//compDiffWithSmoothProperty_1st_decoupled();
	//compDiffWithSmoothProperty_2nd_decoupled(); //2nd order diffusion solver with the advection source terms
	compDiffWithSmoothProperty_2nd_coupled();
	//compDiffWithSmoothProperty_2nd_decoupled_Shuqiang(); //2nd order diffusion solver by Shuqiang
	stop_clock("compDiffWithSmoothProperty");



        start_clock("compSGS");
        //compSGS();	//Subgrid model by Hyunkyun Lim
        stop_clock("compSGS");

	if (debugging("trace"))
	    printf("max_speed after compDiffWithSmoothProperty(): %20.14f\n",
				max_speed);

	// 2) projection step
	accum_dt += m_dt;
	if (accum_dt >= min_dt)
	{
	    start_clock("computeProjection");
	    computeProjection();
	    //computeProjection_Shuqiang(); //Projection step by Shuqiang
	    stop_clock("computeProjection");

	    start_clock("computePressure");
	    computePressure();
	    stop_clock("computePressure");

	    start_clock("computeNewVelocity");
	    computeNewVelocity();
	    stop_clock("computeNewVelocity");
	    accum_dt = 0.0;
	}

	if (debugging("sample_velocity"))
	{
	    sampleVelocity();
	}

	if (debugging("trace"))
	    printf("max_speed after computeNewVelocity(): %20.14f\n",
				max_speed);

	start_clock("copyMeshStates");
	copyMeshStates();
	stop_clock("copyMeshStates");

	setAdvectionDt();
	stop_clock("solve");
}	/* end solve */


double Incompress_Solver_Smooth_3D_Cylindrical::getVorticityX(int i, int j, int k)
{
	int index0,index00,index01,index10,index11;
	double v00,v01,v10,v11;
	double dy,dz;
	double vorticity;

	dy = top_h[1];
	dz = top_h[2];
	index0 = d_index3d(i,j,k,top_gmax);
	index00 = d_index3d(i,j-1,k,top_gmax);
	index01 = d_index3d(i,j+1,k,top_gmax);
	index10 = d_index3d(i,j,k-1,top_gmax);
	index11 = d_index3d(i,j,k+1,top_gmax);
	v00 = -cell_center[index00].m_state.m_U[2];
	v01 =  cell_center[index01].m_state.m_U[2];
	v10 =  cell_center[index10].m_state.m_U[1];
	v11 = -cell_center[index11].m_state.m_U[1];

	vorticity = (v00 + v01)/2.0/dz + (v10 + v11)/2.0/dy;
	return vorticity;
}	/* end getVorticityX */

double Incompress_Solver_Smooth_3D_Cylindrical::getVorticityY(int i, int j, int k)
{
	int index0,index00,index01,index10,index11;
	double v00,v01,v10,v11;
	double dx,dz;
	double vorticity;

	dx = top_h[0];
	dz = top_h[2];
	index0 = d_index3d(i,j,k,top_gmax);
	index00 = d_index3d(i,j,k-1,top_gmax);
	index01 = d_index3d(i,j,k+1,top_gmax);
	index10 = d_index3d(i-1,j,k,top_gmax);
	index11 = d_index3d(i+1,j,k,top_gmax);
	v00 = -cell_center[index00].m_state.m_U[0];
	v01 =  cell_center[index01].m_state.m_U[0];
	v10 =  cell_center[index10].m_state.m_U[2];
	v11 = -cell_center[index11].m_state.m_U[2];

	vorticity = (v00 + v01)/2.0/dx + (v10 + v11)/2.0/dz;
	return vorticity;
}	/* end getVorticityY */

double Incompress_Solver_Smooth_3D_Cylindrical::getVorticityZ(int i, int j, int k)
{
	int index0,index00,index01,index10,index11;
	double v00,v01,v10,v11;
	double dx,dy;
	double vorticity;

	dx = top_h[0];
	dy = top_h[1];
	index0 = d_index3d(i,j,k,top_gmax);
	index00 = d_index3d(i-1,j,k,top_gmax);
	index01 = d_index3d(i+1,j,k,top_gmax);
	index10 = d_index3d(i,j-1,k,top_gmax);
	index11 = d_index3d(i,j+1,k,top_gmax);
	v00 = -cell_center[index00].m_state.m_U[1];
	v01 =  cell_center[index01].m_state.m_U[1];
	v10 =  cell_center[index10].m_state.m_U[0];
	v11 = -cell_center[index11].m_state.m_U[0];

	vorticity = (v00 + v01)/2.0/dy + (v10 + v11)/2.0/dx;
	return vorticity;
}	/* end getVorticityZ */


void Incompress_Solver_Smooth_3D_Cylindrical::copyMeshStates()
{
	int i,j,k,d,index;
	double **vel = field->vel;
	double *pres = field->pres;
	double **vort3d = field->vort3d;

	for (i = imin; i <= imax; ++i)
	for (j = jmin; j <= jmax; ++j)
	for (k = kmin; k <= kmax; ++k)
	{
	    index  = d_index3d(i,j,k,top_gmax);
	    if (ifluid_comp(top_comp[index]))
	    {
		pres[index] = cell_center[index].m_state.m_P;
	    	vel[0][index] = cell_center[index].m_state.m_U[0];
	    	vel[1][index] = cell_center[index].m_state.m_U[1];
	    	vel[2][index] = cell_center[index].m_state.m_U[2];
		vort3d[0][index] = getVorticityX(i,j,k);
		vort3d[1][index] = getVorticityY(i,j,k);
		vort3d[2][index] = getVorticityZ(i,j,k);
	    }
	    else
	    {
	    	pres[index] = 0.0;
		for (d = 0; d < 3; ++d)
		{
		    vel[d][index] = 0.0;
		    vort3d[d][index] = 0.0;
		}
	    }
	}
	FT_ParallelExchGridArrayBuffer(pres,front);
	FT_ParallelExchGridArrayBuffer(vel[0],front);
	FT_ParallelExchGridArrayBuffer(vel[1],front);
	FT_ParallelExchGridArrayBuffer(vel[2],front);
	FT_ParallelExchGridArrayBuffer(vort3d[0],front);
	FT_ParallelExchGridArrayBuffer(vort3d[1],front);
	FT_ParallelExchGridArrayBuffer(vort3d[2],front);
}	/* end copyMeshStates */


void Incompress_Solver_Smooth_3D_Cylindrical::
	compDiffWithSmoothProperty_1st_decoupled_source(void)
{
        COMPONENT comp;
        int index,index_nb[6],size;
        int I,I_nb[6];
	int i,j,k,l,nb,icoords[MAXD];
        L_STATE state;
	double coords[MAXD], crx_coords[MAXD];
	double coeff[6],mu[6],mu0,rho,rhs,U_nb[6];
        double speed;
        double *x;
	GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
	POINTER intfc_state;
	HYPER_SURF *hs;
	int num_iter;
	double rel_residual;

        setIndexMap();

	max_speed = 0.0;

        size = iupper - ilower;
        FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));

	for (l = 0; l < dim; ++l)
	{
            PETSc solver;
            solver.Create(ilower, iupper-1, 7, 7);
	    solver.Reset_A();
	    solver.Reset_b();
	    solver.Reset_x();

            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
            	I  = ijk_to_I[i][j][k];
            	if (I == -1) continue;

            	index  = d_index3d(i,j,k,top_gmax);
            	index_nb[0] = d_index3d(i-1,j,k,top_gmax);
            	index_nb[1] = d_index3d(i+1,j,k,top_gmax);
            	index_nb[2] = d_index3d(i,j-1,k,top_gmax);
            	index_nb[3] = d_index3d(i,j+1,k,top_gmax);
            	index_nb[4] = d_index3d(i,j,k-1,top_gmax);
            	index_nb[5] = d_index3d(i,j,k+1,top_gmax);

		icoords[0] = i;
		icoords[1] = j;
		icoords[2] = k;
		comp = top_comp[index];

            	I_nb[0] = ijk_to_I[i-1][j][k]; //west
            	I_nb[1] = ijk_to_I[i+1][j][k]; //east
            	I_nb[2] = ijk_to_I[i][j-1][k]; //south
            	I_nb[3] = ijk_to_I[i][j+1][k]; //north
            	I_nb[4] = ijk_to_I[i][j][k-1]; //lower
            	I_nb[5] = ijk_to_I[i][j][k+1]; //upper


            	mu0   = cell_center[index].m_state.m_mu;
            	rho   = cell_center[index].m_state.m_rho;

            	for (nb = 0; nb < 6; nb++)
            	{
                    if (FT_StateStructAtGridCrossing(front,icoords,dir[nb],
                                comp,&intfc_state,&hs,crx_coords) &&
                                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
		    {
			U_nb[nb] = getStateVel[l](intfc_state);
			if (wave_type(hs) == DIRICHLET_BOUNDARY || 
			    wave_type(hs) == NEUMANN_BOUNDARY)
			    mu[nb] = mu0;
			else
			    mu[nb] = 1.0/2*(mu0 + 
				cell_center[index_nb[nb]].m_state.m_mu);
		    }
                    else
		    {
                    	U_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[l];
			mu[nb] = 1.0/2*(mu0 + 
				cell_center[index_nb[nb]].m_state.m_mu);
		    }
            	}


                double rr = sqr(cell_center[index].m_coords[2]);
                double r = cell_center[index].m_coords[2];
                coeff[0] = 0.5*m_dt/rho * mu[0]/(rr*top_h[0]*top_h[0]);
                coeff[1] = 0.5*m_dt/rho * mu[1]/(rr*top_h[0]*top_h[0]);
                coeff[2] = 0.5*m_dt/rho * mu[2]/(top_h[1]*top_h[1]);
                coeff[3] = 0.5*m_dt/rho * mu[3]/(top_h[1]*top_h[1]);
                coeff[4] = 0.5*m_dt/rho * ((mu[4]/(top_h[2]*top_h[2]))-
					(mu[4]/(2.0*top_h[2]*r)));
                coeff[5] = 0.5*m_dt/rho * ((mu[5]/(top_h[2]*top_h[2]))+
					(mu[5]/(2.0*top_h[2]*r)));


            	getRectangleCenter(index, coords);
            	computeSourceTerm(coords, state);

                //SOURCE TERM
                double sour[3];
                double temp0, temp1, temp2, temp3, temp4, temp5, temp6, temp7,
                       temp8, temp9, temp10, temp11, temp12, temp13, temp14,
                       temp15, temp16, temp17, temp18, temp19, temp20;

                    temp0 = cell_center[index].m_state.m_U[0];
                    temp1 = cell_center[index].m_state.m_U[1];
                    temp2 = cell_center[index].m_state.m_U[2];
                    temp3 = cell_center[index_nb[0]].m_state.m_U[0];
                    temp4 = cell_center[index_nb[0]].m_state.m_U[1];
                    temp5 = cell_center[index_nb[0]].m_state.m_U[2];
                    temp6 = cell_center[index_nb[1]].m_state.m_U[0];
                    temp7 = cell_center[index_nb[1]].m_state.m_U[1];
                    temp8 = cell_center[index_nb[1]].m_state.m_U[2];
                    temp9 = cell_center[index_nb[2]].m_state.m_U[0];
                    temp10 = cell_center[index_nb[2]].m_state.m_U[1];
                    temp11 = cell_center[index_nb[2]].m_state.m_U[2];
                    temp12 = cell_center[index_nb[3]].m_state.m_U[0];
                    temp13 = cell_center[index_nb[3]].m_state.m_U[1];
                    temp14 = cell_center[index_nb[3]].m_state.m_U[2];
                    temp15 = cell_center[index_nb[4]].m_state.m_U[0];
                    temp16 = cell_center[index_nb[4]].m_state.m_U[1];
                    temp17 = cell_center[index_nb[4]].m_state.m_U[2];
                    temp18 = cell_center[index_nb[5]].m_state.m_U[0];
                    temp19 = cell_center[index_nb[5]].m_state.m_U[1];
                    temp20 = cell_center[index_nb[5]].m_state.m_U[2];

                    if(l == 0)
                    {
                        sour[0] = (m_dt*cell_center[index].m_state.m_mu/cell_center[index].m_state.m_rho
                                                  * (((2.0/sqr(cell_center[index].m_coords[2]))
                                                  * (temp8
                                                  - temp5) / top_h[0])
                                                  - (temp0 / sqr(cell_center[index].m_coords[2]))));
                    }
                    if(l == 1)
                    {
                        sour[1] = 0.0;
                    }
                    if(l == 2)
                    {
                        sour[2] = - (m_dt*cell_center[index].m_state.m_mu/cell_center[index].m_state.m_rho
                                                  * (((2.0/sqr(cell_center[index].m_coords[2]))
                                                  *(temp6
                                                  - temp3) / top_h[0])
                                                  + (temp2 / sqr(cell_center[index].m_coords[2]))));
                    }

        	//first equation
            	solver.Set_A(I,I,1+coeff[0]+coeff[1]+coeff[2]+coeff[3]+
				coeff[4]+coeff[5]);
		rhs = (1-coeff[0]-coeff[1]-coeff[2]-coeff[3]-coeff[4]-coeff[5])*
		      		cell_center[index].m_state.m_U[l];

		for(nb = 0; nb < 6; nb++)
		{
		    if(I_nb[nb] != -1)
		    {
			solver.Set_A(I,I_nb[nb],-coeff[nb]);
			rhs += coeff[nb]*U_nb[nb];
		    }
		    else
			rhs += 2.0*coeff[nb]*U_nb[nb];
		}
		rhs += m_dt*state.m_U[l];
                rhs += sour[l];
		rhs += m_dt*cell_center[index].m_state.f_surf[l];
		rhs -= m_dt*cell_center[index].m_state.grad_q[l]/rho;
 
		solver.Set_b(I, rhs);

            }

            solver.SetMaxIter(40000);
            solver.SetTol(1e-12);

	    start_clock("Befor Petsc solve");
            solver.Solve_GMRES();
            solver.GetNumIterations(&num_iter);
            solver.GetFinalRelativeResidualNorm(&rel_residual);
            // get back the solution
            solver.Get_x(x);

            if (debugging("PETSc"))
                (void) printf("L_CARTESIAN::"
			"compDiffWithSmoothProperty_1st_decoupled_source: "
                        "num_iter = %d, rel_residual = %g. \n",
                        num_iter,rel_residual);

	    for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                I = ijk_to_I[i][j][k];
                index = d_index3d(i,j,k,top_gmax);
                if (I >= 0)
                {
                    cell_center[index].m_state.m_U[l] = x[I-ilower];
                }
                else
                {
                    cell_center[index].m_state.m_U[l] = 0.0;
                }
            }

	    for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index  = d_index3d(i,j,k,top_gmax);
                array[index] = cell_center[index].m_state.m_U[l];
            }
            scatMeshArray();
	    for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index3d(i,j,k,top_gmax);
                cell_center[index].m_state.m_U[l] = array[index];
            }
        }
	for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
            speed = fabs(cell_center[index].m_state.m_U[0]) +
                    fabs(cell_center[index].m_state.m_U[1]) +
                    fabs(cell_center[index].m_state.m_U[2]);
            if (speed > max_speed)
                    max_speed = speed;
	}
        pp_global_max(&max_speed,1);

        free_these(1,x);
}       /* end compDiffWithSmoothProperty3d_decoupled */


void Incompress_Solver_Smooth_3D_Cylindrical::computePressurePmI(void)
{
        int i,j,k,index;

	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
	for (i = imin; i <= imax; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    cell_center[index].m_state.m_P += cell_center[index].m_state.m_phi;
	    array[index] = cell_center[index].m_state.m_P;
	}
	scatMeshArray();

	for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    cell_center[index].m_state.m_P = array[index];
	    cell_center[index].m_state.m_q = array[index];
	}



}        /* end computePressurePmI3d */

void Incompress_Solver_Smooth_3D_Cylindrical::computePressurePmII(void)
{
        int i,j,k,index;
        double mu0;

        for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
	{
            index = d_index3d(i,j,k,top_gmax);
            mu0 = 0.5*cell_center[index].m_state.m_mu;
            cell_center[index].m_state.m_P += 
				cell_center[index].m_state.m_phi -
                        	accum_dt*mu0*cell_center[index].m_state.div_U;
	    cell_center[index].m_state.m_q = cell_center[index].m_state.m_P;
	}
}        /* end computePressurePmII3d */

void Incompress_Solver_Smooth_3D_Cylindrical::computePressurePmIII(void)
{
        int i,j,k,index;
        double mu0;

        for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
	{
            index = d_index3d(i,j,k,top_gmax);
            mu0 = 0.5*cell_center[index].m_state.m_mu;
            cell_center[index].m_state.m_P = 
				cell_center[index].m_state.m_phi -
                        	accum_dt*mu0*cell_center[index].m_state.div_U;
	    cell_center[index].m_state.m_q = 0.0;
	}
}        /* end computePressurePmIII3d */

void Incompress_Solver_Smooth_3D_Cylindrical::computePressure(void)
{
	switch (iFparams->num_scheme.projc_method)
	{
	case BELL_COLELLA:
	    computePressurePmI();
	    break;
	case KIM_MOIN:
	    computePressurePmII();
	    break;
	case SIMPLE:
	case PEROT_BOTELLA:
	    computePressurePmIII();
	    break;
	case ERROR_PROJC_SCHEME:
	default:
	    (void) printf("Unknown competePressure scheme\n");
	    clean_up(ERROR);
	}
	computeGradientQ();
}	/* end computePressure */

void Incompress_Solver_Smooth_3D_Cylindrical::computeGradientQ(void)
{
	int i,j,k,l,index;
	double *grad_q;
	int icoords[MAXD];

	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
	for (i = imin; i <= imax; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    array[index] = cell_center[index].m_state.m_q;

	}
	scatMeshArray();


	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
	for (i = imin; i <= imax; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    grad_q = cell_center[index].m_state.grad_q;
	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;

	    computeFieldPointGrad(icoords,array,grad_q); //grad_q is the Gradient Operator in Cylindrical Coordinate
	}
	for (l = 0; l < dim; ++l)
	{
	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
	    	index = d_index3d(i,j,k,top_gmax);
		array[index] = cell_center[index].m_state.grad_q[l];
	    }
	    scatMeshArray();
	    for (k = 0; k <= top_gmax[2]; k++)
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {
	    	index = d_index3d(i,j,k,top_gmax);
		cell_center[index].m_state.grad_q[l] = array[index];
	    }
	}
}	/* end computeGradientQ3d */

#define		MAX_TRI_FOR_INTEGRAL		100
void Incompress_Solver_Smooth_3D_Cylindrical::surfaceTension(
	double *coords,
	HYPER_SURF_ELEMENT *hse,
	HYPER_SURF *hs,
	double *force,
	double sigma)
{
	int i,j,k,num_tris;
	TRI *tri,*tri_list[MAX_TRI_FOR_INTEGRAL];
	double kappa_tmp,kappa,mag_nor,area,delta;
	double median[MAXD],nor[MAXD];
	POINT *p;

	TriAndFirstRing(hse,hs,&num_tris,tri_list);
	for (i = 0; i < num_tris; ++i)
	{
	    kappa = 0.0;
	    tri = tri_list[i];
	    for (j = 0; j < 3; ++j) median[j] = 0.0;
	    for (j = 0; j < 3; ++j)
	    {
		p = Point_of_tri(tri)[j];
		for (k = 0; k < 3; ++k) 
		    median[k] += Coords(p)[k];
	    	GetFrontCurvature(p,Hyper_surf_element(tri),hs,
				&kappa_tmp,front);
		kappa += kappa_tmp;
		nor[j] = Tri_normal(tri)[j];
	    }
	    kappa /= 3.0;
	    mag_nor = mag_vector(nor,3);
	    area = 0.5*mag_nor;
	    for (j = 0; j < 3; ++j)  
	    {
		nor[j] /= mag_nor;
		median[j] /= 3.0;
	    }
	    delta = smoothedDeltaFunction(coords,median);
	    if (delta == 0.0) continue;
	    for (j = 0; j < dim; ++j) 
	    {
		force[j] += delta*sigma*area*kappa*nor[j];
	    }
	}
}	/* end surfaceTension3d */

void Incompress_Solver_Smooth_3D_Cylindrical::setInitialCondition()
{
	int i,j,k,index,l;
	COMPONENT comp;
	double coords[MAXD];
        double rho_1 = iFparams->rho1;
        double g = iFparams->gravity[1];
	double Omega1 = 0.05;
        double R_1 = 2.538;
        double R_2 = 3.166;
        double A1,B1,term_3;
        double R_1s = sqr(R_1);
        double R_2s = sqr(R_2);
	int size = (int)cell_center.size();

	FT_MakeGridIntfc(front);
	setDomain();
	setComponent();

        m_rho[0] = iFparams->rho1;
        m_rho[1] = iFparams->rho2;
        m_mu[0] = iFparams->mu1;
        m_mu[1] = iFparams->mu2;
	m_comp[0] = iFparams->m_comp1;
	m_comp[1] = iFparams->m_comp2;
	m_smoothing_radius = iFparams->smoothing_radius;
	m_sigma = iFparams->surf_tension;
	mu_min = rho_min = HUGE;
	for (i = 0; i < 2; ++i)
	{
	    if (ifluid_comp(m_comp[i]))
	    {
        	mu_min = std::min(mu_min,m_mu[i]);
        	rho_min = std::min(rho_min,m_rho[i]);
	    }
	}




	// Initialize state at cell_center
        for (i = 0; i < size; i++)
        {
            getRectangleCenter(i, coords);
	    cell_center[i].m_state.setZero();
	    comp = top_comp[i];
	    if (getInitialState != NULL)
	    	(*getInitialState)(comp,coords,cell_center[i].m_state,dim,
						iFparams);

       /* for one phase initialization*/
	    
	    A1 = (Omega1*R_1s)/(R_1s-R_2s);
            B1 = -A1*R_2s;
            term_3 = (sqr(A1)*sqr(cell_center[i].m_coords[2])/2.0)
                          + (2.0*A1*B1*log(cell_center[i].m_coords[2]))
                          - (sqr(B1)/(2.0*sqr(cell_center[i].m_coords[2])))
                          - (sqr(A1)*R_1s/2.0) - (2.0*A1*B1*log(R_1)) + (sqr(B1)/(2.0*R_1s));

            cell_center[i].m_state.m_U[0] = (A1*cell_center[i].m_coords[2])
                                            + (B1/cell_center[i].m_coords[2]);
            cell_center[i].m_state.m_U[1] = 0.0;
            cell_center[i].m_state.m_U[2] = 0.0;        
	    cell_center[i].m_state.m_exactU[0] = cell_center[i].m_state.m_U[0];
	    cell_center[i].m_state.m_exactU[1] = cell_center[i].m_state.m_U[1];
	    cell_center[i].m_state.m_exactU[2] = cell_center[i].m_state.m_U[2];
            cell_center[i].m_state.m_P = (rho_1*term_3) - (rho_1*g*(cell_center[i].m_coords[1])) + 1.0;

        }

	for (l = 0; l < dim; l++) // Scatter the velocity after the initialization
	{
	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		index = d_index3d(i,j,k,top_gmax);
		array[index] = cell_center[index].m_state.m_U[l];
	    }
	    scatMeshArray();
	    for (k = 0; k <= top_gmax[2]; k++)
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {
		index = d_index3d(i,j,k,top_gmax);
		cell_center[index].m_state.m_U[l] = array[index];
	    }
	}

	double totalVolume, cellVolume;
	double averagep, averagep_test;
	double r;
	totalVolume = 0.0;
	averagep = 0.0;
	averagep_test = 0.0;
	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
	for (i = imin; i <= imax; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    r = cell_center[index].m_coords[2];
	    cellVolume = r*top_h[0]*top_h[1]*top_h[2];
	    totalVolume += cellVolume;
	    averagep += cell_center[index].m_state.m_P * cellVolume;
	}

	printf("\ntotalVolume for the processor itself is %.16g\n", totalVolume);

	pp_global_sum(&totalVolume, 1);
	pp_global_sum(&averagep, 1);
	averagep = averagep / totalVolume;
	printf("\nThe averagep is %.16g\n",averagep);
	printf("\nThe totalVolume is %.16g\n",totalVolume);



	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
	for (i = imin; i <= imax; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    cell_center[index].m_state.m_P = cell_center[index].m_state.m_P - averagep;
	    cell_center[index].m_state.m_q = cell_center[index].m_state.m_P;
	    cell_center[index].m_state.m_exactP = cell_center[index].m_state.m_P;
	    array[index] = cell_center[index].m_state.m_P; //For scatter

	    r = cell_center[index].m_coords[2];
	    cellVolume = r*top_h[0]*top_h[1]*top_h[2];
	    averagep_test += cell_center[index].m_state.m_P * cellVolume;
	}

	pp_global_sum(&averagep_test,1);
	averagep_test = averagep_test / totalVolume;
	printf("\nThe average of p after adjusting is %.16g\n", averagep_test);
	printf("\nTotal volume is %.16g\n", totalVolume);

	scatMeshArray();

	for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    cell_center[index].m_state.m_P = array[index];
	    cell_center[index].m_state.m_q = array[index];
	}

	computeGradientQ();
        copyMeshStates();
	setAdvectionDt();
}       /* end setInitialCondition */

double Incompress_Solver_Smooth_3D_Cylindrical::computeFieldPointDiv(
        int *icoords,
        double **field)
{
        int index;
        COMPONENT comp;
        int i, j,k;
	int index_nb[6];
        double div,rur_nb[2],utheta_nb[2],uz_nb[2],ur;
        double crx_coords[MAXD];
        POINTER intfc_state;
        HYPER_SURF *hs;

	i = icoords[0];
	j = icoords[1];
	k = icoords[2];
	index = d_index3d(i,j,k,top_gmax);
        comp = top_comp[index];

	index_nb[0] = d_index3d(i-1,j,k,top_gmax);
	index_nb[1] = d_index3d(i+1,j,k,top_gmax);
	index_nb[2] = d_index3d(i,j-1,k,top_gmax);
	index_nb[3] = d_index3d(i,j+1,k,top_gmax);
	index_nb[4] = d_index3d(i,j,k-1,top_gmax);
	index_nb[5] = d_index3d(i,j,k+1,top_gmax);

	ur = field[2][index];
	double r, r_nb[2];
	r = cell_center[index].m_coords[2];
	r_nb[0] = cell_center[index_nb[4]].m_coords[2];
	r_nb[1] = cell_center[index_nb[5]].m_coords[2];
	

	if (FT_StateStructAtGridCrossing(front,icoords,WEST,
		comp,&intfc_state,&hs,crx_coords) &&
	        wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    utheta_nb[0] = getStateXvel(intfc_state);
	else
	    utheta_nb[0] = (field[0][index] + field[0][index_nb[0]])/2.0;

	if (FT_StateStructAtGridCrossing(front,icoords,EAST,
		comp,&intfc_state,&hs,crx_coords) &&
	        wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    utheta_nb[1] = getStateXvel(intfc_state);
	else
	    utheta_nb[1] = (field[0][index] + field[0][index_nb[1]])/2.0;

	if (FT_StateStructAtGridCrossing(front,icoords,SOUTH,
		comp,&intfc_state,&hs,crx_coords) &&
	        wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    uz_nb[0] = getStateYvel(intfc_state);
	else
	    uz_nb[0] = (field[1][index] + field[1][index_nb[2]])/2.0;

	if (FT_StateStructAtGridCrossing(front,icoords,NORTH,
		comp,&intfc_state,&hs,crx_coords) &&
	        wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    uz_nb[1] = getStateYvel(intfc_state);
	else
	    uz_nb[1] = (field[1][index] + field[1][index_nb[3]])/2.0;

	if (FT_StateStructAtGridCrossing(front,icoords,LOWER,
		comp,&intfc_state,&hs,crx_coords) &&
	        wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    rur_nb[0] = getStateZvel(intfc_state) * ((r + r_nb[0])/2.0);
	else
	    rur_nb[0] = ((field[2][index] + field[2][index_nb[4]])/2.0) * ((r + r_nb[0])/2.0);

	if (FT_StateStructAtGridCrossing(front,icoords,UPPER,
		comp,&intfc_state,&hs,crx_coords) &&
	        wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    rur_nb[1] = getStateZvel(intfc_state) * ((r + r_nb[1])/2.0);
	else
	    rur_nb[1] = ((field[2][index] + field[2][index_nb[5]])/2.0) * ((r + r_nb[1])/2.0);


        div = (utheta_nb[1] - utheta_nb[0])/(r*top_h[0]) + (uz_nb[1] - uz_nb[0])/top_h[1] + (rur_nb[1] - rur_nb[0])/(r*top_h[2]);
        return div;
}       /* end computeFieldPointDiv */

void Incompress_Solver_Smooth_3D_Cylindrical::computeFieldPointGrad(
        int *icoords,
        double *field,
        double *grad_field)
{
        int index;
        COMPONENT comp;
        int i,j,k,nb;
	int index_nb[6];
        double p_nbedge[6],p0;  //the p values on the cell edges and cell center
        double crx_coords[MAXD];
        POINTER intfc_state;
        HYPER_SURF *hs;
        GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};

	i = icoords[0];
	j = icoords[1];
	k = icoords[2];
	
	index = d_index3d(i,j,k,top_gmax);
        comp = top_comp[index];
	p0 = field[index];

	index_nb[0] = d_index3d(i-1,j,k,top_gmax);
	index_nb[1] = d_index3d(i+1,j,k,top_gmax);
	index_nb[2] = d_index3d(i,j-1,k,top_gmax);
	index_nb[3] = d_index3d(i,j+1,k,top_gmax);
	index_nb[4] = d_index3d(i,j,k-1,top_gmax);
	index_nb[5] = d_index3d(i,j,k+1,top_gmax);

	for (nb = 0; nb < 6; nb++)
	{
	    if(FT_StateStructAtGridCrossing(front,icoords,dir[nb],
			comp,&intfc_state,&hs,crx_coords) &&
		        wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    {
		if (wave_type(hs) == DIRICHLET_BOUNDARY &&
		    boundary_state(hs) == NULL)
		{
		    p_nbedge[nb] = 0.0;
		    printf("\n Now using a flow through boundary for pressure!!\n");
		}
		else
		{
		    //p_nbedge[nb] = p0;

		    if (nb == 0)
			p_nbedge[nb] = p0 - 0.5*(field[index_nb[1]] - p0);
		    else if (nb == 1)
			p_nbedge[nb] = p0 + 0.5*(p0 - field[index_nb[0]]);
		    else if (nb == 2)
			p_nbedge[nb] = p0 - 0.5*(field[index_nb[3]] - p0);
		    else if (nb == 3)
			p_nbedge[nb] = p0 + 0.5*(p0 - field[index_nb[2]]);
		    else if (nb == 4)
			p_nbedge[nb] = p0 - 0.5*(field[index_nb[5]] - p0);
		    else if (nb == 5)
			p_nbedge[nb] = p0 + 0.5*(p0 - field[index_nb[4]]);

		}
	    }
	    else
		p_nbedge[nb] = (p0 + field[index_nb[nb]])/2.0;
	}
	double r;
	r = cell_center[index].m_coords[2];

	grad_field[0] = (p_nbedge[1] - p_nbedge[0])/(r*top_h[0]);
	grad_field[1] = (p_nbedge[3] - p_nbedge[2])/top_h[1];
	grad_field[2] = (p_nbedge[5] - p_nbedge[4])/top_h[2];
}



void Incompress_Solver_Smooth_3D_Cylindrical::computeError(void)
{
    int i,j,k,index;
    double error0, error1, error2;
    double L1error_0, L2error_0, Linferror_0, L1error_1, L2error_1, Linferror_1, L1error_2, L2error_2, Linferror_2;
    double errorp;
    double L1error_p, L2error_p, Linferror_p;
    L1error_0 = 0.0; L1error_1 = 0.0; L1error_2 = 0.0;
    L2error_0 = 0.0; L2error_1 = 0.0; L2error_2 = 0.0;
    Linferror_0 = 0.0; Linferror_1 = 0.0; Linferror_2 = 0.0;
    L1error_p = 0.0; L2error_p = 0.0; Linferror_p = 0.0;
    double cellVolume;
    double r0, r1;

    for (k = kmin; k <= kmax; k++)
    for (j = jmin; j <= jmax; j++)
    for (i = imin; i <= imax; i++)
    {
	index = d_index3d(i,j,k,top_gmax);

	error0 = fabs(cell_center[index].m_state.m_U[0] - cell_center[index].m_state.m_exactU[0]);
	error1 = fabs(cell_center[index].m_state.m_U[1] - cell_center[index].m_state.m_exactU[1]);
	error2 = fabs(cell_center[index].m_state.m_U[2] - cell_center[index].m_state.m_exactU[2]);

	errorp = fabs(cell_center[index].m_state.m_P - cell_center[index].m_state.m_exactP);

	r0 = cell_center[index].m_coords[2] - top_h[2]/2.0;
	r1 = cell_center[index].m_coords[2] + top_h[2]/2.0;


	cellVolume = top_h[0]*top_h[1]*(r1*r1 - r0*r0)/2.0;
	L1error_0 = L1error_0 + error0 * cellVolume;
	L1error_1 = L1error_1 + error1 * cellVolume;
	L1error_2 = L1error_2 + error2 * cellVolume;
	L1error_p = L1error_p + errorp * cellVolume;

	L2error_0 = L2error_0 + error0*error0*cellVolume;
	L2error_1 = L2error_1 + error1*error1*cellVolume;
	L2error_2 = L2error_2 + error2*error2*cellVolume;
	L2error_p = L2error_p + errorp*errorp*cellVolume;

	if (error0 >= Linferror_0)
	    Linferror_0 = error0;
	if (error1 >= Linferror_1)
	    Linferror_1 = error1;
	if (error2 >= Linferror_2)
	    Linferror_2 = error2;
	if (errorp >= Linferror_p)
	    Linferror_p = errorp;
    }

    pp_global_max(&Linferror_0,1);
    pp_global_max(&Linferror_1,1);
    pp_global_max(&Linferror_2,1);
    pp_global_max(&Linferror_p,1);

    pp_global_sum(&L1error_0,1);
    pp_global_sum(&L1error_1,1);
    pp_global_sum(&L1error_2,1);
    pp_global_sum(&L1error_p,1);

    pp_global_sum(&L2error_0,1);
    pp_global_sum(&L2error_1,1);
    pp_global_sum(&L2error_2,1);
    pp_global_sum(&L2error_p,1);

    L2error_0 = sqrt(L2error_0);
    L2error_1 = sqrt(L2error_1);
    L2error_2 = sqrt(L2error_2);
    L2error_p = sqrt(L2error_p);

    printf("\nAt time t = %.16g, step = %d\n", front->time, front->step);
    printf("L1   error = { %.16g   %.16g   %.16g  %.16g}\n", L1error_0, L1error_1, L1error_2, L1error_p);
    printf("L2   error = { %.16g   %.16g   %.16g  %.16g}\n", L2error_0, L2error_1, L2error_2, L2error_p);
    printf("Linf error = { %.16g   %.16g   %.16g  %.16g}\n", Linferror_0, Linferror_1, Linferror_2, Linferror_p);
}


void Incompress_Solver_Smooth_3D_Cylindrical::printInteriorVelocity(char *out_name)
{
        int   i,j,k,index,totalpoints;
        double coord_theta,coord_z,coord_r;
        double vel_theta, vel_z, vel_r;
        char filename[200];
        FILE *outfile;
        sprintf(filename,"%s-visual-ts%s",out_name,
                right_flush(front->step,7));
#if defined(__MPI__)
        if (pp_numnodes() > 1)
            sprintf(filename,"%s-nd%s",filename,right_flush(pp_mynode(),4));
#endif
        sprintf(filename,"%s-liquid.vtk",filename);
        i = 0;
        if(iFparams->movie_option->plot_pres || iFparams->movie_option->plot_velo || iFparams->movie_option->plot_vort)
        {
	    outfile = fopen(filename,"w");
            fprintf(outfile,"# vtk DataFile Version 3.0\n");
            fprintf(outfile,"States of the whole computational domain\n");
            fprintf(outfile,"ASCII\n");
            fprintf(outfile,"DATASET STRUCTURED_GRID\n");
            int pointsr, pointsz, pointstheta;
        
	    pointsr = top_gmax[2] + 1;
            pointsz = top_gmax[1] + 1;
            pointstheta = top_gmax[0] + 1;
            totalpoints = pointsr*pointsz*pointstheta;
            fprintf(outfile, "DIMENSIONS %d %d %d\n", pointsr, pointsz, pointstheta);
            fprintf(outfile, "POINTS %d double\n", totalpoints);
            for(i = 0; i <= top_gmax[0]; ++i)
            for(j = 0; j <= top_gmax[1]; ++j)
            for(k = 0; k <= top_gmax[2]; ++k)
            {
                index = d_index3d(i,j,k,top_gmax);
                coord_theta = cell_center[index].m_coords[0];
                coord_z = cell_center[index].m_coords[1];
                coord_r = cell_center[index].m_coords[2];
                fprintf(outfile, "%.16g %.16g %.16g\n", coord_r*cos(coord_theta), coord_r*sin(coord_theta), coord_z);
            }
            fprintf(outfile, "POINT_DATA %i\n", totalpoints);
            if(iFparams->movie_option->plot_velo)
            {
                fprintf(outfile, "VECTORS velocity double\n");
                for(i = 0; i <= top_gmax[0]; ++i)
                for(j = 0; j <= top_gmax[1]; ++j)
                for(k = 0; k <= top_gmax[2]; ++k)
                {
                    index = d_index3d(i,j,k,top_gmax);
                    coord_theta = cell_center[index].m_coords[0];
                    coord_z = cell_center[index].m_coords[1];
                    coord_r = cell_center[index].m_coords[2];
                    //vel_theta = 0.0;
                    vel_theta = cell_center[index].m_state.m_U[0];
                    vel_z = cell_center[index].m_state.m_U[1];
                    vel_r = cell_center[index].m_state.m_U[2];
                    fprintf(outfile, "%.16g %.16g %16g\n",-vel_theta*sin(coord_theta)+vel_r*cos(coord_theta),vel_theta*cos(coord_theta)+vel_r*sin(coord_theta), vel_z);
                }
            }
            if(iFparams->movie_option->plot_pres)
            {
                fprintf(outfile, "SCALARS pressure double\n");
                fprintf(outfile, "LOOKUP_TABLE default\n");
                for(i = 0; i <= top_gmax[0]; ++i)
                for(j = 0; j <= top_gmax[1]; ++j)
                for(k = 0; k <= top_gmax[2]; ++k)
                {
                    index = d_index3d(i,j,k,top_gmax);
                    fprintf(outfile,"%.16g\n",cell_center[index].m_state.m_P);
                } 
            }
            if(iFparams->movie_option->plot_vort)
            {
            }
            fclose(outfile);
        }
}

double Incompress_Solver_Smooth_3D_Cylindrical::getCellVolume(int *icoords)
{

    int index = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
    double r = cell_center[index].m_coords[2];
    double cellVolume = r*top_h[0]*top_h[1]*top_h[2];

    return cellVolume;
}

double Incompress_Solver_Smooth_3D_Cylindrical::getFaceArea(int *icoords, GRID_DIRECTION dir)
{

    int index = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
    double r = cell_center[index].m_coords[2];

    double faceArea;

    switch(dir)
    {
    case WEST:
    case EAST:
	faceArea = top_h[1]*top_h[2];
	break;
    case SOUTH:
    case NORTH:
	faceArea = r*top_h[2]*top_h[0];
	break;
    case LOWER:
	faceArea = (r-top_h[2]/2.0)*top_h[0]*top_h[1];
	break;
    case UPPER:
	faceArea = (r+top_h[2]/2.0)*top_h[0]*top_h[1];
	break;
    default:
	assert(false);
    }
    return faceArea;

}

void Incompress_Solver_Smooth_3D_Cylindrical::getFaceCenter(int *icoords, GRID_DIRECTION dir, double faceCenter[3])
{
    double dx = 0, dy = 0,dz = 0;

    switch(dir)
    {
    case WEST:
	dx = -top_h[0];
	break;
    case EAST:
	dx =  top_h[0];
	break;
    case SOUTH:
	dy = -top_h[1];
	break;
    case NORTH:
	dy =  top_h[1];
	break;
    case LOWER:
	dz = -top_h[2];
	break;
    case UPPER:
	dz =  top_h[2];
	break;
    default:
	assert(false);
    }

    faceCenter[0] = dx/2 + cell_center[d_index3d(icoords[0],icoords[1],icoords[2],top_gmax)].m_coords[0];
    faceCenter[1] = dy/2 + cell_center[d_index3d(icoords[0],icoords[1],icoords[2],top_gmax)].m_coords[1];
    faceCenter[2] = dz/2 + cell_center[d_index3d(icoords[0],icoords[1],icoords[2],top_gmax)].m_coords[2];

}



void Incompress_Solver_Smooth_3D_Cylindrical::computeProjection_Shuqiang(void)
{
//    Incompress_Solver_Smooth_3D_Cylindrical::computeProjection();
//    return;

    int index, index_nb[6], size;
    double rhs, rho[6], rho0;
    int I,I_nb[6];
    int i,j,k,l,icoords[MAXD];
    double P_max,P_min;
    COMPONENT comp;
    double **vel = iFparams->field->vel;

    max_value = 0.0;
    double value;
    double sum_div;
    sum_div = 0.0;
    int num_iter = 0;
    double rel_residual = 0.0;

    PETSc solver;
    solver.Create(ilower, iupper-1, 7, 7);
    solver.Reset_A();
    solver.Reset_b();
    solver.Reset_x();
    size = iupper - ilower;

    setIndexMap();

    for (l = 0; l < dim; ++l)
	for (k = 0; k <= top_gmax[2]; k++)
	    for (j = 0; j <= top_gmax[1]; j++)
		for (i = 0; i <= top_gmax[0]; i++)
		{
		    index  = d_index3d(i,j,k,top_gmax);
		    vel[l][index] = cell_center[index].m_state.m_U[l];
		}

    /* Compute velocity divergence */
    for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		icoords[0] = i;
		icoords[1] = j;
		icoords[2] = k;
		index  = d_index3d(i,j,k,top_gmax);
		array[index] = computeFieldPointDiv(icoords,vel);
	    }
    scatMeshArray();
    for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {
		index  = d_index3d(i,j,k,top_gmax);
		cell_center[index].m_state.div_U = array[index];
	    }


    if(debugging("step_size"))
    {
	for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
		for (i = imin; i <= imax; i++)
		{
		    index = d_index3d(i,j,k,top_gmax);
		    value = fabs(cell_center[index].m_state.div_U);
		    sum_div = sum_div + cell_center[index].m_state.div_U * cell_center[index].m_coords[2] * top_h[0]*top_h[1]*top_h[2];
		    if(value > max_value)
			max_value = value;
		}
	pp_global_sum(&sum_div,1);
	printf("\nThe summation of divergence of U is %.16g\n",sum_div);
	pp_global_max(&max_value,1);
	printf("\nThe max value of divergence of U is %.16g\n",max_value);
	max_value = 0.0;
    }

//    double normal[6][3] = {{-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}};
//    double flux[6];

    double dh[3] = {top_h[0],top_h[1], top_h[2]};

    for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		index  = d_index3d(i,j,k,top_gmax);
		comp = top_comp[index];
		I = ijk_to_I[i][j][k];
		if (I == -1) continue;

		index_nb[0] = d_index3d(i-1,j,k,top_gmax);
		index_nb[1] = d_index3d(i+1,j,k,top_gmax);
		index_nb[2] = d_index3d(i,j-1,k,top_gmax);
		index_nb[3] = d_index3d(i,j+1,k,top_gmax);
		index_nb[4] = d_index3d(i,j,k-1,top_gmax);
		index_nb[5] = d_index3d(i,j,k+1,top_gmax);
		I_nb[0] = ijk_to_I[i-1][j][k];
		I_nb[1] = ijk_to_I[i+1][j][k];
		I_nb[2] = ijk_to_I[i][j-1][k];
		I_nb[3] = ijk_to_I[i][j+1][k];
		I_nb[4] = ijk_to_I[i][j][k-1];
		I_nb[5] = ijk_to_I[i][j][k+1];
		icoords[0] = i;
		icoords[1] = j;
		icoords[2] = k;

		rho0   = cell_center[index].m_state.m_rho;
		for (l = 0; l < 6; ++l)
		{
		    if (I_nb[l] == -1)
			index_nb[l] = index;
		    rho[l] = 1.0/2*(rho0 + cell_center[index_nb[l]].m_state.m_rho);
		    //coeff[l] = 1.0/rho[l]/sqr(top_h[l/2]);
		}

//		double rr = sqr(cell_center[index].m_coords[2]);
		double r = cell_center[index].m_coords[2];
		double cellVolume = getCellVolume(icoords);
		double faceArea;

		rhs = 0;
		// WEST, theta
		faceArea = getFaceArea(icoords,WEST);
		if(I_nb[0]>=0)
		{
		    solver.Add_A(I, I,       -1/r*1/dh[0]*faceArea /rho[0]);
		    solver.Add_A(I, I_nb[0],  1/r*1/dh[0]*faceArea /rho[0]);
		}
		else
		    ;
		// EAST
		faceArea = getFaceArea(icoords,EAST);
		if(I_nb[1]>=0)
		{
		    solver.Add_A(I, I,       -1/r*1/dh[0]*faceArea /rho[1]);
		    solver.Add_A(I, I_nb[1],  1/r*1/dh[0]*faceArea /rho[1]);
		}
		else
		    ;
		// SOUTH, z
		faceArea = getFaceArea(icoords,SOUTH);
		if(I_nb[2]>=0)
		{
		    solver.Add_A(I, I,       -1  *1/dh[1]*faceArea /rho[2]);
		    solver.Add_A(I, I_nb[2],  1  *1/dh[1]*faceArea /rho[2]);
		}
		else
		    ;
		// NORTH
		faceArea = getFaceArea(icoords,NORTH);
		if(I_nb[3]>=0)
		{
		    solver.Add_A(I, I,       -1  *1/dh[1]*faceArea /rho[3]);
		    solver.Add_A(I, I_nb[3],  1  *1/dh[1]*faceArea /rho[3]);
		}
		else
		    ;
		// LOWER, r
		faceArea = getFaceArea(icoords,LOWER);
		if(I_nb[4]>=0)
		{
		    solver.Add_A(I, I,       -1  *1/dh[2]*faceArea /rho[4]);
		    solver.Add_A(I, I_nb[4],  1  *1/dh[2]*faceArea /rho[4]);
		}
		else
		    ;
		// UPPER
		faceArea = getFaceArea(icoords,UPPER);
		if(I_nb[5]>=0)
		{
		    solver.Add_A(I, I,       -1  *1/dh[2]*faceArea /rho[5]);
		    solver.Add_A(I, I_nb[5],  1  *1/dh[2]*faceArea /rho[5]);
		}
		else
		    ;
		rhs = cell_center[index].m_state.div_U/accum_dt * cellVolume;
		solver.Add_b(I, rhs);


	    }

    solver.SetMaxIter(40000);
    solver.SetTol(1e-14);

    start_clock("Before Petsc Solver in Projection step");
    solver.Solve_withPureNeumann();

//    solver.Print_A(NULL);
//    solver.Print_b(NULL);
//    exit(0);


    solver.GetNumIterations(&num_iter);
    solver.GetFinalRelativeResidualNorm(&rel_residual);
    if(rel_residual > 1)
    {
	printf("\n The solution diverges! The residual is %g. Solve again using GMRES!\n",rel_residual);
	solver.Reset_x();
	solver.Solve_withPureNeumann_GMRES();
	solver.GetNumIterations(&num_iter);
	solver.GetFinalRelativeResidualNorm(&rel_residual);
    }
    stop_clock("After Petsc Solver in Projection step");

    double *x;
    FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));
    solver.Get_x(x);

    if (debugging("PETSc"))
	(void) printf("Incompress_Solver_Smooth_3D_Cylindrical::"
		"computeProjection: "
		"num_iter = %d, rel_residual = %g \n",
		num_iter, rel_residual);

    P_max = -HUGE;		P_min = HUGE;
    for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		index = d_index3d(i,j,k,top_gmax);
		I = ijk_to_I[i][j][k];
		array[index] = x[I-ilower];
	    }
    scatMeshArray();
    for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {
		index  = d_index3d(i,j,k,top_gmax);
		cell_center[index].m_state.m_phi = array[index];
	    }

    if(debugging("step_size"))
    {
	for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
		for (i = imin; i <= imax; i++)
		{
		    index = d_index3d(i,j,k,top_gmax);
		    value = fabs(cell_center[index].m_state.m_phi);
		    if (value > max_value)
			max_value = value;
		}
	pp_global_max(&max_value,1);
	printf("\nThe max value of phi is %.16g\n",max_value);
    }

    free_these(1,x);
}

void Incompress_Solver_Smooth_3D_Cylindrical::compDiffWithSmoothProperty_2nd_decoupled_Shuqiang(void)
{
    COMPONENT comp;
    int index,index_nb[6],size;
    int I,I_nb[6];
    double crx_coords[MAXD];
    double mu,rho; //,U0_nb[6],U1_nb[6],U2_nb[6],U0_center,U1_center,U2_center;
    L_STATE source_term, U_nb[6], U_nb_new[6], U_center, rhs;

    int i,j,k,l,nb,icoords[MAXD];
    double speed;
    double *x;
    GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
    POINTER intfc_state;
    HYPER_SURF *hs;
    int num_iter;
    double rel_residual;

    max_speed = 0.0;
    setIndexMap();

    size = iupper - ilower;
    FT_VectorMemoryAlloc((POINTER*)&x,3*size,sizeof(double));

    PETSc solver;
    solver.Create(3*ilower, 3*iupper-1, 10, 10);
    // 7theta + 2r  for the first equation
    // 7z for the second equation
    // 7r + 2theta for the third equation

    solver.Reset_A();
    solver.Reset_b();
    solver.Reset_x();

    printf("\nIn diffusion solver ,m_dt = %.16g\n", m_dt);

    double dtheta, dz, dr; //, faceCenter[3];
    dtheta = top_h[0];
    dz	   = top_h[1];
    dr     = top_h[2];

    double dh[]  = {dtheta,dz,dr};

    for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		I  = ijk_to_I[i][j][k];
		if (I == -1) continue;

		index  = d_index3d(i,j,k,top_gmax);
		//6 neighbours of the center cell
		index_nb[0] = d_index3d(i-1,j,k,top_gmax);
		index_nb[1] = d_index3d(i+1,j,k,top_gmax);
		index_nb[2] = d_index3d(i,j-1,k,top_gmax);
		index_nb[3] = d_index3d(i,j+1,k,top_gmax);
		index_nb[4] = d_index3d(i,j,k-1,top_gmax);
		index_nb[5] = d_index3d(i,j,k+1,top_gmax);


		//6 neighbours of the center cell
		I_nb[0] = ijk_to_I[i-1][j][k];
		I_nb[1] = ijk_to_I[i+1][j][k];
		I_nb[2] = ijk_to_I[i][j-1][k];
		I_nb[3] = ijk_to_I[i][j+1][k];
		I_nb[4] = ijk_to_I[i][j][k-1];
		I_nb[5] = ijk_to_I[i][j][k+1];


		icoords[0] = i;
		icoords[1] = j;
		icoords[2] = k;
		comp = top_comp[index];

		mu = cell_center[index].m_state.m_mu;
		rho = cell_center[index].m_state.m_rho;
		U_center = cell_center[index].m_state;

		for (nb = 0; nb < 6; nb++)
		{
		    if (FT_StateStructAtGridCrossing(front,icoords,dir[nb],
			    comp,&intfc_state,&hs,crx_coords) &&
			    wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
		    {
			// old boundary condition
			U_nb[nb].m_U[0] = getStateVel[0](intfc_state);
			U_nb[nb].m_U[1] = getStateVel[1](intfc_state);
			U_nb[nb].m_U[2] = getStateVel[2](intfc_state);

			FT_StateStructAtGridCrossing(front,icoords,dir[nb],
				comp,&intfc_state,&hs,crx_coords);
			U_nb_new[nb].m_U[0] = getStateVel[0](intfc_state);
			U_nb_new[nb].m_U[1] = getStateVel[1](intfc_state);
			U_nb_new[nb].m_U[2] = getStateVel[2](intfc_state);

			assert(I_nb[nb]<0);
		    }
		    else
		    {
			U_nb[nb] = cell_center[index_nb[nb]].m_state;
		    }
		}

		double r = cell_center[index].m_coords[2];
		double cellVolume = getCellVolume(icoords);
		double faceArea;

		for(nb = 0; nb<6; nb++)
		{
		    faceArea = getFaceArea(icoords,dir[nb]);

		    if(I_nb[nb]>=0)
		    {
			compDiffWithSmoothProperty_cellFace(
				&solver,icoords,I,I_nb[nb],
				dir[nb],
				dh[nb/2],
				faceArea,cellVolume,r,mu,rho,
				U_nb[nb],U_nb_new[nb],U_center);
		    }
		    else
			compDiffWithSmoothProperty_Dirichlet(
				&solver,icoords,I,I_nb[nb],
				dir[nb],
				dh[nb/2],
				faceArea,cellVolume,r,mu,rho,
				U_nb[nb],U_nb_new[nb],U_center);
		}

		compDiffWithSmoothProperty_cellInterior(
			&solver,icoords,I,cellVolume,r,mu,rho,U_center);

	    }

    solver.SetMaxIter(40000);
    solver.SetTol(1e-15);

    start_clock("Before Petsc Solve");
    solver.Solve_GMRES();
    solver.GetNumIterations(&num_iter);
    solver.GetFinalRelativeResidualNorm(&rel_residual);
    stop_clock("After Petsc Solve");

    // get back the solution
    solver.Get_x(x);

    if (debugging("PETSc"))
	(void) printf("Incompress_Solver_Smooth_3D_Cylindrical::compDiffWithSmoothProperty_2nd_decoupled_Shuqiang: "
		"num_iter = %d, rel_residual = %g. \n",
		num_iter,rel_residual);


    for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		I = ijk_to_I[i][j][k];
		index = d_index3d(i,j,k,top_gmax);
		if (I >= 0)
		{
		    cell_center[index].m_state.m_U[0] = x[I*3-ilower*3];
		    cell_center[index].m_state.m_U[1] = x[I*3-ilower*3+1];
		    cell_center[index].m_state.m_U[2] = x[I*3-ilower*3+2];
		    speed = fabs(cell_center[index].m_state.m_U[0]) +
			    fabs(cell_center[index].m_state.m_U[1]) +
			    fabs(cell_center[index].m_state.m_U[2]);
		    if (speed > max_speed)
			max_speed = speed;
		}
		else
		{
		    cell_center[index].m_state.m_U[0] = 0.0;
		    cell_center[index].m_state.m_U[1] = 0.0;
		    cell_center[index].m_state.m_U[2] = 0.0;
		}
	    }
    for (l = 0; l < 3; ++l)
    {
	for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
		for (i = imin; i <= imax; i++)

		{
		    index  = d_index3d(i,j,k,top_gmax);
		    array[index] = cell_center[index].m_state.m_U[l];
		}
	scatMeshArray();
	for (k = 0; k <= top_gmax[2]; k++)
	    for (j = 0; j <= top_gmax[1]; j++)
		for (i = 0; i <= top_gmax[0]; i++)

		{
		    index  = d_index3d(i,j,k,top_gmax);
		    cell_center[index].m_state.m_U[l] = array[index];
		}
    }
    pp_global_max(&max_speed,1);

    free_these(1,x);
}       /* end compDiffWithSmoothProperty3d */

void Incompress_Solver_Smooth_3D_Cylindrical::compDiffWithSmoothProperty_cellFace(
	PETSc *pSolver,
	int *icoords,
	int I, int I_nb,
	GRID_DIRECTION dir,
	double dh,
	double faceArea,
	double cellVolume,
	double r,
	double mu,
	double rho,
	L_STATE &U_nb,
	L_STATE &U_nb_new,
	L_STATE &U_center)
{
    double mymu = 1.0/2*m_dt*mu/rho;

    double coeff = 1;
    if(dir==WEST || dir==EAST)
	coeff *= 1.0/r;

    // used for first order derivative with theta
    int sign = 0;
    if(dir==WEST)
	sign = -1;
    else if(dir==EAST)
	sign = 1;

    // 1st equ
    // u^{*}
    pSolver->Add_A(I*3, I*3,     - coeff/dh*faceArea*(-mymu));
    pSolver->Add_A(I*3, I_nb*3,    coeff/dh*faceArea*(-mymu));
    pSolver->Add_b(I*3,          - coeff/dh*faceArea*( mymu) * U_center.m_U[0]);
    pSolver->Add_b(I*3,            coeff/dh*faceArea*( mymu) * U_nb.m_U[0]);

    if(sign!=0)
    {
	pSolver->Add_A(I*3, I*3+2,    -sign/(r*r)/dh*cellVolume*(-mymu));
	pSolver->Add_A(I*3, I_nb*3+2,  sign/(r*r)/dh*cellVolume*(-mymu));
	pSolver->Add_b(I*3,           -sign/(r*r)/dh*cellVolume*( mymu) * U_center.m_U[2]);
	pSolver->Add_b(I*3,            sign/(r*r)/dh*cellVolume*( mymu) * U_nb.m_U[2]);
    }

    // 2nd equ
    pSolver->Add_A(I*3+1, I*3+1,     - coeff/dh*faceArea*(-mymu));
    pSolver->Add_A(I*3+1, I_nb*3+1,    coeff/dh*faceArea*(-mymu));
    pSolver->Add_b(I*3+1,            - coeff/dh*faceArea*( mymu) * U_center.m_U[1]);
    pSolver->Add_b(I*3+1,              coeff/dh*faceArea*( mymu) * U_nb.m_U[1]);

    // 3rd equ
    pSolver->Add_A(I*3+2, I*3+2,     - coeff/dh*faceArea*(-mymu));
    pSolver->Add_A(I*3+2, I_nb*3+2,    coeff/dh*faceArea*(-mymu));
    pSolver->Add_b(I*3+2,            - coeff/dh*faceArea*( mymu) * U_center.m_U[2]);
    pSolver->Add_b(I*3+2,              coeff/dh*faceArea*( mymu) * U_nb.m_U[2]);

    if(sign!=0)
    {
	sign = -sign;	// this makes the following code similar with 1st equ.
	pSolver->Add_A(I*3+2, I*3,    -sign/(r*r)/dh*cellVolume*(-mymu));
	pSolver->Add_A(I*3+2, I_nb*3,  sign/(r*r)/dh*cellVolume*(-mymu));
	pSolver->Add_b(I*3+2,         -sign/(r*r)/dh*cellVolume*( mymu) * U_center.m_U[0]);
	pSolver->Add_b(I*3+2,          sign/(r*r)/dh*cellVolume*( mymu) * U_nb.m_U[0]);
    }
}
/**
 * TODO: modify the following code for boundary.
 */
void Incompress_Solver_Smooth_3D_Cylindrical::compDiffWithSmoothProperty_Dirichlet(
	PETSc *pSolver,
	int *icoords,
	int I, int I_nb,
	GRID_DIRECTION dir,
	double dh,
	double faceArea,
	double cellVolume,
	double r,
	double mu,
	double rho,
	L_STATE &U_nb,
	L_STATE &U_nb_new,
	L_STATE &U_center)
{
    double mymu = 1.0/2*m_dt*mu/rho;

    double coeff = 1;
    if(dir==WEST || dir==EAST)
	coeff *= 1.0/r;

    // used for first order derivative with theta
    int sign = 0;
    if(dir==WEST)
	sign = -1;
    else if(dir==EAST)
	sign = 1;

    //-----------------------------------------------------
    // 			1st equ
    //-----------------------------------------------------
    // u^{*}
    pSolver->Add_A(I*3, I*3,     - coeff/dh*faceArea*(-mymu) * 2);
//  pSolver->Add_A(I*3, I_nb*3,    coeff/dh*faceArea*(-mymu));
    pSolver->Add_b(I*3,            coeff/dh*faceArea*( mymu) * 2*U_nb_new.m_U[0]);
    pSolver->Add_b(I*3,          - coeff/dh*faceArea*( mymu) * 2*U_center.m_U[0]);
    pSolver->Add_b(I*3,            coeff/dh*faceArea*( mymu) * 2*U_nb.m_U[0]);

    if(sign!=0)
    {
	pSolver->Add_A(I*3, I*3+2,    -sign/(r*r)/dh*cellVolume*(-mymu) * 2);
//	pSolver->Add_A(I*3, I_nb*3+2,  sign/(r*r)/dh*cellVolume*(-mymu));
	pSolver->Add_b(I*3,            sign/(r*r)/dh*cellVolume*( mymu) * 2*U_nb_new.m_U[2]);
	pSolver->Add_b(I*3,           -sign/(r*r)/dh*cellVolume*( mymu) * 2*U_center.m_U[2]);
	pSolver->Add_b(I*3,            sign/(r*r)/dh*cellVolume*( mymu) * 2*U_nb.m_U[2]);
    }

    //-----------------------------------------------------
    // 			2nd equ
    //-----------------------------------------------------
    pSolver->Add_A(I*3+1, I*3+1,     - coeff/dh*faceArea*(-mymu) * 2);
//  pSolver->Add_A(I*3+1, I_nb*3+1,    coeff/dh*faceArea*(-mymu));
    pSolver->Add_b(I*3+1, 	       coeff/dh*faceArea*( mymu) * 2*U_nb_new.m_U[1]);
    pSolver->Add_b(I*3+1,            - coeff/dh*faceArea*( mymu) * 2*U_center.m_U[1]);
    pSolver->Add_b(I*3+1,              coeff/dh*faceArea*( mymu) * 2*U_nb.m_U[1]);

    //-----------------------------------------------------
    // 			3rd equ
    //-----------------------------------------------------
    pSolver->Add_A(I*3+2, I*3+2,     - coeff/dh*faceArea*(-mymu) * 2);
//  pSolver->Add_A(I*3+2, I_nb*3+2,    coeff/dh*faceArea*(-mymu));
    pSolver->Add_b(I*3+2,              coeff/dh*faceArea*( mymu) * 2*U_nb_new.m_U[2]);
    pSolver->Add_b(I*3+2,            - coeff/dh*faceArea*( mymu) * 2*U_center.m_U[2]);
    pSolver->Add_b(I*3+2,              coeff/dh*faceArea*( mymu) * 2*U_nb.m_U[2]);

    if(sign!=0)
    {
	sign = -sign;	// this makes the following code similar with 1st equ.
	pSolver->Add_A(I*3+2, I*3,    -sign/(r*r)/dh*cellVolume*(-mymu) * 2);
//	pSolver->Add_A(I*3+2, I_nb*3,  sign/(r*r)/dh*cellVolume*(-mymu));
	pSolver->Add_b(I*3+2,          sign/(r*r)/dh*cellVolume*( mymu) * 2*U_nb_new.m_U[0]);
	pSolver->Add_b(I*3+2,         -sign/(r*r)/dh*cellVolume*( mymu) * 2*U_center.m_U[0]);
	pSolver->Add_b(I*3+2,          sign/(r*r)/dh*cellVolume*( mymu) * 2*U_nb.m_U[0]);
    }
}
void Incompress_Solver_Smooth_3D_Cylindrical::compDiffWithSmoothProperty_cellInterior(
	PETSc *pSolver,
	int *icoords,
	int I,
	double cellVolume,
	double r,
	double mu,
	double rho,
	L_STATE &U_center)
{
    int index  = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
    double coords[3], rhs;
    double mymu = 1.0/2*m_dt*mu/rho;


    L_STATE source_term;

    getRectangleCenter(index, coords);
    computeSourceTerm(coords, source_term);


    //-----------------------------------------------------
    // 			1st equ
    //-----------------------------------------------------
    pSolver->Add_A(3*I, 3*I, 1*cellVolume);
    pSolver->Add_b(3*I,      1*cellVolume * U_center.m_U[0]);
    pSolver->Add_A(3*I, 3*I, 1/(r*r)*cellVolume*( mymu));
    pSolver->Add_b(3*I,      1/(r*r)*cellVolume*(-mymu) * U_center.m_U[0]);

    rhs  = m_dt*source_term.m_U[0]*cellVolume;
    rhs += m_dt*cell_center[index].m_state.f_surf[0]/rho*cellVolume;
    rhs -= m_dt * cell_center[index].m_state.grad_q[0]/rho*cellVolume;
    rhs -= m_dt * cell_center[index].m_state.m_adv[0] *cellVolume;
    pSolver->Add_b(3*I, rhs);

    //-----------------------------------------------------
    // 			2nd equ
    //-----------------------------------------------------
    pSolver->Add_A(3*I+1, 3*I+1, 1*cellVolume);
    pSolver->Add_b(3*I+1,        1*cellVolume * U_center.m_U[1]);

    rhs  = m_dt*source_term.m_U[1]*cellVolume;
    rhs += m_dt*cell_center[index].m_state.f_surf[1]/rho*cellVolume;
    rhs -= m_dt * cell_center[index].m_state.grad_q[1]/rho*cellVolume;
    rhs -= m_dt * cell_center[index].m_state.m_adv[1] *cellVolume;
    pSolver->Add_b(3*I+1, rhs);

    //-----------------------------------------------------
    // 			3rd equ
    //-----------------------------------------------------
    pSolver->Add_A(3*I+2, 3*I+2, 1*cellVolume);
    pSolver->Add_b(3*I+2,        1*cellVolume * U_center.m_U[2]);
    pSolver->Add_A(3*I+2, 3*I+2, 1/(r*r)*cellVolume*( mymu));
    pSolver->Add_b(3*I+2,        1/(r*r)*cellVolume*(-mymu) * U_center.m_U[2]);
    rhs  = m_dt*source_term.m_U[2]*cellVolume;
    rhs += m_dt*cell_center[index].m_state.f_surf[2]/rho*cellVolume;
    rhs -= m_dt * cell_center[index].m_state.grad_q[2]/rho*cellVolume;
    rhs -= m_dt * cell_center[index].m_state.m_adv[2] *cellVolume;
    pSolver->Add_b(3*I+2, rhs);
}


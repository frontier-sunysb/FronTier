/*******************************************************************
 * 		CARTESIAN.c
 *******************************************************************/
/*
      To-Do: Line 392, the values of ilower and iupper 
      
*/

#include "vcartsn.h"
#include "solver.h"

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
	int i,j,k,index;
	double crds[MAXD];
	int icoords[MAXD];
	int size;
	int cell_index;
	RECTANGLE       rectangle;
	double hmin = HUGE;

	if (debugging("trace")) printf("Entering initMesh()\n");

	FT_MakeGridIntfc(front);
	setDomain();

	size = 1;
	for (i = 0; i < dim; ++i)
	{
	    if (hmin > top_h[i]) hmin = top_h[i];
	    size *= (top_gmax[i] + 1);
	}
	
	switch (dim)
	{
	case 2:
	    FT_ScalarMemoryAlloc((POINTER*)&field,sizeof(VF_FIELD));
	    FT_MatrixMemoryAlloc((POINTER*)&field->vel,2,size,sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&field->eta,size,sizeof(double));
	    imin = (lbuf[0] == 0) ? 1 : lbuf[0];
            jmin = (lbuf[1] == 0) ? 1 : lbuf[1];
            imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
            jmax = (ubuf[1] == 0) ? top_gmax[1] - 1 : top_gmax[1] - ubuf[1];
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
	}

	setComponent();
	FT_FreeGridIntfc(front);
	if (debugging("trace")) printf("Leaving initMesh()\n");
}

void CARTESIAN::setComponent(void)
{
	int i;

        double *coords;
        int *icoords;
        int size = (int)cell_center.size();
        HYPER_SURF_ELEMENT *hse;
        HYPER_SURF *hs;
        double t[MAXD],point[MAXD];
        int n;

        for (i = 0; i < size; i++)
        {
	    cell_center[i].comp = top_comp[i];
        }
}	/* end setComponent */

void CARTESIAN::computeAdvection()
{
}


void CARTESIAN::solve(double dt)
{
	if (debugging("trace")) printf("Entering solve()\n");
	m_dt = dt;
	start_clock("solve");

	setDomain();
        if (debugging("trace")) printf("Passing setDomain()\n");

	setComponent();
	if (debugging("trace")) printf("Passing setComponent()\n");

	computeAdvection();
	if (debugging("trace")) printf("Passing liquid computeAdvection()\n");

	setAdvectionDt();
	if (debugging("trace")) printf("Passing setAdvectionDt()\n");

	stop_clock("solve");
	if (debugging("trace")) printf("Leaving solve()\n");
}


void CARTESIAN::setAdvectionDt()
{
	double D,Dl,Ds;
	static double m_dt_expl,m_dt_impl;  // explicit and implicit time step
	static boolean first = YES;

	if (first)
	{
	    first = NO;
	    min_dt = 0.1*sqr(hmin)/D/(double)dim;
	}

	// For smooth transition to implicit step
	double tstep = (double)front->step;
	double smooth_factor; 
	smooth_factor = 1.0/(1.0 + sqr(tstep/20.0));
	m_dt = m_dt_impl - (m_dt_impl - m_dt_expl)*smooth_factor;
	/* For Crank Nicolson scheme*/
	m_dt = m_dt_expl;

	if (debugging("trace"))
	{
	    printf("In setAdvectionDt: m_dt = %24.18g min_dt = %f\n",
				m_dt,min_dt);
	}
}	/* end setAdvectionDt */

void CARTESIAN::getVelocity(double *p, double *U)
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

void CARTESIAN::getRectangleIndex(int index, int &i, int &j)
{
	i = cell_center[index].icoords[0];
	j = cell_center[index].icoords[1];
}

void CARTESIAN::getRectangleIndex(int index, int &i, int &j, int &k)
{
	i = cell_center[index].icoords[0];
	j = cell_center[index].icoords[1];
	k = cell_center[index].icoords[2];
}


int CARTESIAN::getRectangleComponent(int index)
{	
	return getComponent(cell_center[index].icoords);
}

void CARTESIAN::getRectangleCenter(
	int index, 
	double *coords)
{
	int i;
	for (i = 0; i < dim; ++i)
	    coords[i] = cell_center[index].coords[i];
}

void CARTESIAN::getRectangleCenter(
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

int CARTESIAN::getComponent(
	double *p)
{
	return component(p,front->interf);
}

int CARTESIAN::getComponent(
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

CARTESIAN::CARTESIAN(Front &front):front(&front)
{
}

void CARTESIAN::setDomain()
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
	hmin = top_h[0];
	for (i = 1; i < dim; ++i)
	    if (hmin > top_h[i]) hmin = top_h[i];

}	/* end setDomain */

void CARTESIAN::initDrawVariables()
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
	    CursorAfterString(infile,"Type y to make movie of velocity:");
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'Y' || string[0] == 'y')
	    {
		if (set_bound)
		{
		    CursorAfterString(infile,"Enter min and max velocity:");
                    fscanf(infile,"%lf %lf",&var_min,&var_max);
                    (void) printf("%f %f\n",var_min,var_max);
		}
		FT_AddHdfMovieVariable(front,set_bound,YES,SOLID_COMP,
				"xvel",0,field->vel[0],getStateXvel,
				var_max,var_min);
		FT_AddHdfMovieVariable(front,set_bound,YES,SOLID_COMP,
				"yvel",0,field->vel[1],getStateYvel,
				var_max,var_min);
	    }
	    if (CursorAfterStringOpt(infile,
		"Type y to make movie of eta:"))
	    {
            	fscanf(infile,"%s",string);
            	(void) printf("%s\n",string);
            	if (string[0] == 'Y' || string[0] == 'y')
	    	{
		    if (set_bound)
		    {
		        CursorAfterString(infile,
				"Enter min and max eta:");
                        fscanf(infile,"%lf %lf",&var_min,&var_max);
                        (void) printf("%f %f\n",var_min,var_max);
		    }
		    FT_AddHdfMovieVariable(front,set_bound,YES,SOLID_COMP,
				"eta",0,field->eta,getStateEta,
				var_max,var_min);
		    FT_AddVtkScalarMovieVariable(front,"eta",field->eta);
	    	}
	    }
	    break;
	}
	/* Added for vtk movie of vector field */
	CursorAfterString(infile,"Type y to make vector velocity field movie:");
        fscanf(infile,"%s",string);
        (void) printf("%s\n",string);
        if (string[0] == 'Y' || string[0] == 'y')
	    FT_AddVtkVectorMovieVariable(front,"VELOCITY",field->vel);

	fclose(infile);
}	/* end initMovieVariables */

void CARTESIAN::initFlowState()
{
	int i,j;
	FT_MakeGridIntfc(front);
	setDomain();
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
	}
}	/* end initFlowState */

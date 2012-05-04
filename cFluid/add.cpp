
EXPORT boolean FT_NearestGridPointVar(
	Front *front,
	COMPONENT comp,
	double *coords,
	double *grid_array,
	int range,
	double *ans)
{
        INTERFACE *grid_intfc = front->grid_intfc;
	RECT_GRID *gr = &topological_grid(grid_intfc);
	struct Table *T = table_of_interface(grid_intfc);
        COMPONENT *top_comp = T->components;
	int top_gmax = gr->gmax;
	double *L = gr->L;
	double *h = gr->h;
	double grid_p[MAXD];
	double dist,min_dist = HUGE;
        int i,j,k,dim = gr->dim;
	int imin[MAXD],imax[MAXD];
	int icoords[MAXD];
	boolean grid_pt_found = NO;

	if (!rect_in_which(coords,icoords,gr))
        {
            *ans = 0.0;
            return NO;
        }
	for (i = 0; i < dim; ++i)
	{
	    imin[i] = max(icoords[i]-range,0);
	    imax[i] = min(icoords[i]+range+1,top_gmax[i]+1);
	}
	switch (dim)
	{
	case 1:
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index1d(i,top_gmax);
		if (top_comp[index] != comp) continue;
		grid_p[0] = L[0] + i*h[0];
		dist = distance_between_positions(coords,grid_p,dim);
		if (dist < min_dist)
		{
		    min_dist = dist;
		    icoords[0] = i;    
		    grid_pt_found = YES;
		}
	    }
	    break;
	case 2:
	    for (i = imin[0]; i <= imax[0]; ++i)
	    for (j = imin[1]; j <= imax[1]; ++j)
	    {
		index = d_index2d(i,j,top_gmax);
		if (top_comp[index] != comp) continue;
		grid_p[0] = L[0] + i*h[0];
		grid_p[1] = L[1] + j*h[1];
		dist = distance_between_positions(coords,grid_p,dim);
		if (dist < min_dist)
		{
		    min_dist = dist;
		    icoords[0] = i;    
		    icoords[1] = j;    
		    grid_pt_found = YES;
		}
	    }
	    break;
	case 3:
	    for (i = imin[0]; i <= imax[0]; ++i)
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (k = imin[2]; k <= imax[2]; ++k)
	    {
		index = d_index3d(i,j,k,top_gmax);
		if (top_comp[index] != comp) continue;
		grid_p[0] = L[0] + i*h[0];
		grid_p[1] = L[1] + j*h[1];
		grid_p[2] = L[2] + k*h[2];
		dist = distance_between_positions(coords,grid_p,dim);
		if (dist < min_dist)
		{
		    min_dist = dist;
		    icoords[0] = i;    
		    icoords[1] = j;    
		    icoords[2] = k;    
		    grid_pt_found = YES;
		}
	    }
	}
	if (grid_pt_found == YES)
	{
	    index = d_index(icoords,top_gmax,3);
            *ans = grid_array[index];
            return YES;
	}
	else
	{
            *ans = 0.0;
            return NO;
	}
}	/* end FT_NearestGridPointFieldVar */

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

/*******************************************************************
 * 		G_CARTESIAN.c
 *******************************************************************/
#include "cFluid.h"

static double (*getStateMom[MAXD])(Locstate) =
               {getStateXmom,getStateYmom,getStateZmom};

void G_CARTESIAN::readBaseFront(int i)
{
        char *dir_name = eqn_params->base_dir_name;
        int RestartStep = eqn_params->steps[i];
        F_BASIC_DATA *f_basic = eqn_params->f_basic;
        char restart_state_name[200];

        FT_ScalarMemoryAlloc((POINTER*)&base_front,sizeof(Front));
        f_basic->RestartRun = YES;
        f_basic->size_of_intfc_state = sizeof(STATE);

        sprintf(f_basic->restart_name,"%s/intfc-ts%s",dir_name,
                        right_flush(RestartStep,7));
        sprintf(restart_state_name,"%s/state.ts%s",dir_name,
                        right_flush(RestartStep,7));
	if (debugging("base_files"))
	{
	    (void) printf("Currently only handle one processor:\n");
	    (void) printf("Open base file names:\n");
            (void) printf("base intfc name = %s\n",f_basic->restart_name);
            (void) printf("base state name = %s\n",restart_state_name);
	}

        FT_StartUp(base_front,f_basic);

        readBaseStates(restart_state_name);
}       /* end readBaseFront */

void G_CARTESIAN::readBaseStates(
        char *restart_name)
{
        FILE *infile;
        int i,j,k,l,index;
        char fname[100];
        double *u;
        int *base_gmax;
	int *base_comp;
	Table *T;
        RECT_GRID *base_grid;
	int base_size;
	STATE st_tmp;
	int comp;
	EOS_PARAMS *eos = eqn_params->eos;
	static boolean first = YES;

	base_front->extra1 = (POINTER)eqn_params;
        readFrontStates(base_front,restart_name);

        sprintf(fname,"%s-gas",restart_name);
        infile = fopen(fname,"r");

        /* Initialize states in the interior regions */

        FT_MakeGridIntfc(base_front);
        base_grid = &topological_grid(base_front->grid_intfc);
	T = table_of_interface(base_front->grid_intfc);
	base_comp = T->components;
        base_gmax = base_grid->gmax;
        FT_ScalarMemoryAlloc((POINTER*)&base_field,sizeof(FIELD));
        next_output_line_containing_string(infile,"Interior gas states:");
	if (first)
	{
            switch (dim)
            {
            case 1:
	    	base_size = base_gmax[0] + 1;
            	FT_VectorMemoryAlloc((POINTER*)&base_field->dens,base_size,
				FLOAT);
            	FT_VectorMemoryAlloc((POINTER*)&base_field->pres,base_size,
				FLOAT);
            	FT_VectorMemoryAlloc((POINTER*)&base_field->engy,base_size,
				FLOAT);
            	FT_MatrixMemoryAlloc((POINTER*)&base_field->momn,1,base_size,
				FLOAT);
	    	break;
            case 2:
	    	base_size = (base_gmax[0] + 1)*(base_gmax[1] + 1);
            	FT_VectorMemoryAlloc((POINTER*)&base_field->dens,base_size,
				FLOAT);
            	FT_VectorMemoryAlloc((POINTER*)&base_field->pres,base_size,
				FLOAT);
            	FT_VectorMemoryAlloc((POINTER*)&base_field->engy,base_size,
				FLOAT);
            	FT_MatrixMemoryAlloc((POINTER*)&base_field->momn,2,base_size,
				FLOAT);
	    	break;
            case 3:
	    	base_size = (base_gmax[0] + 1)*(base_gmax[1] + 1)*
				(base_gmax[2] + 1);
            	FT_VectorMemoryAlloc((POINTER*)&base_field->dens,base_size,
				FLOAT);
            	FT_VectorMemoryAlloc((POINTER*)&base_field->pres,base_size,
				FLOAT);
            	FT_VectorMemoryAlloc((POINTER*)&base_field->engy,base_size,
				FLOAT);
            	FT_MatrixMemoryAlloc((POINTER*)&base_field->momn,3,base_size,
				FLOAT);
	    	break;
	    }
	}

        switch (dim)
        {
        case 1:
            for (i = 0; i <= base_gmax[0]; ++i)
            {
                index = d_index1d(i,base_gmax);
		comp = base_comp[index];
                st_tmp.eos = &(eos[comp]);

                fscanf(infile,"%lf",&base_field->dens[index]);
                fscanf(infile,"%lf",&base_field->engy[index]);
                st_tmp.dens = base_field->dens[index];
                st_tmp.engy = base_field->engy[index];
                for (l = 0; l < dim; ++l)
                {
                    fscanf(infile,"%lf",&base_field->momn[l][index]);
                    st_tmp.momn[l] = base_field->momn[l][index];
                }
                base_field->pres[index] = EosPressure(&st_tmp);
            }
            break;
        case 2:
            for (i = 0; i <= base_gmax[0]; ++i)
            for (j = 0; j <= base_gmax[1]; ++j)
            {
                index = d_index2d(i,j,base_gmax);
		comp = base_comp[index];
                st_tmp.eos = &(eos[comp]);

                fscanf(infile,"%lf",&base_field->dens[index]);
                fscanf(infile,"%lf",&base_field->engy[index]);
                st_tmp.dens = base_field->dens[index];
                st_tmp.engy = base_field->engy[index];
                for (l = 0; l < dim; ++l)
                {
                    fscanf(infile,"%lf",&base_field->momn[l][index]);
                    st_tmp.momn[l] = base_field->momn[l][index];
                }
                base_field->pres[index] = EosPressure(&st_tmp);
            }
            break;
        case 3:
            for (i = 0; i <= base_gmax[0]; ++i)
            for (j = 0; j <= base_gmax[1]; ++j)
            for (k = 0; k <= base_gmax[2]; ++k)
            {
                index = d_index3d(i,j,k,base_gmax);
		comp = base_comp[index];
                st_tmp.eos = &(eos[comp]);

                fscanf(infile,"%lf",&base_field->dens[index]);
                fscanf(infile,"%lf",&base_field->engy[index]);
                st_tmp.dens = base_field->dens[index];
                st_tmp.engy = base_field->engy[index];
                for (l = 0; l < dim; ++l)
                {
                    fscanf(infile,"%lf",&base_field->momn[l][index]);
                    st_tmp.momn[l] = base_field->momn[l][index];
                }
                base_field->pres[index] = EosPressure(&st_tmp);
            }
            break;
        }
        fclose(infile);
}       /* end readBaseStates */

void G_CARTESIAN::compareWithBaseData(
        char *out_name)
{
	double L1[6],L2[6],LI[6];
	int ip = front->ip - 2;
	int i,j,k,n,l,index,N;
	double v,v_base;
	double *dens = field.dens;
	double **momn = field.momn;
	double *engy = field.engy;
	int comp;
	double coords[MAXD];
	FILE *dfile,*efile,*mfile,*lfile;
	char fname[200];

	if (ip >= eqn_params->num_step) return;
	readBaseFront(ip);

	for (i = 0; i < 6; ++i)
	    L1[i] = L2[i] = LI[i] = 0.0;
	switch (dim)
	{
	case 1:
	    N = 0;
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index1d(i,top_gmax);
		comp = top_comp[index];
		coords[0] = top_L[0] + i*top_h[0];
		n = 0;
		v = dens[index];
		FT_IntrpStateVarAtCoords(base_front,comp,coords,
			base_field->dens,getStateDens,&v_base,NULL);
		L1[n] += fabs(v - v_base);
		L2[n] += sqr(v - v_base);
		if (LI[n] < fabs(v - v_base)) 
		    LI[n] = fabs(v - v_base);
		n++;
		v = engy[index];
		FT_IntrpStateVarAtCoords(base_front,comp,coords,
			base_field->engy,getStateEngy,&v_base,NULL);
		L1[n] += fabs(v - v_base);
		L2[n] += sqr(v - v_base);
		if (LI[n] < fabs(v - v_base)) 
		    LI[n] = fabs(v - v_base);
		n++;
		for (l = 0; l < dim; ++l)
		{
		    v = momn[l][index];
		    FT_IntrpStateVarAtCoords(base_front,comp,coords,
			    base_field->momn[l],getStateMom[l],&v_base,NULL);
		    L1[n] += fabs(v - v_base);
		    L2[n] += sqr(v - v_base);
		    if (LI[n] < fabs(v - v_base)) 
		        LI[n] = fabs(v - v_base);
		    n++;
		}
		N++;
	    }
	    break;
	case 2:
	    for (i = imin[0]; i <= imax[0]; ++i)
	    for (j = imin[1]; j <= imax[1]; ++j)
	    {
		index = d_index2d(i,j,top_gmax);
		comp = top_comp[index];
		coords[0] = top_L[0] + i*top_h[0];
		n = 0;
		v = dens[index];
		FT_IntrpStateVarAtCoords(base_front,comp,coords,
			base_field->dens,getStateDens,&v_base,NULL);
		L1[n] += fabs(v - v_base);
		L2[n] += sqr(v - v_base);
		if (LI[n] < fabs(v - v_base)) 
		    LI[n] = fabs(v - v_base);
		n++;
		v = engy[index];
		FT_IntrpStateVarAtCoords(base_front,comp,coords,
			base_field->engy,getStateEngy,&v_base,NULL);
		L1[n] += fabs(v - v_base);
		L2[n] += sqr(v - v_base);
		if (LI[n] < fabs(v - v_base)) 
		    LI[n] = fabs(v - v_base);
		n++;
		for (l = 0; l < dim; ++l)
		{
		    v = momn[l][index];
		    FT_IntrpStateVarAtCoords(base_front,comp,coords,
			    base_field->momn[l],getStateMom[l],&v_base,NULL);
		    L1[n] += fabs(v - v_base);
		    L2[n] += sqr(v - v_base);
		    if (LI[n] < fabs(v - v_base)) 
		        LI[n] = fabs(v - v_base);
		    n++;
		}
	    }
	    break;
	case 3:
	    for (i = imin[0]; i <= imax[0]; ++i)
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (k = imin[2]; k <= imax[2]; ++k)
	    {
		index = d_index3d(i,j,k,top_gmax);
		comp = top_comp[index];
		coords[0] = top_L[0] + i*top_h[0];
		n = 0;
		v = dens[index];
		FT_IntrpStateVarAtCoords(base_front,comp,coords,
			base_field->dens,getStateDens,&v_base,NULL);
		L1[n] += fabs(v - v_base);
		L2[n] += sqr(v - v_base);
		if (LI[n] < fabs(v - v_base)) 
		    LI[n] = fabs(v - v_base);
		n++;
		v = engy[index];
		FT_IntrpStateVarAtCoords(base_front,comp,coords,
			base_field->engy,getStateEngy,&v_base,NULL);
		L1[n] += fabs(v - v_base);
		L2[n] += sqr(v - v_base);
		if (LI[n] < fabs(v - v_base)) 
		    LI[n] = fabs(v - v_base);
		n++;
		for (l = 0; l < dim; ++l)
		{
		    v = momn[l][index];
		    FT_IntrpStateVarAtCoords(base_front,comp,coords,
			    base_field->momn[l],getStateMom[l],&v_base,NULL);
		    L1[n] += fabs(v - v_base);
		    L2[n] += sqr(v - v_base);
		    if (LI[n] < fabs(v - v_base)) 
		        LI[n] = fabs(v - v_base);
		    n++;
		}
	    }
	    break;
	}

	for (i = 0; i < dim+2; ++i)
	{
	    L1[i] /= N;
	    L2[i] /= N;		L2[i] = sqrt(L2[i]);
	}
	sprintf(fname,"%s/L-norm-%d",out_name,ip);
	lfile = fopen(fname,"w");
	fprintf(lfile,"Density Norm\n");
	fprintf(lfile,"L1 = %16.12f  L2 = %16.12f  LI = %16.12f\n\n",
			L1[0],L2[0],LI[0]);
	fprintf(lfile,"Energy Norm\n");
	fprintf(lfile,"L1 = %16.12f  L2 = %16.12f  LI = %16.12f\n\n",
			L1[1],L2[1],LI[1]);
	fprintf(lfile,"Momentum Norm\n");
	fprintf(lfile,"L1 = %16.12f  L2 = %16.12f  LI = %16.12f\n\n",
			L1[2],L2[2],LI[2]);
	fclose(lfile);
}	/* end compareWithBaseData */

void G_CARTESIAN::freeBaseFront()
{
	int i = front->ip - 2;
	if (i >= eqn_params->num_step) return;
	FT_FreeFront(base_front);
}	/* end freeBaseFront */

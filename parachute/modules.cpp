/***************************************************************
FronTier is a set of libraries that implements differnt types of 
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

#include <iFluid.h>
#include <airfoil.h>

extern void initParachuteModules(Front *front)
{
	int i,num_canopy;
	FILE *infile = fopen(InName(front),"r");
	SURFACE *surf;

	if (debugging("trace"))
	    (void) printf("Entering initParachuteModules()\n");

	if (debugging("set_module"))
	    gview_plot_interface("module-step-1",front->interf);

	CursorAfterString(infile,"Enter number of canopy surfaces:");
        fscanf(infile,"%d",&num_canopy);
        (void) printf("%d\n",num_canopy);

	for (i = 0; i < num_canopy; ++i)
	{
	    CgalCanopySurface(infile,front,&surf);
	}

	if (debugging("trace"))
	    (void) printf("Leaving initParachuteModules()\n");
}	/* end initParachuteModules */

extern void initParachuteDefault(
	Front *front)
{
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	af_params->is_parachute_system = YES;
	af_params->num_opt_round = 20;
        af_params->spring_model = MODEL1;
	af_params->gore_len_fac = 1.0;
}	/* end initParachuteDefault */

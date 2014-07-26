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

/**********************************************************************
 * 		ifluid_basic.h
 **********************************************************************/

#ifndef _FT_IFLUID_BASIC_H_
#define _FT_IFLUID_BASIC_H_

enum _IF_PROB_TYPE {
        ERROR_TYPE = -1,
        TWO_FLUID_BUBBLE = 1,
        TWO_FLUID_RT,
        TWO_FLUID_KH,
        FLUID_SOLID_CIRCLE,
        FLUID_SOLID_RECT,
        FLUID_SOLID_TRIANGLE,
        BUBBLE_SURFACE,
        FLUID_RIGID_BODY,
        ROTOR_ONE_FLUID,
        ROTOR_TWO_FLUID,
	CHANNEL_FLOW,
	RANDOM_FLOW,
        FLUID_CRYSTAL,
	TAYLOR_GREEN_VORTEX,
        FLUID_SOLID_CYLINDER
};
typedef enum _IF_PROB_TYPE IF_PROB_TYPE;

extern void restart_set_dirichlet_bdry_function(Front*);
extern void iF_flowThroughBoundaryState(double*,HYPER_SURF*,Front*,POINTER,
                        POINTER);
extern void iF_timeDependBoundaryState(double*,HYPER_SURF*,Front*,POINTER,
                        POINTER);
extern void ifluid_point_propagate(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
extern void ifluid_compute_force_and_torque(Front*,CURVE*,double,double*,
                        double*);
extern void setInitialIntfc(Front*,LEVEL_FUNC_PACK*,char*,IF_PROB_TYPE);
extern void init_fluid_state_func(Incompress_Solver_Smooth_Basis*,IF_PROB_TYPE);
extern void read_iFparams(char*,IF_PARAMS*);
extern void read_iF_prob_type(char*,IF_PROB_TYPE*);
extern void recordBdryEnergyFlux(Front*,char*);
#endif

/**********************************************************************
 * 		ifluid_state.h
 **********************************************************************/

#pragma once
struct _STATE {
        double dens;                    /* Density */
        double pres;                    /* Pressure */
        double phi;                     /* Potential */
        double vel[MAXD];               /* Velocities */
        double vort;                    /* Vorticity in 2D */
        double impulse[MAXD];            /* Accum impact from external force */
	double solute;			/* For subsurface problem */
	double temperature;             /* For melting with flow problem */
	double vapor;                   /* For climate problem */
        double supersat;		/* For climate problem */
	double mu;			/* For eddy viscosity */
};
typedef struct _STATE STATE;

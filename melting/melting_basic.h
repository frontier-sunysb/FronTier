/**********************************************************************
 * 		melting_basic.h					      *
 **********************************************************************/

#ifndef _FT_MELTING_BASIC_H_
#define _FT_MELTING_BASIC_H_


struct _STATE {
	double dens;                    /* Density */
        double pres;                    /* Pressure */
        double phi;                     /* Potential */
        double vel[MAXD];               /* Velocities */
        double vort;                    /* Vorticity in 2D */
        double impulse[MAXD];           /* Accum impact from external force */
        double solute;                  /* For subsurface problem */
	double temperature;		/* For melting with flow problem */
};
typedef struct _STATE STATE;

#endif

/*
 * weno.h
 *
 *  Created on: Jun 12, 2011
 *      Author: yli
 */

#ifndef WENO_H_
#define WENO_H_

#include <FronTier.h>

enum _EQN_TYPE {
        UNKNOWN_EQN_TYPE = -1,
        LINEAR_EQN = 1,
        BURGER_EQN
};
typedef enum _EQN_TYPE EQN_TYPE;

enum _INIT_TYPE {
        UNKNOWN_INIT_TYPE = -1,
	SQUARE = 1,
        HAMP,
	COSINE,
        WAVE
};
typedef enum _INIT_TYPE INIT_TYPE;

struct _PARAMS {
	EQN_TYPE eqn_type;
	INIT_TYPE init_type;
	double a;		/* a in u_t+au_x=0 or u_t+a(u^2)_x=0 */
};
typedef struct _PARAMS PARAMS;

void setEquationType(_PARAMS PARAMS);
void Weno5(int mesh_size, double *u_old, double *u_new, double dx, double dt);
void Runge_Kutta_4th(int mesh_size, double *u_old, double *u_new, double dx, double dt);
void TVD_RK_3th(int mesh_size, double *u_old, double *u_new, double dx, double dt);
void Weno5_Get_Flux(double *u_old, double *flux, double lambda, int mesh_size);

#endif /* WENO_H_ */

/*
 * gd.h
 *
 *  Created on: July 1, 2015
 *      Author: xlchen
 */

#ifndef DG_H_
#define DG_H_

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
void Dg5(int mesh_size, double **coef_old, double **coef_new, double dx, double dt);
void TVD_RK_3th(int mesh_size, double **coef_old, double **coef_new, double dx, double dt);
void Dg5_rhs_eval(int mesh_size, double **coef, double **rhs, double dx);

#endif /* DG_H_ */

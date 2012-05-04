/* 
 * File:   EBM2D_ELLIPTIC.h
 * Author: shuqiang (robert) wang
 *
 * Created on Nov 29, 2008.
 * 
 * This class is used for the Embedded Boundary Method for solving either 
 * elliptic boundary value or elliptic interface problem.
 *
 * The EBM2D is in fact an finite colume method with the unknowns defined at the
 * cell center. For the elliptic boundary value problem, only one unknown in
 * each cell is needed. However, for the elliptic interface problem, 1 or 2 or
 * 3 unknowns are needed at the cell center. 
 *
 * There are many different way of storing/referening the unknowns at the cell 
 * center. Therefore, it is necessary to hide the implementation so that other
 * modules are not dependent on the implementation. Even if the implementation
 * is changed, there is no need to change other modules.
 *
 * The equation to be solved is the following
 *      
 *      Div * beta Grad phi = rou
 *
 * References:
 * 1) Hans Johansen and Phillip Colella, A Cartesian Grid Embedded Boundary 
 *    Method for Poisson's Equation on Irregular Domains, JCP 147, 60-85 (1998).
 *
 * See also lpetsc.c.
 *
 *
 * The classes relation:
 *
 * For elliptic solver:
 *      EBM2D_CELL      <--- EBM2D_ELLIPTIC_STATE;
 *      EBM2D_CARTESIAN <--- EBM2D_ELLIPTIC.
 *
 * For Navier-Stokes solver:
 *      EBM2D_CELL      <--- EBM2D_ELLIPTIC_STATE <--- EBM2D_NS_STATE;
 *      EBM2D_CARTESIAN <--- EBM2D_ELLIPTIC       <--- EBM2D_NS.
 *
 * Possibly for Gas solver:
 *      EBM2D_CELL      <--- EBM2D_GAS_STATE;
 *      EBM2D_CARTESIAN <--- EBM2D_GAS.
 *
 */
#ifndef _EBM2D_ELLIPTIC_H
#define _EBM2D_ELLIPTIC_H
#include <EBM2D.h>

class EBM2D_ELLIPTIC_STATE: public EBM2D_CELL {
public:

};

class EBM2D_ELLIPTIC: public EBM2D_CARTESIAN {
public:

};



#endif //_EBM2D_ELLIPTIC_H

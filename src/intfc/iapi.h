/*! \file iapi.h
    
    \brief The iapi.h contains the functions used to operate the interface.
 */

/*! \defgroup CURVE    Curve Functions
/*! \defgroup SURFACE  SURFACE Functions
 **/

#include <intfc/int.h>

                /* Interface Function Prototypes*/

#if defined(c_plusplus) || defined(__cplusplus)
extern "C" {
#endif

/*! \fn void I_MoveNodeToPoint(POINT *pt, CURVE *curve)
    \ingroup CURVE
    \brief For closed curve, move the node to the given point. Do nothing
     if curve is not a closed one or the point is not on curve.
    \param pt @b in Address of point the closed node to be moved to
    \param curve @b inout Address of curve to be operated on
 */
   IMPORT  void I_MoveNodeToPoint(POINT *pt,
			CURVE *curve );

/*! \fn CURVE **I_SplitCurve(POINT *pt, CURVE *curve)
    \ingroup CURVE
    \brief Split a curve into two curves at the given point pt, return
     two curves as a pointer array. Return NULL if the point is not on
     the input curve.
    \param pt @b in Address of point the curve is to split at
    \param curve @b inout Address of curve to be operated on
 */
   IMPORT  CURVE **I_SplitCurve(POINT *pt,
			CURVE *curve );

/*! \fn void I_SmoothSurfColor(SURFACE *surf,int num_rounds)
    \ingroup SURFACE
    \brief A utility function to smooth the color of the surface
     by averaging the color with neighboring triangles.
    \param surf @b inout Address of surface to smooth color
    \param num_round @b in number of round of smoothing
 */
   IMPORT  void I_SmoothSurfColor(SURFACE *surf,
			int num_rounds);

/*! \fn SURFACE *I_CopySurface(SURFACE *surf)
    \ingroup SURFACE
    \brief This function return a copy of the input surface on the
     same interface structure.
    \param surf @b in Input of the surface.
 */
   IMPORT  SURFACE *I_CopySurface(SURFACE *surf);

/*! \fn void I_ShiftSurface(SURFACE *surf,double *displacement)
    \ingroup SURFACE
    \brief This function shift each point of the surface by
     given displacement.
    \param surf @b inout Input of the surface to be operated.
    \param displacement @b in displacement each point is to be shifted.
 */
   IMPORT  void I_ShiftSurface(SURFACE *surf,
			double *displacement);

#if defined(c_plusplus) || defined(__cplusplus)
}
#endif

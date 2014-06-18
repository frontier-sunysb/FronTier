/*! \file iapi.h
    
    \brief The iapi.h contains the functions used to operate the interface.
 */

/*! \defgroup POINT      Point Functions
/*! \defgroup CURVE      Curve Functions
/*! \defgroup SURFACE    SURFACE Functions
/*! \defgroup INTERFACE  INTERFACE Functions
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

/*! \fn SURFACE *I_AddTwoSurfaces(SURFACE *surf1,SURFACE *surf2)
    \ingroup SURFACE
    \brief This function add tris of surf2 to surf1, delete surf2 and
     return surf1, no other operations is performed.
    \param surf1 @b in Input of the surface 1.
    \param surf2 @b in Input of the surface 2.
 */
   IMPORT  SURFACE *I_AddTwoSurfaces(SURFACE *surf1,
			SURFACE *surf2);

/*! \fn void I_TransInteriorIntfcPoints(INTERFACE *intfc,double *displacement)
    \ingroup INTERFACE
    \brief This function makes translation of all interior points in
     the interface.
    \param intfc @b inout Input of the interface.
    \param displacement @b in Displacement of points.
 */
   IMPORT  void I_TransInteriorIntfcPoints(INTERFACE *intfc,
			double *displacement);

/*! \fn void I_SphericalRoratePoint(POINT *p,double *center,double phi,double theta,boolean first)
    \ingroup POINT
    \brief This function rotate the point about the center with spherical angle.
    \param p @b inout Point to be rotated.
    \param center @b in Center of the rotation.
    \param phi @b in Azimuthal angle.
    \param theta @b in Polar angle.
    \param first @b in Flag if rotation parameters is first uesd.
 */
   IMPORT  void I_SphericalRotatePoint(
			POINT *p,
        		double *center,                 /* Rotation center */
        		double phi,                     /* Azimuthal angle */
        		double theta,                   /* Polar angle */
        		boolean first);

/*! \fn void I_PolarRoratePoint(POINT *p,double *center,double phi,boolean first)
    \ingroup POINT
    \brief This function rotates the point about the center with polar angle.
    \param p @b inout Point to be rotated.
    \param center @b in Center of the rotation.
    \param phi @b in Polar angle.
    \param first @b in Flag if rotation parameters is first uesd.
 */
   IMPORT  void I_PolarRotatePoint(
			POINT *p,
        		double *center,                 /* Rotation center */
        		double phi,                     /* Polar angle */
        		boolean first);

/*! \fn void I_RotatePointAboutAxis(POINT *p,double *dir,double *axis,double phi)
    \ingroup POINT
    \brief This function rotate the point about a line as axis.
    \param p @b inout Point to be rotated.
    \param dir @b in Direction vector of rotation line of axis.
    \param axis @b in A point on the line of axis.
    \param phi @b in Angle of rotation in counter-clock wise direction.
 */
   IMPORT  void I_RotatePointAboutAxis(
			POINT *p,
        		double *dir,                
        		double *axis,                
        		double phi);   

/*! \fn void I_SphericalRorateInteriorIntfcPoints(INTERFACE *intfc,double *center,double phi,double theta)
    \ingroup POINT
    \brief This function rotates the interior points of the input interface 
     about the center with spherical angle (bounding curves or surfaces will
     not be rotated.
    \param intfc @b inout Interface whose interior points to be rotated.
    \param center @b in Center of the rotation.
    \param phi @b in Azimuthal angle.
    \param theta @b in Polar angle.
 */
   IMPORT  void I_SphericalRotateInteriorIntfcPoints(
			INTERFACE *intfc,
        		double *center,                 /* Rotation center */
        		double phi,                     /* Azimuthal angle */
        		double theta);                   /* Polar angle */

/*! \fn void I_SphericalRorateInteriorSurfPoints(SURFACE *surf,double *center,double phi,double theta)
    \ingroup POINT
    \brief This function rotates the interior points of the input surface 
     about the center with spherical angle (no curve point is rotated).
    \param surf @b inout Surface whose interior points to be rotated.
    \param center @b in Center of the rotation.
    \param phi @b in Azimuthal angle.
    \param theta @b in Polar angle.
 */
   IMPORT  void I_SphericalRotateInteriorSurfPoints(
			SURFACE *surf,
        		double *center,                 /* Rotation center */
        		double phi,                     /* Azimuthal angle */
        		double theta);                   /* Polar angle */

/*! \fn void I_SphericalRorateInteriorCurvePoints(CURVE *curve,double *center,double phi,double theta)
    \ingroup POINT
    \brief This function rotates the interior points of the input surface 
     about the center with spherical angle (no node point is rotated).
    \param curve @b inout Curve whose interior points to be rotated.
    \param center @b in Center of the rotation.
    \param phi @b in Azimuthal angle.
    \param theta @b in Polar angle.
 */
   IMPORT  void I_SphericalRotateInteriorCurvePoints(
			CURVE *curve,
        		double *center,                 /* Rotation center */
        		double phi,                     /* Azimuthal angle */
        		double theta);                   /* Polar angle */

/*! \fn int I_NumOfSurfInteriorPoints(SURFACE *surf)
 *  \ingroup QUERY
    \brief This function count number of interior points on the surface
     boundary (on curve or node) points are not included.
    \param surf @b in   Pointer to a surface of the front interface.
 */

   IMPORT  int I_NumOfSurfInteriorPoints(SURFACE *surf);

/*! \fn int I_NumOfCurveInteriorPoints(CURVE *curve)
 *  \ingroup QUERY
    \brief This function count number of interior points on the curve
     boundary (on node) points are not included.
    \param curve @b in   Pointer to a curve of the front interface.
 */

   IMPORT  int I_NumOfCurveInteriorPoints(CURVE *curve);

/*! \fn void I_FoldSurface(SURFACE *surf,double *dir,double *axis,double angle,SIDE side,boolean first)
    \ingroup SURFACE
    \brief This function fold the surface around the axis described by 
     dir (direction) and axis (a point on the axis) with input angle.
     Only points on the input side of the axis will be rotated. The
     function does its job of the surface is flat and axis is a line
     on the surface.
    \param surf @b inout Input of the surface to be folded.
    \param dir @b in direction of the folding axis.
    \param axis @b in a point on the folding axis.
    \param angle @b in angle of folding.
    \param side @b in indicating point on the side to be rotated (folded).
    \param first @b in if true, reset all the point of the surface.
 */
   IMPORT  void I_FoldSurface(SURFACE *surf,
			double *dir,
			double *axis,
			double angle,
			SIDE side,
			boolean first);

/*! \fn boolean I_SewSurface(SURFACE *surf,double *crds_start,double *crds_end)
    \ingroup SURFACE
    \brief This function switches the surface along a line segment starting
     from crds_start to crds_end. The sewing line must be along existing 
     curves with tolerance, else return NO.
    \param surf @b inout Input of the surface to be sewed.
    \param crds_start @b in start coordinates of sewing line segment.
    \param crds_end @b in end coordinates of sewing line segment.
 */
   IMPORT  boolean I_SewSurface(SURFACE *surf,
			double *crds_start,
			double *crds_end);

#if defined(c_plusplus) || defined(__cplusplus)
}
#endif

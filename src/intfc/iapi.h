/*! \file iapi.h
    
    \brief The iapi.h contains the functions used to operate the interface.
 */

/*! \defgroup POINT      Point Functions
/*! \defgroup CURVE      Curve Functions
/*! \defgroup SURFACE    SURFACE Functions
/*! \defgroup INTERFACE  INTERFACE Functions
/*! \defgroup QUERY  	 QUERY Functions
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
    \param surf @b in   Pointer to an input surface.
 */

   IMPORT  int I_NumOfSurfInteriorPoints(SURFACE *surf);

/*! \fn int I_NumOfCurveInteriorPoints(CURVE *curve)
 *  \ingroup QUERY
    \brief This function count number of interior points on the curve
     boundary (on node) points are not included.
    \param curve @b in   Pointer to a curve.
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

/*! \fn int I_NumOfIntfcCurves(INTERFACE *intfc)
 *  \ingroup QUERY
    \brief This function count number of curves in the interface
     structure. It returns an integer for total number of curves.
    \param intfc @b in   Pointer to the interface.
 */

   IMPORT  int I_NumOfIntfcCurves(INTERFACE *intfc);

/*! \fn int I_NumOfIntfcSurfaces(INTERFACE *intfc)
 *  \ingroup QUERY
    \brief This function count number of surfaces in the interface
     structure. It returns an integer for total number of surfaces.
    \param intfc @b in   Pointer to the interface.
 */

   IMPORT  int I_NumOfIntfcSurfaces(INTERFACE *intfc);

/*! \fn int I_NumOfIntfcNodes(INTERFACE *intfc)
 *  \ingroup QUERY
    \brief This function count number of nodes in the interface
     structure. It returns an integer for total number of nodes.
    \param intfc @b in   Pointer to the interface.
 */

   IMPORT  int I_NumOfIntfcNodes(INTERFACE *intfc);

/*! \fn int I_NumOfNodeCurves(NODE *node)
 *  \ingroup QUERY
    \brief This function count number of curves attached to the node,
     including both in_curves and out_curves.
    \param node @b in   Pointer to a node of the interface.
 */

   IMPORT  int I_NumOfNodeCurves(NODE *node);

/*! \fn int I_NumOfCurveSurfaces(CURVE *curve)
 *  \ingroup QUERY
    \brief This function count number of surfaces attached to the curve,
     including both pos_surfaces and neg_surfaces.
    \param curve @b in   Pointer to a curve of the interface.
 */

   IMPORT  int I_NumOfCurveSurfaces(CURVE *curve);

/*! \fn int I_ComponentAtCoords(double *coords,INTERFACE *intfc)
 *  \ingroup QUERY
    \brief This function is a topological function which returns
     the component index at the input coordinates.
    \param coords @b in   Double array of coordinates.
    \param intfc @b in   Pointer to the interface.
 */

   IMPORT  int I_ComponentAtCoords(
			double *coords,
			INTERFACE *intfc);

/*! \fn int I_NumOfIntfcBonds(INTERFACE *intfc)
 *  \ingroup QUERY
    \brief This function count number of bonds on the entire interface.
    \param intfc @b in  Pointer to the interface.
 */

   IMPORT  int I_NumOfIntfcBonds(INTERFACE *intfc);

/*! \fn int I_NumOfIntfcPoints(INTERFACE *intfc)
 *  \ingroup QUERY
    \brief This function count number of point on the entire interface.
    \param intfc @b in  Pointer to the interface.
 */

   IMPORT  int I_NumOfIntfcPoints(INTERFACE *intfc);

/*! \fn int I_NumOfCurveBonds(CURVE *curve)
 *  \ingroup QUERY
    \brief This function count number of bond on a curve.
    \param curve @b in  Pointer to curve.
 */

   IMPORT  int I_NumOfCurveBonds(CURVE *curve);

/*! \fn int I_NumOfCurvePoints(CURVE *curve)
 *  \ingroup QUERY
    \brief This function count number of point on the curve including
     its nodes (closed node may be counted twice).
    \param curve @b in  Pointer to a curve.
 */

   IMPORT  int I_NumOfCurvePoints(CURVE *curve);

/*! \fn void I_ArrayOfCurves(INTERFACE *intfc,CURVE **curves)
 *  \ingroup QUERY
    \brief This function put all curves in the interface to an array,
     assuming the memory of the pointer array has already been allocated.
    \param intfc @b in  Pointer to the interface.
    \param curves @b inout  Curve array (memory allocated).
 */

   IMPORT  void I_ArrayOfCurves(INTERFACE *intfc,CURVE **curves);

/*! \fn void I_ArrayOfSurfaces(INTERFACE *intfc,SURFACE **surfs)
 *  \ingroup QUERY
    \brief This function put all surfaces in the interface to an array,
     assuming the memory of the pointer array has already been allocated.
    \param intfc @b in  Pointer to the interface.
    \param surfs @b inout  Surface array (memory allocated).
 */

   IMPORT  void I_ArrayOfSurfaces(INTERFACE *intfc,SURFACE **surfs);

/*! \fn int I_NumOfIntfcTris(INTERFACE *intfc)
 *  \ingroup QUERY
    \brief This function count number of triangles on the entire interface.
    \param intfc @b in  Pointer to the interface.
 */

   IMPORT  int I_NumOfIntfcTris(INTERFACE *intfc);

/*! \fn int I_NumOfSurfTris(SURFACE *surf)
 *  \ingroup QUERY
    \brief This function count number of triangles on the surface.
    \param surf @b in   Pointer to a surface.
 */

   IMPORT  int I_NumOfSurfTris(SURFACE *surf);

/*! \fn int I_NumOfSurfPoints(SURFACE *surf)
 *  \ingroup QUERY
    \brief This function count number of points on the surface.
    \param surf @b in   Pointer to a surface.
 */

   IMPORT  int I_NumOfSurfPoints(SURFACE *surf);

/*! \fn int I_NumOfSurfaces(INTERFACE *intfc)
 *  \ingroup QUERY
    \brief This function returns number of surfaces on the interface.
    \param intfc @b in  Pointer to the interface.
 */

   IMPORT  int I_NumOfSurfaces(INTERFACE *intfc);

/*! \fn void I_ArrayOfIntfcCurves(INTERFACE *intfc, CURVE **curve_array)
 *  \ingroup QUERY
    \brief This function put all the handles (pointers) of curves in the 
     intfc to an array (already allocated with memory) curve_array.
    \param intfc @b in	Pointer to an interface.
    \param curve_array @b inout curve array (with memory allocated).
 */

   IMPORT  void I_ArrayOfIntfcCurves(INTERFACE *intfc ,
				CURVE **curve_array);

/*! \fn void I_ArrayOfCurvePoints(CURVE *curve, POINT **point_array)
 *  \ingroup QUERY
    \brief This function put all the handles (pointers) of points on the 
     curve to an array (already allocated with memory) point_array.
    \param curve @b in	Pointer to a curve of the interface.
    \param point_array @b inout point array (with memory allocated).
 */

   IMPORT  void I_ArrayOfCurvePoints(CURVE *curve,
				POINT **point_array);

/*! \fn int I_NumOfSurfPoints(SURFACE *surf)
 *  \ingroup QUERY
    \brief This function count number of point on the surface.
    \param surf @b in	Pointer to a surface of the interface.
 */

   IMPORT  int I_NumOfSurfPoints(SURFACE *surf);

/*! \fn int I_NumOfIntfcBonds(INTERFACE *intfc)
 *  \ingroup QUERY
    \brief This function count number of bond on the entire interface.
    \param intfc @b in	Pointer to the interface.
 */

   IMPORT  int I_NumOfIntfcBonds(INTERFACE *intfc);

/*! \fn void I_ArrayOfNodeCurves(NODE *node, CURVE **curve_array)
 *  \ingroup QUERY
    \brief This function put all the handles (pointers) of curves on the 
     node to an array (already allocated with memory) curve_array.
    \param node @b in	Pointer to a node of the interface.
    \param curve_array @b inout curve array (with memory allocated).
 */

   IMPORT  void I_ArrayOfNodeCurves(NODE *node,
				CURVE **curve_array);

/*! \fn void I_ArrayOfCurveBonds(CURVE *curve, BOND **bond_array)
 *  \ingroup QUERY
    \brief This function put all the handles (pointers) of bonds on the 
     curve to an array (already allocated with memory) bond_array.
    \param curve @b in	Pointer to a curve of the interface.
    \param bond_array @b inout bond array (with memory allocated).
 */

   IMPORT  void I_ArrayOfCurveBonds(CURVE *curve,
				BOND **bond_array);

/*! \fn void I_ArrayOfSurfPoints(SURFACE *surf, POINT **point_array)
 *  \ingroup QUERY
    \brief This function put all the handles (pointers) of points on the 
     surface to an array (already allocated with memory) point_array.
    \param surf @b in	Pointer to a surface of the interface.
    \param point_array @b inout point array (with memory allocated).
 */

   IMPORT  void I_ArrayOfSurfPoints(SURFACE *surf,
				POINT **point_array);


/*! \fn void I_ArrayOfSurfCurves(SURFACE *surf, CURVE **curve_array)
 *  \ingroup QUERY
    \brief This function put all the handles (pointers) of curves on the 
     surface to an array (already allocated with memory) curve_array.
    \param surf @b in	Pointer to a surface of the interface.
    \param curve_array @b inout curve array (with memory allocated).
 */

   IMPORT  void I_ArrayOfSurfCurves(SURFACE *surf, CURVE **curve_array);

/*! \fn int I_NumOfSurfCurves(SURFACE *surf)
 *  \ingroup QUERY
    \brief This function count and return number of curves on the surface.
    \param surf @b in	Pointer to a surface of the interface.
 */

   IMPORT  int I_NumOfSurfCurves(SURFACE *surf);

/*! \fn int I_FirstRingTrisAroundPoint(POINT *p, TRI *tri, TRI ***tris)
 *  \ingroup QUERY
    \brief This function searches for triangles in the first ring around
     the input point p, and return the number of tris in the first ring.
    \param p @b in	Pointer to a point of the interface.
    \param tri @b in	Pointer to one of the tris around the point.
    \param tris @b out	Pointer to array of tris around the point.
 */

   IMPORT  int I_FirstRingTrisAroundPoint(POINT *p, TRI *tri, TRI ***tris);

/*! \fn CURVE* I_CurveOfPoint(INTERFACE *intfc, POINT *point, BOND **bond)
 *  \ingroup QUERY
    \brief This function looks for the curve on which the point 
     is located. If found, it will return the handle (pointer) of 
     the curve, otherwise it will return NULL.
    \param intfc @b in	Pointer to the interface.
    \param point @b in point on which the curve to be returned.
    \param bond @b out bond on which the point is at.
 */

   IMPORT  CURVE* I_CurveOfPoint(INTERFACE *intfc,POINT *point, BOND **bond);

/*! \fn NODE* I_NodeOfPoint(INTERFACE *intfc, POINT *point)
 *  \ingroup QUERY
    \brief This function looks for the node on which the point 
     is located. If found, it will return the handle (pointer) of 
     the node, otherwise it will return NULL.
    \param intfc @b in	Pointer to the interface.
    \param point @b in point on which the curve to be returned.
 */

   IMPORT  NODE* I_NodeOfPoint(INTERFACE *intfc, POINT *point);
#if defined(c_plusplus) || defined(__cplusplus)
}
#endif

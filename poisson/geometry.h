/*******************************************************************************
 *                geometry.h
 * used for 2D/3D geometrical calculations.
 * 
 * Is namespace a better choice than static function of a class?
 ******************************************************************************/

#ifndef _GEOMETRY_H
#define _GEOMETRY_H
class GEOMETRY {

     static const double m_tol = 1E-12;
public:
     static boolean intersectLineLine2D(double p1[2], double p2[2],
				     double p3[2], double p4[2],
				     double crx[2]);   
     

     void test(void);  // test for all member functions

};

#endif




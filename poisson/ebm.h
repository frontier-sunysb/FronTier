/******************************************************************************
 * File:   ebm.h
 * Author: shuqiang (robert) wang
 * 
 * Created on Dec 12, 2008.
 *
 * This file is used for the Embedded Boundary Method for solving eithe elliptic
 * boundary value or elliptic interface problem.
 *
 * See also: 
 *          EBM2D.h for 2D implementation;
 *          EBM3D.h for 3D implementation.
 *
 * Here are defined some enums & classes shared by EBM2D & EBM3D.
 *  
 ******************************************************************************/
#ifndef _EBM_H
#define _EBM_H

#define EBM_DEBUG   // used to for debugging.

/******************************************************************************
 *            EBM_FACE
 ******************************************************************************/
enum EBM_FACE 
{
     FACE_WEST  = 0, 
     FACE_EAST  = 1, 
     FACE_SOUTH = 2, 
     FACE_NORTH = 3, 
     FACE_DOWN  = 4, 
     FACE_UP    = 5
};

/******************************************************************************
 *             EBM_PLANE
 * Used to represent the coordinate planes.
 * PLANE_X means the plane perpendicular to the x coordinate.
 * PLANE_Y, PLANE_Z has similar meaning.
 *****************************************************************************/
enum EBM_PLANE 
{
     PLANE_X=0, 
     PLANE_Y=2, 
     PLANE_Z=4
};


/****************************************************************************
 *             EBM_COORD
 ***************************************************************************/
enum EBM_COORD
{
     COORD_X = 0,
     COORD_Y = 1,
     COORD_Z = 2
};

#endif

/************************************************************************************
FronTier is a set of libraries that implements differnt types of Front Traking algorithms.
Front Tracking is a numerical method for the solution of partial differential equations 
whose solutions have discontinuities.  


Copyright (C) 1999 by The University at Stony Brook. 
 

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

******************************************************************************/


/*
*				fdecs.h:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains declarations of variables related to the front.  A front
*	consists of an interface (see int.c for documentation) and
*	and the values of the physical state variables on each side of the
*	interface, together with further variables which specify algorithms for
*	the processing of fronts.
*/

#if !defined(_PATRECON_H)
#define _FDECS_H

/* box reconstruction */

enum _TRI_POS {TRIOUT,
	       TRIIN,
	       TRIINSIDE,
	       TRIFLAGA};
typedef enum _TRI_POS TRI_POS;

struct _BBOX {
        double           f[2][3];
        int             side[2][3];     /* yes: the side is constraint */
        struct _BBOX    *prev, *next;
};
typedef struct _BBOX BBOX;

/*recon box */

#define	MAX_RBOX	100
#define	MAX_RBOX_PROC	32
#define RECON_BUFFER	2
#define MAX_DSTORE	100
#define	TABFACTOR	1000

struct _RECON_BOX {
	int	bmin[3], bmax[3];
	double	fmin[3], fmax[3];
	int	flag, proc, np;
	int	number;
	int	procs[MAX_RBOX_PROC], bd_sum[MAX_RBOX_PROC];
	int	shift[MAX_RBOX_PROC][3];
	int	bd_flag[3][2];
};
typedef struct _RECON_BOX RECON_BOX;

struct _GGRID {
	RECON_BOX	rbox;	   /* local */
	int 	gmax[3];	   /* global */
	double	GL[3], GU[3];
	double	h[3];
};
typedef struct _GGRID GGRID;

struct _POINT_BUCKET {
	RECT_GRID	ggr;
	double		btol[3];
	int		***np;	 /* number of points in each bucket */
	POINT		*****pt; /* point pointers in each bucket */
	POINT		**bpt;   /* point store */
};
typedef struct _POINT_BUCKET POINT_BUCKET;

struct _DTABLE{
	int		nbox;
	int		dtab[MAX_RBOX][2];
	int		dtabstore[MAX_DSTORE];
};
typedef struct _DTABLE DTABLE;

#endif /* !defined(_PATRECON_H) */


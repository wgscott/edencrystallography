/************************************************************************

EDEN - Recovery of electron density from X-ray diffraction patterns.
Copyright (C) 2002 Hanna Szoke

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

Hanna Szoke may be reached at szoke2@llnl.gov or at
Lawrence Livermore National Laboratory
7000 East Ave., Livermore, CA 94550-92344. USA.

DISCLAIMER

This software was prepared as an account of work sponsored by an
agency of the United States Government.  Neither the United States Government
nor the University of California nor any of their employees, makes any
warranty, express or implied, or assumes any liability or responsibility for
the accuracy, completeness, or usefulness of any information, apparatus,
product, or process disclosed, or represents that its use would not infringe
privately-owned rights.  Reference herein to any specific commercial products,
process, or service by trade name, trademark, manufacturer, or otherwise, does
not necessarily constitute or imply its endorsement, recommendation, or
favoring by the United States Government or the University of California.  The
views and opinions of authors expressed herein do not necessarily state or
reflect those of the United States Government or the University of California,
and shall not be used for advertising or product endorsement purposes.
*************************************************************************/
/*                                                                      */
/*     Program: ccg.h                                                   */
/*     Date:    2/24/95                                                 */
/*     Author:  Erik M. Johansson, translation of Dennis Goodman's      */
/*              FORTRAN constrained least squares conjugate gradient    */
/*              codes                                                   */
/*                                                                      */
/*     Mods:    none                                                    */
/*                                                                      */
/*     (c) Copyright 1995 the Regents of the University                 */
/*         of California.  All rights reserved.                         */
/*                                                                      */
/*     This software is a result of work performed at Lawrence          */
/*     Livermore National Laboratory.  The United States Government     */
/*     retains certain rights therein.                                  */
/*                                                                      */
/************************************************************************/

extern	real	dfdxpc ;	/* input, frac. change df/dx in getsol */

/* ivec definitions */

#define BOUNDED_NC	0
#define UNCONSTRAINED	1
#define FIXED		2
#define BOUNDED_LC_FREE 3
#define BOUNDED_UC_FREE 4
#define BOUNDED_LC_FIXED 5
#define BOUNDED_UC_FIXED 6

/* istop/iqt definitions */

#define NORMAL_STOP		0
#define FN_LT_FMIN		1
#define MAX_ITER_LIMIT		2
#define MAX_FN_CALLS_EXTRAP	3
#define DX_TOO_SMALL		4
#define MAX_FN_CALLS_INTERP	5
#define AP_EQ_AMX		6
#define	GRADIENT		7
#define	INTERRUPT		8
#define RESSQ_LT_DISCRP		-1
#define TOLERANCE_PROBLEMS	-2
#define INITIAL_AP_LT_0		-3
#define INITIAL_DD_GT_0		-4

#define MAXIT	600		/* max. # iterations in getsol */
#define TOL	1.e-6		/* rel. accuracy of function used by getsol */
#define COSMIN 0.01
#define GAMMA 0.1
#define CCG_RHO   0.0001

				/* 12/22/97 - Care with parentheses! */

#define SIGN(a, b) (((b) >= (0)) ? (fabs(a)) : (-fabs(a)))
#ifndef MIN
#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#endif

#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif


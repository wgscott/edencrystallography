/******************************************************************************

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

                                REMAP.H

  Title:        Include file for code involving close vicinities in Np space.
  Author:       Hanna Szoke
  Date:         3/12/97
  Function:     This file contains data structures used in the 
		generation of blobs (ran.c); in calculating the
		low-resolution surface cost function or Sayre's 
		equation cost function (qshare.c and cost.c);
		and for Tetra and NCS (remap.c and ncsutil.c).  

                This file was modified from an earlier "blobs.h".

******************************************************************************/

#define	CUTOFF	3.e-3	/* exp. radius beyond which Eden ignores contributions */
#define	IRANGE	5	/* Distance in index units to start looking 
			for gexp entries.   */
#define	MARGIN	1	/* extra margin of safety in "small cell" dimensions */
#define	NGMAX	400	/* max # elements in gexp */
#define DINT	8	/* # of subdivisions of a grid cell for interpolation */

extern	int	Nb ;		/* Number of vicinity points */
extern	real	*gexp ;		/* array of vicinity exponential factors */

extern	int	*poff ;		/* Nptotal array of pointers to indexlist  */
extern	int	*indexlist ;	/* array for fast indexing to vicinity pts. */


/*******************************************************************************

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

                                REG.H

  Title:        Include file for postprocessors 
		(regrid.c, shapes.c and count.c)

  Author:       Hanna Szoke

  Date:         Separated from Regrid: Sept, 2001

*******************************************************************************/

typedef struct {
   int	beg, end ;		/* start and end indices */
} LIMITS ;

typedef struct {
   real	beg, end ;	/* start and end fractional coords */
} FLIMITS ;

extern	FLIMITS	flimx, flimy, flimz ;	/* input fractional xyz limits	*/
extern	FLIMITS	flime ;
extern	LIMITS	rlimx, rlimy, rlimz ;	/* coord limits after regrid */
extern	real	*rnump ;		/* electron/voxel initial array */
extern	real	*rho ;			/* for holding electron density */
extern	real	*dout ;			/* rho's that are written out	*/
extern	int	spread_fac ;		/* regridding factor */
extern	char	erange ;		/* T/F ? - apply range (flime) */
extern	int	Nxout, Nyout, Nzout ;
extern	int	Npout ;

extern	real	*drdx ;			/* 3 1st derivatives... */
extern	real	*drdy ;
extern	real	*drdz ;

extern	real	*d2rdxdx ;		/* ... and 6 2nd derivatives. */
extern	real	*d2rdxdy ;
extern	real	*d2rdxdz ;
extern	real	*d2rdydy ;			
extern	real	*d2rdydz ;
extern	real	*d2rdzdz ;

extern	char	vox_filename[] ;
extern	char	list_filename[] ;
extern	char	pdb_filename[] ;


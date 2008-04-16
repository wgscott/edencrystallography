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

				PDB.H

  Title:        Include info re. pdb file
  Author:       Hanna Szoke
  Date:         Feb. 1, 1995 (separated)
  Function:     Contains info needed by all programs using pdb files.

******************************************************************************/

#define MAXNATOMS      500000    /* max # atoms in unit cell      */

#define PROTEIN_DENSITY 0.42    /* changed from 0.48, 2/26/2001	 */
				/* In units of el/cubic Angstrom */
#define SOLVENT_DENSITY 0.34    /* "  "     "  "        "        */

        /*      Input and derived parameters.   */

extern	int     Nin ;                   /* number of atoms in input */
extern	int     Natoms ;                /* number of atoms in unit cell */
extern	int	atno[MAXNATOMS] ;	/* atom number */
extern	char	atid[MAXNATOMS][5] ;	/* atom id */
extern	char	fullatid[MAXNATOMS][6] ;	/* full atom id */
extern	char	groupid[MAXNATOMS][4] ;	/* amino acid id */
extern	int	groupno[MAXNATOMS] ;	/* residue number */
extern	POINT	pos[MAXNATOMS] ;	/* ortho pos of atoms in unit cell */
extern	float	occup[MAXNATOMS] ;	/* occupancies */
extern	real  Bj[MAXNATOMS] ;         /* Debye-Waller factors */
extern	real  fj[MAXNATOMS] ;         /* number of electrons per input atom */
extern	int  Zj[MAXNATOMS] ;         /* atomic number per input atom */
extern	int	marker[MAXNATOMS] ;	/* marker of overlap atoms */


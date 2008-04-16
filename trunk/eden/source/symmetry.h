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

                                SYMMETRY.H

  Title:  Include file for crystallographic symmetry group information. 
  Author: Hanna Szoke
  Date:	  Oct. 1994; rewritten Oct. 1996 for all space groups.

******************************************************************************/


#define	MEL	12			/* # of matrix elements in 1 matop */

		/**************************
		Commonly used functions
		**************************/

void	apply_symop(real, real, real, real *, 
		    real *, real *, real *) ;
void	apply_symop_hkl(int, int, int, real *, int *, int *, int *, real *) ;
void	ortho_to_fract(real, real, real, real *, real *, real *) ;
void	Tr_ortho_to_fract(real, real, real, real *, real *, real *) ;
void	fract_to_ortho(real, real, real, real *, real *, real *) ;
void	index_to_fract(int, int, int, int, real *, real *, real *) ; 
int	fract_to_index(int, real, real, real, int, int, int) ;

extern	char	SGname[] ;		/* name of crystal symmetry group */
extern	char	PGname[] ;		/* name of crystal point group */

extern	char	crystal_system[20] ;	/* 1 of the 7: 
					TRICLINIC, MONOCLINIC, ORTHORHOMBIC, 
					TETRAGONAL, TRIGONAL, HEXAGONAL or 
					CUBIC	*/

extern	int	laue_group ;		/* one of 11: a subdivision of 
					crystal_system used to identify 
					unique data in (h,k,l) */

extern	real	*matop ;		/* to hold matrix operators for doing
					transformations between points related 
					by crystal symmetry */

extern	int	Nprim ;			/* # of primitive RT matrices */

extern	int	maxdenom[3] ;		/* max. denominators of translation 
					terms, needed for ensuring good 
					Nx, Ny and Nz */

extern	int	Ncs ;			/* # of asymmetric units in crystal 
					(= # rotation/translation operations
					for the selected space group)  */

extern	int	Ics ;			/* # of unique grid points in an 
					asymmetric unit */

extern	char	detwin ;		/* is detwinning enabled? T/F */
extern	char	t_type ;		/* A (amplitude) or I (intensity)*/
extern	int	*t_matrix ;		/* twinning matrix transformations */
extern	real	t_frac;			/* twinning fraction */

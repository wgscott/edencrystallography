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

                                DIMS.H

  Title:  Include file for Electron Density Calculation by Holographic 
	  Reconstruction: basic dimensions and related variables.
  Author:  Hanna Szoke
  Date:  Separated off on 8/16/95

******************************************************************************/

real	nusq(int, int, int) ;		
void	fetch_hkl(int, int *, int *, int *) ;

#define	NUMSEGS	5		/* for analyzing shells in (hkl) space) */
#define	SIMPLE		0	/* simple cubic grid type */
#define	BODY_CENTERED	1	/* body-centered cubic grid type */

#define	MIN_DENSITY	0.0	/* "  "     "  "        "        */
#define	MAX_DENSITY	1.e10  	/* "  "     "  "        "        */
 
extern	int     Nx, Ny, Nz ;    /* solution space dimensions */
extern	float	floatNx ;
extern	float	floatNy ;
extern	float	floatNz ;
extern	int     Np ;            /* product of Nx*Ny*Nz */
extern	int     Nptotal ;       /* Np (sc) or 2*Np (bcc) */
extern	int	Nhr ;		/* # highres points	*/
extern	int	Npextended ;	/* Nptotal + Nhr (high resolution)	*/
extern	int	grid_type ;	/* solution grid type (sc or sc+bcc)	*/
extern	int	Nh, Nk, Nl ;	/* reciprocal space dimensions	*/
extern	int     Nhkl ;          /* number of reciprocal lattice points */

extern	real  grid_spacing ;  /* solution resolution in Angstrom */
extern	real	dr ;		/* actual average resolution in Angstrom */
extern	real  eta ;           /* factor describing spread of Gaussians */
extern	real	delfac ;	/* PISQ * eta * dr * dr	*/
extern	real	input_res;	/* resolution input variable	*/
extern	real	low_res_cutoff;	/* low resolution cutoff	*/
extern	real	limit_rat ;	/* 1/(dr*dr)	*/
extern	real	dsq_limit ;	/* effective limit from data or limit_rat */
extern	real	nusqlim[] ;	/* values for 1/8, 1/4,... of (hkl) shells */
extern	int	useSig ;	/* if TRUE, use sigmas in input fobs files */
extern	int	anom_flag ;	/* if TRUE, input fobs/fcalc are anomalous */
extern	int	HighRes ;	/* if TRUE, high-resolution processing on  */
extern	real	hrCutoff ;	/* hrCutoff * <vox> = high res limit */
extern	real	fscale ;	/* scaling factor for fo data */

extern	float	solvent_voxel ;	/* ave. voxel value for the solvent */
extern	float	protein_voxel ;	/* ave. voxel value for the protein */
extern	float	min_voxel ;	/* min. voxel value for the protein */
extern	float	max_voxel ;	/* max. voxel value for the protein */


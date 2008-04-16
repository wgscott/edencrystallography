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

                                MIR.H

  Title:  MIR include file for Electron Density Calculation 
	  by Holographic Reconstruction
  Author:  Hanna Szoke
  Date:  1/27/94

******************************************************************************/

#define	MAXDER		8	/* max # of derivative (isomorphous replacement)
				   structure factor files  */
#define	MAXRES		4	/* max # of resolution sets among data sets */


extern	real	abs_scale_m[] ;
extern	real	relwt_der[] ;   /* relative weights for derivatives  */
extern	real	res_der[] ;   	 /* intrinsic resolution for derivatives  */
extern	int	useset[] ;	 /* which resolution should be used? */
extern	real	*rexpfac[] ;	 /* pointer to expfac array for data sets with
				    resolution r*/
extern	char    *intermask ;    /* intersection of all Fobs masks (MIR or MAD) */
extern	COMPLEX	*rexpfac_bcc[] ;	 /* ditto for expfac_bcc array */
extern	int	Nder;		 /* number of derivatives */
extern	int	NM;		 /* Nder+1 (Native plus MIR files)  */
extern	int	Nres ;		 /* number of resolution sets for expfac */
extern	real	resrat[] ; 	 /* resolution ratios for the resolution sets */

extern	char	fc_heavy_fn[][MAXSTRING] ; /* filenames of MIR 
					heavy atom str factors */
extern	char	fo_der_fn[][MAXSTRING] ; /* filenames of MIR obs str factors */

extern	real	fscale_der[] ;	/* input scale factors for derivative fobs */
extern	int	autoscale ;	/* if TRUE (default) apply rel_scale_m's */

extern	real	fprime[] ;	/* f' for anomalous MAD heavy file */
extern	real	f2prime[] ;	/* f'' for anomalous MAD heavy file */
extern	real	Zheavy[] ;	/* Z for heavy atom - only one per file! */
extern	char	anom_der_flag[] ;	/* flag identifying anomalous data */

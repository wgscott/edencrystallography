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

                                EDEN.H

  Title:  Include file for Electron Density Calculation by Holographic 
	  Reconstruction: main set of data structures, used by Back, 
	  Forth and Solve.
  Author:  Hanna Szoke
  Date:  1/10/92

******************************************************************************/

#include "util.h"		/* general definitions */
#include "cellparams.h"		/* a, b, c, angles, crystal_symmetry */
#include "symmetry.h"		/* symmetry groups supported by Eden */
#include "dims.h"		/* grid_type, Np, Nhkl, and other dimensions */

#define MAXCONSTR	12	/* max. no. of constraints for a problem */
#define DISCRP	1.e-4		/* discrepancy parameter used in 'funct' */
#define  EPS_PHI 1.e-10		/* cutoff angle for gradient exception 
				   processing 	*/

enum ctype {CS, NCS, PHASE_EXT, SAYRE, TARGET, SINGLET, TRIPLET} ; 

        /*      Input and derived parameters.   */

extern	float	R_stop ;	/* R_factor at which code should stop	*/
extern  real	discrp_frac ;   /* fraction of DISCRP set by user (def: 1) */
extern	int	robust ;	/* use robust hkl cost function? (T/F) */
extern	real	rob_factor ;	/* factor to use with robust cost function */

extern	char	fc_filename[] ;/* name of file, calculated structure factors */
extern	char	fo_filename[] ;/* name of file, observed structure factors */
extern	char	md_filename[] ;/* name of file, model electron map */

extern	int	Nconstraints ;	/* number of constraints */
extern	int	Ntargets ;	/* number of target constraints */
extern	int	con_type[] ;	/* target type: 1 of the ctype (see above) */
extern	real	relwt_con[] ;	/* relative weight for a constraint */
extern	real	cost_addend[] ;	/* constant factor in the cost function */
extern	real	con_coef[] ;	/* coefficient for cost function calc. */
extern	char	target_type[][MAXSTRING] ;/* target type, target constraint */
extern	char	ta_filename[][MAXSTRING] ;/* target file, target constraint */
extern	char	wt_filename[][MAXSTRING] ;/* weight file, target constraint */

extern	real	phase_ext_res ;		/* low resolution for phase ext */

	/* Widely-used arrays in (hkl) space.  */

extern	char	*maskfo ;	/* mask for fobs	*/
extern	char	*maskfc ;	/* mask for fcalc	*/
extern	real	*hkl_weight ;	/* weights for (hkl) in half-ellipsoid */
extern	real	*sigma ;	/* for fobs sigma values */
extern	real	*expfac ;       /* exponential factors for sc solutions */
extern	COMPLEX	*expfac_bcc ;   /* exponential factors for bcc solutions */
extern	COMPLEX	*twinned_fc ;	/* array for twinning manipulations 	*/

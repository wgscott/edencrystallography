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

                                EXPANDFC

  Title:        EXPANDFC - Utility to expand an fcalc file to P1 symmetry.
  Author:       Hanna Szoke
  Date:         8/9/95
  Function:     Expand an (hkl) file to the full half-ellipsoid
		Allows "ANOM TRUE" for anomalous dispersion.

*******************************************************************************/

#include "util.h"		/* general definitions */
#include "cellparams.h"		/* a, b, c, angles, crystal_symmetry */
#include "symmetry.h"		/* symmetry groups information */
#include "dims.h"		/* Np, Nhkl, and other dimensions */

void	expandfc_main(int argc, char *argv[])
{
     int	readfcalc1_allh(char *, COMPLEX *, char *) ;
     int	readfcalc2_allh(char *, COMPLEX *, COMPLEX *, char *, char *) ;
     void	writefcalc(char *, COMPLEX *, char *) ;

     COMPLEX	*R1 ;
     COMPLEX	*R2=NULL ;
     char	*maskfc1 ;	
     char	*maskfc2=NULL ;	
     int	triclinic = FALSE ;
     int	N ;

     char	outfilename1[MAXSTRING] ;
     char	outfilename2[MAXSTRING] ;


     /************************************************************
     Check for right ballpark, say hello, identify file names 
     ************************************************************/

     if (argc > optind+3) 	/* too many arguments */
	  ballpark(caller) ; 
     else if (argc == optind+3)
          sprintf(sf_filename, "%s", argv[optind+2]) ;
     else {
	  prompt("\nPlease enter the name of an fcalc file: ") ;
          sscanf(terminp, "%s", sf_filename) ;
	  sprintf(message, "Running %s on fcalc file %s", caller, sf_filename) ;
	  printTwice(message) ;
     }

     hello(caller) ; 

     strip_path(sf_filename, pwd_name) ;

     /************************************************************
     Get input parameters, check whether user has anomalous data.
     and allocate space.                 

     Note outfilename2 is used for anomalous dispersion only if 
     input crystal class is triclinic.
     ************************************************************/

     readBasicInput() ;

     if (anom_flag == FALSE) {
          extend_filename(pwd_name, outfilename1, "_P1") ;
     }
     else {
          extend_filename(pwd_name, outfilename1, "_P1plus") ;
          extend_filename(pwd_name, outfilename2, "_P1minus") ;
     }

     if (laue_group == 1)
	  triclinic = TRUE ;

     R1   = (COMPLEX *) e_malloc(Nhkl*sizeof(COMPLEX), caller) ;
     maskfc1 = (char *) e_malloc(Nhkl*sizeof(char), caller) ;
     
     if ((anom_flag) && (triclinic)) {
          R2   = (COMPLEX *) e_malloc(Nhkl*sizeof(COMPLEX), caller) ;
          maskfc2 = (char *) e_malloc(Nhkl*sizeof(char), caller) ;
     }

     /************************************************************
     Read, expand to P1 (laue_group = 1) and write input file  
     ************************************************************/

     sprintf(message, "Reading file %s...\n", sf_filename) ;
     printTwice(message) ;

     if (anom_flag && triclinic) {
          N = readfcalc2_allh(sf_filename, R1, R2, maskfc1, maskfc2) ;
          sprintf(message, "Writing expanded data (%d lines) to %s and %s", 
	  N, outfilename1, outfilename2) ;
          printTwice(message) ;
          writefcalc(outfilename1, R1, maskfc1);
          writefcalc(outfilename2, R2, maskfc2);
     }
     else {
          N = readfcalc1_allh(sf_filename, R1, maskfc1) ;
          sprintf(message, "Writing expanded data (%d lines) to %s", 
	  N, outfilename1) ;
          printTwice(message) ;
          writefcalc(outfilename1, R1, maskfc1);
     }

     return ;
}

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

=======
EDEN - Recovery of electron density from X-ray diffraction patterns.
Copyright (C) 2002 Hanna Szoke

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,

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
P


>>>>>>> 1.2
                                EXPANDFO

  Title:        EXPANDFO - Utility to expand an fobs file to P1 symmetry.
  Author:       Hanna Szoke
  Date:         8/9/95
  Function:     Expand an (hkl) file to the full half-ellipsoid
		Allows 'anom true' for anomalous dispersion.

*******************************************************************************/

#include "util.h"		/* general definitions */
#include "symmetry.h"		/* symmetry group information */
#include "dims.h"		/* Np, Nhkl, and other dimensions */

void	expandfo_main(int argc, char *argv[])
{
     void	writefobs(char *, real *, char *, real *) ;
     int	readfobs0_allh(char *, real *, char *, real *) ;
     int	readfobs1_allh(char *, real *, char *, real *, char *) ;
     int	readfobs2_allh(char *, real *, real *, char *, char *,
			      real *, real *) ;
     int	checkSigmas(real *, char *) ;     /* in hklinit.c */

     real	*F1 ;
     real	*F2=NULL ;
     real	*sigma1 ;
     real	*sigma2=NULL ;
     char	*maskfo1 ;	
     char	*maskfo2=NULL ;	
     char	*mcentric=NULL ;
     int	triclinic = FALSE ;
     int	Nread = 0 ;

     char	outfilename1[MAXSTRING] ;
     char	outfilename2[MAXSTRING] ;


     /***********************************************************
     Check for right ballpark, say hello, identify file names.
     Note that outfilename2 is used only for anomalous dispersion 
     if input crystal class is triclinic.
     ************************************************************/

     if (argc > optind+3) 	/* too many arguments */
	  ballpark(caller) ; 
     else if (argc == optind+3)
          sprintf(sf_filename, "%s", argv[optind+2]) ;
     else {
	  prompt("\nPlease enter the name of an fobs file: ") ;
          sscanf(terminp, "%s", sf_filename) ;
	  sprintf(message, "Running %s on fobs file %s", caller, sf_filename) ;
	  printTwice(message) ;
     }

     hello(caller) ; 

     strip_path(sf_filename, pwd_name) ;

     /************************************************************
     Get input parameters and allocate space.                
     See whether user has anomalous data.  Default is FALSE.
     ************************************************************/

     readBasicInput() ;

     if (anom_flag == FALSE) 
     {
          extend_filename(pwd_name, outfilename1, "_P1") ;
     }
     else 
     {
          extend_filename(pwd_name, outfilename1, "_P1plus") ;
          extend_filename(pwd_name, outfilename2, "_P1minus") ;
     }

     if (laue_group == 1)
	  triclinic = TRUE ;

     F1     = (real *) e_malloc(Nhkl*sizeof(real), caller) ;
     sigma1 = (real *) e_malloc(Nhkl*sizeof(real), caller) ;
     maskfo1  = (char *) e_malloc(Nhkl*sizeof(char), caller) ;
     
     if (anom_flag) {
          if (triclinic) {
               F2     = (real *) e_malloc(Nhkl*sizeof(real), caller) ;
               sigma2 = (real *) e_malloc(Nhkl*sizeof(real), caller) ;
               maskfo2  = (char *) e_malloc(Nhkl*sizeof(char), caller) ;
          }
	  else {
	       mcentric = (char *) e_malloc(Nhkl*sizeof(char), caller) ;
          }
     }


     /************************************************************
     Read, expand to P1 and write input file.
     ************************************************************/

     sprintf(message, "Reading file %s ...\n", sf_filename) ;
     printTwice(message) ;

     if (anom_flag && triclinic) {
          Nread = readfobs2_allh(sf_filename, F1, F2, maskfo1, maskfo2,
	       			sigma1, sigma2) ;
          useSig = checkSigmas(sigma1, maskfo1) ;
          useSig = checkSigmas(sigma2, maskfo2) ;
	  sprintf(message,"Writing expanded data to %s and %s", 
	  outfilename1, outfilename2) ;
          printTwice(message) ;
          writefobs(outfilename1, F1, maskfo1, sigma1) ;  
          writefobs(outfilename2, F2, maskfo2, sigma2) ;  
     }
     else {
	  if (anom_flag)
               Nread = readfobs1_allh(sf_filename, F1, maskfo1, sigma1, 
					mcentric) ;
	  else
               Nread = readfobs0_allh(sf_filename, F1, maskfo1, sigma1) ;

	  sprintf(message,"Writing expanded data to %s", outfilename1) ;
          printTwice(message) ;
          writefobs(outfilename1, F1, maskfo1, sigma1) ;  
     }

     return ;
}

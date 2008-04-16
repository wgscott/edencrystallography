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

                                CADHKL

  Title:        CADHKL - Utility for adding,subtracting or merging (hkl) files
  Author:       Hanna Szoke
  Date:         10/23/97 (in present form)
  Function:     Add (or subtract) comparable (hkl) structure factors from 
		two input files and place the result in an output file or
		merge 

*******************************************************************************/

#include "util.h"
#include "cellparams.h"
#include "symmetry.h"
#include "dims.h"

static	real	C1, C2 ;	/* coefficients, defaulting to 1. */
static	void	combineFcalcFobsFactors(char *, char *, 
		COMPLEX *, real *, COMPLEX *) ;
static	void	combineFobsFobsFactors(char *, char *) ;
static	void	combineFcalcFactors(char *, char *, 
		COMPLEX *, COMPLEX *, COMPLEX *) ;
static	void	combineFcalcFcalcPhases(char *, char *,
		COMPLEX *, COMPLEX *, COMPLEX *) ;


void	cadhkl_main(int argc, char *argv[])
{
     int      	readfcalc(char *, COMPLEX *, char *) ;
     int      	readfobs(char *, real *, char *, real *) ;
     void	fetch_coefs(real *, real *);
     void	writefcalc(char *, COMPLEX *, char *) ;
     void	writefcalc_non0(char *, COMPLEX *, char *) ;
     void	writefobs(char *, real *, char *, real *) ;

     char	filename1[MAXSTRING] ;
     char	filename2[MAXSTRING] ;
     char	filename3[MAXSTRING] ;
     COMPLEX	*R1=NULL, *R2=NULL, *R3=NULL ;
     real	*F1=NULL, *sig1=NULL ;
     real	*F2=NULL, *sig2=NULL ;
     char	*mask1, *mask2 ;
     char	merge = FALSE ;
     char	elim = FALSE ;
     char	copy = FALSE ;

     /* Check for right ballpark */

     if (argc > optind+5) 	/* too many arguments */
	  ballpark(caller) ; 
     else if (argc == optind+5) {
          sprintf(filename1, "%s", argv[optind+2]) ;
          sprintf(filename2, "%s", argv[optind+3]) ;
          sprintf(filename3, "%s", argv[optind+4]) ;
     }
     else {
	  prompt("\nPlease enter 3 file names - 2 input, 1 output: ") ;
          sscanf(terminp, "%s %s %s", filename1, filename2, filename3) ;
	  sprintf(message, "Running %s on files %s and %s, output to %s", 
	          caller, filename1, filename2, filename3) ;
	  printTwice(message) ;
     }

     hello(caller) ;
     readBasicInput() ;

     if ((strncmp(mode_info, "none", 4) == 0) ||
         (strncmp(mode_info, "add", 3) == 0))  
	 printTwice("\nThis is a run in add mode.\n") ;
     
     else if (strncmp(mode_info, "merg", 4) == 0) {
	 merge = TRUE ;    
	 printTwice("\nThis is a run in merge mode.\n") ;
     }	
     else if (strncmp(mode_info, "elim", 4) == 0) {
	 elim = TRUE ;    
	 printTwice("\nThis is a run in eliminate mode.\n") ;
     }	
     else if (strncmp(mode_info, "copy", 4) == 0) {
	 copy = TRUE ;    
	 printTwice("\nThis is a run in copy (phases) mode.\n") ;
     }	
     else
	  EdenError("Illegal mode_info value!") ;

     fetch_coefs(&C1, &C2) ;

     if (merge) {
         R1 = (COMPLEX *) e_malloc(Nhkl*sizeof(COMPLEX), caller) ;
         F2 =  (real *) e_malloc(Nhkl*sizeof(real),  caller) ;
         sig2 = (real *) e_malloc(Nhkl*sizeof(real),  caller) ;
         R3 = (COMPLEX *) e_malloc(Nhkl*sizeof(COMPLEX), caller) ;
     }
     else if (elim) {
         F1 =  (real *) e_malloc(Nhkl*sizeof(real), caller) ;
         sig1 = (real *) e_malloc(Nhkl*sizeof(real),  caller) ;
         F2 =  (real *) e_malloc(Nhkl*sizeof(real),  caller) ;
         sig2 = (real *) e_malloc(Nhkl*sizeof(real),  caller) ;
     }
     else {
         R1 = (COMPLEX *) e_malloc(Nhkl*sizeof(COMPLEX), caller) ;
         R2 = (COMPLEX *) e_malloc(Nhkl*sizeof(COMPLEX), caller) ;
         R3 = (COMPLEX *) e_malloc(Nhkl*sizeof(COMPLEX), caller) ;
     }

     mask1 = (char *) e_malloc(Nhkl*sizeof(char), caller) ;
     mask2 = (char *) e_malloc(Nhkl*sizeof(char), caller) ;

     if (merge) {
          readfcalc(filename1, R1, mask1) ;
	  readfobs(filename2, F2, mask2, sig2) ;
          combineFcalcFobsFactors(mask1, mask2, R1, F2, R3) ;
          writefcalc(filename3, R3, mask1) ;
     }
     else if (copy) {
          readfcalc(filename1, R1, mask1) ;
	  readfcalc(filename2, R2, mask2) ;
          combineFcalcFcalcPhases(mask1, mask2, R1, R2, R3) ; 
	  writefcalc(filename3, R3, mask1) ;
     }
     else if (elim) {
          readfobs(filename1, F1, mask1, sig1) ;
	  readfobs(filename2, F2, mask2, sig2) ;
          combineFobsFobsFactors(mask1, mask2) ;
          writefobs(filename3, F1, mask1, sig1) ;
     }
     else {
          readfcalc(filename1, R1, mask1) ;
          readfcalc(filename2, R2, mask2) ;
          combineFcalcFactors(mask1, mask2, R1, R2, R3) ;
          writefcalc_non0(filename3, R3, mask1) ;
     }
}

void	combineFcalcFactors(char *mask1, char *mask2, 
			    COMPLEX *R1, COMPLEX *R2, COMPLEX *R3)
{
     int        n ;
     int	sum ;

     for (n = 0; n < Nhkl; n++)
	  *(mask1+n) = *(mask1+n) * *(mask2+n) ;

     for (n = 0, sum = 0; n < Nhkl; n++)
	  sum += *(mask1+n) ;
     sprintf(message, " Number of entries in the output is %d", sum) ;
     printTwice(message) ;

     for (n = 0; n < Nhkl; n++) {
	if (*(mask1+n)) {	/* Both input entries must be present */
          (R3+n)->re = C1 * (R1+n)->re + C2 * (R2+n)->re ;
          (R3+n)->im = C1 * (R1+n)->im + C2 * (R2+n)->im ;
	}
        else {
          (R3+n)->re = 0 ;
          (R3+n)->im = 0 ;
	}
     }
}

void	combineFcalcFobsFactors(char *mask1, char *mask2, 
			    COMPLEX *R1, real *F2, COMPLEX *R3)
{
     int        n ;
     real	phase ;
     real	amp0, amplitude ;

     for (n = 0; n < Nhkl; n++)
	  *(mask1+n) = *(mask1+n) * *(mask2+n) ;

     amp0 = 1. ;	/* we don't use R1->re for scaling amplitudes! */

     for (n = 0; n < Nhkl; n++) {

	/* In order to combine entries meaningfully,
	(1) The F2 amplitude must be > 0;
	(2) The R1 phase must be well-defined; and
	(3) (of course) both input entries must be present */

        if (*(F2+n) == 0)
	    *(mask1+n) = 0 ;

        if (((R1+n)->re == 0) && ((R1+n)->im == 0)) 
	    *(mask1+n) = 0 ;

	if (*(mask1+n)) {	

	    getAmpPhase(R1+n, amp0, &amplitude, &phase) ;
               
            (R3+n)->re = *(F2+n) * cos(phase*DTOR) ;
            (R3+n)->im = *(F2+n) * sin(phase*DTOR) ;

	}
        else {

            (R3+n)->re = 0 ;
            (R3+n)->im = 0 ;

	}
     }
}

void	combineFobsFobsFactors(char *mask1, char *mask2) 
{
     int        n ;

     for (n = 1; n < Nhkl; n++)
	  *(mask1+n) = *(mask1+n) * (!*(mask2+n)) ;

}
void	combineFcalcFcalcPhases(char *mask1, char *mask2, 
			    COMPLEX *R1, COMPLEX *R2, COMPLEX *R3)
{
     int        n ;
     real	phase1, phase2 ;
     real	amp0, amplitude1, amplitude2 ;

     for (n = 0; n < Nhkl; n++) {
	  (R3+n)->re = (R1+n)->re ;
	  (R3+n)->im = (R1+n)->im ;
     }

     amp0 = R1->re ;

     for (n = 1; n < Nhkl; n++) 
	 if (*(mask1+n) && *(mask2+n)) {
	    getAmpPhase(R1+n, amp0, &amplitude1, &phase1) ;
	    getAmpPhase(R2+n, amp0, &amplitude2, &phase2) ;
               
            (R3+n)->re = amplitude1 * cos(phase2*DTOR) ;
            (R3+n)->im = amplitude1 * sin(phase2*DTOR) ;
	}
}

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

                                DPHASE

  Title:        DPHASE - Utility to compare phases in 2 model (hkl) files
  Author:       Hanna Szoke
  Date:         1/10/94
  Function:     This program calculates the phase differences between
		comparable (hkl) structure factors in two input files.
		The phase differences, weighted by the amplitudes, are 
		calculated twice - once based on the first input file and 
		once based on the 2nd.  They are then summed and reported
		over shells in (hkl) space.

		User must supply an Eden-style input file (containing
		the unit cell dimensions, symmetry group and a data 
		resolution) and two (hkl) file names.  

  		Code excludes terms for which the amplitude in either file 
		is 0 and it excludes the (000) term.
*******************************************************************************/

#include "util.h"		/* general definitions */
#include "cellparams.h"		/* a, b, c, angles, crystal_symmetry */
#include "symmetry.h"		/* symmetry group information */
#include "dims.h"		/* Np, Nhkl, and other dimensions */

static	void   	useAP(COMPLEX *, real *, real *) ;
static	void	prepareShells(char *, int *) ;
static	void	comparePhases(real *, real *, real *, int *, char *) ;


void	dphase_main(int argc, char *argv[])
{
     int	readfcalc(char *, COMPLEX *, char *) ;
     real	get1Rfac(char *, real *, COMPLEX *, int) ;
     void	shellHeader_R() ;
     void	set_centrics(char *) ;
     void	Cmaskit(COMPLEX *, char *, int) ;

     COMPLEX	*Rcomp1, *Rcomp2 ;
     char	*maskfc ;	
     char	*mcentric ;
     char	filename1[MAXSTRING] ;
     char	filename2[MAXSTRING] ;
     real	*amp1, *amp2 ;
     real	*pha1, *pha2 ;
     int	*shell ;
     char	*mask ;
     int	n ;

     /************************************************************
     No more switches - -h handled by eden.c
     Check for right ballpark, say hello, identify file names 
     ************************************************************/

     if (argc > optind+4) 	/* too many arguments */
	  ballpark(caller) ; 
     else if (argc == optind+4) {
          sprintf(filename1, "%s", argv[optind+2]) ;
          sprintf(filename2, "%s", argv[optind+3]) ;
     } 
     else {
	  prompt("\nPlease enter 2 filenames: ") ;
          sscanf(terminp, "%s %s", filename1, filename2) ;
	  sprintf(message, "Running %s on files %s and %s", 
		  caller, filename1, filename2) ; 
	  printTwice(message) ;
     }

     hello(caller) ;

     /************************************************************
     Check that both files are available.                     
     ************************************************************/

     check_exist(filename1) ; 
     check_exist(filename2) ;

     /************************************************************
     Get input parameters and allocate space.                 
     ************************************************************/

     readBasicInput() ;

     strcpy(message, "dphase") ;

     amp1 = (real *) e_malloc(Nhkl*sizeof(real),  message) ;
     amp2 = (real *) e_malloc(Nhkl*sizeof(real),  message) ;
     pha1 = (real *) e_malloc(Nhkl*sizeof(real),  message) ;
     pha2 = (real *) e_malloc(Nhkl*sizeof(real),  message) ;
     
     Rcomp1 = (COMPLEX *) e_malloc(Nhkl*sizeof(COMPLEX), message) ;
     Rcomp2 = (COMPLEX *) e_malloc(Nhkl*sizeof(COMPLEX), message) ;
     shell = (int *)     e_malloc(Nhkl*sizeof(int),  message) ;
     
     mask =     (char *) e_malloc(Nhkl*sizeof(char), message) ;
     mcentric = (char *) e_malloc(Nhkl*sizeof(char), message) ;
     maskfc =   (char *) e_malloc(Nhkl*sizeof(char), message) ;
     
     for (n = 0; n < Nhkl; n++)
	  *(mask+n) = 1 ;

     /************************************************************
     Read both input files, determine how to group them by shell. 
     ************************************************************/

     sprintf(message, "Reading %s: ", filename1) ;
     printTwice(message) ;
     readfcalc(filename1, Rcomp1, maskfc) ;

     for (n = 0; n < Nhkl; n++)
	  *(mask+n) *= *(maskfc+n) ;

     sprintf(message, "%s: ", filename2) ; 
     printTwice(message) ; 
     readfcalc(filename2, Rcomp2, maskfc) ; 

     for (n = 0; n < Nhkl; n++) 
	  *(mask+n) *= *(maskfc+n) ; 

     e_free(Nhkl*sizeof(char), maskfc) ;

     Cmaskit(Rcomp1, mask, Nhkl) ;
     Cmaskit(Rcomp2, mask, Nhkl) ;

     useAP(Rcomp1, amp1, pha1) ;
     useAP(Rcomp2, amp2, pha2) ; 

     prepareShells(mask, shell) ;
     set_centrics(mcentric) ;

     /************************************************************
     Do dphase calculation using each (hkl) in turn as basis 
     ************************************************************/

     printTwice(" ") ;
     printTwice( 
     "In the following reports, the average phase difference (dph)") ;
     printTwice(
     "and the average cosine of the phase difference (dcos), both") ;
     printTwice(
     "weighted by amplitudes, are reported over resolution ranges.") ;
     printTwice( 
     "The information is reported twice - first for all phases,") ;
     printTwice(
     "then for restricted (centric) phases only.  Delta phases ") ;
     printTwice(
     "are followed by R factor calculations for comparing amplitudes.") ;

     sprintf(message,"\n\nComparing phases weighted by %s ...", filename1) ;
     printTwice(message) ;
     comparePhases(amp1, pha1, pha2, shell, mcentric) ;

     sprintf(message,"\n\nComparing phases weighted by %s ...", filename2) ;
     printTwice(message) ;
     comparePhases(amp2, pha2, pha1, shell, mcentric) ;

     sprintf(message,
     "\n\nComparing amplitudes weighted (1) by %s, (2) by %s:",
     filename1, filename2) ;
     printTwice(message) ;
     shellHeader_R() ;
     get1Rfac("    (1) ", amp1, Rcomp2, Nhkl) ;
     get1Rfac("    (2) ", amp2, Rcomp1, Nhkl) ;
}

void	useAP(COMPLEX *R, real *amp, real *pha)
{
     int     	n ;
     real	amp0 ;

     amp0 = R->re ;

     for (n = 0; n < Nhkl; n++) 
	  getAmpPhase(R+n, amp0, amp+n, pha+n) ;

}

void	prepareShells(char *mask, int *shell)
{
     void	setDsqLimit(char *, int) ;
     void	hklTabHead() ;
     int	m, n ;
     int	h, k, l ;

     /************************************************************
     Prepare shell indices for all n's for which there were 
     structure factors in both input files, i.e.,  *(mask+n) == 1.     
     WE EXCLUDE THE (000) TERM (m = 0).                   
     ************************************************************/

     setDsqLimit(mask, Nhkl) ;
     hklTabHead() ;

     for (n = 0; n < Nhkl; n++) 
          *(shell+n) = -1 ;

     for (n = 1; n < Nhkl; n++) {
	 if (*(mask+n)) {	

             fetch_hkl(n, &h, &k, &l) ;

             for (m = 1; m < NUMSEGS; m++) 
	          if (nusq(h, k, l) <= nusqlim[m])
		       break ;

	     if (m < NUMSEGS)	/* don't let m spill over array lengths */
                  *(shell+n) = m ;
        }
     }
}

void	comparePhases(real *amp, real *pha1, real *pha2, int *shell, 
		      char *mcentric)
{
     void	shellHeader() ;		/* in hklutil.c */

     real	coef, dph, cosdph ;
     real	totF = 0, totd = 0, totC = 0 ;
     real	ctotF = 0, ctotd = 0, ctotC = 0 ;
     real	sumF[NUMSEGS], denom[NUMSEGS] ;
     real	csumF[NUMSEGS], cdenom[NUMSEGS] ;
     real	sumC[NUMSEGS] ;
     real	csumC[NUMSEGS] ;
     int	cases[NUMSEGS] ;
     int	ccases[NUMSEGS] ;
     int	totalcases = 0 ;
     int	ctotalcases = 0 ;

     int	m, n ;

     /*  Note m = 1; m < NUMSEGS (5) - m=0 is (000) term that is not 
         involved in the dphase calculations. */

     for (m = 1; m < NUMSEGS; m++) {	 /* Initialize summation arrays */
	  sumF[m] = 0 ;
	  sumC[m] = 0 ;
	  denom[m] = 0 ;
	  cases[m] = 0 ;
	  csumF[m] = 0 ;
	  csumC[m] = 0 ;
	  cdenom[m] = 0 ;
	  ccases[m] = 0 ;
     }
     /******************************************************
     Calculate phase differences segment by segment.   	
     All phases, as returned by atan2, are in [-180,+180]; 
     force absolute differences to range [0, 180].                            

     6/27/97: CHANGED sumF to use sum of abs phase diff, dph,  
	      instead of sqrt of sum of dph^2
     ******************************************************/

     for (n = 0; n < Nhkl; n++, mcentric++) {

	if ((m = *(shell+n)) > 0) {

	  coef = *(amp+n) ;
          dph = fabs (*(pha1+n) - *(pha2+n)) ; 
	  if (dph > 180.)
	       dph = 360 - dph ;

          cosdph = cos(DTOR*dph) ;

          sumF[m]    += coef*dph ;
          sumC[m]    += coef*cosdph ;
	  denom[m]   += coef ;
	  cases[m]++ ;

          if (*mcentric) {

               csumF[m]    += coef*dph ;
               csumC[m]    += coef*cosdph ;
	       cdenom[m]   += coef ;
	       ccases[m]++ ;

          }
        }
     }

     for (m = 1; m < NUMSEGS; m++) {

	  totF  += sumF[m];
	  totC  += sumC[m];
	  totd  += denom[m];
	  totalcases += cases[m] ;
	  
          if (denom[m] > 0) {
	       sumF[m]   = sumF[m] / denom[m] ;
	       sumC[m]   = sumC[m] / denom[m] ;
          }
          else {
	       sumF[m]   = 0 ;
	       sumC[m]   = 0 ;
          }

	  ctotF  += csumF[m];
	  ctotC  += csumC[m];
	  ctotd  += cdenom[m];
	  ctotalcases += ccases[m] ;
	  
          if (cdenom[m] > 0) {
	       csumF[m]   = csumF[m] / cdenom[m] ;
	       csumC[m]   = csumC[m] / cdenom[m] ;
          }
          else {
	       csumF[m]   = 0 ;
	       csumC[m]   = 0 ;
          }
          
     }

     shellHeader() ;

     sprintf(message,"dph   %10g %10g %10g %10g %10g", 
	    sumF[1], sumF[2], sumF[3], sumF[4], totF/totd) ;
     printTwice(message) ;
     sprintf(message,"dcos     %7.3g    %7.3g    %7.3g    %7.3g    %7.3g", 
	    sumC[1], sumC[2], sumC[3], sumC[4], totC/totd) ;
     printTwice(message) ;
     sprintf(message,"# cases %8d %10d %10d %10d %10d ", 
	    cases[1], cases[2], cases[3], cases[4], totalcases) ;
     printTwice(message) ;
     sprintf(message,
	  "relwt (%%)  %5.2g      %5.2g      %5.2g      %5.2g       100.", 
	  100*denom[1]/totd, 100*denom[2]/totd, 
	  100*denom[3]/totd, 100*denom[4]/totd) ;
     printTwice(message) ;
     if (ctotalcases > 0) {
          printTwice("\n\t\t Restricted phases only: \n") ;
          sprintf(message,"dph   %10g %10g %10g %10g %10g", 
	    csumF[1], csumF[2], csumF[3], csumF[4], ctotF/ctotd) ;
          printTwice(message) ;
          sprintf(message,"dcos     %7.3g    %7.3g    %7.3g    %7.3g    %7.3g", 
	    csumC[1], csumC[2], csumC[3], csumC[4], ctotC/ctotd) ;
          printTwice(message) ;
          sprintf(message,"# cases %8d %10d %10d %10d %10d ", 
	    ccases[1], ccases[2], ccases[3], ccases[4], ctotalcases) ;
          printTwice(message) ;
          sprintf(message,
	       "relwt (%%)  %5.2g      %5.2g      %5.2g      %5.2g       100.", 
	       100*cdenom[1]/ctotd, 100*cdenom[2]/ctotd, 
	       100*cdenom[3]/ctotd, 100*cdenom[4]/ctotd) ;
          printTwice(message) ;
          printTwice(" ") ;
     }
     else
          printTwice("\n\t\t\t (No restricted phases) ") ;
}

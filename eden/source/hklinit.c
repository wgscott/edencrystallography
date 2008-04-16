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

                                HKLINIT.C

  Title:        Hkl space initialization package for EDEN.
  Author:       Hanna Szoke
  Date:         1/23/92
  Function:     This package contains set-up functions for Solve, Back and
		other utilities.  These include: preparing weights and 
		exponential factors for the solver.

******************************************************************************/

#include	"eden.h"
#include	"mir.h"

  /************************************************************
  declarations of EXTERN variables in eden.h               
  ************************************************************/

real  *hkl_weight;        	/* weights for (hkl) in half-ellipsoid */
real	*sigma ;		/* sigma values in fobs file	*/
real  *expfac ;       	/* exponential factors for sc solutions */
COMPLEX *expfac_bcc ;   	/* exponential factors for bcc solutions */




  	/************************************************************
	Examine sigmas:
	Reset useSig if all sigmas are 0. 
	Exit with error if some sigmas are set and others aren't.
	Then, if useSig is off, replace sigma = 0 by 1  everywhere.
	************************************************************/

int checkSigmas(real *sarray, char *mask)
{
     int	s ;
     real	*sig ;
     int	Nsigpresent = 0, Nsigabsent = 0 ;

     if (useSig) {

          for (s = 0, sig = sarray; s < Nhkl; s++, sig++, mask++)  {
	       if (*mask) {
	            if (*sig == 0) 
	                 Nsigabsent++ ;   
		    else
			 Nsigpresent++ ;
               }
          }
          if (Nsigpresent == 0) {
                printTwice("\nAll sigmas are 0 - turning off useSig flag") ;
		useSig = FALSE ;
          }
	  else if (Nsigabsent > 0) {
	       EdenError(
		     "Bad fobs file: some sigmas are 0, others are not.") ;
          }
     }

     if (!useSig) 
          for (s = 0, sig = sarray; s < Nhkl; s++, sig++)  
	       *sig = 1. ;
     
     return(useSig) ;

}

  	/************************************************************
	Determine range of sigmas; 
	Ignore sigma(000) for purposes of reporting range,
	5/1/01: there was an error and we had NOT ignored sigma(000)
	but now that error is corrected.
	************************************************************/

void	SigmaRange(real *sigmin, real *sigmax, real *sarray, char *mask)
{
     int	s ;
     real	smin = 1000. ;
     real	smax = -1000. ;

     sarray++ ;
     mask++ ;

     for (s = 1; s < Nhkl; s++, sarray++, mask++)  
         if (*mask) {
             if (*sarray < smin)
	         smin = *sarray ;
	     if (*sarray > smax)
	         smax = *sarray ;
          }

     if (smin > 0)  {	/* legal range of sigmas */
          sprintf(message, "Sigma range is %g - %g\n", smin, smax) ;
          printTwice(message) ;
     }
     else {				/* something's wrong */
          sprintf(message, 
           "The range of sigmas: (%g, %g) is illegal!", smin, smax) ;
          EdenError(message) ;
     }

     *sigmin = smin ;
     *sigmax = smax ;

}
  	/************************************************************
        Set up default weights in (hkl) space.
	************************************************************/

void setupWeights() 	
{
     int	m, n ;
     real	*wt ;

     strcpy(message, "setupWeights") ;

     hkl_weight = (real *) e_malloc(NM*Nhkl*sizeof(real), message) ;  

     for (n = 0, wt = hkl_weight; n < NM*Nhkl; n++, wt++) 
          *wt = 1 ;
     
     for (m = 0, wt = hkl_weight; m < NM; m++, wt += Nhkl) 
          *wt = 1./ SQRT2  ;
}

  	/************************************************************
        As of 6/17/97, all of the NM (Nder+1) data sets are 
	normalized (with appropriate relative weights) together.
	As of 9/28/98. sigmas = 1 if useSig is not set.
	************************************************************/

real getNormFac() 	
{
     real	SsqWts ;	
     real	Smasks ;	
     real	denom = 0 ;
     real	numer = 0 ;
     int	m, n ;
     real	*wt ;
     char	*mask ;
     real	*sig ;

     wt = hkl_weight ;
     mask = maskfo ;
     sig = sigma ;

     for (m = 0; m < NM; m++) {
	  
          for (n = 0, SsqWts = 0, Smasks = 0; n < Nhkl; 
	       n++, wt++, mask++, sig++) {

               if (*mask) {
	            SsqWts += (*wt / *sig) * (*wt / *sig) ; 
	            Smasks += 1 ;
               }
          }

	  numer += relwt_der[m] * SsqWts ;
	  denom += relwt_der[m] * Smasks ;
          
     }
     return (numer / denom) ;
}

  	/************************************************************
        Set up expfac(Nhkl), defined as: exp (- nusq * dfac) where
        dfac is	usually 
		PISQ * eta * dr * dr
	aka	Beden/4
	dr is the average grid spacing, input_res * RESFAC (basics.c)
	and nusq is the absolute magnitude squared of the reciprocal 
	space vector (h,k,l) in units of (1/Angstrom) squared.  
	It reduces to
		nusq = (h/a)^2 + (k/b)^2 + (l/c)^2 in orthogonal cells.
	See setupNusq() etc. in crysutil.c.
	The value of dfac is adjusted for high resolution and for
	phase extension.
	************************************************************/

real	*setupexpfac(real dfac)
{
     real	*ef ;
     int        n ;
     int        h, k, l ;

     strcpy(message, "setupexpfac") ;

     ef = (real *) e_malloc(Nhkl*sizeof(real), message) ;  

     for (n = 0; n < Nhkl; n++) {

          fetch_hkl(n, &h, &k, &l) ;
          *(ef+n) = exp(-(dfac * nusq(h, k, l))) ;
     }

     return (ef) ;
}
  	/************************************************************
        Set up COMPLEX expfac(Nhkl) for the body-centered array:
        	 exp[-nusq*dfac] * exp[2pi*i*(h/2Nx + k/2Ny + l/2Nz)]
        where dfac is usually = PISQ * eta * dr * dr
	************************************************************/

COMPLEX *setupexpfac_bcc(real dfac)
{
     COMPLEX 	*cef ;
     real	eval ;
     real	bcc_offset ;
     real	x = (real) Nx ;
     real	y = (real) Ny ;
     real	z = (real) Nz ;
     COMPLEX	*pexp ;
     int        n ;
     int        h, k, l ;

     strcpy(message, "setupexpfac_bcc") ;

     cef = (COMPLEX *) e_malloc(Nhkl*sizeof(COMPLEX), message) ;  

     pexp = cef ;

     for (n = 0; n < Nhkl; n++, pexp++) {

          fetch_hkl(n, &h, &k, &l) ;
          eval = exp(-dfac * nusq(h, k, l)) ;
	  bcc_offset = (h/x + k/y + l/z) / 2. ;
          pexp->re = eval * cos(TWOPI*bcc_offset) ;
          pexp->im = eval * sin(TWOPI*bcc_offset) ;
     }

     return (cef) ;
}

	/*********************************************/
       	/* Calculate F and counts segment by segment.*/	
	/*********************************************/

void 	describeFR(real *F, char *masko, char *maskc)
{
     void	shellHeader_F() ;
     void	hklTabHead() ;

     real	sumF[NUMSEGS+1] ;
     int	ocases[NUMSEGS+1], ccases[NUMSEGS+1], dcases[NUMSEGS+1] ;
     int	valid_cases[NUMSEGS+1] ;
     int	m, n ;
     int	h, k, l ;

     for (m = 0; m <= NUMSEGS; m++) {
	  sumF[m] = 0 ;
	  ocases[m] = 0 ;
	  ccases[m] = 0 ;
	  dcases[m] = 0 ;
	  valid_cases[m] = 0 ;
     }
     hklTabHead() ;

     for (n = 0; n < Nhkl; n++) {

	  if (n == 0)
               m = 0 ;

          else {

               fetch_hkl(n, &h, &k, &l) ;

	       for (m = 1; m < NUMSEGS; m++)
		    if (nusq(h, k, l) < nusqlim[m])
			 break ;
          }
	  if (m < NUMSEGS) {   
               sumF[m]  += *(F+n) ;
	       ocases[m] += *(masko+n) ;
	       ccases[m] += *(maskc+n) && *(masko+n) ;
	       dcases[m] += *(maskc+n) && !*(masko+n) ;
	       valid_cases[m] += 1 ;
          }

     }

     for (m = 0; m < NUMSEGS; m++) {
	  sumF[NUMSEGS] += sumF[m];
	  ocases[NUMSEGS] += ocases[m];
	  ccases[NUMSEGS] += ccases[m];
	  dcases[NUMSEGS] += dcases[m];
	  valid_cases[NUMSEGS] += valid_cases[m];
     }

     for (m = 0; m <= NUMSEGS; m++) {
	  sumF[m] = sumF[m] * 100 / sumF[NUMSEGS] ;
     }

     shellHeader_F() ;

     fprintf(fp_log, 
            "%% F %9g %10g %10g %10g %10g %10g\n", 
	    sumF[0], sumF[1], sumF[2], sumF[3], sumF[4], sumF[5]) ;
     fprintf(fp_log, 
            "# Fo    %5d %10d %10d %10d %10d %10d\n", 
	    ocases[0], ocases[1], ocases[2], ocases[3], ocases[4], ocases[5]) ;
     fprintf(fp_log, 
            "# Fc in %5d %10d %10d %10d %10d %10d\n", 
	    ccases[0], ccases[1], ccases[2], ccases[3], ccases[4], ccases[5]) ;
     fprintf(fp_log, 
            "# Fc out%5d %10d %10d %10d %10d %10d\n", 
	    dcases[0], dcases[1], dcases[2], dcases[3], dcases[4], dcases[5]) ;
     fprintf(fp_log, 
            "max #   %5d %10d %10d %10d %10d %10d\n", 
	    valid_cases[0], valid_cases[1], valid_cases[2], 
	    valid_cases[3], valid_cases[4], valid_cases[5]) ;

     fflush(fp_log) ;
     return ;
}

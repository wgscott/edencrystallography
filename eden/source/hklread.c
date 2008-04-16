/****************************************************************************

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

                                HKLREAD.C

  Title:        Hkl read package for EDEN.
  Author:       Hanna Szoke
  Date:         5/18/2000 separated from hklutil.c
  Function:     This package contains functions that read structure factors:

			readfobs()
			readfobs0_allh()
			readfobs1_allh()
			readfobs2_allh()
			readfcalc()
			readfcalc1_allh()
			readfcalc2_allh()
		and low-level functions called by them.

****************************************************************************/
#include	"util.h"
#include	"symmetry.h"
#include	"dims.h"

#define	COS_EPSI 1.e-6		/* for identifying phases ~= to N*pi */

struct Ntypes {		
     int used ;
     int total ;
} ;

FILE	*fp ;
static	struct	Ntypes	N ;
static	int	nof ;
static	int	fobs_mismatch_count ;	/* counter of errors in fobs */
static	int	fcalc_mismatch_count ;	/* counter of errors in fcalc */
static	char	newform ;
static	char	sigflag ;

static	void	acc_fo(real *,real,char *,real *,real,int,int,int) ;
static	void	average_centrics(real *, real *, char *) ;
static	int	enforce_P1(char *) ;
static	void	fc_conflict(int, COMPLEX, COMPLEX) ;
static	void	finish_expandfo(real *, char *, real *, int *, int *) ; 
static	void	fo_conflict(int, real, real, real, real) ;
static	void	handle_forbidden_fc(COMPLEX *, char *) ; 
static	int	handle_forbidden_fo(real *, char *, real *, real) ; 
static	void	init_fc(COMPLEX *, char *) ;
static	void	init_fo(real *, real *, char *) ;
static	char	read_1_fobs(char *,int *,int *,int *,real *,real *);
static	void	set_fc(COMPLEX *,real,real,char *,int,int,int,real) ;
static	void	set_fc_noFriedel(int, COMPLEX *, real, real,
			         char *, int, int, int, real) ; 
static	void    set_fo(real *,real,char *,real *,real,int,int,int) ;
static	void    set_fo_noFriedel(int, real *, real, char *,
			 real *, real, int, int, int) ;
static	void	set_forbid(char *) ;
static	void	check_centric_phases(COMPLEX *, char *, char *) ;

int	assemble(int, int, int) ;	/* used frequently */

	/*********************************************
	Read experimental (.fo) diffraction pattern.            
	*********************************************/
					
int readfobs(char *name, real *array, char *mask, real *sarray)
{
     int	readfobs0_allh(char *, real *, char *, real *) ;
     int	readfobs1_allh(char *, real *, char *, real *, char *) ;
     int	n ;
     int	Nused ;
     char	*mcentric ;

     /* Deal with "empty" */

     if (strcmp(name, "empty") == 0) {
	  *(array)  = 1. ;
	  *(mask)   = 1 ;
	  *(sarray) = 1. ;

          for (n = 1; n < Nhkl; n++) {
               *(array+n) = 0. ;
               *(mask+n) = 0 ;
               *(sarray+n) = 0 ;	
          }
          printTwice("Fobs(0,0,0) = 1, Sig(0,0,0) = 1 - empty (dummy) file") ;
	  return(0) ;
     }

     /* Be sure it's there; check for 2 distinct formats and presence of sigmas */

     check_exist(name) ;

     fp = fopen(name, "r") ;

     while (fgets(nextline, MAXSTRING, fp) != NULL) {
	  if (((int)strlen(nextline) > 1) &&
              ((strstr(nextline, "IND") != NULL) ||
              ((strstr(nextline, "ind") != NULL)))) {

	      if (((strstr(nextline, "SIG") != NULL) ||
                  ((strstr(nextline, "sig") != NULL)))) {

		  if (sscanf(nextline, "%*s %f %f %f %*s %f %*s %f", 
				    &t1, &t2, &t3, &t4, &t5) == 5)  {
                       newform = FALSE ; 	/* old format, with sigmas */
		       sigflag = TRUE ;
		       break ;
                  }
		  else if (sscanf(nextline, "%*s %f %f %f %*s %f %f %*s %f", 
				    &t1, &t2, &t3, &t4, &t5, &t6) == 6) {
		       newform = TRUE ;		/* new format, with sigmas */
		       sigflag = TRUE ;
		       break ;
                  }
               }
	       else {

		  if (sscanf(nextline, "%*s %f %f %f %*s %f", 
				    &t1, &t2, &t3, &t4) == 4) {
                       newform = FALSE ; 	/* old format, without sigmas */
		       sigflag = FALSE ;
		       break ;
                  }
		  else if (sscanf(nextline, "%*s %f %f %f %*s %f %f",
				    &t1, &t2, &t3, &t4, &t5) == 5) {
		       newform = TRUE ;		/* new format, without sigmas */
		       sigflag = FALSE ;
		       break ;
                  }
               }
           }

/***           else { 	apparently header line
           }				***/
     }
          
     fclose(fp) ;

     /* Do the reading. */

     if (anom_flag) {
          mcentric = (char *) e_malloc(Nhkl*sizeof(char), caller) ;
          Nused = readfobs1_allh(name,array,mask,sarray,mcentric) ;
	  e_free(Nhkl*sizeof(char), mcentric) ;
     }
     else
          Nused = readfobs0_allh(name, array, mask, sarray) ;

     return (Nused) ;
}

	/*********************************************************
	Read and expand experimental (.fo) diffraction pattern.  
	Handle regular diffraction (no anomalous dispersion) 
	with one set of unique points in the half-ellipsoid.
	*********************************************************/
					
int readfobs0_allh(char *name, real *array, char *maskfo, real *sigma)
{
     FILE	*fp ;
     int	Nunique ;
     int	Nforbid ;
     real	amplitude ;
     real	sigvalue ;
     real	off ;		/* returned by apply_symop_hkl(), unused */
     real	*mat ;
     int        n, h, k, l ;
     int        newh, newk, newl ;


     init_fo(array, sigma, maskfo) ;

	/***********************************************************
	Read experimental data, expand to P1.
        Generate (h,k,l) triplets such that the Bragg condition
        is satisfied and at most only one of the two reciprocal
        lattice points -- (h,k,l) and (-h,-k,-l) -- appears.  
	Arbitrarily, we select h >= 0.  We use Friedel's law.
	***********************************************************/
					
     fp = fopen(name, "r") ;

     while (fgets(nextline, MAXSTRING, fp) != NULL)  {

	  if (read_1_fobs(nextline, &h, &k, &l, &amplitude, &sigvalue)){

	     N.total++ ;

             if (nusq(h, k, l) < limit_rat) {

                 N.used++ ;

	         for (n = 0,mat = matop; n < Ncs; n++,mat += MEL) {
			       
		     apply_symop_hkl(h, k, l, mat,
				   &newh, &newk, &newl, &off) ;

	             set_fo(array, amplitude, maskfo, 
				   sigma, sigvalue, newh, newk, newl) ;

                 } 
             }	/* end check for inclusion */
         }	/* end check for something on input line */
     }		/* end check for EOF	*/

     fclose(fp) ;

     finish_expandfo(array, maskfo, sigma, &Nunique, &Nforbid) ;
     Nunique += Nforbid ;

     if (Ncs > 1)
        sprintf(message, 
	"%d (expanded to %d) out of %d Fobs entries read in;",
	    N.used, Nunique, N.total) ;
     else
        sprintf(message, "%d out of %d Fobs entries read in;",
	    N.used, N.total) ;
     printTwice(message) ;

     sprintf(message, "Fobs(0,0,0) = %g, sigma = %g", *array, *sigma) ;
     printTwice(message) ;
     
     if (fobs_mismatch_count > 0) {
          sprintf(message, "There were a total of %d mismatches", 
			   fobs_mismatch_count) ;
          printTwice(message) ;
     }

     if (fobs_mismatch_count > 0.25*N.used) 
          EdenError(
  "readfobs did not work properly - Please check the ANOM flag setting!") ;

     return(Nunique) ;
}

	/*********************************************************
	Read and expand experimental (.fo) diffraction pattern.  
	Handle anomalous dispersion with one set of unique points 
	in the half-ellipsoid.  Centrics are averaged; they should
	all be the same, but if we assume that they are, we may 
	get large numbers of mismatches.
	*********************************************************/
					
int readfobs1_allh(char *name, real *array, char *maskfo, real *sigma, 
		   char *mcentric)
{
     void	set_centrics(char *) ;

     int	Nunique ;
     int	Nforbid ;
     real	amplitude ;
     real	sigvalue ;
     real	off ;		/* returned by apply_symop_hkl(), unused */
     real	*mat ;
     int        n, h, k, l ;
     int        nn, newh, newk, newl ;


     /* Initialize arrays */

     init_fo(array, sigma, maskfo) ;
     set_centrics(mcentric) ;

	/***********************************************************
	Read experimental data, expand to P1.
        Generate (h,k,l) triplets such that the Bragg condition
        is satisfied and at most only one of the two reciprocal
        lattice points -- (h,k,l) and (-h,-k,-l) -- appears.  
	Arbitrarily, we select h >= 0.  We do not use Friedel's law 
	but instead, expand only h >= 0.  We accumulate centrics, 
	afterwards average them.
	***********************************************************/
					
     fp = fopen(name, "r") ;

     while (fgets(nextline, MAXSTRING, fp) != NULL)  {
	      
	  if (read_1_fobs(nextline, &h, &k, &l, &amplitude, &sigvalue)){

	     N.total++ ;

             if (nusq(h, k, l) < limit_rat) {

                 N.used++ ;

	         for (n = 0,mat = matop; n < Ncs; n++,mat += MEL) {
			       
		     apply_symop_hkl(h, k, l, mat,
				   &newh, &newk, &newl, &off) ;

                     nn = assemble(newh, newk, newl) ;
		     if ((nn >=0) && (nn < Nhkl)) {
                     if (*(mcentric+nn) == 0)
		          set_fo_noFriedel(TRUE, array, amplitude, maskfo, 
				      sigma, sigvalue, newh, newk, newl) ;
                     else
		          acc_fo(array, amplitude, maskfo, 
				      sigma, sigvalue, newh, newk, newl) ; 
                     }
                 } 

             }	/* end check for inclusion */
         }	/* end check for something on input line */
     }		/* end check for EOF	*/

     fclose(fp) ;

     average_centrics(array, sigma, mcentric) ; 

     finish_expandfo(array, maskfo, sigma, &Nunique, &Nforbid) ;
     Nunique += Nforbid ;

     if (Ncs > 1)
        sprintf(message, 
	"%d (expanded to %d) out of %d Fobs entries read in;",
	    N.used, Nunique, N.total) ;
     else
        sprintf(message, "%d out of %d Fobs entries read in;",
	    N.used, N.total) ;
     printTwice(message) ;
     sprintf(message, "Fobs(0,0,0) = %g", *array) ;
     printTwice(message) ;
     
     if (fobs_mismatch_count > 0) {
          sprintf(message, "There were a total of %d mismatches", 
			   fobs_mismatch_count) ;
          printTwice(message) ;
     }

     if (fobs_mismatch_count > 0.25 * N.used) 
          EdenError(
  "Expandfo did not work properly - Please check the ANOM flag setting!") ;

     return (Nunique) ;
}

	/*********************************************************
	Read and expand experimental (.fo) diffraction pattern.  
	for anomalous dispersion with two sets of unique points in 
	the half-ellipsoid.
	*********************************************************/

int readfobs2_allh(char *name, real *Fplus, real *Fminus, char *mask_plus, 
		   char *mask_minus, real *splus, real *sminus)
{
     real	amplitude ;
     real	sigvalue ;
     real	off ;		/* returned by apply_symop_hkl(), unused */
     real	*mat ;
     int        n, h, k, l ;
     int        newh, newk, newl ;
     int	Nuniq_plus, Nuniq_minus, Nunique ;
     int	Nforbid ;

     /* Open file */

     fp = fopen(name, "r") ;

     /* Initialize arrays */

     init_fo(Fplus, splus, mask_plus) ;
     init_fo(Fminus, sminus, mask_minus) ;

	/***********************************************************
	Read experimental diffraction pattern.
	Using rotation/translation matrices, expand to P1 using both 
	positive and negative half-ellipsoids, separately.
	***********************************************************/
					
     while (fgets(nextline, MAXSTRING, fp) != NULL)  {
	      
	  if (read_1_fobs(nextline, &h, &k, &l, &amplitude, &sigvalue)){

	     N.total++ ;

             if (nusq(h, k, l) < limit_rat) {

                 N.used++ ;

	         for (n = 0,mat = matop; n < Ncs; n++,mat += MEL) {
			       
		     apply_symop_hkl(h, k, l, mat,
				   &newh, &newk, &newl, &off) ;

		     set_fo_noFriedel(TRUE, Fplus, amplitude, mask_plus, 
				      splus, sigvalue, newh, newk, newl) ;
		     set_fo_noFriedel(FALSE, Fminus, amplitude, mask_minus, 
			              sminus, sigvalue, newh, newk, newl) ;

                 } 

              }		/* end check for octant */
          }		/* end check for something legal on input line */
     }

     fclose(fp) ;

     finish_expandfo(Fplus, mask_plus, splus, &Nuniq_plus, &Nforbid) ;
     Nuniq_plus += Nforbid ;
     finish_expandfo(Fminus, mask_minus, sminus, &Nuniq_minus, &Nforbid) ;
     Nuniq_minus += Nforbid ;

     Nunique = Nuniq_plus + Nuniq_minus ; 

     if (Ncs > 1)
        sprintf(message, "%d (expanded to %d) out of %d Fobs entries read in;",
	    N.used, Nunique, N.total) ;
     else
        sprintf(message, "%d out of %d Fobs entries read in;",
	    N.used, N.total) ;
     printTwice(message) ;
     sprintf(message, "Fobs+(0,0,0) = %g", *Fplus) ;
     printTwice(message) ;
     sprintf(message, "Fobs-(0,0,0) = %g", *Fminus) ;
     printTwice(message) ;
     
     if (fobs_mismatch_count > 0) {
          sprintf(message, "There were a total of %d mismatches", 
			   fobs_mismatch_count) ;
          printTwice(message) ;
     }

     return(Nunique) ;
}

char read_1_fobs(char *line, int *h,int *k,int *l,real *amp, real *sig)
{
    if (((int)strlen(line) > 1) &&
        ((strstr(line, "IND") != NULL) ||
        ((strstr(line, "ind") != NULL)))) {

         if (newform)
              sscanf(line, "%*s %d %d %d %*s %g %*g %*s %g", h, k, l, &t1, &t2) ;
	 else
              sscanf(line, "%*s %d %d %d %*s %g %*s %g", h, k, l, &t1, &t2) ;

	 *amp = t1 ;

         if (sigflag) 
	     *sig = t2 ;
         else
	     *sig = 0 ;

         if (*amp < 0) {
             EdenWarning(
	   "Your input fobs file contains entries with amplitude < 0;") ;
             EdenError("Quitting.") ;
         }
         return(TRUE) ;
     }
     else
         return(FALSE) ;

}
	/***********************************************************
	Use Friedel's law to identify a position in the h >= 0 half
	ellipsoid.  Check that there is no conflict at that point.
	***********************************************************/

void    set_fo(real *F, real value, char *maskfo, real *sigma, 
	               real sigvalue, int h, int k, int l)
{
     int        n ;
	   
     if (h >= 0) {
          n = assemble(h, k, l) ;
				  
          if (n >= 0) {
	      if (*(maskfo + n) == 0) {
                   *(F+n) = value ;
                   *(maskfo + n) = 1 ;
                   *(sigma+n) = sigvalue ;
              }
	      else {
	          fo_conflict(n, *(F+n), value, *(sigma+n), sigvalue) ;
	     }
         }
     }
     if (h <= 0) {
	  n = assemble(-h, -k, -l) ;

          if (n >= 0) {
	      if (*(maskfo + n) == 0) {
                   *(F+n) = value ;
                   *(maskfo + n) = 1 ;
                   *(sigma+n) = sigvalue ;
              }
	      else {
	          fo_conflict(n, *(F+n), value, *(sigma+n), sigvalue) ;
	     }
         }
     }
}

       /*********************************************************
       Determine forbidden reflections and apply them to maskfo;
       report errors (forbidden points at which fobs was non-zero).

       We use the information in maskfo to set up a value for
       min (limit_rat, nusqmax); this should be used
       in place of limit_rat for determining how many forbidden
       points to set!  
       ***********************************************************/

int	handle_forbidden_fo(real *array, char *maskfo, real *sigma, 
			    real sigvalue) 
{
     void       setDsqLimit(char *, int) ;

     char	*mforbid ;
     int	n, h, k, l ;
     int	countm = 0 ;

     mforbid = (char *) e_malloc(Nhkl*sizeof(char), "handle_forbidden_fo") ;

     for (n = 0; n < Nhkl; n++)
	  *(mforbid+n) = 0 ;

     setDsqLimit(maskfo, Nhkl) ;
     set_forbid(mforbid) ;
     n = enforce_P1(mforbid) ;

     for (n = 0; n < Nhkl; n++) 
	  if (*(mforbid+n)) {
	       if (*(maskfo+n)) { 
		    if (*(array+n) > 0) {

			 fetch_hkl(n, &h, &k, &l) ;
			 sprintf(message,
     "Fobs = %g at forbidden reflection! at (%d,%d,%d) - resetting", 
      *(array+n), h, k, l) ;
			 EdenWarning(message) ;

			 *(array+n) = 0. ;
			 *(sigma+n) = sigvalue ;
                    }
               }
               else {
	          countm++ ;
		  *(maskfo+n) = 1 ;
		  *(sigma+n) = sigvalue ;
               }
          }

     e_free(Nhkl*sizeof(char), mforbid) ;

     return (countm) ;
}

	/***********************************************************
	Do NOT use Friedel's law to identify a position in the 
	h >= 0 half ellipsoid or in the h <= 0 half-ellipsoid.  
	Check that there is no conflict at that point.
	***********************************************************/

void    set_fo_noFriedel(int plus, real *F, real value, char *maskfo,
			 real *sigma, real sigvalue, int h, int k, int l)
{
     int        n = 0 ;
     int	doit = FALSE ;
	   
     if (plus && (h >= 0)) {
          n = assemble(h, k, l) ;
          if (n >= 0) 
	       doit = TRUE ;
     }
     else if (!plus && (h <= 0)) {
          n = assemble(-h, -k, -l) ;
          if (n >= 0) 
	       doit = TRUE ;
     }

     if (doit) {
	 if (*(maskfo + n) == 0) {
              *(F+n) = value ;
              *(maskfo + n) = 1 ;
              *(sigma+n) = sigvalue ;
         }
	 else {
	      fo_conflict(n, *(F+n), value, *(sigma+n), sigvalue) ;
	 }
     }
}
	/***********************************************************
	Do NOT use Friedel's law to identify a position in the 
	h >= 0 half ellipsoid or in the h <= 0 half-ellipsoid.  
	Handle centrics by accumulation of SQUARES of sigmas.
	***********************************************************/ 
void    acc_fo(real *F, real value, char *mask, 
	      real *sigma, real sigvalue, int h, int k, int l)
{
     int        n = 0 ;
     int	doit = FALSE ;
	   
     if  (h >= 0) {
          n = assemble(h, k, l) ;
          if (n >= 0) 
	       doit = TRUE ;
     }

     if (doit) {
          *(F+n) += value ;
          *(mask + n) = 1 ;
          *(sigma+n) += sigvalue*sigvalue ;
     }
}

void	average_centrics(real *F, real *sigma, char *mcentric)
{
     int        n ;

	/*********************************************************
     	This looks wrong but it's okay: centric sigmas have their 
	SQUARES accumulated into *sigma in acc_fo() prior to the 
	call to this function, so here we restore the true sigmas.
	*********************************************************/

     for (n = 0; n < Nhkl; n++) {
	  if (*(mcentric+n) > 0) {
	       *(F+n) /= 2 ;
	       *(sigma+n) = sqrt(*(sigma+n)/2) ;
          }
     }
}
	/*********************************************************
	Do some housekeeping: Check that f(000) is set.
        Enforce uniqueness on h=0 plane; set forbidden reflections.
        Give forbidden reflections the best (minimum) sigma around.
	As of 9/28/98, if useSig is false, forbidden sigmas are
	set to 1 (sigmin)  rather than to zero.
        Unnecessarily but harmlessly, apply maskfo to the array. 
	*********************************************************/

void	finish_expandfo(real *array, char *maskfo, real *sigma, 
			        int *Nunique, int *Nforbid) 
{
     void	checkF000(real *, real *, char *) ;
     int	checkSigmas(real *, char *) ;

     int	n ;
     real	sigmin = 1 ;


     checkF000(array, sigma, maskfo) ;
     *Nunique = enforce_P1(maskfo) ;

     useSig = checkSigmas(sigma, maskfo) ;

     *Nforbid = handle_forbidden_fo(array, maskfo, sigma, sigmin) ;

     for (n = 1; n < Nhkl; n++)
          *(array+n) *= *(maskfo+n) ;

}

int readfcalc(char *name, COMPLEX *Fc, char *mask)
{
     int	readfcalc1_allh(char *, COMPLEX *, char *) ;

     int	Nused ; 

     if ((anom_flag) && (laue_group == 1)) 
          EdenError( "Please run expandfc to handle your anomalous data") ;

     init_fc(Fc, mask) ;

     if (strcmp(name, "empty") == 0) {
          printTwice( "Fcalc(0,0,0) = (0, 0) - empty model") ;
	  return(0) ;
     }
     else
        check_exist(name) ;

     Nused = readfcalc1_allh(name, Fc, mask) ;

     return(Nused) ;
}
	/***********************************************************
	Oversee reading and expansion of model structure factors for 
	either regular diffraction or anomalous dispersion with
	one set of unique data points in the half-ellipsoid. 
	***********************************************************/

int	readfcalc1_allh(char *name, COMPLEX *R1, char *maskfc) 
{
     real	cleanPhase(real) ;
     void	set_centrics(char *) ;
     void	Cmaskit(COMPLEX *, char *, int) ;

     int	Nex ;
     real	amp, phase, inphase ;
     real	off ;
     real	*mat ;
     int        n, h, k, l ;
     int	newh, newk, newl ;
     char	*mcentric ;

     init_fc(R1, maskfc) ;

	/***********************************************************
	Read model structure factors.
	Using rotation/translation matrices, expand to P1.
	If the anomalous flag is not set, use Friedel's law;
	otherwise, expand only the non-negative h's.
	***********************************************************/

     fp = fopen(name, "r") ;

     while (fgets(nextline, MAXSTRING, fp) != NULL) {
	  if (((int)strlen(nextline) > 1) &&
              ((strstr(nextline, "IND") != NULL) ||
              ((strstr(nextline, "ind") != NULL)))) {

	     N.total++ ;
             nof = sscanf(nextline, "%*s %d %d %d %*s %g %g", 
			         &h, &k, &l, &t1, &t2) ;
             if (nof != 5)
                  EdenError("Wrong number of fields in your input fcalc file!") ;
             amp=t1 ;
	     inphase = t2 ;
             phase = cleanPhase(inphase) ;

             if (nusq(h, k, l) < limit_rat) {

                 N.used++ ;

                 for (n = 0,mat = matop; n < Ncs; n++,mat += MEL) {

		      apply_symop_hkl(h, k, l, mat, 
				      &newh, &newk, &newl, &off) ;

                      if (!anom_flag)
                         set_fc(R1, amp, phase, maskfc, newh, newk, newl, off) ;
                      else
		         set_fc_noFriedel(TRUE, R1, amp, phase, maskfc, 
				       newh, newk, newl, off) ;

                 }
             }		/* end check for resolution */
         }		/* end check for something on the input line */
     }			/* end check for end-of-file */

     Nex = enforce_P1(maskfc) ;
     handle_forbidden_fc(R1, maskfc) ;

     if (anom_flag) {
          mcentric = (char *) e_malloc(Nhkl*sizeof(char), caller) ;
          set_centrics(mcentric) ;
	  check_centric_phases(R1, mcentric, maskfc) ;
	  e_free(Nhkl*sizeof(char), mcentric) ;
     }

     for (n = 0, Nex = 0; n < Nhkl; n++) 
          Nex += *(maskfc+n) ; 

     Cmaskit(R1, maskfc, Nhkl) ;

     sprintf(message, "%d out of %d Fcalc entries read in;", N.used, N.total) ;
     printTwice(message) ;
     sprintf(message, "%d Fcalc entries are valid,", Nex) ;
     printTwice(message) ;
     if (R1->re > 0) {
          sprintf(message, "Fcalc(0,0,0) = (%g, %g)", R1->re, R1->im) ;
          printTwice(message) ;
     }
     else 
	  EdenError("Missing F(000) in your model!") ;
     
     if (fcalc_mismatch_count > 0) {
          sprintf(message, "There were a total of %d mismatches", 
			   fcalc_mismatch_count) ;
          printTwice(message) ;
     }
     fclose(fp) ;

/**     if (fcalc_mismatch_count > 0.25 * N.used)  **/
     if (fcalc_mismatch_count > 0.5 * N.used) 
          EdenError(
  "Expansion did not work properly - Please check the ANOM flag setting!") ;

     return(Nex) ;
}

	/***********************************************************
	Oversee reading and expansion of model structure factors for 
	anomalous dispersion with two sets of unique data points 
	in the half-ellipsoid. 
	***********************************************************/

int	readfcalc2_allh(char *name, COMPLEX *Rplus, COMPLEX *Rminus, 
			char *maskfc_plus, char *maskfc_minus)
{
     FILE	*fp ;

     real	cleanPhase(real) ;
     void	Cmaskit(COMPLEX *, char *, int) ;

     real	amp, phase, inphase ;
     real	off ;
     real	*mat ;
     int        n, h, k, l ;
     int	newh, newk, newl ;
     int	Nunique ;

     check_exist(name) ;

     init_fc(Rplus, maskfc_plus) ;
     init_fc(Rminus, maskfc_minus) ;

	/***********************************************************
	Read model structure factors.
	Using rotation/translation matrices, expand to P1
	using both positive and negative half-ellipsoids,
	using "mock" Friedel's law.
       ***********************************************************/
	
     fp = fopen(name, "r") ;

     while (fgets(nextline, MAXSTRING, fp) != NULL) {
	  if (((int)strlen(nextline) > 1) &&
              ((strstr(nextline, "IND") != NULL) ||
              ((strstr(nextline, "ind") != NULL)))) {

	     N.total++ ;
             nof = sscanf(nextline, "%*s %d %d %d %*s %g %g", 
			         &h, &k, &l, &t1, &t2) ;
             if (nof != 5)
                  EdenError("Wrong number of fields in your input fcalc file!") ;
             amp=t1 ;
	     inphase = t2 ;
             phase = cleanPhase(inphase) ;

     /* Look for (h,k,l) triplets among the unique points of the ellipsoid */

             if (nusq(h, k, l) < limit_rat) {

                 N.used++ ;

                 for (n = 0,mat = matop; n < Ncs; n++,mat += MEL) {

		      apply_symop_hkl(h, k, l, mat, 
				      &newh, &newk, &newl, &off) ;
		      set_fc_noFriedel(TRUE, Rplus, amp, phase, maskfc_plus, 
		                            newh, newk, newl, off) ;
		      set_fc_noFriedel(FALSE, Rminus, amp, phase, maskfc_minus,
		                            newh, newk, newl, off) ;

		 }
             }		/* end check for non-negative octant of ellipsoid */
        }		/* end check for something on the input line */
     }			/* end check for end-of-file */

     Nunique = enforce_P1(maskfc_plus) ;
     Nunique += enforce_P1(maskfc_minus) ;

     handle_forbidden_fc(Rplus, maskfc_plus) ;
     handle_forbidden_fc(Rminus, maskfc_minus) ;

     Cmaskit(Rplus, maskfc_plus, Nhkl) ;
     Cmaskit(Rminus, maskfc_minus, Nhkl) ;

     if (Ncs > 1)
        sprintf(message, "%d (expanded to %d) out of %d Fcalc entries read;",
	    N.used, Nunique, N.total) ;
     else
        sprintf(message, "%d out of %d Fcalc entries read;",
	    N.used, N.total) ;
     printTwice(message) ;
     if (Rplus->re > 0) {
          sprintf(message, "Fcalc(0,0,0) = (%g, %g)", Rplus->re, Rplus->im) ;
          printTwice(message) ;
     }
     else 
	  EdenError("Missing F(000) in your model!") ;
     
     
     if (fcalc_mismatch_count > 0) {
          sprintf(message, "There were a total of %d mismatches", 
			   fcalc_mismatch_count) ;
          printTwice(message) ;
     }

     fclose(fp) ;

     if (fcalc_mismatch_count > 0.25 * N.used) 
          EdenError(
  "Expandfc did not work properly - Please check the ANOM flag setting!") ;

     return(Nunique) ;
}

	/***********************************************************
	Use Friedel's law to identify a position in the h >= 0 half
	ellipsoid.  Check that there is no conflict at that point.
	***********************************************************/

void	set_fc(COMPLEX *R, real amp, real phase, char *maskfc, 
	               int h, int k, int l, real off) 
{
     int	n ;
     COMPLEX	newR ;

     if (h >= 0) {
     
          newR.re = amp * cos(DTOR*phase - off) ;
          newR.im = amp * sin(DTOR*phase - off) ;

          n = assemble(h, k, l) ;

          if (n >= 0) {
	       if (*(maskfc + n) == 0) {
                    (R+n)->re = newR.re ;
                    (R+n)->im = newR.im ;
                    *(maskfc + n) = 1 ;
               }
	       else 
                    fc_conflict(n, *(R+n), newR) ;
          }
     }
     if (h <= 0) {
     
          newR.re = amp * cos(DTOR*phase - off) ;
          newR.im = -amp * sin(DTOR*phase - off) ;	/* complex conj */

          n = assemble(-h, -k, -l) ;

          if (n >= 0) {
	       if (*(maskfc + n) == 0) {
                    (R+n)->re = newR.re ;
                    (R+n)->im = newR.im ;
                    *(maskfc + n) = 1 ;
               }
	       else 
                    fc_conflict(n, *(R+n), newR) ;
          }
     }
}

	/***********************************************************
	Do NOT use Friedel's law to identify a position in the 
	h >= 0 half ellipsoid or use "mock Friedel" (whatever
	that means) in the h <= 0 half-ellipsoid.  
	Check that there is no conflict at that point.
	***********************************************************/

void	set_fc_noFriedel(int plus, COMPLEX *R, real amp, real phase,
		         char *maskfc, int h, int k, int l, real off) 
{
     int	n = 0 ;
     COMPLEX	newR ;
     int	doit = FALSE ;

     if (plus && (h >= 0)) {

          newR.re = amp * cos(DTOR*phase - off) ;
          newR.im = amp * sin(DTOR*phase - off) ;

          n = assemble(h, k, l) ;
          if (n >= 0) 
	      doit = TRUE ;

     }
     if (!plus && (h <= 0)) {

          newR.re = amp * cos(DTOR*phase - off) ;
          newR.im = -amp * sin(DTOR*phase - off) ;	/* complex conj */

          n = assemble(-h, -k, -l) ;
          if (n >= 0) 
	      doit = TRUE ;

     }

     if (doit) {

	  if (*(maskfc + n) == 0) {
               (R + n)->re = newR.re ;
               (R + n)->im = newR.im ;
               *(maskfc + n) = 1 ;
          }
	  else 
               fc_conflict(n, *(R+n), newR) ;
     }
}

       /*********************************************************
       Determine forbidden reflections and apply them to maskfc;
       report errors (forbidden points at which fcalc was non-zero).

       We use the information in maskfc to set up a value for
       min (limit_rat, nusqmax); this should be used
       in place of limit_rat for determining how many forbidden
       points to set!  
       ***********************************************************/

void	handle_forbidden_fc(COMPLEX *R, char *maskfc) 
{
     void       setDsqLimit(char *, int) ;

     char	*mforbid ;
     int	n, h, k, l ;

     mforbid = (char *) e_malloc(Nhkl*sizeof(char), "expandfc") ;

     for (n = 0; n < Nhkl; n++)
	  *(mforbid+n) = 0 ;

     setDsqLimit(maskfc, Nhkl) ;

     set_forbid(mforbid) ;
     n = enforce_P1(mforbid) ;

     for (n = 0; n < Nhkl; n++) 
	  if (*(mforbid+n)) {
	       if (*(maskfc+n)) {
		    if ((R+n)->re > 0) {

			 fetch_hkl(n, &h, &k, &l) ;
			 sprintf(message,
     "Fcalc = (%g %g) at forbidden reflection! at (%d,%d,%d) - resetting", 
      (R+n)->re, (R+n)->im, h, k, l) ;
			 EdenWarning(message) ;

			 (R+n)->re = 0. ;
			 (R+n)->im = 0. ;
                    }
               }
               else
		  *(maskfc+n) = 1 ;
          }

     e_free(Nhkl*sizeof(char), mforbid) ;
}


	/***********************************************************
	Check centric phases; if all are 0 (SAD) - eliminate them!
	***********************************************************/

void	check_centric_phases(COMPLEX *R, char *mcentric, char *maskfc)
{
     real amp0, Ramp, Rphase ;
     int    n, N0, Ncen ;

     amp0 = R->re ;

     for (n = 1, Ncen = 0, N0 = 0; n < Nhkl; n++) 
	  if (*(maskfc+n) && *(mcentric+n)) {
	       Ncen++ ;
	       getAmpPhase(R+n, amp0, &Ramp, &Rphase) ;
	       if (Rphase == 0)
		    N0++ ;
          }

     if ((N0 == Ncen) && (Ncen > 0)) {
	  EdenWarning(
     "All the centric reflections for your anomalous data have 0 phases!") ;
	  EdenWarning("They will NOT be used...") ;

          for (n = 1; n < Nhkl; n++) 
	       if (*(maskfc+n) && *(mcentric+n)) 
		    *(maskfc+n) = 0 ;
     }
}

	/***********************************************************
	Set mask to show forbidden reflections (systematic absences).
	We look for (hkl)'s which, under a symmetry operation, 
	give back the same reciprocal indices with a phase that 
	differs significantly from 2*pi.
	We use Ncs, not Nprim, for counting.
	***********************************************************/

void set_forbid(char *mforbid)
{
     real	off, coff ;
     real	*mat ;
     int	N = 0 ;
     int        p, n, h, k, l ;
     int	newh, newk, newl ;

     if (dsq_limit == 0)
          EdenError("Dsq_limit is 0! - Please check the log.") ;
		
     for (n = 0; n < Nhkl; n++) 
          *(mforbid+n) = 0 ;

     for (n = 1; n < Nhkl; n++) {

          fetch_hkl(n, &h, &k, &l) ;

          if ((nusq(h, k, l) < dsq_limit) && (*(mforbid+n) == 0)) {

               for (p = 0, mat = matop; p < Ncs; p++, mat += MEL) {

	           apply_symop_hkl(h, k, l, mat, &newh, &newk, &newl, &off) ;

	           coff = cos(off) - 1.;

	           if ((newh == h) && (newk == k) && (newl == l) && 
		       (fabs(coff) > COS_EPSI)) {
		       *(mforbid+n) = 1 ;
		       break ;
                   }		
               }	
          }	
     }			
     for (n = 1; n < Nhkl; n++) 
	  N += *(mforbid+n) ;
}

#define	MAX_MISMATCH	20 	/* # of mismatches reported in non-verbose 
				mode */
#define	DEPS	1.e-2 		/* for judging relative errors in 
				structure factors that should agree */

void fo_conflict(int n, real oldF, real newF, real oldsig, 
			 real newsig)
{
      int     h, k, l ;

      if ((fabs(newF - oldF) > DEPS*(newF + oldF)) ||
          (fabs(newsig - oldsig) > DEPS*(newsig + oldsig))) {

           fobs_mismatch_count++ ;

	   if ((verbose) || (fobs_mismatch_count < MAX_MISMATCH)) {
                fetch_hkl(n, &h, &k, &l) ;
		fprintf(fp_log, "FOBS conflict at h=%d, k=%d, l=%d", h,k,l) ;
		fprintf(fp_log, " old F(sig)= %g (%g), new F(sig) = %g (%g)\n",
		     oldF, oldsig, newF, newsig) ;
           }
      }
}
void	fc_conflict(int n, COMPLEX oldR, COMPLEX newR)
{
     COMPLEX	dR, sR ;
     int	h, k, l ;

     dR.re = newR.re - oldR.re ;
     dR.im = newR.im - oldR.im ;
     sR.re = newR.re + oldR.re ;
     sR.im = newR.im + oldR.im ;

     if (((dR.re * dR.re) + (dR.im * dR.im)) > 
         DEPS*DEPS*((sR.re * sR.re) + (sR.im * sR.im))) {

         fcalc_mismatch_count++ ;

         if ((verbose)||(fcalc_mismatch_count < MAX_MISMATCH)) {
              fetch_hkl(n, &h, &k, &l) ;
              fprintf(fp_log, 
     "FCALC conflict at h=%d, k=%d, l=%d, old=(%g,%g), new=(%g,%g)\n",
              h, k, l, oldR.re, oldR.im, newR.re, newR.im) ;
         }
     }
}

	/*********************************************
        Enforce half-ellipsoid rule;
	count points for which mask is set uniquely.
	*********************************************/

int	enforce_P1(char *mask) 
{
     int	N = 0 ;
     int	n, h, k, l ;

     for (n = 0; n < Nhkl; n++)
	  if (*(mask+n)) {
	       fetch_hkl(n, &h, &k, &l) ;
	       if (((h > 0)) ||
		   ((h == 0) && (k > 0)) ||
		   ((h == 0) && (k == 0) && (l >=0))) 
		     N++ ;
               else
		    *(mask+n) = 0 ;
          }

     return (N) ;
}
void	init_fo(real *array, real *sigma, char *mask)
{

     int	n ;

     for (n = 0; n < Nhkl; n++) {
          *(array+n) = 0. ;
          *(sigma+n) = 1. ;
          *(mask+n) = 0 ;
     }

     fobs_mismatch_count  = 0 ;	
     N.used = N.total = 0 ;
}

void	init_fc(COMPLEX *F, char *mask)
{
     int	n ;

     for (n = 0; n < Nhkl; n++) {
          *(mask+n) = 0 ;
          (F+n)->re = 0. ;
          (F+n)->im = 0. ;
     }

     fcalc_mismatch_count = 0 ;
     N.used = N.total = 0 ;
}

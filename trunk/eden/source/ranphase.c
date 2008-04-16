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

                                RANPHASE

  Title:        RANPHASE - Generate fcalc files with partially random phases.
  Author:       Hanna Szoke
  Date:         6/23/97
  Function:     Using as input an fobs file (Fob) and (optionally) a model fcalc
		file (R), this program generates complex Oran vectors with 
		semi-random phases, such that the vector composition 
				FRan = R + Oran 
		satisfies |FRan| ~ Fob.   The ~ sign indicates that we use
		sigmas (if available) to select a point on each Gaussian.
		Input fobs and fcalc files should be apodized and scaled.
	
*******************************************************************************/

#include "util.h"
#include "cellparams.h"
#include "symmetry.h"
#include "dims.h"

float	ran1() ;		/* packaged in util.c */

static	int	pseed ;			/* starting point for random numbers */
static	char	*mcentric ;
static	real	*mphases ;

static	real	generateRanPhase(int, real, real) ;
static	void	getFRan_noR(real *, real *, char *, COMPLEX *) ; 
static	void	getFRan_withR(COMPLEX *, real *, real *, char *, COMPLEX *) ;
static	void	quadsolve(real, real, real *, real *) ;

void	ranphase_main(int argc, char *argv[])
{
     int      	readfobs(char *, real *, char *, real *) ;
     int	checkSigmas(real *, char *) ;
     void	set_cen_phases(char *, real *) ;
     void	prepare_unique_hklmask(char **) ;
     int      	readfcalc(char *, COMPLEX *, char *) ;
     void	writefcalc(char *, COMPLEX *, char *) ;

     char	fo_filename[MAXSTRING] ;
     char	fc_filename[MAXSTRING] ;

     COMPLEX	*R=NULL, *FRan ;
     real	*Fob, *sigma ;
     char	*maskfo, *maskfc ;
     char	*unique_mask ;
     int	is_model = FALSE ;
     int	n, s = 1 ;

	/***********************************************************
	Check for right ballpark, etc.
	***********************************************************/

     sprintf(fc_filename, "none") ;

     if (argc > optind+5) 	/* too many arguments */
	  ballpark(caller) ; 
     else if (argc == optind+5) {
          sprintf(fo_filename, "%s", argv[optind+2]) ;
	  s = atoi(argv[optind+3]) ;
	  is_model = TRUE ;
	  sprintf(fc_filename, "%s", argv[optind+4]) ;
     }
     else  if (argc == optind+4) {
          sprintf(fo_filename, "%s", argv[optind+2]) ;
	  s = atoi(argv[optind+3]) ;
     }
     else {
	  prompt(
  "\nPlease enter fobs filename, a seed, & optionally, an fcalc filename: ") ;
          sscanf(terminp, "%s %i %s", fo_filename, &s, message) ;
	  if ((int)strlen(message) > 0) {
	       is_model = TRUE ;
	       sprintf(fc_filename, "%s", message) ;
          }

	  sprintf(message, 
		"Running %s on fobs file %s with seed %i and fcalc file %s.", 
	         caller, fo_filename, s, fc_filename) ;
	  printTwice(message) ;
     }

     hello(caller) ;
      
     if (s < 1) 
	 EdenError("Bad seed - stopping.") ;
					    
	/***********************************************************
        Just to ensure that we don't use very small seeds that have 
	the weird property that they occasionally return 1, add a 
	"random" number before complementing.
	***********************************************************/
							    
     pseed = -(s + 137) ;
							  
     sprintf(out_filename, "ran%02d.fcalc", s) ;

     readBasicInput() ;

	/***********************************************************
        Allocations and reads.
	***********************************************************/

     strcpy(message, "ranphase") ;

     Fob   = (real *) e_malloc(Nhkl*sizeof(real), message) ;
     sigma = (real *) e_malloc(Nhkl*sizeof(real), message) ;

     if (is_model) { 
         R = (COMPLEX *) e_malloc(Nhkl*sizeof(COMPLEX), message) ;
     }

     maskfc   = (char *) e_malloc(Nhkl*sizeof(char), message) ;
     maskfo   = (char *) e_malloc(Nhkl*sizeof(char), message) ;
     mcentric = (char *) e_malloc(Nhkl*sizeof(char), message) ;
     mphases = (real *) e_malloc(Nhkl*sizeof(real), message) ;
     
     FRan = (COMPLEX *) e_malloc(Nhkl*sizeof(COMPLEX), message) ;

     readfobs(fo_filename, Fob, maskfo, sigma) ;
     useSig = checkSigmas(sigma, maskfo) ;
     set_cen_phases(mcentric, mphases) ;

	/***********************************************************
	It's desirable that only the unique reflections in the fobs 
	file be used by ranphase.  
	***********************************************************/

     prepare_unique_hklmask(&unique_mask) ;
     for (n = 0; n < Nhkl; n++)
         *(maskfo+n) *= *(unique_mask+n) ;

     if (is_model)
	  readfcalc(fc_filename, R, maskfc) ;

     e_free(Nhkl*sizeof(char), maskfc) ;

	/***********************************************************
        Distinguish between cases with and without a model fcalc.
	***********************************************************/

     if (is_model) 
	 getFRan_withR(R, Fob, sigma, maskfo, FRan) ;
     else
	 getFRan_noR(Fob, sigma, maskfo, FRan) ;

     writefcalc(out_filename, FRan, maskfo) ;
}

real	generateRanPhase(int n, real phamin, real phamax)
{
     real	phase = 0 ;
     real	range ;

	/***********************************************************
	We distinguish between centrics and acentrics; for centrics, 
	we apply one of the permissible centric phases;
	for acentrics, we apply a random phase limited to the range.
	***********************************************************/

     range = phamax - phamin ;

     if (n > 0) {
        if (*(mcentric+n)) {
	   phase = *(mphases+n) * RTOD ;
   	   if (ran1(&pseed) > 0.5) phase += 180. ;
        }
        else  
	   phase = range * ran1(&pseed) + phamin ;
     }

     return (phase) ;
}

				/****************************************
				With no model, phases are genuine
				random (to within centric restrictions).
				The only trickiness involves ensuring 
				that the Gaussian perturbation to Fob 
				doesn't cause the amplitude to be < 0.
				****************************************/

void	 getFRan_noR(real *Fob, real *sigma, char *mask, COMPLEX *FRan) 
{

     COMPLEX	*FR ;
     real	Fphase ;
     real	phamin = 0, phamax = 360 ;
     int	n ;

     for (n = 0, FR = FRan; n < Nhkl; n++, FR++) {
          FR->re = 0 ;
          FR->im = 0 ;
     }

     for (n = 0, FR = FRan; n < Nhkl; n++, FR++, Fob++, sigma++, mask++) {

	  if ((*mask) && (*Fob > 0)) {

               if (n == 0) {		/* handle separately */
	            FR->re = *Fob ;
	            FR->im = 0 ;
               }

	       Fphase = generateRanPhase(n, phamin, phamax) ;

	       FR->re = *Fob * cos(DTOR * Fphase ) ;
	       FR->im = *Fob * sin(DTOR * Fphase ) ;
	  }

     }
}

				/****************************************
				With a model, phases are no longer
				truly random.  The phamin/phamax restrict
				the choice to a range around the direction
				of Rcalc.  Also, we must still ensure
				that the Gaussian perturbation to Fob 
				doesn't cause the amplitude to be < 0.
				****************************************/

void	 getFRan_withR(COMPLEX *Rcalc, real *Fob, real *sigma, 
		       char *mask, COMPLEX *FRan) 
{
     real	normGauss(int *) ;

     COMPLEX	*FR ;
     real	amp0 ;
     real	Famp ;
     real	Oamp, Ophase ;
     real	Ramp, Rphase ;
     real	ranphase, phamin, phamax ;
     real	bcoef, ccoef, discrim ;
     real	x1, x2 ;
     int	n, count_rp = 0 ;

     for (n = 0, FR = FRan; n < Nhkl; n++, FR++) {
          FR->re = 0 ;
          FR->im = 0 ;
     }

     amp0 = Rcalc->re ;

     for (n = 0, FR = FRan; 
	  n < Nhkl; n++, Rcalc++, FR++, Fob++, sigma++, mask++) {

	  if ((*mask) && (*Fob > 0)) {

	       Famp = (useSig) ? *Fob + *sigma * normGauss(&pseed) : *Fob ;

               if (n == 0) {		/* handle separately */
	            FR->re = Famp ;
	            FR->im = 0 ;
               }

               else if (Famp > 0) {

	            getAmpPhase(Rcalc, amp0, &Ramp, &Rphase) ;

	/***********************************************
	We need to find |O| from the complex vector sum:
	   (R.re + O.re)^2 + (R.im + O.im)^2 = Famp^2.  

	Expanding, we must solve the quadratic equation:

	|O|^2 + 2|O||R|cos(ranphase) + (|R|^2 - Famp^2) = 0.
   or:
	|O|^2 + |O|*bcoef + ccoef = 0, where

	ranphase is the >>relative<< phase of O wrt R.
	***********************************************/

                    count_rp++ ;
	            ccoef = Ramp*Ramp - Famp*Famp ;

	            if (Ramp <= Famp) {
		         phamax = 180. ;
		         phamin = -180. ;
		    }
	            else {
		         phamax = 180. + RTOD * asin(Famp/Ramp) ;
		         phamin = 180. - RTOD * asin(Famp/Ramp) ;
		    }

	            ranphase = generateRanPhase(n, phamin, phamax) ;
                    bcoef = 2 * Ramp * cos(DTOR * (ranphase)) ;
	            discrim = bcoef*bcoef - 4*ccoef ;

                    if (discrim < 0) /* Shouldn't get to this */
{
sprintf(message, "Famp = %g, Ramp = %g, rat = %g, phasmax/min = %g, %g",
Famp, Ramp, Famp/Ramp, phamin, phamax) ;
EdenWarning(message) ;
sprintf(message, "n = %d, ranphase=%g, bcoef=%g, discrim=%g", 
n, ranphase, bcoef, discrim) ;
EdenWarning(message) ;
EdenError("Quitting.") ;
}

	            quadsolve(bcoef, ccoef, &x1, &x2) ;

                    if ((x1 < 0) && (x2 < 0)) /* Shouldn't get to this */
{
sprintf(message, "Famp = %g, Ramp = %g, rat = %g, phasmax/min = %g, %g",
Famp, Ramp, Famp/Ramp, phamin, phamax) ;
EdenWarning(message) ;
sprintf(message, "n = %d, ranphase=%g, bcoef=%g, discrim=%g, x1 = %g, x2 = %g",
n, ranphase, bcoef, discrim, x1, x2) ;
EdenWarning(message) ;
EdenError("Quitting.") ;
}

	/***********************************************
        If Famp >= Ramp, there should be one positive 
	solution (which we select); otherwise, we select 
	one of the two positives by random throw. 
	***********************************************/

                    if (Famp >= Ramp) 
		         Oamp = (x1 > x2) ? x1 : x2 ;
	            else		
		         Oamp = (ran1(&pseed) < 0.5) ? x1 : x2 ;

                    if (Oamp < 0)	 /* Shouldn't get to this  */
{
sprintf(message, "Hit a negative amplitude! - n=%d, Famp=%g, Ramp=%g, Oamp=%g",
n, Famp, Ramp, Oamp) ;
EdenWarning(message) ;
EdenError("Quitting.") ;
}

                    Ophase = ranphase + Rphase ;
	            FR->re = Rcalc->re + Oamp * cos(DTOR * Ophase) ;
	            FR->im = Rcalc->im + Oamp * sin(DTOR * Ophase) ;

               	}
          }
     }
     sprintf(message, 
	     " # of random phases applied to SF's = %d\n", count_rp) ;
     printTwice(message) ;
     
     return ;
}

void	quadsolve(real bcoef, real ccoef, real *x1, real *x2) 
{
     int	sgb ;
     real	q ;

		/* Follow Numerical Recipes, p145, replacing a -> 1. */
		
     sgb = (bcoef >= 0) ? 1: -1 ;
     q = -0.5 * (bcoef + sgb * sqrt(bcoef*bcoef - 4 * ccoef)) ;
     *x1 = q ;
     *x2 = ccoef / q ;
}

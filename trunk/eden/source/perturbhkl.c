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

                                PERTURBHKL

  Title:        PERTURBHKL - Perturb fcalc file.
  Author:       Hanna Szoke
  Date:         6/23/97
  Function:     Using as input a model fcalc file, this program generates 
		complex vectors with Gaussian random perturbations to both
		amplitudes and phases.
	
*******************************************************************************/

#include "util.h"
#include "cellparams.h"
#include "symmetry.h"
#include "dims.h"

float	ran1() ;		/* packaged in util.c */

static	int	pseed ;			/* starting point for random numbers */
static	real	pert_frac ;

static	void	 getRpert(COMPLEX *, char *, char *, COMPLEX *) ; 

void	perturbhkl_main(int argc, char *argv[])
{
     void	set_centrics(char *) ;
     void	prepare_unique_hklmask(char **) ;
     int      	readfcalc(char *, COMPLEX *, char *) ;
     void	writefcalc(char *, COMPLEX *, char *) ;

     char	fc_filename[MAXSTRING] ;

     COMPLEX	*R, *Rpert ;
     char	*maskfc ;
     char	*mcentric ;
     char	*unique_mask ;
     int	n, s = 1 ;

	/***********************************************************
	Check for right ballpark, etc.
	***********************************************************/

     if (argc > optind+5) 	/* too many arguments */
	  ballpark(caller) ; 
     else if (argc > optind+3) {
          sprintf(fc_filename, "%s", argv[optind+2]) ;
          pert_frac = atof(argv[optind+3]) ;
          if (argc == optind+5)
	       s = atoi(argv[optind+4]) ;
     } 
     else {
	  prompt("\nPlease enter filename, pert_frac, & optionally, seed: ") ;
          sscanf(terminp, "%s %f %s", fc_filename, &t1, message) ;
          pert_frac = t1 ;
	  if ((int)strlen(message) > 0)
	       s = atoi(message) ;
	  sprintf(message, 
		"Running %s on file %s with pert_frac %f and seed %i", 
	         caller, fc_filename, pert_frac, s) ;
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
							  
     sprintf(out_filename, "pert%02d.fcalc", s) ;

     readBasicInput() ;

	/***********************************************************
        Allocations and reads.
	***********************************************************/

     strcpy(message, "perturbhkl") ;

     R = (COMPLEX *) e_malloc(Nhkl*sizeof(COMPLEX), message) ;
     Rpert = (COMPLEX *) e_malloc(Nhkl*sizeof(COMPLEX), message) ;

     maskfc   = (char *) e_malloc(Nhkl*sizeof(char), message) ;
     mcentric = (char *) e_malloc(Nhkl*sizeof(char), message) ;

     set_centrics(mcentric) ;

	/***********************************************************
	It's desirable that only the unique reflections in the fcalc 
	file be used by perturbhkl.  
	***********************************************************/

     readfcalc(fc_filename, R, maskfc) ;
     prepare_unique_hklmask(&unique_mask) ;

     for (n = 0; n < Nhkl; n++)
         *(maskfc+n) *= *(unique_mask+n) ;

     getRpert(R, maskfc, mcentric, Rpert) ;
     writefcalc(out_filename, Rpert, maskfc) ;
}

				/****************************************
				****************************************/

void	 getRpert(COMPLEX *Rcalc, char *mask, char *mcentric, COMPLEX *Rpert) 
{
     real	normGauss(int *) ;

     COMPLEX	*RP ;
     real	amp0 ;
     real	Ramp, Rphase ;
     real	newamp ;
     int	n ;
     int	p = 0 ;

     for (n = 0, RP = Rpert; n < Nhkl; n++, RP++) {
          RP->re = 0 ;
          RP->im = 0 ;
     }

     amp0 = Rcalc->re ;

     for (n = 0, RP = Rpert; 
	  n < Nhkl; n++, Rcalc++, RP++, mcentric++, mask++) {

	  if (*mask) {

	       getAmpPhase(Rcalc, amp0, &Ramp, &Rphase) ;

	       if (*mcentric) {
	            newamp = Ramp + normGauss(&pseed)*pert_frac*Ramp ;
		    RP->re = newamp*cos(DTOR*Rphase) ;
		    RP->im = newamp*sin(DTOR*Rphase) ;
               }
	       else {
	            RP->re = Rcalc->re + normGauss(&pseed)*pert_frac*Ramp ;
	            RP->im = Rcalc->im + normGauss(&pseed)*pert_frac*Ramp ;
               }
	       p++ ;

          }
     }
     /* Handle F(000) separately */

     Rpert->re = amp0 ;
     Rpert->im = 0 ;
     
     sprintf(message, "%d structure factors were perturbed", p) ;
     printTwice (message) ;
}

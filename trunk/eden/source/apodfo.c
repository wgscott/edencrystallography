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

                                APODFO

  Title:        Apodize a file of fobs data and prepare Wilson-style plots.
  Author:       Hanna Szoke
  Date:         10/05/00: scale_fo() generalized & moved to apodutil.

*******************************************************************************/

#include "util.h"
#include "cellparams.h"
#include "dims.h"	/* for eta and grid_spacing */
#include "apodize.h"

static	real	*Fap ;
static	real	*sigap ;
static	int	Ntot ;			/* # masked-in structure factors */
static  real	*Int ;			/* full set of (amp^2)	*/
static  real	*ksq ;			/* full set of (1/d^2)  */
static  real	*sigsqI ;
static  real	*sigsq ;
static  real	*sratio ;
static 	real	asig, bsig ;		/* coefs for apodization of sigmas */

static	void	analyze_sigsq() ;
static	int	apodfo_startup() ;
static	void	apply_apod_fo() ;	
static	void	prepare_fo() ;		/* make weighted F^2, sig^2 */

void	apodfo_main(int argc, char *argv[])
{
     void	writefobs(char *, real *, char *, real *) ;

	/* in apodutil.c:	*/

     void	allocate_Bin_arrays(int) ;
     void	apply_w0_correction() ;
     real	choose_delta(real, real, real, real, real, real) ;
     void	do_binning(int, real *, real *, real *) ;	
     real	get_delta(real) ;
     void	init_apod_arrays(int, real *, real *, real *) ;
     real	linfit(int, real *, real *, char) ;		
     void	make_wilson(char, real *) ;		
     void	make_selected_wilson(char *) ;		
     void	scale_fo(char *, char *) ;
     void	show_wilson_plots() ;

     char	name_fc[MAXSTRING], name_fo[MAXSTRING] ;
     real	delta1, delta2 ;
     real	dstdev1, dstdev2 ;
     real	B1, B2 ;
     char	corr ;			/* are we applying w0 correction? */

  /*******************************************************
  Check for right ballpark; identify fobs file name.
  *******************************************************/

     if (argc > optind+3) 	/* too many arguments */
	  ballpark(caller) ; 
     else if (argc == optind+3)
          sprintf(sf_filename, "%s", argv[optind+2]) ;
     else {
	  prompt("\nPlease enter the name of an fobs file: ") ;
          sscanf(terminp, "%s", sf_filename) ;
	  sprintf(message, "Running Apodfo on fobs file %s", sf_filename) ;
	  printTwice(message) ;
     }

     strip_path(sf_filename, pwd_name) ;

     hello(caller) ;

  /*******************************************************
  Read the input fobs file using standard read function.
  Set up arrays in (h,k,l) reciprocal lattice triplets,
  Prepare a list of non-zero intensities, sigma squared*
  intensities, and nusq's.  
  *******************************************************/

     Nbin = apodfo_startup() ; 
     
     init_apod_arrays(Ntot, Int, ksq, sigsqI) ;
     prepare_fo() ;

  /*******************************************************
  Bin the information, gathering logs of average binned 
  intensities & weights.  Determine coefficients (a, b) 
  for fitting sigsq as a linear function of Int.
  *******************************************************/

     allocate_Bin_arrays(Nbin) ;
     do_binning(Ntot, Int, ksq, sigsqI) ;

     if (useSig)
         analyze_sigsq() ;

  /*******************************************************
  Do a linear fit to establish slope and intercept.  
  Prepare for display the original binned log(I^2) and 
  the linearly fitted version.
  *******************************************************/

     printTwice("\n\tDoing calculation without correction for proteins ...") ;

     corr = FALSE ; 	/* no correction 1st time */

     B1 = linfit(Nbin, w_lnFsq, apodwt, corr) ; 
     delta1 = get_delta(B1) ; 
     make_wilson(corr, &dstdev1) ;

  /*******************************************************
  If useW0 flag is set,
  repeat calculation for a plot to which the universal W0 
  curve correction has been applied.  Display both sets of 
  plots and let user select the appropriate delta.
  *******************************************************/
   
     if (useW0) {
          printTwice("\n\tRedoing calculation with correction ...") ;

          corr = TRUE ; 	/* correction 2nd time */

          apply_w0_correction() ;
          B2 = linfit(Nbin, w_lnFsq, apodwt, corr) ; 
          delta2 = get_delta(B2) ; 

          make_wilson(corr, &dstdev2) ;
          show_wilson_plots() ;
          delta = choose_delta(delta1, dstdev1, B1, delta2, dstdev2, B2) ;
     }
     else { 
          show_wilson_plots() ;
	  delta = delta1 ;
     }

  /*******************************************************
  Replace Fobs (aka Fap) and sigap by apodized version, 
  insofar as delta > 0.
  *******************************************************/

     if (delta > 0) {
          apply_apod_fo() ;

          extend_filename(pwd_name, out_filename, "_apo") ;

          sprintf(message, "\nWriting %s", out_filename) ;
          printTwice(message) ;

          writefobs(out_filename, Fap, maskap, sigap) ;

          init_apod_arrays(Ntot, Int, ksq, sigsqI) ;
          prepare_fo() ;

          do_binning(Ntot, Int, ksq, sigsqI) ;
     }
     sprintf(name_fo, "%s_wil", pwd_name) ;
     make_selected_wilson(name_fo) ;

  /*******************************************************
  Find out whether user wants to scale fo to fc. 
  *******************************************************/

     prompt("Scale? - y or n: ") ;

     if (terminp[0] == 'n')
	  return ;		/* no scaling  */

    prompt(
    "Enter name of file containing fc Wilson data\nfrom the end of your Apodfc run: ")  ;

     sscanf(terminp, "%s", name_fc) ;
     fprintf(fp_log, "Wilson data file name is %s\n", name_fc) ;

     scale_fo(name_fo, name_fc) ;
}

  /*******************************************************
  Prepare a list of non-zero squared amplitudes (I),
  (sigma*I) squared values, and nusq values.
  We also need sigsq for apodization of sigmas.
  *******************************************************/

void	prepare_fo()
{
     real	amp ;
     real	sig ;
     int	n ;
     int        m ;
     int        h, k, l ;

     /* Handle (h,k,l) triplets except for (000) */

     for (n = 1, m = 0; n < Nhkl; n++) {

          if (*(maskap+n)) { 

               amp = *(Fap+n) ;
	       sig = *(sigap+n) ;
               fetch_hkl(n, &h, &k, &l) ;

               *(Int+m) = amp * amp ;
	       *(ksq+m) = nusq(h, k, l) ;
	       *(sigsq+m) = sig * sig ;
	       *(sigsqI+m) = 4. * amp * amp * sig * sig ;
	       m++ ;
	  } 

     }			
}

  /*******************************************************
  Apply exponential factor to (h,k,l) triplets except for 
  (000). using the standard  PISQ * delta * eta * dr * dr,
  Round off very low amplitudes. Sigmas are apodized using 
  a, b and sratio array.  See analyze_sigsq().
  *******************************************************/

#define F_EPSI	1.e-7

void	apply_apod_fo() 	
{	  
     real	efac ;
     int	n ;
     int        h, k, l ;
     real	Fmin ;
     real	thisF ;

     Fmin = *Fap * F_EPSI ;

     for (n = 1; n < Nhkl; n++) 
	     
          if (*(maskap+n)) {

               fetch_hkl(n, &h, &k, &l) ;

	       efac = exp(-(delta * delfac * nusq(h, k, l))) ;

               *(Fap+n) *= efac ;

	       if (*(Fap+n) < Fmin) 
		    *(Fap+n) = 0 ;

               if (useSig) {
		   thisF = *(Fap+n) ;
	           *(sigap+n) = sqrt(*(sratio+n) * (asig + bsig*thisF*thisF)) ;
	       } 
	  }
}

  /*******************************************************
  Read the input fobs file using standard read function.
  Set up arrays in (h,k,l) reciprocal lattice triplets,
  *******************************************************/

int	apodfo_startup() 
{
     void	a_readInput() ;		/* read input parameters */
     void	setDsqLimit(char *, int) ;
     int	readfobs(char *, real *, char *, real *) ;	
     int	checkSigmas(real *, char *) ;	
     void	checkF000(real *, real *, char *) ;

     int	n ;

     a_readInput() ;

     strcpy(message, "apodfo_startup") ;

     maskap = (char *)  e_malloc(Nhkl*sizeof(char),   message) ; 
     Fap   = (real *) e_malloc(Nhkl*sizeof(real), message) ;
     sigap = (real *) e_malloc(Nhkl*sizeof(real), message) ;

     Ntot = readfobs(sf_filename, Fap, maskap, sigap) ;
     checkF000(Fap, sigap, maskap) ;
     useSig = checkSigmas(sigap, maskap) ;

     strcpy(message, "apodfo_startup") ;

     if (useSig)
        sratio   = (real *) e_malloc(Nhkl*sizeof(real), message) ; 

     Int = (real *) e_malloc(Ntot*sizeof(real), message) ; 
     sigsqI = (real *) e_malloc(Ntot*sizeof(real), message) ; 
     ksq   = (real *) e_malloc(Ntot*sizeof(real), message) ; 
     sigsq = (real *) e_malloc(Ntot*sizeof(real), message) ; 

     if (useSig) 
       fprintf(fp_log, "Sigmas will be used for weighting\n") ;
     else
       fprintf(fp_log, "Sigmas will not be used for weighting\n") ;

     setDsqLimit(maskap, Nhkl) ;

     n = (dsq_limit + 0.5*binwidth)/ binwidth + 1 ;

     sprintf(message, "# of bins set to %d, average bin size is %g.", 
		n, binwidth) ;
     printTwice(message) ;

     return (n) ;
}

  /*******************************************************
  Determine coefficients (a, b) for fitting sigsq as a
  linear function of Int -- i.e. minimize
	Sum [sigsq - (a+b*Int)]^2.
  Letting f to be the above expression, we set partial 
  derivatives of f wrt a and b to 0, then solve.

  We set up sratio to be an array of sigsq / (a + b*Int)
  *******************************************************/

void	analyze_sigsq() 
{
    int	   m, n ;
    real thisI, thissq ;
    real SInt = 0, SIntsq = 0, SsigsqInt = 0, Ssigsq = 0 ;
    real denom, numa, numb ;

    for (m = 0; m < Ntot; m++) {
	 thisI = *(Int+m) ;
	 thissq = *(sigsq+m) ;
	 SInt += thisI ;
	 SIntsq += thisI * thisI ;
	 SsigsqInt += thisI * thissq ;
	 Ssigsq += thissq ;
    }
    denom = SInt * SInt - Ntot * SIntsq ;
    numa = SInt * SsigsqInt - SIntsq * Ssigsq ;
    numb = SInt * Ssigsq - SsigsqInt * Ntot ;

    asig = numa / denom ;
    bsig = numb / denom ;

    sprintf(message, 
             "Sigma apodization coefficients are %g and %g", asig, bsig);
    printTwice(message) ;

    if (asig < 0) {
	 printTwice("Your sigmas are very irregular;") ;
	 printTwice("Please check the fobs input or turn off the useSig flag!") ;
	 EdenError("Quitting...") ;
    }

    if (bsig < 0)  {
	printTwice(
    "Warning: your sigmas are not increasing with reflection intensities;") ;
        printTwice(message) ;
        sprintf(message, "Resetting bsig from %g to 0", bsig) ;
        printTwice(message) ;
        bsig = 0 ;
    }

    for (m = 0, n = 1; n < Nhkl; n++) {
	 if (*(maskap+n)) {
	      *(sratio+n) = *(sigsq+m) / (asig + bsig * *(Int+m)) ; 
	      m++ ;
         }
         else
	      *(sratio+n) = 1. ;
     }
}

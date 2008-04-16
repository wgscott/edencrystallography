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

                                APODFC

  Title:        Apodize a file of fcalc data and prepare Wilson-style plots.
  Author:       Hanna Szoke
  Date:         10/05/00: Minor changes to bring into line with apodfo.c.

*******************************************************************************/

#include "util.h"
#include "cellparams.h"
#include "dims.h"	/* for eta and grid_spacing */
#include "apodize.h"

static	COMPLEX	*Rap ;

static	int	Ntot ;			/* # masked-in structure factors */
static  real	*Int ;			/* full set of (amp^2)	*/
static  real	*ksq ;			/* full set of (1/d^2)  */
static  real	*sigsqI ;

static	int	apodfc_startup() ;
static	void	apply_apod_fc() ;	/* use the delta that was calculated  */
static	void	prepare_fc() ;		/* make weighted R^2 */

void
apodfc_main(int argc, char *argv[])
{
	/* in apodutil.c:       */

     void	allocate_Bin_arrays(int) ;
     void	apply_w0_correction() ;
     real	choose_delta(real, real, real, real, real, real) ;
     void	do_binning(int, real *, real *, real *) ;	
     real	get_delta(real) ;
     void	init_apod_arrays(int, real *, real *, real *) ;
     real	linfit(int, real *, real *, char) ;		
     void	make_wilson(char, real *) ;
     void	make_selected_wilson(char *) ;
     void	show_wilson_plots() ;
     void	writefcalc(char *, COMPLEX *, char *) ;	

     real	delta1, delta2 ;
     real	dstdev1, dstdev2 ;
     real	B1, B2 ;
     char	corr ;			/* are we applying w0 correction? */

  /*******************************************************
  Check for right ballpark; identify Fcalc file name.
  *******************************************************/

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

     strip_path(sf_filename, pwd_name) ;

     hello(caller) ;

  /*******************************************************
  Read the input fcalc file using standard read function.
  Set up arrays in (h,k,l) reciprocal lattice triplets,
  Prepare a list of non-zero intensities, dummy sigma
  squared*intensities, and nusq's.  
  *******************************************************/

     Nbin = apodfc_startup() ;

     init_apod_arrays(Ntot, Int, ksq, sigsqI) ;
     prepare_fc() ;

  /*******************************************************
  Bin the information, gathering logs of average binned 
  intensities & weights, 
  *******************************************************/

     allocate_Bin_arrays(Nbin) ;
     do_binning(Ntot, Int, ksq, sigsqI) ;

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

          corr = TRUE ; 	/* no correction 1st time */

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

  /******************************************************************
  Replace Fcalc (aka Rap) by apodized version, insofar as delta > 0.
  ******************************************************************/

     if (delta > 0) {
          apply_apod_fc() ;

          extend_filename(pwd_name, out_filename, "_apo") ;

          sprintf(message, "\nWriting %s", out_filename) ;
          printTwice(message) ;

          writefcalc(out_filename, Rap, maskap) ;

	  init_apod_arrays(Ntot, Int, ksq, sigsqI) ;
	  prepare_fc() ;

          do_binning(Ntot, Int, ksq, sigsqI) ;
     }
     sprintf(out_filename, "%s_wil", pwd_name) ;
     make_selected_wilson(out_filename) ;

     printTwice("\nPLEASE NOTE:") ;
     sprintf(message, " The Wilson file '%s' should be used for scaling your fobs.",
     out_filename) ;
     printTwice(message) ;

     return ;
}

  /*******************************************************
  Prepare a list of non-zero squared amplitudes (I) and 
  nusq values.  Handle (h,k,l) triplets except for (000).  
  Use dummy sigsqI, for compatibility with fobs.
  *******************************************************/

void prepare_fc()
{
     real	amp, amp0 ;
     real	phase ;
     int	n ;
     int        m ;
     int        h, k, l ;

     amp0 = Rap->re ;

     for (n = 1, m = 0; n < Nhkl; n++) {
	     
          if (*(maskap+n)) {

               fetch_hkl(n, &h, &k, &l) ;
	       getAmpPhase(Rap+n, amp0, &amp, &phase) ;

               *(Int+m) = amp * amp ;
	       *(ksq+m) = nusq(h, k, l) ;
	       *(sigsqI+m) = 4. * amp * amp ;	
	       m++ ;

	  }
     }			

     return ;
}

  /*******************************************************
  Apply exponential factor to (h,k,l) triplets except for 
  (000). using the standard  PISQ * delta * eta * dr * dr,
  *******************************************************/

void	apply_apod_fc() 	
{	  
     real	amp, amp0 ;
     real	phase ;
     int	n ;
     int        h, k, l ;

     amp0 = Rap->re ;

     for (n = 1; n < Nhkl; n++) {
	     
          if (*(maskap+n)) {

               fetch_hkl(n, &h, &k, &l) ;
	       getAmpPhase(Rap+n, amp0, &amp, &phase) ;

               amp *= exp(-(delta * delfac * nusq(h, k, l))) ;
               (Rap+n)->re = amp * cos(DTOR*phase) ;
               (Rap+n)->im = amp * sin(DTOR*phase)  ;
	  }
     }			

     return ;
}

  /*******************************************************
  Read the input fcalc file using standard read function.
  Set up arrays in (h,k,l) reciprocal lattice triplets,
  *******************************************************/

int apodfc_startup()
{
     void	a_readInput() ;		/* read input parameters */
     void	setDsqLimit(char *, int) ;
     int	readfcalc(char *, COMPLEX *, char *) ;	
     int	n ;

     a_readInput() ;

     strcpy(message, "apodfc_startup") ;

     maskap = (char *) e_malloc(Nhkl*sizeof(char), message) ; 
     Rap = (COMPLEX *) e_malloc(Nhkl*sizeof(COMPLEX), message) ;

     Ntot = readfcalc(sf_filename, Rap, maskap) ;
     if (Ntot == 0) 
	  EdenError("Cannot apodize an empty file!") ;

     strcpy(message, "apodfc_startup") ;

     Int = (real *) e_malloc(Ntot*sizeof(real), message) ; 
     sigsqI = (real *) e_malloc(Ntot*sizeof(real), message) ; 
     ksq   = (real *) e_malloc(Ntot*sizeof(real), message) ; 

     setDsqLimit(maskap, Nhkl) ;

     n = (dsq_limit + 0.5*binwidth)/ binwidth + 1 ;

     sprintf(message, "# of bins set to %d, average bin size is %g.", 
		n, binwidth) ;
     printTwice(message) ;

     return (n) ;
}

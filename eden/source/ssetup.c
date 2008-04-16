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

                                SOLVE SETUP

  Title:        Setup Functions for Solve
  Author:       Hanna Szoke
  Date:         Separated from Solve 10/18/01

	ssetup.c works in conjunction with the following source files:

		hklutil.c, hklread.c and hklinit.c handle (hkl) space work;
		init.c and basics.c handle initialization procedures;
		qshare.c handles general functions related to Np space;      
		crysutil.c  handles crystal symmetry;
		util.c has low-level worker routines.

*******************************************************************************/

#include "eden.h"		/* util.h, cellparams.h, symmetry.h, dims.h */
#include "eden_Np.h"
#include "mir.h"
#include "cost.h"

char	*intermask ;	/* intersection of all Fobs masks (MIR or MAD) */
real	abs_scale_m[MAXDER] ;
char	cost_filename[MAXSTRING] ;			/* in cost.h */
real	AveWt ;
real	norm_fac ;	/* 1/sqrt(AveWt), in 'cost.h' */
char	*maskfc ;	/* mask identifying fc input values, in eden.h  */

real	*rexpfac[MAXRES] ;		/* variables in 'mir.h' */
COMPLEX	*rexpfac_bcc[MAXRES] ;

static	void	calc_anomSF(COMPLEX *, real, real, real) ;
static	void	prepare_intersect_mask() ;
static	void	report_mode() ;
static	void	report_memory() ;
static	int	setup_fobs() ;
static	void	setup_fcalc() ;
static	void	allocNhklArrays() ;
static	void	addin_native_calc(COMPLEX *, COMPLEX *, char *, char *) ; 
static	void	compareMasks(char *, char *, int) ; 

void	solver_setup(char *problemName)
{
     void	allocNpArrays() ;		/* in 'qshare.c' */
     void	initNparrays() ;		/* in 'qshare.c' */
     void       echo() ;			/* in 'init.c' */
     void       echo_ignored() ;		/* in 'init.c' */
     void       readInput() ;			/* in 'init.c' */
     void       readInput_Mir() ;		/* in 'init.c' */
     void       readAnomInput() ;		/* in 'init.c' */
     void	prepare_singlets(char *) ;	/* in 'cost.c' */
     void	prepare_triplets(char *) ;	/* in 'cost.c' */

     int	Neq ;		/* # equations = # non-zero F values*/
     int	k ;

	/* set up cost file */

     hello(caller) ;

     if (verbose)
	  costfile_option = TRUE ;

     if (costfile_option) {
	  sprintf(cost_filename, "%s.cost", problemName) ;

          if ((fp_cost = fopen(cost_filename,"w")) == NULL) 
          {
               sprintf(message, "Cannot open %s", cost_filename);
               EdenError(message) ;
          }
     }

  /************************************************************
   Fetch input conditions, report.                         
  ************************************************************/
 
     readInput() ;
     readInput_Mir() ;
     readAnomInput() ;

     report_mode() ;
     echo() ;
     echo_ignored() ;

  /************************************************************
   Do fobs part of the Nhkl-space setups.  Then allocate and 
   initialize arrays for real-space arrays.  Finally, generate
   or read fcalc files.
  ************************************************************/

     report_memory() ;
     Neq = setup_fobs() ;

     allocNpArrays() ;
     initNparrays() ;

     setup_fcalc() ;

     if (detwin)
          twinned_fc = (COMPLEX *) e_malloc(Nhkl*sizeof(COMPLEX), "twinning") ;

     for (k = 0; k < Nconstraints; k++) 
	
          switch (con_type[k]) {

	  case SINGLET :
	       prepare_singlets(ta_filename[k]) ;
               break ;
     
	  case TRIPLET :
	       prepare_triplets(ta_filename[k]) ;
               break ;
     
	  default :
	       break ;
          }

     sprintf(message,
           "\n\t\tThis problem has %d equations, %d unknowns\n", Neq, 
		Npextended) ;
     printTwice(message) ;
}

int	setup_fobs()
{
	/* functions in 'hklinit' */

     void	SigmaRange(real *, real *, real *, char *) ;
     int	checkSigmas(real *, char *) ;
     void       setupWeights() ;
     real	getNormFac() ;

	/* functions in 'hklutil.c' or hklread.c */

     void	setDsqLimit() ;
     void	checkF000(real *, real *, char *) ;
     int	readfobs(char *, real *, char *, real *) ;
     void	scalefobs(real *, real) ;

     /* local variables */

     int	Neq = 0 ;		/* # equations = # non-zero F values*/
     int	m, md, off ;
     real	sigmin, sigmax ;	/* unused here */
     real	*wt ;
     real	*sarray ;
     char	*mask ;
     int	anom_flag_hold ;

  /************************************************************
   Set up arrays in (h,k,l) --  reciprocal lattice triplets   
  ************************************************************/

     printTwice("Setting up initial arrays ...") ;

     allocNhklArrays() ;

  /************************************************************
  Read diffraction pattern, Full, for the full molecule.  
  ************************************************************/

     printTwice("\nReading fobs file ") ;
     printTwice("... for native") ;
     Neq = readfobs(fo_filename, Full, maskfo, sigma) ;
     if (relwt_der[0] != 0) { 
          checkF000(Full, sigma, maskfo) ;
          scalefobs(Full, fscale) ;
          scalefobs(sigma, fscale) ;
     }
     F000 =  *Full ;

  /************************************************************
   Deal with MIR and MAD (which are indistinguishable to Solve)
  ************************************************************/

     anom_flag_hold = anom_flag ;

     for (m = 0, md = 1; m < Nder; m++, md++) {
	  off = md*Nhkl ;
          sprintf(message, "\n...for derivative # %1d", md) ;
          printTwice(message) ;
	  anom_flag = anom_der_flag[md] ;
          Neq += readfobs(fo_der_fn[m], Full+off, maskfo+off, sigma+off) ;
          if (relwt_der[md] != 0) {
               checkF000(Full+off, sigma+off, maskfo+off) ;
               scalefobs(Full+off, fscale_der[md]) ;
               scalefobs(sigma+off, fscale_der[md]) ;
          }
     }
     anom_flag = anom_flag_hold ;
     abs_scale_m[0] = fscale ;

     for (m = 0, md = 1; m < Nder; m++, md++) 
	abs_scale_m[md] = fscale_der[md] ;

     setDsqLimit(maskfo, NM*Nhkl) ;

  /************************************************************
  Check sigma range and reset flag if there were missing sigmas 
  in input.  Set criterion for discrepancy-criterion-satisfied 
  and normalize the weights.
  ************************************************************/

     for (m = 0; m < NM; m++) 
        if (relwt_der[m] > 0)
             useSig = checkSigmas(sigma + m*Nhkl, maskfo + m*Nhkl);

     if (useSig) {
          printTwice("\n\tSigma ranges") ;
	  if (relwt_der[0] > 0) {
               printTwice("\nFor native:") ;
               SigmaRange(&sigmin, &sigmax, sigma, maskfo) ; 
	  }

          for (m = 0, md = 1; m < Nder; m++, md++) {
	      if (relwt_der[md] > 0) {
	           sprintf(message, "For der %1d:", md) ;
	           printTwice(message) ;
	           SigmaRange(&sigmin, &sigmax, 
	     		     	sigma + md*Nhkl, maskfo + md*Nhkl);
	      }
     	  }
     }

  /************************************************************
  Handle weights, including normalization.
  ************************************************************/

     setupWeights() ;
     AveWt = getNormFac() ;
     norm_fac = 1. / sqrt(AveWt) ;

     for (m = 0, wt = hkl_weight, sarray = sigma, mask = maskfo; m < NM*Nhkl; 
	  m++, wt++, sarray++, mask++)

          if (*sarray > 0)
	       *wt *= *mask * norm_fac / *sarray ;

     fprintf(fp_log, "Hkl weight normalization factor = %g\n", norm_fac) ;

     return (Neq) ;
}

void	setup_fcalc()
{
     void	convolve(real *, real *, COMPLEX *, COMPLEX *) ; 
							/* in fftshare.c */
     void	discards(COMPLEX *, char *, char *) ;	/* in hklutil.c */
     real     *setupexpfac(real) ;			/* in hklinit.c */
     COMPLEX    *setupexpfac_bcc(real) ;		/* in hklinit.c */
     void	describeFR(real *, char *, char *) ;	/* in hklinit.c */
     int	readfcalc(char *, COMPLEX *, char *) ;	/* in hklread.c */
     char	*get_P1mask() ;				/* in hklutil.c */
     void	hrPrepareSFT() ;			/* in highres.c */
     void	hrDoSFT(real *, COMPLEX *) ;		/* in highres.c */

     int	m, md, off, r ;
     real	thisdelfac ;
     int	anom_flag_hold ;

  /************************************************************
   Set up arrays in (h,k,l) --  reciprocal lattice triplets   
  ************************************************************/

     thisdelfac = delfac * resrat[0] * resrat[0] ;
     expfac = rexpfac[0] = setupexpfac(thisdelfac) ; 
     if (grid_type == BODY_CENTERED)
	  expfac_bcc = rexpfac_bcc[0] = setupexpfac_bcc(thisdelfac) ;

     if (HighRes)
          hrPrepareSFT() ;

  /************************************************************
  Read a diffraction pattern, R0, corresponding to the known 
  part of the molecule or generate it from knownp; this will be 
  added into an iterated R within the loop. Enforce consistent 
  masking.  Provide some diagnostics on Full and R0. 
  ************************************************************/

     printTwice("\nGenerating fcalc file.") ;

     maskfc = get_P1mask() ;
     convolve(knownp, expfac, expfac_bcc, R0) ;

     if (HighRes) 
          hrDoSFT(knownp + Nptotal, R0) ;

     compareMasks(maskfo, maskfc, Nhkl) ;

  /************************************************************
   Deal with MIR and MAD (which are indistinguishable to Solve)
   Add the native model to the heavy atom fc data (MIR or MAD).
  ************************************************************/

     anom_flag_hold = anom_flag ;

     for (m = 0, md = 1; m < Nder; m++, md++) {
	  off = md*Nhkl ;
	  anom_flag = anom_der_flag[md] ;
          sprintf(message, "\n...for derivative # %1d", md) ;
          printTwice(message) ;
          readfcalc(fc_heavy_fn[m], R0+off, maskfc+off) ;
	  if (anom_der_flag[md])
	       calc_anomSF(R0+off, fprime[md], f2prime[md], Zheavy[m]) ;

          addin_native_calc(R0, R0+off, maskfc, maskfc+off) ;

          compareMasks(maskfo+off, maskfc+off, Nhkl) ;
     }
     anom_flag = anom_flag_hold ;

  /************************************************************
   Deal with assorted intrinsic resolutions of the derivatives;
  ************************************************************/

     for ( r = 1; r < Nres; r++) {
          thisdelfac = delfac * resrat[r] * resrat[r] ;
	  if (resrat[r] != resrat[0]) {
               rexpfac[r] = setupexpfac(thisdelfac) ; 
               if (grid_type == BODY_CENTERED)
	            rexpfac_bcc[r] = setupexpfac_bcc(thisdelfac) ;
          }
	  else {
	       rexpfac[r] = rexpfac[0] ;
               if (grid_type == BODY_CENTERED)
	            rexpfac_bcc[r] = rexpfac_bcc[r] = rexpfac_bcc[0] ;
          }

     }

     if (relwt_der[0] > 0)
	  describeFR(Full, maskfo, maskfc) ;

     for (m = 0, md = 1; m < Nder; m++, md++) {
	 fprintf(fp_log, "\nFor der %1d: ", md) ;
         if (relwt_der[md] > 0)
	      describeFR(Full + md*Nhkl, maskfo + md*Nhkl, maskfc + md*Nhkl);
     }
     if (Nder > 0)
	  prepare_intersect_mask() ;

     if ((relwt_der[0] > 0) && (very_verbose))
	  discards(R0, maskfo, maskfc) ;
}

void allocNhklArrays()
{
	/* functions in 'fftshare.c' or 'cost.c' */

     void	initfft() ; 
     void	initfft_aux() ; 

     int	n ;

  /************************************************************
   Arrays in Nhkl*NM (Native + MIR); Oar is essentially local 
   to cost.c and concosts.c.
   The masks are initialized (usually unneeded) because there
   are coniditions under which nothing is actually read...
  ************************************************************/

     strcpy(message, "allocNhklArrays" );

     maskfc = (char *) e_malloc(Nhkl*NM*sizeof(char), message) ;
     maskfo = (char *) e_malloc(Nhkl*NM*sizeof(char), message) ;

     for (n = 0; n < NM*Nhkl; n++) {
	  *(maskfo+n) = 0 ;
	  *(maskfc+n) = 0 ;
     }
     
     Full  = (real *)e_malloc(Nhkl*NM*sizeof(real), message) ;
     sigma = (real *)e_malloc(Nhkl*NM*sizeof(real), message) ;
     
     R =  (COMPLEX *) e_malloc(Nhkl*NM*sizeof(COMPLEX), message) ;
     R0 = (COMPLEX *) e_malloc(Nhkl*NM*sizeof(COMPLEX), message) ;
     Oar =(COMPLEX *) e_malloc(Nhkl*NM*sizeof(COMPLEX), message) ; 
     
  /************************************************************
  Set up arrays for the FFT package and initialize it.    
  ************************************************************/

     initfft() ;
     initfft_aux() ;
}

void	addin_native_calc(COMPLEX *Rnat, COMPLEX *RD, char *mnat, char *mder)
{
     int        n ;

     for (n = 0; n < Nhkl; n++, Rnat++, RD++, mnat++, mder++) {
          RD->re += Rnat->re ;
          RD->im += Rnat->im ;
	  *mder *= *mnat ;
     }
}

void	compareMasks(char *mfo, char *mfc, int Nitems) 
{
     int	n;
     int	throwoutNc = 0 ;

  /************************************************************
  Determine which (hkl) entries to use in the cost function;
  this is based on maskfo only!  There is no more real use for 
  maskfc, but it is used in a dummy capacity in the loop.
  ************************************************************/


     for (n = 0; n < Nitems; n++, mfo++, mfc++) 
	  throwoutNc += (!(*mfo) && (*mfc)) ;

     sprintf(message, "Fc reflections thrown out: %d\n", throwoutNc) ;
     printTwice(message) ;
}

void report_mode()      /* Summarize Eden mode of operation.         */	
{
     int	k ;

     printTwice(
	"\n\t******************************************************\n") ;

     if (Nder > 0) {
          if (Nder == 1) 
	       sprintf(message, "\tThis is an MIR run with one derivative.") ;
          else
	       sprintf(message, "\tThis is an MIR run with %d derivatives.", 
		    Nder) ;
          printTwice(message) ;
     }
     else {
	  sprintf(message, "\tThis is an Eden run in %s mode.", mode_info) ;
          printTwice(message) ;
     }
     if (detwin) {
          if (t_type == 'A')
               sprintf(message, 
          "\tAmplitude detwinning is enabled with fraction %g.",
                  t_frac) ;
          else
               sprintf(message, 
          "\tIntensity detwinning is enabled with fraction %g.",
                  t_frac) ;
          printTwice(message) ;
     }
     if (HighRes) {
          sprintf(message, "\tHigh Resolution processing is in effect.") ;
          printTwice(message) ;
     }

     if (Nconstraints == 0) {
	  sprintf(message, 
	       "\tThere are no Np constraints in the cost function.") ;
          printTwice(message) ;
     }
     else if (Nconstraints == 1) {
	  sprintf(message, 
	       "\tThere is 1 Np constraint in the cost function:") ;
          printTwice(message) ;
     }
     else {
	  sprintf(message, 
	       "\tThere are %d Np constraints in the cost function:",
	  Nconstraints) ;
          printTwice(message) ;
     }
     for (k = 0; k < Nconstraints; k++) 
	
          switch (con_type[k]) {

	  case TARGET :
	       sprintf(message, 
	            "\t   A %starget with relative weight %g.", target_type[k],
		    relwt_con[k]) ;
               printTwice(message) ;
	       break ;

	  case PHASE_EXT :
	       sprintf(message, 
	            "\t   A phase extension target at a resolution of %g.", 
		    phase_ext_res) ;
               printTwice(message) ;
	       sprintf(message, 
		    "\t   with relative weight %g", relwt_con[k]);
               printTwice(message) ;
	       break ;

	  case CS :
	       sprintf(message, 
                    "\t   Crystal symmetry enforced with relative weight %g.", 
	            relwt_con[k]) ;
               printTwice(message) ;
	       break ;

	  case SAYRE :
	       sprintf(message, 
		    "\t   A high resolution (Sayre's eq.) constraint") ;
               printTwice(message) ;
	       sprintf(message, 
		    "\t   with relative weight %g", relwt_con[k]);
               printTwice(message) ;
               sprintf(message, 
	            "\t   and an optimal cost of %g.", cost_addend[k]) ;
               printTwice(message) ;
	       break ;

	  case SINGLET:
	       sprintf(message, 
	            "\t   Singlet phases with relative weight %g.", 
		    relwt_con[k]) ;
               printTwice(message) ;
               break ;
     
	  case TRIPLET:
	       sprintf(message, 
	            "\t   Triplet phases with relative weight %g.", 
		    relwt_con[k]) ;
               printTwice(message) ;
               break ;
     
	  default :
	       break ;

          }

     printTwice(
	  "\n\t******************************************************\n") ;
}

void	report_memory()
{

     real	mem_Npspace, mem_recspace, mem_tot ;
     int	Nptot_gridPt ;
     int	k, Nex ;

	/*******************************************************
	In Np space, there are 6 Eden real arrays (snump, 
	totp, knownp, nump, minp and maxp) plus 1 Eden int 
	array (typep) and 4 CCG real arrays (xn, gn, go and d),
	each with Nptotal elements.  For each target, there are 
	another 2 real arrays (target, weight) disregarding
	exotic constraints.  There is also a complex FFT array 
	(temp) in Np, NOT Nptotal!
	*******************************************************/
     
     Nptot_gridPt = (10 + 2*Ntargets) * sizeof(real) + 1 * sizeof(int) ; 

     mem_Npspace = Nptot_gridPt * Nptotal + sizeof(COMPLEX) * Np ;

	/*******************************************************
	In reciprocal space, there are 3 Eden complex arrays in
	NM*Nhkl (R, R0 and Oar) and 1 or 2 complex arrays in 
	Nhkl (fftarray, expfac_bcc).  If detwinning is enabled,
	there is another complex array in Nhkl.
	There are 3 real arrays
	in NM*Nhkl (Full, sigma and hkl_weight ) plus 1 real 
	arrays in Nhkl (expfac).  There are 2 byte arrays in 
	NM*Nhkl (maskfc and maskfo) and 2 in Nhkl (intermask and
	unique_hkl).
	We ignore arrays that are 2-dimensional.
	*******************************************************/

     Nex = (grid_type == SIMPLE) ? 0 : 1 ;
     Nex += (detwin) ? 1 : 0 ;

     mem_recspace = 3 * NM * Nhkl * sizeof(COMPLEX) +
		    (1+Nex) * Nhkl * sizeof(COMPLEX) +
		    3 * NM * Nhkl * sizeof(real)  +
		    1 * Nhkl * sizeof(real) +
		    2 * NM * Nhkl * sizeof(char) +
		    2 * Nhkl * sizeof(char) ;

     for (k = 0; k < Nconstraints; k++) 
	
          switch (con_type[k]) {

	  case SINGLET :
	  case TRIPLET :
	       mem_recspace += Nhkl ;
               break ;
     
	  default :
	       break ;
          }

     mem_Npspace  *= BtoMB ;
     mem_recspace *= BtoMB ;
     mem_tot = mem_Npspace + mem_recspace ;

     printTwice("Approximate (underestimated) memory requirements in Mbytes:") ;
     sprintf(message, 
     "Physical space: %4.3g, Reciprocal space: %4.3g, Total: %4.3g.\n",
     mem_Npspace, mem_recspace, mem_tot) ;
     printTwice(message) ;
}

void	calc_anomSF(COMPLEX *array, real fp, real fpp, real Zh)
{
	  /******************************************************
	  Complex multiply hydrogen-like strucure fators by
	  (Z + f') + if''.
	  ******************************************************/

     real	afac ;
     COMPLEX	old ;
     int	n ;

     afac = Zh + fp ;

     for (n = 0; n < Nhkl; n++) {

          old.re = (array+n)->re ;
	  old.im = (array+n)->im ;   
         (array + n)->re = afac*old.re - fpp*old.im ;
         (array + n)->im = afac*old.im + fpp*old.re ;

     }

}
void	prepare_intersect_mask() 
{
	int	m, n ;
	char	*pmask, *pimask ;

     intermask = (char *)e_malloc(Nhkl*sizeof(char), "prepare_intersect_mask") ;

     for (n = 0, pimask = intermask; n < Nhkl; n++, pimask++)
	  *pimask = 1 ;

     for (m = 0, pmask = maskfo; m < NM; m++) {
	  if (relwt_der[m] > 0) {
               for (n = 0, pimask = intermask; n < Nhkl; n++, pmask++, pimask++)
	            *pimask *= *pmask ;
          }
	  else
	       pmask += Nhkl ;
     }
}

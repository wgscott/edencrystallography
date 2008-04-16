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

                                SOLVE

  Title:        Quadratic version of solver with MIR and MAD.
  Author:       Hanna Szoke
  Date:         10/1/94
  Function:     This program estimates electron densities using A. Szoke's 
		holographic restoration method and Dennis Goodman's 
		non-linear complex conjugate gradient optimizer, CCG.

		Output includes a log file plus the solution file in 
		physical space after each outer iteration,  written as a 
		binary file in units of el/voxel.  
		This serves as input to the postprocessor, REGRID.
		Other optional output includes an outlier report, listing all
		reflections that deviate by more than 4*sigma from the fobs
		and a listing of the costs at each step of the CCG.

	Solve.c works in conjunction with the following source files:

		ssetup.c does all the setup functions for Solve.
		ccgutil.c is the interface to ccg.c;
		cost.c, concosts.c have cost function calculations;
		fftwshare.c is the interface to the fftw package;
		highres.c for high-resolution option;
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

static	real	*Rabs ;		/* sqrt(R.re^2 + R.im^2) */
static	char	*P1_mask ;

COMPLEX	*R0 ;		/* Model structure factors for known atoms */
real 	*Full ;	    	/* diffraction pattern amplitudes of full molecule */
real	feden_min ;	/* for comparing cost */
real	F000 ;		/* value of *Full, used in Sayre's term */
real	hkl_cost ;	/* current value of hkl cost function */
real	sing_cost ;	/* current normalized value of singlet cost function */
real	trip_cost ;	/* current normalized value of triplet cost function */

FILE	*fp_cost ;					/* in util.h */
char	costfile_option = FALSE ;			/* in util.h */

	/* variables in 'eden.h' */

char	*maskfo ;	/* mask identifying fo input values */
COMPLEX	*twinned_fc ;

	/* variables in 'cost.h' */

COMPLEX *R ;    	/* Model structure factors for known+found atoms */

static	real	getRfactor() ;
static	void    makeStructureFactor() ; 
static	void	makeRabs() ;
static	void	outliers_summary() ;
static	void	outliers_details() ;
static	void	rescaleFobs() ;
static	void	calcChiSq() ;
static	real	getDiscrep() ;
static	void	fine_tune() ;
static	void	report_costs(char *) ;
static	void	monitor_It_frac(COMPLEX *, real *) ;
static	void	monitor_At_frac(COMPLEX *, real *) ;
static	void	monitor_t_frac(COMPLEX *, real *) ;
static	COMPLEX	*twinner(COMPLEX *) ;
static	real	get_scale_factor(int) ; 

void	solve_main(int argc, char *argv[])
{
     void	solver_setup(char *) ;
     int	checkSolutions() ;		/* in 'qshare.c' */
     void	writeSol(char *) ;	
     void       updateKnownValues() ; 
     void       reinitNump() ;
     int       	HoloMethod() ;			 /* in 'ccgutil.c' */
     int	doStdDev(char *) ;
     void	reportTotFunctCalls() ;
     void	do_inplace_cs_symmetrization(real *, real *) ;
		      			/* in 'crysutil.c' */

     void	hrWriteLists(char *, real *) ;
     void	hrSymmetrize(real *) ;

     real	R_fac ;			/* Current R factor */
     real	R_fac_init ;		/* Initial R factor */
     int        iter ;
     int	stdev_crit ;		/* continue? ret'd by doStdDev */
     int	Holo_crit ;		/* continue? ret'd by HoloMethod */
     int	update_crit ;		/* continue? ret'd by checkSolutions */
     char	problemName[MAXSTRING] = "" ;
     char	out_name[MAXSTRING] ;	/* generic name, used for output */

     if (argc < optind + 2) {	/* wrong number of arguments */
	ballpark(caller) ;
     }

  /******************************************************
  Pull off any prefix defining whereabouts of input file
  and remove the .inp extension.
  Fetch input, prepare for doing optimization.          
  ******************************************************/

     strip_path(input_filename, pwd_name) ;
     strncpy(problemName, pwd_name, strlen(pwd_name) - 4) ;
     sprintf(out_name, "%s", problemName) ;

     solver_setup(problemName) ;
     makeStructureFactor() ;

     if ((R_fac_init = R_fac = getRfactor()) < R_stop) {
          sprintf(message, 
	       "\nQuitting immediately - initial R factor is less than %g",
	       R_stop) ;
          EdenWarning(message) ;
          return ;
     }
        
  /******************************************************
  Prepare an initial outlier report.  
  ******************************************************/ 

     makeRabs() ; 

     outliers_summary() ; 
     if (verbose) 
           outliers_details() ; 

     if (useSig)
           calcChiSq() ;

     e_free(NM*Nhkl*sizeof(real), Rabs) ; 

     if (detwin)
          monitor_t_frac(R, Full);

  /************************************************************
  Iterate; 50 is a magic number to stop code if it doesn't 
  converge, but we've never hit it.  Typically, iter reaches
  2 to 5 before the code quits.
  ************************************************************/

     for (iter = 0;iter < 50;iter++) { 

          printTwice("\n") ;
          printTwice(timestamp()) ;
          printTwice("Applying holographic reconstruction - ") ;

	  sprintf(message, "iteration # %d", iter+1) ; 
	  printTwice(message) ;

	  if (costfile_option)
	       fprintf(fp_cost, "\niteration # %d\n", iter+1) ; 

       /************************************************************
       Do holographic reconstruction, symmetrize.  Add new 
       contributions to snump, write out an electron density map.
       ************************************************************/

          feden_min = getDiscrep() ;
	  stdev_crit = doStdDev("Initial standard deviation") ;  
          if (iter == 0)
	       report_costs("Initial") ;
						/**************/
          Holo_crit = HoloMethod();		/* THIS IS IT */
						/**************/
	  update_crit = checkSolutions(); 
		
          if (Ncs > 1) {
	     stdev_crit = doStdDev("Standard deviation before symmetrization") ;

             do_inplace_cs_symmetrization(nump, snump) ;  
	     if (HighRes)
		  hrSymmetrize(nump+Nptotal) ;

	     stdev_crit = doStdDev("Standard deviation after symmetrization") ; 
          }
	  else
	     stdev_crit = doStdDev("Standard deviation") ; 

          updateKnownValues() ; 

          writeSol(out_name) ; 

          if (HighRes) {

               sprintf(message, "%s.list", problemName) ;
               hrWriteLists(message, totp+Nptotal) ;
	       if (verbose) {
                    sprintf(message, "it%1d.list", iter+1) ;
                    hrWriteLists(message, totp+Nptotal) ;
	       }
          }
          
       /************************************************************
       Do an analysis of the current solution in (hkl) space.
       ************************************************************/

          makeStructureFactor() ;
          R_fac = getRfactor() ;
	  makeRabs() ;

          outliers_summary() ;
          if (verbose) 
               outliers_details() ;
          
	  if (useSig)
               calcChiSq() ;
	  rescaleFobs() ;

          e_free(NM*Nhkl*sizeof(real), Rabs) ; 

          if (detwin)
               monitor_t_frac(R, Full);

          fine_tune() ;

       /************************************************************
       Decide whether to go on, based on various figures of merit.
       ************************************************************/

          if (R_fac < R_stop) {
	       sprintf(message, 
		 "Stopping - Rfac is less than %g\n", R_stop) ;
               printTwice(message) ;
	       break ;
          }

          if ((stdev_crit == STOP) || 
	      (Holo_crit  == STOP) ||
	      (update_crit == STOP)) 
	       break ;

          reinitNump() ;
	  report_costs("Current") ;
     }

    report_costs("Final") ;

    sprintf(message, 
	 "\nOverall R factor changed from %g to %g", R_fac_init, R_fac) ;
    printTwice(message) ;

    reportTotFunctCalls() ;
}

				/****************************************
                 		 Calculate the R factor for each part of
				 the problem (native plus derivatives) 
				****************************************/
real getRfactor()
{
     void	shellHeader_R() ;	/* in hklutil.c */
     real	get1Rfac(char [], real *, COMPLEX *, int) ;		
     real	get1Rfac_weighted( char [], real *, COMPLEX *, int, real *) ;	

     real	thisr ;
     int	m ;
     char	what[MAXSTRING] ;
     COMPLEX	*localR ;

     shellHeader_R() ;
                                                                     
     localR = (detwin) ? twinner(R) : R ;

     if (NM > 1) {
          for (m = 0; m < NM; m++) {
	       if (relwt_der[m] > 0) {
	            if (m == 0)
                         sprintf(what, " Native: ") ;
	            else
                         sprintf(what, " Der. %d: ", m) ;
	            thisr = get1Rfac(what, Full+m*Nhkl, localR+m*Nhkl, Nhkl) ;
               }
          }
	  printTwice(" ") ;			/* extra linefeed */
          sprintf(what, " Total:  ") ;
          thisr = get1Rfac_weighted(what, Full, localR, NM*Nhkl, relwt_der) ;

     }
     else {
          sprintf(what, "         ") ;
          thisr = get1Rfac(what, Full, localR, Nhkl) ;
     }

     return(thisr) ;
}
     			/****************************************************
	  		Report fine-tuned f', f''.
     			****************************************************/
void	fine_tune()
{
     void	recalc_fp_fpp(int, real*, COMPLEX *, real, real) ;
     int	m, md, off ;


     if (NM > 1) {

          for (m = 0, md = 1; m < NM; m++, md++) {
	       if (anom_der_flag[md]) {
		    off = md*Nhkl ;
                    recalc_fp_fpp(md, nump, R0+off, fprime[md], f2prime[md]) ;
               }
          }
     }
}

     		/************************************************************
     		Assemble structure factors, R, from FFT-generated reflections 
     		from found atoms plus R0 (reflections from known atoms).
     		************************************************************/

void makeStructureFactor() 
{
     void	convolve(real *, real *, COMPLEX *, COMPLEX *) ;
     void	hrDoSFT(real *, COMPLEX *) ;
     void	Cmaskit(COMPLEX *, char *, int) ;

     static int first_time = TRUE ;
     int        m, n ;

     /************************************************************
     Make structure factors for the grid(s).  
     Propagate R throughout array.
     ************************************************************/

     convolve(snump, expfac, expfac_bcc, R) ;

     if (first_time)
          first_time = FALSE;
     else
          if (HighRes) 
               hrDoSFT(snump+Nptotal, R) ;

     if (NM > 1) {
	  for (n = Nhkl; n < Nhkl*NM; n++) {
	       m = n % Nhkl ;
	       (R+n)->re = (R+m)->re ;
	       (R+n)->im = (R+m)->im ;
          }
     }

     /************************************************************
     Add in known part of the native structure, write out at the
     end of each iteration (but not at the start of the run). 
     ************************************************************/

     for (n = 0; n < Nhkl*NM; n++) {
          (R+n)->re += (R0+n)->re ;
          (R+n)->im += (R0+n)->im ;
     }
     
     /************************************************************
     In preparation for the next iteration, apply mask to R. 
     ************************************************************/

     Cmaskit(R, maskfo, Nhkl*NM) ;
}

void	makeRabs()
{
     int	n ;
     real	Cabs() ;

     Rabs   = (real *)e_malloc(Nhkl*NM*sizeof(real), "makeRabs") ;

     for (n = 0; n < Nhkl*NM; n++) 
          *(Rabs+n) = Cabs(*(R+n)) ;
}

				/*********************************************
	        		Calculate summary information about the
				outliers for a data set. 
				*********************************************/
#define	NBUCKETS	5 

void outliers_summary()
{
     real	spread ;
     real	buckets[NBUCKETS] ;
     real	p1, p2, p3, p4 ;
     real	sqrt_relwt ;
     real	numer ;

     int	Ntot ;
     int	m, n, nn, b ;
     char	gotit ;
     char	which[10] ;

     if (useSig)
          printTwice(
	  "\n\tPercentages of Fcalcs within sigma intervals of Fobs:\n") ;
     else
          printTwice(
	  "\n\tPercentages of Fcalcs within intervals of Fobs:\n") ;
     printTwice(
  "interval:                1          2          3          4   outliers") ;

     for (m = 0; m < NM; m++) {		/* once for native and each der. */

	if (relwt_der[m] > 0) {

	  Ntot = 0 ;
	  sqrt_relwt = sqrt(relwt_der[m]) ;

          for (b = 0; b < NBUCKETS; b++) 
	       buckets[b] = 0 ;

          if (m == 0)
	       strcpy(which, "native") ;
          else 
	       sprintf(which, "Der %1d ", m) ;
	       
          for (nn = 1; nn < Nhkl; nn++) {

               n = m*Nhkl + nn ;

	       if (*(maskfo+n)) {   

	            numer = fabs(*(Full+n)-*(Rabs+n))*sqrt_relwt ;

                    spread = (useSig) ? numer / *(sigma+n) : numer ;
	            Ntot++ ;

	            for (b = 0; b < NBUCKETS-1; b++) {
		         gotit = FALSE ;
		         if (spread < (real) (b+1)) {
			      buckets[b] += 1. ;
			      gotit = TRUE ;
			      break ;
                         }
                    }
                    if (!gotit) 
		         buckets[NBUCKETS-1] += 1. ;
	       }
          }

          for (b = 0; b < NBUCKETS; b++) 
	       buckets[b] *= 100. / (real) Ntot ;

          sprintf(message, "%% %6s:       %10g %10g %10g %10g %10g", 
          which, buckets[0], buckets[1], buckets[2], buckets[3], buckets[4]) ;
          printTwice(message) ;
        }
     }

     p1 = erfc(1./SQRT2) ;
     p2 = erfc(2./SQRT2) ;
     p3 = erfc(3./SQRT2) ;
     p4 = erfc(4./SQRT2) ;

     sprintf(message, "%% for Gaussian: %10g %10g %10g %10g %10g", 
     100*(1.-p1), 100*(p1-p2), 100*(p2-p3), 100*(p3-p4), 100*p4) ;
     printTwice(message) ;
}

				/*********************************************
	        		Print out details about the far outliers in
				a data set. 
				*********************************************/

void outliers_details()
{
     FILE	*fp ;
     char	filename[12] ;
     real	spread ;
     real	rrr;
     real	sqrt_relwt ;
     real	numer ;

     int	m, n, nn, b ;
     int	h, k, l ;
     char	gotit ;
     char	which[10] ;

     for (m = 0; m < NM; m++) {

	if (relwt_der[m] > 0) {

	  sqrt_relwt = sqrt(relwt_der[m]) ;

	  sprintf(filename, "outlier%1d", m) ;
          if ((fp = fopen(filename, "w")) == NULL)
          EdenError("Cannot open outlier file") ;

          if (m == 0)
	       strcpy(which, "native") ;
          else 
	       sprintf(which, "der %1d ", m) ;

          if (useSig) {
               fprintf(fp, 
  "Outliers whose weighted deviation > 4*sigma for %s\n\n", which) ;

               fprintf(fp, 
"Res (A)         deviation    hkl          Fobs         |Fcalc|      Sig         \n") ;
          }
          else {
               fprintf(fp, 
  "Outliers whose weighted deviation > 4 for %s\n\n", which) ;

               fprintf(fp, 
"Res (A)         deviation    hkl          Fobs         |Fcalc|        \n") ;
          }
          for (nn = 1; nn < Nhkl; nn++) {	/* exclude (000) term */

               n = m*Nhkl + nn ;

	            if (*(maskfo+n)) {   

	            numer = (*(Full+n)-*(Rabs+n))*sqrt_relwt ;
	            spread = (useSig) ? numer / *(sigma+n) : numer ;

	            for (b = 0; b < NBUCKETS-1; b++) {
		         gotit = FALSE ;
		         if (fabs(spread) < (real) (b+1)) {
		              gotit = TRUE ;
			      break ;
                         }
                    }
                    if (!gotit) {

                         fetch_hkl(nn, &h, &k, &l) ;
			 rrr = 1. / sqrt(nusq(h, k, l)) ;
			 if (useSig)
		              fprintf(fp, 
			 " %10g   %10g (%3d %3d %3d) %10g  %10g  %10g\n",
		          rrr, spread, h, k, l, 
			  *(Full+n), *(Rabs+n), *(sigma+n)/sqrt_relwt) ;
			 else
		              fprintf(fp, 
			 " %10g   %10g (%3d %3d %3d) %10g  %10g\n",
		          rrr, spread, h, k, l, 
			  *(Full+n), *(Rabs+n)) ;

	            }
	       }
	    }
	    fclose(fp) ;
        }  
     }
} 

			/************************************************
       			Determine scale factors for native and MIR/MAD
			derivatives, based on a mask that uses the 
			intersection of all the maskfo's.
			************************************************/

void	rescaleFobs()
{
     void	scalefobs(real *, real) ;
     int	m ;
     real	numer = 0, denom = 0, norm ;
     real	rel_scale_m[MAXDER] ;	/* relative scale_m's	*/
     real	scale_m[MAXDER] ;	

     if (Nder == 0)
	  return ;

     for (m = 0; m < NM; m++) {
	  if (relwt_der[m] > 0) 
		scale_m[m] = get_scale_factor(m) ;
	  else 
		scale_m[m] = 0 ;
     }     

     for (m = 0; m < NM; m++) {
	 denom += relwt_der[m] ;
	 numer += relwt_der[m] * scale_m[m] ;
     }
     norm = numer / denom ;
 
     for (m = 0; m < NM; m++) {
	 rel_scale_m[m] = scale_m[m] / norm ;
	 abs_scale_m[m] *= rel_scale_m[m] ;

	 if (autoscale) {
              scalefobs(Full+m*Nhkl, rel_scale_m[m]) ;
              scalefobs(sigma+m*Nhkl, rel_scale_m[m]) ;
	 }
     }

     printTwice("") ;

     for (m = 0; m < NM; m++) {
         sprintf(message, "Relative scale[%d] = %g, absolute scale[%d] = %g",
	         m, rel_scale_m[m], m, abs_scale_m[m]) ;
         printTwice(message) ;
     }
}

			/************************************************
       			Determine ChiSq for native and each derivative.
			************************************************/
void	calcChiSq() 
{
     int	m, n, nn ;
     int	Ntot ;
     real	chi ;
     real	chisq[MAXDER] ;	
     real	wtchisq[MAXDER] ;	
     real	sqrt_relwt ;

    for (m = 0; m < NM; m++) {		/* once for native and each der. */

	Ntot = 0 ;
	chisq[m] = 0 ;
	wtchisq[m] = 0 ;

	if (relwt_der[m] > 0) {

	  sqrt_relwt = sqrt(relwt_der[m]) ;

          for (nn = 1; nn < Nhkl; nn++) {	/* skip (000) term */

               n = m*Nhkl + nn ;

	       if (*(maskfo+n)) {   
	            chi = (*(Full+n) - *(Rabs+n)) * sqrt_relwt / *(sigma+n) ;
	            chisq[m] += chi*chi ;
	            wtchisq[m] += chi*chi / sqrt(1. + chi*chi/16.) ;
	            Ntot++ ;
               }
          }

	  chisq[m] /= Ntot ;
	  wtchisq[m] /= Ntot ;

	  /***relwt_der[m] /= sqrt(chisq[m]) ;***/
	  /***relwt_der[m] *= 0.5*(1/chisq[m] + 1) ; ***/
        }
     }
     printTwice("") ;

     if (NM == 1) {
	  sprintf(message,
"Goodness-of-fit: chisq = %g, robust chisq = %g", chisq[0], wtchisq[0]) ;
	  printTwice(message) ;
     }
     else {
	  printTwice("Goodness-of-fit parameters:") ;
          for (m = 0; m < NM; m++) {
            sprintf(message, 
"Chi sq[%d] = %g, robust Chi sq[%d] = %g, relative weight[%d] = %g",
		    m, chisq[m], m, wtchisq[m], m, relwt_der[m]) ;
	   
	    printTwice(message) ;
          }
     }
}

		/*************************************************************
		We look for a scale factor, s, that minimizes:

			f = (1/2) Sum [ (hkl_weight/s)^2 * (|Fc| - s|Fo|)^2 ],

		where Sum is over Nhkl (for native or derivative).  The
		factor s will multiply both Fo and sigma (which is in the
		denominator of hkl_weight).  We then find s by solving 
		df/ds = 0.
		Note that the mask used for identifying relevant terms is
		the intersection of the maskfo's for native+derivatives;
		see prepare_intersect_mask() in ssetup.c.
		*************************************************************/

real get_scale_factor(int m)
			/* m: 0 for native, 1, 2, ... for derivatives */
{
     int	n ;
     real	wsq ;
     real	numer = 0, denom = 0 ;

     real	*pRabs ;		/* local pointers for convenience */
     real	*pFo, *pwt ;
     char	*pmask ;

     pRabs=Rabs + m*Nhkl ; 
     pFo=Full + m*Nhkl ; 
     pwt=hkl_weight + m*Nhkl ; 
     pmask = intermask ; 

     for (n = 0; n < Nhkl; n++, pRabs++, pFo++, pwt++, pmask++) {
	  if ((*pmask) && (n > 0)) {
               wsq = *pwt * *pwt ;
               numer += wsq * *pRabs * *pRabs ;
	       denom += wsq * *pFo * *pRabs ; 
	  }
     }

     return (numer / denom) ;
}

  	/************************************************************
	Get "discrepancy value" (see paper VI, Appendix D).  
	************************************************************/

real getDiscrep() 	
{
     real	SsqWts ;	
     int	m, n ;
     real	*wt ;
     char	*mask ;
     real	fedmin = 0 ;

     if (!useSig)
	  fedmin = DISCRP ;		/***	DISCRP defined in eden.h */

     else {

	  wt = hkl_weight ;
	  mask = maskfo ;

          for (m = 0; m < NM; m++) {
	  
               for (n = 0, SsqWts = 0; n < Nhkl; n++, wt++, mask++) {

	            SsqWts += *wt * *wt * *mask ;
               }
    
	       fedmin += SsqWts * relwt_der[m] ;

          }
	  fedmin *= 0.5 / AveWt ;

          sprintf(message, "%6.1e", fedmin) ;
          sscanf (message, "%g", &t1) ;
	  fedmin = t1 ;
     }
     fedmin *= discrp_frac ;            /* new 5/17/99 */
     sprintf(message,
       "\nStopping criterion for the (hkl) cost function is %6.1e", fedmin) ;
     printTwice(message) ;

     return(fedmin) ;
}

void report_costs(char *which)
{
        int	k ;
        real	dphase_singlet() ;
        real	dphase_triplet() ;
	real	dph ;

        sprintf(message,
          "%s value of the (hkl) cost function is %6.1e", which, hkl_cost) ;
        printTwice(message) ;

        for (k = 0; k < Nconstraints; k++) 
	
          switch (con_type[k]) {

	  case SINGLET :
             sprintf(message,
               "%s normalized value of the singlet cost function is %6.1e",
	       which, sing_cost) ;
             printTwice(message) ;
             
             dph = dphase_singlet() ;
             sprintf(message,
   "and weighted average phase error of the singlet cost function is %6.1e deg.",
	       dph) ;
             printTwice(message) ;
               break ;
     
	  case TRIPLET :
             sprintf(message,
               "%s normalized value of the triplet cost function is %6.1e",
	       which, trip_cost) ;
             printTwice(message) ;
             
             dph = dphase_triplet() ;
             sprintf(message,
   "and weighted average phase error of the triplet cost function is %6.1e deg.",
	       dph) ;
             printTwice(message) ;
               break ;
     
	  default :
	       break ;
          }
}

COMPLEX	*twinner(COMPLEX *Fc)
{   
     char	*get_P1mask() ;				/* in hklutil.c */
     int	fetch_twin(int, int *, int *) ;		/* in hklutil.c */
     void	Cmaskit(COMPLEX *, char *, int) ;
     int	n, t_n ;
     int	ct ;
     COMPLEX	*pFc, *pTw ;
     COMPLEX	ftw ;
     real	tbar ;
     static	char	first = TRUE ;

     if (first) {
          first = FALSE ;
          P1_mask = get_P1mask() ;
     }

     tbar = 1. - t_frac ;

     if (t_type == 'A') 	/* amplitude detwinning */

          for (n = 0, pTw = twinned_fc, pFc = Fc; n < Nhkl; n++, pTw++, pFc++) {

               pTw->re = tbar * pFc ->re ;
               pTw->im = tbar * pFc ->im ;

               t_n = fetch_twin(n, t_matrix, &ct) ;
               ftw = *(Fc + t_n) ;
               pTw->re += t_frac * ftw.re ;
               pTw->im += ct * t_frac * ftw.im ;
          }

     else 			/* intensity detwinning */

          for (n = 0, pTw = twinned_fc, pFc = Fc; n < Nhkl; n++, pTw++, pFc++) {

               pTw->re = tbar * (pFc->re * pFc->re + pFc->im * pFc->im) ;

               t_n = fetch_twin(n, t_matrix, &ct) ;
               ftw = *(Fc + t_n) ;
               pTw->re += t_frac * (ftw.re * ftw.re + ftw.im * ftw.im) ;
               pTw->re = sqrt(pTw->re) ;
               pTw->im = 0 ; 
          }

     Cmaskit(twinned_fc, P1_mask, Nhkl) ;

     return (twinned_fc) ;
}

  /******************************************************
  We monitor BOTH amplitude AND intensity detwinning, as
  a check on the validity of our choice!
  ******************************************************/

void	monitor_t_frac(COMPLEX *Fclocal, real *Folocal)
{ 
     monitor_At_frac(Fclocal, Folocal) ;
     monitor_It_frac(Fclocal, Folocal) ;
}

void	monitor_It_frac(COMPLEX *Fclocal, real *Folocal)
{ 
     int	fetch_twin(int, int *, int *) ;
     real	Cabs() ;
     int	n, t_n ;
     COMPLEX	*pFc ;
     real 	*pFo, *pwt ;
     int	ct;
     real	denom, numer ;
     real	fcabs, tabs, sqfo, sqfc, sqtw, coef ;
     real	Itf ;
     char	*Pm ;

     numer = denom = 0 ;

     for (n = 0, pFo=Folocal, pFc =Fclocal, pwt = hkl_weight, Pm = P1_mask; 
          n < Nhkl; n++, pFo++, pFc++, pwt++, Pm++) 

       if (n > 0) {

          sqfo = *pFo * *pFo ;
          fcabs = Cabs(*pFc) ;
          sqfc = fcabs * fcabs ;
          t_n = fetch_twin(n, t_matrix, &ct) ;
          tabs = Cabs(*(Fclocal+t_n)) ;
          sqtw = tabs * tabs ;
          coef = *pwt * *pwt * *Pm ;
          numer += coef * (sqfo - sqfc) * (sqfc - sqtw) ;
          denom += coef * (sqfc - sqtw) * (sqfc - sqtw) ;
       }

     if (denom > 0) {
          Itf = -numer / denom ;
          sprintf(message, 
     "\nThe value of t_frac for intensity detwinning is %g\n", Itf) ;
          printTwice(message) ;
     }
     else {
          sprintf(message, 
     "\nThe value of t_frac for intensity detwinning is undefined\n") ;
          printTwice(message) ;
     }
}

void	monitor_At_frac(COMPLEX *Fclocal, real *Folocal)
{ 
     int	fetch_twin(int, int *, int *) ;
     real	Cabs() ;
     void	Cmaskit(COMPLEX *, char *, int) ;
     int	n, t_n ;
     COMPLEX	*pFc, *pTw ;
     real 	*pFo, *pwt ;
     int	ct;
     real	atwin, abar, grad ;
     real	lamw2, resid, fabstw ;
     COMPLEX	dFc, exp_phi ;

     printTwice(
"\nInformation for monitoring the value of t_frac for amplitude detwinning:") ;

     printTwice("t\t\tdf/dt") ;

     for (atwin = 0.8 * t_frac; atwin < 1.21* t_frac; atwin += 0.1 * t_frac) {

          if (atwin > 0.5)
               break;

          abar = 1. - atwin;

	  /* we use twins within the loop  */

          for (n = 0, pTw = twinned_fc, pFc = Fclocal; n < Nhkl; 
               n++, pTw++, pFc++) {
               pTw->re = abar * pFc ->re ;
               pTw->im = abar * pFc ->im ;
          }

          for (n = 0, pTw = twinned_fc; n < Nhkl; n++, pTw++) {
               t_n = fetch_twin(n, t_matrix, &ct) ;
               pTw->re += atwin * (Fclocal+t_n)->re ;
               pTw->im += ct * atwin * (Fclocal+t_n)->im ;
          }

	  Cmaskit(twinned_fc, P1_mask, Nhkl) ;

	/*******************************************************
        Now we do the gradient wrt the twinning fraction, atwin.
	Note that we exclude n = 0. 
        We handle exception processing exactly as in the hkl 
        cost function.  Note: hkl_weight includes maskfo!
	*******************************************************/

          grad = 0. ;

          for (n = 0, pFo=Folocal, pTw = twinned_fc, pwt = hkl_weight;
	       n < Nhkl; n++, pFo++, pTw++, pwt++) {

               fabstw = Cabs(*pTw) ;
               resid = fabstw - *pFo ;
               lamw2 = *pwt * *pwt / 2 ;

               t_n = fetch_twin(n, t_matrix, &ct) ;

               dFc.re = (Fclocal + t_n)->re - (Fclocal + n)->re ;
               dFc.im = ct * (Fclocal + t_n)->im - (Fclocal + n)->im ;

                    if (fabstw > EPS_PHI) {
                         exp_phi.re = pTw->re / fabstw ;
                         exp_phi.im = pTw->im / fabstw ;
                    }
                    else {
                         exp_phi.re = (pTw->re > 0) ? 1. : -1. ;
                         exp_phi.im = 0. ;
                    }
                    grad += lamw2 * resid *
                            (dFc.re * exp_phi.re + dFc.im * exp_phi.im) ;
          }
          sprintf(message, "%f\t%g", atwin, grad) ; 
          printTwice(message);
     }
}

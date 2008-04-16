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

                                COST.C

  Title:        Source code for cost functions w/o constraints.
  Author:       Hanna Szoke
  Date:         6/10/00 (constraint costs separated off)
  Function:     This module contains functions used by Solve and Back
		for calculating the various hkl cost functions.

*******************************************************************************/

#include "eden.h"		
#include "eden_Np.h"		/* Np space arrays also in back.c, qshare.c */
#include "mir.h"		/* Mir arrays also in init.c   */
#include "cost.h"		/* for cost function calculation */

COMPLEX *Oar ;          /* array for holding FFT transform of np's */
static	int	first = TRUE ;

static	real	e_funct_hkl(real *, real *) ;
static	real	prep_grad() ;
static	real	prep_grad_Itw() ;
static	real	e_robust_funct_hkl(real *, real *) ;
static	void	compute_grad(real *) ;
static	void	convolution(real *) ;


void e_funct(real *csol, real *ressq, real *grad, int *iqflag)

/*
real	*csol ;		current solution array	
real	*ressq ;	(returned) residue  
real	*grad ;		(returned) gradient arra
int	*iqflag ;	(returned) flag re success	
*/

{
     real	e_funct_phasext(real *, real, int) ;
     real	e_funct_np(real *, real, int) ;
     real	e_funct_cs(real *, real *, real) ;
     real	e_funct_sayre(real *, real, real) ;
     real	e_funct_singlet(real *, real) ;
     real	e_funct_triplet(real *, real) ;

     real	*pcsol, *ptotp, *psnump, *pknownp ;
     real	con_cost[MAXCONSTR] ;
     real	*pg ;
     real	lamcoef ;
     int	i, n, target_number ;

     if (first) {
	  first = FALSE ;
	  
          if (robust) {
	       sprintf(message, 
	       "Using the robust hkl cost function with coefficient %g.", 
	       rob_factor) ;
	       printTwice(message) ;
          }
          else
	       printTwice("Using the original hkl cost function.") ;
     }

     /************************************************************
     Initialize constraint costs.
     Initialize gradient; each partial cost function will 
     accumulate into the gradient, so that ordering is arbitrary.
     ************************************************************/

     for (i = 0; i < MAXCONSTR; i++)
          con_cost[i] = 0 ;

     for (n = 0, pg = grad; n < Npextended; n++, pg++)
          *pg = 0 ;

     /************************************************************
     Do the basic cost calculation in Fourier space. 
     Then compute physical space cost functions based on the 
     various constraints:
        target(s), crystallographic symmetry, 
        Sayre's equation, and phase extension (in Fourier space).
     In most cases, the space cost function is applied to the total
     electrons/voxel in the starting model (knownp) + those found 
     in previous iterations (snump) + those found in this iteration 
     (csol).  The exception in e_funct_cs(), for which only csol 
     need be considered.
     ************************************************************/

     if (robust)
	  hkl_cost = e_robust_funct_hkl(csol, grad) ;
     else
	  hkl_cost = e_funct_hkl(csol, grad) ;

     if (Nconstraints > 0) {

	for (n=0, ptotp=totp, pcsol=csol, psnump=snump, pknownp=knownp ;
             n < Nptotal; n++, ptotp++, pcsol++, psnump++, pknownp++)

	     *ptotp = *pcsol + *psnump + *pknownp ;


        for (i = 0, target_number = 0; i < Nconstraints; i++) {

	     lamcoef = relwt_con[i] * con_coef[i] ;

	     switch (con_type[i]) {

	     case TARGET :
		  con_cost[i] = e_funct_np(grad, lamcoef, target_number) ;
		  target_number++ ;		/* .. to select next target */
		  break ;

	     case PHASE_EXT :
		  con_cost[i] = e_funct_phasext(grad, lamcoef, target_number) ;
		  target_number++ ;		/* .. to select next target */
		  break ;

	     case CS :
		  con_cost[i] = e_funct_cs(csol, grad, lamcoef) ;
		  break ;

	     case SAYRE :
		  con_cost[i] = e_funct_sayre(grad, lamcoef, cost_addend[i]) ;
		  break ;

             case SINGLET :
		  con_cost[i] = e_funct_singlet(grad, lamcoef) ;
		  sing_cost = con_cost[i] / relwt_con[i] ;
		  break ;

             case TRIPLET :
		  con_cost[i] = e_funct_triplet(grad, lamcoef) ;
		  trip_cost = con_cost[i] / relwt_con[i] ;
		  break ;

	     default  :
		  break ;

             }	
         }
     }
     /************************************************************
     Summarize the full cost function and report.  Separate terms 
     are reported in the appropriate e_funct_*(). 
     ************************************************************/

        *ressq = hkl_cost ;

        for (i = 0; i < Nconstraints; i++)
	     *ressq += con_cost[i] ;

        if (costfile_option) {
             if (Nconstraints > 0) 
	          fprintf(fp_cost, ", Total: %10.4e\n", *ressq) ;
	     else
	          fprintf(fp_cost, "\n") ;
        }

	fflush(fp_cost) ;

	*iqflag = (hkl_cost<= feden_min) ? -1:0 ; 	

	return ;
}

        /***********************************************************************
	This formerly large function has been broken up into three parts:
	convolution() that creates Oar from the input csol; prep_grad() and
	prep_grad_Itw() that prepare gradient arrays; and compute_grad() 
	that does the deconvolutions.
	real	*csol ;		current solution array
	real	*grad ;		(returned) gradient array
        ***********************************************************************/

real	e_funct_hkl(real *csol, real *grad)
{
	void	hrDoSFT(real *, COMPLEX *) ;
        void	hrGradSFT(COMPLEX *, real *) ;
	void	Cmaskit() ;			/* in hklutil.c */

	real	ef_hkl_cost ;

        convolution(csol) ;

        if (HighRes)
             hrDoSFT(csol+Nptotal, Oar) ;

	Cmaskit(Oar, maskfo, NM*Nhkl) ;

        if ((detwin) && (t_type == 'I'))
             ef_hkl_cost = prep_grad_Itw() ;
        else 
             ef_hkl_cost = prep_grad() ;
        
        compute_grad(grad) ;
        if (HighRes) 
             hrGradSFT(Oar, grad+Nptotal) ;
        
	return(ef_hkl_cost) ;
}

        /***********************************************************************
	For each resolution, r, do an FFT on csol and convolve the result with 
	the appropriate exponential factor;  do this once only for each 
	resolution and then propagate the information among the derivatives 
	which share that resolution.  See equation (32). paper V.
	real	*csol ;		current solution array	
        ***********************************************************************/

void	convolution(real *csol)
{
     int	m, n, r ;
     int	mfirst ;
     COMPLEX	*pOar, *thisp, *firstp;	
     void	convolve(real *, real *, COMPLEX *, COMPLEX *) ;

     for (r = 0; r < Nres; r++) {

	mfirst = -1 ;  

        for (m = 0, pOar=Oar; m < NM; m++, pOar += Nhkl) {

	   if (useset[m] == r) {

	      if (mfirst == -1) {

		 mfirst = m ;

	         convolve(csol, rexpfac[r], rexpfac_bcc[r], pOar) ;

              }
	      else {

	          for (n = 0, thisp = pOar, firstp = Oar+mfirst*Nhkl; 
		       n < Nhkl; n++, thisp++, firstp++) {

	               thisp->re = firstp->re ;
	               thisp->im = firstp->im ;

                  } 
              }
           }
        }
     }
     return ;
}

real	prep_grad()
{
	static	char	pg_first = TRUE ;
	static	real	tbar ;
	int	fetch_twin(int, int *, int *) ;

	int	m, n ;
	int	ct, t_n ;
	real	dsq, resid ;
	real	wsq ;
	real	lamw2 ;
        real	subres[MAXRES] ;
	COMPLEX	exp_phi ;
	COMPLEX	*pOar, *pR ;	
	real	*pF, *wt ;
 	COMPLEX	*ptwin ;
	COMPLEX	ftw ;
        real	Cabs();

        /***********************************************************************
         Replace O by R+O in Oar; replace R+O by twinned version, insofar as
	 twinning is enabled;  calculate |R+O| and hold it in dsq;
	 Put |R+O| - Fobs in resid and collect the weighted residuals.
	 Compute residues for native and derivatives separately.
	
         Also, prepare to compute the gradient; exp_phi is (R + O) / |R + O|.
	 (We must check the denominator and take action if |R+O| is too small.)
	 Reuse Oar for kernel in deconvolution, wt^2*(|R+O|-|Fobs|)*exp_phi/2.
	 See Eq. (34) - (39). Paper V.  (We have established that the special 
	 treatment of *pF (Fo) == 0 is unnecessary but harmless.)
        ***********************************************************************/
	 
        for (m = 0, pOar=Oar, pR=R; m < NM; m++) {
	     for (n = 0; n < Nhkl; n++, pOar++, pR++) {
	          pOar->re += pR->re ;
	          pOar->im += pR->im ;
             }
        }

        if (detwin) {
           if (pg_first) {
	      pg_first = FALSE ;
              tbar = 1. - t_frac;
	   }

           for (n=0, ptwin=twinned_fc, pOar=Oar; n < Nhkl; 
                n++, ptwin++, pOar++) {
              ptwin->re = pOar->re ;
              ptwin->im = pOar->im ;
           }

           for (n=0, ptwin=twinned_fc, pOar=Oar; n < Nhkl; 
                n++, ptwin++, pOar++) {

               pOar->re = tbar * ptwin->re ;
               pOar->im = tbar * ptwin->im ;

               t_n = fetch_twin(n, t_matrix, &ct) ;
               ftw = *(twinned_fc + t_n) ;
               pOar->re += t_frac * ftw.re ;
               pOar->im += ct * t_frac * ftw.im ;
          }
        }

        for (m = 0, pOar=Oar, pF=Full, wt=hkl_weight; m < NM; m++) {
	     subres[m] = 0 ;
	     for (n = 0; n < Nhkl; n++, pOar++, pF++, wt++) {
                  dsq = Cabs(*pOar) ; 
                  resid = dsq - *pF ;
	          wsq = *wt * *wt ;
	          subres[m] += wsq * resid * resid ;
		  lamw2 = relwt_der[m] * wsq / 2 ;

	          if (*pF > 0) {
		       if (dsq > EPS_PHI) {
		            exp_phi.re = pOar->re / dsq ;
		            exp_phi.im = pOar->im / dsq ;
                       }
		       else {
		            exp_phi.re = (pOar->re >= 0) ? 1. : -1. ;
		            exp_phi.im = 0. ;
                       }
		       pOar->re = lamw2 * resid * exp_phi.re ;
		       pOar->im = lamw2 * resid * exp_phi.im ;
                  } 
	          else {
	               pOar->re *= lamw2 ;
	               pOar->im *= lamw2 ;
                  } 
             } 
	     subres[m] *= relwt_der[m] / 2 ;
        } 

        if (detwin) {

           for (n=0, ptwin=twinned_fc, pOar=Oar; n < Nhkl; 
                n++, ptwin++, pOar++) {
              ptwin->re = pOar->re ;
              ptwin->im = pOar->im ;
           }

           for (n=0, ptwin=twinned_fc, pOar=Oar; n < Nhkl; 
                n++, ptwin++, pOar++) {

               pOar->re = tbar * ptwin->re ;
               pOar->im = tbar * ptwin->im ;

               t_n = fetch_twin(n, t_matrix, &ct) ;
               ftw = *(twinned_fc + t_n) ;
               pOar->re += t_frac * ftw.re ;
               pOar->im += ct * t_frac * ftw.im ;
          }
        }

        for (m = 0, hkl_cost = 0; m < NM; m++)
	     hkl_cost += subres[m] ;

	if (costfile_option) {
	     fprintf(fp_cost, "Hkl nat: %10.4e", subres[0]) ;
             if (NM > 1) {
	          for (m = 1; m < NM; m++) 
		       fprintf(fp_cost, ", der #%d: %10.4e", m, subres[m]) ;
                  fprintf(fp_cost, ", tot: %10.4e", hkl_cost) ;
             }
        }
	return (hkl_cost) ;
}

real	prep_grad_Itw()
{
	static	char	Ifirst = TRUE ;
	static	real	*At ;
	static	real	tbar ;
	int	fetch_twin(int, int *, int *) ;

	int	m, n;
	int	ct, t_n ;
	real	*pAt ;
	real	resid ;
	real	wsq ;
	real	lamw2 ;
        real	subres[MAXRES] ;
	COMPLEX	*pOar, *pR ;	
 	COMPLEX	*ptwin ;
	COMPLEX	ftw ;
	real	*pF, *wt ;

	if (Ifirst) {
	  Ifirst = FALSE ;
          tbar = 1. - t_frac;
          At = (real *) e_malloc(Nhkl*sizeof(real), "prep_grad_Itw") ;
	}
        /***********************************************************************
         Replace O by R+O in Oar; do NOT replace R+O by twinned version, 
	 but instead keep them both!  We need twinned version for At (aka dsq),
	 twinned_fc for (1-t)*Fc/At (at h) + t* Fc/At (at ht)
	 Finally, we use Oar for the gradient array (as always).
        ***********************************************************************/
	 
        for (m = 0, pOar=Oar, pR=R; m < NM; m++) {
	     for (n = 0; n < Nhkl; n++, pOar++, pR++) {
	          pOar->re += pR->re ;
	          pOar->im += pR->im ;
             }
        }
        for (m = 0, pAt = At, pOar = Oar; m < NM; m++) {
          for (n = 0; n < Nhkl; n++, pAt++, pOar++) {

               *pAt = tbar * (pOar->re * pOar->re + pOar->im * pOar->im) ;

               t_n = fetch_twin(n, t_matrix, &ct) ;
               ftw = *(Oar + t_n) ;
               *pAt += t_frac * (ftw.re * ftw.re + ftw.im * ftw.im) ;
               *pAt = sqrt(*pAt) ;
          }
        }

        for (m = 0, ptwin = twinned_fc, pOar = Oar; m < NM; m++) {
	     for (n = 0; n < Nhkl; n++, ptwin++, pOar++) {
	          ptwin->re = 0. ;
	          ptwin->im = 0. ;
             }
        }

        for (m = 0, ptwin = twinned_fc, pOar = Oar, pAt=At; m < NM; m++) {
	     for (n = 0; n < Nhkl; n++, ptwin++, pOar++, pAt++) {

               if (*pAt > EPS_PHI) {

	          ptwin->re += tbar*pOar->re / *pAt;
	          ptwin->im += tbar*pOar->im / *pAt;

                  t_n = fetch_twin(n, t_matrix, &ct) ;
                  (twinned_fc + t_n)->re += t_frac * pOar->re / *pAt;
                  (twinned_fc + t_n)->im += ct * t_frac * pOar->im / *pAt;
               }
             }
        }

        for (m = 0, pOar=Oar, pF=Full, wt=hkl_weight, ptwin=twinned_fc, pAt=At; 
	     m < NM; m++) {

	     subres[m] = 0 ;

	     for (n = 0; n < Nhkl; n++, pOar++, pF++, wt++, ptwin++, pAt++) {
                  resid = *pAt - *pF ;
	          wsq = *wt * *wt ;
	          subres[m] += wsq * resid * resid ;
		  lamw2 = relwt_der[m] * wsq / 2 ;

		  pOar->re = lamw2 * resid * ptwin->re ;
		  pOar->im = lamw2 * resid * ptwin->im ;
              }
	      subres[m] *= relwt_der[m] / 2 ;
        } 

        for (m = 0, hkl_cost = 0; m < NM; m++)
	     hkl_cost += subres[m] ;

	if (costfile_option) {
	     fprintf(fp_cost, "Hkl nat: %10.4e", subres[0]) ;
             if (NM > 1) {
	          for (m = 1; m < NM; m++) 
		       fprintf(fp_cost, ", der #%d: %10.4e", m, subres[m]) ;
                  fprintf(fp_cost, ", tot: %10.4e", hkl_cost) ;
             }
	     fflush(fp_cost) ;
        } 
	return (hkl_cost) ;
}

void	compute_grad(real *grad)
{
        void	deconvolve(real *, real *, COMPLEX *, COMPLEX *) ;

	int	j, m, n, r ;
	int	mfirst ;
	COMPLEX	*pOar, *thisp, *firstp ;	
	real	*pgrad ;

        /***********************************************************************
        Compute the gradient: g(back) = Re [DFT(-) (Diff * exp(-eta...)) ].
        SUM over the kernel elements that share a resolution before 
	deconvolving.  See Eq. (40), Paper V.
	Note that deconvolution now implies accumulation into the grad array;
	this means that the grad array must be initialized (here and in Back
	and in phase extension) and deconvolve() must accumulate!
        ***********************************************************************/

	for (j = 0, pgrad = grad; j < Nptotal; j++, pgrad++)
	     *pgrad = 0 ;

        for (r = 0; r < Nres; r++) {

	   mfirst = -1 ;  

           for (m = 0, pOar=Oar; m < NM; m++, pOar += Nhkl) {

	      if (useset[m] == r) {

	         if (mfirst == -1) 
		    mfirst = m ;

                 else 
	            for (n = 0, thisp = pOar, firstp = Oar+mfirst*Nhkl; 
		       n < Nhkl; n++, thisp++, firstp++) {

	               firstp->re += thisp->re ;
	               firstp->im += thisp->im ;
                  } 
               }
           }
        deconvolve(grad, rexpfac[r], rexpfac_bcc[r], Oar + mfirst*Nhkl) ;
     }

     for (j = 0, pgrad = grad; j < Nptotal; j++, pgrad++)
          *pgrad *= 2 ;

}

real	e_robust_funct_hkl(real *csol, real *grad)
{
	int	m, n ;
	real	dsq, resid, subres[MAXDER] ;
	real	Aprime, cprime ;
	real	co_factor ;
	real	der_factor ;
/*	real	step_size_correction ;  Didn't seem to help matters... */
	COMPLEX	exp_phi ;

	COMPLEX	*pOar, *pR ;	
	real	*pF, *wt ;

        real	Cabs();

int N1 = 0, N2 = 0, N3 = 0 ;

     convolution(csol) ; 

        /***********************************************************************
         Replace O by R+O in Oar; calculate |R+O| and hold it in dsq;
	 Put |R+O| - Fobs in resid and collect the weighted residuals.

	 Compute residues (for native and derivatives separately), using a
	 robust cost function that is linear beyond the inner well -

	 f = A' (x/c')^2           if |x| < c'   and
	 f = A' ( 2|x|/c' - 1.)    if |x| > c'
	 
	 where x = resid, c' = Ro*sqrt(norm_fac) / *wt and A' = Ro^2*norm_fac.
	 and Ro is the (input) rob_factor.

         Also, prepare to compute the gradient; exp_phi is (R + O) / |R + O|.
	 We use a constant gradient beyond the inner well.
	 (We must check the denominator and take action if |R+O| is too small.)
	 The F(000) is always handled using the quadratic form!

	 Reuse Oar for kernel in deconvolution, wt^2*(|R+O|-|Fobs|)*exp_phi/2.
	 See Eq. (34) - (39). Paper V.
        ***********************************************************************/

        for (m = 0, pOar=Oar, pR=R, pF=Full, wt=hkl_weight; m < NM; m++) {
	     subres[m] = 0 ;
	     for (n = 0; n < Nhkl; n++, pOar++, pR++, pF++, wt++) {
	          pOar->re += pR->re ;
	          pOar->im += pR->im ;
                  dsq = Cabs(*pOar) ; 
                  resid = dsq - *pF ;
		  cprime = rob_factor * sqrt(norm_fac) / *wt ;
		  Aprime = rob_factor * rob_factor * norm_fac ;
	 
	          if ((resid > cprime && n > 0)) {
		       subres[m] += Aprime * (2 * resid / cprime - 1.) ;
		       co_factor = Aprime /cprime ;
		       der_factor = 1. ;
	/*	       step_size_correction = 2 - cprime / resid ; */
N1++ ;
                  }
	          else if ((resid < -cprime && n > 0)) {
		       subres[m] += Aprime * (-2 * resid / cprime  - 1.) ;
		       co_factor = Aprime /cprime ;
		       der_factor = -1. ;
	/*	       step_size_correction = 2 + cprime / resid ; */
N2++ ;
                  }
	          else {
		       subres[m] += Aprime * (resid/cprime) * (resid/cprime) ;
		       co_factor = Aprime / (cprime * cprime) ;
		       der_factor = resid ;
	/*	       step_size_correction = 1. ; */
N3++ ;
                  }
		  co_factor *= relwt_der[m] / 2 ;

	          if (*pF > 0) {
		       if (dsq > EPS_PHI) {
		            exp_phi.re = pOar->re / dsq ;
		            exp_phi.im = pOar->im / dsq ;
                       }
		       else {
		            exp_phi.re = (pOar->re >= 0) ? 1. : -1. ;
		            exp_phi.im = 0. ;
                       }
/*		       pOar->re = step_size_correction * 
				  co_factor * der_factor * exp_phi.re ;
		       pOar->im = step_size_correction * 
				  co_factor * der_factor * exp_phi.im ; 
*/
		       pOar->re = co_factor * der_factor * exp_phi.re ;
		       pOar->im = co_factor * der_factor * exp_phi.im ;
                  } 
	          else {
	               pOar->re *= co_factor ;
	               pOar->im *= co_factor ;
                  } 
             } 
	     subres[m] *= relwt_der[m] / 2 ;
        } 

        for (m = 0, hkl_cost = 0; m < NM; m++)
	     hkl_cost += subres[m] ;

        compute_grad(grad) ;

        /* Reports */

	if (costfile_option) {
	     fprintf(fp_cost, "Hkl nat: %10.4e cats: %6i, %6i, %6i", 
	     subres[0],N1,N2,N3) ;
             if (NM > 1) {

	          for (m = 1; m < NM; m++) 
		       fprintf(fp_cost, ", der #%d: %10.4e", m, subres[m]) ;
                  fprintf(fp_cost, ", tot: %10.4e", hkl_cost) ;
             }
        }

	return (hkl_cost) ;
}
void	b_funct(real *nump, real *fval, real *grad, int *iqflag)
{
	int	n ;
	real	ressq = 0 ;
        void	convolve(real *, real *, COMPLEX *, COMPLEX *) ;
        void	deconvolve(real *, real *, COMPLEX *, COMPLEX *) ;
	COMPLEX	*pOar, *pFc0 ;
	real	*wt ;
	char	*pmaskfc ;

	convolve(nump, expfac, expfac_bcc, Oar) ;

	pOar = Oar ;
	pFc0  = Fc0  ;
	wt = hkl_weight ;
	pmaskfc = maskfc ;

	for (n = 0; n < Nhkl; n++, pOar++, pFc0++, wt++, pmaskfc++) {

	   pOar->re *=  *pmaskfc ;
	   pOar->im *=  *pmaskfc ;

	   pOar->re -=  pFc0->re ;
	   pOar->im -=  pFc0->im ;

	   pOar->re *=  *wt ;
	   pOar->im *=  *wt ;

	   ressq +=  pOar->re * pOar->re + pOar->im * pOar->im ;
        } 

	*fval = .5 * ressq ;

	if (costfile_option) {
	     fprintf(fp_cost, "Residual = %10.4e\n", *fval) ; 
	     fflush(fp_cost) ;
        }

	*iqflag = 0 ;
	if(ressq <= DISCRP*discrp_frac)         /* DISCRP is defined in eden.h,
						discrp_frac is input var. */
	     *iqflag = -1 ;

/*** Compute the gradient: g(back) = Re [DFT(-) (Diff * exp(-eta...)) ] ***/

        for (n = 0; n < Nptotal; n++)
	     *(grad+n) = 0 ;

        deconvolve(grad, expfac, expfac_bcc, Oar) ;

        for (n = 0; n < Nptotal; n++)
	     *(grad+n) *= 2. ;
}


	/************************************************
	minimal_b_funct() skips the gradient calculation; 
	it is called by back_main() to get things going, 
	and after symmetrization, to see how far Back 
	strayed.   
	************************************************/

void	minimal_b_funct(real *fval)
{

	int	n ;
	real	ressq = 0 ;
        void	convolve(real *, real *, COMPLEX *, COMPLEX *) ;
	COMPLEX	*pOar, *pFc0 ;
	real	*wt ;
	char	*pmaskfc ;

	convolve(nump, expfac, expfac_bcc, Oar) ;

	pOar = Oar ;
	pFc0  = Fc0  ;
	wt = hkl_weight ;
	pmaskfc = maskfc ;

	for (n = 0; n < Nhkl; n++, pOar++, pFc0++, wt++, pmaskfc++) {

	   pOar->re *=  *pmaskfc ;
	   pOar->im *=  *pmaskfc ;

	   pOar->re -=  pFc0->re ;
	   pOar->im -=  pFc0->im ;

	   pOar->re *=  *wt ;
	   pOar->im *=  *wt ;

	   ressq +=  pOar->re * pOar->re + pOar->im * pOar->im ;
        } 

	*fval = .5 * ressq ;

	if (costfile_option)
	     fprintf(fp_cost, "Residual = %10.4e\n", *fval) ;

	fflush(fp_cost) ;
	return ;
}

void	recalc_fp_fpp(int md, real *csol, COMPLEX *Fhyd, real fp, real fpp)
{
        void	convolve(real *, real *, COMPLEX *, COMPLEX *) ;

	int	n, r ;
	real	dsq ;
	real	wsq ;
	real	Osign ;
	COMPLEX	Fhop ;
	real	A, B, C, D, E, F ;
	real	fdenom ;
	real	dfprime, df2prime ;

	COMPLEX	*pOar, *pR, *pFhyd ;	
	real	*pF, *wt ;
	char	*pmask ;

        /***********************************************************************
	For the applicable resolution, r, do an FFT on csol and convolve the 
	result with the appropriate exponential factor.
        ***********************************************************************/

        r = useset[md] ;

	convolve(csol, rexpfac[r], rexpfac_bcc[r], Oar) ;

        /* Apply mask    */

        for (n = 0, pmask=maskfo, pOar=Oar; n < Nhkl; n++, pOar++, pmask++) {
             pOar->re *= *pmask ;
             pOar->im *= *pmask ;
        } 

        /***********************************************************************
         Replace O by R+O in Oar; calculate |R+O| and hold it in dsq;
	 Put |R+O| - Fobs in resid and collect the weighted residuals.
	 Re-use the "hydrogen" factor to solve:

		df' * A - df'' * B = C
		df' * B - df'' * E = F

	We use Fhop = Fhyd * (Oar*) (complex variables!) where (Oar*) stands
	for the complex conjugate of Oar.
        ***********************************************************************/
	 

	A = B = C = D = E = F = 0 ;
	pFhyd = Fhyd ;
	pOar = Oar ;
	pR = R + md*Nhkl ;
	pF = Full + md*Nhkl ;
	wt = hkl_weight + md*Nhkl ;

	for (n = 0; n < Nhkl; n++, pOar++, pR++, pF++, wt++, pFhyd++) {
	     pOar->re += pR->re ;
	     pOar->im += pR->im ;
	     dsq = sqrt(pOar->re * pOar->re + pOar->im * pOar->im) ;

	     if (dsq < EPS_PHI) {
		  Osign = (pOar->re >= 0) ? 1. : -1. ;
		  pOar->re = Osign ;
		  pOar->im = 0 ;
		  dsq = 1. ;
             }

	     wsq = *wt * *wt ;

	     Fhop.re =  pFhyd->re * pOar->re + pFhyd->im * pOar->im ;
	     Fhop.im = -pFhyd->re * pOar->im + pFhyd->im * pOar->re ;

	     A += wsq * Fhop.re * Fhop.re / (dsq * dsq) ;
	     B += wsq * Fhop.re * Fhop.im / (dsq * dsq) ;
	     C += wsq * Fhop.re  * (*pF - dsq) / dsq ;
	     D += wsq * Fhop.re * Fhop.im / (dsq * dsq) ;
	     E += wsq * Fhop.im * Fhop.im / (dsq * dsq) ;
	     F += wsq * Fhop.im  * (*pF - dsq) / dsq ;
        } 

	if ((fdenom = A*E - B*D) > 0) {
	     dfprime  = (C*E - B*F) / fdenom ;
	     df2prime = (C*D - A*F) / fdenom ;
        }
	else {
	     EdenWarning("Vanishing denominator in delta f' f''") ;
	     dfprime  = 0 ;
	     df2prime = 0 ;
        }

        /* Reports */

	sprintf(message, 
	"\nder #%d: delta f' f'' = %g %g", md, dfprime, df2prime) ;
        printTwice(message) ;
	sprintf(message, 
	"corresponding f' f'' = %g %g", fp + dfprime, fpp + df2prime) ;
        printTwice(message) ;
}

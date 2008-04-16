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

                                CONCOSTS.C

  Title:        Source code for constraint cost functions 
  Author:       Hanna Szoke
  Date:         6/10/00  (separated from cost.c)
  Function:     This module contains functions used by Solve and
		called from e_funct() in cost.c, for calculating
		the various constraint cost functions.
			e_funct_np()
			e_funct_phasext()
			e_funct_cs()
			e_funct_sayre()
			e_funct_singlets()
			e_funct_triplets()
		and lower-level functions called by them.

*******************************************************************************/

#include "eden.h"		
#include "eden_Np.h"		/* Np space arrays also in back.c, qshare.c */
#include "remap.h"		/* for close vicinity procedures */
#include "cost.h"		/* for cost function calculation */

static	real	*phext_expfac ;       	/* exp factors for phase extension */
static	COMPLEX *phext_expfac_bcc ;   	/* exp factors for phase extension BCC*/
static	char	*phmask ;		/* hkl space mask for phase extension */

          /* data structures for singlet and triplet cost functions */

static	int	Nslist ;		/* length of list of singlet info */
static	int	*singlet_n ;		/* index into (hkl) for a singlet phase */
static	real	*singlet_phase ;	/* singlet target phase */
static	real	*singlet_sig ;		/* sigma associated with singlet target phase */
static	int	Ntlist ;		/* length of list of triplet info */
static	int	*triplet_n ;		/* index into (hkl) for a triplet (hkl) */
static	real	*triplet_phase ;	/* triplet target phase */
static	real	*triplet_sig ;		/* sigma associated with triplet phase */
static	real	*triplet_fac ;		/* +/-1 associated with triplet phase */
static	COMPLEX	*Gst = 0 ;		/* gradient array for both s and t phases */

static	void	prepare_phase_ext(int) ;	
static	real	use_st_info(int, int, int, int *) ;

#define EPSI_AMP	1.e-5

real e_funct_np(real *grad, real lamcoef, int k)

/**	grad 	(returned) gradient array
	lamcoef  relative weight * con_coef
	k 	 target index ***/
{
	int	j ;
	real	res_np = 0 ;
	real	*pwtp, *ptotp, *ptarg, *pgrad ; 

        ptotp = totp ;
	pgrad = grad ;
        pwtp = weightp + Nptotal*k ;
	ptarg = target + Nptotal*k ;

	for (j = 0; j < Nptotal; j++, pgrad++, pwtp++, ptotp++, ptarg++) {

             res_np += *pwtp * *pwtp * (*ptotp - *ptarg) * (*ptotp - *ptarg) ;
	     *pgrad += lamcoef * *pwtp * *pwtp * (*ptotp - *ptarg) ;

        }

        res_np *= lamcoef / 2. ;


        /* Reports */

        if (costfile_option) 
	     fprintf(fp_cost, ", Np: %10.4e", res_np) ;

	return (res_np) ;
}

real e_funct_phasext(real *grad, real lamcoef, int k)
{
     void	apply_weight(real *, real *) ;
     void	convolve(real *, real *, COMPLEX *, COMPLEX *) ;
     void	deconvolve(real *, real *, COMPLEX *, COMPLEX *) ;

     int	j, n ;
     real	pe_cost = 0 ;

     real	*ptotp, *ptarg, *pgrad ; /* ptrs for rapid looping in Np space*/
     char	*pmask ;		 /* ptr for Nhkl space */
     COMPLEX	*pOar, *pFc0 ;		 /* ptrs for Nhkl space */
     real	*wt ;			 /* ptr for hkl_weight */
     static	int	ph_first = TRUE ;

     if (ph_first) {
	  prepare_phase_ext(k) ;
	  ph_first = FALSE ;
     }

     /************************************************************
     Step 1: apply weights to current solution in order to cut out
     whatever is extraneous to the lowres target.  Here, the target 
     array, ptarg = target + Nptotal*k, is used as a spare array.
     Step 2: prepare its low-res (hkl)-space counterpart.
     Step 3: apply a Back-like cost calculation.
     ************************************************************/

     ptarg=target + Nptotal*k ;

     for (n = 0, ptotp = totp; n < Nptotal; n++, ptarg++, ptotp++)
	  *ptarg = *ptotp ;

     ptarg=target + Nptotal*k ;

     apply_weight(ptarg, weightp + Nptotal*k) ;

     convolve(ptarg, phext_expfac, phext_expfac_bcc, Oar) ;

     	/************************************************************
	Now we calculate cost function much as in b_funct().  Again, 
	ptarg is used as a spare array for the gradient contribution.  
	The gradient terms must have the Np-space weights applied 
	before accumulating them into the full gradient.
	************************************************************/

     pOar = Oar ;
     pFc0  = Fc0 ;
     pmask = phmask ;
     wt = hkl_weight ;

     for (n = 0; n < Nhkl; n++, pOar++, pFc0++, pmask++, wt++) {

        pOar->re -=  pFc0->re ;
        pOar->im -=  pFc0->im ;

        pOar->re *=  *pmask  ;
        pOar->im *=  *pmask  ;

        pOar->re *=  *wt  ;
        pOar->im *=  *wt  ;

        pe_cost +=  pOar->re * pOar->re + pOar->im * pOar->im ;

     }

     for (j = 0; j < Nptotal; j++)
	     *(ptarg+j) = 0 ;

     deconvolve(ptarg, expfac, expfac_bcc, Oar) ;

     apply_weight(ptarg, weightp + Nptotal*k) ;

     for (j = 0, pgrad = grad; j < Nptotal; j++, pgrad++, ptarg++) 
	  *pgrad += lamcoef * *ptarg ;

     pe_cost *= 0.5 * lamcoef ;

        /* Reports */

     if (costfile_option) 
          fprintf(fp_cost, ", Phasext: %10.4e", pe_cost) ;

     return (pe_cost) ;
}

real e_funct_cs(real *csol, real *grad, real lamcoef)
{
	void	do_cs_symmetrization(real *, real *) ;   /* in crysutil.c */
	int	j ;
	real	res_cs = 0 ;
	real	*pcsol, *pgrad ;
	real	*pave ;


	/* Compute the cs (crystallographic symmetry) cost function. */

	pave = extra ;
        do_cs_symmetrization(csol, pave) ;
	pcsol = csol ;

	for (j = 0; j < Nptotal; j++, pcsol++, pave++)
             res_cs += (*pcsol - *pave) * (*pcsol - *pave) ;
        res_cs *= lamcoef / 2. ;

      /* Add in the gradient contribution corresponding to the cs constraints */

	pgrad = grad ;
	pcsol = csol ;
	pave = extra ;

	for (j = 0; j < Nptotal; j++, pgrad++, pcsol++, pave++)
	     *pgrad +=  lamcoef * (*pcsol - *pave) ;

        /* Reports */

        if (costfile_option)
	     fprintf(fp_cost, ", Ncs: %10.4e", res_cs) ;
	   
	return (res_cs) ;
}

real e_funct_sayre(real *grad, real lamcoef, real addend) 
/* 	addend ;	cost_addend[k]	*/
{
	real	*pgrad ;
	real	*ptotp ;

        real	res_sayre = 0 ;
	int	n, m ;
	int	close_ind ;
	real	gex ;
	real	contrib ;
	real	mean_evox = 0 ;
	int	*bpoff ;
	int	*nindexlist ;

        /***********************************************************************
        Calculate the global mean, to be subtracted from each e/vox value;
	earlier code had the following "dynamic" mean:

	ptotp = totp ;

	for (n = 0; n < Nptotal; n++, ptotp++) 
	     mean_evox += *ptotp ;
        
	mean_evox /= Nptotal ;
	
	but this would seem to be incorrect, if the # of recovered electrons 
	has risen strongly with respect to the (known) F000.
        ***********************************************************************/

	mean_evox = F000 / Nptotal ;

        /***********************************************************************
	High resolution (Sayre's equation) cost function contribution. 
	We use poff for accessing the appropriate indexlist entries, for  
	direct (fast) indexing into the close vicinity of each n.
        Note that we subtract the global mean from each e/vox value.
        ***********************************************************************/

	ptotp = totp ;
	pgrad = grad ;
	bpoff = poff ;

	for (n = 0; n < Nptotal; n++, bpoff++, ptotp++, pgrad++) {
	   
	     for (m = 0, nindexlist = indexlist + *bpoff, contrib = 0; 
		  m < Nb; m++, nindexlist++) {

		  gex = *(gexp + m) ;
                  close_ind = (n + *nindexlist) ;

		  contrib += gex * (*(totp + close_ind) - mean_evox) ;
                  
             }

	     res_sayre += (*ptotp - mean_evox) * contrib ;
	     *pgrad += -lamcoef * contrib ;

        }

	res_sayre -= addend ;
        res_sayre *= (-0.5) * lamcoef ;

        /* Reports */

	if (costfile_option) 
             fprintf(fp_cost, ", sayre: %10.4e", res_sayre) ;

        return(res_sayre) ; 
}

void prepare_phase_ext(int tnum)
{
     real     *setupexpfac(real) ;
     COMPLEX    *setupexpfac_bcc(real) ;
     void	apply_weight(real *, real *) ;
     void	convolve(real *, real *, COMPLEX *, COMPLEX *) ;

     real	*pwtp, *ptarg ;
     real	delta ;
     real	phdr ;			/* low resolution counterpart to dr */
     int        n, h, k, l ;

  	/************************************************************
        Set up phext_expfac(Nhkl), defined as: exp (-A|k*D^2|) *
	exp[2pi*i* k.dr]
	where D is the low-res counterpart of dr (input res), 
        A is PISQ*eta, k stands for (h,k,l). Set up the comparable 
	COMPLEX array for body-centered grid,
	************************************************************/

     phdr = dr * phase_ext_res / input_res ;
     delta = PISQ * eta * phdr*phdr ;

     phext_expfac = setupexpfac(delta) ;

     if (grid_type == BODY_CENTERED)
          phext_expfac_bcc = setupexpfac_bcc(delta) ;

  	/************************************************************
	Allocate and create Fc0.  Once this is completed, there is no
	further use for the target array in question; it will be used
	in a dummy capacity in e_funct_phasext().
	************************************************************/

     Fc0 = (COMPLEX *) e_malloc(Nhkl*sizeof(COMPLEX), "prepare_phase_ext") ; 

     ptarg = target + Nptotal*tnum ;
     pwtp  = weightp + Nptotal*tnum ;

     apply_weight(ptarg, pwtp) ;

     convolve(ptarg, expfac, expfac_bcc, Fc0) ;

  	/************************************************************
	Create phmask.  (Do we really need it?)
	************************************************************/

     phmask = (char *) e_malloc(Nhkl*sizeof(char), "prepare_phase_ext") ; 

     for (n = 0; n < Nhkl; n++) {

          fetch_hkl(n, &h, &k, &l) ;

	  *(phmask + n) = (nusq(h, k, l) < limit_rat) ? 1 : 0 ;
	  /*** *(phmask + n) = ((Fc0 + n)->re > 0) ? 1 : 0 ; ***/
     }

     return ;
}

				/***********************************************
				Prepare singlet info for e_funct_singlet().
				***********************************************/
void	prepare_singlets(char *filename)
{
     FILE	*fp ;
     real	cleanPhase(real) ;

     real	*checks ;
     int	filelength() ;
     int	s_h, s_k, s_l ;
     float	s_phase, s_sigma ;
     real	off ;
     real	*mat ;
     real	phase, factor ;
     int	newh, newk, newl ;
     int	s, n, p ;
				/***********************************************
				Read a file of singlet phases.  Expected format:
				h1 k1 l1 phase sigma
				Note: no "INDEX" etc.
				***********************************************/

     Nslist = Ncs*filelength(filename) ;
     sprintf(message, "singlet arrays") ;

     singlet_n     = (int *)    e_malloc(Nslist*sizeof(int), message) ; 
     singlet_phase = (real *) e_malloc(Nslist*sizeof(real), message) ; 
     singlet_sig   = (real *) e_malloc(Nslist*sizeof(real), message) ; 
     checks        = (real *) e_malloc(Nhkl  *sizeof(real), message) ;

     if (Gst == 0)
	  Gst = (COMPLEX *) e_malloc(Nhkl*sizeof(COMPLEX), message) ;

     fp = fopen(filename, "r");
     if (fp == NULL)
     {
          sprintf(message, "Cannot open %s", filename);
          EdenError(message) ;
     }

     for (n = 0; n < Nslist; n++)
	  *(singlet_n + n) = 0 ;

     for (n = 0; n < Nhkl; n++)
	  *(checks + n) = -360 ;	/* -360 is definitely illegal! */

     fprintf(fp_log, "\nReading singlet information from %s\n", filename) ;

     s = 0 ;

     while (fgets(nextline, MAXSTRING, fp) != NULL)  {

        s_sigma = 0 ;

        sscanf(nextline, "%d %d %d  %g %g", 
               &s_h, &s_k, &s_l, &s_phase, &s_sigma) ;

        if (s_sigma == 0) {
	     sprintf(message, "Missing sigmas in %s!", filename) ;
	     EdenError(message) ;
        }

        for (p = 0, mat = matop; p < Ncs; p++, mat += MEL) {

              apply_symop_hkl(s_h, s_k, s_l, mat, 
			      &newh, &newk, &newl, &off) ;

	      factor = use_st_info(newh, newk, newl, &n) ;
	      phase = (s_phase - off*RTOD) * factor ;
	      phase = cleanPhase(phase) ;

	      if (!*(maskfo + n)) {
		 sprintf(message, 
"Singlet (%d %d %d) from input (%d %d %d) corresponds to point w/o fobs info!",
		   newh, newk, newl, s_h, s_k, s_l) ;
		 EdenWarning(message) ;
              }
	      else if (*(checks + n) > -360) {
		 if (fabs(*(checks + n) - phase) > EPS_PHI) {
		      printTwice("Mismatching singlet phases!") ;
                      sprintf(message, 
		      "s=%d, p=%d, n=%d, phase=%g, old phase=%g, delta=%g\n", 
                      s, p, n, phase, *(checks+n), phase-*(checks+n)) ;
                      fprintf(fp_log, message) ;
                 }
              } 
		/***********************************************
		Collect legal info.
		***********************************************/

              else {
                 *(singlet_n + s) = n ;
                 *(singlet_phase + s) = *(checks + n) = phase ;
                 *(singlet_sig + s) = s_sigma ;
	         s++ ;
	      }
        }
     }
     Nslist = s ;	/* correction for duplications and empty lines */

     sprintf(message, "Length of singlet data after expansion to P1 is %d", 
	     Nslist);
     printTwice(message) ;
     e_free(Nhkl*sizeof(real), checks) ;
}

				/***********************************************
				Prepare triplet info for e_funct_triplet().
				***********************************************/
void	prepare_triplets(char *filename)
{
     FILE	*fp ;
     real	cleanPhase(real) ;

     real	*checkt ;
     int	filelength() ;
     int	t_h[3], t_k[3], t_l[3] ;
     float	t_phase, t_sigma ;
     real	off ;
     real	*mat ;
     real	factor ;
     int	newh, newk, newl ;
     int	t, n, p, q ;
     int	h0, k0, l0, h1, k1, l1, h2, k2, l2 ;
     int	f0, f1, f2 ;
     int	legal_triplet ;

		/***********************************************
		Read a file of triplet phases.  Expected format:
		h1 k1 l1 h2 k2 l2 h3 k3 l3 phase sigma
		where sum(hi)=sum(ki)=sum(li) = 0.
		***********************************************/

     Ntlist = 3*Ncs*filelength(filename) ;
     sprintf(message, "triplet arrays") ;

     triplet_n     = (int *)    e_malloc(Ntlist*sizeof(int), message) ; 
     triplet_phase = (real *) e_malloc(Ntlist*sizeof(real), message) ; 
     triplet_sig   = (real *) e_malloc(Ntlist*sizeof(real), message) ; 
     triplet_fac   = (real *) e_malloc(Ntlist*sizeof(real), message) ; 
     checkt        = (real *) e_malloc(Nhkl  *sizeof(real), message) ;

     if (Gst == 0)
	  Gst = (COMPLEX *) e_malloc(Nhkl*sizeof(COMPLEX), message) ;

     fp = fopen(filename, "r");
     if (fp == NULL)
     {
          sprintf(message, "Cannot open %s", filename);
          EdenError(message) ;
     }

     for (n = 0; n < Nhkl; n++)
	  *(checkt + n) = -1 ;

     fprintf(fp_log, "\nReading triplet information from %s\n", filename) ;

     t = 0 ;

     while (fgets(nextline, MAXSTRING, fp) != NULL)  
      if ((int)strlen(nextline) > 1) {

        t_sigma = 0 ;

        sscanf(nextline, "%d %d %d %d %d %d %d %d %d  %g %g", 
               t_h, t_k, t_l, t_h+1, t_k+1, t_l+1, 
	       t_h+2, t_k+2, t_l+2, &t_phase, &t_sigma) ;

        if (t_sigma == 0) {
	     sprintf(message, "Missing sigmas in %s!", filename) ;
	     EdenError(message) ;
        }

        for (p = 0, mat = matop; p < Ncs; p++, mat += MEL) { 

           legal_triplet = TRUE ;

           for (q = 0; q < 3; q++) {

              apply_symop_hkl(t_h[q], t_k[q], t_l[q], mat, 
			      &newh, &newk, &newl, &off) ;
	      factor = use_st_info(newh, newk, newl, &n) ;

	      if (!*(maskfo + n)) {
		 sprintf(message, 
"Triplet (%d %d %d) from input (%d %d %d) corresponds to point w/o fobs info!\n",
		   newh, newk, newl, t_h[q], t_k[q], t_l[q]) ;
		 fprintf(fp_log, message) ;
		 legal_triplet = FALSE ;
              }
           }

           if (legal_triplet) {

		/***********************************************
		Collect info; note the following:

		1) t_sigma and t_phase are the same for all
		   symmetry-related (h,k,l)'s.  It may seem that
		   one should apply the sum of the offsets to
		   t_phase, but that sum is 0 (to within e-14).

		2) factor is +1 or -1, -1 if the (hkl) go into
		   the negative half-ellipsoid under the
		   symmetry operation.  

		3) n, of course, differs for each (hkl).
		***********************************************/

              for (q = 0; q < 3; q++) {

                 apply_symop_hkl(t_h[q], t_k[q], t_l[q], mat, 
			         &newh, &newk, &newl, &off) ;

	         factor = use_st_info(newh, newk, newl, &n) ;
	         *(triplet_n + t + q) = n ;
	         *(triplet_fac + t + q) = factor ;
                 *(triplet_sig + t + q) = t_sigma ;
                 *(triplet_phase + t + q) = cleanPhase(t_phase) ;

              }

		/***********************************************
       		Do some checks 
		***********************************************/

              fetch_hkl(*(triplet_n + t + 0), &h0, &k0, &l0) ;
              fetch_hkl(*(triplet_n + t + 1), &h1, &k1, &l1) ;
              fetch_hkl(*(triplet_n + t + 2), &h2, &k2, &l2) ;
	      f0 = (int) *(triplet_fac + t + 0) ;
	      f1 = (int) *(triplet_fac + t + 1) ;
	      f2 = (int) *(triplet_fac + t + 2) ;
              
              if (((h0*f0 + h1*f1 + h2*f2) != 0) || 
	          ((k0*f0 + k1*f1 + k2*f2) != 0) || 
	          ((l0*f0 + l1*f1 + l2*f2) != 0)) {
	          sprintf(message, 
		  "Illegal triplet invariant (%s) for matrix op. # %d  in %s!",
		          nextline, p, filename) ;
	          EdenError(message) ;
              }

              t += 3 ;
	   }
        }
     } 
     Ntlist = t ;	/* correction for empty lines */

     sprintf(message, 
     "Length of valid triplet data after expansion to P1 is %d", Ntlist);
     printTwice(message) ;
}

real	use_st_info(int h, int k, int l, int *nvalue) 
{
      int	assemble(int, int, int) ;
      real	factor ;

      if (h > 0) {
         *nvalue = assemble(h, k, l) ;
         factor = 1. ;
      }
      else if ((h == 0) && (k > 0)){
         *nvalue = assemble(h, k, l) ;
         factor = 1. ;
      }
      else if ((h == 0) && (k == 0) && (l > 0)){
         *nvalue = assemble(h, k, l) ;
         factor = 1. ;
      }
      else {
         *nvalue = assemble(-h, -k, -l) ;
         factor = -1. ;
      }

      return (factor) ;
}

real	e_funct_singlet(real *grad, real lamcoef)
{
     void	convolve(real *, real *, COMPLEX *, COMPLEX *) ;
     void	deconvolve(real *, real *, COMPLEX *, COMPLEX *) ;
     void	Cmaskit(COMPLEX *, char *, int) ;

     int	i, n ;
     int	s ;

     COMPLEX	*pGst ;	

     real	wt ;
     real	amp0, amplitude, phase ;
     real	numer, sdenom ;
     real	this_phase, this_sigma, dphase ;
     real	sigsq, wsq ;
     real	mult, gmult ;
     real	singlet_cost = 0 ;
     real	gnumer ;
     real	s_ratio ;

        /***********************************************************************
	Do an FFT on totp and convolve the result with the exponential factor.
        Apply mask.
        ***********************************************************************/

        convolve(totp, expfac, expfac_bcc, Oar) ;

        Cmaskit(Oar, maskfo, Nhkl) ;

        /***********************************************************************
	The following looks somewhat like the guts of e_funct_hkl().
	It is derived from AS's write-up of Aug 6, 1999 with corrections from 
	Jan 24, 2000.  Where applicable, I include as comments the previous
	version of the code.
        ***********************************************************************/
	 
         amp0 = Oar->re ;
	 numer = sdenom = 0 ;

	 for(i = 0, pGst = Gst; i < Nhkl; i++, pGst++) {
	      pGst->re = 0 ;
	      pGst->im = 0 ;
         }

	 for (s = 0; s < Nslist; s++) 	{

	      this_phase = *(singlet_phase + s) ;
	      this_sigma = *(singlet_sig + s) ;
	      sigsq = this_sigma * this_sigma ;
	      n = *(singlet_n + s) ;
	      wt = *(hkl_weight + n) ;
	      wsq = wt * wt ;
	      s_ratio = wsq / sigsq ;

	      getAmpPhase(Oar+n, amp0, &amplitude, &phase) ;

	      dphase = phase - this_phase ;

	      gmult = *(Full + n) ;
	      mult = gmult * amplitude ;

	      numer += s_ratio * mult * (1. - cos(DTOR*dphase)) ;
	      sdenom += s_ratio ;
	      gnumer = -s_ratio * gmult * sin(DTOR*dphase) ;

	      (Gst + n)->re += gnumer * sin(DTOR*phase) ;
	      (Gst + n)->im += -gnumer * cos(DTOR*phase) ;
         }			

         sdenom /= Nslist ;
	 singlet_cost = lamcoef * numer / sdenom ;

	 for(i = 0, pGst = Gst; i < Nhkl; i++, pGst++) {
	      pGst->re *= lamcoef/sdenom ;
	      pGst->im *= lamcoef/sdenom ;
         }

         deconvolve(grad, expfac, expfac_bcc, Gst) ;

        /* Reports */

	if (costfile_option) 
             fprintf(fp_cost, ", singlet: %10.4e", singlet_cost) ;

        return(singlet_cost) ; 
}

real	e_funct_triplet(real *grad, real lamcoef)
{
     void	convolve(real *, real *, COMPLEX *, COMPLEX *) ;
     void	deconvolve(real *, real *, COMPLEX *, COMPLEX *) ;
     void	Cmaskit(COMPLEX *, char *, int) ;

     int	i, n ;
     int	t, q ;

     COMPLEX	*pGst ;	

     real	wt ;
     real	amp0, amplitude, phase ;
     real	numer, tdenom ;
     real	this_phase, this_sigma, dphase ;
     real	sigsq, wsq, dsq ;
     real	tm1, tm2, mult ;
     real	triplet_cost = 0 ;
     real	gnumer ;
     real	t_ratio ;
     real	psum, factor ;

        /***********************************************************************
	Do an FFT on totp and convolve the result with the exponential factor.
        Apply mask.
        ***********************************************************************/

        convolve(totp, expfac, expfac_bcc, Oar) ;

        Cmaskit(Oar, maskfo, Nhkl) ;

        /***********************************************************************
	The following looks somewhat like the guts of e_funct_hkl().
	It is derived from AS's write-up of Aug. 6, 1999
	with corrections from Jan 24, 2000.
        ***********************************************************************/
	 
         amp0 = Oar->re ;
	 numer = tdenom = 0 ;

	 for(n = 0, pGst = Gst; n < Nhkl; n++, pGst++) {
	      pGst->re = 0 ;
	      pGst->im = 0 ;
         }

	 for (t = 0; t < Ntlist; t += 3) 	{

	      for (q = 0, psum = 0; q < 3; q++)	{

	          this_sigma = *(triplet_sig + t + q) ;
	          sigsq = this_sigma * this_sigma ;
	          n = *(triplet_n + t + q) ;
	          wt = *(hkl_weight + n) ;
	          wsq = wt * wt ;
	          t_ratio = wsq / sigsq ;

	          getAmpPhase(Oar+n, amp0, &amplitude, &phase) ;

	          this_phase = *(triplet_phase + t + q) ;
		  factor = *(triplet_fac + t + q) ;
		  psum += phase * factor ;
              }

	      dphase = psum - this_phase ;

	      for (q = 0, tm1 = tm2 = 1.; q < 3; q++)	{

	          n = *(triplet_n + t + q) ;

	          tm1 *= ((Oar+n)->re*(Oar+n)->re) + ((Oar+n)->im*(Oar+n)->im) ;
	          tm2 *= *(Full+n) * *(Full+n) ;
              }
	      mult = pow(tm1*tm2, (real)1./(real)6.) ;

	      numer += t_ratio * mult * (1. - cos(DTOR*dphase)) ;
	      tdenom += t_ratio ;
	      gnumer = -t_ratio * mult * sin(DTOR*dphase) ;

	      /**************************************************************
	      At this point, we must deal with the 3 components individually,
	      sorting out real and imaginary parts of the gradient - 
	      e.g. for the 1st case:
	      e^[pi/2 - sum(phi) + phi(2) + phi(3)] =
	      e^[pi/2 - sum(phi) + sum(phi) - phi(1)] = e^[pi/2 - phi(1)].
	      **************************************************************/

	      for (q = 0; q < 3; q++)	{

	          n = *(triplet_n + t + q) ;
		  factor = *(triplet_fac + t + q) ;

	          getAmpPhase(Oar+n, amp0, &amplitude, &phase) ;

	          dsq = sqrt(amplitude*amplitude + EPSI_AMP*EPSI_AMP) ;

	          (Gst + n)->re +=  factor * gnumer * sin(DTOR*phase) / dsq ;
	          (Gst + n)->im += -factor * gnumer * cos(DTOR*phase) / dsq ;
              }
         }			

         tdenom /= Ntlist ;
	 triplet_cost = lamcoef * numer / tdenom ;

	 for (i = 0, pGst = Gst; i < Nhkl; i++, pGst++) {
	      pGst->re *= lamcoef/tdenom ;
	      pGst->im *= lamcoef/tdenom ;
         }

         deconvolve(grad, expfac, expfac_bcc, Gst) ;

        /* Reports */

	if (costfile_option) 
             fprintf(fp_cost, ", triplet: %10.4e", triplet_cost) ;

        return(triplet_cost) ; 
}

real	dphase_triplet()
{
     void	convolve(real *, real *, COMPLEX *, COMPLEX *) ;
     real	cleanPhase(real) ;
     void	Cmaskit(COMPLEX *, char *, int) ;

     int	n ;
     int	t, q ;

     FILE	*fp ;

     real	wt ;
     real	amp0, amplitude, phase ;
     real	weighted_numer, weighted_denom ;
     real	this_phase, this_sigma, dphase ;
     real	sigsq, wsq ;
     real	t_ratio ;
     real	psum ;
     real	t_coef ;
     real	factor ;
     int	h0,k0,l0,h1,k1,l1,h2,k2,l2;
     static	int	which = 0 ;
     char	filename[MAXSTRING] ;

        /***********************************************************************
	Prepare to write a new triplet file of Eden's results. 
        ***********************************************************************/
	 
     sprintf(filename, "eden_trip.%1d", which) ;

     fp = fopen(filename, "w");
     if (fp == NULL)
     {
          sprintf(message, "Cannot open %s", filename) ;
          EdenError(message) ;
     }
     which++ ;

        /***********************************************************************
	Do an FFT on totp and convolve the result with the exponential factor.
        Apply mask.
        ***********************************************************************/

        convolve(totp, expfac, expfac_bcc, Oar) ;

        Cmaskit(Oar, maskfo, Nhkl) ;

        /***********************************************************************
	The following is derived from AS's write-up of Jan 11, 2000.
        ***********************************************************************/
	 
         amp0 = Oar->re ;
	 weighted_numer = weighted_denom = 0 ;

	 for (t = 0; t < Ntlist; t += 3) 	{

	      this_sigma = *(triplet_sig + t) ;
	      sigsq = this_sigma * this_sigma ;
	      n = *(triplet_n + t) ;
	      wt = *(hkl_weight + n) ;
	      wsq = wt * wt ;
	      t_ratio = wsq / sigsq ;
	      t_coef = 1. ;

	      for (q = 0, psum = 0; q < 3; q++)	{

	          n = *(triplet_n + t + q) ;
		  factor = *(triplet_fac + t + q) ;

	          getAmpPhase(Oar+n, amp0, &amplitude, &phase) ;

		  psum += phase * factor ;
		  t_coef *= amplitude ;

	          this_phase = *(triplet_phase + t + q) ;
              }

	      dphase = psum - this_phase ;
	      dphase = cleanPhase(dphase) ;

	      if (t_coef > 0) {

	              weighted_numer += t_ratio * t_coef * fabs(dphase) ;
	              weighted_denom += t_ratio * t_coef ;

                      fetch_hkl(*(triplet_n + t + 0), &h0, &k0, &l0) ;
                      fetch_hkl(*(triplet_n + t + 1), &h1, &k1, &l1) ; 
		      fetch_hkl(*(triplet_n + t + 2), &h2, &k2, &l2) ;
	              fprintf(fp, "%d %d %d %d %d %d %d %d %d %g %g %g\n",
	              h0,k0,l0,h1,k1,l1,h2,k2,l2,
	              psum, dphase, t_coef) ;
	      }

         }			
	 fclose(fp) ;

        return(weighted_numer / weighted_denom) ; 
}

real	dphase_singlet()
{
     void	convolve(real *, real *, COMPLEX *, COMPLEX *) ;
     real	cleanPhase(real) ;
     void	Cmaskit(COMPLEX *, char *, int) ;

     FILE	*fp ;
     real	wt ;
     real	amp0, amplitude, phase ;
     real	weighted_numer, weighted_denom ;
     real	this_phase, this_sigma, dphase ;
     real	sigsq, wsq ;
     real	s_ratio ;
     real	s_coef ;
     int	h,k,l;
     static	int	which = 0 ;
     char	filename[MAXSTRING] ;
     int	n, s ;

        /***********************************************************************
	Prepare to write a new singlet file of Eden's results. 
        ***********************************************************************/
	 
     sprintf(filename, "eden_sing.%1d", which) ;

     fp = fopen(filename, "w");
     if (fp == NULL)
     {
          sprintf(message, "Cannot open %s", filename) ;
          EdenError(message) ;
     }
     which++ ;

        /***********************************************************************
	Do an FFT on totp and convolve the result with the exponential factor.
        Apply mask.
        ***********************************************************************/

        convolve(totp, expfac, expfac_bcc, Oar) ;

        Cmaskit(Oar, maskfo, Nhkl) ;

        /***********************************************************************
	The following is derived from AS's write-up of Jan 11, 2000.
        ***********************************************************************/
	 
         amp0 = Oar->re ;
	 weighted_numer = weighted_denom = 0 ;

	 for (s = 0; s < Nslist; s++) 	{

	      this_phase = *(singlet_phase + s) ;
	      this_sigma = *(singlet_sig + s) ;
	      sigsq = this_sigma * this_sigma ;
	      n = *(singlet_n + s) ;
	      wt = *(hkl_weight + n) ;
	      wsq = wt * wt ;
	      s_ratio = wsq / sigsq ;

	      getAmpPhase(Oar+n, amp0, &amplitude, &phase) ;

	      s_coef = amplitude ;
	      dphase = cleanPhase(phase - this_phase) ;
	      dphase = fabs(dphase) ;

	      if (s_coef > 0) {

	              weighted_numer += s_ratio * s_coef * dphase ;
	              weighted_denom += s_ratio * s_coef ;

                      fetch_hkl(*(singlet_n + s), &h, &k, &l) ;
	              fprintf(fp, "%d %d %d %g %g %g\n",
	              h,k,l, phase, dphase, s_coef) ;
	      }

         }			
	 fclose(fp) ;

        return(weighted_numer / weighted_denom) ; 
}
void	apply_weight(real *tar, real *wt)
{
     int	n ;

     for (n = 0; n < Nptotal; n++, tar++, wt++)
	  *tar *= *wt ;
}

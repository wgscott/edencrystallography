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

                                HKLUTIL.C

  Title:        Hkl utility package for EDEN.
  Author:       Hanna Szoke
  Date:         2/21/95
  Function:     This package contains functions that handle (hkl) factors 
		for all codes that need them, except functions that read 
		structure factors.

		They include:

		write functions -
			writefobs()
			writefcalc()
			writefcalc_non0()
			discards()
		low-level functions handling data -
			getAmpPhase()
			cleanPhase()
			fetch_hkl()
			fetch_twin()
			assemble()
			setDsqLimit()
			checkF000()
			scalefobs()
		R factor functions -
			get1Rfac()
			get1Rfac_weighted()
		Table Header functions -
			hklTabHead()
			shellHeader()
			shellHeader_R()
			shellHeader_F()
		Mask-related functions -
			get_P1mask()
			prepare_unique_hklmask()
			check_others()		(internal to hklutil.c)
			set_centrics()
			set_cen_phases()
			Cmaskit()

2/22/2001 -- 	Disallowing +/- Nk/2 and Nl/2 in prepare_unique hklmask().
*******************************************************************************/

#include	"util.h"
#include	"symmetry.h"
#include	"dims.h"
#include	"mir.h"		/* for NM */

#define AMP_EPSI 1.e-7		/* for eliminating round-off in amplitudes */
#define	EPSI_PHASE 1.e-3	/* for rounding phases very close to 0 or 180 */

real	dsq_limit ;
static	char	Arange[NUMSEGS][MAXSTRING] ;
static	char	*unique_hkl ;

static	void	check_others(char *, int, int, int, int) ;

void writefobs(char *filename, real *F, char *mask, real *sarray)          

	/* We write the unique reflections for which the mask is set.  */

{
     FILE	*fpout ;
     int	Nrefl, n, h, k, l ;
     char	*pmask ;
     real	*pF, *psig ;

     fpout = fopen(filename, "w") ;
     if (fpout == NULL)
     {
          sprintf(message, "Cannot open %s", filename);
          EdenError(message) ;
     }
     pmask = mask ;
     for (Nrefl = n = 0; n < Nhkl; n++, pmask++) 
	     Nrefl += *pmask ;

     /***********************************************************/
     /*		write header					*/
     /***********************************************************/
     fprintf(fpout, "NREFlection=   %d\n", Nrefl) ;
     if (anom_flag)
          fprintf(fpout, "ANOMalous=TRUE\n") ;
     else
          fprintf(fpout, "ANOMalous=FALSE\n") ;
     
     fprintf(fpout, 
	     "DECLare NAME=FOBS       DOMAin=RECIprocal,  TYPE=REAL END\n") ;

     pmask = mask ;
     pF = F ;
     psig = sarray ;

     if (useSig) {
          fprintf(fpout, 
		 "DECLare NAME=SIGMA      DOMAin=RECIprocal,  TYPE=REAL END\n") ;
          for (n = 0; n < Nhkl; n++, F++, mask++, psig++) 

	       if (*mask) {
	  
                    fetch_hkl(n, &h, &k, &l) ;
						    
	            fprintf(fpout, 
			 "INDEX   %d   %d   %d   FOBS   %g   SIGMA   %g\n", 
	  	              h, k, l, *F, *psig) ;
               }
     }
     else {
          for (n = 0; n < Nhkl; n++, F++, mask++) 

	       if (*mask) {
	  
                    fetch_hkl(n, &h, &k, &l) ;
						    
	            fprintf(fpout, "INDEX   %d   %d   %d   FOBS   %g\n", 
	  	              h, k, l, *F) ;
               }
     }
     fclose(fpout) ;
}

void writefcalc(char *filename, COMPLEX *R, char *mask)          
{
     FILE	*fpout ;
     real	amplitude ;
     real	phase ;
     real	amp0 ;
     int	Nrefl, n, h, k, l ;
     char	*pmask ;

     fpout = fopen(filename, "w") ;
     if (fpout == NULL)
     {
          sprintf(message, "Cannot open %s", filename);
          EdenError(message) ;
     }

     pmask = mask ;
     for (Nrefl = n = 0; n < Nhkl; n++, pmask++) 
	     Nrefl += *pmask ;

     /***********************************************************/
     /*		write header					*/
     /***********************************************************/
     fprintf(fpout, "NREFlection=   %d\n", Nrefl) ;
     if (anom_flag)
          fprintf(fpout, "ANOMalous=TRUE\n") ;
     else
          fprintf(fpout, "ANOMalous=FALSE\n") ;

     fprintf(fpout, 
	     "DECLare NAME=FCALC       DOMAin=RECIprocal,  TYPE=COMP END\n") ;

     amp0 = R->re ;
     pmask = mask ;

     for (n = 0; n < Nhkl; n++, pmask++) {

	  if (*pmask) {

               fetch_hkl(n, &h, &k, &l) ;
	       getAmpPhase(R+n, amp0, &amplitude, &phase) ;

	       fprintf(fpout, "INDEX   %d   %d   %d   FCALC   %g   %g\n", 
	  	       h, k, l, amplitude, phase) ;
          }
     }
     fclose(fpout) ;
}
void writefcalc_non0(char *filename, COMPLEX *R, char *mask)          
{
     FILE	*fpout ;
     real	amplitude ;
     real	phase ;
     real	amp0 ;
     int	Nrefl, n, h, k, l ;
     char	*pmask ;

     fpout = fopen(filename, "w") ;
     if (fpout == NULL)
     {
          sprintf(message, "Cannot open %s", filename);
          EdenError(message) ;
     }

     pmask = mask ;
     amp0 = R->re ;

     for (Nrefl = n = 0; n < Nhkl; n++, pmask++) 
	  if (*pmask) {

	       getAmpPhase(R+n, amp0, &amplitude, &phase) ;

	       if (amplitude > 0)
	            Nrefl += 1 ;
          }

     /***********************************************************/
     /*		write header					*/
     /***********************************************************/
     fprintf(fpout, "NREFlection=   %d\n", Nrefl) ;
     if (anom_flag)
          fprintf(fpout, "ANOMalous=TRUE\n") ;
     else
          fprintf(fpout, "ANOMalous=FALSE\n") ;

     fprintf(fpout, 
	     "DECLare NAME=FCALC       DOMAin=RECIprocal,  TYPE=COMP END\n") ;

     for (n = 0; n < Nhkl; n++, mask++) {

	  if (*mask) {

               fetch_hkl(n, &h, &k, &l) ;
	       getAmpPhase(R+n, amp0, &amplitude, &phase) ;

	       if (amplitude > 0)
		    fprintf(fpout, "INDEX   %d   %d   %d   FCALC   %g   %g\n", 
	  	            h, k, l, amplitude, phase) ;
          }
     }
     fclose(fpout) ;
}

			/************************************************
       			Print out Fc's that are discarded because they 
			are not in the fobs file.  (Print Fc(000) also,
			to prevent the file from looking illegal.)
			************************************************/

void 	discards(COMPLEX *Rfc, char *masko, char *maskc)
{
     FILE	*fp ;
     real     amplitude ;
     real     phase ;
     real     amp0 ;
     int	n ;
     int	h, k, l ;

     if ((fp = fopen("discards.fcalc", "w")) == NULL)
           EdenError("Cannot open discards.fcalc file") ;

     amp0 = Rfc->re ;

     fprintf(fp, "INDEX   0   0   0   FCALC   1.   0.\n") ;	/* dummy */

     for (n = 1; n < Nhkl; n++) {

          fetch_hkl(n, &h, &k, &l) ;

	  if ( *(maskc+n) && !*(masko+n) ) {

	       getAmpPhase(Rfc+n, amp0, &amplitude, &phase) ;
			     
	       fprintf(fp, "INDEX   %d   %d   %d   FCALC   %g   %g\n",
		            h, k, l, amplitude, phase) ;

          }
     }
     fclose(fp) ;
}

	/**************************************************
	Convert a complex (re, im) R into an amplitude and 
	phase; disallow very small amplitudes and put 
	phases into the expected range. 
	**************************************************/

void getAmpPhase(COMPLEX *R, real amp0, real *amplitude, real *phase) 
{
     real	Cabs();
     real	cleanPhase(real) ;

     *amplitude = Cabs(*R) ; 

     if (*amplitude < amp0*AMP_EPSI)  {
          *amplitude = 0 ;
          *phase = 0 ;
     }
     else 
          *phase = cleanPhase(RTOD * atan2(R->im, R->re)) ;

}
real	cleanPhase(real phase)
{
     while (phase > 180.)
          phase -= 360 ;
     while (phase<= -180.)
          phase += 360 ;
     if (fabs(180 - fabs(phase)) < EPSI_PHASE)
          phase = 180 ;
     else if (fabs(phase) < EPSI_PHASE)
          phase = 0 ;

     return(phase) ;
}
	/**************************************************
	Rearrange negative k and l indices.
	Use position: 0, 1, 2, ... N-1,  N,  N+1, ... 2*N-1  
	for k (or l): 0, 1, 2, ... N-1, -N, -N+1, ... -1                
	**************************************************/

void fetch_hkl(int n, int *h, int *k, int *l) 
{

          *h = n % Nh ;
     	  *k = (n/Nh) % Nk ;
     	  *l = n / (Nh*Nk) ;

	  if (*k >= Nk/2)
	       *k -= Nk ;

	  if (*l >= Nl/2)
	       *l -= Nl ;
}

	/**************************************************
	Apply the input twinning matrix to the index n.
	Return the new index, (n_t). and a constant, ct,
	which is 1 for regular, -1 for Friedel pair.
	**************************************************/

int	fetch_twin(int n, int *tmat, int *ct) 
{
     int h, k, l ;
     int h_t, k_t, l_t, n_t ;

     fetch_hkl(n, &h, &k, &l) ;

     h_t = (*(tmat+0)*h + *(tmat+1)*k + *(tmat+2)*l) % Nh ;
     k_t = (*(tmat+3)*h + *(tmat+4)*k + *(tmat+5)*l) % Nk ;
     l_t = (*(tmat+6)*h + *(tmat+7)*k + *(tmat+8)*l) % Nl ;

     if (((h_t > 0)) ||
         ((h_t == 0) && (k_t > 0)) ||
         ((h_t == 0) && (k_t == 0) && (l_t >=0))) {
          n_t = h_t + ((Nk + k_t) % Nk)*Nh + ((Nl + l_t) % Nl)*Nh*Nk ;
          *ct = 1;
     }
     else {
          n_t = -h_t + ((Nk - k_t) % Nk)*Nh + ((Nl - l_t) % Nl)*Nh*Nk ;
          *ct =-1;
     }
     return (n_t) ;
}
int	assemble(int h, int k, int l) 
{

      if (h < 0 || h >= Nh)
          return (-1) ;
      if ((k < -Nk/2) || (k >= Nk/2))
          return (-1) ;
      if ((l < -Nl/2) || (l >= Nl/2))
          return (-1) ;

      return (h + ((Nk + k) % Nk)*Nh + ((Nl + l) % Nl)*Nh*Nk) ;
}

void	setDsqLimit(char *mask, int N)
{
     real	nusqval = 0 ;
     static	char	DONE_IT = FALSE ;
     int	n ;
     int	h, k, l ;
     real	thisnusq ;


     if (!DONE_IT) {

          /* Fetch max(nusq) for N entries; N may be either Nhkl or NM*Nhkl */

          for (n = 0; n < N; n++, mask++) 
	       if (*mask) {
	            fetch_hkl(n%Nhkl, &h, &k, &l) ;
	            if ((thisnusq = nusq(h, k, l)) > nusqval) 
	                nusqval = thisnusq ; 
               }

          if (nusqval < limit_rat) {
	       dsq_limit = nusqval ;
               fprintf(fp_log, 
	  "\n1/d-squared shell limit is %g, based on data\n\n", eta*dsq_limit) ;
          }
          else {
	       dsq_limit = limit_rat ;
               fprintf(fp_log, 
	  "\n1/d-squared shell limit is %g, based on Eden parameter\n\n", 
	       eta*dsq_limit) ;
          }
          DONE_IT = TRUE ;
     }
}

		/*********************************************
	        Calculate the R factor segment by segment.
		R factor = Sum [ |Full - |R|| ] / Sum [ Full ]
		*********************************************/

real get1Rfac(char *what, real *Flocal, COMPLEX *Rlocal, int Nlocal)
{
     real	sumFR[NUMSEGS] ;
     real	sumF[NUMSEGS]  ;
     real	Rabsq ;
     real	Fabs ;
     double	rat[NUMSEGS] ;
     double	totFR = 0, totF = 0, totrat ;

     int	m, n, nn ;
     int	h, k, l ;

     for (m = 0; m < NUMSEGS; m++) {
	  sumFR[m] = 0 ;
	  sumF[m] = 0 ;
     }

     for (nn = 0; nn < Nlocal; nn++) {

          n = nn % Nhkl ;

	  if (n == 0)
               m = 0 ;

          else {

               fetch_hkl(n, &h, &k, &l) ;

	       for (m = 1; m < NUMSEGS; m++)
		    if (nusq(h, k, l) < nusqlim[m])
			 break ;
          }
	  if (m < NUMSEGS) {   

               Rabsq = (Rlocal+nn)->re * (Rlocal+nn)->re + 
	     	       (Rlocal+nn)->im * (Rlocal+nn)->im ;
	       Fabs  = fabs(*(Flocal+nn)) ;

	       sumFR[m] += fabs(Fabs - sqrt(Rabsq)) ;
               sumF[m]  += Fabs ;
	  }
     }

     for (m = 0; m < NUMSEGS; m++) {
	  totFR += sumFR[m];
	  totF  += sumF[m];
	  if (sumF[m] > 0)
               rat[m] = sumFR[m] / sumF[m] ;
          else
               rat[m] = 0. ;
     }

     totrat = totFR / totF ;

     sprintf(message, 
            "%s %10lg %10lg %10lg %10lg %10lg %10lg", 
	    what, rat[0], rat[1], rat[2], rat[3], rat[4], totrat) ;
     printTwice(message) ;

     return (totrat) ;
}

		/*********************************************
	        Calculate the R factor segment by segment,
		weighted appropriately.
		R factor = Sum [ |Full - |R|| ] / Sum [ Full ]
		*********************************************/

real get1Rfac_weighted(char *what, real *Flocal, COMPLEX *Rlocal, 
			 int Nlocal, real *mirwt)
{
     real	sumFR[NUMSEGS] ;
     real	sumF[NUMSEGS]  ;
     real	Rabsq ;
     real	Fabs ;
     real	rat[NUMSEGS] ;
     real	sqrt_mirwt ;
     real	totFR = 0, totF = 0, totrat ;

     int	m, n, nn ;
     int	h, k, l ;

     for (m = 0; m < NUMSEGS; m++) {
	  sumFR[m] = 0 ;
	  sumF[m] = 0 ;
     }

     for (nn = 0; nn < Nlocal; nn++) {

          n = nn % Nhkl ;
	  sqrt_mirwt = sqrt(mirwt[nn / Nhkl]) ;

	  if (n == 0)
               m = 0 ;

          else {

               fetch_hkl(n, &h, &k, &l) ;

	       for (m = 1; m < NUMSEGS; m++)
		    if (nusq(h, k, l) < nusqlim[m])
			 break ;
          }
	  if (m < NUMSEGS) {   

               Rabsq = (Rlocal+nn)->re * (Rlocal+nn)->re + 
		       (Rlocal+nn)->im * (Rlocal+nn)->im ;
	       Fabs  = fabs(*(Flocal+nn)) ;

	       sumFR[m] += sqrt_mirwt * fabs(Fabs - sqrt(Rabsq)) ;
               sumF[m]  += sqrt_mirwt * Fabs ;
          }
     }

     for (m = 0; m < NUMSEGS; m++) {
	  totFR += sumFR[m];
	  totF  += sumF[m];
	  if (sumF[m] > 0)
               rat[m] = sumFR[m] / sumF[m] ;
          else
               rat[m] = 0. ;
     }

     totrat = totFR / totF ;

     sprintf(message, 
            "%s %10g %10g %10g %10g %10g %10g", 
	    what, rat[0], rat[1], rat[2], rat[3], rat[4], totrat) ;
     printTwice(message) ;

     return (totrat) ;
}

void checkF000(real *fo, real *sig, char *mask)
{
     /* Check that f(000) is set */

     if (*mask == 0) {
	  if (interactive) {

                 prompt("Please enter F(0,0,0) and corresponding sigma: ") ;
		 t1 = 0 ;
		 t2 = 0 ;

                 sscanf(terminp, "%g %g", &t1, &t2) ;
	         fprintf(fp_log, "\nF(0,0,0) = %g, sig(0,0,0) = %g\n",
		 t1, t2) ;
		 *fo = t1 ;
		 *sig = t2 ;
		 *mask = 1 ;
	  }
	  else {
	       EdenWarning( 
	       "Missing F(0,0,0) - total # of electrons in unit cell.") ;
	       EdenWarning( "Please put it in the .fobs file.") ;
	       EdenWarning( 
	       "If you are using sigmas, put in a sigma(0,0,0) as well -") ;
	       EdenWarning( "e.g. sqrt(F(0,0,0))/10.") ;
	       EdenError("No F(0,0,0)") ;
	  }
     }
}
void	scalefobs(real *array, real sfac)
{
     int	n ;

     array++ ;		/* DON'T scale F(000)!	*/
     if (sfac != 1.) 
	  for (n = 1; n < Nhkl; n++, array++) 
	       *array *= sfac ;

}

void hklTabHead()
{
     real	Amin ;
     int	m ;

     Amin = sqrt(1 / dsq_limit) ;

	 /* Initialize segment limits in Angstrom with "nice" numbers. */

     for (m = 1; m < NUMSEGS; m++) 
          sprintf(Arange[m], "%4.1f", pow(2.0, 0.5*(NUMSEGS-m-1))*Amin) ;

   /* Recalculate limits in terms of 1/d^2, but don't fiddle with outermost! */

     for (m = 1; m < NUMSEGS; m++) {	 

          sscanf(Arange[m], "%g", &t1) ;
          nusqlim[m] = 1. / (t1*t1) ;

	  if (m == NUMSEGS-1)		
               nusqlim[m] = dsq_limit ;	

     }
}
void shellHeader()	/* called by dphase; no (000) */
{
     printTwice(
     "\nhkl shell   <1/8    1/8-1/4    1/4-1/2       rest        all") ;

     sprintf(message,
     "resol (A) > %s  %s-%s  %s-%s  %s-%s        all \n",
     Arange[1], Arange[1], Arange[2], Arange[2], Arange[3], Arange[3],
     Arange[4]) ;
     printTwice(message) ;
}
void shellHeader_R()	/* called by Solve & Back; 2 copies (printTwice() */
{
     printTwice("\n\t\tR factors for fractions of 1/d-squared\n") ;
     printTwice(
    "\t       (000)        1/8        1/4        1/2  Remainder    Overall\n") ;

     sprintf(message, 
     "resol (A)                > %s  %s-%s  %s-%s  %s-%s        all \n\n",
     Arange[1], Arange[1], Arange[2], Arange[2], Arange[3], Arange[3],
     Arange[4]) ;
     printTwice(message) ;
}

void shellHeader_F()	/* called by init (for Solve only); goes to log. */
{

     fprintf(fp_log,
     "F Percentages, Fo and Fc counts for fractions of 1/d-squared\n\n") ;
     fprintf(fp_log,
     "        (000)        1/8        1/4        1/2  Remainder    Overall\n\n") ;
     fprintf(fp_log, 
     "res (A)           > %s  %s-%s  %s-%s  %s-%s        all \n\n",
     Arange[1], Arange[1], Arange[2], Arange[2], Arange[3], Arange[3],
     Arange[4]) ;
}
	/*********************************************
        Create a mask in which we enforce the
	half-ellipsoid rule (for expanded output).
	*********************************************/

char	*get_P1mask() 
{
     int	m, n, onen, h, k, l ;
     char	*mask ;

     mask = (char *) e_malloc(NM*Nhkl*sizeof(char), "get_P1mask") ; 
          
     for (m = 0; m < NM; m++) {
       for (n = 0; n < Nhkl; n++) {

          onen = m*Nhkl + n ;
          *(mask+onen) = 0 ;

          fetch_hkl(n, &h, &k, &l) ;

          if (nusq(h, k, l) < limit_rat) {
	  
	       if (((h > 0)) ||
		   ((h == 0) && (k > 0)) ||
		   ((h == 0) && (k == 0) && (l >=0))) 
	            *(mask+onen) = 1 ;
          }
        }
     }

     return (mask) ;
}

		/**************************************************
		Determine unique reflections, for which **pmask = 1.  

		Procedure: Initialize all points in mask to zero.
		Set unique points to 1, related non-unique pts to 2.
		When all mask has been accounted for, reset to 0 
		all points outside the resolution limits and
		all points that were previously set to 2.
		**************************************************/

void prepare_unique_hklmask(char **pmask)
{
     int	n ;
     int	h, k, l ;
     int	kk, ll ;
     int	sum_as1 = 0, sum_notas1 = 0 ;
     char	*umask ;

     unique_hkl = (char *) e_malloc(Nhkl*sizeof(char), "hklutil") ;
     
     *pmask = unique_hkl ;
     
     for (n = 0, umask = unique_hkl; n < Nhkl; n++, umask++)
	  *umask = 0 ;

     umask = unique_hkl ;

     for (h = 0; h < Nh; h++) {

          for (kk = 0; kk < Nk; kk++) {

	       k =  (kk > Nk/2) ? kk - Nk : kk ;

	       if (k != Nk/2) {

                    for (ll = 0; ll < Nl; ll++) {

	                 l =  (ll > Nl/2) ? ll - Nl : ll ;

		         if (l != Nl/2) {

                              if ((n = assemble(h, k, l)) >= 0) {

			           if (*(umask+n) == 1) {
                              sprintf(message, 
	   		      "Illegally set point at hkl = %d %d %d",
     			      h, k, l) ;
                              EdenError(message) ;
			           }
			           else if (*(umask+n) == 0) {
			                *(umask+n) = 1 ;
			                check_others(unique_hkl, h, k, l, n) ;
                                   }
			           /*  else if (*(umask+n) == 2) do nothing */
                              }
                         }
                    }		
               }		
          }	
     }	
     
     /* Reset the mask outside the nusq limit */

     for (n = 0, umask = unique_hkl; n < Nhkl; n++, umask++) {

          fetch_hkl(n, &h, &k, &l) ;
						    
          if (nusq(h, k, l) >= limit_rat)  
	       *umask = 0 ;

     }			

     /* Statistics */

     for (n = 0, umask = unique_hkl; n < Nhkl; n++, umask++) {
	  if (*umask == 1)
	       sum_as1++ ;
	  else if (*umask == 2)
	       sum_notas1++ ;
	  /* else (*umask == 0) not interested */
     }

     sprintf(message, "\nMax. # of unique and non-unique reflections: %d %d\n", 
	     sum_as1, sum_notas1) ;
     printTwice(message) ;


     /* Now replace "others" by "masked-out" */

     for (n = 0, umask = unique_hkl; n < Nhkl; n++, umask++) 
	  if (*umask == 2)
	       *umask = 0 ;
}

void	check_others(char *umask, int h, int k, int l, int n)
{
     int	p ;
     real	off ;
     real	*mat ;
     int	newh, newk, newl, newn ;
     int	nf ;

     for (p = 0, mat = matop; p < Ncs; p++, mat += MEL) {

          apply_symop_hkl(h, k, l, mat, &newh, &newk, &newl, &off) ;

          if (((newn = assemble( newh,  newk,  newl)) >= 0) &&
               (newn != n)) {

		if (*(umask+newn) == 1) {
                     sprintf(message, 
     "Illegally set unique point - hkl = %d %d %d, new h,k,l = %d %d %d", 
                     h, k, l, newh, newk, newl) ;
                     EdenError(message) ;
                }
		else if (*(umask+newn) == 0) 
		     *(umask+newn) = 2 ;
		/* else if (*(umask+newn) == 2) do nothing */
          }
          		
				/* no check for Friedel pairs if anomalous data! */
          if (!anom_flag) {

/***               if (((nf = assemble(-h, -k, -l)) > 0) && (nf != n)) { ***/
               if (((nf = assemble(-newh, -newk, -newl)) > 0) && (nf != n)) {

                   if (*(umask+nf) == 1) {
                        sprintf(message, "Illegally set Friedel point at hkl = %d %d %d",
		           -newh, -newk, -newl) ;
                         EdenError(message) ;
                   }
                   else if (*(umask+nf) == 0) {
                        *(umask+nf) = 2 ;
                   }
                   /*  else if (*(umask+nf) == 2) do nothing */
               }
          }
     }
}

void set_centrics(char *mcentric)
{
     real	off ;
     real	*mat ;
     int        p, n ;
     int        h, k, l ;
     int	newh, newk, newl ;

     for (n = 0; n < Nhkl; n++) 
          *(mcentric+n) = 0 ;

	/***********************************************************
	We look for (hkl)'s which, under a symmetry operation, 
	give back the negatives of the same reciprocal indices.
	See Giacovazzo, Fundamentals of Crystallography, pp 156-157. 
	***********************************************************/

     for (n = 1; n < Nhkl; n++) {

          fetch_hkl(n, &h, &k, &l) ;

          if (nusq(h, k, l) < limit_rat) {

               for (p = 0, mat = matop; p < Ncs; p++, mat += MEL) {

	           apply_symop_hkl(h, k, l, mat, &newh, &newk, &newl, &off) ;

	           if ((newh == -h) && (newk == -k) && (newl == -l)) 
		       *(mcentric+n) += 1 ;
               }	
          }	
     }			
}

void set_cen_phases(char *mcentric, real *mphases)
{
     real	off ;
     real	*mat ;
     int        p, n ;
     int        h, k, l ;
     int	newh, newk, newl ;

     for (n = 0; n < Nhkl; n++) {
          *(mcentric+n) = 0 ;
          *(mphases+n) = 0 ;
     }

	/***********************************************************
	We look for (hkl)'s which, under a symmetry operation, 
	give back the negatives of the same reciprocal indices.
	We set the corresponding phases according to Giacovazzo, 
	3.39 (but without the n*pi factor that is added in ranphase)
	***********************************************************/

     for (n = 1; n < Nhkl; n++) {

          fetch_hkl(n, &h, &k, &l) ;

          if (nusq(h, k, l) < limit_rat) {

               for (p = 0, mat = matop; p < Ncs; p++, mat += MEL) {

	           apply_symop_hkl(h, k, l, mat, &newh, &newk, &newl, &off) ;

	           if ((newh == -h) && (newk == -k) && (newl == -l)) {
		       *(mcentric+n) = 1 ;
		       *(mphases+n) = off / 2. ;
		       break ;
                   }		
               }	
          }	
     }			

     return ;
}
void	Cmaskit(COMPLEX *array, char *mask, int N)
{
      int	n ;

      for (n = 0; n < N; n++, array++, mask++) {
	   array->re *= *mask ;
	   array->im *= *mask ;
      }	   
}

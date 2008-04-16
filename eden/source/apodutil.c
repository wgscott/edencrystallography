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

                                APODUTIL

  Title:        Utilities for running apodfc and apodfo. 
  Author:       Hanna Szoke
  Date:         10/5/00: full comparison of 2 wil files for scaling 
		fobs to fcalc; see scale_fo().

*******************************************************************************/

#include "util.h"
#include "cellparams.h"
#include "dims.h"	
#include "apodize.h"

#define	BINWID 0.002 	/* default width of a bin */

	/* variables declared in apodize.h	*/

int	useW0 ;		/* use W0 correction for proteins? T/F */
int	Nbin ;		/* # bins used in program */
real	delta ;		/* smearing factor */
real	binwidth ;	/* width of a bin */
char	*maskap ;	/* declared here although used only in apodf[co].c */
real	*w_lnFsq ;	/* wilson y coords with or w/o w0 correction */
real	*apodwt ;	/* weights associated with Iave */

	/* variables local to apodutil.c */

static	real	*Ibin ;		/* Binned, not averaged, I's */
static	real	*sigsqbin ;	/* Binned, not averaged, sigma's */
static	real	*Iave ;		/* Binned & averaged I's */
static	real	*sigsqave ;	/* Binned & averaged sigma's */
static	int	*cases ;	/* # of entries in a bin */

static	int	lin_orig ;	/* index of start of linearized signal */
static	int	m0 ;		/* index of first non-empty bin */
static	real	min_x ;		/* min x corresponding to max data res */
static	real	min_x0 ;	/* min x corr. to max data res (corrected) */
static	real	max_x ;		/* max x corresponding to min data res */
static	real	y_intercept ;
static	real	slope ;

static	real	min_res = 3.5 ;		/* default for linearization */
static	real	max_res = 0.05 ;	/* default for linearization */
static	real	min_res0 = 7.0;		/* default for corrected wilson plots */
static	real	apod_res = 0 ;

static	void	send_it_out_asc(char *, real *, int, float, float) ;
static	void	generateBins(int, real *, real *, real *) ;
static	char	checkBins() ;
static	void	averageBins() ;
static	void	newvar() ;
static	int	read_it_in_asc(char *, real **, real **) ; 
static	int	genlinBins(real *, char) ;		



void a_readInput()    /* Read input, identifying data in terms of keywords. 
                         Do requisite reciprocal lattice calculations.  */
{
     int	k ;
     char	good = FALSE ;

     readBasicInput() ;
     binwidth = BINWID ;
     useW0 = TRUE ;

     for (k = 0; k < Nlines; k++) {
          strcpy(nextline, allinp + k*MAXSTRING) ;
          sscanf(nextline, "%s", id) ;

          if (strcmp(id, "APOD_RES") == 0) {
               sscanf(nextline, "%*s %g", &t1) ;
               apod_res = t1 ;
          }
          else if (strcmp(id, "MIN_RES") == 0) {
               sscanf(nextline, "%*s %g", &t1) ;
               min_res = t1 ;
          }
          else if (strcmp(id, "MAX_RES") == 0) {
               sscanf(nextline, "%*s %g", &t1) ;
               max_res = t1 ;
          }
          else if (strcmp(id, "BINWIDTH") == 0) {
               sscanf(nextline, "%*s %g", &t1) ;
               binwidth = t1 ;
          }
     }

     if (apod_res == 0)
	  apod_res = input_res;

     if (min_res < max_res) {
         EdenWarning( "Your resolution limits seem to be inverted!") ;
	 EdenError(" MAX_RES should be a smaller number than MIN_RES") ;
     }

     for (k = 1; k < 100; k++) 
        if (fabs(binwidth - BINWID*k) < 1.e-5) {
	     good = TRUE ;
	     break ;
        }

     if (!good) {
	 sprintf(message,"Binwidth (%g) is not a multiple of %g", binwidth,
		 BINWID) ;
	 printTwice(message) ;
	 EdenWarning("Turning off useW0 flag.\n") ;
	 useW0 = FALSE ;
     }
     
     min_x = 1.0 / (min_res*min_res) ;
     max_x = 1.0 / (max_res*max_res) ;
     min_x0 = 1.0 / (min_res0*min_res0) ;

     sprintf(message, "Linearization limits are (%g, %g) A.", 
	     min_res, max_res) ;
     printTwice(message) ;
     sprintf(message, "Apodization resolution = %g A.", apod_res) ;
     printTwice(message) ;

     return ;
}

void	init_apod_arrays(int Ntot, real *Int, real *ksq, real *sigsqI)
{
     int	n ;

     for (n = 0; n < Ntot; n++) {
          *(Int+n) = 0 ;
          *(ksq+n) = 0 ;
          *(sigsqI+n) = 0 ;
     }
     return ;
}
void	allocate_Bin_arrays(int N) 
{
     cases   = (int *) e_malloc(N*sizeof(int), caller) ;
     Ibin = (real *) e_malloc(N*sizeof(real), caller) ;
     sigsqbin = (real *) e_malloc(N*sizeof(real), caller) ;
     Iave = (real *) e_malloc(N*sizeof(real), caller) ;
     sigsqave = (real *) e_malloc(N*sizeof(real), caller) ;
     w_lnFsq = (real *) e_malloc(N*sizeof(real), caller) ;
     apodwt = (real *) e_malloc(N*sizeof(real), caller) ;
}
void	do_binning(int Ntot, real *Int, real *ksq, real *sigsqI)
{
     char okay ;

     do {
          generateBins(Ntot, Int, ksq, sigsqI) ;
          okay = checkBins() ;
     } while (!okay) ;

     averageBins() ;
     newvar() ;

     return ;

}

			/*********************************************
              		Move Int, ksq, sigsqI from lists to bins.
			*********************************************/

void	generateBins(int Ntot, real *Int, real *ksq, real *sigsqI)
{
     int	m, n ;
     int	overflow = 0 ;
     real	Iover = 0, sover = 0 ;

     /*  Prepare bins */

     for (m = 0; m < Nbin; m++) {
	  *(Ibin+m) = 0 ;
	  *(sigsqbin+m) = 0 ;
	  *(cases+m) = 0 ;
     }

     for (n = 0; n < Ntot; n++) {

	  m = (*(ksq+n) + 0.5*binwidth) / binwidth ;

          if (*(Int+n) > 0) {

               if (m >= Nbin) {
	            overflow++ ;
	            m = Nbin-1 ;
                    Iover += *(Int+n) ;
	            sover += *(sigsqI+n) ; 
	       }

	       *(Ibin+m) += *(Int+n) ;
	       *(sigsqbin+m) += *(sigsqI+n) ;
	       *(cases+m) += 1 ;
	  }

     }

    if (overflow > 0) {
	  sprintf(message, 
	  "%d entries overflowed the highest resolution bin.", overflow) ;
	  printTwice(message) ;
	  sprintf(message, 
  "Their contribution to that bin was I = %g (sig = %g)", Iover, sover) ;
	  printTwice(message) ;
	  sprintf(message, 
  "                            out of I = %g (sig = %g)", *(Ibin+Nbin-1), 
							*(sigsqbin+Nbin-1)) ;
	  printTwice(message) ;
     }

     return ;
}

			/*********************************************
              		Check the binning procedure 
			*********************************************/

char	checkBins()
{
     int	m1, mlast=0 ;
     char	okay = TRUE ;

     /****************************************************
     Check bins: empty ones at the beginning are okay;
     empty ones at the end are also okay, but they cause 
     Nbin to be suitably redefined.   
     NOTE: an empty bin can occur either because no (hkl)
     entries go into it or because the total intensity,
     Iave (not yet averaged) is 0.
     8/23/99 - Corrected code and messages re. empty bins 
     in middle of range.
     *****************************************************/

     if (verbose)  {
	   fprintf(fp_log, "Initial bin information:\n\n") ;

           for (m0 = 0; m0 < Nbin; m0++) 
               fprintf(fp_log," m = %d, Ibin = %g, cases = %d\n",
	       m0, *(Ibin+m0), *(cases+m0)) ;

               fprintf(fp_log, "\n") ;
          }

     for (m0 = 0; m0 < Nbin; m0++) 
          if ((*(cases+m0) > 0) && (*(Ibin+m0) > 0))
              break ;

     for (mlast = Nbin-1; mlast >= m0; mlast--) 
          if ((*(cases+mlast) > 0) && (*(Ibin+mlast) > 0))
              break ;
		
     for (m1 = m0; m1 < mlast; m1++) 
          if ((*(cases+m1) == 0) || (*(Ibin+m1) == 0)){
               sprintf(message,
	       "There is at least one empty bin (# %d) within range (%d, %d)!",
	       m1, m0, mlast) ;
	       printTwice(message) ;
               okay = FALSE ;
               break;
          }	       

     if (!okay) {
          binwidth += BINWID;
          Nbin = (dsq_limit + 0.5*binwidth)/ binwidth + 1 ;
     }
     else
          Nbin = mlast + 1 ;

     fprintf(fp_log, "First and last occupied bins are %d and %d\n", 
	  m0, mlast) ;

     return(okay) ;
}

			/*********************************************
              		Average bin contents based on # of cases.
			Note sigsqave is averaged with 1/sq(# of cases)
			*********************************************/

void	averageBins()
{
     real	d ;
     int	m ;
     char	bad = FALSE ;

     for (m = m0; m < Nbin; m++) {
	  *(Iave+m) = 0. ;
	  *(sigsqave+m) = 0. ;
     }
     
     for (m = m0; m < Nbin; m++) {
	  d = (real) *(cases+m) ;
	  *(Iave+m) = *(Ibin+m) / d ;
	  *(sigsqave+m) = *(sigsqbin+m) / (d*d) ;
     }
     
     if (verbose) {
	  fprintf(fp_log, "\n") ;
	  fprintf(fp_log, " bin      d(A)       1/d*d    count    ave(Fsq)\n") ;
	  fprintf(fp_log, " ---     -----     -------    -----    --------\n") ;

          for (m = m0; m < Nbin; m++)  {
	       d = 1. / sqrt((m+0.5)*binwidth) ;
	       fprintf(fp_log, "%4d  %8.2f  %10.5f  %5d  %10.2f\n", 
			    m,  d, (m+0.5)*binwidth, *(cases+m), *(Iave+m)) ;
               }
     }

     for (m = m0; m < Nbin; m++) {
	  if (*(sigsqave+m)  <= 0) {
	       sprintf(message, "sigsq array is not > 0 at m = %d", m) ;
	       EdenWarning(message) ;
	       bad = TRUE ;
          }
     }
     if (bad)
	  EdenError("Quitting ...") ;
     
     return ;
}
void	newvar() 
{
     int	m ;

     for (m = m0; m < Nbin; m++) {
          *(w_lnFsq+m) = log(*(Iave+m)) ;
	  *(apodwt+m) = *(Iave+m) * *(Iave+m) / *(sigsqave+m) ;
     if (verbose)
	  fprintf(fp_log, 
	  "m = %d, ln Fsq = %g, weight = %g\n", m, *(w_lnFsq+m), *(apodwt+m)) ;
     }
}

int	genlinBins(real *lin_lnFsq, char corr)
{
     int	m ;
     real	x ;
     real	xm ;
     int	m1 ;

     for (m = 0; m < Nbin; m++) {
	  x = (m+0.5)*binwidth ;
	  *(lin_lnFsq+m) = y_intercept + x*slope ;
     } 

     for (m1 = 0; m1 < Nbin; m1++) {
	  x = m1*binwidth ;
	  if (x > max_x) 
	       break ;
     } 

     /* set lin_orig */

     xm = (corr) ? min_x0: min_x ;

     for (lin_orig = 0; lin_orig < Nbin; lin_orig++) 
	  if ((x = lin_orig*binwidth) > xm)
	       break ;
     
     return(m1-lin_orig) ;
}

     /************************************************************************
     Derive best fit of Fsq to y = y0 + slope*x, in range (min_x, max_x).
     See Numerical Recipes, 1986, p.505 (14.2.6).
     *************************************************************************/

real     linfit(int N, real *lnFsq, real *awt, char corr) 
{
     real	Sx = 0, Sy = 0, Sw = 0 ;
     real	Sxx = 0, Sxy = 0 ;
     real	Den = 0 ;
     real	x, y, w ;
     real	xm ;
     real	ksqmin = 100., ksqmax = 0 ;
     real	Bcryst ;
     int	m ;

     xm = (corr) ? min_x0: min_x ;

     for (m = 0; m < N; m++) {

	  x = (m+0.5)*binwidth ;
	  y = *(lnFsq+m) ;
	  w = *(awt+m) ;

	  if ((x > xm) && (x < max_x)) {
	       Sw += w ;
	       Sx += x*w ;
	       Sy += y*w ;
	       Sxx += x*x*w ;
	       Sxy += x*y*w ;
          }
	  if (x < ksqmin)
	       ksqmin = x ;
	  if (x > ksqmax)
	       ksqmax = x ;
     }
     Den = Sw*Sxx - Sx*Sx ;

     if (Den == 0) {
          sprintf(message, "ksq range is (%g, %g), ", ksqmin, ksqmax) ;
          EdenWarning(message) ;
          sprintf(message, 
            "min_x = %g (min_res = %g), max_x = %g (max_res = %g).",
            xm, min_res, max_x, max_res) ;
          EdenWarning(message) ;   
          EdenWarning(
 "If you are apodizing strongly, you should probably reset min_res to around 20.\n") ;
          EdenError("Please set min_res and/or max_res!") ;
     }

     y_intercept = (Sxx * Sy - Sx * Sxy) / Den ;
     slope       = (Sw * Sxy - Sx * Sy) / Den ;

     fprintf(fp_log, 
       "\nCurve fitted over (%g, %g) to (y = y0 + slope*x)\n",  xm, max_x) ;
     fprintf(fp_log, "with y0 = %g, (variance: %g)\n", y_intercept, Sxx / Den) ;

     fprintf(fp_log, "and slope = %g (variance: %g).  Cov(y0, b) = %g.\n", 
		slope, Sw / Den, -Sx / Den) ;

     /************************************************************************
     We use the slope of ln(Fsq) vs. ksq ( = nusq()) to calculate Bcryst.
     *************************************************************************/

     Bcryst = -2. * slope ;
     sprintf(message,
	"\nThe average crystallographic B factor is  %g Asq,", Bcryst) ;
     printTwice(message) ;

     sprintf(message,
	"corresponding to a resolution of %g A.\n", 
	     1./(TWOPI*grid_spacing/input_res) * sqrt(Bcryst/eta)) ;
     printTwice(message) ;


     return (Bcryst) ;
}

real	get_delta(real Bcryst)
{
     real	ratsq ;
     real	Btarget, Beden ;

     /************************************************************************
     There should be NO delta , NO delfac.  See AS's notes, 8/25/02.
     Here is what we should do: 

	Bcryst = -2*slope, where slope is linear fit of ln(Fsq) vs nusq;

     It is the actual average value of |F| vs exp[-Bcryst*nusq/4].
     We want to change B into

	Btarget = (apod_res/input_res)^2 * Beden, where Beden = 4PISQ*dr^2*eta,

     but ONLY if Btarget > Bcryst.  In such a case, each F() has to be corrected:

	F'(h,k,l) = F(h,k,l)*exp[-1/4(Btarget - Bcryst)*nusq].

     (we have corrected the output message.)
     *************************************************************************/

     ratsq = (apod_res / input_res) * (apod_res / input_res) ;

     delta = ratsq - Bcryst / (4 * PISQ * eta * grid_spacing * grid_spacing) ;

     Beden = 4 * PISQ * eta * grid_spacing * grid_spacing ;
     Btarget = ratsq * Beden ;

     if (delta > 0) 
          sprintf(message, 
       "The data should be smeared using a target B factor of %g.\n", Btarget) ;
     else { 
          sprintf(message,  
       "The data should be unsmeared: target B factor (%g) is not greater than Bcryst (%g).\n", 
       Btarget, Bcryst) ; 
	  delta = 0 ;
     }
     printTwice(message) ;

     return (delta) ;
}

void apply_w0_correction() 
{
     FILE 	*fp ;
     real	off ;
     real	epsi ;
     int	n ;

     strcpy(message, getenv("EDENHOME")) ;
     strcat(message, "/source/w0_d");

     fp = fopen(message, "r");
     if (fp == NULL) {
	  strcat(message, " - cannot open this file!") ;
          EdenError(message) ;
     }

     epsi = 0.01 * binwidth ;
     off = m0 * binwidth ;
     n = m0 ;
     while (fgets(nextline, MAXSTRING, fp) != NULL)  {
        sscanf(nextline, "%g %g", &t1, &t2) ;  
        if (t1 >= off - epsi) {
	     *(w_lnFsq + n) += t2 ;
             off += binwidth ;
             n++ ;
        }
        if (n >= Nbin - m0) break ;
     }

     return ;
}

void make_wilson(char corr, real *dstdev)  
{
     real	*l_lnFsq ;
     real	*d_lnFsq ;

     int	Nlin ;			/* length of lin_lnFsq */
     float	origin, interval ;	/* signal origin & interval */

     real	x, dave = 0 ;
     int	count = 0 ;

     int m ;
     real	xm ;
     char	ext[20] ;

  /*******************************************************
  Now prepare for display the original binned log(I^2),
  the linearly fitted version and the difference file.
  *******************************************************/
     
     l_lnFsq = (real *) e_malloc(Nbin*sizeof(real), "apodutil") ;
     d_lnFsq = (real *) e_malloc(Nbin*sizeof(real), "apodutil") ;

     Nlin = genlinBins(l_lnFsq, corr) ;

     if (corr) 
         strcpy(ext, "_w0corr") ;
     else 
         strcpy(ext, "") ;

     interval = binwidth ;
     origin = m0*interval ;

     sprintf(out_filename, "wil%s", ext) ;

     send_it_out_asc(out_filename, w_lnFsq+m0, Nbin-m0, origin, interval) ;
     origin = lin_orig*interval ;

     sprintf(out_filename, "lin_wil%s", ext) ;

     send_it_out_asc(out_filename, l_lnFsq+lin_orig, Nlin, origin, interval) ;
     if (corr) 
	  xm = min_x0 ;
     else 
	  xm = min_x ;

     for (m = 0; m < Nbin; m++)
          *(d_lnFsq+m) = *(l_lnFsq+m) - *(w_lnFsq+m) ;

     for (m = 0; m < Nbin; m++) {

	  x = (m+0.5)*binwidth ;

	  if ((x > xm) && (x < max_x)) {
               dave += *(d_lnFsq+m) ;
	       count++ ;
          }

     }
     dave = dave / count ;
     *dstdev = 0 ;

     for (m = 0; m < Nbin; m++) {

	  x = (m+0.5)*binwidth ;

	  if ((x > xm) && (x < max_x)) 
               *dstdev += (*(d_lnFsq+m)  - dave) * (*(d_lnFsq+m)  - dave) ;

     }
     *dstdev = sqrt(*dstdev / count) ;

     if (corr) 
         strcpy(ext, "corrected") ;
     else 
         strcpy(ext, "uncorrected") ;

     sprintf(message, 
     "Standard deviation of wilson curve (%s) with respect", ext) ;
     printTwice(message) ;
     sprintf(message, "   to its linearized version is %g.", *dstdev) ;
     printTwice(message) ;

     e_free(Nbin*sizeof(real), l_lnFsq) ;
     e_free(Nbin*sizeof(real), d_lnFsq) ;

     return ;
}

void make_selected_wilson(char *name)  
{
     float	origin, interval ;	/* signal origin & interval */

     interval = binwidth ;
     origin = m0*interval ;

     send_it_out_asc(name, w_lnFsq+m0, Nbin-m0, origin, interval) ;

     return ;
}

real	choose_delta(real d1, real crit1, real B1, 
                     real d2, real crit2, real B2)
{
     int	reread_wil() ;	
     real	choice, mychoice, otherchoice ;
     real	mincrit ;
     real	oldBcryst, newBcryst ;
     real	myBcryst, otherBcryst ;
     char	ext[20] ;

     mincrit = (crit2 < crit1) ? crit2 : crit1 ;
   
     if  ((crit1 - crit2) > 0.1*mincrit) {
          mychoice = d2 ;
          otherchoice = d1 ;
	  myBcryst = B2 ;
	  otherBcryst = B1 ;
	  strcpy(ext, "corrected ") ;
     }
     else {
          mychoice = d1 ;
          otherchoice = d2 ;
	  myBcryst = B1 ;
	  otherBcryst = B2 ;
	  strcpy(ext, "uncorrected ") ;
     }
     if (graphics) {
          fprintf(fp_log, "\tProposed plot is %s.", ext) ;
          sprintf(message, "\n\tProposed plot is %s - ok [y/n]? ", ext) ;
          prompt(message) ;

          if (terminp[0] == 'n') {
              choice = otherchoice ;
	      oldBcryst = otherBcryst ;
          }
          else {
              choice = mychoice ;
	      oldBcryst = myBcryst ;
          }
     }
     else {
	  choice = mychoice ;
	  oldBcryst = myBcryst ;
     }

     if (choice == d1) {

          reread_wil() ;

	  if (d1 > 0)
               printTwice("\nUsing original Wilson plot (apodized).") ;
	  else
               printTwice("\nUsing original Wilson plot (unapodized).") ;
     }
     else {
	  if (d2 > 0)
               printTwice("\nUsing corrected Wilson plot (apodized).") ;
	  else
               printTwice("\nUsing corrected Wilson plot (unapodized).") ;
     }

     /************************************************************************
     We want to report the effective average B factor selected.
     See comments in function get_delta().
     *************************************************************************/

     newBcryst = oldBcryst + 
		  (4 * PISQ * eta * choice * grid_spacing * grid_spacing) ;

     sprintf(message,
	"\nThe new average crystallographic B factor is  %g Asq.", newBcryst) ;
     printTwice(message) ;
     sprintf(message,
	"The crystallographic B factor correction is  %g Asq.\n", 
	newBcryst - oldBcryst) ;
     printTwice(message) ;

     return(choice) ;
}
void show_wilson_plots() 
{
     if (graphics) {
          printTwice("\n") ;

          if (useW0)
	       strcpy(message, 
	       "xmgrace wil lin_wil wil_w0corr lin_wil_w0corr -legend load") ;
          else
	       strcpy(message, 
	       "xmgrace wil lin_wil -legend load") ;
          system(message) ;
     }
     else {
	  printTwice(
"You can use a graphics program to plot wil, lin_wil (no w0 correction)") ;
	  printTwice(
"                                   and wil_w0corr, lin_wil_w0corr") ;
     }
}

		/*******************************************************
 		Write an ASCII file of x/y info.
		*******************************************************/

void send_it_out_asc(char *out_name, real *array, int N, float x0, float delx)
{
     FILE 	*fp ;
     int	i ;

     if ((fp = fopen(out_name,"w")) == NULL) 
	  EdenError("Couldn't open ascii file for writing.") ;

     fprintf(fp, "@ title \"Wilson Plot\"\n") ;
     fprintf(fp, "@ xaxis label \"1/resolution_sq (1/A^2)\"\n") ;
     fprintf(fp, "@ yaxis label \"ln(F_sq)\"\n") ;

     for (i = 0; i < N; i++)		/* write pairs of x and y */
          fprintf(fp,"%g\t%g\n", x0+i*delx, *(array+i)) ;
	  
     fclose(fp) ;
}
		/*******************************************************
 		Read an ASCII file of x/y info. Return actual length, 
		Ntotal, which may be <= N, but not > N!
		*******************************************************/

int	read_it_in_asc(char *in_name, real **xarray, real **yarray) 
{
     FILE 	*fp ;
     int	i ;
     int	Ntotal = 0 ;
     real	*xap, *yap;


     if ((fp = fopen(in_name,"r")) == NULL) 
	  EdenError("Couldn't open ascii file for reading.") ;

     while (fgets(nextline, MAXSTRING, fp) != NULL)  
	  Ntotal++ ;
	  
     fclose(fp) ;

     *xarray = (real *) e_malloc(Ntotal*sizeof(real), caller) ;
     *yarray = (real *) e_malloc(Ntotal*sizeof(real), caller) ;
     xap = *xarray ;
     yap = *yarray ;

     for (i = 0; i < Ntotal; i++) {
	  *(xap + i) = 0 ;
	  *(yap + i) = 0 ;
     }

     if ((fp = fopen(in_name,"r")) == NULL) 
	  EdenError("Couldn't open ascii file for reading.") ;

     i = 0 ;

     while (fgets(nextline, MAXSTRING, fp) != NULL)  {
          sscanf(nextline, "%s", id) ;
          if (strcmp(id, "@") != 0) {
               sscanf(nextline, "%g %g", &t1, &t2) ;
	       *(xap+i) = t1 ;
	       *(yap+i) = t2 ;
	       i++ ;
          }
     }
	  
     fclose(fp) ;

     return(i) ;
}
		/*******************************************************
 		Reread the "wil" file and replace w_lnFsq values. 
		*******************************************************/

int	reread_wil() 
{
     FILE 	*fp ;
     int	m ;


     if ((fp = fopen("wil","r")) == NULL) 
	  EdenError("Couldn't open 'wil' file for reading.") ;

     m = m0 ;

     while (fgets(nextline, MAXSTRING, fp) != NULL)  {
          sscanf(nextline, "%s", id) ;
          if (strcmp(id, "@") != 0) {
               sscanf(nextline, "%g %g", &t1, &t2) ;
	       *(w_lnFsq+m) = t2 ;
	       m++ ;
          }
     }
	  
     fclose(fp) ;

     return(m) ;
}
void	scale_fo(char *name_fo, char *name_fc) 
{
     real	*xarr_fc, *yarr_fc ;
     real	*xarr_fo, *yarr_fo ;
     real	fcwidth, fowidth ;
     real	rat_oc, epsi ;
     real	diff = 0 ;
     real	correction ;
     float	origin, interval ;
     char	in_filename[MAXSTRING] ;
     int	Nfc, Nfo ;
     int	nc_skip = 1, no_skip = 1 ;
     int	nc_start = 0, no_start = 0;
     int	nc,no;
     int	count = 0 ;
     int	n ;
     char	good ;

     sprintf(in_filename, "%s", name_fc) ;
     Nfc = read_it_in_asc(in_filename, &xarr_fc, &yarr_fc) ;

     sprintf(in_filename, "%s", name_fo) ;
     Nfo = read_it_in_asc(in_filename, &xarr_fo, &yarr_fo) ;

     /*******************************************************
     Scaling fits two _wil plots which may have different 
     origins, lengths, binwidths ... these must be checked.
     First: do their ranges overlap?
     Then: are their binwidths compatible?
     *******************************************************/

     if ((*xarr_fo >= *(xarr_fc+Nfc-1)) ||
         (*xarr_fc >= *(xarr_fo+Nfo-1))) 
          EdenError("Wilson files do not have overlapping ranges!") ;

     fcwidth = (*(xarr_fc+1) - *xarr_fc) ;
     fowidth = (*(xarr_fo+1) - *xarr_fo) ;
     rat_oc = fowidth /fcwidth ; 
     epsi = 0.01 * fcwidth ;

     if (rat_oc >= 1.) {
         for (n = 1, good = FALSE; n <= 16; n++)
             if (fabs(rat_oc - n) < epsi) {
                 good = TRUE ;
                 nc_skip = n ;
                 break ;			/* widths are compatible */
             }
     }
     else {
         for (n = 1, good = FALSE; n <= 16; n++)
             if (fabs((1./rat_oc) - n) < epsi) {
                 good = TRUE ;
                 no_skip = n ;
                 break ;			/* widths are compatible */
             }
     }

     if (!good) {
	  sprintf(message, "fo binwidth is %g, fc binwidth is %g", 
		  fowidth, fcwidth) ;
	  printTwice(message) ;
	  EdenError("Fo and fc Wilson arrays are incompatible") ;
     }

     /* find first common point in xarr_fc and xarr_fo */

     for (nc_start = 0; nc_start < Nfc; nc_start++) {
	 for (no_start = 0; no_start < Nfo; no_start++) {
	     if ((*(xarr_fo + no_start) >= *(xarr_fc + nc_start))) 
	          break ; 
          }

	  if (fabs((*(xarr_fo + no_start) - *(xarr_fc + nc_start))) < epsi)
	       break ;
      }
      if ((no_start > Nfo/2) || (nc_start > Nfc/2))
	   EdenError(
	   "Couldn't find suitable common starting index for scaling.") ;


     for (nc = nc_start, no = no_start; (nc < Nfc && no < Nfo); 
          nc+= nc_skip,  no+= no_skip) {

          diff += *(yarr_fc + nc) - *(yarr_fo + no) ;
	  count++ ;
     }

     correction = diff / (float) count ;
     origin = *xarr_fo ;
     interval = fowidth ;

     sprintf(message, "\tProposed value for fscale is %g", 
             exp(correction / 2.)) ;
     printTwice(message) ;

     for (n = 0; n < Nfo; n++) 
	  *(yarr_fo + n) += correction ;

     send_it_out_asc("scaled_wil", yarr_fo, Nfo, origin, interval) ;

     if (graphics) {
        sprintf(message, "xmgrace %s %s -legend load", name_fo, name_fc) ;
        system(message) ;

        sprintf(message, "xmgrace %s %s -legend load", "scaled_wil", name_fc) ;
        system(message) ;
     }
     else {
	sprintf(message, "You can use a graphics program to plot %s and %s",
		"scaled_wil", name_fc) ;
        printTwice(message) ;
     }

     return ;
}

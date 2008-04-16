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

                                DISTANCE

  Title:        Calculate distances among files of electron/pixel.    
  Author:       Hanna Szoke
  Date:         29/12/94
  Function:     This program reads M (in range 2 - 9) sets of binary files 
		representing (unregridded) electron/pixels.  For each 
		pair, it calculates and reports 

		* the rms fractional distances among them;
		* the linear absolute fractional distances among them; and
		* the correlation coefficients among them.
		
                Input parameters are read locally. UTIL.C has worker functions.

*******************************************************************************/

#include "util.h"
#include "dims.h"

#define MAXNFILES	9	/* max # of files for distance calculation */
#define MINNFILES	2	/* min # of files for distance calculation */

static	char	dfilename[MAXNFILES][MAXSTRING] ;
static	real	*dnump[MAXNFILES] ;	/* electron densities */
static	int	M ;			/* number of Np file sets to measure */

static  int	dNhr[MAXNFILES] ;
static	real	*hrNump[MAXNFILES] ;
static	real	*hrNx[MAXNFILES] ;
static	real	*hrNy[MAXNFILES] ;
static	real	*hrNz[MAXNFILES] ;

static	void	do_cc() ;
static	void	do_lin() ;
static	void	do_rms() ;
static	real	rms_fracdist(real *, real *, real *, real *) ;
static	real	lin_fracdist(real *, real *, real *, real *) ;
static	real	corr_coef(real *, real *, real *, real *) ;
static	void	printMatrix(real (* )()) ;


void	distance_main(int argc, char *argv[])
{
     char	list_filename[MAXSTRING] ;
     void	readEpixFile(char *, real *) ;
     int	hrGetLength(char *) ;
     void	hrRenameList(real **, real **, real **) ;
     void	hrReadLists(char *, real *) ;		
     void	hrCompareIndices(real *, real *, real *) ;

     int	m ;			/* counters of files */

  /***********************************************************
   Check for right ballpark, say hello, identify file names.
  ***********************************************************/

     M = argc - optind - 2 ;

     if (M > MAXNFILES) 
	  EdenError("Too many files!") ;
     else if (M >= MINNFILES)
          for (m = 0; m < M; m++) 
               sprintf(dfilename[m], "%s", argv[optind+m+2]) ;

     if (M == 0) {
	  prompt("\nPlease enter 2 - 8 file names: ") ;
          sscanf(terminp, "%s %s %s %s %s %s %s %s", dfilename[0], dfilename[1],
	         dfilename[2], dfilename[3], dfilename[4], dfilename[5], 
	         dfilename[6], dfilename[7]) ;

          M = m = 0 ;
          while ((int)strlen(dfilename[m]) > 0)   {
	       M++; 
	       m++;
          }
          if (M < MINNFILES) 
	       EdenError("Too few files!") ;

	  sprintf(message, "Running %s on %d files", caller, M) ; 
	  printTwice(message) ;
     }

     hello(caller) ;

  /************************************************************
  Fetch input conditions, set up arrays in physical space.
  Read the binary el/voxel maps. 
  ************************************************************/

     readBasicInput() ;

     for (m = 0; m < M; m++) 
          dnump[m] = (real *) e_malloc(Nptotal*sizeof(real), "distance") ;

     for (m = 0; m < M; m++) {  
          readEpixFile(dfilename[m], dnump[m]) ;

          if (HighRes) {
               sprintf(list_filename, "%s.list", dfilename[m]) ;
               dNhr[m] = hrGetLength(list_filename) ;	
	       if (m > 0) 
	            if (dNhr[m] != dNhr[0])
	                 EdenError(
               "Your high-res list files do NOT have the same length!") ;

               hrNump[m]= (real *)e_malloc(dNhr[m]*sizeof(real), "distance") ;
               hrReadLists(list_filename, hrNump[m]) ;	
	       hrRenameList(&hrNx[m], &hrNy[m], &hrNz[m]) ;

	       if (m > 0) 
                    hrCompareIndices(hrNx[m], hrNy[m], hrNz[m]) ;
          }
     }
  /************************************************************
  Compute and report distances.
  ************************************************************/

     do_rms() ;
     do_lin() ;
     do_cc() ;

}

void	do_rms() 
{
  /************************************************************
  Compute and report rms distances.
  ************************************************************/

     printTwice(
       "\nThe rms fractional distance between 2 el/pix files - A and B - is\n") ;
     printTwice( "\tsqrt[Sum[(a - b)^2]] / (sqrt[Sum[a^2]] + sqrt[Sum[b^2]])\n ") ;
     printTwice( "where a and b are data values from files A and B") ;
     printTwice( "and Sum[] is a sum over all such elements.\n") ;

     if (M == 2) {
          sprintf(message, 
	  "The rms fractional distance between data in %s and %s is %g\n",
          dfilename[0], dfilename[1], 
	  rms_fracdist(dnump[0], dnump[1], hrNump[0], hrNump[1]));
          printTwice(message) ;
     }
     else {
	  printTwice("Rms fractional distances among the files.\n") ;

          printMatrix(rms_fracdist) ;

     }
}
real	rms_fracdist(real *targ1, real *targ2, real *hr1, real *hr2) 
{
     real 	*ptarg1 ;
     real 	*ptarg2 ;
     int 	n ;
     real	sumsqdiff = 0 ;
     real	sumsqa = 0, sumsqb = 0 ;
     real	a, b,  diff ;

     ptarg1 = targ1 ;
     ptarg2 = targ2 ;

     for (n = 0; n < Nptotal; n++) {
          a = *(ptarg1++) ;
          b = *(ptarg2++) ;
	  diff = a - b ;
	  sumsqdiff += diff * diff ;
	  sumsqa    += a * a ;
	  sumsqb    += b * b ;
     } 
     if (HighRes) {
	  for (n = 0; n < dNhr[0]; n++) {
               a = *(hr1++) ;
               b = *(hr2++) ;
	       diff = a - b ;
	       sumsqdiff += diff * diff ;
	       sumsqa    += a * a ;
	       sumsqb    += b * b ;
          }
     }
     return (sqrt(sumsqdiff) / (sqrt(sumsqa) + sqrt(sumsqb))) ;

}

void	do_lin() 
{
  /************************************************************
  Compute and report linear absolute fractional distances.
  ************************************************************/

     printTwice( 
  "\nThe linear absolute fractional distance between 2 el/pix files - A and B - is\n") ;
     printTwice( "\tSum[|a - b|] /  Sum[a + b]\n ") ;
     printTwice( "where a and b are data values from files A and B") ;
     printTwice( "and Sum[] is a sum over all such elements.\n") ;

     if (M == 2) {
          sprintf(message, 
	  "The linear abs. fractional distance between data in %s and %s is %g\n",
          dfilename[0], dfilename[1], 
	  lin_fracdist(dnump[0], dnump[1], hrNump[0], hrNump[1])) ;
          printTwice(message) ;
     }
     else {
	  printTwice(
        "The linear absolute fractional distances among the files.\n") ;

          printMatrix(lin_fracdist) ;
     }
}
real	lin_fracdist(real *targ1, real *targ2, real *hr1, real *hr2) 
{
     real 	*ptarg1 ;
     real 	*ptarg2 ;
     int 	n ;
     real	sum = 0, diff = 0 ; 

     ptarg1 = targ1 ;
     ptarg2 = targ2 ;

     for (n = 0; n < Nptotal; n++) {
	  diff += fabs(*(ptarg1) - *(ptarg2)) ;
	  sum += (*(ptarg1++) + *(ptarg2++)) ;
     } 
     if (HighRes) {
	  for (n = 0; n < dNhr[0]; n++) {
	       diff += fabs(*(hr1) - *(hr2)) ;
	       sum += (*(hr1++) + *(hr2++)) ;
          }
     }
     return (diff / sum) ; 

}

void	do_cc() 
{
  /************************************************************
  Compute and report correlation coefficients.
  ************************************************************/

     printTwice(
     "\nThe correlation coefficient between 2 el/pix files - A and B - is\n") ;
     printTwice( "\tSum[del(a) * del(b)] / sqrt[Sum[del(a)^2] * Sum[del(b)^2)]]\n") ;

     printTwice( "where del(a) is a - a_ave and del(b) is b - b_ave") ;
     printTwice( "and Sum[] is a sum over all such elements.\n") ;

     if (M == 2) {
          sprintf(message, 
	  "The correlation coefficient between data in %s and %s is %g\n",
          dfilename[0], dfilename[1], 
	  corr_coef(dnump[0], dnump[1], hrNump[0], hrNump[1])) ;
          printTwice(message) ;
     }
     else {
	  printTwice( "The correlation coefficients among the files.\n") ;
          printMatrix(corr_coef) ;
     }
     return ;
}
real	corr_coef(real *targ1, real *targ2, real *hr1, real *hr2) 
{
     real 	*ptarg1 ;
     real 	*ptarg2 ;
     real	mean1 = 0 ;
     real	mean2 = 0 ;
     real	num = 0 ;
     real	denom1 = 0, denom2 = 0 ;
     real	d1, d2 ;
     int 	n ;
     int	N ;

     for (n = 0, ptarg1 = targ1; n < Nptotal; n++, ptarg1++) 
	  mean1 += *ptarg1 ;

     N = Nptotal; 

     if (HighRes) {
	  for (n = 0; n < dNhr[0]; n++) 
	       mean1 += *(hr1+n) ;
          N += dNhr[0] ;
     }
     mean1  *= 1./(float) N ;

     for (n = 0, ptarg2 = targ2; n < Nptotal; n++, ptarg2++) 
	  mean2 += *ptarg2 ;

     if (HighRes) 
	  for (n = 0; n < dNhr[0]; n++) 
	       mean2 += *(hr2+n) ;

     mean2  *= 1./(float) N;

     for (n = 0, ptarg1 = targ1, ptarg2 = targ2; n < Nptotal; 
	  n++, ptarg1++, ptarg2++) {

	  d1 = *ptarg1 - mean1 ;
	  d2 = *ptarg2 - mean2 ;
	  num += d1 * d2 ;
	  denom1 += d1 * d1;
	  denom2 += d2 * d2 ;
     }
     if (HighRes) 
          for (n = 0; n < dNhr[0]; n++) {

	       d1 = *(hr1+n) - mean1 ;
	       d2 = *(hr2+n) - mean2 ;
	       num += d1 * d2 ;
	       denom1 += d1 * d1;
	       denom2 += d2 * d2 ;
          }

     return (num / sqrt(denom1*denom2) ) ;
}

void printMatrix(real (* distfunc)(real *, real *, real *, real *))
{
     int	m, m1, m2 ;		/* counters of files */
     int	l, len ;		/* counters of length message buffer */

     real	average = 0 ;
     real	maver[MAXNFILES] ;
     real	value;			/* value returned by rms_fracdist(() */
     char	alpha[MAXNFILES] ;

  /************************************************************
  Compute and report distances.
  ************************************************************/

	  for (m = 0; m < M; m++) 
	       alpha[m] = 'a' + m ;
	  for (m = 0, l = 0, len = 0; m < M; m++, l+=len) 
	       len = sprintf(message+l, "\t%5s%c", "     ", alpha[m]) ;
          printTwice(message) ;

	  for (m1 = 0; m1 < M; m1++) {
	       len = sprintf(message, "%5s%c", "     ", alpha[m1]) ;
	       for (m2 = 0,l=len, maver[m1] = 0; m2 < M; m2++,l+=len) {
		    value = distfunc
			 (dnump[m1], dnump[m2], hrNump[m1], hrNump[m2]) ;
		    if (m2 > m1) 
			 average += value ;
		    if (m2 != m1) 
			 maver[m1] += value ;
		    len = sprintf(message+l, "\t%6.4f", value) ;
               }
               printTwice(message) ;
          }

          len = sprintf(message, "\n   ave") ;

       	  for (m2 = 0,l=len; m2 < M; m2++,l+=len) {
	       len = sprintf(message+l, "\t%6.4f", maver[m2] / (M-1)) ;
          }
          printTwice(message) ;

	  sprintf(message, "\nGrand average = %g\n", average / (M*(M-1)/2)) ;
          printTwice(message) ;

	  printTwice("Where:") ;
	  for (m = 0; m < M; m++) {
	       sprintf(message, "\t%c%s%s", 
		    alpha[m], " stands for ", dfilename[m]) ;
               printTwice(message) ;
          }
}

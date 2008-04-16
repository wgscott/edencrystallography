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

                                HIGHRES.C

  Title:        HIGHRES 
  Author:       Hanna Szoke
  Date:         Oct. 13, 2000
		Mar 12, 2002 - revised extensively!
  Function:     Deals with the high resolution option in Eden.

void	hrCompareIndices() compares two sets of indices
void	hrDoregSFT does the slow Fourier transform for Regrid.
void	hrDoSFT does the slow Fourier transform for Solve and Forth.
int	hrGetLength() determined Nhr of list to be read.
void	hrGradSFT() computes the gradients.
int	hrInitialize() sets up high-resolution processing in Solve.
void	hrMaskOut() applies hrMask to its input Nptotal array. 
void	hrPrepareSFT() prepares for slow Fourier transform in Solve and Forth.
void	hrReadLists() reads lists.
void	hrRenameList() 
void	hrResetIndices() uses spread factor to change the indices (for Regrid).
void	hrSetup() returns calculated hr values.
void	hrSymmetrize() applies symmetrization.
void	hrWriteLists() writes a given list, hrNx, hrNy, hrNz. 

int	hrCheckPt() adds a new pt to the list (internal).
int	hrCrysSym() determines points related by crystal symmetry (internal).
int	hrMakeLists() sets up hrNx, hrNy, hrNz; determines Nhr (internal).
int	hrMakeMask() sets up hrMask.  returns hrM (internal).
int	hrMatch() checks pt for equality of indices with existing pts (internal).
real	hrSetCutoff() returns a cutoff value (internal).
*******************************************************************************/

#define        ROUNDOFF        1.e-5  /* ... to remove clutter in hrMatch */

#include "util.h"		/* general definitions */
#include "cellparams.h"		/* a, b, c, angles, crystal_symmetry */
#include "symmetry.h"		/* symmetry group information */
#include "dims.h"		/* incl. HighRes, Nhr and Npextended  */

		/* Data Structures */

int	Nhr ; 			/* length (int) of hr array(s)	*/
static  real	average ;	/* used in hrSetCutoff() & hrAnalyze()	*/

static	char	*hrMask ;		/* masks out points > hrCutoff * <vox>	*/
static  int	maxNhr ;
static	real	*hrNx ;			/* non-integer indices in x	*/
static	real	*hrNy ;			/* non-integer indices in y	*/
static	real	*hrNz ;			/* non-integer indices in z	*/
static	real	*hrexpfac ;		/* exponential factors for SF	*/
static	int	*hr_csind ;

static	real	hrSetCutoff(real *) ;
static	int	hrMakeMask(real, real *) ;
static	int	hrMakeLists(real *, real *) ;
static	void	hrAnalyze(real *) ;
static	int	hrCheckPt(real *, real, real, real, int) ;
static	int	hrCrysSym(real *) ;
static	int	hrMatch(real, real, real) ;

void	hpsort(unsigned long, float []) ;

	/************************************************************
	Set up high-resolution processing.  Called from qshare.c,
	initNparrays().  
	************************************************************/

int	hrInitialize(real *array) 
{
     real	cutoff;
     int	hrM;

     cutoff = hrSetCutoff(array) ;
     hrM = hrMakeMask(cutoff, array);

     if ((hrM == 0) || (hrM > 5))
	  hrAnalyze(array) ;

     if (hrM == 0)
          EdenError("Quitting - There are no high-resolution points!") ; 
     else {
          sprintf(message, "There are %d high-resolution points.", hrM) ;
          printTwice(message) ;
     }

	/************************************************************
	Allocate maximum possible space for list arrays (7*hrM).
	In fact, only the first Nhr elements will be used.
	************************************************************/

     maxNhr = 7 * hrM * Ncs ;
     sprintf(message, "hrInitialize") ;
     hrNx =    (real *)e_malloc(maxNhr*sizeof(real), message) ;
     hrNy =    (real *)e_malloc(maxNhr*sizeof(real), message) ;
     hrNz =    (real *)e_malloc(maxNhr*sizeof(real), message) ;

     return(maxNhr) ;
}


void	hrSetup(real *array)
{
     void	hrWriteLists(char *, real *) ;

     int	added ;
     real	*inarray, *outarray ;

     inarray = array ;
     outarray = array + Nptotal ;

     Nhr = hrMakeLists(inarray, outarray) ;
     added = hrCrysSym(outarray) ;

     sprintf(message, "Nhr = %d, added = %d\n", Nhr, added) ;
     printTwice(message) ;

     Nhr += added ;
     Npextended = Nptotal + Nhr ;

     hrWriteLists("initial.list", outarray) ;

     return ;
}

	/************************************************************
	Return value of cutoff, based on user-definable hrCutoff.  
	Function is internal to highres.c.
	************************************************************/
real	hrSetCutoff(real *array) 
{
     int	n ;

     for (n = 0, average = 0; n < Nptotal; n++, array++)
          average += *array ;

     average /= Nptotal ;
     sprintf(message, "average is %g",average);
     printTwice(message) ;

     return (average * hrCutoff) ; 
}

	/************************************************************
	User's input HRCUTOFF was not appropriate; analyze array, 
	printing out 10 top values and the hrCutoff that would 
	correspond to them.  See findCutoff() in maketar.c
	************************************************************/

void	hrAnalyze(real *array) 
{
     float	*copy, *pcopy ;		/* array to use for sorting in place */
     real	*parray ;
     int	n ;
     float	hrval ;
     unsigned long	Nun ;

     /* allocate space, copy model data. */

     copy  = (float *) e_malloc(Nptotal*sizeof(float), "hrAnalyze") ;

     for (n = 0, pcopy = copy, parray = array; n < Nptotal; n++, parray++, pcopy++)
	  *pcopy = *parray ;

     /*  Sort in ascending order */

     Nun = Nptotal ;
     hpsort(Nun, copy) ;

     for (n = 1; n <= 10; n++) {
	  hrval = *(copy + Nptotal - n) / average ;
          sprintf(message, "for %d HighRes points, use HRCUTOFF = %4.1f", 
		  n, hrval) ;
	  printTwice(message) ;
     }

     e_free(Nptotal*sizeof(float), copy) ;
     return ; 
}


	/************************************************************
	Set up hrMask, return hrM. Internal to highres.c.
	************************************************************/

int	hrMakeMask(real cutoff, real *array) 
{
     int	n ;
     int	hrM ;
     char	*phrMask;

     hrMask = (char *)e_malloc( Nptotal*sizeof(char), "hrMakeMask") ;

     for (n = 0, phrMask = hrMask; n < Nptotal; n++, phrMask++)
          *phrMask = 1 ;

     for (n = 0, phrMask = hrMask, hrM = 0; n < Nptotal; 
          n++, phrMask++, array++)
          if (*array > cutoff) {
               *phrMask = 0 ;
               hrM++ ;
          }

     return(hrM);
}

	/************************************************************
	Using the hrMask as a guide, set up hrNx, hrNy, hrNz.  
	Check for redundant entries ; determine Nhr where 
	Nhr <= 7*hrM (allocated lengths of the arrays). 
	>>> This is where the specific algorithm is used!
	hrMakeLists() is internal to highres.c.
	inarray is first Nptotal entries;
	outarray is SAME array, starting at Nptotal.
	************************************************************/

int	hrMakeLists(real *inarray, real *outarray) 
{
     real	fclip(real, real) ;
     int	n ;
     char	*pmask ;
     real	thisx, thisy, thisz ;
     real	thisxmin, thisymin, thiszmin ;
     real	thisxplu, thisyplu, thiszplu ;

     for (n = 0, pmask = hrMask, Nhr = 0; n < Nptotal; n++, pmask++, inarray++) {

          if (*pmask == 0) {			/* hr point	*/
               *(outarray + Nhr) = *inarray ;

               if (n < Np) {			/* simple lattice */
	            *(hrNx+Nhr) = thisx =  n % Nx ;
                    *(hrNy+Nhr) = thisy =  (n/Nx) % Ny ;
                    *(hrNz+Nhr) = thisz =  n/(Nx*Ny) ;
               }
               else {				/* body-centered lattice */
	            *(hrNx+Nhr) = thisx =  n % Nx + 0.5 ;
                    *(hrNy+Nhr) = thisy =  (n/Nx) % Ny + 0.5 ;
                    *(hrNz+Nhr) = thisz =  n/(Nx*Ny) % Nz + 0.5 ;
               }

               thisxmin = fclip(thisx - 0.5, (real) Nx) ;
               thisymin = fclip(thisy - 0.5, (real) Ny) ;
               thiszmin = fclip(thisz - 0.5, (real) Nz) ;
               thisxplu = fclip(thisx + 0.5, (real) Nx) ;
               thisyplu = fclip(thisy + 0.5, (real) Ny) ;
               thiszplu = fclip(thisz + 0.5, (real) Nz) ;

               Nhr ++ ;
               Nhr += hrCheckPt(outarray, thisxplu, thisy,    thisz,    Nhr) ;
               Nhr += hrCheckPt(outarray, thisxmin, thisy,    thisz,    Nhr) ;
               Nhr += hrCheckPt(outarray, thisx,    thisyplu, thisz,    Nhr) ;
               Nhr += hrCheckPt(outarray, thisx,    thisymin, thisz,    Nhr) ;
               Nhr += hrCheckPt(outarray, thisx,    thisy,    thiszplu, Nhr) ;
               Nhr += hrCheckPt(outarray, thisx,    thisy,    thiszmin, Nhr) ;
          }
     }
     return(Nhr) ;
}

	/************************************************************
	Check a new point against all previously established points;
	If it's new, add it to the lists, and return 1 to increment 
	Nhr.  Otherwise, return 0.  Internal to highres.c.
	************************************************************/

int	hrCheckPt(real *outarray, real newx, real newy, real newz, int N)
{
     int	nhr ;
     real	*pNx, *pNy, *pNz ;

     pNx = hrNx ;
     pNy = hrNy ;
     pNz = hrNz ;
 
     for (nhr = 0; nhr < N; nhr++, pNx++, pNy++, pNz++) {

         if ((newx==*pNx) && (newy==*pNy) && (newz == *pNz))
              return (0) ;
     }

     *(outarray + N) = 0 ;
     *(hrNx + N) = newx ;
     *(hrNy + N) = newy ;
     *(hrNz + N) = newz ;

     return (1) ;
}


int	hrMatch(real fx, real fy, real fz)
{
     int	n ;
     real	*pNx, *pNy, *pNz ;

     pNx = hrNx ;
     pNy = hrNy ;
     pNz = hrNz ;

     for (n = 0; n < Nhr; n++, pNx++, pNy++, pNz++) {

	  if ((fabs(*pNx/Nx - fx) < ROUNDOFF) && 
	      (fabs(*pNy/Ny - fy) < ROUNDOFF) && 
	      (fabs(*pNz/Nz - fz) < ROUNDOFF)) 
	           return(n) ;			/* accounted for */
     }

     return(-1) ;				/* new */
}

	/************************************************************
	Find points related by crystal symmetry to the hr list points.
	This function is similar to setup_csind() in crysutil.c.    
	************************************************************/
int	hrCrysSym(real *outarray)
{
     int	i, n ;
     real	*mat ;
     int	newpt ;
     int	added = 0;
     real	*pNx, *pNy, *pNz ;
     real	newfx, newfy, newfz ;
     char	*match ;

     pNx = hrNx ;
     pNy = hrNy ;
     pNz = hrNz ;

     hr_csind = (int *) e_malloc(maxNhr*sizeof(int), "hrCrysSym") ;
     match    = (char *) e_malloc(maxNhr*sizeof(char), "hrCrysSym") ;

     for (n = 0; n < maxNhr; n++)  {
	  *(hr_csind + n) = n ;
	  *(match + n) = 0 ;
     }

     for (n = 0; n < Nhr; n++, pNx++, pNy++, pNz++) 

          if (*(match + n) == 0) {

	       for(i = 0, mat = matop; i < Ncs; i++, mat += MEL) {
		    apply_symop(*pNx/Nx, *pNy/Ny, *pNz/Nz, mat, 
				&newfx, &newfy, &newfz) ;
                    newpt = hrMatch(newfx, newfy, newfz) ;

		    if (newpt < 0) {
			 *(hr_csind + Nhr + added) = n ;
			 *(hrNx + Nhr + added) = newfx*Nx ;
			 *(hrNy + Nhr + added) = newfy*Ny ;
			 *(hrNz + Nhr + added) = newfz*Nz ;
			 *(outarray + Nhr + added) = *(outarray + n) ;
		         *(match + Nhr + added) = 1 ;
			 added++ ;
                    }

		    if (very_verbose) {
		       if (newpt < 0)
            fprintf(fp_log, "old point = %g, %g, %g, newpt = %g, %g, %g\n",
			 *pNx, *pNy, *pNz, newfx*Nx, newfy*Ny, newfz*Nz) ;
		       else
                           fprintf(fp_log, "match already set for %d\n", newpt) ;
                    }
               }
          }

     if (very_verbose) {

          for (n = 0; n < Nhr+added; n++)
             fprintf(fp_log, "n = %d, index = %d\n", n, *(hr_csind + n)) ;

          fprintf(fp_log, "maxNhr=%d, Ncs=%d, Nhr=%d\n", maxNhr, Ncs, Nhr) ;
     }

     e_free(maxNhr*sizeof(char), match) ;
     return(added) ;
}

void	hrSymmetrize(real *array)
{
     int	nhr, index ;
     int	*phrcs ;
     real	*parray ;
     real	*hrsums ;
     int	*hrcounts ;

     sprintf(message, "hrSymmetrize") ;

     hrsums = (real *)e_malloc(Nhr*sizeof(real), message) ;
     hrcounts = (int *)e_malloc(Nhr*sizeof(int), message) ;

     for (nhr = 0; nhr < Nhr; nhr++) {
	  *(hrsums+nhr) = 0 ;
	  *(hrcounts+nhr) = 0 ;
     }
      
     parray = array ;
     phrcs = hr_csind ;

     for (nhr = 0; nhr < Nhr; nhr++, parray++, phrcs++) {
          index = *phrcs ;
	  *(hrsums + index) += *(parray) ;
	  *(hrcounts + index) += 1 ;
     }


     for (nhr = 0; nhr < Nhr; nhr++) 
	  if (*(hrcounts+nhr) > 0)
	       *(hrsums+nhr) /= *(hrcounts+nhr) ;

     parray = array ;
     phrcs = hr_csind ;

     for (nhr = 0; nhr < Nhr; nhr++, parray++, phrcs++) {
          index = *phrcs ;
	  if (very_verbose) {
               sprintf(message, " nhr=%d, array = %g, average = %g", 
               nhr, *(parray), *(hrsums + index)) ;
               printTwice(message) ;
          }  
	  *(parray) = *(hrsums + index) ;
     }
     e_free(Nhr*sizeof(real), hrsums) ;
     e_free(Nhr*sizeof(int), hrcounts) ;

}

	/************************************************************
	Prepare for slow Fourier transform:
	there is no hrSetupexpfac(); we use setupexpfac(d),  
		d = (1/4)*delfac.  
	We can use this for the hr points regardless of whether 
	they are derived from simple or body-centered lattice pts.
	************************************************************/

void	hrPrepareSFT() 
{
     real     *setupexpfac(real) ;

     hrexpfac = setupexpfac(delfac/4.) ;

     return ;
}

	/************************************************************
	Do the slow Fourier Transform: accumulate into O(h, k, l)
        O.re = hrexpfac*S[nump * cos(pix*h+piy*k+piz*l) 
        O.im = hrexpfac*S[nump * sin(pix*h+piy*k+piz*l)
	where S is sum over Nhr.  O is complex, input and output.
	*array values are nump's from Nptotal - Npextended.
	hrDoSFT() is called from solve.c, forth.c, ssetup.c and cost.c.
	************************************************************/

void	hrDoSFT(real *array, COMPLEX *O)
{
     int	n, nhr, h, k, l ;
     real	pix, piy, piz;

     for (nhr = 0; nhr < Nhr; nhr++, array++) {

        if (*array != 0) {

           pix = *(hrNx+nhr) / floatNx ;
           piy = *(hrNy+nhr) / floatNy ;
           piz = *(hrNz+nhr) / floatNz ;

           for (h = 0; h < Nh; h++) {
              for (k = 0; k < Nk; k++) {
                 for (l = 0; l < Nl; l++) {

                    n = h + k*Nh + l*Nh*Nk ;
                    (O + n)->re += *(hrexpfac+n) * *array * 
                                   cos(TWOPI*(h*pix + k*piy + l*piz)) ;
                    (O + n)->im += *(hrexpfac+n) * *array *
                                   sin(TWOPI*(h*pix + k*piy + l*piz)) ;

                 }
              }
           }
        }
     }
}

	/************************************************************
	Do the slow Fourier Transform: accumulate into O(h, k, l).
	This form of SFT is used by regrid and forth.
	************************************************************/

void	hrDoregSFT(COMPLEX *Rhigh, real *array)
{
     int	n, nhr, h, k, l ;
     real	thispt ;
     real	pix, piy, piz;

     for (nhr = 0; nhr < Nhr; nhr++) {

        if (*array != 0) {

           pix = *(hrNx+nhr) / (float) Nx ;
           piy = *(hrNy+nhr) / (float) Ny ;
           piz = *(hrNz+nhr) / (float) Nz ;

           thispt = *(array + nhr) ;

           for (h = 0; h < Nx; h++) {
              for (k = 0; k < Ny; k++) {
                 for (l = 0; l < Nz; l++) {

                    n = h + k*Nx + l*Nx*Ny ;

                    (Rhigh + n)->re += thispt * 
				       cos(TWOPI*(h*pix + k*piy + l*piz)) ;
                    (Rhigh + n)->im += thispt * 
				       sin(TWOPI*(h*pix + k*piy + l*piz)) ;

                 }
              }
           }
        }
     }
}

	/************************************************************
	Compute the gradients from the O which includes hr
	contributions from hrDoSFT().  Calculate, for n = 1,...,  Nhr
	grad(n) = Sum_over_Nhkl [O(n).re*cos() + O(n).im*sin()]
	hrGradSFT() is called from cost.c, e_funct_hkl().
	************************************************************/

void	hrGradSFT(COMPLEX *O, real *grad) 
{
     int	n, nhr, h, k, l ;
     real	pix, piy, piz;
     real	*pgrad, arg;
     real	Cabs() ;

     for (nhr = 0, pgrad = grad; nhr < Nhr; nhr++, pgrad++) 
	  *pgrad = 0.0 ;

     for (n = 0; n < Nhkl; n++) 
          if (Cabs(*(O+n)) != 0) {

               fetch_hkl(n, &h, &k, &l) ;

               for (nhr = 0, pgrad = grad; nhr < Nhr; nhr++, pgrad++) {

                    pix = *(hrNx+nhr) / floatNx ;
                    piy = *(hrNy+nhr) / floatNy ;
                    piz = *(hrNz+nhr) / floatNz ;

	            arg = h*pix + k*piy + l*piz ;

                    *pgrad += ((O + n)->re * cos(TWOPI*arg) +
                               (O + n)->im * sin(TWOPI*arg)) * *(hrexpfac+n) ;

               }
          }
}


	/************************************************************
	Apply hrMask to input array.  Called from cost.c (grad),
	cost.c (maxp) and qshare.c (knownp).
	************************************************************/

void	hrMaskOut(real *array) 
{
      int	n ;
      char	*phrMask ;

      for (n = 0, phrMask = hrMask; n <Nptotal; n++, array++, phrMask++)
         *array *= *phrMask ;

}
	/************************************************************
	Write an hr list incl. array values, hrNx, hrNy, hrNz.  
	Called within highres.c and by solve.c and addmaps.c.
	************************************************************/

void	hrWriteLists(char *hrname, real *array)
{
     FILE	*fpout ;
     int	nhr ;
     real     *phr ;
     real     *pix, *piy, *piz ;

     if ((fpout = fopen(hrname, "w")) == NULL)
     {
          sprintf(message, "hrWriteLists: Cannot open %s", hrname);
          EdenError(message) ;
     }
     fprintf(fpout, "%d\n", Nhr);
     
     for (nhr = 0, phr = array, pix = hrNx, piy = hrNy, piz = hrNz; 
          nhr < Nhr; nhr++, phr++, pix++, piy++, piz++) 
     
          fprintf(fpout, "%g\t %5.2f %5.2f %5.2f\n", 
                          *phr, *pix, *piy, *piz) ;

     fclose(fpout) ;
}

int	hrGetLength(char *hrname) 
{
     FILE	*fpin ;

     if ((fpin = fopen(hrname, "r")) == NULL)
     {
          sprintf(message, "hrGetLength: Cannot open %s", hrname);
          EdenError(message) ;
     }

     while (fgets(nextline, MAXSTRING, fpin) != NULL) {
          sscanf(nextline, "%d", &Nhr) ;
	  if (Nhr > 0)
	       break ;
     }
     
     fclose(fpin) ;
     return (Nhr) ;
}

	/************************************************************
	Read hrNump, hrNx, hrNy, hrNz.  Transfer hrNump to array
	************************************************************/

void	hrReadLists(char *hrname, real *array) 
{
     FILE	*fpin ;
     char	first = TRUE;
     int	n = 0 ;

     if ((fpin = fopen(hrname, "r")) == NULL)
     {
          sprintf(message, "hrReadLists: Cannot open %s", hrname);
          EdenError(message) ;
     }

     while (fgets(nextline, MAXSTRING, fpin) != NULL) {
          if (first) {
              first = FALSE ;
              sscanf(nextline, "%d", &Nhr) ;
              sprintf(message, "hrReadLists") ;
              hrNx =   (real *)e_malloc(Nhr*sizeof(real), message) ;
              hrNy =   (real *)e_malloc(Nhr*sizeof(real), message) ;
              hrNz =   (real *)e_malloc(Nhr*sizeof(real), message) ;
          }
          else {
             sscanf(nextline, "%g %g %g %g", &t0, &t1, &t2, &t3) ;
	     *(array+n) = t0 ;
             *(hrNx+n) = t1 ;
             *(hrNy+n) = t2 ;
             *(hrNz+n) = t3 ;
             n++ ;
          }
     }
     return ;
}

	/************************************************************
	Increase hrNx, hrNy, hrNz by spread_factor, for use in Regrid
	************************************************************/

void	hrResetIndices(int spread_factor) 
{
     int	n;

     for (n = 0; n <Nhr; n++) {
        *(hrNx + n) *= spread_factor ;
        *(hrNy + n) *= spread_factor ;
        *(hrNz + n) *= spread_factor ;
     }
}
	/************************************************************
	Transfer array pointers, for use in Addmaps. 
	************************************************************/

void	hrRenameList(real **hrx, real **hry, real **hrz) 
{
     *hrx   = hrNx ;
     *hry   = hrNy ;
     *hrz   = hrNz ;
}


void	hrCompareIndices(real *hrNx1, real *hrNy1, real *hrNz1)
{
     int	n ;

     for (n = 0; n < Nhr; n++)
	  if ((*(hrNx1 + n) != *(hrNx + n)) ||
	      (*(hrNy1 + n) != *(hrNy + n)) ||
	      (*(hrNz1 + n) != *(hrNz + n))) {

	  fprintf(stderr, "n=%d, nx1=%g, nx2=%g\n", n, *(hrNx+n), *(hrNx1+n)) ;
	  fprintf(stderr, "n=%d, ny1=%g, ny2=%g\n", n, *(hrNy+n), *(hrNy1+n)) ;
	  fprintf(stderr, "n=%d, nz1=%g, nz2=%g\n", n, *(hrNz+n), *(hrNz1+n)) ;
	      EdenError(
      " The high-resolution indices in your files do not agree - quitting.") ;
     }
}
     

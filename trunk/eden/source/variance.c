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

				VARIANCE

  Title:        Calculate stability among files of electron/pixel.    
  Author:       Hanna Szoke
  Date:         13/09/00
  Function:     This program reads M (in range 2 - 50) sets of binary files 
		representing electron/pixels.  It calculates & writes 4 files;

		* the average at each point		(named average.bin)
		* the standard error or sqrt(variance)	(named sterror.bin)
		* error-weighted map          		(named erwm.bin)
		* an "unreliability" map		(named unrely.bin)
		
                Input parameters are read locally. UTIL.C has worker functions.
		High resolution is enabled.

*******************************************************************************/

#include "util.h"
#include "dims.h"

#define MAXNFILES	50	/* max # of files for stability calculations */
#define MINNFILES	 2	/* min # of files for stability calculations */

static	real	*xnump ;
static	real	*average ;
static	real	*sterror ;
static	real	*erwm ;
static	real	*unrely ;

static	char	out_fn[MAXSTRING] ;
static  int	dNhr[MAXNFILES] ;
static	real	*hrNump[MAXNFILES] ;
static	real	*hrNx[MAXNFILES] ;
static	real	*hrNy[MAXNFILES] ;
static	real	*hrNz[MAXNFILES] ;

static	void	add_to_average(real *, int) ;  /* accumulate sum (for average)  */
static	void	add_to_sterror(real *, int) ;  /* accumulate sumsq (for sterror) */
static	void	allocate_var_files() ;	/* ... pointers for Np-space files */
static	void	finish_average(int) ;	/* 1/M, write average	*/
static	void	finish_sterror(int) ;	/* compute, write sterror	*/
static	void	finish_erwm() ;
static	void	finish_unrely(real) ;

void variance_main(int argc, char *argv[])
{
	/* functions called, in alphabetical order */

     int	hrGetLength(char *) ;
     void	hrReadLists(char *, real *) ;		
     void	hrRenameList(real **, real **, real **) ;
     void	hrCompareIndices(real *, real *, real *) ;
     char	list_filename[MAXSTRING] ;
     void	readEpixFile(char *, real *) ;

     char	filename[MAXNFILES][MAXSTRING] ;
     int	M ;			/* number of Np file sets to handle */
     int	m ;			/* counters of files */
     real	pert_frac;

  /***********************************************************
   Check for right ballpark, say hello, identify file names.
  ***********************************************************/

     M = argc - optind - 3 ;

     if (M > MAXNFILES) 
	  EdenError("Too many files!") ;

     else if (M >= MINNFILES)
          pert_frac = atof(argv[optind+2]) ;
          for (m = 0; m < M; m++) 
               sprintf(filename[m], "%s", argv[optind+m+3]) ;

     if (M == 0) {
	  prompt("\nPlease enter pert_frac and 2 - 8 file names: ") ;
          sscanf(terminp, "%f %s %s %s %s %s %s %s %s", &t1,
		 filename[0], filename[1], filename[2], filename[3], 
		 filename[4], filename[5], filename[6], filename[7]) ;

          M = m = 0 ;
	  pert_frac = t1 ;
          while ((int)strlen(filename[m]) > 0)   {
	       M++; 
	       m++;
          }
          if (M < MINNFILES) 
	       EdenError("Too few files!") ;

	  sprintf(message, "Running %s with pert_frac %f on %d files", 
		  caller, pert_frac, M) ; 
	  printTwice(message) ;
     }

     hello(caller) ;

	/************************************************************
	Fetch input conditions, set up arrays in physical space.
	Read each binary el/voxel file, add to an average array.
	************************************************************/

     readBasicInput() ;

     for (m = 0; m < M; m++) {  

          if (HighRes) {
               sprintf(list_filename, "%s.list", filename[m]) ;
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
     Npextended = Nptotal + dNhr[0] ;

     allocate_var_files() ;

     sprintf(message, "\nNumber of files = %d, Pert_frac = %g\n", M, pert_frac);
     printTwice(message) ;

     printTwice("Calculating average among input files.\n") ;

     for (m = 0; m < M;m++)  {
          readEpixFile(filename[m], xnump) ;
          add_to_average(xnump, m) ;
     }
     finish_average(M) ;

	/************************************************************
	Compute and write standard error; this requires that we 
	reread the original files.
	************************************************************/

     printTwice("Calculating standard error among input files.\n") ;

     for (m = 0; m < M; m++)  {
          readEpixFile(filename[m], xnump) ;
          add_to_sterror(xnump, m) ;
     }
     finish_sterror(M) ;

	/************************************************************
	Compute and write the error-weighted map and "unreliability".
	************************************************************/

     printTwice(
     "Calculating the error-weighted average and unreliabiility maps.\n") ;

     finish_erwm() ;
     finish_unrely(pert_frac) ;

     return ;
}

void	allocate_var_files()
{
     int	n ;
     real	*pave, *pste ;

     xnump = (real *) e_malloc(Nptotal*sizeof(real), "variance") ;

     average = (real *) e_malloc(Npextended*sizeof(real), "variance") ;
     sterror = (real *) e_malloc(Npextended*sizeof(real), "variance") ;
     erwm = (real *) e_malloc(Npextended*sizeof(real), "variance") ;
     unrely = (real *) e_malloc(Npextended*sizeof(real), "variance") ;

     for (n = 0, pave = average; n <Npextended; n++, pave++)
          *pave = 0.;

     for (n = 0, pste = sterror; n <Npextended; n++, pste++)
          *pste = 0. ;

     return ;
}
void	add_to_average(real *nump, int m)		/* accumulate nump */  
{
     int 	n, i ;

     real	*pnump, *pave ;

     for (n = 0, pave = average, pnump = nump; n <Nptotal; n++, pave++,pnump++)
          *pave += *pnump ;

     if (HighRes)
          for (n = Nptotal, i = 0; n <Npextended; n++, pave++, i++)
               *pave += *(hrNump[m]+i) ;
     return;

}

void	finish_average(int M)
{
     void       writeEpixFile(char *, real *) ;
     void       hrWriteLists(char *, real *) ;

     int 	n ;
     real	*pave ;
 
     for (n = 0, pave = average; n <Npextended; n++, pave++)
          *pave /= M ;

     sprintf(out_fn, "average") ;
     writeEpixFile(out_fn, average) ;

     if (HighRes) 
          hrWriteLists("average.list", average+Nptotal) ;
     
     return;
}

void	add_to_sterror(real *nump, int m)		/* accumulate nump * nump */  
{
     int 	n, i ;
     real	*pnump, *pste ;

     for (n = 0, pste = sterror, pnump = nump; n <Nptotal; n++, pste++,pnump++)
          *pste += *pnump * *pnump ;

     if (HighRes)
          for (n = Nptotal, i = 0; n <Npextended; n++, pste++, i++)
               *pste += *(hrNump[m]+i) * *(hrNump[m] + i) ;

     return;
}

void	finish_sterror(int M)
{
     void       writeEpixFile() ;
     void       hrWriteLists(char *, real *) ;

     int 	n ;
     real	*pave, *pste ;
 
     for (n = 0, pste = sterror, pave = average; n <Npextended; 
          n++, pste++, pave++) {

          *pste -= M * *pave * *pave ;
          *pste /= (M-1) ;
          *pste = sqrt(fabs(*pste)) ;
     }

     sprintf(out_fn, "sterror") ;
     writeEpixFile(out_fn, sterror) ;

     if (HighRes) 
          hrWriteLists("sterror.list", sterror+Nptotal) ;
     
     return;
}


void	finish_erwm()
{
     void       writeEpixFile() ;
     void       hrWriteLists(char *, real *) ;

     int	n ;
     real	*pave, *pste, *per ;

     for (n = 0, pste = sterror, pave = average, per = erwm; n <Npextended; 
          n++, pste++, pave++, per++) {

          if (*pste > 1.e-5)
              *per = (*pave * *pave) / (*pave + *pste) ;
          else
              *per = *pave ;
     }
     sprintf(out_fn, "erwm") ;
     writeEpixFile(out_fn, erwm) ;

     if (HighRes) 
          hrWriteLists("erwm.list", erwm+Nptotal) ;
     
     return;
}

void	finish_unrely(real pert_frac)
{
     void       writeEpixFile() ;
     void       hrWriteLists(char *, real *) ;

     int	n ;
     real	*pave, *pste, *pun ;

     for (n = 0, pste = sterror, pave = average, pun = unrely; n <Npextended; 
          n++, pste++, pave++, pun++) 

         if (*pste > 1.e-5)
              *pun = *pste / (pert_frac * *pave + *pste) ;
	 else
              *pun = 0 ;
						
     sprintf(out_fn, "unreliability") ;
     writeEpixFile(out_fn, unrely) ;

     if (HighRes) 
          hrWriteLists("unreliability.list", unrely+Nptotal) ;
     
     return;
}

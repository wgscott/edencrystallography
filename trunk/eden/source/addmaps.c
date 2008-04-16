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

                                ADDMAPS & MULTMAPS

  Title:        ADDMAPS - Utility to add or subtract 2 el/voxel files.
                MULTMAPS - Utility to multiply 2 el/voxel files.
  Author:       Hanna Szoke
  Date:         9/29/97
  Function:     This program takes as input two sets of electron/voxel binary
 		files and 
		addmaps: writes out their sum or difference w. coefficients.
		multmaps: writes out their product w. coefficients.

*******************************************************************************/

#include "util.h"
#include "cellparams.h"
#include "dims.h"

static  int	Nhr1 = 0, Nhr2 = 0 ;
static	real	*hrNump1, *hrNump2 ;
static	real	*hrNx1, *hrNx2 ;
static	real	*hrNy1, *hrNy2 ;
static	real	*hrNz1, *hrNz2 ;


void	addmaps_main(int argc, char *argv[])
{
     void	fetch_coefs(real *, real *) ;	/* coefficients */
    void	readEpixFile(char *, real *) ;   	/* Read input file */
     void	analyzEpixFile(real *) ;		/* Report */
     void	writeEpixFile(char *, real *) ;	/* Write output file */
     int	hrGetLength(char *) ;
     void	hrRenameList(real **, real **, real **) ;
     void	hrReadLists(char *, real *) ;		
     void	hrWriteLists(char *, real *) ;		
     void	hrCompareIndices(real *, real *, real *) ;

     char	filename1[MAXSTRING] ;
     char	filename2[MAXSTRING] ;
     char	filename3[MAXSTRING] ;
     char	list_filename[MAXSTRING] ;

     real	*nump1 ;	/* electron densities file # 1*/
     real	*nump2 ;	/* electron densities file # 2 */
     real	*outp ;		/* sum or difference for output */
     real	*pn1, *pn2, *po ;
     real	*addlist ;

     int        n ;
     real	C1, C2 ;	/* coefficients, defaulting to 1. */

  /********************************************************
  Check for right ballpark, say hello, identify file names.
  ********************************************************/

     if (argc > optind+5) 	/* too many arguments */
	  ballpark(caller) ; 
     else if (argc == optind+5) {
          sprintf(filename1, "%s", argv[optind+2]) ;
          sprintf(filename2, "%s", argv[optind+3]) ;
          sprintf(filename3, "%s", argv[optind+4]) ;
     }
     else {
	  prompt("\nPlease enter 3 file names - 2 input, 1 output: ") ;
          sscanf(terminp, "%s %s %s", filename1, filename2, filename3) ;
	  sprintf(message, "Running %s on files %s and %s, output to %s", 
	          caller, filename1, filename2, filename3) ;
	  printTwice(message) ;
     }

     hello(caller) ;

  /********************************************************
  Fetch input conditions.  Set up arrays in physical space.
  Read the binary files.  
  ********************************************************/

     readBasicInput() ;
     fetch_coefs(&C1, &C2) ;

     nump1 = (real *) e_malloc(Nptotal*sizeof(real), caller) ; 
     nump2 = (real *) e_malloc(Nptotal*sizeof(real), caller) ; 
     outp  = (real *) e_malloc(Nptotal*sizeof(real), caller) ; 

     printTwice("Reading solution maps ...") ;

     readEpixFile(filename1, nump1) ;

     if (HighRes) {
          sprintf(list_filename, "%s.list", filename1) ;
	  Nhr1 = hrGetLength(list_filename) ;
          hrNump1 = (real *) e_malloc(Nhr1*sizeof(real), caller) ;
          hrReadLists(list_filename, hrNump1) ;	
	  sprintf(message, "including %d high-resolution points.", Nhr1) ;
	  printTwice(message) ;
	  hrRenameList(&hrNx1, &hrNy1, &hrNz1) ;
     }
     readEpixFile(filename2, nump2) ;

     if (HighRes) {
          sprintf(list_filename, "%s.list", filename2) ;
          Nhr2 = hrGetLength(list_filename) ;	
	  if (Nhr2 != Nhr1)
	       EdenError(
               "Your high-res list files do NOT have the same length!") ;

          hrNump2 = (real *) e_malloc(Nhr1*sizeof(real), caller) ;
          hrReadLists(list_filename, hrNump2) ;	
	  sprintf(message, "including %d high-resolution points.", Nhr2) ;
	  printTwice(message) ;
	  hrCompareIndices(hrNx1, hrNy1, hrNz1) ;
	  hrRenameList(&hrNx2, &hrNy2, &hrNz2) ;
     }

  /********************************************************
  Compute sum with coefficients and write it.
  ********************************************************/

     for (n = 0, pn1=nump1, pn2=nump2, po=outp; n < Nptotal; 
	  n++,   pn1++,     pn2++,     po++)
         *po = C1 * *pn1 + C2 * *pn2 ;

     printTwice("Writing sum ...\n") ;

     analyzEpixFile(outp) ;   
     writeEpixFile(filename3, outp) ;

     if (HighRes) {

	  addlist = (real *) e_malloc(Nhr1*sizeof(real), caller) ;

	  for (n = 0; n < Nhr1; n++)
	       *(addlist+n) = C1 * *(hrNump1+n) + C2 * *(hrNump2+n) ;

          strip_suffix(filename3, ".bin", message) ;
          sprintf(list_filename, "%s.list", message) ;
          hrWriteLists(list_filename, addlist) ;
     }
}


void	multmaps_main(int argc, char *argv[])
{
     void	readEpixFile(char *, real *) ;   	/* Read input file */
     void	analyzEpixFile(real *) ;		/* Report */
     void	writeEpixFile(char *, real *) ;	/* Write output file */

     char	filename1[MAXSTRING] ;
     char	filename2[MAXSTRING] ;
     char	filename3[MAXSTRING] ;

     real	*nump1 ;	/* electron densities file # 1*/
     real	*nump2 ;	/* electron densities file # 2 */
     real	*outp ;		/* sum or difference for output */
     real	*pn1, *pn2, *po ;

     int        n ;

  /********************************************************
  Check for right ballpark, say hello, identify file names.
  ********************************************************/

     if (argc > optind+5) 	/* too many arguments */
	  ballpark(caller) ; 
     else if (argc == optind+5) {
          sprintf(filename1, "%s", argv[optind+2]) ;
          sprintf(filename2, "%s", argv[optind+3]) ;
          sprintf(filename3, "%s", argv[optind+4]) ;
     }
     else {
	  prompt("\nPlease enter 3 file names - 2 input, 1 output: ") ;
          sscanf(terminp, "%s %s %s", filename1, filename2, filename3) ;
	  sprintf(message, "Running %s on files %s and %s, output to %s", 
	          caller, filename1, filename2, filename3) ;
	  printTwice(message) ;
     }

     hello(caller) ;

  /********************************************************
  Fetch input conditions.  Set up arrays in physical space.
  Read the binary files.  
  ********************************************************/

     readBasicInput() ;

     nump1 = (real *) e_malloc(Nptotal*sizeof(real), caller) ; 
     nump2 = (real *) e_malloc(Nptotal*sizeof(real), caller) ; 
     outp  = (real *) e_malloc(Nptotal*sizeof(real), caller) ; 

     printTwice("Reading input maps ...") ;

     readEpixFile(filename1, nump1) ;
     readEpixFile(filename2, nump2) ;

  /********************************************************
  Compute product and write it.
  ********************************************************/

     for (n = 0, pn1=nump1, pn2=nump2, po=outp; n < Nptotal; 
	  n++,   pn1++,     pn2++,     po++)
         *po = *pn1 * *pn2 ;

     printTwice("Writing product ...\n") ;

     analyzEpixFile(outp) ;   
     writeEpixFile(filename3, outp) ;

}

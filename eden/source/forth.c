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

                                FORTH

  Title:        Forward transformation of electron/voxel files.    
  Author:       Hanna Szoke
  Date:         11/11/93
  Function:     This program uses a plain FFT to transform 1 or 2
		files of electrons/voxel to (hkl) space and multiplies
		the transforms by appropriate exponential factors.
		Input is a simple grid or a body-centered grid of given 
		grid spacing. If HighRes is in effect, the non-gridded 
		input list will be read from an ascii file whose name
		corresponds to the .bin file of gridded data.

		Since Solve no longer writes out a "newhkl" file, Forth
		should be invoked whenever the final .fcalc is needed
		(e.g., for XTALview).

		It has been established that forth followed by back
		without apodization produces essentially the identity
		transformation;  the reverse - back followed by forth
		- is not exactly the identity because of possible
		differences in structure factor resolution.

*******************************************************************************/

#include "eden.h"		

static	COMPLEX	*Rforth ;	/* array for structure factors */
static	real	*fnump ;	/* electron per voxel, as read in  */
static	real	*hrlist ;	/* high-res electrons */

static	void f_transform() ;

void	forth_main(int argc, char *argv[])
{
	/* functions packaged here */

     void	f_transform() ;		/* Do the actual forth (FFT's) */

	/* functions packaged elsewhere */

     void	readEpixFile(char *, real *) ;
     int	hrGetLength(char *) ;
     void	hrReadLists(char *, real *) ;		
     void	initfft() ;		/* set up the FFT arrays	*/
     void	initfft_aux() ;		/* set up auxiliary FFT arrays	*/
     void	prepare_unique_hklmask(char **) ;	/* in hklutil.c */
     void 	writefcalc(char *, COMPLEX *, char *) ;          

     char	in_filename[MAXSTRING] ;
     char	list_filename[MAXSTRING] ;
     char	*umask ;

  /***********************************************************
  No switches other than -h (handled in eden.c).
  ***********************************************************/

     if (argc > optind+3) 	/* too many arguments */
	  ballpark(caller) ; 
     else if (argc == optind+3)
          sprintf(in_filename, "%s", argv[optind+2]) ;
     else {
	  prompt("\nPlease enter the name of an electron/voxel file: ") ;
          sscanf(terminp, "%s", in_filename) ;
	  sprintf(message, "Running Forth on el/voxel file %s", in_filename) ;
	  printTwice(message) ;
     }

     hello(caller) ;

  /************************************************************
  Fetch input conditions, read the input electron/voxel file(s)                 
  ************************************************************/

     readBasicInput() ;

     strcpy(message, "forth_main") ;

     fnump = (real *) e_malloc(Nptotal*sizeof(real), message) ; 

     readEpixFile(in_filename, fnump) ;

     if (HighRes) {
          sprintf(list_filename, "%s.list", in_filename) ;
	  Nhr = hrGetLength(list_filename) ;
          hrlist = (real *) e_malloc(Nhr*sizeof(real), message) ; 
	  hrReadLists(list_filename, hrlist) ;
     }

			       
     /************************************************************
     Prepare for FFTs.
     ************************************************************/

     initfft() ;
     initfft_aux() ;

     Rforth = (COMPLEX *) e_malloc(Nhkl*sizeof(COMPLEX), message) ; 

  /************************************************************
  Transform and write out (to pwd) structure factors. 
  ************************************************************/

     f_transform() ;

     printTwice("Writing structure factors ...") ;

     strip_path(in_filename, pwd_name) ;
     extend_filename(pwd_name, sf_filename, "_forth.hkl") ;

     prepare_unique_hklmask(&umask) ; 
     writefcalc(sf_filename, Rforth, umask) ;
}

void	f_transform()          

/* Transform, apply exponential Gaussian factors.   */

{

     real     *setupexpfac(real) ;
     COMPLEX    *setupexpfac_bcc(real) ;
     void	convolve(real *, real *, COMPLEX *, COMPLEX *) ;
     void	hrPrepareSFT() ;
     void	hrDoSFT(real *, COMPLEX *) ;

  /************************************************************
  Prepare exp factors for FFTs using delfac = PISQ*eta*dr*dr.
  ************************************************************/

     expfac = setupexpfac(delfac) ;
     if (grid_type == BODY_CENTERED)
	  expfac_bcc = setupexpfac_bcc(delfac) ;

     convolve(fnump, expfac, expfac_bcc, Rforth) ;

     if (HighRes) {
	  hrPrepareSFT() ;
	  hrDoSFT(hrlist, Rforth) ; 
     }
}


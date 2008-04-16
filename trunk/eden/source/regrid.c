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

                                REGRID 

  Title:        Transform Solve output to sampled electron densities;
		Write a .map file
  Author:       Hanna Szoke
  Date:         10/14/92
  Function:     This program takes as input a solution from Solve or Back
		(# electrons/ voxel, representing Gaussian peaks) and produces 
		a sampled electron density map in el/cubA, at a resolution that 
		is N times greater, where N is an integer (by default, 2).  
		N (spread_fac) is optionally read from the execute line:
		
			eden regrid param_file vox_file [N]

                Program expects to find a parameter file (param_file.inp) and
		a binary file named vox_file.bin. 

                If HighRes is in effect, the non-gridded input list will be 
                read from an ascii file named vox_file.list.

		You may specify in the .inp file that the range in 
		(x,y,z) is to span any fractional part of the unit cell, 
		including a fraction that exceeds 1.  The range (but not 
		necessarily the end points) must be positive.  Alternately,
		you may specify a pdb file name in the input; in this case,
		Regrid determines the output range based on the pdb information.

		Regrid writes vox_file_[N].map.

*******************************************************************************/

#include "util.h"
#include "cellparams.h"
#include "dims.h"
#include "symmetry.h"
#include "reg.h"

/*	variables in reg.h */

FLIMITS	flimx, flimy, flimz ;	/* input fractional xyz limits	*/
FLIMITS	flime ;
LIMITS	rlimx, rlimy, rlimz ;	/* coord limits after regrid */
char	vox_filename[MAXSTRING] ;
char	list_filename[MAXSTRING] ;
char	pdb_filename[MAXSTRING] ;
real	*rho ;		/* regridded electron densities */
real	*rnump ;		/* electron per voxel array read in */
real	*dout ;			/* rho's for output  */
int	spread_fac ;
char	erange ;		/* T/F ? - apply user's range (flime) */
int     Nxout, Nyout, Nzout ;   /* output space dimensions */
int     Npout ;                 /* product of Nxout*Nyout*Nzout */


void	regrid_main(int argc, char *argv[])
{
	/* functions packaged in regutil.c */

     void	pp_readInput() ;	/* Input common to Regrid & Count */
     void	rs_readInput() ;	/* Input specific to Regrid */
     void	redefine_dimensions() ;	/* Internal and output dimensions */
     void	spread(real **) ;	/* Distribute input over larger grid */
     void	transform(real *) ;	/* Do the actual regridding (FFT's) */
     void	prepareOut(real *) ;
     void	analyzeRange() ;	/* Find range for output density map */
     void	printXmap();		/* Write output in XPLOR format */

     void	readEpixFile(char *, real *) ; 
     void	hrRegInput(char *) ;	/* in regutil.c */
     char	strip_name[MAXSTRING] ;

     spread_fac = 2 ;
     
     if (argc > optind+4) 	/* too many arguments */
	  ballpark(caller) ; 
     else if (argc > optind+2) {
          sprintf(vox_filename, "%s", argv[optind+2]) ;
          if (argc == optind+4) 
	       spread_fac = atoi(argv[optind+3]) ;
     } 
     else {
	  prompt("\nPlease enter filename & optionally, spread factor: ") ;
          sscanf(terminp, "%s %s", vox_filename, message) ;
	  if ((int)strlen(message) > 0)
	       spread_fac = atoi(message) ;
	  sprintf(message, "Running %s on file %s with spread factor %d", 
		  caller, vox_filename, spread_fac) ;
	  printTwice(message) ;
     }

     hello(caller) ;

  /************************************************************
  Fetch input conditions, get input data
  ************************************************************/

     pp_readInput() ;
     rs_readInput() ;

     rnump = (real *) e_malloc(Nptotal*sizeof(real), caller) ;
     printTwice("Reading solution map ...") ;

     readEpixFile(vox_filename, rnump) ;

     if (HighRes) 
	  hrRegInput(vox_filename) ;

  /************************************************************
  Spread the gridded data over a finer grid.         
  Transform (FFT followed by FFT-1)       
  ************************************************************/

     redefine_dimensions() ;

     printTwice("Regridding ... ") ;

     spread(&rho) ;
     transform(rho) ;

  /************************************************************
  Output will be written in pwd rather than in directory from
  which the Solve/Back output was taken.
  ************************************************************/

     strip_path(vox_filename, pwd_name) ;
     strip_suffix(pwd_name, ".bin", strip_name) ;

     prepareOut(rho) ; 
     analyzeRange() ;

  /************************************************************
  Prepare filename for electron density map (output).   
  Name includes _N where N is the display resolution.
  ************************************************************/

     printTwice("Writing regridded solution file in XPLOR format.\n") ;
     sprintf(out_filename, "%s_%1d.map", strip_name, spread_fac) ;
     printXmap(dout) ;
}

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

                                TOHU

  Title:        TOHU - Eden preprocessor. 
  Author:       Hanna Szoke
  Date:         8/28/92
  Function:     This program converts .pdb input into a file of calculated
		structure factors.  It is thus a simplified subset of XPLOR.
		It uses a file of parameters and symmetry information plus
		the input .pdb.

                Output consists of an .fcalc file. 

*******************************************************************************/

#include "math.h"
#include "util.h"
#include "pdb.h"
#include "cellparams.h"
#include "symmetry.h"
#include "dims.h"
#include "calc_rho.h"

static	void	initStructureFactor(char *) ;

void	tohu_main(int argc, char *argv[])
{
     char	pdbName[MAXSTRING] ;
     void	fetchPdbInfo(char *) ;
     void	countPdbInfo() ;
     void	expandPdbInfo() ;
     void	estimate_F000() ;

     
  /* Check for right ballpark */

     if (argc > optind+3) 	/* too many arguments */
	  ballpark(caller) ; 
     else if (argc == optind+3)
          sprintf(pdbName, "%s", argv[optind+2]) ;
     else {
	  prompt("\nPlease enter the name of a pdb file: ") ;
          sscanf(terminp, "%s", pdbName) ;
	  sprintf(message, "Running %s on pdb file %s", caller, pdbName) ;
	  printTwice(message) ;
     }

     hello(caller) ;     

     printTwice
      ("\t**********************************************************") ;
     printTwice
      ("\tWARNING: Tohu is a slow substitute for other programs that");
     printTwice
      ("\tcalculate structure factors from PDB information.  It uses") ;
     printTwice
      ("\tB values, but assumes point atoms.  However, it produces  ") ;
     printTwice
      ("\tstructure factors on an absolute scale without any further") ;
     printTwice
      ("\tmanipulation.") ;
     printTwice
      ("\t**********************************************************\n") ;

     readBasicInput() ;

     fetchPdbInfo(pdbName) ;
     countPdbInfo() ;
     expandPdbInfo() ;

     if (!anom_flag) 
          estimate_F000() ;

     if (quick)
	     return ;

     initStructureFactor(pdbName) ;

     return ;
}

void	initStructureFactor(char *name)          

				/* Make structure factors from known atoms */
{
     void	prepare_unique_hklmask(char **) ;
     void	writefcalc_non0(char *, COMPLEX *, char *) ;

     char	output_filename[MAXSTRING] ;
     real     nusqval ;            /* absolute value of K/2pi squared */
     real	fx, fy, fz ;
     real	arg ;
     real	coef ;
     COMPLEX	*R, *bR ;
     char	*umask ;
     char	*unique_hkl ;
     int        j, n, nq ;
     int        h, k, l ;

     /**************************************************
     Initialize structure factors table from 
     CCP4's atomsf.lib.
     ***************************************************/
     fill_ff_tables();

     /**************************************************
     Convert pdb info from orthogonal to fractional 
     coordinates.
     **************************************************/

     for (j = 0; j < Natoms; j++) {
          ortho_to_fract(pos[j].x, pos[j].y, pos[j].z, &fx, &fy, &fz) ;
          pos[j].x = fx ;
          pos[j].y = fy ;
          pos[j].z = fz ;
     }

     /**************************************************
     If the pdb info is for anomalous data (heavy atoms 
     of one species only), we force the atoms to be 
     hydrogens;  the true Z (fj[j]) will be applied 
     with f' and f'' in Solve! 
     **************************************************/

     if (anom_flag) {
          printTwice("Replacing heavy atoms by hydrogens...") ;
          for (j = 0; j < Natoms; j++) 
	       fj[j] = (real) occup[j] ;
     }

     /**************************************************
     Prepare array.
     **************************************************/
     
     R = (COMPLEX *) e_malloc(Nhkl*sizeof(COMPLEX), "tohu") ;

     prepare_unique_hklmask(&unique_hkl) ;

     /* Give user a notion of the time Tohu will take */

     for (n = nq = 0, umask = unique_hkl; n < Nhkl; n++, umask++) 
          nq += *umask ;

     printTwice("Generating unique reflections only.") ;

     /**************************************************
     Now loop over (h,k,l) space and, for each 
     appropriate term, form the structure factor terms 
     (SLOW).
     **************************************************/
    
     for (n = nq = 0, umask = unique_hkl, bR = R; n < Nhkl; 
	  n++, umask++, bR++) {

          fetch_hkl(n, &h, &k, &l) ;
	  nusqval = nusq(h, k, l) ;

          if (*umask) {
		
               nq++ ;
               bR->re = 0. ;
               bR->im = 0. ;

               for (j = 0; j < Natoms; j++) {

		    arg = TWOPI * 
                          (h * pos[j].x + k * pos[j].y + l * pos[j].z) ; 
		    /**************************************************************
                     Calculate scattering factors using 5 parameter approach instead
                     of point atoms.

                     We need to get the sqrt(nusqval) because that that's
		     what the function expects. We also use fj[j]/occup[j] because
		     the function needs Z not Z*occup which is what fj looks to be.
		    ****************************************************************/
		    coef = scatt_f(sqrt(nusqval),(int)((fj[j]/occup[j])+0.5),Bj[j])*occup[j];
		    /*
		      Old point atom code
		      
		      coef = fj[j] * exp(-Bj[j] * nusqval / 4.) ;
		    */

                    bR->re += coef * cos(arg) ;
                    bR->im += coef * sin(arg) ;
               }

	       if (((nq) % 10000) == 0) {
	            sprintf(message,
			 "Finished calculating %d structure factors", nq) ;
                    printTwice(message) ;
               }
          }
     }

     /**************************************************
     Prepare .hkl file for output.
     Output will be written in pwd rather than in 
     directory from which the pdb file was read. 
     **************************************************/

     strip_path(name, pwd_name) ;
     strip_suffix(pwd_name, ".pdb", output_filename) ;
     strcat(output_filename, ".fcalc") ;

     sprintf(message, "Writing %s file\n", output_filename) ;
     printTwice(message) ;

     writefcalc_non0(output_filename, R, unique_hkl) ;

     return ;
}



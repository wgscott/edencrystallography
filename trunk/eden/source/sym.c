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

                                SYM

  Title:        Find symmetry-related points in unit cell.               
  Author:       Hanna Szoke
  Date:         11/11/93
  Function:     This program finds -- for a given symmetry group
		and unit cell data --
		the points that are equivalent to input (x,y,z).  

		In interactive mode, the program sits in an infinite 
		loop until user types "quit".   

		In non-interactive mode, the program reads a pdb file 
		and generates one pdb file.

		Option -c[overlap] (non-interactive mode) allows sym to report 
		and eliminate any atom if it overlaps another atom (in another 
		asymmetric unit or in another monomer).  Eliminated atoms 
		are taken out of all equivalent positions.

  Revisions:	4/20/95: added ncs.
                8/31/95: added option -c.
		12/18/96: report frac. limits in pdb file (non-interactive)
		8/2/01:	 withdrew ncs.

*******************************************************************************/

#include "util.h"
#include "symmetry.h"
#include "cellparams.h"
#include "dims.h"
#include "pdb.h"

static	int	check = FALSE ;
static	real	overlap = 0 ;
static	char	pdbName[MAXSTRING] ;
static	int     Ntot ;		/* total atoms in unit cell */
static	int	Nncs ;	

static	real	*fracx, *fracy, *fracz ;
		/* fract coord version of pdb coords  */

static	void    do_sym_interac() ;
static	void    do_sym_pdb() ;
static	void	s_readInput() ;
static	void	checkOverlaps() ;
static	int	findCsOverlap(int, int, int)  ;


void	sym_main(int argc, char *argv[])
{
     if (argc > optind+3) 	/* too many arguments */
	  ballpark(caller) ; 
     else 
	  if (!interactive) {
               if (argc == optind+3)
                    sprintf(pdbName, "%s", argv[optind+2]) ;
               else {
	            prompt("\nPlease enter the name of a pdb file: ") ;
                    sscanf(terminp, "%s", pdbName) ;
	            sprintf(message, 
		    "Running Sym on el/voxel file %s", pdbName) ;
	            printTwice(message) ;
               }
          }
     hello(caller) ;

     s_readInput() ;

     if (interactive) 
         do_sym_interac() ;
     else
         do_sym_pdb() ; 
}
void s_readInput()    /*   Read Sym input	*/ 
{
     int	k ;

     readBasicInput() ;
     Nncs = 1 ;  

     for (k = 0; k < Nlines; k++) {
          strcpy(nextline, allinp + k*MAXSTRING) ;
          sscanf(nextline, "%s", id) ;

          if (strcmp(id, "OVERLAP") == 0) {
	       sscanf(nextline, "%*s %g", &t1) ;
	       check = TRUE ;
	       overlap = t1 ;
          }
     }

     if (check) {
	  sprintf(message, "Overlap = %g", overlap) ;
          printTwice(message) ;
     }
}

	/***********************************************************
	For each input (x, y, z) entered at the user's terminal, 
	this function reports all symmetry-related (cs) 
	points in the unit cell. It loops until user types 'q'. 
	***********************************************************/
void	do_sym_interac()
{
     int        n, i ;
     real	xin, yin, zin ;
     real	newx, newy, newz ;
     real	fx, fy, fz ;
     real	newfx, newfy, newfz ;
     real	*mat ;

    if (Nncs == 1)
         sprintf(message, "Enter an orthogonal set (x, y, z) in Angstrom: ") ;
    else
         sprintf(message, 
"Enter an orthogonal set (x, y, z) in Angstrom\n within the PRIMARY monomer: ") ;

    prompt(message) ;

    while (terminp[0] != 'q') {

         xin = yin = zin = 0 ;

         sscanf(terminp, "%g %g %g", &t1, &t2, &t3) ;
	 xin = t1 ;
	 yin = t2 ;
	 zin = t3 ;
         fprintf(fp_log, "\nYour Input: %s\n", terminp) ;
         if (verbose) {
	        ortho_to_fract(xin, yin, zin, &fx, &fy, &fz) ;
	        sprintf(message, "... in fractional coords: %g %g %g", 
			fx, fy, fz) ;
	        printTwice(message) ;
         }

         for (n = 0; n < Nncs; n++) {

             ortho_to_fract(xin, yin, zin, &fx, &fy, &fz) ;

             for (i = 0, mat = matop; i < Ncs; i++, mat += MEL) {

                 apply_symop(fx, fy, fz, mat, &newfx, &newfy, &newfz) ;
                 fract_to_ortho(newfx, newfy, newfz, &newx, &newy, &newz) ;

                 sprintf(message, "\tx = %8.4f, y = %8.4f, z = %8.4f\n", 
                                  newx, newy, newz) ;
	         printTwice(message) ;
             }
         }
         prompt("\nEnter next (x, y, z)  or 'q': ") ;
    }
}

	/***********************************************************
	This function reads a pdb file and generates one pdb file 
	containing the coordinates for that monomer.
	***********************************************************/
void	do_sym_pdb()
{
	/* functions packaged in 'pdbutil.c'  */

     void	fetchPdbInfo() ;		/* get fields of pdb file */
     void	countPdbInfo() ;		/* count and report Z's */
     void	expandPartialPdbInfo() ;	/* expand 1 monomer */
     void	writePdbInfo() ;		/* write file for 1 monomer */
     void	estimate_F000() ;
     void	getRange(int, real *, real *, real *, real *, 
			 real *, real *) ;

     int        n ;
     char       label[MAXSTRING] ;
     real	fminx, fmaxx, fminy, fmaxy, fminz, fmaxz ;

   /*******************************************************************
   Storage order is: 
   asu 1 mono 1, asu 2 mono 1, ..., asu [Ncs] mono 1,     (1*Natoms)
   ...  
   ********************************************************************/

     fetchPdbInfo(pdbName) ;
     countPdbInfo() ;

     Natoms = Nin*Ncs ;
     Ntot = Nncs*Natoms ;

     if (Ntot > MAXNATOMS) 
          EdenError( 
             "No room in pdb arrays for all atoms in unit cell!") ;

	/***********************************************************
	Report range of the pdb info (useful for regrid limits).
	***********************************************************/

     printTwice("\nRange of data (fractional coords), based on input .pdb file:") ;
     getRange(Nin, &fminx, &fmaxx, &fminy, &fmaxy, &fminz, &fmaxz) ;

     sprintf(message, "x (%g, %g), y (%g, %g), z (%g, %g)",
	 fminx, fmaxx, fminy, fmaxy, fminz, fmaxz) ;
     printTwice(message) ;

	/* Initially, mark all atoms as non-overlapping */

     for (n = 0; n < Ntot; n++)
          marker[n] = 0 ;

	/* expand first monomer to other asu's (horizontal)  */

     expandPartialPdbInfo(0, Nin, Nin) ;

	/***********************************************************
	Report the full range of the pdb info and if user wants, 
	check the various nc's for overlap. 
	***********************************************************/

     printTwice
     ("\nRange of data (fractional coords) after CS expansion:") ;

     getRange(Ntot, &fminx, &fmaxx, &fminy, &fmaxy, &fminz, &fmaxz) ;

     sprintf(message, "x (%g, %g), y (%g, %g), z (%g, %g)",
	 fminx, fmaxx, fminy, fmaxy, fminz, fmaxz) ;
     printTwice(message) ;

     if (check) 
	  checkOverlaps() ;

     estimate_F000() ;


     /**************************************************
     Now write the monomers to the pdb files 
     Output will be written in pwd rather than in 
     directory from which the pdb file was read. 
     **************************************************/

     strip_path(pdbName, pwd_name) ;
     strip_suffix(pwd_name, ".pdb", message) ;

     for (n = 0; n < Nncs; n++) {
         sprintf(out_filename, "%s.n%d.pdb", message, n+1) ; 
         if (check)
              sprintf(label, "Monomer #%d, overlap distance = %g A", n+1,
			overlap) ;
         else
              sprintf(label, "Monomer #%d", n+1) ;
         writePdbInfo(out_filename, label, n*Natoms) ;
     }
}

	/***********************************************************
	This function determines whether separate asu's and monomers 
	overlap; insofar as they do, overlap points are removed from 
	consideration (marker is set). 
	***********************************************************/
void	checkOverlaps()
{

     real	fclip() ;			/* in util.c      */

     int        i, m, n ;
     int	Nover ;
     real	*fx, *fy, *fz ;			/* dummy pointers */
     real	onefx, onefy, onefz ;		/* single fract. coord. */

	/* prepare fractional coordinates */

     strcpy(message, "sym") ;

     fracx = (real *) e_malloc(Ntot*sizeof(real), message) ;
     fracy = (real *) e_malloc(Ntot*sizeof(real), message) ;
     fracz = (real *) e_malloc(Ntot*sizeof(real), message) ;

         for (i = 0, fx = fracx, fy = fracy, fz = fracz; i < Ntot; 
	      i++, fx++, fy++, fz++) {

                 ortho_to_fract(pos[i].x, pos[i].y, pos[i].z, 
                                &onefx, &onefy, &onefz) ;
                 *fx = fclip(onefx, 1.) ;
                 *fy = fclip(onefy, 1.) ;
                 *fz = fclip(onefz, 1.) ;
         }

         for (i = 0, Nover = 0; i < Nncs; i++)
              for (m = 0; m < Ncs-1; m++) 
                   for (n = m+1; n < Ncs; n++) 
                        Nover += findCsOverlap(m, n, i*Natoms) ;

         if (Nover > 0) {
              sprintf(message,
                "Total # crystal symmetry overlap points = %d", Nover) ;
              printTwice(message) ;
         }

}

   /*******************************************************************
	Check each atom of asu m against each atom of asu n
	(m < n) and eliminate the atoms in m (arbitrary choice) that 
	correspond to points of overlap.  This is done separately for 
	each monomer (value of off).  Storage order is: 

   asu 1 mono 1, asu 2 mono 1, ..., asu [Ncs] mono 1,  (Natoms)
   asu 1 mono 2, asu 2 mono 2, ..., asu [Ncs] mono 2   (2*Natoms) etc. 
   ********************************************************************/

int	findCsOverlap(int m, int n, int off) 

/*	m ;		asymmetric unit (zero origin) to check from */
/*	n ;		asymmetric unit (zero origin) to check to  */
/*	off ;		offset into array to start check */
{
     real	Rsq() ;		/* (crysutil.c) distance in fract coords */
     
     int	i, ii ;		/* id of point in asymmetric unit m */
     int	j ;		/* id of point in asymmetric unit n */
     int	mi, nj ;	/* offsets into arrays of fract coords */
     real	dfx, dfy, dfz ;	/* differences between fract. coords */
     int	Nover = 0 ;	/* counter of overlap points */


     for (i = 0, mi = m*Nin; i < Nin; i++, mi++) 

          if (marker[i] == 0) {

               for (j = 0, nj = n*Nin; j < Nin; j++, nj++) 

                    if (marker[j] == 0) {

                         dfx=fabs (*(fracx + off + mi) - *(fracx + off + nj)) ; 
                         dfy=fabs (*(fracy + off + mi) - *(fracy + off + nj)) ;
                         dfz=fabs (*(fracz + off + mi) - *(fracz + off + nj)) ;

                         if (dfx > 0.5)
                              dfx = 1. - dfx ;
                         if (dfy > 0.5)
                              dfy = 1. - dfy ;
                         if (dfz > 0.5)
                              dfz = 1. - dfz ;

                         if (Rsq(dfx, dfy, dfz) < overlap*overlap)  {
                             sprintf(message,
     "monomer %d: overlap of asym unit %d, # %d and asym unit %d, # %d", 
     off/Natoms + 1, m+1, atno[i], n+1, atno[j]) ;
			     printTwice(message) ;
                             
			     for (ii = i%Nin; ii < Ntot; ii+=Nin)
                                  marker[ii] = 1 ;

                             Nover ++ ;
			     break ;	/* no need to check any more j"s */
                         }
                    }
         }
     return (Nover) ;
}

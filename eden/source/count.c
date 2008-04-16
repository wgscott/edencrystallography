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

                                COUNT

  Title:        Transform Solve output to sampled electron densities;
		count electrons 
  Author:       Hanna Szoke
  Date:         10/14/92
  Function:     This program takes as input a solution from Solve or Back
		(# electrons/ voxel, representing Gaussian peaks) and produces 
		a sampled electron density map in el/cubA, at a resolution that 
		is N times greater, where N is an integer (by default, 2).  
		N (spread_fac) is optionally read from the execute line:
		
			eden count param_file vox_file [N]

                Program expects to find a parameter file (param_file.inp) and
		a binary file named vox_file.bin. 

                If HighRes is in effect, the non-gridded input list will be 
                read from an ascii file named vox_file.list.

		Count expects to find a pdb file named in param_file.inp, 
		It counts electrons in the sampled electron density map within 
		spheres of given radii around each pdb entry, writing the output
		to a file vox_file_{N}.count.

*******************************************************************************/

#include "util.h"
#include "cellparams.h"
#include "dims.h"
#include "pdb.h"
#include "symmetry.h"
#include "reg.h"

        /*      Input and derived parametersi from pdb.h	*/   

int	Natoms ;		/* number of atoms in unit cell	*/
int     Nin ;                	/* number of atoms in input */
int	atno[MAXNATOMS] ;	/* atom number */
char	atid[MAXNATOMS][5] ;	/* atom id */
char	groupid[MAXNATOMS][4] ;	/* amino acid id */
int	groupno[MAXNATOMS] ;	/* residue number */
POINT	pos[MAXNATOMS] ;	/* fractional positions of atoms in unit cell */
float	occup[MAXNATOMS] ;	/* occupancies */
real	Bj[MAXNATOMS] ;         /* Debye-Waller factors */
real	fj[MAXNATOMS] ;         /* number of electrons per input atom */
int	marker[MAXNATOMS] ;	/* marker of overlap atoms */

static	real	lev0, lev1, lev2 ;	/* levels for counting electrons */
static	real	Bcorr ;			/* B factor from apodization */
static	real	*fullx, *fully, *fullz ;	/* arrays used in count */
static	int	*fullp ;

static	void	c_readInput() ;		/* Read input specific to Count. */
static	void	count_electrons() ;	
static	void	print_el_count(char *, real *, real *, real *) ;
static	int	replicate_pdb() ;


void	count_main(int argc, char *argv[])
{
	/* functions called */

     void	pp_readInput() ;
     void	redefine_dimensions() ;	/* Internal and output dimensions */ 
     void	spread(real **) ;	/* Distribute input over larger grid */
     void	transform(real *) ;	/* Do the actual regridding (FFT's) */
     void	hrRegInput(char *) ;	/* input of list  files, in regutil.c */

     void	readEpixFile(char *, real *) ; 
     void	countPdbInfo() ;
     void	expandPdbInfo() ;

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
  Fetch input conditions, get input data.
  ************************************************************/

     pp_readInput() ;
     c_readInput() ;

     rnump = (real *) e_malloc(Nptotal*sizeof(real), caller) ;
     printTwice("Reading solution map ...") ;

     readEpixFile(vox_filename, rnump) ;

     if (HighRes) 
          hrRegInput(vox_filename) ;

     redefine_dimensions() ;

     countPdbInfo() ;
     expandPdbInfo() ;

  /************************************************************
  Spread the solutions over a finer grid.         
  Transform (FFT followed by FFT-1)       
  ************************************************************/

     printTwice("Regridding ... ") ;

     spread(&rho) ;
     transform(rho) ;

  /************************************************************
  Use it to count electrons. 
  Output will be written in pwd rather than in directory from
  which the Solve/Back output was taken.
  ************************************************************/

     strip_path(vox_filename, pwd_name) ;
     strip_suffix(pwd_name, ".bin", strip_name) ;

     sprintf(out_filename, "%s_%1d.count", strip_name, spread_fac) ;
     count_electrons() ;
}

void c_readInput()    /*   Read input specific to Count	*/ 
{
     int	k ;

     /* Defaults for optional and some required input.   */

     lev0 = 1. ;
     lev1 = 1.5 ;
     lev2 = 2. ;
     Bcorr = UNSET ;

     for (k = 0; k < Nlines; k++) {
          strcpy(nextline, allinp + k*MAXSTRING) ;
          sscanf(nextline, "%s", id) ;

          if (strcmp(id, "LEVELS") == 0)  {
	       sscanf(nextline, "%*s %g %g %g", &t0, &t1, &t2) ;
	       lev0 = t0 ;
	       lev1 = t1 ;
	       lev2 = t2 ;
          }
          else if (strcmp(id, "BCORR") == 0) {
               sscanf(nextline, "%*s %g", &t0) ;
               Bcorr = t0 ;
          }
     }

     /************************************************************
     Check for legality of input
     ************************************************************/

     if (Bcorr == UNSET) 
	 EdenError("Missing B correction from apodization!") ;
     if (strcmp(pdb_filename, "none") == 0)
	 EdenError("Missing pdb file name!") ;
     else{
	 sprintf(message,
 "Electrons will be counted to %g, %g and %g fractions of input res.",
	 lev0, lev1, lev2) ;
         printTwice(message) ;
	 sprintf(message, "B correction from apodization is %g.", Bcorr) ;
         printTwice(message) ;
     }
}

        /************************************************************
  	Use data in full regridded array, rho, to count electrons 
	around each pdb atom positions.
	Three radii are used: by default, 1, 1.5 and 2 * input_res.  
	Other choices of fractional limits may be requested using 
	keyword LEVELS in the input.
        ************************************************************/

void	count_electrons()
{
     real	thisnump ;
     real	*count0, *count1, *count2 ;
     real	*holdr ;
     real	*parti ;
     real	*coef ;			/* coefficients for partitioning */
     int	*checkp ;
     int	*holdp ;
     real	x, y, z, x0 ;
     real	rsq ;
     real	rsq0, rsq1, rsq2 ;
     real	Bfac ;
     real	norm0, norm1, norm2 ;
     real	fri, frj, frk ;    
     real	part ;
     int	Nrep_pdb ;
     int 	n, p, h ;
     int	hmax ;
     int	Nused = 0 ;

     printTwice("Counting electrons ...") ;

     /****************************************************
     Move all pdb points into the unit cell.   Replicate 
     those that are close to cell edges.  Return Nrep_pdb
     (Natoms+replicated ones) plus for each, the index p
     and the positions (x, y, z).
     ****************************************************/

     Nrep_pdb = replicate_pdb() ;

     /****************************************************
     Allocate space for identifying characteristics of pdb 
     entries that contribute to each grid point.
     Here, hmax represents an estimate of the total number
     of  pdb's among which a grid point may have to be
     partitioned (surely, <= Natoms).

     Allocate, initialize space for counting electrons.
     ****************************************************/

     sprintf(message, "counting_electrons") ;

     hmax = Natoms ;

     checkp = (int *) e_malloc(Nrep_pdb*sizeof(int), message) ;
     holdr  = (real *) e_malloc(hmax*sizeof(real), message) ;
     parti  = (real *) e_malloc(hmax*sizeof(real), message) ;
     holdp  = (int *) e_malloc(hmax*sizeof(int), message) ;
     coef   = (real *) e_malloc(Natoms*sizeof(real), message) ;

     for (p = 0; p < Natoms; p++) {
          Bfac = Bj[p] + Bcorr ;
          *(coef+p)  = fj[p] * occup[p] *
		       pow((real) (2. * TWOPI) / Bfac, (real) 1.5) ;
     }

     count0 = (real *) e_malloc(3*Natoms*sizeof(real), message) ;
     count1 = count0 + Natoms ;
     count2 = count0 + 2*Natoms ;

     for (p = 0; p < 3*Natoms; p++) 
          *(count0 + p) = 0 ; 

     /****************************************************
     lev0, lev1, and lev2 are radii in units of input_res 
     whose values default to 1, 1.5 and 2. respectively.
     ****************************************************/

     rsq0 = lev0 * lev0 * input_res * input_res ;
     rsq1 = lev1 * lev1 * input_res * input_res ;
     rsq2 = lev2 * lev2 * input_res * input_res ;

     /****************************************************
     Now for the loop over grid points.
     ****************************************************/

     sprintf(message, "There will be %d calculations\n", Np) ;
     printTwice(message) ;

     for (n = 0; n < Np; n++) {

        thisnump = *(rho + n) ;

        /************************************************************
  	Determine, for each position in the unit cell, which pdb 
	atom(s) extend over it and by what fraction (coverage).
	We use a clever scheme to cut down on unnecessary churning
	that uses the fact that the x- and a- axes coincide.
        ************************************************************/

	index_to_fract(n, Nx, Ny, Nz, &fri, &frj, &frk) ;

        if (n%Nx == 0) {

	     fract_to_ortho(fri, frj, frk, &x0, &y, &z) ;

	     for (p = 0; p < Nrep_pdb; p++) {

                  rsq = (
		         (*(fully+p) - y) * (*(fully+p) - y) +
		         (*(fullz+p) - z) * (*(fullz+p) - z)) ;

		  *(checkp+p) = (rsq < rsq2) ? 1 : 0 ;
             }
        }

	x = x0 + fri * aAxis ;
     
	for (p = 0, h = 0; p < Nrep_pdb; p++) {

             if (*(checkp+p) == 1) {

                  rsq = ((*(fullx+p) - x) * (*(fullx+p) - x) +
	                 (*(fully+p) - y) * (*(fully+p) - y) +
	                 (*(fullz+p) - z) * (*(fullz+p) - z)) ;

                  if (rsq < rsq2) {

                       Bfac = Bj[*(fullp+p)] + Bcorr ;
		       *(holdp+h) = *(fullp+p) ;
		       *(holdr+h) = rsq ;
		       *(parti+h) = *(coef+*(fullp+p)) * 
			       exp(-TWOPI*TWOPI*rsq / Bfac) ;
		       h++ ;
                  }
             }
        }
	hmax = h ;		
	
	if (hmax > 0)
	     Nused++ ;

	/****************************************************
	Now hmax is the # of pdb's overlapping grid point n.
	Normalize partition contributions, apply to counts.
	****************************************************/

	for (h = 0, norm0 = norm1 = norm2 = 0; h < hmax; h++) {
	       rsq = *(holdr+h) ;
	       if (rsq < rsq2)
                    norm2 += *(parti+h) ;
	       if (rsq < rsq1)
                    norm1 += *(parti+h) ;
	       if (rsq < rsq0)
                    norm0 += *(parti+h) ;
        }

	for (h = 0; h < hmax; h++) {

	     rsq = *(holdr+h) ;
             p = *(holdp+h) ;
	     part = *(parti+h) ;

	     if ((rsq < rsq2) && (norm2 > 0))
                  *(count2 + p) += thisnump * part / norm2 ;
             if ((rsq < rsq1) && (norm1 > 0))		  
                  *(count1 + p) += thisnump * part / norm1 ;
             if ((rsq < rsq0) && (norm0 > 0))  
                  *(count0 + p) += thisnump * part / norm0 ;
        }

	if ((n > 0) && (n%100000) == 0) {
	     sprintf(message, "Finished calc for n = %d", n) ;
	     printTwice(message) ;
        }
     }

     sprintf(message, 
     "\nCoverage information will be written to %s", out_filename) ;
     printTwice(message) ;

     print_el_count(out_filename, count0, count1, count2) ;

     sprintf(message, "Counting used %d out of %d grid points.", Nused, Np) ;
     printTwice(message) ;
}

        /************************************************************
  	Write out results of counting electrons. 
        ************************************************************/

void	print_el_count(char *filename, 
          real *count0, real *count1, real *count2)
{
     FILE 	*fp ;
     real	totnump ;
     real	remainder0 = 0, remainder1 = 0, remainder2 = 0 ;
     real	volvox ;
     int 	n, p ;

     for (n = 0, totnump = 0; n < Np; n++) 
         totnump += *(rho + n) ;

     remainder0 = totnump ;
     remainder1 = totnump ;
     remainder2 = totnump ;

     for (p = 0; p < Natoms; p++) {
	  remainder0 -= *(count0 + p) ;
	  remainder1 -= *(count1 + p) ;
	  remainder2 -= *(count2 + p) ;
     }
     /****************************************************
     Convert from el/cuA to total electrons, report.
     Note: volvox is NOT the same as cubA_to_vox in the
     rest of Eden because of the spread_fac in Np!
     ****************************************************/

     volvox = crys_vol / Np ;

     if ((fp = fopen(filename, "w")) == NULL)
          EdenError("Couldn't open ,count file.") ;

     fprintf(fp, 
    "\n\tPDB data and electron counts within %3.1f, %3.1f, %3.1f res radii\n\n",
	     lev0, lev1, lev2) ;
     fprintf(fp,
"Ind  At no   At id      x       y      z     occ   B    %5.1f   %5.1f   %5.1f\n\n",
	     lev0, lev1, lev2) ;

     for (p = 0; p < Natoms; p++) {
          fprintf(fp,
          "%4d %5d %4s %3s %7.3f %7.3f %7.3f %5.2f %5.2f %7.3f %7.3f %7.3f\n",
          p, atno[p], atid[p], groupid[p], 
	  pos[p].x, pos[p].y, pos[p].z, occup[p], Bj[p],
	  *(count0+p)*volvox, *(count1+p)*volvox, *(count2+p)*volvox) ;
     }

     sprintf(message, "\nTotal no. of electrons is %g\n", totnump*volvox) ;
     fprintf(fp, message) ;
     printTwice(message) ;

     sprintf(message, "Count of all remaining electrons: %g %g %g\n", 
     remainder0*volvox, remainder1*volvox, remainder2*volvox) ;
     fprintf(fp, message) ;
     printTwice(message) ;

     fclose(fp) ;
}

        /************************************************************
        Ensure that all pdb points are in unit cell and replicate all 
	pdb entries that overlap faces, edges and corners.
	At most, there are 8 such replications.
	In verbose mode, write replication info to log.
        ************************************************************/

int	replicate_pdb()
{
     int 	next, n, p ;
     int	rep_p = 0 ;
     real	fri, frj, frk ;    
     real	newfri, newfrj, newfrk ;    
     real	marginx, marginy, marginz ;
     real	nextfri[8], nextfrj[8], nextfrk[8] ;

     /****************************************************
     Allocate space for identifying replicated pdb entries
     ****************************************************/

     sprintf(message, "replicate_pdb") ;

     fullp = (int *) e_malloc(8*Natoms*sizeof(int), message) ;
     fullx = (real *) e_malloc(8*Natoms*sizeof(real), message) ;
     fully = (real *) e_malloc(8*Natoms*sizeof(real), message) ;
     fullz = (real *) e_malloc(8*Natoms*sizeof(real), message) ;

     /****************************************************
     Prepare margins of error.  1.5 is magic number ...
     so is 0.5 in test of margins ...
     ****************************************************/

     marginx = 1.5 * lev2*input_res / aAxis ;
     marginy = 1.5 * lev2*input_res / bAxis ;
     marginz = 1.5 * lev2*input_res / cAxis ;

     if ((marginx >= 0.5) || (marginy >= 0.5) || (marginz >= 0.5)) 
          EdenError(
          "Level 2 is too large! - margins extend over full cell ...") ;

     for (p = 0; p < Natoms; p++) {

          /****************************************************
          First, ensure all pdb points are in unit cell.
          ****************************************************/

          ortho_to_fract(pos[p].x, pos[p].y, pos[p].z, &fri, &frj, &frk) ;

	  while (fri < 0)
	       fri += 1 ;
	  while (fri > 1)
	       fri -= 1 ;
	  while (frj < 0)
	       frj += 1 ;
	  while (frj > 1)
	       frj -= 1 ;
	  while (frk < 0)
	       frk += 1 ;
	  while (frk > 1)
	       frk -= 1 ;

	  fract_to_ortho(fri, frj, frk, fullx+rep_p, fully+rep_p, fullz+rep_p) ;

          *(fullp+rep_p) = p ;
	  rep_p++ ;

          /****************************************************
          Then, look for pdb's that are close to faces, edges 
          and corners, requiring, respectively, 1, 3 or 7 new 
          entries.
          ****************************************************/

          if ((fri < marginx) || (fri > 1. - marginx) ||
              (frj < marginy) || (frj > 1. - marginy) ||
              (frk < marginz) || (frk > 1. - marginz)) {

              if (verbose) {
                   sprintf(message, 
	           "p = %d, old xyz = %g %g %g", *(fullp+rep_p-1), 
		    *(fullx+rep_p-1), *(fully+rep_p-1), *(fullz+rep_p-1)) ;
                   printTwice(message) ;
              }
	      newfri= fri ;
	      newfrj= frj ;
	      newfrk= frk ;
	      next = 0 ;

	      for (n = 0; n < 8; n++) {
		   nextfri[n] = fri ;
		   nextfrj[n] = frj ;
		   nextfrk[n] = frk ;
              }

	      /************************************
	      Account for faces.
	      ************************************/

              if (fri < marginx) {
		   nextfri[next] = newfri = fri + 1 ;
		   next++ ;
              }
              else if (fri > 1. - marginx) {
		   nextfri[next] = newfri = fri - 1 ;
		   next++ ;
              }
              if (frj < marginy) {
		   nextfrj[next] = newfrj = frj + 1 ;
		   next++ ;
              }
              else if (frj > 1. - marginy) {
		   nextfrj[next] = newfrj = frj - 1 ;
		   next++ ;
              }
              if (frk < marginz) {
		   nextfrk[next] = newfrk = frk + 1 ;
		   next++ ;
              }
              else if (frk > 1. - marginz) {
		   nextfrk[next] = newfrk = frk - 1 ;
		   next++ ;
              }

	      /************************************
	      Account for edges.
	      ************************************/

              if ((newfri != fri) && (newfrj != frj)) {
		   nextfri[next] = newfri ;
		   nextfrj[next] = newfrj ;
		   next++ ;
              }

              if ((newfri != fri) && (newfrk != frk)) {
		   nextfri[next] = newfri ;
		   nextfrk[next] = newfrk ;
		   next++ ;
              }

              if ((newfrj != frj) && (newfrk != frk)) {
		   nextfrj[next] = newfrj ;
		   nextfrk[next] = newfrk ;
		   next++ ;
              }

	      /************************************
	      Finally, account for corners.
	      ************************************/

              if ((newfri != fri) && (newfrj != frj) && (newfrk != frk)) {
		   nextfri[next] = newfri ;
		   nextfrj[next] = newfrj ;
		   nextfrk[next] = newfrk ;
		   next++ ;
              }
	      for (n = 0; n < next; n++) {

	           fract_to_ortho(nextfri[n], nextfrj[n], nextfrk[n], 
				  fullx+rep_p, fully+rep_p, fullz+rep_p) ;
	           *(fullp+rep_p) = p ;

		   if (verbose) {

                       sprintf(message, 
		       "p = %d, new xyz = %g %g %g", *(fullp+rep_p),
		       *(fullx+rep_p), *(fully+rep_p), *(fullz+rep_p)) ;
		       printTwice(message) ;

                   }
	           rep_p++ ;
             }
         }
     }
     sprintf(message, "Number of pdb entries expanded from %d to %d.",
     Natoms, rep_p) ;
     printTwice(message) ;

     return (rep_p) ;
}

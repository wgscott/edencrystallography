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

                                SHAPES 

  Title:        Analyze Solve output in terms of the morphology at each point.
  Author:       Hanna Szoke
  Date:         10/14/92
  Function:     This program takes as input a solution from Solve or Back
		(# electrons/ voxel, representing Gaussian peaks) and produces 
		an array of numbers in range (-1, 9) representing the shape
		of the 3-dim. volume at each point in the requested region.  
		See analyze_shapes() for the interpretation of these numbers.
		N (spread_fac) is optionally read from the execute line:
		
			eden shapes param_file vox_file [N]

		but is 2 by default.
                Program expects to find a parameter file (param_file.inp) and
		a binary file, vox_file.bin. 

                If HighRes is in effect, the non-gridded input list will be 
                read from an ascii file named vox_file.list.

		You may specify in the .inp file that the range in 
		(x,y,z) is to span any fractional part of the unit cell, 
		including a fraction that exceeds 1.  The range (but not 
		necessarily the end points) must be positive.  Alternately,
		you may specify a pdb file name in the input; in this case,
		Shapes determines the output range based on the pdb range.

		Shapes writes vox_file_[N].map.  It also writes files 
		shapes_[f].map containing sshape indices - one per value of 
		f, the density cutoff fraction that defaults to 0.1.

*******************************************************************************/

#include "util.h"
#include "cellparams.h"
#include "dims.h"
#include "symmetry.h"
#include "reg.h"

#define	SMALL		1.e-1		/* for round-off of shapes A, B, C */
#define	RHO_SMALL	1.e-1		/* for round-off of rho */

static	void	sh_readInput() ;
static	void	getAscale() ;
static	void	analyze_shapes() ;
static	void	renorm_derivs() ;
static	void	averageRho() ;
static	void    analyze_one_shape(int, real *, real *, real *) ;
static	real	*dummy ;	/* for holding spread() contents, 2nd time */

real	*drdx ;			/* 3 1st derivatives... */
real	*drdy ;
real	*drdz ;
real	*d2rdxdx ;		/* and 6 2nd derivatives. */
real	*d2rdxdy ;
real	*d2rdxdz ;
real	*d2rdydy ;			
real	*d2rdydz ;
real	*d2rdzdz ;

static	int	Nrlist ;
static	real	*rho_list ;
static	real	rho_min ;
static	real	rho_cutoff ;
static	real	*shape_values ;
static	real	coef1, coef2, coef3, coef4 ;
static 	real	data_res ;
static  real	Ascale1, Ascale2, Ascale3 ;


void	shapes_main(int argc, char *argv[])
{
	/* functions packaged in regutil.c */

     void	pp_readInput() ;	/* Input common to all postprocessors */
     void	rs_readInput() ;	/* Input for Regrid and Shapes */
     void	redefine_dimensions() ;	/* Internal and output dimensions */
     void	allocate_deriv_arrays() ;
     void	spread(real **) ;	/* Distribute input over larger grid */
     void	transform(real *) ;	/* Do the actual regridding (FFT's) */
     void	shapes_transform(real *) ;	/* Derivatives via FFT's */
     void	printXmap();		/* Write output in XPLOR format */
     void	prepareOut(real *) ;	/* Transfer data to output array */

     void	readEpixFile(char *, real *) ; 
     void	hrRegInput(char *) ;		

     int	k ;
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
     sh_readInput() ;

     rnump = (real *) e_malloc(Nptotal*sizeof(real), caller) ;
     printTwice("Reading solution map ...") ;

     readEpixFile(vox_filename, rnump) ;

     if (HighRes) 
	  hrRegInput(vox_filename) ;

  /************************************************************
  Spread the gridded data over a finer grid.         
  Transform (FFT followed by FFT-1)       
  Prepare 1st and 2nd derivatives.
  ************************************************************/

     redefine_dimensions() ;
     allocate_deriv_arrays() ;

     printTwice("Regridding ... ") ;

     spread(&rho) ;
     transform(rho) ;

     spread(&dummy) ;
     shapes_transform(dummy) ;

  /************************************************************
  Output will be written in pwd rather than in directory from
  which the Solve/Back output was taken.
  ************************************************************/

     strip_path(vox_filename, pwd_name) ;
     strip_suffix(pwd_name, ".bin", strip_name) ;

     prepareOut(rho) ; 
     renorm_derivs() ;

  /************************************************************
  Prepare filename for electron density map, write it out.   
  ************************************************************/

     printTwice(
  "Writing regridded solution file and shapes file[s] in XPLOR format.\n") ;
     sprintf(out_filename, "%s_%1d.map", strip_name, spread_fac) ;
     printXmap(dout) ;

  /************************************************************
  At this point the 3 1st and 6 2nd derivative arrays are 
  ready to be analyzed.
  ************************************************************/

     for (k = 0; k < Nrlist; k++) {
	  rho_cutoff = *(rho_list + k) ;
	  sprintf(message, "\nUsing rho_cutoff, fraction of max(density) = %g", 
	       rho_cutoff) ;
	  printTwice(message) ;
          averageRho() ;
          analyze_shapes() ;
	  prepareOut(shape_values) ;
     
  /************************************************************
  The desired part of shape_values has been transferred to dout.
  Prepare filenames for shapes map, write out.
  ************************************************************/

          sprintf(out_filename, "shapes_%g.map", rho_cutoff) ;
          printXmap(dout) ;
     }
}

#define	Nrho	10 

  /************************************************************
  Read levels for rho_cutoff.
  Read "data_res" (defaults to input_res)
  ************************************************************/
void	sh_readInput()
{
     int	k ;
     float	f[Nrho] ;
     float	rlmin = 1. ;
     float	rlmax = 0. ;

     t1 = UNSET ;

     for (k = 0; k < Nrho; k++)
	  f[k] = 0 ;

     for (k = 0; k < Nlines; k++) {
          strcpy(nextline, allinp + k*MAXSTRING) ;
          sscanf(nextline, "%s", id) ;

          if (strcmp(id, "RHO_CUTOFF") == 0)
             sscanf(nextline, "%*s %g %g %g %g %g %g %g %g %g %g", 
			      f, f+1, f+2, f+3, f+4, f+5, f+6, f+7, f+8, f+9) ;
          else if (strcmp(id, "DATA_RES") == 0)
	     sscanf(nextline, "%*s %g", &t1) ;
     }

     data_res = (t1 != UNSET) ? t1 : input_res ;

     if (f[0] == -1) {		/*default = no user input */
	  Nrlist = 1 ;
	  rho_list = (real *) e_malloc(Nrlist*sizeof(real), caller) ;
	  *rho_list = RHO_SMALL ;
     }
     else {
	  Nrlist = 0 ;
	  k = 0 ;
	  while (f[k] >= 0) {
	       k++ ;
	       Nrlist++ ;
          }
	  rho_list = (real *) e_malloc(Nrlist*sizeof(real), caller) ;
	  for (k = 0; k < Nrlist; k++) {
	       *(rho_list+k) = f[k];

               if ((*(rho_list+k) > 1.) || (*(rho_list+k) < 0)) {
	            sprintf(message, 
		         "Rho_list value, %g, is out of range (0, 1) !", 
	                 *(rho_list+k)) ;
                    printTwice(message) ;
	            EdenError("Quitting.") ;
               }
	       else {
		    if (f[k]< rlmin) rlmin = f[k] ;
		    if (f[k]> rlmax) rlmax = f[k] ;
               }
          }
     }

     /*  Report input value to output. */ 

     if (Nrlist == 1)
        sprintf(message, "Rho_cutoff fraction is %g", *rho_list) ;
     else
        sprintf(message, "There are %d Rho_cutoff fractions in range %g - %g",
		  Nrlist, rlmin, rlmax) ;
     printTwice(message) ;
     sprintf(message, "Data resolution (used for scaling derivatives) is %g",
	     data_res) ;
     printTwice(message) ;
}

void	averageRho()	/* Average a solution map  */
{
     int	n, goodn = 0 ;
     real	*this_rho, rho0, rho_max ;

     for (n = 0, rho_max = 0, this_rho = rho; n < Np; n++, this_rho++)
	  if (rho_max < *this_rho)
	       rho_max = *this_rho ;

     rho_min = rho_cutoff * rho_max ;

     for (n = 0, rho0 = 0, this_rho = rho; n < Np; n++, this_rho++)
         if (*this_rho > rho_min) {
	     rho0 += *this_rho ;
             goodn++ ;
         }
     rho0 /= goodn ;

     sprintf(message, "Average density over occupied part of unit cell is %6.3f", 
	     rho0) ;
     printTwice(message) ;
     sprintf(message, "This includes %d out of %d points", goodn, Np) ;
     printTwice(message) ;
}

/*	Ascale1 is |min(Ashape)| / 3.	*/

void	getAscale()
{
     real	Ashape, Bshape, Cshape ;
     real	local_rho ;
     int	n ;

     sprintf(message, "Shapes will ignore density values less than %6.3f", 
	     rho_min) ;  
     printTwice(message) ;

     for (n = 0, Ascale1 = 1000; n < Np; n++) {

          local_rho = *(rho+n) ;

          if (fabs(local_rho) >= rho_min) {
               analyze_one_shape(n, &Ashape, &Bshape, &Cshape) ;
	       if (Ashape < Ascale1)
		    Ascale1 = Ashape ;
          }
     }
     Ascale1 = fabs(Ascale1) / 3. ;
     Ascale2 = Ascale1 * Ascale1 ;
     Ascale3 = Ascale2 * Ascale1 ;

     if (verbose) {
	  sprintf(message, 
   "Shapes will be rescaled for epsilon determination, using factors:") ;
	  printTwice(message) ;
          sprintf(message,"%g (A), %g (B) and %g (C)", Ascale1, Ascale2, Ascale3) ;
          printTwice(message) ;
     }
}

enum {OTHER = -1, UNIFORM, BUBBLE, TUBE, THROAT, NEGLENS,
                  SADDLE,  LENS,   NECK,  SNAKE,    BLOB} ;

void	analyze_shapes()
{
     real	Ashape, Bshape, Cshape ;
     real	Amin = 0, Bmin = 0, Cmin = 0 ;
     real	Amax = 0, Bmax = 0, Cmax = 0 ;
     real	local_rho, this_shape ;
     real	pc ;
     int	n ;
     int	x, y, z ;
     int	uniform = 0, bubble = 0,  tube = 0 ;
     int	throat = 0,  neglens = 0, saddle = 0 ;
     int 	lens = 0,    neck = 0,    snake = 0 ;
     int	blob = 0 ,   other = 0 ;

     shape_values = (real *) e_malloc(Np*sizeof(real), caller) ;
     getAscale() ;

     for (n = 0; n < Np; n++) {

          local_rho = *(rho+n) ;

          if (fabs(local_rho) < rho_min) {
               this_shape = UNIFORM ;
               uniform++ ;
          }
	  else {
               analyze_one_shape(n, &Ashape, &Bshape, &Cshape) ;

/****************************************************************************
			Here are the assignments:
		C = 0			C > 0			C < 0          
	-------------------	-------------------	-------------------
	A<0	A=0	A>0	A<0	A=0	A>0	A<0	A=0	A>0
B<0	 5	 5	 5	 7	 7	 7	 3	 3	 3
B=0	 6	 0	 4	 7	no	no	no	no	 3
B>0	 8	no	 2	 7	no	 1	 9	no	 3
*****************************************************************************/ 

               if (fabs(Ashape) * Ascale1 < SMALL)
                    Ashape = 0 ;
               if (fabs(Bshape) * Ascale2 < SMALL)
                    Bshape = 0 ;
               if (fabs(Cshape) * Ascale3 < SMALL)
                    Cshape = 0 ;

               if ((Ashape == 0) && (Bshape == 0) && (Cshape == 0)) {
                    this_shape = UNIFORM ;
                    uniform++ ;
               }
               else if ((Ashape < 0) && (Bshape > 0) && (Cshape < 0)) {
                    this_shape = BLOB ;
                    blob++ ;
               }
               else if ((Ashape < 0) && (Bshape > 0) && (Cshape == 0)) {
                    this_shape = SNAKE ;
                    snake++ ;
               }
               else if ((Ashape < 0) && (Bshape >= 0) && (Cshape > 0)) {
                    this_shape = NECK ;
                    neck++ ;
               }
               else if ((Ashape < 0) && (Bshape == 0) && (Cshape == 0)) {
                    this_shape = LENS ;
                    lens++ ;
               }
               else if ((Ashape > 0) && (Bshape > 0) && (Cshape > 0)) {
                    this_shape = BUBBLE ;
                    bubble++ ;
               }
               else if ((Ashape > 0) && (Bshape == 0) && (Cshape == 0)) {
                    this_shape = NEGLENS ;
                    neglens++ ;
               }
               else if ((Ashape > 0) && (Bshape > 0) && (Cshape == 0)) {
                    this_shape = TUBE ;
                    tube++ ;
               }
               else if ((Ashape > 0) && (Bshape >= 0) && (Cshape < 0)) {
                    this_shape = THROAT ;
                    throat++ ;
               }
               else if (                (Bshape < 0) && (Cshape < 0)) {
                    this_shape = THROAT ;
                    throat++ ;
               }
               else if (                (Bshape < 0) && (Cshape == 0)) {
                    this_shape = SADDLE ;
                    saddle++ ;
               }
               else if (                (Bshape < 0) && (Cshape > 0)) {
                    this_shape = NECK ;
                    neck++ ;
               }
               else {
                    this_shape = OTHER ;
                    other++ ;
               }
               if (verbose) {
                    x = n % Nx ;
                    y = (n/Nx) % Ny ;
                    z = n / (Nx*Ny) ;
                    fprintf(fp_log, 
    "coefs: %6.3f %6.3f %6.3f %6.3f Type: %g x,y,z: %3d %3d %3d A,B,C: %6.3f %6.3f %6.3f rho: %6.3f\n",
     coef1, coef2, coef3, coef4, this_shape, x, y, z, 
     Ashape*Ascale1, Bshape*Ascale2, Cshape*Ascale3, local_rho) ; 
               }

               if (Amax < Ashape*Ascale1) Amax = Ashape*Ascale1 ;
               if (Bmax < Bshape*Ascale2) Bmax = Bshape*Ascale2 ;
               if (Cmax < Cshape*Ascale3) Cmax = Cshape*Ascale3 ;

               if (Amin > Ashape*Ascale1) Amin = Ashape*Ascale1 ;
               if (Bmin > Bshape*Ascale2) Bmin = Bshape*Ascale2 ;
               if (Cmin > Cshape*Ascale3) Cmin = Cshape*Ascale3 ;
          }
	  *(shape_values+n) = this_shape ;
     }
     pc = 100./Np ;
     printTwice("\nShape Distribution for full unit cell") ;
     printTwice("\nshape name (num)   %   # points\n") ;

    sprintf(message," blobs       (%1d) %4.1f %8d", BLOB, pc*blob, blob) ;
     printTwice(message) ;
    sprintf(message," snakes      (%1d) %4.1f %8d", SNAKE, pc*snake, snake) ;
     printTwice(message) ;
    sprintf(message," necks       (%1d) %4.1f %8d", NECK, pc*neck, neck) ;
     printTwice(message) ;
    sprintf(message," lenses      (%1d) %4.1f %8d", LENS, pc*lens, lens) ;
     printTwice(message) ;
    sprintf(message," saddle pts  (%1d) %4.1f %8d", SADDLE, pc*saddle, saddle) ;
     printTwice(message) ;
    sprintf(message," neg lenses  (%1d) %4.1f %8d", NEGLENS, pc*neglens, neglens);
     printTwice(message) ;
    sprintf(message," throats     (%1d) %4.1f %8d", THROAT, pc*throat, throat) ;
     printTwice(message) ;
    sprintf(message," tubes       (%1d) %4.1f %8d", TUBE, pc*tube, tube) ;
     printTwice(message) ;
    sprintf(message," bubbles     (%1d) %4.1f %8d", BUBBLE, pc*bubble, bubble) ;
     printTwice(message) ;
    sprintf(message," uniform pts (%1d) %4.1f %8d", UNIFORM, pc*uniform, uniform);
     printTwice(message) ;
    sprintf(message," Total                %8d", Np) ;
     printTwice(message) ;
     printTwice("\n") ;

     if (other > 0) {
         sprintf(message,
	 "Note: there are also %2d other (illegal) points with index %2d",
			other, OTHER) ;
         printTwice(message) ;
         printTwice("\n") ;
     } 

     if (verbose && (Ascale1 > 0)) {
          sprintf(message, 
          "Range of A unscaled: (%6.2f to %6.2f) and scaled: (%6.2f to %6.2f)", 
	     Amin, Amax, Amin/Ascale1, Amax/Ascale1) ;
          printTwice(message) ;

          sprintf(message, 
          "Range of B unscaled: (%6.2f to %6.2f) and scaled: (%6.2f to %6.2f)", 
	     Bmin, Bmax, Bmin/Ascale2, Bmax/Ascale2) ;
          printTwice(message) ;

          sprintf(message, 
          "Range of C unscaled: (%6.2f to %6.2f) and scaled: (%6.2f to %6.2f)", 
	     Cmin, Cmax, Cmin/Ascale3, Cmax/Ascale3) ;
          printTwice(message) ;
     }
}

void	renorm_derivs()
{
     real	local_rho ;
     real	Nco, Mco ;
     real	dr_fac ;
     int	n ;

     /*		Renormalize the first and 2nd derivatives;
		Use same criterion for examining them as will be used
		in analyze_shapes() for using them.			*/

     dr_fac = data_res*sqrt(eta / 2.) ;

     for (n = 0; n < Np; n++) {

          local_rho = *(rho+n) ;
          Nco = (fabs(local_rho) > rho_min) ? dr_fac/local_rho : 0 ;

	  *(drdx+n) *= Nco ;
	  *(drdy+n) *= Nco ;
	  *(drdz+n) *= Nco ;
     }

     for (n = 0; n < Np; n++) {

          local_rho = *(rho+n) ;
          Mco = (fabs(local_rho) > rho_min) ? (dr_fac*dr_fac)/local_rho : 0 ;

	  *(d2rdxdx+n) *= Mco ;
	  *(d2rdxdy+n) *= Mco ;
	  *(d2rdxdz+n) *= Mco ;
	  *(d2rdydy+n) *= Mco ;
	  *(d2rdydz+n) *= Mco ;
	  *(d2rdzdz+n) *= Mco ;
      }
}

void	analyze_one_shape(int n, real *Ashape, real *Bshape, real *Cshape)
{
     real	a11, a12, a13, a22, a23, a33 ;
     real	b11, b12, b13, b22, b23, b33 ;

     a11 = 1 + *(drdx+n) * *(drdx+n) ;
     a12 =     *(drdx+n) * *(drdy+n) ;
     a13 =     *(drdx+n) * *(drdz+n) ;
     a22 = 1 + *(drdy+n) * *(drdy+n) ;
     a23 =     *(drdy+n) * *(drdz+n) ;
     a33 = 1 + *(drdz+n) * *(drdz+n) ;

     b11 = *(d2rdxdx+n) ;
     b12 = *(d2rdxdy+n) ;
     b13 = *(d2rdxdz+n) ;
     b22 = *(d2rdydy+n) ;
     b23 = *(d2rdydz+n) ;
     b33 = *(d2rdzdz+n) ;

     coef1 = a11 * (a22*a33 - a23*a23) -
             a12 * (a12*a33 - a13*a23) +
             a13 * (a12*a23 - a13*a22) ;

     coef2 = b11 * (a22*a33 - a23*a23) -
             a12 * (b12*a33 - b13*a23) +
             a13 * (b12*a23 - b13*a22) +
             a11 * (b22*a33 - b23*a23) -
             b12 * (a12*a33 - a13*a23) +
             a13 * (a12*b23 - a13*b22) +
             a11 * (a22*b33 - a23*b23) -
             a12 * (a12*b33 - a13*b23) +
             b13 * (a12*a23 - a13*a22) ;

     coef3 = a11 * (b22*b33 - b23*b23) -
             b12 * (a12*b33 - a13*b23) +
             b13 * (a12*b23 - a13*b22) +
             b11 * (a22*b33 - a23*b23) -
             a12 * (b12*b33 - b13*b23) +
             b13 * (b12*a23 - b13*a22) +
             b11 * (b22*a33 - b23*a23) -
             b12 * (b12*a33 - b13*a23) +
             a13 * (b12*b23 - b13*b22) ;

     coef4 = b11 * (b22*b33 - b23*b23) -
             b12 * (b12*b33 - b13*b23) +
             b13 * (b12*b23 - b13*b22) ;

     *Ashape = coef2 / coef1 ;
     *Bshape = coef3 / coef1 ;
     *Cshape = coef4 / coef1 ;
}
void	allocate_deriv_arrays()
{
     drdx = (real *)e_malloc(Np*sizeof(real), caller) ;
     drdy = (real *)e_malloc(Np*sizeof(real), caller) ;
     drdz = (real *)e_malloc(Np*sizeof(real), caller) ;
     d2rdxdx = (real *)e_malloc(Np*sizeof(real), caller) ;
     d2rdxdy = (real *)e_malloc(Np*sizeof(real), caller) ;
     d2rdxdz = (real *)e_malloc(Np*sizeof(real), caller) ;
     d2rdydy = (real *)e_malloc(Np*sizeof(real), caller) ;
     d2rdydz = (real *)e_malloc(Np*sizeof(real), caller) ;
     d2rdzdz = (real *)e_malloc(Np*sizeof(real), caller) ;
}

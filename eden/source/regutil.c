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

                          UTILITIES FOR POSTPROCESSING

  Title:        Regutil - utilities for postprocessing.
  Author:       Hanna Szoke
  Date:         09/15/01
  Function:     This package contains functions used by 
		Regrid, Shapes and Count.

*******************************************************************************/
#include "util.h"
#include "cellparams.h"
#include "dims.h"
#include "pdb.h"
#include "symmetry.h"
#include "reg.h"

#define	TINY	1.e-4		/* for round-off of output */


static	real	*hrNump ;			/* for hr points */

static	int	Nxin, Nyin, Nzin ;
static	real	normfac ;
static	real	scale_factor ;

static	void	combine_octants(COMPLEX *, real) ;
static	void	sh_combine_octants1(COMPLEX *, real, char, char, char) ;
static	void	sh_combine_octants2(COMPLEX *, real, char, char, char, 
                                    char, char, char) ;
static	void	do_der1(COMPLEX *, COMPLEX *, COMPLEX *, COMPLEX *,
                        char, char, char, real *) ;
static	void	do_der2(COMPLEX *, COMPLEX *, COMPLEX *, COMPLEX *,
                        char, char, char, char, char, char, real *) ;

 	/************************************************************
	Spread elements of rnump onto a finer grid, fine, which will
	contain sampled electrons per cubic Angstrom.  Spread uses the
	given spread_fac.  This function includes assembling the 
	simple and intercalating grids. 
	************************************************************/

void	spread(real **fine)
{
     int 	i, j, k ;
     int 	ii, jj, kk ;

     real 	*prnump ;
     real	*pfine ;
     int	off = spread_fac / 2 ;

     *fine = (real *) e_malloc(Np*sizeof(real), caller) ;

     prnump = rnump ;
     pfine = *fine;

     for (k = 0; k < Np; k++)
	  *(pfine++) = 0 ;

     for (k = 0, pfine = *fine; k < Nzin; k++) {
	kk = spread_fac*k ;
        for (j = 0; j < Nyin; j++){
	   jj = spread_fac*j ;
           for (i = 0; i < Nxin; i++) {
	        ii = spread_fac * i ; 
                *(pfine + ii + jj*Nx + kk*Nx*Ny) = *(prnump++) ;
	   } 
        }
     }

     if (grid_type == BODY_CENTERED)
        for (k = 0, pfine = *fine; k < Nzin; k++) {
	   kk = spread_fac*k + off ;
           for (j = 0; j < Nyin; j++){
	      jj = spread_fac*j + off ;
              for (i = 0; i < Nxin; i++) {
	           ii = spread_fac * i + off ; 
                   *(pfine + ii + jj*Nx + kk*Nx*Ny) = *(prnump++) ;
	      } 
           }
        }

}

        /************************************************************
  	Transfer data from full regridded array, full, to dout, 
	using limits as defined by rlimx, rlimy and rlimz.  
	These limits are not necessarily a proper subset of the full 
	regridded array and may span the edges of the unit cell.   
        ************************************************************/

void	prepareOut(real *full) 
{
     int 	i, j, k ;
     int	n, ii, jj, kk ;
     LIMITS	ilims, jlims, klims ;

     Nxout = (rlimx.end - rlimx.beg) ;
     Nyout = (rlimy.end - rlimy.beg) ;
     Nzout = (rlimz.end - rlimz.beg) ;

     Npout = Nxout * Nyout * Nzout ;
     
     if (dout == NULL)
         dout = (real *) e_malloc(Npout*sizeof(real), message) ;

     ilims = rlimx ;
     jlims = rlimy ;
     klims = rlimz ;
     while (ilims.beg < 0) {
	  ilims.beg += Nx ;
	  ilims.end += Nx ;
     }
     while (jlims.beg < 0) {
	  jlims.beg += Ny ;
	  jlims.end += Ny ;
     }
     while (klims.beg < 0) {
	  klims.beg += Nz ;
	  klims.end += Nz ;
     }
     for (n = 0, k = klims.beg; k < klims.end; k++) {
	kk = k % Nz ;
        for (j = jlims.beg; j < jlims.end; j++){
	   jj = j % Ny ;
           for (i = ilims.beg; i < ilims.end; i++, n++) {
		ii = i % Nx ;
                *(dout+n) = *(full + ii + jj*Nx + kk*Nx*Ny) ;
	   } 
        }
     }
}

#include	<sys/time.h>	/* ... for picking up date & time	*/

	/************************************************************
	Write a file of output electron densities in XPLOR format.
	************************************************************/

void	printXmap(real *array)	
{
     FILE 	*fp ;
     void       get_min_max(real *, int, real *, real *) ; 
     int	z, n, i, j ;
     real	ar_min ;
     real	ar_max ;

     if ((fp = fopen(out_filename, "w")) == NULL)
          EdenError("Couldn't open .map file.") ;

     /* Ensure that no tiny values appear */

     for (n = 0; n < Npout; n++) 
    	  if (fabs(*(array+n)) < TINY)
	       *(array+n) = 0 ;

     get_min_max(array, Npout, &ar_min, &ar_max) ;

     if (ar_min == ar_max) 
	  EdenError("All the data are effectively equal!") ;

	/************************************************************
	May 4, 2004: There was a serious problem, now corrected, 
	in the header: Nxout -> Nx, etc. 
	************************************************************/

     fprintf(fp,"\n       3 !NTITLE\n") ; 
     fprintf(fp," REMARKS Electron Density Information from EDEN\n") ;
     fprintf(fp," REMARKS DATE: %s \n", timestamp()) ;
     fprintf(fp," REMARKS Density range is (%g, %g), scale factor is %g\n", 
             ar_min, ar_max, scale_factor) ;
     fprintf(fp," %7d %7d %7d %7d %7d %7d %7d %7d %7d\n",
	     Nx, rlimx.beg, rlimx.end-1,
	     Ny, rlimy.beg, rlimy.end-1, 
	     Nz, rlimz.beg, rlimz.end-1) ;
     fprintf(fp,"%12.5E%12.5E%12.5E%12.5E%12.5E%12.5E\n",
	     aAxis, bAxis, cAxis, angles[0], angles[1], angles[2]) ;
     fprintf(fp, "ZYX\n") ;

     z = rlimz.beg ;
     for (n = 0; n < Nzout; n++, z++) {
        fprintf(fp, "%8d\n", z) ;
        for (i = 0; i < Nxout*Nyout; i += 6) {
         for (j = 0; j < 6; j++)
	   if (i+j < Nxout*Nyout)
                fprintf(fp, "%12.5E", *(array+n*Nxout*Nyout+i+j)) ;
	 fprintf(fp, "\n") ;
        }
     }
     fprintf(fp, "   -9999\n") ;
     fclose(fp) ;
}

void	analyzeRange()	/* Analyze a solution map  */
{
     void       get_min_max(real *, int, real *, real *) ; 
     real	ar_min ;
     real	ar_max ;
     int	n ;
     char	messout[MAXSTRING] ;

     get_min_max(dout, Npout, &ar_min, &ar_max) ;

     if (erange) {
          scale_factor = (real) (flime.end - flime.beg) / (ar_max - ar_min) ;

          for (n = 0; n < Npout; n++) 
	       *(dout+n) = flime.beg + scale_factor * (*(dout+n) - ar_min) ;
     }

     sprintf(messout,
        "\nRange of output is (%g, %g) el/cubA", ar_min, ar_max) ;
     printTwice(messout) ;
}
	/************************************************************
	Transform, apply exponential Gaussian factors, and then 
	back-transform densities.  
	In contrast to the regular Solve transformations to/from 
	Fourier space, there is no need to expand the COMPLEX R to Nhkl.
	Thus, if there is high resolution, we redefine Nh, Nk and Nl 
	to be Nx, Ny and Nz, respectively so that we can use doSFT to
	accumulate the list contributions into the R array.
	************************************************************/

void transform(real *dnump)          
{
     void	initfft() ;
     void	do_fft(COMPLEX *, int) ;
     void    	hrResetIndices(int);
     void	hrDoregSFT(COMPLEX *, real *);	

     COMPLEX 	*R=NULL, *Rhigh=NULL ;
     COMPLEX	*pR, *pRhigh ;
     real	*pdnump ;
     int        n ;

     initfft() ;

      strcpy(message, "transform") ;
 
      R = (COMPLEX *) e_malloc(Np*sizeof(COMPLEX),  message) ;

     /* Feed the input into R */

     for (n=0, pdnump=dnump, pR=R; n < Np; n++, pR++, pdnump++) {
          pR->re = *pdnump ;
          pR->im = 0 ;
     }

     printTwice("Starting forward FFT...") ;

     do_fft(R, 1) ;			/* FORWARD FFT */

     if (HighRes) {
          hrResetIndices(spread_fac);

          Rhigh = (COMPLEX *) e_malloc(Np*sizeof(COMPLEX),message) ;

          for (n=0, pRhigh=Rhigh; n < Np; n++, pRhigh++) {
               pRhigh->re = 0 ;
               pRhigh->im = 0 ;
          }

          hrDoregSFT(Rhigh, hrNump) ;
     }

     /* Using data from full ellipsoid, apply the spread factor */

     combine_octants(R, delfac) ;

     if (HighRes) {
          combine_octants(Rhigh, delfac/4.) ;

          for (n=0, pR = R, pRhigh=Rhigh; n < Np; n++, pR++, pRhigh++) {
               pR->re += pRhigh->re ;
               pR->im += pRhigh->im ;
          }
     }

     do_fft(R, -1) ;		/* BACK FFT */

     printTwice("...Finished back FFT") ;

  /************************************************************
   Retrieve the real part and renormalize using 1/volume, 
   to convert to electrons per cubic Angstrom. See for example, 
   Crystal Structure Analysis, Glusker & Trueblood, p.204.
  ************************************************************/

     normfac = 1. / crys_vol ;

     for (n=0, pdnump=dnump, pR=R; n < Np; n++, pR++, pdnump++) 
          *pdnump = pR->re * normfac ;
}

	/***********************************************************
	Bring together the exponential function contributions from 
	the 8 octants of the full ellipsoid.  This differs from 
	sum_over_octants, in that R is NOT expanded to the half-
	ellipsoid here;  instead, the factors that would apply if
	it were so expanded are used to "combine" the data.
	***********************************************************/

void	combine_octants(COMPLEX *R, real dfac)
{
     real     one_expfac ;
     int	n, h, k, l ;

     for (n = 0; n < Np; n++) {

          h = n % Nx ;
	  k = (n/Nx) % Ny ;
	  l = n / (Nx*Ny) ;

          one_expfac = exp(-(dfac * nusq(h, k, l))) ;

	  if (k > 0) 
               one_expfac += exp(-(dfac * nusq(h, k-Ny, l))) ;
       
	  if (l > 0)
               one_expfac += exp(-(dfac * nusq(h, k, l-Nz))) ;
           
	  if ((k > 0) && (l > 0)) 
               one_expfac += exp(-(dfac * nusq(h, k-Ny, l-Nz))) ;
           
	  if (h > 0) 
               one_expfac += exp(-(dfac * nusq(h-Nx, k, l))) ;
        
	  if ((h > 0) && (k > 0)) 
               one_expfac += exp(-(dfac * nusq(h-Nx, k-Ny, l))) ;
         
	  if ((h > 0) && (l > 0))
               one_expfac += exp(-(dfac * nusq(h-Nx, k, l-Nz))) ;
          
	  if ((h > 0) && (k > 0) && (l > 0)) 
               one_expfac += exp(-(dfac * nusq(h-Nx, k-Ny, l-Nz))) ;
	   
          (R+n)->re *= one_expfac ;
          (R+n)->im *= one_expfac ;
      }
}

void redefine_dimensions()
{
     /* Redefine internal (regridded) dimensions */

     Nx = spread_fac * Nxin ;
     Ny = spread_fac * Nyin ;
     Nz = spread_fac * Nzin ;
     Np = Nx * Ny * Nz ;

     /* Translate user's limits into limits in the output array */

     rlimx.beg = (floor) (Nx * flimx.beg) ;
     rlimy.beg = (floor) (Ny * flimy.beg) ;
     rlimz.beg = (floor) (Nz * flimz.beg) ;

     rlimx.end = (ceil) (Nx * flimx.end) ;
     rlimy.end = (ceil) (Ny * flimy.end) ;
     rlimz.end = (ceil) (Nz * flimz.end) ;

     if ((flimx.beg != 0) || (flimx.end != 1) ||
         (flimy.beg != 0) || (flimy.end != 1) ||
         (flimz.beg != 0) || (flimz.end != 1)) {

          printTwice(
	    "\nThe actual gridded coordinate ranges in Angstrom follow -") ;
          sprintf(message,
            "\tx limits: (%g, %g)", aAxis*rlimx.beg / Nx, aAxis*rlimx.end / Nx);
          printTwice(message) ;
          sprintf(message,
            "\ty limits: (%g, %g)", bAxis*rlimy.beg / Ny, bAxis*rlimy.end / Ny);
          printTwice(message) ;
          sprintf(message,
            "\tz limits: (%g, %g)", cAxis*rlimz.beg / Nz, cAxis*rlimz.end / Nz);
          printTwice(message) ;

          printTwice("\n") ;
     }
}

void pp_readInput()    /*   Read input common to Regrid, Shapes & Count */ 
{
     void	fetchPdbInfo(char *) ;
     void	getRange(int, real *, real *, real *, real *, 
			 real *, real *) ;
     int	k ;
     int	status = 0 ;

     readBasicInput() ;

     Nxin = Nx ;
     Nyin = Ny ;
     Nzin = Nz ;

     /* Defaults for optional and some required input.   */

     flimx.beg = 0 ;	flimx.end = 1 ;
     flimy.beg = 0 ;	flimy.end = 1 ;
     flimz.beg = 0 ;	flimz.end = 1 ;

     strcpy(pdb_filename, "none") ;

     for (k = 0; k < Nlines; k++) {
          strcpy(nextline, allinp + k*MAXSTRING) ;
          sscanf(nextline, "%s", id) ;

          if (strcmp(id, "X_LIMITS") == 0) {
               sscanf(nextline, "%*s %g %g", &t1, &t2) ;
	       flimx.beg = t1 ; 
	       flimx.end = t2 ;
          }
          else if (strcmp(id, "Y_LIMITS") == 0) {
               sscanf(nextline, "%*s %g %g", &t1, &t2) ;
	       flimy.beg = t1 ; 
	       flimy.end = t2 ;
          }
          else if (strcmp(id, "Z_LIMITS") == 0) {
               sscanf(nextline, "%*s %g %g", &t1, &t2) ;
	       flimz.beg = t1 ; 
	       flimz.end = t2 ;
          }
          else if ((strcmp(id, "PDB_FILENAME") == 0) ||
                   (strcmp(id, "PDB_FN") == 0))
               sscanf(nextline, "%*s %s", pdb_filename) ;
     }

     /************************************************************
     Be sure spread_fac is reasonable & consistent with grid type 
     ************************************************************/

     if (spread_fac > 8)
          EdenError("Your spread factor is too large!") ;

     if ((grid_type == BODY_CENTERED) && ((spread_fac%2) != 0)) 
	  EdenError(
  "You can't have an odd regrid factor with a body-centered grid!") ;

     /************************************************************
     Check for legality of dimensional input
     ************************************************************/

     if ((flimx.end - flimx.beg) < 0)  {
	 EdenWarning("Negative x_limits span in input.") ;
         status = -1 ;
     }
     if ((flimy.end - flimy.beg) < 0)  {
	 EdenWarning("Negative y_limits span in input.") ;
         status = -1 ;
     }
     if ((flimz.end - flimz.beg) < 0)  {
	 EdenWarning("Negative z_limits span in input.") ;
         status = -1 ;
     }

     if (status == -1)
	  EdenError("Quitting.") ;

     /************************************************************
     Fetch Pdb input.
     ************************************************************/

     if (strcmp(pdb_filename, "none") != 0) {
          fetchPdbInfo(pdb_filename) ;
          getRange(Natoms, &flimx.beg, &flimx.end, &flimy.beg, &flimy.end,
                   &flimz.beg, &flimz.end) ;
     }

     /*  Report input values to output. */ 

     sprintf(message, 
          "Unit cell measures %6.2f by %6.2f by %6.2f Angstrom", 
	  aAxis, bAxis, cAxis);
     printTwice(message) ;
     sprintf(message, "Eta, exp. factor, is %g", eta) ;
     printTwice(message) ;
     sprintf(message, "Regridding factor (spread_fac) is %d", spread_fac) ;
     printTwice(message) ;

     printTwice("\nThe requested coordinate ranges in Angstrom follow -") ;
     sprintf(message, 
	       "\tx limits: (%g, %g)", aAxis*flimx.beg, aAxis*flimx.end) ;
     printTwice(message) ;
     sprintf(message, 
	       "\ty limits: (%g, %g)", bAxis*flimy.beg, bAxis*flimy.end) ;
     printTwice(message) ;
     sprintf(message, 
	       "\tz limits: (%g, %g)", cAxis*flimz.beg, cAxis*flimz.end) ;
     printTwice(message) ;
}

void rs_readInput()    /*   Read input for Regrid and Shapes only	*/ 
{
     int	k ;
     int	status = 0 ;

     /* Defaults for optional and some required input.   */

     flime.beg = 0 ;
     flime.end = 1 ;
     erange = FALSE ;

     for (k = 0; k < Nlines; k++) {
          strcpy(nextline, allinp + k*MAXSTRING) ;
          sscanf(nextline, "%s", id) ;

          if (strcmp(id, "E_LIMITS") == 0) {
               sscanf(nextline, "%*s %g %g", &t1, &t2) ;
	       flime.beg = t1 ; 
	       flime.end = t2 ;
          }
     }

     if ((flime.end - flime.beg) < 0)  {
	 EdenWarning("Negative electron limits span in input.") ;
         status = -1 ;
     }

     if (status == -1)
	  EdenError("Quitting.") ;

     if ((flime.beg != 0) || (flime.end != 1))
	  erange = TRUE ;

}

	/************************************************************
	Transform, apply exponential Gaussian factors, and then 
	back-transform densities for 1st and 2nd derivatives.   
	************************************************************/

void shapes_transform(real *dnump)          
{
     void	freefft() ;
     void	do_fft(COMPLEX *, int) ;
     void    	hrResetIndices(int);
     void	hrDoregSFT(COMPLEX *, real *);	

     COMPLEX 	*R=NULL, *Rhold=NULL ;
     COMPLEX 	*Rhigh=NULL, *Rhihold=NULL ;
     COMPLEX	*pR, *pRhigh, *pRhold, *pRhihold ;
     real	*pdnump ;
     int        n ;

     strcpy(message, "shapes_transform") ;

     R = (COMPLEX *) e_malloc(Np*sizeof(COMPLEX),  message) ;
     Rhold = (COMPLEX *) e_malloc(Np*sizeof(COMPLEX),  message) ;
     if (HighRes) {
          Rhigh = (COMPLEX *) e_malloc(Np*sizeof(COMPLEX),message) ;
          Rhihold = (COMPLEX *) e_malloc(Np*sizeof(COMPLEX),message) ;
     }

     /* Feed the input into R */

     for (n=0, pdnump=dnump, pR=R; n < Np; n++, pR++, pdnump++) {
          pR->re = *pdnump ;
          pR->im = 0 ;
     }

     printTwice("Starting forward FFT for shapes...") ;

     do_fft(R, 1) ;			/* FORWARD FFT */

     if (HighRes) {
          hrResetIndices(spread_fac);

          for (n=0, pRhigh=Rhigh; n < Np; n++, pRhigh++) {
               pRhigh->re = 0 ;
               pRhigh->im = 0 ;
          }

          hrDoregSFT(Rhigh, hrNump) ;
     }

     for (n = 0, pR = R, pRhold = Rhold; n < Np; n++, pR++, pRhold++) {
          pRhold->re = pR->re ;
          pRhold->im = pR->im ;
     }
     if (HighRes) {
          for (n = 0, pRhigh = Rhigh, pRhihold = Rhihold; n < Np; 
              n++, pRhigh++, pRhihold++) {
               pRhihold->re = pRhigh->re ;
               pRhihold->im = pRhigh->im ;
          }
     }
     
     	/* First derivatives with respect to x, y, z */

     do_der1(R, Rhold, Rhigh, Rhihold, 1, 0, 0, drdx) ;
     do_der1(R, Rhold, Rhigh, Rhihold, 0, 1, 0, drdy) ;
     do_der1(R, Rhold, Rhigh, Rhihold, 0, 0, 1, drdz) ;

     	/* Second derivatives with respect to xx, xy, xz, yy, yz, zz */

     do_der2(R, Rhold, Rhigh, Rhihold, 1, 0, 0, 0, 0, 0, d2rdxdx) ;
     do_der2(R, Rhold, Rhigh, Rhihold, 0, 1, 0, 0, 0, 0, d2rdxdy) ;
     do_der2(R, Rhold, Rhigh, Rhihold, 0, 0, 1, 0, 0, 0, d2rdxdz) ;
     do_der2(R, Rhold, Rhigh, Rhihold, 0, 0, 0, 1, 0, 0, d2rdydy) ;
     do_der2(R, Rhold, Rhigh, Rhihold, 0, 0, 0, 0, 1, 0, d2rdydz) ;
     do_der2(R, Rhold, Rhigh, Rhihold, 0, 0, 0, 0, 0, 1, d2rdzdz) ;

     printTwice("...Finished back FFT's for shapes.") ;

     e_free(Np*sizeof(COMPLEX), Rhold) ;
     e_free(Np*sizeof(COMPLEX), R) ; 
     if (HighRes) {
          e_free(Np*sizeof(COMPLEX), Rhihold) ;
          e_free(Np*sizeof(COMPLEX), Rhigh) ; 
     }
     freefft() ;
}

void	do_der1(COMPLEX *R, COMPLEX *Rhold, COMPLEX *Rhigh, COMPLEX *Rhihold,
                char hflag, char kflag, char lflag, real *pder)
{
     void	do_fft(COMPLEX *, int) ;
     COMPLEX	*pR, *pRhold, *pRhigh ;
     real	re_part ;
     int	j ;

     for (j = 0, pR = R, pRhold = Rhold; j < Np; j++, pR++, pRhold++) {
          pR->re = pRhold->re ;
          pR->im = pRhold->im ;
     }
     if (HighRes) 
          for (j = 0, pR = Rhigh, pRhold = Rhihold; j < Np; 
               j++, pR++, pRhold++) {
               pR->re = pRhold->re ;
               pR->im = pRhold->im ;
          }
     
     /* Apply the spread factor and factors for 1st derivative */

     sh_combine_octants1(R, delfac, hflag, kflag, lflag) ;

     if (HighRes) {
          sh_combine_octants1(Rhigh, delfac/4., hflag, kflag, lflag) ;

          for (j = 0, pR = R, pRhigh = Rhigh; j < Np; 
               j++, pR++, pRhigh++) {

               pR->re += pRhigh->re ;
               pR->im += pRhigh->im ;
          }
     }

     /* -2pi*i factor */

     for (j = 0, pR = R; j < Np; j++, pR++) {
          re_part = pR->re ;
          pR->re =  TWOPI * pR->im ;
          pR->im = -TWOPI * re_part ;
     }

     /* First derivatives with respect to x, y, z */

     do_fft(R, -1) ;		

     for (j=0, pR=R; j < Np; j++, pR++, pder++) 
          *pder = pR->re *normfac ;

}

void	do_der2(COMPLEX *R, COMPLEX *Rhold, COMPLEX *Rhigh, COMPLEX *Rhihold,
                char hhflag, char hkflag, char hlflag, 
                char kkflag, char klflag, char llflag, real *pder)
{
     void	do_fft(COMPLEX *, int) ;
     COMPLEX	*pR, *pRhigh, *pRhold ;
     real	twopisq = (TWOPI*TWOPI) ;
     int	j ;

     for (j = 0, pR = R, pRhold = Rhold; j < Np; j++, pR++, pRhold++) {
          pR->re = pRhold->re ;
          pR->im = pRhold->im ;
     }
     if (HighRes) 
          for (j = 0, pR = Rhigh, pRhold = Rhihold; j < Np; 
               j++, pR++, pRhold++) {
               pR->re = pRhold->re ;
               pR->im = pRhold->im ;
          }
     
     /* Apply the spread factor and factors for second derivatives  */

     sh_combine_octants2(R, delfac, hhflag, hkflag, hlflag,
                            kkflag, klflag, llflag) ;

     if (HighRes) {
          sh_combine_octants2(Rhigh, delfac/4., hhflag, hkflag, hlflag,
                            kkflag, klflag, llflag) ;

          for (j = 0, pR = R, pRhigh = Rhigh; j < Np; 
               j++, pR++, pRhigh++) {

               pR->re += pRhigh->re ;
               pR->im += pRhigh->im ;
          }
     }

     /* (2pi*i)^2  factor */

     for (j = 0, pR = R; j < Np; j++, pR++) {
          pR->re *= -twopisq ;
          pR->im *= -twopisq ;
     }

     do_fft(R, -1) ;		

     for (j=0, pR=R; j < Np; j++, pR++, pder++) 
          *pder = pR->re * normfac ;
}

	/***********************************************************
	Bring together the exponential function contributions from 
	the 8 octants of the full ellipsoid.  This differs from 
	sum_over_octants, in that R is NOT expanded to the half-
	ellipsoid here;  instead, the factors that would apply if
	it were so expanded are used to "combine" the data.
	***********************************************************/

void	sh_combine_octants1(COMPLEX *Rlocal, real dfac, 
                           char hflag, char kflag, char lflag)
{
     real     one_expfac ;
     int	n, h, k, l ;
     real	fhin, fkin, flin ;
     real	fh, fk, fl ;

     for (n = 0; n < Np; n++) {

          h = n % Nx ;
	  k = (n/Nx) % Ny ;
	  l = n / (Nx*Ny) ;
	  fhin = (real) h ;
	  fkin = (real) k ;
	  flin = (real) l ;

          Tr_ortho_to_fract(fhin, fkin, flin, &fh, &fk, &fl) ;
          one_expfac = exp(-(dfac * nusq(h, k, l))) *
                       (fh * hflag + fk * kflag + fl * lflag) ;

	  if (k > 0) {
	       fhin = (real) h ;
	       fkin = (real) (k-Ny) ;
	       flin = (real) l ;

               Tr_ortho_to_fract(fhin, fkin, flin, &fh, &fk, &fl) ;

               one_expfac += exp(-(dfac * nusq(h, k-Ny, l))) *
                       (fh * hflag + fk * kflag + fl * lflag) ;
          } 
	  if (l > 0){
	       fhin = (real) h ;
	       fkin = (real) k ;
	       flin = (real) (l-Nz) ;

               Tr_ortho_to_fract(fhin, fkin, flin, &fh, &fk, &fl) ;

               one_expfac += exp(-(dfac * nusq(h, k, l-Nz))) *
                       (fh * hflag + fk * kflag + fl * lflag) ;
          } 
	  if ((k > 0) && (l > 0)) {
	       fhin = (real) h ;
	       fkin = (real) (k-Ny) ;
	       flin = (real) (l-Nz) ;

               Tr_ortho_to_fract(fhin, fkin, flin, &fh, &fk, &fl) ;

               one_expfac += exp(-(dfac * nusq(h, k-Ny, l-Nz))) *
                       (fh * hflag + fk * kflag + fl * lflag) ;
          } 
	  if (h > 0) {
	       fhin = (real) (h-Nx) ;
	       fkin = (real) k ;
	       flin = (real) l ;

               Tr_ortho_to_fract(fhin, fkin, flin, &fh, &fk, &fl) ;

               one_expfac += exp(-(dfac * nusq(h-Nx, k, l))) *
                       (fh * hflag + fk * kflag + fl * lflag) ;
          } 
	  if ((h > 0) && (k > 0)) {
	       fhin = (real) (h-Nx) ;
	       fkin = (real) (k-Ny) ;
	       flin = (real) l ;

               Tr_ortho_to_fract(fhin, fkin, flin, &fh, &fk, &fl) ;

               one_expfac += exp(-(dfac * nusq(h-Nx, k-Ny, l))) *
                       (fh * hflag + fk * kflag + fl * lflag) ;
          } 
	  if ((h > 0) && (l > 0)){
	       fhin = (real) (h-Nx) ;
	       fkin = (real) k ;
	       flin = (real) (l-Nz) ;

               Tr_ortho_to_fract(fhin, fkin, flin, &fh, &fk, &fl) ;

               one_expfac += exp(-(dfac * nusq(h-Nx, k, l-Nz))) *
                       (fh * hflag + fk * kflag + fl * lflag) ;
          } 
	  if ((h > 0) && (k > 0) && (l > 0)) {
	       fhin = (real) (h-Nx) ;
	       fkin = (real) (k-Ny) ;
	       flin = (real) (l-Nz) ;

               Tr_ortho_to_fract(fhin, fkin, flin, &fh, &fk, &fl) ;

               one_expfac += exp(-(dfac * nusq(h-Nx, k-Ny, l-Nz))) *
                       (fh * hflag + fk * kflag + fl * lflag) ;
          } 
          (Rlocal+n)->re *= one_expfac ;
          (Rlocal+n)->im *= one_expfac ;

      }
}

void	sh_combine_octants2(COMPLEX *Rlocal, real dfac, 
                           char hhflag, char hkflag, char hlflag,
                           char kkflag, char klflag, char llflag)
{
     real     one_expfac ;
     int	n, h, k, l ;
     real	fhin, fkin, flin ;
     real	fh, fk, fl ;

     for (n = 0; n < Np; n++) {

          h = n % Nx ;
	  k = (n/Nx) % Ny ;
	  l = n / (Nx*Ny) ;
	  fhin = (real) h ;
	  fkin = (real) k ;
	  flin = (real) l ;

          Tr_ortho_to_fract(fhin, fkin, flin, &fh, &fk, &fl) ;

          one_expfac = exp(-(dfac * nusq(h, k, l))) *
                       (fh * fh * hhflag + fh * fk * hkflag + fh * fl * hlflag +
                        fk * fk * kkflag + fk * fl * klflag + 
                        fl * fl * llflag) ;

	  if (k > 0) {
	       fhin = (real) h ;
	       fkin = (real) (k-Ny) ;
	       flin = (real) l ;

               Tr_ortho_to_fract(fhin, fkin, flin, &fh, &fk, &fl) ;

               one_expfac += exp(-(dfac * nusq(h, k-Ny, l))) *
                      (fh * fh * hhflag + fh * fk * hkflag + fh * fl * hlflag +
                       fk * fk * kkflag + fk * fl * klflag + 
                       fl * fl * llflag) ;
          } 
	  if (l > 0){
	       fhin = (real) h ;
	       fkin = (real) k ;
	       flin = (real) (l-Nz) ;

               Tr_ortho_to_fract(fhin, fkin, flin, &fh, &fk, &fl) ;

               one_expfac += exp(-(dfac * nusq(h, k, l-Nz))) *
                      (fh * fh * hhflag + fh * fk * hkflag + fh * fl * hlflag +
                       fk * fk * kkflag + fk * fl * klflag + 
                       fl * fl * llflag) ;
          } 
	  if ((k > 0) && (l > 0)) {
	       fhin = (real) h ;
	       fkin = (real) (k-Ny) ;
	       flin = (real) (l-Nz) ;

               Tr_ortho_to_fract(fhin, fkin, flin, &fh, &fk, &fl) ;

               one_expfac += exp(-(dfac * nusq(h, k-Ny, l-Nz))) *
                      (fh * fh * hhflag + fh * fk * hkflag + fh * fl * hlflag +
                       fk * fk * kkflag + fk * fl * klflag + 
                       fl * fl * llflag) ;
          } 
	  if (h > 0) {
	       fhin = (real) (h-Nx) ;
	       fkin = (real) k ;
	       flin = (real) l ;

               Tr_ortho_to_fract(fhin, fkin, flin, &fh, &fk, &fl) ;

               one_expfac += exp(-(dfac * nusq(h-Nx, k, l))) *
                     (fh * fh * hhflag + fh * fk * hkflag + fh * fl * hlflag +
                      fk * fk * kkflag + fk * fl * klflag + 
                      fl * fl * llflag) ;
          } 
	  if ((h > 0) && (k > 0)) {
	       fhin = (real) (h-Nx) ;
	       fkin = (real) (k-Ny) ;
	       flin = (real) l ;

               Tr_ortho_to_fract(fhin, fkin, flin, &fh, &fk, &fl) ;

               one_expfac += exp(-(dfac * nusq(h-Nx, k-Ny, l))) *
                     (fh * fh * hhflag + fh * fk * hkflag + fh * fl * hlflag +
                      fk * fk * kkflag + fk * fl * klflag + 
                      fl * fl * llflag) ;
          } 
	  if ((h > 0) && (l > 0)){
	       fhin = (real) (h-Nx) ;
	       fkin = (real) k ;
	       flin = (real) (l-Nz) ;

               Tr_ortho_to_fract(fhin, fkin, flin, &fh, &fk, &fl) ;

               one_expfac += exp(-(dfac * nusq(h-Nx, k, l-Nz))) *
                     (fh * fh * hhflag + fh * fk * hkflag + fh * fl * hlflag +
                      fk * fk * kkflag + fk * fl * klflag + 
                      fl * fl * llflag) ;
          } 
	  if ((h > 0) && (k > 0) && (l > 0)) {
	       fhin = (real) (h-Nx) ;
	       fkin = (real) (k-Ny) ;
	       flin = (real) (l-Nz) ;

               Tr_ortho_to_fract(fhin, fkin, flin, &fh, &fk, &fl) ;

               one_expfac += exp(-(dfac * nusq(h-Nx, k-Ny, l-Nz))) *
                     (fh * fh * hhflag + fh * fk * hkflag + fh * fl * hlflag + 
                      fk * fk * kkflag + fk * fl * klflag + 
                      fl * fl * llflag) ;
          } 
          (Rlocal+n)->re *= one_expfac ;
          (Rlocal+n)->im *= one_expfac ;

      }
}
void	hrRegInput(char *vox_filename)
{
     int	hrGetLength(char *) ;
     void	hrReadLists(char *, real *) ;

     sprintf(list_filename, "%s.list", vox_filename) ;
     Nhr = hrGetLength(list_filename) ;
     hrNump = (real *) e_malloc(Nhr*sizeof(real), caller) ;
     hrReadLists(list_filename, hrNump) ;	
     sprintf(message, "including %d high-resolution points.", Nhr) ;
     printTwice(message) ;
}


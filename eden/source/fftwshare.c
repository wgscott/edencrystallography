/*************************************************************************

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

				FFTWSHARE.C

	Title:		FFT-preparing functions (in C) for Eden codes.
	Author:		Hanna Szoke
	Date:		Mar. 31, 1995
	Function:	This package contains initializations,
			low-level routines (extract, insert) and
			the interface to the fft solver fftw.

*********************************************************************/

#include "util.h"
#include "dims.h"


static	COMPLEX	*fftarray ;			
static	COMPLEX	*temp ;

static	fftwnd_plan	planf, planb ;

static	void	convolve0(real *, real *, COMPLEX *) ;
static	void	convolveb(real *, COMPLEX *, COMPLEX *) ;
static	void	spread_over_Nhkl(COMPLEX *, COMPLEX *);
static	void	deconv0(real *, real *, COMPLEX *) ;
static	void	deconvb(real *, COMPLEX *, COMPLEX *) ;
static	void	sum_over_octants(COMPLEX *, COMPLEX *) ;

void initfft()
{
     int	flags ;
     FILE	*fp_wis ;
     char	Exist = FALSE ;

/*     Use this for quick tests: flags = FFTW_ESTIMATE|FFTW_IN_PLACE ; */

     printTwice("Begun making FFT plan ...") ;

     if ((fp_wis = fopen("fft_wis", "r")) != NULL) {
          flags = FFTW_IN_PLACE|FFTW_USE_WISDOM ;
          Exist = TRUE ;
          fftw_import_wisdom_from_file(fp_wis) ;
     }
     else
          flags = FFTW_MEASURE|FFTW_IN_PLACE|FFTW_USE_WISDOM ;

     if ((planf = fftw3d_create_plan(Nz, Ny, Nx, FFTW_FORWARD, flags)) == NULL)
          EdenError("Problem creating FFTW forward plan") ;

     if (very_verbose)
          fftwnd_fprint_plan(fp_log, planf) ;

     if (Exist == FALSE) { 
          fp_wis = fopen("fft_wis", "w") ;
          if (fp_wis == NULL)
               EdenError("Cannot open fft_wis file for writing") ;
          fftw_export_wisdom_to_file(fp_wis) ;
          fflush(fp_wis) ;
          fclose(fp_wis) ; 
     }

     flags = FFTW_IN_PLACE|FFTW_USE_WISDOM ;

     if ((planb = fftw3d_create_plan(Nz, Ny, Nx, FFTW_BACKWARD, flags)) == NULL)
	  EdenError("Problem creating FFTW backward plan") ;

     printTwice(" ... Finished making FFT plan.") ;
} 
void initfft_aux()
{

	/* Initialize auxiliary arrays for fft manipulations.  
	   This is required for Back, Forth and Solve, but not for 
	   Ran or Regrid (which also do fft's).  */

     strcpy(message, "initfft_aux") ;

     fftarray = (COMPLEX *) e_malloc(Nhkl*sizeof(COMPLEX), message) ;
     temp = (COMPLEX *) e_malloc(Np*sizeof(COMPLEX), message) ;
}
void freefft()
{
	/* free up the plans and memory associated with them.	*/

	fftwnd_destroy_plan(planf) ;
	fftwnd_destroy_plan(planb) ;
}

static
int goodval[] =       {   2,   4,   6,   8,  10,  12,  14,  16,  
			 18,  20,  22,  24,  26,  28,  30,  32,  
			 36,  40,  42,  44,  48,  50,  52,  54,  
			 56,  60,  64,  66,  70,  72,  78,  80,  
			 84,  88,  90,  96,  98, 100, 104, 108, 
		        110, 112, 120, 126, 128, 130, 132, 140,
			144, 150, 154, 156, 160, 162, 168, 176,
			180, 182, 192, 196, 198, 200, 208, 210,
			216, 220, 224, 234, 240, 250, 252, 256, 
			260, 264, 270, 280, 286, 288, 294, 300, 
			308, 312, 320, 324, 330, 336, 350, 352,
			360, 364, 378, 384, 390, 392, 396, 400, 
			416, 420, 432, 440, 448, 450, 462, 468,
			480, 486, 490, 500, 512, 1024} ;
#define	JMAX	110	   /* dimension of goodval[] */

#define	JTOP	1024	/* max dimension (Nx, Ny or Nz)	*/

int 	getFftDim(float x, int M)	
			/* Expand x to the closest number that is an even 
			product of powers of 2, 3 and 5  and divisible by M. */

{
   int Mgoodval[JMAX] ;

   int	i = 0, j ;
   int	imax ;
   int	closest_i ;

   for (imax = 0, j = 0; j < JMAX; j++)
	if ((goodval[j] % M) == 0 ) 
             Mgoodval[imax++] = goodval[j] ;

   if (x > JTOP) {
	sprintf(message, " %g is too large a dimension - maximum is %d", 
		x, JTOP) ;
        EdenError(message) ;
   }
   else if (x <= Mgoodval[0]) 
	return (Mgoodval[0]);

   while (x > Mgoodval[i])
      i++ ;

   closest_i = ((Mgoodval[i]-x) < (x - Mgoodval[i-1])) ? i : (i-1) ;

   return (Mgoodval[closest_i]) ;
}

void	convolve(real *array, real *efac, COMPLEX *cefac, COMPLEX *prod)
{
     int	n ;

     for (n = 0; n < Nhkl; n++) {
	  (prod+n)->re = 0. ;
	  (prod+n)->im = 0. ;
     }

     convolve0(array, efac, prod) ;

     if (grid_type == BODY_CENTERED)
          convolveb(array+Np, cefac, prod) ;
}

void convolve0(real *array, real *efac, COMPLEX *prod) 
{
	/***********************************************************
	 Accumulate complex product of DFT(array) and (real) efac.

	 On entry:
		array = Current value of the image.
	        efac = real array to multiply Fourier transform
	 On return:
	 	prod = complex product
	***********************************************************/

        int	i ;
	COMPLEX	*pfft ;
	COMPLEX	*ptemp, *pprod ;
	real	*pefac ;
	real	*pin ;
	void	do_fft(COMPLEX *, int) ;

        for (i = 0, pin=array, ptemp=temp; i < Np; i++, pin++, ptemp++) {
	     ptemp->re = *pin ;
	     ptemp->im = 0 ;
        }

        do_fft(temp, 1) ;
	spread_over_Nhkl(temp, fftarray) ;

	/* We return (not accumulate) the product of fftarray*efac */

        for (i = 0, pfft=fftarray, pprod=prod, pefac=efac; i < Nhkl; 
	     i++, pfft++, pprod++, pefac++) 
        {
	     pprod->re = *pefac * pfft->re ;
	     pprod->im = *pefac * pfft->im ;
        }
}
  
void convolveb(real *array, COMPLEX *efac, COMPLEX *prod) 
{
	/***********************************************************
	 Accumulate complex product of DFT(array) and (complex) efac.

	 On entry:
		array = Current value of the image.
	        efac = Complex array to multiply Fourier transform
	 On return:
	 	prod = complex product

	***********************************************************/

        int	i ;
	COMPLEX	*rfft, *ifft ;
	COMPLEX	*ptemp, *pprod ;
	COMPLEX	*refac, *iefac ;
	real	*pin ;
	void	do_fft(COMPLEX *, int) ;

        for (i = 0, pin=array, ptemp=temp; i < Np; i++, pin++, ptemp++) {
	     ptemp->re = *pin ;
	     ptemp->im = 0 ;
        }

        do_fft(temp, 1) ;
	spread_over_Nhkl(temp, fftarray) ;

	/* We accumulate the product of fftarray*efac */

	rfft = ifft = fftarray ;
	refac = iefac = efac ;

        for (i = 0, pprod = prod; i < Nhkl; 
	     i++, pprod++, rfft++, ifft++, refac++, iefac++) 
	{
	     pprod->re += refac->re * rfft->re - iefac->im * ifft->im ;
	     pprod->im += refac->re * ifft->im + iefac->im * rfft->re ;
        }
}

void	spread_over_Nhkl(COMPLEX *inarray, COMPLEX *outarray)
{
	/***************************************************************
	Rearrange info; inarray contains data packed in the positive 
	octant from (0,0,0) to (Nx,Ny,Nz), excluding the upper limits. 
	Each element (h,k,l) is placed in outarray spanning (0,-Ny,-Nz) 
	to (Nx,Ny,Nz), using in-line coding that effectively does assemble()
	but it faster.

	Then find the elements of inarray corresponding to negative k's 
	and/or l's, and expand the data appropriately.  Note "in" is an 
	index into the packed inarray, "nn" is an index into the full 
	Nhkl outarray.

	June 4, 2004: This code "looks wrong" - also sum_over_octants().
	But it was tested by incrementing a counter(nn) whenever an
	nn is set; then all cases in which counter is not 1 were reported.
	Only (h, -Ny, k) and (h, k, -Nz) were unset - as is correct, and
	in no case was counter(nn) = 2.
	***************************************************************/

      int	n, k, l, nn, in, thisn ;

      for (n = 0; n < Nhkl; n++) {
	   (outarray+n)->re = 0 ;
	   (outarray+n)->im = 0 ;
      }

      for (n = 0; n < Np; n+=Nx) {

           k = (n/Nx) % Ny ;
           l = n / (Nx*Ny) ;

           nn = k * Nh + l * Nh * Nk ;
	   in = n ;
 
	   for (thisn = 0; thisn < Nx; thisn++, nn++, in++) {
 
                (outarray+nn)->re = (inarray+in)->re ;
                (outarray+nn)->im = (inarray+in)->im ;
 
           }

           if ((k > 0)) {

                nn = (Nk-k) * Nh + l * Nh * Nk ;
                in = (Ny-k) * Nx + l * Nx * Ny;

	        for (thisn = 0; thisn < Nx; thisn++, nn++, in++) {
 
                     (outarray+nn)->re = (inarray+in)->re ;
                     (outarray+nn)->im = (inarray+in)->im ;

                }
           }
           if ((l > 0)) {

                nn = k * Nh + (Nl-l) * Nh * Nk ;
                in = k * Nx + (Nz-l) * Nx * Ny ;

	        for (thisn = 0; thisn < Nx; thisn++, nn++, in++) {
 
                     (outarray+nn)->re = (inarray+in)->re ;
                     (outarray+nn)->im = (inarray+in)->im ;

                }
           }
           if ((k > 0) && (l > 0)) {

                nn = (Nk-k) * Nh + (Nl-l) * Nh * Nk ;
                in = (Ny-k) * Nx + (Nz-l) * Nx * Ny;

	        for (thisn = 0; thisn < Nx; thisn++, nn++, in++) {
 
                     (outarray+nn)->re = (inarray+in)->re ;
                     (outarray+nn)->im = (inarray+in)->im ;

                }
           }
      }
}
 
void	deconvolve(real *garray, real *efac, COMPLEX *cefac, COMPLEX *kern)
{
     void deconv0(real *, real *, COMPLEX *) ;
     void deconvb(real *, COMPLEX *, COMPLEX *) ;

	deconv0(garray, efac, kern) ;

	if (grid_type == BODY_CENTERED)
	     deconvb(garray+Np, cefac, kern) ;
}

void	deconv0(real *garray, real *efac, COMPLEX *kern) 
{
	/***********************************************************
	 Calculate the gradient from kern. 
	
	 On entry:
	       efac = real array to multiply Fourier transform
	       kern = r, weighted & with exp_phi
	 On return:
		garray  = gradient 
	***********************************************************/

	int i ;
	real	*pefac, *pout ;
	COMPLEX	*pkern ;
	COMPLEX	*pfft ;
	COMPLEX	*ptemp ;
	void	do_fft(COMPLEX *, int) ;

        for (i = 0, pfft=fftarray, pkern=kern, pefac=efac; i < Nhkl; 
	     i++, pfft++, pkern++, pefac++) 
	{
             pfft->re = *pefac * pkern->re ;
             pfft->im = *pefac * pkern->im ;
        }

        sum_over_octants(fftarray, temp) ;
        do_fft(temp, -1) ;

        for (i = 0, ptemp=temp, pout=garray; i < Np; i++, ptemp++, pout++) 
	     *pout += ptemp->re ;
}

void	deconvb(real *garray, COMPLEX *efac, COMPLEX *kern) 
{
	/***********************************************************
	 Calculate the gradient from kern. 
	
	 On entry:
	       efac = complex array to multiply Fourier transform
		kern = r, weighted & with exp_phi
	 On return:
		garray  = gradient 
	       
	***********************************************************/

	int i ;
	COMPLEX	*pefac ;
	COMPLEX	*pkern ;
	COMPLEX	*pfft ;
	COMPLEX	*ptemp ;
	real	*pout ;
	void	do_fft(COMPLEX *, int) ;

        for (i = 0, pfft=fftarray, pkern=kern, pefac=efac; i < Nhkl; 
	     i++, pfft++, pkern++, pefac++) 
	{
             pfft->re = pefac->re * pkern->re + pefac->im * pkern->im ;
             pfft->im = pefac->re * pkern->im - pefac->im * pkern->re ;
        }

        sum_over_octants(fftarray, temp) ;
        do_fft(temp, -1) ;

        for (i = 0, ptemp=temp, pout=garray; i < Np; i++, ptemp++, pout++) 
	     *pout += ptemp->re ;
}

void	sum_over_octants(COMPLEX *inarray, COMPLEX *outarray)
{
	/***********************************************************
	This function pulls in the contributions from the four octants
	of the h>=0 hemi-ellipsoid -- (h, +/-k, +/-l) -- and stashes 
	the result in an output array measuring Nx by Ny by Nz.
	Note that this is done for the gradient; the +/- sign is
	built in to "inarray" (aka R), so the contributions from the
	octants are simply added together.
	***********************************************************/

      int	n, k, l, nn, in, thisn ;

      for (n = 0; n < Np; n++) {
	   (outarray+n)->re = 0 ;
	   (outarray+n)->im = 0 ;
      }

      for (n = 0; n < Np; n+= Nx) {

	   k = (n/Nx) % Ny ;
	   l = n / (Nx*Ny) ;

	   nn = k * Nh + l * Nh * Nk ;
	   in = n ;

	   for (thisn = 0; thisn < Nx; thisn++, nn++, in++) {

	        (outarray+in)->re += (inarray+nn)->re ;
	        (outarray+in)->im += (inarray+nn)->im ;

           }
	   if ((k > 0)) {

	        nn = (Nk-k) * Nh + l * Nh * Nk ;
                in = (Ny-k) * Nx + l * Nx * Ny;

	        for (thisn = 0; thisn < Nx; thisn++, nn++, in++) {

	             (outarray+in)->re += (inarray+nn)->re ;
	             (outarray+in)->im += (inarray+nn)->im ;

                }
           }
	   if ((l > 0)) {
	        nn = k * Nh + (Nl-l) * Nh * Nk ;
                in = k * Nx + (Nz-l) * Nx * Ny ;

	        for (thisn = 0; thisn < Nx; thisn++, nn++, in++) {

	             (outarray+in)->re += (inarray+nn)->re ;
	             (outarray+in)->im += (inarray+nn)->im ;

                }
           }
	   if ((k > 0) && (l > 0)) { 
	        nn = (Nk-k) * Nh + (Nl-l) * Nh * Nk ;
                in = (Ny-k) * Nx + (Nz-l) * Nx * Ny;

	        for (thisn = 0; thisn < Nx; thisn++, nn++, in++) {

	             (outarray+in)->re += (inarray+nn)->re ;
	             (outarray+in)->im += (inarray+nn)->im ;

                }
           }
      }
}
	 
void	do_fft(COMPLEX *inarray, int isign)
{
	/***************************************************************
	Following usual convention, isign is -1 for forward, 1 for back 
	transform.  Last 3 arguments are ignored since we do in-place
	transforms.
	***************************************************************/

      int	howmany = 1 ;        

      if (isign == -1)
           fftwnd(planf, howmany, (FFTW_COMPLEX *) inarray, 1, 1, 0, 0, 0) ;
      else
           fftwnd(planb, howmany, (FFTW_COMPLEX *) inarray, 1, 1, 0, 0, 0) ;
}


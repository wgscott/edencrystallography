/******************************************************************************

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

                                REMAP.C

  Title:        Utility package for remapping.
  Author:       Hanna Szoke
  Date:         Early parts: 1993, later parts: 1998.
  Function:     Functions for remapping electron/voxel values. 

			OLD PARTS, dealing with gridded info: 

void prepare_relative_index_list()	- lists of close indices, exponents
void set_garrays() 	- worker function to put values into lists

void prepare_off_array()- transfer list to indexlist, for Sayre.
void add_to_indexlist()	- worker function for indexlist preparation
int  searchlist()	- worker function for indexlist preparation
void check_off_array()  - ? does indexlist work ?
	
... plus some trivial worker functions.  
******************************************************************************/

#include "util.h"
#include "cellparams.h"
#include "dims.h"
#include "symmetry.h"
#include "remap.h"		/* for close vicinity procedures */

	/* declarations of remap.h variables */

int	Nb = 0 ;	/* Number of vicinity points */
int	*indexlist ;	/* array for fast indexing to vicinity pts. */
int	*poff ;		/* Nptotal array of pointers to indexlist */
real	*gexp ;		/* array for vicinity exponential factors */

	/* declarations local to remap.c     */

static  int     *gix, *giy, *giz, *gNp ;   /* arrays for vicinity indices */
static	real  denom ;		/* 1. / (dconst*eta*dr*dr) */

static	void	set_garrays(int *, int, int, int, int, real) ;
static	int	searchlist(int, int *) ;
static	void	add_to_indexlist(int *, int *) ;

void prepare_relative_index_list(real dconst)

	/* We use dconst = 1 generally, but dconst = 2 for NCS.
	Derivation gives 2, but 1.5 or 1.0 makes it less smeared out. */
{
     int	i, j, k ;	/* looping indices around (0, 0, 0) */
     real	fx, fy, fz ;	/* fractional coords */
     real	x, y, z ;	/* orthogonal coords of nearby r */
     real	val ;
     real	rsq ;		/* |r - ref|^2 */
     real	rsqlim ;

     /************************************************************
     Set up a table of exponents over a radius r such that 
     exp(-(r/d)^2) >= CUTOFF, defined in remap.h.
     ************************************************************/

     rsqlim = -log(CUTOFF) ;
     denom = 1. / (dconst*eta*dr*dr) ;

     strcpy(message, "prepare_relative_index_list") ;

     gexp = (real *) e_malloc(NGMAX*sizeof(real), message) ;

     gix = (int *) e_malloc(NGMAX*sizeof(int), message) ;
     giy = (int *) e_malloc(NGMAX*sizeof(int), message) ;
     giz = (int *) e_malloc(NGMAX*sizeof(int), message) ;
     gNp = (int *) e_malloc(NGMAX*sizeof(int), message) ;

     /************************************************************
     We take as reference point (0 0 0) and look for nearby points.
     NOTE: PROCEDURE CHANGED 1/8/98 - indices (i, j, k) are stored 
     in the gix etc. rather than their absolute values.
     ************************************************************/

     for (k = -IRANGE; k <= +IRANGE; k++) {

          fz = k / floatNz ;

	  for (j = -IRANGE; j <= +IRANGE; j++) {

               fy = j / floatNy ;

               for (i = -IRANGE; i <= +IRANGE; i++) {

                    fx = i / floatNx ;

                    fract_to_ortho(fx, fy, fz, &x, &y, &z) ;

                    rsq = x*x + y*y + z*z ;

		    if ((val = rsq * denom) <= rsqlim) 
			 set_garrays(&Nb, i, j, k, 0, val) ;
               }
          }
     }

     if (grid_type == BODY_CENTERED) {

          for (k = -(IRANGE+1); k <= +IRANGE; k++) {

               fz = (k+0.5) / floatNz ;

	       for (j = -(IRANGE+1); j <= +IRANGE; j++) {

                    fy = (j+0.5) / floatNy ;

                    for (i = -(IRANGE+1); i <= +IRANGE; i++) {

                         fx = (i+0.5) / floatNx ;

                         fract_to_ortho(fx, fy, fz, &x, &y, &z) ;

                         rsq = x*x + y*y + z*z ;

		         if ((val = rsq * denom) < rsqlim) 
                              set_garrays(&Nb, i, j, k, Np, val) ;
                    }
               }
          }
     }

     sprintf(message, "Length of close vicinity arrays = %d", Nb) ;
     printTwice(message) ;

     /***  if (very_verbose)   ***/
     if (verbose)  {
	  fprintf(fp_log, "Close Vicinity Arrays\n") ;
          for (i = 0; i < Nb; i++)
	       fprintf(fp_log, 
	       "i = %d, (ix, iy, iz, N) = %d %d %d %d, gexp = %g\n",
	       i, *(gix+i), *(giy+i), *(giz+i), *(gNp+i), *(gexp+i)) ;
     }

}

void set_garrays(int *g, int i, int j, int k, int lat, real val) 
{
     *(gix + *g)= i ; 
     *(giy + *g)= j ;
     *(giz + *g)= k ;
     *(gNp + *g)= lat ;
     *(gexp + *g)= exp(-val) ;

     (*g)++ ;
     if (*g >= NGMAX) {
         sprintf(message, "Too many exp factors! - max = %d", NGMAX) ;
         EdenError(message) ;
     }
}

void prepare_off_array()
{
     int        n ;
     int        i, j, k, l ;
     int        m ;
     int        offset ;
     int        *bpoff ;
     int        Nindexlist = 0 ;
     int        nfound ;
     int        *fastind ;
					 
					  
     /************************************************************
     This function prepares indexlist, a list of relative indices
     for a given point to its close vicinity points, and poff, an
     Nptotal array of indices for accessing indexlist.
     For points well in the interior of the unit cell, the relative
     indices to close vicinity points may be used directly, be they
     positive or negative.  But at points close to the edges of the
     unit cell, a more careful indexing is required because of
     wrap-around. These data structures are needed for fast access
     to vicinity points when doing the Sayre cost term.
     ************************************************************/
		
     strcpy(message, "prepare_off_array") ;
	      
     poff = (int *) e_malloc(Nptotal*sizeof(int), message) ;
     indexlist = (int *) e_malloc(Nb*NGMAX*sizeof(int), message) ;
     fastind = (int *) e_malloc(Nb*sizeof(int), message) ;
				      
     for (n = 0, bpoff = poff; n < Nptotal; n++, bpoff++) {
						  
          i = n % Nx ;
          j = (n/Nx) % Ny ;
          k = n / (Nx*Ny) % Nz ;
          l = (n < Np)? 0 : Np ;
	   
          for (m = 0; m < Nb; m++) {

	       offset = *(gix + m) + *(giy + m) * Nx + *(giz + m) * Nx * Ny
		        + *(gNp + m) ;
		
	       if ((i + *(gix + m)) < 0)
		   offset += Nx ;
	       if ((j + *(giy + m)) < 0)
	           offset += Nx*Ny ;
	       if ((k + *(giz + m)) < 0)
	           offset += Nx*Ny*Nz ;

               if ((i + *(gix + m)) > Nx-1)
		   offset -= Nx ;
	       if ((j + *(giy + m)) > Ny-1)
	           offset -= Nx*Ny ;
	       if ((k + *(giz + m)) > Nz-1)
	           offset -= Nx*Ny*Nz ;

	       if ((l + *(gNp + m)) == Nptotal)
	           offset -= Nptotal ;

               *(fastind + m) = offset ;

          }
	  
	  /* Have we already encountered this set of fast indices?
          if not - add it to indexlist.  */
			       
	  if ((nfound = searchlist(Nindexlist, fastind)) == Nindexlist) {
               add_to_indexlist(&Nindexlist, fastind) ;
	       Nindexlist++ ;
          }

          *bpoff = Nb * nfound ;

     }
	          
     if (very_verbose) {
          sprintf(message, "Total # of fast indexing sets: %d", Nindexlist) ;
	  printTwice(message) ;
     }
     
     e_free(NGMAX*sizeof(int), gix) ;
     e_free(NGMAX*sizeof(int), giy) ;
     e_free(NGMAX*sizeof(int), giz) ;
     e_free(NGMAX*sizeof(int), gNp) ;
     e_free(Nb*sizeof(int), fastind) ;

     return ; 
}

void add_to_indexlist(int *Nindexlist, int *fastind)
{
			/* Add one set of direct indices to *indexlist */
     int	*entry ;
     int	m ;

     entry = indexlist + *Nindexlist * Nb ;

     for (m = 0; m < Nb; m++)
          *(entry + m) = *(fastind + m) ;

     if (very_verbose) {
          fprintf(fp_log, "\nIn add_to_indexlist ...\n") ;

          for (m = 0; m < Nb; m++)
     	       fprintf(fp_log, "Entry: %d\n", *(fastind+m)) ;
     }

     return ;
}
int searchlist(int Nindexlist, int *newind)
{
     int	m, n ;
     int	gotit ;
     int	*old ;
     int	*new ;

     for (n = 0, old = indexlist; n < Nindexlist; n++) {
	  for (m = 0, new = newind, gotit = TRUE; m < Nb; m++, new++, old++) {
	       if (*old != *new) 
		    gotit = FALSE ;
          }
	  if (gotit) 
	       return (n) ;
     }
     return(Nindexlist) ;
}
void check_off_array()
{
        int     n, m ;
	int     close_ind ;
	int     *bpoff ;
	int     *nindexlist ;
	int     badstuff = FALSE ;

        /********************************************************************** 
	We use poff for accessing the appropriate indexlist entries, for
	direct (fast) indexing into the close vicinity of each n.
	This is the mode in which poff will be used in the cost functions.
	***********************************************************************/
	bpoff = poff ;

	for (n = 0; n < Nptotal; n++, bpoff++) {
	     for (m = 0, nindexlist = indexlist + *bpoff;
		  m < Nb; m++, nindexlist++) {

		  close_ind = (n + *nindexlist) ;

                  if ((close_ind < 0) || (close_ind >= Nptotal)) {
		       badstuff = TRUE ;
		       break ;
		  }

             }
	     if (badstuff) {
		  sprintf(message,
		  "Bad index! - n=%d, offset into indexlist=%d", n, *bpoff) ;
		  EdenWarning(message) ;

		  for (m = 0, nindexlist = indexlist + *bpoff; m < Nb;
		       m++, nindexlist++) {

                       sprintf(message, "m=%d, offset=%d", m, *nindexlist) ;
		       EdenWarning(message) ;
                  }

                  EdenError("Quitting.") ;

	     }
        }

        return ;
}




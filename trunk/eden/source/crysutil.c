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

                                CRYSUTIL.C

  Title:        Crystal Computational Utilities
  Author:       Hanna Szoke
  Date:         Combined with symutil.c Dec. 2000. 
  Function:	This package contains functions reading, interpreting
		and using symop.lib, for setting up indexing information, 
		for symmetrizing arrays and for propagating values from one 
		asymmetric unit to the rest of the crystal.
                It does low-level computations related to crystal parameters 
		a, b, c, alpha, beta and gamma.

		function			called from ...
		---------			----------------
		fetchSymop() 			basics.c
		prepareMatrix()			internal 
		apply_symop()			sym.c, ncsutil.c, setupncs.c
		apply_symop_nowrap()		pdbutil.c
		apply_symop_hkl() 		cost.c, expandf[co].c, hklutil.c
		setup_csind()			back.c, qshare.c
		get_equiv_list()		internal 
		do_cs_symmetrization()		cost.c
		do_inplace_cs_symmetrization() 	back.c, solve.c
	  	cs_expand()			setupncs.c
		setupNusq()			basics.c, pdbutil.c
		nusq()				(ubiquitous)
		prepare_F_and_Finv()			"     
		volfac()			basics.c
		Rsq()				sym.c
		ortho_to_fract()		(ubiquitous)
		fract_to_ortho()			"
		index_to_fract()			"
		fract_to_index()			"
		Tr_ortho_to_fract()		shapes.c

*******************************************************************************/

#include "util.h"
#include "cellparams.h"
#include "symmetry.h"
#include "dims.h"

/* definitions of cellparams.h declarations (input and derived structures) */

static	real	F[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}} ;		
					/* transformation matrix, 
					orthogonal->fractional coords */
static	real	Finv[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}} ;	
					/* inverse of F */
static	real	Ftr[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}} ;	
					/* transpose of F */

/*	local to this package (coefficients for Nusq)	*/

static	real	cosalpha, cosbeta, cosgamma ;
static	real	cohh, cokk, coll, cokl, colh, cohk ;

static	char	*unique_pmask ;	/* mask of one asymmetric unit in real space */
static	int	*equiv_list ;	/* list of equivalent indices */

#define	COSEPS	1.e-10 
#define	SMALLSTRING	20
#define	ESYM	1.e-1		/* relative measure of asymmetry	*/
#define ABSSYM	1.e-4		/* absolute limit for considering asymmetry */

real	*matop ;		/* to hold matrix operators for doing
				transformations between points related 
				by crystal symmetry */
int	Nprim ;
char	crystal_system[SMALLSTRING] ;	/* 1 of the 7: 
				TRICLINIC, MONOCLINIC, ORTHORHOMBIC, 
				TETRAGONAL, TRIGONAL, HEXAGONAL or CUBIC */

int	laue_group ;		/* one of 11: a subdivision of crystal_system 
				used to identify unique data in (h,k,l) */
int	maxdenom[3] ;		/* max. denominators of translation terms, 
				needed for ensuring good Nx, Ny and Nz */
char	SGname[MAXSTRING] ;
char	PGname[MAXSTRING] ;

int	Ncs ;			/* Number of asymmetric units in crystal */
int	Ics ;			/* # of unique grid points in an asym. unit */

static	void	get_equiv_list(int) ;
static	void	prepareMatrix(real *) ;

	/**********************************************************
	Read the ccp4 file, symop.lib and retrieve the information 
	matching the space group identified by user's input string 
	mySGname.   
	**********************************************************/

int	fetchSymop(char *mySGname) 
{
     FILE	*fp ;
     int	n ;
     int	Nsg ;
     int	gotcha = FALSE ;
     char	*SymGname, *PointGname ;

     for (n = 0; n < 3; n++) 
	  maxdenom[n] = 1. ;	/* initialization, for prepareMatrix() calls */

     strcpy(message, getenv("EDENHOME")) ;
     strcat(message, "/source/symop.lib");
     check_exist(message) ;

     fp = fopen(message, "r") ;

     SymGname = (char *)e_malloc(SMALLSTRING*sizeof(char), caller) ;
     PointGname = (char *)e_malloc(SMALLSTRING*sizeof(char), caller) ;

     while ((fgets(nextline, MAXSTRING, fp) != NULL) && 
	    ((int)strlen(nextline) > 1)) {

         sscanf(nextline, "%d %d %d %s %s %s", 
                &Nsg, &Ncs, &Nprim, SymGname, PointGname, crystal_system) ;

	 if (strcmp(SymGname, mySGname) == 0) {
	     gotcha = TRUE ;
	     strcpy(SGname, SymGname) ;
	     strcpy(PGname, PointGname) ;

             matop = (real *) e_malloc(Ncs*MEL*sizeof(real), "crysutil") ;

	     /* identify the crystal class and the Laue group;
	     see Practical Protein Crystallography, Duncan McRee, p. 65 or 
	     Giacovazzo, pp. 14 and 24 */

	     if (Nsg < 3) 
	          laue_group = 1 ;	/* triclinic */
             else if (Nsg < 16) 
		  laue_group = 2 ;	/* monoclinic */
             else if (Nsg < 75)
		  laue_group = 222 ;	/* orthorhombic */
             else if (Nsg < 89)
		  laue_group = 4 ;	/* 1st class of tetragonals */
             else if (Nsg < 143)
		  laue_group = 422 ;	/* 2nd class of tetragonals */
             else if (Nsg < 149)
		  laue_group = 3 ;	/* trigonals -3 */
             else if (Nsg < 168)
		  laue_group = 32 ;	/* trigonals -32 */
             else if (Nsg < 177)
		  laue_group = 6 ;	/* hexagonals -6 */
             else if (Nsg < 195)
		  laue_group = 622 ;	/* hexagonals -622 */
             else if (Nsg < 207)
		  laue_group = 23 ;	/* 1st class of cubics */
             else 		
		  laue_group = 432 ;	/* 2nd class of cubics */
          }
          if (gotcha)
              break ;
     }

     e_free(SMALLSTRING*sizeof(char), SymGname) ;
     e_free(SMALLSTRING*sizeof(char), PointGname) ;

     if (!gotcha) {
          sprintf(message, " Cannot recognize space group %s", mySGname) ;
          EdenError(message) ;
     }


     for (n = 0; n < Ncs; n++) {
          strcpy(nextline, " ") ;
          if ((fgets(nextline, MAXSTRING, fp) != NULL) && 
              ((int)strlen(nextline) > 1)) 

                prepareMatrix(matop+MEL*n) ; 
     }


     return(Nsg) ;
}

void	prepareMatrix(real *mat)
{
     int	i,  n ;
     int	sign ;
     int	numer ;
     int	denom ;
     char	*next ;
     real	*pmat ;

     for (n = 0, pmat = mat, next = nextline; n < 3; n++, next++, pmat += 4) {

          sign = 1 ;
	  numer = 0 ;
	  denom = 0 ;

          for (i = 0; i < 4; i++)
               *(pmat+i) = 0 ;

          while ((*next != ',') && (*next != '\n')) {

	       switch (*next) {
	       case 'X' :
		    *pmat = sign ;
		    break;
	       case 'Y' :
		    *(pmat+1) = sign ;
		    break;
	       case 'Z' :
		    *(pmat+2) = sign ;
		    break;
	       case '-' :
		    sign = -1 ;
		    break;
	       case '+' :
		    sign = 1 ;
		    break;
              case '1': case '2': case '3': case '4': case '5': 
              case '6': case '7': case '8': case '9':
		    if (numer == 0) {
			 sscanf(next, "%1d", &numer) ;
                    }
		    else {
			 sscanf(next, "%1d", &denom) ;
		         *(pmat+3) = sign * numer / (real) denom ;
			 if (maxdenom[n] < denom)
			      maxdenom[n] = denom ;
                    }
		    break ;
	      default:			/* includes '/'	*/
		    break;
              }
              next++ ;
          }
     }
     if (very_verbose) {  
	for (n = 0, pmat=mat; n < 3; n++, pmat += 4) 
	  fprintf(fp_log, 
	    "mat[%1d][0] - [%1d][3] = %g, %g, %g, %g maxdenom = %d\n",
	    n, n, *(pmat), *(pmat+1), *(pmat+2), *(pmat+3), maxdenom[n]) ;
     } 
     return ;
}

void	apply_symop(real fx, real fy, real fz, real *mat, 
		    real *newfx, real *newfy, real *newfz)

	/**********************************************************
	fx, fy, fz are input fractional coords, not necessarily in
		range [0,1) 
	*mat is	rotation/translation matrix
	*newfx, *newfy, *newfz are output rotated and translated 
		fractional coordinates, in range [0,1]	
	**********************************************************/
{
     real	fclip(real, real) ;
     real	ffx, ffy, ffz ;
     real	rtx, rty, rtz ;
     real	one = 1.000 ;

     ffx = fclip(fx, one) ;
     ffy = fclip(fy, one) ;
     ffz = fclip(fz, one) ;

     rtx = *(mat+0) * ffx + *(mat+1) * ffy + *(mat+2) * ffz + *(mat+3) ;
     mat += 4 ;

     rty = *(mat+0) * ffx + *(mat+1) * ffy + *(mat+2) * ffz + *(mat+3) ;
     mat += 4 ;

     rtz = *(mat+0) * ffx + *(mat+1) * ffy + *(mat+2) * ffz + *(mat+3) ;

     *newfx = fclip(rtx, one) ;
     *newfy = fclip(rty, one) ;
     *newfz = fclip(rtz, one) ;
}

void	apply_symop_nowrap(real fx, real fy, real fz, real *mat, 
		           real *newfx, real *newfy, real *newfz)

	/**********************************************************
	fx, fy, fz are input fractional coords, not necessarily in
		range [0,1) 
	*mat is	rotation/translation matrix
	*newfx, *newfy, *newfz are output rotated and translated 
		fractional coordinates, in range [0,1]	
	**********************************************************/
{
     real	rtx, rty, rtz ;

     rtx = *(mat+0) * fx + *(mat+1) * fy + *(mat+2) * fz + *(mat+3) ;
     mat += 4 ;

     rty = *(mat+0) * fx + *(mat+1) * fy + *(mat+2) * fz + *(mat+3) ;
     mat += 4 ;

     rtz = *(mat+0) * fx + *(mat+1) * fy + *(mat+2) * fz + *(mat+3) ;

     *newfx = rtx ;
     *newfy = rty ;
     *newfz = rtz ;
}

void	apply_symop_hkl(int h, int k, int l, real *mat, 
			int *newh, int *newk, int *newl, real *off)

	/**********************************************************
	h, k, l are input reciprocal space indices	
	*mat is	rotation/translation matrix
	*newh, *newk, *newl are (hkl) after rotation/translation
	*off is offset to apply to imaginary part of str. factors 
	**********************************************************/
{
     real	newfh, newfk, newfl, foff ;

     newfh = *(mat+0) * h + *(mat+4) * k + *(mat+ 8) * l ;
     newfk = *(mat+1) * h + *(mat+5) * k + *(mat+ 9) * l ;
     newfl = *(mat+2) * h + *(mat+6) * k + *(mat+10) * l ;
     foff =  *(mat+3) * h + *(mat+7) * k + *(mat+11) * l ;

     *newh = newfh ;
     *newk = newfk ;
     *newl = newfl ;
     *off  = TWOPI * foff ; 
}

	/**********************************************************
	Set up the index array identifying equivalent positions 
	in asymmetric units. 
	**********************************************************/

void	setup_csind()
{
     int	i ;		/* counter of RT's */
     int	n ;		/* counter of grid points */
     real	fx, fy, fz ;	/* fractional coordinates before RT */
     real	newfx, newfy, newfz ;	/* fractional coordinates after RT */
     int	new_pt ;		/* composite index after RT */
     real	*mat ;
     int	lind ;		/* length of index_array */

      /**********************************************************************
      We set up a mask identifying one asymmetric unit in real space.
      Method: set all points to 0;
      loop over Nptotal; for each pt whose value is 0:
	 set it to 1;
	 find symmetry-related points; for each:
	    if it's 0 set it to 2
	    otherwise leave it alone.

      Then reset all 2's to 0.
      **********************************************************************/

     strcpy(message, "crysutil") ;

     unique_pmask = (char *) e_malloc(Nptotal*sizeof(char), message) ;
     equiv_list = (int *) e_malloc(Ncs*sizeof(int), message) ;

     for (n = 0; n < Nptotal; n++)
	  *(unique_pmask+n) = 0 ;

     for (n = 0; n < Nptotal; n++) {

	  if (*(unique_pmask+n) == 0) {			/* not yet set */

	       *(unique_pmask+n) = 1 ;

               index_to_fract(n, Nx, Ny, Nz, &fx, &fy, &fz) ;

	       for (i = 0, mat = matop; i < Ncs; i++, mat += MEL) {

                    apply_symop(fx, fy, fz, mat, &newfx, &newfy, &newfz) ;

                    new_pt = fract_to_index(n, newfx, newfy, newfz, Nx, Ny, Nz) ;
		    if (*(unique_pmask+new_pt) == 0)
			 *(unique_pmask + new_pt) = 2 ;
               }
          }
     }

     for (n = 0; n < Nptotal; n++) 
	  if (*(unique_pmask+n) == 2) 			/* related point */
	       *(unique_pmask+n) = 0 ;

     for (n = 0, lind = 0; n < Nptotal; n++) 
	  if (*(unique_pmask+n) == 1) 			/* count them */
	       lind++ ;

     Ics = lind ;

     return ;
}


void	get_equiv_list(int n)
	/* n is index into Nptotal for which equiv_list is req'd */
{
     int	i ;		
     real	fx, fy, fz ;	/* fractional coordinates before RT */
     real	newfx, newfy, newfz ;	/* fractional coordinates after RT */
     real	*mat ;

     index_to_fract(n, Nx, Ny, Nz, &fx, &fy, &fz) ;

     for (i = 0, mat = matop; i < Ncs; i++, mat += MEL) {

          apply_symop(fx, fy, fz, mat, &newfx, &newfy, &newfz) ;

	  *(equiv_list + i) = fract_to_index(n, newfx, newfy, newfz, Nx, Ny, Nz) ;

     }
}

	/**********************************************************
	Compute the average among equivalent cs positions in 
	asymmetric units, putting results into array 'average'  
	**********************************************************/

void	do_cs_symmetrization(real *array, real *average)
{
     real	one_ave ;
     int	i, n ;

     if (Ncs == 1) 

          for (n = 0; n<Nptotal; n++) 
	       *(average + n) = *(array + n)  ;

     else

          for (n = 0; n<Nptotal; n++) 

     	       if (*(unique_pmask+n)) {

	            one_ave = 0  ;

	            get_equiv_list(n) ;

	            for (i = 0; i < Ncs; i++)
	                 one_ave += *(array + *(equiv_list + i)) ;

                    one_ave /= Ncs ;

	            for (i = 0; i < Ncs; i++)
	                 *(average + *(equiv_list + i)) = one_ave ;
               }
}

	/**********************************************************
	Compute the average among equivalent cs positions in 
	asymmetric units, replacing original data by averages.
	Report symmetry violations.
	**********************************************************/

void	do_inplace_cs_symmetrization(real *array, real *sarray)
{
     real	one_ave ;
     int	i, n ;
     int	out = 0 ;
     int	bad ;
     int	x, y, z;
     real	sumsqdiff = 0 ;
     real	sumsqsum = 0 ;
     real	sum, diff ;

     if (Ncs == 1)
	  return ;

     for (n = 0; n<Nptotal; n++) {

	  if (*(unique_pmask+n)) {

	       one_ave = 0  ;

	       get_equiv_list(n) ;

	       for (i = 0; i < Ncs; i++)
	            one_ave += *(array + *(equiv_list + i)) ;

               one_ave /= Ncs ;

	       for (i = 0; i < Ncs; i++) {
	            diff = *(array + *(equiv_list + i)) - one_ave ;
	            sumsqdiff += diff * diff ;
	            sum = one_ave + *(sarray+*(equiv_list + i)) ;
	            sumsqsum += sum * sum ;
               }

               if (fabs(one_ave) > ABSSYM) {

	           for (i = 0, bad = FALSE; i < Ncs; i++) {
	                if ((fabs(*(array+*(equiv_list + i)) - one_ave) > 
	            ESYM*fabs(one_ave + *(sarray+*(equiv_list + i)))) &&
	            (ESYM*fabs(one_ave + *(sarray+*(equiv_list + i)))) > ABSSYM)

			     bad = TRUE ;
                   }
                   if (bad) {
	                out++ ;
                        if (very_verbose) {
	                    for (i = 0; i < Ncs; i++) {

			         x = *(equiv_list + i) % Nx ;
			         y = (*(equiv_list + i)/Nx) % Ny ;
			         z = *(equiv_list + i) / (Nx*Ny) ;

                                 fprintf(fp_log, "val(%d,%d,%d) = %6g, ", 
			         x, y, z, *(array+*(equiv_list + i))) ; 
		            }	    
                            fprintf(fp_log,  "and average %g\n\n", one_ave) ; 
                        }
                   }
               }

	       /* Now replace individual values by their average */

	       for (i = 0; i < Ncs; i++) 
	            *(array + *(equiv_list + i)) = one_ave ;
          }
     }

     sprintf(message, 
    "The rms fractional distance between original and symmetrized arrays is %g",
	sqrt(sumsqdiff / sumsqsum)) ;
     printTwice(message) ;
    
     if (out > 0) {
	  sprintf(message,
  "Symmetrization changed %d out of %d asym. unit elements ", out, Ics) ; 
          printTwice(message) ;
	  sprintf(message,
  "\tby more than %6.3f %% of the average.", 100*ESYM) ; 
          printTwice(message) ;
     }
     return ;
}


	/**********************************************************
	Expand 1 cs in any arbitrary position to full unit cell.
	**********************************************************/

void	cs_expand(real *array)
{
     int	i ;		/* counter of RT's */
     int	n ;		/* counter of grid points */
     real	fx, fy, fz ;	/* fractional coordinates before RT */
     real	newfx, newfy, newfz ;	/* fractional coordinates after RT */
     real	val ;
     int	new_pt ;		/* composite index after RT */
     int	*touch ;		/* each pt holds its asym. unit # */
     real	*mat ;


     touch = (int *) e_malloc(Nptotal*sizeof(int), "crysutil") ;

     for (n = 0; n < Nptotal; n++)
	  *(touch+n) = 0 ;

     for (n = 0; n < Nptotal; n++) {

          if (((val = *(array+n)) > 0) && (*(touch+n) == 0)) {

               index_to_fract(n, Nx, Ny, Nz, &fx, &fy, &fz) ;

	       for (i = 0, mat = matop; i < Ncs; i++, mat += MEL) {

                         apply_symop(fx, fy, fz, mat, &newfx, &newfy, &newfz) ;

                         new_pt = 
			    fract_to_index(n, newfx, newfy, newfz, Nx, Ny, Nz) ;
		         *(array + new_pt) = val ;
			 *(touch + new_pt) = i ;
               }
          }
     }

     free (touch) ;
}

void	setupNusq()
/*******************************************************************************

We use the general expression (Giacovazzo, Fundamentals of Crystallography, p66)
to set up the coefficients in the most general (triclinic) expression of 1/dHsq
which we write in the form:

   cohh*h*h + cokk*k*k + coll*l*l + 2*cokl*k*l + 2*colh(l*h + 2*cohk*h*k

*******************************************************************************/
{
     real	denom ;
     real	cosasq ;
     real	cosbsq ;
     real	coscsq ;

     cosalpha = cos(angles[0]*DTOR) ;
     cosbeta  = cos(angles[1]*DTOR) ;
     cosgamma = cos(angles[2]*DTOR) ;	/* cosines of angles */

     /* Clean up roundoff */

     if (fabs(cosalpha) < COSEPS) 
	  cosalpha = 0 ;
     if (fabs(cosbeta) < COSEPS) 
	  cosbeta = 0 ;
     if (fabs(cosgamma) < COSEPS) 
	  cosgamma = 0 ;

     cosasq = cosalpha*cosalpha ;
     cosbsq = cosbeta*cosbeta ;
     coscsq = cosgamma*cosgamma ;

     denom = 1. / (1 - cosasq - cosbsq - coscsq + 2.*cosalpha*cosbeta*cosgamma) ;
     cohh = denom / (aAxis*aAxis) * (1 - cosasq) ;
     cokk = denom / (bAxis*bAxis) * (1 - cosbsq) ;
     coll = denom / (cAxis*cAxis) * (1 - coscsq) ;
     cokl = denom / (bAxis*cAxis) * (cosbeta *cosgamma - cosalpha) ;
     colh = denom / (cAxis*aAxis) * (cosgamma*cosalpha - cosbeta) ;
     cohk = denom / (aAxis*bAxis) * (cosalpha*cosbeta  - cosgamma) ;

     if (fabs(cokl) < COSEPS)
	  cokl = 0 ;
     if (fabs(colh) < COSEPS)
	  colh = 0 ;
     if (fabs(cohk) < COSEPS)
	  cohk = 0 ;

     if (very_verbose) {
          fprintf(fp_log, "Nusq coefficients are:\n") ; 
          fprintf(fp_log, 
          "    cohh, cokk, coll = %e %e %e\nand  cokl, colh, cohk = %e %e %e\n", 
          cohh, cokk, coll, cokl, colh, cohk) ;
     }
}
/*******************************************************************************
nusq is the absolute magnitude squared of the reciprocal space vector (h,k,l)
in units of (1/Angstrom) squared.  It reduces to
	nusq = (h/a)^2 + (k/b)^2 + (l/c)^2 in orthogonal cells.
******************************************************************************/
real	nusq(int h, int k, int l)
{
     return(cohh*h*h + cokk*k*k + coll*l*l + 
	    2*cokl*k*l + 2*colh*l*h + 2*cohk*h*k) ;
}


/*******************************************************************************

	Derive transformation matrices between crystallographic-fractional
	and orthogonal coordinates.

F transforms from orthogonal to crystallographic-fractional;
Finv transforms from crystallographic-fractional to orthogonal.

      T	             -1 T
F =  M  and Finv = (M  ) (Giacovazzo, Fundamentals of Crystallography, p68)

******************************************************************************/

void prepare_F_and_Finv()
{
     real	volfac() ;
     real	singamma ;		/* sine of gamma (only one used) */
     real	D ;

     cosalpha = cos(angles[0]*DTOR) ;
     cosbeta  = cos(angles[1]*DTOR) ;
     cosgamma = cos(angles[2]*DTOR) ;	/* cosines of angles */

     singamma = sin(angles[2]* DTOR) ;

     D = 1. / volfac() ;

     Ftr[0][0] = F[0][0] = 1 / aAxis ;
     Ftr[1][0] = F[0][1] = -cosgamma / (aAxis*singamma) ; 
     Ftr[2][0] = F[0][2] = 
                    D * (cosalpha * cosgamma - cosbeta) / (aAxis*singamma) ;
     Ftr[0][1] = F[1][0] = 0 ;
     Ftr[1][1] = F[1][1] = 1 / (bAxis * singamma) ; 
     Ftr[2][1] = F[1][2] = 
                    D * (cosbeta * cosgamma - cosalpha) / (bAxis * singamma) ;
     Ftr[0][2] = F[2][0] = 0 ;
     Ftr[1][2] = F[2][1] = 0 ;
     Ftr[2][2] = F[2][2] = D * singamma / cAxis ;

     Finv[0][0] = aAxis ;
     Finv[0][1] = bAxis * cosgamma ;
     Finv[0][2] = cAxis * cosbeta ;
     Finv[1][0] = 0 ;
     Finv[1][1] = bAxis * singamma ;
     Finv[1][2] = cAxis * (cosalpha - cosbeta *cosgamma) / singamma ;
     Finv[2][0] = 0 ;
     Finv[2][1] = 0 ;
     Finv[2][2] = cAxis / (D * singamma) ;

}
real     volfac()
{
     return (sqrt(1. - cosalpha*cosalpha - cosbeta*cosbeta - cosgamma*cosgamma
		    + 2. * cosalpha * cosbeta * cosgamma)) ; 
}

real	Rsq(fx, fy, fz)
real	fx, fy, fz ;
{
     return(fx*fx*aAxis*aAxis + fy*fy*bAxis*bAxis + fz*fz*cAxis*cAxis + 
          2*fx*fy*aAxis*bAxis*cosgamma + 
	  2*fx*fz*aAxis*cAxis*cosbeta + 
	  2*fy*fz*bAxis*cAxis*cosalpha) ;
}

/*******************************************************************************
	Transformations between fractional and orthogonal coordinates.
	and from an index into the unit cell array to fractional coordinates.
*******************************************************************************/

void	ortho_to_fract(real xin, real yin, real zin, 
		      real *xout, real *yout, real *zout)
{
     void	matvecmult() ;
     real	vin[3] ;
     real	vout[3] ;

     vin[0] = xin ;
     vin[1] = yin ;
     vin[2] = zin ;

     matvecmult(F, vin, vout) ; 

     *xout = vout[0] ;
     *yout = vout[1] ;
     *zout = vout[2] ;
}
void	Tr_ortho_to_fract(real xin, real yin, real zin, 
		      real *xout, real *yout, real *zout)
{
     void	matvecmult() ;
     real	vin[3] ;
     real	vout[3] ;

     vin[0] = xin ;
     vin[1] = yin ;
     vin[2] = zin ;

     matvecmult(Ftr, vin, vout) ; 

     *xout = vout[0] ;
     *yout = vout[1] ;
     *zout = vout[2] ;
}
void	fract_to_ortho(real xin, real yin, real zin, 
		       real *xout, real *yout, real *zout)
{
     void	matvecmult() ;
     real	vin[3] ;
     real	vout[3] ;

     vin[0] = xin ;
     vin[1] = yin ;
     vin[2] = zin ;

     matvecmult(Finv, vin, vout) ; 

     *xout = vout[0] ;
     *yout = vout[1] ;
     *zout = vout[2] ;
}

void	index_to_fract(int n, int Nx, int Ny, int Nz, 
		       real *fx, real *fy, real *fz) 
{
     int	i, j, k ;

     i = n % Nx ;
     j = (n/Nx) % Ny ;
     k = n / (Nx*Ny) ;

     if (k < Nz) {
           *fx = i / (float) Nx ;
           *fy = j / (float) Ny ;
           *fz = k / (float) Nz ;
     }
     else{
           k -= Nz ;
           *fx = (i + 0.5) / (float) Nx ;
           *fy = (j + 0.5) / (float) Ny ;
           *fz = (k + 0.5) / (float) Nz ;
     }
}
#define ROUND 1.e-4
int	fract_to_index(int nbase, real fx, real fy, real fz, 
			int Nx, int Ny, int Nz)

	/**********************************************************
	nbase is where we started from: determines which grid 
	fx, fy, fz are fractional indices (input) 
	Nx, Ny, Nz are usual # of grid points in each dimension.
	**********************************************************/
{
     int	i, j, k ;
     int	n ;

      i = fx * Nx + ROUND ;
      j = fy * Ny + ROUND ;
      k = fz * Nz + ROUND ;
      n = i + j * Nx + k * Nx * Ny ;

      if (nbase >= Nx*Ny*Nz)
          n += Nx*Ny*Nz ;

      return (n) ;
}

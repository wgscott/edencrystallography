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

                                QSHARE.C

  Title:        Shared utilities in physical space
  Author:       Hanna Szoke
  Date:         Dec. 22, 1993 
  Function:     Worker functions for Solve and Back mainly.

void allocNpArrays() 	- solve		
void  initNparrays() 	- solve
real inv_mean_wt()	- qshare
void updateKnownValues()- solve		
int  checkSolutions()	- solve
void writeSol() 	- solve
void analyzEpixFile()	- back, maketar, qshare
void readEpixFile() 	- addmaps, convert, distance, forth, maketar, qshare, 
			  regrid, variance
void writeEpixFile()	- addmaps, back, convert, maketar, qshare, variance
void analyzeTar()	- qshare
void get_min_max()	- qshare, convert, regutil

******************************************************************************/

#include "eden.h"
#include "eden_Np.h"

	/* declarations of eden_Np.h variables */

real	*nump ;                 /* electrons/voxel, one iteration */
real	*snump ;		/* electrons/voxel, all iterations */
real  *totp ;              	/* running total electrons/voxel   */
real  *knownp ;               /* electrons/voxel corresponding to fc file */
real  *target ;               /* electrons/voxel corresponding to target(s) */
real  *weightp ;              /* weights for Np-space cost function */
real  *extra ;        	/* extra work array for CS constraint */

static	void	analyzeTar(real *, real *, int) ;	
static	real	inv_mean_wt(real *) ;
static  char	solution_filename[MAXSTRING] ;

		/************************************************************
               	Allocate & initialize arrays in Np space needed by solve
		(excluding those needed only by conjugate gradient solver
		and those in Npextended).
		************************************************************/
void allocNpArrays() 
{
     int	NN = 0 ;
     int	k ;
     int	there_is_phase_ext = FALSE ;

     strcpy(message, "allocNpArrays") ;

     knownp = (real *) e_malloc(Nptotal*sizeof(real), message) ;

     /* Allocate arrays for constraints */

     for (k = 0; k < Nconstraints; k++) 
	  if (con_type[k] == PHASE_EXT) {
	       there_is_phase_ext = TRUE ;
	       break ;
          }

     NN = Ntargets*Nptotal ;
     if (there_is_phase_ext)
	  NN += Nptotal ;

     if (NN > 0) {
          target = (real *) e_malloc(NN*sizeof(real), message) ;
          weightp= (real *) e_malloc(NN*sizeof(real), message) ;
     }

     for (k = 0; k < Nconstraints; k++) 
	  if (con_type[k] == CS) {
               extra = (real *) e_malloc(Nptotal*sizeof(real), message) ;
	       break ;
          }
          
     return ;
}

		/************************************************************
   		Decide whether size of solutions warrants a continuation of 
		the program.
		************************************************************/

#define EPS_GO_ON     1.e-3           /* Criterion for convergence */

int checkSolutions()
{
     int        n ;
     real     snumpmax ;
     real     numpmax ;

       /************************************************************
       Find largest solution in the cumulative array. 
       ************************************************************/

     snumpmax = fabs (*snump) ;

     for (n = 1; n < Nptotal; n++)
          if ( fabs (*(snump+n)) >  snumpmax) {
               snumpmax = fabs (*(snump+n)) ;
          }

       /************************************************************
       Find largest solution in the current iteration.  
       ************************************************************/

     numpmax = fabs (*nump) ;

     for (n = 1; n < Nptotal; n++)
          if ( fabs (*(nump+n)) >  numpmax) 
               numpmax = fabs (*(nump+n)) ;

       /************************************************************
       If it's too small, quit (unlikely in real cases).        
       ************************************************************/

     if (numpmax < snumpmax * EPS_GO_ON) {
          sprintf(message, 
		 "Solutions converged to within %g of maximum electrons/voxel", 
		 EPS_GO_ON) ;
          printTwice(message) ;
          return (STOP) ;
     }
     else
          return(GO_ON) ;
}

		/************************************************************
		Initialize Np space arrays, reading files where necessary.
		initMinMaxType() & alloc_ccg_mem() are called by Back as well.
		************************************************************/
		
void	initNparrays() 
{
     void	setup_csind() ;
     void	initMinMaxType() ;
     void	reviseMinMaxp() ;
     void	readEpixFile(char *, real *) ; 
     void	prepare_relative_index_list() ;	
     void	prepare_off_array() ;
     void	check_off_array() ;
     void	alloc_ccg_mem() ;		/* in ccg.c */
     int	hrInitialize() ;		/* in highres.c */
     void	hrSetup() ;			/* in highres.c */
     void	hrMaskOut() ;			/* in highres.c */

     real	volcoef ;	/*  V / (2pi*eta*dr^2)^(3/2)  */
     real	*big ;

     int	n ;
     int	k, NN, maxNhr ;


     for (n = 0; n < Nptotal; n++) 
	  *(knownp+n) = 0 ;

       /************************************************************
       Read model; in correction mode,                        
       revise the lower (and upper) bounds on the optimizer.
       Special treatment of "empty" is not really needed, since
       readEpixFile() deals with it...
	Note that hrInitialize(knownp) DISCARDS all the high points
	in knownp by masking them out!!!
       ************************************************************/

     if (strcmp (md_filename, "empty") != 0) {
          printTwice("Reading model electron map ...") ;
          readEpixFile(md_filename, knownp) ;
     }
     else
	  printTwice("Model electron map file is empty") ;

     setup_csind() ;
     sprintf(message, "initNparrays") ;

     if (HighRes)  {

       /************************************************************
        At this point, we re-allocate the maximum possible size for
        knownp and transfer Nptotal to this "big" array; then we find
        the actual req'd size (Npextended) and eventaully reallocate
        and transfer data back to knownp.
	Maybe use memcpy(big, knownp, Nptotal*sizeof(real))?
       ************************************************************/

          maxNhr = hrInitialize(knownp) ;	

	  big = (real *) e_malloc((Nptotal + maxNhr)*sizeof(real), message);
	 
          for (n = 0; n < Nptotal; n++) 
	       *(big+n) = *(knownp+n) ;

          e_free(Nptotal*sizeof(real), knownp) ;

	  hrSetup(big) ;		/* Npextended is set here */

	  knownp = (real *) e_malloc(Npextended * sizeof(real), message);
	 
          for (n = 0; n < Npextended; n++) 
	       *(knownp+n) = *(big+n) ;

          e_free((Nptotal+maxNhr)*sizeof(real),big) ;

          hrMaskOut(knownp) ;
     }
     else
          Npextended = Nptotal ;
 
     sprintf(message, "Nptotal = %d, Npextended = %d", Nptotal, Npextended) ;
     printTwice(message) ;

     nump = (real *) e_malloc(Npextended*sizeof(real), message) ;
     snump  = (real *) e_malloc(Npextended*sizeof(real), message) ;
     totp   = (real *) e_malloc(Npextended*sizeof(real), message) ;

     /*  TEMPORARILY - use nump for knownp(Nptotal) + hrKnownp */

     for (n = 0; n < Npextended; n++) 
	  *(nump+n) = *(knownp+n) ;

     initMinMaxType(Npextended) ;

     if (strncmp(mode_info, "corr", 4) == 0) 
          reviseMinMaxp(Npextended, nump) ;

     alloc_ccg_mem(Npextended) ;

     for (n = 0; n < Npextended; n++) {
	  *(nump+n) = 0 ;
	  *(snump+n) = 0 ;
     }
       /************************************************************
       Calculate a generic weight-based coefficient, volcoef. 
       ************************************************************/

     volcoef = crys_vol / pow((TWOPI*eta*dr*dr), (real) 1.5) ;
     sprintf(message, 
     "The volume coefficient: V/(2pi*eta*drsq)^3/2 = %g", volcoef) ;
     printTwice(message) ;

     for (k = 0, NN = 0; k < Nconstraints; k++)  {

       /************************************************************
       If there are targets, including phase extension targets,
       read them and their associated weights;  for ordinary
       targets only, factor the inverse-weight-based coefficient 
       into con_coef[k], 

       Note: you may invoke "full" as a keyword to use default 
       weights (1 everywhere), rather than reading in a weight file.

       Report average el/voxel over target. 
       ************************************************************/

          if ((con_type[k] == TARGET) || (con_type[k] == PHASE_EXT)) {

               for (n = 0; n < Nptotal; n++) 
	            *(weightp + NN + n) = 1 ;
	  
               sprintf(message,
		    "\nReading Np target info for constraint # %d...", k+1) ;
               printTwice(message) ;

               readEpixFile(ta_filename[k], target + NN) ;

               if (strcmp (wt_filename[k], "full") != 0) 
	            readEpixFile(wt_filename[k], weightp + NN) ;

	       if (con_type[k] == TARGET) 
	            con_coef[k] *= volcoef * inv_mean_wt(weightp + NN) ;

               sprintf(message, "Initial analysis:") ; 
               printTwice(message) ;

	       analyzeTar(knownp, weightp + NN, Nptotal) ;
               NN += Nptotal ;
          }
     

       /************************************************************ 
       Allocate and prepare arrays for Sayre's equation.  
       ************************************************************/

          else if (con_type[k] == SAYRE) {
               prepare_relative_index_list( (real) 1.) ;	
	       prepare_off_array() ;
	       check_off_array() ;
          }
     }
}

real	inv_mean_wt(real *wtarray)
{
     int	n ;
     real	swsq = 0 ;

       /************************************************************
       We calculate  N / Sum(wt^2) ;
       ************************************************************/

     for (n = 0; n < Nptotal; n++, wtarray++)
	  swsq += *wtarray * *wtarray ;

     if (swsq == 0) 
	  EdenError("All weights are zero! - Quitting.") ;
     
     return (Nptotal / swsq) ; 

}

		/************************************************************
		Transfer the solution from this iteration (nump) to the 
		solution arrays. 
		************************************************************/

void updateKnownValues()	
{
     void	reviseMinMaxp() ;
     int        n ;
     real	sum = 0 ;	/* sum of snump	*/
     real	thisn ;
     real	variance = 0 ;

       /************************************************************
       Add solutions from this iteration to solution set and 
       subtract same amounts from the ccg limiting arrays.
       ************************************************************/

     reviseMinMaxp(Npextended, nump) ;

     for (n = 0; n < Npextended; n++) {
	thisn = *(nump + n) ;
     	*(snump+n) += thisn ;
        variance += thisn*thisn ;
     }
       /************************************************************
       Provide diagnostics                              	  
       ************************************************************/

     for (n = 0; n < Npextended; n++) 
          sum += *(snump + n) ;
          
     sprintf(message,
	    "Cumulative sum of recovered electrons is %g,", sum) ;
     printTwice(message) ;

     variance = sqrt(variance/Npextended) ;

     sprintf(message,
	    "Variance of electrons/voxel for this iteration is %g\n", variance) ;
     printTwice(message) ;
}
		/************************************************************
		Reinitialize nump; used to be in updateKnownValues(), but I
		found a need for nump between that call and the next outer
		iteration.
		************************************************************/

void reinitNump()	
{
     int        n ;

     for (n = 0; n < Npextended; n++) 
	*(nump +n) = 0 ;
}

		/************************************************************
	        Write out the current solution        
		************************************************************/

void writeSol(char *name) 
{
     void	analyzEpixFile(real *) ;
     void	writeEpixFile(char *, real *) ;
     real     sum = 0 ;		/* sum of model + recovered electrons */
     int	k, n, NN ;

     /************************************************************
     The solution consists of the cumulative electrons found. 
     The known starting model is added in.  
     ************************************************************/

     if (strncmp(mode_info, "comp", 4) == 0) {
	  printTwice("Recovered electrons:") ;    
          analyzEpixFile(snump) ;
     }
     for (n = 0; n < Npextended; n++)
          *(totp+n) = *(snump+n) + *(knownp+n) ;


     for (n = 0; n < Npextended; n++) 
          sum += *(totp + n) ;
          
     sprintf(message,
	  "\nTotal electrons (starting model plus recovered): %g\n", sum) ;
     printTwice(message) ;
     
     printTwice("\t Analysis of electron densities:") ;

     analyzEpixFile(totp) ;

     for (k = 0, NN = 0; k < Nconstraints; k++, NN += Nptotal)  {

          if ((con_type[k] == TARGET) || (con_type[k] == PHASE_EXT)) {

               sprintf(message, "\nAnalysis over target # %d:", k+1) ; /* 1-origin */
               printTwice(message) ;

	       analyzeTar(totp, weightp + NN, Nptotal) ;
          }

     }

     writeEpixFile(name, totp) ;
}

#define	ROUNDOFF	1.e-10	/* ... to remove clutter in reports */

void	analyzEpixFile(real *array)
{
     void	get_min_max(real *, int, real *, real *) ;
     int	n, j, jbunch, sj ;
     real	ar_min, ar_max ;
     real	decade[10], sdecade[10] ;
     real	value ;
     real	from, to, jfrom, jto ;

     for (n = 0; n < Nptotal; n++) 
	  if (fabs(*(array+n)) < ROUNDOFF)
	       *(array+n) = 0 ;
    
     get_min_max(array, Nptotal, &ar_min, &ar_max) ;
     
     /* Print message indicating range of solutions */

     sprintf(message,
             "Range is (%g, %g) el/voxel, (%g, %g) el/cubA.",
	     ar_min, ar_max, ar_min/cubA_to_vox, ar_max/cubA_to_vox) ; 
     printTwice(message) ;

     /* Now for distribution of data ... */

     for (j = 0; j < 10; j++)
	  decade[j] = 0 ;

     for (n = 0; n < Nptotal; n++) {
	  value = *(array+n) / cubA_to_vox ;
          j = (int) (value * 10.) ;
	  if (j > 9) j = 9 ;
	  if (j < 0) j = 0 ;
	  decade[j] += 1 ;
     }
     for (j = 0; j < 10; j++)
	  decade[j]  *= 100. / Nptotal ;

     /* Determine degree of bunching */

     for (jbunch = 10, j = 0; j < 10; j++)
          if (decade[j] > 50.)
               jbunch = j ;

     /* Print message indicating distribution of unbunched data */

     if (jbunch == 10) {
          printTwice("\nDistribution of electron densities:") ;
          printTwice("Range (el/cubA)     %\n") ;

          sprintf(message, "    < 0.1:     %#6.2f", decade[0]) ;
          printTwice(message) ;

          from = 0.1 ;
          for (j = 1; j < 9; j++) {
	       to = (j+1) / 10. ;
               sprintf(message, "%3.1f - %3.1f:     %#6.2f", 
	            from, to, decade[j]) ;
	       printTwice(message) ;
	       from = to ;
          }
          sprintf(message, "    > 0.9:     %#6.2f", decade[9]) ;
          printTwice(message) ;
     }

     /* Handle data with bunching.  */

     else {
          for (j = 0; j < 10; j++)
	       sdecade[j] = 0 ;

          for (n = 0; n < Nptotal; n++) {
	       value = *(array+n) / cubA_to_vox ;
               j = (int) (value * 10.) ;
	       if (j > 9) j = 9 ;
	       if (j < 0) j = 0 ;

               if (j == jbunch) {
                    sj = (int) (value*100 - 10*jbunch) ;
                    if (sj > 9) sj = 9 ;
                    if (sj < 0) sj = 0 ;
	            sdecade[sj] += 1 ;
               }
          }
          for (sj = 0; sj < 10; sj++)
	       sdecade[sj]  *= 100. / Nptotal ;

          printTwice(
            "Range (el/cubA)     %    Subrange (el/cubA)     %\n") ;

          jfrom = jbunch / 10. ; 
          jto = jfrom + 0.01 ;

          sprintf(message, 
            "    < 0.1:     %#6.2f     %4.2f - %4.2f:     %#6.2f", 
                  decade[0], jfrom, jto, sdecade[0]) ;
          printTwice(message) ;

          from = 1 / 10. ;
          jfrom = jto;

          for (j = 1; j < 9; j++) {
	       to = (j + 1) / 10. ;
               jto = jbunch/10. + (j + 1) / 100. ;
               sprintf(message, 
             "%3.1f - %3.1f:     %#6.2f     %4.2f - %4.2f:     %#6.2f", 
	            from, to, decade[j], jfrom, jto, sdecade[j]) ;
	       printTwice(message) ;
	       from = to ;
	       jfrom = jto ;
          }

          sprintf(message, 
                "    > 0.9:     %#6.2f     %4.2f - %4.2f:     %#6.2f", 
                  decade[9], jfrom, jfrom + 0.01, sdecade[9]) ;
          printTwice(message) ;
     }
}

		/************************************************************
                Read in / write out binary file.
		************************************************************/

void	readEpixFile(char *name, real *array) 
{
     FILE 	*fp;
     float	*f_array ;	/* for converting real to float */
     int	n ;
     int	Nread ;
     char	*ext = ".bin" ;

     for (n = 0; n < Nptotal; n++) 
          *(array + n) = 0 ;

     if (strcmp (name, "empty") == 0) 
	  return ; 

     f_array = (float *) e_malloc(Nptotal*sizeof(float), "readEpixFile") ;

     sprintf(solution_filename, "%s", name) ;

     if (strstr(solution_filename, ext) == NULL)
	  strcat(solution_filename, ext) ;

     if ((fp = fopen(solution_filename, "rb")) == NULL)
     {
          sprintf(message, "readEpixFile: Cannot open %s", solution_filename); 
	  EdenError(message) ;
     }
     Nread = fread(f_array, sizeof(float), Nptotal, fp);
     if (Nread != Nptotal) {
          sprintf(message, 
     "Expected length : %d, actual length: %d", Nptotal, Nread) ;
          printTwice(message) ;
          sprintf(message, 
     "readEpixFile: wrong number of items in file %s", solution_filename); 
	  EdenError(message) ;
     }

/*** One could also check for a file that's too long: 
     if (n = filelength(solution_filename) > Nptotal) {
          sprintf(message,
     "Expected length : %d, actual length: %d", Nptotal, n) ;
          printTwice(message) ;
          sprintf(message,
     "readEpixFile: file %s is longer than expected", solution_filename); 
          printTwice(message) ;
          printTwice("Was it prepared at a higher resolution?") ;
	  EdenError("Quitting") ; ***/
        
     fclose(fp);

     for (n = 0; n < Nptotal; n++) 
          *(array + n) = *(f_array + n) ;

     e_free(Nptotal*sizeof(float), f_array) ;
}
void	writeEpixFile(char *name, real *array) 
{
     FILE 	*fp;
     float	*f_array ;	/* for converting real to float */
     int	n ;
     char	*ext = ".bin" ;

     f_array = (float *) e_malloc(Nptotal*sizeof(float), "writeEpixFile") ;

     for (n = 0; n < Nptotal; n++) 
          *(f_array + n) = *(array + n) ;

     if (strstr(name, ext) == NULL)
	  strcat(name, ext) ;

     sprintf(solution_filename, "%s", name) ;

     if ((fp = fopen(solution_filename, "wb")) == NULL)
     {
          sprintf(message, "writeEpixFile: Cannot open %s", solution_filename); 
	  EdenError(message) ;
     }
     if (fwrite(f_array, sizeof(float), Nptotal, fp) != Nptotal) {
          sprintf(message, "writeEpixFile: Error writing %s", solution_filename); 
	  EdenError(message) ;
     }
     fclose(fp);

     e_free(Nptotal*sizeof(float), f_array) ;
}

void	analyzeTar(real *array, real *wt, int N)	
{
     int	n ;
     real	Sum = 0 ;
     real 	count = 0 ;
     char	messout[MAXSTRING] ;

     /* Find the average electron level (in el/voxel and el/cubA) over the 
     targetted region (for which *wt > 0).   */

     for (n = 0; n < N; n++, array++, wt++) {
	  Sum += *array * *wt ;
	  count += *wt ;
     }

     sprintf(messout, 
     "\tAverage over targetted region = %7.3f el/voxel. %7.3f el/cubA",
          Sum/count, Sum/(count*cubA_to_vox)) ;
     printTwice(messout) ;
     printTwice("that will be constrained to ~0.34 el/cubA") ;
}

void	get_min_max(real *array, int N, real *a_min, real *a_max)
{   
     int	n ;

     *a_min = 10000. ;
     *a_max = -10000. ;

     for (n = 0; n < N; n++) {
	  if (*(array+n) < *a_min)
	       *a_min = *(array+n) ;
	  if (*(array+n) > *a_max)
	       *a_max = *(array+n) ;
     }
}

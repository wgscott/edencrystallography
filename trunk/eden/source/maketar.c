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

                                MAKETAR

  Title:        MAKETAR - Utility to generate solvent or protein targets
  Author:       Hanna Szoke
  Date:         9/12/94
  Function:     This program takes as input a set of electron/voxel (binary)
 		files, on the basis of which it creates and writes out 2 sets 
		of electron/voxel files -- weight files and target files.

		It creates a set of weight files, either by setting low points 
		to 1 and high ones to zero, or by setting high points 
		to 1 and low ones to zero.  The choice is determined from 
		the input: "TARGET LOW" sets low points to 1 while 
		"TARGET HIGH" set high points to 1.  There is no default.
		
		The fraction of masked-in points is user-definable (e.g., 
		MASK_FRACTION 0.4) and defaults to 50% (0.5).  Alternately, 
		user may set a threshold (in el/cubic A) at which the mask 
		kicks in (e.g., THRESHOLD 0.2).
		
		Maketar also writes out a set of targets containing a 
		constant value wherever the mask is set and the input value
		elsewhere.  The constant value is user-definable and should 
		be entered in el/cubic A (e.g., TARGET_VAL  0.25); the default 
		is SOLVENT_DENSITY (0.34, defined in pdb.h).  
		The target will not be useful in all cases!  
		For example, if the input is a partial model and the intent 
		is to create a stabilizing target, the original input itself 
		should be used as the target.

*******************************************************************************/

#include "util.h"
#include "cellparams.h"
#include "dims.h"
#include "pdb.h"

static	int	use_mf ;	/* use mask_fraction or threshold? */
static	real	mask_fraction ;	/* fraction of points to be masked in */
static	real	threshold ;	/* el. dens. at which masking kicks in */
static	real	target_value ;
static	int	target_in ;	/* what gets targetted	- low or high values? */

enum	targ {NONE, LOW, HIGH} ;

static	void	findCutoff(real *) ;		/* find cutoff corresponding
						to input mask fraction */
static	void    m_readInput() ;		
static	void	makeTargetWeights(real *, real *) ;	/* "do it" */
void	hpsort(unsigned long, float []) ;
					   

void	maketar_main(int argc, char *argv[])
{
	/* non-local functions called */

     void	readEpixFile(char *, real *) ;
     void	analyzEpixFile(real *) ;

     char	filename1[MAXSTRING] ;

     real	*xnump ;	/* electron densities */
     real	*outp ;		/* weight and target arrays for output */

  /********************************************************
  Check for right ballpark, say hello, identify file names.
  ********************************************************/

     if (argc > optind+3) 	/* too many arguments */
	  ballpark(caller) ; 
     else if (argc == optind+3)
          sprintf(filename1, "%s", argv[optind+2]) ;
     else {
	  prompt("\nPlease enter the name of an electron/voxel file: ") ;
          sscanf(terminp, "%s", filename1) ;
	  sprintf(message, "Maketar will use el/voxel file %s", filename1) ;
	  printTwice(message) ;
     }

     hello(caller) ;

  /********************************************************
  Fetch input conditions.  Set up arrays in physical space.
  Read the binary model(s).  
  ********************************************************/

     m_readInput() ;

     if (HighRes) {
	  printTwice("\n\t\tNOTE:") ;
	  printTwice(
     "Maketar does not deal with input lists of high-resolution points.\n") ;
     }

     xnump = (real *) e_malloc(Nptotal*sizeof(real), caller) ;
     outp  = (real *) e_malloc(Nptotal*sizeof(real), caller) ;
          
     printTwice("Reading solution map ...") ;

     readEpixFile(filename1, xnump) ;
     analyzEpixFile(xnump) ;

  /********************************************************
  Compute and write weights and targets.
  ********************************************************/

     findCutoff(xnump) ;

     makeTargetWeights(xnump, outp) ;
}

	/*****************************************************
  	Find the threshold value corresponding to an input
	mask_fraction or find the mask_fraction corresponding 
	to an input threshold.
	******************************************************/

void	findCutoff(real *model) 
{
     float	*copy, *pcopy ;		/* array to use for sorting in place */
     int	n ;
     int	ncutoff ;
     unsigned long	Nun ;

     /* allocate space, copy model data. */

     copy  = (float *) e_malloc(Nptotal*sizeof(float), "maketar") ;

     for (n = 0, pcopy = copy; n < Nptotal; n++, model++, pcopy++)  
	  *pcopy = *model ;

     /*  Sort in ascending order */

     Nun = Nptotal ;
     hpsort(Nun, copy) ;

     if (use_mf) {

          /*  Use mask_fraction to determine threshold */

          if (target_in == HIGH) 
               ncutoff = (1. - mask_fraction) * Nptotal ;	
          else 
               ncutoff = mask_fraction * Nptotal ;	

          if (ncutoff < 1)
	       ncutoff = 1 ;
          else if (ncutoff > Nptotal)
               ncutoff = Nptotal ;

          threshold = *(copy + ncutoff - 1) ;

          sprintf(message, 
               "Requested fraction of targetted points is %g", mask_fraction) ;
          printTwice(message) ;
          sprintf(message, 
	       "corresponding to a threshold of %g el/vox (%g el/cubA).\n", 
	       threshold, threshold / cubA_to_vox) ;
          printTwice(message) ;
     }
     else {

          /*  Use threshold to determine mask_fraction */

	  n = 0 ;

	  if (target_in == HIGH) {
	       while ((n < Nptotal) && (*(copy+n) <= threshold))
                    n++ ;
               ncutoff = n ;
	       mask_fraction = (Nptotal - ncutoff) / (real) Nptotal ;
          }
          else {
	       while ((n < Nptotal) && (*(copy+n) <= threshold))
                    n++ ;
               ncutoff = n ;
	       mask_fraction = ncutoff / (real) Nptotal ;
          }
          sprintf(message, 
	       "Requested threshold is %g el/vox (%g el/cubA),", 
	       threshold, threshold/cubA_to_vox) ;
          printTwice(message) ;
          sprintf(message, 
       "corresponding to a fraction %g of targetted points.\n", mask_fraction) ;
          printTwice(message) ;
     }
     e_free(Nptotal*sizeof(float),copy) ;
}

	/*****************************************************
  	Calculate target weights and values.
	******************************************************/

void	makeTargetWeights(real *model, real *outp) 
{
     void	writeEpixFile(char *, real *) ;
     real 	*pmodel, *poutp ;
     real	avetg = 0, avenottg = 0 ;
     int 	n, Nweights = 0 ;
     char	solname[MAXSTRING] ;

     printTwice("Preparing weights ...\n") ;

     /*  Make weights array */

     if (target_in == HIGH) {
          for (n = 0, pmodel = model, poutp = outp; n < Nptotal; 
	       n++, pmodel++, poutp++)  
       
	       *poutp = (*pmodel >= threshold) ? 1:0 ;
     } 
     else {
          for (n = 0, pmodel = model, poutp = outp; n < Nptotal; 
	       n++, pmodel++, poutp++) 
       
	       *poutp = (*pmodel <= threshold) ? 1:0 ;
     } 
     sprintf(solname, "weight") ; 
     writeEpixFile(solname, outp) ;


     for (n = 0, Nweights = 0, poutp = outp; n < Nptotal; n++, poutp++)  
	  Nweights += *poutp ;

     sprintf(message, "Number of targetted points is %d out of %d.\n",
     Nweights, Nptotal) ;
     printTwice(message) ;

     /* Some statistics.... 	*/     

     if (target_in == HIGH) {
           for (n = 0, pmodel = model; n < Nptotal; n++, pmodel++) {
	        if (*pmodel >= threshold)
		     avetg += *pmodel ;
	        else
		     avenottg += *pmodel ;
           }  
     }
     else {
           for (n = 0, pmodel = model; n < Nptotal; n++, pmodel++) {
	        if (*pmodel <= threshold)
		     avetg += *pmodel ;
	        else
		     avenottg += *pmodel ;
           }  
     }
     if (Nweights > 0)
          avetg /= Nweights ;
     else
	  avetg = 0 ;

     if (Nweights < Nptotal)
          avenottg /= (Nptotal - Nweights) ;
     else
	  avenottg = 0 ;

     sprintf(message, 
     "Average of input over targetted area = %g el/vox (%g el/cub A).", 
	  avetg, avetg/cubA_to_vox) ;
     printTwice(message) ;

     sprintf(message, 
     "Average of input over non-targetted area = %g el/vox (%g el/cub A).\n", 
	  avenottg, avenottg/cubA_to_vox) ;
     printTwice(message) ;

     /*  Make target array; it's not clear what should be put into it.
     No one decision can produce a target that is applicable always. 
     Current procedure is: fill it with target_value, but leave the
     header information for reference..  */

     if (target_in == HIGH) {

          printTwice(
"Preparing target, containing input where value >= threshold, 0 elsewhere.") ;

          for (n = 0, pmodel = model, poutp = outp; n < Nptotal; 
	       n++, pmodel++, poutp++)  
       
	       *poutp *= *pmodel ;
     
     }
     else {

          sprintf(message, 
          "Preparing target, containing %g el/voxel (%g el/cubic A)\n", 
	       target_value*cubA_to_vox, target_value) ;
          printTwice(message) ;

          for (n = 0, poutp = outp; n < Nptotal; n++, poutp++)  
       
	       *poutp *= target_value*cubA_to_vox ;
     }
     sprintf(solname, "target") ; 
     writeEpixFile(solname, outp) ;
}

	/*****************************************************
	Read input, identifying data in terms of keywords.
	******************************************************/

void m_readInput()
{
     char	tval[MAXSTRING] ;
     int	j, k ;

     readBasicInput() ;

     mask_fraction = 0.0 ;
     threshold = UNSET ;
     target_in = NONE ;
     target_value = SOLVENT_DENSITY ;		/* Note - new default! */

     for (k = 0; k < Nlines; k++) {
          strcpy(nextline, allinp + k*MAXSTRING) ;
          sscanf(nextline, "%s", id) ;

          if (strcmp(id, "MASK_FRACTION") == 0) {
               sscanf(nextline, "%*s %g", &t1) ;
               mask_fraction = t1 ;
          }
          if (strcmp(id, "THRESHOLD") == 0) {
               sscanf(nextline, "%*s %g", &t1) ;
               threshold = t1 ;
          }
          if (strcmp(id, "TARGET_VALUE") == 0) {
               sscanf(nextline, "%*s %g", &t1) ;
               target_value = t1 ;
          }
          if (strcmp(id, "TARGET") == 0) {
               sscanf(nextline, "%*s %s", tval) ;

               for (j=0; j<(int)strlen(tval); j++)
                   tval[j] = toupper((int) tval[j]) ;

               if (strcmp(tval, "LOW") == 0)
		   target_in = LOW ;
               else if (strcmp(tval, "HIGH") == 0)
		   target_in = HIGH ;
          }
     }

     /* Checks */

     if (target_in == NONE) 
	 EdenError("Illegal or missing value for TARGET: use LOW or HIGH") ;
     
     if (mask_fraction == 0) {
          use_mf = FALSE ;
	  if (threshold == UNSET) {
	       threshold = 0 ;
	       mask_fraction = 0.5 ;
	       use_mf = TRUE ;
          }
	  else {
	       threshold *= cubA_to_vox ;		/* convert to el/voxel */
          }
     }
     else {
	  use_mf = TRUE ;
	  if (threshold != UNSET) 
	       EdenError("You can't have both THRESHOLD and MASK_FRACTION!") ;
     }

     /* Reports */

     if (target_in == HIGH) 
          fprintf(fp_log, "The target is high.\n") ;
     else 
          fprintf(fp_log, "The target is low.\n") ;
}


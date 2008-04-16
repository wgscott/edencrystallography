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

                                INIT.C

  Title:        General initialization package for EDEN's Solve.
  Author:       Hanna Szoke
  Date:         1/23/92
  Function:     This package contains set-up functions for Solve.
		These include: reading and interpreting .inp file,
		echoing what is in the input to the log file; and
                defining the run conditions for Solve.

		The basic .inp read function (used by ALL Eden's 
		pieces) is packaged in basics.c;
		Input of (hkl) data is in hklread.c.

******************************************************************************/

#include	"eden.h"
#include	"mir.h"
#include	"ccg.h"		/* needed for dfdxpc */

  /************************************************************
  declarations of some EXTERN variables in dims.h               
  ************************************************************/

real	nusqlim[NUMSEGS] ;	/* values for 1/8, 1/4,... of (hkl) shells */
real	fscale ;		/* scaling factor for fo data */
real	hrCutoff ;		/* hrCutoff * <vox> = high res limit	*/

  /************************************************************
  declarations of EXTERN variables in ccg.h               
  ************************************************************/

real	dfdxpc ;		/* frac. change, df/dx passed to getsol */

  /************************************************************
  declarations of EXTERN variables in eden.h               
  ************************************************************/

float	R_stop ;		/* R_factor at which code should stop	*/
real	discrp_frac ;           /* fraction of DISCRP set by user (def: 1) */

char	fo_filename[MAXSTRING] ;/* name of file, observed structure factors */
char	md_filename[MAXSTRING] ;/* name of file, model electron map */

int	Nconstraints ;		/* number of all constraints */
int	Ntargets ;		/* number of target constraints */
int	con_type[MAXCONSTR] ;	/* constraint type */
real	relwt_con[MAXCONSTR] ;	/* relative weight for a constraint */
real	cost_addend[MAXCONSTR] ;/* constant factor in the cost function */
real	con_coef[MAXCONSTR] ;	/* coefficient for cost function calc. */

char	target_type[MAXCONSTR][MAXSTRING] ;/* target type, target constraint */
char	ta_filename[MAXCONSTR][MAXSTRING] ;/* target file, target constraint */
char	wt_filename[MAXCONSTR][MAXSTRING] ;/* weight file, target constraint */

real	phase_ext_res ;		/* low resolution for phase extension */
int	robust ;		/* if TRUE use robust hkl cost function */
real	rob_factor = UNSET ;	/* factor to use with robust cost function */

char	detwin ;		/* is crystal detwinning requested? def: N */
char	t_type ;		/* A (amplitude) or I (intensity)	*/
int	*t_matrix ;		/* twinning matrix transformations */
real	t_frac ;		/* twinning fraction */

 /************************************************************
  declarations of EXTERN variables in util.h               
  ************************************************************/

char	sf_filename[MAXSTRING] ;	/* generic structure factor file name */
char	pwd_name[MAXSTRING] ;	/* name of str. factor file, w/o path */

  /************************************************************
  declarations of EXTERN variables in mir.h                
  NOTE: Solve will reset NM based on Nder.
  ************************************************************/

real	relwt_der[MAXDER] ;	/* relative weights for native & derivatives  */
real	res_der[MAXDER] ;	/* intrinsic res. for native & derivatives  */
int	useset[MAXDER] ;	/* which resolution set should be used? */
int	Nder;			/* number of derivative files  */
int	NM = 1 ;		/* Nder+1 (Native plus derivative files)  */
int	Nres =  1;		 /* number of resolution sets for expfac */
real	resrat[MAXRES] ;	/* resolution ratios for the resolution sets */

char	fc_heavy_fn[MAXDER][MAXSTRING] ;/* file names, heavy atom str factors */
char	fo_der_fn[MAXDER][MAXSTRING] ;	/* file names, deriv obs str factors */
real	fscale_der[MAXDER] ;	/* input scale factors for derivative fobs */
int	autoscale ;		/* if TRUE (default) apply rel_scale_m's */

real	fprime[MAXDER] ;	/* f' for anomalous MAD heavy file */
real	f2prime[MAXDER] ;	/* f'' for anomalous MAD heavy file */
real	Zheavy[MAXDER] ;	/* Z for heavy atom - only one per file! */
char	anom_der_flag[MAXDER] ;	/* flag identifying anomalous data */

  /************************************************************
  declarations of miscellaneous variables local to init.c  
  ************************************************************/

static	char	title[MAXSTRING] ;	/* whatever the user wants ... */

static	real	relwt_native ;	/* rel. weight for native (hkl) set */
static	real	res_native ;	/* intrinsic res. for native (hkl) set */

static	void	setGlobalDefaults() ;
static	void	readConstraintInfo() ;

				/*********************************************
               			Read Solve input 
				*********************************************/
void readInput()    
{
     int	fetchTFStatus(char *, int) ;
     int	k ;

     /************************************************************
     Read basics; set defaults for global variables.
     ************************************************************/

     readBasicInput() ;
     setGlobalDefaults() ;


     for (k = 0; k < Nlines; k++) {
          strcpy(nextline, allinp + k*MAXSTRING) ;
          sscanf(nextline, "%s", id) ;

          if (strcmp(id, "TITLE") == 0) {
               strcpy(title, nextline+6) ;
	       goodline[k] = TRUE ;
          }
          else if (strcmp(id, "DFDX_CRIT") == 0) {
               sscanf(nextline, "%*s %g", &t1) ;
               dfdxpc = t1 ;
	       goodline[k] = TRUE ;
          }

          else if ((strcmp(id, "FO_FILENAME") == 0) ||
              (strcmp(id, "FO_FN") == 0)) {
               sscanf(nextline, "%*s %s", fo_filename) ;
	       goodline[k] = TRUE ;
          }
          else if (strcmp(id, "R_STOP") == 0) {
               sscanf(nextline, "%*s %g", &t1) ;
	       R_stop = t1 ;
	       goodline[k] = TRUE ;
          }
	  else if (strcmp(id, "DISCRP_FRAC") == 0) {
	       sscanf(nextline, "%*s %g", &t1) ;
	       discrp_frac = t1 ;
	       goodline[k] = TRUE ;
	  }
          else if (strcmp(id, "RELWT_NATIVE") == 0) {
               sscanf(nextline, "%*s %g", &t1) ;
               relwt_native = t1 ;
	       goodline[k] = TRUE ;
          }
          else if (strcmp(id, "RES_NATIVE") == 0) {
               sscanf(nextline, "%*s %g", &t1) ;
               res_native = t1 ;
	       goodline[k] = TRUE ;
          }
          else if (strcmp(id, "FSCALE") == 0) {
               sscanf(nextline, "%*s %g", &t1) ;
               fscale = t1 ;
	       goodline[k] = TRUE ;
          }
          else if (strcmp(id, "ROB_FACTOR") == 0) {
               sscanf(nextline, "%*s %g", &t1) ;
               rob_factor = t1 ;
	       goodline[k] = TRUE ;
          }
          if ((strcmp(id, "MD_FILENAME") == 0) ||
              (strcmp(id, "MD_FN")       == 0) ||
              (strcmp(id, "EDENS_FN") == 0)) {
               sscanf(nextline, "%*s %s", md_filename) ;
	       goodline[k] = TRUE ;
          }
          else if (strcmp(id, "NCONSTRAINTS") == 0) {
	       sscanf(nextline, "%*s %d", &Nconstraints) ;
	       goodline[k] = TRUE ;
          }
          else if (strcmp(id, "T_TYPE") == 0) {
               sscanf(nextline, "%*s %s", message) ;
               if ((strncmp(message, "I", 1) == 0) ||
                   (strncmp(message, "i", 1) == 0)) 
                   t_type ='I' ;
               else if ((strncmp(message, "A", 1) == 0) ||
                        (strncmp(message, "a", 1) == 0)) 
                   t_type = 'A' ;
               else
                   EdenError("Illegal detwinning type!") ;

	       goodline[k] = TRUE ;
          }
          else if (strcmp(id, "T_FRAC") == 0) {
	       sscanf(nextline, "%*s %g", &t1) ;
               t_frac = t1 ;
	       goodline[k] = TRUE ;
          }
          else if (strcmp(id, "T_MATRIX") == 0) {
   
               t_matrix = (int *) e_malloc(9*sizeof(int), "init") ;
               sscanf(nextline, "%*s %d %d %d %d %d %d %d %d %d", 
                     (t_matrix+0), (t_matrix+1), (t_matrix+2),
                     (t_matrix+3), (t_matrix+4), (t_matrix+5),
                     (t_matrix+6), (t_matrix+7), (t_matrix+8)) ;
               goodline[k] = TRUE ;
          }
          else if ((HighRes) && (strcmp(id, "HRCUTOFF") == 0)) {
               sscanf(nextline, "%*s %g", &t1) ;
               hrCutoff = t1 ;
               goodline[k] = TRUE ;
          }
     }

     robust = fetchTFStatus("ROBUST", FALSE) ;
     if ((robust) && (rob_factor == UNSET))
         EdenError("Missing rob_factor!") ;

     autoscale = fetchTFStatus("AUTOSCALE", TRUE) ;

     if (strcmp(fo_filename, "none") == 0) 
         EdenError("Missing fo file name!") ;
          
     if (fscale == UNSET)
         EdenError
           ("Missing fscale, scale factor for fo; when in doubt, use 1.") ;

     if (strcmp(mode_info, "none") == 0) 
         strcpy(mode_info, "correction") ;
          
     if ((strncmp(mode_info, "corr", 4) != 0) &&
         (strncmp(mode_info, "comp", 4) != 0)) 
         EdenError("Illegal mode information!") ;
          
     if (strcmp(md_filename, "none") == 0) {
         EdenError("The starting binary file name is missing!");
     }
     if (res_native == UNSET)
         res_native = input_res ;

     detwin = fetchTFStatus("DETWIN", FALSE) ;

     if (detwin) {
         if (t_frac == UNSET)
             EdenError("Missing twinning fraction!") ;
         if (t_matrix == NULL) 
             EdenError("Missing twinning matrix!") ;
         if ((t_frac < 0.) || (t_frac > 0.5)) {
             sprintf(message, 
             "Twinning fraction, %g, is out of range!\n", t_frac) ;
             EdenError(message) ;
         }
     }

     /************************************************************
     Deal with constraints.
     Do NOT free allinp; it may be used by satellite programs.
     ************************************************************/

     readConstraintInfo() ;
}

				/*********************************************
               			Read Solve constraint information
				*********************************************/
void readConstraintInfo()    
{
     char	as_con_type[MAXSTRING][MAXCONSTR] ;
     int	j, k, m ;
     int	status = 0 ;

     for (k = 0; k < MAXCONSTR; k++) {
          strcpy(wt_filename[k], "none") ;
          strcpy(ta_filename[k], "none") ;
          strcpy(as_con_type[k], "") ;
	  relwt_con[k] = 0 ;
	  cost_addend[k] = 0 ;
	  con_coef[k] = 1 ;
     }
     phase_ext_res = 0 ;

     if (Nconstraints == 0)
	  return ;

     for (k = 0; k < Nlines; k++) {
          strcpy(nextline, allinp + k*MAXSTRING) ;
          sscanf(nextline, "%s", id) ;

	/* General constraint input */

          if ((strncmp(id, "RELWT_CON", 9) == 0) && isdigit(id[9]))  {
	       j = atoi(&id[9]) - 1 ;
               sscanf(nextline, "%*s %g", &t1) ;
	       relwt_con[j] = t1 ;
	       goodline[k] = TRUE ;
          }
          else if ((strncmp(id, "COST_ADDEND", 11) == 0) && isdigit(id[11]))  {
	       j = atoi(&id[11]) - 1 ;
               sscanf(nextline, "%*s %g", &t1) ;
	       cost_addend[j] = t1 ;
	       goodline[k] = TRUE ;
	  }
          else if ((strncmp(id, "CON_TYPE", 8) == 0) && isdigit(id[8]))  {
	       j = atoi(&id[8]) - 1 ;
               sscanf(nextline, "%*s %s", as_con_type[j]) ;
	       goodline[k] = TRUE ;
	  }

	/* Target and Phase Extension constraint input */

          else if (((strncmp(id, "TA_FILENAME", 11) == 0) && isdigit(id[11])) ||
                   ((strncmp(id, "TA_FN",        5) == 0) && isdigit(id[ 5]))) {
	       j = atoi(&id[11]) - 1 ;
               sscanf(nextline, "%*s %s", ta_filename[j]) ;
	       goodline[k] = TRUE ;
	  }
          else if (((strncmp(id, "WT_FILENAME", 11) == 0) && isdigit(id[11])) ||
                   ((strncmp(id, "WT_FN",        5) == 0) && isdigit(id[ 5]))) {
	       j = atoi(&id[11]) - 1 ;
               sscanf(nextline, "%*s %s", wt_filename[j]) ;
	       goodline[k] = TRUE ;
	  }

	/* Phase Extension constraint input */

          else if (strcmp(id, "PHASE_EXT_RES") == 0) {
               sscanf(nextline, "%*s %g", &t1) ;
	       phase_ext_res = t1 ;
	       goodline[k] = TRUE ;
          }
     }

     /************************************************************
     Check each constraint type for legality.
     ************************************************************/

     for (k = 0; k < Nconstraints; k++) 
          if (relwt_con[k] == 0) {
               sprintf(message, "Missing relative weight for target # %d!", k+1) ;
	       EdenWarning(message) ;
               status = -1 ; 
           }

     for (k = 0; k < Nconstraints; k++) {
	  
	  for (m=0; m<(int)strlen(as_con_type[k]); m++)
	       as_con_type[k][m] = tolower((int) as_con_type[k][m]) ;

	  if ((strcmp(as_con_type[k], "target") == 0) ||
	      (strcmp(as_con_type[k], "solvent_tar") == 0) ||
	      (strcmp(as_con_type[k], "stabilize_tar") == 0)) {
	       con_type[k] = TARGET ;
	       Ntargets++ ;
               if (strcmp(ta_filename[k], "none") == 0) {
	           sprintf(message, "Missing ta_filename for target # %d!", 
		   k+1) ;
	           EdenWarning(message) ;
	           status = -1 ; 
               }
               if (strcmp(as_con_type[k], "solvent_tar") == 0)
		    strcpy(target_type[k], "solvent ") ;
               else if (strcmp(as_con_type[k], "stabilize_tar") == 0)
		    strcpy(target_type[k], "stabilizing ") ;
               else
		    strcpy(target_type[k], "") ;
          }
	  else if (strcmp(as_con_type[k], "singlet") == 0) {
	       con_type[k] = SINGLET ;
               if (strcmp(ta_filename[k], "none") == 0) {
	           sprintf(message, 
		   "Missing ta_filename for singlet cost function!") ;
	           EdenWarning(message) ;
	           status = -1 ; 
               }
          }
	  else if (strcmp(as_con_type[k], "triplet") == 0) {
	       con_type[k] = TRIPLET ;
               if (strcmp(ta_filename[k], "none") == 0) {
	           sprintf(message, 
		   "Missing ta_filename for triplet cost function!") ;
	           EdenWarning(message) ;
	           status = -1 ; 
               }
          }
	  else if (strcmp(as_con_type[k], "phase_ext") == 0) {
	       con_type[k] = PHASE_EXT ;
               if (strcmp(ta_filename[k], "none") == 0) {
	           sprintf(message, 
		   "Missing filename for phase extension!") ; 
	           EdenWarning(message) ;
	           status = -1 ; 
               }
               if (strcmp(wt_filename[k], "none") == 0) {
	           EdenWarning(  
		   "Missing weight filename for phase extension!") ; 
	           status = -1 ; 
               }
	       if (phase_ext_res == 0) {
		   EdenWarning("Missing phase extension resolution!") ;
	           status = -1 ; 
               }
	       else if (phase_ext_res < input_res) {
		   sprintf(message,
	           "Phase_ext_res (%gA) should be no higher than input_res (%gA)!",
	           phase_ext_res, input_res) ;
	           EdenWarning(message) ;
	           status = -1 ; 
               }
               if (strncmp(mode_info, "corr", 4) != 0) {
	           EdenWarning("Use correction mode for phase extension!") ;
	           status = -1 ; 
               }
          }
	  else if (strcmp(as_con_type[k], "cs") == 0) {
	       con_type[k] = CS ;
	       ;	/* no special checks */
          }
	  else if (strcmp(as_con_type[k], "ncs") == 0) {
	       EdenWarning("NCS has been withdrawn from Eden") ;
	       status = -1 ; 
          }
	  else if (strcmp(as_con_type[k], "sayre") == 0) {
	       con_type[k] = SAYRE ;
	       if (cost_addend[k] == 0) {
	           EdenWarning( 
	  "Missing optimum value for Sayre's equation cost function!") ;
	           status = -1 ; 
               }
          }
	  else  {
	       EdenWarning( "Illegal constraint type") ;	
	       status = -1 ; 
          }
     }
     if (status != 0)
          EdenError("1 or more illegal constraint entries.") ;
}

				/*********************************************
                   		Read Mir input from [name].inp            
				*********************************************/
void readInput_Mir()    
{
     int	j, k, m ;
     int	r, mm ;
     char	gotcha ;

     /* Defaults */

     Nder = 0 ;

     for (m = 0; m < MAXDER; m++) {
          strcpy(fo_der_fn[m], "none") ;
          strcpy(fc_heavy_fn[m], "none") ;
	  fscale_der[m] = 1. ;
	  res_der[m] = input_res ;
     }
     relwt_der[0] = relwt_native ;
     res_der[0] = res_native ;

     for (m = 1; m < MAXDER; m++) 	
	  relwt_der[m] = 1. ;

     for (k = 0; k < Nlines; k++) {
          strcpy(nextline, allinp + k*MAXSTRING) ;
          sscanf(nextline, "%s", id) ;

          if (strcmp(id, "NDER") == 0) {
               sscanf(nextline, "%*s %d", &Nder) ;
	       goodline[k] = TRUE ;
	  }
          else if ((strncmp(id, "RELWT_DER", 9) == 0) && isdigit(id[9]))  {
	       j = atoi(&id[9]) ;
               sscanf(nextline, "%*s %g", &t1) ;
	       relwt_der[j] = t1 ;
	       goodline[k] = TRUE ;
	  }
          else if ((strncmp(id, "RES_DER", 7) == 0) && isdigit(id[7]))  {
	       j = atoi(&id[7]) ;
               sscanf(nextline, "%*s %g", &t1) ;
	       res_der[j] = t1 ;
	       goodline[k] = TRUE ;
	  }
          else if ((strncmp(id, "FSCALE_DER", 10) == 0) && isdigit(id[10]))  {
	       j = atoi(&id[10]) ;
               sscanf(nextline, "%*s %g", &t1) ;
	       fscale_der[j] = t1 ;
	       goodline[k] = TRUE ;
	  }
          else if ((strncmp(id, "FO_DER_FN", 9) == 0) && isdigit(id[9]))  {
	       j = atoi(&id[9]) - 1 ;
               sscanf(nextline, "%*s %s", fo_der_fn[j]) ;
	       goodline[k] = TRUE ;
	  }
          else if ((strncmp(id, "FC_HEAVY_FN", 11) == 0) && isdigit(id[11]))  {
	       j = atoi(&id[11]) - 1 ;
               sscanf(nextline, "%*s %s", fc_heavy_fn[j]) ;
	       goodline[k] = TRUE ;
	  }
     }

     /* Check for fatally flawed input */

     if (Nder > MAXDER) 
	 EdenError("Too many derivative files.") ;

     if (Nder > 0) {

          for (k = 0; k < Nder; k++) {
              if (strcmp(fc_heavy_fn[k], "none") == 0) {
                   sprintf(message, 
			"Missing heavy atom fc file name for k = %d!", k+1) ;
                   EdenError(message) ;
              }
              if (strcmp(fo_der_fn[k], "none") == 0) {
                   sprintf(message, 
		        "Missing derivative fo file name for k = %d!", k+1) ;
                   EdenError(message) ;
              }
          }
     }

     NM = Nder + 1 ;

     /************************************************************
     Do some tedious checks to determine the resolutions for which
     the exp_fac arrays must be set up.  By definition, the native
     uses the 1st resolution set.  Nres counts the number of 
     resolution sets (hopefully, not > 2!)
     ************************************************************/

     useset[0] = 0 ;	
     resrat[0] = res_der[0] / input_res ;

     for (r = 0, m = 1; m < NM; m++) {

	  useset[m] = 0 ;
	  gotcha = FALSE ;

	  for (mm = 0; mm < m; mm++) {
	       if (res_der[m] == res_der[mm]) {  	/* it's old  */
		    gotcha = TRUE ;
		    useset[m] = useset[mm] ;
		    break ;
               }
          }
	  if (!gotcha) {				/* it's new */
	       r++ ;
	       useset[m] = r ;
               resrat[r] = res_der[m] / input_res ;
          }
     }
     Nres = r+1 ;

     if (Nres > MAXRES)
	  EdenError("Too many intrinsic resolutions!") ;
}


				/*********************************************
     				Read special input for anomalous (MAD) data 
				*********************************************/
void readAnomInput()    
{
     int	i, j, k, m ;
     char	bad = FALSE ;

     /* Defaults */


     for (m = 1; m < NM; m++) {
	  fprime[m] = 1000. ;
	  f2prime[m] = 1000. ;
          Zheavy[m] = 0 ;
	  anom_der_flag[m] = FALSE ;
     }

     for (k = 0; k < Nlines; k++) {
          strcpy(nextline, allinp + k*MAXSTRING) ;
          sscanf(nextline, "%s", id) ;

          if (strncmp(id, "Z", 1) == 0) {
	       j = atoi(&id[1]) ;
               sscanf(nextline, "%*s %g", &t1) ;
	       Zheavy[j] = t1 ;
	       goodline[k] = TRUE ;
	  }
          else if ((strncmp(id, "FP_FPP", 6) == 0) && isdigit(id[6]))  {
	       j = atoi(&id[6]) ;
               sscanf(nextline, "%*s %g %g", &t1, &t2) ;
	       fprime[j] = t1 ;
	       f2prime[j] = t2 ;
	       goodline[k] = TRUE ;
	  }
          else if ((strncmp(id, "ANOM_DER", 8) == 0) && isdigit(id[8]))  {
	       j = atoi(&id[8]) ;

               strcpy(message, "") ;

               sscanf(nextline, "%*s %s", message) ;

               for (i=0; i<(int)strlen(message); i++)
                    message[i] = toupper((int) message[i]) ;

	       goodline[k] = TRUE ;
               if ((strcmp(message, "ON") == 0) ||
                   (strcmp(message, "TRUE") == 0))
	            anom_der_flag[j] = TRUE ;
               else if ((strcmp(message, "OFF") == 0) ||
                        (strcmp(message, "FALSE") == 0))
	            anom_der_flag[j] = FALSE ;
               else
	            goodline[k] = FALSE ;

          }
     }

     /* Check for fatally flawed input */

     for (m = 1; m < NM; m++) {

          if (anom_der_flag[m]) {
               if (fprime[m] == 1000.) {
	            sprintf(message, "Missing f prime for derivative %d!", m) ;
	            printTwice(message) ;
                    bad = TRUE ;
               }
               if (f2prime[m] == 1000.) {
	            sprintf(message, "Missing f doubleprime for derivative %d!", m) ;
	            printTwice(message) ;
                    bad = TRUE ;
               }
               if (Zheavy[m] == 0) {
                    sprintf(message, "Missing Z (at. no.) for derivative %d!", m) ;
	            printTwice(message) ;
                    bad = TRUE ;
               }
          }
     }
     if (bad) 
          EdenError("Exiting...") ;
}


#define	DFDX_SOLVE	3.e-2 
				/*********************************************
                       		Set defaults for global variables  
				*********************************************/
void setGlobalDefaults()    
{
     dfdxpc = DFDX_SOLVE ;		

     relwt_native = 1. ;
     res_native = UNSET ;
     R_stop = 0. ;
     discrp_frac = 1. ;
     fscale = UNSET ;
     Nconstraints = 0 ;
     Ntargets = 0 ;
     t_type = ' ' ;
     t_frac = UNSET ;
     t_matrix = NULL ;  
     hrCutoff = 10. ;

     strcpy(fo_filename, "none") ;
     strcpy(md_filename, "none") ;
     strcpy(title, "") ;
}

				/*********************************************
          			Echo information in .inp file to log.     	
				*********************************************/
void echo()      
{
     void	nice_print_to_log(char *, char *) ;
     int	k, m ;

     if ((int)strlen(title) > 0) 
          fprintf(fp_log, "\n%s\n\n", title) ;
     
     /************************************************************
     Information regarding input files.
     ************************************************************/

     fprintf(fp_log, "\n") ;
     nice_print_to_log("Observed structure factors will be read from ",
                            fo_filename) ;
     if (useSig) 
         fprintf(fp_log, "Sigmas will be used for weighting\n") ;
     else
         fprintf(fp_log, "Sigmas will not be used for weighting\n") ;
     fprintf(fp_log, "Data scaling factor is %g\n", fscale) ;

     nice_print_to_log("Structure factors will be calculated from ", md_filename) ;

     fprintf(fp_log, 
		"Stop getsol if df/dx is reduced to %g of its initial value\n",
		dfdxpc) ;

     nice_print_to_log("Starting physical space model will be read from ", 
                        md_filename) ;
     
     if (detwin) {
         fprintf(fp_log, "\n") ;
         if (t_type == 'A')
             fprintf(fp_log, "Amplitude detwinning is enabled \n") ;
         else if (t_type == 'I')
             fprintf(fp_log, "Intensity detwinning is enabled \n") ;

         fprintf(fp_log, "The twinning fraction is %g \n", t_frac) ;
         fprintf(fp_log, "The twinning matrix is \t %3d,%3d,%3d\n", 
                          *(t_matrix+0), *(t_matrix+1), *(t_matrix+2)) ;
         fprintf(fp_log, "                       \t %3d,%3d,%3d\n",
                          *(t_matrix+3), *(t_matrix+4), *(t_matrix+5)) ;
         fprintf(fp_log, "                       \t %3d,%3d,%3d\n",
                          *(t_matrix+6), *(t_matrix+7), *(t_matrix+8)) ;
     }

     /************************************************************
     Constraints, relative weights except for constraints already
     reported in the mode-of-operation banner. 
     ************************************************************/

     fprintf(fp_log, "\n") ;
     fprintf(fp_log, "Relative weight for Nhkl space is %g \n", relwt_native) ;
     if (Nder > 0)
	  fprintf(fp_log, "intrinsic resolution of the native is %g \n",
		  res_native) ;

     fprintf(fp_log, "\n") ;
     if (Ntargets > 0)
	  fprintf(fp_log, "Information on Np space constraints follows:\n\n") ; 

     for (k = 0; k < Nconstraints; k++) {
	  switch (con_type[k]) {
	  case TARGET:
              nice_print_to_log("Target Np values will be read from ", 
                                 ta_filename[k]) ;

	      nice_print_to_log("Target weights will be read from ",
                                 wt_filename[k]) ;
              fprintf(fp_log, 
	         "Target relative weight is %g \n", relwt_con[k]) ;
              break ;
     
	  default :
               break ;

          }
     }

     /************************************************************
     Mir information.  
     ************************************************************/

     if (Nder > 0)

	for (m = 0; m < Nder; m++) {
          fprintf(fp_log, "\n") ;
	  fprintf(fp_log, "For derivative # %d:\n", m+1) ;

          nice_print_to_log("\tfile of observed derivative structure factors is ",
		            fo_der_fn[m]) ;
          
          nice_print_to_log("\tfile of calculated heavy atom structure factors is ",
		            fc_heavy_fn[m]) ;
          
	  if (anom_der_flag[m+1]) {
	       fprintf(fp_log, "\tDerivative %d is anomalous\n", m+1) ;
	       fprintf(fp_log, "\tf' and f'' are %g and %g\n", 
		  fprime[m+1], f2prime[m+1]) ;
	       fprintf(fp_log, "\tZ for the heavy atom is %g\n", Zheavy[m+1]) ;
          }

          fprintf(fp_log, "\tdata scaling factor is %g \n", fscale_der[m+1]) ;
          fprintf(fp_log, "\trelative weight in the cost function is %g \n", 
		  relwt_der[m+1]) ;
	  fprintf(fp_log, "\tintrinsic resolution of the data file is %g \n",
		  res_der[m+1]) ;
     }

     fprintf(fp_log, "\n") ;
     fflush(fp_log) ;
}

void	echo_ignored()
{
     FILE	*fp ;
     int	k = 0 ;
     char 	first = TRUE ;

     /* Was anything ignored? */

     fp = fopen(input_filename, "r") ;

     while (fgets(nextline, MAXSTRING, fp) != NULL)  {
	  if (net_strlen(nextline) > 1) {
	       if (!goodline[k]) {
                    if (first) {
                         printTwice(
                    "The following input lines were ignored:\n") ;
                         first = FALSE ;
                    }
	            sprintf(message, "%s%s", ">>  ", nextline) ;
	            printTwice(message) ;
               }
               k++ ;
          }
     }
     printTwice("") ;
     fclose(fp) ;
}

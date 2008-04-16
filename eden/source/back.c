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

                                BACK

  Title:        Back transformation of calculated diffraction patterns.    
  Author:       Hanna Szoke
  Date:         11/30/92
  Function:     This program estimates electron densities from a set of
		calculated diffraction patterns with phases. It obtains
		a "solution map" - i.e., the amplitudes of a set of
		Gaussian densities of given width, eta*(grid_spacing)^2,
		centered on a simple grid or on a body-centered
		grid of given grid spacing. 

		It is assumed that the input fc file has been correctly
		apodized prior to running Back.

		The main purpose of this calculation is to provide Solve 
		with a map of initial values for the conjugate gradient 
		optimization process which is everywhere non-negative
		and consistent with the reciprocal-space representation.

                Output includes a log of the process, together with a
		binary file of electrons/voxel.
		The structure factors that are compatible with these 
		electron/voxel files are not written out as of 6/7/01.

*******************************************************************************/

#include "eden.h"
#include "eden_Np.h"
#include "cost.h"
#include "ccg.h"		/* needed for dfdxpc */

COMPLEX	*Fc0 ;		/* Calculated structure factor for known atoms */
static	real	*Fhold ;	/* ... to hold |Fc0| for diagnostics 	*/

static	char	b_title[MAXSTRING] ;	/* whatever the user wants ... */

static	void	back_setup() ;
static	void	makeFArray() ;
static	void	b_allocNhklArrays() ;
static	void	report_bmemory() ;
static	void	b_readInput() ;
static	void	b_echo() ;


void	back_main(int argc, char *argv[])
{
     void	minimal_b_funct(real *) ;	
     void	convolve(real *, real *, COMPLEX *, COMPLEX *) ;
     void	do_inplace_cs_symmetrization(real *, real *) ;
     real	get1Rfac(char *, real *, COMPLEX *, int) ;	/* " */
     void	shellHeader_R() ;				/* " */
     void      	b_HoloMethod() ;		/* in ccgutil.c */
     void	freeMinMaxType(int) ;		/* in ccgutil.c	*/
     void	analyzEpixFile(real *) ;	/* in qshare.c */
     void	writeEpixFile(char *, real *) ;	/* in qshare.c */
     void	Cmaskit() ;			/* in hklutil.c */

     real	sum = 0 ;
     real	residue = 0 ;			
     int        n ;
     char	problemName[MAXSTRING] = "" ;
     char	solution_filename[MAXSTRING] ;

  /******************************************************
  Check for right ballpark; 2/25/02: sf_filename may be 
  entered EITHER on execute line OR from input file!
  ******************************************************/

     strcpy(sf_filename, "none") ;

     if (argc > optind+3)  	/* too many arguments */
          ballpark(caller) ;

     if (argc > optind+2)
	  sprintf(sf_filename, argv[optind+2]) ;

  /*********************************************
  Pull off prefix defining whereabouts of input 
  file and remove .inp extension.
  Say hello, prepare cost file.                             
  *********************************************/

     strip_path(input_filename, pwd_name) ;
     strncpy(problemName, pwd_name, strlen(pwd_name) - 4) ;

     hello(caller) ;

     if (verbose)
          costfile_option = TRUE ;
     else
          costfile_option = FALSE ;

     if (costfile_option) {
	  sprintf(cost_filename, "%s_back.cost", problemName) ;

          if ((fp_cost = fopen(cost_filename,"w")) == NULL) 
          {
               sprintf(message, "Cannot open %s", cost_filename);
               EdenError(message) ;
          }
     }

  /*********************************************
  Fetch input conditions, do all initializations.
  Report initial standard deviation 
  *********************************************/

     back_setup() ;

     minimal_b_funct(&residue) ;

     sprintf(message, 
	 "\nInitial standard deviation = %g", sqrt(residue / (Nhkl-1))) ;

     printTwice(message) ;

  /*********************************************
  Do holographic reconstruction.                           
  *********************************************/

     sprintf(message, "\nApplying phased holographic reconstruction \n") ;
     printTwice(message) ;

     b_HoloMethod();

  /*********************************************
  Symmetrize resulting Np arrays.                         
  Report final standard deviation.
  *********************************************/

     if (Ncs > 1)  
          do_inplace_cs_symmetrization(nump, nump) ; 

     minimal_b_funct(&residue) ;

  /*********************************************
  Summarize results, write electron/pixel arrays.
  *********************************************/

     for (n = 0; n < Nptotal; n++) 
          sum += *(nump + n) ;
          
     sprintf(message, "Sum of recovered electrons is %g\n", sum) ;
     printTwice(message) ;

  /*********************************************
  Diagnostics.  Note freeMinMaxType() to 
  provide space for writeEpixFile
  *********************************************/

     freeMinMaxType(Nptotal) ;
     sprintf(solution_filename, "%s_back", problemName) ;
     analyzEpixFile(nump) ;
     writeEpixFile(solution_filename, nump) ; 

     makeFArray() ;
     convolve(nump, expfac, expfac_bcc, Fc0);
     Cmaskit(Fc0, maskfc, Nhkl);
     shellHeader_R() ;
     get1Rfac(" Model:  ", Fhold, Fc0, Nhkl) ;

     return ;
}

void	back_setup()
{
     void	setup_csind() ;			/* in symutil.c */

     void	hklTabHead() ;			/* in hkutil.c */
     int	readfcalc(char *, COMPLEX *, char *) ;	/* in hkutil.c */
     void	setDsqLimit(char *, int) ;		/* in hkutil.c */

     void       setupWeights() ;		/* in hklinit.c */
     real     *setupexpfac(real) ;		/* in hklinit.c */
     COMPLEX    *setupexpfac_bcc(real) ; 	/* in hklinit.c */
     void	echo_ignored() ;		/* in init.c	*/

     void	initMinMaxType(int) ;		/* in ccgutil.c */
     void	alloc_ccg_mem(int) ;		/* in ccg.c */

     int        n ;
     int	Neq ;			/* # equations = # non-zero Fc0 values*/

  /*********************************************
  Read input parameter file.                               
  *********************************************/

     b_readInput() ;
     b_echo() ;
     echo_ignored() ;
     report_bmemory() ;	

  /*********************************************
  Set up arrays in (h,k,l) reciprocal lattice 
  triplets.  Prepare for FFTs.                                    
  *********************************************/

     printTwice("Setting up initial arrays ...") ;

     b_allocNhklArrays() ;

  /*********************************************
  Read a set of diffraction patterns, Fc0.
  It is automatically expanded to P1.
  Prepare exp factors for FFTs using 
  delfac = PISQ*eta*dr*dr.
  *********************************************/

     setupWeights() ; 
     Neq = readfcalc(sf_filename, Fc0, maskfc) ;

     expfac = setupexpfac(delfac) ;
     if (grid_type == BODY_CENTERED)
          expfac_bcc = setupexpfac_bcc(delfac) ;

     setDsqLimit(maskfc, Nhkl) ;
     hklTabHead() ;

  /*********************************************
  Identify scope of problem.                               
  *********************************************/

     sprintf(message,
           "\n\t\tThis problem has %d equations, %d unknowns\n", Neq, Nptotal) ;
     printTwice(message) ;

  /*********************************************
  Prepare arrays for the electron densities and
  for the conjugate gradient solver.  Set up 
  symmetry arrays.                                  
  *********************************************/

     nump = (real *) e_malloc(Nptotal*sizeof(real), "back_setup") ;
          
     for (n = 0; n < Nptotal; n++)
	  *(nump+n) = 0 ;

     initMinMaxType(Nptotal) ;
     alloc_ccg_mem(Nptotal) ;
     if (Ncs > 1)  
          setup_csind() ;

    
     return ;
}

void	b_allocNhklArrays()
{
     void	initfft() ; 
     void	initfft_aux() ; 

  /*******************************************************
  Set up arrays in (h,k,l) reciprocal lattice triplets.
  Prepare for FFTs.                                       
  *******************************************************/

     strcpy(message, "b_allocNhklArrays") ;

     maskfc =   (char *) e_malloc(Nhkl*sizeof(char),    message) ;
     Fc0   = (COMPLEX *) e_malloc(Nhkl*sizeof(COMPLEX), message) ;
     Oar   = (COMPLEX *) e_malloc(Nhkl*sizeof(COMPLEX), message) ;
     
     initfft() ;
     initfft_aux() ;
}

void	makeFArray() 
{
     real	Cabs();
     int        n ;

  /****************************************
   Make Fhold = |Fc0| for diagnostics      
  ****************************************/

     Fhold =  (real *) e_malloc(Nhkl*sizeof(real), "makeFArray") ;

     for (n = 0; n < Nhkl; n++) 
          *(Fhold+n) = Cabs(*(Fc0+n)) ; 
}


void	report_bmemory()
{

     real	mem_Npspace, mem_recspace, mem_tot ;
     int	Nptot_gridPt ;
     int	Nex ;

	/*******************************************************
	In Np space, there are 3 Eden real arrays (nump) 
	minp and maxp) plus 1 Eden int array (typep) and 1 char
	array (unique_pmask), 4 CCG real arrays (xn, gn, go 
	and d), each with Nptotal elements.  There is also a 
	complex FFT array (temp) in Np, NOT Nptotal!
	*******************************************************/
     
     Nptot_gridPt = 7 * sizeof(real) + 1 * sizeof(int) + 1 * sizeof(char) ; 

     mem_Npspace = Nptot_gridPt * Nptotal + sizeof(COMPLEX) * Np ;

	/*******************************************************
	In reciprocal space, there are 3 or 4 Eden complex 
	arrays in Nhkl (Fc0, Oar, fftarray, expfac_bcc).  
	There are 2 real arrays in Nhkl (hkl_weight and expfac). 
	There is 1 byte array (maskfc). 
	We ignore arrays that are 2-dimensional.
	*******************************************************/

     Nex = (grid_type == SIMPLE) ? 0 : 1 ;

     mem_recspace = (3 + Nex) * Nhkl * sizeof(COMPLEX) +
		    2 * Nhkl * sizeof(real) + 
		    1 * Nhkl * sizeof(char) ;

     mem_Npspace *= 1.e-6 ;
     mem_recspace *= 1.e-6 ;
     mem_tot = mem_Npspace + mem_recspace ;

     printTwice("Approximate (underestimated) memory requirements in Mbytes:") ;
     sprintf(message, 
     "Physical space: %4.3g, Reciprocal space: %4.3g, Total: %4.3g.\n",
     mem_Npspace, mem_recspace, mem_tot) ;
     printTwice(message) ;
}


#define	DFDX_BACK	1.e-3 

void	b_readInput()    
{
     int	k ;

     /************************************************************
     Read basics; set defaults for global variables.
     ************************************************************/

     readBasicInput() ;

     dfdxpc = DFDX_BACK ;		/* changed most recently 9/28/98 */
     discrp_frac = 1. ;

     strcpy(b_title, "") ;

     /************************************************************
     Now read input that is specific to Back
     *************************************************************/

     for (k = 0; k < Nlines; k++) {
          strcpy(nextline, allinp + k*MAXSTRING) ;
          sscanf(nextline, "%s", id) ;

          if (strcmp(id, "TITLE") == 0) {
               strcpy(b_title, nextline+6) ;
	       goodline[k] = TRUE ;
          }
          else if ((strcmp(id, "FC_FILENAME") == 0) ||
              (strcmp(id, "FC_FN") == 0)) {

               if (strcmp(sf_filename, "none") == 0) {
                    sscanf(nextline, "%*s %s", sf_filename) ;
	            goodline[k] = TRUE ;
	       }
          }
          else if (strcmp(id, "DFDX_CRIT") == 0) {
               sscanf(nextline, "%*s %g", &t1) ;
               dfdxpc = t1 ;
	       goodline[k] = TRUE ;
          }
     }

     /************************************************************
     Check for fatally flawed input.
     ************************************************************/

     if (strcmp(sf_filename, "none") == 0)
          EdenError("Missing fc file name!") ;

}

				/*********************************************
          			Echo information in .inp file to log.     	
				*********************************************/
void b_echo()      
{
     void nice_print_to_log(char *, char *) ;

     if ((int)strlen(b_title) > 0) 
          fprintf(fp_log, "\n%s\n\n", b_title) ;
     
     /************************************************************
     Information regarding input files.
     ************************************************************/

     fprintf(fp_log, "\n") ;

     nice_print_to_log(
          "Calculated structure factors will be read from ", sf_filename) ;

     fprintf(fp_log, 
		"Stop getsol if df/dx is reduced to %g of its initial value\n",
		dfdxpc) ;

     fprintf(fp_log, "\n") ;
     fflush(fp_log) ;
}

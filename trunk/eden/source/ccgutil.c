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

                                CCGUTIL.C

  Title:        Interface with conjugate gradient solver package.
  Author:       Hanna Szoke
  Date:         Aug. 27, 1997 (separated from qshare.c)
  Function:     Provide interface between Solve or Back and the conjugate
		gradient solver package, ccg.c.

void funct()		- ccg	
int  doStdDev()		- solve
int  HoloMethod() 	- solve
void b_HoloMethod() 	- back
int  DGsolve() 		- (local)
void reportTotFunctCalls()	- solve
void initMinMaxType() 	- back, qshare
void reviseMinMaxp() 	- qshare

******************************************************************************/

#include "util.h"		
#include "dims.h"		
#include "ccg.h"
#include "eden_Np.h"		/* needed for nump */
#include "mir.h"		/* needed for NM in reporting standard dev. */

#define	COUNTER	3		/* # iterations over which we track std. dev. */

	/* declarations of local variables */

static	real  *minp ;                 /* min. electrons/voxel */
static	real  *maxp ;                 /* max. electrons/voxel */
static	int	*typep ;		/* flag for type of solution point */
static	real	std_dev[COUNTER] ;
static	int	totFunctCalls = 0 ;	/* total number of calls to "funct" */
static	int	Nclsrch = 0 ;		/* total number of calls to "clsrch" */

static	int   	DGsolve(real *) ;

	/*************************************************************
	nc ;		 number of variables (Nptotal)
	csol ;		 current solution array (nump or totp)
	ressq ;		 (returned) residue  		
	grad ;		 (returned) gradient array
	iqflag ;	 (returned) flag re success	
	*************************************************************/
void funct(int nc, real *csol, real *ressq, real *grad, int *iqflag)
{
     void	e_funct(real *, real *, real *, int *) ;
     void	b_funct(real *, real *, real *, int *) ;

     if (strcmp(caller, "solve") == 0)
	 e_funct(csol, ressq, grad, iqflag) ;  
     else
         b_funct(csol, ressq, grad, iqflag) ;

}


#define DEL_STD_DEV	0.03	/* min. fract. change in std. dev. over COUNTER 
				iterations before code quits */

				/*********************************************
				Report standard deviation  
				*********************************************/
int doStdDev(char *descrip) 
{
     void	e_funct(real *, real *, real *, int *) ;
     static 	int	count  = -1;
     real	*junk ;
     real	residue ;
     int	iqflag ;	/* (unused) */
     int	j ;

       /* First time in - initialize deviations array  */

     if (count == -1) 
	for (j = 0; j < COUNTER; j++)
          std_dev[j] = 1.e10 ;


     count = (count+1) % COUNTER ;

	/* call funct(), using junk as a spare array for unneeded gradients */

     junk = (real *) e_malloc(Npextended*sizeof(real), "DoStdDev - temp.") ;

     e_funct(nump, &residue, junk, &iqflag) ;

     e_free(Npextended*sizeof(real), junk) ;

     std_dev[count] = sqrt(residue / (Nhkl*NM)) ;

     sprintf(message, "%s = %g", descrip, std_dev[count]) ;
     printTwice(message) ;

       /************************************************************
       The following criterion is often invoked to end an eden run. 
       If after 3 calls (initial, before symmetrization, after 
       symmetrization), the standard deviation does not decrease by 
       more than DEL_STD_DEV (3%), quit.     
       *************************************************************/

     if ((count == 2) && 
	 ((std_dev[0] - std_dev[count]) < DEL_STD_DEV * std_dev[count])){
	printTwice("\nStopping - standard deviation is not decreasing ...\n") ;
	return(STOP) ;
     }

     return(GO_ON) ;
}

		/************************************************************
       		Solve holographic reconstruction equation in Eden
		************************************************************/
int HoloMethod()            
{
     int	DGval ;
     int        n ;
     real     sum = 0 ;		/* sum of solutions, this iteration */
     real     residue ;
     real	hrT ;
     real	*array ;

       /************************************************************
       DGsolve is where it's all done. 
       ************************************************************/

     DGval = DGsolve(&residue) ;

     for (n = 0; n < Nptotal; n++) 
          sum += *(nump + n) ;

     if (HighRes) {
          for (n = 0, hrT = 0, array = nump+Nptotal; n <Nhr; n++, array++)
              hrT  += *array ;

          sprintf(message,
	    "Sum of recovered gridded electrons in this iteration is %g", sum) ;
          printTwice(message) ;
          sprintf(message, "Sum of high-resolution electrons is %g", hrT) ;
          printTwice(message) ;
          sprintf(message, "Total is %g", sum + hrT) ;
          printTwice(message) ;
     }
     else {
          sprintf(message,
	    "Sum of recovered electrons in this iteration is %g", sum) ;
          printTwice(message) ;
     }

     if (DGval == STOP)
	return(STOP);
     else
	return(GO_ON) ;
}
		/*********************************************
       		Solve holographic reconstruction equation in Back
		**********************************************/

void b_HoloMethod()            
{

     int	DGval ;
     real     residue ;

     DGval = DGsolve(&residue) ;

     sprintf(message, 
	    "Standard deviation of solution = %g", sqrt(residue / (Nhkl-1))) ;
     printTwice(message) ;

     return ;
}

		/************************************************************
                Solve Mat*nump = Hdiff by the conjugate gradient method, 
		using FFTs in place of matrix multiplications.      
		residue is the value of 'funct' at minimum.
		************************************************************/

int DGsolve(real *residue)  
{
     void getsol(int, real *, real *, real *, int *, int, int *, int *, 
		 int *, real *, real, real *, real, int *);

     int	iter = MAXIT ;	/* max # allowed iterations	*/
     int	itn ;		/* actual # iterations performed	*/
     int	ifn ;		/* # times 'funct' was invoked	*/
     int	istop = 0 ;	/* stopping condition		*/
     real	fn ;		/* value of function at minimum */
     real	fmin = 0 ;	/* (input) min 'funct' value	*/
     real	tol = TOL ;	/* (input) relative tolerance	*/
     static	real	df0 = 0;	
				/* Initial df (used & returned by getsol)*/
     int	status ;	/* return status */


     getsol(Npextended, nump, minp, maxp, typep, 
	    iter, &itn, &ifn, &istop, &fn, fmin, &df0, tol, &Nclsrch) ; 

     totFunctCalls += ifn ;		
     *residue = fn ;

     switch (istop)
     {
     case NORMAL_STOP: 
          sprintf(message, 
		 "\ngetsol worked, %d funct calls", ifn) ;
          status = GO_ON ;
	  break ;
     case RESSQ_LT_DISCRP:
          sprintf(message, 
	    "\nStopping - discrepancy principle satisfied; %d funct calls", 
	    ifn) ;
          status = STOP ;
	  break ;
     case FN_LT_FMIN: 
          sprintf(message, 
	    "\nCost function is negative, %d funct calls.", ifn-1) ;
          status = STOP ;
	  break ;
     case MAX_ITER_LIMIT:
          sprintf(message, 
	 "\nToo many iterations in Getsol. The solver is STUCK - please re-examine your input!") ;
          status = STOP ;
	  break ;
     case MAX_FN_CALLS_EXTRAP:
     case MAX_FN_CALLS_INTERP:
          sprintf(message, 
	 "\nToo many iterations in Clsrch. The solver is STUCK - please re-examine your input!") ;
          status = STOP ;
	  break ;
     case DX_TOO_SMALL:
          sprintf(message, 
		 "\nDead in the water - making no progress ... Quitting.") ;
          status = STOP ;
	  break ;
     case INTERRUPT: 
          sprintf(message, 
		 "\nQuitting - interrupt received") ;
          status = STOP ;
	  break ;
     case GRADIENT: 
          sprintf(message, 
		 "\ndf/dx went down enough, %d funct calls", ifn) ;
          status = GO_ON ;
	  break ;
     default:
	  sprintf(message,
		 "\nOther reason for ending getsol; istop = %d", istop) ;
          status = GO_ON ;
	  break ;
     }
     printTwice(message) ;

     return (status) ;
}

void reportTotFunctCalls()
{
     sprintf(message, 
	  "Total number of solver search directions was %d", Nclsrch) ;
     printTwice(message) ;
     sprintf(message, 
	  "Total number of cost function calls was %d\n", totFunctCalls) ;
     printTwice(message) ;
     return ;
}

		/************************************************************
		Allocate and initialize the arrays specific to "getsol".
		>>> This is where the non-negativity constraint 
		    on electron density is established!
		Note: Back calls this with N=Nptotal, Solve with N=Npextended.
		************************************************************/
void	initMinMaxType(int N) 
{
     int	n ;

     strcpy(message, "ccgutil") ;

     minp = (real *) e_malloc(N*sizeof(real), message) ;
     maxp = (real *) e_malloc(N*sizeof(real), message) ;
     typep = (int *)   e_malloc(N*sizeof(int), message) ;

     for (n = 0; n < N; n++) {
	  *(minp+n) = min_voxel ;	/* defined in dims.h ... */
	  *(maxp+n) = max_voxel ;	/* ... but user-redefinable */
	  *(typep+n) = BOUNDED_NC ;	/* defined in ccg.h  */
	  }
}
		/************************************************************
		Revise minp and maxp to account for known information, 
		either in correction mode at start of solve process. 
		or between outer solver iterations, to account for 
		el/voxel values just found.
		Obviously, if the upper bound on the solver is very large,
		there is no particular point in pushing it down a little; 
		however, things can change...
		************************************************************/

void reviseMinMaxp(int N, real *array)	
{
     int        n ;

     for (n = 0; n < N; n++, array++) {
	*(minp +n) -= *array ;
	*(maxp +n) -= *array ;
     }

     return ;
}
		/************************************************************
		Free up the memory associated with xmin, xmax & typep arrays.  
		************************************************************/

void freeMinMaxType(N)	
{

     e_free(N*sizeof(real), minp) ;
     e_free(N*sizeof(real), maxp) ;
     e_free(N*sizeof(int), typep) ;

     return ;
}

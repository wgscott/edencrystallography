/************************************************************************

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
*************************************************************************/
/*                                                                      */
/*                                                                      */
/*     Program: getsol.c                                                */
/*     Date:    2/24/95                                                 */
/*     Author:  Erik M. Johansson, translation of Dennis Goodman's      */
/*              FORTRAN constrained least squares conjugate gradient    */
/*              codes                                                   */
/*                                                                      */
/*     (c) Copyright 1995 the Regents of the University                 */
/*         of California.  All rights reserved.                         */
/*                                                                      */
/*     This software is a result of work performed at Lawrence          */
/*     Livermore National Laboratory.  The United States Government     */
/*     retains certain rights therein.                                  */
/*                                                                      */
/************************************************************************/

#include <stdio.h>
#include <math.h>
#include "util.h"
#include "ccg.h"
#include <signal.h>
#include <setjmp.h>

static	real	*xn, *gn, *go, *d ;
  
static
void clsrch(int, real *, real *, real *, real *, real *, real *, 
	    int *, real *, real *, real *, real, real, real, int *,
	    int *, int *);
static	void update(int, real, int *, real *, real *, real *, real *,
	    real *, int *);
static	void update_d(int, real *, int *, real *, real *, int *);
static	real dot_product_ivec(int, real *, real *, int *);
static	void countem(int, int *);
static	void	handle_ctrl_c() ;

static	jmp_buf	sjbuf ;

/************************************************************************

 This is the general nonlinear version of getsol.

 Algorithm to minimize a function of multiple variables that are
 constrained to be nonnegative using a modified version of the
 conjugate gradient algorithm. The algorithm is described in
 "On applying the conjugate-gradient algorithm to image processing
 problems," D.M. Goodman, E.M. Johansson, and T.W. Lawrence. A
 book chapter in "Multivariate Analysis: Future Directions,"
 C.R. Rao, editor, 1993.

 The algorithm was implemented in a Fortran version of getsol and clsrch
 exactly as described in that article; since then, Goodman has made a 
 number of changes that are reflected in the current C version, the most
 important change being the following:

 In the original code, the new search direction was orthogonalized to xn - xo, 
 where xo is the original position and xn is the position returned by clsrch.
 In the new code, the new search direction is orthogonal to the last search
 direction that clsrch used before it got to its minimum.  The difference 
 between these two methods lies in the variables that were bent because they
 hit a boundary.  The new approach assures that the variables that were stuck
 on a boundary but became free in the new search direction always move in a
 descent direction.  But note that orthogonality to all the previous directions
 can never be assured.

 Dennis M. Goodman
 L-495 PO Box 808
 Lawrence Livermore National Lab.
 Livermore, CA 94550
 Tel. (510)-423-7893

************************************************************************/

void getsol(int nc, real *xo, real *xmin, real *xmax, int *ivec, 
	    int iter, int *itn, int *ifn, int *istop, real *fn, 
	    real fmin, real *df0, real tol, int *Nclsrch)

/************************************************************************
	 on entry: 

	nc	 the number of variables; unchanged by getsol		
	xo	 vector of variables.  Set to initial estimate.
  		   On return, it holds the location of the minimum  	
	xmin	 lower bounds on the variables			
	*xmax	 upper bounds on the variables			
	ivec	 array of integers providing information about the
		variables.  See ccg.h					
	iter	 maximum number of iterations allowed in getsol	
	fmin	 smallest possible value for function being minimized	
	df0	 initial value for df for this call to getsol.	
                   On return, same for next call to getsol		
	tol	 relative tolerance for function			

	 on return:

	itn	 the number of iterations in getsol			
	ifn	 the number of times funct was evaluated		
	istop	 the reason for stopping.  See ccg.h			
	fn	 the value of funct at the minimum			

	 on entry and return:

	Nclsrch 	cumulative # calls to clsrch() 		
************************************************************************/

{
  void	funct() ;

  int j, restrt, istuck, itr, tstdf, tstdx, tstgr, tstga, icatch, debug;
  int	iqflag ;
  static real dfdx0 = 0 ;
  real eps, sqteps, tol2, tol3, snc, df, dg, fo, fold, tola, amx,
    dx, xnm, gfre, glofre, ghifre, gogo, top, bottom, xdiff, gnj, gnj2,
    gdiff, rgnrm, beta, beta1, beta2, dfdx, dnm, gnm, ap, bitst, gjsqr;
  real *p_xmin, *p_xmax, *p_xo, *p_xn, *p_d, *p_go,
    *p_gn;
  int *p_ivec;

  /* initialization phase ***************************************************/

  /* create the temporary arrays */

  /* machine precision, eps, is defined as the smallest representable */
  /* number such that 1.0 + eps != 1.0. For 64 bit IEEE floating point, */
  /* eps = 2^(-52)  */
  eps = EPS ;

  /* square root of machine precision */
  sqteps = sqrt(eps);

  /* tolerances */
  tol2 = sqrt(tol);
  tol3 = pow(tol, 1./3.);

  /* square root of the number of variables */
  snc = sqrt((real) nc);
  
  /* controls output from code */
  debug = FALSE;

  /* check to be certain that the number of variables is > 0 */
  if (nc < 1)
      EdenError("Error in getsol, nc <= 0, stopping.");

  /* Force bounded variables in starting estimate to be within their bounds. */
  p_ivec = ivec;
  p_xmin = xmin;
  p_xmax = xmax;
  p_xo = xo;
  for (j = 0; j < nc; j++, p_ivec++, p_xo++, p_xmin++, p_xmax++)
    {
      switch (*p_ivec)
	{
	case BOUNDED_NC:
	case BOUNDED_LC_FREE:
	case BOUNDED_UC_FREE:
	case BOUNDED_LC_FIXED:
	case BOUNDED_UC_FIXED:
	    if (*p_xo <= *p_xmin)
	      {
		*p_xo = *p_xmin;
		*p_ivec = BOUNDED_LC_FREE;
	      }
	    else if (*p_xo >= *p_xmax)
	     {
	       *p_xo = *p_xmax;
	       *p_ivec = BOUNDED_UC_FREE;
	     }
	    else
	      *p_ivec = BOUNDED_NC;

	    break;
	case UNCONSTRAINED:
	case FIXED:
	    break;
	default:
	  sprintf(message, "Error: illegal value in ivec: %d", *p_ivec);
	  EdenError(message) ;
	}
    }

    /* set xn=xo for all cases ***/

    memcpy(xn, xo, nc * sizeof(real));

  /* Initialize restrt, ifn, res, fn, gn, df and go */
  restrt = FALSE;
  *ifn = 1;

  funct(nc, xo, fn, gn, &iqflag) ;

  if (*fn < fmin)
    {
      *istop = FN_LT_FMIN;
      return;
    }
  else if (iqflag == -1)    
    {
      *istop = RESSQ_LT_DISCRP;
      return;
    }

  if (*df0 == 0)
      df = fmin - *fn;
  else
      df = *df0;	

  /* Put gn into go, calculate initial d = -g+ and dg = |g+|^2, reset ivec.
  istuck counts the number of bounded variables at upper or lower 
  bounds that will not move in next linesearch iteration (i.e., 
  those variables with ivec[j] = BOUNDED_LC_FIXED or BOUNDED_UC_FIXED) */

  memcpy(go, gn, nc * sizeof(real));
  update_d(nc, &dg, &istuck, gn, d, ivec);

  /* check if initial reduced gradient is small enough to stop here */
  if (sqrt(fabs((double)dg)) <= snc * eps)
    {
      if (debug)
	{
	  fprintf(stdout, "GETSOL: Initial gradient small enough to stop\n");
	  fprintf(stdout, "Norm = %g\n", sqrt(fabs(dg)));
	}
      *istop = NORMAL_STOP;
      return;
    }



  if (debug)
    fprintf(stdout, "Begin optimization, fn = %g dg = %g istuck = %d\n",
	    *fn, dg, istuck);

  fold = *fn ;

  signal(SIGINT,handle_ctrl_c) ;
  setjmp(sjbuf) ;

  /* main iteration loop *************************************************/

  for (*itn = 1; *itn <= iter; (*itn)++)
    {
      fo = *fn;

      /* absolute tolerance for linesearch. */
      tola = tol * (1.0 + fabs(*fn));

      /* set maximum size of the steplength and estimated steplength     */

      amx = (fmin - *fn) / (CCG_RHO*dg) ;
      if ((-df) > 10*eps)
	   ap = 2 * df / dg ;
      else
	   ap = -2 * 10 * eps / dg ;

      if (debug)
	fprintf(stdout, "iteration = %d fn = %g ap = %g dg = %g is = %d\n",
		*itn, *fn, ap, dg, istuck);

      /* perform the line search */

      (*Nclsrch)++ ;
      if (costfile_option)
	   fprintf(fp_cost, "%d: New search direction.\n", *Nclsrch) ;

      clsrch(nc, xo, xn, gn, d, xmin, xmax, ivec, fn, &dg, &ap, 
	     amx, fmin, tola, &itr, istop, &istuck);
    
      *ifn += itr;
      *df0 = *fn - fold ;
      fold = *fn ;

      switch (*istop)
	{
	case RESSQ_LT_DISCRP:
	case FN_LT_FMIN:
	case INTERRUPT:
	  memcpy(xo, xn, nc * sizeof(real));
	  return;
	case INITIAL_AP_LT_0:
	case INITIAL_DD_GT_0:
	  return;
	default:
	  break;
	}
      
      /* Compute norms of solution (xnm) and change in solution (dx). */
      /* Compute the three parts of reduced gradient: gfre, glofre, */
      /* and ghifre. */
      
      dx = xnm = gfre = glofre = ghifre = 0.0;
      p_ivec = ivec;
      p_xo = xo;
      p_xn = xn;
      p_gn = gn;
      for (j = 0; j < nc; j++, p_ivec++, p_xo++, p_xn++, p_gn++)
	{
	  gnj = *p_gn;
	  gnj2 = gnj * gnj;
	  xdiff = *p_xn - *p_xo;
	  dx += xdiff * xdiff;
	  xnm += *p_xn * *p_xn;
	  switch (*p_ivec)
	    {
	    case BOUNDED_NC:
	    case UNCONSTRAINED:
	      gfre +=  gnj2;
	      break;
	    case BOUNDED_LC_FREE:
	    case BOUNDED_LC_FIXED:
	      if (gnj < 0.0)
		glofre += gnj2;
	      break;
	    case BOUNDED_UC_FREE:
	    case BOUNDED_UC_FIXED:
	      if (gnj > 0.0)
		ghifre +=  gnj2;
	      break;
	    case FIXED:
	      break;
	    }
	}

      dx = sqrt(dx);
      xnm = sqrt(xnm);

      /* Check stopping tests from pp 306-7 of Gill Murray and Wright. */
      rgnrm = sqrt(gfre + glofre + ghifre);
      df = *fn - fo;

      tstdf = - df < (tol *(1.0 + fabs(*fn)));
      tstdx = dx < (tol2 * (1.0 + xnm));
      tstgr = rgnrm <= (tol3 * snc * (1.0 + fabs(*fn)));
      tstga = rgnrm <= snc * eps;

      if ((tstdf && tstdx && tstgr) || tstga)
	{
	  if (debug)
	    fprintf(stdout, "Normal stop in getsol.\n");
	  *istop = NORMAL_STOP;
	  return;
	}

       /*----------------------------------------------------------
       EDEN CONDITION: if df/dx goes down enough, get out (normal 
       stop).  As of 11/18/98, we do a memcpy of xn to xo, in order
       to get the best possible result and one for which a direct
       call to e_funct() from doStdDev() will give the same result 
       as the most recent regular cost function calculation.  
       ----------------------------------------------------------*/

      dfdx = fabs(df / dx) ;
      if (*itn == 1) 
	   dfdx0 = dfdx ;

      else if (dfdx < dfdxpc * dfdx0)
      {
         if (debug)
           fprintf(stdout, "Gradient went down enough in getsol.\n");

         memcpy(xo, xn, nc * sizeof(real));

         *istop = GRADIENT;
         return ;
      }
      if (dx <= sqteps * snc * xnm)
	{
	  if (restrt)
	    {
	      if (debug)
		fprintf(stdout, "No progress in getsol: stop.\n");
	      *istop = DX_TOO_SMALL;
	      return;
	    }
	  else
	    restrt = TRUE;
	}
      else
	restrt = FALSE;

       /*----------------------------------------------------------
       Determine which variables that are on their boundaries will 
       be allowed to move away from them. This is done by selecting
       those whose negative gradients point away from the boundary.
       However, the absolute value of these gradients must be greater 
       than bitst. This may prevent "chattering." Also, compute the 
       numerator (top) and denominators (bottom and gogo) used to 
       calculate beta.  Because we used xn - xo rather than d in 
       calculating the denominator for the Hestenes-Stiefel version, 
       we have to scale with ap.
       In our applications, beta1 should be used. 
       ----------------------------------------------------------*/

      bitst = (GAMMA * sqrt(gfre)) / snc;
      top = bottom = gogo = 0.0;
      icatch = 0;
      p_ivec = ivec;
      p_xo = xo;
      p_xn = xn;
      p_gn = gn;
      p_go = go;
      for (j = 0; j < nc; j++, p_ivec++, p_xn++, p_xo++, p_gn++, p_go++)
	{
	  xdiff = *p_xn - *p_xo;
	  gdiff = *p_gn - *p_go;

	  switch (*p_ivec)
	    {
	    case BOUNDED_LC_FREE:
	    case BOUNDED_LC_FIXED:
	      if (*p_gn >= - bitst) 
		  *p_ivec = BOUNDED_LC_FIXED;
              else
		  {
		  top += *p_gn * gdiff;
		  icatch++;
		  *p_ivec = BOUNDED_LC_FREE;
		  }
	      break;
	    case BOUNDED_UC_FREE:
	    case BOUNDED_UC_FIXED:
	      if (*p_gn <= bitst) 
		  *p_ivec = BOUNDED_UC_FIXED;
              else
	          {
		  top += *p_gn * gdiff;
		  icatch++;
		  *p_ivec = BOUNDED_UC_FREE;
		  }
	      break; 	
	    case BOUNDED_NC:
	    case UNCONSTRAINED:
	      top += *p_gn * gdiff;
	      bottom += xdiff * gdiff; 
	      gogo += *p_go * *p_go;
	      break;
	    case FIXED:
	      break;
	    }

	}

      beta1 = ap * (top / bottom);  /* Hestenes-Stiefel version */
      beta2 = top / gogo;           /* Polak-Ribiere version */

      beta = (bottom != 0) ? beta1 : 0 ; 

      if (beta < 0) {
	   printTwice("New gradient search");
	   if (costfile_option)
		fprintf(fp_cost, 
		"New gradient search: resetting negative beta to 0\n");
	   beta = 0;
      }
				/* .. negative beta causes istop=-3 
				if we get here often, something more
				sophisticated is needed.  (5/23/95) */

      if (debug)
	fprintf(stdout, "H-S BETA = %g P-R BETA = %g\n", beta1, beta2);

      if (debug)
	{
	  dfdx = 1.0e20;
	  if( fabs(dx) >= 1.0e-20)
	    dfdx = df / dx;
	  fprintf(stdout, 
	  "it = %d fn = %g ap = %g dg = %g is = %d\n icatch = %d df = %g dx = %g dfdx = %g\n",
		  itr, *fn, ap, dg, istuck, icatch, df, dx, dfdx);
	}

      /* Calculate the new direction, calculate the new directional 
      derivative, and calculate the norms of the direction vector and 
      the gradient vector for the current set of unconstrained variables. */

      dg = dnm = gnm = 0.0;
      istuck = 0;
      p_ivec = ivec;
      p_d = d;
      p_gn = gn;
      for (j = 0; j < nc; j++, p_ivec++, p_d++, p_gn++)

	switch (*p_ivec)
	  {
	  case BOUNDED_NC:
	  case UNCONSTRAINED:
	    *p_d = - *p_gn + beta * *p_d;   
	    dnm += *p_d * *p_d;
	    gnm += *p_gn * *p_gn;
	    dg += *p_d * *p_gn;
	    break;

	  case BOUNDED_LC_FREE:
	  case BOUNDED_UC_FREE:
	    *p_d = - *p_gn;
	    gjsqr = *p_gn * *p_gn;
	    dnm += gjsqr;
	    gnm += gjsqr;
	    dg -= gjsqr;
	    *p_ivec = BOUNDED_NC;
	    break;

	  case BOUNDED_LC_FIXED:
	  case BOUNDED_UC_FIXED:
	    *p_d = 0.0;
	    istuck++;
	    break;

	  case FIXED:
	    break;
	  }

      dnm = sqrt(dnm);
      gnm = sqrt(gnm);
      memcpy(go, gn, nc * sizeof(real));
      memcpy(xo, xn, nc * sizeof(real));
      
      /* Check conditions for restarting the iteration with a gradient step.
      This is done when the angle between the direction vector 
      and the negative gradient is so large that the the direction 
      vector is not "sufficiently downward". cosmin = .01 corresponds
      to an angle between the direction vector and the negative gradient 
      of approximately 89.5 degrees. A restart is also tried if the 
      algorithm is making no progress as measured by the norm of the
      change in x. If no progress is made after this restart, then getsol 
      returns with istop = DX_TOO_SMALL. */

      if((fabs(dg) <  COSMIN * dnm * gnm) || restrt)
	{
	  if (debug)
	    fprintf(stdout, "Restarting\n");

          update_d(nc, &dg, &istuck, gn, d, ivec);

	}
    }

  if (debug)
    fprintf(stdout, "Exceeded iteration limit in getsol: stop.\n");
  *istop = MAX_ITER_LIMIT;
  return ;

}

/***************************************************************************/

void alloc_ccg_mem(int nc)
{
     strcpy(message, "alloc_ccg_mem") ;

     xn = (real *) e_malloc(nc*sizeof(real), message) ;
     gn = (real *) e_malloc(nc*sizeof(real), message) ;
     go = (real *) e_malloc(nc*sizeof(real), message) ;
     d  = (real *) e_malloc(nc*sizeof(real), message) ;
}
  
/***************************************************************************/


/************************************************************************/
/*                                                                      */
/*     Program: clsrch.c                                                */
/*     Date:    2/24/95                                                 */
/*     Author:  Erik M. Johansson, translation of Dennis Goodman's      */
/*              FORTRAN constrained least squares conjugate gradient    */
/*              codes                                                   */
/*                                                                      */
/************************************************************************/

/*-------------------------------------------------------------------------
The search for a value of ap that minimizes the function is done using
safeguarded cubic extrapolation/interpolation.  A bracket is established 
within which a minimizing ap must lie.  The next guess of the minimum is 
obtained by fitting a cubic.  The next guess is the minimum of the cubic
if it lies in a subinterval of the bracket.  Otherwise, the next guess is
one of the endpoints of the subinterval.  To speed convergence, if one end 
of the bracket has stayed the same too long, evaluation at the appropriate
endpoint of the subinterval is forced.  The optimum ap must also be such
that (ap,fvap) lies below the "upper rho line" defined by (ap,fv + rho*dg*ap)
where fv and dg are the initial function value and directional derivatives.
If sigma > rho, the bracket contains acceptable points that lie below the
upper rho line (Armijo-Goldstein condition), and have a directional
derivative that has been reduced by a factor of sigma.  These are the normal
conditions for terminating the linesearch; they are called the Wolfe-Powell
conditions.  Unfortunately the minimum may occurr at a value of alpha where
a variable encounters a bound.  The derivative is discontinuous here, and it
may be impossible to reduce its magnitude by sigma.  To deal with this, after
NGOLD iterations the Goldstein conditions which require that (ap,fvap) lie
below the upper rho line and above the lower rho line defined by 
(ap,fv + (1-rho)*dg*ap) are tested.  If they are satisfied the linesearch 
terminates.  Parameters rho, sigma and tau are as defined in Fletcher's book.
-----------------------------------------------------------------------------*/

#define MAXITN 50
#define NGOLD 5
#define NSAME 3
#define RHO 1.e-4
#define SIGMA 0.1
#define TAU1 9.0
#define TAU2 0.1
#define TAU3 0.5

void clsrch(int nc, real *xo, real *xn, real *gn, real *d, 
	    real *xmin, real *xmax, int *ivec, real *fn, real *dg, 
	    real *ap, real amx, real fmin, real tol, int *itn, 
	    int *iqt, int *istuck)
/************************************************************************

	 on entry: 

	nc	 the number of variables; unchanged by clsrch		
	xo	 current value of the variables; unchanged by clsrch	
	d	 search direction; unchanged by clsrch		
	xmin	 lower bounds on the variables			
	xmax	 upper bounds on the variables			
	ivec	 array of integers providing information about the
		variables.  See ccg.h					
	fn	 function value at xo; On return, value at xn		
	dg	 magnitude of directional derivative at xo; on return,
		its magnitude at xn					
	ap	 initial guess for steplength; on return, steplength	
	amx	 maximum value for steplength				
	fmin	 smallest possible value for function being minimized	
	tol	 relative tolerance for function			

	 on return:

	xn	 new value of the variables; xn = xo + ap*d     
	gn	 the gradient at xn				
	itn	 the number of iterations in getsol	
	iqt	 status of linesearch on return.  See ccg.h
	istuck	 number of stuck positions          	
************************************************************************/

{
  
  int braktd, btr, ingold, belohi, abovlo, debug, debug1;
  real apbtr, apwrs, fnbtr, fnwrs, dgbtr, dgwrs, 
    apdif, fnap, dlt, tp1, atp1, sdgbtr, sdgwrs, dgbw, tp2, scale, top,
    bottom, aplow, aphigh, dgap;
  int izwrs, izbtr, nswrs, nsbtr, iqflag, izap, izdif;

  void funct();

  debug = FALSE;
  debug1 = FALSE;

  /* Check to ensure that the starting guess is positive and that
     the starting directional derivative is negative. */

/* SWITCHING ORDER OF CHECK - 9/18/95 */

  if (*dg >= 0.0)
    {
      if (debug1)
	fprintf(stdout, "LINESEARCH ERROR: starting dg >= 0\n");
      *iqt = INITIAL_DD_GT_0;
      return;
    }
  if (*ap <= 0.0)
    {
      if (debug1)
	fprintf(stdout, "LINESEARCH ERROR: starting ap <= 0\n");
      *iqt = INITIAL_AP_LT_0;
      return;
    }

  /* initialize things. */
  apbtr = apwrs = 0.0;
  braktd = btr = FALSE;
  izwrs = izbtr = *istuck;
  nswrs = nsbtr = -1;
  fnbtr = fnwrs = *fn;
  dgbtr = dgwrs = *dg;
  iqflag = 0;

  /* Main loop. */

  *itn = 0;
  for (;;)
    {
      (*itn)++;

      if (setjmp(sjbuf) != 0) {
	  *iqt = INTERRUPT ;
	  return ;
      }
    
      if ((*itn > MAXITN) || (fabs((apbtr - *ap) * dgbtr) < tol))
	{
	  *ap = apbtr;
	  *fn = fnbtr;
	  *dg = dgbtr;
	  *istuck = izbtr;
	  if (*itn > MAXITN)
	    {
	      /* Maximum iteration count exceeded. */
	      if (debug)
		fprintf(stdout, "LINESEARCH EXIT: iterations exceeded MAXITN\n");
	      if (braktd)
		*iqt = MAX_FN_CALLS_INTERP;
	      else
		*iqt = MAX_FN_CALLS_EXTRAP;
	    }
	  else
	    {
	      /* Tolerance problems, so quit. */
	      *iqt = TOLERANCE_PROBLEMS;
	      if (debug1)
		fprintf(stdout,"LINESEARCH EXIT: tolerance problems\n");
	    }

	  /* For certain abnormal exits, the current xn is not the best.
	     Hence, call update one last time to get xn at the
	     best point.If calling funct is particularly expensive, this
	     extra call could be avoided by an extra vector to keep track
	     of g. However, typically the algorithm won't get here very
	     often.  */
	  if (!btr)
	    {
	      update(nc, *ap, &izap, xo, xn, d, xmin, xmax, ivec); 

              funct(nc, xn, &fnap, gn, &iqflag) ; 

	      *fn = fnap;
	      dgap = dot_product_ivec(nc, gn, d, ivec);
	      *dg = dgap;
	      *istuck = izap;
	      (*itn)++;
	      if (debug)
		fprintf(stdout,
			"LINESEARCH ABNORMAL EXIT: CALL UPDATE AGAIN\n");
	    }
	  return;
	}

      /* Update the function value, and the directional
	 derivative at the new value of alpha. */
      update(nc, *ap, &izap, xo, xn, d, xmin, xmax, ivec); 

      funct(nc, xn, &fnap, gn, &iqflag) ; 

      dgap = dot_product_ivec(nc, gn, d, ivec);
      apdif = apbtr - apwrs;
      izdif = izbtr - izwrs;
      if (debug) {
	fprintf(stdout, "linesearch iteration = %d bracketed is %d apdif = %g izdif = %d\n",
		*itn, braktd, apdif, izdif) ;
	fprintf(stdout, "apw = %g fnw = %g dgw = %g izw = %d\n", 
		apwrs, fnwrs, dgwrs, izwrs) ;
	fprintf(stdout, "apb = %g fnb = %g dgb = %g izb = %d\n", 
		apbtr, fnbtr, dgbtr, izbtr);
	fprintf(stdout, "ap = %g fn = %g dg = %g iz = %d\n",
		*ap, fnap, dgap, izap);
      }
      if (fnap < fmin)
	{
	  /*----------------------------------------------------------
	    fnap is smaller then the expected minimum, so quit.
	    ----------------------------------------------------------*/
	  *fn = fnap;
	  *dg = dgap;
	  *istuck = izap;
	  *iqt = FN_LT_FMIN;
	  if (debug1)
	    fprintf(stdout, "LINESEARCH EXIT: function < fmin\n");
	  return;
	}
      else if (iqflag < 0)
	{
	  /*----------------------------------------------------------
	    funct returned a negative iqflag, so quit.
	    ----------------------------------------------------------*/
	  *fn = fnap;
	  *dg = dgap;
	  *istuck = izap;
	  *iqt = iqflag;
	  if (debug1)
	    fprintf(stdout, "LINESEARCH EXIT: negative iqflag\n");
	  return;
	}

      /* Determine the current endpoints of the bracket and if an
	 acceptable point has been bracketed. */
      ingold = (*itn < NGOLD);
      belohi = fnap <= *fn + RHO * *ap * *dg;
      abovlo = fnap >= *fn + (1.0 - RHO) * *ap * *dg;
      if (ingold && ((!belohi) || (fnap >= fnbtr)))
	{
	  /*----------------------------------------------------------
	    The value of ap was not better than apbtr, so leave apbtr
	    alone and replace apwrs with ap. An acceptable point has
	    been bracketed because ap either caused the graph to cross
	    the upper RHO line or caused an increase.
	    ----------------------------------------------------------*/
	  nswrs = 0;
	  nsbtr++;
	  braktd = TRUE;
	  btr = FALSE;
	  apwrs = *ap;
	  fnwrs = fnap;
	  dgwrs = dgap;
	  izwrs = izap;
	}
      else if (belohi && (fabs(dgap) <= -SIGMA * *dg))
	{
	  /*----------------------------------------------------------
	    ap reduced the directional derivative sufficiently to
	    terminate the linesearch, so return.
	    ----------------------------------------------------------*/
	  apbtr = *ap;
	  *fn = fnap;
	  *dg = dgap;
	  *istuck = izap;
	  if (debug1)
	    fprintf(stdout, "LINESEARCH EXIT: standard conditions\n");
	  iqt = NORMAL_STOP;
	  return;
	}
      else if ((!ingold) && belohi && abovlo)
	{
	  /*----------------------------------------------------------
	    After NGOLD or more function evaluations, we still have no
	    satisfactory point. This is probably due to the fact that
	    we have bracketed a point with a derivative discontinuity
	    at which the derivative changes sign. Therefore, check to
	    see if the current best point at least is above the
	    lower line in the Goldstein conditions, and return if
	    this is indeed the case.
	    ----------------------------------------------------------*/
	  if (debug1)	
	    fprintf(stdout, "LINESEARCH EXIT: satisfied Goldstein\n"); 
	  apbtr = *ap;
	  *fn = fnap;
	  *dg = dgap;
	  *istuck = izap;
	  iqt = NORMAL_STOP;
	  return;
	}
      else
	{
	  /*----------------------------------------------------------
	    ap wasn't good enough to stop the linesearch, but it is 
	    still the current best point. First see what should be
	    done with the other end of the bracket, then replace apbtr
	    with ap.
	    ----------------------------------------------------------*/
	  if ((!braktd) || (((apwrs - apbtr) * dgap) >= 0.0))
	    {
	      /*-------------------------------------------------------
		Replace apwrs with the current apbtr. If have not yet 
		bracketed, check if the derivative has changed sign
		(i.e. become positive) in which case desirable points
		have been bracketed.
		-------------------------------------------------------*/
	      apwrs = apbtr;
	      fnwrs = fnbtr;
	      dgwrs = dgbtr;
	      izwrs = izbtr;
	      nswrs = -1;
	      if (dgap >= 0.0)
		braktd = TRUE;
	    }
	  btr = TRUE;
	  apbtr = *ap;
	  fnbtr = fnap;
	  dgbtr = dgap;
	  izbtr = izap;
	  nswrs++;
	  nsbtr = 0;
	}

      if (debug)
	fprintf(stdout, "nswrs = %d nsbtr = %d\n", nswrs, nsbtr);

      /* Perform Hermite cubic interpolation, calculate minimum of the cubic */
      dlt = apbtr - apwrs;
      tp1 = dgwrs + dgbtr - 3.0 * (fnwrs - fnbtr) / (apwrs - apbtr);

      /* tp2 is sqrt(tp1**2 - dgwrs*dgbtr) if tp1**2 - dgwrs*dgbtr > 0.
	 Otherwise, tp2 = 0. This is the square root of a sum of squares,
	 so do it carefully to avoid underflow and overflow. */
      atp1 = fabs(tp1);
      sdgbtr = SIGN(sqrt(fabs(dgbtr)), dgbtr);
      sdgwrs = SIGN(sqrt(fabs(dgwrs)), dgwrs);
      dgbw = sdgbtr * sdgwrs;
      if (dgbw >= 0.0)
	{
	  /* subtracting two squares: */
	  if (atp1 <= dgbw)
	    tp2 = 0.0;
	  else
	    tp2 = sqrt(atp1 + dgbw) * sqrt(atp1 - dgbw);
	}
      else
	{
	  /* adding two squares:*/
	  scale = atp1 - dgbw;
	  tp2 = scale * sqrt((atp1 * atp1) / (scale * scale) +
			     (dgbw * dgbw) / (scale * scale));
	}
	
      if (apwrs > apbtr)
	   tp2 = -tp2 ; 		/* Correction 12/23/97 */

      /* the cubic interpolation estimate is:
	 apcube = apbtr-(apbtr-apwrs)*(dgbtr+tp2-tp1)/(dgbtr-dgwrs+two*tp2)
	 but we don't do it this way because unpleasant things can happen
	 if the denominator is small. Instead, we compute the numerator
	 and denominator separately, and check first if the estimate falls
	 within an acceptible region. If not, we don't calculate it. We
	 multiply inequalities by bottom, so we need bottom >=  0. */

      top = (apbtr - apwrs) * (dgbtr + tp2 - tp1);
      bottom = dgbtr - dgwrs + 2.0 * tp2;
      if (bottom < 0.0)
	{
	  top *= -1;
	  bottom *= -1;
	}

      if (!braktd)
	{
	  /*----------------------------------------------------------
	    Extrapolate.
	    ----------------------------------------------------------*/
	  if (apbtr == amx)
	    {
	      /*-------------------------------------------------------
		Return because at the maximum allowable ap without
		bracketing any acceptable points.
		-------------------------------------------------------*/
	      *fn = fnap;
	      *dg = dgap;
	      *istuck = izap;
	      *iqt = AP_EQ_AMX;
	      if (debug1)
		fprintf(stdout, "LINESEARCH EXIT: at amx\n");
	      return;
	    }
	  else if ((apbtr + dlt) >= amx)
	    /*-------------------------------------------------------
	      Maintain the upper limit amx.
	      -------------------------------------------------------*/
	    *ap = amx;
	  else
	    /*-------------------------------------------------------
	      Make certain that the next ap is moved to the right by
	      at least the length of the current bracket, but by no
	      more than the TAU1 times this length. Hence, apcube
	      is accepted if and only if:
	      apbtr + dlt < apcube < min(amx,apbtr + TAU1*dlt)
	      -------------------------------------------------------*/
	    if ((dlt * bottom) >=  -top)
	      *ap = apbtr + dlt;
	    else if (-top >= (bottom * MIN(amx - apbtr, TAU1 * dlt)))
	      *ap = MIN(amx, apbtr + TAU1 * dlt);
	    else
	      *ap = apbtr - top / bottom;
	}
      else
	{
	  /*----------------------------------------------------------
	    Apply safeguarded interpolation. The new ap is apcube if
	    apcube lies between apbtr-TAU2*dlt and apwrs+TAU3*dlt.
	    Otherwise ap is set to the one of these points closest to
	    apcube.
	    ----------------------------------------------------------*/
	  if (apbtr < apwrs)
	    {
	      aplow = apbtr - TAU2 * dlt;
	      aphigh = apwrs + TAU3 * dlt;
	    }
	  else
	    {
	      aplow = apwrs + TAU3 * dlt;
	      aphigh = apbtr - TAU2 * dlt;
	    }
	  if (((nsbtr / NSAME) * NSAME == nsbtr) && (nsbtr != 0))
	    {
	      if (debug1)
		fprintf(stdout, "FORCING EVALUATION NEAR APBTR\n");
	      *ap = apbtr - TAU2 * dlt;
	    }
	  else if (((nswrs / NSAME) * NSAME == nswrs) && (nswrs != 0))
	    {
	      if (debug1)
		fprintf(stdout, "FORCING EVALUATION NEAR APWRS\n");
	      *ap = apwrs + TAU3 * dlt;
	    }
	  else if (bottom * (aplow - apbtr) >= -top)
	    *ap = aplow;
	  else if (bottom * (aphigh - apbtr)<= -top)
	    *ap = aphigh;
	  else
	    *ap = apbtr - top / bottom;
	}
    } /* end of for (;;) loop */
}

void update(int nc, real ap, int *itoss, real *xo, real *xn, real *d, 
	    real *xmin, real *xmax, int *ivec)
{
  int debug, n, *p_ivec;
  real *p_xn, *p_xo, *p_d, *p_xmin, *p_xmax;

  debug = FALSE;
  *itoss = 0;

  /* ----------------------------------------------------------------
     Update the residual at the new value of alpha 
     and calculate the function and directional derivative. 
     ----------------------------------------------------------------*/

  if (debug)
    countem(nc, ivec);

  p_ivec = ivec;
  p_xn = xn;
  p_xo = xo;
  p_d = d;
  p_xmin = xmin;
  p_xmax = xmax;
  for (n = 0; n < nc; n++, p_ivec++, p_xn++, p_xo++, p_d++, p_xmin++, p_xmax++)
    switch (*p_ivec)
      {
      case BOUNDED_NC:
      case BOUNDED_LC_FREE:
      case BOUNDED_UC_FREE:
	   
	    /*-------------------------------------------------------------
	      path for bounded variables that are free to move but are
	      either at their bounds or may hit
	      their bounds during this iteration. they are moved the
	      appropriate distance and forced to be within their bounds.
	      ivec[n] is set to BOUNDED_NC if the nth variable is strictly
	      within its bounds, to BOUNDED_LC_FREE if it is at the lower 
	      bound, and to BOUNDED_UC_FREE if it is at the upper bound.
	      -------------------------------------------------------------*/
	    *p_xn = *p_xo + ap * *p_d;

	    if (*p_xn <= *p_xmin)
	      {
		*p_xn = *p_xmin;
		*p_ivec = BOUNDED_LC_FREE;
		(*itoss)++;
	      }
	    else if (*p_xn >= *p_xmax)
	      {
		*p_xn = *p_xmax;
		*p_ivec = BOUNDED_UC_FREE;
		(*itoss)++;
	      }
	    else
		*p_ivec = BOUNDED_NC;
	 
        break;

      case UNCONSTRAINED:
	/*-------------------------------------------------------------
	  path for the unbounded variables. they are moved the
	  appropriate distance. Also path for bounded variables
	  that are not at bounds and will not hit them during
	  this iteration.
	  -------------------------------------------------------------*/
	*p_xn = *p_xo + ap * *p_d;
	break;
	
      case BOUNDED_LC_FIXED:
      case BOUNDED_UC_FIXED:
	    /*-------------------------------------------------------------
	      path for bounded variables that are fixed for this
	      iteration of the algorithm.  Nothing is done to them.
	      The xn = xo statement should have no effect.
	      -------------------------------------------------------------*/
	    *p_xn = *p_xo;
	    (*itoss)++;
	    break;

      case FIXED:
	    /*-------------------------------------------------------------
	      path for the fixed variables. nothing is done to them.
	      The xn = xo statement should have no effect.
	      -------------------------------------------------------------*/
	    *p_xn = *p_xo;
	    break;

      }	

  if (debug)
    countem(nc, ivec);

}

void update_d(int nc, real *dg, int *istuck, real *gn, real *d, int *ivec)
{

  int	j, *p_ivec ;
  real *p_gn, *p_d ;

  /* Calculate initial d = -g+ and dg = |g+|^2, initialize ivec 
  istuck counts the number of bounded variables at upper or lower 
  bounds that will not move in next linesearch iteration (i.e., 
  those variables with ivec[j] = BOUNDED_LC_FIXED or BOUNDED_UC_FIXED) */

  *dg = 0.0;
  *istuck = 0;
  p_gn = gn;
  p_d = d;
  p_ivec = ivec;
  for (j = 0; j < nc; j++, p_gn++, p_d++, p_ivec++)
    {
      switch (*p_ivec)
	{
	case BOUNDED_NC:
	case UNCONSTRAINED:
	  *p_d = - *p_gn;
	  *dg -= *p_gn * *p_gn;
	  break;
	case BOUNDED_LC_FREE:
	case BOUNDED_LC_FIXED:
	  if (*p_gn < 0.0)
	    {
	      *p_d = - *p_gn;
	      *dg -= *p_gn * *p_gn;
	      *p_ivec = BOUNDED_NC;
	    }
	   else
	     {
	       *p_d = 0.0;
	       *p_ivec = BOUNDED_LC_FIXED;
	       (*istuck)++;
	     }
	  break;
	case BOUNDED_UC_FREE:
	case BOUNDED_UC_FIXED:
	  if (*p_gn > 0.0)
	    {
	      *p_d = - *p_gn;
	      *dg -= *p_gn * *p_gn;
	      *p_ivec = BOUNDED_NC;
	    }
	  else
	    {
	      *p_d = 0.0;
	      *p_ivec = BOUNDED_UC_FIXED;
	      (*istuck)++;
	    }
	  break;
	case FIXED:
	  *p_d = 0.0;
	  break;
	}
    }
}



real dot_product_ivec(int nv, real *g, real *d, int *ivec)
{
  real *p_g, *p_d, sum;
  int *p_ivec, i;
  
  sum = 0.0;
  p_ivec = ivec;
  p_g = g;
  p_d = d;
  for (i = 0; i < nv; i++, p_g++, p_d++, p_ivec++)
    switch (*p_ivec)
      {
      case BOUNDED_NC:
      case UNCONSTRAINED:
	sum += *p_g * *p_d;
      default:
	break;
      }
  return sum;
}

/****************************************************************************/

void countem(int nv, int *ivec)
{

  int itype[5], *p_ivec, i;

  for (i = 0; i < 5; i++)
    itype[i] = 0;

  p_ivec = ivec;
  for (i = 0; i < nv; i++, p_ivec++)
    itype[*p_ivec]++;

  fprintf(stdout, "i1 = %-5d i2 = %-5d i3 = %-5d i4 = %-5d i5 = %-5d\n",
	  itype[0], itype[1], itype[2], itype[3], itype[4]);
}

void handle_ctrl_c(int dummy)   
{
    signal(SIGINT,handle_ctrl_c) ;
    longjmp(sjbuf, 0) ;
}


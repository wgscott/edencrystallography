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

                                UTIL.H

  Title:  Include file needed by the full suite of Eden programs
  Author:  Hanna Szoke
  Date:  3/3/93

******************************************************************************/

#include <stdlib.h>	/* declarations of malloc() et al. */
#include <stdio.h>
#include <math.h>
#include <ctype.h>	/* declarations of isalpha() etc. */
#include <string.h>	/* declarations of strcpy etc. */
#include <unistd.h>	/* Added 3/2/99, for Linux getopt() */
#include <sys/types.h>
#include <sys/file.h>
#include <sys/stat.h>

#define MAXSTRING       400     /* max length for reading character strings */
#define	UNSET		-1	/* initialization for reading input	*/

#define FALSE            0
#define TRUE             1

#define	SQRT2	sqrt((double)2.) 
#define RTOD	180. / M_PI
#define DTOR	M_PI / 180.
#define PISQ	M_PI*M_PI
#define TWOPI	2*M_PI
#define EPS	pow(2.0, -52.)
#define	BtoMB	pow(2.0, -20.) 

#define	COMMENT_CHAR	'#'	/* can be changed if you like! */

#ifdef DOUBLE
     typedef	double	real;
     #include "dfftw.h"
     #define VERSION "EDEN V5.3 - Double Precision:"
#else
     typedef	float	real ;
     #include "sfftw.h"
     #define VERSION "EDEN V5.3 - Single Precision:"
#endif

typedef struct {
  real re, im;
} COMPLEX;

typedef struct {
   real       x ;
   real       y ;
   real       z ;
} POINT ;

/***********************************************************
 * GLOBAL VARIABLES
 ***********************************************************/
extern	FILE    *fp_log ;	/* pointer to log file */
extern	FILE    *fp_cost ;	/* pointer to cost file */
extern	char	log_filename[] ;	/* name of log file */
extern	char	caller[] ;	/* name of calling program (eden, back, etc.) */
extern	char	command_line[] ;	/* command line to eden */
extern	char	message[] ;	/* general-purpose message buffer */
extern	char	*allinp ;	/* buffer for all input */
extern	char	nextline[] ;	/* ...for reading ASCII files	*/
extern	char	goodline[] ;	/* for identifying useful lines of input */
extern	char	terminp[] ;	/* used by prompt() to capture user input */
extern	char	id[] ;		/* first symbol on a line	*/
extern	int	Nlines ;	/* actual # of lines in .inp file */
extern	char	input_filename[] ;	/* parameter file name */
extern	char	sf_filename[] ;	/* generic structure factor file name */
extern	char	pwd_name[] ;	/* name of str. factor file, w/o path */
extern	char	out_filename[] ;	/* generic output file name	*/
extern	char	mode_info[] ;	/* operating mode info  */
extern	int	verbose ;	/* execute-line option -v */
extern	int	very_verbose ;	/* execute-line option -V */
extern	int	interactive ;	/* execute-line option -i */
extern	int	graphics ;	/* execute-line option -g */
extern	int	quick ;		/* execute-line option -q */
extern	int	silent ;	/* execute-line option -s */
extern	char	costfile_option ;	/* save costs? (goes with -v)	*/
extern	float	t0, t1, t2, t3, t4, t5, t6 ;	/* for reading w. sscanf() */


/***********************************************************
 * Commonly-used utility functions
 ***********************************************************/
extern	void	ballpark(char *) ;
extern	void	check_exist(char *) ;
extern	void	hello(char *) ;
extern	void	printTwice(char *) ;
extern	void	EdenError(char *) ;
extern	void	EdenWarning(char *) ;
extern	void	prompt(char *) ;
extern	void	readBasicInput() ;
extern	void	strip_path(char *, char *) ;
extern	void	strip_suffix(char *, char *, char *) ;
extern	void	extend_filename(char *, char *, char *) ;
extern	int	net_strlen(char *) ;		
extern	char	*timestamp() ;	
extern	void	*e_malloc(unsigned int, char *) ;
extern	void	e_free(unsigned int, void *) ;
extern	void    getAmpPhase(COMPLEX *, real, real *, real *) ;

void matvecmult(real mat[3][3], real vec[3], real *prod);
void matinv(real m[3][3], real minv[3][3]);

/***********************************************************
 * The EDEN methods
 ***********************************************************/
typedef void eden_method(int,char *[]);

eden_method
  addmaps_main,
  apodfc_main,
  apodfo_main,
  back_main,
  bin2map_main,
  bin2view_main,
  cadhkl_main,
  count_main,
  distance_main,
  dphase_main,
  expandfc_main,
  expandfo_main,
  forth_main,
  map2bin_main,
  maketar_main,
  multmaps_main,
  perturbhkl_main,
  ranphase_main,
  refine_main,
  regrid_main,
  shapes_main,
  solve_main,
  sym_main,
  tohu_main,
  variance_main,
  view2bin_main;


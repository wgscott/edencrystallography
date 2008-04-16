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

				UTIL.C

  Title:        Utility function package for EDEN 
  Author:       Hanna Szoke
  Date:         10/19/92
  Function:     This package contains functions used by Solve and 
		by all its satellite executable programs.  They include 

		String manipulation functions -
			net_strlen()
			strip_path()
			strip_suffix()
			first3chars()
		Message traffic functions -
			help()
			ballpark()
			hello()
			nice_print_to_log()
			printTwice()
			prompt()
			EdenWarning()
			EdenError()
		File manipulation functions -
			extend_filename()
			check_exist()
			filelength()
		Miscellaneous functions -
			timestamp()
			clip()
			fclip()
			Cabs()
			hpsort()
		Matrix manipulation functions -
			matmatmult()
			matvecmult() 
			matinv() 
			printmat()
			printvec()
		Random number functions -
			ran1()
			normGauss()
		Memory management functions -
			e_malloc()
			e_free()
		Record keeping functions -
			setup_records() 
			fatal_error_signal ()
			start_record() 
			end_record() 

*******************************************************************************/
#include  <sys/time.h>	/* ... for picking up date & time	*/
#include  <signal.h>	/* ... for capturing kill		*/
#include  "util.h"

char	message[MAXSTRING] ;
char	terminp[MAXSTRING] ;	/* used by prompt() to capture user input */
char	out_filename[MAXSTRING] ;

static	char	e_message[MAXSTRING] ;	/* used when there's a conflict
					with message[] */
static	char	record_filename[MAXSTRING] ="" ;
static	char	*record ;		/* array holding record info */
static	int	rec_length = 0 ;	/* current length of record */

static	real	ememory = 0 ;
static	real	maxmem = 0 ;

static	void	fatal_error_signal(int);
static	volatile sig_atomic_t fatal_error_in_progress = 0;
     

				/****************************************
				String manipulation functions
				****************************************/

int	net_strlen(char text[])	/*  strip comments, '=' and ',' */
{
     int	 i, m, jb ;
     char	*cp ;

     if (((m = (int)strlen(text)) > 0) && 
	 ((cp = strchr(text, COMMENT_CHAR)) != NULL)) {

	       while ((text+m) >= cp) {
		    text[m] = '\0' ;
		    m-- ;
               }
     }
     if ((m = (int)strlen(text)) > 0) { 	/* no trailing commas ... */
               if (text[m] == ',') {
		    text[m] = '\0' ;
		    m-- ;
               }
     }
     if ((m = (int)strlen(text)) > 0) { 	/* ... and no equal sign */
          for (i = 0; i < m; i++) {
               if (text[i] == '=') 
		    text[i] = ' ' ;
          }
     }
                    
     if ((m = (int)strlen(text)) > 0) { 	/* ... and no leading blanks. */
	  for (i = jb = 0; i < m; i++) {
	       if (isalpha(text[i]))
		    break ;
               else
		    jb++ ;
          }
	  for (i = 0; i < m-jb; i++) 
               text[i] = text[i+jb] ;

	  for (i = m-jb; i < m; i++) 
               text[i] = '\0' ;

     }

     m = (int)strlen(text) ;

     if (very_verbose)
	  fprintf(stderr, "line = %s, length = %d\n", text, m) ;

     return (m) ;
}

void	strip_path(char *fullpath, char *filename)
{
	/* Strip file name from potentially full path name */

     char	c = '/' ;
     char	*p ;
     int	i ;

     if ((int)strlen(fullpath) == 0)
	  EdenError("Missing filename for strip_path!") ;

     for (i = 0; i < MAXSTRING; i++)
          *(filename+i) = '\0' ;
	   
     if ((p = strrchr(fullpath, c)) == NULL) 
          strcpy(filename, fullpath) ;
     else
          strncpy(filename, p+1, (int)strlen(p) - 1);
}

void	strip_suffix(char *oldname, char *suffix, char *newname)
{
     int	i, l ;

	/* Strip suffix from potentially full name */

     strcpy(newname, oldname) ;

     if (strstr(oldname, suffix) != NULL) {

	  l = strlen(newname) ;
	  for (i = 1; i < (int)strlen(suffix)+1; i++)
	       *(newname+l-i) = '\0' ;
     }
}

int first3chars(char *line, char *symbol) 
{
     int	k0, k ;
     int	kmax ;

     id[0] = '\0' ;

     sprintf(id, "%s", line) ;
     kmax = (int)strlen(id) ;

     /*****************************************************
     Access first 3 non-blank characters in id, converted 
     to upper-case.  Compare them to "symbol".
     *****************************************************/

     for (k0=0; k0<kmax; k0++)
	  if (isalpha(id[k0]))
	       break ;

     for (k=0; k<kmax; k++)
	 id[k] = toupper((int) id[k+k0]) ;

     return (strncmp(id, symbol, 3)) ;
}
     
				/****************************************
				Message Traffic functions
				****************************************/

void	help(char *name) 
{
     sprintf(message, "$EDENHOME/help/lookat $EDENHOME/help/%s.hlp %s &", 
	     name, name) ;
     if (system(message) != 0) {
             sprintf(message, "more $EDENHOME/help/%s.hlp", name) ; 
             system(message) ;
          } 
}

void	ballpark(char *name) 
{
     sprintf(message, "more %s/%s/%s.inf", "$EDENHOME", "help", name) ;
     if (system(message) != 0)
	  ballpark("eden") ;
     exit(-1) ;
}

void	hello(char *caller)
{
     printTwice(" ") ;
     printTwice(timestamp()) ;
     sprintf(message, "\n\n\t%s %s\n", VERSION, caller) ;
     printTwice(message) ;
}

void	printTwice(char *mess)	/* send info to terminal (stdout) and to log */
{

     fprintf(fp_log, "%s\n", mess) ;
     fflush(fp_log) ;
     if (!silent) {
          printf("%s\n", mess) ;
          fflush(stdout) ;
     }
}
 
void	prompt(char *mess)	/* Send message to user, await reply */
{
    fprintf(stdout, mess) ;

    while (fgets(terminp, MAXSTRING, stdin) != NULL)  {
	if ((int)strlen(terminp) > 0) 
            break ;
    }
}

void	nice_print_to_log(char *string1, char *string2)
{
    int p1, p2 ;

    p1 = (int)strlen(string1) ;
    p2 = (int)strlen(string2) ;

    if ((p1+p2) < 80)
          fprintf(fp_log, "%s %s\n", string1, string2) ;
     else {
          fprintf(fp_log, "%s\n", string1) ;
          fprintf(fp_log, "%s\n", string2) ;
     }
}

void    EdenWarning(char *mess)         /* send info to stderr, to log, 
					DON'T quit. */
{
     fprintf(fp_log, "%s\n", mess) ;
     fflush(fp_log) ;
     fprintf(stderr, "%s\n", mess) ;
     fflush(stderr) ;
}
 
void    EdenError(char *mess)		/* send info to stderr, to log, and quit. */
{
     void	end_record(int);
  
     fprintf(fp_log, "ERROR: %s\n", mess) ;
     fflush(fp_log) ;
     fprintf(stderr, "ERROR: %s\n", mess) ;
     fflush(stderr) ;

     printTwice("\n") ;
     printTwice(timestamp()) ;
     printTwice("\n") ;

     if (strcmp(record_filename, "") != 0)
          end_record(1) ;
     else
          exit (-1);
}
				/****************************************
				File manipulation functions
				****************************************/

void extend_filename(char *oldname, char *newname, char *addstring)
{
     int	n ;

     if (strrchr(oldname, '.') == NULL) 
	  sprintf(newname, "%s%s", oldname, addstring) ;
     else {
          n = (int)strlen(oldname) - (int)strlen(strrchr(oldname, '.')) ;
          strncpy (newname, oldname, n) ;
          sprintf (newname+n, "%s%s", addstring, oldname+n) ;
     }

     return ;
}
     
void check_exist(char *name)
{
     FILE	*fp ;

     if ((fp = fopen(name, "r")) == NULL)
     {
          sprintf(e_message, "Cannot open %s", name);
          EdenError(e_message) ;

     /* note use of e_message to prevent conflict in 
        cases where 'name' is encoded in 'message'. */
     }
     fclose(fp) ;
}

int	filelength(char *name)
{
     int	length = 0 ;
     FILE	*fpin ;

     sprintf(message, "wc %s >tmp", name) ;

     system(message) ;
     if ((fpin = fopen("tmp", "r")) != NULL) {
          fscanf(fpin, "%d", &length) ;
          fclose (fpin) ;
          system("rm tmp") ;
     }
     return (length) ;

}
				/****************************************
				Miscellaneous functions 
				****************************************/

char	*timestamp()		/* Put time into string 		*/
{
     struct	tm *localtime(const time_t *) ;
     char	*asctime(const struct tm *) ;
/***    time_t	time(char *) ;***/
     time_t	t ;
     char	*c ;

     t = time(NULL) ;
     c = (char *)asctime(localtime(&t)) ;

     *(c+(int)strlen(c)-1) = '\0' ;
     return (c) ;
}	     

int	clip(int value, int range)	/* Clip integers to (0, range) */
{
     while (value < 0)
          value += range ;
     while (value >= range)
          value -= range ;

     return (value) ;
}

real	fclip(real value, real range)
{
     while (value < 0)
          value += range ; 
     while (value >= range)
          value -= range ;

     return (value) ;
}
real	Cabs(COMPLEX cx)		/* abs (modulus) of complex argument */
{ 
     return (sqrt(cx.re*cx.re + cx.im*cx.im)) ;
}

void hpsort(unsigned long n,float ra[])	
			/* Numerical Recipes, 2nd edition, K & R version,
			painfully rewritten for zero-origin arrays by HS. */
{
	unsigned long i,ir,j,l;
	float rra;

        if (n < 2) return ;
	l=(n >> 1)+1;
	ir=n-1;
	for (;;) {
		if (l > 1) {
			rra=ra[--l-1];
		} else {
			rra=ra[ir];
			ra[ir]=ra[0];
			if (--ir == 0) {
				ra[0]=rra;
				break;
			}
		}
		i=l-1;
		j=l+l-1;
		while (j <= ir) {
			if (j < ir && ra[j] < ra[j+1]) j++;
			if (rra < ra[j]) {
				ra[i]=ra[j];
				i=j;
				j = (2*j+1);
			} else j=ir+1;
		}
		ra[i]=rra;
	}
}

				/****************************************
				Matrix manipulation functions
				****************************************/

/*********************************************************************
	Worker functions for 3 by 3 matrices ONLY! - 
		matmatmult() multiplies 2 matrices,
	   	matvecmult() multiplies a matrix and a vector,
		matinv() inverts a matrix, and
		printmat() printvec() provide diagnostics.
*********************************************************************/
void matmatmult(real mat1[3][3], real mat2[3][3], real *prod)
{
     int	i, j, k ;

     for (i = 0; i < 3; i++) {
	  for (j = 0; j < 3; j++, prod++) {
	       *prod = 0 ;
	       for (k = 0; k < 3; k++)
		    *prod += mat1[i][k]*mat2[k][j] ;
          }
     }
}
void matvecmult(real mat[3][3], real vec[3], real *prod)
{
     int	j, k ;

     for (j = 0; j < 3; j++, prod++) {
          *prod = 0 ;
          for (k = 0; k < 3; k++) 
   	    *prod += mat[j][k]*vec[k] ;
     }
}
void matinv(real m[3][3], real minv[3][3])
{
     real	determ ;
     real	*p_minv ;
     int	j ;

     /* calculate cofactors */

     minv[0][0] = m[1][1]*m[2][2] - m[1][2]*m[2][1] ;
     minv[1][0] = m[1][2]*m[2][0] - m[1][0]*m[2][2] ;
     minv[2][0] = m[1][0]*m[2][1] - m[1][1]*m[2][0] ;

     minv[0][1] = m[2][1]*m[0][2] - m[2][2]*m[0][1] ;
     minv[1][1] = m[2][2]*m[0][0] - m[2][0]*m[0][2] ;
     minv[2][1] = m[2][0]*m[0][1] - m[2][1]*m[0][0] ;

     minv[0][2] = m[0][1]*m[1][2] - m[0][2]*m[1][1] ;
     minv[1][2] = m[0][2]*m[1][0] - m[0][0]*m[1][2] ;
     minv[2][2] = m[0][0]*m[1][1] - m[0][1]*m[1][0] ;

     /* calculate determinant */

     determ = m[0][0]*minv[0][0] + m[0][1]*minv[1][0] + m[0][2]*minv[2][0] ;

     /* normalize */

     for (j = 0, p_minv = &minv[0][0]; j < 9; j++, p_minv++)
          *p_minv /= determ ;
}
void	printmat(char *string, real mat[3][3]) 
{
     int	i ;

     for (i = 0; i < 3; i++) {
	  sprintf (message, "%s, row %d: %f \t%f \t%f", string, i, 
                  mat[i][0], mat[i][1], mat[i][2]) ;
          printTwice(message) ;
     }

     printTwice("") ;
}
void	printvec(char *string, real vec[3]) 
{
     sprintf(message, "%s:\t%g\t%g\t%g\n", string, vec[0], vec[1], vec[2]) ;
     printTwice(message) ;
}

				/****************************************
				Random number functions
				****************************************/

/*******************************************************************************

  "Minimal" random number generator of Park and Miller with Bays-Durham shuffle 
  and added safeguards.  Returns a uniform random deviate between 0.0 and 1.0 
  (exclusive of the endpoint values).  Call with 'idum' a negative integer to 
  initialize; thereafter, do not alter 'idum' between successive deviates in a 
  sequence.  rnmx should approximate the largest floating value that is less 
  than 1.

  Ref: Press & Teukolsky, "Portable Random Number Generators", Computers in 
  Physics, Vol. 6, #5, p. 522  (Sept/Oct 1992).

*******************************************************************************/

#define	NTAB	32	/* table length */

float ran1(int *idum)
{
    int	 	ia = 16807, im = 2147483647, iq = 127773 ;
    int		ir = 2836, ndiv = 1 + (im-1)/NTAB ;
    float	eps = 1.2e-7, rnmx = 1. - eps ;
    float	am = 1. / im ;
    float	val ;
    static int 	iv[NTAB], iy = 0 ;
    int 	j, k ;

/*  Initialize; ensure idum != 0; load shuffle table after 8 warm-ups. */

      if (*idum <= 0 || iy == 0) {
	   *idum = -*idum ;
	   if (*idum < 1) *idum = 1;
	   for (j = NTAB+7; j >= 0; j--) {
		k = *idum/iq ;
		*idum = ia*(*idum-k*iq)-ir*k ;	
		if (*idum < 0) *idum += im ;
		if (j < NTAB) iv[j] = *idum ;
           }
	   iy = iv[0] ;
      }

/*  Start here when not initializing; compute idum = mod(ia*idum, im)
    without overflow, by Schrage's method.  j will be in range 0:NTAB-1. */

      k = *idum/iq ;
      *idum = ia*(*idum-k*iq)-ir*k ;
      if (*idum < 0) *idum += im ;
      j = iy/ndiv ;

/*  Output previously stored value and refill shuffle table. */

      iy = iv[j] ;
      iv[j] = *idum ;
      val = am*iy ;
      if (val == 1.) 	/*  Ensure that endpoint value is not reached. */
      {		
	   printTwice("Hit an endpoint in the random number generator!") ;
	   val = rnmx ;
      }

      return (val) ;
}
real	normGauss(int *pseed)
{
		/* Function adapted from Numerical Recipes (gasdev) 
		returns a normally distributed deviate with zero
		mean and unit variance.  */

	float	u1, u2 ;
	real	v1, v2, s ;

	do {
	     u1 = ran1(pseed) ;
	     u2 = ran1(pseed) ;
	     v1 = 2 * u1 - 1 ;
	     v2 = 2 * u2 - 1 ;
	     s = v1*v1 + v2*v2 ;
        } while (s >= 1) ;

	return(v1 * sqrt(-2*log(s)/s)) ;
}

				/****************************************
				Memory management functions
				****************************************/

void	*e_malloc(unsigned int length, char *func_name)
{
        void *p ;

        if ((p = (void *)calloc(length, 1)) == NULL) { 
	     sprintf(e_message, "Allocated %g Mbytes; new request: %g Mbytes.", 
	             ememory*BtoMB, length*BtoMB) ;
	     EdenWarning(e_message) ;
             sprintf(e_message, "Problem allocating space from %s", func_name) ;
             EdenError(e_message) ;
        } 
	ememory += length ;
	if (ememory > maxmem)
		maxmem = ememory ;

	if (very_verbose) {
             sprintf(e_message, 
		     "Allocated memory: %d bytes from %s", length, func_name) ;
	     printTwice(e_message) ;
             sprintf(e_message, 
		     "Total: %g Mbytes, Maximum: %g Mbytes", 
		     BtoMB*ememory, BtoMB*maxmem) ;
	     printTwice(e_message) ;
        } 

	return(p) ;
}
void	e_free(unsigned int length, void *p) 
{
        free(p) ;
        p = NULL ;     

	ememory -= length ;

	if (very_verbose) {
             sprintf(e_message, "Deallocated memory: %d", length) ;
	     printTwice(e_message) ;
        } 

	return ;
}    

/*******************************************************************************
	Record keeping functions: Note the record file name and allocate space
	for collecting information.  Enable capturing of ctrl-c and various 
	"kill" commands.  As of 6/30/00, this has been tested and works under
	Linux; in particular, control-c is SIGINT.  The others may all be
	activated by running "kill" with suitable argument.  Note: on another 
	system one should first invoke "stty -a" to determine what's what.  
	Under Linux, for example, control-z causes suspension, 
*******************************************************************************/

void	setup_records(char *record_info) 
{
     int	rec_lines = 5 ;		/* # lines written for bookkeeping */

     strcpy(record_filename, record_info) ;
     record = (char *) e_malloc(rec_lines*MAXSTRING*sizeof(char), 
                                "setup_records") ;

     signal(SIGHUP,fatal_error_signal);
     signal(SIGQUIT,fatal_error_signal);
     signal(SIGKILL,fatal_error_signal);
     signal(SIGTERM,fatal_error_signal);
     signal(SIGINT,fatal_error_signal);

     return;
}


void	fatal_error_signal (int sig)
{
     void	end_record(int) ;

	  /********************************************************************
          Since this handler is established for more than one kind of signal, 
          it might still get invoked recursively by delivery of some other 
          kind of signal.  Use a static variable to keep track of that. 
	  ********************************************************************/

       if (fatal_error_in_progress)
         raise (sig);
       fatal_error_in_progress = 1;
     
       /* Now do the clean up actions: */

          end_record(2) ;
     
	  /********************************************************************
          Now reraise the signal.  We reactivate the signal's default handling, 
	  which is to terminate the process.  We could just call `exit' or 
	  `abort', but reraising the signal sets the return status from the 
	  process correctly. 
	  ********************************************************************/

       signal (sig, SIG_DFL);
       raise (sig);
    }

					/**********************************
					If bookkeeping is invoked, write
					the start-up message to it. 
					***********************************/
void	start_record() 
{
     char	*cwd ;
     int	k ;

     rec_length = k = sprintf(record, "\n") ;
     rec_length += sprintf(record+k, timestamp()) ;
     k = rec_length ;

     if ((cwd = getcwd(NULL, 120)) == NULL) 
	  EdenError("Can't get pwd") ;

     rec_length += sprintf(record+k, "\nRunning EDEN from %s\n", cwd) ;
     k = rec_length ;
     free(cwd) ;
     rec_length += sprintf(record+k, "Command Line: %s\n", command_line) ;
}
					/**********************************
					If bookkeeping is invoked, write 
					the record messages (incl. how Eden 
					ended) to the record file and exit.
					***********************************/
void	end_record(int flag) 
{
     FILE	*fprec ;
     int	k ;

     k = rec_length ;
     if (flag == 0)
          rec_length += sprintf(record+k, "Eden ended successfully.\n") ;
     else if (flag == 1)
          rec_length += sprintf(record+k, 
                        "Eden ended with an error - check %s.\n", log_filename) ;
     else
          rec_length += sprintf(record+k, "Eden killed.\n") ;

     if ((fprec = fopen(record_filename, "a+")) == NULL)
     {
          sprintf(message, "Cannot open %s", record_filename);
          printTwice(message) ;
          exit(-1) ;
     }
     if (fwrite(record, sizeof(char), rec_length, fprec) != rec_length) {
	  sprintf(message, "end_record: Error writing %s", record_filename);
	  EdenError(message) ;
     }

     fclose(fprec) ;
     exit(flag) ;
}

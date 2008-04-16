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

                                EDEN

  Title:        Controller for Eden - the solver and its utilities.
  Author:       Hanna Szoke
  Date:         5/25/95
  Function:     This controller calls one of the following:

		addmaps		apodfc		apodfo 		back
		bin2map		cadhkl		count		
		distance	dphase		expandfc 	expandfo 	
		forth		maketar		map2bin		multmaps	
		perturbhkl 	ranphase	regrid		shapes		
		solve		sym 		tohu 		variance	

*******************************************************************************/


#include "util.h"

FILE    *fp_log ; 	

char	log_filename[MAXSTRING] ;
char	caller[MAXSTRING] ;	
char	command_line[MAXSTRING] ;

int	batch ;		/* note: NOT in util.h! it's local to eden.c */
int	silent ;
int	verbose ;
int	very_verbose ;
int	interactive ;
int	graphics ;	
int	quick ;
float	t0, t1, t2, t3, t4, t5, t6 ;	/* for reading w. sscanf() */

static	int	check_switches(int, char *[]) ;
static	void   	get_unique_logname() ;
static	void	go_do_it(int, char *[]) ;
static	int	valid_name(char *) ;

typedef struct {
  char name[32];
  int  nwarn;
  eden_method *em;
  char *warning;
} method;

method methods[] = {
  { "addmaps",	    0, addmaps_main,	NULL},	   
  { "apodfc",	    0, apodfc_main,	NULL},	   
  { "apodfo",	    0, apodfo_main,	NULL},	   
  { "back",	    0, back_main,	NULL},	   
  { "bin2map",	    0, bin2map_main,	NULL},	   
  { "bin2view",     1, NULL,		"bin2view is no longer a part of Eden"},
  { "cadhkl",	    0, cadhkl_main,	NULL},	   
  { "count",	    0, count_main,	NULL},	   
  { "distance",     0, distance_main,	NULL},   
  { "dphase",	    0, dphase_main,	NULL},	   
  { "expandfc",     0, expandfc_main,	NULL},   
  { "expandfo",     0, expandfo_main,	NULL},   
  { "forth",	    0, forth_main,	NULL},	   
  { "map2bin",	    0, map2bin_main,	NULL},	   
  { "maketar",	    0, maketar_main,	NULL},	   
  { "multmaps",     0, multmaps_main,	NULL},   
  { "perturbhkl",   0, perturbhkl_main,	NULL}, 
  { "ranphase",     0, ranphase_main,	NULL},   
  { "refine",       0, refine_main,	NULL},   
  { "regrid",	    0, regrid_main,	NULL},	   
  { "shapes",	    0, shapes_main,	NULL},	   
  { "solve",	    0, solve_main,	NULL},	   
  { "sym",	    0, sym_main,	NULL},	   
  { "tohu",	    0, tohu_main,	NULL},	   
  { "variance",     0, variance_main,	NULL},   
  { "view2bin",     1, NULL,		"view2bin is no longer a part of Eden"},
  { "drho",         1, distance_main,  "'drho' has been replaced by 'distance'" },
  { "setupncs",     2, NULL,            "'setupncs' (for NCS option) has been disabled" },
  { "tetra",        2, NULL,            "'tetra' has been disabled.  Use 'tohu' + 'back'" },
  { "keywords",     0, NULL, NULL},
  { "news",         0, NULL, NULL},
  { "switches",     0, NULL, NULL}
};
#define NMETHOD (sizeof(methods)/sizeof(methods[0]))

int	check_switches(int argc, char *argv[])
{
  /******************************************************
  Check switches: everybody can ask for help and many 
  programs have a meaningful 'verbose' at 2 levels.
  Flag 'graphics' enables call to xmgr in apodfo/c.
  ******************************************************/

     int	ch ;
     int	send_help =FALSE ;

     batch = FALSE ;
     silent = FALSE ;
     interactive = FALSE ;
     verbose = FALSE ;
     very_verbose = FALSE ;
     graphics = FALSE ;
     quick = FALSE ;

     /******************************************************
     Note: getopt() provides EDEN with variable 'optind'
     which appears throughout the code.  See man getopt.
     ******************************************************/

     while ((ch = getopt(argc, argv, "bghiqsvV")) != EOF) 
     {
       switch(ch)
       {
       case 'b':
            batch = TRUE ;
	    break ;
       case 'h':
            send_help = TRUE ;
	    break ;
       case 's':
            silent = TRUE ;
	    break ;
       case 'i':
            interactive = TRUE ;
	    break ;
       case 'v':
            verbose = TRUE ;
	    break ;
       case 'V':
            very_verbose = TRUE ;
            verbose = TRUE ;
	    break ;
       case 'g':
            graphics = TRUE ;
	    break ;
       case 'q':
            quick = TRUE ;
	    break ;
       default:
	    break ;
       }
    }
    return(send_help) ;
}

int     valid_name(char *caller)
{
  int i;
  
  for(i=0; (i<NMETHOD); i++) {
    if (strcmp(caller,methods[i].name) == 0) {
      if (methods[i].nwarn > 0)
	printf("WARNING: %s\n",methods[i].warning);
      if (methods[i].nwarn < 2)
	return TRUE;
    }
  }
  printf("\n\t%s: No such program in Eden\n", caller) ;
  ballpark("eden");
  return FALSE;
}

void     go_do_it(int argc, char *argv[])
{
  int i;
  
  for(i=0; (i<NMETHOD); i++) {
    if (strcmp(caller,methods[i].name) == 0) {
      methods[i].em(argc,argv);
      return;
    }
  }
}

void	get_unique_logname()
{   
     int	j ;
     FILE	*fpin ;
     char	*cwd ;

     /* Try the basic name - [caller].log */

     sprintf(log_filename, "%s.log", caller) ;

     if ((!batch) &&((fpin = fopen(log_filename, "r")) != NULL)) {
     
     /* Basic name is in use; find first free name of form [caller][j].log */

	 for (j = 1; j < 10; j++) {

             sprintf(log_filename, "%s%1d.log", caller, j) ;

             if ((fpin = fopen(log_filename, "r")) == NULL) 
	          break ;
             else
                  fclose(fpin) ;
         }

         /* Fall-through; all names are in use - exit with error */

         if (j == 9) {
              fprintf(stderr, "\nThere are no more free names for a log!\n") ;
	      exit(-1) ;	/* Can't call EdenError yet*/
         }

     }
     if ((fp_log = fopen(log_filename,"w")) == NULL) 
     {
          fprintf(stderr, "Cannot open %s\n", log_filename);
          exit(-1);
     }
     else
          fprintf(fp_log, "\n\t\t\t%s\n\n", log_filename) ;

     if ((cwd = getcwd(NULL, 120)) == NULL) 
	  EdenError("Can't get pwd") ;

     fprintf(fp_log, "\nRunning EDEN from %s\n", cwd) ;
     free(cwd) ;

     return;    
}

int	main(int argc, char *argv[])
{
     void	help(char *) ;
     void	end_record(int);

     int	send_help ;
     int	k ;

     send_help = check_switches(argc, argv) ;

  /******************************************************
  Check existence of environmental variable EDENHOME.
  If user types only "eden" or "eden -h", provide a brief 
  overview and quit. Otherwise, decide what's called.
  Note that exits - good or bad - cannot be put into
  the record file because it hasn't yet been established!
  ******************************************************/
 
     if (getenv("EDENHOME")  == NULL) {
	  fprintf(stderr, 
	    "Please set environment variable EDENHOME before running Eden!\n") ;
	  exit(-1) ;
     }

     if ( (argc == 1) || ((argc == 2) && (optind == 2))) 
	  ballpark("eden") ;  
     
     strcpy(caller, argv[optind]) ;

     for (k=0; k<(int) strlen(caller); k++)
          caller[k] = tolower((int) caller[k]) ;

     if (!valid_name(caller))
	  exit(-1) ;

  /******************************************************
  Deal with a request for help.  Note: users may get help 
  for "keywords", "news" and "switches" even though 
  they're not legal programs.
  ******************************************************/

     if (send_help) {
          help(caller) ;
	  exit(0) ;
     }

  /******************************************************
  Deal with e.g.  "eden regrid" (no arg list).  
  Calls to ballpark() are terminal.
  ******************************************************/

     if ((argc == 2) && (optind == 1)) 
	  ballpark(caller) ;  

  /******************************************************
  Prepare a log file, echo pwd & input line to the log.
  Identify parameter file name (which may or may not 
  have the .inp extension).
  ******************************************************/

     get_unique_logname() ;

     sprintf(command_line, argv[0]) ;

     for (k = 1; k < argc; k++) { 
         strcat(command_line, " ") ;
         strcat(command_line, argv[k]) ;
     }

     fprintf(fp_log, "Command line: %s\n", command_line) ;
     fflush(fp_log) ;

     sprintf(input_filename, "%s", argv[optind+1]) ;
     if (strstr(input_filename, ".inp") == NULL)
	  strcat(input_filename, ".inp") ;

  /******************************************************
  Finally: go do it!
  ******************************************************/

     go_do_it(argc, argv) ;

     if (!silent) 
	  printf("A log of this run has been written to %s\n", log_filename) ;
	
     printTwice("\n") ;
     printTwice(timestamp()) ;
     printTwice("\n") ;

     end_record(0) ;		/* end_record exits. */
     exit(0) ;			/* to make compiler happy */
}

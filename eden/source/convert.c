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

                           BIN2MAP, MAP2BIN 

  Title:        Bin2map, Map2bin.
  Author:       Hanna Szoke
  Date:         May 7, 2003   
  Function:     Convert .bin to .map  or .map to .bin files.
*******************************************************************************/

#include "util.h"
#include "dims.h"
#include "cellparams.h"

static  char	*mapex = ".map" ;
static	int	headread ;


void bin2map_main(int argc, char *argv[])
{
     void	readEpixFile(char *, real *) ;	/* Read input files */
     void	get_min_max(real *, int, real *, real *) ;

     FILE	*fp ;
     char	filename[MAXSTRING] ;
     char	pwd_name[MAXSTRING], out_name[MAXSTRING] ;
     real	*array ;		/* electron densities */
     void	strip_suffix(char *, char *, char *) ;
     int        n, i, j ;
     real     ar_min ;
     real     ar_max ;
     real        scale_factor = 1. ;     /* scaling factor */


  /***********************************************************
   Check for right ballpark, say hello, identify file names.
  ***********************************************************/

     if (argc > optind+3) 	/* too many arguments */
	  ballpark(caller) ; 
     else if (argc == optind+3)
          sprintf(filename, "%s", argv[optind+2]) ;
     else {
	  prompt("\nPlease enter the name of an electron/voxel file: ") ;
          sscanf(terminp, "%s", filename) ;
	  sprintf(message, "Running %s on el/voxel file %s", caller, filename) ;
	  printTwice(message) ;
     }

     hello(caller) ;

  /************************************************************
  Fetch input conditions, set up array in physical space.
  Read the bin file and write the .bin solution maps. 
  Output will be written in pwd rather than in directory from
  which the .bin files were read. 
  ************************************************************/

     readBasicInput() ;

     if (grid_type == BODY_CENTERED) {
         printTwice("You can't create a .map from a body-centered .bin file.") ;
         EdenError("Use Regrid instead!") ;
     }

     array = (real *) e_malloc(Nptotal*sizeof(real), caller) ;
          
     readEpixFile(filename, array) ;

     strip_path(filename, pwd_name) ;
     strip_suffix(pwd_name, "bin", out_name) ;
     strcat(out_name, mapex) ;

     if ((fp = fopen(out_name, "w")) == NULL)
          EdenError("Couldn't open .map file.") ;

     get_min_max(array, Np, &ar_min, &ar_max) ;

     if (ar_min == ar_max)
          EdenError("All the data are effectively equal!") ;

     fprintf(fp,"\n       3 !NTITLE\n") ;
     fprintf(fp," REMARKS Electron Density Information from EDEN\n") ;
     fprintf(fp," REMARKS DATE: %s \n", timestamp()) ;
     fprintf(fp," REMARKS Density range is (%g, %g), scale factor is %g\n",
             ar_min, ar_max, scale_factor) ;
     fprintf(fp," %7d %7d %7d %7d %7d %7d %7d %7d %7d\n",
             Nx, 0, Nx-1, Ny, 0, Ny-1, Nz, 0, Nz-1) ;
     fprintf(fp,"%12.5E%12.5E%12.5E%12.5E%12.5E%12.5E\n",
             aAxis, bAxis, cAxis, angles[0], angles[1], angles[2]) ;
     fprintf(fp, "ZYX\n") ;

     for (n = 0; n < Nz; n++) {

        fprintf(fp, "%8d\n", n) ;
        for (i = 0; i < Nx*Ny; i += 6) {
         for (j = 0; j < 6; j++)
           if (i+j < Nx*Ny)
                fprintf(fp, "%12.5E", *(array+n*Nx*Ny+i+j)) ;
         fprintf(fp, "\n") ;
        }
     }
     fprintf(fp, "   -9999\n") ;
     fclose(fp) ;
}


void map2bin_main(int argc, char *argv[])
{
     void	readXmap_header(char *, int *, int *, int *) ;  /* in maputil.c */
     void	readXmap(char *, int, int, int, real *) ;     /* in maputil.c */
     void	writeEpixFile(char *, real *) ;	        /* in qshare.c */
     void	strip_suffix(char *, char *, char *) ;

     char	map_filename[MAXSTRING] ;
     char	pwd_name[MAXSTRING], out_name[MAXSTRING] ;

     int	Nxmap, Nymap, Nzmap ;
     real	*xnump ;		/* electron densities */

  /***********************************************************
   Check for right ballpark, say hello, identify file names.
  ***********************************************************/

     if (argc > optind+3) 	/* too many arguments */
	  ballpark(caller) ; 
     else if (argc == optind+3)
          sprintf(map_filename, "%s", argv[optind+2]) ;
     else {
	  prompt("\nPlease enter the name of a .map file: ") ;
          sscanf(terminp, "%s", map_filename) ;
     }

     hello(caller) ;

  /************************************************************
  Fetch input conditions, set up array in physical space.
  Read the bin file and write the .bin solution maps. 
  ************************************************************/

     readBasicInput() ;

     sprintf(message, "Reading map file %s", map_filename) ;
     printTwice(message) ;

     readXmap_header(map_filename, &Nxmap, &Nymap, &Nzmap) ;

     Nx = Nxmap ;
     Ny = Nymap ;
     Nz = Nzmap ;
     Nptotal = Nx*Ny*Nz ;

     xnump = (real *) e_malloc(Nptotal*sizeof(real), caller) ;

     readXmap(map_filename, Nx, Ny, Nz, xnump) ;

  /************************************************************
  Output will be written in pwd rather than in directory from
  which the .bin file was read. 
  ************************************************************/

     strip_path(map_filename, pwd_name) ;
     strip_suffix(pwd_name, mapex, out_name) ;
     strcat(out_name, ".bin") ;

     sprintf(message, "Writing binary file to %s", out_name) ;
     printTwice(message) ;

     sprintf(message, "Note - this will be a simple, not body-centered file.") ;
     printTwice(message) ;

     writeEpixFile(out_name, xnump) ;
}

        /************************************************************
        Read header of file of electron densities in XPLOR format.
        ************************************************************/

void    readXmap_header(char *in_filename, int *Nx, int *Ny, int *Nz)
{
     int        first3chars() ;
     FILE       *fp ;
     int	Nxmap, Nymap, Nzmap ;
     int	rlimx[2], rlimy[2], rlimz[2] ;	/* coord limits */

     if (strstr(in_filename, mapex) == NULL)
	  strcat(in_filename, mapex) ;

     if ((fp = fopen(in_filename, "r")) == NULL)
     {
          sprintf(message,"Cannot open %s", in_filename) ;
          EdenError(message) ;
     }

     /************************************************************
     			Read header stuff  
     ************************************************************/

     headread = 0 ;

     do {
          fgets(nextline, MAXSTRING, fp) ;
          fprintf(fp_log, "%s", nextline) ;
          headread++ ; 
     } while (first3chars(nextline, "REMARKS") != 0) ;

     do {
          fgets(nextline, MAXSTRING, fp) ;
          fprintf(fp_log, "%s", nextline) ;
          headread++ ;
     } while (first3chars(nextline, "REMARKS") == 0) ;

     sscanf(nextline, "%8d %8d %8d %8d %8d %8d %8d %8d %8d", 
                      &Nxmap, &rlimx[0], &rlimx[1], 
                      &Nymap, &rlimy[0], &rlimy[1], 
                      &Nzmap, &rlimz[0], &rlimz[1]) ;
     *Nx = Nxmap ;
     *Ny = Nymap ;
     *Nz = Nzmap ;

     fgets(nextline, MAXSTRING, fp) ;   /* cell dimensions */
     headread++ ;
     fprintf(fp_log, "%s", nextline) ;
     fgets(nextline, MAXSTRING, fp) ;   /* ordering */
     headread++ ;
     fprintf(fp_log, "%s", nextline) ;
     fprintf(fp_log, "\n") ;

     fclose(fp);
     return ;
}

        /************************************************************
        Read data in file of electron densities in XPLOR format.
        ************************************************************/

void    readXmap(char *in_filename, int Nx, int Ny, int Nz, real *array)
{
     FILE       *fp ;
     int        n, i, j ;
     int	nval ;
     real     *part ;   
     float     p[6] ;    

     if ((fp = fopen(in_filename, "r")) == NULL)
     {
          sprintf(message,"Cannot open %s", in_filename) ;
          EdenError(message) ;
     }

     /************************************************************
     			Skip header stuff  
     ************************************************************/

     for  (n = 0; n < headread; n++) 
          fgets(nextline, MAXSTRING, fp) ;

     /************************************************************
     			Now read data	
     ************************************************************/

     for (n = 0; n < Nz; n++) {

        fgets(nextline, MAXSTRING, fp) ;
        part = array + n*Nx*Ny ;

        for (i = 0; i < Nx*Ny; i += 6) {
             fgets(nextline, MAXSTRING, fp) ;
             nval = sscanf(nextline, "%g %g %g %g %g %g",
                    &p[0], &p[1], &p[2], &p[3], &p[4], &p[5]) ;
             for (j = 0; j < nval; j++) 
                  *(part+i+j) = p[j] ;
        }
     }

     fclose(fp);
     return ;
}

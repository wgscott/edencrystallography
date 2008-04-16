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

				PDBUTIL.C

  Title:        PDB Utilities for Count, Regrid, Sym, Refine and Tohu.
  Author:       Hanna Szoke
  Date:         Aug. 13, 2004 (most recent changes)
  Function:     These functions handle input of pdb info and its
		propagation throughout the unit cell.

  Assumptions:  point electrons 

*******************************************************************************/

#include "util.h"
#include "cellparams.h"
#include "symmetry.h"
#include "pdb.h"

static	void	propagate_info(int, int) ;

static	float	zatoms ;	
static	float	zatomsSq ;

static	void	printPdbInfo(int, float, char[]) ;

void fetchPdbInfo(char *filename)  	/*   Read and process input .pdb data file */
{
     FILE       *fp ;
     int	first3chars() ;
     char	pdb_filename[MAXSTRING] ;
     char	gid[10], aid[10] ;
     char       tmp_char;
     char       *next_field_start;

     /**************************************************
     Examine filename; if it does NOT have extension
     .pdb, add that to the pdb_filename; otherwise, use
     pdb_filename as is but truncate filename (remove
     the 4 chars ".pdb") since the rest of the name is 
     used by Tohu, Tetra and Sym to form other names. 

     Open .pdb input file. 
     **************************************************/

     strcpy(pdb_filename, filename) ;
     if (strstr(filename, ".pdb") == NULL)
	  strcat(pdb_filename, ".pdb") ;

     check_exist(pdb_filename) ;
     fp = fopen(pdb_filename, "r") ;

     /**************************************************
     Read .pdb input. Because fields tend to coalesce, I 
     have used a fussy procedure, based on the PDB specs 
     for ATOM; my guide is pageless, so I can't refer to 
     chapter and verse.
     **************************************************/

     Nin = 0 ;

     while (fgets(nextline, MAXSTRING, fp) != NULL)  {
	  if (((int)strlen(nextline) > 1) &&
              ((strstr(nextline, "ATOM") == nextline) ||
	       (strstr(nextline, "HETATM") == nextline))) {

               if (Nin >= MAXNATOMS) {
                    sprintf(message, 
			 "Please increase MAXNATOMS (%d)", MAXNATOMS) ;
                    EdenWarning(message) ;
                    EdenError("Too many atoms!") ;
               }

     		/**************************************************
     		Original code had:
     		sscanf(nextline, "%*s %d %s %s %d %g %g %g %g %g", 
            		&atno[Nin], aid, gid, 
            		&groupno[Nin], &t1, &t2, &t3, &t4, &t5) ;

     		Filipe Maia: Extract the fixed width fields from 
	     	the pdb, according to 
     		http://www.ccp4.ac.uk/dist/html/pdbformat.html
     		(which he believes is exactly the same as the formal 
     		format description).  
     		The atof/atoi should be changed to strtof/strtol to 
     		detect input errors.
     		**************************************************/

	       tmp_char = nextline[11];		/* limit the field*/
	       nextline[11] = 0;
	       next_field_start = &(nextline[6]);
	       atno[Nin] = atoi(next_field_start); /* restore original line*/
	       nextline[11] = tmp_char;		/* advance field pointer */
	       next_field_start+= 6;
	       tmp_char = nextline[16];
	       nextline[16] = 0;

	       /* Next copy Remoteness indicator and Branch 
		  designator besides the Chemical symbol */
	       strncpy(aid,next_field_start,5);	       
	       aid[4] = 0;
	       nextline[16] = tmp_char;

	       /* Next: skip over alternate location indicator */
	       next_field_start+= 5;
	       tmp_char = nextline[20];
	       nextline[20] = 0;
	       strncpy(gid,next_field_start,4);
	       gid[4] = 0;

	       /* Next: skip over chain identifier */
	       next_field_start+= 5;
	       tmp_char = nextline[26];
	       nextline[26] = 0;
	       groupno[Nin] = atoi(next_field_start);
	       nextline[26] = tmp_char;

	       /* Start retrieving coordinates*/
	       next_field_start += 8;
	       tmp_char = nextline[38];
	       nextline[38] = 0;
	       t1 = (float)atof(next_field_start);
	       nextline[38] = tmp_char;
	       next_field_start += 8;
	       tmp_char = nextline[38];
	       nextline[46] = 0;
	       t2 = (float)atof(next_field_start);
	       nextline[46] = tmp_char;
	       next_field_start += 8;
	       tmp_char = nextline[54];
	       nextline[54] = 0;
	       t3 = (float)atof(next_field_start);

	       /* Now retrieve Occupancy and B */
	       nextline[54] = tmp_char;
	       next_field_start += 8;
	       tmp_char = nextline[60];
	       nextline[60] = 0;
	       t4 = (float)atof(next_field_start);
	       nextline[60] = tmp_char;
	       next_field_start += 6;
	       tmp_char = nextline[66];
	       nextline[66] = 0;
	       t5 = (float)atof(next_field_start);
	       nextline[66] = tmp_char;

               strncpy(atid[Nin], aid, 2) ; 
               strncpy(fullatid[Nin], aid, 5) ; 
	       fullatid[Nin][4] = 0;
	       atid[Nin][2] = 0;
	       if(isspace(atid[Nin][0]) ||
	          isdigit(atid[Nin][0])){
		 atid[Nin][0] = atid[Nin][1];		/* left justify */
		 atid[Nin][1] = ' ';
	       }
               strncpy(groupid[Nin], gid, 3) ; 
	       groupid[Nin][3] = 0;

               /* Following check is based on our sad experience */

               if (t5 < 0) {
                    sprintf(message, "ERROR -- Bj[%d] = %g < 0!\n", Nin, t5) ;
                    EdenError(message) ;
               }
               if (t4 > 1.) {
                    sprintf(message, "occupancy[%d] = %g > 1!\n", Nin, t4) ;
                    EdenWarning(message) ;
               }

               pos[Nin].x = t1 ;
               pos[Nin].y = t2 ;
               pos[Nin].z = t3 ;
               occup[Nin] = t4 ;
               Bj[Nin] = t5 ; 

               if (verbose)
                    fprintf(fp_log, 
                "# %4d, at %4s, group %s, g# %d, pos %g %g %g, occ %g, B %g\n",
                    atno[Nin], atid[Nin], groupid[Nin], groupno[Nin], 
                    pos[Nin].x, pos[Nin].y, pos[Nin].z, occup[Nin], Bj[Nin]) ;

	       Nin++ ;
          }			/* skip anything else */
     }
     fclose(fp);

     Natoms = Nin ;
}

#define NPER  106
char	legal_atom_names[NPER][2] ;

void	get_legal_atom_names()
{
     /* The full periodic table as 2 character identifiers, left-justified */
     
     strcpy(legal_atom_names [0], "H HELIBEB C N O F NE") ;
     strcpy(legal_atom_names[10], "NAMGALSIP S CLARK CA") ; 
     strcpy(legal_atom_names[20], "SCTIV CRMNFECONICUZN") ; 
     strcpy(legal_atom_names[30], "GAGEASSEBRKRRBSRY ZR") ; 
     strcpy(legal_atom_names[40], "NBMOTCRURHPDAGCDINSN") ; 
     strcpy(legal_atom_names[50], "SBTEI XECSBALACEPRND") ; 
     strcpy(legal_atom_names[60], "PMSMEUGDTBDYHOERTMYB") ; 
     strcpy(legal_atom_names[70], "LUHFTAW REOSIRPTAUHG") ; 
     strcpy(legal_atom_names[80], "TLPBBIPOATRNFRRAACTH") ; 
     strcpy(legal_atom_names[90], "PAU NPPUAMCMBKCFESFM") ; 
     strcpy(legal_atom_names[100], "MDNOLRRFHA") ;

     return;
}
void countPdbInfo()  	/*   Count atoms in a .pdb data file */
{
     void	get_legal_atom_names() ;
     int	n, z ;
     int	natoms[NPER] ;
     float	nel[NPER] ;
                                                                        
     get_legal_atom_names() ;

     for (z = 0; z < NPER; z++) {
	  natoms[z] = 0 ;
	  nel[z] = 0 ;
     }

     zatoms = 0 ;
     zatomsSq = 0 ;

     for (n = 0; n < Natoms; n++) {

	  /* Get rid of numerical junk in 2nd position - e.g. O1, N3 */

          if (!isupper(atid[n][1]))
	       atid[n][1] = ' ' ;

	  for (z = 0; z < NPER; z++) 
	       if ((atid[n][0] == legal_atom_names[z][0]) &&
	           (atid[n][1] == legal_atom_names[z][1])) {
		    fj[n] = (z+1) * occup[n] ;
		    natoms[z] ++ ;    ;
		    nel[z] += fj[n] ;
		    zatoms += fj[n] ;
		    zatomsSq += (z+1) * fj[n] ;
		    break ;
               }

          if (z == NPER) {
	      sprintf(message, "Couldn't identify '%2s'!", atid[n]) ;
	      printTwice(message) ;
	      printTwice("Probably you should reformat your pdb file.") ;
	      EdenError("Please refer to 'awk_pdb' in $EDENHOME/tools") ;
          }
     }

     printTwice("\nPdb file contains:") ;

     for (z = 0; z < NPER; z++)
          printPdbInfo(natoms[z], nel[z], legal_atom_names[z]) ;
     
     printTwice("  ") ;

     if (Natoms > 1) {
	  sprintf(message,
          "Total: %6d atoms. %8g electrons", Natoms, zatoms) ;
          printTwice(message) ;
	  sprintf(message,
          "Total zsq : %g \n", zatomsSq) ;
          printTwice(message) ;
     }
}


/* Return the atomic number of a 2 char symbol or 0 on error */

int getZfromSymbol(char * symbol)
{
  void	get_legal_atom_names() ;
  int z;

  get_legal_atom_names() ;

  if(isalpha(symbol[0])){
    symbol[0] = toupper(symbol[0]);
  }else if(isalpha(symbol[1])){
    symbol[0] = toupper(symbol[1]);
    symbol[1] = ' ';
  }else{
    return 0;
  }
  if(isalpha(symbol[1]))
    symbol[1] = toupper(symbol[1]);
  else
    symbol[1] = ' ';
  symbol[2] = 0;
  for (z = 0; z < NPER; z++) 
    if ((symbol[0] == legal_atom_names[z][0]) &&
	(symbol[1] == legal_atom_names[z][1])) 
      return z+1;     
  return 0;
}

void	repositionPdb()
{
     int	n;
     real	sumx = 0, sumy = 0, sumz = 0, cases = 0 ;
     void	writePdbInfo(char *, char *, int) ;

     for (n = 0; n < Nin; n++) {

          sumx += pos[n].x * occup[n] ;
          sumy += pos[n].y * occup[n] ;
          sumz += pos[n].z * occup[n] ;
	  cases += occup[n] ;

     }
     sumx /= cases ;
     sumy /= cases ;
     sumz /= cases ;

     for (n = 0; n < Nin; n++) {

          pos[n].x -= sumx ;
          pos[n].y -= sumy ;
          pos[n].z -= sumz ;

     }
     printTwice("Writing repositioned Pdb file to 'new.pdb'") ;
     writePdbInfo("new.pdb", " ", 0) ;
}

void	randomize(int *seed)
{
     float      ran1(int *) ;                /* packaged in util.c */
     void	writePdbInfo(char *, char *, int) ;

     int	n ;
     real	x, y, z ;
     real	halfa = aAxis / 2. ;
     real	halfb = bAxis / 2. ;
     real	halfc = cAxis / 2. ;


     for (n = 0; n < Nin; n++) {
	  x = ran1(seed) ;
	  y = ran1(seed) ;
	  z = ran1(seed) ;
          pos[n].x = aAxis*x - halfa ;
          pos[n].y = bAxis*y - halfb ;
          pos[n].z = cAxis*z - halfc ;
     }
     printTwice("Writing random-positioned Pdb file to 'ran.pdb'") ;
     writePdbInfo("ran.pdb", " ", 0) ;
}

void	printPdbInfo(int totatoms, float totel, char name[2])
{
     if (totatoms > 0) {
         sprintf(message,
     "    %6d %c%c atoms, %8g electrons   ", 
	       totatoms, name[0], name[1], totel) ;
         printTwice(message) ;
     }
}


void expandPdbInfo()  	/*   expand input PDB data to unit cell, using 
			crystal symmetry.  */
{
	 void	apply_symop_nowrap(real, real, real, real *, real *,
				   real *, real *) ;
	 int	j, jnew, i ;
	 real	*mat ;
	 real	fx, fy, fz ;
	 real	newfx, newfy, newfz ;
     
        if ((Natoms = Nin*Ncs) >= MAXNATOMS) {
             sprintf(message, "Please increase MAXNATOMS (%d)", MAXNATOMS) ;
	     EdenWarning(message) ;
             EdenError("Too many atoms!") ;
        }

	for (j = 0; j < Nin; j++) {

              ortho_to_fract(pos[j].x, pos[j].y, pos[j].z, &fx, &fy, &fz) ;

              for (i = 1, mat = matop+MEL; i < Ncs; i++, mat += MEL) {

                   jnew = i*Nin + j ;

                   apply_symop_nowrap(fx, fy, fz, mat, &newfx, &newfy, &newfz) ;
                   fract_to_ortho(newfx, newfy, newfz, 
				  &pos[jnew].x, &pos[jnew].y, &pos[jnew].z) ;

                   propagate_info(j, jnew) ;
 	     }
 	}
}
void expandPartialPdbInfo(int N1, int N2, int skip)  	

	/*   expand partial PDB data to unit cell, using crystal symmetry. 
	N1, N2 : index range of input 
	skip :	 How far (in index space) to skip between asu's   */ 
{
	 void	apply_symop_nowrap(real, real, real, real *, real *,
				   real *, real *) ;
	 int	j, jnew, i ;
	 real	*mat ;
	 real	fx, fy, fz ;
	 real	newfx, newfy, newfz ;
     
        if ((N2 - N1)*Ncs >= MAXNATOMS) {
             sprintf(message, "Please increase MAXNATOMS (%d)", MAXNATOMS) ;
	     EdenWarning(message) ;
             EdenError("Too many atoms!") ;
        }

	for (j = N1; j < N2; j++) {

              ortho_to_fract(pos[j].x, pos[j].y, pos[j].z, &fx, &fy, &fz) ;

              for (i = 1, mat = matop+MEL; i < Ncs; i++, mat += MEL) {

                   jnew = i*skip + j ;

                   apply_symop_nowrap(fx, fy, fz, mat, &newfx, &newfy, &newfz) ;
                   fract_to_ortho(newfx, newfy, newfz, 
				  &pos[jnew].x, &pos[jnew].y, &pos[jnew].z) ;

                   propagate_info(j, jnew) ;
              }	     
 	}

}
void	propagate_info(int j, int jnew) 
{
         strcpy(atid[jnew], atid[j]) ;
         strcpy(groupid[jnew], groupid[j]) ;

         atno[jnew] = atno[j] ;
         groupno[jnew] = groupno[j] ;
         occup[jnew] = occup[j] ;
         Bj[jnew] = Bj[j] ;
         fj[jnew] = fj[j] ;
}

#include <sys/time.h>
void	writePdbInfo(char *filename, char *label, int N)
{
       FILE	*fpout, *fpin ;
       int	i ;
       char	myself[MAXSTRING] ;

       if ((fpout = fopen(filename, "w")) == NULL)
       {
          sprintf(message, "Cannot open %s", filename);
          EdenError(message) ;
       }

     /* find out who is the user    */

     system("whoami >tmp") ;
     if ((fpin = fopen("tmp", "r")) != NULL) {
          fscanf(fpin, "%s", myself) ;
          fclose (fpin) ;
          system("rm tmp") ;
     }
     else
          strcpy(myself, "unknown ") ;

       fprintf(fpout, "REMARK FILENAME=%s, %s \n", filename, label) ;
       fprintf(fpout, "REMARK DATE: %s created by user: %s \n", 
               timestamp(), myself) ;

       /*  Note %-4s below for left justification - we try to copy Fortran,
	but it still is not quite the same as an XPLOR-generated file! */

       for (i = N; i < N + Natoms; i++) {
           if (marker[i] == 0)
              fprintf(fpout, 
/*** Fixes for nasty pdb files 3/30/98 and 5/15/99
                 "ATOM %6d %-4s %3s %5d     %7.3f %7.3f %7.3f %5.2f %5.2f\n",
                 "ATOM %6d  %-3s %3s %5d     %7.3f %7.3f %7.3f %5.2f %5.2f\n",
*****/
                 "ATOM %5d  %-3s %3s %5d     %7.3f %7.3f %7.3f %5.2f %5.2f\n",
                  atno[i], atid[i], groupid[i], groupno[i], pos[i].x, 
	          pos[i].y, pos[i].z, occup[i], Bj[i]) ;
           else
              fprintf(fpout, 
                 "%%ATOM %6d %-4s %3s %5d     %7.3f %7.3f %7.3f %5.2f %5.2f\n",
                 atno[i], atid[i], groupid[i], groupno[i], pos[i].x, 
		 pos[i].y, pos[i].z, occup[i], Bj[i]) ;
       }
       fprintf(fpout, "END\n") ;
}
void	estimate_F000()
{
     int	Fest ;
     real	Matthews ;
     real	Nel_prot  ;
     real	protein_frac ; 

     /************************************************
     Generate and report an estimate of F000, using

     F000 = [1 - rho_solvent/rho_protein]*Ncs*zatoms + 
	    rho_solvent*Vunit_cell.

     For rho_solvent = 0.34, rho_protein = 0.42
     in electrons/cuA.    Also, 

     Matthews' coefficient = Vunit_cell / (Ncs*zatoms)
     and  protein_fraction = Vprot / Vunit_cell.
     *************************************************/

     Nel_prot = Ncs*zatoms ;
     Fest = (int) ((1. - SOLVENT_DENSITY/PROTEIN_DENSITY)*Nel_prot +
		   SOLVENT_DENSITY * crys_vol) ;
     if (Fest > 10000)
	  Fest = 1000 * (int) ((Fest + 500.)/1000.) ;
     else if (Fest > 100)
	  Fest = 10 * (int) ((Fest + 5.)/10.) ;

     Matthews = crys_vol / Nel_prot ;
     protein_frac = Nel_prot / (crys_vol * PROTEIN_DENSITY) ;

     sprintf(message, 
     "total pdb electrons in unit cell %g corresponds to Fobs(000) ~ %d", 
	    Nel_prot, Fest);
     printTwice(message) ;

     sprintf(message, "Matthews' coefficient is %g", Matthews) ;
     printTwice(message) ;

     sprintf(message, "and protein fraction is %g.", protein_frac) ;
     printTwice(message) ;
}

void	getRange(int N, real *fminx, real *fmaxx, real *fminy,
		 real *fmaxy, real *fminz, real *fmaxz) 
{
     int	i ;
     real	onefx, onefy, onefz ;		/* single fract. coord. */

     *fminx = +10., *fminy = +10., *fminz = +10. ;
     *fmaxx = -10., *fmaxy = -10., *fmaxz = -10. ;

     /* Determine fractional limits of pdb file */

     for (i = 0; i < N; i++) {
            ortho_to_fract(pos[i].x, pos[i].y, pos[i].z, 
                                &onefx, &onefy, &onefz) ;

            if (onefx < *fminx) *fminx = onefx ;
            if (onefx > *fmaxx) *fmaxx = onefx ;
            if (onefy < *fminy) *fminy = onefy ;
            if (onefy > *fmaxy) *fmaxy = onefy ;
            if (onefz < *fminz) *fminz = onefz ;
            if (onefz > *fmaxz) *fmaxz = onefz ;
     }
}

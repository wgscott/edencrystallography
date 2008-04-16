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

                                BASICS.C

  Title:        General initialization package for EDEN.
  Author:       Hanna Szoke
  Date:         1/23/92
  Function:     This package contains set-up functions for all EDEN programs,
		involving the reading and interpretation of simple .inp file,
		and the echoing of its contents to the log.

******************************************************************************/

#include  "util.h"
#include  "cellparams.h"
#include  "dims.h"
#include  "symmetry.h"
#include  "pdb.h"

#define   MAXLINES 200
#define	  RIGHT	90.
#define	  DEVANG  15.	/* max deviation from 90 deg for valid bcc grid_type*/

  /************************************************************
  declarations of EXTERN variables in cellparams.h        
  ************************************************************/

real  aAxis, bAxis, cAxis ;	/* dimensions (Angstrom) of unit cell	*/
real	angles[3] ;		/* 3 angles (degrees) of axes of unit cell */
real	crys_vol ; 		/* volfac() * aAxis * bAxis * cAxis	*/
float	cubA_to_vox ;		/* conversion factor, el/cubic A to el/voxel */
int	crystal_symmetry ;

  /************************************************************
  declarations of most EXTERN variables in dims.h               
  ************************************************************/

int     Nx, Ny, Nz ;            /* solution space dimensions */
float	floatNx ;
float	floatNy ;
float	floatNz ;
int     Np ;                    /* product of Nx*Ny*Nz */
int     Nptotal ;               /* Np (sc) or 2*Np (bcc) */
int	Npextended ;		/* Nptotal + Nhr (high resolution)	*/
int	grid_type;		/* SIMPLE or BODY_CENTERED */
int	Nh, Nk, Nl ;		/* reciprocal space dimensions		*/
int     Nhkl ;                  /* number of reciprocal lattice points */
real	grid_spacing ;          /* gridding resolution in Angstrom */
real	dr ;			/* actual average resolution in Angstrom */
real	eta ;                   /* factor used in setting up exp factors */
real	delfac ;		/* PISQ * eta * dr * dr	*/
real	input_res;		/* resolution input variable	*/
real	low_res_cutoff = 1000;	/* low resolution cutoff */
real	limit_rat ;    		/* 1/(dr*dr)	*/

int     useSig ;        	/* if TRUE, use sigmas in input fobs files */
int	anom_flag = FALSE ;	/* if TRUE, input fobs/fcalc are anomalous */
float   solvent_voxel ;		/* ave. voxel value for the solvent */
float   protein_voxel ; 	/* ave. voxel value for the protein */
int	HighRes ;		/* if TRUE, high-resolution processing on  */
float   min_voxel ; 		/* min. voxel value for a crystal (for ccg) */
float   max_voxel ; 		/* max. voxel value for a crystal (for ccg) */

 /************************************************************
  declarations of EXTERN variables in util.h               
  ************************************************************/

char	input_filename[MAXSTRING] ;
char	*allinp ;
char	nextline[MAXSTRING] ;
char	goodline[MAXLINES] ;
char	id[MAXSTRING] ;
char	mode_info[MAXSTRING] ;  /* information about run mode */
int	Nlines ;

int	fetchTFStatus(char *, int) ;

static	real	Na_error, Nb_error, Nc_error ;
static	float	sdens, pdens, mindens, maxdens ; 
static	char	record_info[MAXSTRING] ;  

static	void	dump_to_log() ;
static	int	fetchAllInput() ;
static	void	readCellInput() ;
static	void	get_res_params() ;
static	int	interpret_sym_group() ;
static	void	intercheck() ;
static	int	fetchNpdims() ;
static  void	fetchDensities() ;
static	void	basic_echo() ;
static	int	get_grid_type() ;
static	void	fetch_mode_info() ;
static	void	fetch_record_info() ;

     /************************************************************
     Read basic input: this is called by ALL Eden programs     
     ************************************************************/

void readBasicInput()    
{
     void	prepare_F_and_Finv() ;
     void	setupNusq() ;
     real	volfac() ;
     void	start_record() ;

     /************************************************************
     Dump the .inp file as-is to the .log file for reference.
     Fetch all the input, place it in a buffer, allinp.
     ************************************************************/

     check_exist(input_filename) ;
     sprintf(message, "\tReading file %s\n", input_filename) ;
     printTwice(message) ;

     dump_to_log() ;
     Nlines = fetchAllInput() ;

     /************************************************************
     Read input, identifying data in terms of keywords;
     a, b, c, angles and input_res are always read.
     Optionally, grid_type and eta are also read.
     ************************************************************/

     readCellInput() ;

     grid_type = get_grid_type() ;
     get_res_params() ;

     /************************************************************
     Identify the symmetry group. Check for consistency.        
     ************************************************************/

     crystal_symmetry = interpret_sym_group() ;
     intercheck() ;
		
     /************************************************************
     Define dimensions.   
     ************************************************************/

     Np = fetchNpdims() ;

     floatNx = (float) Nx ;
     floatNy = (float) Ny ;
     floatNz = (float) Nz ;

     Nptotal = (grid_type == SIMPLE) ? Np : 2*Np ;
     Npextended = Nptotal ;
     Nh = Nx ;
     Nk = 2 * Ny ;
     Nl = 2 * Nz ;
     Nhkl = Nh * Nk * Nl ;

     dr = sqrt(((aAxis/Nx)*(aAxis/Nx) + (bAxis/Ny)*(bAxis/Ny) + 
		(cAxis/Nz)*(cAxis/Nz))/3) ;

     Na_error = (Nx - (aAxis/dr)) / (aAxis/dr) * 100. ;
     Nb_error = (Ny - (bAxis/dr)) / (bAxis/dr) * 100. ;
     Nc_error = (Nz - (cAxis/dr)) / (cAxis/dr) * 100. ;

     limit_rat = 1 / (dr*dr) ;
     delfac = PISQ * eta * dr * dr ;

     /************************************************************
     Use sigmas? Average solvent and protein voxel values?
     Minimum. maximum densities for ccg?
     Note: we call prepare_F_and_Finv "just to be on the safe 
     side" although it's needed only in connection with pdb files.
     ************************************************************/

     setupNusq() ;
     useSig = fetchTFStatus("USESIG", TRUE) ;

     crys_vol = volfac() * aAxis * bAxis * cAxis ;
     cubA_to_vox = crys_vol / Nptotal ;

     fetchDensities() ;
     solvent_voxel = sdens * cubA_to_vox ;
     protein_voxel = pdens * cubA_to_vox ;
     min_voxel = mindens * cubA_to_vox ;
     max_voxel = maxdens * cubA_to_vox ;

     HighRes = fetchTFStatus("HIGHRES", FALSE) ;
     fetch_mode_info() ;
     prepare_F_and_Finv() ;

     anom_flag = fetchTFStatus("ANOM", FALSE) ;
     
     fetch_record_info() ;
     start_record();

     basic_echo() ;
}

     /************************************************************
     Dump all the .inp file as-is to the .log file.    
     ************************************************************/

void	dump_to_log()
{
     FILE       *fp ;

     fp = fopen(input_filename, "r") ;
     
     fprintf(fp_log,"Dump of %s follows:\n", input_filename) ;

     while (fgets(nextline, MAXSTRING, fp) != NULL)  {
          fprintf(fp_log, "\t%s", nextline) ;
     }

     fprintf(fp_log,"\nEnd of %s dump.\n", input_filename) ;
     fclose(fp);
}

int	fetchAllInput()
{
     FILE       *fp ;
     int	N ;
     int	k, n ;
     int	nid ;
     int	filelength() ;


     check_exist(input_filename) ;

               /***********************************************
               Count lines of information in input file.
               ***********************************************/

     N = filelength(input_filename) ;

     if (N > MAXLINES)
         EdenError("Not enough space for input!  Please check your input file!") ;

     allinp   = (char *) e_malloc(sizeof(char)*MAXSTRING*N, 
		"fetchAllInput") ;

     for (n = 0; n < MAXSTRING*N; n++)
	  *(allinp+n) = '\b';

     for (n = 0; n < MAXLINES; n++)
	  goodline[n] = FALSE ;

     fp = fopen(input_filename, "r") ;

     N = 0 ;

     while (fgets(nextline, MAXSTRING, fp) != NULL)  {
          id[0] = '\0' ;
          if (net_strlen(nextline) > 1) {
               sscanf(nextline, "%s", id) ;
	       nid = (int)strlen(id) ;

               /***********************************************
               Identify data in terms of upper-case keywords.
               ***********************************************/

	       for (k=0; k<nid; k++)
		    nextline[k] = toupper((int) id[k]) ;

               sprintf((allinp+N*MAXSTRING), "%s", nextline) ;
	       N++ ;
          }
     }
     fclose(fp);

     return (N) ;
}

void readCellInput()    
{
     int	k ;

     /**************
     Defaults 
     **************/

     t1 = t2 = t3 = t4 = t5 = t6 = 0 ;
     aAxis = bAxis = cAxis = 0 ;
     angles[0] = angles[1] = angles[2] = 0. ;

     /* Note use of &t1 to capture single precision input! */

     for (k = 0; k < Nlines; k++) {
          strcpy(nextline, allinp + k*MAXSTRING) ;
          sscanf(nextline, "%s", id) ;

          if (strcmp(id, "CELL") == 0) {
               sscanf(nextline, "%*s %g %g %g %g %g %g", 
			   &t1, &t2, &t3, &t4, &t5, &t6) ;
               aAxis = t1 ;
               bAxis = t2 ;
               cAxis = t3 ;
               angles[0] = t4 ;
               angles[1] = t5 ;
               angles[2] = t6 ;
	       goodline[k] = TRUE ;
	       break ;
          }
     }

     /*********************************************
     Check for potentially missing information. 
     Check for fatally flawed input. 
     By convention, all angles should be either 
     <= 90 degrees or >= 90 degrees.
     *********************************************/

     if ((aAxis == 0) || (bAxis == 0) || (cAxis == 0)) 
	 EdenError("Missing dimension of unit cell!") ; 
     
     if ((angles[0] == 0) || (angles[1] == 0) || (angles[2] == 0)) 
	 EdenError("Missing angles of unit cell!") ; 


     if (!(((angles[0] <= RIGHT) && (angles[1] <= RIGHT) && (angles[2] <= RIGHT)) ||
         ((angles[0] >= RIGHT) && (angles[1] >= RIGHT) && (angles[2] >= RIGHT)))) 

	  EdenError(
     "By convention, crystallographic angles should all be <= 90 or >= 90.") ;
		 
}     

int	interpret_sym_group()
{
     int	k ;
     char       symmetry_group[MAXSTRING] ;
     int	fetchSymop(char *) ;		/* packaged in crysutil.c */

     strcpy(symmetry_group, "none") ;		

     for (k = 0; k < Nlines; k++) {
          strcpy(nextline, allinp + k*MAXSTRING) ;
          sscanf(nextline, "%s", id) ;

          if (strcmp(id, "SYMMETRY") == 0) {
               sscanf(nextline, "%*s %s", symmetry_group) ;
	       goodline[k] = TRUE ;
               break ;
          }
     }

     if (strcmp(symmetry_group, "none") == 0) 
	 EdenError("Missing symmetry group information!") ;

     for (k=0; k<(int)strlen(symmetry_group); k++)
         symmetry_group[k] = toupper((int) symmetry_group[k]) ;

     /*********************************************
     Compare input against the 230 symmetry groups, 
     retrieve the symmetry system and rotation-
     translation matrices.
     *********************************************/

     return (fetchSymop(symmetry_group)) ;
}

     /*****************************************************************
     Prior to 6/23/98, "lattice" was required input;  since then, it is 
     allowed, but the default is set on the basis of the cell angles.
     As of 2/22/99, name changed to "grid_type".
     *****************************************************************/

int get_grid_type()  
{
     char       input_grid_type[MAXSTRING] ;
     int	k ;
     int	default_value ;
     int	value ;

     strcpy(input_grid_type, "none") ;

     for (k = 0; k < Nlines; k++) {
          strcpy(nextline, allinp + k*MAXSTRING) ;
          sscanf(nextline, "%s", id) ;

          if (strcmp(id, "GRID_TYPE") == 0) {
               sscanf(nextline, "%*s %s", input_grid_type) ;
	       goodline[k] = TRUE ;
               break ;
          }
     }

     /****************************************************
     Identify the grid_type.  By default - use angle info.
     ****************************************************/

     if ((fabs(angles[0] - RIGHT) <= DEVANG) &&
         (fabs(angles[1] - RIGHT) <= DEVANG) &&
         (fabs(angles[2] - RIGHT) <= DEVANG))
 
         default_value = value = BODY_CENTERED ;
     else 
         default_value = value = SIMPLE ;
      

     if (strcmp(input_grid_type, "none") != 0) {

         for (k=0; k<(int)strlen(input_grid_type); k++)
             input_grid_type[k] = toupper((int) input_grid_type[k]) ;

         if (strcmp(input_grid_type, "SIMPLE") == 0)

    	    value = SIMPLE ;

         else if ((strcmp(input_grid_type, "BODY_CENTERED") == 0) ||
                  (strcmp(input_grid_type, "BODY-CENTERED") == 0))

	       value = BODY_CENTERED ;

         else {
	     sprintf(message, 
		  "Cannot interpret grid type %s", input_grid_type) ;
	     EdenError(message) ;
         }

     /**************************
     Somewhat superfluous check!
     **************************/

	 if ((default_value == SIMPLE) && (value == BODY_CENTERED)) {
             for (k = 0; k < 3; k++) {		 
	         if (fabs(angles[k] - RIGHT) > DEVANG) {
                     sprintf(message, 
		 "WARNING: angle (%g) deviates significantly from 90 degrees",
		 angles[k]) ;
	             EdenWarning(message) ;
	             EdenWarning(   
	         "It's preferable to use a simple grid for this crystal!") ;
                 }
             } 
         }
     }

     return (value) ;
}

void intercheck() 
{

     /************************************************************
     The following checks are taken from International Tables for
     Crystallography, Volume A (1983), p. 13.
     There are no restrictions on triclinic crystals.

     Note that the b-unique setting is assumed throughout.  Note 
     too the assumption of hexagonal axes for all trigonal crystals.
     ************************************************************/

     if (strcmp(crystal_system, "MONOCLINIC") == 0) {
	if ((angles[0] != RIGHT) || (angles[2] != RIGHT)) 
	   EdenError("EDEN requires monoclinic crystals with b axis unique.") ;
     }		 
     else if (strcmp(crystal_system, "ORTHORHOMBIC") == 0) {
	if ((angles[0] != RIGHT) || (angles[1] != RIGHT) || (angles[2] != RIGHT)) 
	   EdenError(
	     "All angles must be 90 degrees for orthorhombic crystals.") ;
     }		 
     else if (strcmp(crystal_system, "TETRAGONAL") == 0) {
	if ((angles[0] != RIGHT) || (angles[1] != RIGHT) || (angles[2] != RIGHT))
	   EdenError(
             "All angles must be 90 degrees for tetragonal crystals.") ;
	if (aAxis != bAxis) 
	   EdenError(
             "Cell dimensions a and b must be equal for tetragonal crystals.") ;
     }		 
     else if (strcmp(crystal_system, "TRIGONAL") == 0) {
	if ((angles[0] != RIGHT) || (angles[1] != RIGHT) || (angles[2] != 120.)) {
	   EdenWarning(    
      "Angles must be 90. 90 and 120 degrees for trigonal crystals.") ;
	   EdenError("Eden does not handle rhombohedral axes.") ;
         }
	if (aAxis != bAxis) 
	   EdenError(
              "Cell dimensions a and b must be equal for trigonal crystals.") ;
     }		 
     else if (strcmp(crystal_system, "HEXAGONAL") == 0) {
	if ((angles[0] != RIGHT) || (angles[1] != RIGHT) || (angles[2] != 120.)) 
	   EdenError(
            "Angles must be 90, 90 and 120 degrees for hexagonal crystals.") ;
         
	if (aAxis != bAxis) 
	   EdenError(
             "Cell dimensions a and b must be equal for hexagonal crystals.") ;
     }		 
     else if (strcmp(crystal_system, "CUBIC") == 0) {
	if ((angles[0] != RIGHT) || (angles[1] != RIGHT) || (angles[2] != RIGHT)) 
	   EdenError("All angles must be 90 degrees for cubic crystals.") ;
         
	if ((aAxis != bAxis) || (aAxis != cAxis)) 
	   EdenError("All cell dimensions must be equal for cubic crystals.") ;
     }		 
     else if (strcmp(crystal_system, "TRICLINIC") != 0) {
	sprintf(message, "Unknown crystal system - %s", crystal_system) ;
        EdenError(message) ;
     }		 
}

#define	SC_RESFAC	0.6 	/* see paper V	*/
#define	BCC_RESFAC	0.7 	/* (2/sqrt(3)) * SC_RESFAC	*/

     /************************************************************
     The following were CHANGED from long-standing old values on 
     6/23/00.  User input of grid_spacing was also dropped.   
     NOTE: Default res factors and etas lead to a pdb B value of
     about 11-1/2.
     ************************************************************/

#define	NEW_SC_ETA	0.8 	/* used to be 1.0	*/
#define	NEW_BCC_ETA	0.6 	/* used to be 0.75	*/

void	get_res_params()
{
     int	k ;

     /* initializations */

     input_res = UNSET ;
     eta = UNSET ;

     for (k = 0; k < Nlines; k++) {
          strcpy(nextline, allinp + k*MAXSTRING) ;
          sscanf(nextline, "%s", id) ;

          if (strcmp(id, "ETA") == 0) {
               sscanf(nextline, "%*s %g", &t1) ;
               eta = t1 ;
	       goodline[k] = TRUE ;
          }
          else if ((strcmp(id, "INPUT_RES") == 0) || 
                   (strcmp(id, "RESOLUTION") == 0)) {
               sscanf(nextline, "%*s %g", &t1) ;
               input_res = t1 ;
	       goodline[k] = TRUE ;
          }
          else if ((strcmp(id, "LOWRES") == 0)){
	    sscanf(nextline, "%*s %g", &t1) ;
               low_res_cutoff = t1 ;
	       goodline[k] = TRUE ;
          }
     }

     /* put it all together ... 
     no more old-style input of grid_spacing;
     new style is obligatory input_res only  */

     if (input_res == UNSET) 
	  EdenError("Missing input resolution!") ;
     
     if (grid_type == SIMPLE) { 
	  if (eta == UNSET) 
               eta = NEW_SC_ETA ;

	  grid_spacing = input_res * SC_RESFAC ;
     }
     else {
	  if (eta == UNSET) 
               eta = NEW_BCC_ETA ;

	  grid_spacing = input_res * BCC_RESFAC ;
     }

     if ((eta < 0.5) || (eta > 0.9)) {
	  sprintf(message, "Eta (%g) is out of range (0.5 - 0.9).", eta) ;
	  EdenError(message) ;
     }
}

void	fetchDensities()
{
     int	k ;

     sdens = SOLVENT_DENSITY ;  /* default: 0.34 - disordered water */
     pdens = PROTEIN_DENSITY ;  /* default: 0.42 - average protein */
     mindens = MIN_DENSITY ;	 /* default: 0 */
     maxdens = MAX_DENSITY ;	 /* default: 1.e10 */

     for (k = 0; k < Nlines; k++) {
          strcpy(nextline, allinp + k*MAXSTRING) ;
          sscanf(nextline, "%s", id) ;

          if (strcmp(id, "SOLVENT") == 0) {
               sscanf(nextline, "%*s %g", &sdens) ;
               goodline[k] = TRUE ;
          }
          else if (strcmp(id, "PROTEIN") == 0) {
               sscanf(nextline, "%*s %g", &pdens) ;
               goodline[k] = TRUE ;
          }
          else if (strcmp(id, "MIN_DENS") == 0) {
               sscanf(nextline, "%*s %g", &mindens) ;
               goodline[k] = TRUE ;
          }
          else if (strcmp(id, "MAX_DENS") == 0) {
               sscanf(nextline, "%*s %g", &maxdens) ;
               goodline[k] = TRUE ;
          }
     }
}

void fetch_coefs(real *C1, real *C2)    /* Read C1 and C2 */
{
     int	k ;

     *C1 = *C2 = 1. ;

     for (k = 0; k < Nlines; k++) {
          strcpy(nextline, allinp + k*MAXSTRING) ;
          sscanf(nextline, "%s", id) ;

          if (strcmp(id, "C1") == 0) {
               sscanf(nextline, "%*s %g", &t1) ;
               *C1 = t1 ;
          }
          else if (strcmp(id, "C2") == 0) {
               sscanf(nextline, "%*s %g", &t1) ;
               *C2 = t1 ;
          }
     }

     sprintf(message, "Using coefficients (%g, %g).", *C1, *C2) ;
     printTwice(message) ;
}

void	fetch_mode_info()
{
     int	k ;

     strcpy(mode_info, "none") ;

     for (k = 0; k < Nlines; k++) {
          strcpy(nextline, allinp + k*MAXSTRING) ;
          sscanf(nextline, "%s", id) ;

          if (strcmp(id, "MODE") == 0) {
               sscanf(nextline, "%*s %s", mode_info) ;
	       goodline[k] = TRUE ;
          }
     }

     for (k=0; k<(int)strlen(mode_info); k++)
          mode_info[k] = tolower((int) mode_info[k]) ;

}
void	fetch_record_info()
{
     void	setup_records(char *);
     int	k ;

     strcpy(record_info, "history") ;

     for (k = 0; k < Nlines; k++) {
          strcpy(nextline, allinp + k*MAXSTRING) ;
          sscanf(nextline, "%s", id) ;

          if (strcmp(id, "RECORD") == 0) {
               sscanf(nextline, "%*s %s", record_info) ;
	       goodline[k] = TRUE ;
          }
     }

     setup_records(record_info) ;
}

int	fetchNpdims() 

     /**********************************************************************
     fetchNpdims is called by all Eden programs that deal with gridded 
     arrays of electrons/voxel and Fourier transforms.  It bumps dx, dy and 
     dz to the closest product of powers allowed by the FFTs.  For certain 
     symmetries, we must impose further conditions on the dimensions.  
     These are expressed in maxdenom[i], i = 0, 2 and reflect the need to 
     have integer indices when applying a translational symmetry operation 
     on a given (x, y, z).   For most space groups, maxdenom[i] = 1 or 2, 
     but for P41, e.g., maxdenom[2] = 4 and other possible values are 3 & 6.

     The actual effective resolution is determined by averaging 
     aAxis/Nx, bAxis/Ny and cAxis/Nz.         
     **********************************************************************/

{
     int	getFftDim(float, int) ;
     float	dx, dy, dz ;

     dx = aAxis / grid_spacing ;
     dy = bAxis / grid_spacing ;
     dz = cAxis / grid_spacing ;

     Nx = getFftDim(dx, maxdenom[0]) ;
     Ny = getFftDim(dy, maxdenom[1]) ;
     Nz = getFftDim(dz, maxdenom[2]) ;

     return (Nx * Ny * Nz) ;
}

int fetchTFStatus(char *key, int dfault)  
{
     char	setting[MAXSTRING] ;
     int	value ;
     int	k ;

     strcpy(setting, "") ;

     for (k = 0; k < Nlines; k++) {
          strcpy(nextline, allinp + k*MAXSTRING) ;
          sscanf(nextline, "%s", id) ;

          if (strcmp(id, key) == 0) {
               sscanf(nextline, "%*s %s", setting) ;
	       goodline[k] = TRUE ;
               break ;
          }
     }

     value = dfault ;

     if ((int)strlen(setting) > 0) {

          for (k=0; k<(int)strlen(setting); k++)
               setting[k] = toupper((int) setting[k]) ;

          if ((strcmp(setting, "ON") == 0) ||
              (strcmp(setting, "TRUE") == 0)) 
	       value = TRUE ;
          else if ((strcmp(setting, "OFF") == 0) ||
              (strcmp(setting, "FALSE") == 0)) 
	       value = FALSE ;
	  else
	       goodline[k] = FALSE ;
     }
     return (value) ;
}


void	basic_echo() 
{
     fprintf(fp_log, "\n") ;
     fprintf(fp_log, "Unit cell measures %6.2f by %6.2f by %6.2f Angstrom\n",
             aAxis, bAxis, cAxis) ;

     fprintf(fp_log, "Alpha = %6.2f, beta = %6.2f, gamma = %6.2f degrees\n",
             angles[0], angles[1], angles[2]) ;

     fprintf(fp_log, "Scale factor for converting el/A^3 to el/grid pt is %g\n",
		     crys_vol / Nptotal) ;
     fprintf(fp_log, "Symmetry is %s.\n", SGname) ;
     fprintf(fp_log, "Input resolution is %g Angstrom.\n", input_res) ;

     /************************************************************
     Display and solution dimensions.
     ************************************************************/

     fprintf(fp_log, "\n") ;
     fprintf(fp_log, "Gridding resolution is %g Angstrom, eta is %g\n", 
		      grid_spacing, eta) ;
     fprintf(fp_log, "Input unit cell partitions are %4.1f by %4.1f by %4.1f\n",
		      aAxis/grid_spacing, bAxis/grid_spacing, 
		      cAxis/grid_spacing) ;
     fprintf(fp_log, "Actual unit cell partitions are %d by %d by %d\n", 
		Nx, Ny, Nz) ;
     fprintf(fp_log, "Average resolution in Angstrom is dr = %g\n", dr) ;
     fprintf(fp_log, 
     	"Partition errors (%%) with respect to the average resolution\n") ;
     fprintf(fp_log, "  in a, b and c are %g, %g, and %g.\n", 
		Na_error, Nb_error, Nc_error) ;
     fprintf(fp_log, "\n") ;
     if (grid_type == SIMPLE)
    	  fprintf(fp_log, "Eden grid is simple.\n") ; 
     else
    	  fprintf(fp_log, "Eden grid is body-centered.\n") ; 
     if (anom_flag)
	  fprintf(fp_log, "The anomalous data flag is set.\n") ;
     else
	  fprintf(fp_log, "The anomalous data flag is not set.\n") ;

     fprintf(fp_log, "\n") ;
     fprintf(fp_log, "Run summaries will be written to %s.\n", record_info) ;

}

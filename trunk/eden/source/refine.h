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

				REFINE.H  

  Title:        Include file for 'refine'.
  Author:       Filipe Maia
  Date:         July, 2004
  Function:     

**********************************************************************************/

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_cblas.h>		/* try to circumvent MacOSX bug */
#define CBLAS_INDEX size_t


#define MAX_ATOM_RADIUS 2.5	/* cutoff beyond which the atom should have
			           no effect on the electron density (in A) */
#define MAX_ATOM_RADIUS2 6.2	/* MAX_ATOM_RADIUS squared */
#define MAX_B 99.0

/* These define the precision used for in the electron density tables */
#define BSTEP 2
#define RSTEP 0.005

#define MAX(A,B) ((A)>(B)?(A):(B))
#define MIN(A,B) ((A)>(B)?(B):(A))

typedef struct{
  int position;
  real density_contribution;
  real drhodo;
  real drhodx;
  real drhody;
  real drhodz;
  real drhodb;
}AtomPoint;

typedef struct{
  int Nx;
  int Ny;
  int Nz;
  int dmin[3];
  int dmax[3];
  real cell[6];
  /* matrices to transform between map and real space */
  real ortho2frac[3][3];
  real frac2ortho[3][3];  
  real * ed;
  /* this will hold the real_space coordinates of the grid points */
  real * rs_grid;
}XplorMap;

typedef struct{
  XplorMap xmap;
  /* 
     the rest of the data
     allows us to say which grid points
     are influenced by each atom 
  */
  int natoms;
  /* this is a triple pointer which can be better understood as:
     real * ed_p[natoms][max_grid_pointer_per_atom]
   */
  AtomPoint ** ed_p;
  /* This is used to define the cube around the atoms */
  float grid_cube_dimensions[3];
}CalcMap;

#define REFINE_X 1
#define REFINE_Y 2
#define REFINE_Z 4
#define REFINE_B 8
#define REFINE_O 16

typedef struct{
  XplorMap * xmap;
  CalcMap * cmap;
  int flags;
  int np;
}Params;


void    refine_main(int argc, char *argv[]);
FILE * refine_init(float * tolerance);
int read_xplor(FILE * mapfile,XplorMap * map);
int calc_ed_from_pdb(CalcMap * map,XplorMap * xplor);
void map2ortho(XplorMap * map,real x,real y,real z,real * res);
void ortho2map(XplorMap * map,real x,real y,real z,real * res);
int getZfromSymbol(char * symbol);
int write_xplor(XplorMap * map,char * file);

real get_residual(CalcMap * map,XplorMap * xplor);
void get_all_derivatives(Params * p,gsl_vector *df);
double gsl_residual(const gsl_vector *v, void *params);
void gsl_derivatives(const gsl_vector *v, void *params, gsl_vector *df);
void gsl_residual_and_derivatives(const gsl_vector *v, void *params,double *f, gsl_vector *df);
void run_gsl_minimization(CalcMap * map,XplorMap * xplor,int flags,int atomp,float tolerance);
void free_cmap(CalcMap * map);
void write_pdb(XplorMap * xplor);
real get_RSRfactor(CalcMap * map,XplorMap * xplor);
void correct_pdb_ed(Params * p,real * old_pos,real old_Bj,real old_occup,int n);
void update_atoms(Params * p,const gsl_vector *v,int i);



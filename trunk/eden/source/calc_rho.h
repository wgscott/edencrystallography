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

				CALC_RHO.H

  Title:        Calculate electron density, include file 
  Author:       Filipe Maia
  Date:         July, 2004

**********************************************************************************/

#define SHELLS 300
/* number of total elements */
#define ELEMENTS 118

#define RHO 1
#define DRHODR 2
#define DRHODB 4

typedef struct{
  /* pointers to the actual density */
  real * rho[ELEMENTS];
  real rstep;
  real bstep;
  real hres;
  real lres;
  int nrsteps;
  int nbsteps;
}RhoTable;

/* compute the root mean square of 2 numbers */
real rms(real a,real b);

/* read the form factors coefficients originally from $CLIB/data/atomsf.lib */
void fill_ff_tables(); 

/* return the atomic scattering factor of the atom with atomic number Z */
/* Ions not supported (specially because the PDB parser doesn't recognize them)*/

real scatt_f(real d,int Z,real b);

  /* fourier transform of a solid sphere */

double G(double r, double d);

real g_di(real r,real b,int Z,real d1,real d0,real * result,int flags);

real g_d(real r,real b,int Z,real res_high,real res_low,real * result,int flags);

void init_rho_table(real rstep,real bstep);

int build_rho_table(int Z,real res_high, real res_low);

real lookup_rho(int Z,real b, real r);

void lookup_rho_and_derivatives(int Z,real b, real r, real * result);

/* Decompose a radial derivative into its x y and z components.
   The origin for the position must be the same as the one for the radius */

void decompose(real pos[3],real r, real derivative,real res[3]);



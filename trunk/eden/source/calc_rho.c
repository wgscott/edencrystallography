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

				CALC_RHO.C

  Title:        Calculate electron density; functions used by 'refine'.
  Author:       Filipe Maia
  Date:         July, 2004
  Function:     

void fill_ff_tables()		called by init_rho_table()
real rms()
real scatt_f()
double G()
real g_di()
real g_d()
void init_rho_table()
int build_rho_table()
real lookup_rho()
void lookup_rho_and_derivatives()
void decompose()

******************************************************************************/

#include "util.h"
#include "refine.h"
#include "calc_rho.h"

#define MIN_RADIUS 0.00001
#define MIN_B 0.00001

/* Syntax sugar to access the rho_table */

#define RHOTAB(Z,B,R) (rho_table.rho[(Z)][((B)*rho_table.nrsteps+(R))*3])
#define DRHODRTAB(Z,B,R) (rho_table.rho[(Z)][((B)*rho_table.nrsteps+(R))*3+1])
#define DRHODBTAB(Z,B,R) (rho_table.rho[(Z)][((B)*rho_table.nrsteps+(R))*3+2])

/* Do not change to real, otherwise you have to change the format in the sscanf of fill_ff_tabless */
float atomsf[ELEMENTS][9];

static RhoTable rho_table;

		/**********************************************************
		Use atomsf.lib to obtain details for all atoms in a problem.
		**********************************************************/
void fill_ff_tables(){
  char * edenhome; 
  char path[1024];
  FILE * ff;
  char line[1024];
  char * p;
  int z;
  int i;
  edenhome= getenv("EDENHOME");
  strncpy(path,edenhome,1000);
  strcat(path,"/source/atomsf.lib");
  ff = fopen(path,"r");
  if(!ff){
    EdenError("Could not open atomsf.lib!");
  }
  
  /* Init atomsf array */
  for(i = 0;i<ELEMENTS;i++){
    atomsf[i][0] = -1;
  }
  /* skip the comments */
  for(fgets(line,1024,ff);strstr(line,"AD") == line;fgets(line,1024,ff));

  /* we're at the beginning of an atom record */
  while(line[0]){
    fgets(line,1024,ff);
    /* get Z */
    sscanf(line,"%d",&z);
    if(atomsf[z][0] == -1){
      p = line+23;
      sscanf(p,"%f",&atomsf[z][8]);
      fgets(line,1024,ff);
      p = line;
      sscanf(p,"%f",&atomsf[z][0]);
      p+= 16; 
      sscanf(p,"%f",&atomsf[z][1]);
      p+= 16; 
      sscanf(p,"%f",&atomsf[z][2]);
      p+= 16; 
      sscanf(p,"%f",&atomsf[z][3]);
      fgets(line,1024,ff);
      p = line;
      sscanf(p,"%f",&atomsf[z][4]);
      p+= 16; 
      sscanf(p,"%f",&atomsf[z][5]);
      p+= 16; 
      sscanf(p,"%f",&atomsf[z][6]);
      p+= 16; 
      sscanf(p,"%f",&atomsf[z][7]);
      /* get the last line of the atom record */
      fgets(line,1024,ff);
      /* get the first line of the next atom record */
      line[0] = 0;
      fgets(line,1024,ff);      
    }else{
      /* skip record */
      fgets(line,1024,ff);      
      fgets(line,1024,ff);      
      fgets(line,1024,ff);      
      line[0] = 0;
      fgets(line,1024,ff);            
    }
  }
}
    
      
real rms(real a,real b){
  return sqrt((a*a+b*b)/2);
}

real scatt_f(real d,int Z,real B){
  real res = 0;
  int i;
  for(i = 0;i<4;i++){                        /* --- check this */
    res+= atomsf[Z][i]*exp((-atomsf[Z][i+4]-B)*d*d/4.0);
  }                
  res += atomsf[Z][8]*exp(-B*d*d/4.0);
  return res;    
}

double G(double r,double d){
  return d*d*d*((sin(TWOPI*d*r)-TWOPI*d*r*cos(TWOPI*d*r))/((TWOPI*d*r)*(TWOPI*d*r)*(TWOPI*d*r)));
}

/* Calculate g(d_bar^*_i)*(G(r,d1)-G(r,d0) and its derivatives, 
   according to the flags, and place the result in result, 
   for the shell between d1 and d0, See Acta Crys A, Chapman(1995) 
   equation 10 for the meaning of the symbols.  
    r - radius from atom center 
    b - temperature factor
    Z - atomic number
    d1 - radius of outer sphere
    d0 - radius of inner sphere
    result - array where result is placed
    flag - tells us what to calculate */

real g_di(real r,real b,int Z,real d1,real d0,real * result,int flags)
{
  double G1 = G(r,d1);
  double G0 = G(r,d0);
  real g = scatt_f(rms(d1,d0),Z,b);

  if(flags & RHO)
    result[0] = g*(G1-G0);
  if(flags & DRHODR){		/* equation 25 and 29 of the paper */
                                /* drhodr = g * (dg1dr - dg0dr) */

    result[1] = (((d1*d1*d1 * sin(TWOPI*d1*r)) / (TWOPI*d1*r)) - 3.0 * G1) / r; 
    result[1] -=(((d0*d0*d0 * sin(TWOPI*d0*r)) / (TWOPI*d0*r)) - 3.0 * G0) / r;
    result[1] *= g;
  }
  if(flags & DRHODB){		/* equation 17 and 19 of the paper */
    result[2] = -(((d1*d1 + d0*d0) / 2.0)/4.0)*g*(G1-G0);    
                      /* accuracy critical region ^^^^^ */
  }
  return result[0];
}


/* Do the sum of g_di over all the shells */
real g_d(real r,real b,int Z,real res_high,real res_low,real * result,int flags){
  real d1,d0;
  real sum[3];
  real deltad;
  int i;
  if(r<MIN_RADIUS)
    r = MIN_RADIUS;
  if(b<MIN_B)
    b = MIN_B;
  /* transform the resolution to reciprocal space */
  res_low = 1.0/res_low;
  if(res_high)
    res_high = 1.0/res_high;
  else
    res_high = 1000;
  deltad = (res_high-res_low)/SHELLS;
  /* d0 is the radius of the inner sphere */
  /* d1 is the radius of the outer sphere */
  d0=res_low;
  d1=res_low+deltad;
    if(flags & RHO)
      result[0] = 0;
    if(flags & DRHODR)
      result[1] = 0;
    if(flags & DRHODB)
      result[2]  = 0;
  for ( i=0; i<SHELLS; i++ ) {
    g_di(r,b,Z,d1,d0,sum,flags);
    if(flags & RHO)
      result[0] += sum[0];
    if(flags & DRHODR)
      result[1] += sum[1];
    if(flags & DRHODB)
      result[2] += sum[2];
    d0 += deltad;
    d1 += deltad;
  }
  result[0] *= 2*TWOPI;
  result[1] *= 2*TWOPI;
  result[2] *= 2*TWOPI;
  return result[0];
}

void init_rho_table(real rstep,real bstep){
  int i;
  fill_ff_tables();
  for(i = 0;i<ELEMENTS;i++){
    rho_table.rho[i] = NULL;
  }
  rho_table.rstep = rstep;
  rho_table.bstep = bstep;
}

int build_rho_table(int Z,real res_high, real res_low){
  /*  real b,r;
      int i,j;*/
  /* build tables for atom Z */
  if(Z > ELEMENTS)
    return -1;
  if(rho_table.rho[Z])
    return 0;
  rho_table.nrsteps = (int)(MAX_ATOM_RADIUS/rho_table.rstep)+1;
  rho_table.nbsteps = (int)(MAX_B/rho_table.bstep)+1;
  rho_table.hres = res_high;
  rho_table.lres = res_low;
  /* we have to allocate space for the density and for the 
   derivative with respect to B and r.
   The derivative with respect to Occ is the same as the density
   (because to get the real density we still have to multiply
   by the occupancy)
  */

  rho_table.rho[Z] = (real *)e_malloc(
     sizeof(real)*3*rho_table.nrsteps*rho_table.nbsteps, "build_rho_table");
  /*initialize table*/
  memset(rho_table.rho[Z],0,sizeof(real)*3*rho_table.nrsteps*rho_table.nbsteps);
  
  /* We will fill in the table during the lookups so that it starts more quickly */
  /*

  for(j = 0;j<rho_table.nbsteps;j++){
    b = j*rho_table.bstep;
    for(i =0;i<rho_table.nrsteps;i++){
      r = i*rho_table.rstep;
      g_d(r,b,Z,res_high,res_low, &(rho_table.rho[Z][(j*rho_table.nrsteps+i)*3]),RHO|DRHODR|DRHODB);
    }
  }
  */
  return 0;
}
 
/* There is no error checking here for speed's sake,
  so when you call this function make sure that the
  table is properly built or you'll get a nice seg fault.
  We'll use a simple bilinear interpolation to get the value */

real lookup_rho(int Z,real b, real r)
{
  int r0,b0;
  real t,u;

  b/=rho_table.bstep;
  r/=rho_table.rstep;
  b0 = (int)b;
  r0 = (int)r;
  if(!RHOTAB(Z,b0,r0))
    g_d(r0*rho_table.rstep,b0*rho_table.bstep,Z,rho_table.hres,rho_table.lres, &(rho_table.rho[Z][(b0*rho_table.nrsteps+r0)*3]),RHO|DRHODR|DRHODB);
  if(!RHOTAB(Z,b0+1,r0))
    g_d(r0*rho_table.rstep,(b0+1)*rho_table.bstep,Z,rho_table.hres,rho_table.lres, &(rho_table.rho[Z][((b0+1)*rho_table.nrsteps+r0)*3]),RHO|DRHODR|DRHODB);
  if(!RHOTAB(Z,b0,r0+1))
    g_d((r0+1)*rho_table.rstep,b0*rho_table.bstep,Z,rho_table.hres,rho_table.lres, &(rho_table.rho[Z][(b0*rho_table.nrsteps+r0+1)*3]),RHO|DRHODR|DRHODB);
  if(!RHOTAB(Z,b0+1,r0+1))
    g_d((r0+1)*rho_table.rstep,(b0+1)*rho_table.bstep,Z,rho_table.hres,rho_table.lres, &(rho_table.rho[Z][((b0+1)*rho_table.nrsteps+r0+1)*3]),RHO|DRHODR|DRHODB);

  t = (b-b0);
  u = (r-r0);
  return  (1.0-t) * (1.0-u) * RHOTAB(Z,b0,r0)+
                t * (1.0-u) * RHOTAB(Z,b0+1,r0)+
                      t * u * RHOTAB(Z,b0+1,r0+1)+
                 (1.0-t)* u * RHOTAB(Z,b0,r0+1);
}
  

void lookup_rho_and_derivatives(int Z,real b, real r, real * result)
{
  int r0,b0;
  real t,u;

  b/=rho_table.bstep;
  r/=rho_table.rstep;
  b0 = (int)b;
  r0 = (int)r;
  if(!RHOTAB(Z,b0,r0))
    g_d(r0*rho_table.rstep,b0*rho_table.bstep,Z,rho_table.hres,rho_table.lres, &(rho_table.rho[Z][(b0*rho_table.nrsteps+r0)*3]),RHO|DRHODR|DRHODB);
  if(!RHOTAB(Z,b0+1,r0))
    g_d(r0*rho_table.rstep,(b0+1)*rho_table.bstep,Z,rho_table.hres,rho_table.lres, &(rho_table.rho[Z][((b0+1)*rho_table.nrsteps+r0)*3]),RHO|DRHODR|DRHODB);
  if(!RHOTAB(Z,b0,r0+1))
    g_d((r0+1)*rho_table.rstep,b0*rho_table.bstep,Z,rho_table.hres,rho_table.lres, &(rho_table.rho[Z][(b0*rho_table.nrsteps+r0+1)*3]),RHO|DRHODR|DRHODB);
  if(!RHOTAB(Z,b0+1,r0+1))
    g_d((r0+1)*rho_table.rstep,(b0+1)*rho_table.bstep,Z,rho_table.hres,rho_table.lres, &(rho_table.rho[Z][((b0+1)*rho_table.nrsteps+r0+1)*3]),RHO|DRHODR|DRHODB);

  t = (b-b0);
  u = (r-r0);
  *(result) = (1.0-t)*(1.0-u)*RHOTAB(Z,b0,r0)+
                    t*(1.0-u)*RHOTAB(Z,b0+1,r0)+
                          t*u*RHOTAB(Z,b0+1,r0+1)+
                    (1.0-t)*u*RHOTAB(Z,b0,r0+1);
  *(result+1) = (1.0-t)*(1.0-u)*DRHODRTAB(Z,b0,r0)+
                      t*(1.0-u)*DRHODRTAB(Z,b0+1,r0)+
                            t*u*DRHODRTAB(Z,b0+1,r0+1)+
                      (1.0-t)*u*DRHODRTAB(Z,b0,r0+1);
  *(result+2) = (1.0-t)*(1.0-u)*DRHODBTAB(Z,b0,r0)+
                      t*(1.0-u)*DRHODBTAB(Z,b0+1,r0)+
                            t*u*DRHODBTAB(Z,b0+1,r0+1)+
                      (1.0-t)*u*DRHODBTAB(Z,b0,r0+1);
}

void decompose(real pos[3],real r, real derivative,real res[3])
{
  if(r< MIN_RADIUS)
    r = MIN_RADIUS;
  res[0] = (pos[0]/r)*derivative;
  res[1] = (pos[1]/r)*derivative;
  res[2] = (pos[2]/r)*derivative;
}



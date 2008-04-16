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

				REFINE.C

  Title:        Do real-space refinement
  Author:       Filipe Maia
  Date:         July, 2004
  Function:     Find the best pdb file that reflects the electron density.

**********************************************************************************/

#include "eden.h"
#include "pdb.h"
#include "reg.h"
#include "calc_rho.h"
#include "refine.h"
#ifdef DEBUG
#include <sys/time.h>
#endif 

int Zj[MAXNATOMS];
char fullatid[MAXNATOMS][6];
int max_pp_atom;		/*pp  stands for points_per  */


void  refine_main(int argc, char *argv[]){
  FILE * mapname;   
  XplorMap map;
  CalcMap c_map;
  /* This controls how hard the minimizer will try 
     keyword RS_TOL
  */
  float tolerance = 0.15;

  c_map.ed_p = NULL;
  c_map.xmap.ed = NULL;
  c_map.xmap.rs_grid = NULL;

  hello(caller) ;

  mapname = refine_init(&tolerance) ;
    
  read_xplor(mapname,&map);
  
  run_gsl_minimization(&c_map,&map,REFINE_X|REFINE_Y|REFINE_Z|REFINE_O|REFINE_B,5,tolerance);
      
  sprintf(message,"Final RS R factor - %f\n",get_RSRfactor(&c_map,&map));
  printTwice(message) ;

  write_pdb(&map);
  return ;
}

FILE * refine_init(float * tolerance)
{
  void fetchPdbInfo(char *filename);
  char  xplor_filename[1024] = "";
  FILE * map = NULL;
  int k,i;
  pdb_filename[0] = 0;

  readBasicInput();
  if(low_res_cutoff == 1000){
    printTwice("Low resolution cutoff not specified; using 1000 A.");
  }
  for (k = 0; k < Nlines; k++) {
    strcpy(nextline, allinp + k*MAXSTRING) ;
    sscanf(nextline, "%s", id) ;
    
    if (strcmp(id, "XPLOR_FILENAME") == 0)  {
      sscanf(nextline, "%*s %s", xplor_filename) ;
    } else if ((strcmp(id, "PDB_FILENAME") == 0) ||
	       (strcmp(id, "PDB_FN") == 0)){
      sscanf(nextline, "%*s %s", pdb_filename) ;
    } else if ((strcmp(id, "RS_TOL") == 0)){
      sscanf(nextline, "%*s %f", tolerance);
    }  
    
  }
  if(!xplor_filename[0])
    EdenError("Missing xplor file name!") ;
 
  if(!pdb_filename[0])
    EdenError("Missing pdb file name!") ;
  
  init_rho_table(RSTEP,BSTEP);
  fetchPdbInfo(pdb_filename);

  for(i = 0;i<Nin;i++){
    Zj[i] = getZfromSymbol(atid[i]);
    if(Zj[i])
      build_rho_table(Zj[i],input_res,low_res_cutoff);
    else
      EdenError("Couldn't identify input atom chemical symbol.\nBear in mind that eden still doesn't correctly support ions or other strange entities");
  }
  map = fopen(xplor_filename,"r");
  if(!map) {
    sprintf(message, "Error opening %s", xplor_filename) ;
    EdenError(message) ;
  }
  

  return map;
}



	/* Read_xplor should return non 0 on error */

int read_xplor(FILE * mapfile,XplorMap * map)
{
  char line[1024];
  int imult;
  int jmult;
  int kmult;
  int imax,jmax,kmax;
  char * p = NULL;
  int n,i,j,k;
  real a4, a5 ;
  char tmp;

  /* read title, skip remarks */

  fgets(line,1024,mapfile);
  fgets(line,1024,mapfile);
  for(fgets(line,1024,mapfile);
  strstr(line,"REMARKS");
  fgets(line,1024,mapfile));

  sscanf(line," %7d %7d %7d %7d %7d %7d %7d %7d %7d",
	 &(map->Nx),&(map->dmin[0]),&(map->dmax[0]),
	 &(map->Ny),&(map->dmin[1]),&(map->dmax[1]),
	 &(map->Nz),&(map->dmin[2]),&(map->dmax[2]));
  fgets(line,1024,mapfile);

  sscanf(line,"%f %f %f %f %f %f",&t1, &t2, &t3, &t4, &t5, &t6) ;
  map->cell[0] = t1 ;
  map->cell[1] = t2 ;
  map->cell[2] = t3 ;
  map->cell[3] = t4 ;
  map->cell[4] = t5 ;
  map->cell[5] = t6 ;

  a4 = cos(DTOR * map->cell[4]) ;
  a5 = cos(DTOR * map->cell[5]) ;

  /* The column are the cell edge vectors */
    map->frac2ortho[0][0] = map->cell[0];
    map->frac2ortho[1][0] = 0;
    map->frac2ortho[2][0] = 0;
    map->frac2ortho[0][1] = cos(DTOR * map->cell[5])*map->cell[1];
    map->frac2ortho[1][1] = sin(DTOR * map->cell[5])*map->cell[1];
    map->frac2ortho[2][1] = 0;
    map->frac2ortho[0][2] = map->cell[2]*a4;
    map->frac2ortho[1][2] = map->cell[2]* (cos(DTOR * map->cell[3])- a4*a5)/
                                           sin(DTOR * map->cell[5]);
    map->frac2ortho[2][2] = sqrt(map->cell[2]*map->cell[2]-
                                 map->frac2ortho[0][2]*map->frac2ortho[0][2]-
                                 map->frac2ortho[1][2]*map->frac2ortho[1][2]);

    matinv(map->frac2ortho,map->ortho2frac);

  /* Internally maps will always be stored in ZYX format, meaning 
     that X is the index that changes fastest.  However, this is 
     not necessarily the case for the input.
     Because of that we need to do a little bit of index trickery. */

  fgets(line,1024,mapfile);
  if(line[0] == 'Z'){		 /* The usual*/
    imult = map->Nx*map->Ny;
    imax = map->Nz;
  }else if(line[0] == 'Y'){
    imult = map->Nx;
    imax = map->Ny;
  }else if(line[0] == 'X'){
    imult = 1;
    imax = map->Nx;
  }else{
    return 1;
  }

  if(line[1] == 'Z'){
    jmult = map->Nx*map->Ny;
    jmax = map->Nz;
  }else if(line[1] == 'Y'){	 /* The usual*/
    jmult = map->Nx;
    jmax = map->Ny;
  }else if(line[1] == 'X'){
    jmult = 1;
    jmax = map->Nx;
  }else{
    return 1;
  }

  if(line[2] == 'Z'){
    kmult = map->Nx*map->Ny;
    kmax = map->Nz;
  }else if(line[2] == 'Y'){
    kmult = map->Nx;
    kmax = map->Ny;
  }else if(line[2] == 'X'){	 /* The usual*/
    kmult = 1;
    kmax = map->Nx;
  }else{
    return 1;
  }
  n = 0;
  map->ed = (real *)e_malloc(sizeof(real)*map->Nx*map->Ny*map->Nz, 
              "read_xplor");
  map->rs_grid = (real *)e_malloc(sizeof(real)*map->Nx*map->Ny*map->Nz*3, 
              "read_xplor");
  for(i = 0;i<imax;i++){
    /* get the layer number line*/
    fgets(line,1024,mapfile);      
    n = 0;
    for(j = 0;j<jmax;j++){
      for(k = 0;k<kmax;k++){
	if(n % 6 == 0){
	  /*we're done with this line.
	   Read the next one and position pointer */
	  fgets(line,1024,mapfile);      
	  p = line;
	}
	tmp = p[12];
	p[12] = 0;

	sscanf(p,"%g",&t1) ;
	map->ed[i*imult+j*jmult+k*kmult] = t1 ;

	p[12] = tmp;
	p+=12;
	n++;
      }
    }
  }
  for(k = 0;k<map->Nz;k++){
    for(j = 0;j<map->Ny;j++){
      for(i = 0;i<map->Nx;i++){
	map2ortho(map,i,j,k,&(map->rs_grid[(k*map->Ny*map->Nx+j*map->Nx+i)*3]));
      }
    }
  }
  return 0;
}
  
/* Calculate a density map with the same properties
   as the xplor map from the PDB information */

int calc_ed_from_pdb(CalcMap *map, XplorMap *xplor)
{
  int i,j,k,n;
  int ti,tj,tk;
  real ntreal_pos[3];
  unsigned int wrap_flag;

  int imult,jmult,kmult;
  int mini,minj,mink;
  int maxi,maxj,maxk; 
  real grid_point_density;
  real grid_pos[3];
  real *real_pos;
  real rel_pos[3];
  real dist;
  int grid_pp_atom;		/*pp  stands for points_per  */
  real tmp[3];
  
#ifdef DEBUG
  struct timezone * tz = NULL;
  struct timeval  tvi;
  struct timeval  tvf;
  gettimeofday(&tvi, tz);
  fprintf(stderr,"calc_ed called\n");
#endif
  map->xmap = *xplor;

  for(i = 0;i<3;i++){
    map->xmap.dmin[i] = xplor->dmin[i];
    map->xmap.dmax[i] = xplor->dmax[i];
    for(j = 0;j<3;j++){
      map->xmap.frac2ortho[i][j] = xplor->frac2ortho[i][j];
      map->xmap.ortho2frac[i][j] = xplor->ortho2frac[i][j];
    }
  }
  for(i = 0;i<6;i++){
    map->xmap.cell[i]= xplor->cell[i];
  }
  imult = map->xmap.Ny*map->xmap.Nx;
  jmult = map->xmap.Nx;
  kmult = 1;
  grid_point_density = map->xmap.Nx/xplor->cell[0]*
                       map->xmap.Ny/xplor->cell[1]*
                       map->xmap.Nz/xplor->cell[2];

  /*product of the number of points in the x, y and z directions */

  max_pp_atom =  ceil(MAX_ATOM_RADIUS*2*map->xmap.Nx/
                              xplor->frac2ortho[0][0]+1);
  max_pp_atom *= ceil(MAX_ATOM_RADIUS*2*map->xmap.Ny/
                              xplor->frac2ortho[1][1]+1);
  max_pp_atom *= ceil(MAX_ATOM_RADIUS*2*map->xmap.Nz/
                              xplor->frac2ortho[2][2]+1);

  map->xmap.ed = (real *)e_malloc(sizeof(real)*map->xmap.Nx
                        *map->xmap.Ny*map->xmap.Nz, "calc_ed_from_pdb");
  /* make sure the density is all 0
     Lets hope that 0 mean a float point 0 in all architectures 
  */
  memset(map->xmap.ed,0,sizeof(real)*map->xmap.Nx*map->xmap.Ny*map->xmap.Nz);

  sprintf(message, "calc_ed_from_pdb");
  map->natoms = Nin;
  map->ed_p = (AtomPoint **)e_malloc(sizeof(AtomPoint *)*map->natoms, message) ;
  map->grid_cube_dimensions[0] = xplor->frac2ortho[0][0]/map->xmap.Nx;
  map->grid_cube_dimensions[1] = xplor->frac2ortho[1][1]/map->xmap.Ny;
  map->grid_cube_dimensions[2] = xplor->frac2ortho[2][2]/map->xmap.Nz;

  /* We can't compute the derivatives until we have the complete density
     so we'll need to do 2 passes */

  for(n = 0;n<map->natoms;n++){
    ortho2map(&map->xmap,(real)pos[n].x,(real)pos[n].y,(real)pos[n].z,grid_pos);
    grid_pp_atom = 0;

    map->ed_p[n] = (AtomPoint *)
                   e_malloc(sizeof(AtomPoint)*max_pp_atom,message) ;
    mini = floor(grid_pos[2]-(MAX_ATOM_RADIUS/map->grid_cube_dimensions[2]+1));
    maxi = ceil(grid_pos[2]+(MAX_ATOM_RADIUS/map->grid_cube_dimensions[2]))+1;
    minj = floor(grid_pos[1]-(MAX_ATOM_RADIUS/map->grid_cube_dimensions[1]+1));
    maxj = ceil(grid_pos[1]+(MAX_ATOM_RADIUS/map->grid_cube_dimensions[1]))+1;
    mink = floor(grid_pos[0]-(MAX_ATOM_RADIUS/map->grid_cube_dimensions[0]+1));
    maxk = ceil(grid_pos[0]+(MAX_ATOM_RADIUS/map->grid_cube_dimensions[0]))+1;
    for(i = mini;i<maxi;i++){
      ti = i;
      wrap_flag = 0;
      while(ti<0){
	ti+= map->xmap.Nz;
	wrap_flag |= 1;
      }
      while(ti>=map->xmap.Nz){
	ti-= map->xmap.Nz;	
	wrap_flag |= 1;
      }
      for(j = minj;j<maxj;j++){
	tj = j;
	wrap_flag &= ~2;
	while(tj<0){
	  tj+= map->xmap.Ny;
	  wrap_flag |= 2;
	}
	while(tj>=map->xmap.Ny){
	  tj-= map->xmap.Ny;
	  wrap_flag |= 2;
	}
	for(k = mink;k<maxk;k++){
	  tk = k;
	  wrap_flag &= ~4;
	  while(tk<0){
	    tk+= map->xmap.Nx;
	    wrap_flag |= 4;
	  }
	  while(tk>=map->xmap.Nx){
	    tk-= map->xmap.Nx;
	    wrap_flag |= 4;
	  }
	  if(wrap_flag){
	    map2ortho(&map->xmap,k,j,i,ntreal_pos);
	    real_pos = ntreal_pos;
	  }else{
	    real_pos = &(map->xmap.rs_grid[(i*map->xmap.Ny*map->xmap.Nx+j*map->xmap.Nx+k)*3]);
	  }
	  dist = (real_pos[0]-pos[n].x)*(real_pos[0]-pos[n].x)+
	         (real_pos[1]-pos[n].y)*(real_pos[1]-pos[n].y)+
	         (real_pos[2]-pos[n].z)*(real_pos[2]-pos[n].z);

	  if(dist < MAX_ATOM_RADIUS2){
	    dist = sqrt(dist);
	    map->ed_p[n][grid_pp_atom].position = ti*imult+tj*jmult+tk*kmult;

	    if(grid_pp_atom == max_pp_atom)
	      EdenError("Not enough space for the points per atom lists!");
	    
	    lookup_rho_and_derivatives(Zj[n],Bj[n],dist,tmp);
	    map->ed_p[n][grid_pp_atom].drhodo = tmp[0];
	    tmp[0]*= occup[n];
	    tmp[1]*= occup[n];
	    tmp[2]*= occup[n];
	    map->xmap.ed[ti*imult+tj*jmult+tk*kmult] += tmp[0];
	    map->ed_p[n][grid_pp_atom].density_contribution = tmp[0];
	    map->ed_p[n][grid_pp_atom].drhodb = tmp[2];
	    
	    rel_pos[0] = real_pos[0]-pos[n].x;
	    rel_pos[1] = real_pos[1]-pos[n].y;
	    rel_pos[2] = real_pos[2]-pos[n].z;
	    decompose(rel_pos,dist,tmp[1],tmp);
	    map->ed_p[n][grid_pp_atom].drhodx = tmp[0];
	    map->ed_p[n][grid_pp_atom].drhody = tmp[1];
	    map->ed_p[n][grid_pp_atom].drhodz = tmp[2];
	    grid_pp_atom++;
	  }
	}
      }
    }
    map->ed_p[n][grid_pp_atom].position = -1;
  }
  /* second pass to compute the derivatives.  at least this time 
     we already have the points which are part of the atom */

#ifdef DEBUG
  gettimeofday(&tvf, tz);
  fprintf(stderr,"Time used in calc_ed - %.1f ms\n",(tvf.tv_sec-tvi.tv_sec)*1000+(tvf.tv_usec-tvi.tv_usec)/1000.0);
#endif
  return 0;
}

void correct_pdb_ed(Params * p,real * old_pos,real old_Bj,real old_occup,int n){
  int i,j,k;
  int ti,tj,tk;
  real ntreal_pos[3];
  unsigned int wrap_flag;
  int grid_pp_atom;		/*pp  stands for points_per  */
  int imult,jmult,kmult;
  int mini,minj,mink;
  int maxi,maxj,maxk; 
  real grid_pos[3];
  real *real_pos;
  real rel_pos[3];
  real dist;
  real tmp[3];

  i = 0;
  /* First remove the old density */
  while(p->cmap->ed_p[n] && p->cmap->ed_p[n][i].position!= -1) {
    p->cmap->xmap.ed[p->cmap->ed_p[n][i].position] -= p->cmap->ed_p[n][i].density_contribution;
    i++;
  }
  
  imult = p->cmap->xmap.Ny*p->cmap->xmap.Nx;
  jmult = p->cmap->xmap.Nx;
  kmult = 1;
    

  /* Now add the new density */
  ortho2map(&p->cmap->xmap,(real)pos[n].x,(real)pos[n].y,(real)pos[n].z,grid_pos);
  grid_pp_atom = 0;
  
  mini = floor(grid_pos[2]-(MAX_ATOM_RADIUS/p->cmap->grid_cube_dimensions[2]+1));
  maxi = ceil(grid_pos[2]+(MAX_ATOM_RADIUS/p->cmap->grid_cube_dimensions[2]))+1;
  minj = floor(grid_pos[1]-(MAX_ATOM_RADIUS/p->cmap->grid_cube_dimensions[1]+1));
  maxj = ceil(grid_pos[1]+(MAX_ATOM_RADIUS/p->cmap->grid_cube_dimensions[1]))+1;
  mink = floor(grid_pos[0]-(MAX_ATOM_RADIUS/p->cmap->grid_cube_dimensions[0]+1));
  maxk = ceil(grid_pos[0]+(MAX_ATOM_RADIUS/p->cmap->grid_cube_dimensions[0]))+1;
  for(i = mini;i<maxi;i++){
    ti = i;
    wrap_flag = 0;
    while(ti<0){
      ti+= p->cmap->xmap.Nz;
      wrap_flag |= 1;
    }
    while(ti>=p->cmap->xmap.Nz){
      ti-= p->cmap->xmap.Nz;	
      wrap_flag |= 1;
    }
    for(j = minj;j<maxj;j++){
      tj = j;
      wrap_flag &= ~2;
      while(tj<0){
	tj+= p->cmap->xmap.Ny;
	wrap_flag |= 2;
      }
      while(tj>=p->cmap->xmap.Ny){
	tj-= p->cmap->xmap.Ny;
	wrap_flag |= 2;
      }
      for(k = mink;k<maxk;k++){
	tk = k;
	wrap_flag &= ~4;
	while(tk<0){
	  tk+= p->cmap->xmap.Nx;
	  wrap_flag |= 4;
	}
	while(tk>=p->cmap->xmap.Nx){
	  tk-= p->cmap->xmap.Nx;
	  wrap_flag |= 4;
	}
	if(wrap_flag){
	  map2ortho(&p->cmap->xmap,k,j,i,ntreal_pos);
	  real_pos = ntreal_pos;
	}else{
	  real_pos = &(p->cmap->xmap.rs_grid[(i*p->cmap->xmap.Ny*p->cmap->xmap.Nx+j*p->cmap->xmap.Nx+k)*3]);
	}
	dist = (real_pos[0]-pos[n].x)*(real_pos[0]-pos[n].x)+
	  (real_pos[1]-pos[n].y)*(real_pos[1]-pos[n].y)+
	  (real_pos[2]-pos[n].z)*(real_pos[2]-pos[n].z);
	
	if(dist < MAX_ATOM_RADIUS2){
	  dist = sqrt(dist);
	  p->cmap->ed_p[n][grid_pp_atom].position = ti*imult+tj*jmult+tk*kmult;
	  
	  if(grid_pp_atom == max_pp_atom)
	    EdenError("Not enough space for the points per atom lists!");
	  
	  lookup_rho_and_derivatives(Zj[n],Bj[n],dist,tmp);
	  p->cmap->ed_p[n][grid_pp_atom].drhodo = tmp[0];
	  tmp[0]*= occup[n];
	  tmp[1]*= occup[n];
	  tmp[2]*= occup[n];
	  p->cmap->xmap.ed[ti*imult+tj*jmult+tk*kmult] += tmp[0];
	  p->cmap->ed_p[n][grid_pp_atom].density_contribution = tmp[0];
	  p->cmap->ed_p[n][grid_pp_atom].drhodb = tmp[2];
	  
	  rel_pos[0] = real_pos[0]-pos[n].x;
	  rel_pos[1] = real_pos[1]-pos[n].y;
	  rel_pos[2] = real_pos[2]-pos[n].z;
	  decompose(rel_pos,dist,tmp[1],tmp);
	  p->cmap->ed_p[n][grid_pp_atom].drhodx = tmp[0];
	  p->cmap->ed_p[n][grid_pp_atom].drhody = tmp[1];
	  p->cmap->ed_p[n][grid_pp_atom].drhodz = tmp[2];
	  grid_pp_atom++;
	}
      }
    }
  }
  p->cmap->ed_p[n][grid_pp_atom].position = -1;
}

real get_residual(CalcMap * map,XplorMap * xplor){
  real cost = 0;
  int i;

  for(i = 0;i<xplor->Nz*xplor->Ny*xplor->Nx;i++){
    cost+= (xplor->ed[i]-map->xmap.ed[i])*(xplor->ed[i]-map->xmap.ed[i]);
  }
  return cost;
} 

/* The following function is used ONLY for reporting progress in refinement.
   Originally, it computed sums of squares in both numerator and denominator.*/

real get_RSRfactor(CalcMap * map,XplorMap * xplor)
{
  real cost = 0;
  real denominator=0;
  int i;

  for(i = 0;i<xplor->Nz*xplor->Ny*xplor->Nx;i++){
    cost       +=fabs(xplor->ed[i]-map->xmap.ed[i]);
    denominator+=fabs(xplor->ed[i]+map->xmap.ed[i]);
  }  
  return cost/(denominator);
} 

void get_all_derivatives(Params * p,gsl_vector *df){
  /*  real grid_pos[3];*/

  int n,i;
  real dd,dx,dy,dz,db,docc;
  for(n = 0;n<p->cmap->natoms;n++){

    /*    ortho2map(&p->cmap->xmap,
                    (real)pos[n].x,(real)pos[n].y,(real)pos[n].z,grid_pos);*/
    i = 0;
    dx=dy=dz=db=docc= 0;
    while(p->cmap->ed_p[n] && p->cmap->ed_p[n][i].position!= -1)
    {
	dd = p->xmap->ed[p->cmap->ed_p[n][i].position]
       -p->cmap->xmap.ed[p->cmap->ed_p[n][i].position] ;

	dx -= 2 * dd * p->cmap->ed_p[n][i].drhodx;
	dy -= 2 * dd * p->cmap->ed_p[n][i].drhody;
	dz -= 2 * dd * p->cmap->ed_p[n][i].drhodz;
	db -= 2 * dd * p->cmap->ed_p[n][i].drhodb;
	docc -= 2 * dd * p->cmap->ed_p[n][i].drhodo;
        i++;
    }
    i = 0;
    if(p->flags & REFINE_X){
      gsl_vector_set(df, n*p->np+i, -dx);
      i++;
    }
    if(p->flags & REFINE_Y){
      gsl_vector_set(df, n*p->np+i, -dy);
      i++;
    }
    if(p->flags & REFINE_Z){
      gsl_vector_set(df, n*p->np+i, -dz);
      i++;
    }
    if(p->flags & REFINE_B){
      gsl_vector_set(df, n*p->np+i, db);    
      i++;
    }
    if(p->flags & REFINE_O){
      if(docc>0 && occup[n] >= 1)
	docc = 0;
      gsl_vector_set(df, n*p->np+i, docc);    
      i++;
    }
  }  
}


void get_derivative(Params * p,gsl_vector *df,int n){
  /*  real grid_pos[3];*/

  int i;
  real dd,dx,dy,dz,db,docc;
  /*    ortho2map(&p->cmap->xmap,
	(real)pos[n].x,(real)pos[n].y,(real)pos[n].z,grid_pos);*/
  i = 0;
  dx=dy=dz=db=docc= 0;
  while(p->cmap->ed_p[n] && p->cmap->ed_p[n][i].position!= -1)
    {
      dd = p->xmap->ed[p->cmap->ed_p[n][i].position]
	-p->cmap->xmap.ed[p->cmap->ed_p[n][i].position] ;
      
      dx -= 2 * dd * p->cmap->ed_p[n][i].drhodx;
      dy -= 2 * dd * p->cmap->ed_p[n][i].drhody;
      dz -= 2 * dd * p->cmap->ed_p[n][i].drhodz;
      db -= 2 * dd * p->cmap->ed_p[n][i].drhodb;
      docc -= 2 * dd * p->cmap->ed_p[n][i].drhodo;
      i++;
    }
  i = 0;
  if(p->flags & REFINE_X){
    gsl_vector_set(df, n*p->np+i, -dx);
    i++;
  }
  if(p->flags & REFINE_Y){
    gsl_vector_set(df, n*p->np+i, -dy);
    i++;
  }
  if(p->flags & REFINE_Z){
    gsl_vector_set(df, n*p->np+i, -dz);
    i++;
  }
  if(p->flags & REFINE_B){
    gsl_vector_set(df, n*p->np+i, db);    
    i++;
  }
  if(p->flags & REFINE_O){
    if(docc>0 && occup[n] >= 1)
      docc = 0;
    gsl_vector_set(df, n*p->np+i, docc);    
    i++;
  }
}

double gsl_residual(const gsl_vector *v, void *params){
  int i;

  Params * p = (Params *)params;
  /* 
     the gsl_vector is a list of the position bfactors of the atoms x
     in the form of x1,y1,z1,b1,x2,y2,z2,b2...xn,yn,zn,bn 
  */

  /* params contains a pointer to CalcMap and to XplorMap*/

  /* update atom properties */
  for(i= 0;i<Nin;i++){
    update_atoms(p,v,i);
  }    
  
  return get_residual(p->cmap,p->xmap);
}

void gsl_derivatives(const gsl_vector *v, void *params, gsl_vector *df){
  int i;
  Params * p = (Params *)params;
  /* 
     the gsl_vector is a list of the position bfactors of the atoms x
     in the form of x1,y1,z1,b1,x2,y2,z2,b2...xn,yn,zn,bn 
  */

  /* params contains a pointer to CalcMap and to XplorMap*/

  /* update atom properties */
  for(i= 0;i<Nin;i++){
    update_atoms(p,v,i);
  }
  get_all_derivatives(p,df);    
}

void gsl_residual_and_derivatives(const gsl_vector *v, void *params,double *f, gsl_vector *df){

  int i;
  Params * p = (Params *)params;
  /* 
     the gsl_vector is a list of the position bfactors of the atoms x
     in the form of x1,y1,z1,b1,x2,y2,z2,b2...xn,yn,zn,bn 
  */

  /* params contains a pointer to CalcMap and to XplorMap*/

  /* update atom properties */
  for(i= 0;i<Nin;i++){
    update_atoms(p,v,i);
  }

  get_all_derivatives(p,df);  
  *f = get_residual(p->cmap,p->xmap);
}

void map2ortho(XplorMap * map,real x,real y,real z,real * res){
  int a;
  real frac[3];
  frac[0] = x/map->Nx;
  frac[1] = y/map->Ny;
  frac[2] = z/map->Nz;
  /* Lets hope this comparison with floats works */
  if(map->cell[3] == 90 && map->cell[4] == 90 && map->cell[5] == 90){
    /* simple orthorhombic case */
    for(a = 0;a<3;a++){
      res[a] = frac[a]*map->cell[a];
    }
  }else{
    matvecmult(map->frac2ortho,frac,res);
  }
}

void ortho2map(XplorMap * map,real x,real y,real z,real * res){
  real xyz[3];
  /* Lets hope this comparison with floats works */
  if(map->cell[3] == 90 && map->cell[4] == 90 && map->cell[5] == 90){
    /* simple orthorhombic case */
    res[0] = x*map->Nx/map->cell[0];
    res[1] = y*map->Ny/map->cell[1];
    res[2] = z*map->Nz/map->cell[2];
  }else{
    xyz[0] = x;
    xyz[1] = y;
    xyz[2] = z;
    matvecmult(map->ortho2frac,xyz,res);
    res[0] *= map->Nx;
    res[1] *= map->Ny;
    res[2] *= map->Nz;
  }
}
/* write_xplor() is essentially the same as printXmap() in regutil.c,
   except for the data structures.                              ***/

int write_xplor(XplorMap * map,char * file){
  FILE * fp;
  int z,i,j,n;
  fp = fopen(file,"w");
  if(!fp)
    EdenError("Unable to write xplor file.");
  
  fprintf(fp,"\n       2 !NTITLE\n") ; 
  fprintf(fp," REMARKS Electron Density Information from EDEN\n") ;
  fprintf(fp," REMARKS DATE: %s \n", timestamp()) ;
  fprintf(fp," %7d %7d %7d %7d %7d %7d %7d %7d %7d\n",
	  map->Nx, map->dmin[0], map->dmax[0],
	  map->Ny, map->dmin[1], map->dmax[1],
	  map->Nz, map->dmin[2], map->dmax[2]);
     fprintf(fp,"%12.5E%12.5E%12.5E%12.5E%12.5E%12.5E\n",
	     map->cell[0],map->cell[1],map->cell[2],map->cell[3],map->cell[4],map->cell[5]);
     fprintf(fp, "ZYX\n") ;

     z = map->dmin[2];
     for (n = 0; n < map->Nz; n++, z++) {
        fprintf(fp, "%8d\n", z) ;
        for (i = 0; i < map->Nx*map->Ny; i += 6) {
         for (j = 0; j < 6; j++)
	   if (i+j < map->Nx*map->Ny)
                fprintf(fp, "%12.5E", map->ed[n*map->Nx*map->Ny+i+j]);
	 fprintf(fp, "\n") ;
        }
     }
     fprintf(fp, "   -9999\n") ;
     fclose(fp) ;
  return 0;
}

void run_gsl_minimization(CalcMap * map,XplorMap * xplor,int flags,int atomp,float tolerance){

  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *s = NULL;
  Params param;
  gsl_vector *x;
  gsl_vector *grad;
  real grad_norm = 0;
  gsl_multimin_function_fdf my_func;
  int status;
  int np = Nin*atomp;
  int n,i, Nstop;
  int iter = 0;

  param.cmap = map;
  param.xmap = xplor;
  param.flags = flags;
  param.np = atomp;

  /* calculate initial density */
  calc_ed_from_pdb(map,xplor);

  my_func.f = &gsl_residual;
  my_func.df = &gsl_derivatives;
  my_func.fdf = &gsl_residual_and_derivatives;
  my_func.n = np;
  my_func.params = &param;
  x = gsl_vector_alloc (np);
  for(n = 0;n<Nin;n++){
    i = 0;
    if(flags & REFINE_X){
      gsl_vector_set (x, n*atomp+i, pos[n].x);
      i++;
    }
    if(flags & REFINE_Y){
      gsl_vector_set (x,n*atomp+i, pos[n].y);
      i++;
    }
    if(flags & REFINE_Z){
      gsl_vector_set (x, n*atomp+i, pos[n].z);
      i++;
    }
    if(flags & REFINE_B){
      gsl_vector_set (x, n*atomp+i, Bj[n]);
      i++;
    }
    if(flags & REFINE_O){
      gsl_vector_set (x, n*atomp+i, MIN(occup[n],1));
      i++;
    }
  }

  /*  T = gsl_multimin_fdfminimizer_conjugate_fr;*/

  T = gsl_multimin_fdfminimizer_vector_bfgs;
  s = gsl_multimin_fdfminimizer_alloc (T, np);

  gsl_multimin_fdfminimizer_set (s, &my_func, x, 0.2, 1e-4);
  sprintf(message,"Initial RS R factor - %f\n",get_RSRfactor(map,xplor));
  printTwice(message) ;
  Nstop = (int)(100*(log(Nin)+1)) ;
  sprintf(message, 
    "At most, there will be %d iterations\n", Nstop) ;
  printTwice(message) ;

  fprintf (fp_log, " it # |rho_obs-rho_pdb|  gradient \n") ;
  fprintf (fp_log, "----- ----------------- ----------\n") ;
  do{
    iter++;
    if(iter % 10 == 0){
      write_pdb(xplor);
      sprintf(message,
        "Iter # = %d, RS R factor = %f",iter, get_RSRfactor(map,xplor));      
      printTwice(message) ;
    }
    status = gsl_multimin_fdfminimizer_iterate(s);
    
/* gsl_errno.h in /sw/include/gsl (on MacOSX)
   defines GSL_SUCCESS = 0; DSL_FAILURE = -1. GSL_CONTINUE = -2
   and various positive values (1 - 32) for various problems.    */

    if (status)
      break;
    grad = gsl_multimin_fdfminimizer_gradient(s);
    for(i = 0;i<np;i++){
      grad_norm += gsl_vector_get(grad, i)*gsl_vector_get(grad, i);
    }
    grad_norm = sqrt(grad_norm);
    fprintf (fp_log,"%5d %17.5f %10.5f\n", iter, s->f,grad_norm);
    fflush (fp_log) ;
    status = gsl_multimin_test_gradient (s->gradient, np*tolerance);
  }
  while (status == GSL_CONTINUE && iter < Nstop) ;
  if(Nstop <= iter){
    sprintf(message,
            "Maximum number of iterations reached!\nFinal gradient = %f",grad_norm);
    printTwice(message) ;
  } else if(status == GSL_FAILURE){
    sprintf(message,
	    "Minimization failled! Final gradient = %f",grad_norm);
    printTwice(message) ;
  } else if(status == GSL_SUCCESS){
    sprintf (message, "Minimum found! Final gradient = %f",grad_norm);
    printTwice(message) ;
  }else if(status == GSL_ENOPROG){
    sprintf (message, "Minimization finished. No longer making progress towards solution.\nFinal gradient = %f",
	     grad_norm);
    printTwice(message) ;
  }else{
    sprintf (message, "GSL error: %d\nPlease check gsl_errno.h for the meaning.",status);
    printTwice(message) ;
  }

  write_pdb(xplor);
  gsl_multimin_fdfminimizer_free (s);
  gsl_vector_free (x);

}

void free_cmap(CalcMap * map){
  int i;
  if(map){
    if(map->ed_p){
      i = 0;
      while(i<Nin){
	if(map->ed_p[i])
	  e_free(sizeof(AtomPoint *)*max_pp_atom, map->ed_p[i]) ;
	i++;
      }
      e_free(sizeof(AtomPoint *)*map->natoms, map->ed_p) ;
    }
    if(map->xmap.ed)
      e_free(sizeof(real)*map->xmap.Nx *map->xmap.Ny*map->xmap.Nz, map->xmap.ed) ;

    /* do not free rs_grid, as this is just a pointer to xplor rs_grid*/
  }
}
/** The following is slightly different from "writePdbInfo()" in pdbutil.c ... */

void write_pdb(XplorMap * map){
  int i;
  char buffer[1024];
  FILE * pdb;
  /*  printf("%s%9.3f%9.3f%9.3f%9.3f%9.3f%9.3f P 1           1\n","CRYST1",map->cell[0],map->cell[1],map->cell[2],map->cell[3],map->cell[4],map->cell[5]);*/
  /* this needs to be changed in the future*/
  sprintf(buffer,"grep CRYST1 %s > out.pdb",pdb_filename);
  system(buffer);
  sprintf(buffer,"grep SCALE %s >> out.pdb",pdb_filename);
  system(buffer);
  pdb = fopen("out.pdb","a");
  for(i = 0;i<Nin;i++){
    fprintf(pdb,"ATOM  %5d %4s%4s%6d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",atno[i],fullatid[i],groupid[i],groupno[i],pos[i].x,pos[i].y,pos[i].z,occup[i],Bj[i]);
  }
  fclose(pdb);
}

void update_atoms(Params * p,const gsl_vector *v,int i){
  int flag = 0;
  int j = 0;
  real old_pos[3];
  real old_Bj;
  real old_occup;
  old_pos[0] = pos[i].x;
  old_pos[1] = pos[i].y;
  old_pos[2] = pos[i].z;
  old_Bj = Bj[i];
  old_occup = occup[i];
  
  
  /* Only update the atom properties if the difference is significative */
  if(p->flags & REFINE_X){
    if(fabs(pos[i].x - gsl_vector_get(v, i*p->np+j)) > RSTEP/2){
      pos[i].x = gsl_vector_get(v, i*p->np+j);
      flag = 1;	
    }
    j++;
  }
  if(p->flags & REFINE_Y){
    if(fabs(pos[i].y - gsl_vector_get(v, i*p->np+j)) > RSTEP/2){
      pos[i].y = gsl_vector_get(v, i*p->np+j);
      flag = 1;
    }
    j++;
  }
  if(p->flags & REFINE_Z){
    if(fabs(pos[i].z - gsl_vector_get(v, i*p->np+j)) > RSTEP/2){
      pos[i].z = gsl_vector_get(v, i*p->np+j);
      flag = 1;
    }
    j++;
  }
  if(p->flags & REFINE_B){
    if(fabs(Bj[i] - gsl_vector_get(v, i*p->np+j)) > BSTEP/2){
      Bj[i] = gsl_vector_get(v, i*p->np+j);
      flag = 1;
    }
    j++;
  }
  if(p->flags & REFINE_O){
    if(fabs(occup[i] - gsl_vector_get(v, i*p->np+j)) > 0.01){
      occup[i] = MIN(1,gsl_vector_get(v, i*p->np+j));
      flag = 1;
    }
    j++;
  }
  if(flag){	  
    correct_pdb_ed(p,old_pos,old_Bj,old_occup,i);
  }
}




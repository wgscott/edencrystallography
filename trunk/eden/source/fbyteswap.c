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
*******************************************************************************/

#include <stdio.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
void byteswap();

#define ARRSIZE 512
#define BUFSIZE ARRSIZE*4

/*************************************************************************
The problem with moving binary data files between the DECs and Suns or SGI
is that the DEC is a lit4le-endian machine while most standard unix boxes
are big-endian. What that means is that when you write out a float on
a Sun, the 4 bytes are written out like this:
        A B C D 
but on the DECs, they are written out like this:
        D C B A
 
If you run "dd conv=swab", you get
        B A D C 
 
which is not what you want, either.
 
Greg Tomaschke wrote this little program that will read in 
a file of floats, and do the conversion
        D C B A  ->  A B C D 
 
Run it like this
        fbyteswap binfile1 binfile2
 
If binfile1 is a SUN-format file, then binfile2 will be a DEC-format file
and vice-versa.
 
This assumes that binfile1 was written using the "write" system call,
and contains only floats.
 
**************************************************************************/

main (int argc, char **argv) {

  int ifile, ofile;
  int i, idx;
  float x[ARRSIZE], y[ARRSIZE];

  int n_read, n_wrtn;
  
  if (argc != 3 ) {
    printf("Usage: byteswap [infile] [outfile]\n");
    return 0;
  }
  
  ifile = open(argv[1], O_RDONLY);
  ofile = open(argv[2], 
		  O_WRONLY | O_CREAT | O_TRUNC , 0600);

  idx = 0;
  n_read=read(ifile, x, BUFSIZE) ;
  if (n_read < 0) {
    perror("Error reading from input file");
  }
  while (n_read == BUFSIZE) {
    byteswap(x, ARRSIZE);
    n_wrtn=write(ofile, x, BUFSIZE);
    if (n_wrtn < 0) {
      perror("Error writing to output file");
    }
    
    n_read=read(ifile, x, BUFSIZE) ;
    if (n_read < 0) {
      perror("Error reading from input file");
    }
  }
  if (n_read > 0) {
    byteswap(x, n_read/4);
    n_wrtn=write(ofile, x, n_read);
    if (n_wrtn < 0) {
      perror("Error writing to output file");
    }
  }
  close(ofile);
  return 0;
}


void byteswap(char *a, int len) {

  int idx, n, i;
  char tmp;
  char *cptr;
 
  cptr = a;
  idx = 0;

  for (n=0; n<len; n++) {
    /*    for (i=0; i<4; i++) printf(" %3o",a[idx+i]); */
    for (i=0; i<2; i++) {
      tmp = a[idx+i];
      a[idx+i] = a[idx+3-i];
      a[idx+3-i] = tmp;
    }
/*    printf (" --> ");
    for (i=0; i<4; i++) printf(" %3o",a[idx+i]);
    printf("\n"); */
    idx += 4;
  }
  return;
}
     



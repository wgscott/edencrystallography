########################################################################
#
#			Makefile for Eden 
#
#	eden.c is a little driver program, directing control to one of 
#		the branches of Eden:
#
#	addmaps.c adds terms in two Np-space files 
#	apodfc.c does a Wilson-like calculation to apodize fcalc files        
#	apodfo.c does a Wilson-like calculation to apodize fobs files        
#	back.c transforms an hkl file onto a grid in physical space 
#	cadhkl.c combines terms in two hkl files 
#	convert.c reformats from .bin to .map or from .map to .bin 
#	count.c counts electrons in a .bin file using a pdb file info 
#	distance.c measures a metric distance among Np-space files 
#	dphase.c compares two hkl files 
#	expandfc.c expands an fcalc file to P1 symmetry 
#	expandfo.c expands an fobs file to P1 symmetry 
#	forth.c transforms a grid in physical space to an hkl file 
#	maketar.c creates targets and weights from a model  
#       perturbhkl.c perturbs an fcalc file 
#	ranphase.c creates an fcalc with random phases from an fobs 
#	refine.c does real-space ref't of a pdb file wrt a map file
#	regrid.c converts el/voxel data to sampled densities 
#	shapes.c determines the morphology of sampled densities 
#	solve.c is the main solver module with or w/o Mir, MAD 
#	sym.c uses symmetry & unit cell info to generate full pdb 
#	tohu.c converts .pdb info to hkl info (preprocessor) 
#	variance.c measures variance among el/voxelfiles. 
#
#	Other lower-level files:
#
#	apodutil.c contains functions used by both Apodfc and -fo 
#	basics.c handles input parameters for all Eden programs 
#	calc_rho.c converts pdb info to an xplor-style map.
#	ccg.c contains the least squares conjugate gradient solver 
#	ccgutil.c provides the interface between Back/Solve and ccg 
#	cost.c oversees the cost function calcs for Solve and Back   
#	concost.c does the constraint cost function calculations 
#	crysutil.c contains crystallographic computational functions 
#	fftwshare.c has functions shared among codes that call FFT's 
#	hklread.c handles read functions in (hkl) space 
#	hklinit.c handles initializations in (hkl) space 
#	hklutil.c handles write functions & utilities in (hkl) space 
#	highres.c handles high-resolution processing for Solve 
#	init.c handles initialization functions for Solve and Back 
#	pdbutil.c manipulates pdb info used by Sym and Tohu 
#       qshare.c handles Np-space functions 
#       remap.c handles functions for remapping in close vicinity 
#       regutil.c handles functions common to Regrid and Count  
#	ssetup.c contains the set-up functions for Solve 
#	util.c handles general utilities needed by all programs 
#	
#	fbyteswap.c deals with low- <-> high-endian binary switching.
#
#######################################################################
# 	Directory where the FFTW include and lib are
#######################################################################
FFTW		= /sw
LIB  		= $(FFTW)/lib
INC  		= $(FFTW)/include

#######################################################################
# 	Directory to install the binary, version
#######################################################################
BIN 		= $(HOME)/bin
VERSION		=V5.3

#######################################################################
# 	Compiler flags for optimization and/or debugging
#
# -O2 -ieee 	for no core dump
# -O2  -g 		for debug
# -O3 on HP, -O2 on SGI or DEC.
#
#######################################################################
CFLAGS		=   -Wall -pipe -g3
#CFLAGS		=   -Wall -pipe -marchi686 -pedantic  -gdwarf-2 -g3
#CFLAGS		=   -Wall -O3 -march=i686 
ifneq ($(FFTW),/usr/local)
CFLAGS		+= -I$(INC)
endif

#			NOTE:
# DOUBLESWITCH = ON	for double precision;
# DOUBLESWITCH = OFF	for single precision.
#
# Use different libraries and include files for single and
# double precision versions of EDEN.
#######################################################################
DOUBLESWITCH 	= ON

ifeq ($(DOUBLESWITCH),ON)
LDFLAGS	 = -L$(LIB) -ldfftw -lm  
CFLAGS	+= -DDOUBLE 
TARGET	 = deden
else
LDFLAGS	 = -L$(LIB) -lsfftw -lm  
TARGET	 = seden
endif

LDFLAGS += -g -lgsl -lgslcblas 

# GSL_PREFIX:= /sw/lib
# LDFLAGS += -l$(GSL_PREFIX)/libgsl.a -l$(GSL_PREFIX)/libgslcblas.a 

.c.o:
	$(CC) $(CFLAGS) -c $<

OBJ	= 	addmaps.o   apodfc.o   apodfo.o     back.o     cadhkl.o  \
		convert.o   count.o    distance.o   dphase.o  expandfc.o \
		expandfo.o  forth.o   maketar.o  perturbhkl.o ranphase.o \
		refine.o    regrid.o   shapes.o     solve.o    sym.o     \
                tohu.o      variance.o                                   \
		apodutil.o  basics.o   calc_rho.o  ccg.o       ccgutil.o \
                concosts.o  cost.o     crysutil.o  fftwshare.o highres.o \
                hklinit.o   hklread.o  hklutil.o   init.o      pdbutil.o \
                qshare.o    regutil.o  remap.o     ssetup.o    util.o     

all:		eden	fbyteswap 

eden:  		eden.o   $(OBJ)
		$(CC)  $^ $(LDFLAGS) -o $@

fbyteswap:	fbyteswap.c		
		$(CC) -o $@ $^

install:	all
		mv fbyteswap $(BIN)
		mv eden $(BIN)/$(TARGET)
		chmod 751 $(BIN)/$(TARGET)
		chmod +x mkeden.csh
		./mkeden.csh $(BIN) $(BIN)/$(TARGET)

dist:		clean
		cd ../..; tar -cvf eden_$(VERSION).tar --exclude CVS --exclude \#\*\# --exclude \*~ --exclude Misc --exclude .\* eden; gzip -f eden_$(VERSION).tar

clean:
		rm -f eden eden.o $(OBJ)

#######################################################################
# 	Dependencies for all the files, entered manually
#	Note: eden.h includes util.h, cellparams.h, symmetry.h, dims.h
#######################################################################
eden.o:  	eden.c util.h 
addmaps.o:	addmaps.c util.h cellparams.h dims.h
apodfc.o:	apodfc.c util.h cellparams.h dims.h apodize.h
apodfo.o:	apodfo.c util.h cellparams.h dims.h apodize.h
apodutil.o:	apodutil.c util.h cellparams.h dims.h apodize.h
back.o:		back.c ccg.h cost.h eden.h eden_Np.h 
basics.o:	basics.c cellparams.h dims.h symmetry.h util.h pdb.h
cadhkl.o:	cadhkl.c util.h cellparams.h dims.h symmetry.h
calc_rho.o:	calc_rho.c calc_rho.h util.h refine.h
convert.o:	convert.c util.h cellparams.h dims.h 
ccg.o:		ccg.c ccg.h util.h 
ccgutil.o:	ccgutil.c ccg.h dims.h eden_Np.h mir.h util.h
concosts.o:	concosts.c eden.h eden_Np.h cost.h remap.h
cost.o:		cost.c eden.h eden_Np.h cost.h remap.h 
count.o:	count.c util.h cellparams.h dims.h symmetry.h reg.h pdb.h 
crysutil.o:	crysutil.c util.h cellparams.h symmetry.h dims.h
distance.o:	distance.c dims.h util.h 
dphase.o:	dphase.c dims.h symmetry.h util.h cellparams.h
expandfc.o:	expandfc.c cellparams.h dims.h symmetry.h util.h 
expandfo.o:	expandfo.c dims.h symmetry.h util.h 
fftwshare.o:	fftwshare.c dims.h util.h
forth.o:	forth.c  eden.h 
highres.o:	highres.c  cellparams.h dims.h symmetry.h util.h
hklinit.o:	hklinit.c eden.h mir.h 
hklread.o:	hklread.c dims.h symmetry.h util.h 
hklutil.o:	hklutil.c dims.h symmetry.h util.h mir.h
init.o:		init.c eden.h mir.h ccg.h
maketar.o:	maketar.c cellparams.h dims.h util.h pdb.h
perturbhkl.o:	perturbhkl.c cellparams.h dims.h symmetry.h util.h
pdbutil.o:	pdbutil.c cellparams.h symmetry.h pdb.h util.h
qshare.o:	qshare.c eden.h eden_Np.h 
ranphase.o:	ranphase.c  util.h cellparams.h dims.h symmetry.h
refine.o:	refine.c refine.h eden.h pdb.h reg.h
regrid.o:	regrid.c cellparams.h dims.h symmetry.h util.h reg.h
regutil.o:	regutil.c cellparams.h dims.h symmetry.h util.h reg.h pdb.h
remap.o:	remap.c cellparams.h dims.h remap.h symmetry.h util.h
shapes.o:	shapes.c cellparams.h dims.h symmetry.h util.h reg.h
solve.o:	solve.c eden.h eden_Np.h cost.h mir.h 
ssetup.o:	ssetup.c eden.h eden_Np.h cost.h mir.h 
sym.o:  	sym.c cellparams.h dims.h symmetry.h pdb.h util.h
tohu.o:		tohu.c cellparams.h dims.h pdb.h symmetry.h util.h
util.o:		util.c util.h
variance.o:   	variance.c dims.h util.h

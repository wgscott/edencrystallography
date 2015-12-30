# EDEN: Electron DENsity holographic refinement #

# Authors #

Eden is written by Hanna and Abraham Szöke.

# Current Stable Version #

eden-5.3

# What is Eden? #
Eden is crystallographic real-space electron-density refinement and optimization program that produces electron density maps with minimal model bias in a robust manner.

This program improves macromolecular crystallographic electron density maps in a maximally unbiased manner.  There are now two versions of eden, "seden" and "deden." seden is compiled with single-precision fftw libs, and is faster and less memory-intensive. deden is double- precision and is possibly more accurate (but in practice the differences appear insignificant). ieden is a new python Tkinter-based GUI. Type "eden" or "seden" to invoke single-precision eden, "ieden" for the GUI, and "deden" for double-precision eden. GUI users need to define the environment variable EDITOR, eg: export EDITOR='see' or setenv EDITOR vim.

# PDF Manual #

[EdenUserManual.pdf](http://edencrystallography.googlecode.com/files/EdenUserManual.pdf)

# Quick Links #

The [Eden Home Page](http://www.edencrystallography.org/). (This page is currently off-line).

I have an [Eden Quick Start](http://chemistry.ucsc.edu/%7Ewgscott/eden/eden_quickstart.html) page. This gives an example of how to use Eden to make a minimally biased electron density map given a set of coordinates and a set of Fobs.

My [Eden on OS X](http://xanana.ucsc.edu/xtal/eden.html) page.


# How to Install Eden #

Eden should compile and run on any unix or linux platform (linux, sgi, osx, etc). It has no configure (yet) so you have to edit the makefile.  It needs gsl and fftw2 to compile, and prefers to have grace (xmgr), tk, and python present at runtime. It is primarily a command-line program, but now has a primitive tkinter gui.


## Linux Packages ##


### Debian Linux package (i386-ubuntu) ###

http://diablo.ucsc.edu/~wgscott/debian/deb/eden-5.3-2_i386.deb

### Gentoo Linux ###

To install Eden on Gentoo Linux, do this:

emerge eden


### RH/Fedora rpm ###
(converted from above with alien)

http://diablo.ucsc.edu/~wgscott/debian/rpm/eden-5.3-2.i386.rpm




## Mac OS X Package ##

### Using the Fink package manager on OS X ###

Eden is [available as a package in Fink](http://pdb.finkproject.org/pdb/package.php/eden) and by far the easiest thing to do is to issue the command

> fink install eden

and it takes about two minutes to compile, and everything is set up for you.  Issue

> fink describe eden

for additional details.

## Manual Install ##

This is OS X specific, but completely generalizable.

1. [l Download Eden from the Eden GoogleCode website](http://code.google.com/p/edencrystallography/downloads/list), and put the tar file into /usr/local/eden (or your favorite equivalent).

2. Unpack the tar file and read the README\_FIRST file. Install fftw2, gsl and grace if you don't have these already (which is why doing all this with fink makes more sense.)

3. Edit the Makefile to include the following:

> libs=-L/sw/lib -lfftw -lm
> opt=-O3
> CC=cc

(Use the -L/path(s) to wherever your fftw2 and gsl libraries are.)

4. Make other suggested edits to the Makefile, including:

> mv eden /usr/local/bin/eden

5. Issue the command

> export EDENHOME=/usr/local/eden

if you use a Bourne-like shell, or

> setenv EDENHOME /usr/local/eden

if you use tcsh.  Put the relevant command into the relevant startup file. Then to build the program, issue the commands

> cd $EDENHOME/source
> sudo make  (or run make after becoming root).

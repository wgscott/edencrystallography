#
# Makefile for EDEN manual
#
# In CVS we have the tex file and the pdf file (not gzipped,
# because then CVS can not keep track of the history and changes)
# Although the pdf can be generated in most cases we keep it in CVS as well
# in case a user does not have latex installed.
#
TEX	= pdflatex
IDX	= makeindex

all:		UserManual.pdf

UserManual.pdf:	UserManual.tex
		$(TEX) $^
		$(IDX) -s hfill.ist UserManual
		$(TEX) $^

clean:
	/bin/rm -r UserManual.pdf UserManual.aux UserManual.log UserManual.toc

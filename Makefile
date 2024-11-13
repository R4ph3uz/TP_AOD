# Makefile qui genere l'executable distanceEdition et fait des tests de verification
#
#
CC=gcc
LATEXC=pdflatex
DOCC=doxygen
CFLAGS=-g -Wall 
OPT=$(CFLAGS)

REFDIR=.
SRCDIR=$(REFDIR)/src
BINDIR=$(REFDIR)/bin
DOCDIR=$(REFDIR)/doc
TESTDIR=$(REFDIR)/tests
REPORTDIR=$(REFDIR)/report

LATEXSOURCE=$(wildcard $(REPORTDIR)/*.tex)
CSOURCE=$(wildcard $(SRCDIR)/*.c)
PDF=$(LATEXSOURCE:.tex=.pdf)

all: binary binary_debug
#report binary_perf

binary: $(BINDIR)/distanceEdition

binary_debug: $(BINDIR)/distanceEditiondebug

binary_perf: $(BINDIR)/distanceEdition-perf

$(BINDIR)/distanceEdition-perf: $(SRCDIR)/distanceEdition.c $(BINDIR)/aware.o
	$(CC) $(OPT) -D__PERF_MESURE__ -I$(SRCDIR) -o $(BINDIR)/distanceEdition-perf $(BINDIR)/aware.o $(SRCDIR)/distanceEdition.c 

report: $(PDF) 

doc: $(DOCDIR)/index.html


$(BINDIR)/distanceEdition: $(SRCDIR)/distanceEdition.c $(BINDIR)/aware.o $(SRCDIR)/characters_to_base.h
	$(CC) $(OPT) -I$(SRCDIR) -o $(BINDIR)/distanceEdition $(BINDIR)/aware.o $(SRCDIR)/distanceEdition.c 

$(BINDIR)/Needleman-Wunsch-recmemo.o: $(SRCDIR)/Needleman-Wunsch-recmemo.h $(SRCDIR)/Needleman-Wunsch-recmemo.c 
	$(CC) $(OPT) -I$(SRCDIR) -c  -o $(BINDIR)/Needleman-Wunsch-recmemo.o $(SRCDIR)/Needleman-Wunsch-recmemo.c

$(BINDIR)/aware.o: $(SRCDIR)/aware.h $(SRCDIR)/aware.c 
	$(CC) $(OPT) -I$(SRCDIR) -c  -o $(BINDIR)/aware.o $(SRCDIR)/aware.c

$(BINDIR)/extract-fasta-sequences-size: $(SRCDIR)/extract-fasta-sequences-size.c
	$(CC) $(OPT) -I$(SRCDIR) -o $(BINDIR)/extract-fasta-sequences-size $(SRCDIR)/extract-fasta-sequences-size.c

clean:
	rm -rf $(DOCDIR) $(BINDIR)/* $(REPORTDIR)/*.aux $(REPORTDIR)/*.log  $(REPORTDIR)/rapport.pdf 

#$(BINDIR)/distanceEdition: $(CSOURCE)
#	$(CC) $(CFLAGS)  $^ -o $@ 

$(BINDIR)/distanceEditiondebug: $(SRCDIR)/distanceEdition.c $(BINDIR)/aware.o $(SRCDIR)/characters_to_base.h
	$(CC) $(CFLAGS)  $^ -o $@ -DDEBUG

%.pdf: $(LATEXSOURCE)
	$(LATEXC) -output-directory $(REPORTDIR) $^ 

$(DOCDIR)/index.html: $(SRCDIR)/Doxyfile $(CSOURCE)
	$(DOCC) $(SRCDIR)/Doxyfile


test: $(BINDIR)/distanceEdition $(TESTDIR)/Makefile-test
	cd $(TESTDIR) ; make -f Makefile-test all 
	
test-valgrind: $(BINDIR)/distanceEdition $(TESTDIR)/Makefile-test
	make -f $(TESTDIR)/Makefile-test all-valgrind
	
.PHONY: all doc bin report 


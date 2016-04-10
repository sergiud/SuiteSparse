#===============================================================================
# Mongoose/Tcov/Makefile
#===============================================================================

# run statement-coverage test without Valgrind
default: $(METIS) go

# to run with Valgrind as well
valgrind: $(METIS) vgo

ccode: all

include ../../SuiteSparse_config/SuiteSparse_config.mk

#-------------------------------------------------------------------------------

LIBS = \
    ../../SuiteSparse_config/libsuitesparseconfig.a \
    -lrt

INC = -I./Include -I../../SuiteSparse_config

OBJ = \
  cs.o \
  BoundaryHeap.o \
  Coarsening.o \
  Conditioning.o \
  EdgeSeparator.o \
  FidduciaMattheyes.o \
  GuessCut.o \
  Interop.o \
  Matching.o \
  Nuri_qpGradProj.o \
  QPboundary.o \
  QPgradproj.o \
  QPlinks.o \
  QPmaxheap.o \
  QPminheap.o \
  QPnapdown.o \
  QPnapsack.o \
  QPnapsack2.o \
  QPnapup.o \
  Refinement.o \
  Waterdance.o

CF =  -pg \
	-O0 -g -fprofile-arcs -ftest-coverage \
	-Wall -W -Wshadow -Winline -Wno-unused-parameter \
	-Wredundant-decls -Wdisabled-optimization \
	-ansi -fexceptions

C = $(CXX) $(CF) $(INC)

all: test1

library: test1

purge: distclean

distclean: clean
	- $(RM) test1 *.mtx.txt pfile tfile cov.out *.dot *.hpp
	- $(RM) -r $(PURGE)

clean:
	- $(RM) -r $(CLEAN)

#------------------------------------------------------------------------------

copy_resources:
	cp ../Source/*.cpp ./
	cp ../Include/*.hpp ./

test1: test1.cpp $(OBJ)
	$(C) test1.cpp -o test1 $(OBJ) $(LIBS) -lm
	
go: copy_resources test1
	- rm cov.out	
	- ./test1 ~/PROBLEMS/139.mtx > test1_139.mtx.txt
	- ./test1 ~/PROBLEMS/515.mtx > test1_515.mtx.txt
	- ./cov

#------------------------------------------------------------------------------

cs.o: cs.cpp
	$(C) -c $(INC) $<
	
BoundaryHeap.o: BoundaryHeap.cpp
	$(C) -c $(INC) $<

Coarsening.o: Coarsening.cpp
	$(C) -c $(INC) $<

Conditioning.o: Conditioning.cpp
	$(C) -c $(INC) $<

EdgeSeparator.o: EdgeSeparator.cpp
	$(C) -c $(INC) $<

FidduciaMattheyes.o: FidduciaMattheyes.cpp
	$(C) -c $(INC) $<

GuessCut.o: GuessCut.cpp
	$(C) -c $(INC) $<

Interop.o: Interop.cpp
	$(C) -c $(INC) $<

Matching.o: Matching.cpp
	$(C) -c $(INC) $<

Nuri_qpGradProj.o: Nuri_qpGradProj.cpp
	$(C) -c $(INC) $<

QPboundary.o: QPboundary.cpp
	$(C) -c $(INC) $<

QPgradproj.o: QPgradproj.cpp
	$(C) -c $(INC) $<

QPlinks.o: QPlinks.cpp
	$(C) -c $(INC) $<

QPmaxheap.o: QPmaxheap.cpp
	$(C) -c $(INC) $<

QPminheap.o: QPminheap.cpp
	$(C) -c $(INC) $<

QPnapdown.o: QPnapdown.cpp
	$(C) -c $(INC) $<

QPnapsack.o: QPnapsack.cpp
	$(C) -c $(INC) $<

QPnapsack2.o: QPnapsack2.cpp
	$(C) -c $(INC) $<

QPnapup.o: QPnapup.cpp
	$(C) -c $(INC) $<

Refinement.o: Refinement.cpp
	$(C) -c $(INC) $<

Waterdance.o: Waterdance.cpp
	$(C) -c $(INC) $<

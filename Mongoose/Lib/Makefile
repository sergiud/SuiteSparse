#-------------------------------------------------------------------------------
# Mongoose/Lib/Makefile
#-------------------------------------------------------------------------------

default: all

include ../../SuiteSparse_config/SuiteSparse_config.mk

OBJS = \
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

CLIB = ../../SuiteSparse_config/libsuitesparseconfig.a

INC = -I../Include -I../../SuiteSparse_config

#-------------------------------------------------------------------------------

all: runtime

library: runtime

lib: runtime

distclean: purge

clean:
	$(RM) -r $(CLEAN)

purge: clean
	$(RM) -r $(PURGE)

runtime: $(OBJS)
	$(AR) -r "libMongoose.a" $(CLIB) $(OBJS)

#-------------------------------------------------------------------------------

cs.o: ../Source/cs.cpp
	$(CXX) -c $(INC) $<
	
BoundaryHeap.o: ../Source/BoundaryHeap.cpp
	$(CXX) -c $(INC) $<

Coarsening.o: ../Source/Coarsening.cpp
	$(CXX) -c $(INC) $<

Conditioning.o: ../Source/Conditioning.cpp
	$(CXX) -c $(INC) $<

EdgeSeparator.o: ../Source/EdgeSeparator.cpp
	$(CXX) -c $(INC) $<

FidduciaMattheyes.o: ../Source/FidduciaMattheyes.cpp
	$(CXX) -c $(INC) $<

GuessCut.o: ../Source/GuessCut.cpp
	$(CXX) -c $(INC) $<

Interop.o: ../Source/Interop.cpp
	$(CXX) -c $(INC) $<

Matching.o: ../Source/Matching.cpp
	$(CXX) -c $(INC) $<

Nuri_qpGradProj.o: ../Source/Nuri_qpGradProj.cpp
	$(CXX) -c $(INC) $<

QPboundary.o: ../Source/QPboundary.cpp
	$(CXX) -c $(INC) $<

QPgradproj.o: ../Source/QPgradproj.cpp
	$(CXX) -c $(INC) $<

QPlinks.o: ../Source/QPlinks.cpp
	$(CXX) -c $(INC) $<

QPmaxheap.o: ../Source/QPmaxheap.cpp
	$(CXX) -c $(INC) $<

QPminheap.o: ../Source/QPminheap.cpp
	$(CXX) -c $(INC) $<

QPnapdown.o: ../Source/QPnapdown.cpp
	$(CXX) -c $(INC) $<

QPnapsack.o: ../Source/QPnapsack.cpp
	$(CXX) -c $(INC) $<

QPnapsack2.o: ../Source/QPnapsack2.cpp
	$(CXX) -c $(INC) $<

QPnapup.o: ../Source/QPnapup.cpp
	$(CXX) -c $(INC) $<

Refinement.o: ../Source/Refinement.cpp
	$(CXX) -c $(INC) $<

Waterdance.o: ../Source/Waterdance.cpp
	$(CXX) -c $(INC) $<

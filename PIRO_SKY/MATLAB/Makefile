#===============================================================================
# PIRO_SKY/MATLAB/Makefile: for compiling the MATLAB interfaces for the sparse
# SVD library. 
#
# TODO: delete this file.  use a *_make.m file exclusively.
#===============================================================================

default: all

include ../../SuiteSparse_config/SuiteSparse_config.mk

#MATLAB_FLAGS = -DNDEBUG -DMATLAB_MEX_FILE
#C = $(CC) $(CFLAGS) $(PERF) $(MATLAB_FLAGS) -DLONG -DMATLAB_MEX_FILE
C = $(CC) $(CFLAGS) $(PERF) $(MATLAB_FLAGS) -DLONG 
MX = $(MEX) $(PERF) $(MXFLAGS) -DLONG

purge: distclean

distclean: clean
	- $(RM) *.mex* *.dll *.o *.obj

clean:
	- $(RM) *.mex* *.dll *.o *.obj

I = -I../Include -I../../SuiteSparse_config \
    -I/usr/local/matlab/R2008b/extern/include/ -I../../PIRO_BAND/Include

MEXLIB = blksky blksky_symbolic

OBJ = blksky.o blksky_symbolic.o blksky_utils.o genband_givens.o

all: $(OBJ) $(MEXLIB)

blksky_symbolic:
	$(MX) $(I) blksky_symbolic_mex.c blksky_symbolic.o \
	blksky_utils.o -output blksky_symbolic

blksky:
	$(MX) $(I) blksky_mex.c genband_givens.o \
	blksky_symbolic.o blksky_utils.o blksky.o \
	../../PIRO_BAND/Source/libpiro_band.a -output blksky

blksky.o: ../Source/blksky.c
	$(C) -c $(I) $< -o $@

blksky_symbolic.o: ../Source/blksky_symbolic.c
	$(C) -c $(I) $< -o $@

blksky_utils.o: ../Source/blksky_utils.c
	$(C) -c $(I) $< -o $@

genband_givens.o: ../Source/genband_givens.c
	$(C) -c $(I) $< -o $@


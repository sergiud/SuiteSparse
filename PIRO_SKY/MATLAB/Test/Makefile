#===============================================================================
# PIRO_SKY/MATLAB/Makefile: for compiling the MATLAB interfaces for the sparse
# SVD library. 
#
# TODO: delete this and use a *.m file instead
#===============================================================================

default: all

include ../../../SuiteSparse_config/SuiteSparse_config.mk

CFLAGS = -g -DLONG
C = $(CC) $(CFLAGS)

purge: distclean

distclean: clean
	- $(RM) .dll *.o *.obj ./test_blksky

clean:
	- $(RM) *.dll *.o *.obj 

I = -I../Include -I../../SuiteSparse_config -I../../CXSparse/Include

LIB = ../../../CXSparse/Lib/libcsparse.a ../Source/libpiro_sky.a -lm ../../PIRO_BAND/Source/libpiro_band.a

all: test_blksky

test_blksky: ../../Tcov/test_blksky.c
	$(C) $(I) $< -o $@ $(LIB)


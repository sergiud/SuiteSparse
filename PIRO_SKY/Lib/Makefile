#===============================================================================
# SkylineSVD/Lib/Makefile: for compiling the C library.
#===============================================================================

default: all

ccode: all

include ../../SuiteSparse_config/SuiteSparse_config.mk

C = $(CC) $(CFLAGS) $(PIRO_BAND_REDUCE_CONFIG) -DLONG
#CFLAGS = -g
#C = $(CC) $(CFLAGS)

all: libpiro_sky.a

library: libpiro_sky.a

purge: distclean

distclean: clean
	- $(RM) -r $(PURGE)

clean:
	- $(RM) -r $(CLEAN)

DI = blksky.o blksky_symbolic.o blksky_utils.o genband_givens.o


OBJ = $(DI)

libpiro_sky.a: $(OBJ)
	$(ARCHIVE) libpiro_sky.a $(OBJ)
	$(RANLIB) libpiro_sky.a

I = -I../Include -I../../SuiteSparse_config -I../../PIRO_BAND/Include/

#-------------------------------------------------------------------------------
blksky.o: ../Source/blksky.c
	$(C) $(PERF) -c $(I) $< -o $@

blksky_utils.o: ../Source/blksky_utils.c
	$(C) $(PERF) -c $(I) $< -o $@

genband_givens.o: ../Source/genband_givens.c
	$(C) $(PERF) -c $(I) $< -o $@

blksky_symbolic.o: ../Source/blksky_symbolic.c
	$(C) $(PERF) -c $(I) $< -o $@

#-------------------------------------------------------------------------------

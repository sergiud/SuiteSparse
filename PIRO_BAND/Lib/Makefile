#-------------------------------------------------------------------------------
# PIRO_BAND/Lib/Makefile: for compiling the band reduction library
#-------------------------------------------------------------------------------

default: all

ccode: all

include ../../SuiteSparse_config/SuiteSparse_config.mk

#C = $(CC) $(CFLAGS) $(PIRO_BAND_REDUCE_CONFIG)
C = $(CC) $(CF)

all: libpiroband.a

library: libpiroband.a

purge: distclean

distclean: clean
	- $(RM) -r $(PURGE)

clean:
	- $(RM) -r $(CLEAN)

PERF =
# BENCHMARK is for internal use. 
#PERF =  -DBENCHMARK

INC =   ../Include/piro_band.h \
        ../Include/piro_band_lapack.h \
        ../Source/piro_band_blocksize.h \
        ../Source/piro_band_memory.h \
        ../Source/piro_band_util.h \
        ../Source/piro_band_uv_update.h \
        ../Source/piro_band_internal.h \
        ../Source/piro_band_lapack_internal.h

#-------------------------------------------------------------------------------
# piro_band_reduce (double, int, real)
#-------------------------------------------------------------------------------

DI = piro_band_main.o piro_band_lapack_main.o piro_band_blocksize_main.o \
	piro_band_uv_update_main.o

#-------------------------------------------------------------------------------

# to compile just the double/int/real version, use OBJ = $(DI) $(BS_I)
OBJ = $(DI) 

libpiroband.a: $(OBJ)
	$(ARCHIVE)  libpiroband.a $(OBJ)
	- $(RANLIB) libpiroband.a

$(OBJ): $(INC)

I = -I../Include -I../../SuiteSparse_config -I../Source

#-------------------------------------------------------------------------------
piro_band_main.o: ../Source/piro_band_main.c $(INC) ../Source/piro_band.c \
    ../Source/piro_band_givens.c
	$(C) $(PERF) -c $(I) $< -o $@

#-------------------------------------------------------------------------------
piro_band_lapack_main.o: ../Source/piro_band_lapack_main.c $(INC) \
    ../Source/piro_band_lapack.c
	$(C) $(PERF) -c $(I) $< -o $@

#-------------------------------------------------------------------------------
piro_band_blocksize_main.o: ../Source/piro_band_blocksize_main.c $(INC) \
    ../Source/piro_band_blocksize.c
	$(C) $(PERF) -c $(I) $< -o $@

#-------------------------------------------------------------------------------
piro_band_uv_update_main.o: ../Source/piro_band_uv_update_main.c $(INC) \
    ../Source/piro_band_uv_update.c
	$(C) $(PERF) -c $(I) $< -o $@

#-------------------------------------------------------------------------------

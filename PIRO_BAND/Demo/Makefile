#-------------------------------------------------------------------------------
# PIRO_BAND/Demo/Makefile
#-------------------------------------------------------------------------------

default: all

include ../../SuiteSparse_config/SuiteSparse_config.mk

CC = gcc 

C = $(CC) $(CF)

all: run

I = -I../Source -I../Include -I../../SuiteSparse_config/ -I. 

LIB2 = ../Lib/libpiroband.a ../../SuiteSparse_config/libsuitesparseconfig.a \
        $(LDLIBS)

run: test_piro_band
	./test_piro_band

INC =   ../Include/piro_band.h \
        ../Include/piro_band_lapack.h \
        ../Source/piro_band_blocksize.h \
        ../Source/piro_band_internal.h \
        ../Source/piro_band_lapack_internal.h \
        ../Source/piro_band_memory.h \
        ../Source/piro_band_util.h \
        ../Source/piro_band_uv_update.h \
        ../Demo/test_piro_band.h \
        ../Demo/test_rand.h \
        ../../SuiteSparse_config/SuiteSparse_config.h

test_piro_band: test_piro_band.c test_piro_band_main.c test_rand.c Makefile \
    test_piro_band.c piro_band_testutils.c  $(INC) library
	$(C) $(I) test_rand.c test_piro_band_main.c -o test_piro_band $(LIB2)

library:
	( cd ../../SuiteSparse_config ; $(MAKE) library )
	( cd ../Lib ; $(MAKE) )

purge: distclean

distclean: clean
	- $(RM) -r test_piro_band *.out $(PURGE)

clean:
	- $(RM) -r $(CLEAN)

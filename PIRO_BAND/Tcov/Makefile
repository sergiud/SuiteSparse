#-------------------------------------------------------------------------------
# PIRO_BAND/Tcov/Makefile
#-------------------------------------------------------------------------------

default: all

include ../../SuiteSparse_config/SuiteSparse_config.mk

CC = gcc 

I = -I../Source -I../Include -I../../SuiteSparse_config/ -I. -I../Demo/

CF = -Wall -Wshadow -Wmissing-prototypes -Wstrict-prototypes \
    -Wredundant-decls -Wnested-externs -Wdisabled-optimization \
    -g -O0 -ftest-coverage -fprofile-arcs -ansi

#for profiling
#CF = -g -ftest-coverage -fprofile-arcs

C = $(CC) $(CF)

# PRETTY = grep -v "^\#" | indent -bl -nce -bli0 -i4 -sob 
# PRETTY = grep -v "^\#"
# PRETTY = cat

V = valgrind

GEN_SRC = cov_piro_band_main.c cov_piro_band_lapack_main.c \
	cov_piro_band_blocksize_main.c cov_piro_band_uv_update_main.c

SRC = piro_band_cover_main.c memory.c ../Demo/test_rand.c \
    ../../SuiteSparse_config/SuiteSparse_config.c

DEP_SRC = \
    ../Demo/piro_band_testutils.c \
    ../Source/piro_band.c \
    ../Source/piro_band_blocksize.c \
    ../Source/piro_band_blocksize_main.c \
    ../Source/piro_band_givens.c \
    ../Source/piro_band_lapack.c \
    ../Source/piro_band_lapack_main.c \
    ../Source/piro_band_main.c \
    ../Source/piro_band_uv_update.c \
    ../Source/piro_band_uv_update_main.c \
    ../Tcov/piro_band_cover.c \
    ../Demo/piro_band_testutils.c

INC = \
    ../Include/piro_band.h \
    ../Include/piro_band_lapack.h \
    ../Source/piro_band_blocksize.h \
    ../Source/piro_band_internal.h \
    ../Source/piro_band_lapack_internal.h \
    ../Source/piro_band_memory.h \
    ../Source/piro_band_util.h \
    ../Source/piro_band_uv_update.h \
    ../Tcov/mem.h \
    ../Tcov/piro_band_cover.h \
    ../Demo/test_piro_band.h \
    ../Demo/test_rand.h

all: run

cov: all

run: cov_piro_band
	./cov_piro_band
	./covall

memcheck: all
	$(V) ./cov_piro_band
	./covall

cov_piro_band: piro_band_main piro_band_lapack_main piro_band_blocksize_main \
        piro_band_uv_update_main piro_band_cover_main.c memory.c $(DEP_SRC) \
        Makefile $(SRC) $(INC)
	$(C) -g $(I) $(GEN_SRC) $(SRC) -o cov_piro_band -lm

# was
#       $(C) -g $(I) -DMALLOC=my_malloc2 -DFREE=my_free2 $(GEN_SRC) $(SRC) \
#           -o cov_piro_band -lm

purge: distclean

distclean: clean
	- $(RM) cov_piro_band *.out tmp/*.output *.sort cov_*.c
	- $(RM) -r $(PURGE)

clean:
	- $(RM) -r $(CLEAN)


.c.o:
	$(C) -c $(I) $*.c

#-------------------------------------------------------------------------------

piro_band_main: $(DEP_SRC) $(SRC) Makefile
	$(C) -DTCOV -DMALLOC=my_malloc2 -DFREE=my_free2 -E $(I) \
            ../Source/piro_band_main.c \
            | $(PRETTY) > cov_piro_band_main.c

piro_band_lapack_main: $(DEP_SRC) $(SRC) Makefile
	$(C) -DTCOV -DMALLOC=my_malloc2 -DFREE=my_free2 -E $(I) \
            ../Source/piro_band_lapack_main.c \
            | $(PRETTY) > cov_piro_band_lapack_main.c

piro_band_blocksize_main: $(DEP_SRC) $(SRC) Makefile
	$(C) -DTCOV -DMALLOC=my_malloc2 -DFREE=my_free2 -E $(I) \
            ../Source/piro_band_blocksize_main.c \
            | $(PRETTY) > cov_piro_band_blocksize_main.c

piro_band_uv_update_main: $(DEP_SRC) $(SRC) Makefile
	$(C) -DTCOV -DMALLOC=my_malloc2 -DFREE=my_free2 -E $(I) \
            ../Source/piro_band_uv_update_main.c \
            | $(PRETTY) > cov_piro_band_uv_update_main.c

#-------------------------------------------------------------------------------

#===============================================================================
# PIRO_SKY/Tcov/Makefile
#===============================================================================

default: run

include ../../SuiteSparse_config/SuiteSparse_config.mk

# select one of the following (valgrind is for Linux only)
  V = valgrind --leak-check=full
# V = valgrind
# V = 

CC = gcc
CFLAGS = -O -g -fprofile-arcs -ftest-coverage \
    	-Wall -W -Wshadow -Wmissing-prototypes -Wstrict-prototypes \
    	-Wredundant-decls -Wnested-externs -Wdisabled-optimization -ansi \
	-Wno-unused-parameter -Werror
C = $(CC) $(CFLAGS)

purge: distclean

I2 = -I../Include -I../../SuiteSparse_config -I../../CXSparse/Include -I../../PIRO_BAND/Include

LIB = -lm

# ../../CSparse/Lib/libcsparse.a ../Source/libpiro_sky.a -lm ../../PIRO_BAND/Source/libpiro_band.a

#===============================================================================# CXSparse
#===============================================================================

I1 = -I../../SuiteSparse_config -I../../CXSparse/Include

# AR = ar cr
# RANLIB = ranlib

CS_SOURCE =  cs_dupl.c cs_entry.c cs_util.c \
        cs_cumsum.c cs_norm.c cs_malloc.c \
	cs_transpose.c cs_compress.c \
	cs_print.c cs_load.c

CS_DL_OBJ =  cs_dupl_dl.o cs_entry_dl.o cs_util_dl.o \
        cs_cumsum_dl.o cs_norm_dl.o cs_malloc_dl.o \
	cs_transpose_dl.o cs_compress_dl.o \
	cs_print_dl.o cs_load_dl.o

# cs_%.c:
# 	- ln -s ../../CXSparse/Source/cs_$*.c

CS = $(CS_DL_OBJ)

$(CS): ../../CXSparse/Include/cs.h Makefile

cs_%_dl.o : ../../CXSparse/Source/cs_%.c
	$(C) $(I1) -DCS_LONG -c $< -o $@

# cs_%_dl.o : cs_%.c
# 	$(C) $(I1) -DCS_LONG -c $< -o $@

# .PRECIOUS: cs_%.c

#===============================================================================


INC =  \
    ../Include/blksky.h \
    ../Include/genband_internal.h \
    ../Include/genband_util.h \
    ../../PIRO_BAND/Include/piro_band.h \
    ../../PIRO_BAND/Include/piro_band_blocksize.h \
    ../../PIRO_BAND/Include/piro_band_internal.h \
    ../../PIRO_BAND/Include/piro_band_lapack.h \
    ../../PIRO_BAND/Include/piro_band_lapack_internal.h \
    ../../PIRO_BAND/Include/piro_band_memory.h \
    ../../PIRO_BAND/Include/piro_band_util.h \
    ../../PIRO_BAND/Include/piro_band_uv_update.h \
    ../../PIRO_BAND/Source/piro_band.c \
    ../../PIRO_BAND/Source/piro_band_blocksize.c \
    ../../PIRO_BAND/Source/piro_band_givens.c \
    ../../PIRO_BAND/Source/piro_band_lapack.c \
    ../../PIRO_BAND/Source/piro_band_uv_update.c

#===============================================================================
# PIRO_SKY
#===============================================================================

blksky.o: ../Source/blksky.c $(INC)
	$(C) $(I2) -DLONG -c $< -o $@

blksky_symbolic.o: ../Source/blksky_symbolic.c $(INC)
	$(C) $(I2) -DLONG -c $< -o $@

blksky_utils.o: ../Source/blksky_utils.c $(INC)
	$(C) $(I2) -DLONG -c $< -o $@

genband_givens.o: ../Source/genband_givens.c $(INC)
	$(C) $(I2) -DLONG -c $< -o $@

PIROSKY = blksky.o blksky_symbolic.o blksky_utils.o genband_givens.o


#===============================================================================
# PIRO_BAND
#===============================================================================

PIROBAND = piro_band_blocksize_main.o piro_band_lapack_main.o \
    piro_band_main.o piro_band_uv_update_main.o

piro_%.o : ../../PIRO_BAND/Source/piro_%.c $(INC)
	$(C) $(I2) -c $< -o $@

#===============================================================================

test_blksky.o: test_blksky.c $(INC)
	$(C) $(I2) -DLONG -c $< -o $@

all: $(CS) $(PIROSKY) $(PIROBAND) test_blksky

run: all
	- $(V) ./test_blksky < Matrices/R1.txt
	- $(V) ./test_blksky < Matrices/R1b.txt
	- $(V) ./test_blksky < Matrices/R3.txt
	- $(V) ./test_blksky < Matrices/Rmat
	- $(V) ./test_blksky < Matrices/Tina_DisCal_R.txt
	- $(V) ./test_blksky < Matrices/Tina_DisCal_R2.txt
	./covall

test_blksky: test_blksky.o  $(CS) $(PIROSKY) $(PIROBAND)
	$(C) $(I2) $< -o $@ $(CS) $(PIROSKY) $(PIROBAND) $(LIB)

distclean: clean
	- $(RM) -r $(PURGE) ./test_blksky

clean:
	- $(RM) -r $(CLEAN)


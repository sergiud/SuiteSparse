#
#  This file is part of MUMPS 5.7.3, released
#  on Mon Jul 15 11:44:21 UTC 2024
#
topdir = ..
libdir = $(topdir)/lib

default: d

.PHONY: default all s d c z multi clean
.SECONDEXPANSION:

all:	c z s d multi

c:	csimpletest csimpletest_save_restore
z:	zsimpletest zsimpletest_save_restore
s:	ssimpletest ssimpletest_save_restore
d:	c_mumps_paru dsimpletest dsimpletest_save_restore c_example_save_restore c_example
multi:	multiple_arithmetics_example


include $(topdir)/Makefile.inc



LIBSMUMPS = -L$(libdir) -lsmumps$(PLAT) -lmumps_common$(PLAT)
LIBDMUMPS = -L$(libdir) -ldmumps$(PLAT) -lmumps_common$(PLAT)
LIBCMUMPS = -L$(libdir) -lcmumps$(PLAT) -lmumps_common$(PLAT)
LIBZMUMPS = -L$(libdir) -lzmumps$(PLAT) -lmumps_common$(PLAT)

c_mumps_paru: $$@.o
	$(FC) $(OPTF) -I. -I$(topdir)/include -I$(topdir)/src $(INCS) -c  mmio.c mmio.o
	$(FL) -o $@ $(OPTL) mmio.o c_mumps_paru.o $(LIBDMUMPS) $(LORDERINGS) $(LIBS) $(LIBBLAS) $(LIBOTHERS)

ssimpletest:   $$@.o
	$(FL) -o $@ $(OPTL) ssimpletest.o  $(LIBSMUMPS) $(LORDERINGS) $(LIBS) $(RPATH_OPT) $(LIBBLAS) $(LIBOTHERS)

dsimpletest:  $$@.o 
	$(FL) -o $@ $(OPTL) dsimpletest.o  $(LIBDMUMPS) $(LORDERINGS) $(LIBS) $(RPATH_OPT) $(LIBBLAS) $(LIBOTHERS)

csimpletest:  $$@.o
	$(FL) -o $@ $(OPTL) csimpletest.o  $(LIBCMUMPS) $(LORDERINGS) $(LIBS) $(RPATH_OPT) $(LIBBLAS) $(LIBOTHERS)

zsimpletest:  $$@.o
	$(FL) -o $@ $(OPTL) zsimpletest.o  $(LIBZMUMPS) $(LORDERINGS) $(LIBS) $(RPATH_OPT) $(LIBBLAS) $(LIBOTHERS)

c_example: $$@.o
	$(FL) -o $@ $(OPTL) $@.o $(LIBDMUMPS) $(LORDERINGS) $(LIBS) $(RPATH_OPT) $(LIBBLAS) $(LIBOTHERS)


multiple_arithmetics_example: $$@.o
	$(FL) -o $@ $(OPTL) $@.o $(LIBSMUMPS) $(LIBDMUMPS) $(LIBCMUMPS) $(LIBZMUMPS) $(LORDERINGS) $(LIBS) $(RPATH_OPT) $(LIBBLAS) $(LIBOTHERS)

ssimpletest_save_restore:   $$@.o
	$(FL) -o $@ $(OPTL) ssimpletest_save_restore.o  $(LIBSMUMPS) $(LORDERINGS) $(LIBS) $(LIBBLAS) $(LIBOTHERS)

dsimpletest_save_restore:  $$@.o 
	$(FL) -o $@ $(OPTL) dsimpletest_save_restore.o  $(LIBDMUMPS) $(LORDERINGS) $(LIBS) $(LIBBLAS) $(LIBOTHERS)

csimpletest_save_restore:  $$@.o
	$(FL) -o $@ $(OPTL) csimpletest_save_restore.o  $(LIBCMUMPS) $(LORDERINGS) $(LIBS) $(LIBBLAS) $(LIBOTHERS)

zsimpletest_save_restore:  $$@.o
	$(FL) -o $@ $(OPTL) zsimpletest_save_restore.o  $(LIBZMUMPS) $(LORDERINGS) $(LIBS) $(LIBBLAS) $(LIBOTHERS)

c_example_save_restore: $$@.o
	$(FL) -o $@ $(OPTL) $@.o $(LIBDMUMPS) $(LORDERINGS) $(LIBS) $(LIBBLAS) $(LIBOTHERS)

.SUFFIXES: .c .F .o
.F.o:
	$(FC) $(OPTF) -I. -I$(topdir)/include -I$(topdir)/src $(INCS) -c $*.F $(OUTF)$*.o
.c.o:
	$(CC) $(OPTC) $(CDEFS) -I. -I$(topdir)/include -I$(topdir)/src $(INCS) -c $*.c $(OUTC)$*.o

clean:
	$(RM) *.o c_mumps_paru [sdcz]simpletest c_example multiple_arithmetics_example ssimpletest_save_restore dsimpletest_save_restore csimpletest_save_restore zsimpletest_save_restore c_example_save_restore

#-------------------------------------------------------------------------------
# Makefile for PIRO_BAND.
#-------------------------------------------------------------------------------

SUITESPARSE ?= $(realpath $(CURDIR)/..)
export SUITESPARSE

VERSION = 1.0.0

default: all

include ../SuiteSparse_config/SuiteSparse_config.mk

# Compile all C library, but not the MATLAB interface, and run the test
all:
	( cd Demo ; $(MAKE) )

# Exhaustive statement coverage test
cov:
	( cd Tcov ; $(MAKE) )

# Compile the C library
library:
	( cd Lib ; $(MAKE) )

# remove object files, but keep library *.a file and MATLAB mexFunctions
clean:
	( cd Lib ; $(MAKE) clean )
	$(RM) MATLAB/*.o
	( cd Demo ; $(MAKE) clean )
	( cd Tcov ; $(MAKE) clean )

# remove all except original files in the distribution
distclean:
	( cd Lib ; $(MAKE) distclean )
	$(RM) MATLAB/*.mex* MATLAB/*.o MATLAB/Test/*.mex* MATLAB/Test/*.o
	( cd Demo ; $(MAKE) distclean )
	( cd Tcov ; $(MAKE) distclean )

purge: distclean

# install PIRO_BAND
install:
	$(CP) Lib/libpiroband.a $(INSTALL_LIB)/libpiroband.$(VERSION).a
	( cd $(INSTALL_LIB) ; ln -sf libpiroband.$(VERSION).a libpiroband.a )
	$(CP) Include/piro_band*.h $(INSTALL_INCLUDE)
	chmod 744 $(INSTALL_LIB)/libpiroband*.a
	chmod 644 $(INSTALL_INCLUDE)/piro_band*.h

# uninstall PIRO_BAND
uninstall:
	$(RM) $(INSTALL_LIB)/libpiroband*.a
	$(RM) $(INSTALL_INCLUDE)/piro_band*.h


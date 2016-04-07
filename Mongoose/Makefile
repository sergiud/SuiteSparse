#-------------------------------------------------------------------------------
# Mongoose Makefile
#-------------------------------------------------------------------------------

SUITESPARSE ?= $(realpath $(CURDIR)/..)
export SUITESPARSE

VERSION = 1.0.0

default: all

include ../SuiteSparse_config/SuiteSparse_config.mk

# Compile the C-callable libraries only (there is no Demo for this package)
all: library demo

demo: library
	( cd Demo ; $(MAKE) )

# Compile the C-callable libraries only.
library:
	( cd Lib ; $(MAKE) )

# Remove all files not in the original distribution
purge:
	( cd Lib ; $(MAKE) purge )
	( cd Demo ; $(MAKE) purge )
	( cd Tcov ; $(MAKE) purge )
	( cd MATLAB ; $(RM) *.mex* *.o *.obj )

# Remove all files not in the original distribution, except keep the 
# compiled libraries.
clean:
	( cd Lib ; $(MAKE) clean )
	( cd Demo ; $(MAKE) clean )
	( cd Tcov ; $(MAKE) clean )

distclean: purge

ccode: all

# install Mongoose
install:
	$(CP) Lib/libMongoose.a $(INSTALL_LIB)/libMongoose.$(VERSION).a
	( cd $(INSTALL_LIB) ; ln -sf libMongoose.$(VERSION).a libMongoose.a )

# cs.hpp Graph.hpp Interop.hpp Mongoose.hpp Options.hpp

# uninstall Mongoose
uninstall:
	$(RM) $(INSTALL_LIB)/libMongoose*.a

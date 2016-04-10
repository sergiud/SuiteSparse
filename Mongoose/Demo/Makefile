#===============================================================================
# Mongoose/Demo/Makefile
#===============================================================================

default: all

include ../../SuiteSparse_config/SuiteSparse_config.mk

#-------------------------------------------------------------------------------

LIBS = \
    ../Lib/libMongoose.a \
    ../../SuiteSparse_config/libsuitesparseconfig.a \
    -lrt

INC = -I../Include -I../../SuiteSparse_config

C = $(CXX) $(CF) $(INC)

#------------------------------------------------------------------------------

all: demo

demo: $(LIBS)
	$(C) Demo.cpp -o Demo $(LIBS)

purge: distclean

distclean: clean
	- $(RM) -r Demo $(PURGE)

clean:
	- $(RM) -r $(CLEAN)

#-------------------------------------------------------------------------------
# PIRO_BAND/Doc/Makefile
#-------------------------------------------------------------------------------

default: all

include ../../SuiteSparse_config/SuiteSparse_config.mk

all: UserGuide.pdf

I = \
	../Include/piro_band.h \
	../Include/piro_band_lapack.h

C = \
	../Source/piro_band.c

M = \
	../MATLAB/piro_band.m \
	../MATLAB/piro_band_svd.m

UserGuide.pdf: UserGuide.tex $(I) $(C) $(M) Makefile rule.awk header.tex footer.tex getmproto mfooter.tex mheader.tex mfile.awk
	./getmproto ../MATLAB/piro_band.m > _piro_band_m.tex
	./getmproto ../MATLAB/piro_band_svd.m > _piro_band_svd_m.tex
	./getproto '/int piro_band\_reduce\_dri/, /\) ;/' ../Include/piro_band.h > _piro_band_reduce_1.tex
	./getproto '/int piro_band\_reduce\_drl/, /\) ;/' ../Include/piro_band.h > _piro_band_reduce_2.tex
	./getproto '/int piro_band\_reduce\_sri/, /\) ;/' ../Include/piro_band.h > _piro_band_reduce_3.tex
	./getproto '/int piro_band\_reduce\_srl/, /\) ;/' ../Include/piro_band.h > _piro_band_reduce_4.tex
	./getproto '/int piro_band\_reduce\_dci/, /\) ;/' ../Include/piro_band.h > _piro_band_reduce_5.tex
	./getproto '/int piro_band\_reduce\_dcl/, /\) ;/' ../Include/piro_band.h > _piro_band_reduce_6.tex
	./getproto '/int piro_band\_reduce\_sci/, /\) ;/' ../Include/piro_band.h > _piro_band_reduce_7.tex
	./getproto '/int piro_band\_reduce\_scl/, /\) ;/' ../Include/piro_band.h > _piro_band_reduce_8.tex
	./getproto '/void piro_band\_get\_blocksize /, /\) ;/' ../Include/piro_band.h > _piro_band_blk_1.tex
	./getproto '/void piro_band\_get\_blocksize_l/, /\) ;/' ../Include/piro_band.h > _piro_band_blk_2.tex
	./getproto '/void piro_band\_dgbbrd /, /\) ;/' ../Include/piro_band_lapack.h > _piro_band_lapack_1.tex
	./getproto '/void piro_band\_dgbbrd\_l/, /\) ;/' ../Include/piro_band_lapack.h > _piro_band_lapack_2.tex
	./getproto '/void piro_band\_zgbbrd /, /\) ;/' ../Include/piro_band_lapack.h > _piro_band_lapack_3.tex
	./getproto '/void piro_band\_zgbbrd\_l/, /\) ;/' ../Include/piro_band_lapack.h > _piro_band_lapack_4.tex
	./getproto '/void piro_band\_sgbbrd /, /\) ;/' ../Include/piro_band_lapack.h > _piro_band_lapack_5.tex
	./getproto '/void piro_band\_sgbbrd\_l/, /\) ;/' ../Include/piro_band_lapack.h > _piro_band_lapack_6.tex
	./getproto '/void piro_band\_cgbbrd /, /\) ;/' ../Include/piro_band_lapack.h > _piro_band_lapack_7.tex
	./getproto '/void piro_band\_cgbbrd\_l/, /\) ;/' ../Include/piro_band_lapack.h > _piro_band_lapack_8.tex
	./getproto '/void piro_band\_dsbtrd /, /\) ;/' ../Include/piro_band_lapack.h > _piro_band_lapack_9.tex
	./getproto '/void piro_band\_dsbtrd\_l/, /\) ;/' ../Include/piro_band_lapack.h > _piro_band_lapack_10.tex
	./getproto '/void piro_band\_zhbtrd /, /\) ;/' ../Include/piro_band_lapack.h > _piro_band_lapack_11.tex
	./getproto '/void piro_band\_zhbtrd\_l/, /\) ;/' ../Include/piro_band_lapack.h > _piro_band_lapack_12.tex
	./getproto '/void piro_band\_ssbtrd /, /\) ;/' ../Include/piro_band_lapack.h > _piro_band_lapack_13.tex
	./getproto '/void piro_band\_ssbtrd\_l/, /\) ;/' ../Include/piro_band_lapack.h > _piro_band_lapack_14.tex
	./getproto '/void piro_band\_chbtrd /, /\) ;/' ../Include/piro_band_lapack.h > _piro_band_lapack_15.tex
	./getproto '/void piro_band\_chbtrd\_l/, /\) ;/' ../Include/piro_band_lapack.h > _piro_band_lapack_16.tex
	pdflatex UserGuide
	#bibtex UserGuide
	pdflatex UserGuide
	pdflatex UserGuide

distclean: purge

purge: clean
	- $(RM) -r _temp.awk _*.tex $(PURGE)

clean:
	- $(RM) -r $(CLEAN)

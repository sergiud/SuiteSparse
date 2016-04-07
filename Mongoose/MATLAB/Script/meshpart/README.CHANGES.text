
The meshpart toolbox was originally written for Matlab 4.2.

This is a list of the original files, from John Gilbert, that Tim Davis
modified to work with Matlab 5.1, and a summary of the changes required.
The renaming of functions could have been fixed with making sure the Matlab
"path" was set properly, but I prefer not to deal with conflicting names
(you never know for sure which version you're getting).  The functions
etreeplot, radon, and treeplot are all functions in Matlab 5.1.

July 7, 1998.


Contents.m	etreeplot renamed Gilbert_etreeplot

analyze.m	etreeplot renamed Gilbert_etreeplot

chaco.m		changed X==[] to isempty(X)

centerpoint.m	using Gilbert_radon instead of radon

etreeplot.m	renamed to Gilbert_etreeplot.  Using Gilbert_treeplot instead
		of treeplot.

geond.m		using Gilbert_etreeplot instead of etreeplot

gspart.m	changed X==[] to isempty(X)

meshdemo.m	using Gilbert_etreeplot instead of etreeplot

radon.m		radon renamed Gilbert_radon

sepline.m	removed spurious "end;" statement

treeplot.m	renamed to Gilbert_treeplot.  Lots of changes.
		X==[] to isempty(X), and X~=[] to ~isempty(X).
		Removed labeling of nodes with node numbers (rather than
		fix the now-incompatible usage of the "text" function, which
		changed in Matlab 5.1.

./chaco:
Makefile	changed MATLAB, CHACO locations.  Using cc instead of gcc.
		changed CFLAGS and OFLAGS.  Using mex -V4 instead of the
		obsolete cmex.  Note that the *.c files in ./chaco were
		not modified (except user_params.c) to be Matlab 5.1-compliant.

user_params.c	deleted duplicate declarations of global variables:
		SIMULATOR, SIMULATION_ITNS, PERCENTAGE_OUTPUT, CUT_COST,
		HOP_COST, BDY_COST, BDY_HOP_COST, STARTUP_COST.
		The cc compiler complains about this, and aborts the
		compilation.


Jul 2001 and Feb 2002 (John Gilbert):

Further modified the toolbox to:
- add support for the Metis partitioning package 
  (via Robert Bridson's mexfile interface)
- clean up some details for Matlab 6
- rename Gilbert_treeplot to treeplotg
- rename Gilbert_etreeplot to etreeplotg
- rename Gilbert_radon to radong
- rename treelayout to treelayoutg (since treelayout is also a routine in Matlab)
- restore the node number labeling to treeplotg
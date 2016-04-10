% METISMEX : Mex-file interface to METIS
% 
% /****************************************************************************
% * metismex.c
% * Public domain MATLAB CMEX-file to let you use METIS-4.0 from MATLAB.
% * Usage:
% * [part,edgecut] = metismex('PartGraphRecursive',A,nparts,wgtflag,options)
% * [part,edgecut] = metismex('PartGraphKway',A,nparts,wgtflag,options)
% * [perm,iperm] = metismex('EdgeND',A,options)
% * [perm,iperm] = metismex('NodeND',A,options)
% * sep = metismex('NodeBisect',A,wgtflag,options)
% *
% * Output arguments, along with the wgtflag and options input arguments,
% * are optional. See the METIS manual for the meanings of the arguments.
% *
% * Note that error checking is not done: make sure A is structurally
% * symmetric or it will crash.
% *
% * To compile, you need to have Metis 4, and do something like
% *   mex -I<metis.h directory> -L<libmetis.a directory> -lmetis metismex.c
% *
% * Robert Bridson
% *****************************************************************************/
% 
% /* Modified 2 July 2001 by JRG to compile under Windows with lcc */
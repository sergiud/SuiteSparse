function [map,t]=mlchaco(A,Adiag,vwgtsP,ewgtsP,xy,assignment, ...
    architecture,ndims_tot,mesh_dims,goal,global_method,local_method, ...
    rqi_flag,vmax,ndims,eigtol,seed);
% MLCHACO : Hendrickson/Leland's graph partitioner.
%
%  MLCHACO is the mex-file version of the Chaco graph partitioner.
%  The user will normally call CHACO, which calls MLCHACO.
%  See CHACO for a description of parameters.
%
% This is part of a Matlab interface to the software described in
% B. Hendrickson and R. Leland, "The Chaco User's Guide (Version 2.0)",
% Sandia National Laboratories report SAND94-2692, October 1994.
% Interface written by John Gilbert, Xerox PARC, February 1995, and
% Copyright (c) 1994-1996 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.

% Above is the text for HELP MLCHACO; the following will be executed 
% only if the mex-file appropriate for the machine can't be found.

disp('Warning:  Executable mlchaco.mex not found for this architecture.');
disp('          See meshpart/chaco/README for installation instructions.');
map = zeros(1,size(A,2));

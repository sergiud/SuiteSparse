function umfpack_install
%UMFPACK_INSTALL to compile and install umfpack and amd2 for use in MATLAB
%   Your current directory must be UMFPACK/MATLAB for this function to work.
%
% Example:
%   umfpack_install
%
% See also umfpack, amd2.

% UMFPACK, Copyright (c) 2005-2022, Timothy A. Davis, All Rights Reserved.
% SPDX-License-Identifier: GPL-2.0+

% compile and install UMFPACK
umfpack_path = pwd ;
addpath (umfpack_path) ;
umfpack_make

% compile and install AMD
cd ../../AMD/MATLAB
amd_path = pwd ;
addpath (amd_path) ;
amd_make ;

cd (umfpack_path)

fprintf ('Now trying the umfpack_simple demo.\n');
umfpack_simple

fprintf ('Added the following directories to the path.  You may wish to add\n');
fprintf ('these permanently with the MATLAB pathtool command:\n') ;
fprintf ('%s\n', umfpack_path) ;
fprintf ('%s\n', amd_path) ;


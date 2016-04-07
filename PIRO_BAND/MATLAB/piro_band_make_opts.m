function [flags lapack obj_ext include lib] = piro_band_make_opts
%PIRO_BAND_MAKE_OPTS determines compiler flags for compiling PIRO_BAND.
% This function is not meant to be user-callable.
%
% Example:
%   [flags lapack obj_ext include lib] = piro_band_make_opts
%
% See also piro_band_make.

% Copyright 2012, Sivasankaran Rajamanickam, Timothy A. Davis
% http://www.cise.ufl.edu/research/sparse

fprintf (' for MATLAB version %s\n', version) ;

flags = '' ;
is64 = ~isempty (strfind (computer, '64')) ;
if (is64)
    flags = '-largeArrayDims' ;
end

try
    % ispc does not appear in MATLAB 5.3
    pc = ispc ;
    mac = ismac ;
catch                                                                       %#ok
    % if ispc fails, assume we are on a Windows PC if it's not unix
    pc = ~isunix ;
    mac = 0 ;
end

%-------------------------------------------------------------------------------
% BLAS option
%-------------------------------------------------------------------------------

% This is exceedingly ugly.  The MATLAB mex command needs to be told where to
% fine the LAPACK and BLAS libraries, which is a real portability nightmare.

if (pc)
    if (verLessThan ('matlab', '6.5'))
        % MATLAB 6.1 and earlier: use the version supplied here
        lapack = 'lcc_lib/libmwlapack.lib' ;
    elseif (verLessThan ('matlab', '7.5'))
        lapack = 'libmwlapack.lib' ;
    else
        lapack = 'libmwlapack.lib libmwblas.lib' ;
    end
else
    if (verLessThan ('matlab', '7.5'))
        lapack = '-lmwlapack' ;
    else
        lapack = '-lmwlapack -lmwblas' ;
    end
end

if (is64 && ~verLessThan ('matlab', '7.8'))
    % versions 7.8 and later on 64-bit platforms use a 64-bit BLAS
    fprintf ('with 64-bit BLAS\n') ;
    flags = [flags ' -DBLAS64'] ;
end

if (pc)
    obj_ext = '.obj' ;
    include = ' -IWindows' ;
else
    obj_ext = '.o' ;
    include = '' ;
end

lib = ''  ;
if (~(pc || mac))
    % for POSIX timing routine
    lib = [lib ' -lrt'] ;
end

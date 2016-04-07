function piro_sky_make
% piro_sky_make compiles PIRO_SKY package for MATLAB.
% Your current directory must be PIRO_SKY/MATLAB for this function to work.
% PIRO_BAND is also compiled.
%
% Usage :
%   piro_sky_make

%   Copyright 2009, Sivasankaran Rajamanickam, Timothy A. Davis
%   http://www.cise.ufl.edu/research/sparse

details = 0 ;           % 1 if details of each command are to be printed

%-------------------------------------------------------------------------------
% install piro_band, for qr_unsqueeze and piro_band_make_opts
%-------------------------------------------------------------------------------

cd ../../PIRO_BAND/MATLAB
piro_band_make
addpath (pwd)
cd ../../PIRO_SKY/MATLAB

%-------------------------------------------------------------------------------
% compile piro_sky
%-------------------------------------------------------------------------------

fprintf('Compiling PIRO_SKY ...') ;
[flags lapack obj_ext include lib] = piro_band_make_opts ;

include = '-I. -I../Include -I../../SuiteSparse_config -I../../PIRO_BAND/Include -I../../PIRO_BAND/MATLAB -I../../PIRO_BAND/Source' ;

piro_sky_src = { ...
    '../Source/blksky', ...
    '../Source/blksky_symbolic', ...
    '../Source/blksky_utils', ...
    '../Source/genband_givens', ...
    '../../PIRO_BAND/MATLAB/piro_band_mexutil', ...
    '../../PIRO_BAND/Source/piro_band_main', ...
    '../../PIRO_BAND/Source/piro_band_blocksize_main', ...
    '../../PIRO_BAND/Source/piro_band_lapack_main', ...
    '../../PIRO_BAND/Source/piro_band_uv_update_main'} ;

piro_sky_mex_src = { ...
    'blksky', ...
    'blksky_symbolic' } ;

if (ispc)
    obj_extension = '.obj' ;
else
    obj_extension = '.o' ;
end

obj = '' ;
for f = piro_sky_src
    f = strrep (f {1}, '/', filesep) ;
    slash = strfind (f, filesep) ;
    if (isempty (slash))
        slash = 1 ;
    else
        slash = slash (end) + 1 ;
    end
    o = f (slash:end) ;
    obj = [obj  ' ' o obj_extension] ;
    s = sprintf ('mex %s %s -c %s.c', flags, include, f);
    do_cmd (s, details) ;
end

 % compile each mexFunction
 for f = piro_sky_mex_src
    s = sprintf('mex %s %s %s_mex.c -output %s', flags, include, f{1}, f{1}) ;
    s = [s ' ' obj ' ' lapack ' ' lib] ;
    do_cmd (s, details) ;
 end

fprintf ('PIRO_SKY successfully compiled.\n') ;

addpath (pwd) ;

fprintf ('\nThe following path has been added.  You may wish to add \n') ;
fprintf ('it permanently, using the MATLAB pathtool command:\n') ;
fprintf ('%s\n', pwd) ;

fprintf('\nTo run the basic tests cut-and-paste this statement:\n');
fprintf('cd Test ; test_blksky ; test_blksky_symbolic\n') ;

end

%-------------------------------------------------------------------------------

function do_cmd (s, details)
if (details)
    fprintf ('%s\n', s) ;
else
    fprintf('.')
end
eval (s) ;
end

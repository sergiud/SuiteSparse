function piro_band_test_make
%PIRO_BAND_TEST_MAKE compiles the PIRO_BAND test functions for MATLAB.
% These functions are for extensive testing only.  They are not required for
% general usage.  Your current directory must be PIRO_BAND/MATLAB/Test.
%
% Example:
%   piro_band_test_make
%
% See also piro_band, piro_band_svd, piro_band_qr.

% Copyright 2012, Sivasankaran Rajamanickam, Timothy A. Davis
% http://www.cise.ufl.edu/research/sparse

% get compile flags
fprintf ('Compiling PIRO_BAND tests') ;
[flags lapack obj_ext include] = piro_band_make_opts ;

%-------------------------------------------------------------------------------

include = [include ...
    ' -I. -I.. -I../../Include -I../../Source -I../../../SuiteSparse_config'] ;

piro_band_src = { ...
    '../../../SuiteSparse_config/SuiteSparse_config.c', ...
    '../../Source/piro_band_main', ...
    '../../Source/piro_band_blocksize_main', ...
    '../../Source/piro_band_lapack_main', ...
    '../../Source/piro_band_uv_update_main', ...
    '../piro_band_qr_main', ...
    '../piro_band_mexutil' } ;

piro_band_mex_src = { ...
    'piro_band_qr', ...
    'piro_band_lapack' } ;

obj = '' ;
for k = 1:length (piro_band_src)
    f = piro_band_src {k} ;
    f = strrep (f, '/', filesep) ;
    slash = strfind (f, filesep) ;
    slash = slash (end) + 1 ;
    o = f (slash:end) ;
    obj = [obj  ' ' o obj_ext] ; %#ok
    s = sprintf ('mex %s -O %s -c %s.c', flags, include, f);
    fprintf('.')
    eval (s) ;
 end

% compile each mexFunction
for k = 1: length (piro_band_mex_src)
    f = piro_band_mex_src {k} ;
    s = sprintf ('mex -O %s %smex.c -output %s', include, f, f) ;
    s = [s ' ' obj ' ' lapack ' ' flags] ;  %#ok
    fprintf('.')
    eval (s) ;
end

fprintf ('\nPIRO_BAND test functions successfully compiled.\n') ;

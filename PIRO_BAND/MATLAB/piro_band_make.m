function piro_band_make
%PIRO_BAND_MAKE compiles and tests the PIRO_BAND package for MATLAB.
% Your current directory must be PIRO_BAND/MATLAB to use this function.
%
% Example:
%   piro_band_make
%
% See also piro_band, piro_band_svd.

% Copyright 2012, Sivasankaran Rajamanickam, Timothy A. Davis
% http://www.cise.ufl.edu/research/sparse

fprintf ('Compiling PIRO_BAND') ;
[flags lapack obj_ext include lib] = piro_band_make_opts ;

%-------------------------------------------------------------------------------

include = [include ' -I. -I../Include -I../Source -I../../SuiteSparse_config'] ;

piro_band_src = { ...
    '../../SuiteSparse_config/SuiteSparse_config', ...
    '../Source/piro_band_main', ...
    '../Source/piro_band_blocksize_main', ...
    '../Source/piro_band_lapack_main', ...
    '../Source/piro_band_uv_update_main', ...
    'piro_band_qr_main', ...
    'piro_band_mexutil' } ;

piro_band_mex_src = { ...
    'piro_band', ...
    'piro_band_svdkernel' } ;

obj = '' ;
for k = 1:length (piro_band_src)
    f = piro_band_src {k} ;
    f = strrep (f, '/', filesep) ;
    slash = strfind (f, filesep) ;
    if (isempty (slash))
        slash = 1 ;
    else
        slash = slash (end) + 1 ;
    end
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
    s = [s ' ' obj ' ' lapack ' ' flags ' ' lib] ;  %#ok
    fprintf('.')
    eval (s) ;
end

% compile the qr_unsqueeze mexFunction
s = sprintf ('mex %s -O qr_unsqueezemex.c  -output qr_unsqueeze', flags) ;
eval (s) ;

fprintf ('\nPIRO_BAND successfully compiled.\n') ;
addpath (pwd) ;
fprintf ('\nThe following path has been added for this session:\n') ;
fprintf ('    %s\n', pwd) ;
fprintf ('Use the pathtool or savepath commands to add it permanently.\n') ;

% run a quick test
fprintf ('Testing ') ;
load west0479
A = west0479 ;
s1 = piro_band_svd (A) ;
fprintf ('.') ;
s2 = svd (full (A)) ;
fprintf ('.') ;
anorm = s2 (1) ;
err0 = norm (s1-s2) / anorm ;
[U1,S1,V1] = piro_band_svd (A) ;
fprintf ('.') ;
[U2,S2,V2] = svd (full (A)) ;
fprintf ('.') ;
err1 = norm (full (U1*S1*V1' - A)) / anorm ;
err2 = norm (full (U2*S2*V2' - A)) / anorm ;
err = max ([err0 err1 err2]) ;

if (err > 1e-12)
    fprintf ('error %g\n', err) ;
    error ('test failed!') ;
end

fprintf (' all tests passed.\n') ;


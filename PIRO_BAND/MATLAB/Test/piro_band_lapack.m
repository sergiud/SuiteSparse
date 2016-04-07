function [B, U, V] = piro_band_lapack (A,sym,uplo)                          %#ok
%PIRO_BAND_LAPACK reduces a matrix to bidiagonal form, for testing only.
% A must be sparse.  Tests PIRO_BAND's LAPACK-style interface.
%
% Example:
%   [B, U, V] = piro_band_lapack(A) ;
%   norm (A - U*B*V')                   % norm will be small
%   [B, U] = piro_band_lapack (A, 'sym')
%   norm (A - U*B*U')                   % norm will be small
% 
% For the symmetric case, only the upper trianglar part of A is accessed.
% With a struct:
%
% ...  = piro_band_lapack(A, opts)
%
% opts.sym and opts.uplo are optional arguments. If opts.sym=1 then the matrix
% is assumed to be symmetric.   If opts.uplo is 'U', or not present, then only
% the upper triangular part is considered, and the lower triangular part is
% assumed to be the transpose of the upper part.  If opts.uplo is 'L', then the
% lower triangular part of A is used instead.
%
% See also piro_band.

% Copyright 2012, Sivasankaran Rajamanickam, Timothy A. Davis
% http://www.cise.ufl.edu/research/sparse

error('piro_band_lapack mexFunction not found') ;


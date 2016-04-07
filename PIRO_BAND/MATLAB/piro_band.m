function [B, U, V]  = piro_band (A, vargin)                                 %#ok
%PIRO_BAND reduce a band matrix to upper bidiagonal form.
%
% piro_band reduces a band matrix (stored in sparse or full format) to an upper
% bidiagonal matrix using blocked and interleaved Givens rotations.
%
% Usage:
%
% [B, U, V] = piro_band (A)
% [B, U, V] = piro_band (A, opts)
%
% A is m-by-n and can be real or complex.  n must be 2 or more.  B is real, and
% is the same size as A.  It is upper bidiagonal for the unsymmetric case, and
% symmetric tridiagonal for the symmetric case.  A is equal to U*B*V', where U
% and V are full orthogonal matrices.  If A is m-by-n, B is the same size, U is
% m-by-m, and V is n-by-n.  A and B have the same singular values.
%
% opts can be a single parameter containing the following fields:
%
%   opts.sym:  0 if A is unsymmetric, 1 if symmetric.  Default is 0.
%   opts.blks: a vector of size 4
%
% alternatively, you can specify a list of individual arguments in any order.
% For example:
%
%   ... = piro_band (A, 'sym', [32 32 64 64]) ;
%
% where 'sym' is the same as opts.sym=1, and a vector of size 4 is used for
% the opts.blks parameter.
%
% If opts.sym=1, then A must be square and is assumed to be symmetric.  Only
% the entries in the upper triangular part are considered, and the lower
% triangular part is assumed to be the transpose of the upper part.  The
% symmetric matrix C = A + tril(A',-1) is operated on via symmetric reductions,
% and B is a symmetric tridiagonal matrix.  B and C have the same singular
% values, and the same eigenvalues.
%
% blks is an array of size 4 that determines the block sizes used in the block
% reduction for the upper and lower block sizes respectively.  The block size
% for the upper band is blks(1)-by-blks(2), and the lower band is reduced with
% blocks of size blks(3)-by-blks(4).  The opts.blks option is primarily meant
% for performance experiments, since the default block sizes usually give the
% best performance.
%
% In contrast to piro_band_svd, no fill reducing ordering is used.  Using
% symrcm and permuting the matrix prior to calling piro_band is a good option
% for reducing the bandwidth and thus the total work.
%
% Example:
%
%   A = rand (4) ;
%   [B,U,V] = piro_band (A) ;
%   A - U*B*V'
%   svd (A) - svds (B)
%
% See also piro_band_svd, svd, symrcm.

% Copyright 2012, Sivasankaran Rajamanickam, Timothy A. Davis
% http://www.cise.ufl.edu/research/sparse

error ('piro_band mexFunction not found')


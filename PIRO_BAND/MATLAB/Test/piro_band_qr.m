function [Q, R] = piro_band_qr (A)                                          %#ok
%PIRO_BAND_QR computes the QR factorization of a band matrix A.
% A is a band matrix in sparse or full format.
% This function is meant for testing only.
%
% Example:
%
% [Q,R] = piro_band_qr(A) ;
%
%      If A is m-by-n with m <= n, then
%          Q is m-by-m full matrix, and
%          R is m-by-n sparse or full matrix.
%      If m > n, then
%          Q is m-by-n full matrix, and
%          R is n-by-n sparse or full matrix.
%      This is the same as [Q,R] = qr(A,0) in MATLAB.
%      The full or sparse storage format of R is the same as A.
%
%  QR = piro_band_qr(A) ;
%
%      returns QR as a struct with fields V, Beta, R, bl and bu.
%      V        Householder vectors stored in band format.
%      beta     beta parameters for each Householder transformation.
%      R        upper triangular matrix in band format.
%      bl       lower bandwidth of A.
%      bu       upper bandwidth of A.
%
% See also qr, piro_band.

% Copyright 2012, Sivasankaran Rajamanickam, Timothy A. Davis
% http://www.cise.ufl.edu/research/sparse

error ('piro_band_qr mexFunction not found')

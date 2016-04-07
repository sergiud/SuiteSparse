% PIRO_BAND: a package for computing the singular value decomposition and QR
% factorization of band matrices, and for reducing band matrices to bidiagonal
% or symmetric tridiagonal form via orthogonal operations.
%
% Copyright 2012, Sivasankaran Rajamanickam and Timothy A. Davis,
% University of Florida.
%
% Files
%   piro_band      - reduce a band matrix to upper bidiagonal form.
%   piro_band_svd  - singular value decomposition of a sparse or full matrix.
%   piro_band_make - compiles and tests the PIRO_BAND package for MATLAB.

% Files not meant to be called by the end-user:
%   piro_band_make_opts - determines compiler flags for compiling PIRO_BAND.
%   piro_band_svdkernel - this is not meant to be called directry by the end-user.
%   qr_unsqueeze        - returns a row permutation that 'unsqueezes' a rank-deficient R.

% TODO run MATLAB source code tools

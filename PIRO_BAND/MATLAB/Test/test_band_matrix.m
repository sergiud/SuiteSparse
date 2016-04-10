function A = test_band_matrix (m, n, hi, lo, kind)
%TEST_BAND_MATRIX generates a small sparse banded matrix, with random values.
%
% A = test_band_matrix (m, n, hi, lo)
%
% generate a band matrix with both upper bandwidth hi and lower bandwidth lo.
% hi=1 and lo=1 gives a tridiagonal matrix.
%
% Example:
%
%   A = test_band_matrix (10, 5, 2, 3)
%
% See also spdiags, randn, tril, triu.

% Copyright 2012, Sivasankaran Rajamanickam, Timothy A. Davis
% http://www.cise.ufl.edu/research/sparse

if (nargin < 5)
    kind = 'real' ;
end

if (strcmp (kind , 'complex'))
    Areal = test_band_matrix (m, n, hi, lo) ;
    Aimag = test_band_matrix (m, n, hi, lo) ;
    A = Areal + 1i * Aimag ;
else
    lo = abs (lo) ;
    hi = abs (hi) ;
    A = spdiags (randn (min (m,n), lo+1+hi), -lo:hi, m, n) ;
end

% assert (nnz (tril (A, -lo-1)) == 0) ;
% assert (nnz (triu (A,  hi+1)) == 0) ;


function err = test_piro_band_special (pr)
%TEST_PIRO_BAND_SPECIAL tests piro_band functions on a set of matrices.
% Example:
%   test_piro_band_special
%
% Tests the band reduction method with random input matrices.
%
% See also test_piro_band.

% Copyright 2012, Sivasankaran Rajamanickam, Timothy A. Davis
% http://www.cise.ufl.edu/research/sparse

if (nargin == 0)
    pr = 0 ;
end

err = 0 ;

fprintf ('\ntest piro_band: ') ;

if (pr)
    fprintf ('\n--------------Symmetric case  ---------------- \n\n') ;
end

% Only the upper half is stored.

% tridiagonal case
A = test_band_matrix (10, 10, 1, 0) ;
err = max (err, piro_band_sym_verify (A, [0 0 0 0], pr, 0)) ;

A = test_band_matrix (10, 10, 2, 0) ;
err = max (err, piro_band_sym_verify (A, [1 1 0 0], pr, 1)) ;

A = test_band_matrix (10, 10, 5, 0) ;
err = max (err, piro_band_sym_verify (A, [2 2 0 0], pr, 0)) ;

A = test_band_matrix (10, 10, 5, 0) ;
err = max (err, piro_band_sym_verify (A, [2 3 0 0], pr, 1)) ;

if (pr)
    fprintf ('\n--------------Complex case  ---------------- \n\n') ;
end

% Only the upper half is stored.

% tridiagonal case
A = test_band_matrix (10, 10, 1, 0) ;
A1 = test_band_matrix (10, 10, 1, 0) ;
d = diag (A1) ;
A1 = A1 - diag (d) ;
A = A + 1i * A1 ;
err = max (err, piro_band_sym_verify (A, [0 0 0 0], pr, 0)) ;

A = test_band_matrix (10, 10, 2, 0) ;
A1 = test_band_matrix (10, 10, 2, 0) ;
d = diag (A1) ;
A1 = A1 - diag (d) ;
A = A + 1i * A1 ;
err = max (err, piro_band_sym_verify (A, [1 1 0 0], pr, 1)) ;

A = test_band_matrix (10, 10, 5, 0) ;
A1 = test_band_matrix (10, 10, 5, 0) ;
d = diag (A1) ;
A1 = A1 - diag (d) ;
A = A + 1i * A1 ;
err = max (err, piro_band_sym_verify (A, [2 2 0 0], pr, 0)) ;

A = test_band_matrix (10, 10, 5, 0) ;
A1 = test_band_matrix (10, 10, 5, 0) ;
d = diag (A1) ;
A1 = A1 - diag (d) ;
A = A + 1i * A1 ;
err = max (err, piro_band_sym_verify (A, [2 3 0 0], pr, 1)) ;

if (pr)
    fprintf ('\n-------------- LAPACK style interface ------------- \n\n') ;
    fprintf ('\n--------------Symmetric case  ---------------- \n\n') ;
end

% Only the upper half is stored.
A = test_band_matrix (10, 10, 0, 0) ;
err = max (err, piro_band_lapack_sym_verify (A, 0, pr, 1)) ;
err = max (err, piro_band_lapack_sym_verify (A, 0, pr, 0)) ;

% tridiagonal case
A = test_band_matrix (10, 10, 1, 0) ;
err = max (err, piro_band_lapack_sym_verify (A, 0, pr, 1)) ;
err = max (err, piro_band_lapack_sym_verify (A, 0, pr, 0)) ;

A = test_band_matrix (10, 10, 2, 0) ;
err = max (err, piro_band_lapack_sym_verify (A, 0, pr, 1)) ;
err = max (err, piro_band_lapack_sym_verify (A, 0, pr, 0)) ;

A = test_band_matrix (10, 10, 5, 0) ;
err = max (err, piro_band_lapack_sym_verify (A, 0, pr, 1)) ;
err = max (err, piro_band_lapack_sym_verify (A, 0, pr, 0)) ;

A = test_band_matrix (10, 10, 6, 0) ;
err = max (err, piro_band_lapack_sym_verify (A, 0, pr, 1)) ;
err = max (err, piro_band_lapack_sym_verify (A, 0, pr, 0)) ;

% Only the lower half is stored.
% tridiagonal case
A = test_band_matrix (10, 10, 0, 1) ;
err = max (err, piro_band_lapack_sym_verify (A, 1, pr, 1)) ;
err = max (err, piro_band_lapack_sym_verify (A, 1, pr, 0)) ;

A = test_band_matrix (10, 10, 0, 2) ;
err = max (err, piro_band_lapack_sym_verify (A, 1, pr, 1)) ;
err = max (err, piro_band_lapack_sym_verify (A, 1, pr, 0)) ;

A = test_band_matrix (10, 10, 0, 5) ;
err = max (err, piro_band_lapack_sym_verify (A, 1, pr, 1)) ;
err = max (err, piro_band_lapack_sym_verify (A, 1, pr, 0)) ;

A = test_band_matrix (10, 10, 0, 6) ;
err = max (err, piro_band_lapack_sym_verify (A, 1, pr, 1)) ;
err = max (err, piro_band_lapack_sym_verify (A, 1, pr, 0)) ;

if (pr)
    fprintf ('\n--------------Complex case  ---------------- \n\n') ;
end

% Only the upper half is stored.
A = test_band_matrix (10, 10, 0, 0) ;
A1 = test_band_matrix (10, 10, 0, 0) ;
d = diag (A1) ;
A1 = A1 - diag (d) ;
A = A + 1i * A1 ;
err = max (err, piro_band_lapack_sym_verify (A, 0, pr, 1)) ;
err = max (err, piro_band_lapack_sym_verify (A, 0, pr, 0)) ;

% tridiagonal case
A = test_band_matrix (10, 10, 1, 0) ;
A1 = test_band_matrix (10, 10, 1, 0) ;
d = diag (A1) ;
A1 = A1 - diag (d) ;
A = A + 1i * A1 ;
err = max (err, piro_band_lapack_sym_verify (A, 0, pr, 1)) ;
err = max (err, piro_band_lapack_sym_verify (A, 0, pr, 0)) ;

A = test_band_matrix (10, 10, 2, 0) ;
A1 = test_band_matrix (10, 10, 2, 0) ;
d = diag (A1) ;
A1 = A1 - diag (d) ;
A = A + 1i * A1 ;
err = max (err, piro_band_lapack_sym_verify (A, 0, pr, 1)) ;
err = max (err, piro_band_lapack_sym_verify (A, 0, pr, 0)) ;

A = test_band_matrix (10, 10, 5, 0) ;
A1 = test_band_matrix (10, 10, 5, 0) ;
d = diag (A1) ;
A1 = A1 - diag (d) ;
A = A + 1i * A1 ;
err = max (err, piro_band_lapack_sym_verify (A, 0, pr, 1)) ;
err = max (err, piro_band_lapack_sym_verify (A, 0, pr, 0)) ;

A = test_band_matrix (10, 10, 5, 0) ;
A1 = test_band_matrix (10, 10, 5, 0) ;
d = diag (A1) ;
A1 = A1 - diag (d) ;
A = A + 1i * A1 ;
err = max (err, piro_band_lapack_sym_verify (A, 0, pr, 1)) ;
err = max (err, piro_band_lapack_sym_verify (A, 0, pr, 0)) ;

% Only the lower half is stored.
% tridiagonal case
A = test_band_matrix (10, 10, 0, 1) ;
A1 = test_band_matrix (10, 10, 0, 1) ;
d = diag (A1) ;
A1 = A1 - diag (d) ;
A = A + 1i * A1 ;
err = max (err, piro_band_lapack_sym_verify (A, 1, pr, 1)) ;
err = max (err, piro_band_lapack_sym_verify (A, 1, pr, 0)) ;

A = test_band_matrix (10, 10, 0, 2) ;
A1 = test_band_matrix (10, 10, 0, 2) ;
d = diag (A1) ;
A1 = A1 - diag (d) ;
A = A + 1i * A1 ;
err = max (err, piro_band_lapack_sym_verify (A, 1, pr, 1)) ;
err = max (err, piro_band_lapack_sym_verify (A, 1, pr, 0)) ;

A = test_band_matrix (10, 10, 0, 5) ;
A1 = test_band_matrix (10, 10, 0, 5) ;
d = diag (A1) ;
A1 = A1 - diag (d) ;
A = A + 1i * A1 ;
err = max (err, piro_band_lapack_sym_verify (A, 1, pr, 1)) ;
err = max (err, piro_band_lapack_sym_verify (A, 1, pr, 0)) ;

A = test_band_matrix (10, 10, 0, 6) ;
A1 = test_band_matrix (10, 10, 0, 6) ;
d = diag (A1) ;
A1 = A1 - diag (d) ;
A = A + 1i * A1 ;
err = max (err, piro_band_lapack_sym_verify (A, 1, pr, 1)) ;
err = max (err, piro_band_lapack_sym_verify (A, 1, pr, 0)) ;

fprintf ('\nerror: %g OK', err) ;


%-------------------------------------------------------------------------------
% Find the SVD for symmetric case and verify the result.
function err = piro_band_sym_verify (A, blks, pr, usestruct) 
[m, n] = size (A) ;
d = diag (A) ;
A2 = A + A' - diag (d) ;
sgood = svd (full (A2)) ;

if (usestruct)
    % test the opts struct method
    opts.sym = 1 ;
    opts.blks = blks ;
    [B, U] = piro_band (A, opts) ;
else
    % test the non-struct method
    [B, U] = piro_band (A, 'sym', blks) ;
end
if (~issparse (B))
    error ('B not sparse?') ;
end
if (exist ('spok') == 3)    %#ok
    spok (B) ;
end

s = svd (full (B)) ;
err1 = norm (sgood-s) / sgood (1) ;

A1 = U * B * U' ;
err2 = norm (A1 - A2) / sgood (1) ;

err = test_piro_band_error (m, n, err1, err2, pr) ;

%-------------------------------------------------------------------------------
% Find the symmetric SVD using LAPACK style call and verify the result.
function err = piro_band_lapack_sym_verify (A, uplo, pr, usestruct) 

[m, n] = size (A) ;
d = diag (A) ;
A2 = A + A' - diag (d) ;
sgood = svd (full (A2)) ;

if (uplo == 0)
    uplo = 'upper' ;
else
    uplo = 'lower' ;
end

if (usestruct)
    opts.sym = 1 ;
    opts.uplo = uplo ;
    [B, U] = piro_band_lapack (A, opts) ;
else
    [B, U] = piro_band_lapack (A, 'sym', uplo) ;
end
if (~issparse (B))
    error ('B not sparse?') ;
end
if (exist ('spok') == 3)    %#ok
    spok (B) ;
end

s = svd (full (B)) ;
err1 = norm (sgood-s) / sgood (1) ;

A1 = U * B * U' ;
err2 = norm (A1 - A2) / sgood (1) ;

% repeat, but just get B
if (usestruct)
    B2 = piro_band_lapack (A, opts) ;
else
    B2 = piro_band_lapack (A, 'sym', uplo) ;
end
if (~issparse (B2))
    error ('B2 not sparse?') ;
end
if (exist ('spok') == 3)    %#ok
    spok (B2) ;
end
if (~isequal (B,B2))
    error ('B vs B2 mismatch') ;
end

err = test_piro_band_error (m, n, err1, err2, pr) ;

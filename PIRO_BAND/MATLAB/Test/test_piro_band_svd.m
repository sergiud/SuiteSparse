function maxerr = test_piro_band_svd (pr)
%TEST_PIRO_BAND_SVD tests piro_band_svd on a set of matrices.
% Example:
%   test_piro_band_svd
%
% Finds the SVD using the band reduction algorithm for random input matrices
% and compares the result with the MATLAB svd.  Uses random input matrices.
%
% See also test_piro_band.

% Copyright 2012, Sivasankaran Rajamanickam, Timothy A. Davis
% http://www.cise.ufl.edu/research/sparse

if (nargin < 1)
    pr = 0 ;
end

maxerr = 0 ;
kinds = {'real', 'complex'} ;

for t = 1:2
    err = 0 ;
    kind = kinds {t} ;
    fprintf ('\ntest piro_band_svd, %s case:\n', kind) ;

    % diagonal case
    A = test_band_matrix (10, 10, 0, 0, kind) ;
    err = max (err, piro_band_verify (A, [0 0 0 0], pr)) ;

    A = test_band_matrix (10, 20, 0, 0, kind) ;
    err = max (err, piro_band_verify (A, [0 0 0 0], pr)) ;

    A = test_band_matrix (20, 10, 0, 0, kind) ;
    err = max (err, piro_band_verify (A, [0 0 0 0], pr)) ;

    % bidiagonal case
    A = test_band_matrix (10, 10, 1, 0, kind) ;
    err = max (err, piro_band_verify (A, [0 0 0 0], pr)) ;

    A = test_band_matrix (10, 20, 1, 0, kind) ;
    err = max (err, piro_band_verify (A, [0 0 0 0], pr)) ;

    A = test_band_matrix (20, 10, 1, 0, kind) ;
    err = max (err, piro_band_verify (A, [0 0 0 0], pr)) ;

    % tridiagonal case
    A = test_band_matrix (10, 10, 1, 1, kind) ;
    err = max (err, piro_band_verify (A, [0 0 1 1], pr)) ;

    A = test_band_matrix (10, 20, 1, 1, kind) ;
    err = max (err, piro_band_verify (A, [0 0 1 1], pr)) ;

    A = test_band_matrix (20, 10, 1, 1, kind) ;
    err = max (err, piro_band_verify (A, [0 0 1 1], pr)) ;

    fprintf ('\n') ;

    % upper banded case
    A = test_band_matrix (10, 10, 5, 0, kind) ;
    err = max (err, piro_band_verify (A, [2 2 0 0], pr)) ;

    A = test_band_matrix (10, 20, 6, 0, kind) ;
    err = max (err, piro_band_verify (A, [3 3 0 0], pr)) ;

    A = test_band_matrix (20, 10, 5, 0, kind) ;
    err = max (err, piro_band_verify (A, [4 1 0 0], pr)) ;

    % lower Hessenberg case
    A = test_band_matrix (10, 10, 1, 5, kind) ;
    err = max (err, piro_band_verify (A, [0 0 1 5], pr)) ;

    A = test_band_matrix (10, 20, 1, 7, kind) ;
    err = max (err, piro_band_verify (A, [0 0 7 1], pr)) ;

    A = test_band_matrix (20, 10, 1, 6, kind) ;
    err = max (err, piro_band_verify (A, [0 0 3 3], pr)) ;

    % lower banded case
    A = test_band_matrix (10, 10, 0, 5, kind) ;
    err = max (err, piro_band_verify (A, [0 0 1 5], pr)) ;

    A = test_band_matrix (10, 20, 0, 7, kind) ;
    err = max (err, piro_band_verify (A, [0 0 7 1], pr)) ;

    A = test_band_matrix (20, 10, 0, 6, kind) ;
    err = max (err, piro_band_verify (A, [0 0 3 3], pr)) ;

    fprintf ('\n') ;

    % scalar
    A = test_band_matrix (1, 1, 0, 0, kind) ;
    err = max (err, piro_band_verify (A, [ ], pr)) ;

    % 2-by-2diagonal
    A = test_band_matrix (2, 2, 0, 0, kind) ;
    err = max (err, piro_band_verify (A, [ ], pr)) ;

    % row vector
    A = test_band_matrix (1, 10, 1, 1, kind) ;
    err = max (err, piro_band_verify (A, [ ], pr)) ;

    % column vector
    A = A' ;
    err = max (err, piro_band_verify (A, [ ], pr)) ;

    % small banded matrices
    A = test_band_matrix (10, 10, 5, 5, kind) ;
    err = max (err, piro_band_verify (A, [2 2 2 2], pr)) ;

    err = max (err, piro_band_verify (full (A), [ ], pr)) ;

    A = test_band_matrix (10, 10, 6, 5, kind) ;
    err = max (err, piro_band_verify (A, [3 2 2 2], pr)) ;

    A = test_band_matrix (10, 20, 6, 9, kind) ;
    err = max (err, piro_band_verify (A, [3 2 3 2], pr)) ;

    A = test_band_matrix (10, 20, 5, 5, kind) ;
    err = max (err, piro_band_verify (A, [ ], pr)) ;

    fprintf ('\n') ;

    % with an empty row
    state = warning ;       % save existing warning state
    warning ('off', 'PIRO_BAND_SVD:noSPQR') ;
    lastwarn ('') ;
    A (3,:) = 0 ;
    err = max (err, piro_band_verify (A, [ ], pr)) ;
    [~, id] = lastwarn ;
    err = max (err, double (~strcmp (id, 'PIRO_BAND_SVD:noSPQR'))) ;
    warning (state) ;       % restore warning state
    err = test_piro_band_error (size (A,1), size(A,2), err, err, pr) ;

    A = test_band_matrix (20, 10, 5, 5, kind) ;
    err = max (err, piro_band_verify (A, [ ], pr)) ;

    A = test_band_matrix (20, 10, 5, 6, kind) ;
    err = max (err, piro_band_verify (A, [ ], pr)) ;

    A = test_band_matrix (20, 10, 7, 8, kind) ;
    err = max (err, piro_band_verify (A, [3 3 3 3], pr)) ;

    A = test_band_matrix (20, 30, 11, 11, kind) ;
    err = max (err, piro_band_verify (A, [3 5 3 5], pr)) ;

    A = test_band_matrix (30, 20, 11, 11, kind) ;
    err = max (err, piro_band_verify (A, [3 5 3 5], pr)) ;

    A = test_band_matrix (50, 60, 6, 9, kind) ;
    err = max (err, piro_band_verify (A, [ ], pr)) ;

    A = test_band_matrix (60, 50, 6, 9, kind) ;
    err = max (err, piro_band_verify (A, [ ], pr)) ;

    A = test_band_matrix (50, 50, 6, 9, kind) ;
    err = max (err, piro_band_verify (A, [ ], pr)) ;

    A = test_band_matrix (100, 100, 30, 30, kind) ;
    err = max (err, piro_band_verify (A, [10 10 10 10], pr)) ;

    A = test_band_matrix (100, 100, 30, 40, kind) ;
    err = max (err, piro_band_verify (A, [15 15 17 9], pr)) ;

    fprintf ('\n[') ;

    % using the west0479 test matrix, distributed with MATLAB
    load west0479 ;
    A = west0479 ;
    if (t == 2)
        A = A + 1i * sprandn (A) ;
    end
    err = max (err, piro_band_verify (A, [ ], pr)) ;

    fprintf (']\n[') ;

    % with a dense row
    A (1,:) = 42 ;
    err = max (err, piro_band_verify (A, [ ], pr)) ;

    fprintf (']\n%s case error: %g OK\n', kind, err) ;
    maxerr = max (maxerr, err) ;

end

%-------------------------------------------------------------------------------
function err = piro_band_verify (A, blks, pr)
[m, n] = size (A) ;
k = min (m,n) ;
Smatlab = svd (full (A)) ;
err = 0 ;

% using an opts struct
for ordering = {'rcm', 'colamd', 'amd', 'none'}
    for qropt = 0:2
        for econ = 0:1
            opts = struct ('econ', econ, 'qr', qropt, 'ordering', ordering) ;
            S = piro_band_svd (A, opts) ;
            err = max (err, check_svd (A, [ ], [ ], S, [ ], Smatlab, pr)) ;
            [U, S, V] = piro_band_svd (A, opts);
            err = max (err, check_svd (A, [ ], U, S, V, Smatlab, pr)) ;
            [U1, S] = piro_band_svd (A, opts);
            err = max (err, check_svd (A, U1, U, S, [ ], Smatlab, pr)) ;
        end
    end
end

% economy, with a string option

S = piro_band_svd (A, 'econ') ;
err = max (err, check_svd (A, [ ], [ ], S, [ ], Smatlab, pr)) ;
[U, S, V] = piro_band_svd (A, 'econ');
err = max (err, check_svd (A, [ ], U, S, V, Smatlab, pr)) ;
[U1, S] = piro_band_svd (A, 'econ');
err = max (err, check_svd (A, U1, U, S, [ ], Smatlab, pr)) ;

% non-economy, no opts argument

[U, S, V] = piro_band_svd (A) ;
err = max (err, check_svd (A, [ ], U, S, V, Smatlab, pr)) ;
S = piro_band_svd (A) ;
err = max (err, check_svd (A, [ ], [ ], S, [ ], Smatlab, pr)) ;
[U1, S] = piro_band_svd (A) ;
err = max (err, check_svd (A, U1, U, S, [ ], Smatlab, pr)) ;

% error-handling (error is expected)
if (k >= 2)
    try
        S = piro_band_svd (A, struct ('ordering', 'gunk')) ;
        err = 1 ;
    catch me                                                                %#ok
    end
end
err = test_piro_band_error (m, n, err, err, pr) ;

if (n > 1)
    % test piro_band
    if (~isempty (blks))
        [B, U, V] = piro_band (A, struct ('blks', blks)) ;
        err = max (err, check_band (A, U, B, V, Smatlab, pr)) ;
        [B, U, V] = piro_band (A, blks) ;
        err = max (err, check_band (A, U, B, V, Smatlab, pr)) ;
    end
    [B, U, V] = piro_band (A) ;
    err = max (err, check_band (A, U, B, V, Smatlab, pr)) ;
    [B, U] = piro_band (A) ;
    err = max (err, check_band (A, U, B, [ ], Smatlab, pr)) ;
    B = piro_band (A) ;
    err = max (err, check_band (A, [ ], B, [ ], Smatlab, pr)) ;
    % test the LAPACK interface
    [B, U, V] = piro_band_lapack (sparse (A)) ;
    err = max (err, check_band (A, U, B, V, Smatlab, pr)) ;
    [B, U] = piro_band_lapack (sparse (A)) ;
    err = max (err, check_band (A, U, B, [ ], Smatlab, pr)) ;
    B = piro_band_lapack (sparse (A)) ;
    err = max (err, check_band (A, [ ], B, [ ], Smatlab, pr)) ;
end

% test the band QR
[Q, R] = piro_band_qr (A) ;
err = max (err, check_band (A, Q, R, speye (n), Smatlab, pr, issparse (A))) ;

QR = piro_band_qr (A) ;
[Q, R] = piro_band_extract_qr (QR) ;
err = max (err, check_band (A, Q, R, speye (n), Smatlab, pr, issparse (A))) ;

if (~pr)
    fprintf ('.') ;
end

%-------------------------------------------------------------------------------

function err = check_svd (A, U1, U, S, V, Smatlab, pr)
% err = check_svd (U, S, V, Smatlab, pr)
[m, n] = size (A) ;
k = min (m,n) ;
if (issparse (S))
    % A = U*S*V', where S is a sparse diagonal matrix
    if (exist ('spok') == 3)    %#ok
        spok (S) ;
    end
    if (isempty (U1))
        err2 = norm (U*S*V' - A) / Smatlab (1) ;
    else
        err2 = norm (U - U1) ;
    end
    if (k > 1)
        S = diag (S) ;
    elseif (k == 1)
        S = S (1) ;
    end
else
    % S = svd (A) where S is a dense vector
    err2 = 0 ;
end
err1 = norm (Smatlab - S) / Smatlab (1) ;
err = test_piro_band_error (m, n, err1, err2, pr) ;

%-------------------------------------------------------------------------------
function err = check_band (A, U, B, V, Smatlab, pr, Bsparse)
% err = check_band (A, U, B, V, Smatlab, pr)
% A should equal U*B*V', and the singular values of A and B should match.
% Smatlab on input are the singular values of A.  If V is empty, A == U*B*V'
% is not checked.
[m, n] = size (A) ;
k = min (m,n) ;
if (nargin < 7)
    Bsparse = false ;
end
if (Bsparse)
    % B should be sparse
    if (~issparse (B))
        error ('B should be sparse') ;
    end
    if (exist ('spok') == 3)    %#ok
        spok (B) ;
    end
end
s = svd (full (B)) ;
err1 = norm (Smatlab-s) / Smatlab (1) ;
err2 = 0 ;
if (~isempty (V))
    A1 = U * B * V' ;
    if (issparse (A1))
        err2 = norm (A1 - A, 1) / Smatlab (1) ;
    else
        err2 = norm (A1 - A) / Smatlab (1) ;
    end
end
err = test_piro_band_error (m, n, err1, err2, pr) ;

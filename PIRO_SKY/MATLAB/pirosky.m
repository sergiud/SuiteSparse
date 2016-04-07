function [s, stats, ss] = pirosky (A, s)
%PIROSKY s = pirosky(A) computes all the singular values of a sparse matrix A.
% Default parameter settings are used.
%
% s = pirosky (A, opts) uses non-default parameter settings:
%
%   opts.ordering:  The sparse QR factoriztion of A(:,q) is computed, where
%       q is selected according to the following options:
%       'natural':  no preordering of A.
%       'symrcm':   q = symrcm(A'*A).  This is the default.
%       'colamd':   q = colamd(A), reduces the fill-in in R better than symrcm,
%                   but often takes more work in the skyline reduction of R to
%                   bidiagonal.
%
%   opts.demo:  if true, display a demo and print diagnostic info
%   opts.tol:   same as tol for svd.  If < 0, then use the default for svd.
%
% [s,stats] = pirosky  (...) returns optional statistics
%
%   stats.rnz       nnz(R) for the QR factorization of A
%   stats.rank_est  rank as esimated by qr or spqr
%   stats.rank      numerical rank, computed from singular values
%
% Requires the PIRO_BAND/MATLAB/qr_unsqueeze mexFunction

% TODO also return U and V

if (~isreal (A))
    error ('complex sparse matrices not yet supported') ;
end

if (~issparse (A))
    % force A to be sparse
    A = sparse (A) ;
end

% default options
if (nargin < 2)
    s = 'symrcm' ;
end
if (ischar (s))
    opts.ordering = s ;
elseif (isstruct (s))
    opts = s ;
else
    error ('opts must be a string or a struct') ;
end

if (~isfield (opts, 'ordering'))
    opts.ordering = 'symrcm' ;
end
if (~isfield (opts, 'demo'))
    opts.demo = 0 ;
end
if (~isfield (opts, 'tol'))
    opts.tol = -1 ;
end
if (~isfield (opts, 'name'))
    opts.name = 'A' ;
end

[m n] = size (A) ;
m1 = m ;
n1 = n ;

if (opts.demo)
    fprintf ('pirosky:  A is %d-by-%d with %d nonzeros\n', m, n, nnz (A)) ;
    subplot (2,2,1) ;
    spy2 (A, opts.name) ;
    drawnow ;
end

ss {1} = svd (full (A)) ;       % hack test

% make sure A is tall and thin, not short and fat
if (m < n)
    A = A' ;
    [m n] = size (A) ;
end

ss {2} = svd (full (A)) ;       % hack test

% order the columns of A
switch opts.ordering
    case 'natural'
        ;
    case 'symrcm'
        S = spones (A) ;
        q = symrcm (S'*S) ;
        clear S
        A = A (:,q) ;
    case 'colamd'
        q = colamd (A) ;
        A = A (:,q) ;
    default
        error ('unknown ordering option') ;
end

% TODO need the permutation q for U and/or V

ss {3} = svd (full (A)) ;       % hack test

if (opts.demo)
    subplot (2,2,2) ;
    spy2 (A, 'permuted A') ;
    drawnow ;
end

% QR factorization of A, or permuted version of A.  A is not permuted,
% except that R(1:r,1:r) should have a zero-free diagonal.
spqr_opts.ordering = 'natural' ;
spqr_opts.Q = 'discard' ;                   % do not keep Q (TODO: needed for U)
spqr_opts.tol = 0 ;                         % rank-detecting tolerance
[Q,R,P,info] = spqr (A, spqr_opts) ;
info
r = info.rank_A_estimate ;
clear Q P                                   % TODO keep Q and P for U and V

ss {4} = svd (full (R)) ;       % hack test

if (r < (4/5)*n)    % TODO determine this threshold

    % If the rank is very small, then blksky will be very inefficient.
    % So R is transposed and factorized again.

    % TODO need to handle U and V here too.  The old Q is V.

    fprintf ('          transposing R and doing second QR\n') ;

    % TODO use try/catch, and punt to MATLAB qr
    R = R (1:r,:)' ;
    [Q,R,P,info] = spqr (R, spqr_opts) ;
    r2 = info.rank_A_estimate ;

ss {5} = svd (full (R)) ;       % hack test

    if (r2 < r)
        % an unusual case
        [p,r] = qr_unsqueeze (R) ;
        R = R (p,:) ;
        % TODO need to keep p for U and V
    end
    [m n] = size (R) ;

elseif (r < n)

ss {5} = svd (full (R)) ;       % hack test

    % the matrix has substantial rank, so just unsqueeze it
    fprintf ('          just unsqueeze\n') ;
    [p,r] = qr_unsqueeze (R) ;
    R = R (p,:) ;

end

ss {6} = svd (full (R)) ;       % hack test

stats.rank_est = r ;

% ensure R is square for blksky
size (R)
n
R = R (1:n, 1:n) ;
% svd (full (R)) - svd (full (A))
% save R R

stats.rnz = nnz (R) ;

if (opts.demo)
    subplot (2,2,3) ;
    spy2 (R, 'R') ;
    fprintf ('          R is %d-by-%d nnz(R): %d nnz(diag(R)): %d\n', ...
        n, n, stats.rnz, nnz(diag(R))) ;
    fprintf ('          estimated rank: %d\n', stats.rank_est) ;
    drawnow ;
end

save R R
% reduce the upper triangular sparse R to upper bidiagonal
B = blksky (R) ;
clear R

ss {7} = svd (full (B)) ;       % hack test

% find the singular values of the bidiagonal matrix
if (n < 2 || nnz (diag (B,1)) == 0)
    s = abs (full (diag (B))) ;
else
    s = piro_band_svd (B) ;
end

% pad s with zeros for the rank-deficient case
if (length (s) < min (m1,n1))
    s (min (m1,n1)) = 0 ;
end

ss {8} = s ;       % hack test

% ensure s is a column vector
if (size (s,2) > 1)
    s = s' ;
end

s = sort (s, 'descend') ;

ss {9} = s ;       % hack test

if (opts.tol < 0)
    opts.tol = max (m1,n1) * eps (max (s)) ;
end

stats.rank = sum (s > opts.tol) ;

if (opts.demo)
    subplot (2,2,4) ;
    semilogy (1:length(s), s, 'o', [1 length(s)], [opts.tol opts.tol], 'r-') ;
    title (sprintf ('singular values (rank %d)', stats.rank)) ;
    fprintf ('          numerical rank: %d\n', stats.rank) ;
    drawnow ;
end

%-------------------------------------------------------------------------------
function spy2 (A,t)
%SPY2:  draw a sparse matrix.
% Use CSparse if available, or the MATLAB spy otherwise
try
    cspy (A)
catch
    spy (A)
end
title (t)

function [U, S, V] = piro_band_svd (A, varargin)
%PIRO_BAND_SVD singular value decomposition of a sparse or full matrix.
% A may be sparse or full, and real or complex.  This function can be many
% times faster than SVD if A is banded or can be permuted into a form with a
% small band.
%
% Usage:
%
% s = piro_band_svd (A) ;
%
%   The singular values of A are returned in the vector s.  This does the same
%   the same thing as s = svd (full (A)), but exploiting the band of A.
%
% [U, S, V] = piro_band_svd (A) ;
%
%   Identical to [U,S,V] = svd (full (A)), except that S is returned as sparse.
%
% [U, S, V] = piro_band_svd (A, 'econ') ;
%
%   Identical to [U,S,V] = svd (full (A), 'econ'), except that S is sparse.
%
% [ ... ] = piro_band_svd (A, opts),  where opts is a struct:
%
%       opts.econ: if true (nonzero), this is the same as 'econ', above.
%           The default value is false.
%       opts.qr: if true, perform a QR factorization first.  if false: skip QR.
%           The default is to do the QR.
%       opts.ordering: a string describing the fill-reducing ordering to use.
%           'rcm': band-reducing ordering (SYMRCM).  This is the default.
%           'amd': AMD ordering
%           'colamd': COLAMD ordering
%           'none': no ordering
%
% For most accurate results for rank-deficient matrices, download and install
% SPQR from http://www.cise.ufl.edu/research/sparse .
%
% Example
%
%   load west0479
%   A = west0479 ;
%   s = piro_band_svd (A) ;
%   [U,S,V] = piro_band_svd (A) ;
%
% See also SVD, PIRO_BAND, SYMRCM, AMD, COLAMD, SPQR.

% Copyright 2012, Sivasankaran Rajamanickam, Timothy A. Davis
% http://www.cise.ufl.edu/research/sparse

% Details: For most accurate results for rank-deficient matrices (and when
% opts.qr is true), the SPQR function from SuiteSparse should be used with a
% non-default drop tolerance of zero.  Otherwise, the built-in QR will be used
% instead.  The built-in QR uses SPQR, but it does not allow the drop tolerance
% to be modified, and thus small singular values will not be computed
% accurately.  A warning is issued if this case occurs.  To obtain SPQR and all
% of SuiteSparse, see http://www.cise.ufl.edu/research/sparse .
%
% To set the blocksize used internally in the mexFunction, use
% piro_band_svd(A,opts) where opts.blks is a vector of length 4.  This option
% is meant for development and performance evaluation only.
%
% For testing purposes, opts.qr can be 0: false, 1: use SPQR if available or
% MATLAB QR if not, 2: use MATLAB QR.

%-------------------------------------------------------------------------------
% get input options
%-------------------------------------------------------------------------------

opts = struct ;
if (nargin > 1)
    arg = varargin {1} ;
    if (isstruct (arg))
        opts = arg ;
    else
        opts.econ = (ischar (arg) && strcmp (arg, 'econ')) ;
    end
end
if (~isfield (opts, 'econ'))
    opts.econ = false ;
end
if (~isfield (opts, 'qr'))
    opts.qr = true ;
end
if (~isfield (opts, 'ordering'))
    opts.ordering = 'rcm' ;
end

%-------------------------------------------------------------------------------
% compute the SVD
%-------------------------------------------------------------------------------

[m n] = size (A) ;

if (min (m,n) < 2)

    %---------------------------------------------------------------------------
    % special case:  A is a scalar, vector, or empty matrix.  Use MATLAB svd.
    %---------------------------------------------------------------------------

    A = full (A) ;
    if (nargout > 1)
        if (opts.econ)
            [U,S,V] = svd (A, 'econ') ;
        else
            [U,S,V] = svd (A) ;
        end
        S = sparse (S) ;
    else
        U = svd (A) ;
    end

else

    %---------------------------------------------------------------------------
    % general case: A is a matrix.  Use piro_band_svdkernel
    %---------------------------------------------------------------------------

    %---------------------------------------------------------------------------
    % transpose the matrix if it is short and fat
    %---------------------------------------------------------------------------

    trans = (m < n) ;
    if (trans)
        B = A' ;
    else
        B = A ;
    end
    [m n] = size (B) ;

    %---------------------------------------------------------------------------
    % order the matrix, if requested
    %---------------------------------------------------------------------------

    % find the ordering
    symmetric_ordering = false ;
    switch (opts.ordering)
        case 'rcm'
            S = spones (B) ;
            T = S' ;                % T is n-by-m
            if (opts.qr || m ~= n)
                % symrcm (B'*B), but remove dense rows first
                keep = full (sum (T)) < 16 * sqrt (n) ;
                if (~all (keep))
                    T = T (:,keep) ;
                end
                p = symrcm (T*T') ;
            else
                symmetric_ordering = true ;
                p = symrcm (S+T) ;
            end
        case 'colamd'
            p = colamd (B) ;
        case 'amd'
            if (m == n)
                symmetric_ordering = true ;
                p = amd (B) ;
            else
                p = colamd (B) ;
            end
        case 'none'
            p = [ ] ;
        otherwise
            error ('PIROBAND:unrecognizedOrdering', ...
                'unrecognized ordering option') ;
    end

    % apply the ordering
    if (symmetric_ordering)
        C = B (p,p) ;
    elseif (~isempty (p))
        C = B (:,p) ;
    else
        C = B ;
    end

    %---------------------------------------------------------------------------
    % do a QR factorization first, if requested
    %---------------------------------------------------------------------------

    Q = speye (m) ;
    if (opts.qr)
        if (~issparse (C))
            C = sparse (C) ;
        end
        qropts.tol = 0 ;
        caution = false ;
        if ((opts.qr ~= 2) && (exist ('spqr', 'file') == 3))
            if (nargout > 1)
                [Q,R] = spqr (C, qropts) ;
            else
                R = spqr (C, qropts) ;
            end
        else
            % The built-in qr uses spqr, but it does not take a opts.tol
            % input parameter.  Results may be less accurate if A is
            % rank-deficient.
            caution = true ;
            if (nargout > 1)
                [Q,R] = qr (C) ;
            else
                R = qr (C) ;
            end
        end
        % unsqueeze R if rank deficient
        [p2, rnk] = qr_unsqueeze (R) ;
        if (rnk < n)
            R = R (p2,:) ;
            Q = Q (:,p2) ;
            % now Q*R = C, so p2 is no longer needed
            clear p2
            if (caution)
                warning ('PIRO_BAND_SVD:noSPQR', ...
                    ['Matrix is rank-deficient and SPQR is not available.' ...
                    '  Results may be inaccurate.']) ;
            end
        end
    else
        R = C ;
    end

    % now Q*R = C, where C is equal to B, B(:p), or B(p,p)
    % fprintf ('QR-C %g\n', norm (Q*R-C,1)) ;

    %---------------------------------------------------------------------------
    % compute the SVD of R using the PIRO_BAND SVD kernel
    %---------------------------------------------------------------------------

    if (nargout > 1)

        [U,S,V] = piro_band_svdkernel (R, opts) ;
        % U*S*V' = R
        % (Q*U)*S*V' = Q*R = C
        U = Q*U ;
        % now U*S*V' = C
        % fprintf ('USV''-C: %g\n', norm (U*S*V'-C,1)) ;

        % merge the fill-reducing into U and V
        if (symmetric_ordering)
            % U*S*V' = B(p,p)
            U (p,:) = U ;
            V (p,:) = V ;
            % now U*S*V' = B
            % fprintf ('USV''-B sym: %g\n', norm (U*S*V'-B,1)) ;
        elseif (~isempty (p))
            % U*S*V' = B(:,p)
            V (p,:) = V ;
            % now U*S*V' = B
            % fprintf ('USV''-B unsym: %g\n', norm (U*S*V'-B,1)) ;
        end

        % undo the transpose
        if (trans)
            % U*S*V' = B = A'
            % fprintf ('USV''-A'' trans: %g\n', norm (U*S*V'-A',1)) ;
            S = S' ;
            t = V ;
            V = U ;
            U = t ;
            % now U*S*V' = A
            % fprintf ('USV''-A trans: %g\n', norm (U*S*V'-A,1)) ;
        end

    else

        % s is a vector with just the singular values
        s = piro_band_svdkernel (R, opts) ;
        % return s as the first return value; it happens to be called 'U'
        U = s ;

    end
end


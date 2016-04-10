function [p, rnk] = qr_unsqueeze (R)                                        %#ok
%QR_UNSQUEEZE returns a row permutation that 'unsqueezes' a rank-deficient R.
%
%   [p, r] = qr_unsqueeze (R)
%
% If a sparse matrix A is rank-deficient, the MATLAB statement [Q,R]=qr(A)
% returns a 'squeezed R' matrix that is upper trapezoidal.  This function
% returns a row permutation that moves the squeezed rows so that they are back
% on the diagonal.  The second output is the rank of A as estimated by the
% sparse QR.  R must be a 'squeezed' upper triangular/trapezoidal matrix of
% size m-by-n with m >= n.  As long as A has m >= n, then R=qr(A) or
% [Q,R]=qr(A) will return the required R for this function.
%
% In a squeezed R, a column k is declared 'dead' if it is found to be linearly
% dependent on columns 1:k-1 in A has no corresponding row in R.  For example,
% if the 2nd column is 'dead' then R for a 3-by-3 matrix will look like this:
%
%   1 x x
%   . . 3
%   . . .
%
% where 1 and 3 denote the two 'live' entries on the 'diagonal'.  This
% diagonal has been squeezed, however.  This function moves the live rows
% downwards, placing these entries back on the diagonal:
%
%   1 x x
%   . . .
%   . . 3
%
% Example:
%
%   A = sparse ([1 1 1 ; 1 1 2 ; 1 1 3]) ;
%   [Q,R] = qr (A) ;
%   [p,r] = qr_unsqueeze (R)
%   disp (full (R))
%   disp (full (R (p,:)))
%
% See also qr.

% Copyright 2012, Timothy A. Davis
% http://www.cise.ufl.edu/research/sparse

error ('qr_unsqueeze mexFunction not found') ;

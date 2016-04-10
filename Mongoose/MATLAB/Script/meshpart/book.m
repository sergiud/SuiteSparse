function [A,xy] = book(k,n)
% BOOK : k-page graph
%
% [A,xy] = book(k,n):
% Generate a book of k pages, each with n*2 points.
%
% Outputs:  A is the Laplacian; xy is coordinates for a planar drawing.
%
% John Gilbert, 1994.
% Copyright (c) 1990-1996 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.

A = grid5(n);
A = blockdiags(A,0,k,k);
blocks = 1:(k*n^2);
for i = n^2 : n^2 : (k-1)*n^2
    blocks(i+1 : i+n) = blocks(1:n);
end;
A = contract(A,blocks);

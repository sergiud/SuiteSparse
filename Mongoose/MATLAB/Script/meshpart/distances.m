function D = distances(A,xy)
% DISTANCES : Matrix of distances between adjacent mesh points
%
% Input:  A is the (sparse) adjacency matrix of a directed or undirected graph.
%         xy is the vertex coordinates, one per row, in any dimensional space.
%
% Output: D = distances(A,xy) is the matrix of distances between endpoints
%         of edges.  If A(i,j) is nonzero, then D(i,j) is the Euclidean
%         distance between xy(i,:) and xy(j,:); otherwise D(i,j) = 0.
%
% John Gilbert, 1993.
% Copyright (c) 1990-1996 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.

[m,n] = size(A);
if m~=n, error('A must be square'), end;
[m,k] = size(xy);
if m~=n, error('xy must have n rows'), end;

[i,j] = find(A);
dxy = (xy(i,:)-xy(j,:))';
d = sqrt(sum(dxy .^ 2));

D = sparse(i,j,d,n,n);

function [A,xy] = cockroach(k)
% COCKROACH : Planar graph for which spectral partitioning works poorly.
%
% [A,xy] = cockroach(k):
% Generate a mesh (with 6*k points) whose best edge separator has size 2,
% but for which the spectral algorithm gives a separator of size O(k).
% (From Guattery and Miller, "On the performance of spectral graph
% partitioning methods," SODA 1995.)
%
% Outputs:  A is the Laplacian; xy is coordinates for a planar drawing.
%
% John Gilbert, 1994.
% Copyright (c) 1990-1996 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.

n = 6*k;
A = blockdiags([-1 -1 -1], -1:1, n, n);
B = fliplr(speye(2*k,2*k));
middle = (2*k+1):(4*k);
A(middle,middle) = A(middle,middle)-B;

x = [((3*k-1):-1:0)' ; (0:(3*k-1))'];
y = [zeros(3*k,1) ; ones(3*k,1)];
xy = [x y];



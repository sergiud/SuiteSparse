function [A,xyz] = grid3dt(k)
% GRID3DT : Generate 3-dimensional tetrahedral finite element mesh.
%
% [A,xyz] = GRID3DT(k) returns a k^3-by-k^3 symmetric positive definite 
%        matrix A of the k-by-k-by-k grid, with cells divided into tetrahedra, 
%        and an array xyz of coordinates for the grid points.
%
% John Gilbert, Xerox PARC, 1992.
% Copyright (c) 1990-1996 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.

% 1-d mesh
a = blockdiags ([-1 14 -1], -1:1, k, k);

% glue to laminate 1-d meshes into 2-d meshes
b = blockdiags ([-1 -1], [0 1], k, k);

% 2-d mesh
aa = blockdiags ([b a b'], -1:1, k, k);

% glue to laminate 2-d meshes into 3-d meshes
bb = blockdiags ([b b'], [0 1], k, k);

% 3-d mesh
A = blockdiags ([bb aa bb'], -1:1, k, k);

% xyz coordinates of nodes
xyz = zeros(k^3,3);
x = ones(k,1) * (1:k);
y = x';
j = 1;
for i = 1:k
    slab = j:j+k^2-1;
    xyz(slab,1) = x(:);
    xyz(slab,2) = y(:);
    xyz(slab,3) = i * ones(k^2,1);
    j = j + k^2;
end;

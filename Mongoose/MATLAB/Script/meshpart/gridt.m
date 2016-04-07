function [A,xy] = gridt(k)
% GRIDT : Generate triangular mesh.
%
% [A,xy] = GRIDT(k) returns a k*(k+1)/2-square symmetric positive 
% definite matrix A with the structure of an equilateral triangular 
% mesh with side k, and an array xy of coordinates of the mesh points.
%
% John Gilbert, Xerox PARC, 1992.
% Copyright (c) 1990-1996 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.

% Start with a square 7-point grid.

[A,xy] = grid7(k);

% Mask off one triangle of it.

xy = xy-1;
f = find(xy(:,1)+xy(:,2)<k);
A = A(f,f);
xy = xy(f,:);

% Make the other triangle equilateral.

T = [ 1 0 ; 1/2 2/sqrt(5)];
xy = xy*T;

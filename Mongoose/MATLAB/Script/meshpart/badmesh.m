function [A,xy] = badmesh(k,alpha)
% BADMESH : Mesh that can't be partitioned well with a straight line.
%
% [A,xy] = badmesh(k,alpha):
% Generate a k-level mesh (with 4*k points) with no straight-line separator,
% with alpha<1 the ratio between shell sizes.  (alpha defaults to 4/5).
% alpha=0 means linearly spaced shells.
%
% See also GEOPART.
%
%       John Gilbert and Shanghua Teng, 1992-1993.
% Copyright (c) 1990-1996 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.


if nargin < 2
    alpha = 4/5;
end;

if alpha
    x = alpha .^ (1:k);
else
    x = [1:k];
end

x = x';

xy = [x x ; -x x ; -x -x ; x -x];

square = [1 1 0 1 ; 1 1 1 0 ; 0 1 1 1 ; 1 0 1 1];
I = speye(4,4);
A = blockdiags ([I square I], [-1:1], k, k);

kk = 0:k-1;
shuffle = [4*kk 4*kk+1 4*kk+2 4*kk+3] + ones(1,4*k);

A = A(shuffle,shuffle);

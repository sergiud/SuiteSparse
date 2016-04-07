function p = coordpart(A,xy,picture);
% COORDPART : Coordinate bisection partition of a mesh.
%
% p = coordpart(A,xy) returns a list of the vertices on one side of a partition
%     obtained by bisection perpendicular to a coordinate axis.  We try every
%     coordinate axis and return the best cut.
%     Input A is the adjacency matrix of the mesh; 
%     each row of xy is the coordinates of a point in d-space.
%
% coordpart(A,xy,1) also draws a picture.
%
% See also GEOPART, SPECPART, GSPART
%
% John Gilbert  24 May 1994
% Copyright (c) 1990-1996 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.

if nargin < 3
    picture = (nargout == 0);
end;

d = size(xy,2);
best_cut = inf;
for dim = 1:d
    v = zeros(d,1);
    v(dim) = 1;
    pd = partition(xy,v);
    this_cut = cutsize(A,pd);
    if this_cut < best_cut
        best_cut = this_cut;
        p = pd;
    end;
end;

if picture
    clf reset
    gplotpart(A,xy,p);
    title('Coordinate Bisection')
end;

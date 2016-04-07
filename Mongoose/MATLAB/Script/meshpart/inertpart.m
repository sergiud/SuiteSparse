function p = inertpart(A,xy,picture);
% INERTPART : Inertial partition of a graph.
%
% p = inertpart(A,xy) returns a list of the vertices on one side of a partition
%     obtained by bisection with a line or plane normal to a moment of inertia
%     of the vertices, considered as points in Euclidean space.
%     Input A is the adjacency matrix of the mesh (used only for the picture!);
%     each row of xy is the coordinates of a point in d-space.
%
% inertpart(A,xy,1) also draws a picture.
%
% See also GEOPART, SPECPART, ISPART
%
% John Gilbert  11 Nov 1994
% Copyright (c) 1990-1996 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.

if nargin < 3
    picture = (nargout == 0);
end;
[n,d] = size(xy);

% We shift the points to have mean zero,
% and then cut normal to the singular vector.

xymean = mean(xy);
[U,S,V] = svd(xy - ones(n,1)*xymean,0);
v = V(:,1)';
p = partition(xy,v);

if picture
    clf reset
    gplotpart(A,xy,p);
    title('Inertial Partition')
end;


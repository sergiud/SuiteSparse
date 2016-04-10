function xyz = conmapinv(c,xyzmap)
% CONMAPINV : Inverse conformal map for geometric partitioning.
%
% [xyz,xy] = conmapinv(c,xyzmap):
% Input is points xyzmap on a unit sphere in d-space
% (one point per row of matrix xyz),
% and a point c in d-space represented as a row vector.
%
% We compute the inverse of the conformal mapping in conmap.m.
%
% The output is the new points on the sphere.
%
% See GEOPART.
%
% John Gilbert and Shanghua Teng, 1992-1993.
% Copyright (c) 1990-1996 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.


% Compute the reflection and stretch.
[Q,r] = reflector(c);
alpha = sqrt((1+r)/(1-r));

% Unstretch.
xymap = stereodown(xyzmap);
xyref = alpha * xymap;
xyzref = stereoup(xyref);

% Unreflect.
xyz = xyzref * Q';

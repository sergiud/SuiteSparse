function xy = stereodown(xyz)
% STEREODOWN : Stereographic projection from sphere to plane.
%
% xy = stereodown(xyz):  
% Project points xyz stereographically from the unit sphere in dim dimensions
% to the dim-1 dimensional space.
% The number of dimensions of the mesh is the number of columns of xyz minus 1.
% This is the inverse of stereoup.
%
% See GEOPART.
%
% John Gilbert and Shanghua Teng, 1992-1993.
% Copyright (c) 1990-1996 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.

[n,dim] = size(xyz);   % dim is the number of dimensions of the sphere.

xy = xyz(:,1:dim-1) ./ ((ones(n,1)-xyz(:,dim)) * ones(1,dim-1));

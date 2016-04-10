function xyz = stereoup(xy)
% STEREOUP : Stereographic projection from plane to sphere.
%
% xyz = stereoup(xy):  
% Project points xy stereographically from the dim-1 dimensional space 
% onto the unit sphere in dim dimensions.
% The number of dimensions of the mesh is the number of columns of xy.
% This is the inverse of stereodown.
%
% See GEOPART.
%
% John Gilbert and Shanghua Teng, 1992-1993.
% Copyright (c) 1990-1996 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.

[n,dim] = size(xy);
dim = dim+1;       % dim is the number of dimensions of the sphere.

edim = [zeros(n,dim-1) ones(n,1)];
normsquared = [(xy.^2) ones(n,1)] * ones(dim,dim);

% Each row of "edim" is the dim-dimensional unit vector [0 0 ... 0 1].
% "normsquared" is the square of the norm of (xy - edim), if we think
% of xy as embedded in dim dimensions.
%
% Thus xyz = edim + 2*(xy - edim) ./ normsquared.
 
xyz = edim + 2 * [xy -ones(n,1)] ./ normsquared;




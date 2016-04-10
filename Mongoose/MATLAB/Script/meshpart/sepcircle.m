function [gc,gcquality] = sepcircle(A,xyz,ntries)
% SEPCIRCLE : Good great circle for geometric partitioning.
%
% [gc,gcquality] = sepcircle(A,xyz,ntries)
% Generate a great circle on the unit sphere in d dimensions,
% trying to find one that gives a good edge separator of the graph A.
% "xyz" is the input points, conformally mapped on the unit sphere
%       in d-space so that the approximate centerpoint is at the origin.
%
% We return the best of ntries (default 40) attempts.
%
% The circle is returned only as a direction, a d-vector normal to its plane.
% The second output is the partition quality as measured by sepquality.
%
% See GEOPART.
%
% John Gilbert and Shanghua Teng, 1992-1993.
% Copyright (c) 1990-1996 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.

if nargin < 3
    ntries = 40;
end;

[npoints,d] = size(xyz);

% For inertial weighting, we weight the randomly chosen great circles
% according to a power (say 2) of the inertial matrix of the points.

M = (xyz' * xyz) ^ 2;

% Now choose all ntries random great circles.  
% Each circle is represented by a vector normal to its plane.
% Normally distributed points give vectors with uniformly distributed 
% directions, which we then weight by the matrix M from above.

vv = randn(ntries,d) * M;

quality = Inf * ones(ntries,1);
for i = 1 : size(vv,1)
    v = vv(i,:);
    if norm(v) ~= 0
        quality(i) = sepquality(v,A,xyz);
    end;
end;

[gcquality,i] = min(quality);
gc = vv(i,:)';

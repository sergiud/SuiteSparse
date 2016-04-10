function c = centerpoint(xyz,n)
% CENTERPOINT : Approximate centerpoint in any number of dimensions.
%
% c = centerpoint(xyz,n):
% Input: xyz: matrix of points, one per row.
%        n: sample size (optional, defaults to number of points, 
%                        rounded to a multiple of d+2 if necessary).
%
% See GEOPART.
%
% John Gilbert and Shanghua Teng, 1992-1994.
% Copyright (c) 1990-1996 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.
%
% Modified by Tim Davis, for Matlab 5.1.  July 6, 1998
% Modified by John Gilbert for Matlab 6, Feb 2002

% Determine sample size n, and round it to be congruent to 1 mod d+1.
[npoints,d] = size(xyz);
if nargin < 2
    n = npoints;
else
    n = min(n,npoints);
end;
n = (d+1) * floor((n-1)/(d+1)) + 1;

% Sample points without replacement.

xyzsize = ceil(n*(1+1/(d+1)));
xyzs = zeros(xyzsize,d);
sample = randperm(npoints);
sample = sample(1:n);
xyzs(1:n,:) = xyz(sample,:);

% Perform Radon reduction according to a full d+2-ary tree.

queuehead = 1;
queuetail = n;
while queuehead < queuetail
    queuetail = queuetail+1;
    xyzs(queuetail,:) = radong(xyzs(queuehead:queuehead+d+1,:));
    queuehead = queuehead+d+2;
end;

c = xyzs(queuetail,:);

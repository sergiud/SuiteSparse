function [A,xyz] = tree3xpath(dep,len);
% TREE3XPATH : Graph for which spectral partitioning works poorly.
%
%   Generate a triple-tree version of the tree-cross-path mesh 
%   from Guattery and Miller,  "On the performance of spectral 
%   graph partitioning methods," SODA 1995.
%
% [A,xyz] = tree3xpath(n)
%   The generated mesh has about n vertices.
% [A,xyz] = tree3xpath(dep,len)
%   The tree has depth "dep" and the path has length "len",
%   so the whole mesh has (3*2^dep - 2)*len vertices.
%
% Outputs:  A is the Laplacian; xyz is coordinates for a drawing.
%           If no output arguments, we draw the picture.
%
% See also TREEXPATH.
%
% John Gilbert, 1994.
% Copyright (c) 1990-1996 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.

if nargin <= 1
    n = dep;
    dep = floor(2/3 * log(n)/log(2));
end;
ntree = 3*2^dep - 2;
if nargin <= 1
    len = round(n/ntree);
end;

% The tree is three complete binary trees hanging off a common center.

if dep == 0
    parent = [1];
else
    p = ntree/2 - 1;
    parent = [0:p ; 0:p];
    parent = parent(:);
    parent(1) = 1;
    parent(2) = 1;
end;

Tree = sparse(1:ntree,parent,1,ntree,ntree);
Tree = Tree|Tree'|speye(ntree);

xy = zeros(ntree,2);
firstpoint = 2;
radius = 1;
firstangle = 0;
levelsize = 3;
for d = 1:dep
    range = (0:levelsize-1)'+firstpoint;
    theta = ((0:levelsize-1)' * 2 * pi / levelsize) + firstangle;
    xy(range,1) = radius * cos(theta);
    xy(range,2) = radius * sin(theta);
    firstangle = firstangle-pi/2/levelsize;
    firstpoint = firstpoint+levelsize;
    levelsize = 2*levelsize;
    alpha = pi/levelsize;
    radius = radius * (cos(alpha) + sin(alpha)/tan(pi/3-alpha));
end;

A = blockdiags([Tree Tree Tree], -1:1, len, len);
A = laplacian(A);

index = ones(ntree,1) * (1:len);
zindex = index(:);
index = (1:ntree)' * ones(1,len);
xyindex = index(:);
xyz = [xy(xyindex,:) zindex];

if nargout == 0
    if len == 1
        gplotg(A,xy);
    else
        gplotg(A,xyz);
    end;
end;

function [A,xyz] = treexpath(dep,len);
% TREEXPATH : Graph for which spectral partitioning works poorly.
%
%   Generate the tree-cross-path mesh from Guattery and Miller, 
%   "On the performance of spectral graph partitioning methods," SODA 1995.
%
% [A,xyz] = treexpath(n)
%   The generated mesh has about n vertices.
% [A,xyz] = treexpath(dep,len)
%   Each tree has depth "dep" and the path has length "len",
%   so the whole mesh has (2^(dep+2) - 2)*len vertices.
% Outputs:  A is the Laplacian; xyz is coordinates for a drawing.
%           If no output arguments, we draw the picture.
%
% Try [A,xyz] = treexpath(700), and then use various partitioners on the result.
%
% John Gilbert, 1994.
% Copyright (c) 1990-1996 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.

if nargin <= 1
    n = dep;
    ntree = (n/4)^(2/3);
    dep = floor(log(ntree+2)/log(2))-2;
    dep = max(dep,0);
end;
ntree = 2^(dep+2) - 2;
if nargin <= 1
    len = round(n/ntree);
end;
fprintf('TreeDepth=%d, TreeSize=%d, PathSize=%d, TotalSize=%d\n', ...
         dep, ntree, len, len*ntree);

% The tree is two complete binary trees joined by a edge.

p = ntree/2 - 1;
parent = [0:p ; 0:p];
parent = parent(:);
parent(1:2) = [1;1];

Tree = sparse(1:ntree,parent,1,ntree,ntree);
Tree = Tree|Tree'|speye(ntree);

xy = zeros(ntree,2);
firstpoint = 1;
radius = .5;
firstangle = 0;
levelsize = 2;
for d = 0:dep
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

% The whole graph is the cross product of the tree and a path.

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

function [x,y,h,s] = treelayoutg(parent,post)
% TREELAYOUTG : Lay out a tree or forest.
%	[x,y,h,s] = treelayoutg(parent,post)
%	    parent is the vector of parent pointers, with 0 for a root.
%	    post is a postorder permutation on the the tree nodes.
%	    (If post is omitted we compute it here.)
%	    x and y are vectors of coordinates in the unit square at which 
%	    to lay out the nodes of the tree to make a nice picture.
%	    Optionally, h is the height of the tree and s is the 
%	    number of vertices in the top-level separator.
%
% See also TREEPERM, ETREE, TREEPLOTG, ETREEPLOTG.
% Copyright (c) 1990-1996 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.
% Copyright (c) 1984-92 by The MathWorks, Inc.
%
% Written 1991 by John Gilbert.
% Modified 21 May 1993 by JRG to use depth not height for y.
% Modified 23 Jun 1993 by JRG to call treeperm.
% Modified  4 Aug 1994 by JRG to return height in vertices.
% Modified  6 Feb 2002 by JRG for Matlab 6.

% This is based on the C code in sptrees.c by John Gilbert.
% Leaves are spaced evenly on the x axis, and internal
% nodes are centered over their descendant leaves with
% y coordinate proportional to 1 - depth in the tree.

if nargin < 2,
    post = treeperm(parent);
end;

n = max(size(parent));

% Add a dummy root node #n+1, and identify the leaves.

parent = rem (parent+n, n+1) + 1;  % change all 0s to n+1s
isaleaf = ones(1,n+1);
isaleaf(parent) = zeros(n,1);

% In postorder, compute heights and descendant leaf intervals.
% Space leaves evenly in x (in postorder).

xmin = n*ones(1,n+1);
xmax = zeros(1,n+1);
height = zeros(1,n+1);
depth = zeros(1,n+1);
nkids = zeros(1,n+1);
nleaves = 0;

for i = 1:n,
    node = post(i);
    if isaleaf(node),
        nleaves = nleaves+1;
        xmin(node) = nleaves;
        xmax(node) = nleaves;
    end;
    dad = parent(node);
    height(dad) = max (height(dad), height(node)+1);
    xmin(dad)   = min (xmin(dad),   xmin(node));
    xmax(dad)   = max (xmax(dad),   xmax(node));
    nkids(dad)  = nkids(dad)+1;
end;

% In reverse postorder, compute depths.
 
for i = n:-1:1
    node = post(i);
    dad = parent(node);
    depth(node) = depth(dad)+1;
end;

% Compute coordinates, leaving a little space on all sides.

treeht = height(n+1) - 1;
deltax = 1/(nleaves+1);
deltay = 1/(treeht+2);
x = deltax * (xmin+xmax)/2;
y = deltay * (treeht-depth+2);

% Omit the dummy node.

x = x(1:n);
y = y(1:n);

% Return the height and top separator size.
% (We return the height in vertices, standard sparse matrix practice.)

h = treeht+1;
s = n+1 - max(find(nkids~=1));

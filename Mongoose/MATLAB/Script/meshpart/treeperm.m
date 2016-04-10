function post = treeperm(parent)
% TREEPERM :  Compute postorder permutation for a tree or forest.
%	post = treeperm(parent)
%	    parent is the vector of parent pointers, with 0 for a root.
%	    post is a postorder permutation on the the tree nodes.
%
%	See also ETREE, TREELAYOUTG, TREEPLOTG, ETREEPLOTG.
%
% Written 23 Jun 1993 by John Gilbert.
% Copyright (c) 1990-1996 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.


% Create the adjacency matrix A of the given tree,
% and get the postorder with a call to etree.

[nr,n] = size(parent);
if nr > 1
    parent = parent';
    n = nr;
end;
kids = find(parent);
if size(parent,1) ~= 1 | any(parent(kids) <= kids) | any(parent(kids) > n)
    error('Input must be a tree in topological order');
end;
i = [kids parent(kids) 1:n];
j = [parent(kids) kids 1:n];
A = sparse (i, j, 1, n, n);
[ignore, post] = etree(A);

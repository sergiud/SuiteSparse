function [R,xyR] = contract(A,blocks,xyA)
% CONTRACT : Condense a graph according to a given block structure.
%
% R = contract(A,blocks):
%
% Input:   A is an n by n matrix representing a directed or undirected graph.
%          blocks is an n-vector of vertex labels.
%
% Output:  R is the k-vertex "contracted graph" obtained by contracting 
%          each block of A (with like-labeled vertices) into a single node.  
%          Edge weights in R are the number of edges joining the blocks,
%          and vertex weights are the number of vertices in the block.
%
% If a third input argument "xyA" is present, a second output argument "xyR"
% is returned, with xy coordinates for the vertices of R:
%
%                 [R,xyR] = contract(A,blocks,xyA)
%
% John Gilbert and Shanghua Teng, 1992-1993.
% Copyright (c) 1990-1996 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.

[n,m] = size(A);
if m ~= n
    error('A should be square');
end;
if max(size(blocks)) ~= n | min(size(blocks)) ~= 1
    error('Vertex labels should be an n-vector');
end;

A = spones(A) | speye(n);

% Normalize "blocks" to be contiguous integers 1:k.
blocks = ranks(blocks);
k = max(blocks);

sizes = full(sparse(blocks,1,1,k,1));  % block sizes.

% "blocks" maps vertices of A to vertices of R.
% Map each edge of A to an edge of R.

[i,j] = find(A);
R = sparse(blocks(i),blocks(j),1,k,k);
R = R - diag(diag(R)) + diag(sparse(sizes));

% If coordinates for A were given, get coords for each vertex
% of R by averaging the coords of the corresponding block of A.

if nargin >= 3
    xR = full(sparse(blocks,1,xyA(:,1),k,1)) ./ sizes;
    yR = full(sparse(blocks,1,xyA(:,2),k,1)) ./ sizes;
    xyR = [xR yR];
end;

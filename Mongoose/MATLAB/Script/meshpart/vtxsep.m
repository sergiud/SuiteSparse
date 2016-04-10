function [sep,parta,partb] = vtxsep(A,a,b)
% VTXSEP : Convert an edge separator (or node partition) to a vertex separator.
%
% sep = vtxsep(A,a,b):
% A is a symmetric 0-1 matrix representing an undirected graph.
% a and b partition the vertices of A.  (b is optional.)
% This function returns a vertex separator sep of minimum size
% such that a-sep and b-sep are in different components of A-sep.
%
% Optional outputs: [sep,parta,partb], where sep is as above
% and parta = a - sep, partb = b - sep.
%
% John Gilbert, Xerox PARC, 1993
% Copyright (c) 1990-1996 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.
 
% We use dmperm to find a maximum matching on the bipartite subgraph 
% of A induced by the edges joining a to b, and then use the matching
% to construct a minimum vertex cover of that subgraph, which is the
% desired vertex separator.

% aborder is the points of a adjacent to b, and similarly bborder

A = spones(A);
if nargin < 3
    b = other(a,A);
end;

% The following "if"s are because "max" behaves differently
% on a 1-row matrix than on a multirow matrix.
if length(b) > 1
    aborder = a(find(max(A(b,a))));
else
    aborder = a(find(A(b,a)));
end;
if length(a) > 1
    bborder = b(find(max(A(a,b))));
else
    bborder = b(find(A(a,b)));
end;

if length(aborder) == 0

    % The parts are disconnected, so the separator is empty.
    sepa = [];
    sepb = [];

else

    % Use dmperm to find a matching of the bipartite graph.
    % The separator is points of a in the horizontal subgraph,
    % plus points of b in the vertical subgraph.

    [p,q,r,s] = dmperm(A(aborder,bborder));
    sepa = aborder(p(r(1):(r(2)-1)));
    sepb = bborder(q(s(2):(s(length(s))-1)));

end;

sep = [sepa sepb];

t = zeros(1,length(A));
t(a) = ones(size(a));
t(sepa) = zeros(size(sepa));
parta = find(t);

t = zeros(1,length(A));
t(b) = ones(size(b));
t(sepb) = zeros(size(sepb));
partb = find(t);

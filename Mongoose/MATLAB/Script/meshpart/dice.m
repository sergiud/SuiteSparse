function [map,sepij,sepA] = dice(method,levels,A,arg2,arg3,arg4)
% DICE : Separate a graph recursively.
%
% [map,sepij,sepA] = dice(method,levels,A,arg2,arg3,arg4)
% partitions the mesh or graph A recursively by a specified method.
%
% Input:
%   'method' is the name of the 2-way edge separator function to call.
%   levels   is the number of levels of partitioning.
%   A        is the adjacency matrix of the graph.
%   arg2, arg3, arg4 are optional additional arguments to the 2-way function.
%            arg2 is special : If it has the same number of rows as A,
%            then the recursive calls use the appropriate subset of arg2
%            as well as of A.  This is useful if the partitioner uses xy coords.
%
% Output:
%   map      is a vector of integers from 0 to 2^levels-1, indexed by vertex,
%            giving the partition number for each vertex.
%   sepij    is a list of separating edges; each row is an edge [i j].
%   sepA     is the adjacency matrix of the graph less the separating edges.
%
% For example,  dice('geopart',7,A,xy,0,ntries) 
%   splits A into 2^7=128 parts, using a call at each stage that looks like
%   geopart(A,xy,0,ntries).
% 
% John Gilbert, 1994.
% Copyright (c) 1990-1996 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.

minpoints = 8;   % Don't separate pieces smaller than this.

nargcall = nargin - 2;

n = size(A,1);

if n < minpoints | levels < 1
    map = zeros(1,n);
else
    % Call the partitioner to split the graph, giving one part.
    if nargcall == 4
        a = feval(method, A, arg2, arg3, arg4);
    elseif nargcall == 3
        a = feval(method, A, arg2, arg3);
    elseif nargcall == 2
        a = feval(method, A, arg2);
    else
        a = feval(method, A);
    end;

    % Get the complement of that part, giving the other part.
    b = other(a,n);

    % Call recursively.
    if nargcall >= 2
        if size(arg2,1) == n
            arg2a = arg2(a,:);
            arg2b = arg2(b,:);
        else
            arg2a = arg2;
            arg2b = arg2;
        end;
    end;
    if nargcall == 4
        mapa = dice(method,levels-1,A(a,a),arg2a,arg3,arg4);
        mapb = dice(method,levels-1,A(b,b),arg2b,arg3,arg4);
    elseif nargcall == 3
        mapa = dice(method,levels-1,A(a,a),arg2a,arg3);
        mapb = dice(method,levels-1,A(b,b),arg2b,arg3);
    elseif nargcall == 2
        mapa = dice(method,levels-1,A(a,a),arg2a);
        mapb = dice(method,levels-1,A(b,b),arg2b);
    else
        mapa = dice(method,levels-1,A(a,a));
        mapb = dice(method,levels-1,A(b,b));
    end;

    % Set up the whole map.
    map = zeros(1,n);
    mapb = mapb + max(mapa) + 1;
    map(a) = mapa;
    map(b) = mapb;
end;

% Set up the separating edge list and separated graph.
if nargout >= 2
    [i,j] = find(A);
    f = find(map(i) > map(j));
    sepij = [i(f) j(f)];
    if nargout >= 3
        f = find(map(i) == map(j));
        sepA = sparse(i(f), j(f), 1, n, n);
    end;
end;

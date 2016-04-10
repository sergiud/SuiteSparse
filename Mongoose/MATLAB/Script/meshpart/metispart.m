function [p1,p2] = metispart(A,xy);
% METISPART : Partition a graph using Metis default method.
%
% p = metispart(A) returns a list of the vertices on one side of
%     a partition obtained by Metis 4.0 applied to the graph of A.
%     
% Optional arguments:
%   metispart(A,xy) draws a picture of the partitioned graph,
%                   using the rows of xy as vertex coordinates.
%   [p1,p2] = metispart(...)   also returns the list of vertices
%                              on the other side of the partition.
%
% See also METISMEX (which accepts all the Metis options), 
% GEOPART, GSPART, SPECPART, METISDICE, METISND.
%
% John Gilbert  3 Jul 01
% Copyright (c) 1990-2001 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.

if nargin < 2
    xy = 0;
end;
picture = max(size(xy)) > 1;

map = metismex('PartGraphRecursive',A,2);
[p1,p2] = other(map);

if picture
    gplotpart(A,xy,p1);
    title('Metis Partition')
end;
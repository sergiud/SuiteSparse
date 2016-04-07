function map = specdice(A,a,b);
% SPECDICE : Spectral multiway partition.
%
% map = specdice(A,nlevels)
% A is the adjacency matrix of a graph.
% This uses spectral partitioning to divide A into 2^nlevels pieces
% of equal size (within one node), with relatively small connections.
%
% If xy coordinates are given as a second or third argument,
% specdice draws a picture of the result.
%
% See also SPECPART, DICE, GPLOTMAP, GSDICE, GEODICE.
%
% John Gilbert, 1994.
% Copyright (c) 1990-1996 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.

if nargin >= 3
    if length(a) == 1
        nlevels = a; xy = b; picture = 1;
    else
        nlevels = b; xy = a; picture = 1;
    end;
elseif nargin == 2
    if length(a) == 1
        nlevels = a; picture = 0;
    else
        nlevels = 4; xy = a; picture = 1;
    end;
else
        nlevels = 4; picture = 0;
end;

map = dice('specpart',nlevels,A);

if picture
    gplotmap(A,xy,map);
    title('Spectral Partition');
end;

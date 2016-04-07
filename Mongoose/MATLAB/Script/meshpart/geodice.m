function map = geodice(A,a,b,c);
% GEODICE : Geometric multiway partition.
%
% map = geodice(A,xy,nlevels,ntries)
% A is the adjacency matrix of a graph.
% xy is its vertex coordinates, one vertex per row.
% This uses geometric partitioning to divide A into 2^nlevels pieces
% of equal size (within one node), with relatively small connections.
% ntries (optional, default 30) is the number of trial circles.
%
% If no output argument is given, geodice draws a picture of the result.
%
% See also GEOPART, DICE, GPLOTMAP, SPECDICE, GSDICE.
%
% John Gilbert, 1994.
% Copyright (c) 1990-1996 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.

picture = (nargout == 0);

% Sort out the ntries.
ntries = 30;
if nargin == 4
    ntries = c;
elseif nargin == 3
    if length(a) == 1 & length(b) == 1
        error('No coordinates given?');
    end;
end;

% Sort out the coordinates and nlevels.
if nargin >= 3
    if length(a) == 1
        nlevels = a;
        if length(b) > 1
            xy = b; 
        else
            error('No coordinates given?');
        end;
    else
        nlevels = b; xy = a; 
    end;
elseif nargin == 2
    if length(a) == 1
        error('No coordinates given?');
    else
        nlevels = 4; xy = a; 
    end;
else
        error('No coordinates given?');
end;

% Enough with the arguments!

map = dice('geopart',nlevels,A,xy,0,ntries);

if picture
    gplotmap(A,xy,map);
    title('Geometric Partition');
end;

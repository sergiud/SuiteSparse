function map = gsdice(A,a,b,c);
% GSDICE : Geometric spectral multiway partition.
%
% map = gsdice(A,nlevels)
% A is the adjacency matrix of a graph.
% This uses geometric spectral partitioning to divide A into 2^nlevels 
% pieces of equal size (within one node), with relatively small connections.
%
%       gsdice(A,nlevels,xy) or
%       gsdice(A,xy,nlevels):  Draw a picture of the result, as well.
%       gsdice( ... ,ntries):  Use "ntries" trials in the geometric routine.
%
% See also GSPART, DICE, GPLOTMAP, SPECDICE, GEODICE.
%
% John Gilbert, 1994.
% Copyright (c) 1990-1996 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.

% Sort out the ntries
ntries = 50;
if nargin == 4
    ntries = c;
elseif nargin == 3
    if length(a) == 1 & length(b) == 1
        ntries = b;
    end;
end;

% Sort out the coordinates and nlevels.
if nargin >= 3
    if length(a) == 1
        nlevels = a;
        if length(b) > 1
            xy = b; picture = 1;
        else
            picture = 0;
        end;
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

% Enough with the arguments!

map = dice('gspart',nlevels,A,0,ntries);

if picture
    gplotmap(A,xy,map);
    title('Geometric Spectral Partition');
end;

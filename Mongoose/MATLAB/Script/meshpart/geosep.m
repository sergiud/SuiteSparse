function [sep,a,b] = geosep(A,xy,picture,ntries)
% GEOSEP : 2-D geometric separator for a finite element mesh.
%
% s = geosep(A,xy) returns the geometric separator for
% the mesh whose structure is A and whose coordinates are xy.
%
% The number of dimensions is the number of columns of xy.
%
% [s,a,b] = geosep(...) also returns the partition [a,b] of A-s.
%
% There are two optional arguments:
% geosep(A,xy,picture,ntries)
%    picture (default 0) : Draw a picture (2 dimensions only).
%    ntries (default 30) : Number of random choices to make.
%
% See also GEOPART, GEODICE, GEOND, SHOWGEOSEP, SHOWDICE.
%
% John Gilbert and Shanghua Teng, 1992-1993.
% Copyright (c) 1990-1996 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.

if nargin < 3 | size(xy,2) > 2
    picture = 0;
end;
if nargin < 4
    ntries = 30;
end;

[parta,partb] = geopart(A,xy,-abs(picture),ntries);
[sep,a,b] = vtxsep(A,parta,partb);

if picture
    na_nb_nc = [length(a) length(b) length(sep)]
    highlight(A,xy,sep);
    axis('off');
    drawnow;
end;

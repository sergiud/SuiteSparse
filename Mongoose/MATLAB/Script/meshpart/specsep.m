function [sep,a,b] = specsep(A,xy)
% SPECSEP : Spectral separators for a finite element mesh.
%
% s = specsep(A) returns a separator based on the Fiedler vector
% for the graph whose structure is A.
%
% If coords are given as a second argument 'xy', it draws a picture.
%
% [s,a,b] = specsep(...) also returns the partition [a,b] of A-s.
%
% John Gilbert 11 Jun 93.
% Copyright (c) 1990-1996 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.

if nargin < 2
  xy = 0;
end;

picture = (max(size(xy)) >= 2);

[parta,partb] = specpart(A);
[sep,a,b] = vtxsep(A,parta,partb);

if picture
    na_nb_nc = [length(a) length(b) length(sep)]
    highlight(A,xy,sep);
    axis('off');
    drawnow;
end;


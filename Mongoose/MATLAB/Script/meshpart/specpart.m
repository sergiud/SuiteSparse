function [part1,part2,sep1,sep2] = specpart(A,xy,ignore);
% SPECPART : Spectral partition of a graph.
%
% [part1,part2] = specpart(A) returns a partition of the n vertices
%                 of A into two lists part1 and part2 according to the
%                 spectral bisection algorithm of Simon et al:  
%                 Label the vertices with the components of the Fiedler vector
%                 (the second eigenvector of the Laplacian matrix) and partition
%                 them about the median value.
%
% [part1,part2,sep1,sep2] = specpart(.) also returns the separating edges.
%
% If vertex coordinates are given as a second argument,
% specpart(A,xy) draws a picture of the result.
%
% See also LAPLACIAN, FIEDLER.
%
% John Gilbert  9 March 1993
% Copyright (c) 1990-1996 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.

if nargin < 2
    xy = 0;
end;
picture = max(size(xy)) > 1;

x = fiedler(A);
x = (x(:))';
mid = median(x);
part1 = find(x < mid);
part2 = find(x >= mid);

if nargout > 2 | picture
    [sep1,sep2] = find(A(part1,part2));
    sep1 = part1(sep1);
    sep2 = part2(sep2);
    if picture
        na_nb_ne = [length(part1) length(part2) length(sep1)];
        gplotpart(A,xy,part1);
        title('Spectral Partition');
    end;
end;

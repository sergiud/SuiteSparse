function [out1,out2] = other(in1,in2);
% OTHER : Find the other part of a partition, or
%         convert a partition to the other representation.
%
% part2 = other(part1,n)     : part1 is a list of subscripts in 1:n;
%                              part2 becomes the complementary list.
% part2 = other(part1,A)     : Same, except n is taken as max(size(A)).
% [part2,p] = other(part1,n) : Also returns 0/1 partition vector of length n.
% [part2,p] = other(part1,A) : Same.
% [part1,part2] = other(p)   : Converts 0/1 vector to lists of subscripts.
%
% John Gilbert, 1993.
% Copyright (c) 1990-1996 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.

% Modified 3 Jul 01 by JRG to fix a bug in the "map" case

if nargin == 1

    out1 = find(in1 == 0);
    out2 = find(in1 ~= 0);

else

    if max(size(in2)) > 1
        in2 = max(size(in2));
    end;
    out2 = full(spones(sparse(1,in1,1,1,in2)));
    out1 = find(out2 == 0);

end;

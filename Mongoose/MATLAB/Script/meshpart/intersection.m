function C = intersection(A,B)
% INTERSECTION : Set intersection
%
% C = intersection(A,B) returns a row vector with one occurrence
%     of each element common to A and B, in sorted order.
%
% Based on code by Asko Huuskonen, Finnish Meteorological Institute,
% Asko.Huuskonen@fmi.fi, posted to comp.soft-sys.matlab.
% This version by John Gilbert, Xerox PARC, 27 Sep 1993
% Copyright (c) 1990-1996 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.

if nargin < 2
    B = A;
end;
if length(A) == 0
    C = [];
    return;
end;
if length(B) == 0 
    C = [];
    return;
end;

A = setfilter(A);
B = setfilter(B);
AB = sort([A B]);
C = AB(find(diff(AB)==0));

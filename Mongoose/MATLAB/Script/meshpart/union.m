function C = union(A,B)
% UNION : Set union
%
% C = union(A,B) returns a row vector with one occurrence
%     of each different element of A and B, in sorted order.
%
% Based on code by Asko Huuskonen, Finnish Meteorological Institute,
% Asko.Huuskonen@fmi.fi, posted to comp.soft-sys.matlab.
% This version by John Gilbert, Xerox PARC, 27 Sep 1993
% Copyright (c) 1990-1996 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.

if nargin < 2
    B = [];
end;
if length(A) == 0
    C = setfilter(B);
    return;
end;
if length(B) == 0 
    C = setfilter(A);
    return;
end;

C = setfilter([A(:);B(:)]);

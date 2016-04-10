function B = setfilter(A)
% SETFILTER : Sort and remove duplicates.
%
% B = setfilter(A) returns a row vector with one occurrence
%     of each different element of A, in sorted order.
%
% Based on code by Asko Huuskonen, Finnish Meteorological Institute,
% Asko.Huuskonen@fmi.fi, posted to comp.soft-sys.matlab.
% This version by John Gilbert, Xerox PARC, 27 Sep 1993
% Copyright (c) 1990-1996 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.

if length(A) == 0
    B = []; 
    return;
end;

B = sort(A(:));
B(find(diff(B)==0)) = [];
B = B.';

function [Q,r] = reflector(c)
% REFLECTOR : Orthogonal transformation to put point on last axis.
%
% [Q,r] = reflector(c):
% Given a point c in d-space (as a row vector),
% return an orthogonal Q such that c*Q is on the last axis,
% and has magnitude r.  Q is a Householder reflection.
%
% See GEOPART.
%
% John Gilbert and Shanghua Teng, 1992-1993.
% Copyright (c) 1990-1996 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.

d = length(c);
p = d:-1:1;
[Q,r] = qr(c(p)');
Q = Q(p,p);
r = r(1);

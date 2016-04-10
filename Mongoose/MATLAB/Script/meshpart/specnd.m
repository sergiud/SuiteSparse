function p = specnd(S,minsep);
% SPECND : Spectral nested dissection ordering.
%
% p = specnd(S,minsep).  Nested dissection ordering of S.
% For a symmetric positive definite matrix S, this returns
% a nested dissection permutation p, so that S(p,p) tends to 
% have a sparser Cholesky factor than S.  
%
% minsep   (optional, default 3) is the smallest subgraph that will
%          be separated recursively.
%
% With no output argument, specnd reports statistics
% and draws a picture of the elimination tree.
%
% See also SPECPART, NDPERM, GSND, GEOND.
%
% John Gilbert, 1993.
% Copyright (c) 1990-1996 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.
%
% Modified by John Gilbert for Matlab 6, Feb 2002

if nargin < 2, minsep = 3; end  % This is the minsep default.
if minsep < 3, minsep = 3; end  % No point in separating less than 3 vertices.

p = ndperm('specpart',minsep,S);

if nargout == 0
    Sp = S(p,p);
    analyze(Sp);
    etreeplotg(Sp);
    title('Spectral Nested Dissection');
end;

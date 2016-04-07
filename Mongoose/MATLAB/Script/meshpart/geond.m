function p = geond(S,Sxy,minsep,ntries);
% GEOND : Geometric nested dissection ordering.
%
% p = geond(S,Sxy,minsep,ntries).  Nested dissection ordering of S.
% For a symmetric positive definite matrix S, this returns
% a nested dissection permutation p, so that S(p,p) tends to 
% have a sparser Cholesky factor than S.  
%
% Sxy      coordinates for the vertices of S.
% minsep   (optional, default 3) is the smallest subgraph that will
%          be separated recursively.
% ntries   (optional, default 20) is the number of random choices
%          to try for each separator.  
%
% With no output argument, geond reports statistics 
% and draws a picture of the elimination tree.
%
% See also GEOPART, NDPERM, SPECND, GSND.
%
% John Gilbert and Shanghua Teng, 1992-1993.
% Copyright (c) 1990-1996 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.
%
% Modified by John Gilbert for Matlab 6, Feb 2002

if nargin < 3, minsep = 3; end  % This is the minsep default.
if minsep < 3, minsep = 3; end  % No point in separating less than 3 vertices.

if nargin < 4
    ntries = 20;
end;

p = ndperm('geopart',minsep,S,Sxy,0,ntries);

if nargout == 0
    Sp = S(p,p);
    analyze(Sp);
    etreeplotg(Sp);
    title('Geometric Nested Dissection');
end;

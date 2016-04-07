function p = metisnd(S);
% METISND : Metis nested dissection ordering.
%
% p = metisnd(S):  Nested dissection ordering of S.
% For a symmetric positive definite matrix S, this returns
% a nested dissection permutation p, so that S(p,p) tends to 
% have a sparser Cholesky factor than S.  
%
% This just calls Metis's nested dissection routine and then
% postorders the elimination tree (same fill, nicer pictures).
% With no output argument, metisnd reports statistics 
% and draws a picture of the elimination tree.
%
% See also METISMEX (which accepts all the Metis options), 
% GEOND, GSND, SPECND, METISPART, METISDICE.
%
% John Gilbert  3 Jul 01
% Copyright (c) 1990-2001 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.
%
% Modified by John Gilbert for Matlab 6, Feb 2002

p = metismex('NodeND',S);
[t,q] = etree(S(p,p));
p = p(q);

if nargout == 0
    Sp = S(p,p);
    analyze(Sp);
    etreeplotg(Sp);
    title('Metis Nested Dissection');
end;

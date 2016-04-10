function etreeplotg(A,b,c,d)
% ETREEPLOTG : Plot the elimination tree.
%	etreeplotg(A):  Plot the elimination tree of A (or A+A', if asymmetric).
%	etreeplotg(A,1):  Also label the nodes with their numbers.
%	etreeplotg(A,c,d):  See treeplotg for optional parameters c and d.
%	etreeplotg(A,1,c,d):  See treeplotg for optional parameters c and d.
%
%	See also TREEPLOTG, ETREE, TREELAYOUTG.
% Copyright (c) 1990-1996 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.
%	Copyright (c) 1984-92 by The MathWorks, Inc.
%
% Modified by Tim Davis, for Matlab 5.1.  July 1998
% Modified by John Gilbert for Matlab 6 Feb 2002.

B = spones(A);
if nargin == 1,
    treeplotg(etree(B+B'));
elseif nargin == 2,
    treeplotg(etree(B+B'),b);
elseif nargin == 3,
    treeplotg(etree(B+B'),b,c);
else
    treeplotg(etree(B+B'),b,c,d);
end

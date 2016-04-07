function cutsize = sepquality(v,A,xyz)
% SEPQUALITY : Separator quality.
%
% cutsize = sepquality(v,A,xyz)
% Return the number of edges crossing a partition of the vertices of A,
% at positions xyz, by the plane described by v.
%
% See GEOPART.
%
% John Gilbert and Shanghua Teng, 1992-1993.
% Copyright (c) 1990-1996 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.

[a,b] = partition(xyz,v);
if min(length(a),length(b))
    cutsize = nnz(A(a,b));
else
    cutsize = 0;
end

function map = metisdice(A,a,b);
% METISDICE : Metis multiway partition.
%
% map = metisdice(A,nparts)
% A is the adjacency matrix of a graph.
% This uses Metis to divide A into nparts pieces of approximately
% equal size, with relatively small connections.
%
%       gsdice(A,nparts,xy) or
%       gsdice(A,xy,nparts):  Draw a picture of the result, as well.
%
% See also METISMEX (which accepts all the Metis options), 
% GEODICE, GSDICE, SPECDICE, METISPART, METISND.
%
% John Gilbert  3 Jul 01
% Copyright (c) 1990-2001 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.

% Sort out the input arguments.
picture = (nargin >= 3);
if picture
    if length(a) == 1
        nparts = a;
        xy = b;
    else
        nparts = b;
        xy = a;
    end;
else
    nparts = a;
end;
        
map = metismex('PartGraphKway',A,nparts);

if picture
    gplotmap(A,xy,map);
    title('Metis Multiway Partition');
end;

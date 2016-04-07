function  [part1,part2,chaco_time] = chaco(A,xy,method,nparts,goal)
% CHACO : Hendrickson/Leland's graph partitioner.
% 
%  [part1,part2,chaco_time] = chaco(A) returns a 50/50 vertex partition of the mesh
%  whose symmetric adjacency matrix is A.  
%
%  Optional arguments:
%  map = chaco(A, xy, method, nparts, goal);
%  map = chaco(A, xy, method, inmap, goal);
%
%  A:          Depending on "method", A may contain vertex and edge weights.
%
%  xy:         Each row of xy is the coordinates of a vertex.
%              If xy is non-null and there is no output, draw a picture.
%
%  method:     Scalar or vector describing the desired method.  
%              Default is multilevel Kernighan-Lin; other possibilities below.
%
%  nparts      Number of parts to divide into.  Default is 2.  If nparts is 
%    or        present, the output is a "map vector", see below.  (If method(5) 
%  inmap:      is specified, nparts is interpreted differently; see below.  In
%              any case, the default is to divide into two parts.)
%              If method(1) = 7 (see below), this argument is a map vector
%              specifying an initial 2-way partition, and Chaco refines it.
%
%  goal:       Optionally, a vector of desired sizes (or total vertex weights)
%              for each of the nparts parts.  Default is all sizes equal.
%
%  map:        If nparts and inmap are not present, the output is a vector of 
%              the n/2 vertex numbers in one part of the 2-way partition, for
%              compatibility with geopart and specpart.
%              If nparts or imap is present, the output is a vector of the
%              n part numbers, from 0 to nparts-1, assigned to the vertices.
%
% This is a Matlab interface to the graph partitioning software described
% in B. Hendrickson and R. Leland, "The Chaco User's Guide (Version 2.0)",
% Sandia National Laboratories report SAND94-2692, October 1994.
% This interface was written by John Gilbert, Xerox PARC, and is
% Copyright (c) 1994-1996 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.
%
% Modified by Tim Davis, for Matlab 5.1.  July 6, 1998.
%
% See also GEOPART, SPECPART.
%
% "method" is a vector of flags as follows.  Not all combinations are supported.
% See Section 6.10 of the Chaco manual for more details on all the arguments.
% If "method" is shorter than 10, we use the defaults for unspecified entries.
%
% method(1):  Global partitioning method  ("global_method" in the Chaco manual).
%             1 Multilevel Kernighan-Lin (default)
%             2 Spectral
%             3 Inertial
%             4 Linear
%             5 Random
%             6 Scattered
%             7 Use "inmap" as the global (2-way) partition
%
% method(2):  Local refinement method  ("local_method" in the Chaco manual).
%             1 Kernighan-Lin (default)
%             2 None
%
% method(3):  Vertex weighting.
%             0 No weights (default)
%             1 Use diag(A) as (positive integer) vertex weights
%
% method(4):  Edge weighting.
%             0 No weights (default)
%             1 Use off-diagonals of A as (positive integer) edge weights
%
% method(5):  Target architecture  ("architecture" in the Chaco manual).
%             If method(5) = 0, the target is a hypercube, "nparts" is the 
%             number of dimensions, and the partition is into 2^nparts parts.  
%             If method(5) = 1, 2, or 3, the target is a 1-, 2-, or 3-D grid,
%             "nparts" is a vector of the sizes of the grid in each dimension,
%             and the partition is into prod(nparts) parts.
%             Default is method(5) = 1, so nparts is the number of parts.
%
% method(6):  Partitioning dimension  ("ndims" in the Chaco manual).
%             1 Bisection (default)
%             2 Quadrisection
%             3 Octasection
%
% method(7):  Number of vertices to coarsen to  ("vmax" in the Chaco manual).
%             Default is 50.
%
% method(8):  Eigensolver  ("rqi_flag" in the Chaco manual).
%             0 RQI/Symmlq (default)
%             1 Lanczos 
%
% method(9):  Eigensolver convergence tolerance  ("eigtol" in the Chaco manual).
%             Default is .001
%
% method(10): Seed for random number generator  ("seed" in the Chaco manual).
%             Default is 7654321.
%
% Many esoteric details of Chaco's behavior can be changed by placing a file
% called "User_Params" in the same directory as the executable mlchaco.mex.
% As always, see the manual for details.

DefaultMethod = [1 1 0 0 1 1 50 0 .001 7654321];

% Fill in default arguments.
if nargin < 2, xy = []; end;
if nargin < 3, method = DefaultMethod; end;
if nargin < 4, nparts = []; end;
if nargin < 5, goal = []; end;
if length(method) < length(DefaultMethod)
    method = [method DefaultMethod(length(method)+1 : length(DefaultMethod))];
end;

% Decide on output and graphics.
if (isempty (nparts))
    mapvector = 0;
else
    mapvector = 1;
end;
picture = (nargout == 0) & (size(xy,2) >= 2);

% Chaco numbers vertices from 1 and the Matlab sparse data structure 
% numbers rows from 0, so we add an empty first row to make things line up.
% This code also makes sure the arg to Chaco will be sparse.
[n,n] = size(A);
Adiag = diag(diag(A));
Aout = [sparse(1,n) ; A-Adiag];

% Make sure all args except the adj matrix are full;
if issparse(xy)
    xy = full(xy);
end;
if issparse(method)
    method = full(method);
end;
if issparse(goal)
    goal = full(goal);
end;

% Decode "method" to get the actual args to Chaco.
% Note that "nparts" may correspond to any of several Chaco
% parameters, depending on the method.

global_method = method(1);
local_method = method(2);
if method(3)
    vwgts = full(Adiag);
    totalvwgt = sum(vwgts);
else
    vwgts = [];
    totalvwgt = size(A,2);
end;
ewgtsP = method(4);  % This is just true or false; the weights are in Aout.
architecture = method(5);
if global_method == 7
    % Refine an input partition: "nparts" is the input partition.
    % This seems to work only for hypercube architecture, 
    % so we force a 1-D hypercube with 2-way partitioning.
    architecture = 0;
    assignment = nparts;
    ndims_tot = 1;
    mesh_dims = [];
    ndims = 1;
    nsets = 2;
elseif architecture == 0
    % Partition for hypercube: "nparts" is # of dimensions (default 1).
    assignment = [];
    if nparts == []
        ndims_tot = 1;
    else
        ndims_tot = nparts;
    end;
    mesh_dims = [];
    ndims = method(6);
    nsets = 2^ndims_tot;
else
    % Partition for mesh: "nparts" is vector of mesh sizes in each
    % dimension, default [2 1 ... 1] with "architecture" dimensions.
    assignment = [];
    ndims_tot = [];
    if (isempty (nparts))
        mesh_dims = ones(1,architecture);
        mesh_dims(1) = 2;
    else
        mesh_dims = nparts;
    end;
    ndims = method(6);
    nsets = prod(mesh_dims);
end;
if length(goal) ~= nsets
    goal = totalvwgt/nsets * ones(1,nsets);
end;
vmax = method(7);
rqi_flag = method(8);
eigtol = method(9);
seed = method(10);

% The args to the mex-file interface to Chaco are almost the same as
% the args to the Chaco "interface" routine as described in the manual.

% For debugging, save the arguments:
%
%save mlchaco.mat Aout vwgts ewgtsP xy assignment architecture ...
%    ndims_tot mesh_dims goal global_method local_method rqi_flag ...
%    vmax ndims eigtol seed

[map,chaco_time]=mlchaco(Aout, vwgts, ewgtsP, xy, assignment, architecture, ...
    ndims_tot, mesh_dims, goal, global_method, local_method, rqi_flag, ...
    vmax, ndims, eigtol, seed);

% Draw the picture.
if picture
    if nsets == 2
        gplotpart(A,xy,find(map==0));
    else
        gplotmap(A,xy,map);
    end;
    if     method(1)==1, heading = 'Multilevel Kernighan-Lin';
    elseif method(1)==2, heading = 'Spectral';
    elseif method(1)==3, heading = 'Inertial';
    elseif method(1)==4, heading = 'Linear';
    elseif method(1)==5, heading = 'Random';
    elseif method(1)==6, heading = 'Scattered';
    elseif method(1)==7, heading = 'Input'; 
    end;
    heading = [heading ' Partition'];
    if method(2)==1 & method(1) ~= 1 
        heading =[heading ' Refined by KL'];
    end;
    title(heading);
end;

% Put output in the right form.
if mapvector
    part1 = map;
else
    part1 = find(map==0);
    part2 = find(map==1);
end;

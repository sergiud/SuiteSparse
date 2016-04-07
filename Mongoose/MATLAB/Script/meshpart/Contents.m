% Mesh Partitioning Toolbox.  Distribution version of 8 Feb 2002. 
% 
% Demonstrations.
%   Say "meshdemo" for some examples.
%   Also, many of the individual routines can draw pictures of what they do.
% 
% Partitioning methods.
%   geopart        - Geometric. 
%   specpart       - Spectral.
%   gspart         - Geometric spectral.
%   coordpart      - Coordinate bisection.
%   inertpart      - Inertial bisection.
%   chaco          - Multilevel Kernighan-Lin, and other options.
%   metispart      - Multilevel method from Metis.
%   metismex       - Interface to more options of Metis.
% 
% Multiway partitions.
%   dice           - Use any 2-way partitioner to get a multiway partition.
%   geodice        - Recursive geometric partitioning.
%   specdice       - Recursive spectral partitioning.
%   gsdice         - Recursive geometric spectral partitioning.
%   chaco          - Can also produce multiway partitions directly.
%   metisdice      - Multiway partitioning from Metis.
% 
% Vertex separators.
%   vtxsep         - Convert a 2-way partition to a vertex separator.
%   geosep         - Vertex separator from geometric partitioning.
%   specsep        - Vertex separator from spectral partitioning.
% 
% Nested dissection.
%   ndperm         - Use any 2-way partitioner for nested dissection.
%   geond          - Geometric nested dissection ordering.
%   specnd         - Spectral nested dissection ordering.
%   gsnd           - Geometric spectral nested dissection ordering.
%   metisnd        - Nested dissection ordering from Metis.
%   analyze        - Predict fill, opcount, etc. for an elimination ordering.
% 
% Meshes and graph generators.
%   meshes.mat     - Three sample meshes with coordinates:
%     Eppstein     - A 2D finite-element mesh with 547 nodes.
%     Smallmesh    - A 2D finite-element mesh with 136 nodes.
%     Tapir        - A 2D finite-element mesh with 1024 nodes.
%   grid5          - 2D square 5-point mesh.
%   grid7          - 2D square 7-point mesh.
%   grid9          - 2D square 9-point mesh.
%   gridt          - 2D triangular mesh.
%   grid3d         - 3D cubical mesh.
%   grid3dt        - 3D cubical simplicial mesh.
%   badmesh        - A mesh that has no good straight-line cut.
%   cockroach      - A mesh for which spectral bisection does poorly.
%   treexpath      - A mesh for which spectral bisection does poorly.
% 
% Visualization and graphics.
%   (All the partitioners can also draw pictures of what they do.)
%   gplotpart      - Draw a 2-way partition.
%   gplotmap       - Draw a multiway partition.
%   highlight      - Draw a mesh with some vertices highlighted.
%   gplotg         - Draw a 2D or 3D mesh (replaces Matlab's gplot).
%   etreeplotg     - Draw an elimination tree (replaces Matlab's etreeplot).
%   spypart        - Matrix spy plot with partition boundaries.
%   dmspy          - Spy plot of matrix in block triangular form.
%  
% Utilities.
%   cutsize        - Find or count edges cut by a partition.
%   other          - Other side of a partition, or change representations.
%   intersection   - Intersection of two sets.
%   union          - Union of two sets.
%   fiedler        - Fiedler vector of a graph.
%   laplacian      - Laplacian matrix of a graph.
%   components     - Connected components of a graph.
%   contract       - Condense a graph according to a given block structure.
%   distances      - Distances between adjacent mesh points.
%   blockdiags     - Create matrix with specified block diagonals.
%   resetrandoms   - Reset random number generators to startup values.
%
% See meshpart.html for references.

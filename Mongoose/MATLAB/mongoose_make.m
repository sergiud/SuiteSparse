function mongoose_make ( )
%MONGOOSE_MAKE compiles the GP mexFunctions

details = 0 ;	    % 1 if details of each command are to be printed

v = getversion ;
try
    % ispc does not appear in MATLAB 5.3
    pc = ispc ;
catch                                                                       %#ok
    % if ispc fails, assume we are on a Windows PC if it's not unix
    pc = ~isunix ;
end

d = '' ;
is64 = (~isempty (strfind (computer, '64'))) ;
if (is64)
    % 64-bit MATLAB
    d = '-largeArrayDims' ;
end

include = '-I. -I../Include -I../../SuiteSparse_config' ;

% Linux/Unix require these flags for large file support
include = [include ' -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE'] ;

% We're compiling this from within a mex function.
include = [include ' -DGP_MEX_FUNCTION'] ;

lib = '-L../Lib -L../../SuiteSparse_config -lsuitesparseconfig -lrt' ;

% Fix the include & library path.
include = strrep (include, '/', filesep) ;
lib = strrep (lib, '/', filesep) ;

%-------------------------------------------------------------------------------

gp_src = {
    '../Source/cs', ...
    '../Source/Coarsening', ...
    '../Source/Conditioning', ...
    '../Source/EdgeSeparator', ...
    '../Source/GuessCut', ...
    '../Source/Matching', ...
    '../Source/Refinement', ...
    '../Source/Waterdance', ...
    '../Source/BoundaryHeap', ...
    '../Source/FidduciaMattheyes', ...
    '../Source/QPboundary', ...
    '../Source/QPgradproj', ...
    '../Source/QPlinks', ...
    '../Source/QPminheap', ...
    '../Source/QPmaxheap', ...
    '../Source/QPnapdown', ...
    '../Source/QPnapup', ...
    '../Source/QPnapsack', ...
    '../Source/Nuri_qpGradProj', ...
    '../Source/Interop', ...
    './mex_get_graph', ...
    './mex_get_options', ...
    './mex_getput_vector', ...
    './mex_put_options', ...
    './mex_struct_util' } ;

gp_mex_src = { 'gp_getDefaultOptions', 'gp_computeEdgeSeparator', 'gp_exportGraph' } ;

obj_extension = '.o' ;

% compile each library source file
obj = '' ;

kk = 0 ;

fprintf('\n\nBuilding Mongoose\n') ;
for f = gp_src
    ff = strrep (f{1}, '/', filesep) ;
    slash = strfind (ff, filesep) ;
    if (isempty (slash))
        slash = 1 ;
    else
        slash = slash (end) + 1 ;
    end
    o = ff (slash:end) ;
    obj = [obj ' ' o obj_extension] ;                      %#ok
    s = sprintf ('mex %s -DDLONG -O %s -c %s.cpp', d, include, ff) ;
    kk = do_cmd (s, kk, details) ;
end

fprintf('\n\nBuilding MEX\n');
for f = gp_mex_src
    s = sprintf ('mex %s -DDLONG -O %s %s.cpp', d, include, f{1}) ;
    s = [s obj ' ' lib] ;                                  %#ok
    kk = do_cmd (s, kk, details) ;
end

% clean up
s = ['delete ' obj] ;
do_cmd (s, kk, details) ;
fprintf ('\nMongoose successfully compiled\n') ;

%-------------------------------------------------------------------------------
function kk = do_cmd (s, kk, details)
%DO_CMD: evaluate a command, and either print it or print a "."

if (details)
    fprintf ('%s\n', s) ;
else
    if (mod (kk, 60) == 0)
	fprintf ('\n') ;
    end
    kk = kk + 1 ;
    fprintf ('.') ;
end
eval (s) ;

%-------------------------------------------------------------------------------
function v = getversion
% determine the MATLAB version, and return it as a double.
v = sscanf (version, '%d.%d.%d') ;
v = 10.^(0:-1:-(length(v)-1)) * v ;

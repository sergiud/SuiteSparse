%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function filter = gp_getFilter(inclusion, exclusion, forceRebuild)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INCLUSION
% 1) Everything
% 2) Small Problems
% 3) Medium Problems
% 4) Large Problems
% 5) Problems with max degree 0.1*n
% 6) Problems with max degree 0.01*n
% 7) Problems with max degree 0.001*n
%
% EXCLUSION
% 1) Not 2D/3D Mesh Problems      (set if you want ONLY 2D/3D problems)
% 2) 2D/3D Mesh Problems          (set if you want NO 2D/3D problems)
% 3) Not "graph" Problems         (set if you want ONLY graph problems)
% 4) "graph" Problems             (set if you want NO graph problems)
% 5) Not "optimization problem"s  (set if you want ONLY these)
% 6) "optimization problem"s      (set if you want NO these)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Configure and handle variable inputs.
if nargin < 1, return; end
if nargin < 2, exclusion = 0; end
if nargin < 3, forceRebuild = false; end

% Build file IO
outputDir = '~/workspace/GraphPartitioning/Script/filter' ;
filename = sprintf('%s/filter-i%d-e%d', outputDir, inclusion, exclusion) ;

% If we're not forcing a rebuild then try to get filter from cache.
if ~forceRebuild
    try
        load(filename);
        return;
    catch
    end
end

% Build a new one from scratch
index = UFget() ;
filter = find(index.nnz)' ;
filter(:,2) = 0 ;

% Depending on what we want to include:
switch inclusion

    % Include everything
    case 1
        filter(:,2) = 1 ;

    % Include Small Problems (nz < 10000)
    case 2
        filter(find(index.nnz < 1e4), 2) = 1 ;

    % Include Medium Problems (nz >= 10000 && nz < 1000000)
    case 3
        filter(find(1e4 <= index.nnz < 1e6), 2) = 1 ;

    % Include Large Problems (nz >= 1000000)
    case 4
        filter(find(index.nnz >= 1e6), 2) = 1 ;

    % Include problems with dense nodes.
    case {5, 6, 7}
        nMin = 10000;
        nMax = 1000000;
        switch(inclusion)
            case 5, density = 0.100;
            case 6, density = 0.010;
            case 7, density = 0.001;
        end

        for i = 1:length(index.nnz)
            if(index.ncols(i) < nMin || index.ncols(i) > nMax), continue; end;
            if index.nrows ~= index.ncols, continue; end;

            problem = UFget(i);
            if ~isreal(problem.A), continue; end;

            degrees = sum(logical(problem.A));
            maxDegree = max(degrees);

            if maxDegree > density*index.ncols(i), filter(i,2) = 1; end
        end
end

% Depending on the exclusion parameter, exclude things from the return set.
switch exclusion

    % Exclude non-2D/3D mesh problems.
    case 1
        filter(find(~index.isND), 2) = 0 ;

    % Exclude 2D/3D mesh problems.
    case 2
        filter(find(index.isND), 2) = 0 ;

    % Exclude non-graph problems.
    case 3
        kinds = UFkinds;
        for n = 1:length(kinds)
            if isempty(strfind(char(kinds(n)), 'graph')), filter(n, 2) = 0; end
        end

    % Exclude graph problems.
    case 4
        kinds = UFkinds;
        for n = 1:length(kinds)
            if ~isempty(strfind(char(kinds(n)), 'graph')) filter(n, 2) = 0; end
        end

    % Exclude non-"optimization problem"s.
    case 5
        kinds = UFkinds;
        for n = 1:length(kinds)
            if isempty(strfind(char(kinds(n)), 'optimization problem')), filter(n, 2) = 0; end
        end

    % Exclude "optimization problem"s.
    case 6
        kinds = UFkinds;
        for n = 1:length(kinds)
            if ~isempty(strfind(char(kinds(n)), 'optimization problem')) filter(n, 2) = 0; end
        end

end

save(filename, 'filter') ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [config] = gp_getConfig(inConfig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT VALUES
% config.kstart = 1 ;
% config.kend = 2530 ;                  (k=2537 id=2455 too big for mp machine)
% config.experimentList = [1:3] ;
% config.ewList = [0:1] ;
% config.methodList = [] ;
% config.dmperm = true ;
% config.recompile = true ;
% config.dofast = false ;
% config.maxtime = 1 ;
% config.outputDir = '~/Desktop/results' ;
% config.exportGraph = false ;
% config.exportPath = '~/workspace/GraphPartitioning/graphs/troll.mtx' ;
% config.wSeed1 = 15761 ;
% config.wSeed2 = 79059 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% experimentList - contains the id of the problem type we want to run:
%   1) Square Symmetric (only applicable if the matrix is square)
%   2) Bipartite Graph
%   3) AA' or A'A (depending on the larger of the dimensions)
%
% methodList - contains the various options and methods to test:
%   1) fprintf('METIS - Default');
%   2) fprintf('GP - Default');
%   3) fprintf('GP - Consider Quads');
%   4) fprintf('GP - Ignore Trinary & Quads');
%   5) fprintf('GP - Ignore Gemini & Trinary & Quads');
%   6) fprintf('GP - 2 scale');
%   7) fprintf('GP - 3 scale');
%   8) fprintf('GP - 4 scale');
%   9) fprintf('GP - 0.010 tolerance');
%  10) fprintf('GP - 0.015 tolerance');
%  11) fprintf('GP - 900 coarsen limit');
%  12) fprintf('GP - 800 coarsen limit');
%  13) fprintf('GP - 700 coarsen limit');
%  14) fprintf('GP - 600 coarsen limit');
%  15) fprintf('GP - 500 coarsen limit');
%  16) fprintf('GP - 400 coarsen limit');
%  17) fprintf('GP - 300 coarsen limit');
%  18) fprintf('GP - 200 coarsen limit');
%  19) fprintf('GP - 100 coarsen limit');
%  20) fprintf('GP - 50 coarsen limit');
%  21) fprintf('GP - 0.95 coreThreshold');
%  22) fprintf('GP - 0.90 coreThreshold');
%  23) fprintf('GP - 0.85 coreThreshold');
%  24) fprintf('GP - 0.80 coreThreshold');
%
% ewList - whether we want to consider edge weights in the problem:
%   0) Ignore Edge Weights
%   1) Consider Edge Weights
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Allow the user to pass a partially constructed configuration structure.
config = [] ;

if nargin == 1
    % If the input is not a structure then return the default config.
    if ~isstruct(inConfig)
        config = gp_getConfig() ;
        return ;
    % Else sanitize the input.
    else
        config = inConfig ;
    end
end

% Fill any unspecified fields with these defaults.
if ~isfield(config, 'kstart'), config.kstart = 1; end
if ~isfield(config, 'kend'), config.kend = 2200; end
if ~isfield(config, 'experimentList'), config.experimentList = [1]; end
if ~isfield(config, 'ewList'), config.ewList = [0:1]; end
if ~isfield(config, 'methodList'), config.methodList = []; end
if ~isfield(config, 'dmperm'), config.dmperm = true; end
if ~isfield(config, 'recompile'), config.recompile = true; end
if ~isfield(config, 'dofast'), config.dofast = true; end
if ~isfield(config, 'maxtime'), config.maxtime = 1; end
if ~isfield(config, 'outputDir'), config.outputDir = '~/Desktop/results'; end
if ~isfield(config, 'exportGraph'), config.exportGraph = false; end
if ~isfield(config, 'exportPath'), config.exportPath = '~/workspace/Mongoose/troll.mtx'; end
if ~isfield(config, 'wSeed1'), config.wSeed1 = 15761; end
if ~isfield(config, 'wSeed2'), config.wSeed2 = 97059; end

% If not provided, configure problem ids sorted by number of edges ascending.
if ~isfield(config, 'ids')
    index = UFget ;
    ids = find(index.nnz) ;
    [ignore, i] = sort(index.nnz(ids)) ;
    config.ids = ids(i) ;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gp_runExperiment()
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

config.recompile = true ;
config.dmperm = true ;
config.dofast = true ;

config.experimentList = [1] ;
config.ewList = [0] ;

config.kstart = 1 ;
config.kend   = 2000 ;

% Configure METIS experiment list.
%config.dmperm = true ;
%config.methodList = [1] ;
%config.outputDir = '~/Desktop/results/FRODO/dmperm/METIS' ;
%gp_execfight(config) ;
%config.dmperm = false ;
%config.outputDir = '~/Desktop/results/FRODO/nodmperm/METIS' ;
%gp_execfight(config) ;

% Configure primary experiment list.
config.dmperm = true ;
config.methodList = [70,72] ;
config.outputDir = '~/Desktop/results/FRODO/dmperm/GP' ;
gp_execfight(config) ;
%config.dmperm = false ;
%config.outputDir = '~/Desktop/results/FRODO/nodmperm/GP' ;
%gp_execfight(config) ;


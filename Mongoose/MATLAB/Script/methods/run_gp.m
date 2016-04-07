function [time,cost,bal,P] = run_gp(A, ew, O, config)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A          - square symmetric matrix
% ew         - 1 or 0 flag indicating whether edge weights should be considered
% O          - the matlab struct that represents user options
% config     - the input configuration structure for experiments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find the GP edge cut.
runs = 0 ; tElapsed = 0 ; failure = 0 ; maxtime = config.maxtime ;
partition = [];

tStart = tic ;
while tElapsed < maxtime && failure == 0
    try
        partition = gp_computeEdgeSeparator(A,O);
        runs = runs + 1;
        tElapsed = toc(tStart);
    catch exception
        fprintf('GP threw an exception\n') ;
        failure = 1 ;
    end
    if config.dofast == true, break; end
end

% If we didn't fail then save the details to the record
if failure == 0
    time = tElapsed / runs ;

    if nargout >= 4
        [cost,bal,P] = gp_interpretPart(A, partition) ;
    else
        [cost,bal] = gp_interpretPart(A, partition) ;
    end
else
    time = Inf ;
    cost = Inf ;
    bal = Inf ;
    if nargout >= 3, P = 0; end
end


function [time,cost,bal,P] = run_metis(A, ew, config)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A  - square symmetric matrix
% ew - 1 or 0 flag indicating whether edge weights should be considered
% config - the input configuration structure for experiments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find the METIS edge cut.
runs = 0 ; tElapsed = 0 ; failure = 0 ; maxtime = config.maxtime ;
tStart = tic;
while tElapsed < maxtime && failure == 0
%    try
        map = metismex4('PartGraphRecursive', A, 2, 1);
        [pL,pR] = other(map);

        runs = runs + 1;
        tElapsed = toc(tStart);

%    catch exception
%        fprintf('METIS threw an exception\n') ;
%        failure = 1 ;
%    end
    if config.dofast, break; end
end

% If we didn't fail then save the details to the record
if failure == 0
    % Configure the METIS partition choices.
    partition(1,pR) = 0 ;
    partition(1,pL) = 1 ;
    partition = partition' ;

    % Save data to the record
    time = tElapsed / runs ;

    if nargout >= 3
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


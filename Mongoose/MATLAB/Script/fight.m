%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T,C,B] = fight(config, method, experiment, ew, fail)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% config must be a matlab structure at this point. Terminate if it isn't.
if ~isstruct(config), return; end

% Read values from the configuration structure.
kstart = config.kstart ;
kend = config.kend ;
do_dmperm = config.dmperm ;
ids = config.ids ;
exportGraph = config.exportGraph ;
exportPath = config.exportPath ;

SQSYM  = 1 ;
BIPART = 2 ;
AAT    = 3 ;

% Initialize the outputs as all empty.
T = zeros(size(ids)) ;
C = zeros(size(ids)) ;
B = zeros(size(ids)) ;

k = 1 ;
for id = ids

    % Find out if the k is a failure case.
    isFailure = size(find(fail == id), 2) ;

    % If we're within our run bounds and we aren't a failure case then:
    if k >= kstart && isFailure == 0

        % Log the k and id as well as a timestamp.
        gp_printDetails(k, id) ;

        % Get the conditioned problem and skip it if we don't have enough memory
        [A,ex] = gp_getConditionedProblem(id,ew) ;
        if ex == 1
            k = k + 1 ;
            if k > kend, break; end
            continue;
        end

        % If we're square, form the symmetric graph.
        dosq = 0 ;
        if experiment == SQSYM
            try
                if size(A,1) == size(A,2)
                    sqA = 0.5 * (A + A') ;
                    A = sqA ;
                    dosq = 1 ;
                end
            catch exception
                gp_log('Error: Out of Memory when forming square symmetric graph.\n') ;
            end
        end

        % Form the bipartite graph.
        dobi = 0 ;
        if experiment == BIPART
            try
                biA = spaugment(A, 0) ;
                A = biA ;
                dobi = 1 ;
            catch
                gp_log('Error: Out of Memory when forming bipartite graph.\n') ;
            end
        end

        % Form AA' or A'A
        doaat = 0 ;
        if experiment == AAT
            try
                if size(A,1) >= size(A,2)
                    aatA = A'*A ;
                else
                    aatA = A*A' ;
                end
                A = aatA ;
                doaat = 1 ;
            catch
                gp_log('Error: Out of Memory when forming AAt (or AtA).\n') ;
            end
        end

        % Select the largest component to cut.
        if dosq || dobi || doaat
            [A,ex] = gp_getLargestComponent(A);
            if ex == 1
                k = k + 1 ;
                if k > kend, break; end
                continue;
            end
        end

        % If we want to export the graph then export it to the troll graph.
        if exportGraph, gp_export(A) ; end

        % Default the time and cost for both experiments to Inf
        time = Inf; cost = Inf; bal = Inf;

        % Partition using METIS
        switch method
            case 1
                if dosq,  [time,cost,bal] = run_metis4(A, ew, config); end
                if dobi,  [time,cost,bal] = run_metis4(A, ew, config); end
                if doaat, [time,cost,bal] = run_metis4(A, ew, config); end

            case 2
                if dosq,  [time,cost,bal] = run_metis5(A, ew, config); end
                if dobi,  [time,cost,bal] = run_metis5(A, ew, config); end
                if doaat, [time,cost,bal] = run_metis5(A, ew, config); end

            % Partition using GP
            otherwise

                % Generate the UserParams based on the method.
                O = gp_generateOptions(method) ;

                % Run the appropriate experiment.
                if dosq,  [time,cost,bal] = run_gp(A, ew, O, config); end
                if dobi,  [time,cost,bal] = run_gp(A, ew, O, config); end
                if doaat, [time,cost,bal] = run_gp(A, ew, O, config); end
        end

        % Save results
        T(k) = time ;
        C(k) = cost ;
        B(k) = bal ;

        % Log results
        if dosq || dobi || doaat
            % Extra logging for debugging purposes.
            gp_log(sprintf('   Time = %f\n', time), false) ;
            gp_log(sprintf('   Cost = %f\n', cost), false) ;
            gp_log(sprintf('   Imba = %f\n', bal), false) ;
        else
            gp_log(sprintf('   Not suitable for experiment.\n'), false) ;
        end
    else
        if isFailure == 1
            T(k) = Inf ;
            C(k) = Inf ;
            B(k) = Inf ;
        end
    end

    k = k + 1 ;
    if k > kend, break; end
end

% Prepend the list of ids and transpose the results so that they're by column.
T = vertcat(ids, T)';
C = vertcat(ids, C)';
B = vertcat(ids, B)';


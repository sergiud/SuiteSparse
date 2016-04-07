%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,ex] = gp_getConditionedProblem(probID, ew)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 1
    gp_log('Warning: gp_getConditionedProblem missing param ew. Using 0 as default.\n') ;
    ew = 0;
end

A = [] ;
ex = 0 ;

try
    % Get the problem by ID.
    Prob = UFget(probID) ;
    A = Prob.A ;

    % Condition the edge weights.
    if ew == 0
        A = spones(A) ;
    else
        A = abs(A) ;
    end

    % Remove the diagonal.
    A = tril(A,-1) + triu(A,1) ;

catch exception
    ex = 1 ;
    gp_log('Error: Out of Memory when preconditioning graph for experiment.\n') ;
end

if nnz(A) == 0
    ex = 1;
    gp_log('Error: nnz(A) == 0.\n');
    return ;
end


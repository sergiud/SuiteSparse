%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gp_runSingle(probID, ew, useDMPerm, methodID, recompile)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make sure we're running the latest GP.
if nargin == 5
    if recompile, make_mongoose; end;
end

% Get the conditioned problem and skip it if we don't have enough memory
[A,ex] = gp_getConditionedProblem(probID,ew) ;
if ex == 1, return; end

dosq = 0 ;
try
    if size(A,1) == size(A,2)
        sqA = 0.5 * (A + A') ;
        A = sqA ;
        dosq = 1 ;
    end
catch exception
    gp_log('Error: Out of Memory when forming square symmetric graph.\n') ;
end

if dosq
    [A,ex] = gp_getLargestComponent(A);
    if ex == 1, return; end

    gp_export(A);

    for mid=methodID

%        try
            % Do the cut.    
            if(mid ~= 1)
                gp_showCut(A, ew, gp_generateOptions(mid)) ;
            else
                metis_showCut(A, ew) ;
            end
%        catch
%            gp_log('Error: Unhandled exception encountered.\n') ;
%        end

        fprintf('Press any key to continue...\n');
        pause ;
    end
end


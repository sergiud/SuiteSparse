%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gp_execfight(config)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sanitize the incoming config structure.
config = gp_getConfig(config) ;

% Read parameters from the configuration structure.
outputDir = config.outputDir ;
recompile = config.recompile ;

% If we're supposed to recompile then do so now.
if recompile, make_mongoose; end

% This is the bulk of this script. We run through every combination of tests,
% and we save them to files so that they can be accessed later.

for method = config.methodList

    for experiment = config.experimentList

        for ew = config.ewList

            % Log the executive details.
            gp_printExecDetails(method, experiment, ew) ;

            % Get the fail vector, which depends on method, experiment, and ew.
            fail = gp_getFailVector(method, experiment, ew) ;

            [T,C,B] = fight(config, method, experiment, ew, fail) ;

            if config.dofast
                filename = sprintf('%s/fast-E%d-EW%d-M%d.mat', outputDir, experiment, ew, method) ;
                save(filename, 'C', 'B') ;
            else
                filename = sprintf('%s/E%d-EW%d-M%d.mat', outputDir, experiment, ew, method) ;
                save(filename, 'T', 'C', 'B') ;
            end
        end
    end
end


function gp_showCut(A, inEW)
format compact;

% Get the default config
config = gp_getConfig();
config.dofast = true;

% Configure whether we should consider edge weights.
ew = 0;
if nargin >= 2, ew = inEW; end

% Do the cut and save the time, cost, balance, and Partition.
[T,C,B,P] = run_metis(A, ew, config);
fprintf('  Exec Time: %f\n', T);
fprintf('   Cut Cost: %f\n', C);
fprintf('  Imbalance: %f\n', B);

% Visualize the cut.
cspy(P);


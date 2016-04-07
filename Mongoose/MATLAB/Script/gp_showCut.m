function gp_showCut(A, inEW, O)
format compact;

% Get default options if none were provided on input
if nargin < 3, O = gp_getDefaultOptions(); end

% Get the default config
config = gp_getConfig();
config.dofast = true;

% Configure whether we should consider edge weights.
ew = 0;
if nargin >= 2, ew = inEW; end

% Do the cut and save the time, cost, balance, and Partition.
[T,C,B,P] = run_gp(A, ew, O, config);
fprintf('  Exec Time: %f\n', T);
fprintf('   Cut Cost: %f\n', C);
fprintf('  Imbalance: %f\n', B);

% Visualize the cut.
cspy(P);


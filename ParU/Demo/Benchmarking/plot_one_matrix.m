function ttmax = plot_one_matrix (fignum, name, what_result, package_names, package_timings, threads, base_package)
%
% inputs:
% fignum: figure number to use
% name: name of the matrix, as a string
% what_result: what the times are (factor, solve, sym, all, etc)
% package_names: cell array of strings, name of each package
%   (umf+amd, paru+met, etc)
% package_timings: cell array of run times, one per package
% threads: a list of the thread counts used for each run time
% base_package: an integer in range 1 to # of packages

figure (fignum)
clf (fignum)
hold on
name (find (name == '_')) = '-' ;
title (sprintf ('matrix: %s, what: %s', name, what_result)) ;

npackages = length (package_names) ;
ntimes = length (threads) ;
assert (isequal (npackages, length (package_timings))) ;

base_times = package_timings {base_package} ;

s = { 'ko-', 'bo-', 'go-', 'ro-', 'co-', 'mo-', 'yo-' } ;

for k = 1:ntimes
    thread_counts {ntimes-k+1} = sprintf ('%d', threads (k)) ;
end

ttmax = 1 ;

for k = 1:npackages

    times = package_timings {k} ;
    assert (isequal (size (times), [1 ntimes])) ;

    t = base_times ./ times ;
    t = t (end:-1:1) ;
    plot (1:ntimes, t, s {k}, 'LineWidth', 2) ;
    xticks (1:ntimes) ;
    xticklabels (thread_counts) ;

    tmax = max (t) ;
    ttmax = max (ttmax, tmax) ;

end

legend (package_names) ;

fprintf ('max speedup vs base package %s: %g\n', ...
    package_names {base_package}, ttmax) ;

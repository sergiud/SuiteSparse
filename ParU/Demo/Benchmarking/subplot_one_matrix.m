function ttmax = subplot_one_matrix (name, what_result, package_names, package_timings, threads, base_package, s)
%
% inputs:
% name: name of the matrix, as a string
% what_result: what the times are (factor, solve, sym, all, etc)
% package_names: cell array of strings, name of each package
%   (umf+amd, paru+met, etc)
% package_timings: cell array of run times, one per package
% threads: a list of the thread counts used for each run time
% base_package: an integer in range 1 to # of packages

name (find (name == '_')) = '-' ;

if (isempty (what_result))
    title (sprintf ('%s', name), 'FontSize', 14) ;
else
    title (sprintf ('matrix: %s, what: %s', name, what_result)) ;
end

npackages = length (package_names) ;
ntimes = length (threads) ;
assert (isequal (npackages, length (package_timings))) ;

% use the 1-threaded base method (typically UMFPACK, with AMD/COLAMD), as the base time
base_times = package_timings {base_package} ;
base_time = base_times (end) ;

if (nargin < 7)
    s = { 'ko-', 'bo-', 'go-', 'ro-', 'co-', 'mo-', 'yo-' } ;
end

for k = 1:ntimes
    thread_counts {ntimes-k+1} = sprintf ('%d', threads (k)) ;
end

ttmax = 1 ;

for k = 1:npackages

    times = package_timings {k} ;
    assert (isequal (size (times), [1 ntimes])) ;

%   t = base_times ./ times ;
    t = base_time ./ times ;

    t = t (end:-1:1) ;
    p = plot (1:ntimes, t, s {k}, 'LineWidth', 2) ;
    xlim ([1 ntimes]) ;
    xticks (1:ntimes) ;
    xticklabels (thread_counts) ;
    ax = gca ;
    ax.XAxis.FontSize = 14 ;
    ax.YAxis.FontSize = 14 ;

    tmax = max (t) ;
    ttmax = max (ttmax, tmax) ;

end

fprintf ('max speedup vs base package %s: %g\n', ...
    package_names {base_package}, ttmax) ;

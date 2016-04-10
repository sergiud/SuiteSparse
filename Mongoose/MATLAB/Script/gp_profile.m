function gp_profile(R, figNum, ids)

% Read the optional filter parameter.
if nargin < 3 || size(ids,1) == 0
    hasFilter = false;
    filter = [];    
else
    hasFilter = true;
    filter = ids;
end

% Read the number of entries and number of methods.
numEntries = size(R,1) ;
numMethods = size(R,2)-1 ;                 % Column 1 is a list of problem IDs.

% Build the Rf (R filtered).
Rf = [] ;
if hasFilter
    goodIDs = filter(find(filter(:,2) == 1), 1) ;
    first = true ;    
    for n = 1:length(goodIDs)
        row = R(find(R(:,1) == n), 2:(numMethods+1)) ;
        if first
            Rf = row ;
            first = false;
        else
            Rf = vertcat(Rf, row) ; 
        end
    end
else
    Rf = R(:,2:(numMethods+1)) ;
end

% Do the trick with the min
Rmin = min (Rf')' ;
Rel = Rf ./ repmat (Rmin, 1, numMethods) ;
Rel = sort (Rel) ;

% Configure colors
c = {'b', 'g', 'r', 'c', 'k', 'bo', 'go', 'ro', 'co', 'ko' } ;

% Do the plot
numFilteredEntries = size(Rel,1);
subplot(1, 3, figNum) ;
plot (Rel, 1:numFilteredEntries) ;
axis ([0.8 2.5 0 numFilteredEntries]) ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T,C,B] = buildResults2(experiment, ew, methodList, fast, ids)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T = [] ; C = [] ; B = [] ;
first = 1 ;
drawT = [] ; drawC = [] ; drawB = [] ;

if nargin < 4, fast = true; end
if nargin < 5, ids = []; end

for method = methodList
    if fast
        filename = sprintf('~/Desktop/results/fast-E%d-EW%d-M%d.mat', experiment, ew, method) ;
    else
        filename = sprintf('~/Desktop/results/E%d-EW%d-M%d.mat', experiment, ew, method) ;
    end
    fprintf('Comparing %s\n', filename) ;

    load(filename) ;
    filename = '';

    badBalanceCount = 0 ;
    for i=1:size(B,1)
        if B(i,2) >= 0.40
            badBalanceCount = badBalanceCount + 1 ;            
%            fprintf('Problem %d is out of balance.\n', B(i,1)) ;
            C(i,2) = inf ;
        end
    end
    fprintf('Bad Balance Count = %d\n', badBalanceCount);

    if first == 1
        if ~fast, drawT = T; end
        drawC = C;
        drawB = B;
        first = 0 ;
    else
        if ~fast, drawT = horzcat(drawT, T(:,2)); end
        drawC = horzcat(drawC, C(:,2));
        drawB = horzcat(drawB, B(:,2));
    end
end

% Draw the result plots.
if ~fast, gp_profile(drawT, 1, ids) ; end
gp_profile(drawC, 2, ids) ;
gp_profile(drawB, 3, ids) ;

T = drawT ;
C = drawC ;
B = drawB ;

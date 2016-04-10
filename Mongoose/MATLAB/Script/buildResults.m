%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T,C,B] = buildResults(experiment, ew, methodList, fast, ids)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T = [] ; C = [] ; B = [] ;
first = 1 ;
drawT = [] ; drawC = [] ; drawB = [] ;

doFast = false ;
if nargin >= 4, doFast = fast; end

useIDS = false;
if nargin >= 5, useIDS = true; end

for method = methodList
    if doFast
        filename = sprintf('~/Desktop/results/fast-E%d-EW%d-M%d.mat', experiment, ew, method) ;
    else
        filename = sprintf('~/Desktop/results/E%d-EW%d-M%d.mat', experiment, ew, method) ;
    end
    fprintf('Comparing %s\n', filename) ;

    load(filename) ;
    filename = '';

    if first == 1
        if ~doFast, drawT = T; end
        drawC = C;
        drawB = B;
        first = 0 ;
    else
        if ~doFast, drawT = horzcat(drawT, T(:,2)); end
        drawC = horzcat(drawC, C(:,2));
        drawB = horzcat(drawB, B(:,2));
    end
end

% Draw the result plots.
if useIDS
    if ~doFast, gp_profile(drawT, 1, ids) ; end
    gp_profile(drawC, 2, ids) ;
    gp_profile(drawB, 3, ids) ;
else
    if ~doFast, gp_profile(drawT, 1) ; end
    gp_profile(drawC, 2) ;
    gp_profile(drawB, 3) ;
end

T = drawT ;
C = drawC ;
B = drawB ;

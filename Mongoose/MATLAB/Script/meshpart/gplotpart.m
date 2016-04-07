function handle = gplotpart(A,xy,part1,color1,color2,color3)
% GPLOTPART : Plot a partitioned graph in 2 or 3 dimensions.
%
%	gplotpart(A,xy,part1) plots the n-vertex graph specified 
%       by the n by n adjacency (or Laplacian) matrix A 
%       and the n by 2 or 3 matrix of coordinates xy.
%       Argument part1 is a vector of vertex names (integers 1:n);
%       the subgraphs induced by part1 and its complement are plotted
%       in different colors, with the edges joining them in a third color.
%       Three more optional arguments give the three colors.
%
% See also GPLOTMAP.
%
% John Gilbert  9 March 1993
% Copyright (c) 1990-1996 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.

if nargin < 3, 
    part1 = [];
end;
if nargin < 4,
    color1 = 'yellow';
end;
if nargin < 5,
    color2 = 'cyan';
end;
if nargin < 6,
    color3 = 'red';
end;

[n,n] = size(A);
part1 = part1(:);
part2 = 1:n;
part2(part1)=zeros(size(part1));
part2 = find(part2);
part2 = part2(:);
cut = spaugment(A(part1,part2),1);
cutxy = xy([part1; part2],:);

clf reset
colordef(gcf,'black')
if size(xy,2) == 2
    axis([min(xy(:,1)) max(xy(:,1)) min(xy(:,2)) max(xy(:,2))]);
else
    axis([min(xy(:,1)) max(xy(:,1)) min(xy(:,2)) max(xy(:,2)) ...
        min(xy(:,3)) max(xy(:,3))]); 
end;

axis equal;
axis off;

hold on
if length(part1) > 0 & length(part2) > 0
    set(gplotg(cut,cutxy,'-'),'color',color3);
    xlabel([int2str(cutsize(A,part1)) ' cut edges'],'visible','on');
else
    xlabel('0 cut edges','visible','on');
end;
if length(part1) > 0 
    set(gplotg(A(part1,part1),xy(part1,:),'-'),'color',color1);
    if n < 500,
        if size(xy,2) == 2
            set(plot(xy(part1,1),xy(part1,2),'o'),'color',color1);
        else
            set(plot3(xy(part1,1),xy(part1,2),xy(part1,3),'o'),'color',color1);
        end;
    end;
end;
if length(part2) > 0 
    set(gplotg(A(part2,part2),xy(part2,:),'-'),'color',color2);
    if n < 500,
        if size(xy,2) == 2
            set(plot(xy(part2,1),xy(part2,2),'o'),'color',color2);
        else
            set(plot3(xy(part2,1),xy(part2,2),xy(part2,3),'o'),'color',color2);
        end;
    end;
end;
hold off

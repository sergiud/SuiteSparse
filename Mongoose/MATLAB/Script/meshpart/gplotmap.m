function handle = gplotmap(A,xy,map,pcolor,ecolor)
% GPLOTMAP : Plot a partitioned graph in 2 or 3 dimensions.
%
%	gplotmap(A,xy,map) plots the n-vertex graph specified 
%       by the n by n adjacency (or Laplacian) matrix A 
%       and the n by 2 or 3 matrix of coordinates xy.
%       Argument map is an n-vector of part numbers, which
%       assigns each vertex to a part.
%
%       By default, edges that join different parts are omitted, and
%       the picture shows each part in a different color.  The call
%
%       gplotmap(A,xy,map,pcolor,ecolor)
%
%       uses color "pcolor" for the parts and "ecolor" for the
%       edges between the parts.  If "pcolor" has multiple rows,
%       each part gets the color of one row (in cyclic order).
%
% See also GPLOTPART.
%
% John Gilbert  2 August 1994
% Copyright (c) 1990-1996 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.

[n,n] = size(A);
if nargin < 3, 
    map = zeros(1,n);
end;
if nargin < 4 
    pcolor = [];
end;
if nargin < 5,
    ecolor = [];
end;

parts = setfilter(map);
nparts = length(parts);
if length(pcolor) == 0,
    pcolor = hsv(nparts);
    pcolor = pcolor(randperm(nparts),:);
end;

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

% Count and plot the separating edges.
[i,j] = find(A);
f = find(map(i) > map(j));
if length(f) 
    xlabel([int2str(length(f)) ' cut edges'],'visible','on');
    if length(ecolor)
        cut = sparse(i(f),j(f),1,n,n);
        set(gplotg(cut,xy,'-'),'color',ecolor);
    end;
else
    xlabel('0 cut edges','visible','on');
end;

% Plot each piece.
ncolor = 1;
for partnumber = parts;
    c = pcolor(ncolor,:);
    ncolor = rem(ncolor, size(pcolor,1)) + 1;
    part = find(map == partnumber);
    set(gplotg(A(part,part),xy(part,:),'-'),'color',c);
    if n < 500,
        if size(xy,2) == 2
            set(plot(xy(part,1),xy(part,2),'o'),'color',c);
        else
            set(plot3(xy(part,1),xy(part,2),xy(part,3),'o'),'color',c);
        end;
    end;
end;

hold off

function highlight(A,xy,sep,highcolor,meshcolor,dotsize)
% HIGHLIGHT : Plot a mesh with subgraph highlighted.
%
%  highlight(A,xy,sep) plots a picture of the mesh A with coordinates xy,
%  highlighting the subgraph induced by the vertices in sep.
%  Optional fourth and fifth arguments are color/linetype of highlight and mesh.
%  Optional sixth argument is size of highlighting dot.
%
% John Gilbert and Shanghua Teng, 1992-1993.
% Copyright (c) 1990-1996 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.


if nargin < 4
    highcolor = 'w-';
end
if nargin < 5
    meshcolor = 'r-';
end
if nargin < 6
    dotsize = 15;
end

[n,n] = size(A);
[i,j] = find(A);

% Plot main graph with vertices before edges.

[ignore, p] = sort(max(i,j));
i = i(p);
j = j(p);

% Create a long, NaN-seperated list of line segments,
% rather than individual segments.

X = [ xy(i,1) xy(j,1) NaN*ones(size(i))]';
Y = [ xy(i,2) xy(j,2) NaN*ones(size(i))]';
X = X(:);
Y = Y(:);

xymin = min(xy);
xymax = max(xy);
plot (X, Y, meshcolor,'erasemode','none');
axis([xymin(1) xymax(1) xymin(2) xymax(2)]);
axis('equal');
axis('off');
hold on;

% Highlight sep set.

B = A(sep,sep);
xB = xy(sep,1);
yB = xy(sep,2);
[i,j] = find(B);
X = [xB(i) xB(j)]';
Y = [yB(i) yB(j)]';
plot (X, Y, highcolor,'erasemode','none');
if n < 1200
   handle = plot(xB,yB,[highcolor(1) '.'],'erasemode','none');
   set(handle,'markersize',dotsize)
end

hold off;

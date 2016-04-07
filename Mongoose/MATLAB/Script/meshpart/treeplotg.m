function treeplotg(p,b,c,d)
% TREEPLOTG : Plot a picture of a tree.
%	TREEPLOTG(p,b,c,d)
%	TREEPLOTG(p,c,d)
%	    p is a vector of parent pointers, with p(i) == 0 for a root.
%	    b > 0 labels nodes with numbers in the given font size. 
%	    c is a color and character for nodes, or '' to not plot nodes.
%	    d is a color and character for edges, or '' to not plot edges.
%	    b, c, and d may be omitted, and reasonable defaults are used.
%
%	See also TREEPERM, ETREE, TREELAYOUTG, ETREEPLOTG.
%
% Copyright (c) 1990-1996 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.
% Copyright (c) 1984-92 by The MathWorks, Inc.
%
% Written 1991 by John Gilbert.
% Modified 23 Jun 93 by JRG to turn the axes off.
% Modified  1 Aug 94 by JRG to label the nodes.
% Modified 29 Dec 94 by JRG to fix a marker-symbol bug.
% Modified by Tim Davis, for Matlab 5.1.  July 6, 1998.
% Modified 14 Jun 01 by JRG to restore node numbering
% Modified by John Gilbert for Matlab 6, Feb 2002

n = max(size(p));

% Figure out which argument is which, and set defaults.
if nargin == 1
    b = 0;
    c = 'r';
    d = 'r';
elseif nargin == 2
    if isstr(b)
        c = b;
        d = c; 
        b = 0;
    else
        c = 'r';
        d = 'r';
    end;
elseif nargin == 3
    if isstr(b)
        d = c;
        c = b;
        b = 0;
    else
        d = c;
    end;
end;

% Place the nodes.
[x,y,h]=treelayoutg(p);

% Set up NaN-separated lists of edges.
f = find(p~=0);
pp = p(f);
X = [x(f); x(pp); NaN*ones(size(f))];
Y = [y(f); y(pp); NaN*ones(size(f))];
X = X(:);
Y = Y(:);

% Make sure the color and symbol arguments are clean.
ValidMarkers = '.xo*+adv^<>ph';
[NodeStyle,NodeColor] = colstyle(c);
if (isempty (findstr(ValidMarkers,NodeStyle)))
    NodeStyle = 'o';
    DrawNodes = n < 500  &  ~isempty (c) ;
else
    DrawNodes = 1;
end;
if (isempty (NodeColor))
    NodeColor = 'r';
end;
[EdgeStyle,EdgeColor] = colstyle(d);
if (isempty (EdgeStyle) | ~isempty (findstr(ValidMarkers,EdgeStyle)))
    EdgeStyle = '-';
    DrawEdges = ~isempty (d) ;
else
    DrawEdges = 1;
end;
if (isempty (EdgeColor))
    EdgeColor = NodeColor;
end;
NodeSpec = [NodeColor NodeStyle];
EdgeSpec = [EdgeColor EdgeStyle];

% Draw the plot.
%if DrawNodes & DrawEdges
%    plot (x, y, NodeSpec, X, Y, EdgeSpec);
%elseif DrawNodes
 %   plot (x, y, NodeSpec);
 %elseif DrawEdges
%else
%end;
if DrawEdges & DrawNodes
    plot(X, Y, EdgeSpec, x, y, NodeSpec, 'MarkerSize', 4);
elseif DrawEdges
    plot (X, Y, EdgeSpec);
elseif DrawNodes
    plot (x, y, NodeSpec, 'MarkerSize', 4);
else
end;

if b 
% Label nodes with node numbers, in font size b (or default to 6 if b==1).
    if b == 1
        b = 6;
    end;
    % Don't label the nodes in the middle of a vertical chain.
    f = find((x(1:n-2) ~= x(2:n-1)) | (x(2:n-1) ~= x(3:n)));
    f = [1 f+1 n];
    nf = length(f);
    numbers = zeros(nf,length(int2str(n)));
    for k = 1:nf;
        nk = int2str(f(k));
        numbers(k,1:length(nk)) = nk;
    end;
    th = text(x(f)+.01, y(f), char(numbers));
    set(th,'fontsize',b);
 end;


axis([0 1 0 1]);
axis('off');
xlabel(['height = ' int2str(h)],'visible','on');

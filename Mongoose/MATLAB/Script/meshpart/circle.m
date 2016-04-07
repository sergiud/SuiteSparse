function [xx,yy] = circle(x,y,radius,arg4)
% CIRCLE : Plot a circle.
%
% circle(x, y, radius,'c') plots a circle with specified center and radius, 
%       in color 'c' (optional).
% [xx,yy] = circle(x, y, radius) returns coords instead of plotting;
%       in this case the optional fourth argument is # of points to use.
%
% If x, y, and radius are vectors then one circle is plotted for each element.
%
% John Gilbert, Xerox PARC, 1992.
% Copyright (c) 1990-1996 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.

c = 'r-';
npoints = 100;
if nargin >= 4
    if nargout < 1
        c = arg4;
    else
        npoints = arg4;
    end
end

% Make the arguments column vectors
x = x(:);
y = y(:);
radius = radius(:);

theta = [0:(npoints-1) 0] * 2 * pi / npoints;
xx = radius * cos(theta) + x * ones(1,npoints+1);
yy = radius * sin(theta) + y * ones(1,npoints+1);
xx = xx';
yy = yy';

if nargout < 1
    plot(xx,yy,c)
end

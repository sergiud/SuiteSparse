function [center,radius] = mscircle(v,cpt,xyz)
% MSCIRCLE : Describe a separating circle in mesh space coordinates.
%
% [center,radius] = mscircle(v,cpt,xyz)
% This routine, which is only used to draw a picture and for optional return
% values from geopart, converts the description of the separator in 
% dim+1-space to a description of a circle in dim-space in mesh coordinates.
%
%   v: the normal vector to the separating plane in dim+1 conformal map space.
% cpt: the centerpoint for the conformal mapping.
% xyz: the points themselves.
%
% See GEOPART.
%
% John Gilbert and Shanghua Teng, 1992-1993.
% Copyright (c) 1990-1996 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.

dim = length(v)-1;

% Make the normal vector a unit length column vector.
v = v(:)/norm(v);

% First get dim+1 points on a unit circle in conformally mapped space,
% with center at the origin, normal to the direction of v.
% We will take d-1 of the points to be orthogonal,
% namely a basis for the null space of v.
 
circ = null(v');
circ = circ';
s = -sum(circ);
circ = [circ ; s/norm(s)];

% Move the great circle over on the unit sphere to balance the cut.
    
% In conformally mapped space, 
% find out how far from the origin the separating plane is.

d = median(conmap(cpt,xyz) * v);

% In conformally mapped space, 
% find a circle through a plane normal to v, at distance d from the origin.

s = sqrt(1-d^2);
circ = s * circ + d * ones(size(v))*v';

% Now each of the dim+1 rows of circ is a point on the separating circle
% in conformally mapped space.  Map them back to mesh space.

circ = conmapinv(cpt,circ);
circ = stereodown(circ);

% Now each of the dim+1 rows of circ is a point on the separating circle
% (or sphere) in mesh space.  Figure the center and radius.

% The points x equidistant from points a and b are those in
% the plane of the perpendicular bisector of a and b, which
% is to say (x - (a+b)/2) is orthogonal to a-b.
% We can rewrite that as 2*(a-b)' * x = a'*a - b'*b.
% This gives us dim equations for the center x, taking a to
% be any of the first dim points on the circle and b to be the last.

M = 2 * ( circ(1:dim,:) - ones(dim,1)*circ(dim+1,:) );
circnorms = sum((circ .* circ)')';
q = circnorms(1:dim) - ones(dim,1)*circnorms(dim+1);

% Now the center is the solution to M*x = q:

center = M \ q;

% The radius is the distance from x to the first point.

radius = norm(center'-circ(1,:));

function [part1,part2,sep1,sep2,center,radius] = ...
            geopart(A,xy,picture,ntries)
% GEOPART : Geometric separators for a finite element mesh.
%
% [part1,part2] = geopart(A,xy) returns the vertex partition from geometric
% partitioning for the mesh whose structure is A and whose coordinates are xy.
%
% The number of dimensions is the number of columns of xy.
%
% [part1,part2,sep1,sep2,center,radius] = geopart(...) also returns:
%        sep1, sep2 are the vectors of endpoints of the separating edges.
%        center, radius describe the separating circle in the mesh space.
%        (If the separating circle is actually a line/plane, then radius=Inf
%        and the line's equation is center(1:d) * t = center(d+1). )
%
% There are two optional arguments:
% geopart(A,xy,picture,ntries)
%    picture (default 0) : 1 means draw a picture.
%                          2 means draw a picture and print statistics.
%                         -1 means just print statistics.
%    ntries (default 30) : Number of separating great circles to try.
%
% If there is no output argument, geopart always draws a picture.
%
% See also GEOSEP, GEODICE, GEOND.
%
% John Gilbert and Shanghua Teng, 1992-1993.
% Copyright (c) 1990-1996 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.

[npoints,dim] = size(xy);

if nargin < 3 
    picture = 0;
end;
stats = (picture < 0) | (picture > 1);
picture = (picture > 0) | (nargout == 0);
if nargin < 4
    ntries = 30;
end;

% The "ntries" tries will include "nouter" centerpoint computations,
% each with "ninner" hyperplanes; 
% and also "nlines" cut planes in the dim-dimensional mesh space.
% The following division is ad hoc but seems about right.

nlines = floor((ntries/2) ^ (dim/(dim+1)));
nouter = ceil(log(ntries-nlines+1)/log(20));
ninner = floor((ntries-nlines) / nouter);
nlines = ntries - nouter*ninner;
if stats
    lines_outer_inner = [nlines nouter ninner]
end;

% Size of sample for approximate centerpoint computation.
csample = min(npoints,(dim+3)^4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here's the beginning of the partitioning computation.  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Scale points to have mean zero and max coordinate magnitude 1.
xyshift = mean(xy);
xy = xy - ones(npoints,1)*xyshift;
xyscale = max(max(abs(xy)));
if xyscale == 0
    error('All points have the same coordinates.');
end;
xy = xy / xyscale;

% Project points stereographically up to the sphere.
xyz = stereoup(xy);

% Compute "nouter" approximate centerpoints.
circlequality = Inf;
for i = 1:nouter

    if stats, disp(sprintf('Starting outer iteration %d.',i)); end;

    cpt = centerpoint(xyz,csample);
    if 1-norm(cpt) < sqrt(eps)
        % Centerpoint is on sphere (probably due to duplicate input points);
        % move it in a little.
        dc = randn(size(cpt));
        cpt = cpt * .9 + dc/(20*norm(dc));
    end;
    xyzmap = conmap(cpt,xyz);
    [greatcircle,gcquality] = sepcircle(A,xyzmap,ninner);
    if gcquality < circlequality
        circlequality = gcquality;
        bestcircle = greatcircle;
        bestcpt = cpt;
    end;

end;

if stats, disp('Finished last outer iteration.'); end;

% Also try separating with a straight line in mesh space.
if nlines >= 1
    [bestline,linequality] = sepline(A,xy,nlines);
else 
    linequality = Inf;
end;

if linequality < circlequality
    [part1,part2] = partition(xy,bestline);
else
    [part1,part2] = partition(conmap(bestcpt,xyz),bestcircle);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here's the end of the partitioning computation.        %
% The rest of the mess is just for the optional return   %
% values, and to draw the picture.                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargout > 2 | stats
    % Find the separating edges.
    [sep1,sep2] = find(A(part1,part2));
    sep1 = part1(sep1);
    sep2 = part2(sep2);
end;

if nargout > 3 | picture 
    if linequality < circlequality
        % Describe the separating line in the mesh space.
        radius = Inf;
        v = bestline/norm(bestline);
        b = median(xy * bestline);
        % Now the line's equation is v'*p == b.
        % Put it back in the original coordinate system.
        b = xyscale*b + xyshift*v;
        center = [v' b];
    else
        % Describe the separating circle in the mesh space.
        [center, radius] = mscircle(bestcircle,bestcpt,xyz);
        % Put it back in the original coordinate system.
        center = center*xyscale + xyshift';
        radius = radius*xyscale;
    end;
end;

if stats
    bestline_bestcircle = [linequality circlequality]
    na_nb_ne = [length(part1) length(part2) length(sep1)]
end;
if picture 
    % Plot the points in the original coordinate system.
    xy = xy * xyscale + ones(npoints,1)*xyshift;
    gplotpart(A,xy,part1);
    title('Geometric Partition');
end;
if picture & dim == 2
    % Draw the separating line or circle.
    hold on;
    if radius == Inf
        v = center(1:2)';
        b = center(3);
        p = xyshift' + (b-v'*xyshift')*v;
        q = 2*xyscale*null(v');
        t = [p + q , p - q];
        plot(t(1,:),t(2,:),'w');
    else
        circle(center(1),center(2),radius,'w');
    end;
    hold off;
end;
if picture 
    drawnow;
end;

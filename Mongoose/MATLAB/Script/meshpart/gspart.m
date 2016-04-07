function [p,uv] = gspart(A,picture,ntries,uv);
% GSPART : Geometric/spectral partition of a graph.
%
% p = gspart(A) returns a list of the vertices on one side of a partition
%     obtained by the geometric separator algorithm applied to the graph
%     of A with coordinates chosen from the eigenvectors u2, u3, ...
%     of the Laplacian.
%
% There are three optional arguements:
% gspart(A,picture,ntries,uv)
%   picture (default 0): If nonzero scalar, draw a picture in spectral coords.
%                        If array, use as xy coords for another picture.
%   ntries (default 50): Number of trials for geometric routine.
%   uv (default 2):      If scalar, use this many eigenvectors.
%                        If array, use uv as the matrix of eigenvectors.
% [p,uv] = gspart(...)   also returns the two eigenvectors.
%
% If there is no output argument, gspart always draws a picture.
%                 
% See also SPECPART 
%
% John Gilbert  31 Jan 1994
% Copyright (c) 1990-1996 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.
%
% Modified by Tim Davis, for Matlab 5.1.  July 6, 1998.

% Figure out what the arguments mean...

if nargin < 4
    [uv,lambda] = fiedler(A,2);
end;
if length(uv) == 1
    k = uv;
    [uv,lambda] = fiedler(A,k);
else
    k = size(uv,2);
end;

if nargin < 3
    ntries = [];
end;
if (isempty (ntries))
    ntries = 50;
end;

if nargin < 2
    picture = [];
end;
if (isempty (picture))
    picture = 0;
end;
if length(picture) == 1
    xypicture = 0;
    uvpicture = picture | (nargout == 0);
else
    xypicture = 1;
    uvpicture = 1;
    xy = picture;
end;

% Enough with the arguments, already!

if uvpicture
    uvfigure = gcf;
    xyfigure = gcf+1;
    clf reset
end;
p = geopart(A,uv,uvpicture,ntries);
if uvpicture
    title('Geometric Spectral Partition')
end;
if xypicture
    figure(xyfigure)
    clf reset
    gplotpart(A,xy,p);
    title('Geometric Spectral Partition')
end;
%if uvpicture & xypicture & k == 2
%    % Highlight a couple of corresponding vertices in each picture.
%    uvmean = mean(uv);
%    uvdist = (uv(:,1)-uvmean(1)).^2 + (uv(:,2)-uvmean(2)).^2;
%    [ignore,i] = max(uvdist);
%    xymean = mean(xy);
%    xydist = (xy(:,1)-xymean(1)).^2 + (xy(:,2)-xymean(2)).^2;
%    [ignore,j] = max(xydist);
%    figure(uvfigure)
%    hold on
%    drawpoints(uv(i,:),'m.',20);
%    drawpoints(uv(j,:),'w.',20);
%    hold off
%    figure(xyfigure)
%    hold on
%    drawpoints(xy(i,:),'m.',20);
%    drawpoints(xy(j,:),'w.',20);
%    hold off
%end;
if uvpicture
    figure(uvfigure)
end;

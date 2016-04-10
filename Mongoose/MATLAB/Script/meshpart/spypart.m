function spypart(S, rp, cp);
%SPYPART Spy plot with partitioning.
%   SPYPART(S,rp,cp) plots the sparsity pattern of a matrix S,
%   with lines marking a block partition described by 
%   rp (rows) and cp (columns).
%   If S is square, cp may be omitted and defaults to rp.
%
%   Partitions are specified as in the output of DMPERM:
%   There are length(rp)-1 row blocks, of sizes diff(rp), with rp(1)=1.
%
% See also DMSPY, DMPERM.

% John Gilbert, 1990-2002.
%   Copyright 1984-2000 The MathWorks, Inc. 
%   $Revision: 5.5 $  $Date: 2000/06/01 03:46:26 $

if (nargin < 3), cp = rp; end

clf
colordef(gcf,'black')
[m,n] = size(S);
if max(m,n) > 100 
    spy(S,3)
else
    spy(S,5)
end
hold on

if length(rp) > 2
   k = length(rp)-2;
   X = [zeros(1,k); n+ones(1,k)];
   Y = rp(2:k+1) - 0.5;
   Y = [Y; Y];
   plot(X,Y,'w-')
end
if length(cp) > 2
   k = length(cp)-2;
   X = cp(2:k+1) - .5;
   X = [X; X];
   Y = [zeros(1,k); m+ones(1,k)];
   plot(X,Y,'w-') 
end
axis('ij')
hold off
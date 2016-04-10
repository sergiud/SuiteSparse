%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [q,b,Cout] = gp_interpretPart(A, part)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sanity Checks
%try
%if size(A,1) ~= size(A,2)
%    gp_log('Error: The given matrix is rectangular!\n') ;
%    size(A)
%    return
%end
%if norm(A-A') > 1e-10
%    gp_log('Error: The given matrix is not symmetric!\n') ;
%    A-A'
%    return
%end
%catch exception
%    gp_log('Warning: Out of Memory when interpreting partition result (sanity checks).\n') ;
%end

% Get the permutation after doing the partition on the vertices
[ignore, p] = sort(part) ;

% Permute the matrix into the fancy  [l q]
%                                    [q l]
C = A(p,p) ;
if nargout >= 3, Cout = C; end

% Sum the edge weights that fall in the edge cut ('q' in the above sketch).
nleft = length(find(part == 1)) ;
nright = length(find(part == 0)) ;

q = full(sum(sum(C(1:nleft, nleft+1:end)))) ;

% The imbalance is the difference between the largest component's ratio and 0.50,
% the perfect cut. For example, a 60:40 cut results in an imbalance of 0.10.
b = (max(nleft,nright) / (nleft+nright)) - 0.50 ;

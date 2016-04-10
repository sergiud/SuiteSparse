%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,ex] = gp_getLargestComponent(A)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ex = 0 ;

% Cut the largest component using dmperm:
try
    A = A + speye(size(A)) ;
    [p,q,r,s] = dmperm(A) ;
    A = A(p,q) ;
    [~, g] = max(diff(r)) ;
    st = r(g) ;
    en = r(g+1)-1 ;
    A = A(st:en, st:en) ;
catch
    ex = 1 ;
    gp_log('Error: Out of Memory when selecting largest component for experiment.\n') ;
    return ;
end

% Remove the diagonal.
A = tril(A,-1) + triu(A,1) ;

if nnz(A) == 0
    ex = 1;
    gp_log('Error: nnz(A) == 0.\n');
    return ;
end

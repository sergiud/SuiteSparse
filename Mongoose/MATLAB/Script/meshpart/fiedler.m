function [x,lambda] = fiedler(G,k,verbose);
% FIEDLER : Fiedler vector of a graph.
%
% x = fiedler(G)  returns the eigenvector corresponding to 
%                 the second-smallest eigenvalue of the laplacian matrix 
%                 of the graph whose adjacency matrix is G|G'.
%                 (If G is connected this is the smallest positive eval.)
% [x,lambda] = fiedler(G)  also returns the eigenvalue.
% [x,lambda] = fiedler(G,k) returns the k smallest nonzero eigenvalues
%                 and their eigenvectors.
% [x,lambda] = fiedler(G,k,verbose) prints some trace information:
%                 verbose = 0 for none, 1 for some, 2 for more.
%
% If G is full or smaller than 100 vertices we use Matlab's full "eig";
% if G is sparse we use an iterative method.
% This is *not* an efficient implementation, it's just for playing with
% small examples.
%
% See also LAPLACIAN, SPECPART.
%
% John Gilbert  2 February 1994
% Copyright (c) 1990-1996 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.

if nargin < 3, verbose = 0; end;
babble = verbose > 1;
if nargin < 2, k = 1; end;

[n,n] = size(G);
if n < 100, 
    G = full(G); 
end;

comp = components(G);
if max(comp) > 1  % Add edges to connect the graph.
    if verbose
        disp('Fiedler:  Graph is disconnected.  Adding edges.');
        [x,ignore] = hist(comp,max(comp));
        Sizes = x
    end;
    for c = 1:max(comp)
        i(c) = min(find(comp==c));
    end;
    S = sparse(i(1),i(2:max(comp)),1,n,n);
    G = G | (S | S');
end;

L = laplacian(G);

if babble
    if issparse(L), ss = 'sparse'; else ss = 'full'; end;
    fprintf(1,'Fiedler: %d vertices, %d edges, %s adjacency matrix.\n',...
               n, (nnz(L)-n)/2, ss);
end;

if issparse(G)
    x = ones(n,1)/sqrt(n);
    lambda = 1;
    shift = -eps^.25;
    I = speye(n);
    p = symmmd(L);
    L = L(p,p);
    R = chol(L-shift*I);
    Rt = R';
    for j = 1:k
        if verbose
            disp(['Eigenvector u' int2str(j+1) '...']);
        end;
        % Append the j+1'st eigenvector to x.
    
        % Use a few power iterations with a small negative shift first.
        if babble, fprintf(1,'Starting %d power iterations.\n',20*j);end;
        y = rand(n,1);
        y = y - x*(x'*y);
        for i = 1:20*j
            y = R\(Rt\y);
            y = y - x*(x'*y);
            y = y/norm(y);
        end;

        % Now use Rayleigh quotient iteration with a variable shift.
        if babble, disp('Starting Rayleigh quotient iterations.');end;
        tol = eps^.25;
        maxiter = 10;
        iteration = 0;
        done = 0;
        while ~done,
            iteration = iteration+1;
            estimate = y'*L*y;
            if babble, iteration, end;
            if babble, estimate, end;
            newy = (L - estimate*I) \ y;
            if any(isnan(newy)|isinf(newy))
                if verbose
                    disp(['Rayleigh iteration halted at iteration ' ...
                    int2str(iteration)]);
                end;
                done = 1;
            else
                newy = newy/norm(newy);
                measure = norm(newy-y)*norm(newy+y);
                done = measure < tol;
                y = newy;
            end;
            % Project y to be orthogonal to earlier eigenvectors.
            y = y - x*(x'*y);
            y = y/norm(y);
            if iteration == maxiter
                done = 1;
                disp('Rayleigh iteration reached maxiter.');
            end
        end;
        x = [x y];
        lambda = [lambda estimate];
    end;
    x = x(:,2:k+1);
    x(p,:) = x;
    lambda = lambda(2:k+1);

else

    % Full graph -- just use eig.
    [V,D] = eig(L);
    [d,i] = sort(diag(D));
    lambda = d(2:k+1);
    x = V(:,i(2:k+1));

end;

if babble, lambda, end;

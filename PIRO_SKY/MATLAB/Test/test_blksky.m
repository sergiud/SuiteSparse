function test_blksky(N, cond)
% Usage :
%           test_blksky
%           test_blksky(N, cond)
% Get first N square matrices that satisfy the condition cond from the UF sparse
% matrix collection. The default is to get
% 20 square matrices with dimension between 100 and 500. 
% Compute the numerical factorization of the R-bidiagoalization, after a QR
% factorization of A. We currently use a static data structure for R, 
% and use MATLAB full svd to compare the results.
%
% TBD : Rectangular matrices, fix the bug in the dynamic data
% structure case for R (opt == 3).

doplot = 1 ;
debug = 1 ;

if (nargin < 1)
    N = 40 ;
end
if (nargin < 2)
    cond =  'find ((index.nrows == index.ncols) & (index.nrows < 10))'
    %cond =  'find ((index.nrows == index.ncols) & (index.nrows >= 10) & (index.nrows <=100) )'
    %cond =  'find ((index.nrows == index.ncols) & (index.nrows >= 100) & (index.nrows <=500) )'
    %cond =  'find ((index.nrows == index.ncols) & (index.nrows > 500) & (index.nrows <=1000) )'
    %cond =  'find ((index.nrows == index.ncols) & (index.nrows > 1000) & (index.nrows <=2000) )'
    % cond =  'find ((index.nrows == index.ncols) & (index.nrows > 2000) & (index.nrows <=3500) )'
end

index = UFget ;
ids = eval(cond) ;
[ignore, i] = sort (index.nnz (ids)) ;
ids = ids (i) ;

fprintf('relative error is (sgood-snow/sgood(1)) where sgood is from \n')
fprintf('MATLAB full SVD, snow from the SVD of the bidiagonalized matrix\n\n')
fprintf('ID \t Size \t (relative error)\n\n') 
[junk, idsize] = size(ids) ;
least = min (idsize, N) ;

for id = ids(1:least)
    Prob = UFget (id) ;   % Prob is a struct (matrix, name, meta-data, ...)
    A = Prob.A ;          % A is a symmetric sparse matrix
    %pause

    if (~isreal(A))
        continue ;
    end

    [m, n] = size(A) ;

    % What is the best ordering ?
    if (debug)
    fprintf('Find ordering \n')
    end
    tic 
    %q = colamd(A) ;
    q = symrcm(A'*A) ;
    t1 = toc ;
    %A1 = A(:, q) ;
    A1 = A(q, q) ;


    % Find "Q-less QR decomposition" for now.
    if (debug)
    fprintf('Find qr \n')
    end
    tic
    R = qr(A1) ;
    t2 = toc ;

    if (nnz(diag(R)) < n)
        continue ;
    end

    if (doplot)
        clf ;
        figure(1) ;
        hold on ;

        subplot(2, 2, 1) ;
        spy(A) 
        title('A') 

        subplot(2, 2, 2) ;
        spy(A1) 
        title('A1') 

        subplot(2, 2, 3) ;
        spy(R) 
        title('R') 
        drawnow
    end

    % fprintf('Calling skyline\n') ;
    %Prob.id
    %pause
    if (debug)
    fprintf('Find bandwidth \n')
    end
    tic
    bw = find_bw(R, 1) ;
    t_findbw = toc ;

    if (bw == 0)
        continue ;
    end

    if (debug)
    fprintf('Find R-bidiag, bw %g \n', bw)
    end
    tic
    [b1, b2] = blksky(R, 1, bw) ;
    t3 = toc ;

    newtime =  t1 + t2 + t3 + t_findbw ;

    % find the svd of B
    if (debug)
    fprintf('Find bidiag SVD \n')
    end
    B = zeros(n) ;
    B = spdiags([b1, b2], 0:1, B) ;
    snow = svd(full(B)) ;


    % fprintf('Calling matlab svd\n') ;
    if (debug)
    fprintf('Find SVD \n')
    end
    tic
    sgood = svd(full(A1)) ;
    t4 = toc ;

    maxerr = 0.0 ;
    for i = 1 : n
        err = abs(sgood(i) - snow(i)) ;
        maxerr = max(maxerr, err) ;
    end

    fprintf('%d %d-by-%d %g\n', Prob.id, m, n, maxerr/sgood(1)) 

end


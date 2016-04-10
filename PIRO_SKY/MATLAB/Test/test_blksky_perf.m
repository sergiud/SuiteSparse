function test_blksky_perf(N, cond)
% Usage :
%           test_blksky_perf
%           test_blksky_perf(N, cond)
% Get first N square matrices that satisfy the condition cond from the UF sparse
% matrix collection. The default is to get
% 20 square matrices with dimension between 100 and 500. 
% Compute the numerical factorization of the R-bidiagoalization, after a QR
% factorization of A. We currently use a static data structure for R, 
% and compare the time to MATLAB svds.
%
% TBD : Rectangular matrices, fix the bug in the dynamic data
% structure case for R (opt == 3).

doplot = 0 ;

if (nargin < 1)
    N = 100 ;
end
if (nargin < 2)
     % cond =  'find ((index.nrows == index.ncols) & (index.nrows >= 100) & (index.nrows <=500) )'
     % cond =  'find ((index.nrows == index.ncols) & (index.nrows >= 400) & (index.nrows <=600) )'
      cond =  'find ((index.nrows == index.ncols) & (index.nrows >= 1000) & (index.nrows <=1500) )'
     %cond =  'find ((index.nrows == index.ncols) & (index.nrows >= 3000) & (index.nrows <=3500) )'
end

index = UFget ;
ids = eval(cond) ;
[ignore, i] = sort (index.nnz (ids)) ;
ids = ids (i) ;

fprintf('Relative time is (svds time)/(colamd + qr + R-bidiag time) \n')
fprintf('R-bidiag will use a static data structure for this tests \n') 
fprintf('These tests will be slow because of svds(). \n\n\n') 
fprintf('ID \t Size \t (relative time)\n\n') 

for id = ids(1:N)
    Prob = UFget (id) ;   % Prob is a struct (matrix, name, meta-data, ...)
    A = Prob.A ;          % A is a symmetric sparse matrix
    %pause

    if (Prob.id == 1466)
        continue ;
    end

    [m, n] = size(A) ;

    % What is the best ordering ?
    tic 
    % q = colamd(A) ;
    q = symrcm(A'*A) ;
    t1 = toc ;
    %A1 = A(:, q) ;
    A1 = A(q, q) ;


    % Find "Q-less QR decomposition" for now.
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
        pause
    end

    % fprintf('Calling skyline\n') ;
    tic
    [b1, b2] = blksky(R) ;
    t3 = toc ;
    newtime =  t1 + t2 + t3 ;

    % fprintf('Calling matlab svd\n') ;
    tic
    sgood = svds(A1, floor(n/8)) ;
    t4 = toc ;

    % fprintf('Calling matlab svd\n') ;
    tic
    sgood = svd(full(A1)) ;
    t5 = toc ;

    fprintf('%d \t %d-by-%d \t %g %g \n',Prob.id, m, ...
            n, t4/newtime, t5/newtime) 


end


function test_blksky_symbolic(N, cond)
% Usage :
%           test_blksky_symbolic
%           test_blksky_symbolic(N, cond)
% Get first N square matrices that satisfy the condition cond from the UF sparse
% matrix collection. The default is to get
% 20 square matrices with dimension between 100 and 500. Run symbolic 
% factorization for R-bidiagonalization using the blocked rightmost corner
% method.

doplot = 0 ;

if (nargin < 1)
    N = 20 ;
end
if (nargin < 2)
     cond =  'find ((index.nrows == index.ncols) & (index.nrows >= 100) & (index.nrows <=500) )'
     % cond =  'find ((index.nrows == index.ncols) & (index.nrows >= 400) & (index.nrows <=600) )'
     % cond =  'find ((index.nrows == index.ncols) & (index.nrows >= 3000) & (index.nrows <=3500) )'
end

index = UFget ;
ids = eval(cond) ;
[ignore, i] = sort (index.nnz (ids)) ;
ids = ids (i) ;

fprintf('ID \t Size \t nnz(A) \t nnz(R) \t (skyline size/nnz(R)) \t #givens \t #swaps \t flops \n\n') 

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
    q = colamd(A) ;
    t1 = toc ;
    A1 = A(:, q) ;


    % Find "Q-less QR decomposition" for now.
    tic
    R = qr(A1) ;
    t2 = toc ;

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

    [rR, cR] = size(R)  ;
    [toprow, ng, nswap, fl] = blksky_symbolic(R) ;

    csize = 0 ;
    for k = 1 : cR
        csize = csize + k - toprow(k) + 1 ;
    end

    if (nnz(R) ~= 0)
        fprintf('%d \t %d-by-%d \t %d \t\t %d \t\t %0.4f \t %g \t\t %g \t\t %g\n',...
        Prob.id, m, n, nnz(A), nnz(R), csize/nnz(R), ng, nswap, fl) 
    else
        fprintf('%d \t %d-by-%d \t %d \t\t %d \t\t %g \t\t %g \t\t %g\n',Prob.id, m, ...
                n, nnz(A), nnz(R), ng, nswap, fl) 
    end
end


load sky_perf_compare_mp_aug23_2.mat
tid = v_id' ;
tn = v_n' ;
t_sky_1 = v_sky_1' ;
t_mat_n_1000 = v_mat_n_1000' ;
t_mat_n_500 = v_mat_n_500' ;
t_mat_n_100 = v_mat_n_100' ;
t_mat_n_50 = v_mat_n_50' ;
t_mat_n_25 = v_mat_n_25' ;
t_mat_n_10 = v_mat_n_10' ;
t_mat_n_5 = v_mat_n_5' ;
data  = [ tid tn t_sky_1 t_mat_n_1000 t_mat_n_500 t_mat_n_100 t_mat_n_50 t_mat_n_25 t_mat_n_10 t_mat_n_5] ;
[row, junk] = size(data) ;


N = 1000 ;
     % cond =  'find ((index.nrows == index.ncols) & (index.nrows >= 100) & (index.nrows <=500) )'
     % cond =  'find ((index.nrows == index.ncols) & (index.nrows >= 400) & (index.nrows <=600) )'
      cond =  'find ((index.nrows == index.ncols) & (index.nrows >= 1000) & (index.nrows <= 5000) )'
     %cond =  'find ((index.nrows == index.ncols) & (index.nrows >= 3000) & (index.nrows <=3500) )'

index = UFget ;
ids = eval(cond) ;
[ignore, i] = sort (index.nnz (ids)) ;
ids = ids (i) ;
[junk, idsize] = size(ids) ;
least = min (idsize, N) ;

fprintf('Relative time is (svds time)/(colamd + qr + R-bidiag time) \n')
fprintf('R-bidiag will use a static data structure for this tests \n') 
fprintf('These tests will be slow because of svds(). \n\n\n') 
fprintf('ID \t Size \t Sky Time \t n/1000 \t n/500\t  n/100\t  n/50\t  n/25 \t n/10\t  n/5\n\n') 

j1 = [ 1000 500 100 50 10 5 ] ;
invalid1 = 0 ;
invalid = -1 ;

for id = ids(1:least)
    if ((row < 1 | isempty(find (data(:, 1) == id ))) & id ~= 371)
    Prob = UFget (id) ;   % Prob is a struct (matrix, name, meta-data, ...)
    A = Prob.A ;          % A is a symmetric sparse matrix
    [m, n] = size(A) ;

    v_id = [v_id id] ;
    v_n = [v_n n] ;

    fprintf(' %d ', id) 
    fprintf(' %d ', n) 

    if (~isreal(A))
        v_sky_1 = [v_sky_1 invalid1 ] ;
	v_mat_n_1000 = [v_mat_n_1000 invalid1] ;
	v_mat_n_500 = [v_mat_n_500 invalid1] ;
	v_mat_n_100 = [v_mat_n_100 invalid1] ;
	v_mat_n_50 = [v_mat_n_50 invalid1] ;
	v_mat_n_25 = [v_mat_n_25 invalid1] ;
	v_mat_n_10 = [v_mat_n_10 invalid1] ;
	v_mat_n_5 = [v_mat_n_5 invalid1] ;
	save sky_perf_compare_mp_aug23_2.mat v_*
        fprintf(' \n ' ) 
        continue ;
    end

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
        v_sky_1 = [v_sky_1 invalid1 ] ;
	v_mat_n_1000 = [v_mat_n_1000 invalid1] ;
	v_mat_n_500 = [v_mat_n_500 invalid1] ;
	v_mat_n_100 = [v_mat_n_100 invalid1] ;
	v_mat_n_50 = [v_mat_n_50 invalid1] ;
	v_mat_n_25 = [v_mat_n_25 invalid1] ;
	v_mat_n_10 = [v_mat_n_10 invalid1] ;
	v_mat_n_5 = [v_mat_n_5 invalid1] ;
	save sky_perf_compare_mp_aug23_2.mat v_*
        fprintf(' \n ' ) 
        continue ;
    end

    tic
    bw = find_bw(R, 1) ;
    t_findbw = toc ;

    if (bw == 0)
        v_sky_1 = [v_sky_1 invalid1 ] ;
	v_mat_n_1000 = [v_mat_n_1000 invalid1] ;
	v_mat_n_500 = [v_mat_n_500 invalid1] ;
	v_mat_n_100 = [v_mat_n_100 invalid1] ;
	v_mat_n_50 = [v_mat_n_50 invalid1] ;
	v_mat_n_25 = [v_mat_n_25 invalid1] ;
	v_mat_n_10 = [v_mat_n_10 invalid1] ;
	v_mat_n_5 = [v_mat_n_5 invalid1] ;
	save sky_perf_compare_mp_aug23_2.mat v_*
        fprintf(' \n ' ) 
        continue ;
    end

    % fprintf('Calling skyline\n') ;
    tic
    [b1, b2] = blksky(R, 1, bw) ;
    t3 = toc ;
    newtime =  t1 + t2 + t3 + t_findbw ;
    fprintf(' %0.4f ', newtime) 
    v_sky_1 = [v_sky_1 newtime ] ;

    if (n <= 1000)
	fprintf(' %0.4f ', invalid) ;
	v_mat_n_1000 = [v_mat_n_1000 invalid] ;
    else
        tic
        sgood = svds(A1, floor(n/1000)) ;
        t4 = toc ;
	fprintf(' %0.4f ', t4) 
	v_mat_n_1000 = [v_mat_n_1000 t4] ;
    end

    if (n <= 500)
	fprintf(' %0.4f ', invalid) ;
	v_mat_n_500 = [v_mat_n_500 invalid] ;
    else
        tic
        sgood = svds(A1, floor(n/500)) ;
        t4 = toc ;
	fprintf(' %0.4f ', t4) 
	v_mat_n_500 = [v_mat_n_500 t4] ;
    end

    if (n <= 100)
	fprintf(' %0.4f ', invalid) ;
	v_mat_n_100 = [v_mat_n_100 invalid] ;
    else
        tic
        sgood = svds(A1, floor(n/100)) ;
        t4 = toc ;
	fprintf(' %0.4f ', t4) 
	v_mat_n_100 = [v_mat_n_100 t4] ;
    end

    if (n <= 50)
	fprintf(' %0.4f ', invalid) ;
	v_mat_n_50 = [v_mat_n_50 invalid] ;
    else
        tic
        sgood = svds(A1, floor(n/50)) ;
        t4 = toc ;
	fprintf(' %0.4f ', t4) 
	v_mat_n_50 = [v_mat_n_50 t4] ;
    end

    if (n <= 25)
	fprintf(' %0.4f ', invalid) ;
	v_mat_n_25 = [v_mat_n_25 invalid] ;
    else
        tic
        sgood = svds(A1, floor(n/25)) ;
        t4 = toc ;
	fprintf(' %0.4f ', t4) 
	v_mat_n_25 = [v_mat_n_25 t4] ;
    end

    if (n <= 10)
	fprintf(' %0.4f ', invalid) ;
	v_mat_n_10 = [v_mat_n_10 invalid] ;
    else
        tic
        sgood = svds(A1, floor(n/10)) ;
        t4 = toc ;
	fprintf(' %0.4f ', t4) 
	v_mat_n_10 = [v_mat_n_10 t4] ;
    end

    if (n <= 5)
	fprintf(' %0.4f ', invalid) ;
	v_mat_n_5 = [v_mat_n_5 invalid] ;
    else
        tic
        sgood = svds(A1, floor(n/5)) ;
        t4 = toc ;
	fprintf(' %0.4f ', t4) 
	v_mat_n_5 = [v_mat_n_5 t4] ;
    end

     fprintf(' \n ' ) 
    save sky_perf_compare_mp_aug23_2.mat v_*

    end
end



clear
index = UFget ;
f = find (index.isReal) ;
[ignore i] = sort (index.nnz (f)) ;
f = f (i) ;
nmat = length (f) ;

% just do the first 100
nmat = 100 ;

figure (1) ;
clf
opts.demo = 1 ;

% S = cell (length (index.nnz), 1) ;

for k = 9 % 1:nmat
    id = f (k) ;
    Prob = UFget (id, index) ;
    opts.name = Prob.name ;
    A = Prob.A ;

    if (0)
        % hack
        [m n] = size (A) ;
        if (m < n)
            A = A' ;
        end
        [m n] = size (A) ;
        Ahack = qr (A) ;
        Ahack = Ahack (1:n, 1:n) ;
    else
        Ahack = A ;
    end

    % clear Prob ;
    try
        [s, stats, ss] = pirosky (Ahack,opts) ;
        % S {id} = s ;
        % save S S
    catch
        disp (lasterr)
    end

    % compare with SVD
    s2 = svd (full (Ahack)) ;
    err = norm (s-s2) / max (s2) ;
    if (err > 1e-12)
        [s s2 s-s2]
        error ('!')
    end

end


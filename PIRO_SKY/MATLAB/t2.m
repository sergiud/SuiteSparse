clear
index = UFget ;
f = find (index.isReal) ;
[ignore i] = sort (index.nnz (f)) ;
f = f (i) ;
nmat = length (f) ;

for k = 1:nmat
    id = f (k) ;
    Prob = UFget (id, index) ;
    A = Prob.A ;

    s1 = pirosky (A) ;
    s2 = svd (full (A)) ;

    err = norm (s1-s2) / s2 (1) ;
    fprintf ('%4d : %4d : %10.3e', k, id, err) ;

    if (err > 1e-8)
        fprintf ('   *********') ;
    end
    fprintf ('\n') ;;

end


load R
full (R)
B = blksky (R) ;

s1 = svd (full (B)) ;
s2 = svd (full (R)) ;

[s1 s2 s1-s2]

err = norm (s1-s2)

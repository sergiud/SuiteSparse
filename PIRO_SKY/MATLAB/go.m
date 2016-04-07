blksky_make
% A = sprandn (10,10,0.1) ;
% R = qr (A) ;
% C = qr_unsqueeze (R) ;
load C
full (C)
B = blksky (C)

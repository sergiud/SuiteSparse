function B = blksky(R, opt)
%
% Reduces sparse upper triangular R into bidiagonal matrix B, using 
% rightmost block, leftmost corner method.
% 
% Usage :
% B = blksky(R) 
% B = blksky(R, opt) 
% 
% TODO I don't like 1,2,3.  This is a hack.  Make it a string?
% opt = 1 => Fixed block size, vanilla symbolic factorization and a static 
% skyline. This is the default option if opt is ignored.
% opt = 2 => Dynamic block size, vanilla symbolic factorization and a static 
% skyline.
% opt = 3 => Fixed block size, no symbolic factorization and a dynamically
% allocated skyline. FAILS. NEED TO BE FIXED. TODO
% 
% TODO : Rectangular matrices, fix the bug in the dynamic data
% structure case for R (opt == 3).
%
% TODO return U and V, optionally
% [U,B,V] = blksky(R) 
% [U,B,V] = blksky(R,opt) 
%
% TODO use an opts struct?  or a string?

error('blksky mexFunction not found')

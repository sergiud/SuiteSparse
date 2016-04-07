function [trow, no_givens, no_swaps, flop_count]  = blksky_symbolic(R)
% Usage :
% [toprow] = bllsky_symbolic_mex(R) ;
% [toprow, no_givens] = bllsky_symbolic_mex(R) ;
% [toprow, no_givens, no_swaps] = bllsky_symbolic_mex(R) ;
% [toprow, no_givens, no_swaps, flop_count] = bllsky_symbolic_mex(R) ;
%
% trow is the minimum row index for every column, while reducing upper
% triangular R, using the blocked rightmost corner method. The maximum required
% space for R csize is
%   csize = 0 ;
%   for k = 1 : #cols of R
%       csize = csize + k - trow(k) + 1 ;
%   end
%
% TBD : Handle rectangular matrices. 

error('blksky_symbolic mexFunction not found')



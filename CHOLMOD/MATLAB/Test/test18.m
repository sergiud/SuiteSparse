function test18
%TEST18 test cholmod2 on a few large matrices
% Example:
%   test18
% See also cholmod_test

% Copyright 2006-2022, Timothy A. Davis, All Rights Reserved.
% SPDX-License-Identifier: GPL-2.0+

fprintf ('=================================================================\n');
fprintf ('test18: test cholmod2 on a few large matrices\n') ;

rand ('state',1) ;
randn ('state',1) ;

Prob = ssget (887)							    %#ok
A = Prob.A ;
n = size (A,1) ;
b = rand (n,1) ;
x = cholmod2 (A,b) ;
norm (A*x-b,1)

clear all

Prob = ssget (936)							    %#ok
A = Prob.A ;
n = size (A,1) ;
b = rand (n,1) ;
x = cholmod2 (A,b) ;
norm (A*x-b,1)

clear all

Prob = ssget (887)							    %#ok
A = Prob.A ;
n = size (A,1) ;
b = rand (n,1) ;
x = cholmod2 (A,b) ;
norm (A*x-b,1)

fprintf ('test18 passed\n') ;

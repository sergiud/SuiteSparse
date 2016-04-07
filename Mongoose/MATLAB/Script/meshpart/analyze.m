function [nonzeros,ops,height,frontsize] = analyze(A,figno)
% ANALYZE : Determine complexity of Cholesky factorization.
%
% [nonzeros,ops,height,frontsize] = analyze(A) 
% predicts the complexity of performing Cholesky factorization on a
% symmetric, positive definite matrix with the nonzero structure of A,
% without actually performing the factorization.
%
% Output:
% nonzeros    Number of nonzeros in the Cholesky factor.
% ops         Number of floating point arithmetic operations to factor A.
% height      Parallel depth of factorization, or height of elimination tree.
% frontsize   Largest frontal matrix (or largest clique size) in factorization.
%
% If no output arguments are given, analyze prints its results.
% If a second argument "figno" is given, analyze draws the elimination tree.
%
% See also SYMBFACT (which gives more detailed information), CHOL.
%
% John Gilbert, 6 Aug 1994.
% Copyright (c) 1990-1996 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.
%
% Modified by Tim Davis, for Matlab 5.1.  July 1998
% Modified by John Gilbert for Matlab 6, Feb 2002

[count,height] = symbfact(A);
nonzeros = sum(count);
ops = sum(count .^ 2);
frontsize = max(count);

if nargout == 0
    fprintf(1,'Nonzeros:   %d\n', nonzeros);
    fprintf(1,'Flops:      %d\n', ops);
    fprintf(1,'Height:     %d\n', height);
    fprintf(1,'Front size: %d\n', frontsize);
end;
if nargin >= 2 
    figure(figno);
    clf reset;
    colordef(gcf,'black')
    if size(A,1) <= 250
        etreeplotg(A,1);
    else
        etreeplotg(A);
    end;
end;


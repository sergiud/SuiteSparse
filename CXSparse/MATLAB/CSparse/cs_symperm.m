function C = cs_symperm (A,p)                                               %#ok
%CS_SYMPERM symmetric permutation of a symmetric matrix.
%   C = cs_symperm(A,p) computes C = A(p,p), but accesses only the
%   upper triangular part of A, and returns C upper triangular (A and C are
%   symmetric with just their upper triangular parts stored).  A must be square.
%
%   Example:
%       Prob = ssget ('HB/bcsstk01') ; A = Prob.A ;
%       p = cs_amd (A) ;
%       C = cs_symperm (A, p) ;
%       cspy (A (p,p)) ;
%       cspy (C) ;
%       C - triu (A (p,p))
%
%   See also CS_PERMUTE, SUBSREF, TRIU.

% CXSparse, Copyright (c) 2006-2022, Timothy A. Davis. All Rights Reserved.
% SPDX-License-Identifier: LGPL-2.1+

error ('cs_symperm mexFunction not found') ;

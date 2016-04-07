function p = ndperm(method,minsep,A,arg2,arg3,arg4)
% NDPERM : Nested dissection ordering for a symmetric matrix.
%
% p = ndperm(method,minsep,A)
% returns a nested dissection permutation of symmetric positive definite A,
% such that A(p,p) tends to have a sparser Cholesky factor than A.
%
% Input:
%   'method' is the name of the 2-way edge separator function to call.
%   minsep   is the size of piece not to partition further.  Default is 8.
%   A        is the symmetric matrix.
%   arg2, arg3, arg4 are optional additional arguments to the 2-way function.
%            arg2 is special : If it has the same number of rows as A,
%            then the recursive calls use the appropriate subset of arg2
%            as well as of A.  This is useful if the partitioner uses xy coords.
%
% Output:
%     p      is a permutation vector.
%
% For example,  ndperm('geopart',10,A,xy,0,ntries)
%   orders A by geometric nested dissection down to parts of size 10, 
%   using a call at each stage that looks like geopart(A,xy,0,ntries).
%
% John Gilbert, 1994.
% Copyright (c) 1990-1996 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.
%
% Modified 6 Feb 02 by JRG to use symamd if present on the pieces.

if length(minsep) == 0
    minsep = 8; 
end;
nargcall = nargin - 2;

[mA,nA] = size(A);
if mA~=nA, error('argument must be square'); end;
if nA == 0, p = []; return; end;
p = zeros(1,nA);

% Find the connected components of A
comps = components(A);

% For each connected component,
numbered = 0;
for i = min(comps) : max(comps)
    comp = find(comps == i);
    nS = length(comp);
    if nS < minsep

        % Bottom of recursion, don't separate further
        % (use minimum degree instead)

        if exist('symamd')
            pS = symamd(A(comp,comp)); 
        else
            pS = symmmd(A(comp,comp));
        end

    else

        % Find a separator for the current component.
        S = A(comp,comp);
        if nargcall >= 2
            if size(arg2,1) == nA
                Sarg2 = arg2(comp,:);
            else
                Sarg2 = arg2;
            end;
        end;
        if nargcall == 4
            a = feval(method, S, Sarg2, arg3, arg4);
        elseif nargcall == 3
            a = feval(method, S, Sarg2, arg3);
        elseif nargcall == 2
            a = feval(method, S, Sarg2);
        else
            a = feval(method, S);
        end;

        % Turn the separator into a vertex separator.
        sep = vtxsep(S,a);

        % Recursive call to number the separated graph.
        unsep = other(sep,S);
        U = S(unsep,unsep);
        if nargcall >= 2
            if size(arg2,1) == nA
                Uarg2 = Sarg2(unsep,:);
            else
                Uarg2 = Sarg2;
            end;
        end;
        if nargcall == 4
            q = ndperm(method,minsep,U,Uarg2,arg3,arg4);
        elseif nargcall == 3
            q = ndperm(method,minsep,U,Uarg2,arg3);
        elseif nargcall == 2
            q = ndperm(method,minsep,U,Uarg2);
        else
            q = ndperm(method,minsep,U);
        end;

        if length(unsep)
            pS = [unsep(q) sep];
        else
            pS = sep;
        end
    end

    % Plug the component permutation into the main permutation
    p(numbered+1:numbered+nS) = comp(pS);
    numbered = numbered + nS;
end

function p = colndperm(method,minsep,A,arg2,arg3,arg4)
% 0-0 unfinished COLNDPERM : Column nested dissection ordering for a nonsymmetric matrix.
%
% p = colndperm(method,minsep,A)
% returns a nested dissection permutation for the columns of A, such
% that A(:,p) tends to have a sparser LU or QR factorization than A.
% This is to be used for nonsymmetric matrices.
%
% Input:
%   'method' is the name of the 2-way edge separator function to call.
%   minsep   is the size of piece not to partition further.  Default is 8.
%   A        is the matrix.
%   arg2, arg3, arg4 are optional additional arguments to the 2-way function.
%            arg2 is special : If it has the same number of rows as A has cols,
%            then the recursive calls use the appropriate subset of arg2
%            as well as of A.  This is useful if the partitioner uses xy coords.
%
% Output:
%     p      is a permutation vector.
%
% For example,  ndperm('geopart',10,A,xy,0,ntries)
%   orders A by geometric nested dissection down to parts of size 10, 
%   using a call at each stage that looks like geopart(A,xy,0,ntries),
%   where H is the augmented matrix of A.
%
% John Gilbert, 9 July 2001
% Copyright (c) 1990-2001 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.
%
% Modified 6 Feb 02 by JRG to use colamd if present on the pieces.

if length(minsep) == 0
    minsep = 8; 
end;
nargcall = nargin - 2;

[nrA,ncA] = size(A);
if ncA == 0, p = []; return; end;
p = zeros(1,ncA);

% Form the augmented matrix H,
% whose graph is the bipartite graph of A
H = [speye(nrA,nrA) A ; A' speye(ncA,ncA)];

% Find the connected components of H
comps = components(H);

% For each connected component,
numbered = 0;
for i = min(comps) : max(comps)
    rowcomp = find(comps(1:nrA) == i);
    colcomp = find(comps(nrA+1 : nrA+ncA) == i);
    nS = length(colcomp);
    if nS < minsep

        % Bottom of recursion, don't separate further
        % (use minimum degree instead)

        if exist('colamd')
            pS = colamd(A(comp,comp)); 
        else
            pS = colmmd(A(comp,comp));
        end


    else

        % Find a separator for the current component.
        comp = [rowcomp colcomp];
        S = H(comp,comp);
        if nargcall >= 2
            if size(arg2,1) == ncA
                Sarg2 = arg2(colcomp,:);
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
        sep = vtxsep(S,a) 0-0;

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

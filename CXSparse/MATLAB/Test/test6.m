function test6
%TEST6 test cs_reach, cs_reachr, cs_lsolve, cs_usolve
%
% Example:
%   test6
% See also: testall

% CXSparse, Copyright (c) 2006-2022, Timothy A. Davis. All Rights Reserved.
% SPDX-License-Identifier: LGPL-2.1+

rand ('state', 0)
maxerr = 0 ;
clf
for trial = 1:201
    n = fix (100 * rand (1)) ;
    d = 0.1 * rand (1) ;
    L = tril (sprandn (n,n,d),-1) + sprand (speye (n)) ;
    b = sprandn (n,1,d) ;

    if (~ispc)
        if (rand ( ) > .5)
            L = L + 1i*sprand(L) ;
        end
        if (rand ( ) > .5)
            b = b + 1i*sprand(b) ;
        end
    end

    for uplo = 0:1

        if (uplo == 1)
            % solve Ux=b instead ;
            L = L' ;
        end

        x = L\b ;
        sr = 1 + cs_reachr (L,b) ;
        sz = 1 + cs_reachr (L,b) ;

        check_if_same (sr,sz) ;

        s2 = 1 + cs_reach (L,b) ;

        try
            if (uplo == 0)
                x3 = cs_lsolve (L,b) ;
            else
                x3 = cs_usolve (L,b) ;
            end
        catch
            if (isreal (L) & isreal (b))                                    %#ok
                lasterr
                error ('!') ;
            end
            % punt: sparse(L)\sparse(b) not handled by cs_lsolve or cs_usolve
            x3 = L\b ;
        end

        % x3 is NOT returned in sorted order, so it is not a valid MATLAB
        % sparse matrix.  It might also have explicit zeros.  Double-transpose
        % it to sort it, and multiply one to remove explicit zeros.
        x3 = 1 * (x3')' ;

        spy ([L b x x3])
        drawnow

        s = sort (sr) ;
        [i j xx] = find (x) ;                                               %#ok
        [i3 j3 xx3] = find (x3) ;                                           %#ok

        if (isempty (i))
            if (~isempty (s))
                i       %#ok
                s       %#ok
                error ('!') ;
            end
        elseif (any (s ~= i))
            i       %#ok
            s       %#ok
            error ('!') ;
        end

        if (isempty (i3))
            if (~isempty (s))
                i3      %#ok
                s       %#ok
                error ('!') ;
            end
        elseif (any (s ~= sort (i3)))
            s       %#ok
            i3      %#ok
            error ('!') ;
        end

        if (any (s2 ~= sr))
            s2      %#ok
            sr      %#ok
            error ('!') ;
        end

        err = norm (x-x3,1) ;
        if (err > 1e-12)
            x       %#ok
            x3      %#ok
            uplo    %#ok
            err     %#ok
            error ('!') 
        end

        maxerr = max (maxerr, err) ;

    end

    drawnow
end
fprintf ('maxerr = %g\n', maxerr) ;

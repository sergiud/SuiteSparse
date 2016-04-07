function bw = find_bw(R, space_fact)
fprintf ('function find_bw:\n') ;
R
if (max(size(R) < 16))
    full (R)
end
space_fact

[m, n] = size(R) ;
[B, d] = spdiags(R) ;
maxd = max(d) ;
cumfill = 0 ;
for i1 = 0 : maxd
    diagid = find ( d == i1) ;
    avail = n - i1 ;
    if (isempty(diagid))
        fill = avail ;
    else
        nz = nnz(B(:, diagid)) ;
        fill = avail - nz ;
    end
    cumfill = cumfill + fill ;
    if (cumfill > (space_fact * nnz(R)))
        i1 ;
        break
    end
end
if (cumfill > (space_fact * nnz(R)))
    bw = i1 - 1 ;
else
    bw = maxd ;
end

bw

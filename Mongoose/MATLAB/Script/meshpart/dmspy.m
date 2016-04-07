function dmspy(A)
% DMSPY		Show Dulmage-Mendelsohn decomposition.
% dmspy(A):  
% A is permuted to block triangular form, and plotted with
% the boundaries between diagonal blocks shown.  Adjacent
% blocks of size 1 are merged, so the resulting diagonal
% blocks are either irreducible (strong Hall) or triangular.
%
% See also SPYPART, DMPERM.

% John Gilbert, 1990-2002.

B = spones(A);
[p,q,r,s] = dmperm(B);

% Collapse 1-by-1 blocks.
dr1 = diff([r inf]);
dr2 = diff([-inf r]);
fr = find(dr1>1 | dr2>1);
r = r(fr);

ds1 = diff([s inf]);
ds2 = diff([-inf s]);
fs = find(ds1>1 | ds2>1);
s = s(fs);

spypart(B(p,q),r,s);
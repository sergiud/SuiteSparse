function [Q,R] = piro_band_extract_qr (QR)
% PIRO_BAND_EXTRACT_QR extracts the Q and R factorization from a QR struct
% returned by QR = piro_band_qr (A).  This function is for testing only.
%
% Example
%   A = test_band_matrix (10, 5, 2, 3)
%   QR = piro_band_qr (A)
%   [Q,R] = piro_band_extract_qr (QR)
%   norm (A-Q*R,1)
%
% See also piro_band_qr.

% Copyright 2012, Sivasankaran Rajamanickam, Timothy A. Davis
% http://www.cise.ufl.edu/research/sparse

m = QR.m ;
n = QR.n ;

% extract R
Rband = QR.R ;
t = size (Rband,1) ;
Ri = zeros (t*n, 1) ;
Rj = zeros (t*n, 1) ;
Rx = zeros (t*n, 1) ;
rnz = 0 ;
for j = 1:n
    % Rband (t,j) becomes R (j,j)
    % Rband (t-1,j) becomes R (j-1,j) ...
    % ...
    % Rband (1,j) becomes R (j-t+1,j) ...
    for k = 1:t
        i = j - t + k ;
        if (i >= 1 && i <= m)
            % R (i,j) = Rband (k,j) ;
            rnz = rnz + 1 ;
            Ri (rnz) = i ;
            Rj (rnz) = j ;
            Rx (rnz) = Rband (k,j) ;
        end
    end
end
R = sparse (Ri (1:rnz), Rj (1:rnz), Rx (1:rnz), m, n) ;

% extract Q as a product of Householder vectors
I = speye (m) ;
Q = I ;
V = QR.V ;
Beta = QR.beta ;
nhouse = size (V,2) ;
t = size (V,1);
for k = 1:nhouse
    beta = Beta (k) ;
    if (beta ~= 0)
        Vi = (k:(t+k-1))' ;
        Vx = V (:,k) ;
        ok = (Vi <= m) ;
        v = sparse (Vi (ok), 1, Vx (ok), m, 1) ;
        H = I - beta * v * v' ;
        Q = H'*Q ;
    end
end
Q = Q' ;

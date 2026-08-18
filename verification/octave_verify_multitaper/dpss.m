function [S,lam] = dpss(N,NW,K)
%DPSS_REF  Independent dense reference for DPSS tapers/eigenvalues, NOT
% using tridieig/tridisolve/fftfilt (i.e. not sharing machinery with the
% library under test). Forms the tridiagonal commuting matrix T densely,
% diagonalizes with Octave's built-in eig(), then gets each concentration
% eigenvalue via the direct quadratic form e'*B*e with B the Toeplitz
% sinc (prolate) matrix -- the textbook definition, not the fftfilt trick
% used elsewhere in this codebase.
W = NW/N;
n = (0:N-1)';
d = ((N-1-2*n).^2)*0.25*cos(2*pi*W);
i1 = (1:N-1)';
ee = i1.*(N-i1)/2;
T = diag(d) + diag(ee,1) + diag(ee,-1);
[V,D] = eig(T);
[~,order] = sort(diag(D),'descend');
V = V(:,order);
S = zeros(N,K);
lam = zeros(K,1);
lags = (0:N-1)';
sincvec = 2*W*sinc(2*W*lags);
sincvec(1) = 2*W;
B = toeplitz(sincvec);
for k=1:K
    e = V(:,k);
    e = e/norm(e);
    S(:,k) = e;
    lam(k) = e'*B*e;
end
end

function x = tridisolve(e,d,b,N)
%TRIDISOLVE  Real implementation (the shipped .m is a MEX-stub with no
% compiled binary in this codebase). Thomas algorithm per Golub & Van Loan,
% "Matrix Computations" 2nd ed., p.156, per the original stub's own doc
% comment. Independent hand-transcription from that reference, not derived
% from the Python port, for a genuine cross-check.
x = b(:);
d = d(:);
n = length(d);
for k = 2:n
    mu = e(k-1)/d(k-1);
    d(k) = d(k) - mu*e(k-1);
    x(k) = x(k) - mu*x(k-1);
end
x(n) = x(n)/d(n);
for k = n-1:-1:1
    x(k) = (x(k) - e(k)*x(k+1))/d(k);
end
end

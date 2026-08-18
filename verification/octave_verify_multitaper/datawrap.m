function y = datawrap(x,M)
%DATAWRAP  Not a MATLAB core builtin here (Signal Processing Toolbox);
% independent shim: alias N>M samples into M by summing wrapped copies.
N = length(x);
pad = mod(-N, M);
if pad > 0
    x = [x(:); zeros(pad,1)];
end
x = reshape(x, M, []);
y = sum(x, 2);
end

function [lambda_all] = firstNlambdaDPSS_fixed(N,W,firstN)
% Same as the original firstNlambdaDPSS.m, with the two bug fixes applied
% in the Python translation: (1) v = v(end:-1:1) before use, (2) parity
% check uses the same index as the seed vector.

    K = ceil(2*N*W-1/2);
    d = ((N-1-2*(0:N-1)').^2)*.25*cos(2*pi*W);
    ee = (1:N-1)'.*(N-1:-1:1)'/2;

    t = (0:N-1)'/(N-1)*pi;
    s = [4*W*sinc(2*W*((N-1):(-1):1)');2*W];
    k = [1,firstN];

    S1 = zeros(N,k(2)-k(1)+1);
    lambda_all = zeros(k(2)-k(1)+1,1);
    idx = 0;

    v = tridieig(d,[0;ee],N-k(2)+1,N-k(1)+1);
    v = v(end:-1:1);  % FIX 1

    for j = 1:(k(2)-k(1)+1)
       seed_index = j+k(1)-1;
       e = sin(seed_index*t);
       e = tridisolve(ee,d-v(j),e,N);
       e = tridisolve(ee,d-v(j),e/norm(e),N);
       e = tridisolve(ee,d-v(j),e/norm(e),N);

       parity_index = seed_index;  % FIX 2 (was k(2)+1-j)
       if(mod(parity_index,2) == 0)
           if(e(2) > 0)
               e = e-e(end:-1:1);
           else
               e = e(end:-1:1)-e;
           end
       else
           if(sum(e) > 0)
               e = e+e(end:-1:1);
           else
               e = -e-e(end:-1:1);
           end
       end

       S1(:,idx+j) = e/norm(e);
       lambda_all(idx+j) = fftfilt(S1(end:-1:1,idx+j),S1(:,idx+j))'*s;

       if(k(1) > 1)
        k = [max(1,k(1)-ceil(log(N))),k(1)-1];
        idx = idx + j;
       end
    end
end

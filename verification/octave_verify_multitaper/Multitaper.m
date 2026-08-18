function MTSE = Multitaper(N,W,cutoff)
%% Multitaper.m 
% Function implementing the exact multitaper spectral estimator
% 
% Inputs:
% -N: Length of DPSSs
% -W: Half Bandwidth for DPSSs
% -cutoff: Projection onto the first K = cutoff DPSSs 
%          OR Projection onto the first K = #{ell : lambda_{N,W}^{(ell)} >= cutoff} DPSSs.
% Outputs:
% -MTSE an object containing an exact Multitaper Spectral Estimate function
%
% Given an Nx1 vector x, performing
% y = MTSE.SpectralEstimate(x,M);
% results in y[m] = (1/K)x^*E_(m/M)S_KS_K^*E_(m/M)^*x
% i.e. y[m] is the multitaper spectral estimate at f = m/M
%
% Notes: 
%
% Most recent change - 01/30/2019
%
% Copyright 2019, Santhosh Karnik

%% Exact Multitaper Function
MTSE.SpectralEstimate = @SpectralEstimate;

%% Store signal length, bandwidth, and number of DPSSs
MTSE.N = N;
MTSE.W = W;
MTSE.cutoff = cutoff;

%% Compute and store DPSSs
if(cutoff >= 1)
    K = floor(cutoff+eps);
    [MTSE.S,MTSE.lambda] = dpss(N,N*W,K);
elseif(cutoff > eps)
    if(cutoff > 0.5)
        Kest = ceil(2*N*W);
    else
        Kest = median([1,ceil(2*N*W+(2/pi^2)*log(N)*log(1./cutoff))+10,N]);
    end
    [S,lambda] = dpss(N,N*W,Kest);
    K = nnz(lambda > cutoff);
    K = max(K,1);
    MTSE.S = S(:,1:K);
    MTSE.lambda = lambda(1:K);
end
MTSE.K = K;

%% Exact Multitaper Spectral Estimation
function y = SpectralEstimate(x,M)
    if(nargin < 2)
        M = MTSE.N;
    end
    
    % Compute tapered signals
    if(MTSE.N > M)
       Sx = bsxfun(@times,MTSE.S,x);
       % Datawrap Sx
       A = zeros(M,MTSE.K);
       for k = 1:MTSE.K,
           A(:,k) = datawrap(Sx(:,k),M);
       end
    else
        A = bsxfun(@times,MTSE.S,x);
    end
    
    % Compute FFTs of tapered signals and sum squared magnitudes
    A = fft(A,M,1);
    A = real(A).^2+imag(A).^2;
    y = sum(A,2)/MTSE.K;
end

end
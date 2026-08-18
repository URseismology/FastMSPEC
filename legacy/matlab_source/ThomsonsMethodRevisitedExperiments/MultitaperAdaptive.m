function MTA = MultitaperAdaptive(N,W,cutoff)
%% MultitaperAdaptive.m 
% Function implementing the exact multitaper spectral estimator with an adaptive weighting scheme
% 
% Inputs:
% -N: Length of DPSSs
% -W: Half Bandwidth for DPSSs
% -cutoff: Projection onto the first K = cutoff DPSSs 
%          OR Projection onto the first K = #{ell : lambda_{N,W}^{(ell)} >= cutoff} DPSSs.
%
% Outputs:
% -MTA an object containing an exact Multitaper Spectral Estimate function with adaptive weights
%
% Given an Nx1 vector x, performing
% y = MTA.SpectralEstimate(x,M);
% results in y[m] as the adaptively weighted multitaper spectral estimate at f = m/M
%
% Notes: 
%
% Most recent change - 01/30/2019
%
% Copyright 2019, Santhosh Karnik

%% Multitaper Adaptive Function
MTA.SpectralEstimate = @SpectralEstimate;

%% Store signal length, bandwidth, and number of DPSSs
MTA.N = N;
MTA.W = W;
MTA.cutoff = cutoff;

%% Compute and store DPSSs
if(cutoff >= 1)
    K = floor(cutoff+eps);
    [MTA.S,MTA.lambda] = dpss(N,N*W,K);
elseif(cutoff > eps)
    if(cutoff > 0.5)
        Kest = ceil(2*N*W);
    else
        Kest = median(1,ceil(2*N*W+(2/pi^2)*log(N)*log(1./cutoff))+10,N);
    end
    [S,lambda] = dpss(N,N*W,Kest);
    K = nnz(lambda > cutoff);
    MTA.S = S(:,1:K);
    MTA.lambda = lambda(1:K);
end
MTA.K = K;

%% Multitaper Spectral Estimation with Adaptive Weights
function [Sf,Skf,wkf,num_iter] = SpectralEstimate(x,M)
    if(nargin < 2)
        M = MTA.N;
    end
    
    % Compute single taper estimates for each DPSS taper
    if(MTA.N > M)
       Sx = bsxfun(@times,MTA.S,x);
       % Datawrap Sx
       Skf = zeros(M,MTA.K);
       for k = 1:MTA.K,
           Skf(:,k) = datawrap(Sx(:,k),M);
       end
    else
        Skf = bsxfun(@times,MTA.S,x);
    end    
    Skf = fft(Skf,M,1);
    Skf = real(Skf).^2+imag(Skf).^2;
    
    % Iterative scheme to compute frequency dependent weights and weighted spectral estimate
    sigma2 = x'*x/MTA.N;
    oneminuslambdasigma2 = repmat((1-MTA.lambda')*sigma2,M,1);
    
    wkf = repmat(MTA.lambda',M,1);
    Sf = (Skf(:,1)+Skf(:,2))/2;
    Sf_old = zeros(M,1);
    num_iter = 0;
    while(sum(abs(Sf-Sf_old)) > 5e-4*sigma2/M)
        num_iter = num_iter + 1;
        Sf_old = Sf;
        wkf = ((Sf.^2)*MTA.lambda')./(Sf*MTA.lambda'+oneminuslambdasigma2).^2;
        Sf = sum(wkf.*Skf,2)./sum(wkf,2);
    end
    
end

end
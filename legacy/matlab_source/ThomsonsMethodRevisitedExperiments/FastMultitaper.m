function FMTSE = FastMultitaper(N,W,cutoff,epsilon)
%% FastMultitaper.m 
% Function implementing the fast multitaper spectral estimator
% 
% Inputs:
% -N: Length of DPSSs
% -W: Half Bandwidth for DPSSs
% -cutoff: Projection onto the first K = cutoff DPSSs 
%          OR Projection onto the first K = #{ell : lambda_{N,W}^{(ell)} >= cutoff} DPSSs.
% -epsilon: Relative error in approximation
%
% Outputs:
% -FMTSE an object containing a Fast Multitaper Spectral Estimate function
%
% Given an Nx1 vector x, performing
% y = FMTSE.SpectralEstimate(x,M);
% results in y[m] \approx (1/K)x^*E_(m/M)S_KS_K^*E_(m/M)^*x
% i.e. y[m] is approximately the multitaper spectral estimate at f = m/M
%
% Notes: 
%
% Most recent change - 01/30/2019
%
% Copyright 2019, Santhosh Karnik

%% Fast Multitaper Function
FMTSE.SpectralEstimate = @SpectralEstimate;

%% Store signal length, bandwidth, tolerance, and number of DPSSs
FMTSE.N = N; % Signal Length
FMTSE.W = W; % Half-Bandwidth
FMTSE.cutoff = cutoff; % Number of DPSSs
FMTSE.epsilon = epsilon; % Tolerance

%% Compute DPSSs and store transition region DPSSs
[S,lambda,lambda_all,lowerindex,upperindex] = transitionDPSS(N,W,epsilon,cutoff);
if(cutoff >= 1)
    K = cutoff;
    eig_weights = [ones(K-lowerindex+1,1); zeros(upperindex-K,1)] - lambda;
elseif(cutoff > eps)
    eig_weights = (lambda > cutoff) - lambda; %positive and negative weighing of weights to eigs
    K = nnz(lambda > cutoff) + lowerindex-1;
end

FMTSE.S = bsxfun(@times,S,sqrt(abs(eig_weights'))); %the abolsute of the eigen weigh takes care of the negative values
FMTSE.index_plus = (eig_weights > 0);
FMTSE.r = length(eig_weights);
FMTSE.K = K;
FMTSE.lambda = lambda;
FMTSE.lambda_all=lambda_all;

%% Precompute entries of the prolate matrix B_{N,W}.
FMTSE.vecHalfSinc = 2*W*sinc(2*W*(1:(N-1))');


%% Fast Multitaper Spectral Estimation
function z = SpectralEstimate(x,M)
    if(FMTSE.N ~= size(x,1))
        disp('ERROR: x has incorrect length')
    elseif(size(x,2) > 1)
        disp('x must be a single column vector')
    else        
        % Compute eigenvalue-weighted multitaper spectral estimate using FFTs
        if(nargin < 2)
            M = FMTSE.N; 
        end
        L = M*ceil(2*FMTSE.N/M);
        vecSinc = [2*FMTSE.W;FMTSE.vecHalfSinc;zeros(L-2*FMTSE.N+1,1);flipud(FMTSE.vecHalfSinc)];
        Fx = fft(x,L,1);
        z0 = ifft(fft(real(Fx).^2+imag(Fx).^2,L,1).*vecSinc,L,1);
        if(L > M)
            z0 = z0(1:(L/M):end);
        end

        % Compute tapered spectral estimates using transition-region tapers
        if(FMTSE.N > M)
           Sx = bsxfun(@times,FMTSE.S,x);
           % Datawrap Sx
           A = zeros(M,FMTSE.r);
           for k = 1:FMTSE.r
               A(:,k) = datawrap(Sx(:,k),M);
           end
        else
            A = bsxfun(@times,FMTSE.S,x);
        end
        A = fft(A,M,1);
        A = A.*conj(A); %real(A).^2+imag(A).^2;
        z1 = sum(A(:,FMTSE.index_plus),2)-sum(A(:,~FMTSE.index_plus),2);
        
        % Add eigenvalue-weighted estimate to transition-region correction
        z = (z0+z1)/FMTSE.K;
        
        % Set any entries <=0 to machine precision times maximum
        maxz = max(z);
        z = max(z,eps*maxz);
    end
end

end
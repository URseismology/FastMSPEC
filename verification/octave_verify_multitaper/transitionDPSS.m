function [S,lambda,lambda_all,lowerindex,upperindex] = transitionDPSS(N,W,epsilon,cutoff)
%% transitionDPSS.m 
% Function calculating the DPSS eigenvectors and eigenvalues in the
% transition region (epsilon(1),1-epsilon(2)). 
%
% This function is used by FastSlepianTransform.m,
% FastSlepianTransformMem.m, FastSlepianProjection.m,
% FastProlatePseudoinverse.m, and FastProlateTikhonov.m
%
% This function requires tridiageig.mexw64 and tridisolve.mexw64.
%
% This function reuses some code from the DPSSCALC function writen by 
% T. Krauss, C. Moler, E. Breitenberger
%
% Inputs:
% -N: Length of DPSSs
% -W: Half Bandwidth for DPSSs
% -epsilon: Define transition region as (epsilon(1),1-epsilon(2))
% -cutoff: If cutoff is an integer, then the transition region is extended
% either leftwards to include the (K+1)-th DPSS or rightwards to include the
% K-th DPSS. If cutoff is a real number in (0,1), then the transition
% region is extended as (min(epsilon(1),cutoff),max(1-epsilon(2),cutoff))
%
% Outputs:
% -S: an Nxr matrix whose columns contain the transition region DPSSs
% -lambda: an rx1 vector containing the transition region eigenvlues
% -lowerindex: index of the first transition region DPSS
% -upperindex: index of the last transition region DPSS
%
% Note: All DPSS indices are MATLAB indices (i.e. starting from 1 not 0)
%
% Most recent change - 08/07/2017
%
% Copyright 2017, Santhosh Karnik
%
% This file is part of Fast Slepian Transform (FST) Toolbox version 1.0.
%
%    FST is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FST is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FST.  If not, see <http://www.gnu.org/licenses/>.
%
% If you use this code in an academic paper, please cite our paper:
% S. Karnik, Z. Zhu, M. B. Wakin, J. K. Romberg, and M. A. Davenport., 
% The fast Slepian transform, Appl. Comput. Harmon. Anal. (2017), 
% http://dx.doi.org/10.1016/j.acha.2017.07.005

% Process input arguments
if(nargin < 3)
    % Use a default tolerance if none is specified
    epsilon = 4*eps*sqrt(N)*[1,1]; 
end
if(length(epsilon) == 1)
    epsilon0 = epsilon;
    epsilon1 = epsilon;
else
    epsilon0 = epsilon(1);
    epsilon1 = epsilon(2);
end
if(nargin < 4)
    cutoff = ceil(2*N*W-1/2);
end
if(cutoff >= 1)    
    K = median([1,floor(cutoff),N-1]);
elseif(cutoff > eps)
    K = ceil(2*N*W-1/2);
    epsilon0 = min(epsilon0,cutoff);
    epsilon1 = min(epsilon1,1-cutoff);
end

% Generate the diagonals
d = ((N-1-2*(0:N-1)').^2)*.25*cos(2*pi*W);  % diagonal of T
ee = (1:N-1)'.*(N-1:-1:1)'/2;               % super diagonal of T with a leading zero

% Stuff needed in eigenvector/eigenvalue computation
t = (0:N-1)'/(N-1)*pi;
s = [4*W*sinc(2*W*((N-1):(-1):1)');2*W];

% Overestimate Region for which epsilon0 < lambda < 1/2
halfgap0 = log(max(N*sqrt(sin(2*pi*W)),2))*log(1./epsilon0)/pi^2+3;
k = [K+1,min(N,ceil(2*N*W-1/2)+ceil(halfgap0))];
if(k(1) > k(2))
    S0 = zeros(N,0);
    lambda0 = zeros(0,1);
    index0 = 0;
    flag = false;
else
    S0 = zeros(N,k(2)-k(1)+1);
    lambda0 = zeros(k(2)-k(1)+1,1); 
    idx = 0;
    flag = true;
end

% Repeat procedure if overestimate isn't large enough (this shouldn't happen)
while(flag)
    % Get the eigenvalues of T.
    v = tridieig(d,[0;ee],N-k(2)+1,N-k(1)+1);
    v = v(end:-1:1);

    lastwarn(''); msg = '';
    % Compute eigenvectors and eigenvalues of T one by one until lambda < epsilon0
    for j = 1:(k(2)-k(1)+1),
       e = sin((j+k(1)-1)*t);
       e = tridisolve(ee,d-v(j),e,N);
       e = tridisolve(ee,d-v(j),e/norm(e),N);
       e = tridisolve(ee,d-v(j),e/norm(e),N);  
       [msg2,id2] = lastwarn('');
       if ~isempty(msg2)
           % warning(message('transitionDPSS:DPSS', j+k(1)-1));
           lastwarn('');
           msg = msg2;
       end
       if(mod(j+k(1)-1,2) == 0)
           % Polarize and symmetrize anti-symmetric dpss
           if(e(2) > 0)
               e = e-e(end:-1:1);
           else
               e = e(end:-1:1)-e;
           end           
       else
           % Polarize and symmetrize symmetric dpss
           if(sum(e) > 0)
               e = e+e(end:-1:1);
           else
               e = -e-e(end:-1:1);
           end
       end
       S0(:,idx+j) = e/norm(e);       
       lambda0(idx+j) = fftfilt(S0(end:-1:1,idx+j),S0(:,idx+j))'*s;
       if(lambda0(idx+j) <= epsilon0)
           % If eigenvalue is less than epsilon0,
           % stop computing eigenvectors and eigenvalues
           flag = false;
           index0 = idx+j-1;
           break;
       end
    end

    if(k(2) == N)
        index0 = idx+j;
        flag = false;
    end
    if((k(2) < N) && flag)
        % If not enough eigenvalues computed, compute more
        disp('Computing More Eigenpairs: 0')
        k = [k(2)+1,min(N,k(2)+ceil(log(N)))];
        idx = idx + j;
    end
end

% Overestimate Region for which 1/2 <= lambda < 1-epsilon1
halfgap1 = log(max(N*sqrt(sin(2*pi*W)),2))*log(1./epsilon1)/pi^2+3;
k = [max(1,ceil(2*N*W+1/2)-ceil(halfgap1)),K];
if(k(1) > k(2))
    S1 = zeros(N,0);
    lambda1 = zeros(0,1);
    index1 = 0;
    flag = false;
else
    S1 = zeros(N,k(2)-k(1)+1);
    lambda1 = zeros(k(2)-k(1)+1,1); 
    idx = 0;
    flag = true;
end

% Repeat procedure if overestimate isn't large enough (this shouldn't happen)
while(flag)
    % Get the eigenvalues of T.
    v = tridieig(d,[0;ee],N-k(2)+1,N-k(1)+1);

    lastwarn(''); msg = '';    
    % Compute eigenvectors and eigenvalues one by one until lambda < epsilon0
    for j = 1:(k(2)-k(1)+1),
       e = sin((k(2)+1-j)*t);
       e = tridisolve(ee,d-v(j),e,N);
       e = tridisolve(ee,d-v(j),e/norm(e),N);
       e = tridisolve(ee,d-v(j),e/norm(e),N);  
       [msg2,id2] = lastwarn('');
       if ~isempty(msg2)
     % warning(message('transitionDPSS:DPSS', j+k(1)-1));
           lastwarn('');
           msg = msg2;
       end
       if(mod(k(2)+1-j,2) == 0)
           % Polarize and symmetrize anti-symmetric dpss
           if(e(2) > 0)
               e = e-e(end:-1:1);
           else
               e = e(end:-1:1)-e;
           end           
       else
           % Polarize and symmetrize symmetric dpss
           if(sum(e) > 0)
               e = e+e(end:-1:1);
           else
               e = -e-e(end:-1:1);
           end
       end
       S1(:,idx+j) = e/norm(e);  
       lambda1(idx+j) = fftfilt(S1(end:-1:1,idx+j),S1(:,idx+j))'*s;
       if(lambda1(idx+j) >= 1-epsilon1)
           % If eigenvalue is greater than 1-epsilon1,
           % stop computing eigenvectors and eigenvalues
           flag = false;
           index1 = idx+j-1;
           break;
       end
    end
    
    if(k(1) == 1)
        index1 = idx+j;
        flag = false;
    end
    if((k(1) > 1) && flag)
        % If not enough eigenvalues computed, compute more
        disp('Computing More Eigenpairs: 1')
        k = [max(1,k(1)-ceil(log(N))),k(1)-1];
        idx = idx + j;
    end
end

% Merge eigenvalues and eigenvectors
lambda = [lambda1(index1:-1:1);lambda0(1:index0)];
%lambda_all = [lambda1;lambda0];
lambda_all.lambda1 = lambda1;
lambda_all.lambda0 = lambda0;

S = [S1(:,index1:-1:1),S0(:,1:index0)];
upperindex = K+index0;
lowerindex = K-index1+1;
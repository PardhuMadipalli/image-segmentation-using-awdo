function [rec] = purecmaes_wdo(counteval,rec,npop,pressure)
%   coeff -- CMAES optimized coefficients returned to WDO.
%   counteval -- Iteration counter.
%   rec -- Record of prior values used in CMAES.
%   npop -- number of population members from WDO, each member gets their
%          own set of coefficients determined by the CMAES.
%   pressure -- pressure(cost function) computed by WDO for the set of
%               coefficients that CMEAS picked last iteration
%---------------------------------------------------------
%   This script is modified version of the publicly available purecmaes.m
%   Original "purecmaes.m" can be accessed from:
%   URL: http://www.lri.fr/~hansen/purecmaes.m
%
%   This function takes in pressure values from the WDO and computes new
%   set of coefficients for the WDO and returns them. In simpler terms, it
%   optimizes the WDO coefficients of "alpha, g, RT, c"
%
%   This function is not optimized, only proof-of-concept code to
%   illustrate that the WDO coefficients can be optimized by CMAES, 
%   which in turn creates and adaptive WDO algorithm.
%---------------------------------------------------------

if counteval==2    %initialization only happens when the CMAES called for the first time.
    rec.N = 4;  % problem dimension is fixed to 4: alp, g, RT, c 
    rec.xmean = rand(rec.N,1);  % objective variables initial point
    rec.sigma = 0.5;        % coordinate wise standard deviation (step size)
  
    % Strategy parameter setting: Selection  
%     lambda = 4+floor(3*log(N));  % population size, offspring number
    rec.lambda = npop; %this is defined by the wdo population size.
    rec.mu = rec.lambda/2;               % number of parents/points for recombination
    rec.weights = log(rec.mu+1/2)-log(1:rec.mu)';
    if counteval==2
        rec.weights;
    end
        % muXone array for weighted recombination
    rec.mu = floor(rec.mu);        
    rec.weights = rec.weights/sum(rec.weights);     % normalize recombination weights array
    rec.mueff=sum(rec.weights)^2/sum(rec.weights.^2); % variance-effectiveness of sum w_i x_i

    % Strategy parameter setting: Adaptation
    rec.cc = (4 + rec.mueff/rec.N) / (rec.N+4 + 2*rec.mueff/rec.N); % time constant for cumulation for C
    rec.cs = (rec.mueff+2) / (rec.N+rec.mueff+5);  % t-const for cumulation for sigma control
    rec.c1 = 2 / ((rec.N+1.3)^2+rec.mueff);    % learning rate for rank-one update of C
    rec.cmu = min(1-rec.c1, 2 * (rec.mueff-2+1/rec.mueff) / ((rec.N+2)^2+rec.mueff));  % and for rank-mu update
    rec.damps = 1 + 2*max(0, sqrt((rec.mueff-1)/(rec.N+1))-1) + rec.cs; % damping for sigma 
                                                      % usually close to 1
    % Initialize dynamic (internal) strategy parameters and constants
    rec.pc = zeros(rec.N,1); 
    rec.ps = zeros(rec.N,1);   % pc, ps: evolution paths for C and sigma
    rec.B = eye(rec.N,rec.N);                       % B defines the coordinate system
    rec.D = ones(rec.N,1);                      % diagonal D defines the scaling
    rec.C = rec.B * diag(rec.D.^2) * rec.B';            % covariance matrix C
    rec.invsqrtC = rec.B * diag(rec.D.^-1) * rec.B';    % C^-1/2 
    rec.eigeneval = 0;                      % track update of B and D
    rec.chiN=rec.N^0.5*(1-1/(4*rec.N)+1/(21*rec.N^2));  % expectation of  ||N(0,I)|| == norm(randn(N,1))
end

    % get the fitness from WDO pressure:
    rec.arfitness = pressure';
    
    % Sort by fitness and compute weighted mean into xmean
    [rec.arfitness, rec.arindex] = sort(rec.arfitness, 1, 'descend');  % maximization
    rec.xold = rec.xmean;
    rec.xmean = rec.arx(:,rec.arindex(1:rec.mu)) * rec.weights;  % recombination, new mean value
    
    % Cumulation: Update evolution paths
    rec.ps = (1-rec.cs) * rec.ps ... 
          + sqrt(rec.cs*(2-rec.cs)*rec.mueff) * rec.invsqrtC * (rec.xmean-rec.xold) / rec.sigma; 
    rec.hsig = sum(rec.ps.^2)/(1-(1-rec.cs)^(2*counteval/rec.lambda))/rec.N < 2 + 4/(rec.N+1);
    rec.pc = (1-rec.cc) * rec.pc ...
          + rec.hsig * sqrt(rec.cc*(2-rec.cc)*rec.mueff) * (rec.xmean-rec.xold) / rec.sigma; 

    % Adapt covariance matrix C
    rec.artmp = (1/rec.sigma) * (rec.arx(:,rec.arindex(1:rec.mu)) - repmat(rec.xold,1,rec.mu));  % mu difference vectors
    rec.C = (1-rec.c1-rec.cmu) * rec.C ...                   % regard old matrix  
         + rec.c1 * (rec.pc * rec.pc' ...                % plus rank one update
                 + (1-rec.hsig) * rec.cc*(2-rec.cc) * rec.C) ... % minor correction if hsig==0
         + rec.cmu * rec.artmp * diag(rec.weights) * rec.artmp'; % plus rank mu update 

    % Adapt step size sigma
    rec.sigma = rec.sigma * exp((rec.cs/rec.damps)*(norm(rec.ps)/rec.chiN - 1)); 
    
    % Update B and D from C
    if counteval - rec.eigeneval > rec.lambda/(rec.c1+rec.cmu)/rec.N/10  % to achieve O(N^2)
      rec.eigeneval = counteval;
      rec.C = triu(rec.C) + triu(rec.C,1)'; % enforce symmetry
      rec.C=min(100000, rec.C);
      [rec.B,rec.D] = eig(rec.C);           % eigen decomposition, B==normalized eigenvectors
      rec.D = sqrt(diag(rec.D));        % D contains standard deviations now
      rec.invsqrtC = rec.B * diag(rec.D.^-1) * rec.B';
    end
    
    % Generate and evaluate lambda offspring
    for k=1:rec.lambda,
        rec.arx(:,k) = rec.xmean + rec.sigma * rec.B * (rec.D .* randn(rec.N,1)); % m + sig * Normal(0,C) 
    end

end
% end-of-CMAES.
%----------------------------------------------------------------------
%----------------------------------------------------------------------
%----------------------------------------------------------------------

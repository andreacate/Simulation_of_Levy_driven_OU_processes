function [prices, conf_int] = priceAmericanLS (S, T, K, r)
% Computation of American Call price with algorithm Longstaff-Schwartz [2]
%
% INPUT:
% S:            simulated underlying
% T:            time to maturity
% K:            strikes
% r:            rate
%
% OUTPUT:
% prices:       prices with Longstaff-Schwartz algorithm 
% conf_int:     confindence intervals 

  %% Quantities of interest
    % Calculate number of simulations (N), time steps (M), and time step size (dt)
    N = size(S, 1);
    M = size(S, 2)-1;
    dt = T/M;
    
    % Initialize exercise time for all simulations 
    exerciseTime = M*ones(N, 1);
    
    % Initialize price and confidence interval storage
    prices = zeros(size(K));
    conf_int = zeros(2, size(K,2));

    %% LS Algorithm
    % Loop through each strike price
    for jj = 1:length(prices)
        % Recursive computation of option value at each time step
        %   1. Start with payoff at maturity (max(underlying - strike))
        payoffAM = max(0, S(:,end)-K(jj));
        for ii = M:-1:1
            % Identify simulations where option is In-The-Money (positive payoff)
            IntrValue = max(0, S(:,ii)-K(jj));
            indexITM = find(IntrValue > 0);
            
            % Perform linear regression on ITM payoffs at previous time step to predict continuation value
            %   a) Gather ITM underlying prices and payoffs
            %   b) Build regression matrix (ones, underlying price, underlying price squared)
            %   c) Solve for regression weights
            b = payoffAM(indexITM) .* exp(-r*dt*(exerciseTime(indexITM)-ii+1));

            % Canonical basis
            A = [ones(length(indexITM),1), S(indexITM, ii), S(indexITM, ii).^2];

            % Legendre basis (uncomment for use)
            % A = [ones(length(indexITM),1), S(indexITM, ii), 0.5*(3*S(indexITM, ii).^2-1)];

            % Solving the linear system
            weights = A\b;
            
            % Calculate Continuation Value (CV) and Intrinsic Value (IV) for ITM simulations
            CV = A * weights;
            IV = IntrValue(indexITM);
            
            % Identify simulations for early exercise (where intrinsic value is greater than continuation value)
            indexEarlyEx = find(IV>CV);
            
            % Update payoff and exercise time for early exercise simulations
            payoffAM(indexITM(indexEarlyEx)) = IV(indexEarlyEx);
            exerciseTime(indexITM(indexEarlyEx)) = ii-1;
        end
        
        % Discounted price and confidence interval for current strike
        discPrice = payoffAM .* exp(-r*dt*(exerciseTime));
        [prices(jj), ~, conf_int(:, jj)] = normfit(discPrice);
    end
end % function priceAmericanLS

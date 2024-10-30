function [prices, conf_intervals] = priceAmericanGaussian(S0, b, sigma, T, r, moneyness, M, Nsim)
% Computation of American Call price for Gaussian OU process
%
% INPUT:
% S0:        initial underlying condition
% b:         mean reverting parameter
% sigma:     volatility
% T:         time to maturity
% r:         rate
% moneyness: vector of moneyness
% Nsim:      number of simulations
%
% OUTPUT:
% prices:          prices of European call
% std_dev:         standard deviation of prices
% conf_intervals:  confidence intervals of prices
%
% USES:
% function sim_Gaussian_OU_ED(x0, b, sigma, T, M, Nsim)
% Simulation of the process

    [~, ~, ~, ~, ~,logFwd] = sim_Gaussian_OU_ED(0, b,sigma, T, M, Nsim);
    
     % Forward computation
    fwd = S0 * exp(r*T) ;
    
    % Strikes (expressed as a proportion of the forward prices
    K = fwd * exp(-moneyness);
    
    % Calculate forward price
    F = fwd * exp(logFwd);
    
    % Underlying prices at each time step
    S = F .* exp(-r*T/M*[M:-1:0]);
    
    % Price computation using algorithm in [2] Longstaff & Schwartz
    [prices, conf_intervals] = priceAmericanLS(S, T, K, r);


end % function priceAmericanGaussian



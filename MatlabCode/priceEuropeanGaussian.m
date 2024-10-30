function [prices, std_dev, conf_intervals] = priceEuropeanGaussian(S0, b, sigma, T, r, moneyness, Nsim)
% Computation of European Call price for Gaussian OU process
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

    [~, ~, ~, ~, ~,logFwd] = sim_Gaussian_OU_ED(0, b,sigma, T, 1, Nsim);
    
    % Forward price using risk-free interest rate
    fwd = S0 * exp(r*T) ; 

    % Strikes 
    K = fwd * exp(-moneyness);

    % Underlying price at maturity from simulations
    underlying_sim = fwd * exp(logFwd(:, end));

    % Discounted payoff calculation
    payoff = max(repmat(underlying_sim, 1, length(K)) - K, 0);  
    disc_payoff = payoff * exp(-r*T);  
    
    % Price estimation using normal distribution fitting
    [prices, std_dev, conf_intervals] = normfit(disc_payoff);
    
end % function priceEuropeanGaussian
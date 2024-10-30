function [b_opt, alpha_opt, beta_p_opt, beta_n_opt, c_p_opt, c_n_opt, gamma_c_opt] = ...
                 calibration(S0, T, moneyness, r, flag, M_fft, marketPrices)
% Parameters calibration for TS-OU and OU-TS processes
%
% INPUT:
% S0:            initial underlying condition
% T:             time to maturity
% moneyness:     vector of moneyness
% r:             rate
% flag:          1 -> OU-CTS Finite Activity
%                2 -> CTS-OU Finite Variation
%                3 -> OU-CTS Finite Variation 
% M_fft:         parameter of FFT
% marketPrices:  quoted prices
%
% OUTPUT:
% b_opt:       optimal value for mean reverting parameter b
% alpha_opt:   optimal value for stability parameter alpha
% beta_p_opt:  optimal value for parameter beta positive
% beta_n_opt:  optimal value for parameter beta negative
% c_p_opt:     optimal value for parameter c positive
% c_n_opt:     optimal value for parameter c negative
% gamma_c_op:  optimal value for drift parameter 
%
% USES:
% function priceEuropeanLewis_FFT(S0, b, alpha, beta_p, beta_n, c_p, c_n, gamma_c, T, moneyness, r, flag, M_fft)

% NOTE: in the variable "params" (1, 2, 3, 4, 5, 6, 7) = (b, alpha, beta_p,
% beta_n, c_p, c_n, gamma_c)

   
    %% Price calculation with Lewis formula
    % Create a function handle to the pricing function using Lewis formula
    priceFunHandle = @(params) priceEuropeanLewis_FFT(S0, params(1), params(2), params(3), ...
                     params(4), params(5), params(6), params(7), T, moneyness, r, flag, M_fft);
    
    %% Optimization process
    % Initial guess for the parameters (b, alpha, beta_p, beta_n, c_p, c_n, gamma_c)
    initialGuess = [1, 0.5-(flag==1), 1, 1, 1, 1, 0];
    % Lower bounds for the parameters 
    lowerBounds = [0, -1.5*(flag==1), 0, 0, 0, 0, -10];
    % Upper bounds for the parameters
    upperBounds = [10, 1*((flag==3)+(flag==2)), 10, 10, 10, 10, 10];
    
    options = optimoptions('lsqnonlin', 'Display', 'off');
    
    % Perform non-linear least squares optimization using lsqnonlin
    optParams = lsqnonlin(@(params) priceFunHandle(params)-marketPrices, ...
                           initialGuess, lowerBounds, upperBounds, options);
    
    %% Extracting optimized parameters
    
    b_opt = optParams(1);
    alpha_opt = optParams(2);
    beta_p_opt = optParams(3);
    beta_n_opt = optParams(4);
    c_p_opt = optParams(5);
    c_n_opt = optParams(6);
    gamma_c_opt = optParams(7);
    
end % function calibration

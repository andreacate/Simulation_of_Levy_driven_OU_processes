function [price, conf_int] = priceAmerican (x0, alpha, b, beta_p, beta_n, ...
    c_p, c_n, gamma_c, T, M, M_fft, Nsim, K, r, S0, scale, model, method)
% Computation of American Call price using the simulation procedure via
% Exact Decomposition (Sabino [4]) or FGMC (Baviera [1]) for the different
% model
%
% INPUT:
% x0:       initial condition
% alpha:    stability parameter
% b:        mean reverting parameter
% beta_p:   positive beta
% beta_n:   negative beta
% c_p:      positive c
% c_n:      negative c
% gamma_c:  drift
% T:        time to maturity
% M:        number of steps
% M_fft:    parameter of FFT
% Nsim:     number of simulations
% K:        strikes
% r:        rate
% scale:    scaling factor
% model:    1 -> OU-CTS Finite Activity
%           2 -> CTS-OU Finite Variation
%           3 -> OU-CTS Finite Variation
% method:   1 -> Exact Decomposition
%           2 -> FGMC
%
% OUTPUT:
% price:    price of American call
% conf_int: confindence intervals
%
% USES:
% function sim_OU_TS_FinAct_ED(x0, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c, T, M, Nsim)
% function sim_OU_TS_FinAct_FGMC (x0, b, alpha, beta_p, beta_n, c_p, c_n, gamma_c, T, Nsim, M, M_fft, scale)
% function sim_OU_TS_FinVar_ED(x0, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c, T, M, Nsim)
% function sim_OU_TS_FinVar_FGMC (x0, b, alpha, beta_p, beta_n, c_p, c_n, gamma_c, T, Nsim, M, M_fft, scale)
% function sim_TS_OU_FinVar_ED(x0, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c, T, M, Nsim)
% function sim_TS_OU_FinVar_FGMC (x0, b, alpha, beta_p, beta_n, c_p, c_n, gamma_c, T, Nsim, M, M_fft, scale)
% function priceAmericanLS (S, T, K, r)
    
    
    % Quantities of interest
    
    dt = T/M;
    
    % Simulations based on model and method used
    switch model
        case 1
            if (method == 1)
                [~, ~, ~, ~, ~, logFwd] = sim_OU_TS_FinAct_ED(x0, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c, T, M, Nsim);
            else
                [~, ~, ~, ~, ~, logFwd] = sim_OU_TS_FinAct_FGMC(x0, b, alpha, beta_p, beta_n, c_p, c_n, gamma_c, T, Nsim, M, M_fft, scale);
            end
        case 2
            if (method == 1)
                [~, ~, ~, ~, ~, logFwd] = sim_TS_OU_FinVar_ED(x0, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c, T, M, Nsim);
            else
                [~, ~, ~, ~, ~, logFwd]  = sim_TS_OU_FinVar_FGMC(x0, b, alpha, beta_p, beta_n, c_p, c_n, gamma_c, T, Nsim, M, M_fft, scale);
            end
        case 3
            if (method == 1)
                [~, ~, ~, ~, ~, logFwd] = sim_OU_TS_FinVar_ED(x0, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c, T, M, Nsim);
            else
                [~, ~, ~, ~, ~, logFwd]  = sim_OU_TS_FinVar_FGMC(x0, b, alpha, beta_p, beta_n, c_p, c_n, gamma_c, T, Nsim, M, M_fft, scale);
            end
    end

    % Forward computation
    fwd = S0 * exp(r*T) ;
    
    % Calculate forward price
    F = fwd * exp(logFwd);

    % Underlying prices at each time step
    S = F .* exp(-r*dt*[M:-1:0]); 

    % Price computation using algorithm in [2] Longstaff & Schwartz
    [price, conf_int] = priceAmericanLS(S, T, K, r);
    
end % function priceAmerican
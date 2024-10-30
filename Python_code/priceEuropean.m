function [prices, std_dev, conf_intervals] = priceEuropean (x0, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c, T, M, ...
                                             M_fft, Nsim, r, S0, moneyness, scale, model, method)
% Computation of European Call price using the simulation procedure via
% Exact Decomposition (Sabino [4]) or FGMC (Baviera [1]) for the different
% model
%
% INPUT:
% x0:        initial condition
% alpha:     stability parameter    
% b:         mean reverting parameter
% beta_p:    positive beta
% beta_n:    negative beta
% c_p:       positive c
% c_n:       negative c
% gamma_c:   drift
% T:         time to maturity
% M:         number of steps
% M_fft:     parameter of FFT
% Nsim:      number of simulations
% r:         rate
% S0:        initial underlying condition
% moneyness: vector of moneyness
% scale:     scaling factor
% model:     1 -> OU-CTS Finite Activity
%            2 -> CTS-OU Finite Variation
%            3 -> OU-CTS Finite Variation 
% method:    1 -> Exact Decomposition
%            2 -> FGMC
%
% OUTPUT:
% prices:          prices of European call
% std_dev:         standard deviation of prices
% conf_intervals:  confidence intervals of prices
%
% USES:
% function sim_OU_TS_FinAct_ED(x0, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c, T, M, Nsim)
% function sim_OU_TS_FinAct_FGMC (x0, b, alpha, beta_p, beta_n, c_p, c_n, gamma_c, T, Nsim, M, M_fft, scale)
% function sim_OU_TS_FinVar_ED(x0, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c, T, M, Nsim)
% function sim_OU_TS_FinVar_FGMC (x0, b, alpha, beta_p, beta_n, c_p, c_n, gamma_c, T, Nsim, M, M_fft, scale)
% function sim_TS_OU_FinVar_ED(x0, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c, T, M, Nsim)
% function sim_TS_OU_FinVar_FGMC (x0, b, alpha, beta_p, beta_n, c_p, c_n, gamma_c, T, Nsim, M, M_fft, scale)


        % Simulations based on model and method used
        
        switch model
            case 1
                if (method == 1)
                    [~, ~, ~, ~, ~, logFwd] = sim_OU_TS_FinAct_ED(x0, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c, T, M, Nsim);
                else
                    [~, ~, ~, ~, ~, logFwd] = sim_OU_TS_FinAct_FGMC (x0, b, alpha, beta_p, beta_n, c_p, c_n, gamma_c, T, Nsim, M, M_fft, scale);
                end
            case 2
                if (method == 1)
                    [~, ~, ~, ~, ~, logFwd] = sim_TS_OU_FinVar_ED(x0, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c, T, M, Nsim);
                else
                    [~, ~, ~, ~, ~, logFwd]  = sim_TS_OU_FinVar_FGMC (x0, b, alpha, beta_p, beta_n, c_p, c_n, gamma_c, T, Nsim, M, M_fft, scale);
                end
            case 3
                if (method == 1)
                    [~, ~, ~, ~, ~, logFwd] = sim_OU_TS_FinVar_ED(x0, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c, T, M, Nsim);
                else
                    [~, ~, ~, ~, ~, logFwd]  = sim_OU_TS_FinVar_FGMC (x0, b, alpha, beta_p, beta_n, c_p, c_n, gamma_c, T, Nsim, M, M_fft, scale);
                end
        end

       

        % Forward computation
        fwd = S0 * exp(r*T) ;

        % Strikes computation from moneyness
        K = fwd * exp(-moneyness);

        % Underlying
        underlying_sim = fwd * exp(logFwd(:, end));

        % Discounted payoff
        payoff = max(repmat(underlying_sim, 1, length(K)) - K, 0);
        disc_payoff = payoff * exp(-r*T);

        % Prices and confidence interval
        [prices, std_dev, conf_intervals] = normfit(disc_payoff);

end % function priceEuropean
function [X, theorCumulantsT, simCumulantsT, theorCumulants_dt, simCumulants_dt,logFwd] = ...
                                                        sim_OU_TS_FinVar_FGMC (x0, b, alpha, beta_p, ...
                                                        beta_n, c_p, c_n, gamma_c, T, Nsim, M, M_fft, scale)
% Simulation of OU-TS Finite Variation with Fast and General Monte Carlo
% following the algorithm in [1] Baviera & Manzoni

% INPUT
% x0:       initial condition  
% b:        mean reverting parameter
% alpha:    stability parameter  
% beta_p:   positive beta
% beta_n:   negative beta
% c_p:      positive c
% c_n:      negative c
% gamma_c:  drift
% T:        time to maturity
% Nsim:     number of simulations
% M:        number of steps
% M_fft:    parameter for FFT
% scale:    scale to be applied to the CDF
% 
% OUTPUT
% X:                   matrix with simulations
% theorCumulantsT:     theoretical cumulants on the total maturity
% simCumulantsT:       cumulants of the simulations on the total maturity
% theorCumulants_dt:   theoretical cumulants on the first time step
% simCumulants_dt:     cumulants of the simulations on the first time step
% logFwd:              log of fwd prices
% 
% USES
% function optimal_du_formula(alpha, a, b, c_p, c_n, dt, M, scale, flag)
% function integralViaFFT_CDF(M_fft, du, dt, alpha, a, b, beta_p, beta_n, c_p, c_n, gamma_c, scale, flag)
% function selectCDF(x, discreteCDF, toll)
% function samplingFromDiscreteCDF(discrete_x, discrete_CDF, N_sample)
% function bctsCumulants(X0, alpha, beta_p, beta_n, c_p, c_n, gamma_c, dt, b, flag)
% function compCumulants(vec)
% function psi_Z_OU_TS(u, dt, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c)

    %% Quantities of interest

    flagOU_TS_FV = 3;  % Flag for OU-TS Finite Variation

    % Calculate the shifting parameter 'a' based on beta_p and beta_n
    a = -1/2 * max(beta_p, beta_n) * ((beta_p > beta_n) - (beta_n > beta_p));
    % Computing R necessary in the formula for the CDF
    R = (a < 0);  
    % Time step size
    dt = T / M;  
      
    
    %% Search of optimal du with formula

    % Calculate the optimal du for FFT using a formula
    du_optimal = optimal_du_formula(alpha, a, b, c_p, c_n, dt, M_fft, scale, flagOU_TS_FV);

    
    %% CDF computation with FFT
 
    % Compute the integral and grid for FFT
    [integral_FFT_u, x_grid] = integralViaFFT_CDF(M_fft, du_optimal, dt, alpha, ...
                               a, b, beta_p, beta_n, c_p, c_n, gamma_c, scale, flagOU_TS_FV);
 
    % Compute the discrete CDF from the FFT results
    discreteCDF_FFT = R - exp(a .* x_grid) / pi .* 1/2 .* integral_FFT_u;
 

    %% Plotting the CDF before selection
    % Uncomment for CDF plot before selection
    %   figure()
    % 
    %   plot(x_grid, discreteCDF_FFT)
    %   legend('discrete CDF')
    %   title('discrete CDF before selection')
    %   xlim([-1 1])

    %% Selecting the discrete CDF
    
    toll = 1e-8;  % Tolerance for selecting CDF points
    [x_CDF, values_CDF] = selectCDF(x_grid, discreteCDF_FFT, toll);

    %% Plotting the CDF after selection
    % Uncomment for CDF plot after selection
    %     figure()
    %     plot(x_CDF, values_CDF)
    %     legend('discrete CDF')
    %     title('discrete CDF after selection')

    %% Sampling from the CDF

    X = zeros(Nsim, M+1);  % Initialize simulation matrix
    X(:, 1) = x0 * ones(Nsim, 1);  % Set initial condition
    drift = zeros(1, M);  % Initialize drift vector

    for ii = 1:M
        % Sample values from the selected discrete CDF
        sampleValues = samplingFromDiscreteCDF(x_CDF, values_CDF, Nsim / 2);
        
        % Update process values with the sampled increments and exponential decay
        X(:, ii+1) = X(:, ii) * exp(-b * dt) + sampleValues / scale;
        
        % Compute the drift term
        drift(ii) = -psi_Z_OU_TS(-1i, dt * ii, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c);
    end

    %% Computing the log forward in the risk neutral measure

    % Compute the log forward prices by adding the drift term
    logFwd = X + repmat([0, drift], Nsim, 1);

    %% Computing cumulants

    % Checking the last time step
    simCumulantsT = compCumulants(X(:, M+1)) * 1000;  
    theorCumulantsT = bctsCumulants(x0, alpha, beta_p, beta_n, c_p, c_n, gamma_c, T, b, flagOU_TS_FV) * 1000;  

    % Checking the first time step
    simCumulants_dt = compCumulants(X(:, 2)) * 1000;  
    theorCumulants_dt = bctsCumulants(x0, alpha, beta_p, beta_n, c_p, c_n, gamma_c, dt, b, flagOU_TS_FV) * 1000;  % Theoretical cumulants for the first time step

end % function sim_OU_TS_FinVar_FGMC
function [X, theorCumulantsT, simCumulantsT, theorCumulants_dt, simCumulants_dt, logFwd] ...
                                                     = sim_OU_TS_FinAct_FGMC (x0, b, alpha, ...
                                                       beta_p, beta_n, c_p, c_n, gamma_c, T, Nsim, M, M_fft, scale)
% Simulation of OU-TS Finite Activity with Fast and General Monte Carlo 
% following the algorithm in [1] Baviera & Manzoni
%
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
% function integralViaFFT_CDF(M_fft, du, dt, alpha, a, b, beta_p, beta_n, c_p, c_n, gamma_c, scale, flag)
% function selectCDF(x, discreteCDF, toll)
% function samplingFromDiscreteCDF(discrete_x, discrete_CDF, N_sample)
% function bctsCumulants(X0, alpha, beta_p, beta_n, c_p, c_n, gamma_c, dt, b, flag)
% function compCumulants(vec)
% function psi_Z_OU_TS(u, dt, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c)


    %% Quantities of interest

    flagOU_TS_FA = 1;  % Flag for OU-TS Finite Activity model

    % Calculate the shifting parameter 'a' based on beta_p and beta_n
    a = -1/4 * max(beta_p, beta_n) * ((beta_p > beta_n) - (beta_n > beta_p));  
    % Computing R necessary in the formula for the CDF
    R = (a < 0);  
    % Time step size
    dt = T / M; 

    % Lambda computation
    Gamma = gamma(-alpha);  
    lambda_p = c_p * Gamma * beta_p^alpha;  
    lambda_n = c_n * Gamma * beta_n^alpha; 
    lambda_tot = lambda_p + lambda_n; 
    
    %% Sampling of the Bernoulli

    U = rand(Nsim, M);  % Uniform random variables for Bernoulli sampling
    Bern = U < (1 - exp(-lambda_tot * dt));  % Bernoulli process for jump occurrences

    %% CDF computation with FFT 

    % Parameters of FFT
    % Optimal step size for the FFT integral
    du_optimal = 0.4;  
    
    % Computing the integral via FFT in the formula for the CDF
    [integral_FFT_u, x_grid] = integralViaFFT_CDF(M_fft, du_optimal, dt, alpha, a, ...
                               b, beta_p, beta_n, c_p, c_n, gamma_c, scale, flagOU_TS_FA);

    % Computation of CDF
    discreteCDF_FFT = R - exp(a .* x_grid) / pi .* 1/2 .* integral_FFT_u;

    %% Plotting the CDF before selection
    % Uncomment for CDF plot before selection
    % figure()
    % plot(x_grid, discreteCDF_FFT)
    % legend('discrete CDF')
    % title('discrete CDF before selection')
    % xlim([-1 1])

    %% Selecting the discrete CDF
    
    toll = 0;  % Tolerance for selecting points from the CDF
    [x_CDF, values_CDF] = selectCDF(x_grid, discreteCDF_FFT, toll);

    %% Plotting the CDF after selection
     % Uncomment for CDF plot after selection
    % figure()
    % plot(x_CDF, values_CDF)
    % legend('discrete CDF')
    % title('discrete CDF after selection')
    % xlim([-1 1])

    %% Sampling from the CDF

    X = zeros(Nsim, M+1);  % Initialize the matrix to store the simulations
    X(:,1) = x0 * ones(Nsim, 1);  % Set the initial condition for all simulations
    drift = zeros(1, M);  % Initialize the drift vector

    for ii = 1:M
        [sampleValues] = samplingFromDiscreteCDF(x_CDF, values_CDF, Nsim/2);  % Sample from the discrete CDF
        % Update process just if the value of the bernoulli is 1
        X(:, ii+1) = X(:, ii) * exp(-b * dt) + sampleValues .* Bern(:,ii) / scale;  % Update process just if 
        drift(ii) = -psi_Z_OU_TS(-1i, dt * ii, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c);  % Compute drift term
    end

    %% Computing the log forward in the risk neutral measure

    logFwd = X + repmat([0, drift], Nsim, 1);  % Compute log forward prices

    %% Computing cumulants

    % Checking the last time step
    simCumulantsT = compCumulants(X(:, M+1)) * 1000;  
    theorCumulantsT = bctsCumulants(x0, alpha, beta_p, beta_n, c_p, c_n, gamma_c, T, b, flagOU_TS_FA) * 1000; 

    % Checking the first time step
    simCumulants_dt = compCumulants(X(:, 2)) * 1000;  
    theorCumulants_dt = bctsCumulants(x0, alpha, beta_p, beta_n, c_p, c_n, gamma_c, dt, b, flagOU_TS_FA) * 1000;  

end % function sim_OU_TS_FinAct_FGMC
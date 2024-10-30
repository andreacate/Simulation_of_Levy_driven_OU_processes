function [X, theorCumulantsT, simCumulantsT, theorCumulants_dt, simCumulants_dt, logFwd] = sim_TS_OU_FinVar_ED(x0, alpha, b, beta_p, ...
                                                               beta_n, c_p, c_n, gamma_c, T, M, Nsim)
% Simulation of TS-OU Finite Variation with Exact Decomposition following
% the algorithm 1 in [4] Sabino
%
% INPUT
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
% Nsim:     number of simulations
%
% OUTPUT
% X:                   matrix with simulations
% cumulants_process:   theoretical cumulants of the process
% cumulants_sim:       cumulants of the simulations
% logFwd:              log of fwd prices
%
% USES
% function simulateTS(alpha, beta, c)
% function theorCumulants(X0, alpha, beta, c, dt, b, flag)
% function bctsCumulants(X0, alpha, beta_p, beta_n, c_p, c_n, gamma_c, dt, b, flag)
% function compCumulants(vec)


    %% Quantities of interest

    flag_TS_OU = 2;
    X = zeros(Nsim, M+1);  % Initialize simulation matrix
    X(:, 1) = x0 * ones(Nsim, 1);  % Set initial condition
    drift = zeros(1, M);  % Initialize drift vector
    dt = T / M;  % Time step size

    % Lambda computation
    Gamma = gamma(1 - alpha);  
    a = exp(-b * dt); 
    lambda_p = c_p * Gamma * beta_p^(alpha) / alpha * (1 - a^(alpha));  
    lambda_n = c_n * Gamma * beta_n^(alpha) / alpha * (1 - a^(alpha));  

    %% Simulations

    % Poisson simulation
    poisson_p = random('Poisson', lambda_p, Nsim, M);  % Positive Poisson process
    poisson_n = random('Poisson', lambda_n, Nsim, M);  % Negative Poisson process

    % Process simulation
    for ii = 1:M
        % TS simulation
        X1_p = zeros(Nsim, 1);  % Positive increments
        X1_n = zeros(Nsim, 1);  % Negative increments

        [X1_p, ~] = simulateTS(alpha, beta_p, c_p * (1 - a^(alpha)),Nsim);  % Positive tempered stable increments
        [X1_n, ~] = simulateTS(alpha, beta_n, c_n * (1 - a^(alpha)),Nsim);  % Negative tempered stable increments


        % Compound Poisson simulation
        max_p = max(poisson_p(:, ii));  % Maximum positive jumps
        max_n = max(poisson_n(:, ii));  % Maximum negative jumps

        uniform_p = rand(Nsim, max_p);  % Uniform random variables for positive jumps
        uniform_n = rand(Nsim, max_n);  % Uniform random variables for negative jumps

        V_p = (1 + (a^(-alpha) - 1) * uniform_p).^(1 / alpha);  % Scaled positive jumps
        V_n = (1 + (a^(-alpha) - 1) * uniform_n).^(1 / alpha);  % Scaled negative jumps

        beta_hat_p = beta_p * V_p;  % Adjusted positive beta
        beta_hat_n = beta_n * V_n;  % Adjusted negative beta

        j_hat_p = random('Gamma', (1 - alpha) * ones(size(uniform_p)), 1 ./ beta_hat_p);  % Gamma distributed positive jumps
        j_hat_n = random('Gamma', (1 - alpha) * ones(size(uniform_n)), 1 ./ beta_hat_n);  % Gamma distributed negative jumps

        count_matrix_p = repmat(1:max_p, Nsim, 1);  % Count matrix for positive jumps
        count_matrix_n = repmat(1:max_n, Nsim, 1);  % Count matrix for negative jumps

        rep_matrix_p = repmat(poisson_p(:, ii), 1, max_p);  % Replicated positive Poisson counts
        rep_matrix_n = repmat(poisson_n(:, ii), 1, max_n);  % Replicated negative Poisson counts

        find_p = rep_matrix_p >= count_matrix_p;  % Indicator for positive jumps
        find_n = rep_matrix_n >= count_matrix_n;  % Indicator for negative jumps

        % Computation of the theoretical cumulants of the increment process to compensate it
        % Theoretical positive cumulants of the increment process
        cumulants_p = ctsCumulants(0, alpha, beta_p, c_p, dt, b, flag_TS_OU);
        % Theoretical negative cumulants of the increment process
        cumulants_n = ctsCumulants(0, alpha, beta_n, c_n, dt, b, flag_TS_OU);

        X2_p = sum(j_hat_p .* find_p, 2);  % Summed positive jumps
        X2_n = sum(j_hat_n .* find_n, 2);  % Summed negative jumps

        % Updating the process
        X(:, ii+1) = gamma_c * dt + a * X(:, ii) + X1_p + X2_p - cumulants_p(1) ...
                    - (X1_n + X2_n - cumulants_n(1));

        % Computing the drift term
        drift(ii) = -(psi_X_TS_OU(-1i, alpha, beta_p, beta_n, c_p, c_n, gamma_c) - ...
                      psi_X_TS_OU(-1i * exp(-b * dt * ii), alpha, beta_p, beta_n, c_p, c_n, gamma_c));
    end

    %% Computing cumulants

    % Checking the last time step
    simCumulantsT = compCumulants(X(:, M+1)) * 1000;  
    theorCumulantsT = bctsCumulants(x0, alpha, beta_p, beta_n, c_p, c_n, gamma_c, T, b, flag_TS_OU) * 1000;  

    % Checking the first time step
    simCumulants_dt = compCumulants(X(:, 2)) * 1000;  
    theorCumulants_dt = bctsCumulants(x0, alpha, beta_p, beta_n, c_p, c_n, gamma_c, dt, b, flag_TS_OU) * 1000;  % Theoretical cumulants for the first time step

    %% Computing the log forward in the risk neutral measure

    % Compute the log forward prices by adding the drift term
    logFwd = X + repmat([0, drift], Nsim, 1);

end % function sim_TS_OU_FinVar_ED